//! Bonded interactions: harmonic bonds, harmonic angles, and periodic proper/
//! improper dihedrals.
//!
//! M3 introduces bonded terms alongside the pair (non-bonded) model. To compose
//! them, the force-accumulation contract changes from *overwrite* (M7, single
//! model) to *accumulate*: each model adds into `system.forces` and the virial,
//! and the integrator zeroes both once before the composite force pass. See
//! [`crate::composite::CompositeForceModel`].
//!
//! Every term satisfies the conservation contract dU/dr = −F: bonds and angles
//! use closed-form gradients; the dihedral uses an exact analytic gradient of
//! the dihedral angle (derived so the four atomic gradients sum to zero —
//! translation invariance), verified by the finite-difference property test in
//! `tests/bonded.rs`.

use crate::error::Result;
use crate::force_model::{ForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::system::System;
use crate::vec3::Vec3;

/// One harmonic bond: atoms `i`, `j`, equilibrium length `r0`, stiffness `k`.
/// `U = ½ k (r − r0)²`.
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub r0: f64,
    pub k: f64,
}

/// One harmonic angle: atoms `i`–`j`–`k` (`j` is the vertex), equilibrium angle
/// `theta0` (radians), stiffness `ktheta`. `U = ½ k (θ − θ0)²`.
#[derive(Debug, Clone, Copy)]
pub struct Angle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub theta0: f64,
    pub ktheta: f64,
}

/// One periodic dihedral term: atoms `i`–`j`–`k`–`l` (the dihedral is the angle
/// between the i–j–k and j–k–l planes). The potential is a sum of Fourier terms
/// `U = Σ_n k_n (1 + cos(n·φ − δ_n))` (the OPLS/CHARMM proper-dihedral form;
/// impropers use the same form, typically a single n=2 term for planarity).
#[derive(Debug, Clone, Copy)]
pub struct DihedralTerm {
    pub n: i32,
    pub k: f64,
    pub delta: f64,
}

/// One dihedral: the four atoms plus its Fourier terms.
#[derive(Debug, Clone)]
pub struct Dihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub terms: Vec<DihedralTerm>,
}

/// The full bonded topology.
#[derive(Debug, Clone, Default)]
pub struct BondedTopology {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
}

impl BondedTopology {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty() && self.angles.is_empty() && self.dihedrals.is_empty()
    }
}

/// Force model evaluating all bonded terms. **Accumulates** into `system.forces`
/// and `virial` (does not zero them — the integrator / composite model does).
pub struct BondedForceModel {
    pub topology: BondedTopology,
}

impl BondedForceModel {
    pub fn new(topology: BondedTopology) -> Self {
        Self { topology }
    }
}

impl ForceModel for BondedForceModel {
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        _neighbors: &mut NeighborList,
    ) -> Result<f64> {
        let box_size = system.box_size;
        let mut pe = 0.0;
        let mut w = 0.0;

        for &b in &self.topology.bonds {
            let (u, fi, fj, vw) = bond_term(system.pos(b.i), system.pos(b.j), b.r0, b.k, box_size);
            system.add_force(b.i, fi);
            system.add_force(b.j, fj);
            pe += u;
            w += vw;
        }
        for &a in &self.topology.angles {
            let (u, fi, fj, fk, vw) = angle_term(
                system.pos(a.i),
                system.pos(a.j),
                system.pos(a.k),
                a.theta0,
                a.ktheta,
            );
            system.add_force(a.i, fi);
            system.add_force(a.j, fj);
            system.add_force(a.k, fk);
            pe += u;
            w += vw;
        }
        for d in &self.topology.dihedrals {
            let (u, fi, fj, fk, fl, vw) = dihedral_term(
                system.pos(d.i),
                system.pos(d.j),
                system.pos(d.k),
                system.pos(d.l),
                &d.terms,
                box_size,
            );
            system.add_force(d.i, fi);
            system.add_force(d.j, fj);
            system.add_force(d.k, fk);
            system.add_force(d.l, fl);
            pe += u;
            w += vw;
        }

        virial.w += w;
        Ok(pe)
    }
}

// ---------------------------------------------------------------------------
// Bond (harmonic). U = ½ k (r − r0)², F_i = −k(r−r0) r̂, F_j = −F_i.
// ---------------------------------------------------------------------------

#[inline]
fn bond_term(
    ri: Vec3,
    rj: Vec3,
    r0: f64,
    k: f64,
    box_size: crate::system::SimBox,
) -> (f64, Vec3, Vec3, f64) {
    let rvec = box_size.min_image(ri - rj); // j → i
    let r = rvec.norm();
    let dr = r - r0;
    let u = 0.5 * k * dr * dr;
    // F_i = -dU/dr * (rvec/r) = -k*dr * rvec/r ; F_j = -F_i.
    let mag = -k * dr / r;
    let fi = rvec * mag;
    let fj = fi * -1.0;
    // Virial contribution: r_ij · F_i.
    let vw = rvec.dot(fi);
    (u, fi, fj, vw)
}

// ---------------------------------------------------------------------------
// Angle (harmonic). U = ½ k (θ − θ0)², vertex j in the middle.
// r1 = r_i − r_j (j→i), r2 = r_k − r_j (j→k), cosθ = (r1·r2)/(|r1||r2|).
// ---------------------------------------------------------------------------

#[inline]
fn angle_term(ri: Vec3, rj: Vec3, rk: Vec3, theta0: f64, k: f64) -> (f64, Vec3, Vec3, Vec3, f64) {
    let r1 = ri - rj; // j → i
    let r2 = rk - rj; // j → k
    let a = r1.norm();
    let b = r2.norm();
    let c = (r1.dot(r2) / (a * b)).clamp(-1.0, 1.0);
    let s = (1.0 - c * c).sqrt();
    let theta = c.acos();
    let dtheta = theta - theta0;
    let u = 0.5 * k * dtheta * dtheta;
    let du_dtheta = k * dtheta;

    // Exact gradients ∂θ/∂r (translation-invariant: g_i + g_j + g_k = 0).
    // g_i = -(1/s)[ r2/(ab) − c·r1/a² ],  g_k = -(1/s)[ r1/(ab) − c·r2/b² ].
    let inv_s = if s > 1e-12 { 1.0 / s } else { 0.0 };
    let gi = (r2 * (1.0 / (a * b)) - r1 * (c / (a * a))) * (-inv_s);
    let gk = (r1 * (1.0 / (a * b)) - r2 * (c / (b * b))) * (-inv_s);
    let gj = (gi + gk) * -1.0;

    // F = -dU/dθ · g.
    let fi = gi * (-du_dtheta);
    let fj = gj * (-du_dtheta);
    let fk = gk * (-du_dtheta);
    // Virial: Σ r_{a,vertex}·F_a — use the two bond vectors at the vertex.
    let vw = r1.dot(fi) + r2.dot(fk);
    (u, fi, fj, fk, vw)
}

// ---------------------------------------------------------------------------
// Dihedral (periodic). Atoms i–j–k–l; φ is the angle between planes (i,j,k)
// and (j,k,l). U = Σ_n k_n (1 + cos(n φ − δ_n)).
//
// Conventions (GROMACS): r_ij = r_j − r_i, r_kj = r_j − r_k, r_kl = r_l − r_k,
// m = r_ij × r_kj, n = r_kj × r_kl. cosφ = (m·n)/(|m||n|); the sign of sinφ
// comes from the triple product with r_kj so φ ∈ (−π, π].
//
// The gradient is the exact analytic dφ (derived by hand; the four atomic
// gradients sum to zero). F_a = -(dU/dφ)(∂φ/∂r_a).
// ---------------------------------------------------------------------------

#[inline]
fn dihedral_term(
    ri: Vec3,
    rj: Vec3,
    rk: Vec3,
    rl: Vec3,
    terms: &[DihedralTerm],
    box_size: crate::system::SimBox,
) -> (f64, Vec3, Vec3, Vec3, Vec3, f64) {
    // Bond vectors through the minimum image (the four atoms are covalently
    // linked, so across a periodic boundary the bonded vector is the min image).
    let r_ij = box_size.min_image(rj - ri);
    let r_kj = box_size.min_image(rj - rk);
    let r_kl = box_size.min_image(rl - rk);

    let m = r_ij.cross(r_kj);
    let n = r_kj.cross(r_kl);
    let mm = m.norm();
    let nn = n.norm();
    // Degenerate (collinear) → no torque, zero contribution.
    if mm < 1e-12 || nn < 1e-12 {
        return (
            0.0,
            Vec3::zero(),
            Vec3::zero(),
            Vec3::zero(),
            Vec3::zero(),
            0.0,
        );
    }
    let r = r_kj.norm();
    let rm2 = r * r;
    let cosphi = (m.dot(n) / (mm * nn)).clamp(-1.0, 1.0);
    let sinphi = (m.cross(n)).dot(r_kj) / (mm * nn * r);
    let phi = sinphi.atan2(cosphi);

    // Convention: the raw cross-product angle (`phi` above) is offset by π from
    // the standard IUPAC/GROMACS dihedral (it gives 0 for trans, π for cis).
    // Shift to the standard convention (trans = π, cis = 0). The gradient of φ
    // is unchanged by adding a constant, so the analytic gradient below (derived
    // for `phi`) is also the gradient of `phi_std`.
    let phi_std = phi + std::f64::consts::PI;

    // Energy and dU/dφ (standard convention).
    let mut u = 0.0;
    let mut du_dphi = 0.0;
    for t in terms {
        let nf = t.n as f64;
        u += t.k * (1.0 + (nf * phi_std - t.delta).cos());
        du_dphi += -t.k * nf * (nf * phi_std - t.delta).sin();
    }

    // --- exact gradient of φ w.r.t. the bond vectors ---
    let cp = cosphi;
    let sp = sinphi;
    let mnr = mm * nn * r;
    let coef_s = cp / mnr; // multiplies ∂S
    let coef_c = sp / (mm * nn); // multiplies ∂C
    let coef_r = cp * sp / rm2; // multiplies r_kj (the dR term)

    // ∂C/∂r_ij = r_kj × n ; ∂S/∂r_ij = R² n − (r_kj·n) r_kj.
    let g_ij = (n * rm2 - r_kj * n.dot(r_kj)) * coef_s - r_kj.cross(n) * coef_c;
    // ∂C/∂r_kl = m × r_kj ; ∂S/∂r_kl = R² m − (m·r_kj) r_kj.
    let g_kl = (m * rm2 - r_kj * m.dot(r_kj)) * coef_s - m.cross(r_kj) * coef_c;
    // ∂C/∂r_kj = n×r_ij + r_kl×m
    let dc_kj = n.cross(r_ij) + r_kl.cross(m);
    // ∂S/∂r_kj = ((r_ij·n)+(m·r_kl)) r_kj − (r_ij·r_kj) n − (r_kl·r_kj) m + (m×n)
    let ds_kj =
        r_kj * (r_ij.dot(n) + m.dot(r_kl)) - n * r_ij.dot(r_kj) - m * r_kl.dot(r_kj) + m.cross(n);
    let g_kj = ds_kj * coef_s - dc_kj * coef_c - r_kj * coef_r;

    // Atomic gradients of φ.
    let g_i = g_ij * -1.0;
    let g_l = g_kl;
    let g_j = g_ij + g_kj;
    let g_k = (g_kj + g_kl) * -1.0;

    // F = -dU/dφ · g.
    let fi = g_i * (-du_dphi);
    let fj = g_j * (-du_dphi);
    let fk = g_k * (-du_dphi);
    let fl = g_l * (-du_dphi);
    // Virial: use the three bonded vectors' dot with the end-atom forces.
    let vw = r_ij.dot(fi) + r_kl.dot(fl) + r_kj.dot(fj);
    (u, fi, fj, fk, fl, vw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::SimBox;

    fn box10() -> SimBox {
        SimBox::cubic(100.0)
    }

    #[test]
    fn bond_force_signs_and_energy() {
        // r = 1.5, r0 = 1.0, k = 4 → stretched → attractive (F_i toward j → -x).
        let ri = Vec3::new(1.5, 0.0, 0.0);
        let rj = Vec3::new(0.0, 0.0, 0.0);
        let (u, fi, fj, _) = bond_term(ri, rj, 1.0, 4.0, box10());
        assert!((u - 0.5 * 4.0 * 0.25).abs() < 1e-12);
        assert!(
            fi.0[0] < 0.0,
            "stretched bond: F_i should pull toward j (-x)"
        );
        assert!((fi + fj).norm() < 1e-12, "Newton's 3rd");
    }

    #[test]
    fn angle_force_sums_to_zero() {
        let (u, fi, fj, fk, _) = angle_term(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::zero(),
            Vec3::new(0.0, 1.0, 0.0),
            1.7,
            5.0,
        );
        assert!(u >= 0.0);
        assert!(
            (fi + fj + fk).norm() < 1e-9,
            "angle forces must sum to zero"
        );
    }

    #[test]
    fn dihedral_forces_sum_to_zero() {
        // A non-planar quadrilateral.
        let ri = Vec3::new(1.0, 0.0, 0.0);
        let rj = Vec3::new(0.0, 0.0, 0.0);
        let rk = Vec3::new(0.0, 1.0, 0.0);
        let rl = Vec3::new(0.0, 1.0, 1.0);
        let terms = vec![DihedralTerm {
            n: 1,
            k: 2.0,
            delta: 0.0,
        }];
        let (u, fi, fj, fk, fl, _) = dihedral_term(ri, rj, rk, rl, &terms, box10());
        assert!(u.is_finite());
        assert!(
            (fi + fj + fk + fl).norm() < 1e-9,
            "dihedral forces must sum to zero"
        );
    }
}
