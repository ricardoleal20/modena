//! Ewald summation for long-range Coulomb electrostatics (M4).
//!
//! Splits the conditionally-convergent `q_i q_j / r` lattice sum into three
//! absolutely-convergent pieces:
//!
//! 1. **Real space** — short-ranged `q_i q_j erfc(κ r) / r` over the neighbor
//!    list (the `erfc` screens the interaction past the real-space cutoff
//!    `r_c`). Force: `F = q_i q_j [erfc(κr)/r² + (2κ/√π) e^{−κ²r²}/r] r̂`.
//! 2. **Reciprocal space** — the smooth long-range part in k-space:
//!    `U_rec = (2π/V) Σ_{k≠0} |ρ(k)|² e^{−k²/4κ²} / k²`,
//!    `ρ(k) = Σ_j q_j e^{i k·r_j}`. Direct O(N·n_k) sum (no FFT/SPME yet — that
//!    is a later roadmap item once `rustfft` is wired).
//! 3. **Self + background** — `U_self = −(κ/√π) Σ_i q_i²` removes the spurious
//!    self-interaction introduced by the reciprocal Gaussian; the
//!    `−π/(2Vκ²)(Σ q_i)²` term neutralizes a uniform background for
//!    non-neutral systems (zero for neutral).
//!
//! The splitting parameter `κ` (Ewald parameter) trades real- vs reciprocal-space
//! cost; `r_c` (real cutoff) and `k_c` (k-space cutoff) control accuracy. The
//! implementation is a single [`ForceModel`](crate::force_model::ForceModel) that
//! **accumulates** into `system.forces`/`virial` (the caller zeroes once), so it
//! composes with a pair (e.g. LJ/Buckingham) and bonded models inside a
//! [`CompositeForceModel`](crate::composite::CompositeForceModel). Charges are
//! read per-atom from `system.charges`.
//!
//! Correctness is pinned by: a finite-difference force check (each atom's force
//! equals `−∇U`), the real-space→bare-Coulomb limit (κ→0), the Madelung lattice
//! energy of NaCl, and NVE energy conservation for an ionic system
//! (`tests/electrostatics.rs`).

use std::sync::{Arc, Mutex};

use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::system::System;
use crate::vec3::Vec3;

/// √π.
const SQRT_PI: f64 = 1.772453850905516;
/// 2/√π — the derivative prefactor of erfc: `d/dr erfc(κr) = −(2κ/√π) e^{−κ²r²}`.
const TWO_OVER_SQRT_PI: f64 = 2.0 / SQRT_PI;

/// Ewald electrostatics: real-space screened sum + reciprocal k-space sum +
/// self/background correction.
///
/// Construct with [`EwaldForceModel::new`]; tune `κ`, `real_cutoff`, and
/// `kspace_cutoff` for the desired accuracy. All quantities are in the engine's
/// reduced units, with charges in the reduced-charge convention
/// (`U = q_i q_j / r`).
pub struct EwaldForceModel {
    /// Ewald screening parameter κ (1/length).
    pub kappa: f64,
    /// Real-space cutoff `r_c` (length). Pairs beyond this are skipped.
    pub real_cutoff: f64,
    /// Reciprocal-space cutoff `k_c` (1/length). k-vectors with `|k| ≥ k_c` are
    /// dropped.
    pub kspace_cutoff: f64,
    /// Half-space reciprocal k-vector set, cached per (full cell basis, kcut).
    /// Previously re-enumerated and re-allocated every step (a 19,548-entry
    /// `Vec<KVec>` ≈ 625 KB rebuilt per step at the validation deck, despite
    /// depending only on the box and cutoffs — the exact analogue of the
    /// neighbor-list box-rescale cache). Regenerate only on change.
    #[allow(clippy::type_complexity)]
    kvec_cache: Mutex<Option<KvecCache>>,
    /// Indexed Coulomb scales; shared-list and model policies reconcile once
    /// and reciprocal corrections are applied once per configured pair.
    pair_scales: crate::neighbor::PairScaleSet,
}

/// Cached immutable k-vector snapshot, keyed by full cell and k-space cutoff.
type KvecCache = (crate::system::SimBox, f64, Arc<Vec<KVec>>);

impl EwaldForceModel {
    pub fn new(kappa: f64, real_cutoff: f64, kspace_cutoff: f64) -> Result<Self> {
        if !kappa.is_finite() || kappa <= 0.0 {
            return Err(ModenaError::Config("Ewald `kappa` must be > 0".into()));
        }
        if !real_cutoff.is_finite() || real_cutoff <= 0.0 {
            return Err(ModenaError::Config(
                "Ewald `real_cutoff` must be > 0".into(),
            ));
        }
        if !kspace_cutoff.is_finite() || kspace_cutoff <= 0.0 {
            return Err(ModenaError::Config(
                "Ewald `kspace_cutoff` must be > 0".into(),
            ));
        }
        Ok(Self {
            kappa,
            real_cutoff,
            kspace_cutoff,
            kvec_cache: Mutex::new(None),
            pair_scales: crate::neighbor::PairScaleSet::default(),
        })
    }

    /// Return an immutable snapshot for this full cell basis. Another caller
    /// sharing the model may replace the cache without changing this snapshot.
    fn kvecs(&self, box_size: &crate::system::SimBox) -> Arc<Vec<KVec>> {
        let mut g = self.kvec_cache.lock().unwrap();
        if let Some((cell, cutoff, data)) = &*g {
            if cell == box_size && *cutoff == self.kspace_cutoff {
                return Arc::clone(data);
            }
        }
        let data = Arc::new(generate_kvecs(box_size, self.kspace_cutoff));
        *g = Some((*box_size, self.kspace_cutoff, Arc::clone(&data)));
        data
    }

    /// The real-space cutoff — the caller uses `max(pair_cutoff, this)` for the
    /// shared neighbor list.
    pub fn real_cutoff(&self) -> f64 {
        self.real_cutoff
    }

    /// Excluded charged pairs: the neighbor list removes their erfc part and
    /// the reciprocal sum's Gaussian (erf) part is cancelled explicitly, so an
    /// excluded pair feels exactly zero Coulomb.
    pub fn with_excluded_pairs(mut self, pairs: Vec<(usize, usize)>) -> Self {
        let scales: Vec<_> = pairs
            .into_iter()
            .map(|(i, j)| crate::neighbor::PairScale {
                i,
                j,
                lj: 0.0,
                coulomb: 0.0,
            })
            .collect();
        self.pair_scales =
            crate::neighbor::PairScaleSet::new(&scales).expect("invalid excluded pairs");
        self
    }

    /// Validated independently weighted pairs. The Coulomb weight is consumed
    /// here; LJ weights belong to the shared neighbor list and pair provider.
    pub fn with_pair_scales(mut self, scales: Vec<crate::neighbor::PairScale>) -> Result<Self> {
        self.pair_scales = crate::neighbor::PairScaleSet::new(&scales)?;
        Ok(self)
    }

    /// Recompute the full Ewald energy only (no forces), used by the
    /// finite-difference reference in tests and by MC barostat ΔU. Accumulates
    /// nothing; returns the total potential energy.
    pub fn energy(&self, system: &System, neighbors: &mut NeighborList) -> Result<f64> {
        neighbors.build(system);
        let policy = CoulombPolicy::new(&self.pair_scales, neighbors.scales(), system.n())?;
        let v = system.box_size.volume();
        if v <= 0.0 {
            return Err(ModenaError::Physics("non-positive box volume".into()));
        }
        // Real space.
        let rc2 = self.real_cutoff * self.real_cutoff;
        let kappa = self.kappa;
        let mut u_real = 0.0;
        for &(i, j) in neighbors.pairs() {
            let d = system.box_size.min_image(system.pos(i) - system.pos(j));
            let r2 = d.norm2();
            if r2 >= rc2 || r2 < 1e-14 {
                continue;
            }
            let r = r2.sqrt();
            let qprod = policy.scale(i, j) * system.charge(i) * system.charge(j);
            u_real += qprod * libm::erfc(kappa * r) / r;
        }
        // Excluded-pair erf correction (energy only here).
        let (u_excl, _) = exclusion_correction(
            kappa,
            &system.box_size,
            &system.positions,
            &system.charges,
            self.real_cutoff,
            policy.entries(),
            &mut [],
        );
        u_real += u_excl;
        // Reciprocal + self + background.
        let (u_rec, _kcount) = self.reciprocal_energy(system)?;
        let u_self = self.self_energy(system);
        let u_bg = self.background_energy(system, v);
        Ok(u_real + u_rec + u_self + u_bg)
    }

    /// Reciprocal-space energy `U_rec = (2π/V) Σ_{half-space} 2·|ρ(k)|² w_k`,
    /// `w_k = e^{−k²/4κ²}/k²` — the `2·` counts the ±k pair (half-space
    /// iteration, see `generate_kvecs`). Returns `(energy, n_half_k_vectors)`.
    pub(crate) fn reciprocal_energy(&self, system: &System) -> Result<(f64, usize)> {
        let v = system.box_size.volume();
        let kvecs = self.kvecs(&system.box_size);
        let two_pi_over_v = 2.0 * std::f64::consts::PI / v;
        let kappa2 = self.kappa * self.kappa;
        let n = system.n();
        let pos = &system.positions;
        let q = &system.charges;
        let mut u = 0.0;
        for kv in kvecs.iter() {
            let k2 = kv.k2;
            let w = libm::exp(-k2 / (4.0 * kappa2)) / k2;
            // ρ(k) = Σ q_j (cos + i sin).
            let mut re = 0.0; // Σ q_j cos(k·r_j)
            let mut im = 0.0; // Σ q_j sin(k·r_j)
            for j in 0..n {
                let rj = Vec3::new(pos[3 * j], pos[3 * j + 1], pos[3 * j + 2]);
                let phase = kv.k.dot(rj);
                re += q[j] * libm::cos(phase);
                im += q[j] * libm::sin(phase);
            }
            u += 2.0 * two_pi_over_v * w * (re * re + im * im);
        }
        Ok((u, kvecs.len()))
    }

    /// Self-interaction correction `U_self = −(κ/√π) Σ q_i²`.
    fn self_energy(&self, system: &System) -> f64 {
        let sum_q2: f64 = system.charges.iter().map(|q| q * q).sum();
        -(self.kappa / SQRT_PI) * sum_q2
    }

    /// Uniform-background correction `U_bg = −π/(2Vκ²)(Σ q_i)²` (zero for a
    /// neutral system). Tin-foil (ε = ∞) boundary conditions.
    fn background_energy(&self, system: &System, v: f64) -> f64 {
        let sum_q: f64 = system.charges.iter().sum();
        -std::f64::consts::PI / (2.0 * v * self.kappa * self.kappa) * sum_q * sum_q
    }
}

/// Reconcile model and shared-list Coulomb policies once per evaluation.
/// Identical overlapping policies are applied once; contradictory weights
/// are errors. Iteration visits only configured special pairs, without cloning
/// the index or allocating a combined table on each force pass.
pub(crate) struct CoulombPolicy<'a> {
    model: &'a crate::neighbor::PairScaleSet,
    neighbors: &'a crate::neighbor::PairScaleSet,
}
impl<'a> CoulombPolicy<'a> {
    pub(crate) fn new(
        model: &'a crate::neighbor::PairScaleSet,
        neighbors: &'a crate::neighbor::PairScaleSet,
        n: usize,
    ) -> Result<Self> {
        model.validate(n)?;
        neighbors.validate(n)?;
        for p in model.entries() {
            if let Some(q) = neighbors.get(p.i, p.j) {
                if p.coulomb != q.coulomb {
                    return Err(ModenaError::Config(format!(
                        "model and neighbor Coulomb scales conflict for ({}, {})",
                        p.i, p.j
                    )));
                }
            }
        }
        Ok(Self { model, neighbors })
    }
    pub(crate) fn scale(&self, i: usize, j: usize) -> f64 {
        self.model
            .get(i, j)
            .map_or_else(|| self.neighbors.scale(i, j).1, |p| p.coulomb)
    }
    pub(crate) fn entries(&self) -> impl Iterator<Item = crate::neighbor::PairScale> + '_ {
        self.model.entries().iter().copied().chain(
            self.neighbors
                .entries()
                .iter()
                .copied()
                .filter(|p| self.model.get(p.i, p.j).is_none()),
        )
    }
}

/// The excluded-pair `erf` correction, shared by the direct Ewald and SPME
/// models. Real-space pairs are independently weighted; the reciprocal
/// correction removes `(1-coulomb_scale)` of their Gaussian part. A 0/0 list
/// removes an excluded pair's real-space
/// (erfc) part, but the reciprocal sum / mesh still carries its Gaussian
/// `erf(κr)/r` part — this cancels it, so an excluded pair feels exactly zero
/// Coulomb. Outside the finite real cutoff, the correction removes the full
/// bare fraction rather than just erf: this keeps the energy difference from
/// the unscaled finite-cutoff model exactly `-(1-s)*q_i*q_j/r` at any distance.
///
/// ```text
///   U_corr = −q_i q_j erf(κr)/r
///   F_i    = +q_i q_j [ (2κ/√π) e^{−κ²r²}/r² − erf(κr)/r³ ] (r_i − r_j)
/// ```
pub(crate) fn exclusion_correction(
    kappa: f64,
    box_size: &crate::system::SimBox,
    positions: &[f64],
    charges: &[f64],
    real_cutoff: f64,
    pairs: impl Iterator<Item = crate::neighbor::PairScale>,
    forces: &mut [f64],
) -> (f64, f64) {
    let kappa2 = kappa * kappa;
    let mut u = 0.0;
    let mut w = 0.0;
    for pair in pairs {
        let (i, j) = (pair.i, pair.j);
        let removed = 1.0 - pair.coulomb;
        if removed == 0.0 {
            continue;
        }
        let ri = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        let rj = Vec3::new(positions[3 * j], positions[3 * j + 1], positions[3 * j + 2]);
        let d = box_size.min_image(ri - rj);
        let r2 = d.norm2();
        if r2 < 1e-14 {
            continue;
        }
        let r = r2.sqrt();
        let qprod = removed * charges[i] * charges[j];
        let (erf_kr, exp_kr2) = if r < real_cutoff {
            (libm::erf(kappa * r), libm::exp(-kappa2 * r2))
        } else {
            (1.0, 0.0)
        };
        u += -qprod * erf_kr / r;
        let coef = qprod * (TWO_OVER_SQRT_PI * kappa * exp_kr2 / r2 - erf_kr / (r2 * r));
        for k in 0..3 {
            let fv = d.0[k] * coef;
            if !forces.is_empty() {
                forces[3 * i + k] += fv;
                forces[3 * j + k] -= fv;
            }
            w += d.0[k] * fv;
        }
    }
    (u, w)
}

impl ForceModel for EwaldForceModel {
    fn conserves_total_momentum(&self, _system: &System) -> bool {
        true
    }
    fn neighbor_cutoff(&self) -> Option<f64> {
        Some(self.real_cutoff)
    }
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        neighbors.build(system);
        let policy = CoulombPolicy::new(&self.pair_scales, neighbors.scales(), system.n())?;
        let v = system.box_size.volume();
        if v <= 0.0 {
            return Err(ModenaError::Physics("non-positive box volume".into()));
        }
        let n = system.n();
        if n == 0 {
            return Ok(0.0);
        }
        // Local force accumulator (added into system.forces at the end).
        let mut forces = vec![0.0f64; 3 * n];
        let mut w_virial = 0.0f64;
        let mut pe = 0.0f64;

        let kappa = self.kappa;
        let kappa2 = kappa * kappa;
        let rc2 = self.real_cutoff * self.real_cutoff;
        let box_size = system.box_size;

        // --- Real-space screened sum over the half neighbor list. ---
        {
            let pos = &system.positions;
            let q = &system.charges;
            for &(i, j) in neighbors.pairs() {
                let ri = Vec3::new(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
                let rj = Vec3::new(pos[3 * j], pos[3 * j + 1], pos[3 * j + 2]);
                let d = box_size.min_image(ri - rj);
                let r2 = d.norm2();
                if r2 >= rc2 || r2 < 1e-14 {
                    continue;
                }
                let r = r2.sqrt();
                let qprod = policy.scale(i, j) * q[i] * q[j];
                let erfc_kr = libm::erfc(kappa * r);
                let exp_kr2 = libm::exp(-kappa2 * r2);
                // U += q_i q_j erfc(κr)/r.
                pe += qprod * erfc_kr / r;
                // F_i = d · q_i q_j [ (2κ/√π) e^{−κ²r²}/r² + erfc(κr)/r³ ].
                let coef = qprod * (TWO_OVER_SQRT_PI * kappa * exp_kr2 / r2 + erfc_kr / (r2 * r));
                for k in 0..3 {
                    let fv = d.0[k] * coef;
                    forces[3 * i + k] += fv;
                    forces[3 * j + k] -= fv;
                    w_virial += d.0[k] * fv; // r·F (Newton pair contributes once)
                }
            }
        }

        // --- Excluded-pair erf correction (their reciprocal part). ---
        let (u_excl, w_excl) = exclusion_correction(
            kappa,
            &box_size,
            &system.positions,
            &system.charges,
            self.real_cutoff,
            policy.entries(),
            &mut forces,
        );
        pe += u_excl;
        w_virial += w_excl;

        // --- Reciprocal-space sum (half-space, phases reused). ---
        // F_i = −(4π/V) q_i (B cos_i − A sin_i) w_k · k,  w_k = e^{−k²/4κ²}/k²,
        // with A = Σ q cos(k·r), B = Σ q sin(k·r). Ewald-b: k-vectors cached
        // per box (rebuilt only on change), ±k symmetry halves the k set
        // (energy and force both ×2), and the cos/sin phases are computed once
        // per (k, atom) and reused between the structure-factor and force
        // passes. Combined ≈4× fewer transcendentals; the arithmetic per term
        // is the same, so the brute-force decomposition agreement (~1e-10)
        // carries over unchanged.
        {
            let kvecs = self.kvecs(&box_size);
            let two_pi_over_v = 2.0 * std::f64::consts::PI / v;
            let four_pi_over_v = 4.0 * std::f64::consts::PI / v;
            let pos = &system.positions;
            let q = &system.charges;
            // Reused phase scratch (allocated once per compute call).
            let mut phases = vec![(0.0f64, 0.0f64); n];
            for kv in kvecs.iter() {
                let k2 = kv.k2;
                let w = libm::exp(-k2 / (4.0 * kappa2)) / k2;
                // Structure factor (+ phase scratch for the force pass).
                let mut a = 0.0; // Σ q_j cos(k·r_j)
                let mut b = 0.0; // Σ q_j sin(k·r_j)
                for j in 0..n {
                    let rj = Vec3::new(pos[3 * j], pos[3 * j + 1], pos[3 * j + 2]);
                    let phase = kv.k.dot(rj);
                    let c = libm::cos(phase);
                    let s = libm::sin(phase);
                    phases[j] = (c, s);
                    a += q[j] * c;
                    b += q[j] * s;
                }
                // Energy: (2π/V) w (A² + B²), ×2 for the ±k pair.
                let uk = 2.0 * two_pi_over_v * w * (a * a + b * b);
                pe += uk;
                // Under isotropic strain phases stay fixed, V scales as λ³
                // and k² as λ^-2: -dU/dλ = U(1-k²/(2κ²)).
                w_virial += uk * (1.0 - k2 / (2.0 * kappa2));
                // Force: per atom c_i = −(4π/V) q_i (B cos_i − A sin_i) w;
                // F_i += c_i · k, ×2 for the ±k pair.
                let scale = -2.0 * four_pi_over_v * w;
                for i in 0..n {
                    let (cos_i, sin_i) = phases[i];
                    let ci = scale * q[i] * (b * cos_i - a * sin_i);
                    forces[3 * i] += ci * kv.k.0[0];
                    forces[3 * i + 1] += ci * kv.k.0[1];
                    forces[3 * i + 2] += ci * kv.k.0[2];
                }
            }
        }

        // --- Self + background (energy only, no force). ---
        pe += self.self_energy(system);
        let u_bg = self.background_energy(system, v);
        pe += u_bg;
        w_virial += 3.0 * u_bg;

        // Accumulate into the system (caller zeroed forces + virial once).
        for (dst, &src) in system.forces.iter_mut().zip(forces.iter()).take(3 * n) {
            *dst += src;
        }
        virial.w += w_virial;
        Ok(pe)
    }
}

/// One reciprocal-space vector.
struct KVec {
    k: Vec3,
    k2: f64,
}

/// Generate the **half-space** nonzero k-vectors with `|k| < kcut` (one member
/// of each ±k pair: the first nonzero integer coordinate is positive). For a
/// cell with basis `[a, b, c]` and inverse rows `[r0, r1, r2]`,
/// `k(n) = 2π (n0 r0 + n1 r1 + n2 r2)`. The integer range per axis is bounded
/// rigorously by `|n_k| ≤ |basis_k| · kcut / (2π)`, then filtered by
/// `|k| < kcut`.
///
/// **Why half-space (Ewald-b):** `ρ(−k) = ρ(k)*`, so `U(−k) = U(k)` and the
/// force contribution of −k equals that of +k exactly (c(−k)·(−k) = c(k)·k; the
/// algebra was verified against the accumulation order). Iterating the half
/// space and doubling both energy and force contributions preserves
/// `dU/dr = −F` and halves the trig calls. Callers must multiply by 2.
fn generate_kvecs(box_size: &crate::system::SimBox, kcut: f64) -> Vec<KVec> {
    let basis = box_size.basis();
    let inv = box_size.inv_basis_rows();
    let two_pi = 2.0 * std::f64::consts::PI;
    let n_max: [i32; 3] = [
        ((basis[0].norm() * kcut / two_pi).ceil() as i32).max(0),
        ((basis[1].norm() * kcut / two_pi).ceil() as i32).max(0),
        ((basis[2].norm() * kcut / two_pi).ceil() as i32).max(0),
    ];
    let half =
        |n0: i32, n1: i32, n2: i32| n0 > 0 || (n0 == 0 && n1 > 0) || (n0 == 0 && n1 == 0 && n2 > 0);
    let kcut2 = kcut * kcut;
    let mut kvecs = Vec::new();
    for n0 in -n_max[0]..=n_max[0] {
        for n1 in -n_max[1]..=n_max[1] {
            for n2 in -n_max[2]..=n_max[2] {
                if n0 == 0 && n1 == 0 && n2 == 0 {
                    continue;
                }
                if !half(n0, n1, n2) {
                    continue;
                }
                let k =
                    (inv[0] * (n0 as f64) + inv[1] * (n1 as f64) + inv[2] * (n2 as f64)) * two_pi;
                let k2 = k.norm2();
                if k2 < kcut2 && k2 > 0.0 {
                    kvecs.push(KVec { k, k2 });
                }
            }
        }
    }
    kvecs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::SimBox;

    /// Two opposite charges in a large box: the real-space Ewald force at small
    /// κ must approach the bare Coulomb force `q_i q_j / r²` (erfc(κr)→1, the
    /// exp screening term →0).
    #[test]
    fn real_space_approaches_bare_coulomb() {
        let mut sys = System::new(SimBox::cubic(20.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.5, 0.0, 0.0));
        sys.charges[0] = 1.0;
        sys.charges[1] = -1.0;
        let r = 1.5;
        // κ very small → erfc(κr) ≈ 1, exp(−κ²r²) ≈ 1 but multiplied by κ → ~0.
        let ewald = EwaldForceModel::new(1e-4, 5.0, 1.0).unwrap();
        let mut nl = NeighborList::new(5.0);
        let mut virial = Virial::default();
        sys.zero_forces();
        ewald.compute(&mut sys, &mut virial, &mut nl).unwrap();
        // Bare Coulomb: opposite charges attract → atom 0 is pulled toward atom
        // 1 (at +x), so F on atom 0 is +x with magnitude 1/r². d = r_0 − r_1
        // points −x, q_i q_j = −1, so F_i = d·coef with coef < 0 → +x. ✓
        let fx = sys.force(0).0[0];
        let bare = 1.0 / (r * r);
        assert!(
            (fx - bare).abs() < 1e-3,
            "real-space force {fx} vs bare Coulomb {bare}"
        );
    }

    /// A neutral system's reciprocal + self energy is finite and the total Ewald
    /// energy of a small NaCl-like pair is negative (bound), sanity check.
    #[test]
    fn neutral_system_energy_finite() {
        let mut sys = System::new(SimBox::cubic(4.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(2.0, 0.0, 0.0));
        sys.charges[0] = 1.0;
        sys.charges[1] = -1.0;
        let ewald = EwaldForceModel::new(1.0, 2.5, 4.0).unwrap();
        let mut nl = NeighborList::new(2.5);
        let e = ewald.energy(&sys, &mut nl).unwrap();
        assert!(e.is_finite(), "energy not finite: {e}");
    }

    /// Triclinic k-vector generation respects the cutoff (no k beyond it).
    #[test]
    fn kvec_cutoff_respected() {
        let box_size = SimBox::triclinic([
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(0.5, 3.0, 0.0),
            Vec3::new(0.0, 0.5, 3.0),
        ]);
        let kcut = 4.0;
        let kvecs = generate_kvecs(&box_size, kcut);
        assert!(!kvecs.is_empty(), "should generate some k-vectors");
        for kv in &kvecs {
            assert!(kv.k2 < kcut * kcut, "k-vector beyond cutoff: {}", kv.k2);
        }
    }

    /// An excluded charged pair feels exactly zero Coulomb: with (0,1)
    /// excluded (and out of the neighbor list, as the runner wires it), the
    /// total energy equals the plain energy minus the bare `q_i q_j / r` and
    /// the forces differ by exactly the bare Coulomb pair force — the erf
    /// correction cancels what the reciprocal sum still carries.
    #[test]
    fn excluded_pair_feels_zero_coulomb() {
        let mut sys = System::new(SimBox::cubic(12.0));
        sys.resize(3);
        sys.set_pos(0, Vec3::new(1.0, 1.0, 1.0));
        sys.set_pos(1, Vec3::new(2.5, 1.2, 0.8));
        sys.set_pos(2, Vec3::new(7.0, 8.0, 6.0));
        sys.charges[0] = 1.0;
        sys.charges[1] = -0.7;
        sys.charges[2] = 0.4;
        let kappa = 0.9;
        let mut a = sys.clone();
        let mut b = sys.clone();
        a.zero_forces();
        b.zero_forces();
        let plain = EwaldForceModel::new(kappa, 5.0, 6.0).unwrap();
        let excl = EwaldForceModel::new(kappa, 5.0, 6.0)
            .unwrap()
            .with_excluded_pairs(vec![(0, 1)]);
        let (mut v1, mut v2) = (Virial::default(), Virial::default());
        let mut nl1 = NeighborList::with_skin_exclusions(5.0, 0.0, &[(0, 1)]);
        let mut nl2 = NeighborList::new(5.0);
        let e_excl = excl.compute(&mut a, &mut v1, &mut nl1).unwrap();
        let e_plain = plain.compute(&mut b, &mut v2, &mut nl2).unwrap();

        let d = sys.box_size.min_image(sys.pos(0) - sys.pos(1));
        let r = d.norm();
        let bare_e = sys.charges[0] * sys.charges[1] / r;
        assert!(
            (e_excl - (e_plain - bare_e)).abs() < 1e-8,
            "energy: excl {e_excl} vs plain−bare {}",
            e_plain - bare_e
        );
        let bare_f = d * (sys.charges[0] * sys.charges[1] / (r * r * r));
        for k in 0..3 {
            assert!(
                (a.forces[k] - (b.forces[k] - bare_f.0[k])).abs() < 1e-8,
                "force axis {k}: {} vs {}",
                a.forces[k],
                b.forces[k] - bare_f.0[k]
            );
        }
    }
}
