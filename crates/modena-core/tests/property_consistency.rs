//! Property-based tests enforcing the conservation contract (Feature B).
//!
//! These verify, for the registered `PairPotential`(s), that:
//!   1. The per-atom force equals −∇U (finite difference over the full system).
//!   2. The pair kernel is self-consistent: dU/dr = −F(r) (finite difference on
//!      `energy_force` directly).
//!   3. Outputs are finite (no NaN/Inf) for random in-range distances.
//!   4. The shifted-force form is continuous at the cutoff (F → 0, U continuous).
//!   5. A short NVE run conserves energy within a symplectic-order bound.
//!
//! If any of these break, the engine's headline correctness guarantee — that
//! the reported PE is the exact integral of the reported force — is violated,
//! and the NVE drift signal becomes untrustworthy.

use modena_core::{
    lattice, ForceModel, Integrator, LennardJones, NeighborList, PairForceModel, PairPotential,
    System, Vec3, VelocityVerlet, Virial,
};
use proptest::prelude::*;

fn lj() -> LennardJones {
    LennardJones::new(1.0, 1.0, 2.5).unwrap()
}

/// Build a non-singular configuration: a 3×3×3 FCC lattice (108 atoms) with a
/// per-axis jitter in `[-jitter, jitter]`. Jitter is small enough that no two
/// atoms approach the r → 0 singularity.
///
/// **Why 3×3×3, not 2×2×2:** the box must satisfy `cutoff < L/2` so that every
/// interacting pair has a *unique, perturbation-stable* minimum image. With the
/// standard LJ cutoff 2.5 and lattice constant 1.6796, a 2×2×2 box has L/2 =
/// 1.68 < 2.5 — ill-posed: pairs sit on the L/2 symmetry boundary where an
/// infinitesimal perturbation flips the min-image vector, making the potential
/// energy discontinuous and force ≠ −∇U *by finite difference* even though the
/// kernel is correct. 3×3×3 gives L/2 = 2.52 > 2.5, so no interacting pair
/// (distance ≤ cutoff < L/2) can cross the min-image boundary; PE is smooth and
/// the analytic force is its exact gradient.
fn jittered_system(jitters: &[(f64, f64, f64)], jitter: f64) -> System {
    let mut sys = lattice::fcc(1.6796, [3, 3, 3], "Ar");
    assert_eq!(sys.n(), 108);
    for (i, &(dx, dy, dz)) in jitters.iter().enumerate() {
        // Clamp to the promised bound defensively.
        sys.positions[i] += Vec3::new(
            dx.clamp(-jitter, jitter),
            dy.clamp(-jitter, jitter),
            dz.clamp(-jitter, jitter),
        );
    }
    sys.wrap_positions();
    sys
}

fn total_pe(sys: &System, model: &dyn ForceModel, neighbors: &mut NeighborList) -> f64 {
    let mut v = Virial::default();
    // compute overwrites forces; we only need the returned PE here.
    model
        .compute(&mut sys.clone(), &mut v, neighbors)
        .expect("force eval")
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    /// Per-atom force == −∇U, by central finite difference on the total PE.
    #[test]
    fn force_is_negative_gradient_of_energy(
        jitters in proptest::collection::vec(
            (any::<f64>(), any::<f64>(), any::<f64>()),
            108,
        )
    ) {
        let pot = lj();
        let cutoff = pot.cutoff();
        let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
        let sys = jittered_system(&jitters, 0.10);
        let mut neighbors = NeighborList::new(cutoff);

        // Reference forces.
        let mut virial = Virial::default();
        let mut sys_for_forces = sys.clone();
        model.compute(&mut sys_for_forces, &mut virial, &mut neighbors).unwrap();
        let forces: Vec<Vec3> = sys_for_forces.forces.clone();

        let h = 1e-5_f64;
        let mut max_err = 0.0_f64;
        for (i, _) in forces.iter().enumerate() {
            for k in 0..3 {
                let mut sp = sys.clone();
                let mut sm = sys.clone();
                sp.positions[i].0[k] += h;
                sm.positions[i].0[k] -= h;
                let pe_p = total_pe(&sp, &model, &mut neighbors);
                let pe_m = total_pe(&sm, &model, &mut neighbors);
                // F = -dU/dx  =>  dU/dx = (pe_p - pe_m)/(2h); force = -that.
                let fd_force = -(pe_p - pe_m) / (2.0 * h);
                let err = (forces[i].0[k] - fd_force).abs();
                let scale = 1.0 + forces[i].0[k].abs();
                prop_assert!(
                    err < 1e-3 * scale,
                    "atom {} axis {}: force={}, fd={}, err={}",
                    i, k, forces[i].0[k], fd_force, err
                );
                max_err = max_err.max(err / scale);
            }
        }
        let _ = max_err;
    }

    /// Pair kernel self-consistency: dU/dr = -F(r) by finite difference.
    #[test]
    fn pair_energy_derivative_equals_negative_force(r in 0.9_f64..2.49_f64) {
        let pot = lj();
        let dr = 1e-6_f64;
        let u_at = |rr: f64| -> f64 { pot.energy_force(rr * rr, 0, 0).0 };
        let (u, f_over_r) = pot.energy_force(r * r, 0, 0);
        let f_radial = f_over_r * r;
        let du_dr = (u_at(r + dr) - u_at(r - dr)) / (2.0 * dr);
        let _ = u;
        // Contract: dU/dr = -F(r).
        prop_assert!(
            (du_dr - (-f_radial)).abs() < 1e-4 * (1.0 + f_radial.abs()),
            "r={}: dU/dr={}, -F={}", r, du_dr, -f_radial
        );
    }

    /// Outputs are finite for all in-range distances (no NaN/Inf leak).
    #[test]
    fn outputs_are_finite(r2 in 1e-8_f64..(2.5f64 * 2.5)) {
        let pot = lj();
        let (u, f_over_r) = pot.energy_force(r2, 0, 0);
        prop_assert!(u.is_finite(), "U not finite at r2={r2}: {u}");
        prop_assert!(f_over_r.is_finite(), "f_over_r not finite at r2={r2}: {f_over_r}");
    }

    /// Shifted-force continuity at the cutoff: F → 0 and U is continuous.
    ///
    /// The shifted-force form sets `F_sf(rc) = 0` and `U_sf(rc) = 0` exactly.
    /// Just inside, `F_sf(rc − eps) = F(rc − eps) − F(rc) = −eps·F'(rc) + O(eps²)`,
    /// so the residual is **O(eps)**, not machine-zero — the bound must scale
    /// with `eps`. `U_sf` is C1 at the cutoff (`U_sf'(rc) = −F_sf(rc) = 0`), so
    /// its residual is O(eps²). The bounds below are tight enough to catch a
    /// real discontinuity (an unshifted potential jumps by |U(rc)| ≈ 0.016 and
    /// |F(rc)| ≈ 0.039 at the cutoff), yet accommodate the legitimate residual.
    #[test]
    fn cutoff_continuity(eps in 1e-9_f64..1e-3_f64) {
        let pot = lj();
        let rc = pot.cutoff();
        // Force just inside the cutoff scales as O(eps) (shifted-force residual).
        let (_, f_over_r_in) = pot.energy_force((rc - eps).powi(2), 0, 0);
        let f_in = f_over_r_in * (rc - eps);
        prop_assert!(
            f_in.abs() < 10.0 * eps,
            "F just inside cutoff = {f_in}, expected O(eps) ~ {}",
            10.0 * eps
        );
        // U at the cutoff is ~0 by construction.
        let (u_rc, _) = pot.energy_force(rc * rc, 0, 0);
        prop_assert!(u_rc.abs() < 1e-9, "U(rc) = {u_rc}, expected 0");
        // U continuous (C1): residual is O(eps²), well under a real jump (~0.016).
        let u_just_below = pot.energy_force((rc - eps).powi(2), 0, 0).0;
        prop_assert!(
            (u_just_below - u_rc).abs() < 1e-4,
            "U discontinuous at cutoff: below={u_just_below}, at={u_rc}"
        );
    }
}

/// A short NVE run conserves energy within a symplectic bound. Parametrized over
/// timestep: velocity-Verlet is 4th-order in the shadow Hamiltonian, so drift
/// stays well under 1e-3 for these stable configs.
#[test]
fn nve_drift_within_symplectic_bound() {
    for &(dt, bound) in &[(0.002_f64, 1e-4), (0.005, 1e-3)] {
        let pot = lj();
        let cutoff = pot.cutoff();
        let model = PairForceModel::new(pot);
        let mut system = lattice::fcc(1.6796, [2, 2, 2], "Ar");
        for i in 0..system.n() {
            let fi = i as f64;
            system.velocities[i] = Vec3::new(
                0.5 * fi.sin(),
                0.5 * (1.3 * fi).cos(),
                0.5 * (0.7 * fi).sin(),
            );
        }
        let mut com = Vec3::zero();
        for v in &system.velocities {
            com += *v;
        }
        com *= 1.0 / system.n() as f64;
        for v in system.velocities.iter_mut() {
            *v -= com;
        }

        let mut neighbors = NeighborList::new(cutoff);
        let mut integrator = VelocityVerlet::new();
        let mut virial = Virial::default();
        let _ = model.compute(&mut system, &mut virial, &mut neighbors);

        let mut e0 = None;
        let mut e_last = 0.0;
        for _ in 0..1000 {
            let e = integrator
                .step(&mut system, &model, &mut neighbors, dt)
                .expect("step");
            if e0.is_none() {
                e0 = Some(e.total);
            }
            e_last = e.total;
        }
        let e0 = e0.unwrap();
        let drift = ((e_last - e0) / e0.abs().max(1e-12)).abs();
        assert!(
            drift < bound,
            "dt={dt}: drift {drift:.3e} exceeds bound {bound:.0e}"
        );
        // The conservation contract's verdict should agree.
        let report = modena_core::ConservationReport::new(e0, e_last);
        assert!(
            report.verdict != modena_core::ConservationVerdict::Drifting,
            "dt={dt}: verdict should not be Drifting, got {:?} (drift {drift:.3e})",
            report.verdict
        );
    }
}
