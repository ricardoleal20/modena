//! Long-range electrostatics + triclinic box integration tests (M4).
//!
//! Four correctness anchors:
//! 1. **Force = −∇U by finite differences** — every atom's Ewald force matches
//!    the central-difference gradient of the Ewald energy. This is the
//!    sign/correctness check for the real-space, reciprocal, and (force-free)
//!    self/background terms together; a sign error or a missing factor in any
//!    term shows up here.
//! 2. **Madelung lattice energy of NaCl** — the Ewald energy per ion pair of a
//!    rocksalt crystal converges to the analytic Madelung value
//!    `−α q²/r₀` (α = 1.74756), the quantitative gold standard for Ewald.
//! 3. **Triclinic minimum image** — `SimBox::min_image`/`wrap` are correct for a
//!    general cell tensor, and a diagonal triclinic cell matches the
//!    orthorhombic path exactly.
//! 4. **Ionic NVE energy conservation** — a rocksalt crystal (Coulomb + a
//!    short-range repulsive core) run under NVE conserves total energy, the
//!    integration-level check that the Ewald forces are conservative and wired
//!    through the runner/composite without clobbering the pair model.

use modena_core::{
    lattice, EwaldForceModel, ForceModel, NeighborList, SimBox, Simulation, System, Vec3, Virial,
};
use std::path::Path;

/// Build a small neutral ion cluster with non-symmetric positions (so no force
/// vanishes by symmetry and the finite-difference check exercises every term).
fn ion_cluster() -> System {
    let mut sys = System::new(SimBox::cubic(6.0));
    sys.resize(6);
    let pos = [
        Vec3::new(0.9, 1.1, 1.0),
        Vec3::new(2.3, 0.7, 1.4),
        Vec3::new(1.2, 2.6, 0.8),
        Vec3::new(2.8, 2.9, 2.1),
        Vec3::new(4.1, 1.3, 2.7),
        Vec3::new(3.6, 3.8, 1.2),
    ];
    let q = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0];
    for (i, p) in pos.iter().enumerate() {
        sys.set_pos(i, *p);
        sys.charges[i] = q[i];
    }
    sys
}

fn ewald() -> EwaldForceModel {
    EwaldForceModel::new(1.0, 3.0, 5.0).unwrap()
}

/// Ewald energy of a system (force-free evaluation).
fn ewald_energy(ewald: &EwaldForceModel, system: &System) -> f64 {
    let mut nl = NeighborList::new(3.0);
    ewald.energy(system, &mut nl).unwrap()
}

#[test]
fn ewald_force_is_negative_gradient() {
    let ewald = ewald();
    let mut sys = ion_cluster();

    // Computed forces.
    sys.zero_forces();
    let mut virial = Virial::default();
    let mut nl = NeighborList::new(3.0);
    ewald.compute(&mut sys, &mut virial, &mut nl).unwrap();

    let h = 1e-5;
    let n = sys.n();
    for i in 0..n {
        for k in 0..3 {
            let mut sp = sys.clone();
            {
                let p = sp.pos(i);
                let mut a = p.0;
                a[k] += h;
                sp.set_pos(i, Vec3(a));
            }
            let mut sm = sys.clone();
            {
                let p = sm.pos(i);
                let mut a = p.0;
                a[k] -= h;
                sm.set_pos(i, Vec3(a));
            }
            let ep = ewald_energy(&ewald, &sp);
            let em = ewald_energy(&ewald, &sm);
            let fd = -(ep - em) / (2.0 * h); // central FD of -dU/dx
            let f = sys.force(i).0[k];
            if fd.abs() > 1e-6 {
                let rel = (f - fd).abs() / fd.abs().max(1e-9);
                assert!(
                    rel < 1e-3,
                    "atom {i} axis {k}: force {f:.6e} vs FD {fd:.6e} (rel {rel:.3e})"
                );
            } else {
                assert!(
                    (f - fd).abs() < 1e-4,
                    "atom {i} axis {k}: force {f:.6e} vs FD {fd:.6e}"
                );
            }
        }
    }
}

#[test]
fn ewald_madelung_nacl() {
    // NaCl: a = 2.0 → nearest cation–anion distance r0 = a/2 = 1.0, q = ±1.
    // Madelung energy per ion pair = -α q² / r0, α = 1.747558775.
    let alpha = 1.747558775;
    let r0 = 1.0_f64;
    let madelung_per_pair = -alpha * r0.recip(); // q² = 1

    let sys = lattice::rocksalt(2.0, [3, 3, 3], "Na", "Cl", 1.0, 1.0, 1.0, -1.0);
    let n_pairs = sys.n() / 2;
    // κ = 1.6, real cutoff 3.0 (< L/2 = 3.0), k-space cutoff 6.0 → well
    // converged (the erfc tail and the k-space truncation are both < 1e-3).
    let ewald = EwaldForceModel::new(1.6, 2.9, 6.0).unwrap();
    let e = ewald_energy(&ewald, &sys);
    let per_pair = e / n_pairs as f64;
    let err = (per_pair - madelung_per_pair).abs();
    assert!(
        err < 0.02,
        "Madelung: E/pair = {per_pair:.5}, expected {madelung_per_pair:.5} (err {err:.4})"
    );
}

#[test]
fn triclinic_min_image_matches_orthorhombic_for_diagonal_cell() {
    // A triclinic cell whose vectors are axis-aligned must reproduce the
    // orthorhombic minimum image exactly.
    let ortho = SimBox::cubic(5.0);
    let tri = SimBox::triclinic([
        Vec3::new(5.0, 0.0, 0.0),
        Vec3::new(0.0, 5.0, 0.0),
        Vec3::new(0.0, 0.0, 5.0),
    ]);
    for d in [
        Vec3::new(2.9, -2.6, 1.3),
        Vec3::new(4.5, 4.5, 4.5),
        Vec3::new(-3.7, 2.1, -4.9),
        Vec3::new(0.1, -0.1, 0.2),
    ] {
        let a = ortho.min_image(d);
        let b = tri.min_image(d);
        assert!(
            (a - b).norm() < 1e-12,
            "diagonal triclinic != orthorhombic for {d:?}: {a:?} vs {b:?}"
        );
    }
}

#[test]
fn triclinic_min_image_wraps_lattice_vectors() {
    // A displacement equal to a box vector is a lattice translation → wraps to
    // ~0. General (sheared) cell.
    let a = Vec3::new(3.0, 0.0, 0.0);
    let b = Vec3::new(0.5, 3.0, 0.0);
    let c = Vec3::new(0.0, 0.5, 3.0);
    let box_size = SimBox::triclinic([a, b, c]);
    for v in [a, b, c, a + b, b + c, a + b + c, a * 2.0 - c] {
        let w = box_size.min_image(v);
        assert!(
            w.norm() < 1e-10,
            "lattice vector {v:?} did not wrap to ~0: {w:?} (|w|={})",
            w.norm()
        );
    }
    // A sub-period displacement is returned unwrapped (it's already the min
    // image): d = 0.4·a should come back ≈ 0.4·a.
    let d = a * 0.4;
    let w = box_size.min_image(d);
    assert!(
        (w - d).norm() < 1e-12,
        "sub-period displacement changed: {w:?}"
    );
}

#[test]
fn triclinic_wrap_brings_points_into_cell() {
    let a = Vec3::new(3.0, 0.0, 0.0);
    let b = Vec3::new(0.5, 3.0, 0.0);
    let c = Vec3::new(0.0, 0.5, 3.0);
    let box_size = SimBox::triclinic([a, b, c]);
    let inv = box_size.inv_basis_rows();
    // A point clearly outside the cell wraps to one whose fractional coords are
    // in [0,1).
    let outside = Vec3::new(7.5, -2.0, 4.0);
    let inside = box_size.wrap(outside);
    let s0 = inv[0].dot(inside);
    let s1 = inv[1].dot(inside);
    let s2 = inv[2].dot(inside);
    for (name, s) in [("s0", s0), ("s1", s1), ("s2", s2)] {
        assert!(
            (-1e-10..1.0 + 1e-10).contains(&s),
            "wrapped fractional {name} = {s} not in [0,1)"
        );
    }
    // Wrapping is a pure lattice translation: difference is an integer combo.
    let diff = outside - inside;
    let f0 = inv[0].dot(diff).round();
    let f1 = inv[1].dot(diff).round();
    let f2 = inv[2].dot(diff).round();
    let reconstructed = a * f0 + b * f1 + c * f2;
    assert!(
        (diff - reconstructed).norm() < 1e-10,
        "wrap translation not a lattice vector"
    );
}

#[test]
fn ewald_ionic_nve_energy_conservation() {
    // NaCl with a repulsive core + Ewald, integrated under NVE. The Ewald
    // forces are conservative → total energy is conserved (the integration-level
    // check that Ewald composes with the pair model and the runner without
    // clobbering forces/virial).
    let toml = r#"
units = "reduced"

[structure]
type = "rock_salt"
lattice_constant = 2.2
replicas = [2, 2, 2]
cation = "Na"
anion = "Cl"
cation_charge = 1.0
anion_charge = -1.0

[model]
type = "soft_sphere"
cutoff = 2.0
[[model.types]]
sigma = 0.9
epsilon = 0.5
n = 12
[[model.types]]
sigma = 0.9
epsilon = 0.5
n = 12

[electrostatics]
kappa = 1.0
real_cutoff = 2.0
kspace_cutoff = 4.0

[integrator]
type = "velocity_verlet"
dt = 0.004

[run]
steps = 800
initial_temperature = 0.01
seed = 7
"#;
    let config = modena_core::InputConfig::from_str(toml).expect("config must parse");
    let mut sim = Simulation::from_config(&config, Path::new(".")).expect("sim must build");

    // Run NVE and read the energy range from the conservation report.
    let summary = sim.run().expect("run must complete");
    let drift = summary.conservation.relative_drift;
    assert!(
        drift < 1e-2,
        "ionic NVE relative energy drift too large: {drift:.3e}"
    );
    // The run must not have diverged (finite forces/energies throughout).
    assert!(summary.final_energies.total.is_finite());
}

#[test]
fn spme_ionic_nve_energy_conservation() {
    // Same ionic deck as the direct-Ewald conservation test, with
    // `method = "spme"` — the mesh path must conserve energy just as well
    // (the ik forces are the exact derivative of the mesh energy functional).
    let toml = r#"
units = "reduced"

[structure]
type = "rock_salt"
lattice_constant = 2.2
replicas = [2, 2, 2]
cation = "Na"
anion = "Cl"
cation_charge = 1.0
anion_charge = -1.0

[model]
type = "soft_sphere"
cutoff = 2.0
[[model.types]]
sigma = 0.9
epsilon = 0.5
n = 12
[[model.types]]
sigma = 0.9
epsilon = 0.5
n = 12

[electrostatics]
method = "spme"
kappa = 1.0
real_cutoff = 2.0
kspace_cutoff = 4.0

[integrator]
type = "velocity_verlet"
dt = 0.004

[run]
steps = 800
initial_temperature = 0.01
seed = 7
"#;
    let config = modena_core::InputConfig::from_str(toml).expect("config must parse");
    let mut sim = Simulation::from_config(&config, Path::new(".")).expect("sim must build");

    let summary = sim.run().expect("run must complete");
    let drift = summary.conservation.relative_drift;
    assert!(
        drift < 1e-2,
        "SPME ionic NVE relative energy drift too large: {drift:.3e}"
    );
    assert!(summary.final_energies.total.is_finite());
}

#[test]
fn spme_matches_ewald_deck_energy() {
    // The same deck under both electrostatics methods: the SPME mesh (6^3 for
    // this box/cutoff) must reproduce the direct Ewald total energy within
    // mesh error on the initial (undisplaced, T=0) configuration.
    let base = r#"
units = "reduced"

[structure]
type = "rock_salt"
lattice_constant = 2.2
replicas = [2, 2, 2]
cation = "Na"
anion = "Cl"
cation_charge = 1.0
anion_charge = -1.0

[model]
type = "soft_sphere"
cutoff = 2.0
[[model.types]]
sigma = 0.9
epsilon = 0.5
n = 12
[[model.types]]
sigma = 0.9
epsilon = 0.5
n = 12

[electrostatics]
kappa = 1.0
real_cutoff = 2.0
kspace_cutoff = 4.0

[integrator]
type = "velocity_verlet"
dt = 0.004

[run]
steps = 0
initial_temperature = 0.0
seed = 7
"#;
    let run = |method: &str| -> f64 {
        let toml = if method.is_empty() {
            base.to_string()
        } else {
            base.replacen(
                "[electrostatics]",
                &format!("[electrostatics]\nmethod = \"{method}\""),
                1,
            )
        };
        let config = modena_core::InputConfig::from_str(&toml).expect("config must parse");
        let mut sim = Simulation::from_config(&config, Path::new(".")).expect("sim must build");
        let summary = sim.run().expect("run must complete");
        summary.final_energies.total
    };
    let e_ewald = run("");
    let e_spme = run("spme");
    assert!(
        (e_spme - e_ewald).abs() < 1e-2 * e_ewald.abs(),
        "spme {e_spme} vs ewald {e_ewald}"
    );
}

#[test]
fn ewald_real_space_is_bare_coulomb_for_small_kappa() {
    // As κ → 0, erfc(κr) → 1 and the screening term ∝ κ → 0, so the real-space
    // force approaches bare Coulomb q_i q_j / r². Two ions, one pair.
    let mut sys = System::new(SimBox::cubic(30.0));
    sys.resize(2);
    sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
    sys.set_pos(1, Vec3::new(2.0, 0.0, 0.0));
    sys.charges[0] = 1.0;
    sys.charges[1] = 1.0;
    let r = 2.0_f64;
    let ewald = EwaldForceModel::new(1e-4, 6.0, 0.5).unwrap();
    // Tiny k-space cutoff + tiny κ → reciprocal and self terms ~0; real-space
    // dominates. Like charges repel → atom 0 pushed -x (away from atom 1 at +x).
    sys.zero_forces();
    let mut virial = Virial::default();
    let mut nl = NeighborList::new(6.0);
    ewald.compute(&mut sys, &mut virial, &mut nl).unwrap();
    let fx = sys.force(0).0[0];
    let bare = -1.0 / (r * r); // repulsion: atom 0 driven -x
    assert!(
        (fx - bare).abs() / bare.abs() < 1e-3,
        "real-space force {fx} vs bare Coulomb {bare}"
    );
}
