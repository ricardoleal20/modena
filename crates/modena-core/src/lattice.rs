//! Crystal lattice builders.

use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;

use crate::error::{ModenaError, Result};
use crate::system::{SimBox, System};
use crate::vec3::Vec3;

/// One species entry for a [`random_alloy`] builder.
#[derive(Debug, Clone)]
pub struct AlloySpecies {
    pub element: String,
    pub mass: f64,
    /// Fraction of atoms (0..1); must sum to 1 across all species.
    pub concentration: f64,
    /// Per-atom charge for this species (M4 electrostatics; reduced units).
    /// Defaults to 0 (uncharged).
    pub charge: f64,
}

impl AlloySpecies {
    /// Convenience constructor for an uncharged species.
    pub fn neutral(element: &str, mass: f64, concentration: f64) -> Self {
        Self {
            element: element.to_string(),
            mass,
            concentration,
            charge: 0.0,
        }
    }
}

/// Build an FCC lattice with `replicas` unit cells per axis and lattice
/// constant `a`. 4 atoms per cell. Box is orthorhombic `a·replicas`.
pub fn fcc(a: f64, replicas: [u32; 3], element: &str) -> System {
    let basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];
    let [nx, ny, nz] = replicas;
    let box_size = SimBox::new(Vec3::new(a * nx as f64, a * ny as f64, a * nz as f64));
    let mut system = System::new(box_size);
    let total = (nx * ny * nz * 4) as usize;
    system.resize(total);
    let mut idx = 0;
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                for b in &basis {
                    system.set_pos(
                        idx,
                        Vec3::new(
                            a * (ix as f64 + b.0[0]),
                            a * (iy as f64 + b.0[1]),
                            a * (iz as f64 + b.0[2]),
                        ),
                    );
                    system.elements[idx] = element.to_string();
                    idx += 1;
                }
            }
        }
    }
    system
}

/// Build an FCC lattice and assign each atom a species by concentration, using a
/// seeded shuffle so the assignment is deterministic and reproducible. Per-atom
/// masses and element labels are set from the chosen species. Concentrations
/// must sum to 1 (within tolerance).
pub fn random_alloy(
    a: f64,
    replicas: [u32; 3],
    species: &[AlloySpecies],
    seed: u64,
) -> Result<System> {
    if species.is_empty() {
        return Err(ModenaError::Config(
            "random_alloy needs at least one species".into(),
        ));
    }
    let total_conc: f64 = species.iter().map(|s| s.concentration).sum();
    if (total_conc - 1.0).abs() > 1e-6 {
        return Err(ModenaError::Config(format!(
            "random_alloy concentrations must sum to 1.0, got {total_conc}"
        )));
    }
    let mut system = fcc(a, replicas, &species[0].element);
    let n = system.n();

    // Build a species label per atom by drawing `count` slots for each species.
    // Use a deterministic shuffle: fill an array of species indices in
    // concentration order, then shuffle with the seeded RNG. The last species
    // absorbs rounding so the total is exactly `n`.
    let mut labels: Vec<u32> = Vec::with_capacity(n);
    for (i, s) in species.iter().enumerate() {
        // Last species takes the remainder to guarantee the count is exact.
        let count = if i == species.len() - 1 {
            n - labels.len()
        } else {
            (s.concentration * n as f64).round() as usize
        };
        labels.extend(std::iter::repeat_n(i as u32, count));
    }
    // Guard against the rounding producing the wrong length.
    labels.resize(n, (species.len() - 1) as u32);

    let mut rng = StdRng::seed_from_u64(seed);
    labels.shuffle(&mut rng);

    for (i, &sp) in labels.iter().enumerate() {
        let s = &species[sp as usize];
        system.species[i] = sp;
        system.masses[i] = s.mass;
        system.charges[i] = s.charge;
        system.elements[i] = s.element.clone();
    }
    Ok(system)
}

/// Build a rocksalt (NaCl) lattice: two interpenetrating FCC sublattices. The
/// cation occupies the FCC basis; the anion occupies the same basis shifted by
/// `(a/2, 0, 0)` (the octahedral hole). Nearest cation–anion distance is `a/2`.
/// Cations get `cation_charge`, anions `anion_charge` (M4 electrostatics). The
/// conventional cell holds 4 cations + 4 anions (8 atoms); `replicas` tiles it.
#[allow(clippy::too_many_arguments)]
pub fn rocksalt(
    a: f64,
    replicas: [u32; 3],
    cation_element: &str,
    anion_element: &str,
    cation_mass: f64,
    anion_mass: f64,
    cation_charge: f64,
    anion_charge: f64,
) -> System {
    let fcc_basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];
    let shift = Vec3::new(0.5, 0.0, 0.0);
    let [nx, ny, nz] = replicas;
    let box_size = SimBox::new(Vec3::new(a * nx as f64, a * ny as f64, a * nz as f64));
    let mut system = System::new(box_size);
    let total = (nx * ny * nz * 8) as usize;
    system.resize(total);
    let mut idx = 0;
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let origin = Vec3::new(a * ix as f64, a * iy as f64, a * iz as f64);
                for b in &fcc_basis {
                    // Cation.
                    let p = origin + *b * a;
                    system.set_pos(idx, p);
                    system.masses[idx] = cation_mass;
                    system.charges[idx] = cation_charge;
                    system.species[idx] = 0;
                    system.elements[idx] = cation_element.to_string();
                    idx += 1;
                    // Anion (cation position + (a/2,0,0)).
                    let p = origin + (*b + shift) * a;
                    system.set_pos(idx, p);
                    system.masses[idx] = anion_mass;
                    system.charges[idx] = anion_charge;
                    system.species[idx] = 1;
                    system.elements[idx] = anion_element.to_string();
                    idx += 1;
                }
            }
        }
    }
    system
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fcc_atom_count() {
        let sys = fcc(1.0, [2, 2, 2], "Ar");
        assert_eq!(sys.n(), 2 * 2 * 2 * 4);
        assert_eq!(sys.elements[0], "Ar");
        assert_eq!(sys.species[0], 0);
        assert_eq!(sys.masses[0], 1.0);
    }

    #[test]
    fn random_alloy_counts_and_masses() {
        let sys = random_alloy(
            1.0,
            [3, 3, 3],
            &[
                AlloySpecies::neutral("Ar", 39.948, 0.8),
                AlloySpecies::neutral("Ne", 20.180, 0.2),
            ],
            42,
        )
        .unwrap();
        let n = sys.n();
        assert_eq!(n, 3 * 3 * 3 * 4);
        let n_ar = sys.elements.iter().filter(|e| *e == "Ar").count();
        let n_ne = sys.elements.iter().filter(|e| *e == "Ne").count();
        assert_eq!(n_ar + n_ne, n);
        // ~80/20 split (3*3*3*4 = 108 atoms → 86 Ar, 22 Ne after rounding).
        assert!((n_ar as f64 / n as f64 - 0.8).abs() < 0.05);
        // Masses follow species.
        for i in 0..n {
            let expected = if sys.elements[i] == "Ar" {
                39.948
            } else {
                20.180
            };
            assert!((sys.masses[i] - expected).abs() < 1e-9);
        }
        // Deterministic: same seed → same assignment.
        let sys2 = random_alloy(
            1.0,
            [3, 3, 3],
            &[
                AlloySpecies::neutral("Ar", 39.948, 0.8),
                AlloySpecies::neutral("Ne", 20.180, 0.2),
            ],
            42,
        )
        .unwrap();
        assert_eq!(sys.species, sys2.species);
    }

    #[test]
    fn random_alloy_rejects_bad_concentration() {
        let err = random_alloy(1.0, [2, 2, 2], &[AlloySpecies::neutral("Ar", 1.0, 0.5)], 42);
        assert!(err.is_err());
    }

    #[test]
    fn rocksalt_structure_and_charges() {
        let sys = rocksalt(2.0, [2, 2, 2], "Na", "Cl", 1.0, 1.0, 1.0, -1.0);
        // 2×2×2 conventional cells × 8 atoms = 64 atoms; 32 cations, 32 anions.
        assert_eq!(sys.n(), 64);
        let n_cat = sys.charges.iter().filter(|&&q| q == 1.0).count();
        let n_an = sys.charges.iter().filter(|&&q| q == -1.0).count();
        assert_eq!(n_cat, 32);
        assert_eq!(n_an, 32);
        // Neutral.
        let total_q: f64 = sys.charges.iter().sum();
        assert!(total_q.abs() < 1e-12);
        // Nearest cation–anion distance is a/2 = 1.0 (atom 0 cation at origin,
        // atom 1 anion at (a/2,0,0)).
        let d = sys.pos(0) - sys.pos(1);
        assert!((d.norm() - 1.0).abs() < 1e-12);
    }
}
