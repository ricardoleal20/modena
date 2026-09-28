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

/// Rigid SPC/E water on a simple-cubic lattice at `density` (g/cm³). The box
/// edge is computed in Å from the SPC/E molecular mass (15.9994 + 2·1.008
/// amu) and then scaled to the simulation's length units, so the requested
/// density is exact; use perfect-cube `n_waters` (64, 216, 512, 1000) for a
/// full lattice fill.
///
/// Each molecule is a **contiguous (O, H, H) triple** — the layout the
/// rigid-water constraint projection chunks over. Orientation is a balanced
/// alternating-dipole pattern (site-parity flip + a deterministic per-site
/// twist): a cold *random*-orientation lattice collapses into a deep artifact
/// minimum (≈ −200 kcal/mol·water on both MODENA and LAMMPS), while the
/// balanced pattern starts around −5 kcal/mol·water. Species: 0 = O, 1 = H.
/// `d_oh`/`q_*`/`m_*` arrive already converted to the simulation's units
/// (reduced under `units = "real"`); `hoh_angle` is unitless degrees. `seed`
/// is accepted for API symmetry (the orientation is RNG-free).
#[allow(clippy::too_many_arguments)]
pub fn water_lattice(
    n_waters: usize,
    density_g_cm3: f64,
    d_oh: f64,
    hoh_angle_deg: f64,
    q_o: f64,
    q_h: f64,
    m_o: f64,
    m_h: f64,
    seed: u64,
    to_reduced_length: f64,
) -> Result<System> {
    if n_waters == 0 {
        return Err(ModenaError::Config("n_waters must be >= 1".into()));
    }
    if !(density_g_cm3 > 0.0) {
        return Err(ModenaError::Config("density must be > 0 (g/cm³)".into()));
    }
    if !(d_oh > 0.0) {
        return Err(ModenaError::Config("d_oh must be > 0".into()));
    }
    // Box edge: Å from the physical mass and density, then to sim units.
    let m_amu = 15.9994 + 2.0 * 1.008;
    let vol_a3 = n_waters as f64 * m_amu * 1.66053906660 / density_g_cm3;
    let edge = vol_a3.cbrt() * to_reduced_length;

    let n_side = (n_waters as f64).cbrt().ceil() as usize;
    let spacing = edge / n_side as f64;

    // Local-frame geometry: O at the origin, H's symmetric about +x.
    let half = hoh_angle_deg.to_radians() * 0.5;
    let h1l = Vec3::new(d_oh * half.cos(), d_oh * half.sin(), 0.0);
    let h2l = Vec3::new(d_oh * half.cos(), -d_oh * half.sin(), 0.0);

    let mut system = System::new(SimBox::cubic(edge));
    system.resize(3 * n_waters);
    let _ = seed;
    let mut idx = 0usize;
    'outer: for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                if idx >= n_waters {
                    break 'outer;
                }
                let site = Vec3::new(
                    spacing * ix as f64,
                    spacing * iy as f64,
                    spacing * iz as f64,
                );
                // Balanced alternating-dipole orientation: dipole along ±x by
                // site parity, H's scattered in the yz-plane by a deterministic
                // per-site twist.
                let parity = (ix + iy + iz) & 1;
                let flip = if parity == 0 { 1.0 } else { -1.0 };
                let twist = ((ix * 7 + iy * 11 + iz * 13) % 24) as f64 * 15.0_f64.to_radians();
                let (ct, st) = (twist.cos(), twist.sin());
                let u = Vec3::new(flip, 0.0, 0.0);
                let v = Vec3::new(0.0, ct, st);
                let w = Vec3::new(0.0, -st, ct);
                let rot = |p: Vec3| u * p.0[0] + v * p.0[1] + w * p.0[2];
                system.set_pos(3 * idx, site);
                system.set_pos(3 * idx + 1, site + rot(h1l));
                system.set_pos(3 * idx + 2, site + rot(h2l));
                system.masses[3 * idx] = m_o;
                system.masses[3 * idx + 1] = m_h;
                system.masses[3 * idx + 2] = m_h;
                system.charges[3 * idx] = q_o;
                system.charges[3 * idx + 1] = q_h;
                system.charges[3 * idx + 2] = q_h;
                system.species[3 * idx] = 0;
                system.species[3 * idx + 1] = 1;
                system.species[3 * idx + 2] = 1;
                system.elements[3 * idx] = "O".into();
                system.elements[3 * idx + 1] = "H".into();
                system.elements[3 * idx + 2] = "H".into();
                idx += 1;
            }
        }
    }
    system.wrap_positions();
    Ok(system)
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
