use modena_core::{lattice, ForceModel, LennardJones, PairForceModel, SpmeForceModel, System};
/// Static f64 orthorhombic fixtures; molecular electrostatics is not a complete force field.
pub fn workloads() -> Vec<(String, System, Box<dyn ForceModel>, f64)> {
    let mut result: Vec<(String, System, Box<dyn ForceModel>, f64)> = Vec::new();
    for n in [4_u32, 6, 8] {
        result.push((
            format!("lj_N{}", 4 * n * n * n),
            lattice::fcc(1.6796, [n; 3], "Ar"),
            Box::new(PairForceModel::new(
                LennardJones::new(1.0, 1.0, 2.5).unwrap(),
            )),
            2.5,
        ));
    }
    for n in [2_u32, 3, 4] {
        let ions = lattice::rocksalt(5.0, [n; 3], "Na", "Cl", 1.0, 1.0, 1.0, -1.0);
        result.push((
            format!("ionic_spme_N{}", ions.n()),
            ions,
            Box::new(SpmeForceModel::new(0.8, 4.2, 6.4).unwrap()),
            4.2,
        ));
    }
    for count in [27, 64, 125] {
        let water = lattice::water_lattice(
            count, 1.0, 1.0, 109.47, -0.8476, 0.4238, 15.9994, 1.008, 87287, 1.0,
        )
        .unwrap();
        let excluded = (0..count)
            .flat_map(|i| {
                [
                    (3 * i, 3 * i + 1),
                    (3 * i, 3 * i + 2),
                    (3 * i + 1, 3 * i + 2),
                ]
            })
            .collect();
        result.push((
            format!("molecular_spme_N{}", water.n()),
            water,
            Box::new(
                SpmeForceModel::new(0.34, 4.2, 3.2)
                    .unwrap()
                    .with_excluded_pairs(excluded),
            ),
            4.2,
        ));
    }
    result
}
