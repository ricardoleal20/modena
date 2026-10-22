use modena_core::{InputConfig, Simulation};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
static SEQ: AtomicUsize = AtomicUsize::new(0);
fn fixture(body: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "modena-chemical-{}-{}",
        std::process::id(),
        SEQ.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("atoms.xyz"), body).unwrap();
    dir
}
fn deck(units: &str, extra: &str) -> String {
    format!(
        r#"units = "{units}"
[structure]
type="xyz"
file="atoms.xyz"
box_size=[24.0,24.0,24.0]
[model]
type="lennard_jones"
cutoff=8.0
[[model.types]]
sigma=2.0
epsilon=0.5
[[model.types]]
sigma=2.0
epsilon=0.5
[integrator]
type="velocity_verlet"
dt=0.00001
[run]
steps=1
initial_temperature=0.0
{extra}"#
    )
}
#[test]
fn typed_real_xyz_preserves_heterogeneous_identity_and_converts_once() {
    let dir=fixture("2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nNa 2 4 6 0 23 1\nCl 6 4 6 1 35 -1\n");
    let c = InputConfig::from_str(&deck("real", "")).unwrap();
    let s = Simulation::from_config_without_outputs(&c, &dir).unwrap();
    assert_eq!(s.system().species, vec![0, 1]);
    assert_eq!(s.system().masses, vec![23.0, 35.0]);
    assert_eq!(s.system().positions, vec![1.0, 2.0, 3.0, 3.0, 2.0, 3.0]);
    let q = (332.06371_f64 / (2.0 * 0.5)).sqrt();
    assert!((s.system().charges[0] - q).abs() < 1e-12);
    assert!((s.system().charges[1] + q).abs() < 1e-12);
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn malformed_typed_xyz_and_unknown_type_are_rejected() {
    for atom in [
        "Na 2 4 6 0 -23 1",
        "Na 2 4 6 8 23 1",
        "Na 2 4 6 0 23 NaN",
        "Na 2 4 6 0 23",
    ] {
        let dir = fixture(&format!(
            "1\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\n{atom}\n"
        ));
        let c = InputConfig::from_str(&deck("reduced", "")).unwrap();
        assert!(
            Simulation::from_config_without_outputs(&c, &dir).is_err(),
            "accepted {atom}"
        );
        std::fs::remove_dir_all(dir).unwrap();
    }
}
#[test]
fn graph_shortest_distances_and_explicit_precedence_scale_lj_and_coulomb_once() {
    let dir=fixture("4\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nA 2 2 2 0 1 1\nA 5 2 2 0 1 -1\nA 8 2 2 0 1 1\nA 11 2 2 0 1 -1\n");
    for method in ["ewald", "spme"] {
        let electro=format!("[electrostatics]\nmethod=\"{method}\"\nkappa=0.4\nreal_cutoff=4.0\nkspace_cutoff=2.0\n");
        let unscaled = InputConfig::from_str(&deck("reduced", &electro)).unwrap();
        let base = Simulation::from_config_without_outputs(&unscaled, &dir).unwrap();
        let topo = r#"[topology]
lj_scales=[0.0,0.2,0.5]
coulomb_scales=[0.1,0.4,0.8]
[[topology.bonds]]
atoms=[0,1]
r0=3.0
k=0.0
[[topology.bonds]]
atoms=[1,2]
r0=3.0
k=0.0
[[topology.bonds]]
atoms=[2,3]
r0=3.0
k=0.0
[[topology.exceptions]]
atoms=[3,0]
lj=0.7
coulomb=0.6
"#;
        let parsed = InputConfig::from_str(&deck("reduced", &(electro + topo)));
        assert!(
            parsed.is_ok(),
            "graph scale config rejected: {:?}",
            parsed.err()
        );
        let scaled = Simulation::from_config_without_outputs(&parsed.unwrap(), &dir).unwrap();
        let pot = modena_core::LennardJones::new(2.0, 0.5, 8.0).unwrap();
        use modena_core::PairPotential;
        let mut expected = pe(&base, &dir);
        for i in 0..4 {
            for j in i + 1..4 {
                let d = j - i;
                let r = 3.0 * d as f64;
                let (lj, c) = if d == 3 {
                    (0.7, 0.6)
                } else {
                    ([0.0, 0.2][d - 1], [0.1, 0.4][d - 1])
                };
                let q = if d % 2 == 0 { 1.0 } else { -1.0 };
                expected += (lj - 1.0)
                    * if r < 8.0 {
                        pot.energy_force(r * r, 0, 0).0
                    } else {
                        0.0
                    }
                    + (c - 1.0) * q / r;
            }
        }
        assert!(
            (pe(&scaled, &dir) - expected).abs() < 1e-10,
            "{method}: {} vs {expected}",
            pe(&scaled, &dir)
        );
    }
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn legacy_exclude_bonded_uses_graph_even_without_angles() {
    let dir = fixture("4\nlegacy\nA 2 2 2\nA 5 2 2\nA 8 2 2\nA 11 2 2\n");
    let topo = "[topology]\nexclude_bonded=true\n".to_string()
        + &(0..3)
            .map(|i| format!("[[topology.bonds]]\natoms=[{i},{}]\nr0=3.0\nk=0.0\n", i + 1))
            .collect::<String>();
    let c = InputConfig::from_str(&deck("reduced", &topo)).unwrap();
    let s = Simulation::from_config_without_outputs(&c, Path::new(&dir)).unwrap();
    assert!(
        pe(&s, &dir).abs() < 1e-14,
        "nonbonded graph neighbors remain {}",
        pe(&s, &dir)
    );
    std::fs::remove_dir_all(dir).unwrap();
}

fn pe(s: &Simulation, dir: &Path) -> f64 {
    let path = dir.join("probe.chk");
    s.write_checkpoint(&path).unwrap();
    modena_core::io::restart::load(path)
        .unwrap()
        .continuation
        .unwrap()
        .potential_energy
}

#[test]
fn real_coulomb_energy_and_force_are_reference_invariant() {
    // Independent physical expectations at r=4 Angstrom for charges +1/-2 e.
    let energy = -2.0 * 332.06371 / 4.0;
    let radial_force = -2.0 * 332.06371 / (4.0 * 4.0);
    for (sigma, epsilon, mass) in [(2.0, 0.5, 1.0), (3.4, 0.238, 39.948), (1.0, 1.0, 23.0)] {
        let f = modena_core::UnitFactors::new(sigma, epsilon, mass);
        let qprod = f.to_reduced_charge(1.0) * f.to_reduced_charge(-2.0);
        let r = f.to_reduced_length(4.0);
        assert!(
            (qprod / r * epsilon - energy).abs() < 1e-10,
            "sigma={sigma}: incorrect physical Coulomb energy"
        );
        assert!(
            (qprod / (r * r) * epsilon / sigma - radial_force).abs() < 1e-10,
            "sigma={sigma}: incorrect physical Coulomb force"
        );
    }
}

#[test]
fn shortest_ring_path_wins_and_exception_conflicts_fail() {
    use modena_core::io::input::{BondSpec, BondedTopologySpec, PairExceptionSpec};
    let mut t = BondedTopologySpec {
        lj_scales: Some([0.1, 0.2, 0.3]),
        coulomb_scales: Some([0.4, 0.5, 0.6]),
        ..Default::default()
    };
    for atoms in [[0, 1], [1, 2], [2, 3], [3, 0]] {
        t.bonds.push(BondSpec {
            atoms,
            r0: 1.0,
            k: 0.0,
        });
    }
    let p = modena_core::config::topology_pair_scales(&t, 4).unwrap();
    assert_eq!(
        p.iter()
            .find(|p| p.i == 0 && p.j == 3)
            .map(|p| (p.lj, p.coulomb)),
        Some((0.1, 0.4))
    );
    assert_eq!(
        p.iter()
            .find(|p| p.i == 0 && p.j == 2)
            .map(|p| (p.lj, p.coulomb)),
        Some((0.2, 0.5))
    );
    t.exceptions = vec![
        PairExceptionSpec {
            atoms: [0, 3],
            lj: 0.8,
            coulomb: 0.9,
        },
        PairExceptionSpec {
            atoms: [3, 0],
            lj: 0.7,
            coulomb: 0.9,
        },
    ];
    assert!(modena_core::config::topology_pair_scales(&t, 4).is_err());
    t.exceptions = vec![PairExceptionSpec {
        atoms: [0, 8],
        lj: 0.0,
        coulomb: 0.0,
    }];
    assert!(modena_core::config::topology_pair_scales(&t, 4).is_err());
}
#[test]
fn typed_import_rejects_global_identity_overrides_and_invalid_declarations() {
    for (comment, extra) in [
        (
            "Properties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1",
            "mass=10.0",
        ),
        ("Properties=species:S:1:pos:R:3:unknown:I:1", ""),
    ] {
        let dir = fixture(&format!("1\n{comment}\nA 2 3 4 0 1 0\n"));
        let text = deck("real", "").replace(
            "file=\"atoms.xyz\"",
            &format!("file=\"atoms.xyz\"\n{extra}"),
        );
        let c = InputConfig::from_str(&text).unwrap();
        assert!(Simulation::from_config_without_outputs(&c, &dir).is_err());
        std::fs::remove_dir_all(dir).unwrap();
    }
}
#[test]
fn typed_identity_and_graph_policy_participate_in_continuation() {
    let dir=fixture("2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nA 2 2 2 0 1 0\nB 5 2 2 1 2 0\n");
    let topo =
        "[topology]\nlj_scales=[0.5,1.0,1.0]\n[[topology.bonds]]\natoms=[0,1]\nr0=3.0\nk=0.0\n";
    let c = InputConfig::from_str(&deck("reduced", topo)).unwrap();
    let s = Simulation::from_config_without_outputs(&c, &dir).unwrap();
    s.write_checkpoint(&dir.join("continue.chk")).unwrap();
    let restored =
        deck("reduced", topo).replace("steps=1", "steps=1\nrestore_file=\"continue.chk\"");
    assert!(Simulation::from_config_without_outputs(
        &InputConfig::from_str(&restored).unwrap(),
        &dir
    )
    .is_ok());
    assert!(Simulation::from_config_without_outputs(
        &InputConfig::from_str(&restored.replace("lj_scales=[0.5", "lj_scales=[0.6")).unwrap(),
        &dir
    )
    .is_err());
    std::fs::write(dir.join("atoms.xyz"),"2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nA 2 2 2 0 1 0\nB 5 2 2 1 3 0\n").unwrap();
    assert!(Simulation::from_config_without_outputs(
        &InputConfig::from_str(&restored).unwrap(),
        &dir
    )
    .is_err());
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn real_runner_electrostatic_pair_energy_and_force_match_physical_bare_pair() {
    let dir=fixture("2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nNa 2 4 6 0 23 1\nCl 6 4 6 1 35 -2\n");
    for method in ["ewald", "spme"] {
        let electro=format!("[electrostatics]\nmethod=\"{method}\"\nkappa=0.4\nreal_cutoff=8.0\nkspace_cutoff=2.0\n");
        let bare = Simulation::from_config_without_outputs(
            &InputConfig::from_str(&deck("real", &electro)).unwrap(),
            &dir,
        )
        .unwrap();
        let excluded = Simulation::from_config_without_outputs(
            &InputConfig::from_str(&deck(
                "real",
                &(electro + "[[topology.exceptions]]\natoms=[0,1]\nlj=1.0\ncoulomb=0.0\n"),
            ))
            .unwrap(),
            &dir,
        )
        .unwrap();
        // Difference isolates exact bare-pair removal, independent of mesh/self/background terms.
        assert!(
            ((pe(&bare, &dir) - pe(&excluded, &dir)) * 0.5 - (-2.0 * 332.06371 / 4.0)).abs()
                < 1e-10
        );
        assert!(
            ((bare.system().forces[0] - excluded.system().forces[0]) * 0.5 / 2.0
                - (2.0 * 332.06371 / 16.0))
                .abs()
                < 1e-10
        );
    }
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn rust_builder_real_typed_input_uses_shared_normalization() {
    let dir=fixture("2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nA 2 4 6 0 23 0\nB 6 4 6 0 35 0\n");
    let config = Simulation::builder()
        .units(modena_core::UnitSystem::Real)
        .xyz("atoms.xyz", [24.0; 3])
        .lj(2.0, 0.5, 8.0)
        .velocity_verlet(0.0001)
        .temperature(0.0)
        .to_config()
        .unwrap();
    let sim = Simulation::from_config_without_outputs(&config, &dir).unwrap();
    assert_eq!(sim.system().masses, vec![23.0, 35.0]);
    assert_eq!(sim.system().pos(1).0, [3.0, 2.0, 3.0]);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn legacy_angle_and_dihedral_endpoint_exclusions_remain_compatible() {
    use modena_core::io::input::{AngleSpec, BondedTopologySpec};
    let t = BondedTopologySpec {
        exclude_bonded: true,
        angles: vec![AngleSpec {
            atoms: [0, 1, 2],
            theta0: 1.2,
            k: 0.0,
        }],
        ..Default::default()
    };
    let scales = modena_core::config::topology_pair_scales(&t, 3).unwrap();
    assert_eq!(
        scales
            .iter()
            .find(|p| p.i == 0 && p.j == 2)
            .map(|p| (p.lj, p.coulomb)),
        Some((0.0, 0.0))
    );
}
#[test]
fn explicit_water_exceptions_must_be_compatible_with_mandatory_exclusions() {
    let txt = deck(
        "real",
        "[[topology.exceptions]]\natoms=[0,1]\nlj=0.5\ncoulomb=0.5\n",
    )
    .replace(
        "type=\"xyz\"\nfile=\"atoms.xyz\"\nbox_size=[24.0,24.0,24.0]",
        "type=\"water\"\nn_waters=8\ndensity=0.3",
    );
    let c = InputConfig::from_str(&txt).unwrap();
    assert!(
        Simulation::from_config_without_outputs(&c, Path::new(".")).is_err(),
        "conflicting explicit water exception was silently discarded"
    );
}

#[test]
fn absurd_typed_atom_counts_are_errors_without_capacity_panics() {
    let dir = fixture(&format!(
        "{}\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\n",
        usize::MAX
    ));
    let result = std::panic::catch_unwind(|| {
        modena_core::read_xyz(dir.join("atoms.xyz"), modena_core::SimBox::cubic(10.0))
    });
    assert!(
        matches!(result, Ok(Err(_))),
        "invalid count panicked instead of returning input error"
    );
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn previous_real_charge_convention_checkpoint_cannot_continue_silently() {
    let dir=fixture("2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nNa 2 4 6 0 23 1\nCl 6 4 6 1 35 -1\n");
    let c = InputConfig::from_str(&deck("real", "")).unwrap();
    let sim = Simulation::from_config_without_outputs(&c, &dir).unwrap();
    sim.write_checkpoint(&dir.join("old.chk")).unwrap();
    let mut cp = modena_core::io::restart::load(dir.join("old.chk")).unwrap();
    let state = cp.continuation.as_mut().unwrap();
    for q in &mut state.charges {
        *q *= 2.0;
    }
    state.physical_signature["charges"] = serde_json::json!(state.charges);
    let factor = state.physical_signature["output_unit_factors"][6]
        .as_f64()
        .unwrap();
    state.physical_signature["output_unit_factors"][6] = serde_json::json!(factor * 2.0);
    modena_core::io::restart::save(dir.join("old.chk"), &cp).unwrap();
    let restored = deck("real", "").replace("steps=1", "steps=1\nrestore_file=\"old.chk\"");
    assert!(Simulation::from_config_without_outputs(
        &InputConfig::from_str(&restored).unwrap(),
        &dir
    )
    .is_err());
    std::fs::remove_dir_all(dir).unwrap();
}
