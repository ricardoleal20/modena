use std::path::{Path, PathBuf};

use modena_core::{
    ForceModel, InputConfig, LennardJones, NeighborList, PairForceModel, Simulation, Virial,
};

const BASE: &str = r#"
[structure]
type = "fcc"
lattice_constant = 2.4
replicas = [1,1,1]
mass = 1.0
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.003
[run]
steps = 5
initial_temperature = 1.0
seed = 51
log_every = 1
"#;

fn dir(label: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "modena_hardening_{label}_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}
fn sim(text: &str, dir: &Path) -> Simulation {
    Simulation::from_config(&InputConfig::from_str(text).unwrap(), dir).unwrap()
}

#[test]
fn initialization_uses_conserved_com_dof() {
    let s = sim(BASE, Path::new("."));
    assert_eq!(s.system().dof(), 9);
    assert_eq!(s.system().conserved_com_dof(), 3);
    assert!((s.system().kinetic_energy() - 4.5).abs() < 1e-12);
}

#[test]
fn initial_constraints_are_projected_before_sampling() {
    let text = format!("{BASE}\n[[constraints.bonds]]\natoms = [0,1]\nlength = 1.6\n");
    let s = sim(&text, Path::new("."));
    let sys = s.system();
    let d = sys.box_size.min_image(sys.pos(0) - sys.pos(1));
    assert!((d.norm() - 1.6).abs() < 1e-8);
    assert!(d.dot(sys.vel(0) - sys.vel(1)).abs() < 1e-8);
    assert_eq!(sys.dof(), 8);
    assert!((sys.temperature() - 1.0).abs() < 1e-12);
}

#[test]
fn langevin_initialization_retains_stochastic_com_modes() {
    let text = BASE.replace(
        "type = \"velocity_verlet\"",
        "type = \"langevin_baoab\"\ngamma = 1.0\ntarget_temperature = 1.0",
    );
    let s = sim(&text, Path::new("."));
    assert_eq!(s.system().dof(), 12);
    assert_eq!(s.system().conserved_com_dof(), 0);
    assert!(s.system().total_momentum().norm() > 1e-6);
}

#[test]
fn box_scaling_refreshes_force_energy_and_virial_caches() {
    let text = format!("{BASE}\n[barostat]\ntype = \"berendsen\"\ntarget_pressure = -2.0\ntau = 0.1\ncompressibility = 1.0\n");
    let mut s = sim(&text, Path::new("."));
    let e = s.run().unwrap().final_energies;
    let mut sys = s.system().clone();
    sys.zero_forces();
    let mut v = Virial::default();
    let m = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
    let pe = m
        .compute(&mut sys, &mut v, &mut NeighborList::new(2.5))
        .unwrap();
    assert!(
        (e.potential - pe).abs() < 1e-12,
        "cached PE={} fresh={pe}",
        e.potential
    );
    assert!((e.virial - v.w).abs() < 1e-12);
    assert_eq!(s.system().forces, sys.forces);
}

#[test]
fn strict_config_rejects_invalid_numeric_and_scientific_combinations() {
    let cases = [
        BASE.replace("dt = 0.003","dt = inf"),
        BASE.replace("dt = 0.003","dt = 0.0"),
        BASE.replace("mass = 1.0","mass = -1.0"),
        BASE.replace("initial_temperature = 1.0","initial_temperature = nan"),
        BASE.replace("lattice_constant = 2.4","lattice_constant = inf"),
        BASE.replace("cutoff = 2.5","cutoff = inf"),
        BASE.replace("type = \"velocity_verlet\"","type = \"langevin_baoab\"\ngamma = -1.0\ntarget_temperature = 1.0"),
        format!("{BASE}\n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 1.0\ntau = 0.0"),
        format!("{BASE}\n[barostat]\ntype = \"berendsen\"\ntarget_pressure = 1.0\nanisotropic = true"),
        format!("{BASE}\n[barostat]\ntype = \"parrinello_rahman\"\ntarget_pressure = 1.0\ntarget_temperature = 1.0"),
        format!("{BASE}\n[[constraints.bonds]]\natoms = [0,1]\nlength = 1.0\n[[constraints.bonds]]\natoms = [1,0]\nlength = 1.1"),
        format!("{BASE}\n[[topology.bonds]]\natoms = [0,999]\nk = 1.0\nr0 = 1.0"),
        format!("{BASE}\n[barostat]\ntype = \"berendsen\"\ntarget_pressure = 1.0\n[[constraints.bonds]]\natoms = [0,1]\nlength = 1.0"),
        format!("{}\n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 1.0", BASE.replace("initial_temperature = 1.0","initial_temperature = 0.0")),
        format!("{}\n[thermostat]\ntype = \"berendsen\"\ntarget_temperature = 1.0", BASE.replace("type = \"velocity_verlet\"","type = \"langevin_baoab\"\ngamma = 1.0\ntarget_temperature = 1.0")),
    ];
    let mut accepted = Vec::new();
    for (i, text) in cases.iter().enumerate() {
        // Parsing directly exercises the JSON/programmatic common boundary.
        let c: InputConfig = toml::from_str(text).unwrap();
        if Simulation::from_config(&c, Path::new(".")).is_ok() {
            accepted.push(i);
        }
    }
    assert!(accepted.is_empty(), "invalid cases accepted: {accepted:?}");
}

#[test]
fn unknown_keys_are_errors() {
    for text in [
        BASE.replace("dt = 0.003", "dt = 0.003\ndtt = 0.003"),
        BASE.replace("steps = 5", "steps = 5\nstepps = 5"),
    ] {
        assert!(InputConfig::from_str(&text).is_err());
    }
}

#[test]
fn real_units_toml_and_direct_config_have_identical_physics() {
    let text = format!(
        "units = \"real\"\n{}",
        BASE.replace("lattice_constant = 2.4", "lattice_constant = 8.16")
            .replace("mass = 1.0", "mass = 39.948")
            .replace("sigma = 1.0", "sigma = 3.4")
            .replace("epsilon = 1.0", "epsilon = 0.238")
            .replace("cutoff = 2.5", "cutoff = 8.5")
            .replace("initial_temperature = 1.0", "initial_temperature = 120.0")
    );
    let normalized = InputConfig::from_str(&text).unwrap();
    let raw: InputConfig = toml::from_str(&text).unwrap();
    let json: InputConfig = serde_json::from_value(serde_json::to_value(&raw).unwrap()).unwrap();
    let mut a = Simulation::from_config(&normalized, Path::new(".")).unwrap();
    let mut b = Simulation::from_config(&raw, Path::new(".")).unwrap();
    let mut c = Simulation::from_config(&json, Path::new(".")).unwrap();
    assert_eq!(a.system().positions, b.system().positions);
    assert_eq!(a.system().masses, b.system().masses);
    let e = a.run().unwrap().final_energies.total;
    assert_eq!(e, b.run().unwrap().final_energies.total);
    assert_eq!(e, c.run().unwrap().final_energies.total);
    assert_eq!(a.system().positions, b.system().positions);
}

#[test]
fn stochastic_and_extended_ensemble_split_continuation_is_exact() {
    let d = dir("seams");
    let variants=[
        BASE.to_string(),
        BASE.replace("type = \"velocity_verlet\"","type = \"langevin_baoab\"\ngamma = 0.7\ntarget_temperature = 1.0\nseed = 771"),
        format!("{BASE}\n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 1.0\ntau = 0.1\nseed = 231"),
        BASE.replace("type = \"velocity_verlet\"","type = \"nose_hoover_chain\"\ntarget_temperature = 1.0\ntau = 0.1\nchain_len = 3"),
        format!("{BASE}\n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 1.0\ntau = 0.1\nseed = 231\n[barostat]\ntype = \"monte_carlo\"\ntarget_pressure = 0.1\ntarget_temperature = 1.0\nmax_log_volume = 0.25\ninterval = 3\nseed = 93"),
    ];
    let mut failed = Vec::new();
    for (i, base) in variants.iter().enumerate() {
        let mut whole = sim(&base.replace("steps = 5", "steps = 17"), &d);
        whole.run().unwrap();
        let mut first = sim(&base.replace("steps = 5", "steps = 7"), &d);
        first.run().unwrap();
        first.write_checkpoint(&d.join("seam.chk")).unwrap();
        let resume = base.replace("steps = 5", "steps = 10\nrestore_file = \"seam.chk\"");
        let mut second = sim(&resume, &d);
        second.run().unwrap();
        if second.system().positions != whole.system().positions
            || second.system().velocities != whole.system().velocities
            || second.system().forces != whole.system().forces
        {
            failed.push(i);
        }
    }
    std::fs::remove_dir_all(d).unwrap();
    assert!(
        failed.is_empty(),
        "split continuation differs for variants: {failed:?}"
    );
}

#[test]
fn constraint_continuation_preserves_tangent_state_and_metadata() {
    let d = dir("constraint_seam");
    std::fs::write(
        d.join("atoms.xyz"),
        "3\nfixture\nX 1 1 1\nX 2.2 1 1\nX 5 5 5\n",
    )
    .unwrap();
    let base = BASE.replace(
        "type = \"fcc\"\nlattice_constant = 2.4\nreplicas = [1,1,1]",
        "type = \"xyz\"\nfile = \"atoms.xyz\"\nbox_size = [8.0,8.0,8.0]",
    );
    for integrator in [
        "type = \"velocity_verlet\"",
        "type = \"langevin_baoab\"\ngamma = 0.7\ntarget_temperature = 1.0",
        "type = \"nose_hoover_chain\"\ntarget_temperature = 1.0\ntau = 0.1",
    ] {
        let base = format!(
            "{}\n[[constraints.bonds]]\natoms = [0,1]\nlength = 1.2",
            base.replace("type = \"velocity_verlet\"", integrator)
        );
        let mut whole = sim(&base.replace("steps = 5", "steps = 25"), &d);
        let summary = whole.run().unwrap();
        assert_eq!(summary.metadata.constraint_count, 1);
        assert_eq!(
            summary.metadata.pressure_status,
            modena_core::PressureStatus::DiagnosticWithoutConstraintImpulseStress
        );
        let mut first = sim(&base.replace("steps = 5", "steps = 11"), &d);
        first.run().unwrap();
        first.write_checkpoint(&d.join("cp.chk")).unwrap();
        let mut second = sim(
            &base.replace("steps = 5", "steps = 14\nrestore_file = \"cp.chk\""),
            &d,
        );
        second.run().unwrap();
        assert_eq!(whole.system().positions, second.system().positions);
        assert_eq!(whole.system().velocities, second.system().velocities);
        assert_eq!(whole.system().forces, second.system().forces);
    }
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn intentional_protocol_switch_and_legacy_import_require_state_only() {
    let d = dir("state_only");
    let mut first = sim(BASE, &d);
    first.run().unwrap();
    first.write_checkpoint(&d.join("cp.chk")).unwrap();
    let changed = BASE
        .replace("steps = 5", "steps = 0\nrestore_file = \"cp.chk\"")
        .replace("dt = 0.003", "dt = 0.006");
    assert!(Simulation::from_config(&InputConfig::from_str(&changed).unwrap(), &d).is_err());
    let state_only = changed.replace(
        "restore_file = \"cp.chk\"",
        "restore_file = \"cp.chk\"\nrestart_mode = \"state_only\"",
    );
    let second = sim(&state_only, &d);
    assert_eq!(first.system().positions, second.system().positions);
    assert_eq!(first.system().velocities, second.system().velocities);
    let mut legacy = modena_core::io::restart::load(d.join("cp.chk")).unwrap();
    legacy.continuation = None;
    // Genuine v1 imports lack atom/ensemble state and must opt in explicitly.
    let mut text = format!(
        "MODENA_RESTART 1\natoms {}\nstep {}\ntime {}\nbox 2.4 2.4 2.4\npositions\n",
        legacy.atoms, legacy.step, legacy.time
    );
    for x in &legacy.positions {
        text.push_str(&format!("{x:.16e}\n"));
    }
    text.push_str("velocities\n");
    for x in &legacy.velocities {
        text.push_str(&format!("{x:.16e}\n"));
    }
    std::fs::write(d.join("cp.chk"), text).unwrap();
    assert!(Simulation::from_config(&InputConfig::from_str(&changed).unwrap(), &d).is_err());
    let second = sim(&state_only, &d);
    assert_eq!(first.system().positions, second.system().positions);
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn version_two_carries_atom_identity_and_rejects_corrupt_payloads() {
    let d = dir("v2");
    let s = sim(BASE, &d);
    let p = d.join("cp.chk");
    s.write_checkpoint(&p).unwrap();
    let text = std::fs::read_to_string(&p).unwrap();
    assert!(text.starts_with("MODENA_RESTART 2\n"));
    let cp = modena_core::io::restart::load(&p).unwrap();
    let state = cp.continuation.unwrap();
    assert_eq!(state.masses, s.system().masses);
    assert_eq!(state.species, s.system().species);
    assert_eq!(state.charges, s.system().charges);
    assert_eq!(state.elements, s.system().elements);
    assert_eq!(state.forces, s.system().forces);
    assert_eq!(state.conserved_com_dof, 3);
    let value: serde_json::Value = serde_json::from_str(text.split_once('\n').unwrap().1).unwrap();
    let mut length = value.clone();
    length["positions"] = serde_json::json!([0.0]);
    let mut overflow = value.clone();
    overflow["atoms"] = serde_json::json!(u64::MAX);
    let mut kind = value.clone();
    kind["continuation"]["integrator"] = serde_json::json!({"kind":"wrong"});
    for json in [length, overflow] {
        std::fs::write(&p, format!("MODENA_RESTART 2\n{json}\n")).unwrap();
        assert!(modena_core::io::restart::load(&p).is_err());
    }
    std::fs::write(&p, format!("MODENA_RESTART 2\n{kind}\n")).unwrap();
    assert!(Simulation::from_config(
        &InputConfig::from_str(&BASE.replace("steps = 5", "steps = 1\nrestore_file = \"cp.chk\""))
            .unwrap(),
        &d
    )
    .is_err());
    std::fs::write(&p, format!("{text}corruption\n")).unwrap();
    assert!(modena_core::io::restart::load(&p).is_err());
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn mc_config_requires_consistent_canonical_temperature_and_positive_proposal_controls() {
    let canonical=format!("{BASE}\n[thermostat]\ntype=\"velocity_rescaling\"\ntarget_temperature=1.0\n[barostat]\ntype=\"monte_carlo\"\ntarget_temperature=1.0\ntarget_pressure=0.1\nmax_log_volume=0.1\ninterval=3");
    assert!(InputConfig::from_str(&canonical).is_ok());
    for text in [
        canonical.replace(
            "target_temperature=1.0\ntarget_pressure",
            "target_temperature=2.0\ntarget_pressure",
        ),
        canonical.replace("interval=3", "interval=0"),
        canonical.replace("max_log_volume=0.1", "max_log_volume=-0.1"),
        canonical.replace("type=\"velocity_rescaling\"", "type=\"berendsen\""),
    ] {
        assert!(InputConfig::from_str(&text).is_err());
    }
}

#[test]
fn finite_input_that_overflows_kinetic_initialization_is_rejected() {
    let text = BASE.replace("initial_temperature = 1.0", "initial_temperature = 1e308");
    let c: InputConfig = toml::from_str(&text).unwrap();
    assert!(Simulation::from_config(&c, Path::new(".")).is_err());
}

#[test]
fn full_u64_seed_range_is_supported() {
    let text = BASE.replace("seed = 51", "seed = 18446744073709551615");
    // JSON supports the complete u64 range even though TOML integers do not.
    let mut value = serde_json::to_value(toml::from_str::<InputConfig>(BASE).unwrap()).unwrap();
    value["run"]["seed"] = serde_json::json!(u64::MAX);
    let c: InputConfig = serde_json::from_value(value).unwrap();
    assert!(
        Simulation::from_config(&c, Path::new(".")).is_ok(),
        "JSON seed must preserve u64 range; TOML test source: {text}"
    );
}

#[test]
fn continuation_preserves_com_modes_after_a_state_only_protocol_switch() {
    let d = dir("com_protocol");
    let stochastic = BASE.replace(
        "type = \"velocity_verlet\"",
        "type = \"langevin_baoab\"\ngamma = 1.0\ntarget_temperature = 1.0",
    );
    let mut first = sim(&stochastic, &d);
    first.run().unwrap();
    first.write_checkpoint(&d.join("stochastic.chk")).unwrap();
    let protocol = BASE.replace(
        "steps = 5",
        "steps = 3\nrestore_file = \"stochastic.chk\"\nrestart_mode = \"state_only\"",
    );
    let mut second = sim(&protocol, &d);
    second.run().unwrap();
    second.write_checkpoint(&d.join("protocol.chk")).unwrap();
    assert_eq!(second.system().conserved_com_dof(), 0);
    let continuation = BASE.replace("steps = 5", "steps = 0\nrestore_file = \"protocol.chk\"");
    let continued = Simulation::from_config(&InputConfig::from_str(&continuation).unwrap(), &d);
    assert!(
        continued.is_ok(),
        "continuation must retain the protocol's initialized COM modes"
    );
    assert_eq!(continued.unwrap().system().dof(), 12);
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn lattice_generator_count_overflow_is_a_config_error() {
    let text = BASE.replace("replicas = [1,1,1]", "replicas = [65536,65536,1]");
    let c: InputConfig = toml::from_str(&text).unwrap();
    let result = std::panic::catch_unwind(|| Simulation::from_config(&c, Path::new(".")));
    assert!(
        result.is_ok_and(|r| r.is_err()),
        "generator arithmetic overflow must fail before construction without a panic"
    );
}

#[test]
fn legacy_parser_rejects_nonfinite_overflow_and_trailing_corruption() {
    let d = dir("malformed");
    let path = d.join("legacy.chk");
    let valid="MODENA_RESTART 1\natoms 1\nstep 0\ntime 0\nbox 3 3 3\npositions\n0\n0\n0\nvelocities\n0\n0\n0\n";
    let cases = [
        valid.replace("time 0", "time nan"),
        valid.replace("box 3 3 3", "box -3 3 3"),
        format!("{valid}corrupt\n"),
        valid.replace("MODENA_RESTART 1", "MODENA_RESTART 1 extra"),
        valid.replace("atoms 1", "atoms 18446744073709551615"),
    ];
    let mut accepted = Vec::new();
    for (i, text) in cases.iter().enumerate() {
        std::fs::write(&path, text).unwrap();
        if std::panic::catch_unwind(|| modena_core::io::restart::load(&path))
            .is_ok_and(|r| r.is_ok())
        {
            accepted.push(i);
        }
    }
    assert!(accepted.is_empty(), "corrupt files accepted: {accepted:?}");
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn invalid_checkpoint_write_preserves_existing_file() {
    let d = dir("atomic");
    let path = d.join("cp.chk");
    let text="MODENA_RESTART 1\natoms 1\nstep 0\ntime 0\nbox 3 3 3\npositions\n0\n0\n0\nvelocities\n0\n0\n0\n";
    std::fs::write(&path, text).unwrap();
    let mut cp = modena_core::io::restart::load(&path).unwrap();
    cp.time = f64::NAN;
    assert!(modena_core::io::restart::save(&path, &cp).is_err());
    assert_eq!(std::fs::read_to_string(&path).unwrap(), text);
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn continuation_rejects_changed_physical_identity() {
    let d = dir("identity");
    let mut first = sim(BASE, &d);
    first.run().unwrap();
    first.write_checkpoint(&d.join("cp.chk")).unwrap();
    let restore = BASE.replace("steps = 5", "steps = 1\nrestore_file = \"cp.chk\"");
    for text in [
        restore.replace("dt = 0.003", "dt = 0.004"),
        restore.replace("mass = 1.0", "mass = 2.0"),
        restore.replace("epsilon = 1.0", "epsilon = 2.0"),
        restore.replace("seed = 51", "seed = 52"),
        restore.replace("mass = 1.0", "mass = 1.0\ncharge = 0.5"),
        restore.replace("mass = 1.0", "mass = 1.0\nelement = \"Ar\""),
        restore.replace(
            "type = \"velocity_verlet\"",
            "type = \"langevin_baoab\"\ngamma = 1.0\ntarget_temperature = 1.0",
        ),
        format!("{restore}\n[[constraints.bonds]]\natoms = [0,1]\nlength = 1.6"),
        format!("{restore}\n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 1.0"),
    ] {
        assert!(
            Simulation::from_config(&InputConfig::from_str(&text).unwrap(), &d).is_err(),
            "changed physical identity accepted"
        );
    }
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn continuation_binds_electrostatic_presence_method_and_parameters() {
    let d = dir("electrostatic_identity");
    let charged = BASE.replace("mass = 1.0", "mass = 1.0\ncharge = 0.5");
    let electrostatics = "\n[electrostatics]\nmethod = \"ewald\"\nkappa = 0.8\nreal_cutoff = 2.5\nkspace_cutoff = 4.0\n";
    let original = format!("{charged}{electrostatics}");
    let first = sim(&original, &d);
    first
        .write_checkpoint(&d.join("electrostatic.chk"))
        .unwrap();
    let restore = original.replace(
        "steps = 5",
        "steps = 0\nrestore_file = \"electrostatic.chk\"",
    );
    assert!(Simulation::from_config(&InputConfig::from_str(&restore).unwrap(), &d).is_ok());
    let variants = [
        restore.replace(electrostatics, ""),
        restore.replace("method = \"ewald\"", "method = \"spme\""),
        restore.replace("kappa = 0.8", "kappa = 0.9"),
        restore.replace("real_cutoff = 2.5", "real_cutoff = 2.7"),
        restore.replace("kspace_cutoff = 4.0", "kspace_cutoff = 5.0"),
    ];
    let mut accepted = Vec::new();
    for (i, text) in variants.iter().enumerate() {
        if Simulation::from_config(&InputConfig::from_str(text).unwrap(), &d).is_ok() {
            accepted.push(i);
        }
    }
    sim(&charged, &d)
        .write_checkpoint(&d.join("pair.chk"))
        .unwrap();
    let added = original.replace("steps = 5", "steps = 0\nrestore_file = \"pair.chk\"");
    if Simulation::from_config(&InputConfig::from_str(&added).unwrap(), &d).is_ok() {
        accepted.push(5);
    }
    std::fs::remove_dir_all(d).unwrap();
    assert!(
        accepted.is_empty(),
        "electrostatic identity changes accepted: {accepted:?}"
    );
}

#[test]
fn generated_water_rejects_overlapping_explicit_distance_pairs() {
    let water = r#"
units = "real"
[structure]
type = "water"
n_waters = 8
[model]
type = "lennard_jones"
cutoff = 9.0
[[model.types]]
sigma = 3.166
epsilon = 0.1553
[[model.types]]
sigma = 3.166
epsilon = 0.0
[integrator]
type = "velocity_verlet"
dt = 0.0005
[run]
steps = 0
initial_temperature = 298.15
"#;
    assert_eq!(sim(water, Path::new(".")).system().constraint_count(), 24);
    let mut not_duplicate_errors = Vec::new();
    for (i, (pair, length)) in [
        ([0, 1], 1.0),
        ([1, 0], 1.1),
        ([0, 2], 1.0),
        ([2, 0], 1.1),
        ([1, 2], 1.633),
        ([2, 1], 1.7),
    ]
    .into_iter()
    .enumerate()
    {
        let text = format!("{water}\n[[constraints.bonds]]\natoms = {pair:?}\nlength = {length}\n");
        let raw: InputConfig = toml::from_str(&text).unwrap();
        if !matches!(Simulation::from_config(&raw, Path::new(".")),
            Err(modena_core::ModenaError::Config(message)) if message.contains("duplicate constraint pair"))
        {
            not_duplicate_errors.push(i);
        }
    }
    assert!(
        not_duplicate_errors.is_empty(),
        "water overlaps did not fail as duplicate config pairs: {not_duplicate_errors:?}"
    );
}

#[test]
fn state_only_mass_changes_reassess_com_and_preserve_subsequent_continuation() {
    let d = dir("mass_protocol");
    let multicomponent = BASE
        .replace(
            "type = \"fcc\"\nlattice_constant = 2.4\nreplicas = [1,1,1]\nmass = 1.0",
            "type = \"rock_salt\"\nlattice_constant = 4.8\nreplicas = [1,1,1]\ncation = \"Na\"\nanion = \"Cl\"\ncation_mass = 1.0\nanion_mass = 1.0\ncation_charge = 0.0\nanion_charge = 0.0",
        )
        .replace(
            "sigma = 1.0\nepsilon = 1.0\ncutoff = 2.5",
            "cutoff = 2.5\n[[model.types]]\nsigma = 1.0\nepsilon = 1.0\n[[model.types]]\nsigma = 1.0\nepsilon = 1.0",
        );
    let first = sim(&multicomponent, &d);
    assert_eq!(first.system().conserved_com_dof(), 3);
    assert!(first.system().total_momentum().norm() < 1e-13);
    first.write_checkpoint(&d.join("original.chk")).unwrap();
    let import = multicomponent.replace(
        "steps = 5",
        "steps = 0\nrestore_file = \"original.chk\"\nrestart_mode = \"state_only\"",
    );
    assert_eq!(sim(&import, &d).system().conserved_com_dof(), 3);
    let changed = import.replace("cation_mass = 1.0", "cation_mass = 2.0");
    let imported = sim(&changed, &d);
    assert_eq!(imported.system().velocities, first.system().velocities);
    assert!(imported.system().total_momentum().norm() > 1e-5);
    assert_eq!(imported.system().conserved_com_dof(), 0);
    assert_eq!(imported.system().dof(), 24);

    let mut whole = sim(&changed.replace("steps = 0", "steps = 7"), &d);
    whole.run().unwrap();
    let mut split = sim(&changed.replace("steps = 0", "steps = 3"), &d);
    split.run().unwrap();
    split.write_checkpoint(&d.join("changed.chk")).unwrap();
    let continuation = multicomponent
        .replace("cation_mass = 1.0", "cation_mass = 2.0")
        .replace("steps = 5", "steps = 4\nrestore_file = \"changed.chk\"");
    let mut continued = sim(&continuation, &d);
    assert_eq!(continued.system().conserved_com_dof(), 0);
    continued.run().unwrap();
    assert_eq!(continued.system().positions, whole.system().positions);
    assert_eq!(continued.system().velocities, whole.system().velocities);
    assert_eq!(continued.system().forces, whole.system().forces);
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn mc_temperature_matching_requires_positive_controller_and_relative_small_scale_agreement() {
    let barostat = "\n[barostat]\ntype = \"monte_carlo\"\ntarget_pressure = 0.1\ntarget_temperature = 1e-13\nmax_log_volume = 0.1\ninterval = 3\n";
    let mut accepted = Vec::new();
    for (i, controller) in [
        "type = \"velocity_verlet\"\ndt = 0.003\n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 1e-13\ntau = 0.1",
        "type = \"langevin_baoab\"\ndt = 0.003\ngamma = 0.7\ntarget_temperature = 1e-13",
        "type = \"nose_hoover_chain\"\ndt = 0.003\ntau = 0.1\ntarget_temperature = 1e-13",
    ]
    .into_iter()
    .enumerate()
    {
        let text = BASE
            .replace("type = \"velocity_verlet\"\ndt = 0.003", controller)
            .replace("initial_temperature = 1.0", "initial_temperature = 1e-13");
        assert!(Simulation::from_config(&InputConfig::from_str(&format!("{text}{barostat}")).unwrap(), Path::new(".")).is_ok());
        for (j, target) in ["0.0", "2e-13"].into_iter().enumerate() {
            let mismatched = format!("{}{barostat}", text.replace("target_temperature = 1e-13", &format!("target_temperature = {target}")));
            if InputConfig::from_str(&mismatched).is_ok() {
                accepted.push((i, j));
            }
        }
    }
    assert!(
        accepted.is_empty(),
        "invalid near-zero canonical pairings accepted: {accepted:?}"
    );
}
