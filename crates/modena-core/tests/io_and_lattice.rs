//! Unit tests for lattice builders, I/O round-trips, config parsing, and
//! thermodynamic observables.

use std::io::Write;

use modena_core::{
    lattice, ForceModel, LennardJones, LogWriter, SimBox, System, TrajectoryWriter, Vec3,
};

mod lattice_tests {
    use super::*;

    #[test]
    fn fcc_atom_count_and_density() {
        // 4 atoms/cell × replicas³.
        let sys = lattice::fcc(1.5444, [3, 3, 3], "Ar");
        assert_eq!(sys.n(), 4 * 27);
        // ρ = N / V = 4 / a³.
        let a = 1.5444_f64;
        let rho = sys.n() as f64 / sys.box_size.volume();
        assert!((rho - 4.0 / a.powi(3)).abs() < 1e-9);
    }

    #[test]
    fn fcc_elements_populated() {
        let sys = lattice::fcc(1.0, [1, 1, 1], "Ar");
        assert_eq!(sys.elements.len(), 4);
        assert!(sys.elements.iter().all(|e| e == "Ar"));
    }

    #[test]
    fn fcc_nearest_neighbor_distance() {
        // FCC nearest-neighbor = a/√2. With a=√2, NN = 1.0.
        let sys = lattice::fcc(2f64.sqrt(), [2, 2, 2], "Ar");
        let p0 = sys.pos(0); // (0,0,0)
                             // Find the closest other atom.
        let mut dmin = f64::INFINITY;
        for i in 1..sys.n() {
            let d = (sys.pos(i) - p0).norm();
            if d < dmin {
                dmin = d;
            }
        }
        assert!((dmin - 1.0).abs() < 1e-9, "FCC NN = {dmin}, expected 1.0");
    }
}

mod io_tests {
    use super::*;

    #[test]
    fn xyz_round_trip() {
        let dir = tempdir();
        let path = dir.join("frame.xyz");
        // Write a 3-atom frame.
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(3);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.0, 2.0, 3.0));
        sys.set_pos(2, Vec3::new(-1.0, -2.0, -3.0));
        sys.elements = vec!["Ar".into(), "Ar".into(), "Ne".into()];

        {
            let mut w = TrajectoryWriter::create(&path).unwrap();
            w.write_frame(&sys, 0, 0.0, "test").unwrap();
        }

        // Re-read.
        let sys2 = modena_core::read_xyz(&path, SimBox::cubic(10.0)).unwrap();
        assert_eq!(sys2.n(), 3);
        assert_eq!(sys2.elements[2], "Ne");
        for i in 0..3 {
            assert!((sys2.pos(i) - sys.pos(i)).norm() < 1e-6);
        }
    }

    #[test]
    fn trajectory_writer_appends_multiple_frames() {
        let dir = tempdir();
        let path = dir.join("traj.xyz");
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(2);
        {
            let mut w = TrajectoryWriter::create(&path).unwrap();
            for step in 0..5 {
                sys.set_pos(0, Vec3::new(step as f64, 0.0, 0.0));
                w.write_frame(&sys, step, step as f64 * 0.1, "").unwrap();
            }
        }
        let text = std::fs::read_to_string(&path).unwrap();
        // 5 frames × (2 header lines + 2 atoms) = 20 lines.
        assert_eq!(text.lines().count(), 5 * 4);
        assert!(text.contains("step=0"));
        assert!(text.contains("step=4"));
    }

    #[test]
    fn log_writer_csv_format() {
        let dir = tempdir();
        let path = dir.join("energy.csv");
        {
            let mut w = LogWriter::create(&path).unwrap();
            w.log(0, 0.0, 1.0, -2.0, -1.0, 0.5, 1.5, 3.0).unwrap();
            w.log(10, 0.05, 1.1, -2.1, -1.0, 0.55, 1.6, 3.1).unwrap();
            w.flush().unwrap();
        }
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 rows
        assert!(lines[0].starts_with("step,time,kinetic"));
        let cols: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(cols.len(), 8);
        assert_eq!(cols[0], "0");
    }

    /// Minimal tempdir (avoid pulling a dev-dependency just for this).
    fn tempdir() -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "modena_test_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&p).unwrap();
        // Auto-clean on drop via a panic-free guard isn't worth the complexity
        // for tests; the OS tmpdir is periodically reaped.
        p
    }
}

mod config_tests {
    use super::*;

    #[test]
    fn reduced_config_parses() {
        let toml = r#"
units = "reduced"

[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [2, 2, 2]
element = "Ar"

[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

[integrator]
type = "velocity_verlet"
dt = 0.005

[run]
steps = 100
initial_temperature = 1.0
seed = 42
"#;
        let path = {
            let dir = std::env::temp_dir().join(format!("modena_cfg_{}.toml", nanos()));
            let mut f = std::fs::File::create(&dir).unwrap();
            f.write_all(toml.as_bytes()).unwrap();
            dir
        };
        let cfg = modena_core::InputConfig::from_path(&path).unwrap();
        assert_eq!(cfg.run.steps, 100);
        assert!(matches!(
            cfg.structure,
            modena_core::StructureSpec::Fcc { .. }
        ));
    }

    #[test]
    fn real_units_converted_to_reduced() {
        // M2: `units = "real"` is now accepted. The config is parsed and every
        // numeric field is converted to reduced units using factors derived from
        // the model's (σ, ε) and the structure mass. Here σ=ε=1 (so length and
        // energy factors are 1) but the lattice constant and cutoff are in Å and
        // must come through unchanged on the reduced side.
        let toml = r#"
units = "real"
[structure]
type = "fcc"
lattice_constant = 5.0
replicas = [1, 1, 1]
element = "Ar"
mass = 39.948
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 1.0
[run]
steps = 1
initial_temperature = 120.0
"#;
        let cfg = modena_core::InputConfig::from_str(toml).unwrap();
        // Factors were computed and stored for output conversion.
        assert!(cfg.unit_factors.is_some(), "real units must record factors");
        let f = cfg.unit_factors.unwrap();
        // σ=ε=1 → reduced == real for length/energy; T_factor = 1/kB ≈ 503.2 K.
        assert!((f.temperature - (1.0 / 0.001987204259)).abs() < 1e-3);
        // initial_temperature 120 K → reduced 120 / 503.2 ≈ 0.2385.
        assert!((cfg.run.initial_temperature - 120.0 / f.temperature).abs() < 1e-9);
        // Lattice constant 5.0 Å / σ(=1) = 5.0 reduced.
        match &cfg.structure {
            modena_core::StructureSpec::Fcc {
                lattice_constant, ..
            } => {
                assert!((lattice_constant - 5.0).abs() < 1e-12);
            }
            _ => panic!("expected fcc"),
        }
        // dt 1 ps → reduced 1 / τ. With σ=ε=1, m=39.948: τ = 0.0488882·√39.948 ≈ 0.3089 ps.
        let tau = 0.0488882 * 39.948_f64.sqrt();
        match &cfg.integrator {
            modena_core::IntegratorSpec::VelocityVerlet { dt } => {
                assert!(
                    (dt - 1.0 / tau).abs() < 1e-9,
                    "dt = {dt}, expected {}",
                    1.0 / tau
                );
            }
            _ => panic!("expected velocity_verlet"),
        }
    }

    fn nanos() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    }

    #[test]
    fn topology_and_constraints_parse() {
        // M3: a TOML config with bonded topology + SHAKE constraints.
        let toml = r#"
units = "reduced"

[structure]
type = "fcc"
lattice_constant = 5.0
replicas = [1, 1, 1]
element = "Ar"

[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

[[topology.bonds]]
atoms = [0, 1]
r0 = 1.0
k = 200.0

[[topology.angles]]
atoms = [0, 1, 2]
theta0 = 1.9
k = 50.0

[[topology.dihedrals]]
atoms = [0, 1, 2, 3]
[[topology.dihedrals.terms]]
n = 2
k = 3.0
delta = 0.0

[[constraints.bonds]]
atoms = [0, 1]
length = 1.0

[integrator]
type = "velocity_verlet"
dt = 0.005

[run]
steps = 10
"#;
        let cfg = modena_core::InputConfig::from_str(toml).unwrap();
        let topo = cfg.topology.as_ref().expect("topology present");
        assert_eq!(topo.bonds.len(), 1);
        assert_eq!(topo.bonds[0].atoms, [0, 1]);
        assert!((topo.bonds[0].r0 - 1.0).abs() < 1e-12);
        assert_eq!(topo.angles.len(), 1);
        assert_eq!(topo.angles[0].atoms, [0, 1, 2]);
        assert_eq!(topo.dihedrals.len(), 1);
        assert_eq!(topo.dihedrals[0].terms.len(), 1);
        assert_eq!(topo.dihedrals[0].terms[0].n, 2);
        let cons = cfg.constraints.as_ref().expect("constraints present");
        assert_eq!(cons.bonds.len(), 1);
        assert_eq!(cons.bonds[0].atoms, [0, 1]);
        assert!((cons.bonds[0].length - 1.0).abs() < 1e-12);
    }
}

mod thermo_tests {
    use super::*;

    #[test]
    fn kinetic_energy_and_temperature() {
        // 1 atom, m=2, v=(1,1,1) → KE = ½·2·3 = 3. dof=3 → T = 2·3/3 = 2.
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(1);
        sys.masses[0] = 2.0;
        sys.set_vel(0, Vec3::new(1.0, 1.0, 1.0));
        assert!((sys.kinetic_energy() - 3.0).abs() < 1e-12);
        assert!((sys.temperature() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn total_momentum_tracks_mass_weighted_velocity() {
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(2);
        sys.masses[0] = 1.0;
        sys.masses[1] = 3.0;
        sys.set_vel(0, Vec3::new(2.0, 0.0, 0.0));
        sys.set_vel(1, Vec3::new(1.0, 0.0, 0.0));
        let p = sys.total_momentum();
        // p_x = 1·2 + 3·1 = 5.
        assert!((p.0[0] - 5.0).abs() < 1e-12);
        assert!((p.0[1]).abs() < 1e-12);
    }

    #[test]
    fn pressure_ideal_gas_limit() {
        // Two non-interacting atoms (place far apart so LJ ≈ 0): P = 2KE/(3V).
        let mut sys = System::new(SimBox::cubic(20.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(10.0, 10.0, 10.0));
        sys.set_vel(0, Vec3::new(1.0, 0.0, 0.0));
        sys.set_vel(1, Vec3::new(-1.0, 0.0, 0.0));
        // Force model with no pairs in range → PE=0, virial=0.
        let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
        let model = modena_core::PairForceModel::new(lj);
        let mut nl = modena_core::NeighborList::new(2.5);
        let mut virial = modena_core::Virial::default();
        let pe = model.compute(&mut sys, &mut virial, &mut nl).unwrap();
        assert!(pe.abs() < 1e-9);
        assert!(virial.w.abs() < 1e-9);
        let e = modena_core::Energies::from(&sys, pe, &virial);
        let ke = sys.kinetic_energy();
        let v = sys.box_size.volume();
        let p_ideal = 2.0 * ke / (3.0 * v);
        assert!(
            (e.pressure - p_ideal).abs() < 1e-9,
            "P={}, ideal={}",
            e.pressure,
            p_ideal
        );
    }
}

mod box_tests {
    use super::*;

    #[test]
    fn wrap_into_primary_image() {
        let b = SimBox::cubic(10.0);
        let r = b.wrap(Vec3::new(12.5, -3.0, 25.0));
        assert!((r - Vec3::new(2.5, 7.0, 5.0)).norm() < 1e-12);
    }

    #[test]
    fn cubic_volume() {
        assert!((SimBox::cubic(3.0).volume() - 27.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod xyz_real_units {
    use super::*;
    use std::io::Write;

    fn nanos() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    }

    #[test]
    fn xyz_positions_are_converted_under_real_units() {
        // Regression (audit critical C2): under units = "real", an XYZ file's
        // positions were never converted to reduced units — box_size/mass/
        // charge were scaled but the coordinates come from the file and were
        // left in Å, silently scattering atoms relative to the reduced box.
        // Write an Å-scale XYZ (1 atom at (4, 0, 0) Å in a 10 Å box) with
        // σ = 2.0 (length factor 1/σ = 0.5) and require the built system's
        // position to be (2, 0, 0) reduced.
        let xyz_path = {
            let dir = std::env::temp_dir().join(format!("modena_xyz_{}.xyz", nanos()));
            let mut f = std::fs::File::create(&dir).unwrap();
            f.write_all(
                b"1\nMODENA real-units xyz regression\nAr  4.0  0.0  0.0\n",
            )
            .unwrap();
            dir
        };
        let toml = format!(
            r#"
units = "real"
[structure]
type = "xyz"
file = "{}"
box_size = [10.0, 10.0, 10.0]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 2.0
epsilon = 1.0
cutoff = 5.0
[integrator]
type = "velocity_verlet"
dt = 1.0
[run]
steps = 0
initial_temperature = 1.0
"#,
            xyz_path.display()
        );
        let cfg = modena_core::InputConfig::from_str(&toml).unwrap();
        let base = xyz_path.parent().unwrap();
        let sim = modena_core::runner::Simulation::from_config(&cfg, base).unwrap();
        let sys = sim.system();
        assert_eq!(sys.n(), 1);
        // 10 Å box / σ(=2) = 5 reduced; atom at 4 Å → 2 reduced.
        assert!(sys.box_size.size.0[0] - 5.0 < 1e-12);
        let p = sys.pos(0);
        assert!(
            (p - modena_core::Vec3::new(2.0, 0.0, 0.0)).norm() < 1e-12,
            "atom must be at the reduced (2,0,0), got {p:?}"
        );
        // And under reduced units (no factor), positions pass through untouched.
        let toml_reduced = format!(
            r#"
units = "reduced"
[structure]
type = "xyz"
file = "{}"
box_size = [10.0, 10.0, 10.0]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 2.0
epsilon = 1.0
cutoff = 5.0
[integrator]
type = "velocity_verlet"
dt = 0.005
[run]
steps = 0
initial_temperature = 1.0
"#,
            xyz_path.display()
        );
        let cfg_r = modena_core::InputConfig::from_str(&toml_reduced).unwrap();
        let sim_r = modena_core::runner::Simulation::from_config(&cfg_r, base).unwrap();
        let pr = sim_r.system().pos(0);
        assert!(
            (pr - modena_core::Vec3::new(4.0, 0.0, 0.0)).norm() < 1e-12,
            "reduced units must pass positions through, got {pr:?}"
        );
        let _ = std::fs::remove_file(&xyz_path);
    }
}
