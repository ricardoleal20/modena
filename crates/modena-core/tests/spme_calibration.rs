#[path = "../examples/spme_calibration/mod.rs"]
mod calibration;
use calibration::*;

fn protocol() -> Protocol {
    frozen_protocol("cubic", 2, 1, false).unwrap()
}
fn fixture() -> Snapshot {
    materialize(&protocol()).unwrap().remove(0)
}
fn eval(forces: Vec<[f64; 3]>) -> Evaluation {
    Evaluation {
        energy: 0.,
        strain_w: 0.,
        forces,
    }
}

#[test]
fn vector_rms_denominator_and_raw_net_force_are_explicit() {
    let r = eval(vec![[3., 4., 0.], [0., 0., 0.]]);
    let c = eval(vec![[3., 4., 0.], [0., 0., 1.]]);
    let m = metrics(&r, &c, 2, 1e-8).unwrap();
    assert!((m.reference_force_rms - 12.5_f64.sqrt()).abs() < 1e-14);
    assert!((m.force_rms_abs - 0.5_f64.sqrt()).abs() < 1e-14);
    assert!((m.force_rms_rel.unwrap() - 0.2).abs() < 1e-14);
    assert_eq!(m.max_atom_force_error, 1.);
    assert_eq!(m.net_force, [3., 4., 1.]);
}
#[test]
fn relative_zero_denominator_is_undefined() {
    let m = metrics(&eval(vec![[0.; 3]]), &eval(vec![[1., 0., 0.]]), 1, 1e-8).unwrap();
    assert_eq!(m.force_rms_rel, None);
    assert!(!passes_target(&m, 1e-3, &protocol(), true).unwrap());
}
#[test]
fn force_mean_is_not_removed() {
    let m = metrics(
        &eval(vec![[1., 0., 0.]; 2]),
        &eval(vec![[2., 0., 0.]; 2]),
        2,
        1e-8,
    )
    .unwrap();
    assert_eq!(m.force_rms_abs, 1.);
    assert_eq!(m.net_force, [4., 0., 0.]);
}
#[test]
fn all_candidate_side_gates_are_required() {
    let p = protocol();
    let base = metrics(
        &eval(vec![[1., 0., 0.]]),
        &eval(vec![[1., 0., 0.]]),
        1,
        1e-8,
    )
    .unwrap();
    assert!(passes_target(&base, 1e-3, &p, true).unwrap());
    assert!(!passes_target(&base, 1e-3, &p, false).unwrap());
    for i in 0..5 {
        let mut m = base.clone();
        match i {
            0 => m.force_rms_rel = Some(0.002),
            1 => m.force_rms_abs = 0.002,
            2 => m.max_atom_force_error = 0.02,
            3 => m.energy_per_atom_abs = 0.002,
            _ => m.w_per_atom_abs = 0.02,
        };
        assert!(!passes_target(&m, 1e-3, &p, true).unwrap());
    }
}
#[test]
fn finite_dimensions_charges_and_scales_are_validated() {
    for i in 0..8 {
        let mut s = fixture();
        match i {
            0 => s.cell[0] = f64::NAN,
            1 => s.cell[1] = 0.,
            2 => s.charges.pop().map(|_| ()).unwrap(),
            3 => s.charges.fill(0.),
            4 => s.positions[1] = s.positions[0],
            5 => s.charges[0] = f64::INFINITY,
            6 => s.positions[0][0] = f64::NAN,
            _ => s.scales.push(Scale {
                i: 0,
                j: 1,
                coulomb: 1.1,
            }),
        };
        assert!(validate_snapshot(&mut s, &protocol()).is_err(), "case {i}");
    }
    let mut s = fixture();
    s.scales = vec![
        Scale {
            i: 1,
            j: 0,
            coulomb: 0.4,
        },
        Scale {
            i: 0,
            j: 1,
            coulomb: 0.4,
        },
    ];
    validate_snapshot(&mut s, &protocol()).unwrap();
    assert_eq!(s.scales.len(), 1);
    assert_eq!(s.scales[0].i, 0);
    s.scales.push(Scale {
        i: 0,
        j: 1,
        coulomb: 0.7,
    });
    assert!(validate_snapshot(&mut s, &protocol()).is_err());
    for pair in [(0, 0), (0, 999)] {
        let mut s = fixture();
        s.scales.push(Scale {
            i: pair.0,
            j: pair.1,
            coulomb: 0.,
        });
        assert!(validate_snapshot(&mut s, &protocol()).is_err());
    }
    assert!(metrics(&eval(vec![]), &eval(vec![]), 0, 1e-8).is_err());
    assert!(metrics(
        &eval(vec![[f64::NAN, 0., 0.]]),
        &eval(vec![[0.; 3]]),
        1,
        1e-8
    )
    .is_err());
}
#[test]
fn unsupported_units_pbc_or_triclinic_rejected() {
    let mut s = fixture();
    s.units = "real".into();
    assert!(validate_snapshot(&mut s, &protocol()).is_err());
    s.units = "reduced".into();
    s.pbc[1] = false;
    assert!(validate_snapshot(&mut s, &protocol()).is_err());
    let mut v = serde_json::to_value(fixture()).unwrap();
    v["cell"] = serde_json::json!([[4.4, 0., 0.], [0., 4.4, 0.], [0., 0., 4.4]]);
    assert!(serde_json::from_value::<Snapshot>(v).is_err());
}
#[test]
fn allocation_and_term_caps_are_checked_before_work() {
    let l = protocol().limits;
    assert!(check_work(
        [usize::MAX, 2, 2],
        &WorkCounts {
            real_terms: 0,
            atom_k_terms: 0
        },
        0,
        &l
    )
    .is_err());
    assert!(check_work(
        [1024, 1024, 2],
        &WorkCounts {
            real_terms: 0,
            atom_k_terms: 0
        },
        0,
        &l
    )
    .is_err());
    assert!(check_work(
        [4; 3],
        &WorkCounts {
            real_terms: 100_000_001,
            atom_k_terms: 0
        },
        0,
        &l
    )
    .is_err());
    assert!(check_work(
        [4; 3],
        &WorkCounts {
            real_terms: 0,
            atom_k_terms: 100_000_001
        },
        0,
        &l
    )
    .is_err());
    assert!(check_work(
        [4; 3],
        &WorkCounts {
            real_terms: 0,
            atom_k_terms: 0
        },
        268_435_457,
        &l
    )
    .is_err());
    let mut p = protocol();
    p.limits.max_atoms = 63;
    assert!(materialize(&p).is_err());
    assert!(frozen_protocol("special", 3, 1, false).is_err());
    assert!(frozen_protocol("cubic", 2, 2, false).is_err());
}
#[test]
fn fixture_materialization_is_byte_stable() {
    let p = protocol();
    let a = materialize(&p).unwrap();
    let b = materialize(&p).unwrap();
    assert_eq!(
        serde_json::to_vec(&a).unwrap(),
        serde_json::to_vec(&b).unwrap()
    );
    assert_eq!(a[0].positions.len(), 64);
    assert_eq!(a[0].seed, Some(7));
    assert_eq!(a[0].amplitude, Some(0.02));
    let z = (1_664_525_u64 + 1_013_904_223 + 7 * 69_069) % 65_536;
    let expected = (0.02 * 2.2 * (2. * z as f64 / 65535. - 1.)).rem_euclid(4.4);
    assert_eq!(a[0].positions[0][0], expected);
    for s in &a {
        for r in &s.positions {
            for (d, x) in r.iter().enumerate() {
                assert!(*x >= 0. && *x < s.cell[d]);
            }
        }
    }
}
#[test]
fn holdout_ids_are_disjoint() {
    for family in ["cubic", "ortho", "special"] {
        for replicas in [2, 3] {
            if family == "special" && replicas == 3 {
                continue;
            }
            let p = frozen_protocol(family, replicas, 1, false).unwrap();
            let s = materialize(&p).unwrap();
            let selection: Vec<_> = s
                .iter()
                .filter(|s| s.role == FixtureRole::Selection)
                .map(|s| &s.id)
                .collect();
            let confirmation: Vec<_> = s
                .iter()
                .filter(|s| s.role == FixtureRole::Confirmation)
                .map(|s| &s.id)
                .collect();
            assert!(!selection.is_empty() && !confirmation.is_empty());
            assert!(selection.iter().all(|id| !confirmation.contains(id)));
            validate_frozen(&p, &s).unwrap();
        }
    }
    let p = protocol();
    let mut s = materialize(&p).unwrap();
    s[0].positions[0][0] += 0.001;
    assert!(validate_frozen(&p, &s).is_err());
}
#[test]
fn prepare_is_pure_and_rechecks_hashes_without_force_evaluator() {
    // The pure module has no evaluator callback/API. Preparation can run with no model instantiated.
    let temp = tempfile::tempdir().unwrap();
    let out = temp.path().join("frozen");
    prepare(&out, protocol()).unwrap();
    validate_prepared(&out).unwrap();
    assert!(!out.join("selected.toml").exists());
    assert!(prepare(&out, protocol()).is_err());
    std::fs::write(out.join("snapshots.json"), b"[]").unwrap();
    assert!(validate_prepared(&out).is_err());
}

#[test]
fn strict_records_reject_nonfinite_and_unknown_fields() {
    let mut s = fixture();
    s.positions[0][1] = f64::NAN;
    assert!(serde_json::to_vec(&s).is_err());
    let e = eval(vec![[f64::INFINITY, 0., 0.]]);
    assert!(serde_json::to_vec(&e).is_err());
    let mut p = protocol();
    p.targets[0] = f64::INFINITY;
    assert!(serde_json::to_vec(&p).is_err());
    assert!(toml::from_str::<ReferenceSettings>(
        "kappa = inf\nreal_cutoff = 8.0\nkspace_cutoff = 9.0"
    )
    .is_err());
    let mut value = serde_json::to_value(fixture()).unwrap();
    value["triclinic"] = serde_json::json!(true);
    assert!(serde_json::from_value::<Snapshot>(value).is_err());
}
#[test]
fn candidate_cutoff_does_not_restrict_independent_reference_images() {
    let p = Parameters {
        kappa: 0.5,
        real_cutoff: 8.8,
        kspace_cutoff: 4.,
    };
    assert!(validate_parameters(&p, [4.4; 3]).is_err());
    assert!(validate_reference_settings(&ReferenceSettings {
        kappa: p.kappa,
        real_cutoff: p.real_cutoff,
        kspace_cutoff: p.kspace_cutoff
    })
    .is_ok());
    assert!(validate_reference_settings(&ReferenceSettings {
        kappa: 0.,
        real_cutoff: 8.8,
        kspace_cutoff: 4.
    })
    .is_err());
}
#[test]
fn wrapping_background_special_policy_and_protocol_options_are_exact() {
    let p = protocol();
    let mut s = fixture();
    let expected = s.clone();
    for axis in 0..3 {
        s.positions[0][axis] += s.cell[axis];
    }
    validate_snapshot(&mut s, &p).unwrap();
    for axis in 0..3 {
        assert!((s.positions[0][axis] - expected.positions[0][axis]).abs() < 1e-14);
    }
    let all = materialize(&p).unwrap();
    for s in all.iter().filter(|s| s.id.contains("background")) {
        assert_eq!(s.charges[1], -0.5);
        assert_eq!(s.charges.iter().sum::<f64>(), 0.5);
    }
    let p = frozen_protocol("special", 2, 4, true).unwrap();
    assert_eq!(
        p.reference_ladder,
        vec![
            [4., 8.],
            [5., 10.],
            [6., 12.],
            [7., 14.],
            [8., 16.],
            [9., 18.]
        ]
    );
    let all = materialize(&p).unwrap();
    assert_eq!(all[0].positions.len(), 3);
    let scaled = all
        .iter()
        .find(|s| s.id == "special-3-selection-neutral-scaled")
        .unwrap();
    assert_eq!(
        scaled.scales,
        vec![
            Scale {
                i: 0,
                j: 1,
                coulomb: 0.
            },
            Scale {
                i: 0,
                j: 2,
                coulomb: 0.4
            }
        ]
    );
    let charged = all
        .iter()
        .find(|s| s.id == "special-3-selection-charged-bare")
        .unwrap();
    assert!((charged.charges.iter().sum::<f64>() - 0.7).abs() < 1e-15);
    let p = frozen_protocol("ortho", 3, 1, false).unwrap();
    let all = materialize(&p).unwrap();
    assert_eq!(all[0].positions.len(), 216);
    assert_eq!(all[0].cell, [2.2 * 3., 2.2 * 3. * 1.15, 2.2 * 3. * 1.3]);
}
#[test]
fn altered_caps_roles_and_frozen_protocols_are_rejected() {
    let mut p = protocol();
    p.limits.max_atoms = 257;
    assert!(materialize(&p).is_err());
    let mut p = protocol();
    p.limits.max_candidates = 125;
    assert!(materialize(&p).is_err());
    let mut p = protocol();
    p.limits.max_accuracy_calls = 1;
    assert!(materialize(&p).is_err());
    let mut p = protocol();
    p.tolerances.candidate_force_max_factor = 11.;
    assert!(materialize(&p).is_err());
    let p = protocol();
    let mut s = materialize(&p).unwrap();
    s[0].role = FixtureRole::Confirmation;
    assert!(validate_frozen(&p, &s).is_err());
}

#[test]
fn ionic_rosters_and_selection_confirmation_roles_are_fixed() {
    for family in ["cubic", "ortho"] {
        for replicas in [2, 3] {
            let p = frozen_protocol(family, replicas, 1, false).unwrap();
            let all = materialize(&p).unwrap();
            assert_eq!(all.len(), 14);
            let observed: Vec<_> = all
                .iter()
                .filter(|s| s.role != FixtureRole::RegressionOnly)
                .map(|s| {
                    (
                        s.seed.unwrap(),
                        s.amplitude.unwrap(),
                        s.role.clone(),
                        s.charges[1],
                    )
                })
                .collect();
            assert_eq!(
                observed,
                vec![
                    (7, 0.02, FixtureRole::Selection, -1.),
                    (19, 0.05, FixtureRole::Selection, -1.),
                    (31, 0.08, FixtureRole::Selection, -1.),
                    (43, 0.035, FixtureRole::Confirmation, -1.),
                    (59, 0.065, FixtureRole::Confirmation, -1.),
                    (19, 0.05, FixtureRole::Selection, -0.5),
                    (59, 0.065, FixtureRole::Confirmation, -0.5)
                ]
            );
            let crystal = all.iter().find(|s| s.id.ends_with("-crystal")).unwrap();
            assert_eq!(crystal.role, FixtureRole::RegressionOnly);
            assert_eq!(crystal.positions[0], [0.; 3]);
            assert_eq!(crystal.charges[0..2], [1., -1.]);
        }
    }
}

#[test]
fn frozen_snapshot_hashes_match_the_pre_force_protocol_record() {
    let temp = tempfile::tempdir().unwrap();
    let p = frozen_protocol("cubic", 2, 1, false).unwrap();
    let mut bytes = serde_json::to_vec_pretty(&materialize(&p).unwrap()).unwrap();
    bytes.push(b'\n');
    let path = temp.path().join("cubic-64.json");
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(
        sha256(&path).unwrap(),
        "b3c8833c77185a93d799b0f2eee83df95695e9bf007cbd1b9033b3f943ece514"
    );
    let p = frozen_protocol("cubic", 3, 1, false).unwrap();
    let mut bytes = serde_json::to_vec_pretty(&materialize(&p).unwrap()).unwrap();
    bytes.push(b'\n');
    let path = temp.path().join("cubic-216.json");
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(
        sha256(&path).unwrap(),
        "13de8ab7275ef2c317cbe73ebfac405c9030f5dacd6d3b28bf29b81cb622f3b9"
    );
    let p = frozen_protocol("ortho", 2, 1, false).unwrap();
    let mut bytes = serde_json::to_vec_pretty(&materialize(&p).unwrap()).unwrap();
    bytes.push(b'\n');
    let path = temp.path().join("ortho-64.json");
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(
        sha256(&path).unwrap(),
        "9a1f626df4aabe9338d6f31f5a41f938de24fee65826c10ea633fce9c347eefe"
    );
    let p = frozen_protocol("ortho", 3, 1, false).unwrap();
    let mut bytes = serde_json::to_vec_pretty(&materialize(&p).unwrap()).unwrap();
    bytes.push(b'\n');
    let path = temp.path().join("ortho-216.json");
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(
        sha256(&path).unwrap(),
        "9c43eab24ee2e13a035772d30f487923c66c7d8673bfaba51a1fc9ce2b5e3dc6"
    );
    let p = frozen_protocol("special", 2, 1, false).unwrap();
    let mut bytes = serde_json::to_vec_pretty(&materialize(&p).unwrap()).unwrap();
    bytes.push(b'\n');
    let path = temp.path().join("special-3.json");
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(
        sha256(&path).unwrap(),
        "43c697afef02cf27c9f22f2acb69fa3208bbbda564c6d345485928488593cc62"
    );
}

const ORACLE_A: f64 = 2.2;
fn reference_settings() -> ReferenceSettings {
    // Registered before evaluation, strictly away from the tested spherical shells.
    ReferenceSettings {
        kappa: 1.1 / ORACLE_A,
        real_cutoff: 8.73,
        kspace_cutoff: 3.75,
    }
}
fn special_cases() -> Vec<Snapshot> {
    materialize(&frozen_protocol("special", 2, 1, false).unwrap())
        .unwrap()
        .into_iter()
        .filter(|s| {
            s.role == FixtureRole::Selection
                && (s.id.ends_with("-bare") || s.id.ends_with("-scaled"))
        })
        .collect()
}
fn oracle(s: &Snapshot, p: &ReferenceSettings) -> Evaluation {
    reference(s, p.clone(), &protocol().limits).unwrap()
}
fn close(x: f64, y: f64, tol: f64, label: &str) {
    assert!(
        (x - y).abs() <= tol,
        "{label}: {x:.17e} versus {y:.17e}, residual {:.9e}, cap {tol:.9e}",
        (x - y).abs()
    );
}
fn evaluations_close(a: &Evaluation, b: &Evaluation, tol: f64) {
    close(a.energy, b.energy, tol / ORACLE_A, "energy");
    close(a.strain_w, b.strain_w, tol / ORACLE_A, "W");
    for (fa, fb) in a.forces.iter().zip(&b.forces) {
        for axis in 0..3 {
            close(fa[axis], fb[axis], tol / ORACLE_A.powi(2), "force");
        }
    }
}
// Independent exact tuple sets, only for these small controlled derivative fixtures.
// These detect a cutoff crossing even when a tiny term would hide its energy jump.
type CutoffMembers = (
    std::collections::BTreeSet<(usize, usize, [i32; 3])>,
    std::collections::BTreeSet<[i32; 3]>,
);
fn cutoff_members(s: &Snapshot, p: &ReferenceSettings) -> CutoffMembers {
    let mut real = std::collections::BTreeSet::new();
    let mut reciprocal = std::collections::BTreeSet::new();
    let b = s.cell.map(|l| (p.real_cutoff / l).ceil() as i32 + 1);
    for i in 0..s.positions.len() {
        for j in 0..s.positions.len() {
            for x in -b[0]..=b[0] {
                for y in -b[1]..=b[1] {
                    for z in -b[2]..=b[2] {
                        let n = [x, y, z];
                        if i == j && n == [0; 3] {
                            continue;
                        }
                        let r2: f64 = (0..3)
                            .map(|a| {
                                (s.positions[i][a] - s.positions[j][a]
                                    + f64::from(n[a]) * s.cell[a])
                                    .powi(2)
                            })
                            .sum();
                        if r2 < p.real_cutoff.powi(2) {
                            real.insert((i, j, n));
                        }
                    }
                }
            }
        }
    }
    let b = s
        .cell
        .map(|l| (p.kspace_cutoff * l / (2. * std::f64::consts::PI)).ceil() as i32);
    for x in -b[0]..=b[0] {
        for y in -b[1]..=b[1] {
            for z in -b[2]..=b[2] {
                let m = [x, y, z];
                let k2: f64 = (0..3)
                    .map(|a| (2. * std::f64::consts::PI * f64::from(m[a]) / s.cell[a]).powi(2))
                    .sum();
                if k2 > 0. && k2 < p.kspace_cutoff.powi(2) {
                    reciprocal.insert(m);
                }
            }
        }
    }
    (real, reciprocal)
}
fn derivative_oracle(
    s: &Snapshot,
    p: &ReferenceSettings,
    expected: &CutoffMembers,
) -> (Evaluation, ReferenceAudit) {
    // Exact independent geometry gate runs before the energy-bearing oracle call.
    assert_eq!(
        expected,
        &cutoff_members(s, p),
        "derivative cutoff membership changed"
    );
    let mut budget = ReferenceBudget::new(&protocol().limits);
    let (evaluation, audit) =
        reference_with_budget(s, p.clone(), &protocol().limits, &mut budget).unwrap();
    assert_eq!(audit.real_included_terms, expected.0.len() as u64);
    assert_eq!(audit.reciprocal_vectors, expected.1.len() as u64);
    (evaluation, audit)
}
fn same_trace(base: &ReferenceAudit, other: &ReferenceAudit) {
    assert_eq!(base.real_membership_hash, other.real_membership_hash);
    assert_eq!(
        base.reciprocal_membership_hash,
        other.reciprocal_membership_hash
    );
}
#[test]
fn reference_force_is_negative_energy_gradient() {
    let p = reference_settings();
    for s in special_cases() {
        let members = cutoff_members(&s, &p);
        println!(
            "FORCE_MEMBERSHIP_BASE {} real={:?} reciprocal={:?}",
            s.id, members.0, members.1
        );
        let (e, base_trace) = derivative_oracle(&s, &p, &members);
        for h in [1e-5 * ORACLE_A, 0.5e-5 * ORACLE_A] {
            for atom in 0..3 {
                for axis in 0..3 {
                    let mut plus = s.clone();
                    plus.positions[atom][axis] += h;
                    let mut minus = s.clone();
                    minus.positions[atom][axis] -= h;
                    assert_eq!(members, cutoff_members(&plus, &p), "plus cutoff crossing");
                    assert_eq!(members, cutoff_members(&minus, &p), "minus cutoff crossing");
                    let (ep, plus_trace) = derivative_oracle(&plus, &p, &members);
                    let (em, minus_trace) = derivative_oracle(&minus, &p, &members);
                    same_trace(&base_trace, &plus_trace);
                    same_trace(&base_trace, &minus_trace);
                    println!("FORCE_TRACE {} atom={atom} axis={axis} h={h:.17e} plus={plus_trace:?} minus={minus_trace:?}",s.id);
                    let fd = -(ep.energy - em.energy) / (2. * h);
                    let analytic = e.forces[atom][axis];
                    println!("FORCE_FD {} h={h:.17e} atom={atom} axis={axis} analytic={analytic:.17e} fd={fd:.17e} residual={:.17e} real_members={} k_members={}",s.id,(analytic-fd).abs(),members.0.len(),members.1.len());
                    close(
                        analytic,
                        fd,
                        1e-7 / ORACLE_A.powi(2) + 1e-6 * analytic.abs().max(fd.abs()),
                        "force FD",
                    );
                }
            }
        }
    }
}
#[test]
fn reference_w_is_isotropic_strain_derivative() {
    let p = reference_settings();
    for s in special_cases() {
        let members = cutoff_members(&s, &p);
        println!(
            "STRAIN_MEMBERSHIP_BASE {} real={:?} reciprocal={:?}",
            s.id, members.0, members.1
        );
        let (e, base_trace) = derivative_oracle(&s, &p, &members);
        for h in [1e-5, 0.5e-5] {
            let strain = |lambda: f64| {
                let mut t = s.clone();
                t.cell = t.cell.map(|x| x * lambda);
                for r in &mut t.positions {
                    *r = r.map(|x| x * lambda);
                }
                t
            };
            let plus = strain(1. + h);
            let minus = strain(1. - h);
            assert_eq!(
                members,
                cutoff_members(&plus, &p),
                "strain plus cutoff crossing"
            );
            assert_eq!(
                members,
                cutoff_members(&minus, &p),
                "strain minus cutoff crossing"
            );
            let (ep, plus_trace) = derivative_oracle(&plus, &p, &members);
            let (em, minus_trace) = derivative_oracle(&minus, &p, &members);
            same_trace(&base_trace, &plus_trace);
            same_trace(&base_trace, &minus_trace);
            println!(
                "STRAIN_TRACE {} h={h:.17e} plus={plus_trace:?} minus={minus_trace:?}",
                s.id
            );
            let fd = -(ep.energy - em.energy) / (2. * h);
            println!("STRAIN_FD {} h={h:.17e} analytic={:.17e} fd={fd:.17e} residual={:.17e} real_members={} k_members={}",s.id,e.strain_w,(e.strain_w-fd).abs(),members.0.len(),members.1.len());
            close(
                e.strain_w,
                fd,
                1e-7 / ORACLE_A + 1e-6 * e.strain_w.abs().max(fd.abs()),
                "W FD",
            );
        }
    }
}
// Half-space reciprocal sum uses pair differences, never the reference rho or k loop.
fn independent_half_reciprocal(s: &Snapshot, p: &ReferenceSettings) -> Evaluation {
    let pi = std::f64::consts::PI;
    let c = 4. * pi / s.cell.iter().product::<f64>();
    let mut out = eval(vec![[0.; 3]; s.positions.len()]);
    let b = s
        .cell
        .map(|l| (p.kspace_cutoff * l / (2. * pi)).ceil() as i32);
    for x in 0..=b[0] {
        for y in -b[1]..=b[1] {
            for z in -b[2]..=b[2] {
                if !(x > 0 || x == 0 && y > 0 || x == 0 && y == 0 && z > 0) {
                    continue;
                }
                let k = [
                    2. * pi * f64::from(x) / s.cell[0],
                    2. * pi * f64::from(y) / s.cell[1],
                    2. * pi * f64::from(z) / s.cell[2],
                ];
                let k2 = k.iter().map(|v| v * v).sum::<f64>();
                if k2 >= p.kspace_cutoff.powi(2) {
                    continue;
                }
                let pref = c * (-k2 / (4. * p.kappa.powi(2))).exp() / k2;
                let mut density = s.charges.iter().map(|q| q * q).sum::<f64>();
                for i in 0..s.positions.len() {
                    for j in 0..i {
                        let phase = (0..3)
                            .map(|a| k[a] * (s.positions[i][a] - s.positions[j][a]))
                            .sum::<f64>();
                        density += 2. * s.charges[i] * s.charges[j] * phase.cos();
                        let f = 2. * pref * s.charges[i] * s.charges[j] * phase.sin();
                        for (a, ka) in k.iter().enumerate() {
                            out.forces[i][a] += f * ka;
                            out.forces[j][a] -= f * ka;
                        }
                    }
                }
                out.energy += pref * density;
                out.strain_w += pref * density * (1. - k2 / (2. * p.kappa.powi(2)));
            }
        }
    }
    out
}
#[test]
fn reference_full_and_half_k_oracles_agree() {
    for s in special_cases() {
        let p = reference_settings();
        let full = oracle(&s, &p);
        // K=.1 contains no modes: this independently isolates real/self/background/correction.
        let mut no_k = p.clone();
        no_k.kspace_cutoff = 0.1;
        let base = oracle(&s, &no_k);
        let mut half = independent_half_reciprocal(&s, &p);
        half.energy += base.energy;
        half.strain_w += base.strain_w;
        for (f, b) in half.forces.iter_mut().zip(base.forces) {
            for a in 0..3 {
                f[a] += b[a];
            }
        }
        evaluations_close(&full, &half, 1e-10);
        println!(
            "HALF_K {} energy_delta={:.17e} w_delta={:.17e}",
            s.id,
            full.energy - half.energy,
            full.strain_w - half.strain_w
        );
    }
}
#[test]
fn reference_self_image_energy_is_included() {
    let mut s = special_cases().remove(0);
    s.cell = [ORACLE_A; 3];
    s.positions = vec![[0.4, 0.5, 0.6]];
    s.charges = vec![1.];
    s.scales.clear();
    let p = ReferenceSettings {
        kappa: 0.5,
        real_cutoff: 2.42,
        kspace_cutoff: 0.1,
    };
    let e = oracle(&s, &p);
    let pi = std::f64::consts::PI;
    let u_image = 3. * libm::erfc(p.kappa * ORACLE_A) / ORACLE_A;
    let w_image = 3.
        * (libm::erfc(p.kappa * ORACLE_A) / ORACLE_A
            + 2. * p.kappa / pi.sqrt() * (-(p.kappa * ORACLE_A).powi(2)).exp());
    let bg = -pi / (2. * ORACLE_A.powi(3) * p.kappa.powi(2));
    close(
        e.energy,
        u_image - p.kappa / pi.sqrt() + bg,
        1e-13 / ORACLE_A,
        "six periodic self energies",
    );
    close(
        e.strain_w,
        w_image + 3. * bg,
        1e-13 / ORACLE_A,
        "six periodic self W",
    );
    for f in e.forces[0] {
        close(
            f,
            0.,
            1e-13 / ORACLE_A.powi(2),
            "self-image force cancellation",
        );
    }
}
#[test]
fn reference_charged_background_w_is_three_u() {
    for s in special_cases() {
        let p = reference_settings();
        let mut budget = ReferenceBudget::new(&protocol().limits);
        let (_, audit) =
            reference_with_budget(&s, p.clone(), &protocol().limits, &mut budget).unwrap();
        let q = s.charges.iter().sum::<f64>();
        let expected = -std::f64::consts::PI * q * q
            / (2. * s.cell.iter().product::<f64>() * p.kappa * p.kappa);
        close(
            audit.background_energy,
            expected,
            1e-13 / ORACLE_A,
            "background coefficient",
        );
        close(
            audit.background_w,
            3. * expected,
            1e-13 / ORACLE_A,
            "background W",
        );
        println!(
            "BACKGROUND {} Q={q:.17e} U={:.17e} W={:.17e} audit={audit:?}",
            s.id, audit.background_energy, audit.background_w
        );
    }
}
#[test]
fn reference_bare_pair_removal_inside_and_outside_cutoff() {
    let p = ReferenceSettings {
        kappa: 0.5,
        real_cutoff: 3.7,
        kspace_cutoff: 3.71,
    };
    for s in special_cases().into_iter().filter(|s| s.scales.is_empty()) {
        let base = oracle(&s, &p);
        for j in [1, 2] {
            for c in [0., 0.4, 0.7, 1.] {
                let mut scaled = s.clone();
                scaled.scales = vec![Scale {
                    i: 0,
                    j,
                    coulomb: c,
                }];
                let got = oracle(&scaled, &p);
                let d: [f64; 3] = std::array::from_fn(|a| {
                    let v = s.positions[0][a] - s.positions[j][a];
                    v - s.cell[a] * (v / s.cell[a]).round()
                });
                let r = d.iter().map(|x| x * x).sum::<f64>().sqrt();
                assert_eq!(r < p.real_cutoff, j == 1);
                let removed = (1. - c) * s.charges[0] * s.charges[j] / r;
                let mut expected = base.clone();
                expected.energy -= removed;
                expected.strain_w -= removed;
                for (a, da) in d.iter().enumerate() {
                    let f = removed * da / (r * r);
                    expected.forces[0][a] -= f;
                    expected.forces[j][a] += f;
                }
                evaluations_close(&got, &expected, 1e-11);
            }
        }
    }
}
#[test]
fn reference_rigid_translation_invariant() {
    let p = reference_settings();
    for s in special_cases() {
        let base = oracle(&s, &p);
        for shift in protocol().translations {
            let mut moved = s.clone();
            for r in &mut moved.positions {
                for a in 0..3 {
                    r[a] += shift[a] * s.cell[a] / 16.;
                }
            }
            evaluations_close(&base, &oracle(&moved, &p), 1e-9);
        }
    }
}
#[test]
fn reference_full_box_wrap_invariant() {
    let p = reference_settings();
    for s in special_cases() {
        let base = oracle(&s, &p);
        for atom in 0..3 {
            for a in 0..3 {
                for sign in [-1., 1.] {
                    let mut moved = s.clone();
                    moved.positions[atom][a] += sign * s.cell[a];
                    // Raw noncanonical coordinates reach the oracle, not materialize/validate_snapshot.
                    evaluations_close(&base, &oracle(&moved, &p), 1e-9);
                }
            }
        }
    }
}

#[test]
fn reference_prechecks_bounds_heap_work_and_retained_deadline() {
    let s = special_cases().remove(0);
    let p = reference_settings();
    let limits = protocol().limits;
    let mut budget = ReferenceBudget::new(&limits);
    let (_, audit) = reference_with_budget(&s, p.clone(), &limits, &mut budget).unwrap();
    assert_eq!(audit.work.real_terms, 1122); // 3*3*5^3 - 3 zero-self tuples.
    assert_eq!(audit.prechecked_work.real_terms, audit.work.real_terms);
    assert_eq!(audit.work.atom_k_terms, 6 * audit.reciprocal_vectors); // two actual N-atom passes.
    assert_eq!(
        audit.prechecked_work.atom_k_terms,
        6 * audit.reciprocal_enumerated_vectors
    );
    assert!(
        audit.deadline_checks
            >= (audit.work.real_terms
                + audit.work.atom_k_terms
                + audit.reciprocal_enumerated_vectors)
                / 1024
    );
    assert_eq!(budget.work, audit.work);
    let mut reduced = limits.clone();
    reduced.max_real_terms_per_family = audit.prechecked_work.real_terms;
    let mut budget = ReferenceBudget::new(&reduced);
    reference_with_budget(&s, p.clone(), &reduced, &mut budget).unwrap();
    let old = budget.work.clone();
    assert!(matches!(
        reference_with_budget(&s, p.clone(), &reduced, &mut budget),
        Err(CalibrationFailure::BudgetExceeded(_))
    ));
    assert_eq!(budget.work, old);
    assert_eq!(budget.calls, 1);
    for kind in 0..6 {
        let mut reduced = limits.clone();
        match kind {
            0 => reduced.max_real_terms_per_call = audit.prechecked_work.real_terms - 1,
            1 => reduced.max_atom_k_terms_per_call = audit.prechecked_work.atom_k_terms - 1,
            2 => reduced.max_owned_heap_bytes = audit.owned_heap_bytes - 1,
            3 => reduced.max_real_terms_per_family = audit.prechecked_work.real_terms - 1,
            4 => reduced.max_atom_k_terms_per_family = audit.prechecked_work.atom_k_terms - 1,
            _ => reduced.max_family_wall_seconds = 0,
        }
        let mut b = ReferenceBudget::new(&reduced);
        assert!(
            matches!(
                reference_with_budget(&s, p.clone(), &reduced, &mut b),
                Err(CalibrationFailure::BudgetExceeded(_))
            ),
            "cap {kind}"
        );
        assert_eq!(b.calls, 0);
        assert_eq!(
            b.work,
            WorkCounts {
                real_terms: 0,
                atom_k_terms: 0
            }
        );
    }
    assert!(matches!(
        ReferenceBudget::resume(
            &limits,
            audit.work.clone(),
            1,
            std::time::Duration::from_secs(limits.max_family_wall_seconds),
            0
        ),
        Err(CalibrationFailure::BudgetExceeded(_))
    ));
    let mut resumed = ReferenceBudget::resume(
        &limits,
        audit.work.clone(),
        1,
        std::time::Duration::ZERO,
        32,
    )
    .unwrap();
    let (_, next) = reference_with_budget(&s, p.clone(), &limits, &mut resumed).unwrap();
    assert_eq!(resumed.work.real_terms, 2 * audit.work.real_terms);
    assert_eq!(resumed.work.atom_k_terms, 2 * audit.work.atom_k_terms);
    assert_eq!(next.owned_heap_bytes, audit.owned_heap_bytes + 32);
    for (r, k) in [(f64::MAX, 1.), (1., f64::MAX)] {
        let overflow = ReferenceSettings {
            kappa: 0.5,
            real_cutoff: r,
            kspace_cutoff: k,
        };
        assert!(matches!(
            reference(&s, overflow, &limits),
            Err(CalibrationFailure::BudgetExceeded(_))
        ));
    }
    println!(
        "BUDGET actual={:?} upper={:?} heap={} checks={} resumed={:?}",
        audit.work,
        audit.prechecked_work,
        audit.owned_heap_bytes,
        audit.deadline_checks,
        resumed.work
    );
}
#[test]
fn reference_rejects_invalid_physics_and_nonfinite_arithmetic() {
    let base = special_cases().remove(0);
    let p = reference_settings();
    let limits = protocol().limits;
    for kind in 0..8 {
        let mut s = base.clone();
        match kind {
            0 => s.positions[1] = s.positions[0],
            1 => s.positions[0][0] = f64::NAN,
            2 => s.cell[0] = 0.,
            3 => s.charges.fill(0.),
            4 => s.pbc[0] = false,
            5 => s.units = "real".into(),
            6 => {
                s.scales = vec![
                    Scale {
                        i: 0,
                        j: 1,
                        coulomb: 0.4
                    };
                    2
                ]
            }
            _ => {
                s.scales = vec![Scale {
                    i: 1,
                    j: 0,
                    coulomb: 0.4,
                }]
            }
        }
        assert!(
            reference(&s, p.clone(), &limits).is_err(),
            "invalid physics {kind}"
        );
    }
    let mut s = base.clone();
    s.charges[0] = f64::MAX;
    assert!(matches!(
        reference(&s, p, &limits),
        Err(CalibrationFailure::ReferencePrecisionLimited(_))
    ));
}
fn reference_delta_score(
    left: &Evaluation,
    right: &Evaluation,
    final_value: &Evaluation,
    label: &str,
) -> f64 {
    let m = metrics(right, left, left.forces.len(), 1e-8 / ORACLE_A.powi(2)).unwrap();
    let denominator = metrics(
        final_value,
        final_value,
        final_value.forces.len(),
        1e-8 / ORACLE_A.powi(2),
    )
    .unwrap()
    .reference_force_rms;
    let absolute = 1e-8 / ORACLE_A.powi(2);
    let rms_budget = if denominator > absolute {
        absolute.min(0.01 * 1e-4 * denominator)
    } else {
        absolute
    };
    let score = (m.force_rms_abs / rms_budget)
        .max(m.max_atom_force_error / (1e-7 / ORACLE_A.powi(2)))
        .max(m.energy_per_atom_abs / (1e-9 / ORACLE_A))
        .max(m.w_per_atom_abs / (1e-9 / ORACLE_A));
    println!("STRICT_REFERENCE {label} final_D={denominator:.17e} rms_budget={rms_budget:.17e} score={score:.17e} metrics={m:?}");
    assert!(score.is_finite());
    score
}
fn strict_reference_delta(final_value: &Evaluation, value: &Evaluation, label: &str) {
    assert!(
        reference_delta_score(value, final_value, final_value, label) <= 1.,
        "{label}: reference budgets"
    );
}
fn crystal_levels(p: &Protocol, s: &Snapshot, budget: &mut ReferenceBudget) -> Vec<Evaluation> {
    let mut levels = Vec::new();
    for [r, k] in &p.reference_ladder {
        let settings = ReferenceSettings {
            kappa: p.reference_kappa_a / ORACLE_A,
            real_cutoff: r * ORACLE_A,
            kspace_cutoff: k / ORACLE_A,
        };
        let (value, audit) = reference_with_budget(s, settings.clone(), &p.limits, budget).unwrap();
        println!("MADELUNG_LEVEL replicas={} extended={} settings={settings:?} energy={:.17e} W={:.17e} audit={audit:?}",p.replicas,p.reference_ladder.len()==6,value.energy,value.strain_w);
        levels.push(value);
    }
    levels
}
#[test]
fn reference_default_crystal_remains_unconverged() {
    let p = frozen_protocol("cubic", 2, 1, false).unwrap();
    let s = materialize(&p)
        .unwrap()
        .into_iter()
        .find(|s| s.id.ends_with("-crystal"))
        .unwrap();
    let mut budget = ReferenceBudget::new(&p.limits);
    let outcome = converged_reference_with_budget(&s, &p, &mut budget).unwrap();
    assert_eq!(outcome.record.status, Status::ReferenceUnconverged);
    assert_eq!(outcome.record.chosen_level, None);
    assert_eq!(outcome.record.levels.len(), 4);
    assert_eq!(outcome.record.crosschecks.len(), 2);
    assert!(outcome.refinement_scores[1] > 1. && outcome.refinement_scores[2] <= 1.);
    println!(
        "MADELUNG_DEFAULT {}",
        serde_json::to_string(&outcome).unwrap()
    );
}
#[test]
fn reference_madelung_extended_refinement_sign_and_normalization() {
    let mut per_pair_values = Vec::new();
    for replicas in [2, 3] {
        // Separately registered, finite preset; never automatic growth after a failed test.
        let p = frozen_protocol("cubic", replicas, 1, true).unwrap();
        println!(
            "MADELUNG_EXTENDED_PROTOCOL {}",
            serde_json::to_string(&p).unwrap()
        );
        let s = materialize(&p)
            .unwrap()
            .into_iter()
            .find(|s| s.id.ends_with("-crystal"))
            .unwrap();
        let mut budget = ReferenceBudget::new(&p.limits);
        let levels = crystal_levels(&p, &s, &mut budget);
        let final_value = levels.last().unwrap();
        let scores: Vec<_> = levels
            .windows(2)
            .enumerate()
            .map(|(i, pair)| {
                reference_delta_score(
                    &pair[0],
                    &pair[1],
                    final_value,
                    &format!("extended crystal replicas={replicas} refinement={i}"),
                )
            })
            .collect();
        assert!(
            scores[3] <= 1. && scores[4] <= 1.,
            "extended final two refinements"
        );
        for ka in p.crosscheck_kappa_a {
            let settings = ReferenceSettings {
                kappa: ka / ORACLE_A,
                real_cutoff: 9. * ORACLE_A,
                kspace_cutoff: 18. / ORACLE_A,
            };
            let (cross, audit) =
                reference_with_budget(&s, settings.clone(), &p.limits, &mut budget).unwrap();
            println!("MADELUNG_CROSS replicas={replicas} extended=true settings={settings:?} audit={audit:?}");
            strict_reference_delta(final_value, &cross, "extended Madelung cross kappa");
        }
        let per_pair = final_value.energy / (s.positions.len() as f64 / 2.);
        close(
            per_pair,
            -1.747558775 / (ORACLE_A / 2.),
            2e-5 / ORACLE_A,
            "extended Madelung sign/normalization sanity only",
        );
        println!("MADELUNG_EXTENDED_SANITY replicas={replicas} energy_per_pair={per_pair:.17e} total_work={:?}",budget.work);
        per_pair_values.push(per_pair);
    }
    close(
        per_pair_values[0],
        per_pair_values[1],
        1e-9 / ORACLE_A,
        "extended finite-cell crystal replication",
    );
}
#[test]
fn reference_secondary_core_ewald_only_after_image_and_tail_checks() {
    use modena_core::{
        EwaldForceModel, ForceModel, NeighborList, PairScale, SimBox, System, Vec3, Virial,
    };
    for s in special_cases() {
        let limits = protocol().limits;
        let mut budget = ReferenceBudget::new(&limits);
        let mut levels = Vec::new();
        for (r, k) in [(8.8, 22.), (11., 24.), (13.2, 26.)] {
            let settings = ReferenceSettings {
                kappa: 2.,
                real_cutoff: r,
                kspace_cutoff: k,
            };
            let (value, audit) =
                reference_with_budget(&s, settings.clone(), &limits, &mut budget).unwrap();
            println!(
                "SECONDARY_LEVEL {} settings={settings:?} energy={:.17e} W={:.17e} audit={audit:?}",
                s.id, value.energy, value.strain_w
            );
            levels.push(value);
        }
        let converged = &levels[2];
        strict_reference_delta(converged, &levels[1], "secondary final refinement");
        // Both consecutive differences use the final denominator.
        let mut previous = converged.clone();
        previous.energy += levels[0].energy - levels[1].energy;
        previous.strain_w += levels[0].strain_w - levels[1].strain_w;
        for (i, f) in previous.forces.iter_mut().enumerate() {
            for (a, component) in f.iter_mut().enumerate() {
                *component += levels[0].forces[i][a] - levels[1].forces[i][a];
            }
        }
        strict_reference_delta(converged, &previous, "secondary earlier refinement");
        for kappa in [1.8, 2.2] {
            let settings = ReferenceSettings {
                kappa,
                real_cutoff: 13.2,
                kspace_cutoff: 26.,
            };
            let (cross, audit) = reference_with_budget(&s, settings, &limits, &mut budget).unwrap();
            println!("SECONDARY_CROSS {} kappa={kappa} audit={audit:?}", s.id);
            strict_reference_delta(converged, &cross, "secondary independent kappa");
        }
        let short = ReferenceSettings {
            kappa: 2.,
            real_cutoff: 4.41,
            kspace_cutoff: 26.,
        };
        let (omitted, audit) =
            reference_with_budget(&s, short.clone(), &limits, &mut budget).unwrap();
        println!("SECONDARY_OMITTED {} audit={audit:?}", s.id);
        strict_reference_delta(
            converged,
            &omitted,
            "secondary omitted images and real tail",
        );
        let mut system = System::new(SimBox::from_array(s.cell));
        system.resize(s.positions.len());
        for (i, r) in s.positions.iter().enumerate() {
            system.set_pos(i, Vec3(*r));
        }
        system.charges = s.charges.clone();
        let scales: Vec<_> = s
            .scales
            .iter()
            .map(|c| PairScale {
                i: c.i,
                j: c.j,
                lj: 1.,
                coulomb: c.coulomb,
            })
            .collect();
        let model = EwaldForceModel::new(short.kappa, short.real_cutoff, short.kspace_cutoff)
            .unwrap()
            .with_pair_scales(scales.clone())
            .unwrap();
        let mut neighbors =
            NeighborList::with_pair_scales(short.real_cutoff, 0.1, &scales).unwrap();
        let mut virial = Virial::default();
        system.zero_forces();
        let energy = model
            .compute(&mut system, &mut virial, &mut neighbors)
            .unwrap();
        let core = Evaluation {
            energy,
            strain_w: virial.w,
            forces: system
                .forces
                .chunks_exact(3)
                .map(|f| [f[0], f[1], f[2]])
                .collect(),
        };
        strict_reference_delta(converged, &core, &format!("secondary core {}", s.id));
    }
}

// These synthetic evaluations test convergence control, never oracle physics.
fn ladder_fixture(extended: bool) -> (Protocol, Snapshot) {
    let p = frozen_protocol("special", 2, 1, extended).unwrap();
    let s = materialize(&p)
        .unwrap()
        .into_iter()
        .find(|s| s.id == "special-3-selection-neutral-bare")
        .unwrap();
    (p, s)
}
fn synthetic_ladder(p: &Protocol, s: &Snapshot, energies: &[f64]) -> ReferenceOutcome {
    let mut index = 0;
    converged_reference_test(s, p, &mut ReferenceBudget::new(&p.limits), |_, _| {
        let mut value = eval(vec![[1., 0., 0.]; s.positions.len()]);
        value.energy = energies[index];
        index += 1;
        Ok(value)
    })
    .unwrap()
}
#[test]
fn two_consecutive_reference_passes_required() {
    let (p, s) = ladder_fixture(false);
    let unit = 1e-9 / p.a * 3.;
    let failed = synthetic_ladder(&p, &s, &[10. * unit, 2. * unit, 0., 0., 0., 0.]);
    assert_eq!(failed.record.status, Status::ReferenceUnconverged);
    assert_eq!(failed.record.levels.len(), 4);
    assert_eq!(failed.record.crosschecks.len(), 2);
    assert_eq!(failed.record.chosen_level, None);
    let passed = synthetic_ladder(&p, &s, &[10. * unit, 0., 0., 0., 0., 0.]);
    assert_eq!(passed.record.status, Status::Passed);
    assert_eq!(passed.record.chosen_level, Some(3));
    assert_eq!(passed.refinement_scores.len(), 3);
}
#[test]
fn cross_kappa_disagreement_blocks_calibration() {
    let (p, s) = ladder_fixture(false);
    let result = synthetic_ladder(&p, &s, &[0., 0., 0., 0., 1., 0.]);
    assert_eq!(result.record.status, Status::ReferenceUnconverged);
    assert_eq!(result.record.crosschecks.len(), 2);
    assert!(result.crosscheck_scores[0] > 1.);
    assert_eq!(result.record.chosen_level, None);
}
#[test]
fn stagnation_above_budget_is_explicit() {
    let (p, s) = ladder_fixture(false);
    let result = synthetic_ladder(&p, &s, &[3., 2., 1., 0., 0., 0.]);
    assert_eq!(result.record.status, Status::ReferenceStagnated);
    assert!(matches!(
        result.failure,
        Some(CalibrationFailure::ReferenceStagnated(_))
    ));
    assert!(serde_json::to_string(&result)
        .unwrap()
        .contains("ReferenceStagnated"));
    assert_eq!(
        synthetic_ladder(&p, &s, &[0.; 6]).record.status,
        Status::Passed
    );
}
#[test]
fn nonfinite_reference_is_precision_limited() {
    let (p, s) = ladder_fixture(false);
    let result = synthetic_ladder(&p, &s, &[0., f64::NAN]);
    assert_eq!(result.record.status, Status::ReferencePrecisionLimited);
    assert_eq!(result.record.levels.len(), 1);
    assert!(result.failed_attempt.is_some());
    assert!(serde_json::to_string(&result)
        .unwrap()
        .contains("ReferencePrecisionLimited"));
}
#[test]
fn extended_reference_is_chosen_before_run() {
    let (p, s) = ladder_fixture(true);
    let result = synthetic_ladder(&p, &s, &[0.; 8]);
    assert_eq!(result.record.levels.len(), 6);
    assert_eq!(result.record.levels[5].settings.real_cutoff, 9. * p.a);
    assert_eq!(result.record.crosschecks.len(), 2);
    let mut illegal = p.clone();
    illegal.reference_ladder.push([10., 20.]);
    assert!(converged_reference_test(
        &s,
        &illegal,
        &mut ReferenceBudget::new(&p.limits),
        |_, _| panic!("invalid ladder evaluated")
    )
    .is_err());
}
#[test]
fn reference_deadline_preserves_failure_record() {
    let (p, s) = ladder_fixture(false);
    let mut budget = ReferenceBudget::new(&p.limits);
    let mut index = 0;
    let result = converged_reference_test(&s, &p, &mut budget, |_, b| {
        index += 1;
        if index == 3 {
            b.calls += 1;
            b.work.real_terms += 5;
            return Err(CalibrationFailure::BudgetExceeded(
                "injected deadline".into(),
            ));
        }
        b.calls += 1;
        b.work.real_terms += 17;
        Ok(eval(vec![[1., 0., 0.]; 3]))
    })
    .unwrap();
    assert_eq!(result.record.status, Status::BudgetExceeded);
    assert_eq!(result.record.levels.len(), 2);
    assert_eq!(result.usage.calls, 3);
    assert_eq!(result.usage.work.real_terms, 39);
    assert_eq!(result.failed_attempt.as_ref().unwrap().usage, result.usage);
    assert!(serde_json::to_string(&result)
        .unwrap()
        .contains("injected deadline"));
}

#[test]
fn ladder_final_denominator_and_absolute_floor_are_used() {
    let (p, s) = ladder_fixture(false);
    let mut i = 0;
    let result = converged_reference_test(&s, &p, &mut ReferenceBudget::new(&p.limits), |_, _| {
        let force = [1., 0.5, 0.25, 0.125, 0.125, 0.125][i];
        i += 1;
        Ok(eval(vec![[force, 0., 0.]; 3]))
    })
    .unwrap();
    for level in result.record.levels.iter().skip(1) {
        assert_eq!(
            level
                .delta_from_previous
                .as_ref()
                .unwrap()
                .reference_force_rms,
            metrics(
                &result.record.levels[3].evaluation,
                &result.record.levels[3].evaluation,
                3,
                1e-8 / p.a.powi(2)
            )
            .unwrap()
            .reference_force_rms
        );
    }
    let zero = converged_reference_test(&s, &p, &mut ReferenceBudget::new(&p.limits), |_, _| {
        Ok(eval(vec![[0.; 3]; 3]))
    })
    .unwrap();
    assert_eq!(zero.record.status, Status::RelativeUndefined);
    assert_eq!(zero.record.chosen_level, Some(3));
    assert!(zero.record.crosschecks[0]
        .delta_from_previous
        .as_ref()
        .unwrap()
        .force_rms_rel
        .is_none());
    let m = zero.record.crosschecks[0]
        .delta_from_previous
        .as_ref()
        .unwrap();
    assert!(!passes_target(m, 1e-4, &p, true).unwrap());
}
#[test]
fn ladder_retained_capacity_and_failed_partial_deltas_are_bounded() {
    let (mut p, s) = ladder_fixture(false);
    p.limits.max_owned_heap_bytes = 100_000;
    let mut b = ReferenceBudget::new(&p.limits);
    let mut i = 0;
    let result = converged_reference_test(&s, &p, &mut b, |_, _| {
        i += 1;
        if i == 3 {
            return Err(CalibrationFailure::BudgetExceeded("stop after two".into()));
        }
        Ok(eval(vec![[1., 0., 0.]; 3]))
    })
    .unwrap();
    assert!(result.record.levels[1].delta_from_previous.is_some());
    assert!(result.owned_heap_bytes().unwrap() > 0);
    assert_eq!(b.owned_heap_baseline, 0);
    let mut b = ReferenceBudget::new(&p.limits);
    b.owned_heap_baseline = p.limits.max_owned_heap_bytes;
    let result =
        converged_reference_test(&s, &p, &mut b, |_, _| panic!("overfull heap evaluated")).unwrap();
    assert_eq!(result.record.status, Status::BudgetExceeded);
    assert_eq!(b.owned_heap_baseline, p.limits.max_owned_heap_bytes);
}
#[test]
fn ladder_real_three_ion_records_and_shared_family_budget() {
    // Preregistered in Task 3 evidence before this test's first physics run.
    let (mut p, s) = ladder_fixture(false);
    p.limits.max_accuracy_calls = 7;
    let mut b = ReferenceBudget::new(&p.limits);
    let first = converged_reference_with_budget(&s, &p, &mut b).unwrap();
    println!(
        "TASK3_REAL_FIRST {}",
        serde_json::to_string(&first).unwrap()
    );
    assert_eq!(first.record.levels.len(), 4);
    assert_eq!(first.record.crosschecks.len(), 2);
    assert_eq!(first.audits.len(), 6);
    assert_eq!(first.usage.calls, 6);
    let expected = if first.refinement_scores[1..].iter().all(|x| *x <= 1.)
        && first.crosscheck_scores.iter().all(|x| *x <= 1.)
    {
        Status::Passed
    } else {
        Status::ReferenceUnconverged
    };
    assert_eq!(first.record.status, expected);
    // Current input/protocol are charged by runner; this retained first result is external.
    b.owned_heap_baseline = first.owned_heap_bytes().unwrap();
    let second = converged_reference_with_budget(&s, &p, &mut b).unwrap();
    println!(
        "TASK3_REAL_SECOND {}",
        serde_json::to_string(&second).unwrap()
    );
    assert_eq!(second.record.status, Status::BudgetExceeded);
    assert_eq!(second.record.levels.len(), 1);
    assert_eq!(second.usage.calls, 7);
    assert!(second.usage.work.real_terms > first.usage.work.real_terms);
    assert!(second.usage.elapsed_ns >= first.usage.elapsed_ns);
    assert_eq!(second.failed_attempt.as_ref().unwrap().usage, second.usage);
    assert!(second.audits[0].owned_heap_bytes > first.audits[0].owned_heap_bytes);
    let resumed = ReferenceBudget::resume(
        &p.limits,
        second.usage.work.clone(),
        second.usage.calls,
        std::time::Duration::from_nanos(second.usage.elapsed_ns),
        first.owned_heap_bytes().unwrap() + second.owned_heap_bytes().unwrap(),
    )
    .unwrap();
    assert_eq!(resumed.calls, 7);
}

#[test]
fn ladder_actual_deadline_returns_a_failed_record() {
    let (mut p, s) = ladder_fixture(false);
    p.limits.max_family_wall_seconds = 1;
    let mut b = ReferenceBudget::new(&p.limits);
    let result = converged_reference_test(&s, &p, &mut b, |_, _| {
        std::thread::sleep(std::time::Duration::from_millis(1100));
        Ok(eval(vec![[1., 0., 0.]; 3]))
    })
    .unwrap();
    assert_eq!(result.record.status, Status::BudgetExceeded);
    assert_eq!(result.record.levels.len(), 1);
    assert!(result.failed_attempt.is_some());
    assert!(result.usage.elapsed_ns >= 1_000_000_000);
}

#[test]
fn ladder_metric_overflow_keeps_finite_raw_evaluation() {
    let (p, s) = ladder_fixture(false);
    let mut i = 0;
    let out = converged_reference_test(&s, &p, &mut ReferenceBudget::new(&p.limits), |_, _| {
        i += 1;
        Ok(eval(vec![[if i == 1 { 0. } else { f64::MAX }, 0., 0.]; 3]))
    })
    .unwrap();
    assert_eq!(out.record.status, Status::ReferencePrecisionLimited);
    assert_eq!(out.record.levels.len(), 2);
    assert!(out.record.levels[1].evaluation.forces[0][0].is_finite());
    assert!(serde_json::to_string(&out).is_ok());
}

#[test]
fn ladder_convenience_returns_failed_record_and_rejects_invalid_setup() {
    let (mut p, s) = ladder_fixture(false);
    p.limits.max_owned_heap_bytes = 1;
    let record = converged_reference(&s, &p).unwrap();
    assert_eq!(record.status, Status::BudgetExceeded);
    assert!(record.levels.is_empty());
    assert!(record.chosen_level.is_none());
    p.reference_ladder.push([8., 16.]);
    assert!(matches!(
        converged_reference(&s, &p),
        Err(CalibrationFailure::InvalidInput(_))
    ));
}

#[test]
fn actual_mesh_deduplication_preserves_aliases() {
    let p = protocol();
    let rows = candidate_grid(&p, [4.4; 3]).unwrap();
    assert_eq!(rows.iter().map(|r| r.aliases.len()).sum::<usize>(), 126);
    let mut keys = std::collections::BTreeSet::new();
    for row in &rows {
        assert!(keys.insert((
            row.parameters.kappa.to_bits(),
            row.parameters.real_cutoff.to_bits(),
            row.actual_mesh
        )));
        assert!(row.timing.is_none() && row.rows.is_empty());
        for alias in &row.aliases {
            let b = modena_core::system::SimBox::from_array([4.4; 3]);
            assert_eq!(row.actual_mesh, modena_core::spme::mesh_for_box(&b, *alias));
        }
    }
    // Small cells force several raw kcut values onto the same production mesh.
    let tiny = candidate_grid(&p, [0.25; 3]).unwrap();
    assert!(tiny.len() < 126);
    assert!(tiny.iter().any(|r| r.aliases.len() > 1));
    assert!(candidate_grid(&p, [f64::MAX; 3]).is_err());
}

#[test]
fn one_shift_failure_rejects_candidate() {
    let p = protocol();
    let s = fixture();
    let mut c = candidate_grid(&p, s.cell).unwrap().remove(0);
    let good = metrics(
        &eval(vec![[1., 0., 0.]]),
        &eval(vec![[1., 0., 0.]]),
        1,
        1e-8,
    )
    .unwrap();
    c.rows = (0..7)
        .map(|shift_index| MetricRow {
            snapshot_id: s.id.clone(),
            shift_index,
            metrics: good.clone(),
            passing_targets: p.targets.clone(),
            status: Status::Passed,
        })
        .collect();
    c.status = Status::Passed;
    assert_eq!(
        feasible_targets(&c, &[s.id.as_str()], &p).unwrap(),
        p.targets
    );
    c.rows[6].metrics.w_per_atom_abs = 1.;
    c.rows[6].passing_targets.clear();
    c.rows[6].status = Status::Failed;
    assert!(feasible_targets(&c, &[s.id.as_str()], &p)
        .unwrap()
        .is_empty());
    c.rows.pop();
    assert!(feasible_targets(&c, &[s.id.as_str()], &p).is_err());
}

#[test]
fn relative_undefined_does_not_certify_target() {
    let p = protocol();
    let m = metrics(&eval(vec![[0.; 3]]), &eval(vec![[0.; 3]]), 1, 1e-8).unwrap();
    assert_eq!(m.force_rms_rel, None);
    assert!(!passes_target(&m, 1e-3, &p, true).unwrap());
}

#[test]
fn different_families_do_not_share_denominators() {
    let a = metrics(
        &eval(vec![[1., 0., 0.]]),
        &eval(vec![[1.001, 0., 0.]]),
        1,
        1e-8,
    )
    .unwrap();
    let b = metrics(
        &eval(vec![[0.01, 0., 0.]]),
        &eval(vec![[0.011, 0., 0.]]),
        1,
        1e-8,
    )
    .unwrap();
    assert!((a.force_rms_rel.unwrap() - 0.001).abs() < 1e-14);
    assert!((b.force_rms_rel.unwrap() - 0.1).abs() < 1e-14);
}

#[test]
fn allocator_arithmetic_and_rollback_are_checked() {
    use calibration::allocator::HeapCounter;
    let counter = HeapCounter::new(64);
    assert!(counter.reserve(32));
    assert!(counter.reserve(32));
    assert!(!counter.reserve(1));
    assert_eq!(counter.live(), 64);
    counter.release(32);
    assert_eq!(counter.live(), 32);
    assert!(!counter.reserve(usize::MAX));
    counter.release(32);
    assert_eq!(counter.live(), 0);
    assert_eq!(counter.peak(), 64);
    assert!(!counter.release(1));
    assert_eq!(counter.live(), 0);
}

#[test]
fn allocator_child_probe() {
    let Ok(mode) = std::env::var("MODENA_B05_ALLOCATOR_PROBE") else {
        return;
    };
    use calibration::allocator::{self, NoAllocationGuard};
    match mode.as_str() {
        "cap" => {
            let n = allocator::PROCESS_HEAP_CAP + 1;
            std::hint::black_box(vec![0_u8; std::hint::black_box(n)]);
        }
        "alloc" => {
            let guard = NoAllocationGuard::enter().unwrap();
            assert!(NoAllocationGuard::enter().is_none());
            std::hint::black_box(vec![1_u8; std::hint::black_box(128)]);
            drop(guard);
        }
        "dealloc" => {
            let value = vec![1_u8; 128];
            let guard = NoAllocationGuard::enter().unwrap();
            drop(std::hint::black_box(value));
            drop(guard);
        }
        "realloc" => {
            let mut value = Vec::<u8>::with_capacity(1);
            value.push(1);
            let guard = NoAllocationGuard::enter().unwrap();
            value.reserve(std::hint::black_box(128));
            drop(guard);
        }
        "warm1" | "warm4" => {
            let threads = if mode == "warm1" { 1 } else { 4 };
            let (mut protocol, s) = special_selection(true, true);
            protocol.threads = threads;
            let params = Parameters {
                kspace_cutoff: 32. / 2.2,
                ..special_parameters()
            };
            let mut evaluator = SpmeEvaluator::with_protocol(&s, params, &protocol).unwrap();
            for _ in 0..5 {
                std::hint::black_box(evaluator.compute_without_output_copy().unwrap());
            }
            let guard = NoAllocationGuard::enter().unwrap();
            for _ in 0..3 {
                // No error formatting, cloning, clocks or returned force allocation.
                match evaluator.compute_without_output_copy() {
                    Ok(value) => {
                        std::hint::black_box(value);
                    }
                    Err(_) => std::process::exit(89),
                }
            }
            drop(guard);
        }
        _ => panic!("unknown child mode"),
    }
}

#[test]
fn allocator_rejection_is_isolated_and_typed() {
    use calibration::allocator::{HEAP_EXIT, HEAP_MARKER, TIMING_EXIT, TIMING_MARKER};
    for (mode, exit, marker) in [
        ("cap", HEAP_EXIT, HEAP_MARKER),
        ("alloc", TIMING_EXIT, TIMING_MARKER),
        ("dealloc", TIMING_EXIT, TIMING_MARKER),
        ("realloc", TIMING_EXIT, TIMING_MARKER),
    ] {
        let result = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "allocator_child_probe",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("MODENA_B05_ALLOCATOR_PROBE", mode)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(exit), "{mode}: {:?}", result);
        assert!(
            result.stderr.windows(marker.len()).any(|w| w == marker),
            "{mode}"
        );
    }
}

fn special_selection(charged: bool, scaled: bool) -> (Protocol, Snapshot) {
    let p = frozen_protocol("special", 2, 1, true).unwrap();
    let snapshots = materialize(&p).unwrap();
    let s = snapshots
        .into_iter()
        .find(|s| {
            s.role == FixtureRole::Selection
                && (s.charges.iter().sum::<f64>().abs() > 0.1) == charged
                && (s.scales.len() == 2) == scaled
                && (scaled || s.scales.is_empty())
        })
        .unwrap();
    (p, s)
}
fn special_parameters() -> Parameters {
    Parameters {
        kappa: 1.6 / 2.2,
        real_cutoff: 0.4 * 9.,
        kspace_cutoff: 24. / 2.2,
    }
}
#[test]
fn scales_are_applied_once() {
    let (_, bare) = special_selection(true, false);
    let (_, scaled) = special_selection(true, true);
    let p = special_parameters();
    let b = evaluate_spme(&bare, p.clone()).unwrap();
    let c = evaluate_spme(&scaled, p).unwrap();
    let mut expected = b.clone();
    for scale in &scaled.scales {
        let d: [f64; 3] = std::array::from_fn(|axis| {
            let x = bare.positions[scale.i][axis] - bare.positions[scale.j][axis];
            x - bare.cell[axis] * (x / bare.cell[axis]).round()
        });
        let r = d.iter().map(|x| x * x).sum::<f64>().sqrt();
        let q = (1. - scale.coulomb) * bare.charges[scale.i] * bare.charges[scale.j];
        expected.energy -= q / r;
        expected.strain_w -= q / r;
        for (axis, x) in d.iter().enumerate() {
            expected.forces[scale.i][axis] -= q * x / r.powi(3);
            expected.forces[scale.j][axis] += q * x / r.powi(3);
        }
    }
    let m = metrics(&expected, &c, 3, 1e-8).unwrap();
    assert!(
        m.force_rms_abs < 1e-12 && m.energy_per_atom_abs < 1e-12 && m.w_per_atom_abs < 1e-12,
        "{m:?}"
    );
}
#[test]
fn existing_explicit_model_parameters_unchanged() {
    let (_, s) = special_selection(false, true);
    let params = special_parameters();
    let mut evaluator = SpmeEvaluator::new(&s, params.clone()).unwrap();
    assert_eq!(evaluator.parameters(), params);
    let first = evaluator.evaluate().unwrap();
    let second = evaluator.evaluate().unwrap();
    assert_eq!(first, second); // no stale force/W accumulation
    assert_eq!(evaluator.parameters(), params);
    let model = modena_core::spme::SpmeForceModel::new(
        params.kappa,
        params.real_cutoff,
        params.kspace_cutoff,
    )
    .unwrap();
    assert_eq!(
        [model.kappa, model.real_cutoff, model.kspace_cutoff],
        [params.kappa, params.real_cutoff, params.kspace_cutoff]
    );
}
#[test]
fn mesh_lifecycle_uses_existing_cache_keys() {
    let (_, s) = special_selection(true, true);
    let p = special_parameters();
    let mut e = SpmeEvaluator::new(&s, p.clone()).unwrap();
    let original = e.evaluate().unwrap();
    assert_eq!(e.neighbor_rebuilds(), 1);
    assert_eq!(e.evaluate().unwrap(), original);
    assert_eq!(e.neighbor_rebuilds(), 1);
    for factor in [1.001, 1.2, 1.] {
        let mut changed = s.clone();
        changed.cell = s.cell.map(|x| x * factor);
        for pos in &mut changed.positions {
            for x in pos {
                *x *= factor;
            }
        }
        e.set_snapshot(&changed).unwrap();
        let cached = e.evaluate().unwrap();
        let fresh = evaluate_spme(&changed, p.clone()).unwrap();
        assert_eq!(cached, fresh);
        assert_eq!(
            e.actual_mesh(),
            modena_core::spme::mesh_for_box(
                &modena_core::system::SimBox::from_array(changed.cell),
                p.kspace_cutoff
            )
        );
    }
    assert_eq!(e.neighbor_rebuilds(), 4);
    assert_eq!(e.evaluate().unwrap(), original);
}

#[test]
fn candidate_fixed_mesh_force_and_isotropic_strain_derivatives() {
    let (_, s) = special_selection(true, true);
    let p = special_parameters();
    let central = evaluate_spme(&s, p.clone()).unwrap();
    let mesh = modena_core::spme::mesh_for_box(
        &modena_core::system::SimBox::from_array(s.cell),
        p.kspace_cutoff,
    );
    let membership = |s: &Snapshot| -> Vec<(usize, usize)> {
        let mut pairs = vec![];
        for i in 0..s.positions.len() {
            for j in 0..i {
                let r2 = (0..3)
                    .map(|axis| {
                        let d = s.positions[i][axis] - s.positions[j][axis];
                        let d = d - s.cell[axis] * (d / s.cell[axis]).round();
                        d * d
                    })
                    .sum::<f64>();
                if r2 < p.real_cutoff * p.real_cutoff {
                    pairs.push((j, i));
                }
            }
        }
        pairs
    };
    let central_pairs = membership(&s);
    for h in [1e-5, 5e-6] {
        for atom in 0..3 {
            for axis in 0..3 {
                let mut plus = s.clone();
                let mut minus = s.clone();
                plus.positions[atom][axis] += h;
                minus.positions[atom][axis] -= h;
                assert_eq!(membership(&plus), central_pairs);
                assert_eq!(membership(&minus), central_pairs);
                let ep = evaluate_spme(&plus, p.clone()).unwrap().energy;
                let em = evaluate_spme(&minus, p.clone()).unwrap().energy;
                assert!((-(ep - em) / (2. * h) - central.forces[atom][axis]).abs() < 2e-8);
            }
        }
        let mut energies = vec![];
        for lambda in [1. + h, 1. - h] {
            let mut perturbed = s.clone();
            perturbed.cell = s.cell.map(|x| x * lambda);
            for pos in &mut perturbed.positions {
                for x in pos {
                    *x *= lambda;
                }
            }
            assert_eq!(membership(&perturbed), central_pairs);
            assert_eq!(
                modena_core::spme::mesh_for_box(
                    &modena_core::system::SimBox::from_array(perturbed.cell),
                    p.kspace_cutoff
                ),
                mesh,
                "MeshBoundary cannot certify production W"
            );
            energies.push(evaluate_spme(&perturbed, p.clone()).unwrap().energy);
        }
        assert!((-(energies[0] - energies[1]) / (2. * h) - central.strain_w).abs() < 2e-8);
    }
}

#[test]
fn warmed_production_guard_retains_serial_success_and_parallel_rejection() {
    for (mode, threads) in [("warm1", "1"), ("warm4", "4")] {
        let result = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "allocator_child_probe",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("MODENA_B05_ALLOCATOR_PROBE", mode)
            .env("RAYON_NUM_THREADS", threads)
            .output()
            .unwrap();
        if mode == "warm1" {
            assert!(result.status.success(), "{mode}: {result:?}");
        } else {
            assert_eq!(
                result.status.code(),
                Some(calibration::allocator::TIMING_EXIT),
                "{result:?}"
            );
            assert!(result
                .stderr
                .windows(calibration::allocator::TIMING_MARKER.len())
                .any(|w| w == calibration::allocator::TIMING_MARKER));
            println!("WARM4_INVALID {}", String::from_utf8_lossy(&result.stderr));
        }
    }
}
#[test]
fn actual_selection_candidate_retains_all_shifts_and_reference_certificates() {
    let prepared = if let Some(dir) = std::env::var_os("MODENA_B05_TASK4_PREPARED") {
        std::path::PathBuf::from(dir)
    } else {
        let dir = std::env::temp_dir().join(format!(
            "modena-b05-task4-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        prepare(&dir, frozen_protocol("special", 2, 1, true).unwrap()).unwrap();
        dir
    };
    let (p, snapshots, _) = validate_prepared(&prepared).unwrap();
    assert_eq!(p.family, "special-3");
    assert_eq!(p.reference_ladder.len(), 6);
    let mut reference_budget = ReferenceBudget::new(&p.limits);
    let mut references = Vec::new();
    for s in snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
    {
        // Charge OTHER snapshots, protocol and the retained outer/results vector;
        // current snapshot and protocol are charged by the reference runner itself.
        let others = snapshots
            .iter()
            .filter(|other| other.id != s.id)
            .map(snapshot_test_heap)
            .sum::<u64>();
        reference_budget.owned_heap_baseline = others
            + (snapshots.capacity() * std::mem::size_of::<Snapshot>()) as u64
            + (references.capacity() * std::mem::size_of::<ReferenceOutcome>()) as u64
            + references
                .iter()
                .map(|r: &ReferenceOutcome| r.owned_heap_bytes().unwrap())
                .sum::<u64>();
        let outcome = converged_reference_with_budget(s, &p, &mut reference_budget).unwrap();
        if let Some(root) = std::env::var_os("MODENA_B05_TASK4_OUTPUT") {
            let root = std::path::PathBuf::from(root);
            std::fs::create_dir_all(&root).unwrap();
            let path = root.join(format!("reference-{}.json", s.id));
            std::fs::write(path, serde_json::to_vec_pretty(&outcome).unwrap()).unwrap();
        }
        assert!(
            outcome.record.chosen_level.is_some(),
            "{}: {:?}",
            s.id,
            outcome.record.status
        );
        references.push(outcome);
    }
    reference_budget.owned_heap_baseline = 0;
    let mut budget = CandidateBudget::new(&p.limits, reference_budget);
    let mut certificates = vec![];
    let candidate = candidate_grid(&p, snapshots[0].cell)
        .unwrap()
        .into_iter()
        .find(|c| {
            c.parameters.kappa == 1.6 / 2.2
                && c.parameters.real_cutoff == 0.4 * 9.
                && c.parameters.kspace_cutoff == 24. / 2.2
        })
        .unwrap();
    let external = (snapshots.capacity() * std::mem::size_of::<Snapshot>()) as u64
        + snapshots.iter().map(snapshot_test_heap).sum::<u64>()
        + (p.schema.capacity()
            + p.fixture_set.capacity()
            + p.family.capacity()
            + p.targets.capacity() * 8
            + p.kappa_a.capacity() * 8
            + p.cutoff_lmin.capacity() * 8
            + p.kcut_a.capacity() * 8
            + p.reference_ladder.capacity() * 16
            + p.translations.capacity() * 24
            + candidate.family.capacity()
            + candidate.aliases.capacity() * 8
            + prepared.capacity()) as u64;
    budget.external_owned_heap_bytes = external;
    let calls_before = budget.reference.calls;
    let outcome = evaluate_prepared_candidate(
        &prepared,
        &candidate,
        &references,
        &mut budget,
        &mut certificates,
    )
    .unwrap();
    assert!(outcome.failure.is_none(), "{:?}", outcome.failure);
    let selection = snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .count();
    assert_eq!(outcome.raw.len(), selection * 7);
    assert_eq!(outcome.report.rows.len(), selection * 7);
    assert_eq!(outcome.translation_variation.len(), selection);
    assert_eq!(certificates.len(), selection);
    assert!(certificates
        .iter()
        .all(|c| c.status == Status::Passed && c.shifts.len() == 6));
    assert_eq!(budget.reference.calls, calls_before + 6 * selection as u64);
    assert_eq!(budget.accuracy_calls, 7 * selection as u64);
    assert!(outcome.report.timing.is_none());
    for s in snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
    {
        let r = references
            .iter()
            .find(|r| r.record.snapshot_id == s.id)
            .unwrap();
        let re = &r.record.levels[r.record.chosen_level.unwrap()].evaluation;
        for shift in 0..7 {
            let row = outcome
                .report
                .rows
                .iter()
                .find(|row| row.snapshot_id == s.id && row.shift_index == shift)
                .unwrap();
            let raw = outcome
                .raw
                .iter()
                .find(|row| row.snapshot_id == s.id && row.shift_index == shift)
                .unwrap();
            assert_eq!(
                row.metrics,
                metrics(re, &raw.evaluation, 3, 1e-8 / (p.a * p.a)).unwrap()
            );
        }
    }
    let before = budget.reference.calls;
    budget.external_owned_heap_bytes = external + outcome.owned_heap_bytes().unwrap();
    let repeat = evaluate_prepared_candidate(
        &prepared,
        &candidate,
        &references,
        &mut budget,
        &mut certificates,
    )
    .unwrap();
    assert!(repeat.failure.is_none());
    assert_eq!(budget.reference.calls, before);
    assert_eq!(repeat.report.rows, outcome.report.rows);
    let mut observed = 0;
    let stopped = evaluate_prepared_candidate_with_observer(
        &prepared,
        &candidate,
        &references,
        &mut budget,
        &mut certificates,
        |event| {
            if let CandidateEvent::CandidateRow { row, raw, .. } = event {
                assert_eq!(row.shift_index, raw.shift_index);
                observed += 1;
                if observed == 3 {
                    return Err(CalibrationFailure::Io("journal stop control".into()));
                }
            }
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(observed, 3);
    assert_eq!(stopped.raw.len(), 3);
    assert_eq!(stopped.report.rows.len(), 3);
    assert!(matches!(stopped.failure, Some(CalibrationFailure::Io(_))));
    assert_eq!(budget.reference.calls, before);
    budget.external_owned_heap_bytes +=
        repeat.owned_heap_bytes().unwrap() + stopped.owned_heap_bytes().unwrap();
    let old = budget.accuracy_calls;
    budget.accuracy_calls = p.limits.max_accuracy_calls;
    let capped = evaluate_prepared_candidate(
        &prepared,
        &candidate,
        &references,
        &mut budget,
        &mut certificates,
    )
    .unwrap();
    assert_eq!(capped.report.status, Status::BudgetExceeded);
    assert!(capped.raw.is_empty());
    budget.accuracy_calls = old;
    let mut wrong = references.clone();
    wrong[0].record.snapshot_sha256 = "0".repeat(64);
    budget.external_owned_heap_bytes +=
        repeat.owned_heap_bytes().unwrap() + capped.owned_heap_bytes().unwrap();
    let rejected = evaluate_prepared_candidate(
        &prepared,
        &candidate,
        &wrong,
        &mut budget,
        &mut certificates,
    )
    .unwrap();
    assert!(rejected.failure.is_some());
    assert!(rejected.raw.is_empty());
    if let Some(root) = std::env::var_os("MODENA_B05_TASK4_OUTPUT") {
        let root = std::path::PathBuf::from(root);
        std::fs::create_dir_all(&root).unwrap();
        for (name, value) in [
            (
                "references.json",
                serde_json::to_value(&references).unwrap(),
            ),
            ("candidate.json", serde_json::to_value(&outcome).unwrap()),
            (
                "certificates.json",
                serde_json::to_value(&certificates).unwrap(),
            ),
            (
                "budget-failure.json",
                serde_json::to_value(&capped).unwrap(),
            ),
        ] {
            std::fs::write(root.join(name), serde_json::to_vec_pretty(&value).unwrap()).unwrap();
        }
    }
    println!(
        "TASK4 rows={} refs={} calls={} status={:?} heap={} peak={} counters={:?}",
        outcome.raw.len(),
        references.len(),
        budget.reference.calls,
        outcome.report.status,
        outcome.usage.process_requested_rust_heap_bytes,
        outcome.usage.peak_process_requested_rust_heap_bytes,
        outcome.usage.reference.work
    );
}
fn snapshot_test_heap(s: &Snapshot) -> u64 {
    (s.positions.capacity() * std::mem::size_of::<[f64; 3]>()
        + s.charges.capacity() * 8
        + s.scales.capacity() * std::mem::size_of::<Scale>()
        + s.id.capacity()
        + s.family.capacity()
        + s.units.capacity()) as u64
}
#[test]
fn full_sweep_preflight_keeps_caps_and_confirmation_roster() {
    let p = protocol();
    let snapshots = materialize(&p).unwrap();
    let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    b.accuracy_calls = p.limits.max_accuracy_calls;
    assert!(matches!(
        preflight_selection_sweep(&p, &snapshots, &b),
        Err(CalibrationFailure::BudgetExceeded(_))
    ));
    let mut p = protocol();
    p.limits.max_real_terms_per_family = 1;
    let b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    assert!(matches!(
        preflight_selection_sweep(&p, &materialize(&p).unwrap(), &b),
        Err(CalibrationFailure::BudgetExceeded(_))
    ));
}

#[test]
fn timing_guard_rejects_an_inflight_callback_before_clocks() {
    use calibration::allocator::CallbackState;
    let state = CallbackState::new();
    assert!(state.enter_callback());
    assert!(!state.begin_guard());
    state.leave_callback();
    assert!(state.begin_guard());
    assert!(!state.begin_guard()); // no nested guard can clear an active guard
    assert!(!state.enter_callback());
    state.leave_callback();
    state.end_guard();
    assert!(state.begin_guard());
    state.end_guard();
}

#[test]
fn candidate_two_pass_whole_sweep_projection_rejects_insufficient_budget() {
    let mut p = frozen_protocol("special", 2, 1, true).unwrap();
    let snapshots = materialize(&p).unwrap();
    let b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    let projected = preflight_selection_sweep(&p, &snapshots, &b).unwrap();
    // Independent rectangular count: (2*12+1)*(2*14+1)*(2*15+1)-1.
    // Each nonzero mode visits all3 atoms for rho AND again for forces.
    let per_shift = 2 * 3 * (25 * 29 * 31 - 1);
    assert_eq!(per_shift, 134_844);
    let distinct_meshes = candidate_grid(&p, snapshots[0].cell)
        .unwrap()
        .into_iter()
        .map(|c| c.actual_mesh)
        .collect::<std::collections::BTreeSet<_>>()
        .len() as u64;
    let selection = snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .count() as u64;
    assert_eq!(
        projected.atom_k_terms,
        6 * distinct_meshes * selection * per_shift
    );
    // Between the old one-pass projection and the required two-pass reservation.
    p.limits.max_atom_k_terms_per_family = projected.atom_k_terms * 3 / 4;
    let b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    assert!(matches!(
        preflight_selection_sweep(&p, &materialize(&p).unwrap(), &b),
        Err(CalibrationFailure::BudgetExceeded(_))
    ));
}

// Task5 pure controls: removing any physical gate, ordinary timing validation,
// deterministic tie field, or frozen-ID constraint must break these tests.
fn cost_candidate(id: usize, ns: u64) -> CandidateReport {
    let p = protocol();
    let s = fixture();
    let mut c = candidate_grid(&p, s.cell).unwrap().remove(0);
    c.id = id;
    c.status = Status::Passed;
    let good = metrics(
        &eval(vec![[1., 0., 0.]]),
        &eval(vec![[1., 0., 0.]]),
        1,
        1e-8,
    )
    .unwrap();
    c.rows = (0..7)
        .map(|shift_index| MetricRow {
            snapshot_id: s.id.clone(),
            shift_index,
            metrics: good.clone(),
            passing_targets: p.targets.clone(),
            status: Status::Passed,
        })
        .collect();
    c.timing = Some(TimingRecord {
        reset_included: true,
        cold_ns: 1,
        rebuild_ns: 1,
        batch_total_ns: vec![ns * 20; 5],
        calls_per_batch: 20,
        median_ns_per_call: ns as f64,
        q1_ns_per_call: ns as f64,
        q3_ns_per_call: ns as f64,
    });
    c
}
#[test]
fn selector_requires_all_gates_and_finite_timing() {
    let p = protocol();
    let s = fixture();
    let ids = [s.id.as_str()];
    assert_eq!(
        select(&[cost_candidate(7, 100)], &ids, 1e-3, &p)
            .unwrap()
            .selected_id,
        Some(7)
    );
    for gate in 0..5 {
        let mut c = cost_candidate(7, 100);
        match gate {
            0 => c.rows[6].metrics.force_rms_rel = Some(1.),
            1 => c.rows[6].metrics.force_rms_abs = 1.,
            2 => c.rows[6].metrics.max_atom_force_error = 1.,
            3 => c.rows[6].metrics.energy_per_atom_abs = 1.,
            _ => c.rows[6].metrics.w_per_atom_abs = 1.,
        }
        assert_eq!(select(&[c], &ids, 1e-3, &p).unwrap().selected_id, None);
    }
    for kind in 0..6 {
        let mut c = cost_candidate(7, 100);
        let t = c.timing.as_mut().unwrap();
        match kind {
            0 => t.median_ns_per_call = f64::NAN,
            1 => t.reset_included = false,
            2 => {
                t.batch_total_ns.pop();
            }
            3 => t.calls_per_batch = 0,
            4 => t.q3_ns_per_call = 99.,
            _ => t.batch_total_ns[0] = 0,
        }
        assert!(select(&[c], &ids, 1e-3, &p).is_err());
    }
    let mut missing = cost_candidate(8, 100);
    missing.timing = None;
    assert!(select(&[cost_candidate(7, 100), missing], &ids, 1e-3, &p).is_err());
}
#[test]
fn five_percent_tie_rule_is_deterministic() {
    let p = protocol();
    let s = fixture();
    let ids = [s.id.as_str()];
    let mut fast = cost_candidate(9, 100);
    fast.actual_mesh = [12; 3];
    let mut tied = cost_candidate(8, 105);
    tied.actual_mesh = [10; 3];
    let mut outside = cost_candidate(1, 106);
    outside.actual_mesh = [4; 3];
    for rows in [
        vec![fast.clone(), tied.clone(), outside.clone()],
        vec![outside, tied.clone(), fast],
    ] {
        assert_eq!(select(&rows, &ids, 1e-3, &p).unwrap().selected_id, Some(8));
    }
    let mut lower_rc = tied.clone();
    lower_rc.id = 6;
    lower_rc.parameters.real_cutoff *= 0.9;
    let mut lower_kappa = lower_rc.clone();
    lower_kappa.id = 5;
    lower_kappa.parameters.kappa *= 0.9;
    let mut lower_kcut = lower_kappa.clone();
    lower_kcut.id = 4;
    lower_kcut.parameters.kspace_cutoff *= 0.9;
    assert_eq!(
        select(&[tied, lower_rc, lower_kappa, lower_kcut], &ids, 1e-3, &p)
            .unwrap()
            .selected_id,
        Some(4)
    );
}
#[test]
fn holdout_failure_does_not_select_another_candidate() {
    let p = protocol();
    let s = fixture();
    let ids = [s.id.as_str()];
    let chosen = cost_candidate(7, 100);
    let mut decision = select(&[chosen.clone(), cost_candidate(8, 200)], &ids, 1e-3, &p).unwrap();
    let mut heldout = chosen;
    heldout.rows[6].metrics.w_per_atom_abs = 1.;
    finish_confirmation(&mut decision, &heldout, &ids, &p).unwrap();
    assert_eq!(decision.selected_id, Some(7));
    assert_eq!(decision.confirmation_status, Status::ConfirmationFailed);
    assert!(selected_fragment(&decision, &heldout, &p)
        .unwrap()
        .is_none());
    assert!(finish_confirmation(&mut decision, &cost_candidate(8, 1), &ids, &p).is_err());
}
#[test]
fn target_not_reached_emits_no_config() {
    let p = protocol();
    let s = fixture();
    let ids = [s.id.as_str()];
    let mut c = cost_candidate(7, 100);
    c.rows[0].metrics.force_rms_rel = Some(1.);
    let decision = select(&[c.clone()], &ids, 1e-3, &p).unwrap();
    assert_eq!(decision.selection_status, Status::TargetNotReached);
    assert_eq!(decision.selected_id, None);
    assert!(selected_fragment(&decision, &c, &p).unwrap().is_none());
    let mut passed = select(&[cost_candidate(7, 100)], &ids, 1e-3, &p).unwrap();
    let good = cost_candidate(7, 100);
    finish_confirmation(&mut passed, &good, &ids, &p).unwrap();
    let fragment = selected_fragment(&passed, &good, &p).unwrap().unwrap();
    assert!(fragment.contains("empirical helper calibration; phase profiling pending"));
    assert!(fragment.contains("units = \"reduced\""));
}
#[test]
fn ordinary_batch_child_probe() {
    if std::env::var("MODENA_B05_ORDINARY_CHILD").is_err() {
        return;
    }
    let (p, s) = special_selection(true, true);
    let mut e = SpmeEvaluator::with_protocol(&s, special_parameters(), &p).unwrap();
    for _ in 0..5 {
        e.compute_without_output_copy().unwrap();
    }
    let before = e.neighbor_rebuilds();
    let duration = e.ordinary_batch(20).unwrap();
    assert!(duration > 0);
    assert_eq!(e.neighbor_rebuilds(), before);
    let raw = e.evaluate().unwrap();
    let expected = evaluate_spme(&s, special_parameters()).unwrap();
    evaluations_close(&raw, &expected, 1e-12);
    assert!(e.ordinary_batch(0).is_err());
}
#[test]
fn ordinary_batch_includes_reset_and_reuses_persistent_model() {
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "ordinary_batch_child_probe",
            "--nocapture",
            "--test-threads=1",
        ])
        .env("MODENA_B05_ORDINARY_CHILD", "1")
        .env("RAYON_NUM_THREADS", "1")
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
}

#[test]
fn timing_summary_retains_raw_batches_and_quartiles() {
    let p = protocol();
    let t = ordinary_timing_summary(11, 13, vec![100, 500, 300, 400, 200], &p).unwrap();
    assert_eq!(t.batch_total_ns, vec![100, 500, 300, 400, 200]);
    assert_eq!(
        (t.median_ns_per_call, t.q1_ns_per_call, t.q3_ns_per_call),
        (15., 10., 20.)
    );
    assert!(ordinary_timing_summary(1, 1, vec![1; 4], &p).is_err());
}
#[test]
fn freeze_requires_complete_grid_before_any_holdout() {
    let dir = std::env::temp_dir().join(format!("b05-task5-freeze-{}", std::process::id()));
    if dir.exists() {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    let p = frozen_protocol("special", 2, 1, true).unwrap();
    prepare(&dir, p.clone()).unwrap();
    let snapshots = materialize(&p).unwrap();
    let ids: Vec<_> = snapshots
        .iter()
        .filter(|s| s.role == FixtureRole::Selection)
        .map(|s| s.id.as_str())
        .collect();
    let mut grid = candidate_grid(&p, snapshots[0].cell).unwrap();
    for c in &mut grid {
        let sample = cost_candidate(c.id, 100);
        c.rows = ids
            .iter()
            .flat_map(|id| {
                sample.rows.iter().map(move |row| MetricRow {
                    snapshot_id: (*id).into(),
                    ..row.clone()
                })
            })
            .collect();
        c.timing = sample.timing;
        c.status = Status::Passed;
    }
    assert!(freeze_prepared_selection(&dir, &grid[..grid.len() - 1]).is_err());
    let frozen = freeze_prepared_selection(&dir, &grid).unwrap();
    assert_eq!(frozen.decisions().len(), 2);
    assert!(frozen
        .decisions()
        .iter()
        .all(|d| d.selected_id.is_some() && d.confirmation_status == Status::Pending));
    let frozen_path = dir.join("frozen.json");
    let frozen_sha = frozen.write_new(&frozen_path).unwrap();
    let reloaded = FrozenSelection::load(&dir, &frozen_path, &frozen_sha).unwrap();
    assert_eq!(frozen.decisions(), reloaded.decisions());
    assert!(FrozenSelection::load(&dir, &frozen_path, &"0".repeat(64)).is_err());
    let mut timing = TimingBudget {
        calls: p.limits.max_timing_calls,
    };
    let mut budget = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    assert!(matches!(
        measure_feasible_candidates(&dir, &mut grid, &mut timing, &mut budget),
        Err(CalibrationFailure::BudgetExceeded(_))
    ));
    assert_eq!(budget.reference.calls, 0);
    assert_eq!(timing.calls, 20_000);
    grid[0].timing = None;
    assert!(freeze_prepared_selection(&dir, &grid).is_err());
    for c in &mut grid {
        c.timing = None;
        for row in &mut c.rows {
            row.metrics.force_rms_rel = Some(1.);
        }
    }
    let not_reached = freeze_prepared_selection(&dir, &grid).unwrap();
    let result = not_reached
        .confirm(&dir, &mut budget, &mut Vec::new())
        .unwrap();
    assert_eq!(budget.reference.calls, 0);
    assert_eq!(budget.accuracy_calls, 0);
    assert!(result.references.is_none() && result.candidates.is_empty());
    assert!(result
        .decisions
        .iter()
        .all(|d| d.selection_status == Status::TargetNotReached));

    // A persisted reference file cannot be trusted using its own supplied digest.
    let path = dir.join("untrusted-reference.json");
    std::fs::write(&path, b"{}").unwrap();
    let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
    assert!(PreparedReferences::load_selection(&dir, &path, &"0".repeat(64), &mut b).is_err());
    assert_eq!(b.accuracy_calls, 0);
    assert_eq!(b.reference.calls, 0);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn ordinary_cost_observer_preserves_completed_batches() {
    if cfg!(debug_assertions) {
        return;
    }
    // Isolated release-only child: the global guard must never observe another
    // harness thread. This is exercised by the bounded smoke, not CI speed.
    if std::env::var("MODENA_B05_COST_OBSERVER_CHILD").is_ok() {
        let prepared = if let Some(path) = std::env::var_os("MODENA_B05_TASK5_PREPARED") {
            std::path::PathBuf::from(path)
        } else {
            let dir =
                std::env::temp_dir().join(format!("b05-task5-observer-{}", std::process::id()));
            prepare(&dir, frozen_protocol("special", 2, 1, true).unwrap()).unwrap();
            dir
        };
        let (p, snapshots, _) = validate_prepared(&prepared).unwrap();
        let s = snapshots
            .iter()
            .find(|s| s.id == "special-3-selection-charged-scaled")
            .unwrap();
        let mut b = CandidateBudget::new(&p.limits, ReferenceBudget::new(&p.limits));
        let mut timing = TimingBudget::default();
        let out = measure_ordinary_with_observer(
            s,
            special_parameters(),
            &p,
            &mut timing,
            &mut b,
            |event| {
                if event.stage == CostStage::Batch && event.outcome.batch_total_ns.len() == 2 {
                    return Err(CalibrationFailure::Io(
                        "induced cost journal failure".into(),
                    ));
                }
                Ok(())
            },
        )
        .unwrap();
        assert!(matches!(out.failure, Some(CalibrationFailure::Io(_))));
        assert!(out.cold_ns.is_some() && out.rebuild_ns.is_some());
        assert_eq!(out.warmups_completed, 5);
        assert_eq!(out.batch_total_ns.len(), 2);
        assert_eq!(timing.calls, 47);
        assert!(out.timing.is_none());
        return;
    }
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "ordinary_cost_observer_preserves_completed_batches",
            "--nocapture",
            "--test-threads=1",
        ])
        .env("MODENA_B05_COST_OBSERVER_CHILD", "1")
        .env("RAYON_NUM_THREADS", "1")
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
}
