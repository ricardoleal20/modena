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
