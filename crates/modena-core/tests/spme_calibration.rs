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
                &format!("default crystal refinement={i}"),
            )
        })
        .collect();
    assert!(!scores.windows(2).any(|pair| pair[0] <= 1. && pair[1] <= 1.));
    assert!(scores[1] > 1. && scores[2] <= 1.);
    println!("MADELUNG_DEFAULT status=ReferenceUnconverged scores={scores:?}; no Madelung sanity acceptance");
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
