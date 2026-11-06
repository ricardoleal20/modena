use super::*;
const A: f64 = 2.2;
pub fn frozen_protocol(
    family: &str,
    replicas: u32,
    threads: usize,
    extended: bool,
) -> Result<Protocol, CalibrationFailure> {
    if ![2, 3].contains(&replicas)
        || ![1, 4].contains(&threads)
        || !matches!(family, "cubic" | "ortho" | "special")
        || family == "special" && replicas != 2
    {
        return Err(CalibrationFailure::UnsupportedDomain(
            "family cubic|ortho with replicas 2|3, special with replicas 2; threads 1|4".into(),
        ));
    }
    let n = if family == "special" {
        3
    } else {
        8 * replicas.pow(3)
    };
    let mut ladder = vec![[4., 8.], [5., 10.], [6., 12.], [7., 14.]];
    if extended {
        ladder.extend([[8., 16.], [9., 18.]]);
    }
    Ok(Protocol {
        schema: SCHEMA.into(),
        fixture_set: "c2-spme-v1".into(),
        family: format!("{family}-{n}"),
        replicas,
        a: A,
        q_unit: 1.,
        targets: vec![1e-3, 1e-4],
        kappa_a: vec![1.2, 1.6, 2., 2.4, 3., 3.6],
        cutoff_lmin: vec![0.30, 0.40, 0.49],
        kcut_a: vec![8., 12., 16., 24., 32., 48., 64.],
        reference_ladder: ladder,
        reference_kappa_a: 1.1,
        crosscheck_kappa_a: [0.9, 1.3],
        translations: vec![
            [0., 0., 0.],
            [0.125, 0.375, 0.625],
            [0.25, 0.5, 0.75],
            [0.5, 0.5, 0.5],
            [0.75, 0.25, 0.5],
            [0.37, 0.61, 0.83],
            [0.93, 0.17, 0.41],
        ],
        limits: Limits {
            max_atoms: 256,
            max_mesh_points: 1_048_576,
            max_candidates: 126,
            max_accuracy_calls: 25_000,
            max_timing_calls: 20_000,
            max_profile_calls: 1_000,
            max_real_terms_per_call: 100_000_000,
            max_atom_k_terms_per_call: 100_000_000,
            max_real_terms_per_family: 5_000_000_000,
            max_atom_k_terms_per_family: 500_000_000,
            max_owned_heap_bytes: 268_435_456,
            max_peak_rss_bytes: 536_870_912,
            max_family_wall_seconds: 1800,
            max_threads: 4,
        },
        tolerances: Tolerances {
            reference_force_rms_abs_units: 1e-8,
            reference_force_target_fraction: 0.01,
            reference_force_max_units: 1e-7,
            reference_energy_per_atom_units: 1e-9,
            reference_w_per_atom_units: 1e-9,
            relative_denominator_floor: 1e-8,
            candidate_force_rms_abs_factor: 1.,
            candidate_force_max_factor: 10.,
            candidate_energy_factor: 1.,
            candidate_w_factor: 10.,
            cost_tie_fraction: 0.05,
        },
        threads,
        warmups: 5,
        batches: 5,
        batch_calls: 20,
    })
}
fn shell(
    p: &Protocol,
    id: String,
    role: FixtureRole,
    cell: [f64; 3],
    positions: Vec<[f64; 3]>,
    charges: Vec<f64>,
) -> Snapshot {
    Snapshot {
        id,
        role,
        family: p.family.clone(),
        units: "reduced".into(),
        pbc: [true; 3],
        seed: None,
        amplitude: None,
        cell,
        positions,
        charges,
        scales: vec![],
    }
}
fn ionic(
    p: &Protocol,
    seed: Option<u64>,
    amplitude: f64,
    role: FixtureRole,
    label: &str,
) -> Snapshot {
    let sys = modena_core::lattice::rocksalt(A, [p.replicas; 3], "Na", "Cl", 1., 1., 1., -1.);
    let edge = A * p.replicas as f64;
    let factors = if p.family.starts_with("ortho-") {
        [1., 1.15, 1.3]
    } else {
        [1.; 3]
    };
    let cell = factors.map(|f| edge * f);
    let positions = sys
        .positions
        .chunks_exact(3)
        .enumerate()
        .map(|(i, r)| {
            std::array::from_fn(|axis| {
                let mut v = r[axis].rem_euclid(edge);
                if let Some(seed) = seed {
                    let z = ((i as u64 + 1) * 1_664_525
                        + (axis as u64 + 1) * 1_013_904_223
                        + seed * 69_069)
                        % 65_536;
                    v += amplitude * A * (2. * z as f64 / 65_535. - 1.);
                }
                // Perturb cubic coordinates, wrap, then remap fractional coordinates.
                if factors == [1.; 3] {
                    v.rem_euclid(edge)
                } else {
                    (v.rem_euclid(edge) / edge) * cell[axis]
                }
            })
        })
        .collect();
    let mut s = shell(
        p,
        format!("{}-{label}", p.family),
        role,
        cell,
        positions,
        sys.charges,
    );
    s.seed = seed;
    s.amplitude = Some(amplitude);
    s
}
fn wrap_regressions(snapshots: &mut Vec<Snapshot>, base: &Snapshot) {
    // Exercise each signed full-box wrap before canonical materialization/hashing.
    // Later physical wrap checks must also apply the raw shift at the evaluator boundary.
    for axis in 0..3 {
        for sign in [-1., 1.] {
            let mut s = base.clone();
            s.id = format!(
                "{}-wrap-{axis}-{}",
                base.id,
                if sign > 0. { "plus" } else { "minus" }
            );
            s.role = FixtureRole::RegressionOnly;
            s.positions[0][axis] += sign * s.cell[axis];
            snapshots.push(s);
        }
    }
}
pub fn materialize(p: &Protocol) -> Result<Vec<Snapshot>, CalibrationFailure> {
    let family = p
        .family
        .split('-')
        .next()
        .ok_or_else(|| invalid("family"))?;
    let expected = frozen_protocol(family, p.replicas, p.threads, p.reference_ladder.len() == 6)?;
    // Limits may be reduced for preflight, never increased. Other frozen physical fields are exact.
    let mut physical = p.clone();
    physical.limits = expected.limits.clone();
    if physical != expected {
        return Err(invalid("protocol differs from frozen physical protocol"));
    }
    let caps = serde_json::to_value(&p.limits).map_err(io)?;
    let maximum = serde_json::to_value(&expected.limits).map_err(io)?;
    for (key, value) in caps.as_object().ok_or_else(|| invalid("limits record"))? {
        let cap = value.as_u64().ok_or_else(|| invalid("limit integer"))?;
        if cap == 0
            || cap
                > maximum[key]
                    .as_u64()
                    .ok_or_else(|| invalid("limit maximum"))?
        {
            return Err(CalibrationFailure::BudgetExceeded(format!(
                "invalid/increased limit {key}"
            )));
        }
    }
    if p.threads > p.limits.max_threads {
        return Err(CalibrationFailure::BudgetExceeded("thread cap".into()));
    }
    let candidates = (p.kappa_a.len() as u64)
        .checked_mul(p.cutoff_lmin.len() as u64)
        .and_then(|n| n.checked_mul(p.kcut_a.len() as u64))
        .ok_or_else(|| CalibrationFailure::BudgetExceeded("candidate product overflow".into()))?;
    if candidates > p.limits.max_candidates as u64 {
        return Err(CalibrationFailure::BudgetExceeded("candidate cap".into()));
    }
    let n = if family == "special" {
        3
    } else {
        8 * p.replicas.pow(3) as usize
    };
    if n > p.limits.max_atoms {
        return Err(CalibrationFailure::BudgetExceeded(
            "fixture atoms before lattice allocation".into(),
        ));
    }
    let heap = (n as u64)
        .checked_mul(128 * 32)
        .ok_or_else(|| CalibrationFailure::BudgetExceeded("fixture bytes overflow".into()))?;
    check_work(
        [1; 3],
        &WorkCounts {
            real_terms: 0,
            atom_k_terms: 0,
        },
        heap,
        &p.limits,
    )?;
    let mut result = vec![];
    if family == "special" {
        for (confirmation, role) in [
            (false, FixtureRole::Selection),
            (true, FixtureRole::Confirmation),
        ] {
            for charged in [false, true] {
                for (label, scales) in [
                    ("bare", vec![]),
                    (
                        "scaled",
                        vec![
                            Scale {
                                i: 0,
                                j: 1,
                                coulomb: 0.,
                            },
                            Scale {
                                i: 0,
                                j: 2,
                                coulomb: 0.4,
                            },
                        ],
                    ),
                    (
                        "c1",
                        vec![Scale {
                            i: 0,
                            j: 1,
                            coulomb: 1.,
                        }],
                    ),
                    (
                        "c07",
                        vec![Scale {
                            i: 0,
                            j: 1,
                            coulomb: 0.7,
                        }],
                    ),
                ] {
                    let mut s = shell(
                        p,
                        format!(
                            "special-3-{}-{}-{label}",
                            if confirmation {
                                "confirmation"
                            } else {
                                "selection"
                            },
                            if charged { "charged" } else { "neutral" }
                        ),
                        role.clone(),
                        [9., 10., 11.],
                        vec![[1.1, 2.3, 3.4], [2.7, 3.2, 3.8], [6.2, 5.8, 7.]],
                        vec![1., -0.7, if charged { 0.4 } else { -0.3 }],
                    );
                    if confirmation {
                        for (k, d) in [0.031, -0.017, 0.023].into_iter().enumerate() {
                            s.positions[1][k] += d;
                        }
                    }
                    s.scales = scales;
                    result.push(s);
                }
            }
        }
    } else {
        let cases = [
            (7, 0.02, FixtureRole::Selection),
            (19, 0.05, FixtureRole::Selection),
            (31, 0.08, FixtureRole::Selection),
            (43, 0.035, FixtureRole::Confirmation),
            (59, 0.065, FixtureRole::Confirmation),
        ];
        for (seed, amp, role) in cases {
            result.push(ionic(
                p,
                Some(seed),
                amp,
                role,
                &format!("neutral-seed-{seed}"),
            ));
        }
        for (seed, amp, role) in [
            (19, 0.05, FixtureRole::Selection),
            (59, 0.065, FixtureRole::Confirmation),
        ] {
            let mut s = ionic(p, Some(seed), amp, role, &format!("background-seed-{seed}"));
            s.charges[1] = -0.5;
            result.push(s);
        }
        result.push(ionic(p, None, 0., FixtureRole::RegressionOnly, "crystal"));
    }
    let base = result[0].clone();
    wrap_regressions(&mut result, &base);
    // Snapshot arrays are the exact wrapped materialization. Wrap identities retain distinct roles/IDs.
    for s in &mut result {
        validate_snapshot(s, p)?;
    }
    let calls = result
        .iter()
        .filter(|s| s.role != FixtureRole::RegressionOnly)
        .count() as u64
        * 7
        * 126;
    if calls > p.limits.max_accuracy_calls {
        return Err(CalibrationFailure::BudgetExceeded(
            "family selection+confirmation call cap".into(),
        ));
    }
    Ok(result)
}
pub fn validate_snapshot(s: &mut Snapshot, p: &Protocol) -> Result<(), CalibrationFailure> {
    let fail = |reason: &str| invalid(format!("{}: {reason}", s.id));
    if s.units != "reduced" || s.pbc != [true; 3] || s.family != p.family {
        return Err(CalibrationFailure::UnsupportedDomain(format!(
            "{}: reduced units, full PBC and exact family required",
            s.id
        )));
    }
    let n = s.positions.len();
    if n == 0 || n != s.charges.len() {
        return Err(fail("nonempty matching positions/charges required"));
    }
    if n > p.limits.max_atoms {
        return Err(CalibrationFailure::BudgetExceeded(format!(
            "{}: atom cap",
            s.id
        )));
    }
    if !p.a.is_finite()
        || p.a <= 0.
        || s.cell.iter().any(|x| !x.is_finite() || *x <= 0.)
        || s.positions
            .iter()
            .flatten()
            .chain(s.charges.iter())
            .any(|x| !x.is_finite())
        || s.charges.iter().all(|x| *x == 0.)
        || s.amplitude.is_some_and(|a| !a.is_finite() || a < 0.)
    {
        return Err(fail(
            "finite positive cell/scales and nonzero charges required",
        ));
    }
    let mut scales = BTreeMap::new();
    for c in &s.scales {
        if c.i == c.j
            || c.i >= n
            || c.j >= n
            || !c.coulomb.is_finite()
            || !(0.0..=1.0).contains(&c.coulomb)
        {
            return Err(fail("invalid pair scale/index"));
        }
        let key = (c.i.min(c.j), c.i.max(c.j));
        if scales
            .insert(key, c.coulomb)
            .is_some_and(|old| old != c.coulomb)
        {
            return Err(fail("conflicting duplicate pair policy"));
        }
    }
    for i in 0..n {
        for j in 0..i {
            let r2 = (0..3)
                .map(|d| {
                    let delta = s.positions[i][d] - s.positions[j][d];
                    let delta = delta - s.cell[d] * (delta / s.cell[d]).round();
                    delta * delta
                })
                .sum::<f64>();
            if !r2.is_finite() || r2.sqrt() < 0.2 * p.a {
                return Err(fail("minimum separation below 0.2*a"));
            }
        }
    }
    s.scales = scales
        .into_iter()
        .map(|((i, j), coulomb)| Scale { i, j, coulomb })
        .collect();
    for r in &mut s.positions {
        for (d, x) in r.iter_mut().enumerate() {
            *x = x.rem_euclid(s.cell[d]);
            if *x == s.cell[d] {
                *x = 0.;
            }
        }
    }
    Ok(())
}
pub fn validate_frozen(p: &Protocol, snapshots: &[Snapshot]) -> Result<(), CalibrationFailure> {
    let canonical = materialize(p)?;
    if snapshots != canonical {
        return Err(CalibrationFailure::UnsupportedDomain(
            "custom fixtures/roles are NotCalibrated; exact frozen arrays required".into(),
        ));
    }
    Ok(())
}
