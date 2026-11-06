use modena_core::{
    GraphBudget, GraphBuildError, PeriodicEdge, PeriodicGraph, PeriodicGraphOptions, SimBox,
    System, Vec3,
};

fn options(cutoff: f64, include_self: bool) -> PeriodicGraphOptions {
    PeriodicGraphOptions {
        cutoff,
        pbc: [true; 3],
        include_periodic_self_images: include_self,
        budget: GraphBudget {
            max_atom_pairs: 1_000_000,
            max_image_candidates: 10_000_000,
            max_edges: 1_000_000,
            max_bytes: 134_217_728,
            max_abs_shift: 1_000_000,
        },
    }
}
fn system(cell: SimBox, positions: &[Vec3]) -> System {
    let mut s = System::new(cell);
    s.resize(positions.len());
    for (i, &r) in positions.iter().enumerate() {
        s.set_pos(i, r);
    }
    s
}
fn tuples(g: &PeriodicGraph<'_>) -> Vec<(usize, usize, [i64; 3])> {
    g.edges()
        .iter()
        .map(|e| (e.source, e.destination, e.edge_cell_shift))
        .collect()
}
fn assert_budget(
    result: Result<PeriodicGraph<'_>, GraphBuildError>,
    resource: &'static str,
    limit: u64,
) {
    assert!(
        matches!(result, Err(GraphBuildError::BudgetExceeded { resource: r, limit: l }) if r == resource && l == limit)
    );
}
#[test]
fn multi_image_graph_has_exactly_four_directed_tuples() {
    let s = system(
        SimBox::cubic(2.0),
        &[Vec3::zero(), Vec3::new(0.6, 0.0, 0.0)],
    );
    let g = PeriodicGraph::build(&s, options(1.5, true)).unwrap();
    assert_eq!(
        tuples(&g),
        vec![
            (0, 1, [-1, 0, 0]),
            (0, 1, [0, 0, 0]),
            (1, 0, [0, 0, 0]),
            (1, 0, [1, 0, 0])
        ]
    );
    assert_eq!(g.stats().atom_pair_visits, 8);
    assert_eq!(g.stats().image_candidate_visits, 12);
    assert_eq!(g.stats().edges, 4);
    assert!(std::ptr::eq(g.system(), &s));
    assert_eq!(g.options().cutoff, 1.5);
    assert!(g.edges().windows(2).all(|w| w[0] < w[1]));
    for e in g.edges() {
        assert!(g.edges().contains(&PeriodicEdge {
            source: e.destination,
            destination: e.source,
            edge_cell_shift: e.edge_cell_shift.map(|x| -x)
        }));
    }
}
#[test]
fn periodic_self_policy_and_zero_self() {
    let s = system(SimBox::cubic(1.0), &[Vec3::zero()]);
    let g = PeriodicGraph::build(&s, options(1.01, true)).unwrap();
    assert_eq!(
        tuples(&g),
        vec![
            (0, 0, [-1, 0, 0]),
            (0, 0, [0, -1, 0]),
            (0, 0, [0, 0, -1]),
            (0, 0, [0, 0, 1]),
            (0, 0, [0, 1, 0]),
            (0, 0, [1, 0, 0])
        ]
    );
    assert_eq!(g.stats().image_candidate_visits, 54);
    let excluded = PeriodicGraph::build(&s, options(1.01, false)).unwrap();
    assert!(excluded.edges().is_empty());
    assert_eq!(excluded.stats().image_candidate_visits, 54);
    // Distinct coincident indices remain edges, regardless of chemical arrays.
    let mut overlap = system(SimBox::cubic(5.0), &[Vec3::zero(), Vec3::zero()]);
    overlap.atomic_numbers.clear();
    overlap.species.clear();
    overlap.forces.clear();
    assert_eq!(
        tuples(&PeriodicGraph::build(&overlap, options(0.5, false)).unwrap()),
        vec![(0, 1, [0, 0, 0]), (1, 0, [0, 0, 0])]
    );
}
#[test]
fn all_pbc_masks() {
    let s = system(SimBox::cubic(1.0), &[Vec3::zero()]);
    for mask in 0u32..8 {
        let mut o = options(1.01, true);
        o.pbc = std::array::from_fn(|k| mask & (1 << k) != 0);
        let g = PeriodicGraph::build(&s, o).unwrap();
        let mut expected = Vec::new();
        for k in 0..3 {
            if o.pbc[k] {
                for sign in [-1, 1] {
                    let mut shift = [0; 3];
                    shift[k] = sign;
                    expected.push((0, 0, shift));
                }
            }
        }
        expected.sort_unstable();
        assert_eq!(tuples(&g), expected);
        assert_eq!(
            g.stats().image_candidate_visits,
            2 * 3u64.pow(mask.count_ones())
        );
    }
}
#[test]
fn strict_cutoff_next_float() {
    for (x, count) in [(1.0f64.next_down(), 2), (1.0, 0), (1.0f64.next_up(), 0)] {
        let s = system(SimBox::cubic(8.0), &[Vec3::zero(), Vec3::new(x, 0.0, 0.0)]);
        let mut o = options(1.0, false);
        o.pbc = [false; 3];
        assert_eq!(PeriodicGraph::build(&s, o).unwrap().edges().len(), count);
    }
}
#[test]
fn malformed_geometry_fails_empty_and_nonempty() {
    for n in [0, 1] {
        for cutoff in [
            0.0,
            -1.0,
            f64::NAN,
            f64::INFINITY,
            f64::MAX,
            f64::MIN_POSITIVE,
            1e-160,
        ] {
            let s = system(SimBox::cubic(2.0), &vec![Vec3::zero(); n]);
            assert!(matches!(
                PeriodicGraph::build(&s, options(cutoff, true)),
                Err(GraphBuildError::InvalidInput(_))
            ));
        }
        for edge in [0.0, f64::NAN, f64::INFINITY] {
            let s = system(SimBox::cubic(edge), &vec![Vec3::zero(); n]);
            assert!(matches!(
                PeriodicGraph::build(&s, options(1.0, true)),
                Err(GraphBuildError::InvalidInput(_))
            ));
        }
        let cell = SimBox::triclinic([
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0f64.next_up(), 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]);
        let s = system(cell, &vec![Vec3::zero(); n]);
        assert!(matches!(
            PeriodicGraph::build(&s, options(1.0, true)),
            Err(GraphBuildError::NumericalGeometry(_))
        ));
        for shift in [-1, (1i64 << 52) + 1] {
            let s = system(SimBox::cubic(2.0), &vec![Vec3::zero(); n]);
            let mut o = options(1.0, true);
            o.budget.max_abs_shift = shift;
            assert!(matches!(
                PeriodicGraph::build(&s, o),
                Err(GraphBuildError::ShiftRange)
            ));
        }
    }
    let mut s = system(SimBox::cubic(2.0), &[Vec3::zero()]);
    s.positions.pop();
    assert!(matches!(
        PeriodicGraph::build(&s, options(1.0, true)),
        Err(GraphBuildError::InvalidInput(_))
    ));
    s.positions.push(f64::NAN);
    assert!(matches!(
        PeriodicGraph::build(&s, options(1.0, true)),
        Err(GraphBuildError::InvalidInput(_))
    ));
    let s = system(
        SimBox::cubic(2.0),
        &[
            Vec3::new(f64::MAX, 0.0, 0.0),
            Vec3::new(-f64::MAX, 0.0, 0.0),
        ],
    );
    let mut o = options(1.0, true);
    o.pbc = [false; 3];
    assert!(matches!(
        PeriodicGraph::build(&s, o),
        Err(GraphBuildError::NumericalGeometry(_))
    ));
}
#[test]
fn budget_failure_preserves_source_and_prior_graph() {
    let mut s = system(
        SimBox::cubic(2.0),
        &[Vec3::zero(), Vec3::new(0.6, 0.0, 0.0)],
    );
    s.forces.fill(7.0);
    s.atomic_numbers.copy_from_slice(&[6, 8]);
    let prior = PeriodicGraph::build(&s, options(1.5, true)).unwrap();
    let before = tuples(&prior);
    let saved = (
        s.positions.clone(),
        s.forces.clone(),
        s.atomic_numbers.clone(),
    );
    for (resource, limit) in [
        ("atom_pairs", 7),
        ("image_candidates", 11),
        ("edges", 3),
        ("bytes", 0),
    ] {
        let mut o = options(1.5, true);
        match resource {
            "atom_pairs" => o.budget.max_atom_pairs = limit,
            "image_candidates" => o.budget.max_image_candidates = limit,
            "edges" => o.budget.max_edges = limit as usize,
            _ => o.budget.max_bytes = limit as usize,
        }
        assert_budget(PeriodicGraph::build(&s, o), resource, limit);
        assert_eq!(
            (&s.positions, &s.forces, &s.atomic_numbers),
            (&saved.0, &saved.1, &saved.2)
        );
        assert_eq!(tuples(&prior), before);
    }
    let stats = prior.stats();
    let mut o = options(1.5, true);
    o.budget.max_atom_pairs = 8;
    o.budget.max_image_candidates = 12;
    o.budget.max_edges = 4;
    o.budget.max_bytes = stats.requested_bytes;
    assert_eq!(PeriodicGraph::build(&s, o).unwrap().stats(), stats);
    o.budget.max_bytes -= 1;
    assert_budget(
        PeriodicGraph::build(&s, o),
        "bytes",
        o.budget.max_bytes as u64,
    );
}
#[test]
fn shift_budget_is_not_top_k() {
    let s = system(SimBox::cubic(1.0), &[Vec3::zero()]);
    let mut o = options(1.01, true);
    o.budget.max_abs_shift = 1;
    assert!(matches!(
        PeriodicGraph::build(&s, o),
        Err(GraphBuildError::ShiftRange)
    ));
}
#[test]
fn global_range_failure_is_atomic() {
    let s = system(
        SimBox::cubic(1.0),
        &[Vec3::zero(), Vec3::new(8.0, 0.0, 0.0)],
    );
    let prior = PeriodicGraph::build(&s, options(0.25, true)).unwrap();
    let saved = tuples(&prior);
    let mut o = options(0.25, true);
    o.budget.max_abs_shift = 7;
    o.budget.max_image_candidates = 2;
    assert!(matches!(
        PeriodicGraph::build(&s, o),
        Err(GraphBuildError::ShiftRange)
    ));
    assert_eq!(tuples(&prior), saved);
}
#[test]
#[allow(non_snake_case)]
fn configured_B_never_defines_completeness() {
    let s = system(
        SimBox::cubic(1.0),
        &[Vec3::zero(), Vec3::new(8.0, 0.0, 0.0)],
    );
    let large = PeriodicGraph::build(&s, options(0.25, true)).unwrap();
    assert_eq!(tuples(&large), vec![(0, 1, [-8, 0, 0]), (1, 0, [8, 0, 0])]);
    let mut o = options(0.25, true);
    o.budget.max_abs_shift = 9;
    assert_eq!(
        tuples(&PeriodicGraph::build(&s, o).unwrap()),
        tuples(&large)
    );
}
#[test]
fn fresh_borrow_after_identity_and_box_mutations() {
    let mut s = system(
        SimBox::cubic(2.0),
        &[Vec3::zero(), Vec3::new(0.6, 0.0, 0.0)],
    );
    {
        let g = PeriodicGraph::build(&s, options(1.5, true)).unwrap();
        assert_eq!(g.edges().len(), 4);
    }
    s.atomic_numbers[0] = 8;
    s.positions.swap(0, 3);
    s.box_size = SimBox::cubic(8.0);
    let g = PeriodicGraph::build(&s, options(1.5, true)).unwrap();
    assert_eq!(tuples(&g), vec![(0, 1, [0, 0, 0]), (1, 0, [0, 0, 0])]);
    assert_eq!(g.system().atomic_numbers[0], 8);
}
#[test]
fn valid_empty_graph_has_fixed_workspace_and_zero_counts() {
    let s = system(SimBox::cubic(2.0), &[]);
    let g = PeriodicGraph::build(&s, options(1.0, true)).unwrap();
    assert!(g.edges().is_empty());
    assert_eq!(g.stats().atom_pair_visits, 0);
    assert_eq!(g.stats().image_candidate_visits, 0);
    assert_eq!(g.stats().edges, 0);
    assert!(g.stats().requested_bytes > 0);
    let mut o = options(1.0, true);
    o.budget.max_bytes = g.stats().requested_bytes - 1;
    assert_budget(
        PeriodicGraph::build(&s, o),
        "bytes",
        o.budget.max_bytes as u64,
    );
}

#[test]
fn error_precedence_validates_before_work_and_range_before_candidates() {
    let s = system(SimBox::cubic(1.0), &[Vec3::zero(), Vec3::zero()]);
    let mut o = options(f64::NAN, true);
    o.budget.max_atom_pairs = 0;
    assert!(matches!(
        PeriodicGraph::build(&s, o),
        Err(GraphBuildError::InvalidInput(_))
    ));
    let invalid = system(SimBox::cubic(0.0), &[Vec3::zero()]);
    o.cutoff = 1.01;
    assert!(matches!(
        PeriodicGraph::build(&invalid, o),
        Err(GraphBuildError::InvalidInput(_))
    ));
    o.budget.max_atom_pairs = 7;
    o.budget.max_abs_shift = 1;
    assert_budget(PeriodicGraph::build(&s, o), "atom_pairs", 7);
    o.budget.max_atom_pairs = 8;
    o.budget.max_image_candidates = 0;
    assert!(matches!(
        PeriodicGraph::build(&s, o),
        Err(GraphBuildError::ShiftRange)
    ));
    let empty = system(SimBox::cubic(1.0), &[]);
    o.budget.max_atom_pairs = 0;
    o.budget.max_edges = 0;
    let g = PeriodicGraph::build(&empty, o).unwrap();
    assert_eq!(
        (
            g.stats().atom_pair_visits,
            g.stats().image_candidate_visits,
            g.stats().edges
        ),
        (0, 0, 0)
    );
}

#[path = "support/periodic_oracle.rs"]
#[allow(dead_code)]
mod periodic_oracle;

#[test]
fn fixed_cube_exact_tuples_for_every_declared_fixture() {
    use periodic_oracle::*;
    for f in fixtures() {
        let expected = fixed_cube(&f);
        assert_eq!(actual(&f), expected, "{}", f.id);
        assert!(expected.windows(2).all(|w| w[0] < w[1]));
        for &(i, j, s) in &expected {
            assert!(
                expected.binary_search(&(j, i, s.map(|v| -v))).is_ok(),
                "{}",
                f.id
            );
            assert!((0..3).all(|k| f.pbc[k] || s[k] == 0));
        }
        if f.official {
            assert!(cutoff_margin(&f) >= 1e-8, "{}", f.id);
        }
    }
}
#[test]
fn independent_fixture_cardinalities_and_adjacent_f64_boundary() {
    use periodic_oracle::*;
    for f in fixtures() {
        let count = fixed_cube(&f).len();
        let expected = match f.id.as_str() {
            "empty" | "self_disabled" | "self_short" | "boundary" | "boundary_outside" => 0,
            "self6" => 6,
            "multi4" => 4,
            "boundary_inside" => 2,
            _ => continue,
        };
        assert_eq!(count, expected, "{}", f.id);
    }
}
#[test]
fn gauge_tuples_and_vectors_transform_by_atom_image() {
    use periodic_oracle::*;
    let all = fixtures();
    let base = all.iter().find(|f| f.id == "mace3").unwrap();
    let gauge = all.iter().find(|f| f.id == "gauge").unwrap();
    let t = [[0, 0, 0], [3, -2, 0], [0, 0, 0]];
    let mut transformed = Vec::new();
    for e @ (i, j, s) in actual(base) {
        let moved = (i, j, std::array::from_fn(|k| s[k] + t[i][k] - t[j][k]));
        let a = displacement(base, e);
        let b = displacement(gauge, moved);
        for k in 0..3 {
            assert!((a[k] - b[k]).abs() <= 2e-12 + 2e-12 * a[k].abs());
        }
        transformed.push(moved);
    }
    transformed.sort_unstable();
    assert_eq!(actual(gauge), transformed);
}
#[test]
fn independent_oracle_detects_incomplete_and_wrong_image_graphs() {
    use periodic_oracle::*;
    let all = fixtures();
    let f = all.iter().find(|f| f.id == "unreduced").unwrap();
    let expected = fixed_cube(f);
    assert!(expected.contains(&(0, 1, [-2, 0, 0])));
    let incomplete: Vec<_> = actual(f)
        .into_iter()
        .filter(|e| e.2.iter().all(|x| x.abs() <= 1))
        .collect();
    assert_ne!(
        incomplete, expected,
        "a closest/+-1-only path must be rejected"
    );
    let mut wrong = actual(f);
    wrong[0].2[0] += 1;
    assert_ne!(
        wrong, expected,
        "distance/count-only checks would miss tuple changes"
    );
}
#[test]
fn fresh_snapshots_and_options_match_independent_cube() {
    use periodic_oracle::*;
    let mut f = fixtures().into_iter().find(|f| f.id == "multi4").unwrap();
    assert_eq!(actual(&f), fixed_cube(&f));
    f.positions[1][0] = 1.1;
    assert_eq!(actual(&f), fixed_cube(&f));
    f.positions.swap(0, 1);
    f.atomic_numbers.swap(0, 1);
    assert_eq!(actual(&f), fixed_cube(&f));
    f.atomic_numbers[0] = 1;
    let s = periodic_oracle::system(&f);
    let g = PeriodicGraph::build(&s, periodic_oracle::options(&f)).unwrap();
    assert_eq!(g.system().atomic_numbers, f.atomic_numbers);
    drop(g);
    f.cutoff = 2.1;
    assert_eq!(actual(&f), fixed_cube(&f));
    f.pbc = [false, true, false];
    assert_eq!(actual(&f), fixed_cube(&f));
    f.periodic_self = false;
    assert_eq!(actual(&f), fixed_cube(&f));
    f.positions.pop();
    f.atomic_numbers.pop();
    assert_eq!(actual(&f), fixed_cube(&f));
}

#[test]
fn independent_triangular_fixture_certificates_fit_fixed_cube() {
    use periodic_oracle::*;
    for f in fixtures() {
        // Direct triangular equations only: no inverse or production bound.
        let h = f.cell;
        assert_eq!([h[0][1], h[0][2], h[1][2]], [0.0; 3]);
        let span: [f64; 3] = std::array::from_fn(|k| {
            f.positions
                .iter()
                .flat_map(|a| f.positions.iter().map(move |b| (a[k] - b[k]).abs()))
                .fold(0.0, f64::max)
        });
        // Search-only cushion; membership uses the unchanged strict R^2.
        let q = f.cutoff + 1e-8;
        let sz = ((q + span[2]) / h[2][2].abs()).floor();
        let sy = ((q + span[1] + h[2][1].abs() * sz) / h[1][1].abs()).floor();
        let sx = ((q + span[0] + h[1][0].abs() * sy + h[2][0].abs() * sz) / h[0][0].abs()).floor();
        let certificate = [sx, sy, sz];
        let expected = if f.id == "gauge" {
            [5.0, 3.0, 1.0]
        } else if f.id == "unreduced" {
            [2.0, 0.0, 0.0]
        } else if f.id == "oblique" {
            [1.0, 0.0, 0.0]
        } else if f.id.starts_with("boundary") {
            [0.0; 3]
        } else if f.id == "mace3"
            || f.id.starts_with("mace3_mask")
            || f.id == "deform"
            || f.id == "displaced"
        {
            [2.0, 2.0, 1.0]
        } else if f.id == "self6" || f.id == "self_disabled" {
            [1.0; 3]
        } else {
            [1.0, 0.0, 0.0]
        };
        assert!(
            (0..3).all(|k| certificate[k] <= expected[k]),
            "{}: {:?}",
            f.id,
            certificate
        );
        assert!(expected.iter().all(|&v| v < 16.0));
        for (_, _, shift) in fixed_cube(&f) {
            assert!(
                (0..3).all(|k| (shift[k].abs() as f64) <= expected[k]),
                "{}",
                f.id
            );
        }
    }
}
