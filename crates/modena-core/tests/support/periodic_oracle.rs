//! Test-only fixed-cube oracle. No production geometric search is imported.
use modena_core::{GraphBudget, PeriodicGraph, PeriodicGraphOptions, SimBox, System, Vec3};
use serde::{Deserialize, Serialize};

pub const CASE_BYTES: &str = include_str!("../fixtures/periodic_graph/cases.json");
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Fixture {
    pub id: String,
    pub cell: [[f64; 3]; 3],
    pub positions: Vec<[f64; 3]>,
    pub atomic_numbers: Vec<u8>,
    pub pbc: [bool; 3],
    pub cutoff: f64,
    pub periodic_self: bool,
    pub cube: i64,
    pub official: bool,
}
pub type Tuple = (usize, usize, [i64; 3]);
#[derive(Debug, Serialize)]
pub struct WireEdge {
    pub source: usize,
    pub destination: usize,
    pub edge_cell_shift: [i64; 3],
}
pub fn fixtures() -> Vec<Fixture> {
    serde_json::from_str(CASE_BYTES).expect("fixed fixture JSON")
}
pub fn options(f: &Fixture) -> PeriodicGraphOptions {
    PeriodicGraphOptions {
        cutoff: f.cutoff,
        pbc: f.pbc,
        include_periodic_self_images: f.periodic_self,
        budget: GraphBudget {
            max_atom_pairs: 1_000_000,
            max_image_candidates: 10_000_000,
            max_edges: 1_000_000,
            max_bytes: 134_217_728,
            max_abs_shift: 1_000_000,
        },
    }
}
pub fn system(f: &Fixture) -> System {
    assert_eq!(f.positions.len(), f.atomic_numbers.len());
    let mut s = System::new(SimBox::triclinic(f.cell.map(Vec3)));
    s.resize(f.positions.len());
    for (i, &r) in f.positions.iter().enumerate() {
        s.set_pos(i, Vec3(r));
        s.atomic_numbers[i] = f.atomic_numbers[i];
    }
    s
}
pub fn actual(f: &Fixture) -> Vec<Tuple> {
    let s = system(f);
    PeriodicGraph::build(&s, options(f))
        .unwrap_or_else(|e| panic!("{}: {e}", f.id))
        .edges()
        .iter()
        .map(|e| (e.source, e.destination, e.edge_cell_shift))
        .collect()
}
pub fn wire_edges(edges: &[Tuple]) -> Vec<WireEdge> {
    edges
        .iter()
        .map(|&(source, destination, edge_cell_shift)| WireEdge {
            source,
            destination,
            edge_cell_shift,
        })
        .collect()
}
pub fn displacement(f: &Fixture, edge: Tuple) -> [f64; 3] {
    let (i, j, shift) = edge;
    std::array::from_fn(|k| {
        let x0 = f.positions[j][k] - f.positions[i][k];
        let x1 = x0 + f.cell[0][k] * shift[0] as f64;
        let x2 = x1 + f.cell[1][k] * shift[1] as f64;
        x2 + f.cell[2][k] * shift[2] as f64
    })
}
fn norm2(d: [f64; 3]) -> f64 {
    (d[0] * d[0] + d[1] * d[1]) + d[2] * d[2]
}
fn visit_cube(f: &Fixture, mut visit: impl FnMut(Tuple, f64)) {
    assert_eq!(f.cube, 16, "the oracle cube is fixed, never search-derived");
    for i in 0..f.positions.len() {
        for j in 0..f.positions.len() {
            for sx in -16..=16 {
                for sy in -16..=16 {
                    for sz in -16..=16 {
                        let shift = [sx, sy, sz];
                        if (0..3).any(|k| !f.pbc[k] && shift[k] != 0)
                            || (i == j && (shift == [0; 3] || !f.periodic_self))
                        {
                            continue;
                        }
                        let edge = (i, j, shift);
                        let n2 = norm2(displacement(f, edge));
                        assert!(n2.is_finite(), "{}: oracle arithmetic", f.id);
                        visit(edge, n2);
                    }
                }
            }
        }
    }
}
pub fn fixed_cube(f: &Fixture) -> Vec<Tuple> {
    let mut edges = Vec::new();
    visit_cube(f, |edge, n2| {
        if n2 < f.cutoff * f.cutoff {
            edges.push(edge);
        }
    });
    edges.sort_unstable();
    edges
}
pub fn cutoff_margin(f: &Fixture) -> f64 {
    let mut margin = f64::INFINITY;
    visit_cube(f, |_, n2| margin = margin.min((n2.sqrt() - f.cutoff).abs()));
    margin
}
