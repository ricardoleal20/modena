//! Directed periodic image graphs borrowed from an immutable system snapshot.

use crate::periodic_geometry::{
    certified_cell, checked_shift_endpoint, cutoff_squared, edge_displacement, image_bounds,
    strict_radius, CertifiedCell, GeometryError, GeometryOptions,
};
use crate::System;
use std::mem::size_of;

/// A directed image: `r_destination - r_source + H * edge_cell_shift`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct PeriodicEdge {
    pub source: usize,
    pub destination: usize,
    pub edge_cell_shift: [i64; 3],
}

/// Explicit work, storage, and supported image-domain limits; no truncation occurs.
#[derive(Clone, Copy, Debug)]
pub struct GraphBudget {
    pub max_atom_pairs: u64,
    pub max_image_candidates: u64,
    pub max_edges: usize,
    pub max_bytes: usize,
    pub max_abs_shift: i64,
}

#[derive(Clone, Copy, Debug)]
pub struct PeriodicGraphOptions {
    pub cutoff: f64,
    pub pbc: [bool; 3],
    pub include_periodic_self_images: bool,
    pub budget: GraphBudget,
}

/// Both enumeration passes are included in the visit counters.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GraphBuildStats {
    pub atom_pair_visits: u64,
    pub image_candidate_visits: u64,
    pub edges: usize,
    /// Requested edge storage plus graph header and certified-cell workspace.
    /// This is a logical storage charge, not allocator RSS.
    pub requested_bytes: usize,
}

#[derive(Debug, thiserror::Error)]
pub enum GraphBuildError {
    #[error("invalid periodic graph input: {0}")]
    InvalidInput(String),
    #[error("uncertifiable periodic geometry: {0}")]
    NumericalGeometry(String),
    #[error("periodic graph image shift outside supported range")]
    ShiftRange,
    #[error("periodic graph budget {resource} exceeded: {limit}")]
    BudgetExceeded { resource: &'static str, limit: u64 },
    #[error("periodic graph allocation failed")]
    Allocation,
}

/// An all-image graph tied to the exact positions, atom order, identity, and cell.
///
/// End the borrow and rebuild after changing the snapshot or graph options.
/// Chemical arrays are not validated by this geometry interface.
///
/// The snapshot cannot be mutated while its graph remains in use:
///
/// ```compile_fail,E0502
/// use modena_core::{GraphBudget, PeriodicGraph, PeriodicGraphOptions, SimBox, System};
/// fn graph_borrow_prevents_mutation() {
/// let mut system = System::new(SimBox::cubic(2.0));
/// system.resize(1);
/// let graph = PeriodicGraph::build(&system, PeriodicGraphOptions {
///     cutoff: 1.0, pbc: [true; 3], include_periodic_self_images: true,
///     budget: GraphBudget { max_atom_pairs: 2, max_image_candidates: 100,
///         max_edges: 100, max_bytes: 10000, max_abs_shift: 100 },
/// }).unwrap();
/// system.atomic_numbers[0] = 8; // graph_borrow_prevents_mutation
/// println!("{:?}", graph.edges());
/// }
/// ```
pub struct PeriodicGraph<'a> {
    system: &'a System,
    options: PeriodicGraphOptions,
    edges: Vec<PeriodicEdge>,
    stats: GraphBuildStats,
}

impl<'a> PeriodicGraph<'a> {
    /// Build atomically, or return an explicit input, geometry, or resource error.
    pub fn build(
        system: &'a System,
        options: PeriodicGraphOptions,
    ) -> Result<Self, GraphBuildError> {
        Self::build_allocating(system, options, allocate_edges)
    }

    fn build_allocating(
        system: &'a System,
        options: PeriodicGraphOptions,
        allocate: impl FnOnce(usize, GraphBudget) -> Result<Vec<PeriodicEdge>, GraphBuildError>,
    ) -> Result<Self, GraphBuildError> {
        let cutoff2 = cutoff_squared(options.cutoff).map_err(geometry_error)?;
        checked_shift_endpoint(0.0, options.budget.max_abs_shift).map_err(geometry_error)?;
        let n = system.n();
        if n.checked_mul(3) != Some(system.positions.len())
            || system.positions.iter().any(|x| !x.is_finite())
        {
            return Err(GraphBuildError::InvalidInput(
                "positions must contain exactly 3*N finite coordinates".into(),
            ));
        }
        // Cell validation is mandatory even when the snapshot has no atoms.
        let cell = certified_cell(system.box_size).map_err(geometry_error)?;
        let pair_visits = charged_pairs(n, options.budget)?;
        let counted = enumerate(system, &cell, options, cutoff2, |_| {})?;
        let requested_bytes = storage_bytes(counted.edges, options.budget)?;
        let mut edges = allocate(counted.edges, options.budget)?;
        let emitted = enumerate(system, &cell, options, cutoff2, |edge| edges.push(edge))?;
        if emitted != counted || edges.len() != counted.edges {
            return Err(GraphBuildError::NumericalGeometry(
                "immutable enumeration passes disagree".into(),
            ));
        }
        edges.sort_unstable();
        // Capacity is checked after allocation and again at publication. Pushes
        // cannot grow the vector: pass one certified its exact accepted count.
        storage_bytes(edges.capacity(), options.budget)?;
        let stats = GraphBuildStats {
            atom_pair_visits: pair_visits,
            image_candidate_visits: counted.candidates,
            edges: counted.edges,
            requested_bytes,
        };
        Ok(Self {
            system,
            options,
            edges,
            stats,
        })
    }

    pub fn edges(&self) -> &[PeriodicEdge] {
        &self.edges
    }
    pub fn system(&self) -> &'a System {
        self.system
    }
    pub fn options(&self) -> PeriodicGraphOptions {
        self.options
    }
    pub fn stats(&self) -> GraphBuildStats {
        self.stats
    }
}

fn geometry_error(error: GeometryError) -> GraphBuildError {
    match error {
        GeometryError::Invalid => {
            GraphBuildError::InvalidInput("invalid cell, cutoff, or coordinates".into())
        }
        GeometryError::Uncertified => GraphBuildError::NumericalGeometry(
            "geometry bounds or Cartesian evaluation could not be certified".into(),
        ),
        GeometryError::ShiftRange => GraphBuildError::ShiftRange,
    }
}

fn budget_error(resource: &'static str, limit: u64) -> GraphBuildError {
    GraphBuildError::BudgetExceeded { resource, limit }
}
fn storage_limit(limit: usize) -> u64 {
    // Current Rust targets have usize <= u64. Preserve a bounded error value
    // if a wider target is introduced; work/capacity arithmetic remains checked.
    u64::try_from(limit).unwrap_or(u64::MAX)
}
fn charged_pairs(n: usize, budget: GraphBudget) -> Result<u64, GraphBuildError> {
    let error = || budget_error("atom_pairs", budget.max_atom_pairs);
    let n = u64::try_from(n).map_err(|_| error())?;
    let visits = u128::from(n)
        .checked_mul(u128::from(n))
        .and_then(|pairs| pairs.checked_mul(2))
        .ok_or_else(error)?;
    if visits > u128::from(budget.max_atom_pairs) {
        return Err(error());
    }
    u64::try_from(visits).map_err(|_| error())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Counts {
    candidates: u64,
    edges: usize,
}

fn span_product(bounds: [(i64, i64); 3], limit: u64) -> Result<u128, GraphBuildError> {
    if bounds.iter().any(|&(lo, hi)| lo > hi) {
        return Ok(0);
    }
    let mut product = 1u128;
    for (lo, hi) in bounds {
        let length = u128::try_from(i128::from(hi) - i128::from(lo) + 1)
            .map_err(|_| budget_error("image_candidates", limit))?;
        product = product
            .checked_mul(length)
            .ok_or_else(|| budget_error("image_candidates", limit))?;
    }
    Ok(product)
}

fn enumerate(
    system: &System,
    cell: &CertifiedCell,
    options: PeriodicGraphOptions,
    cutoff2: f64,
    mut emit: impl FnMut(PeriodicEdge),
) -> Result<Counts, GraphBuildError> {
    let budget = options.budget;
    let geometry = GeometryOptions {
        cutoff: options.cutoff,
        pbc: options.pbc,
        max_abs_shift: budget.max_abs_shift,
    };
    let mut counts = Counts {
        candidates: 0,
        edges: 0,
    };
    for source in 0..system.n() {
        let ri = system.pos(source);
        for destination in 0..system.n() {
            let rj = system.pos(destination);
            // image_bounds first runs the B-independent global certificate;
            // no local envelope or candidate enumeration can precede it.
            let bounds = image_bounds(cell, ri, rj, geometry).map_err(geometry_error)?;
            let candidates = span_product(bounds, budget.max_image_candidates)?;
            let charged = candidates
                .checked_mul(2)
                .and_then(|n| n.checked_add(u128::from(counts.candidates)))
                .ok_or_else(|| budget_error("image_candidates", budget.max_image_candidates))?;
            if charged > u128::from(budget.max_image_candidates) {
                return Err(budget_error(
                    "image_candidates",
                    budget.max_image_candidates,
                ));
            }
            counts.candidates = u64::try_from(charged)
                .map_err(|_| budget_error("image_candidates", budget.max_image_candidates))?;
            if candidates == 0 {
                continue;
            }
            for sx in bounds[0].0..=bounds[0].1 {
                for sy in bounds[1].0..=bounds[1].1 {
                    for sz in bounds[2].0..=bounds[2].1 {
                        let shift = [sx, sy, sz];
                        if source == destination
                            && (shift == [0; 3] || !options.include_periodic_self_images)
                        {
                            continue;
                        }
                        let d = edge_displacement(cell, ri, rj, shift).map_err(geometry_error)?;
                        if strict_radius(d, cutoff2).map_err(geometry_error)? {
                            counts.edges = counts.edges.checked_add(1).ok_or_else(|| {
                                budget_error("edges", storage_limit(budget.max_edges))
                            })?;
                            if counts.edges > budget.max_edges {
                                return Err(budget_error("edges", storage_limit(budget.max_edges)));
                            }
                            emit(PeriodicEdge {
                                source,
                                destination,
                                edge_cell_shift: shift,
                            });
                        }
                    }
                }
            }
        }
    }
    Ok(counts)
}

fn storage_bytes(capacity: usize, budget: GraphBudget) -> Result<usize, GraphBuildError> {
    let error = || budget_error("bytes", storage_limit(budget.max_bytes));
    let bytes = capacity
        .checked_mul(size_of::<PeriodicEdge>())
        .and_then(|n| n.checked_add(size_of::<PeriodicGraph<'_>>()))
        .and_then(|n| n.checked_add(size_of::<CertifiedCell>()))
        .ok_or_else(error)?;
    if bytes > budget.max_bytes {
        return Err(error());
    }
    Ok(bytes)
}
fn reserve_edges(count: usize) -> Result<Vec<PeriodicEdge>, ()> {
    let mut edges = Vec::new();
    edges.try_reserve_exact(count).map_err(|_| ())?;
    Ok(edges)
}
fn allocation_result(
    allocation: Result<Vec<PeriodicEdge>, ()>,
    count: usize,
    budget: GraphBudget,
) -> Result<Vec<PeriodicEdge>, GraphBuildError> {
    let edges = allocation.map_err(|_| GraphBuildError::Allocation)?;
    if !edges.is_empty() || edges.capacity() < count {
        return Err(GraphBuildError::Allocation);
    }
    storage_bytes(edges.capacity(), budget)?;
    Ok(edges)
}
fn allocate_edges(count: usize, budget: GraphBudget) -> Result<Vec<PeriodicEdge>, GraphBuildError> {
    allocation_result(reserve_edges(count), count, budget)
}
/// Inject only allocation failure, without requesting huge actual memory.
#[cfg(test)]
fn allocate_edges_with<'a>(
    system: &'a System,
    options: PeriodicGraphOptions,
    allocate: impl FnOnce(usize) -> Result<Vec<PeriodicEdge>, ()>,
) -> Result<PeriodicGraph<'a>, GraphBuildError> {
    PeriodicGraph::build_allocating(system, options, |count, budget| {
        allocation_result(allocate(count), count, budget)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SimBox;
    fn options() -> PeriodicGraphOptions {
        PeriodicGraphOptions {
            cutoff: 1.01,
            pbc: [true; 3],
            include_periodic_self_images: true,
            budget: GraphBudget {
                max_atom_pairs: 100,
                max_image_candidates: 1000,
                max_edges: 100,
                max_bytes: 10000,
                max_abs_shift: 100,
            },
        }
    }
    #[test]
    fn injected_allocation_failure_preserves_source_and_prior_graph() {
        let mut s = System::new(SimBox::cubic(1.0));
        s.resize(1);
        s.forces.fill(7.0);
        s.atomic_numbers[0] = 8;
        let prior = PeriodicGraph::build(&s, options()).unwrap();
        let edges = prior.edges().to_vec();
        let saved = (
            s.positions.clone(),
            s.forces.clone(),
            s.atomic_numbers.clone(),
        );
        assert!(matches!(
            allocate_edges_with(&s, options(), |_| Err(())),
            Err(GraphBuildError::Allocation)
        ));
        assert_eq!(prior.edges(), edges);
        assert_eq!(
            (&s.positions, &s.forces, &s.atomic_numbers),
            (&saved.0, &saved.1, &saved.2)
        );
    }
    #[test]
    fn allocator_capacity_is_charged_before_publication() {
        let mut s = System::new(SimBox::cubic(1.0));
        s.resize(1);
        let mut options = options();
        options.budget.max_bytes = storage_bytes(6, options.budget).unwrap();
        let result = allocate_edges_with(&s, options, |_| Ok(Vec::with_capacity(7)));
        assert!(matches!(
            result,
            Err(GraphBuildError::BudgetExceeded {
                resource: "bytes",
                ..
            })
        ));
    }
    #[test]
    fn empty_snapshot_charges_exact_header_and_workspace() {
        let s = System::new(SimBox::cubic(1.0));
        let graph = PeriodicGraph::build(&s, options()).unwrap();
        assert_eq!(
            graph.stats(),
            GraphBuildStats {
                atom_pair_visits: 0,
                image_candidate_visits: 0,
                edges: 0,
                requested_bytes: size_of::<PeriodicGraph<'_>>() + size_of::<CertifiedCell>()
            }
        );
    }
    #[test]
    fn empty_spans_and_extreme_integer_counts_are_checked() {
        assert_eq!(
            span_product([(3, 2), (i64::MIN, i64::MAX), (0, 0)], 100).unwrap(),
            0
        );
        assert_eq!(
            span_product([(i64::MIN, i64::MAX), (i64::MIN, i64::MAX), (3, 2)], 100).unwrap(),
            0
        );
        assert!(matches!(
            span_product([(i64::MIN, i64::MAX); 3], 100),
            Err(GraphBuildError::BudgetExceeded {
                resource: "image_candidates",
                limit: 100
            })
        ));
        assert!(matches!(
            charged_pairs(usize::MAX, options().budget),
            Err(GraphBuildError::BudgetExceeded {
                resource: "atom_pairs",
                limit: 100
            })
        ));
        assert!(matches!(
            storage_bytes(usize::MAX, options().budget),
            Err(GraphBuildError::BudgetExceeded {
                resource: "bytes",
                limit: 10000
            })
        ));
    }
}
