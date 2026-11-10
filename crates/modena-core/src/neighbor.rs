//! Cell-list neighbor search with a Verlet skin, for orthorhombic periodic
//! boxes; triclinic cells use guaranteed-correct all-pair enumeration with
//! the lattice minimum image, since a Cartesian bounding stencil is unsafe.
//!
//! Builds a *half* neighbor list (each pair `(i, j)` with `i < j` stored once)
//! and the force kernel applies Newton's third law, halving the pair work vs
//! the Phase-1 full list. A Verlet `skin` avoids rebuilding every step: the list
//! is rebuilt only when the largest per-atom displacement since the last build
//! exceeds `skin / 2`, so a well-chosen skin (`~0.3σ`) amortizes the grid build
//! over many steps. The force kernel still gates each pair by the true `cutoff`
//! (pairs within the skin but beyond the cutoff are skipped).
//!
//! M7: positions read from the flat SoA buffer; minimum-image dispatch goes
//! through `SimBox::min_image` (the triclinic seam).

use rayon::prelude::*;

use crate::system::System;
use crate::vec3::Vec3;

mod triclinic;
use triclinic::{build_triclinic, TriclinicBuildStrategy};

/// Requested triclinic builder, available only for explicit diagnostic builds.
#[cfg(feature = "neighbor-diagnostics")]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NeighborBuildStrategy {
    Auto,
    Scalar,
    FractionalStencil,
}

/// Safety/cost reason for falling back from a triclinic stencil.
#[cfg(feature = "neighbor-diagnostics")]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NeighborBuildFallbackReason {
    UncertifiedGeometry,
    FractionalBoundary,
    ImageDomain,
    GridCap,
    StencilCap,
    WorkCap,
    MemoryCap,
    AllBins,
    CandidateCost,
}

/// Actual successful triclinic rebuild, retained when the skin reuses a list.
#[cfg(feature = "neighbor-diagnostics")]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct NeighborBuildDiagnostics {
    pub selected: NeighborBuildStrategy,
    pub fallback: Option<NeighborBuildFallbackReason>,
    /// Number of minimum-image predicates evaluated, after 0/0 exclusions.
    pub candidate_checks: u64,
    /// Requested temporary workspace, excluding the emitted pair buffer.
    pub requested_workspace_bytes: usize,
}

#[cfg(feature = "neighbor-diagnostics")]
impl From<triclinic::BuildReport> for NeighborBuildDiagnostics {
    fn from(report: triclinic::BuildReport) -> Self {
        use triclinic::FallbackReason as F;
        use NeighborBuildFallbackReason as D;
        Self {
            selected: match report.selected {
                TriclinicBuildStrategy::Auto => NeighborBuildStrategy::Auto,
                TriclinicBuildStrategy::Scalar => NeighborBuildStrategy::Scalar,
                TriclinicBuildStrategy::FractionalStencil => {
                    NeighborBuildStrategy::FractionalStencil
                }
            },
            fallback: report.fallback.map(|reason| match reason {
                F::UncertifiedGeometry => D::UncertifiedGeometry,
                F::FractionalBoundary => D::FractionalBoundary,
                F::ImageDomain => D::ImageDomain,
                F::GridCap => D::GridCap,
                F::StencilCap => D::StencilCap,
                F::WorkCap => D::WorkCap,
                F::MemoryCap => D::MemoryCap,
                F::AllBins => D::AllBins,
                F::CandidateCost => D::CandidateCost,
            }),
            candidate_checks: report.candidate_checks,
            requested_workspace_bytes: report.requested_workspace_bytes,
        }
    }
}

pub struct NeighborList {
    cutoff: f64,
    skin: f64,
    /// List-building cutoff = cutoff + skin.
    list_cutoff: f64,
    list_cutoff2: f64,
    /// Half pairs `(i, j)` with `i < j` within `list_cutoff`.
    pairs: Vec<(usize, usize)>,
    /// Reference positions at the last rebuild (for the skin displacement test).
    ref_positions: Vec<f64>,
    /// Full box at the last rebuild. A barostat can rescale the box without
    /// any atom moving: the grid dimensions were computed from the old box, and
    /// after a shrink `cell_size < list_cutoff` would silently miss pairs. A
    /// box change therefore forces a rebuild.
    ref_box: Option<crate::system::SimBox>,
    built: bool,
    /// Number of grid rebuilds performed (for tests + diagnostics).
    rebuild_count: usize,
    /// Indexed independently weighted special pairs; 0/0 pairs are omitted.
    exclusions: PairScaleSet,
    #[cfg(feature = "neighbor-diagnostics")]
    triclinic_strategy: TriclinicBuildStrategy,
    #[cfg(feature = "neighbor-diagnostics")]
    triclinic_diagnostics: Option<NeighborBuildDiagnostics>,
    // Grid state (rebuilt as needed).
    dims: [usize; 3],
    cell_size: [f64; 3],
    head: Vec<i32>, // first atom index in cell, -1 if empty
    next: Vec<i32>, // linked list of atoms within a cell
}

/// Independently weighted nonbonded pair. Missing pairs use (1, 1).
/// Indices are normalized by constructors. Finite nonnegative weights are
/// accepted (including amplification > 1); self pairs and conflicting duplicates
/// are errors. Atom bounds are validated against the live System by providers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairScale {
    pub i: usize,
    pub j: usize,
    pub lj: f64,
    pub coulomb: f64,
}

/// Indexed normalized policy; an empty policy has no per-pair hashing cost.
#[derive(Default, Clone)]
pub(crate) struct PairScaleSet {
    entries: Vec<PairScale>,
    index: std::collections::HashMap<(usize, usize), usize>,
}

impl PairScaleSet {
    pub(crate) fn new(scales: &[PairScale]) -> crate::error::Result<Self> {
        use crate::error::ModenaError;
        let mut set = Self::default();
        for &scale in scales {
            let (i, j) = if scale.i < scale.j {
                (scale.i, scale.j)
            } else {
                (scale.j, scale.i)
            };
            if i == j
                || !scale.lj.is_finite()
                || !scale.coulomb.is_finite()
                || scale.lj < 0.0
                || scale.coulomb < 0.0
            {
                return Err(ModenaError::Config(
                    "pair scales need distinct indices and finite nonnegative LJ/Coulomb weights"
                        .into(),
                ));
            }
            let scale = PairScale { i, j, ..scale };
            if let Some(old) = set.get(i, j) {
                if old != scale {
                    return Err(ModenaError::Config(format!(
                        "conflicting scales for pair ({i}, {j})"
                    )));
                }
            } else {
                set.index.insert((i, j), set.entries.len());
                set.entries.push(scale);
            }
        }
        Ok(set)
    }
    pub(crate) fn get(&self, i: usize, j: usize) -> Option<PairScale> {
        if self.entries.is_empty() {
            return None;
        }
        let key = if i < j { (i, j) } else { (j, i) };
        self.index.get(&key).map(|&idx| self.entries[idx])
    }
    pub(crate) fn scale(&self, i: usize, j: usize) -> (f64, f64) {
        self.get(i, j).map_or((1.0, 1.0), |p| (p.lj, p.coulomb))
    }
    fn contains(&self, i: usize, j: usize) -> bool {
        self.get(i, j)
            .is_some_and(|p| p.lj == 0.0 && p.coulomb == 0.0)
    }
    pub(crate) fn entries(&self) -> &[PairScale] {
        &self.entries
    }
    pub(crate) fn validate(&self, n: usize) -> crate::error::Result<()> {
        if self.entries.iter().any(|p| p.j >= n) {
            return Err(crate::error::ModenaError::Config(
                "pair-scale index outside System atom range".into(),
            ));
        }
        Ok(())
    }
}

/// Atom counts at which the cell-list pair scan parallelizes over `i`-ranges.
/// Below this the rayon overhead exceeds the scan (mirrors the pair-pass
/// `SERIAL_PAIR_THRESHOLD` pattern in `force_model.rs`).
const NEIGHBOR_PAR_MIN: usize = 4096;

impl NeighborList {
    /// New list with force cutoff `cutoff` and Verlet `skin` (skin = 0 disables
    /// the Verlet optimization — rebuild every `build` call, matching Phase 1).
    pub fn with_skin(cutoff: f64, skin: f64) -> Self {
        Self::with_skin_exclusions(cutoff, skin, &[])
    }

    /// Like [`with_skin`](Self::with_skin), with bonded-exclusion pairs that
    /// the pair kernel must never see (`exclusions` are normalized to
    /// `(min, max)` internally).
    pub fn with_skin_exclusions(cutoff: f64, skin: f64, exclusions: &[(u32, u32)]) -> Self {
        let scales: Vec<PairScale> = exclusions
            .iter()
            .map(|&(i, j)| PairScale {
                i: i as usize,
                j: j as usize,
                lj: 0.0,
                coulomb: 0.0,
            })
            .collect();
        Self::with_pair_scales(cutoff, skin, &scales)
            .expect("invalid neighbor/exclusion parameters")
    }

    /// Construct a validated independent LJ/Coulomb policy. Only 0/0 pairs
    /// are removed from the list, so one disabled interaction cannot hide the other.
    pub fn with_pair_scales(
        cutoff: f64,
        skin: f64,
        scales: &[PairScale],
    ) -> crate::error::Result<Self> {
        if !cutoff.is_finite()
            || cutoff <= 0.0
            || !skin.is_finite()
            || skin < 0.0
            || !(cutoff + skin).is_finite()
        {
            return Err(crate::error::ModenaError::Config(
                "neighbor cutoff must be finite positive and skin finite nonnegative".into(),
            ));
        }
        let list_cutoff = cutoff + skin;
        let excl = PairScaleSet::new(scales)?;
        Ok(Self {
            cutoff,
            skin,
            list_cutoff,
            list_cutoff2: list_cutoff * list_cutoff,
            pairs: Vec::new(),
            ref_positions: Vec::new(),
            ref_box: None,
            built: false,
            rebuild_count: 0,
            exclusions: excl,
            #[cfg(feature = "neighbor-diagnostics")]
            triclinic_strategy: TriclinicBuildStrategy::Auto,
            #[cfg(feature = "neighbor-diagnostics")]
            triclinic_diagnostics: None,
            dims: [1, 1, 1],
            cell_size: [cutoff; 3],
            head: Vec::new(),
            next: Vec::new(),
        })
    }

    /// Phase-1-compatible constructor: no skin (rebuild every call).
    pub fn new(cutoff: f64) -> Self {
        Self::with_skin(cutoff, 0.0)
    }

    /// Select a diagnostic request. Auto remains Scalar. Forced FractionalStencil
    /// attempts the certified stencil, with safety/cost fallbacks to Scalar;
    /// reports always name the actual builder. Changing the request invalidates
    /// the cached list so the next build honors it.
    #[cfg(feature = "neighbor-diagnostics")]
    pub fn with_triclinic_strategy_for_diagnostics(
        mut self,
        strategy: NeighborBuildStrategy,
    ) -> Self {
        self.triclinic_strategy = match strategy {
            NeighborBuildStrategy::Auto => TriclinicBuildStrategy::Auto,
            NeighborBuildStrategy::Scalar => TriclinicBuildStrategy::Scalar,
            NeighborBuildStrategy::FractionalStencil => TriclinicBuildStrategy::FractionalStencil,
        };
        self.built = false;
        self.triclinic_diagnostics = None;
        self
    }

    /// Last successful triclinic rebuild; None before a build, after an empty
    /// build or failed bounded rebuild, or after an orthorhombic rebuild.
    #[cfg(feature = "neighbor-diagnostics")]
    pub fn triclinic_build_diagnostics(&self) -> Option<NeighborBuildDiagnostics> {
        self.triclinic_diagnostics
    }

    pub fn pair_scale(&self, i: usize, j: usize) -> (f64, f64) {
        self.exclusions.scale(i, j)
    }

    pub(crate) fn scales(&self) -> &PairScaleSet {
        &self.exclusions
    }

    pub fn pair_scales(&self) -> &[PairScale] {
        self.exclusions.entries()
    }

    pub fn validate_pair_scales(&self, n: usize) -> crate::error::Result<()> {
        self.exclusions.validate(n)
    }

    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    pub fn skin(&self) -> f64 {
        self.skin
    }

    /// The half pair list `(i, j)` with `i < j`.
    pub fn pairs(&self) -> &[(usize, usize)] {
        &self.pairs
    }

    /// Number of stored half-pairs.
    pub fn n_pairs(&self) -> usize {
        self.pairs.len()
    }

    /// Skin-aware rebuild. Rebuilds the grid + pair list only when the largest
    /// per-atom displacement since the last build exceeds `skin / 2`; otherwise
    /// reuses the current list (the force kernel recomputes distances from the
    /// live positions each step, so a stale *list* is correct as long as no pair
    /// crossed the `list_cutoff` boundary).
    pub fn build(&mut self, system: &System) {
        let n = system.n();
        if n == 0 {
            self.built = false;
            self.pairs.clear();
            #[cfg(feature = "neighbor-diagnostics")]
            {
                self.triclinic_diagnostics = None;
            }
            return;
        }
        if self.built
            && self.skin > 0.0
            && self.state_is_current(system)
            && !self.displacement_exceeds(system)
        {
            return;
        }
        self.rebuild(system);
    }

    /// Whether the cached grid still matches the system shape: same atom
    /// count and same box lengths. A barostat rescale changes the box without
    /// moving atoms (the displacement test would pass), but the grid
    /// dimensions and cell sizes were computed from the old box — after a
    /// shrink `cell_size < list_cutoff` and the ±1 stencil silently misses
    /// pairs. A changed atom count would index `ref_positions` out of bounds.
    fn state_is_current(&self, system: &System) -> bool {
        self.ref_positions.len() == 3 * system.n() && self.ref_box == Some(system.box_size)
    }

    /// Force a rebuild regardless of skin (used by tests and the initial build).
    pub fn force_rebuild(&mut self, system: &System) {
        self.rebuild(system);
    }

    /// Total grid rebuilds performed since construction (tests + diagnostics).
    pub fn rebuild_count(&self) -> usize {
        self.rebuild_count
    }

    fn displacement_exceeds(&self, system: &System) -> bool {
        let half_skin = 0.5 * self.skin;
        let half_skin2 = half_skin * half_skin;
        let n = system.n();
        let box_size = system.box_size;
        for i in 0..n {
            // Minimum-image displacement: the integrator wraps positions into
            // the box each step, so an atom that crossed a periodic face shows
            // a raw displacement of ~L. Comparing raw positions would flag it
            // as "moved" forever after, forcing a rebuild every step and
            // defeating the skin. The minimum-image displacement is the true
            // physical displacement.
            let d = box_size.min_image(Vec3::new(
                system.positions[3 * i] - self.ref_positions[3 * i],
                system.positions[3 * i + 1] - self.ref_positions[3 * i + 1],
                system.positions[3 * i + 2] - self.ref_positions[3 * i + 2],
            ));
            if d.norm2() > half_skin2 {
                return true;
            }
        }
        false
    }

    fn rebuild(&mut self, system: &System) {
        self.rebuild_impl::<false>(system, usize::MAX)
            .expect("unbounded pair scan");
    }

    /// Snapshot-only: stream distance candidates into one capped pair
    /// vector. No speculative N*27 reserve or per-thread staging vectors.
    pub(crate) fn rebuild_bounded(
        &mut self,
        system: &System,
        max_pairs: usize,
    ) -> crate::Result<()> {
        self.built = false;
        if self.pairs.capacity() > max_pairs {
            self.pairs = Vec::new();
        }
        self.rebuild_impl::<true>(system, max_pairs)
    }

    fn rebuild_impl<const BOUNDED: bool>(
        &mut self,
        system: &System,
        max_pairs: usize,
    ) -> crate::Result<()> {
        let n = system.n();
        let box_size = system.box_size;
        self.rebuild_count += 1;
        #[cfg(feature = "neighbor-diagnostics")]
        {
            self.triclinic_diagnostics = None;
        }

        // Grid dimensions: cell edge >= list_cutoff. Clamp to 1 cell per axis if
        // the box is smaller than the list cutoff (valid via minimum-image).
        let mut dims = [1usize; 3];
        let mut cell_size = [self.list_cutoff; 3];
        for k in 0..3 {
            let b = box_size.size.0[k];
            let d = (b / self.list_cutoff).floor() as usize;
            let d = d.max(1);
            dims[k] = d;
            cell_size[k] = b / d as f64;
        }
        self.dims = dims;
        self.cell_size = cell_size;

        // Triclinic cell faces need a geometry-safe search: Cartesian wrapped
        // bins can miss partners in nonadjacent bounding-box cells. Enumerate
        // all pairs for these cells. The orthorhombic path remains accelerated.
        // Small-box fast path (committee finding K): when every axis has ≤ 3
        // cells the ±1 stencil (with periodic wrap) covers the entire box, so
        // the binning buys nothing — it became a ~11× overscan of distance
        // checks plus a cache- hostile linked-list chase (measured: the grid
        // rebuild was ~32% of a 1-thread step at benchmark sizes). A plain
        // contiguous j>i double loop visits the same candidate set
        // (∏min(d,3)·N/… ≈ N²/2 here) with dense sequential access and no
        // per-cell indirection. Same pair set, same min-image predicate.
        let mut pairs = std::mem::take(&mut self.pairs);
        pairs.clear();
        if box_size.is_triclinic() || dims.iter().all(|&d| d <= 3) {
            #[cfg(feature = "neighbor-diagnostics")]
            let strategy = self.triclinic_strategy;
            #[cfg(not(feature = "neighbor-diagnostics"))]
            let strategy = TriclinicBuildStrategy::Auto;
            let report = build_triclinic::<BOUNDED>(
                system,
                self.list_cutoff2,
                &self.exclusions,
                max_pairs,
                strategy,
                &mut pairs,
            )?;
            #[cfg(feature = "neighbor-diagnostics")]
            if box_size.is_triclinic() {
                self.triclinic_diagnostics = Some(report.into());
            }
            #[cfg(not(feature = "neighbor-diagnostics"))]
            let _ = report;
        } else {
            let ncells = dims[0] * dims[1] * dims[2];
            self.head.clear();
            self.head.resize(ncells, -1);
            self.next.clear();
            self.next.resize(n, -1);

            // Bin atoms into cells.
            for i in 0..n {
                let c = self.cell_of(system.pos(i));
                self.next[i] = self.head[c];
                self.head[c] = i as i32;
            }

            // Deduplicated per-axis stencil (see Phase-1 comment: with < 3 cells per
            // axis the raw {-1,0,1} stencil wraps onto the same cell and overcounts).
            let dxs = stencil_deltas(dims[0]);
            let dys = stencil_deltas(dims[1]);
            let dzs = stencil_deltas(dims[2]);

            // The env knob forces serial MD scans for benchmarking. Snapshots
            // always use one buffer to enforce the cap across the whole build.
            let force_serial = std::env::var("MODENA_SERIAL_NEIGHBOR").is_ok();
            if !BOUNDED && n >= NEIGHBOR_PAR_MIN && !force_serial {
                // Parallel over `i`-ranges, concatenated in order: rayon's indexed
                // `par_chunks().map().collect()` preserves per-chunk order and the
                // chunk sequence, and the per-`i` scan is unchanged, so the pair
                // list is byte-identical to the serial build — which keeps the
                // force-accumulation order (and hence runnable determinism)
                // unchanged. The division follows the serial threshold so small
                // rebuilds (NVE at benchmark sizes) keep today's exact path.
                let idx: Vec<usize> = (0..n).collect();
                let chunk = n.div_ceil(rayon::current_num_threads().max(2));
                let buffers: Vec<Vec<(usize, usize)>> = idx
                    .par_chunks(chunk)
                    .map(|is| {
                        let mut buf: Vec<(usize, usize)> = Vec::new();
                        for &i in is {
                            self.scan_atom::<false>(
                                system,
                                i,
                                [&dxs, &dys, &dzs],
                                &mut buf,
                                usize::MAX,
                            )
                            .expect("unbounded pair scan");
                        }
                        buf
                    })
                    .collect();
                for buf in buffers {
                    pairs.extend(buf);
                }
            } else {
                for i in 0..n {
                    self.scan_atom::<BOUNDED>(
                        system,
                        i,
                        [&dxs, &dys, &dzs],
                        &mut pairs,
                        max_pairs,
                    )?;
                }
            }
        }
        self.pairs = pairs;

        // Snapshot reference positions and box for the skin test.
        self.ref_positions.clear();
        self.ref_positions.extend_from_slice(system.positions_raw());
        self.ref_box = Some(box_size);
        self.built = true;
        Ok(())
    }

    /// Scan the ±1-stencil neighbor cells of atom `i` and push every qualifying
    /// half-list pair `(i, j)` with `j > i`, skipping exclusions, into `out`.
    /// The shared core of the serial and parallel builds — the pair order (and
    /// thus the force-accumulation order) is identical in both.
    #[inline]
    fn scan_atom<const BOUNDED: bool>(
        &self,
        system: &System,
        i: usize,
        [dxs, dys, dzs]: [&[i32]; 3],
        out: &mut Vec<(usize, usize)>,
        max_pairs: usize,
    ) -> crate::Result<()> {
        let dims = self.dims;
        let lc2 = self.list_cutoff2;
        let box_size = system.box_size;
        let ri = system.pos(i);
        let ci = self.cell_of(ri);
        let (ix, iy, iz) = self.unravel(ci);
        for &dz in dzs {
            for &dy in dys {
                for &dx in dxs {
                    let nx = wrap_cell(ix as i32 + dx, dims[0]);
                    let ny = wrap_cell(iy as i32 + dy, dims[1]);
                    let nz = wrap_cell(iz as i32 + dz, dims[2]);
                    let c = nx * dims[1] * dims[2] + ny * dims[2] + nz;
                    let mut p = self.head[c];
                    while p >= 0 {
                        let j = p as usize;
                        // Half list: only store i < j, skipping exclusions.
                        if j > i && !self.exclusions.contains(i, j) {
                            let rj = system.pos(j);
                            let d = box_size.min_image(ri - rj);
                            if d.norm2() < lc2 {
                                push_pair::<BOUNDED>(out, (i, j), max_pairs)?;
                            }
                        }
                        p = self.next[j];
                    }
                }
            }
        }
        Ok(())
    }
    #[inline]
    fn cell_of(&self, r: Vec3) -> usize {
        let dims = self.dims;
        let cs = self.cell_size;
        let ix = ((r.0[0] / cs[0]).floor() as i32).rem_euclid(dims[0] as i32) as usize;
        let iy = ((r.0[1] / cs[1]).floor() as i32).rem_euclid(dims[1] as i32) as usize;
        let iz = ((r.0[2] / cs[2]).floor() as i32).rem_euclid(dims[2] as i32) as usize;
        ix * dims[1] * dims[2] + iy * dims[2] + iz
    }

    #[inline]
    fn unravel(&self, c: usize) -> (usize, usize, usize) {
        let dims = self.dims;
        let iz = c % dims[2];
        let iy = (c / dims[2]) % dims[1];
        let ix = c / (dims[1] * dims[2]);
        (ix, iy, iz)
    }
}

#[inline]
fn push_pair<const BOUNDED: bool>(
    pairs: &mut Vec<(usize, usize)>,
    pair: (usize, usize),
    max_pairs: usize,
) -> crate::Result<()> {
    if BOUNDED && pairs.len() >= max_pairs {
        return Err(crate::ModenaError::Input(format!(
            "snapshot exceeds {max_pairs} neighbor pairs"
        )));
    }
    if BOUNDED && pairs.len() == pairs.capacity() {
        // Never let Vec::push choose an uncapped geometric growth request.
        let next = pairs.len().saturating_mul(2).max(4).min(max_pairs);
        pairs.try_reserve_exact(next - pairs.len()).map_err(|_| {
            crate::ModenaError::Input("cannot allocate snapshot neighbor pairs".into())
        })?;
    }
    pairs.push(pair);
    Ok(())
}

#[inline]
fn wrap_cell(idx: i32, dim: usize) -> usize {
    let dim = dim as i32;
    (((idx % dim) + dim) % dim) as usize
}

/// Distinct neighbor-cell offsets for one axis of length `dim`, deduplicated
/// against periodic wraparound.
#[inline]
fn stencil_deltas(dim: usize) -> Vec<i32> {
    match dim {
        0 | 1 => vec![0],
        2 => vec![0, 1],
        _ => vec![-1, 0, 1],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::SimBox;

    fn sys_with(positions: &[Vec3], box_edge: f64) -> System {
        let mut sys = System::new(SimBox::cubic(box_edge));
        sys.resize(positions.len());
        for (i, &p) in positions.iter().enumerate() {
            sys.set_pos(i, p);
        }
        sys
    }

    #[test]
    fn triclinic_bounded_zero_skin_matches_scalar() {
        let mut sys = sys_with(
            &[Vec3::zero(), Vec3::new(1., 0., 0.), Vec3::new(9.5, 0., 0.)],
            10.,
        );
        sys.box_size = SimBox::triclinic([
            Vec3::new(10., 0., 0.),
            Vec3::new(3., 10., 0.),
            Vec3::new(0., 0., 10.),
        ]);
        let mut bounded = NeighborList::new(2.);
        let mut scalar = NeighborList::new(2.);
        for x in [1., 4.] {
            sys.set_pos(1, Vec3::new(x, 0., 0.));
            bounded.rebuild_bounded(&sys, 3).unwrap();
            scalar.build(&sys);
            assert_eq!(bounded.pairs(), scalar.pairs());
            assert!(bounded.pairs.capacity() <= 3);
        }
        assert_eq!(bounded.pairs(), &[(0, 2)]);
        assert_eq!(bounded.rebuild_count(), 2);
    }

    #[test]
    fn triclinic_bounded_overflow_invalidates_and_recovers() {
        let mut sys = sys_with(
            &[Vec3::zero(), Vec3::new(1., 0., 0.), Vec3::new(9.5, 0., 0.)],
            10.,
        );
        sys.box_size = SimBox::triclinic([
            Vec3::new(10., 0., 0.),
            Vec3::new(3., 10., 0.),
            Vec3::new(0., 0., 10.),
        ]);
        for skin in [0., 0.3] {
            let mut nl = NeighborList::with_skin(2., skin);
            nl.rebuild_bounded(&sys, 3).unwrap();
            assert_eq!(nl.pairs(), &[(0, 1), (0, 2), (1, 2)]);
            #[cfg(feature = "neighbor-diagnostics")]
            assert!(nl.triclinic_build_diagnostics().is_some());
            let reference_positions = nl.ref_positions.clone();
            let error = nl.rebuild_bounded(&sys, 2).unwrap_err();
            assert!(error.to_string().contains("2 neighbor pairs"));
            assert!(!nl.built);
            #[cfg(feature = "neighbor-diagnostics")]
            assert!(nl.triclinic_build_diagnostics().is_none());
            assert!(nl.pairs().is_empty());
            assert!(nl.pairs.capacity() <= 2);
            assert_eq!(nl.ref_positions, reference_positions);
            assert_eq!(nl.ref_box, Some(sys.box_size));
            assert_eq!(nl.rebuild_count(), 2);
            nl.build(&sys);
            assert_eq!(nl.rebuild_count(), 3);
            assert_eq!(nl.pairs(), &[(0, 1), (0, 2), (1, 2)]);
        }
    }

    #[test]
    fn finds_neighbors_in_cubic_box() {
        let sys = sys_with(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)], 10.0);
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        // Half list: one pair (0,1).
        assert_eq!(nl.pairs(), &[(0, 1)]);
    }

    #[test]
    fn bounded_build_limits_capacity_and_recovers_after_dense_failure() {
        // Exercise both the dense-loop and binned stencil, including PBC
        // and an excluded contact. The three remaining pairs fit exactly.
        for edge in [6., 20.] {
            let sys = sys_with(
                &[
                    Vec3::new(edge - 0.1, 0., 0.),
                    Vec3::new(0.1, 0., 0.),
                    Vec3::new(0.4, 0., 0.),
                    Vec3::new(0.7, 0., 0.),
                ],
                edge,
            );
            let mut nl = NeighborList::with_skin_exclusions(2., 0.1, &[(0, 2), (0, 3), (1, 3)]);
            nl.rebuild_bounded(&sys, 3).unwrap();
            let expected = nl.pairs().to_vec();
            let mut sorted = expected.clone();
            sorted.sort_unstable();
            assert_eq!(sorted, [(0, 1), (1, 2), (2, 3)]);
            assert!(nl.pairs.capacity() <= 3);
            assert!(nl.rebuild_bounded(&sys, 2).is_err());
            assert!(nl.pairs().is_empty());
            assert!(nl.pairs.capacity() <= 2);
            nl.build(&sys); // same box/positions and positive skin: cannot reuse partial data
            assert_eq!(nl.pairs(), expected);
        }
    }

    #[test]
    fn bounded_build_accepts_sparse_system_with_zero_pair_budget() {
        let points: Vec<_> = (0..100).map(|i| Vec3::new(3. * i as f64, 0., 0.)).collect();
        let sys = sys_with(&points, 300.);
        let mut nl = NeighborList::with_skin(2., 0.5);
        nl.rebuild_bounded(&sys, 0).unwrap();
        assert!(nl.pairs().is_empty());
        assert_eq!(nl.pairs.capacity(), 0);
    }

    #[test]
    fn respects_cutoff() {
        let sys = sys_with(
            &[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(5.0, 0.0, 0.0),
            ],
            10.0,
        );
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        assert_eq!(nl.pairs(), &[(0, 1)]);
    }

    #[test]
    fn finds_periodic_image_neighbor() {
        let sys = sys_with(&[Vec3::new(0.1, 0.0, 0.0), Vec3::new(9.9, 0.0, 0.0)], 10.0);
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        assert_eq!(nl.pairs(), &[(0, 1)]);
    }

    #[test]
    fn no_duplicate_pairs_with_coarse_grid() {
        let sys = sys_with(
            &[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(2.0, 0.0, 0.0),
            ],
            5.0,
        );
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        // Pairs: (0,1), (0,2), (1,2) — each once.
        let mut sorted = nl.pairs().to_vec();
        sorted.sort();
        assert_eq!(sorted, vec![(0, 1), (0, 2), (1, 2)]);
    }

    #[test]
    fn skin_skips_rebuild_until_threshold() {
        // skin = 1.0 → rebuild only when displacement > 0.5.
        let mut sys = sys_with(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)], 10.0);
        let mut nl = NeighborList::with_skin(2.5, 1.0);
        nl.build(&sys);
        assert_eq!(nl.pairs(), &[(0, 1)]);

        // Small move (< 0.5): list reused, no rebuild. Nudge atom 1.
        sys.set_pos(1, Vec3::new(1.2, 0.0, 0.0));
        nl.build(&sys);
        assert_eq!(nl.pairs(), &[(0, 1)]);

        // Large move (> 0.5): rebuild. Atom 1 now at 6.0 → out of cutoff.
        sys.set_pos(1, Vec3::new(6.0, 0.0, 0.0));
        nl.build(&sys);
        assert!(nl.pairs().is_empty());
    }

    #[test]
    fn skin_survives_periodic_wrapping() {
        // Regression: the displacement test must use minimum-image deltas.
        // The integrator wraps positions into the box every step, so an atom
        // crossing a periodic face has a raw displacement of ~L but a physical
        // displacement of ~0. With raw comparison the skin would trigger a
        // rebuild every step forever after — defeating the Verlet skin.
        let mut sys = sys_with(&[Vec3::new(0.05, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)], 10.0);
        let mut nl = NeighborList::with_skin(2.5, 1.0); // rebuild only > 0.5
        nl.build(&sys);
        assert_eq!(nl.rebuild_count(), 1);
        assert_eq!(nl.pairs(), &[(0, 1)]);

        // Atom 0 drifts to -0.05 and is wrapped to 9.95: physical move 0.1,
        // raw jump 9.9. The list must be reused (no rebuild), and the pair
        // must survive via the periodic image.
        sys.set_pos(0, Vec3::new(-0.05, 0.0, 0.0));
        sys.wrap_positions();
        assert_eq!(sys.pos(0).0[0], 9.95);
        nl.build(&sys);
        assert_eq!(
            nl.rebuild_count(),
            1,
            "periodic wrapping must not trigger a rebuild"
        );
        assert_eq!(nl.pairs(), &[(0, 1)]);

        // A real physical move > skin/2 across the boundary still rebuilds.
        sys.set_pos(0, Vec3::new(-0.7, 0.0, 0.0));
        sys.wrap_positions();
        nl.build(&sys);
        assert_eq!(nl.rebuild_count(), 2, "physical move > skin/2 must rebuild");
    }

    #[test]
    fn skin_rebuilds_on_box_rescale_and_count_change() {
        // A barostat rescale changes the box without moving atoms. The
        // displacement test would pass, but the stale grid was computed from
        // the old box — after a shrink the ±1 stencil silently misses pairs,
        // so a box change must force a rebuild. Likewise an atom-count change
        // (which would otherwise index ref_positions out of bounds).
        let mut sys = sys_with(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)], 10.0);
        let mut nl = NeighborList::with_skin(2.5, 1.0);
        nl.build(&sys);
        assert_eq!(nl.rebuild_count(), 1);

        // Isotropic shrink 10 -> 9, atoms unmoved: must rebuild.
        sys.box_size = SimBox::cubic(9.0);
        nl.build(&sys);
        assert_eq!(
            nl.rebuild_count(),
            2,
            "box rescale without atom motion must rebuild"
        );
        assert_eq!(nl.pairs(), &[(0, 1)]);

        // Same box, same positions: no rebuild.
        nl.build(&sys);
        assert_eq!(nl.rebuild_count(), 2);

        // Atom count change: must rebuild (no ref_positions OOB).
        sys.resize(3);
        nl.build(&sys);
        assert_eq!(nl.rebuild_count(), 3, "atom-count change must rebuild");
    }

    #[test]
    fn exclusions_remove_pairs_in_both_paths() {
        // Three atoms in a 10-box: (0,1) excluded. Both the dims<=3 brute path
        // and the binned path must skip it.
        let sys = sys_with(
            &[
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.5, 0.0, 0.0),
                Vec3::new(2.0, 0.0, 0.0),
            ],
            10.0,
        );
        // dims = floor(10 / 2.5) = 4 -> binned path.
        let mut nl = NeighborList::with_skin_exclusions(2.5, 0.0, &[(0, 1)]);
        nl.build(&sys);
        let mut sorted = nl.pairs().to_vec();
        sorted.sort();
        assert_eq!(sorted, vec![(0, 2), (1, 2)], "binned path exclusion");
        // Small box -> dims <= 3 brute path; reversed declaration must normalize.
        let mut nl2 = NeighborList::with_skin_exclusions(2.5, 0.0, &[(1, 0)]);
        let mut sys2 = sys.clone();
        sys2.box_size = SimBox::cubic(6.0);
        nl2.build(&sys2);
        assert!(
            !nl2.pairs()
                .iter()
                .any(|&(a, b)| (a == 0 && b == 1) || (a == 1 && b == 0)),
            "brute-path exclusion must hold regardless of normalization; got {:?}",
            nl2.pairs()
        );
    }

    #[test]
    fn skin_list_yields_identical_forces() {
        // Physics equivalence: for as long as no atom moves more than skin/2
        // between rebuilds, a skinned (possibly stale) list must give the same
        // forces/PE as a fresh skinless list — the pair sets within the force
        // cutoff coincide; the extra skin pairs are gated out by the kernel.
        // Not bit-identical (a stale list can enumerate pairs in a different
        // order → different FP summation order), so compare with a tight
        // tolerance.
        use crate::{ForceModel, LennardJones, PairForceModel, Virial};

        let mut sys0 = crate::lattice::fcc(1.6796, [3, 3, 3], "Ar");
        // Deterministic small "velocities" (well below skin/2 per step).
        let vels: Vec<Vec3> = (0..sys0.n())
            .map(|i| {
                let f = i as f64;
                Vec3::new(
                    0.02 * f.sin(),
                    0.02 * (1.7 * f).cos(),
                    0.02 * (0.9 * f).sin(),
                )
            })
            .collect();
        let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());

        let mut bare = NeighborList::new(2.5); // skin 0: rebuild every call
        let mut skinned = NeighborList::with_skin(2.5, 0.6); // rebuild > 0.3

        for step in 0..10 {
            for (i, &v) in vels.iter().enumerate() {
                sys0.set_pos(i, sys0.pos(i) + v);
            }
            sys0.wrap_positions();
            let mut s_bare = sys0.clone();
            let mut s_skin = sys0.clone();
            s_bare.zero_forces();
            s_skin.zero_forces();
            let mut w1 = Virial::default();
            let mut w2 = Virial::default();
            let pe1 = model.compute(&mut s_bare, &mut w1, &mut bare).unwrap();
            let pe2 = model.compute(&mut s_skin, &mut w2, &mut skinned).unwrap();
            assert!(
                (pe1 - pe2).abs() < 1e-9 * (1.0 + pe1.abs()),
                "step {step}: PE {pe1} vs {pe2}"
            );
            for i in 0..sys0.n() {
                let df = s_bare.force(i) - s_skin.force(i);
                assert!(
                    df.norm2() < 1e-16,
                    "step {step} atom {i}: force delta {:?}",
                    df.norm2()
                );
            }
        }
        // The skinned list must actually have amortized: after 10 steps of
        // ~0.02/step (cumulative ~0.2 < 0.3), at most the initial build + one
        // rebuild — far fewer than the skinless 11 builds.
        assert!(skinned.rebuild_count() <= 2, "{}", skinned.rebuild_count());
        assert!(bare.rebuild_count() >= 10);
    }

    /// The parallel cell-list rebuild (n ≥ NEIGHBOR_PAR_MIN) must be
    /// deterministic (two rebuilds byte-identical — this is what keeps the
    /// force-accumulation order fixed) and produce exactly the brute-force
    /// half-list with the same exclusions applied.
    #[test]
    fn parallel_rebuild_deterministic_and_complete() {
        use crate::system::SimBox;
        use rand::rngs::StdRng;
        use rand::RngExt;
        use rand::SeedableRng;

        let n = 4096usize; // exactly the parallel boundary
        assert!(
            n >= NEIGHBOR_PAR_MIN,
            "test must exercise the parallel path"
        );
        let mut sys = System::new(SimBox::cubic(60.0));
        sys.resize(n);
        let mut rng = StdRng::seed_from_u64(11);
        for i in 0..n {
            sys.set_pos(
                i,
                Vec3::new(
                    rng.random_range(0.5..59.5),
                    rng.random_range(0.5..59.5),
                    rng.random_range(0.5..59.5),
                ),
            );
        }
        // Exclusion set: a few thousand tagged pairs (some definitely in-range).
        let excl: Vec<(u32, u32)> = (0..n)
            .step_by(11)
            .flat_map(|i| [(i as u32, (i + 1) as u32), (i as u32, (i + 7) as u32)])
            .collect();
        let mut nl = NeighborList::with_skin_exclusions(2.5, 0.0, &excl);
        nl.force_rebuild(&sys);
        let first = nl.pairs().to_vec();
        // Determinism.
        nl.force_rebuild(&sys);
        nl.force_rebuild(&sys);
        assert_eq!(nl.pairs(), &first[..], "rebuild must be byte-deterministic");
        // Completeness: same set as the naive min-image half-list.
        let lc2 = 2.5 * 2.5;
        let box_size = sys.box_size;
        let mut naive: Vec<(usize, usize)> = Vec::new();
        for i in 0..n {
            let ri = sys.pos(i);
            for j in (i + 1)..n {
                if excl.contains(&(i as u32, j as u32)) || excl.contains(&(j as u32, i as u32)) {
                    continue;
                }
                if box_size.min_image(ri - sys.pos(j)).norm2() < lc2 {
                    naive.push((i, j));
                }
            }
        }
        let mut a = first.clone();
        let mut b = naive;
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "parallel rebuild pair set must equal the brute force");
    }

    #[test]
    fn bounded_stencil_overflow_never_succeeds() {
        let mut sys = sys_with(&[], 32.);
        sys.box_size = SimBox::triclinic([
            Vec3::new(32., 0., 0.),
            Vec3::new(8., 32., 0.),
            Vec3::new(0., 0., 32.),
        ]);
        sys.resize(8);
        for (cluster, base) in [
            [1., 1.5, 2.5],
            [17., 1.5, 2.5],
            [5., 17.5, 18.5],
            [21., 17.5, 18.5],
        ]
        .into_iter()
        .enumerate()
        {
            for dx in 0..2 {
                sys.set_pos(
                    2 * cluster + dx,
                    Vec3::new(base[0] + dx as f64, base[1], base[2]),
                );
            }
        }
        for skin in [0., 0.5] {
            let mut nl = NeighborList::with_skin(2.5, skin);
            // Exercise the actual bounded stencil even in the default build;
            // public diagnostics remain feature gated, bounded API private.
            let mut out = Vec::new();
            let report = build_triclinic::<true>(
                &sys,
                (2.5 + skin) * (2.5 + skin),
                &nl.exclusions,
                4,
                TriclinicBuildStrategy::FractionalStencil,
                &mut out,
            )
            .unwrap();
            #[cfg(feature = "neighbor-diagnostics")]
            assert_eq!(report.selected, TriclinicBuildStrategy::FractionalStencil);
            #[cfg(not(feature = "neighbor-diagnostics"))]
            let _ = report;
            assert_eq!(out, [(0, 1), (2, 3), (4, 5), (6, 7)]);
            out.clear();
            assert!(build_triclinic::<true>(
                &sys,
                (2.5 + skin) * (2.5 + skin),
                &nl.exclusions,
                3,
                TriclinicBuildStrategy::FractionalStencil,
                &mut out
            )
            .is_err());
            assert!(out.len() <= 3 && out.capacity() <= 4);
            #[cfg(feature = "neighbor-diagnostics")]
            {
                nl.triclinic_strategy = TriclinicBuildStrategy::FractionalStencil;
            }
            nl.rebuild_bounded(&sys, 4).unwrap();
            let reference = nl.ref_positions.clone();
            assert!(nl.rebuild_bounded(&sys, 3).is_err());
            assert!(!nl.built);
            assert!(nl.pairs().is_empty());
            assert_eq!(nl.ref_positions, reference);
            assert_eq!(nl.ref_box, Some(sys.box_size));
            #[cfg(feature = "neighbor-diagnostics")]
            assert!(nl.triclinic_build_diagnostics().is_none());
            nl.build(&sys);
            assert_eq!(nl.pairs(), &[(0, 1), (2, 3), (4, 5), (6, 7)]);
        }
    }
}
