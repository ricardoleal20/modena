//! Certified fractional-bin candidates with the unchanged scalar predicate.

use super::{push_pair, PairScaleSet};
use crate::periodic_geometry::{
    certified_cell, fractional_bin, fractional_position, scalar_min_image_search_radius,
    CertifiedCell, GeometryError, Interval,
};
use crate::System;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TriclinicBuildStrategy {
    Auto,
    Scalar,
    FractionalStencil,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum FallbackReason {
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

pub(super) struct BuildReport {
    #[cfg(feature = "neighbor-diagnostics")]
    pub(super) selected: TriclinicBuildStrategy,
    #[cfg(feature = "neighbor-diagnostics")]
    pub(super) fallback: Option<FallbackReason>,
    #[cfg(feature = "neighbor-diagnostics")]
    pub(super) candidate_checks: u64,
    #[cfg(feature = "neighbor-diagnostics")]
    pub(super) requested_workspace_bytes: usize,
}

impl BuildReport {
    fn new(
        _selected: TriclinicBuildStrategy,
        _fallback: Option<FallbackReason>,
        _candidate_checks: u64,
        _requested_workspace_bytes: usize,
    ) -> Self {
        Self {
            #[cfg(feature = "neighbor-diagnostics")]
            selected: _selected,
            #[cfg(feature = "neighbor-diagnostics")]
            fallback: _fallback,
            #[cfg(feature = "neighbor-diagnostics")]
            candidate_checks: _candidate_checks,
            #[cfg(feature = "neighbor-diagnostics")]
            requested_workspace_bytes: _requested_workspace_bytes,
        }
    }
}

const GRID_CAP: usize = 1_000_000;
const STENCIL_CAP: usize = 4096;
const WORK_CAP: usize = 50_000_000;
const MEMORY_CAP: usize = 134_217_728;
type BinAtom = (usize, usize);

struct Prepared {
    dims: [usize; 3],
    offsets: Vec<[usize; 3]>,
    atoms: Vec<BinAtom>,
}

fn geometry(error: GeometryError) -> FallbackReason {
    match error {
        GeometryError::ShiftRange => FallbackReason::ImageDomain,
        _ => FallbackReason::UncertifiedGeometry,
    }
}

fn fallback_for_cost(candidate_estimate: u128, n: usize) -> bool {
    let n128 = n as u128;
    let all_pairs = match n128.checked_mul(n128.saturating_sub(1)) {
        Some(value) => value / 2,
        None => return true,
    };
    match candidate_estimate
        .checked_mul(5)
        .zip(all_pairs.checked_mul(4))
    {
        Some((candidate_scaled, scalar_scaled)) => candidate_scaled >= scalar_scaled,
        None => true,
    }
}

// These checks operate on counts, before any corresponding heap request. The
// sorted (bin, atom) list replaces a million-bin dense table and never inserts
// an atom twice. The output pairs have their separate existing bounded cap.
fn check_workspace(bytes: usize) -> Result<(), FallbackReason> {
    if bytes > MEMORY_CAP {
        Err(FallbackReason::MemoryCap)
    } else {
        Ok(())
    }
}

fn preflight(n: usize, dims: [usize; 3], offsets: usize) -> Result<usize, FallbackReason> {
    let grid = dims
        .into_iter()
        .try_fold(1usize, |p, d| p.checked_mul(d))
        .filter(|&p| p > 0 && p <= GRID_CAP)
        .ok_or(FallbackReason::GridCap)?;
    if offsets > STENCIL_CAP {
        return Err(FallbackReason::StencilCap);
    }
    if n > i32::MAX as usize || n.checked_mul(offsets).is_none_or(|v| v > WORK_CAP) {
        return Err(FallbackReason::WorkCap);
    }
    if offsets == grid {
        return Err(FallbackReason::AllBins);
    }
    workspace_bytes(n, offsets)
}

fn workspace_bytes(n: usize, offsets: usize) -> Result<usize, FallbackReason> {
    let bytes = n
        .checked_mul(std::mem::size_of::<BinAtom>())
        .and_then(|b| {
            offsets
                .checked_mul(std::mem::size_of::<[usize; 3]>())
                .and_then(|s| b.checked_add(s))
        })
        .and_then(|b| {
            b.checked_add(
                std::mem::size_of::<Prepared>()
                    + std::mem::size_of::<CertifiedCell>()
                    + std::mem::size_of::<[AxisOffsets; 3]>(),
            )
        })
        .ok_or(FallbackReason::MemoryCap)?;
    Ok(bytes)
}

#[derive(Clone, Copy)]
struct AxisOffsets {
    values: [usize; 7],
    len: usize,
}

fn axis_offsets(m: usize, w: f64) -> Result<AxisOffsets, FallbackReason> {
    let mut result = AxisOffsets {
        values: [0; 7],
        len: 0,
    };
    let p = Interval::point;
    let width = p(m as f64).and_then(|m| m.mul(p(w)?)).map_err(geometry)?;
    // m=floor(1/w) (rounded conservatively down) implies m*w<=1
    // whenever m>1; outward multiplication can add one endpoint ulp,
    // so ceil(width)+1<=3 and at most seven pre-wrap offsets suffice.
    // For w>=1, m=max(1,floor(1/w))=1 and every integer offset is
    // residue zero. Deduplicate that case BEFORE enumerating any range.
    let end = if m == 1 {
        0
    } else {
        let e = width.hi.ceil() + 1.0;
        if e > 3.0 {
            return Err(FallbackReason::UncertifiedGeometry);
        }
        e as i64
    };
    let threshold = p(w)
        .and_then(|w| w.add(p(1.)?.div(p(m as f64)?)?))
        .map_err(geometry)?;
    for q in -end..=end {
        let wrapped = q.rem_euclid(m as i64) as usize;
        if result.values[..result.len].contains(&wrapped) {
            continue;
        }
        // dist(q/m,Z)=min(residue,m-residue)/m. Reject only when the
        // LOWER enclosure exceeds the UPPER necessary bound w+1/m.
        // Uncertain equality is included. All integer operands are <=1e6.
        let distance = p(wrapped.min(m - wrapped) as f64)
            .and_then(|v| v.div(p(m as f64)?))
            .map_err(geometry)?;
        if distance.lo > threshold.hi {
            continue;
        }
        if result.len == result.values.len() {
            return Err(FallbackReason::StencilCap);
        }
        result.values[result.len] = wrapped;
        result.len += 1;
    }
    Ok(result)
}

impl Prepared {
    // Each residue occurs once; translation modulo dims is a bijection.
    // source<=target therefore visits each unordered occupied-bin pair once,
    // independent of aliases for dimensions 1, 2 or 3. No global pair dedup.
    fn visit<E>(
        &self,
        mut visit: impl FnMut(&[BinAtom], &[BinAtom]) -> Result<(), E>,
    ) -> Result<(), E> {
        for source in self.atoms.chunk_by(|a, b| a.0 == b.0) {
            let bin = source[0].0;
            let origin = [
                bin % self.dims[0],
                (bin / self.dims[0]) % self.dims[1],
                bin / (self.dims[0] * self.dims[1]),
            ];
            for offset in &self.offsets {
                let t = std::array::from_fn::<_, 3, _>(|k| (origin[k] + offset[k]) % self.dims[k]);
                let target = t[0] + self.dims[0] * (t[1] + self.dims[1] * t[2]);
                if target < bin {
                    continue;
                }
                let first = self.atoms.partition_point(|a| a.0 < target);
                let last = self.atoms.partition_point(|a| a.0 <= target);
                if first != last {
                    visit(source, &self.atoms[first..last])?;
                }
            }
        }
        Ok(())
    }

    fn candidate_estimate(&self) -> Result<u128, FallbackReason> {
        let mut total = 0u128;
        self.visit(|a, b| {
            let na = a.len() as u128;
            let count = if a[0].0 == b[0].0 {
                na.checked_mul(na.saturating_sub(1)).map(|x| x / 2)
            } else {
                na.checked_mul(b.len() as u128)
            };
            total = count
                .and_then(|x| total.checked_add(x))
                .ok_or(FallbackReason::CandidateCost)?;
            Ok(())
        })?;
        Ok(total)
    }
}

fn prepare(
    system: &System,
    cutoff2: f64,
    requested: &mut usize,
) -> Result<Prepared, FallbackReason> {
    let cert = certified_cell(system.box_size).map_err(geometry)?;
    // The caller stores R2. An upward sqrt supplies a radius whose square
    // covers that exact stored threshold; Q is search-only and final pruning
    // still uses the original R2, with no reassociation or epsilon shell.
    let q = scalar_min_image_search_radius(&cert, &system.positions, cutoff2.sqrt().next_up())
        .map_err(geometry)?;
    let mut dims = [1; 3];
    let mut axes = [AxisOffsets {
        values: [0; 7],
        len: 0,
    }; 3];
    for k in 0..3 {
        let w = Interval::point(q)
            .and_then(|q| q.mul(Interval::point(cert.inverse_row_norm_upper[k])?))
            .map_err(geometry)?
            .hi;
        let reciprocal = Interval::point(1.)
            .and_then(|one| one.div(Interval::point(w)?))
            .map_err(geometry)?;
        // Taking the LOWER endpoint cannot overstate subdivisions at a
        // reciprocal integer boundary. Any coarser positive grid preserves
        // the bin-offset proof; reject an over-cap grid before integer casts.
        let m = reciprocal.lo.floor().max(1.);
        if m > GRID_CAP as f64 {
            return Err(FallbackReason::GridCap);
        }
        dims[k] = m as usize;
        axes[k] = axis_offsets(dims[k], w)?;
    }
    let count = axes
        .iter()
        .try_fold(1usize, |p, a| p.checked_mul(a.len))
        .ok_or(FallbackReason::StencilCap)?;
    *requested = preflight(system.n(), dims, count)?;
    check_workspace(*requested)?;
    let mut plan = Prepared {
        dims,
        offsets: Vec::new(),
        atoms: Vec::new(),
    };
    plan.offsets
        .try_reserve_exact(count)
        .map_err(|_| FallbackReason::MemoryCap)?;
    plan.atoms
        .try_reserve_exact(system.n())
        .map_err(|_| FallbackReason::MemoryCap)?;
    let actual = workspace_bytes(plan.atoms.capacity(), plan.offsets.capacity())?;
    *requested = (*requested).max(actual);
    check_workspace(*requested)?;
    for &x in &axes[0].values[..axes[0].len] {
        for &y in &axes[1].values[..axes[1].len] {
            for &z in &axes[2].values[..axes[2].len] {
                plan.offsets.push([x, y, z]);
            }
        }
    }
    // For exact f_i=U*r_i and an accepted physical image, each periodic
    // difference is <=w. Certified bins give f_i modulo Z in [b_i/m,(b_i+1)/m).
    // Thus dist((b_i-b_j)/m,Z)<=w+1/m and |q|<=ceil(m*w)+1 suffices.
    // scalar_min_image_search_radius includes the scalar floating discrepancy;
    // fractional_position covers the true inverse; fractional_bin certifies a
    // single bin. Together these establish real AND finite-precision coverage.
    for i in 0..system.n() {
        let fractions = fractional_position(&cert, system.pos(i)).map_err(geometry)?;
        let mut b = [0; 3];
        for k in 0..3 {
            b[k] = fractional_bin(fractions[k], dims[k]).map_err(|e| match e {
                GeometryError::ShiftRange => FallbackReason::ImageDomain,
                _ => FallbackReason::FractionalBoundary,
            })?;
        }
        plan.atoms
            .push((b[0] + dims[0] * (b[1] + dims[1] * b[2]), i));
    }
    plan.atoms.sort_unstable();
    if fallback_for_cost(plan.candidate_estimate()?, system.n()) {
        return Err(FallbackReason::CandidateCost);
    }
    Ok(plan)
}

/// Append the original lexicographic i<j enumeration to the caller's buffer.
/// Output allocation/cap handling stays shared with the orthorhombic builder.
pub(super) fn build_triclinic<const BOUNDED: bool>(
    system: &System,
    list_cutoff2: f64,
    exclusions: &PairScaleSet,
    max_pairs: usize,
    strategy: TriclinicBuildStrategy,
    out: &mut Vec<(usize, usize)>,
) -> crate::Result<BuildReport> {
    // Auto stays scalar until G4. Preparation cannot append pairs: a safety
    // fallback starts the unchanged scalar scan with no partial stencil data.
    let mut requested = 0;
    let mut fallback = None;
    if strategy == TriclinicBuildStrategy::FractionalStencil && system.box_size.is_triclinic() {
        match prepare(system, list_cutoff2, &mut requested) {
            Ok(plan) => {
                let mut checks = 0u64;
                plan.visit::<crate::ModenaError>(|a, b| {
                    for &(_, ai) in a {
                        for &(_, bj) in b {
                            if a[0].0 == b[0].0 && ai >= bj {
                                continue;
                            }
                            let (i, j) = (ai.min(bj), ai.max(bj));
                            if exclusions.contains(i, j) {
                                continue;
                            }
                            checks += 1; // <=WORK_CAP*N within the admitted i32 domain.
                            let d = system.box_size.min_image(system.pos(i) - system.pos(j));
                            if d.norm2() < list_cutoff2 {
                                // A bounded failure propagates; it NEVER becomes a scalar
                                // retry or a successful partial list.
                                push_pair::<BOUNDED>(out, (i, j), max_pairs)?;
                            }
                        }
                    }
                    Ok(())
                })?;
                out.sort_unstable();
                return Ok(BuildReport::new(
                    TriclinicBuildStrategy::FractionalStencil,
                    None,
                    checks,
                    requested,
                ));
            }
            Err(reason) => fallback = Some(reason),
        }
    }
    let n = system.n();
    if !BOUNDED {
        out.reserve(n * 27);
    }
    let mut candidate_checks = 0u64;
    for i in 0..n {
        let ri = system.pos(i);
        for j in (i + 1)..n {
            if exclusions.contains(i, j) {
                continue;
            }
            candidate_checks = candidate_checks.saturating_add(1);
            let d = system.box_size.min_image(ri - system.pos(j));
            if d.norm2() < list_cutoff2 {
                push_pair::<BOUNDED>(out, (i, j), max_pairs)?;
            }
        }
    }
    Ok(BuildReport::new(
        TriclinicBuildStrategy::Scalar,
        fallback,
        candidate_checks,
        requested,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrapped_offsets_are_unique_for_dims_1_2_3() {
        for m in 1..=3 {
            let axis = axis_offsets(m, 0.9 / m as f64).unwrap();
            let mut values = axis.values[..axis.len].to_vec();
            values.sort_unstable();
            assert_eq!(values, (0..m).collect::<Vec<_>>());
            for source in 0..m {
                let mut targets: Vec<_> = values.iter().map(|q| (source + q) % m).collect();
                targets.sort_unstable();
                assert_eq!(targets, (0..m).collect::<Vec<_>>());
            }
        }
        // Mixed 1/2/3 dimensions exercise real traversal and its estimate.
        let axes = [1, 2, 3].map(|m| axis_offsets(m, 0.9 / m as f64).unwrap());
        let mut plan = Prepared {
            dims: [1, 2, 3],
            offsets: Vec::new(),
            atoms: (0..6).map(|b| (b, b)).collect(),
        };
        for &x in &axes[0].values[..axes[0].len] {
            for &y in &axes[1].values[..axes[1].len] {
                for &z in &axes[2].values[..axes[2].len] {
                    plan.offsets.push([x, y, z]);
                }
            }
        }
        let mut visits = Vec::new();
        plan.visit::<()>(|a, b| {
            visits.push((a[0].0, b[0].0));
            Ok(())
        })
        .unwrap();
        visits.sort_unstable();
        let expected: Vec<_> = (0..6).flat_map(|i| (i..6).map(move |j| (i, j))).collect();
        assert_eq!(visits, expected);
        assert_eq!(plan.candidate_estimate().unwrap(), 15);
    }

    #[test]
    fn resource_caps_precede_allocation_and_cost_uses_checked_arithmetic() {
        assert_eq!(preflight(8, [101; 3], 27), Err(FallbackReason::GridCap));
        assert_eq!(
            preflight(8, [usize::MAX; 3], 27),
            Err(FallbackReason::GridCap)
        );
        // Normal m=floor(1/w) admits <=7^3 offsets, so the 4096 limit is
        // defensive. Inject COUNTS here, not huge atom vectors or altered caps.
        assert_eq!(
            preflight(8, [100; 3], 4097),
            Err(FallbackReason::StencilCap)
        );
        assert_eq!(
            preflight(2_000_000, [100; 3], 27),
            Err(FallbackReason::WorkCap)
        );
        assert_eq!(
            preflight(usize::MAX, [100; 3], 27),
            Err(FallbackReason::WorkCap)
        );
        assert_eq!(preflight(8, [3; 3], 27), Err(FallbackReason::AllBins));
        assert_eq!(
            check_workspace(preflight(9_000_000, [100; 3], 1).unwrap()),
            Err(FallbackReason::MemoryCap)
        );
        assert!(!fallback_for_cost(7, 5)); // 7 < .8*10
        assert!(fallback_for_cost(8, 5));
        assert!(fallback_for_cost(u128::MAX, 5));
        assert!(fallback_for_cost(0, 0));
        assert!(fallback_for_cost(0, 1));
        assert!(fallback_for_cost(0, usize::MAX)); // 4*all_pairs overflows u128 on 64-bit.
    }

    #[test]
    fn occupancy_cost_counts_same_and_cross_bin_pairs_once() {
        let plan = Prepared {
            dims: [10, 1, 1],
            offsets: vec![[0, 0, 0], [1, 0, 0], [9, 0, 0]],
            atoms: vec![(0, 0), (0, 1), (1, 2), (1, 3), (1, 4), (8, 5)],
        };
        // C(2,2)+C(3,2)+2*3=10; isolated atom contributes zero.
        assert_eq!(plan.candidate_estimate().unwrap(), 10);
    }
}
