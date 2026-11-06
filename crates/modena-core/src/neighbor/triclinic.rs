//! Scalar triclinic reference and the staged strategy-selection seam.

use super::{push_pair, PairScaleSet};
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
    pub(super) selected: TriclinicBuildStrategy,
    pub(super) fallback: Option<FallbackReason>,
    pub(super) candidate_checks: u64,
    pub(super) requested_workspace_bytes: usize,
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
    // Auto deliberately stays scalar until the cost gate. FractionalStencil
    // is only a staged request in Task 5; no stencil or safety test runs yet.
    let selected = match strategy {
        TriclinicBuildStrategy::Auto
        | TriclinicBuildStrategy::Scalar
        | TriclinicBuildStrategy::FractionalStencil => TriclinicBuildStrategy::Scalar,
    };
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
    Ok(BuildReport {
        selected,
        fallback: None,
        candidate_checks,
        requested_workspace_bytes: 0,
    })
}
