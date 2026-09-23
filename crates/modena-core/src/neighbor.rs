//! Cell-list neighbor search with a Verlet skin, for orthorhombic periodic
//! boxes (triclinic routes through `SimBox::min_image` — M4 swaps the body).
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

use crate::system::System;
use crate::vec3::Vec3;

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
    /// Box lengths at the last rebuild. A barostat can rescale the box without
    /// any atom moving: the grid dimensions were computed from the old box, and
    /// after a shrink `cell_size < list_cutoff` would silently miss pairs. A
    /// box change therefore forces a rebuild.
    ref_box: [f64; 3],
    built: bool,
    /// Number of grid rebuilds performed (for tests + diagnostics).
    rebuild_count: usize,
    // Grid state (rebuilt as needed).
    dims: [usize; 3],
    cell_size: [f64; 3],
    head: Vec<i32>, // first atom index in cell, -1 if empty
    next: Vec<i32>, // linked list of atoms within a cell
}

impl NeighborList {
    /// New list with force cutoff `cutoff` and Verlet `skin` (skin = 0 disables
    /// the Verlet optimization — rebuild every `build` call, matching Phase 1).
    pub fn with_skin(cutoff: f64, skin: f64) -> Self {
        assert!(cutoff > 0.0, "neighbor cutoff must be > 0");
        assert!(skin >= 0.0, "skin must be >= 0");
        let list_cutoff = cutoff + skin;
        Self {
            cutoff,
            skin,
            list_cutoff,
            list_cutoff2: list_cutoff * list_cutoff,
            pairs: Vec::new(),
            ref_positions: Vec::new(),
            ref_box: [0.0; 3],
            built: false,
            rebuild_count: 0,
            dims: [1, 1, 1],
            cell_size: [cutoff; 3],
            head: Vec::new(),
            next: Vec::new(),
        }
    }

    /// Phase-1-compatible constructor: no skin (rebuild every call).
    pub fn new(cutoff: f64) -> Self {
        Self::with_skin(cutoff, 0.0)
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
            self.pairs.clear();
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
        self.ref_positions.len() == 3 * system.n()
            && self.ref_box
                == [
                    system.box_size.size.0[0],
                    system.box_size.size.0[1],
                    system.box_size.size.0[2],
                ]
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
        let n = system.n();
        let box_size = system.box_size;
        self.rebuild_count += 1;

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

        self.pairs.clear();
        let lc2 = self.list_cutoff2;
        for i in 0..n {
            let ri = system.pos(i);
            let ci = self.cell_of(ri);
            let (ix, iy, iz) = self.unravel(ci);
            for &dz in &dzs {
                for &dy in &dys {
                    for &dx in &dxs {
                        let nx = wrap_cell(ix as i32 + dx, dims[0]);
                        let ny = wrap_cell(iy as i32 + dy, dims[1]);
                        let nz = wrap_cell(iz as i32 + dz, dims[2]);
                        let c = nx * dims[1] * dims[2] + ny * dims[2] + nz;
                        let mut p = self.head[c];
                        while p >= 0 {
                            let j = p as usize;
                            // Half list: only store i < j.
                            if j > i {
                                let rj = system.pos(j);
                                let d = box_size.min_image(ri - rj);
                                if d.norm2() < lc2 {
                                    self.pairs.push((i, j));
                                }
                            }
                            p = self.next[j];
                        }
                    }
                }
            }
        }

        // Snapshot reference positions and box for the skin test.
        self.ref_positions.clear();
        self.ref_positions.extend_from_slice(system.positions_raw());
        self.ref_box = [box_size.size.0[0], box_size.size.0[1], box_size.size.0[2]];
        self.built = true;
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
    fn finds_neighbors_in_cubic_box() {
        let sys = sys_with(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)], 10.0);
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        // Half list: one pair (0,1).
        assert_eq!(nl.pairs(), &[(0, 1)]);
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
}
