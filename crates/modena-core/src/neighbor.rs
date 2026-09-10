//! Cell-list neighbor search for orthorhombic periodic boxes.
//!
//! Builds a *full* neighbor list (each pair appears in both atoms' lists) by
//! binning atoms into a uniform grid with cell edge ≥ `cutoff`, then scanning
//! the 27-cell stencil per atom. Rebuilt every step in Phase 1 (no skin/Verlet
//! list yet — correctness over cleverness; see roadmap M7).

use crate::system::System;
use crate::vec3::Vec3;

pub struct NeighborList {
    cutoff: f64,
    cutoff2: f64,
    /// `neighbors[i]` = all atom indices within `cutoff` of atom `i` (i≠i).
    neighbors: Vec<Vec<usize>>,
    // Grid state (rebuilt each call).
    dims: [usize; 3],
    cell_size: Vec3,
    head: Vec<i32>, // first atom index in cell, -1 if empty
    next: Vec<i32>, // linked list of atoms within a cell
}

impl NeighborList {
    pub fn new(cutoff: f64) -> Self {
        assert!(cutoff > 0.0, "neighbor cutoff must be > 0");
        Self {
            cutoff,
            cutoff2: cutoff * cutoff,
            neighbors: Vec::new(),
            dims: [1, 1, 1],
            cell_size: Vec3::splat(cutoff),
            head: Vec::new(),
            next: Vec::new(),
        }
    }

    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    pub fn neighbors_of(&self, i: usize) -> &[usize] {
        &self.neighbors[i]
    }

    /// Rebuild the list for the current positions in `system`.
    pub fn build(&mut self, system: &System) {
        let n = system.n();
        self.neighbors.clear();
        self.neighbors.resize(n, Vec::new());
        if n == 0 {
            return;
        }

        let box_size = system.box_size.size;
        // Grid dimensions: at least 1 per axis; cell edge ≥ cutoff. If the box
        // is smaller than the cutoff along an axis we clamp to 1 cell and rely
        // on minimum-image distances (valid for the force, though it means an
        // atom may interact with a periodic image of itself — avoided by the
        // r2 < 1e-14 guard in the force loop for genuine self-pairs).
        let mut dims = [1usize; 3];
        let mut cell_size = [self.cutoff; 3];
        for k in 0..3 {
            let b = box_size.0[k];
            let d = (b / self.cutoff).floor() as usize;
            let d = d.max(1);
            dims[k] = d;
            cell_size[k] = b / d as f64;
        }
        self.dims = dims;
        self.cell_size = Vec3::new(cell_size[0], cell_size[1], cell_size[2]);

        let ncells = dims[0] * dims[1] * dims[2];
        self.head.clear();
        self.head.resize(ncells, -1);
        self.next.clear();
        self.next.resize(n, -1);

        // Bin atoms into cells via a linked list.
        for i in 0..n {
            let c = self.cell_of(system.positions[i]);
            self.next[i] = self.head[c];
            self.head[c] = i as i32;
        }

        // Per-axis neighbor-cell offsets, *deduplicated*. With a coarse grid
        // (fewer than 3 cells per axis) the raw {-1,0,1} stencil wraps multiple
        // offsets onto the same cell; visiting it more than once overcounts
        // pairs by a position-dependent factor and breaks energy conservation.
        //   dim 1 → {0}, dim 2 → {0,1}, dim ≥ 3 → {-1,0,1}.
        let dxs = stencil_deltas(dims[0]);
        let dys = stencil_deltas(dims[1]);
        let dzs = stencil_deltas(dims[2]);

        // For each atom, scan the deduplicated stencil and collect neighbors.
        let cutoff2 = self.cutoff2;
        for i in 0..n {
            let ri = system.positions[i];
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
                            if j != i {
                                let rj = system.positions[j];
                                let d = (ri - rj).min_image(box_size);
                                if d.norm2() < cutoff2 {
                                    self.neighbors[i].push(j);
                                }
                            }
                            p = self.next[j];
                        }
                    }
                }
            }
        }
    }

    #[inline]
    fn cell_of(&self, r: Vec3) -> usize {
        let dims = self.dims;
        let cs = self.cell_size;
        let ix = ((r.0[0] / cs.0[0]).floor() as i32).rem_euclid(dims[0] as i32) as usize;
        let iy = ((r.0[1] / cs.0[1]).floor() as i32).rem_euclid(dims[1] as i32) as usize;
        let iz = ((r.0[2] / cs.0[2]).floor() as i32).rem_euclid(dims[2] as i32) as usize;
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
/// against periodic wraparound. See `NeighborList::build`.
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
    use crate::system::{SimBox, System};

    #[test]
    fn finds_neighbors_in_cubic_box() {
        // Two atoms 1.0 apart in a 10³ box, cutoff 2.5 → mutual neighbors.
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(2);
        sys.positions[0] = Vec3::new(0.0, 0.0, 0.0);
        sys.positions[1] = Vec3::new(1.0, 0.0, 0.0);
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        assert_eq!(nl.neighbors_of(0), &[1]);
        assert_eq!(nl.neighbors_of(1), &[0]);
    }

    #[test]
    fn respects_cutoff() {
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(3);
        sys.positions[0] = Vec3::new(0.0, 0.0, 0.0);
        sys.positions[1] = Vec3::new(1.0, 0.0, 0.0); // within 2.5
        sys.positions[2] = Vec3::new(5.0, 0.0, 0.0); // outside 2.5
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        assert_eq!(nl.neighbors_of(0), &[1]);
        assert_eq!(nl.neighbors_of(2), &[] as &[usize]);
    }

    #[test]
    fn finds_periodic_image_neighbor() {
        // Atom 0 at 0.1, atom 1 at 9.9 in a 10-box → min-image distance 0.2.
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(2);
        sys.positions[0] = Vec3::new(0.1, 0.0, 0.0);
        sys.positions[1] = Vec3::new(9.9, 0.0, 0.0);
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        assert_eq!(nl.neighbors_of(0), &[1]);
    }

    #[test]
    fn no_duplicate_neighbors_with_coarse_grid() {
        // Box edge 5.0, cutoff 2.5 → only 2 cells per axis. The raw {-1,0,1}
        // stencil would wrap -1 and +1 onto the same cell and double-count.
        // Each atom must list each neighbor exactly once.
        let mut sys = System::new(SimBox::cubic(5.0));
        sys.resize(3);
        sys.positions[0] = Vec3::new(0.0, 0.0, 0.0);
        sys.positions[1] = Vec3::new(1.0, 0.0, 0.0);
        sys.positions[2] = Vec3::new(2.0, 0.0, 0.0);
        let mut nl = NeighborList::new(2.5);
        nl.build(&sys);
        // Atom 0 within 2.5 of atom 1 (dist 1) and atom 2 (dist 2).
        let n0 = nl.neighbors_of(0);
        assert_eq!(n0.len(), 2, "atom 0 neighbors: {n0:?}");
        assert!(n0.contains(&1) && n0.contains(&2));
        // No duplicates.
        let mut sorted = n0.to_vec();
        sorted.sort();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            n0.len(),
            "duplicate neighbors detected: {n0:?}"
        );
    }
}
