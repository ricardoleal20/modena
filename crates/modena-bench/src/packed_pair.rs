//! Scalar packed-layout experiment for reduced, shifted-force LJ.
//!
//! Interleaved delegates to the unmodified production PairForceModel. Packed
//! modes deliberately use serial scalar lane arithmetic, not SIMD or Rayon.
//! Cluster masks follow neighbor epochs, with exact pair/scale comparisons to
//! handle replacement lists whose epoch counters happen to coincide.
use std::{collections::BTreeMap, sync::Mutex};

use modena_core::{
    ForceModel, LennardJones, ModenaError, NeighborList, PairForceModel, PairPotential, Result,
    System, Vec3, Virial,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    Interleaved,
    Soa,
    Aosoa4,
    Cluster4,
}

impl Layout {
    pub fn name(self) -> &'static str {
        match self {
            Self::Interleaved => "interleaved",
            Self::Soa => "soa",
            Self::Aosoa4 => "aosoa4",
            Self::Cluster4 => "cluster4",
        }
    }
}

fn invalid(message: &str) -> ModenaError {
    ModenaError::Input(message.into())
}

fn validate_geometry(s: &System) -> Result<()> {
    if s.box_size.is_triclinic() || s.box_size.size.0.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return Err(invalid(
            "B03 requires a finite positive orthorhombic periodic box",
        ));
    }
    if s.positions.len() != 3 * s.n() || s.positions.iter().any(|x| !x.is_finite()) {
        return Err(invalid("B03 requires exactly 3N finite coordinates"));
    }
    Ok(())
}

/// Conservative experiment domain. The potential is fixed by construction:
/// shifted-force LJ sigma=epsilon=1, cutoff=2.5; arbitrary potentials are not accepted.
pub fn validate_domain(s: &System, nl: &NeighborList) -> Result<()> {
    validate_geometry(s)?;
    if s.forces.len() != 3 * s.n()
        || s.velocities.len() != 3 * s.n()
        || s.species.len() != s.n()
        || s.species.iter().any(|&x| x != 0)
        || s.masses.iter().any(|&m| m != 1.)
        || s.forces.iter().chain(&s.velocities).any(|x| !x.is_finite())
    {
        return Err(invalid(
            "B03 requires finite 3N force/velocity buffers, species=0 and mass=1",
        ));
    }
    if nl.cutoff() != 2.5 || nl.skin() != 0.3 {
        return Err(invalid("B03 requires cutoff=2.5 and skin=0.3"));
    }
    nl.validate_pair_scales(s.n())
}

/// Explicit slot-to-logical-ID bijection; None marks padding, never an atom.
#[derive(Default, Debug)]
pub struct Permutation {
    ids: Vec<Option<usize>>,
    inverse: Vec<usize>,
}

impl Permutation {
    pub fn spatial(s: &System) -> Result<Self> {
        validate_geometry(s)?;
        let mut ids: Vec<_> = (0..s.n()).collect();
        // Cell-major ordering, then coordinates, then logical ID for stable ties.
        ids.sort_unstable_by(|&a, &b| {
            let pa = s.box_size.wrap(s.pos(a));
            let pb = s.box_size.wrap(s.pos(b));
            for k in 0..3 {
                let c = (pa.0[k] / 2.8).floor().total_cmp(&(pb.0[k] / 2.8).floor());
                if !c.is_eq() {
                    return c;
                }
            }
            for k in 0..3 {
                let c = pa.0[k].total_cmp(&pb.0[k]);
                if !c.is_eq() {
                    return c;
                }
            }
            a.cmp(&b)
        });
        Ok(Self::from_ids(ids))
    }
    fn from_ids(ids: Vec<usize>) -> Self {
        let mut inverse = vec![0; ids.len()];
        for (slot, &id) in ids.iter().enumerate() {
            inverse[id] = slot;
        }
        let mut ids: Vec<_> = ids.into_iter().map(Some).collect();
        ids.resize(ids.len().div_ceil(4) * 4, None);
        Self { ids, inverse }
    }
    pub fn logical_ids(&self) -> &[Option<usize>] {
        &self.ids
    }
    pub fn pack<T: Clone>(&self, logical: &[T]) -> Result<Vec<Option<T>>> {
        if logical.len() != self.inverse.len() {
            return Err(invalid("logical buffer length mismatch"));
        }
        Ok(self
            .ids
            .iter()
            .map(|id| id.map(|i| logical[i].clone()))
            .collect())
    }
    pub fn unpack<T: Clone>(&self, packed: &[Option<T>]) -> Result<Vec<T>> {
        if packed.len() != self.ids.len() {
            return Err(invalid("packed buffer length mismatch"));
        }
        self.inverse
            .iter()
            .map(|&slot| {
                packed[slot]
                    .clone()
                    .ok_or_else(|| invalid("missing live lane"))
            })
            .collect()
    }
    fn bytes(&self) -> usize {
        self.ids.capacity() * size_of::<Option<usize>>()
            + self.inverse.capacity() * size_of::<usize>()
    }
}

#[derive(Debug)]
struct Mask {
    blocks: (usize, usize),
    bits: u16,
    scales: [f64; 16],
}

/// Width-four spatial clusters. Only bits derived from exact neighbor pairs are
/// evaluated; zero-LJ/nonzero-Coulomb pairs remain represented with zero scale.
#[derive(Default, Debug)]
pub struct ClusterView {
    permutation: Permutation,
    masks: Vec<Mask>,
}

impl ClusterView {
    /// Caller supplies a neighbor list already built for this system.
    pub fn new(s: &System, nl: &NeighborList) -> Result<Self> {
        nl.validate_pair_scales(s.n())?;
        let permutation = Permutation::spatial(s)?;
        let mut masks: BTreeMap<(usize, usize), Mask> = BTreeMap::new();
        for &(i, j) in nl.pairs() {
            if i >= j || j >= s.n() {
                return Err(invalid("invalid logical half pair"));
            }
            let (mut a, mut b) = (permutation.inverse[i], permutation.inverse[j]);
            if a > b {
                std::mem::swap(&mut a, &mut b);
            }
            let key = (a / 4, b / 4);
            let bit = (a % 4) * 4 + b % 4;
            let mask = masks.entry(key).or_insert(Mask {
                blocks: key,
                bits: 0,
                scales: [0.; 16],
            });
            if mask.bits & (1 << bit) != 0 {
                return Err(invalid("duplicate logical half pair"));
            }
            mask.bits |= 1 << bit;
            mask.scales[bit] = nl.pair_scale(i, j).0;
        }
        Ok(Self {
            permutation,
            masks: masks.into_values().collect(),
        })
    }
    pub fn logical_pairs(&self) -> Vec<(usize, usize, u64)> {
        let mut pairs = Vec::new();
        for mask in &self.masks {
            for bit in 0..16 {
                if mask.bits & (1 << bit) == 0 {
                    continue;
                }
                let i = self.permutation.ids[mask.blocks.0 * 4 + bit / 4].expect("live mask lane");
                let j = self.permutation.ids[mask.blocks.1 * 4 + bit % 4].expect("live mask lane");
                pairs.push((i.min(j), i.max(j), mask.scales[bit].to_bits()));
            }
        }
        pairs
    }
    fn bytes(&self) -> usize {
        self.permutation.bytes() + self.masks.capacity() * size_of::<Mask>()
    }
}

#[derive(Default)]
struct Workspace {
    coordinates: Vec<f64>,
    forces: Vec<f64>,
    cluster: ClusterView,
    epoch: Option<usize>,
    pairs: Vec<(usize, usize, u64)>,
    reconstructions: usize,
    evaluations: usize,
    evaluated_pairs: u64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct LayoutStats {
    pub evaluations: usize,
    pub evaluated_pairs: u64,
    pub cluster_reconstructions: usize,
    /// Owned packed Vec capacities only; excludes core neighbor/model allocations
    /// and transient BTreeMap nodes during cluster reconstruction.
    pub retained_packed_bytes: usize,
}

pub struct PackedForceModel {
    layout: Layout,
    lj: LennardJones,
    baseline: PairForceModel,
    baseline_serial: bool,
    scratch: Mutex<Workspace>,
}

impl PackedForceModel {
    pub fn new(layout: Layout) -> Self {
        let lj = LennardJones::reduced(2.5).expect("fixed valid LJ");
        Self {
            layout,
            baseline: PairForceModel::new(lj.clone()),
            baseline_serial: false,
            lj,
            scratch: Mutex::default(),
        }
    }
    /// Equivalent serial pair-execution control. Neighbor construction still uses
    /// the selected Rayon pool, exactly as it does in all packed candidates.
    pub fn serial_baseline() -> Self {
        let mut model = Self::new(Layout::Interleaved);
        model.baseline_serial = true;
        model
    }
    pub fn stats(&self) -> LayoutStats {
        let w = self.scratch.lock().expect("benchmark workspace poisoned");
        LayoutStats {
            evaluations: w.evaluations,
            evaluated_pairs: w.evaluated_pairs,
            cluster_reconstructions: w.reconstructions,
            retained_packed_bytes: (w.coordinates.capacity() + w.forces.capacity())
                * size_of::<f64>()
                + w.pairs.capacity() * size_of::<(usize, usize, u64)>()
                + w.cluster.bytes(),
        }
    }
    /// Diagnostic expansion of the masks actually used by the most recent
    /// cluster evaluation. Never called in the measured step loop.
    pub fn cached_cluster_pairs(&self) -> Option<Vec<(usize, usize, u64)>> {
        (self.layout == Layout::Cluster4).then(|| {
            self.scratch
                .lock()
                .expect("benchmark workspace poisoned")
                .cluster
                .logical_pairs()
        })
    }
}

#[inline]
fn index(layout: Layout, slot: usize, axis: usize, n: usize) -> usize {
    match layout {
        Layout::Soa => axis * n + slot,
        _ => (slot / 4) * 12 + axis * 4 + slot % 4,
    }
}

impl ForceModel for PackedForceModel {
    fn conserves_total_momentum(&self, _: &System) -> bool {
        true
    }
    fn neighbor_cutoff(&self) -> Option<f64> {
        Some(2.5)
    }
    fn compute(&self, s: &mut System, v: &mut Virial, nl: &mut NeighborList) -> Result<f64> {
        validate_domain(s, nl)?;
        if !v.w.is_finite() {
            return Err(invalid("nonfinite incoming virial"));
        }
        if self.layout == Layout::Interleaved {
            return if self.baseline_serial {
                self.baseline.compute_deterministic(s, v, nl)
            } else {
                self.baseline.compute(s, v, nl)
            };
        }
        nl.build(s);
        let mut guard = self
            .scratch
            .lock()
            .map_err(|_| invalid("benchmark workspace poisoned"))?;
        let w = &mut *guard;
        if self.layout == Layout::Cluster4 {
            let same = w.epoch == Some(nl.rebuild_count())
                && w.cluster.permutation.inverse.len() == s.n()
                && w.pairs.len() == nl.n_pairs()
                && w.pairs
                    .iter()
                    .zip(nl.pairs())
                    .all(|(&(a, b, scale), &(i, j))| {
                        (a, b) == (i, j) && scale == nl.pair_scale(i, j).0.to_bits()
                    });
            if !same {
                w.cluster = ClusterView::new(s, nl)?;
                w.pairs.clear();
                w.pairs.extend(
                    nl.pairs()
                        .iter()
                        .map(|&(i, j)| (i, j, nl.pair_scale(i, j).0.to_bits())),
                );
                w.epoch = Some(nl.rebuild_count());
                w.reconstructions += 1;
            }
        }
        let n = s.n();
        let slots = if self.layout == Layout::Soa {
            n
        } else {
            n.div_ceil(4) * 4
        };
        w.coordinates.resize(slots * 3, 0.);
        w.forces.resize(slots * 3, 0.);
        w.forces.fill(0.);
        for slot in 0..slots {
            let id = if self.layout == Layout::Cluster4 {
                w.cluster.permutation.ids[slot]
            } else {
                (slot < n).then_some(slot)
            };
            for k in 0..3 {
                w.coordinates[index(self.layout, slot, k, n)] =
                    id.map_or(0., |i| s.positions[3 * i + k]);
            }
        }
        let (mut energy, mut virial, mut count) = (0., 0., 0_u64);
        let mut pair = |a: usize, b: usize, scale: f64| {
            let ri = Vec3::new(
                w.coordinates[index(self.layout, a, 0, n)],
                w.coordinates[index(self.layout, a, 1, n)],
                w.coordinates[index(self.layout, a, 2, n)],
            );
            let rj = Vec3::new(
                w.coordinates[index(self.layout, b, 0, n)],
                w.coordinates[index(self.layout, b, 1, n)],
                w.coordinates[index(self.layout, b, 2, n)],
            );
            let d = s.box_size.min_image(ri - rj);
            let r2 = d.norm2();
            if !(1e-14..6.25).contains(&r2) || scale == 0. {
                return;
            }
            let (u, fr) = self.lj.energy_force(r2, 0, 0);
            let f = d * (fr * scale);
            for k in 0..3 {
                w.forces[index(self.layout, a, k, n)] += f.0[k];
                w.forces[index(self.layout, b, k, n)] -= f.0[k];
            }
            energy += u * scale;
            virial += d.dot(f);
            count += 1;
        };
        if self.layout == Layout::Cluster4 {
            for mask in &w.cluster.masks {
                for bit in 0..16 {
                    if mask.bits & (1 << bit) != 0 {
                        pair(
                            mask.blocks.0 * 4 + bit / 4,
                            mask.blocks.1 * 4 + bit % 4,
                            mask.scales[bit],
                        );
                    }
                }
            }
        } else {
            for &(i, j) in nl.pairs() {
                pair(i, j, nl.pair_scale(i, j).0);
            }
        }
        if !energy.is_finite() || !virial.is_finite() || w.forces.iter().any(|x| !x.is_finite()) {
            return Err(invalid("nonfinite packed pair result"));
        }
        for slot in 0..slots {
            let id = if self.layout == Layout::Cluster4 {
                w.cluster.permutation.ids[slot]
            } else {
                (slot < n).then_some(slot)
            };
            if let Some(i) = id {
                for k in 0..3 {
                    s.forces[3 * i + k] += w.forces[index(self.layout, slot, k, n)];
                }
            }
        }
        v.w += virial;
        w.evaluations += 1;
        w.evaluated_pairs += count;
        Ok(energy)
    }
}
