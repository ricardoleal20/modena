//! Composite force model — runs several `ForceModel`s in sequence and sums
//! their potential-energy and virial contributions.
//!
//! This is the M3 glue that lets a pair (non-bonded) model coexist with a
//! [`crate::bonded::BondedForceModel`] (and, later, a Coulomb/Ewald model). It
//! owns the **zero-then-accumulate** contract: `compute` zeroes `system.forces`
//! and `virial` once, then each sub-model adds into them.

use crate::bonded::BondedForceModel;
use crate::error::Result;
use crate::ewald::EwaldForceModel;
use crate::force_model::{ForceModel, PairForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::spme::SpmeForceModel;
use crate::system::System;

/// The long-range Coulomb sub-model: direct reciprocal-sum Ewald or SPME
/// (mesh-based, O(N log N)). Selected via `[electrostatics] method` in the
/// input; both share the same real-space cutoff contract.
pub enum LongRangeForceModel {
    Ewald(EwaldForceModel),
    Spme(SpmeForceModel),
}

impl LongRangeForceModel {
    /// The real-space erfc cutoff — the shared-neighbor-list maximum.
    pub fn real_cutoff(&self) -> f64 {
        match self {
            Self::Ewald(m) => m.real_cutoff(),
            Self::Spme(m) => m.real_cutoff(),
        }
    }
}

impl ForceModel for LongRangeForceModel {
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        match self {
            Self::Ewald(m) => m.compute(system, virial, neighbors),
            Self::Spme(m) => m.compute(system, virial, neighbors),
        }
    }
}

/// A pair model plus optional bonded and long-range Coulomb sub-models.
/// Additional sub-models (many-body, ML surrogate) slot in as further
/// `Option<Box<dyn ForceModel>>` fields without changing the contract.
pub struct CompositeForceModel {
    pub pair: PairForceModel,
    pub bonded: Option<BondedForceModel>,
    pub long_range: Option<LongRangeForceModel>,
}

impl CompositeForceModel {
    pub fn new(pair: PairForceModel) -> Self {
        Self {
            pair,
            bonded: None,
            long_range: None,
        }
    }

    pub fn with_bonded(mut self, bonded: BondedForceModel) -> Self {
        self.bonded = Some(bonded);
        self
    }

    /// Add direct reciprocal-sum Ewald electrostatics to the composition.
    pub fn with_ewald(mut self, ewald: EwaldForceModel) -> Self {
        self.long_range = Some(LongRangeForceModel::Ewald(ewald));
        self
    }

    /// Add SPME (mesh) electrostatics to the composition.
    pub fn with_spme(mut self, spme: SpmeForceModel) -> Self {
        self.long_range = Some(LongRangeForceModel::Spme(spme));
        self
    }
}

impl ForceModel for CompositeForceModel {
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        // Zero once; sub-models accumulate.
        system.forces.fill(0.0);
        virial.w = 0.0;

        let mut pe = self.pair.compute(system, virial, neighbors)?;
        if let Some(b) = &self.bonded {
            pe += b.compute(system, virial, neighbors)?;
        }
        if let Some(e) = &self.long_range {
            pe += e.compute(system, virial, neighbors)?;
        }
        Ok(pe)
    }
}
