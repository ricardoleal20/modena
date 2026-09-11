//! Composite force model — runs several `ForceModel`s in sequence and sums
//! their potential-energy and virial contributions.
//!
//! This is the M3 glue that lets a pair (non-bonded) model coexist with a
//! [`crate::bonded::BondedForceModel`] (and, later, a Coulomb/Ewald model). It
//! owns the **zero-then-accumulate** contract: `compute` zeroes `system.forces`
//! and `virial` once, then each sub-model adds into them.

use crate::bonded::BondedForceModel;
use crate::error::Result;
use crate::force_model::{ForceModel, PairForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::system::System;

/// A pair model plus an optional bonded model. (The struct is intentionally
/// small now; additional sub-models — Coulomb, many-body — slot in as further
/// `Option<Box<dyn ForceModel>>` fields without changing the contract.)
pub struct CompositeForceModel {
    pub pair: PairForceModel,
    pub bonded: Option<BondedForceModel>,
}

impl CompositeForceModel {
    pub fn new(pair: PairForceModel) -> Self {
        Self { pair, bonded: None }
    }

    pub fn with_bonded(mut self, bonded: BondedForceModel) -> Self {
        self.bonded = Some(bonded);
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
        Ok(pe)
    }
}
