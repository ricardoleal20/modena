use super::*;
use modena_core::{
    force_model::{ForceModel, Virial},
    neighbor::{NeighborList, PairScale},
    spme::SpmeForceModel,
    system::{SimBox, System},
};

/// Persistent production model/list. Changing pair policy reconstructs both together;
/// other snapshot changes use existing production mesh/influence/list cache keys.
pub struct SpmeEvaluator {
    pub(super) model: SpmeForceModel,
    pub(super) system: System,
    pub(super) neighbors: NeighborList,
    parameters: Parameters,
    protocol: Protocol,
    scales: Vec<Scale>,
}
impl SpmeEvaluator {
    pub fn new(snapshot: &Snapshot, parameters: Parameters) -> Result<Self, CalibrationFailure> {
        let family = snapshot
            .family
            .split('-')
            .next()
            .ok_or_else(|| invalid("family"))?;
        let replicas = if snapshot.positions.len() == 216 {
            3
        } else {
            2
        };
        let protocol = frozen_protocol(family, replicas, 1, true)?;
        Self::with_protocol(snapshot, parameters, &protocol)
    }
    pub fn with_protocol(
        snapshot: &Snapshot,
        parameters: Parameters,
        p: &Protocol,
    ) -> Result<Self, CalibrationFailure> {
        if p.limits.max_owned_heap_bytes != allocator::PROCESS_HEAP_CAP as u64
            || p.limits.max_threads > 4
            || p.threads > p.limits.max_threads
            || rayon::current_num_threads() != p.threads
        {
            return Err(invalid(
                "SPME evaluator requires registered process heap cap and exact Rayon thread count",
            ));
        }
        let mut checked = snapshot.clone();
        validate_snapshot(&mut checked, p)?;
        validate_parameters(&parameters, checked.cell)?;
        super::candidates::checked_mesh(checked.cell, parameters.kspace_cutoff, &p.limits)?;
        let scales: Vec<_> = checked
            .scales
            .iter()
            .map(|s| PairScale {
                i: s.i,
                j: s.j,
                lj: 1.,
                coulomb: s.coulomb,
            })
            .collect();
        let model = SpmeForceModel::new(
            parameters.kappa,
            parameters.real_cutoff,
            parameters.kspace_cutoff,
        )
        .and_then(|m| m.with_pair_scales(scales.clone()))
        .map_err(|e| invalid(e.to_string()))?;
        let neighbors = NeighborList::with_pair_scales(parameters.real_cutoff, 0.3, &scales)
            .map_err(|e| invalid(e.to_string()))?;
        let mut system = System::new(SimBox::from_array(checked.cell));
        system.resize(checked.charges.len());
        system.positions = checked.positions.iter().flatten().copied().collect();
        system.charges.copy_from_slice(&checked.charges);
        Ok(Self {
            model,
            system,
            neighbors,
            parameters,
            protocol: p.clone(),
            scales: checked.scales,
        })
    }
    pub fn set_snapshot(&mut self, snapshot: &Snapshot) -> Result<(), CalibrationFailure> {
        let mut checked = snapshot.clone();
        validate_snapshot(&mut checked, &self.protocol)?;
        validate_parameters(&self.parameters, checked.cell)?;
        super::candidates::checked_mesh(
            checked.cell,
            self.parameters.kspace_cutoff,
            &self.protocol.limits,
        )?;
        if checked.scales != self.scales || checked.charges.len() != self.system.n() {
            *self = Self::with_protocol(&checked, self.parameters.clone(), &self.protocol)?;
        } else {
            self.system.box_size = SimBox::from_array(checked.cell);
            for (dst, src) in self
                .system
                .positions
                .iter_mut()
                .zip(checked.positions.iter().flatten())
            {
                *dst = *src;
            }
            self.system.charges.copy_from_slice(&checked.charges);
        }
        Ok(())
    }
    /// Raw additive production computation, resetting force/W each call.
    pub fn evaluate(&mut self) -> Result<Evaluation, CalibrationFailure> {
        self.system.zero_forces();
        let mut virial = Virial::default();
        let energy = self
            .model
            .compute(&mut self.system, &mut virial, &mut self.neighbors)
            .map_err(|e| invalid(e.to_string()))?;
        let output = Evaluation {
            energy,
            strain_w: virial.w,
            forces: self
                .system
                .forces
                .chunks_exact(3)
                .map(|f| [f[0], f[1], f[2]])
                .collect(),
        };
        if !output.energy.is_finite()
            || !output.strain_w.is_finite()
            || output.forces.iter().flatten().any(|f| !f.is_finite())
        {
            return Err(CalibrationFailure::ReferencePrecisionLimited(
                "nonfinite production SPME evaluation".into(),
            ));
        }
        Ok(output)
    }
    /// Ordinary persistent production cost, including force/W initialization and
    /// the normal neighbor-cache check. Guard admission/drop and output checking
    /// are outside both clock reads. No profiling, result copies or budget calls.
    #[cfg(test)]
    pub fn post_batch_forces(&self) -> &[f64] {
        &self.system.forces
    }
    pub fn ordinary_batch(&mut self, calls: usize) -> Result<u64, CalibrationFailure> {
        if calls == 0 {
            return Err(invalid("ordinary batch requires nonzero calls"));
        }
        let guard = allocator::NoAllocationGuard::enter().ok_or_else(|| {
            CalibrationFailure::BudgetExceeded(
                "primary timing allocator guard not quiescent".into(),
            )
        })?;
        let mut last = (0., 0.);
        let started = std::time::Instant::now();
        for _ in 0..calls {
            self.system.zero_forces();
            let mut virial = Virial::default();
            let result = self
                .model
                .compute(&mut self.system, &mut virial, &mut self.neighbors);
            let energy = match result {
                Ok(value) => value,
                Err(e) => {
                    drop(guard);
                    return Err(invalid(e.to_string()));
                }
            };
            std::hint::black_box((energy, virial.w, self.system.forces.as_slice()));
            last = (energy, virial.w);
        }
        let elapsed = started.elapsed();
        drop(guard);
        if !last.0.is_finite()
            || !last.1.is_finite()
            || self.system.forces.iter().any(|f| !f.is_finite())
        {
            return Err(CalibrationFailure::ReferencePrecisionLimited(
                "nonfinite timed SPME force".into(),
            ));
        }
        u64::try_from(elapsed.as_nanos())
            .map_err(|_| CalibrationFailure::BudgetExceeded("timing duration overflow".into()))
    }
    /// Separate diagnostic: forced list construction and one normal evaluation.
    /// This intentionally permits setup allocations and never enters primary cost.
    pub fn forced_rebuild_call(&mut self) -> Result<u64, CalibrationFailure> {
        let started = std::time::Instant::now();
        self.neighbors.force_rebuild(&self.system);
        self.compute_without_output_copy()
            .map_err(|e| invalid(e.to_string()))?;
        u64::try_from(started.elapsed().as_nanos())
            .map_err(|_| CalibrationFailure::BudgetExceeded("rebuild duration overflow".into()))
    }
    pub fn actual_mesh(&self) -> [usize; 3] {
        modena_core::spme::mesh_for_box(&self.system.box_size, self.parameters.kspace_cutoff)
    }
    pub fn parameters(&self) -> Parameters {
        self.parameters.clone()
    }
    pub fn neighbor_rebuilds(&self) -> usize {
        self.neighbors.rebuild_count()
    }
    /// Physics/lifecycle diagnostic only: no clock, result allocation or claimed primary timing.
    pub fn compute_without_output_copy(&mut self) -> modena_core::error::Result<f64> {
        self.system.zero_forces();
        let mut virial = Virial::default();
        let energy = self
            .model
            .compute(&mut self.system, &mut virial, &mut self.neighbors)?;
        std::hint::black_box((&self.system.forces, virial.w));
        Ok(energy)
    }
}
pub fn evaluate_spme(snapshot: &Snapshot, p: Parameters) -> Result<Evaluation, CalibrationFailure> {
    SpmeEvaluator::new(snapshot, p)?.evaluate()
}
