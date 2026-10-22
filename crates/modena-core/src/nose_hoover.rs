//! Nosé–Hoover chain (NHC) thermostat — a deterministic integrator that samples
//! the exact canonical (NVT) ensemble.
//!
//! Unlike the velocity-only [`Thermostat`](crate::thermostat::Thermostat)
//! implementations (Berendsen, velocity-rescaling), the NHC is an
//! [`Integrator`](crate::integrator::Integrator): it interleaves the thermostat
//! chain propagation with the velocity-Verlet kick/drift, because the chain
//! variables couple to (and rescale) the particle velocities mid-step.
//!
//! The chain (length `M`, default 4) is integrated with the 3-fold
//! Suzuki–Yoshida decomposition (Martyna, Tuckerman, Tobias, Klein, 1992). The
//! extended Hamiltonian
//!
//! ```text
//!   H* = KE + U + Σ_k p_ξk²/(2 Q_k) + N_f kT ξ_1 + Σ_{k≥2} kT ξ_k
//! ```
//!
//! is conserved by the integrator (a shadow Hamiltonian) — the conservation
//! test in `tests/ensembles.rs` is the correctness check. Chain masses couple to
//! the time constant `tau`: `Q_1 = N_f kT τ²`, `Q_k = kT τ²` for k ≥ 2.

use crate::constraints::Constraints;
use crate::error::Result;
use crate::force_model::{ForceModel, Virial};
use crate::integrator::{prepare_constraints, project_drift, Integrator};
use crate::neighbor::NeighborList;
use crate::system::System;
use crate::thermo::Energies;

/// Suzuki–Yoshida (3-fold) weights. Sum = 1; corrects the integrator to 4th
/// order in the chain propagation.
const SY3: [f64; 3] = [1.3512071919596578, -1.7024143839193153, 1.3512071919596578];

const SHAKE_TOL: f64 = 1e-10;
const SHAKE_MAX_ITER: usize = 100;

/// Nosé–Hoover chain NVT integrator (velocity-Verlet + NHC).
pub struct NoseHooverChain {
    pub target_temperature: f64,
    pub tau: f64,
    chain_len: usize,
    /// Chain momenta p_ξk.
    p_xi: Vec<f64>,
    /// Chain positions ξk.
    xi: Vec<f64>,
    pub constraints: Option<Constraints>,
}

impl NoseHooverChain {
    pub fn new(target_temperature: f64, tau: f64, chain_len: usize) -> Self {
        let chain_len = chain_len.max(1);
        Self {
            target_temperature,
            tau,
            chain_len,
            p_xi: vec![0.0; chain_len],
            xi: vec![0.0; chain_len],
            constraints: None,
        }
    }

    pub fn with_constraints(mut self, constraints: Constraints) -> Self {
        self.constraints = Some(constraints);
        self
    }

    /// Piston masses: Q_1 = N_f kT τ², Q_k = kT τ² (k ≥ 2).
    fn masses(&self, dof: f64) -> Vec<f64> {
        let kt = self.target_temperature;
        let mut q = vec![0.0; self.chain_len];
        q[0] = (dof * kt).max(1e-12) * self.tau * self.tau;
        for q_k in q.iter_mut().take(self.chain_len).skip(1) {
            *q_k = kt.max(1e-12) * self.tau * self.tau;
        }
        q
    }

    /// The conserved extended quantity H* (for testing / logging).
    pub fn extended_hamiltonian(&self, system: &System, pe: f64) -> f64 {
        let dof = system.dof() as f64;
        let kt = self.target_temperature;
        let q = self.masses(dof);
        let ke = system.kinetic_energy();
        let mut h = ke + pe;
        for (p_xi_k, q_k) in self.p_xi.iter().zip(q.iter()).take(self.chain_len) {
            h += p_xi_k * p_xi_k / (2.0 * q_k);
        }
        // Potential of the chain: N_f kT ξ_1 + Σ_{k≥2} kT ξ_k.
        h += dof * kt * self.xi[0];
        for k in 1..self.chain_len {
            h += kt * self.xi[k];
        }
        h
    }
}

impl Integrator for NoseHooverChain {
    fn checkpoint_state(&self) -> Result<serde_json::Value> {
        if let Some(cons) = &self.constraints {
            cons.validate_parameters()?;
        }
        if !self.target_temperature.is_finite()
            || !self.tau.is_finite()
            || self.p_xi.iter().chain(&self.xi).any(|x| !x.is_finite())
        {
            return Err(crate::thermostat::state_error("nonfinite NHC state"));
        }
        Ok(serde_json::json!({"kind":"nose_hoover_chain","version":1,
            "target_temperature":self.target_temperature,"tau":self.tau,"chain_len":self.chain_len,
            "constraints":self.constraints,"p_xi":self.p_xi,"xi":self.xi}))
    }
    fn restore_state(&mut self, state: &serde_json::Value) -> Result<()> {
        crate::thermostat::validate_state_parameters(
            state,
            &self.checkpoint_state()?,
            &[
                "kind",
                "version",
                "target_temperature",
                "tau",
                "chain_len",
                "constraints",
            ],
        )?;
        let parse = |key: &str| -> Result<Vec<f64>> {
            let v: Vec<f64> = serde_json::from_value(state[key].clone())
                .map_err(|err| crate::thermostat::state_error(&err.to_string()))?;
            if v.len() != self.chain_len || v.iter().any(|x| !x.is_finite()) {
                return Err(crate::thermostat::state_error(
                    "invalid NHC extended-state array",
                ));
            }
            Ok(v)
        };
        let p_xi = parse("p_xi")?;
        let xi = parse("xi")?;
        self.p_xi = p_xi;
        self.xi = xi;
        Ok(())
    }
    fn step(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        dt: f64,
    ) -> Result<Energies> {
        prepare_constraints(system, self.constraints.as_ref())?;
        let n = system.n();
        let dof = system.dof() as f64;
        let kt = self.target_temperature;
        let q = self.masses(dof);

        // --- NHC half-step (advances the chain by dt/2, rescales v) ---
        propagate_chain(system, &mut self.p_xi, &mut self.xi, &q, kt, dof, 0.5 * dt);

        // First half-kick: v += dt/2 · F/m  (using cached F(t)).
        kick(system, n, 0.5 * dt);

        project_drift(system, self.constraints.as_ref(), dt, false)?;

        // New forces F(t+dt).
        system.forces.fill(0.0);
        let mut virial = Virial::default();
        let pe = model.compute(system, &mut virial, neighbors)?;

        // Second half-kick: v += dt/2 · F/m.
        kick(system, n, 0.5 * dt);

        // RATTLE.
        if let Some(cons) = &self.constraints {
            cons.rattle(system, SHAKE_TOL, SHAKE_MAX_ITER)?;
        }

        // --- NHC half-step again (completes the full dt) ---
        propagate_chain(system, &mut self.p_xi, &mut self.xi, &q, kt, dof, 0.5 * dt);

        Ok(Energies::from(system, pe, &virial))
    }
}

/// One NHC propagation of length `dt`: Suzuki–Yoshida (3 substeps), advancing
/// the chain momenta `p_xi`, positions `xi`, and rescaling the particle
/// velocities.
fn propagate_chain(
    system: &mut System,
    p_xi: &mut [f64],
    xi: &mut [f64],
    q: &[f64],
    kt: f64,
    dof: f64,
    dt: f64,
) {
    let m = p_xi.len();
    for &w in &SY3 {
        // Suzuki–Yoshida 3rd order: the composition S₂(wᵢ·dt) advances the
        // chain by Σ wᵢ·dt = dt in total. The weights already sum to 1, so an
        // extra `/ SY3.len()` divided the whole propagate by 3 — every NHC run
        // silently used τ_eff = 3·τ (the user-facing response time was 3× too
        // slow). The extended-Hamiltonian conservation test cannot catch this
        // (its bookkeeping balances per substep regardless of the time scale).
        let ddt = w * dt;
        if ddt == 0.0 {
            continue;
        }
        let half = 0.5 * ddt;

        // Chain *forces* G_k (energy units — the force ON p_ξk, so the kick is
        // p_ξk += half·G_k, with NO extra /Q_k). `p_xi` stores the chain momentum
        // p_ξ (consistent with `extended_hamiltonian`'s p_ξ²/(2Q) term and the
        // dξ/dt = p_ξ/Q position update). Only the p²/Q inside G_k divides by a
        // mass (it is the kinetic energy of the previous link).
        //
        //   G_1 = 2·KE − N_f·kT
        //   G_k = p_{k-1}²/Q_{k-1} − kT        (k ≥ 2)
        let force = |k: usize, p_xi: &[f64], ke: f64| -> f64 {
            if k == 0 {
                2.0 * ke - dof * kt
            } else {
                p_xi[k - 1] * p_xi[k - 1] / q[k - 1] - kt
            }
        };

        // Backward half-kick of the chain momenta (M → 1). Each G_k uses the
        // not-yet-kicked p_{k-1} (lower index), which is correct for the backward
        // sweep.
        let ke = system.kinetic_energy();
        for k in (0..m).rev() {
            let gk = force(k, p_xi, ke);
            p_xi[k] += half * gk;
        }

        // Half-scale (the friction operator iL_scale). The NHC is *hierarchical*:
        // link k's momentum is damped by link k+1 exactly as the particle
        // velocities are damped by link 1. So the scaling is applied end-to-start
        // (p_ξ{M-1} by p_ξM, …, p_ξ1 by p_ξ2), then the particle velocities by
        // the just-scaled p_ξ1. p_ξM (the chain end) is unscaled. Omitting the
        // chain-link scaling makes any chain of length ≥ 2 run away — the
        // downstream links decouple and pump energy into link 1.
        scale_friction(system, p_xi, q, half);

        // Drift the chain positions: dξ_k/dt = p_ξk / Q_k.
        for k in 0..m {
            xi[k] += ddt * p_xi[k] / q[k];
        }

        // Second half-scale.
        scale_friction(system, p_xi, q, half);

        // Forward half-kick (1 → M). Each G_k uses the just-kicked p_{k-1} (the
        // forward sweep visits lower indices first); G_1 is recomputed with the
        // post-scaling KE.
        let ke = system.kinetic_energy();
        for k in 0..m {
            let gk = if k == 0 {
                2.0 * ke - dof * kt
            } else {
                p_xi[k - 1] * p_xi[k - 1] / q[k - 1] - kt
            };
            p_xi[k] += half * gk;
        }
    }
}

/// The NHC friction half-step `iL_scale(half)`: scale each chain link's momentum
/// by the next link (end → start), then the particle velocities by link 1.
fn scale_friction(system: &mut System, p_xi: &mut [f64], q: &[f64], half: f64) {
    let m = p_xi.len();
    // Chain links: p_ξk *= exp(−half · p_ξ{k+1}/Q_{k+1}), k = M-2 … 0.
    for k in (0..m.saturating_sub(1)).rev() {
        let f = (-half * p_xi[k + 1] / q[k + 1]).exp();
        p_xi[k] *= f;
    }
    // Particle velocities: v *= exp(−half · p_ξ1/Q_1), using the just-scaled p_ξ1.
    let fv = (-half * p_xi[0] / q[0]).exp();
    for v in system.velocities.iter_mut() {
        *v *= fv;
    }
}

#[inline]
fn kick(system: &mut System, n: usize, factor: f64) {
    let v = &mut system.velocities;
    let f = &system.forces;
    let m = &system.masses;
    for i in 0..n {
        let inv_m = factor / m[i];
        for k in 0..3 {
            v[3 * i + k] += f[3 * i + k] * inv_m;
        }
    }
}
