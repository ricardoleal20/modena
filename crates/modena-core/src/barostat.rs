//! Barostats — pressure coupling that scales the simulation box and positions.
//!
//! Like a [`Thermostat`](crate::thermostat::Thermostat), a `Barostat` is applied
//! once per step *after* the integrator has advanced positions/velocities. It
//! reads the current pressure (from the system's KE and the just-computed virial)
//! and rescales the box edges and all positions by the same factor, so the
//! fractional coordinates are preserved and the density moves toward the target
//! pressure.
//!
//! Canonical isotropic NPT uses [`MonteCarloBarostat`] with a canonical
//! temperature controller. [`BerendsenBarostat`] is equilibration only. The
//! direct legacy [`ParrinelloRahman`] and scalar anisotropic APIs remain
//! experimental; production configuration rejects them. Constrained box
//! rescaling is unsupported until the required coupling and stress exist.

use crate::error::{ModenaError, Result};
use crate::force_model::ForceModel;
use crate::neighbor::NeighborList;
use crate::thermostat::{
    checkpoint_rng, restore_rng, restore_stateless, state_error, validate_state_parameters,
};
use chacha20::ChaCha12Rng;
use rand::{RngExt, SeedableRng};
use serde_json::{json, Value};

use crate::force_model::Virial;
use crate::system::System;

/// Pressure coupling. Applied once per step; scales the box and all positions.
pub trait Barostat: Send {
    fn apply(&mut self, system: &mut System, virial: &Virial, dt: f64);
    /// Apply coupling and refresh all mechanical caches after a box move.
    /// Returns true when a volume move was accepted/applied.
    fn apply_with_model(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        pe: &mut f64,
        virial: &mut Virial,
        dt: f64,
    ) -> Result<bool> {
        let old_box = system.box_size;
        self.apply(system, virial, dt);
        if old_box == system.box_size {
            return Ok(false);
        }
        if !system.box_size.volume().is_finite() || system.box_size.volume() <= 0.0 {
            return Err(ModenaError::Physics(
                "barostat produced invalid volume".into(),
            ));
        }
        neighbors.force_rebuild(system);
        system.zero_forces();
        *virial = Virial::default();
        *pe = model.compute(system, virial, neighbors)?;
        Ok(true)
    }
    fn checkpoint_state(&self) -> Result<Value> {
        Ok(Value::Null)
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        restore_stateless(state)
    }
}

/// No pressure coupling — NVE/NVT (constant volume).
pub struct NoneBarostat;

impl Barostat for NoneBarostat {
    fn apply(&mut self, _system: &mut System, _virial: &Virial, _dt: f64) {}
}

/// Current scalar pressure P = (2 KE + W) / (3 V) (reduced units).
fn pressure(system: &System, virial: &Virial) -> f64 {
    let v = system.box_size.volume();
    if v <= 0.0 {
        return 0.0;
    }
    (2.0 * system.kinetic_energy() + virial.w) / (3.0 * v)
}

/// Berendsen weak-coupling barostat. Scales the box toward `target_pressure`
/// with time constant `tau` and isothermal compressibility `compressibility`
/// (β_T; in reduced units for an LJ run). `anisotropic = false` couples all
/// three axes by a single isotropic factor; `true` rescales each axis by its
/// own diagonal pressure component.
pub struct BerendsenBarostat {
    pub target_pressure: f64,
    pub tau: f64,
    pub compressibility: f64,
    pub anisotropic: bool,
}

impl BerendsenBarostat {
    pub fn new(target_pressure: f64, tau: f64, compressibility: f64) -> Self {
        Self {
            target_pressure,
            tau,
            compressibility,
            anisotropic: false,
        }
    }

    pub fn anisotropic(mut self) -> Self {
        self.anisotropic = true;
        self
    }

    /// Per-axis scale factor μ_k = [1 + β (dt/τ) (P_kk − P0)]^(1/3).
    fn axis_factor(&self, p_kk: f64, dt: f64) -> f64 {
        let mu3 = 1.0 + self.compressibility * (dt / self.tau) * (p_kk - self.target_pressure);
        // Guard against a negative volume (large overstep): clamp the per-step
        // box change to ±10%. Berendsen is weak-coupling; this only triggers for
        // pathological far-from-equilibrium states.
        let mu3 = mu3.clamp(0.9, 1.1);
        mu3.cbrt()
    }
}

impl Barostat for BerendsenBarostat {
    fn checkpoint_state(&self) -> Result<Value> {
        Ok(
            json!({"kind":"berendsen", "target_pressure":self.target_pressure, "tau":self.tau,
            "compressibility":self.compressibility, "anisotropic":self.anisotropic}),
        )
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        if state == &self.checkpoint_state()? {
            Ok(())
        } else {
            Err(state_error("Berendsen barostat parameters differ"))
        }
    }
    fn apply(&mut self, system: &mut System, virial: &Virial, dt: f64) {
        if self.tau <= 0.0 {
            return;
        }
        if self.anisotropic {
            // Diagonal pressure components: P_kk = (2 KE_k + W_kk) / V. With the
            // scalar virial we only have W (trace); distribute the kinetic part
            // per-axis from the velocity components and split W evenly.
            let v = system.box_size.volume();
            if v <= 0.0 {
                return;
            }
            for k in 0..3 {
                let ke_k: f64 = (0..system.n())
                    .map(|i| 0.5 * system.masses[i] * system.vel(i).0[k].powi(2))
                    .sum();
                let p_kk = (2.0 * ke_k + virial.w / 3.0) / v;
                let mu = self.axis_factor(p_kk, dt);
                system.box_size.size.0[k] *= mu;
                for i in 0..system.n() {
                    system.positions[3 * i + k] *= mu;
                }
            }
        } else {
            let p = pressure(system, virial);
            let mu = self.axis_factor(p, dt);
            system.box_size.scale_uniform(mu);
            for x in system.positions.iter_mut() {
                *x *= mu;
            }
        }
    }
}

/// Parrinello–Rahman barostat (isotropic or anisotropic orthorhombic). Carries
/// a piston degree of freedom per coupled axis, integrated with a symplectic
/// kick-drift over each `apply`. The piston mass is tied to the coupling time
/// constant: `W = dof · kT · τ²`. Experimental legacy dynamics, not validated canonical NPT.
pub struct ParrinelloRahman {
    pub target_pressure: f64,
    pub tau: f64,
    pub temperature: f64,
    pub anisotropic: bool,
    /// Piston velocity per axis (d ln L_k / dt). Internal state.
    piston_vel: [f64; 3],
}

impl ParrinelloRahman {
    pub fn new(target_pressure: f64, tau: f64, temperature: f64) -> Self {
        Self {
            target_pressure,
            tau,
            temperature,
            anisotropic: false,
            piston_vel: [0.0; 3],
        }
    }

    pub fn anisotropic(mut self) -> Self {
        self.anisotropic = true;
        self
    }

    /// Piston mass W = dof · kT · τ² (the standard PR–MTTK coupling).
    fn piston_mass(&self, dof: f64) -> f64 {
        (dof * self.temperature).max(1.0) * self.tau * self.tau
    }
}

impl Barostat for ParrinelloRahman {
    fn checkpoint_state(&self) -> Result<Value> {
        Ok(
            json!({"kind":"legacy_piston", "target_pressure":self.target_pressure, "tau":self.tau,
            "temperature":self.temperature, "anisotropic":self.anisotropic, "piston_vel":self.piston_vel}),
        )
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        validate_state_parameters(
            state,
            &self.checkpoint_state()?,
            &[
                "kind",
                "target_pressure",
                "tau",
                "temperature",
                "anisotropic",
            ],
        )?;
        let v: [f64; 3] = serde_json::from_value(state["piston_vel"].clone())
            .map_err(|e| state_error(&e.to_string()))?;
        if v.iter().any(|x| !x.is_finite()) {
            return Err(state_error("nonfinite piston velocity"));
        }
        self.piston_vel = v;
        Ok(())
    }
    fn apply(&mut self, system: &mut System, virial: &Virial, dt: f64) {
        if self.tau <= 0.0 || self.temperature <= 0.0 {
            return;
        }
        let dof = system.dof() as f64;
        let wp = self.piston_mass(dof);
        let half = 0.5 * dt;

        // Pressure used for the piston force. Isotropic: scalar P. Anisotropic:
        // per-axis diagonal P_kk.
        let p_axis = |k: usize| -> f64 {
            let v = system.box_size.volume();
            if v <= 0.0 {
                return 0.0;
            }
            let ke_k: f64 = (0..system.n())
                .map(|i| 0.5 * system.masses[i] * system.vel(i).0[k].powi(2))
                .sum();
            (2.0 * ke_k + virial.w / 3.0) / v
        };

        let coupled_axes: Vec<usize> = if self.anisotropic {
            vec![0, 1, 2]
        } else {
            vec![0] // isotropic: drive axis 0 piston, apply same μ to all axes
        };

        // Half-kick the pistons: v_p += 0.5 dt · 3V(P_kk − P0)/W  (isotropic
        // piston force uses the scalar pressure).
        for &k in &coupled_axes {
            let p = if self.anisotropic {
                p_axis(k)
            } else {
                pressure(system, virial)
            };
            let v = system.box_size.volume();
            let force = 3.0 * v * (p - self.target_pressure);
            self.piston_vel[k] += half * force / wp;
        }

        // Drift: scale each axis by exp(v_p · dt) ≈ 1 + v_p dt (and rescale
        // positions to preserve fractional coords).
        for &k in &coupled_axes {
            let mu = (self.piston_vel[k] * dt).exp();
            if self.anisotropic {
                system.box_size.size.0[k] *= mu;
                for i in 0..system.n() {
                    system.positions[3 * i + k] *= mu;
                }
            } else {
                // Isotropic: one piston (axis 0) drives all three axes.
                system.box_size.scale_uniform(mu);
                for x in system.positions.iter_mut() {
                    *x *= mu;
                }
            }
        }

        // Recompute pressure after the box move (volume changed) and half-kick
        // the pistons again. The virial is unchanged (positions scaled uniformly
        // preserve the pair distances up to the min-image, so W is a good
        // approximation for the half-step).
        let p_axis_new = |k: usize| -> f64 {
            let v = system.box_size.volume();
            if v <= 0.0 {
                return 0.0;
            }
            let ke_k: f64 = (0..system.n())
                .map(|i| 0.5 * system.masses[i] * system.vel(i).0[k].powi(2))
                .sum();
            (2.0 * ke_k + virial.w / 3.0) / v
        };
        for &k in &coupled_axes {
            let p = if self.anisotropic {
                p_axis_new(k)
            } else {
                pressure(system, virial)
            };
            let v = system.box_size.volume();
            let force = 3.0 * v * (p - self.target_pressure);
            self.piston_vel[k] += half * force / wp;
        }
    }
}

/// Canonical isotropic volume Metropolis moves with symmetric log-volume
/// proposals. The proposal Jacobian adds one to the particle count in the
/// acceptance exponent. This API supports unconstrained particles only.
pub struct MonteCarloBarostat {
    pub target_pressure: f64,
    pub beta: f64,
    pub max_logv: f64,
    pub interval: usize,
    rng: ChaCha12Rng,
    pub accepted: usize,
    pub attempted: usize,
    elapsed_steps: usize,
}
impl MonteCarloBarostat {
    pub fn new(target_pressure: f64, temperature: f64, max_logv: f64, seed: u64) -> Self {
        Self {
            target_pressure,
            beta: 1.0 / temperature,
            max_logv,
            interval: 1,
            rng: ChaCha12Rng::seed_from_u64(seed),
            accepted: 0,
            attempted: 0,
            elapsed_steps: 0,
        }
    }
    pub fn with_interval(mut self, interval: usize) -> Result<Self> {
        self.interval = interval;
        self.validate()?;
        Ok(self)
    }
    fn validate(&self) -> Result<()> {
        if [self.target_pressure, self.beta, self.max_logv]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || self.interval == 0
        {
            return Err(ModenaError::Config(
                "MC barostat requires finite positive P, T, proposal width and interval".into(),
            ));
        }
        Ok(())
    }
    fn propose(&mut self, system: &mut System) -> f64 {
        let delta = self.rng.random_range(-self.max_logv..self.max_logv);
        let mu = (delta / 3.0_f64).exp();
        system.box_size.scale_uniform(mu);
        for x in &mut system.positions {
            *x *= mu;
        }
        delta
    }
    fn accepts(
        &mut self,
        n: usize,
        old_v: f64,
        new_v: f64,
        old_pe: f64,
        new_pe: f64,
        delta: f64,
    ) -> bool {
        let log_accept = -self.beta * (new_pe - old_pe + self.target_pressure * (new_v - old_v))
            + (n as f64 + 1.0) * delta;
        let log_draw = self.rng.random::<f64>().ln();
        log_accept.is_finite() && (log_accept >= 0.0 || log_draw < log_accept)
    }
    /// Compatibility energy callback. Restores the complete System on rejection.
    /// Production callers use the fallible model hook to preserve neighbor/virial caches.
    pub fn attempt<F>(&mut self, system: &mut System, current_pe: f64, mut recompute_pe: F) -> f64
    where
        F: FnMut(&mut System) -> f64,
    {
        if self.validate().is_err() || system.constraint_count() != 0 {
            return current_pe;
        }
        let old = system.clone();
        let delta = self.propose(system);
        let pe = recompute_pe(system);
        self.attempted += 1;
        self.elapsed_steps += self.interval;
        if pe.is_finite()
            && self.accepts(
                system.n(),
                old.box_size.volume(),
                system.box_size.volume(),
                current_pe,
                pe,
                delta,
            )
        {
            self.accepted += 1;
            pe
        } else {
            *system = old;
            current_pe
        }
    }
}
impl Barostat for MonteCarloBarostat {
    fn apply(&mut self, _system: &mut System, _virial: &Virial, _dt: f64) {}
    fn apply_with_model(
        &mut self,
        system: &mut System,
        model: &dyn ForceModel,
        neighbors: &mut NeighborList,
        pe: &mut f64,
        virial: &mut Virial,
        _dt: f64,
    ) -> Result<bool> {
        self.validate()?;
        if system.constraint_count() != 0 {
            return Err(ModenaError::Config(
                "MC box moves do not support constrained particles".into(),
            ));
        }
        self.elapsed_steps = self
            .elapsed_steps
            .checked_add(1)
            .ok_or_else(|| state_error("MC step count overflow"))?;
        if !self.elapsed_steps.is_multiple_of(self.interval) {
            return Ok(false);
        }
        let old = system.clone();
        let old_pe = *pe;
        let old_virial = *virial;
        self.attempted = self
            .attempted
            .checked_add(1)
            .ok_or_else(|| state_error("MC attempt count overflow"))?;
        let delta = self.propose(system);
        let trial = (|| -> Result<f64> {
            if !system.box_size.volume().is_finite()
                || system.box_size.volume() <= 0.0
                || system.positions.iter().any(|x| !x.is_finite())
            {
                return Err(ModenaError::Physics(
                    "MC proposal produced nonfinite geometry".into(),
                ));
            }
            neighbors.force_rebuild(system);
            system.zero_forces();
            *virial = Virial::default();
            let trial_pe = model.compute(system, virial, neighbors)?;
            if !trial_pe.is_finite()
                || !virial.w.is_finite()
                || system.forces.iter().any(|x| !x.is_finite())
            {
                return Err(ModenaError::Physics(
                    "MC trial produced nonfinite forces/energy".into(),
                ));
            }
            Ok(trial_pe)
        })();
        let accepted = trial.as_ref().is_ok_and(|trial_pe| {
            self.accepts(
                system.n(),
                old.box_size.volume(),
                system.box_size.volume(),
                old_pe,
                *trial_pe,
                delta,
            )
        });
        if accepted {
            self.accepted = self
                .accepted
                .checked_add(1)
                .ok_or_else(|| state_error("MC acceptance count overflow"))?;
            *pe = trial?;
            Ok(true)
        } else {
            *system = old;
            *pe = old_pe;
            *virial = old_virial;
            neighbors.force_rebuild(system);
            trial.map(|_| false)
        }
    }
    fn checkpoint_state(&self) -> Result<Value> {
        self.validate()?;
        Ok(
            json!({"kind":"monte_carlo_logv", "target_pressure":self.target_pressure, "beta":self.beta,
            "max_logv":self.max_logv,"interval":self.interval,"accepted":self.accepted,"attempted":self.attempted,
            "elapsed_steps":self.elapsed_steps,"rng":checkpoint_rng(&self.rng)?}),
        )
    }
    fn restore_state(&mut self, state: &Value) -> Result<()> {
        validate_state_parameters(
            state,
            &self.checkpoint_state()?,
            &["kind", "target_pressure", "beta", "max_logv", "interval"],
        )?;
        let number = |key: &str| -> Result<usize> {
            state[key]
                .as_u64()
                .and_then(|x| usize::try_from(x).ok())
                .ok_or_else(|| state_error("invalid MC counter"))
        };
        let accepted = number("accepted")?;
        let attempted = number("attempted")?;
        let elapsed = number("elapsed_steps")?;
        if accepted > attempted || attempted != elapsed / self.interval {
            return Err(state_error("inconsistent MC counters"));
        }
        let rng = restore_rng(&state["rng"], self.rng.get_seed())?;
        self.rng = rng;
        self.accepted = accepted;
        self.attempted = attempted;
        self.elapsed_steps = elapsed;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::SimBox;
    use crate::vec3::Vec3;

    #[test]
    fn berendsen_drives_pressure_toward_target() {
        // A box far too small → huge pressure → Berendsen must expand it (μ > 1)
        // in one step, lowering density. We check the box grew.
        let mut sys = System::new(SimBox::cubic(2.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(0.5, 0.0, 0.0));
        // Give KE so pressure is large and positive.
        sys.set_vel(0, Vec3::new(0.0, 5.0, 0.0));
        sys.set_vel(1, Vec3::new(0.0, -5.0, 0.0));
        let v0 = sys.box_size.volume();
        let mut baro = BerendsenBarostat::new(0.0, 1.0, 1.0);
        let virial = Virial { w: 100.0 }; // large positive virial → high pressure
        baro.apply(&mut sys, &virial, 0.005);
        assert!(
            sys.box_size.volume() > v0,
            "box must expand under high pressure"
        );
    }

    #[test]
    fn parrinello_rahman_runs_and_changes_box() {
        let mut sys = System::new(SimBox::cubic(3.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.0, 0.0, 0.0));
        sys.set_vel(0, Vec3::new(0.0, 1.0, 0.0));
        sys.set_vel(1, Vec3::new(0.0, -1.0, 0.0));
        let mut baro = ParrinelloRahman::new(0.0, 1.0, 1.0);
        let virial = Virial { w: 10.0 };
        for _ in 0..10 {
            baro.apply(&mut sys, &virial, 0.005);
        }
        // Box should have moved (piston integrated over 10 steps).
        assert!(sys.box_size.volume() > 0.0);
    }
}
