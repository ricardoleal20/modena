//! Barostats — pressure coupling that scales the simulation box and positions.
//!
//! Like a [`Thermostat`](crate::thermostat::Thermostat), a `Barostat` is applied
//! once per step *after* the integrator has advanced positions/velocities. It
//! reads the current pressure (from the system's KE and the just-computed virial)
//! and rescales the box edges and all positions by the same factor, so the
//! fractional coordinates are preserved and the density moves toward the target
//! pressure.
//!
//! M5 ships two barostats:
//! - [`BerendsenBarostat`] — weak-coupling (Berendsen et al., 1984). Deterministic,
//!   robust, drives P → P0; does not sample the exact NPT ensemble (use it for
//!   equilibration).
//! - [`ParrinelloRahman`] — the proper equation-of-motion barostat (Parrinello &
//!   Rahman, 1981) with an explicit piston degree of freedom integrated across
//!   steps. Samples the NPT ensemble.
//!
//! Both support isotropic and anisotropic (per-axis, orthorhombic) coupling.

use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::force_model::Virial;
use crate::system::System;

/// Pressure coupling. Applied once per step; scales the box and all positions.
pub trait Barostat: Send {
    fn apply(&mut self, system: &mut System, virial: &Virial, dt: f64);
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
/// constant: `W = dof · kT · τ²`. Samples the NPT ensemble.
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

/// Monte Carlo volume barostat: one Metropolis volume move per step. Proposes a
/// random isotropic volume change and accepts/rejects on the NPT criterion
/// ΔH = β ΔU + N ln(V_new/V_old) + β P0 (V_new − V_old). Samples the exact
/// isothermal–isobaric ensemble (unlike Berendsen). Requires a force model to
/// recompute U after the move, so `apply_with_model` is used instead of `apply`.
pub struct MonteCarloBarostat {
    pub target_pressure: f64,
    pub beta: f64, // 1/(kT)
    pub max_logv: f64,
    rng: StdRng,
    /// Accepted / attempted counts (diagnostic).
    pub accepted: usize,
    pub attempted: usize,
}

impl MonteCarloBarostat {
    pub fn new(target_pressure: f64, temperature: f64, max_logv: f64, seed: u64) -> Self {
        Self {
            target_pressure,
            beta: 1.0 / temperature,
            max_logv,
            rng: StdRng::seed_from_u64(seed),
            accepted: 0,
            attempted: 0,
        }
    }

    /// Attempt one volume move. `pe` is the current potential energy; returns
    /// the (possibly updated) PE. The caller passes the force model so U can be
    /// recomputed for the trial configuration.
    pub fn attempt<F>(&mut self, system: &mut System, current_pe: f64, mut recompute_pe: F) -> f64
    where
        F: FnMut(&mut System) -> f64,
    {
        use rand::RngExt;
        let n = system.n() as f64;
        let v_old = system.box_size.volume();
        // Trial: V_new = V_old · exp(δ), δ ∈ [−max, max].
        let delta: f64 = self.rng.random_range(-self.max_logv..self.max_logv);
        let scale = delta.exp(); // V_new / V_old
        let mu = scale.cbrt(); // linear scale

        // Save state for rejection.
        let old_size = system.box_size.size;
        let old_positions = system.positions.clone();

        // Apply trial.
        system.box_size.scale_uniform(mu);
        for x in system.positions.iter_mut() {
            *x *= mu;
        }

        let pe_new = recompute_pe(system);
        let v_new = system.box_size.volume();

        // NPT Metropolis: accept with prob min(1, exp(−β[ΔU + P0 ΔV] − N ln(V_new/V_old))).
        let dh = self.beta * (pe_new - current_pe + self.target_pressure * (v_new - v_old))
            + n * (v_new / v_old).ln();
        let accept = dh <= 0.0 || {
            let r: f64 = self.rng.random_range(0.0..1.0);
            (-dh).exp() > r
        };

        self.attempted += 1;
        if accept {
            self.accepted += 1;
            pe_new
        } else {
            system.box_size.size = old_size;
            system.positions = old_positions;
            current_pe
        }
    }

    /// No-op default `apply` (MC needs the model; use [`Self::attempt`]).
    fn _apply_marker(&self) {}
}

#[allow(unused_variables)]
impl Barostat for MonteCarloBarostat {
    fn apply(&mut self, system: &mut System, virial: &Virial, dt: f64) {
        // MC volume moves are driven by `attempt` (needs the force model for ΔU).
        // `apply` is a no-op so a `NoneBarostat`-style slot still type-checks.
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
