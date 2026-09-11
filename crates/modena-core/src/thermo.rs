//! Thermodynamic observables computed from a `System` + virial.

use crate::force_model::Virial;
use crate::system::System;

#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct Energies {
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
    pub pressure: f64,
    pub virial: f64,
}

impl Energies {
    pub fn from(system: &System, potential: f64, virial: &Virial) -> Self {
        let kinetic = system.kinetic_energy();
        let total = kinetic + potential;
        let temperature = system.temperature();
        let volume = system.box_size.volume();
        let pressure = if volume > 0.0 {
            (2.0 * kinetic + virial.w) / (3.0 * volume)
        } else {
            0.0
        };
        Self {
            kinetic,
            potential,
            total,
            temperature,
            pressure,
            virial: virial.w,
        }
    }
}
