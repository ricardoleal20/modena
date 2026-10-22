//! One scientific validation boundary for every configuration entry point.
use std::collections::HashMap;

use crate::error::{ModenaError, Result};
use crate::io::input::*;
use crate::system::System;
use crate::units::UnitSystem;

fn error(message: impl Into<String>) -> ModenaError {
    ModenaError::Config(message.into())
}
fn positive(name: &str, v: f64) -> Result<()> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(error(format!("{name} must be finite and positive")))
    }
}
fn nonnegative(name: &str, v: f64) -> Result<()> {
    if v.is_finite() && v >= 0.0 {
        Ok(())
    } else {
        Err(error(format!("{name} must be finite and nonnegative")))
    }
}
fn finite_values(value: &toml::Value, path: &str) -> Result<()> {
    match value {
        toml::Value::Float(v) if !v.is_finite() => {
            return Err(error(format!("{path} must be finite")))
        }
        toml::Value::Array(v) => {
            for (i, x) in v.iter().enumerate() {
                finite_values(x, &format!("{path}[{i}]"))?;
            }
        }
        toml::Value::Table(v) => {
            for (key, x) in v {
                finite_values(x, &format!("{path}.{key}"))?;
            }
        }
        _ => {}
    }
    Ok(())
}
fn lattice_count(replicas: [u32; 3], per_cell: usize) -> Result<()> {
    let mut n = per_cell;
    for r in replicas {
        if r == 0 {
            return Err(error("structure replicas must be >= 1"));
        }
        n = n
            .checked_mul(r as usize)
            .ok_or_else(|| error("structure atom count overflow"))?;
    }
    // FCC/rocksalt generators form their total in u32 before allocating.
    if matches!(per_cell, 4 | 8) && n > u32::MAX as usize {
        return Err(error(
            "structure atom count exceeds lattice generator u32 capacity",
        ));
    }
    n.checked_mul(3)
        .and_then(|x| x.checked_mul(std::mem::size_of::<f64>()))
        .filter(|x| *x <= isize::MAX as usize)
        .ok_or_else(|| error("structure buffer size overflow"))?;
    Ok(())
}

pub(crate) fn validate(c: &InputConfig) -> Result<()> {
    // TOML's inspection value has signed integers; stochastic seeds remain
    // full u64 on JSON/Rust paths. They carry no floats needing this scan.
    let mut numeric_view = c.clone();
    numeric_view.run.seed = 0;
    match &mut numeric_view.structure {
        StructureSpec::RandomAlloy { seed, .. } | StructureSpec::Water { seed, .. } => *seed = 0,
        _ => {}
    }
    if let IntegratorSpec::LangevinBaoab { seed, .. } = &mut numeric_view.integrator {
        *seed = 0;
    }
    if let Some(ThermostatSpec::VelocityRescaling { seed, .. }) = &mut numeric_view.thermostat {
        *seed = 0;
    }
    if let Some(BarostatSpec::MonteCarlo { seed, .. }) = &mut numeric_view.barostat {
        *seed = 0;
    }
    let value = toml::Value::try_from(&numeric_view).map_err(|e| error(e.to_string()))?;
    finite_values(&value, "config")?;
    if let Some(f) = c.unit_factors {
        for x in [
            f.sigma,
            f.epsilon,
            f.mass_ref,
            f.time,
            f.temperature,
            f.pressure,
            f.charge,
        ] {
            positive("unit conversion factor", x)?;
        }
        if c.units != UnitSystem::Real {
            return Err(error("unit_factors require real units"));
        }
    }
    nonnegative("run.initial_temperature", c.run.initial_temperature)?;
    nonnegative("run.neighbor_skin", c.run.neighbor_skin)?;
    for (name, value) in [
        ("log_every", c.run.log_every),
        ("trajectory_every", c.run.trajectory_every),
        ("restart_every", c.run.restart_every),
    ] {
        if value == 0 {
            return Err(error(format!("run.{name} must be >= 1")));
        }
    }
    match &c.structure {
        StructureSpec::Fcc {
            lattice_constant,
            replicas,
            mass,
            ..
        } => {
            positive("lattice_constant", *lattice_constant)?;
            lattice_count(*replicas, 4)?;
            if let Some(m) = mass {
                positive("mass", *m)?;
            }
        }
        StructureSpec::Xyz { box_size, mass, .. } => {
            for edge in box_size {
                positive("box_size", *edge)?;
            }
            if let Some(m) = mass {
                positive("mass", *m)?;
            }
        }
        StructureSpec::RandomAlloy {
            lattice_constant,
            replicas,
            species,
            ..
        } => {
            positive("lattice_constant", *lattice_constant)?;
            lattice_count(*replicas, 4)?;
            if species.is_empty() {
                return Err(error("structure species must not be empty"));
            }
            for s in species {
                positive("species.mass", s.mass)?;
                nonnegative("species.concentration", s.concentration)?;
            }
            if (species.iter().map(|s| s.concentration).sum::<f64>() - 1.0).abs() > 1e-10 {
                return Err(error("species concentrations must sum to 1"));
            }
        }
        StructureSpec::RockSalt {
            lattice_constant,
            replicas,
            cation_mass,
            anion_mass,
            ..
        } => {
            positive("lattice_constant", *lattice_constant)?;
            lattice_count(*replicas, 8)?;
            for m in [cation_mass, anion_mass].into_iter().flatten() {
                positive("mass", *m)?;
            }
        }
        StructureSpec::Water {
            n_waters,
            density,
            d_oh,
            hoh_angle,
            m_o,
            m_h,
            ..
        } => {
            if c.units != UnitSystem::Real {
                return Err(error("water requires units = real"));
            }
            if *n_waters == 0 {
                return Err(error("n_waters must be >= 1"));
            }
            lattice_count([*n_waters, 1, 1], 3)?;
            for (name, x) in [
                ("density", *density),
                ("d_oh", *d_oh),
                ("m_o", *m_o),
                ("m_h", *m_h),
            ] {
                positive(name, x)?;
            }
            if *hoh_angle <= 0.0 || *hoh_angle >= 180.0 {
                return Err(error("hoh_angle must be in (0,180)"));
            }
        }
    }
    match &c.integrator {
        IntegratorSpec::VelocityVerlet { dt } => positive("dt", *dt)?,
        IntegratorSpec::LangevinBaoab {
            dt,
            gamma,
            target_temperature,
            ..
        } => {
            positive("dt", *dt)?;
            nonnegative("gamma", *gamma)?;
            nonnegative("target_temperature", *target_temperature)?;
        }
        IntegratorSpec::NoseHooverChain {
            dt,
            tau,
            target_temperature,
            chain_len,
        } => {
            positive("dt", *dt)?;
            positive("tau", *tau)?;
            positive("target_temperature", *target_temperature)?;
            if *chain_len == 0 {
                return Err(error("chain_len must be >= 1"));
            }
        }
    }
    if let Some(t) = &c.thermostat {
        match t {
            ThermostatSpec::None => {}
            ThermostatSpec::Berendsen {
                target_temperature,
                tau,
            }
            | ThermostatSpec::VelocityRescaling {
                target_temperature,
                tau,
                ..
            } => {
                if !matches!(c.integrator, IntegratorSpec::VelocityVerlet { .. }) {
                    return Err(error("Langevin and Nose-Hoover integrators cannot be combined with a separate thermostat"));
                }
                nonnegative("thermostat target_temperature", *target_temperature)?;
                positive("thermostat tau", *tau)?;
                if matches!(t, ThermostatSpec::VelocityRescaling { .. })
                    && *target_temperature > 0.0
                    && c.run.initial_temperature == 0.0
                    && c.run.restore_file.is_none()
                {
                    return Err(error(
                        "positive-temperature CSVR requires nonzero initial kinetic energy",
                    ));
                }
            }
        }
    }
    if let Some(b) = &c.barostat {
        if c.constraints.as_ref().is_some_and(|x| !x.bonds.is_empty())
            || matches!(c.structure, StructureSpec::Water { .. })
        {
            return Err(error("constrained box rescaling is unsupported: constraint impulse stress is not validated"));
        }
        match b {
            BarostatSpec::ParrinelloRahman { .. } => {
                return Err(error(
                    "legacy Parrinello-Rahman piston is unvalidated; use isotropic monte_carlo NPT",
                ))
            }
            BarostatSpec::Berendsen {
                tau,
                compressibility,
                anisotropic,
                ..
            } => {
                if *anisotropic {
                    return Err(error("anisotropic coupling requires a full stress tensor; scalar virial is insufficient"));
                }
                positive("barostat tau", *tau)?;
                positive("compressibility", *compressibility)?;
            }
            BarostatSpec::MonteCarlo {
                target_pressure,
                target_temperature,
                max_log_volume,
                interval,
                ..
            } => {
                positive("MC target_pressure", *target_pressure)?;
                positive("MC target_temperature", *target_temperature)?;
                positive("max_log_volume", *max_log_volume)?;
                if *interval == 0 {
                    return Err(error("MC interval must be >= 1"));
                }
                let canonical = match &c.integrator {
                    IntegratorSpec::LangevinBaoab {
                        target_temperature,
                        gamma,
                        ..
                    } if *gamma > 0.0 => Some(*target_temperature),
                    IntegratorSpec::NoseHooverChain {
                        target_temperature, ..
                    } => Some(*target_temperature),
                    IntegratorSpec::VelocityVerlet { .. } => match &c.thermostat {
                        Some(ThermostatSpec::VelocityRescaling {
                            target_temperature, ..
                        }) => Some(*target_temperature),
                        _ => None,
                    },
                    _ => None,
                };
                if canonical.is_none_or(|t| {
                    t <= 0.0
                        || (t - target_temperature).abs()
                            > 1e-12 * t.abs().max(target_temperature.abs())
                }) {
                    return Err(error("MC NPT requires a canonical thermostat/integrator at the same strictly positive target_temperature"));
                }
            }
        }
    }
    if let Some(e) = &c.electrostatics {
        for x in [e.kappa, e.real_cutoff, e.kspace_cutoff] {
            positive("electrostatics parameter", x)?;
        }
    }
    if let Some(cons) = &c.constraints {
        let mut pairs = HashMap::new();
        for b in &cons.bonds {
            positive("constraint length", b.length)?;
            if b.atoms[0] == b.atoms[1] {
                return Err(error("constraint atoms must differ"));
            }
            let key = (b.atoms[0].min(b.atoms[1]), b.atoms[0].max(b.atoms[1]));
            if pairs.insert(key, b.length).is_some() {
                return Err(error(
                    "duplicate constraint pair (constraints must be independent)",
                ));
            }
        }
    }
    // Constructor tables must not silently overwrite a conflicting override.
    if let Some(overrides) = value
        .get("model")
        .and_then(|v| v.get("overrides"))
        .and_then(|v| v.as_array())
    {
        let mut seen = HashMap::new();
        for o in overrides {
            let mut params = o.as_table().unwrap().clone();
            let types = params.remove("types").unwrap();
            let ids = types.as_array().unwrap();
            let i = ids[0].as_integer().unwrap();
            let j = ids[1].as_integer().unwrap();
            let key = (i.min(j), i.max(j));
            if seen
                .insert(key, params.clone())
                .is_some_and(|old| old != params)
            {
                return Err(error("contradictory pair overrides"));
            }
        }
    }
    Ok(())
}

pub(crate) fn validate_system(c: &InputConfig, s: &System) -> Result<()> {
    let n = s.n();
    if n == 0 {
        return Err(error("structure must contain atoms"));
    }
    for x in s.box_size.size.0 {
        positive("box edge", x)?;
    }
    positive("box volume", s.box_size.volume())?;
    for x in &s.masses {
        positive("atom mass", *x)?;
    }
    if s.positions
        .iter()
        .chain(&s.charges)
        .chain(&s.velocities)
        .any(|x| !x.is_finite())
        || !s.kinetic_energy().is_finite()
    {
        return Err(error(
            "atom positions, charges, velocities and kinetic energy must be finite",
        ));
    }
    let indices = |ids: &[usize]| -> Result<()> {
        if ids.iter().any(|i| *i >= n) || ids.iter().enumerate().any(|(i, a)| ids[..i].contains(a))
        {
            Err(error(
                "topology/constraint indices must be distinct and within atom count",
            ))
        } else {
            Ok(())
        }
    };
    if let Some(t) = &c.topology {
        let mut bonds = HashMap::new();
        for b in &t.bonds {
            indices(&b.atoms)?;
            positive("bond r0", b.r0)?;
            nonnegative("bond k", b.k)?;
            let key = (b.atoms[0].min(b.atoms[1]), b.atoms[0].max(b.atoms[1]));
            let p = (b.r0, b.k);
            if bonds.insert(key, p).is_some_and(|old| old != p) {
                return Err(error("contradictory bonded pair"));
            }
        }
        for a in &t.angles {
            indices(&a.atoms)?;
            nonnegative("angle k", a.k)?;
            if a.theta0 <= 0.0 || a.theta0 >= std::f64::consts::PI {
                return Err(error("angle theta0 must be in (0,pi)"));
            }
        }
        for d in &t.dihedrals {
            indices(&d.atoms)?;
            if d.terms.is_empty() {
                return Err(error("dihedral terms must not be empty"));
            }
            for p in &d.terms {
                if p.n <= 0 {
                    return Err(error("dihedral n must be positive"));
                }
                nonnegative("dihedral k", p.k)?;
            }
        }
    }
    if let Some(c) = &c.constraints {
        for b in &c.bonds {
            indices(&b.atoms)?;
        }
    }
    Ok(())
}
