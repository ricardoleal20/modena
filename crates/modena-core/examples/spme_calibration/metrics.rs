use super::*;
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
pub fn metrics(
    reference: &Evaluation,
    candidate: &Evaluation,
    n: usize,
    force_floor: f64,
) -> Result<ErrorMetrics, CalibrationFailure> {
    if n == 0
        || reference.forces.len() != n
        || candidate.forces.len() != n
        || !force_floor.is_finite()
        || force_floor < 0.
    {
        return Err(invalid("metrics shape/floor mismatch"));
    }
    for e in [reference, candidate] {
        if !e.energy.is_finite()
            || !e.strain_w.is_finite()
            || e.forces.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(invalid("metrics nonfinite evaluation"));
        }
    }
    let mut ref_norm = 0_f64;
    let mut err_norm = 0_f64;
    let mut maximum = 0_f64;
    let mut component = 0_f64;
    let mut net = [0.; 3];
    for (r, c) in reference.forces.iter().zip(&candidate.forces) {
        let d = std::array::from_fn(|k| c[k] - r[k]);
        ref_norm = ref_norm.hypot(norm(*r));
        let error = norm(d);
        err_norm = err_norm.hypot(error);
        maximum = maximum.max(error);
        for k in 0..3 {
            component = component.max(d[k].abs());
            net[k] += c[k];
        }
    }
    let denominator = ref_norm / (n as f64).sqrt();
    let absolute = err_norm / (n as f64).sqrt();
    let energy = (candidate.energy - reference.energy) / n as f64;
    let w = (candidate.strain_w - reference.strain_w) / n as f64;
    let result = ErrorMetrics {
        force_rms_abs: absolute,
        force_rms_rel: if denominator > force_floor {
            Some(absolute / denominator)
        } else {
            None
        },
        reference_force_rms: denominator,
        max_atom_force_error: maximum,
        max_component_force_error: component,
        energy_per_atom_signed: energy,
        energy_per_atom_abs: energy.abs(),
        w_per_atom_signed: w,
        w_per_atom_abs: w.abs(),
        net_force: net,
        net_force_norm: norm(net),
    };
    validate_metrics(&result)?;
    Ok(result)
}
fn validate_metrics(m: &ErrorMetrics) -> Result<(), CalibrationFailure> {
    let values = [
        m.force_rms_abs,
        m.reference_force_rms,
        m.max_atom_force_error,
        m.max_component_force_error,
        m.energy_per_atom_abs,
        m.w_per_atom_abs,
        m.net_force_norm,
    ];
    if values.iter().any(|x| !x.is_finite() || *x < 0.)
        || [m.energy_per_atom_signed, m.w_per_atom_signed]
            .iter()
            .chain(m.net_force.iter())
            .any(|x| !x.is_finite())
        || m.force_rms_rel.is_some_and(|x| !x.is_finite() || x < 0.)
    {
        return Err(invalid("nonfinite/negative metric"));
    }
    Ok(())
}
pub fn passes_target(
    m: &ErrorMetrics,
    target: f64,
    p: &Protocol,
    reference_converged: bool,
) -> Result<bool, CalibrationFailure> {
    validate_metrics(m)?;
    if !target.is_finite()
        || target <= 0.
        || !p.a.is_finite()
        || p.a <= 0.
        || !p.q_unit.is_finite()
        || p.q_unit <= 0.
    {
        return Err(invalid("invalid target/units"));
    }
    let f = p.q_unit * p.q_unit / (p.a * p.a);
    let e = p.q_unit * p.q_unit / p.a;
    let t = &p.tolerances;
    Ok(reference_converged
        && m.force_rms_rel.is_some_and(|r| r <= target)
        && m.force_rms_abs <= t.candidate_force_rms_abs_factor * target * f
        && m.max_atom_force_error <= t.candidate_force_max_factor * target * f
        && m.energy_per_atom_abs <= t.candidate_energy_factor * target * e
        && m.w_per_atom_abs <= t.candidate_w_factor * target * e)
}
