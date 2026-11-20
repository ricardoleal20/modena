//! Read a new raw lib-test observation of the private production projection operator.
use super::protocol::*;
use serde::Deserialize;
use std::path::Path;
type V = [f64; 9];
type M = [V; 9];
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Conditional {
    positions: V,
    masses: [f64; 3],
    velocities: V,
    decay: f64,
    #[serde(rename = "kT")]
    kt: f64,
    draws_by_seed: [usize; 8],
    total: usize,
    raw_mean: V,
    expected_raw_mean: V,
    reaction_mean: V,
    projected_mean: V,
    expected_projected_mean: V,
    raw_covariance: M,
    reaction_covariance: M,
    projected_covariance: M,
    expected_projected_covariance: M,
    normalized_covariance_error: M,
    covariance_se: M,
}
pub fn run_ou(out: &Path) -> R<Vec<GateResult>> {
    let text = std::fs::read_to_string(out.join("jobs/ou-draws/stdout.log")).map_err(err)?;
    let mut matches = text
        .lines()
        .filter_map(|l| l.split_once("conditional_statistics ").map(|(_, v)| v));
    let row: Conditional =
        serde_json::from_str(matches.next().ok_or("missing actual OU measurement")?)
            .map_err(err)?;
    if matches.next().is_some()
        || row.draws_by_seed != [8192; 8]
        || row.total != 65536
        || row.masses != [16., 1., 2.]
        || row.kt != 1.
        || row.decay != (-0.2_f64).exp()
    {
        return Err("conditional sample identity".into());
    }
    for v in [
        &row.positions,
        &row.velocities,
        &row.raw_mean,
        &row.expected_raw_mean,
        &row.reaction_mean,
        &row.projected_mean,
        &row.expected_projected_mean,
    ] {
        finite(v)?
    }
    for m in [
        &row.raw_covariance,
        &row.reaction_covariance,
        &row.projected_covariance,
        &row.expected_projected_covariance,
        &row.normalized_covariance_error,
        &row.covariance_se,
    ] {
        for r in m {
            finite(r)?
        }
    }
    let fixture: ValidationFixture = read(&out.join("inputs/water-nonprincipal-heavy.json"))?;
    let state = fixture.system()?;
    if state.positions != row.positions
        || state.velocities != row.velocities
        || state.masses != row.masses
    {
        return Err("OU snapshot differs from frozen nonprincipal water".into());
    }
    let n = row.total as f64;
    let mut result = vec![];
    for i in 0..9 {
        let rv = (1. - row.decay.powi(2)) * row.masses[i / 3];
        let pv = row.expected_projected_covariance[i][i];
        for (name, x, exp, var) in [
            ("raw", row.raw_mean[i], row.expected_raw_mean[i], rv),
            (
                "projected",
                row.projected_mean[i],
                row.expected_projected_mean[i],
                pv,
            ),
            ("reaction", row.reaction_mean[i], 0., rv - pv),
        ] {
            result.push(gate(
                &format!("ou_{name}_mean_{i}"),
                &[x],
                &[exp],
                &[4. * (var / n).sqrt() + 1e-3],
            )?);
        }
        for j in 0..9 {
            let scale = (pv * row.expected_projected_covariance[j][j]).sqrt();
            let se = row.covariance_se[i][j];
            finite(&[scale, se])?;
            result.push(gate(
                &format!("ou_cov_{i}_{j}"),
                &[
                    row.projected_covariance[i][j],
                    row.normalized_covariance_error[i][j],
                ],
                &[row.expected_projected_covariance[i][j], 0.],
                &[4. * se, 0.02],
            )?);
        }
    }
    Ok(result)
}
