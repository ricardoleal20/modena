//! `modena` — command-line interface for the MODENA MD engine.

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

use modena_core::runner::run_from_file;
use modena_core::{Ensemble, InputConfig, PressureStatus};

#[derive(Parser)]
#[command(
    name = "modena",
    version,
    about = "MODENA — Molecular Dynamics Engine for Numerical Atomistics"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Validate an input file without running (parses + builds the system).
    Check {
        #[arg(short, long)]
        input: PathBuf,
    },
    /// Run a simulation from an input file.
    Run {
        #[arg(short, long)]
        input: PathBuf,
        /// Emit a machine-readable energy diagnostics line (reduced energy/time).
        #[arg(long)]
        diagnostics_json: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::Check { input } => {
            let config = InputConfig::from_path(&input)
                .with_context(|| format!("failed to read input file {}", input.display()))?;
            let base_dir = input.parent().unwrap_or_else(|| std::path::Path::new("."));
            let sim = modena_core::Simulation::from_config_without_outputs(&config, base_dir)
                .context("failed to build simulation from config")?;
            println!(
                "OK: {} atoms, box = {:?}, ready to run {} steps.",
                sim.system().n(),
                sim.system().box_size.size.into_inner(),
                config.run.steps,
            );
        }
        Command::Run {
            input,
            diagnostics_json,
        } => {
            println!("MODENA: loading {}", input.display());
            let start = Instant::now();
            let summary = run_from_file(&input).context("simulation failed")?;
            let elapsed = start.elapsed();
            let e = &summary.final_energies;
            println!(
                "MODENA: finished {} steps in {:.3?}",
                summary.steps, elapsed
            );
            println!(
                "  final  KE = {:.6}   PE = {:.6}   E_total = {:.6}",
                e.kinetic, e.potential, e.total
            );
            println!(
                "  final  T = {:.6}   P = {:.6}   W = {:.6}",
                e.temperature, e.pressure, e.virial
            );
            println!(
                "  thermo units = {:?} (real: kcal/mol, K, atm; reduced: engine units)",
                summary.metadata.output_units
            );
            if summary.metadata.pressure_status
                == PressureStatus::DiagnosticWithoutConstraintImpulseStress
            {
                println!("  pressure diagnostic: omits constraint impulse stress; constrained pressure is not validated");
            }
            let c = &summary.conservation;
            if summary.metadata.ensemble == Ensemble::Nve {
                println!(
                    "  conservation  drift = {:.3e}  (legacy relative change; not the verdict basis)",
                    c.relative_drift
                );
            } else {
                println!(
                    "  energy exchanged  relative change = {:.3e}  (legacy; NVE verdict not applicable)",
                    c.relative_drift
                );
            }
            if let Some(d) = &c.diagnostics {
                let normalized = d
                    .normalized_max_abs_deviation
                    .map(|x| format!("{x:.3e}"))
                    .unwrap_or_else(|| "unavailable (no finite normalized metric)".into());
                println!(
                    "  energy-error advisory  max|delta E| / initial KE = {normalized}  ({:?})",
                    c.verdict
                );
                let slope = d
                    .slope
                    .map(|x| format!("{x:.3e}"))
                    .unwrap_or_else(|| "unavailable".into());
                println!("  energy statistics (reduced energy/time): delta/N={:.3e}  range={:.3e}  RMS={:.3e}  slope={slope}  samples={}",
                    d.delta_per_atom, d.energy_range, d.rms_fluctuation, d.samples);
            }
            if diagnostics_json {
                println!(
                    "  energy diagnostics = {}",
                    serde_json::to_string(c).context("failed to serialize energy diagnostics")?
                );
            }
            let per_step_us = elapsed.as_secs_f64() / summary.steps.max(1) as f64 * 1e6;
            println!(
                "  perf   {:.2} µs/step  ({:.0} steps/s)",
                per_step_us,
                1.0 / per_step_us * 1e6
            );
        }
    }
    Ok(())
}
