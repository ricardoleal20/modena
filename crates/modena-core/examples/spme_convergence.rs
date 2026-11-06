//! Frozen SPME calibration preparation. Later numerical phases are separate.
mod spme_calibration;
use spme_calibration::{frozen_protocol, prepare, CalibrationFailure};
fn run() -> Result<(), CalibrationFailure> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args == ["--help"] || args == ["-h"] {
        println!("spme_convergence prepare --out NEW_DIR --family cubic|ortho|special --replicas 2|3 --threads 1|4 [--reference-extended]");
        return Ok(());
    }
    if args.first().map(String::as_str) != Some("prepare") {
        return Err(CalibrationFailure::InvalidInput(
            "expected prepare command".into(),
        ));
    }
    let mut values = std::collections::BTreeMap::new();
    let mut extended = false;
    let mut i = 1;
    while i < args.len() {
        let key = &args[i];
        if key == "--reference-extended" {
            if extended {
                return Err(CalibrationFailure::InvalidInput(
                    "duplicate extended flag".into(),
                ));
            }
            extended = true;
            i += 1;
            continue;
        }
        if !["--out", "--family", "--replicas", "--threads"].contains(&key.as_str())
            || i + 1 >= args.len()
            || values.insert(key.as_str(), args[i + 1].as_str()).is_some()
        {
            return Err(CalibrationFailure::InvalidInput(format!(
                "unknown, duplicate or missing argument: {key}"
            )));
        }
        i += 2;
    }
    let required = |k| {
        values
            .get(k)
            .copied()
            .ok_or_else(|| CalibrationFailure::InvalidInput(format!("missing {k}")))
    };
    let replicas = required("--replicas")?
        .parse()
        .map_err(|_| CalibrationFailure::InvalidInput("invalid replicas".into()))?;
    let threads = required("--threads")?
        .parse()
        .map_err(|_| CalibrationFailure::InvalidInput("invalid threads".into()))?;
    let p = frozen_protocol(required("--family")?, replicas, threads, extended)?;
    let m = prepare(std::path::Path::new(required("--out")?), p)?;
    println!("Prepared {} {}", m.schema, m.snapshots_sha256);
    Ok(())
}
fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(e.exit_code());
    }
}
