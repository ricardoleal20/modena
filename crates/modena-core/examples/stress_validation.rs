//! Frozen B07 measurement child. Public phases enter through the host campaign controller.
#[path = "stress_validation/mod.rs"]
mod validation;
use std::{
    path::{Path, PathBuf},
    process::{Command, ExitCode},
};
use validation::{costs, ensemble, fixtures, matched, ou, protocol::*};
fn child(args: &[String]) -> R<i32> {
    let [phase, out, rest @ ..] = args else {
        return Err("expected phase and campaign output".into());
    };
    let dir = Path::new(out);
    if phase == "freeze-inputs" && rest.is_empty() {
        for f in matched::fixtures()? {
            f.system()?;
            write_new(&dir.join(format!("{}.json", f.id)), &f)?;
        }
        write_new(&dir.join("roster.json"), &matched::roster())?;
        return Ok(0);
    }
    let protocol = dir.join("protocol.json");
    let sum = Command::new("shasum")
        .args(["-a", "256"])
        .arg(&protocol)
        .output()
        .map_err(err)?;
    let identity = std::env::var("B07_PROTOCOL_SHA256").map_err(err)?;
    if !sum.status.success()
        || String::from_utf8(sum.stdout)
            .map_err(err)?
            .split_whitespace()
            .next()
            != Some(&identity)
    {
        return Err("campaign protocol differs from verified launch".into());
    }
    let p: ValidationProtocol = read(&protocol)?;
    p.validate()?;
    let fixture = |id: &str| {
        p.fixtures
            .iter()
            .find(|f| f.id == id)
            .ok_or_else(|| format!("unregistered fixture {id}"))
    };
    let g = match (phase.as_str(), rest) {
        ("deterministic", []) => fixtures::run_deterministic(dir)?,
        ("ou", []) => ou::run_ou(dir)?,
        ("ensemble", []) => ensemble::run_ensemble(dir)?,
        ("costs", []) => {
            costs::measure(dir, &p)?;
            return Ok(0);
        }
        ("matched-evaluate", [id, b, input, output]) => {
            matched::evaluate(fixture(id)?, b, Path::new(input), Path::new(output))?;
            return Ok(0);
        }
        ("matched-trajectory", [id, output]) => {
            let roster = matched::roster();
            let j = roster
                .iter()
                .find(|j| &j.id == id)
                .ok_or("unregistered trajectory")?;
            matched::trajectory(fixture(&j.fixture)?, j, Path::new(output))?;
            return Ok(0);
        }
        _ => return Err(format!("unsupported phase or arguments {args:?}")),
    };
    write_new(&dir.join(format!("{phase}-gates.json")), &g)?;
    Ok(if g.iter().any(|r| r.status == CaseStatus::Failed) {
        3
    } else if g.iter().any(|r| r.status != CaseStatus::Passed) {
        4
    } else {
        0
    })
}
fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = if args.first().map(String::as_str) == Some("--child") {
        child(&args[1..])
    } else {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        Command::new("python3")
            .arg(root.join("bench/b07_validation/campaign.py"))
            .args(&args)
            .status()
            .map_err(err)
            .and_then(|s| s.code().ok_or("host launcher terminated".into()))
    };
    match result {
        Ok(status) => ExitCode::from(status as u8),
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
