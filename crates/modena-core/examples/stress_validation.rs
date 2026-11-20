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
    let sum = Command::new("/usr/bin/shasum")
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
        ("ase-six-strain", []) => {
            validation::ase::measure(dir)?;
            return Ok(0);
        }
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
fn host_python(args: &[String]) -> R<PathBuf> {
    if args.first().is_some_and(|p| p == "prepare") {
        return Ok(PathBuf::from("/usr/bin/python3"));
    }
    let out = args
        .windows(2)
        .find(|a| a[0] == "--out")
        .ok_or("missing campaign output")?;
    let dir = Path::new(&out[1]);
    let settings: serde_json::Value = read(&dir.join("inputs/settings.json"))?;
    let exe = PathBuf::from(
        settings["host_python"]
            .as_str()
            .ok_or("missing prepared interpreter")?,
    );
    if !exe.is_absolute() {
        return Err("interpreter path is not absolute".into());
    }
    let protocol: ValidationProtocol = read(&dir.join("protocol.json"))?;
    let path = exe.to_str().ok_or("unrepresentable interpreter")?;
    let expect = protocol
        .source_sha256
        .get(path)
        .ok_or("interpreter absent from source inventory")?;
    let hash = Command::new("/usr/bin/shasum")
        .args(["-a", "256", path])
        .output()
        .map_err(err)?;
    if !hash.status.success()
        || String::from_utf8(hash.stdout)
            .map_err(err)?
            .split_whitespace()
            .next()
            != Some(expect.as_str())
    {
        return Err("pinned interpreter changed".into());
    }
    Ok(exe)
}
fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = if args.first().map(String::as_str) == Some("--child") {
        child(&args[1..])
    } else {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        host_python(&args).and_then(|python| {
            Command::new("/bin/sh")
                .arg(
                    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                        .join("../../bench/b07_validation/launch.sh"),
                )
                .arg(python)
                .arg(root.join("bench/b07_validation/campaign.py"))
                .args(&args)
                .status()
                .map_err(err)
                .and_then(|s| s.code().ok_or("host launcher terminated".into()))
        })
    };
    match result {
        Ok(status) => ExitCode::from(status as u8),
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
