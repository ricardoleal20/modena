//! Restart/checkpoint files.
//!
//! A checkpoint carries the **full mechanical state** of a run — box, per-atom
//! positions and velocities (flat SoA `3N`), and the global step counter — in
//! a simple, human-inspectable text format. Everything is in the engine's
//! reduced units; species/masses/charges are deliberately *not* carried (they
//! come from the structure config, whose `convert_to_reduced` already
//! normalizes them), so a restore is only valid against a config that builds
//! the same atom count (checked at load time in the runner).
//!
//! Format (`{:0.16e}` floats — 17 significant digits, exact f64 round-trip):
//!
//! ```text
//! MODENA_RESTART 1
//! atoms 192
//! step 200
//! time 1.0000000000000000e-1
//! box 1.234…e1 1.234…e1 1.234…e1
//! triclinic ax ay az bx by bz cx cy cz      (only for a triclinic cell)
//! positions
//! …3N lines, one float per line…
//! velocities
//! …3N lines, one float per line…
//! ```
//!
//! `time` is informational (`step · dt` at save time); the step is the
//! authoritative counter. 10k atoms → ~60k lines (~1.5 MB) — fine for text.

use std::fmt::Write as _;
use std::path::Path;

use crate::error::{ModenaError, Result};
use crate::system::SimBox;
use crate::vec3::Vec3;

/// Format tag (first token of the header line).
const MAGIC: &str = "MODENA_RESTART";
/// Current format version (second token of the header line).
const VERSION: u32 = 1;

/// A loaded/saved checkpoint. All values in reduced units.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    /// Atom count (sanity-checked against the restoring config's structure).
    pub atoms: usize,
    /// Global step counter at save time.
    pub step: usize,
    /// Reduced time `step · dt` at save time (informational; the step is
    /// authoritative — `dt` belongs to the config, not the checkpoint).
    pub time: f64,
    /// Simulation box (orthorhombic edges, or a triclinic cell).
    pub box_size: SimBox,
    /// Flat SoA positions, length `3 · atoms`.
    pub positions: Vec<f64>,
    /// Flat SoA velocities, length `3 · atoms`.
    pub velocities: Vec<f64>,
}

/// Serialize a checkpoint to `path` (atomically enough for our purposes:
/// truncated + written in one `fs::write`).
pub fn save<P: AsRef<Path>>(path: P, cp: &Checkpoint) -> Result<()> {
    debug_assert_eq!(cp.positions.len(), 3 * cp.atoms);
    debug_assert_eq!(cp.velocities.len(), 3 * cp.atoms);
    let mut out = String::with_capacity(64 * (6 * cp.atoms + 8));
    // `writeln!` into a String is infallible; `.ok()` discards the never-
    // produced `fmt::Error` without an unwrap in every call.
    writeln!(out, "{MAGIC} {VERSION}").ok();
    writeln!(out, "atoms {}", cp.atoms).ok();
    writeln!(out, "step {}", cp.step).ok();
    writeln!(out, "time {:0.16e}", cp.time).ok();
    let s = cp.box_size.size;
    writeln!(
        out,
        "box {:0.16e} {:0.16e} {:0.16e}",
        s.0[0], s.0[1], s.0[2]
    )
    .ok();
    if let Some([a, b, c]) = cp.box_size.vectors() {
        writeln!(
            out,
            "triclinic {:0.16e} {:0.16e} {:0.16e} {:0.16e} {:0.16e} {:0.16e} {:0.16e} {:0.16e} {:0.16e}",
            a.0[0], a.0[1], a.0[2], b.0[0], b.0[1], b.0[2], c.0[0], c.0[1], c.0[2]
        )
        .ok();
    }
    out.push_str("positions\n");
    for &p in &cp.positions {
        writeln!(out, "{p:0.16e}").ok();
    }
    out.push_str("velocities\n");
    for &v in &cp.velocities {
        writeln!(out, "{v:0.16e}").ok();
    }
    if let Some(parent) = path.as_ref().parent() {
        // The log/trajectory writers assume the target dir exists; be a bit
        // kinder for hand-written decks that point restart_file at results/.
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(path, out)?;
    Ok(())
}

/// Parse a checkpoint from `path`.
pub fn load<P: AsRef<Path>>(path: P) -> Result<Checkpoint> {
    let text = std::fs::read_to_string(path.as_ref())?;
    let mut lines = text.lines();

    let bad = |msg: &str| {
        ModenaError::Input(format!("restart file `{}`: {msg}", path.as_ref().display()))
    };

    // Header line: MAGIC VERSION.
    let header = lines.next().ok_or_else(|| bad("file is empty"))?;
    let mut header_parts = header.split_whitespace();
    if header_parts.next() != Some(MAGIC) {
        return Err(bad("not a MODENA restart file (bad magic line)"));
    }
    let version: u32 = header_parts
        .next()
        .ok_or_else(|| bad("missing format version"))?
        .parse()
        .map_err(|e| bad(&format!("bad format version: {e}")))?;
    if version != VERSION {
        return Err(bad(&format!(
            "unsupported version {version} (expected {VERSION})"
        )));
    }

    // Key/value header lines until the `positions` marker.
    let mut atoms: Option<usize> = None;
    let mut step = 0usize;
    let mut time = 0.0f64;
    let mut box_edges: Option<[f64; 3]> = None;
    let mut triclinic: Option<[Vec3; 3]> = None;
    loop {
        let line = lines
            .next()
            .ok_or_else(|| bad("missing `positions` block"))?;
        let mut it = line.split_whitespace();
        match it.next() {
            Some("atoms") => {
                atoms = Some(parse_usize(it.next(), &bad)?);
            }
            Some("step") => step = parse_usize(it.next(), &bad)?,
            Some("time") => time = parse_f64_tokens(&mut it, 1, &bad)?[0],
            Some("box") => {
                let e = parse_f64_tokens(&mut it, 3, &bad)?;
                box_edges = Some([e[0], e[1], e[2]]);
            }
            Some("triclinic") => {
                let v = parse_f64_tokens(&mut it, 9, &bad)?;
                triclinic = Some([
                    Vec3::new(v[0], v[1], v[2]),
                    Vec3::new(v[3], v[4], v[5]),
                    Vec3::new(v[6], v[7], v[8]),
                ]);
            }
            Some("positions") => break,
            _ => return Err(bad(&format!("malformed header line: {line}"))),
        }
    }
    let atoms = atoms.ok_or_else(|| bad("missing `atoms` header"))?;
    let edges = box_edges.ok_or_else(|| bad("missing `box` header"))?;
    let box_size = match triclinic {
        Some(v) => SimBox::triclinic(v),
        None => SimBox::from_array(edges),
    };

    let positions = parse_block(&mut lines, 3 * atoms, &bad)?;
    expect_marker(&mut lines, "velocities", &bad)?;
    let velocities = parse_block(&mut lines, 3 * atoms, &bad)?;

    Ok(Checkpoint {
        atoms,
        step,
        time,
        box_size,
        positions,
        velocities,
    })
}

/// Parse one integer (usize) header token.
fn parse_usize(tok: Option<&str>, bad: &impl Fn(&str) -> ModenaError) -> Result<usize> {
    tok.and_then(|t| t.parse().ok())
        .ok_or_else(|| bad("expected an unsigned integer header value"))
}

/// Parse exactly `n` whitespace-separated f64 tokens from a header line.
fn parse_f64_tokens(
    it: &mut std::str::SplitWhitespace<'_>,
    n: usize,
    bad: &impl Fn(&str) -> ModenaError,
) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let tok = it.next().ok_or_else(|| bad("truncated header line"))?;
        out.push(
            tok.parse()
                .map_err(|e| bad(&format!("bad float `{tok}`: {e}")))?,
        );
    }
    Ok(out)
}

/// Parse `n` one-float-per-line values.
fn parse_block<'a>(
    lines: &mut std::str::Lines<'a>,
    n: usize,
    bad: &impl Fn(&str) -> ModenaError,
) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let line = lines.next().ok_or_else(|| bad("truncated float block"))?;
        out.push(
            line.trim()
                .parse()
                .map_err(|e| bad(&format!("bad float line `{line}`: {e}")))?,
        );
    }
    Ok(out)
}

fn expect_marker<'a>(
    lines: &mut std::str::Lines<'a>,
    marker: &str,
    bad: &impl Fn(&str) -> ModenaError,
) -> Result<()> {
    match lines.next() {
        Some(l) if l.trim() == marker => Ok(()),
        _ => Err(bad(&format!("missing `{marker}` marker"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 17 significant digits round-trip exactly, including the awkward ones.
    #[test]
    fn float_roundtrip_is_exact() {
        let awkward = [
            0.1_f64 + 0.2,
            f64::MIN_POSITIVE,
            -1.0e-300,
            1234567890123456.7,
        ];
        for &v in &awkward {
            let mut s = String::new();
            writeln!(s, "{v:0.16e}").unwrap();
            let back: f64 = s.trim().parse().unwrap();
            assert_eq!(back.to_bits(), v.to_bits(), "{v} -> {s}");
        }
    }

    #[test]
    fn save_load_roundtrip() {
        let dir = std::env::temp_dir().join(format!(
            "modena_restart_ut_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cp.chk");
        let cp = Checkpoint {
            atoms: 2,
            step: 91,
            time: 0.455,
            box_size: SimBox::from_array([3.0, 4.0, 5.0]),
            positions: vec![0.1 + 0.2, 1.0, 2.0, 3.0, 4.0, 5.0],
            velocities: vec![-0.5, 0.25, 0.125, 1e-300, 6.0, 7.0],
        };
        save(&path, &cp).unwrap();
        let back = load(&path).unwrap();
        assert_eq!(back.atoms, cp.atoms);
        assert_eq!(back.step, cp.step);
        assert_eq!(back.time.to_bits(), cp.time.to_bits());
        assert_eq!(back.box_size, cp.box_size);
        assert_eq!(back.positions, cp.positions);
        assert_eq!(back.velocities, cp.velocities);
        // Triclinic cell round-trips through its vectors.
        let tri = Checkpoint {
            box_size: SimBox::triclinic([
                Vec3::new(2.0, 0.1, 0.0),
                Vec3::new(0.0, 2.0, 0.2),
                Vec3::new(0.05, 0.0, 2.0),
            ]),
            ..cp
        };
        save(&path, &tri).unwrap();
        let back = load(&path).unwrap();
        assert!(back.box_size.is_triclinic());
        assert_eq!(back.box_size, tri.box_size);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn load_rejects_garbage() {
        let dir = std::env::temp_dir().join(format!(
            "modena_restart_bad_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bad.chk");
        std::fs::write(&path, "hello world\n").unwrap();
        assert!(load(&path).is_err());
        std::fs::write(&path, "MODENA_RESTART 2\natoms 1\n").unwrap();
        assert!(load(&path).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
