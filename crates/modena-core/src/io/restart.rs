//! Versioned, validated and atomic scientific continuation checkpoints.
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::error::{ModenaError, Result};
use crate::system::SimBox;
use crate::vec3::Vec3;

/// Format tag (first token of the header line).
const MAGIC: &str = "MODENA_RESTART";
/// Current format version (second token of the header line).
const VERSION: u32 = 2;
static TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);

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
    /// None for legacy v1 or mechanical-only imports. Full continuation needs v2.
    pub continuation: Option<ContinuationState>,
}

/// All identity, cached mechanics, stochastic and extended variables needed
/// for exact continuation. Physical signatures are normalized configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContinuationState {
    pub forces: Vec<f64>,
    pub masses: Vec<f64>,
    pub species: Vec<u32>,
    pub charges: Vec<f64>,
    pub elements: Vec<String>,
    pub constraint_count: usize,
    pub conserved_com_dof: usize,
    pub potential_energy: f64,
    pub virial: f64,
    pub physical_signature: Value,
    pub integrator: Value,
    pub thermostat: Value,
    pub barostat: Value,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Version2 {
    atoms: usize,
    step: usize,
    time: f64,
    box_edges: [f64; 3],
    box_vectors: Option<[[f64; 3]; 3]>,
    positions: Vec<f64>,
    velocities: Vec<f64>,
    continuation: Option<ContinuationState>,
}

/// Write to a unique sibling, sync the contents, atomically rename, then sync
/// the containing directory. Validation/serialization failures preserve the old file.
pub fn save<P: AsRef<Path>>(path: P, cp: &Checkpoint) -> Result<()> {
    validate(cp)?;
    let payload = Version2 {
        atoms: cp.atoms,
        step: cp.step,
        time: cp.time,
        box_edges: cp.box_size.size.0,
        box_vectors: cp.box_size.vectors().map(|v| v.map(|x| x.0)),
        positions: cp.positions.clone(),
        velocities: cp.velocities.clone(),
        continuation: cp.continuation.clone(),
    };
    let json = serde_json::to_vec(&payload)
        .map_err(|e| ModenaError::Input(format!("checkpoint serialization: {e}")))?;
    let path = path.as_ref();
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    std::fs::create_dir_all(parent)?;
    let sequence = TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let file_name = path
        .file_name()
        .ok_or_else(|| ModenaError::Input("checkpoint path has no filename".into()))?;
    let temp = parent.join(format!(
        ".{}.{}.{}.tmp",
        file_name.to_string_lossy(),
        std::process::id(),
        sequence
    ));
    let mut created = false;
    let result = (|| -> Result<()> {
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp)?;
        created = true;
        writeln!(file, "{MAGIC} {VERSION}")?;
        file.write_all(&json)?;
        file.write_all(b"\n")?;
        file.sync_all()?;
        drop(file);
        std::fs::rename(&temp, path)?;
        std::fs::File::open(parent)?.sync_all()?;
        Ok(())
    })();
    if result.is_err() && created {
        let _ = std::fs::remove_file(&temp);
    }
    result
}

fn validate(cp: &Checkpoint) -> Result<()> {
    let bad = |s: &str| ModenaError::Input(format!("checkpoint: {s}"));
    let count = cp
        .atoms
        .checked_mul(3)
        .ok_or_else(|| bad("atom buffer size overflow"))?;
    if cp.atoms == 0 || cp.positions.len() != count || cp.velocities.len() != count {
        return Err(bad("invalid atom count or buffer lengths"));
    }
    if !cp.time.is_finite()
        || cp.time < 0.0
        || cp
            .positions
            .iter()
            .chain(&cp.velocities)
            .any(|x| !x.is_finite())
    {
        return Err(bad("nonfinite or invalid mechanical state"));
    }
    if cp
        .box_size
        .size
        .0
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !cp.box_size.volume().is_finite()
        || cp.box_size.volume() <= 0.0
        || cp
            .box_size
            .inv_basis_rows()
            .iter()
            .flat_map(|x| x.0)
            .any(|x| !x.is_finite())
    {
        return Err(bad("box must have finite positive edges and volume"));
    }
    if let Some(c) = &cp.continuation {
        if c.forces.len() != count
            || c.masses.len() != cp.atoms
            || c.species.len() != cp.atoms
            || c.charges.len() != cp.atoms
            || c.elements.len() != cp.atoms
        {
            return Err(bad("continuation atom buffers have inconsistent lengths"));
        }
        if c.forces.iter().chain(&c.charges).any(|x| !x.is_finite())
            || c.masses.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || !c.potential_energy.is_finite()
            || !c.virial.is_finite()
        {
            return Err(bad(
                "invalid continuation masses, energy, forces or charges",
            ));
        }
        if c.constraint_count
            .checked_add(c.conserved_com_dof)
            .is_none_or(|x| x > count)
            || ![0, 3].contains(&c.conserved_com_dof)
            || !c.physical_signature.is_object()
        {
            return Err(bad("invalid physical DOF or signature"));
        }
    }
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
    if header_parts.next().is_some() {
        return Err(bad("trailing header tokens"));
    }
    if version == VERSION {
        let payload: Version2 = serde_json::from_str(&lines.collect::<Vec<_>>().join("\n"))
            .map_err(|e| bad(&format!("malformed v2 payload: {e}")))?;
        let box_size = match payload.box_vectors {
            Some(v) => {
                let vectors = v.map(|a| Vec3::new(a[0], a[1], a[2]));
                if v.iter().flatten().any(|x| !x.is_finite())
                    || !vectors[0].dot(vectors[1].cross(vectors[2])).is_finite()
                    || vectors[0].dot(vectors[1].cross(vectors[2])).abs() == 0.0
                {
                    return Err(bad("invalid triclinic basis"));
                }
                let box_size = SimBox::triclinic(vectors);
                if box_size.size.0 != payload.box_edges {
                    return Err(bad("inconsistent triclinic box edges"));
                }
                box_size
            }
            None => SimBox::from_array(payload.box_edges),
        };
        let cp = Checkpoint {
            atoms: payload.atoms,
            step: payload.step,
            time: payload.time,
            box_size,
            positions: payload.positions,
            velocities: payload.velocities,
            continuation: payload.continuation,
        };
        validate(&cp)?;
        return Ok(cp);
    }
    if version != 1 {
        return Err(bad(&format!(
            "unsupported version {version} (expected {VERSION})"
        )));
    }

    // Key/value header lines until the `positions` marker.
    let mut atoms: Option<usize> = None;
    let mut step = None;
    let mut time = None;
    let mut keys = std::collections::HashSet::new();
    let mut box_edges: Option<[f64; 3]> = None;
    let mut triclinic: Option<[Vec3; 3]> = None;
    loop {
        let line = lines
            .next()
            .ok_or_else(|| bad("missing `positions` block"))?;
        let mut it = line.split_whitespace();
        let key = it.next();
        if !keys.insert(key) {
            return Err(bad("duplicate header key"));
        }
        match key {
            Some("atoms") => {
                atoms = Some(parse_usize(it.next(), &bad)?);
            }
            Some("step") => step = Some(parse_usize(it.next(), &bad)?),
            Some("time") => time = Some(parse_f64_tokens(&mut it, 1, &bad)?[0]),
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
            Some("positions") => {
                if it.next().is_some() {
                    return Err(bad("trailing marker tokens"));
                }
                break;
            }
            _ => return Err(bad(&format!("malformed header line: {line}"))),
        }
        if it.next().is_some() {
            return Err(bad("trailing header values"));
        }
    }
    let atoms = atoms.ok_or_else(|| bad("missing `atoms` header"))?;
    let edges = box_edges.ok_or_else(|| bad("missing `box` header"))?;
    let box_size = match triclinic {
        Some(v) => {
            if v.iter().flat_map(|x| x.0).any(|x| !x.is_finite())
                || !v[0].dot(v[1].cross(v[2])).is_finite()
                || v[0].dot(v[1].cross(v[2])).abs() == 0.0
            {
                return Err(bad("invalid triclinic basis"));
            }
            let b = SimBox::triclinic(v);
            if b.size.0 != edges {
                return Err(bad("inconsistent triclinic box edges"));
            }
            b
        }
        None => SimBox::from_array(edges),
    };

    let count = atoms
        .checked_mul(3)
        .ok_or_else(|| bad("atom buffer size overflow"))?;
    if count > lines.clone().count() {
        return Err(bad("atom count exceeds available data"));
    }
    let positions = parse_block(&mut lines, count, &bad)?;
    expect_marker(&mut lines, "velocities", &bad)?;
    let velocities = parse_block(&mut lines, count, &bad)?;
    if lines.any(|line| !line.trim().is_empty()) {
        return Err(bad("trailing corrupt data"));
    }

    let cp = Checkpoint {
        atoms,
        step: step.ok_or_else(|| bad("missing step"))?,
        time: time.ok_or_else(|| bad("missing time"))?,
        box_size,
        positions,
        velocities,
        continuation: None,
    };
    validate(&cp)?;
    Ok(cp)
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
    use std::fmt::Write as _;

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
            continuation: None,
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
