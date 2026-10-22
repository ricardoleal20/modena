//! XYZ structure/trajectory I/O.

use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::error::{ModenaError, Result};
use crate::system::{SimBox, System};
use crate::vec3::Vec3;

/// A parsed XYZ frame: element symbols + positions.
pub struct XyzFrame {
    pub elements: Vec<String>,
    pub positions: Vec<Vec3>,
    pub comment: String,
}

/// Read a single-frame XYZ file into a `System`. The simulation box must be
/// supplied by the caller (XYZ carries no box info).
pub fn read_xyz<P: AsRef<Path>>(path: P, box_size: SimBox) -> Result<System> {
    Ok(read_xyz_with_identity(path, box_size)?.0)
}

/// Typed XYZ uses exactly `Properties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1`.
/// Type indices are zero-based model indices; mass/charge/positions use input units.
/// Plain XYZ retains legacy mass=1, charge=0, type=0 behavior.
pub fn read_xyz_with_identity<P: AsRef<Path>>(path: P, box_size: SimBox) -> Result<(System, bool)> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();

    reader.read_line(&mut line)?;
    let n: usize = line
        .trim()
        .parse()
        .map_err(|e| ModenaError::Input(format!("XYZ atom count not an integer: {e}")))?;

    if n > u32::MAX as usize || n.checked_mul(3).is_none() {
        return Err(ModenaError::Input(
            "XYZ atom count exceeds supported indexed atom capacity".into(),
        ));
    }
    // Plain comment is ignored; an explicit Properties declaration is strict.
    line.clear();
    reader.read_line(&mut line)?;

    let typed = line
        .split_whitespace()
        .find(|s| s.starts_with("Properties="));
    let typed = match typed {
        None => false,
        Some("Properties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1") => true,
        Some(_) => return Err(ModenaError::Input("unsupported XYZ Properties declaration; expected species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1".into())),
    };
    let mut identities = Vec::new();
    let mut elements = Vec::new();
    let mut positions = Vec::new();
    for _ in 0..n {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(ModenaError::AtomCount {
                expected: n,
                found: positions.len(),
            });
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 || (typed && parts.len() != 7) {
            return Err(ModenaError::Input(format!("malformed XYZ line: {line}")));
        }
        let el = parts[0].to_string();
        let x: f64 = parts[1]
            .parse()
            .map_err(|e| ModenaError::Input(format!("bad x coordinate: {e}")))?;
        let y: f64 = parts[2]
            .parse()
            .map_err(|e| ModenaError::Input(format!("bad y coordinate: {e}")))?;
        let z: f64 = parts[3]
            .parse()
            .map_err(|e| ModenaError::Input(format!("bad z coordinate: {e}")))?;
        if ![x, y, z].iter().all(|v| v.is_finite()) {
            return Err(ModenaError::Input("XYZ coordinates must be finite".into()));
        }
        if typed {
            let species: u32 = parts[4].parse().map_err(|_| {
                ModenaError::Input("XYZ type must be a zero-based integer model index".into())
            })?;
            let mass: f64 = parts[5]
                .parse()
                .map_err(|_| ModenaError::Input("bad XYZ mass".into()))?;
            let charge: f64 = parts[6]
                .parse()
                .map_err(|_| ModenaError::Input("bad XYZ charge".into()))?;
            if !mass.is_finite() || mass <= 0.0 || !charge.is_finite() {
                return Err(ModenaError::Input(
                    "typed XYZ requires positive finite masses and finite charges".into(),
                ));
            }
            identities.push((species, mass, charge));
        }
        elements.push(el);
        positions.push(Vec3::new(x, y, z));
    }

    let mut system = System::new(box_size);
    system.resize(n);
    // Flatten the per-atom Vec3 positions into the SoA buffer.
    for (i, p) in positions.iter().enumerate() {
        system.set_pos(i, *p);
    }
    system.elements = elements;
    for (i, (species, mass, charge)) in identities.into_iter().enumerate() {
        system.species[i] = species;
        system.masses[i] = mass;
        system.charges[i] = charge;
    }
    Ok((system, typed))
}

/// Streaming XYZ trajectory writer (append mode).
pub struct TrajectoryWriter {
    writer: BufWriter<File>,
}

impl TrajectoryWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        if let Some(parent) = path.as_ref().parent() {
            // Decks point trajectory_file at results/…; create it on first use.
            std::fs::create_dir_all(parent)?;
        }
        let file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(path)?;
        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    pub fn write_frame(
        &mut self,
        system: &System,
        step: usize,
        time: f64,
        comment: &str,
    ) -> Result<()> {
        self.write_frame_scaled(system, step, time, comment, 1.0)
    }

    /// Like [`write_frame`](Self::write_frame) but multiplies each position
    /// component by `length_scale`. Used to render real-unit (Å) coordinates in
    /// the trajectory when the core runs in reduced units (`length_scale = σ`).
    pub fn write_frame_scaled(
        &mut self,
        system: &System,
        step: usize,
        time: f64,
        comment: &str,
        length_scale: f64,
    ) -> Result<()> {
        let n = system.n();
        writeln!(self.writer, "{n}")?;
        writeln!(self.writer, "step={step} time={time:.6} {comment}")?;
        for i in 0..n {
            let el = system.elements.get(i).map(|s| s.as_str()).unwrap_or("X");
            let p = system.pos(i);
            writeln!(
                self.writer,
                "{el} {:.8} {:.8} {:.8}",
                p.0[0] * length_scale,
                p.0[1] * length_scale,
                p.0[2] * length_scale,
            )?;
        }
        self.writer.flush()?;
        Ok(())
    }
}
