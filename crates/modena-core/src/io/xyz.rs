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
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();

    reader.read_line(&mut line)?;
    let n: usize = line
        .trim()
        .parse()
        .map_err(|e| ModenaError::Input(format!("XYZ atom count not an integer: {e}")))?;

    // Comment line (ignored on read).
    line.clear();
    reader.read_line(&mut line)?;

    let mut elements = Vec::with_capacity(n);
    let mut positions = Vec::with_capacity(n);
    for _ in 0..n {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(ModenaError::AtomCount {
                expected: n,
                found: positions.len(),
            });
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
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
    Ok(system)
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
