//! Tabulated pair potential: `(r, U, F)` sampled on a grid, interpolated.
//!
//! For potentials without a convenient closed form (e.g. derived from ab-initio
//! data or a coarse-grained map). The table supplies both `U(r)` and the radial
//! force `F(r)` directly, so the conservation contract `dU/dr = −F` is whatever
//! the table says — the caller is responsible for a self-consistent table.
//! Linear interpolation between nodes; shifted-force at the table's last point
//! so both `U` and `F` are continuous at the cutoff.

#![allow(clippy::neg_cmp_op_on_partial_ord)]

use crate::error::{ModenaError, Result};
use crate::potential::check_cutoff;

use super::shifted_force;
use super::PairPotential;

/// Interpolation scheme between table nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InterpKind {
    /// Piecewise-linear between nodes.
    Linear,
}

/// A tabulated pair potential. The grid must be strictly increasing in `r`,
/// start at `r > 0`, and the last node is the cutoff.
#[derive(Debug, Clone)]
pub struct Tabulated {
    /// Strictly increasing radial nodes.
    r: Vec<f64>,
    /// Energy at each node.
    u: Vec<f64>,
    /// Radial force magnitude `F` at each node (positive = repulsive).
    f: Vec<f64>,
    interp: InterpKind,
    cutoff: f64,
    cutoff2: f64,
}

impl Tabulated {
    /// Build from explicit node arrays. `r` must be strictly increasing and
    /// positive; the last `r` is the cutoff.
    pub fn new(r: Vec<f64>, u: Vec<f64>, f: Vec<f64>, interp: InterpKind) -> Result<Self> {
        if r.len() < 2 {
            return Err(ModenaError::Config(
                "tabulated potential needs at least 2 nodes".into(),
            ));
        }
        if r.len() != u.len() || r.len() != f.len() {
            return Err(ModenaError::Config(
                "tabulated r/U/F arrays must have equal length".into(),
            ));
        }
        for win in r.windows(2) {
            if !(win[1] > win[0]) {
                return Err(ModenaError::Config(
                    "tabulated r nodes must be strictly increasing".into(),
                ));
            }
        }
        if !(r[0] > 0.0) {
            return Err(ModenaError::Config("tabulated r nodes must be > 0".into()));
        }
        let cutoff = *r.last().unwrap();
        check_cutoff(cutoff, "tabulated")?;
        Ok(Self {
            r,
            u,
            f,
            interp,
            cutoff,
            cutoff2: cutoff * cutoff,
        })
    }

    /// Parse a whitespace-delimited table: one `r U F` triple per line, `#`
    /// comments and blank lines ignored. Nodes must be in increasing `r` order.
    pub fn from_text(text: &str, interp: InterpKind) -> Result<Self> {
        let mut r = Vec::new();
        let mut u = Vec::new();
        let mut f = Vec::new();
        for (lineno, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 {
                return Err(ModenaError::Input(format!(
                    "tabulated line {}: expected `r U F`, got {parts:?}",
                    lineno + 1
                )));
            }
            let parse = |s: &str, name: &str| -> Result<f64> {
                s.parse::<f64>()
                    .map_err(|e| ModenaError::Input(format!("tabulated {name} parse: {e}")))
            };
            r.push(parse(parts[0], "r")?);
            u.push(parse(parts[1], "U")?);
            f.push(parse(parts[2], "F")?);
        }
        Self::new(r, u, f, interp)
    }

    /// Linear interpolation of `(U, F)` at distance `r`. Returns the raw
    /// (unshifted) values; the caller applies the shifted-force form.
    #[inline]
    fn interp(&self, r: f64) -> (f64, f64) {
        // Clamp to the table range (the caller already gates r2 < cutoff2, so
        // r <= cutoff; values below r[0] are extrapolated flat to the first node).
        match self.interp {
            InterpKind::Linear => {
                if r <= self.r[0] {
                    return (self.u[0], self.f[0]);
                }
                // Binary search for the bracketing interval.
                let idx = match self
                    .r
                    .binary_search_by(|v| v.partial_cmp(&r).unwrap_or(std::cmp::Ordering::Equal))
                {
                    Ok(i) => return (self.u[i], self.f[i]),
                    Err(i) => i,
                };
                // idx is the insertion point; bracket is [idx-1, idx].
                let i0 = idx - 1;
                let i1 = idx.min(self.r.len() - 1);
                let r0 = self.r[i0];
                let r1 = self.r[i1];
                let t = if r1 > r0 { (r - r0) / (r1 - r0) } else { 0.0 };
                let u = self.u[i0] + t * (self.u[i1] - self.u[i0]);
                let f = self.f[i0] + t * (self.f[i1] - self.f[i0]);
                (u, f)
            }
        }
    }
}

impl PairPotential for Tabulated {
    #[inline]
    fn cutoff(&self) -> f64 {
        self.cutoff
    }

    #[inline]
    fn cutoff2(&self) -> f64 {
        self.cutoff2
    }

    #[inline]
    fn energy_force(&self, r2: f64, _ti: u32, _tj: u32) -> (f64, f64) {
        let r = r2.sqrt();
        let (u_raw, f_raw) = self.interp(r);
        let (u, f_sf) = shifted_force::shift(
            u_raw,
            f_raw,
            r,
            self.cutoff,
            *self.u.last().unwrap(),
            *self.f.last().unwrap(),
        );
        (u, f_sf / r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> Tabulated {
        // A small monotonic table.
        let text = "
        # r      U       F
        0.9     10.0    50.0
        1.0      5.0    20.0
        1.5      0.5     2.0
        2.0      0.05    0.2
        2.5      0.0     0.0
        ";
        Tabulated::from_text(text, InterpKind::Linear).unwrap()
    }

    #[test]
    fn node_values_match_exactly() {
        let t = table();
        for (i, &r) in t.r.iter().enumerate() {
            let (u, f_over_r) = t.energy_force(r * r, 0, 0);
            // Shifted-force: at the last node U=F=0; at interior nodes the shift
            // subtracts U(rc) and F(rc) (both 0 here), so raw == shifted.
            let f_radial = f_over_r * r;
            assert!((u - t.u[i]).abs() < 1e-9, "U at node {i}");
            assert!((f_radial - t.f[i]).abs() < 1e-9, "F at node {i}");
        }
    }

    #[test]
    fn energy_and_force_zero_at_cutoff() {
        let t = table();
        let (u, f_over_r) = t.energy_force(t.cutoff2, 0, 0);
        assert!(u.abs() < 1e-9, "shifted U at cutoff = {u}");
        assert!((f_over_r * t.cutoff()).abs() < 1e-9, "shifted F at cutoff");
    }

    #[test]
    fn linear_interp_midpoint() {
        let t = table();
        // Midway between r=1.0 (U=5) and r=1.5 (U=0.5): U=2.75.
        let r = 1.25;
        let (u, _) = t.interp(r);
        assert!((u - 2.75).abs() < 1e-9);
    }

    #[test]
    fn rejects_non_increasing_nodes() {
        let err = Tabulated::new(
            vec![1.0, 0.9, 2.0],
            vec![1.0, 2.0, 3.0],
            vec![1.0, 2.0, 3.0],
            InterpKind::Linear,
        );
        assert!(err.is_err());
    }
}
