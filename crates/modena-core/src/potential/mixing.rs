//! Per-type-pair parameter tables and combining (mixing) rules.
//!
//! A multi-species potential holds one entry per ordered-or-symmetric type pair.
//! The pair index for `(ti, tj)` is `ti * n_types + tj` (the table is stored
//! densely and symmetric: `(i,j)` and `(j,i)` share the same value, so callers
//! may use either ordering). Lorentz–Berthelot is the default combining rule
//! for Lennard-Jones; explicit per-pair overrides take precedence.

use crate::error::{ModenaError, Result};

/// How cross-species pair parameters are derived from single-species parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MixingRule {
    /// Lorentz–Berthelot: `σ_ij = (σ_i + σ_j)/2`, `ε_ij = √(ε_i · ε_j)`.
    LorentzBerthelot,
}

/// A dense symmetric table of one `f64` parameter per type pair, indexed
/// `ti * n_types + tj`. Built once at construction; read in the force hot loop.
#[derive(Debug, Clone)]
pub struct PairParamTable {
    /// `n_types * n_types` entries.
    pub values: Vec<f64>,
    pub n_types: usize,
}

impl PairParamTable {
    pub fn new(n_types: usize) -> Self {
        Self {
            values: vec![0.0; n_types * n_types],
            n_types,
        }
    }

    /// Set the symmetric pair `(i, j)` (and `(j, i)`).
    pub fn set_symmetric(&mut self, i: usize, j: usize, v: f64) {
        self.values[i * self.n_types + j] = v;
        self.values[j * self.n_types + i] = v;
    }

    /// Set the diagonal `(i, i)`.
    pub fn set_diagonal(&mut self, i: usize, v: f64) {
        self.values[i * self.n_types + i] = v;
    }

    /// Get the pair value for `(ti, tj)` (symmetric; order-independent).
    #[inline]
    pub fn get(&self, ti: u32, tj: u32) -> f64 {
        // Order-independent lookup: the table is symmetric, so ti/tj order does
        // not matter, but we keep the canonical index for cache friendliness.
        let idx = ti as usize * self.n_types + tj as usize;
        // SAFETY: callers guarantee ti, tj < n_types (the potential validated it
        // at construction). The force kernel never sees out-of-range species.
        self.values[idx]
    }

    /// Fill a table from per-type parameters using a combining rule.
    pub fn from_types(
        per_type: &[f64],
        rule: MixingRule,
        overrides: &[(usize, usize, f64)],
    ) -> Result<Self> {
        let n = per_type.len();
        if n == 0 {
            return Err(ModenaError::Config(
                "at least one species type is required".into(),
            ));
        }
        let mut table = Self::new(n);
        match rule {
            MixingRule::LorentzBerthelot => {
                for i in 0..n {
                    for j in 0..n {
                        // Geometric mean — the LB rule for an energy-like
                        // parameter (ε). Length-like parameters (σ) use the
                        // arithmetic mean; callers build separate tables.
                        let v = (per_type[i] * per_type[j]).sqrt();
                        table.set_symmetric(i, j, v);
                    }
                }
            }
        }
        for &(i, j, v) in overrides {
            if i >= n || j >= n {
                return Err(ModenaError::Config(format!(
                    "override type index {i}/{j} out of range (n_types={n})"
                )));
            }
            table.set_symmetric(i, j, v);
        }
        Ok(table)
    }
}

/// Lorentz–Berthelot combining rules as free functions (used during table
/// construction by potentials that store σ and ε in separate tables).
pub mod lorentz_berthelot {
    /// σ_ij = (σ_i + σ_j) / 2  (arithmetic mean — the LB length rule).
    #[inline]
    pub fn sigma(sigma_i: f64, sigma_j: f64) -> f64 {
        0.5 * (sigma_i + sigma_j)
    }

    /// ε_ij = √(ε_i · ε_j)  (geometric mean — the LB energy rule).
    #[inline]
    pub fn epsilon(epsilon_i: f64, epsilon_j: f64) -> f64 {
        (epsilon_i * epsilon_j).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lb_table_is_symmetric() {
        let table = PairParamTable::from_types(&[1.0, 4.0, 9.0], MixingRule::LorentzBerthelot, &[])
            .unwrap();
        // Geometric mean √(ε_i·ε_j): the diagonal (i,i) = √(ε_i²) = ε_i.
        assert!((table.get(0, 0) - 1.0).abs() < 1e-12);
        assert!((table.get(1, 1) - 4.0).abs() < 1e-12);
        assert!((table.get(2, 2) - 9.0).abs() < 1e-12);
        // Off-diagonal: √(1·4)=2, √(1·9)=3, √(4·9)=6, symmetric.
        assert!((table.get(0, 1) - 2.0).abs() < 1e-12);
        assert!((table.get(1, 0) - 2.0).abs() < 1e-12);
        assert!((table.get(0, 2) - 3.0).abs() < 1e-12);
        assert!((table.get(2, 0) - 3.0).abs() < 1e-12);
        assert!((table.get(1, 2) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn override_replaces_mixing_default() {
        let table =
            PairParamTable::from_types(&[1.0, 4.0], MixingRule::LorentzBerthelot, &[(0, 1, 99.0)])
                .unwrap();
        assert!((table.get(0, 1) - 99.0).abs() < 1e-12);
        assert!((table.get(1, 0) - 99.0).abs() < 1e-12);
        // Unmodified diagonal still uses mixing.
        assert!((table.get(0, 0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn out_of_range_override_errors() {
        let err = PairParamTable::from_types(&[1.0], MixingRule::LorentzBerthelot, &[(0, 5, 1.0)]);
        assert!(err.is_err());
    }

    #[test]
    fn empty_types_errors() {
        assert!(PairParamTable::from_types(&[], MixingRule::LorentzBerthelot, &[]).is_err());
    }
}
