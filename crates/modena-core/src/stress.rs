//! Symmetric configurational virials in energy units, ordered xx, yy, zz, xy, xz, yz.
//!
//! Positive diagonal components are compressive. The scalar virial is the
//! trace; unavailable tensor data must never be inferred from that scalar.

use serde::{Deserialize, Deserializer, Serialize};

use crate::vec3::Vec3;

/// Symmetric tensor with six independent components in energy units.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "[f64; 6]")]
pub struct SymTensor6(pub [f64; 6]);

impl SymTensor6 {
    pub const ZERO: Self = Self([0.0; 6]);

    /// Symmetric dyad `(a ⊗ b + b ⊗ a) / 2`.
    #[inline]
    pub fn dyad(a: Vec3, b: Vec3) -> Self {
        let [ax, ay, az] = a.0;
        let [bx, by, bz] = b.0;
        Self([
            ax * bx,
            ay * by,
            az * bz,
            0.5 * (ax * by + ay * bx),
            0.5 * (ax * bz + az * bx),
            0.5 * (ay * bz + az * by),
        ])
    }

    pub fn trace(self) -> f64 {
        self.0[0] + self.0[1] + self.0[2]
    }

    pub fn scaled(self, x: f64) -> Self {
        Self(self.0.map(|v| v * x))
    }

    #[inline]
    pub fn plus(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] + rhs.0[i]))
    }

    pub fn finite(self) -> bool {
        self.0.iter().all(|x| x.is_finite())
    }
}

impl TryFrom<[f64; 6]> for SymTensor6 {
    type Error = &'static str;

    fn try_from(value: [f64; 6]) -> std::result::Result<Self, Self::Error> {
        let tensor = Self(value);
        if tensor.finite() {
            Ok(tensor)
        } else {
            Err("tensor components must be finite")
        }
    }
}

/// Why a provider cannot supply a complete configurational tensor.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    content = "value",
    rename_all = "snake_case",
    deny_unknown_fields
)]
pub enum MissingTensorReason {
    LegacyForceModel,
    UnsupportedBackend,
    MissingInitialEvaluation,
    PostStepTransform,
    LegacyRestart,
    ProviderPropertyUnavailable,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MissingTensor {
    pub component: String,
    pub reason: MissingTensorReason,
}

/// Missing contributions are explicit; partial tensor sums are not complete.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    content = "value",
    rename_all = "snake_case",
    deny_unknown_fields
)]
pub enum TensorAvailability {
    Complete(SymTensor6),
    Missing(Vec<MissingTensor>),
}

/// One provider's additive energy, scalar virial, and tensor contribution.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelEvaluation {
    #[serde(deserialize_with = "finite_scalar")]
    pub energy: f64,
    #[serde(deserialize_with = "finite_scalar")]
    pub scalar_virial: f64,
    pub tensor: TensorAvailability,
}

impl ModelEvaluation {
    /// A known empty contribution, with a complete zero tensor.
    pub fn zero() -> Self {
        Self {
            energy: 0.0,
            scalar_virial: 0.0,
            tensor: TensorAvailability::Complete(SymTensor6::ZERO),
        }
    }

    /// Sum scalar contributions and tensors, retaining every missing reason.
    pub fn combine(self, other: Self) -> Self {
        let tensor = match (self.tensor, other.tensor) {
            (TensorAvailability::Complete(a), TensorAvailability::Complete(b)) => {
                TensorAvailability::Complete(a.plus(b))
            }
            (TensorAvailability::Missing(mut a), TensorAvailability::Missing(b)) => {
                a.extend(b);
                TensorAvailability::Missing(a)
            }
            (TensorAvailability::Missing(a), _) | (_, TensorAvailability::Missing(a)) => {
                TensorAvailability::Missing(a)
            }
        };
        Self {
            energy: self.energy + other.energy,
            scalar_virial: self.scalar_virial + other.scalar_virial,
            tensor,
        }
    }
}

fn finite_scalar<'de, D: Deserializer<'de>>(deserializer: D) -> Result<f64, D::Error> {
    let value = f64::deserialize(deserializer)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(serde::de::Error::custom(
            "evaluation scalars must be finite",
        ))
    }
}
