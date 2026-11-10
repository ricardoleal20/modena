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

/// Scientific domain status. Reading a claimed record never authenticates it.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    content = "value",
    rename_all = "snake_case",
    deny_unknown_fields
)]
pub enum ValidationStatus {
    #[default]
    Diagnostic,
    Validated {
        pub_record: String,
    },
}

/// Mechanical endpoint observation. Kinetic, model and constraint tensors are
/// energy numerators; total and scalar are in reduced pressure units.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PressureObservation {
    pub kinetic: SymTensor6,
    pub model: TensorAvailability,
    pub constraint: SymTensor6,
    pub total: TensorAvailability,
    pub scalar: Option<f64>,
    pub volume: f64,
    pub validation: ValidationStatus,
}

impl<'de> Deserialize<'de> for PressureObservation {
    fn deserialize<D: Deserializer<'de>>(d: D) -> std::result::Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Record {
            kinetic: SymTensor6,
            model: TensorAvailability,
            constraint: SymTensor6,
            total: TensorAvailability,
            scalar: Option<f64>,
            volume: f64,
            #[serde(rename = "validation")]
            _claimed_status: ValidationStatus,
        }
        let record = Record::deserialize(d)?;
        if !record.volume.is_finite()
            || record.volume <= 0.
            || record.scalar.is_some_and(|p| !p.is_finite())
        {
            return Err(serde::de::Error::custom(
                "invalid endpoint volume or scalar pressure",
            ));
        }
        match (&record.model, &record.total, record.scalar) {
            (TensorAvailability::Missing(model), TensorAvailability::Missing(total), None)
                if model == total => {}
            (TensorAvailability::Complete(_), TensorAvailability::Complete(_), Some(_)) => (),
            _ => {
                return Err(serde::de::Error::custom(
                    "inconsistent endpoint tensor availability or scalar",
                ))
            }
        }
        Ok(Self {
            kinetic: record.kinetic,
            model: record.model,
            constraint: record.constraint,
            total: record.total,
            scalar: record.scalar,
            volume: record.volume,
            validation: ValidationStatus::Diagnostic,
        })
    }
}

/// Mechanical kinetic numerator Σ mᵢ vᵢ⊗vᵢ (twice kinetic energy in trace).
pub fn kinetic_tensor(system: &crate::System) -> SymTensor6 {
    (0..system.n()).fold(SymTensor6::ZERO, |sum, i| {
        let v = system.vel(i);
        sum.plus(SymTensor6::dyad(v, v).scaled(system.masses[i]))
    })
}

/// Assemble pressure at one snapshot, using its cached conservative forces.
pub fn endpoint_pressure(
    system: &crate::System,
    model: &ModelEvaluation,
    constraints: Option<&crate::Constraints>,
) -> crate::Result<PressureObservation> {
    use crate::error::ModenaError;
    let volume = system.box_size.volume();
    if !volume.is_finite() || volume <= 0. {
        return Err(ModenaError::Physics(format!(
            "invalid endpoint volume {volume:e}"
        )));
    }
    if !model.energy.is_finite() || !model.scalar_virial.is_finite() {
        return Err(ModenaError::Physics(
            "nonfinite endpoint model evaluation".into(),
        ));
    }
    let empty = crate::Constraints::new();
    let constraint =
        crate::constraint_stress::reaction_from_constraints(system, constraints.unwrap_or(&empty))?
            .tensor;
    let kinetic = kinetic_tensor(system);
    if !kinetic.finite() {
        return Err(ModenaError::Physics(
            "nonfinite endpoint kinetic tensor".into(),
        ));
    }
    let (total, scalar) = match &model.tensor {
        TensorAvailability::Missing(why) => (TensorAvailability::Missing(why.clone()), None),
        TensorAvailability::Complete(tensor) => {
            let numerator = kinetic.plus(*tensor).plus(constraint);
            let pressure = SymTensor6(numerator.0.map(|x| x / volume));
            let scalar = pressure.trace() / 3.;
            if !pressure.finite() || !scalar.is_finite() {
                return Err(ModenaError::Physics("nonfinite endpoint pressure".into()));
            }
            (TensorAvailability::Complete(pressure), Some(scalar))
        }
    };
    Ok(PressureObservation {
        kinetic,
        model: model.tensor.clone(),
        constraint,
        total,
        scalar,
        volume,
        validation: ValidationStatus::Diagnostic,
    })
}
