//! Named scalar design parameters. Values are SI, and updates retain the model's
//! normal validation and replace atomically. No solver-specific state is stored.
use crate::units::Quantity;
use crate::{Elastic, Error, Model, Result, SectionKind, Spring};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaterialConstant {
    Young,
    Poisson,
    Density,
    E1,
    E2,
    E3,
    Nu12,
    Nu13,
    Nu23,
    G12,
    G13,
    G23,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Parameter {
    Material { name: String, constant: MaterialConstant },
    ShellThickness { name: String },
    SpringStiffness { name: String },
}

fn missing() -> Error {
    Error::Invalid("design parameter does not name an existing compatible property".into())
}

fn material_slot(material: &mut crate::Material, constant: MaterialConstant) -> Result<&mut f64> {
    use MaterialConstant as C;
    if constant == C::Density {
        return material.density.as_mut().ok_or_else(missing);
    }
    match (material.elastic.as_mut(), constant) {
        (Some(Elastic::Isotropic { young, .. }), C::Young) => Ok(young),
        (Some(Elastic::Isotropic { poisson, .. }), C::Poisson) => Ok(poisson),
        (Some(Elastic::EngineeringConstants { young, poisson, shear }), c) => {
            let (array, i) = match c {
                C::E1 => (young, 0),
                C::E2 => (young, 1),
                C::E3 => (young, 2),
                C::Nu12 => (poisson, 0),
                C::Nu13 => (poisson, 1),
                C::Nu23 => (poisson, 2),
                C::G12 => (shear, 0),
                C::G13 => (shear, 1),
                C::G23 => (shear, 2),
                _ => return Err(missing()),
            };
            array.get_mut(i).ok_or_else(missing)
        }
        _ => Err(missing()),
    }
}

impl Parameter {
    pub fn value(&self, model: &Model) -> Result<f64> {
        match self {
            Self::Material { name, constant } => {
                let mut material = model.material(name).cloned().ok_or_else(missing)?;
                Ok(*material_slot(&mut material, *constant)?)
            }
            Self::ShellThickness { name } => match model.section(name).map(|s| &s.kind) {
                Some(SectionKind::Shell { thickness, .. }) => Ok(*thickness),
                _ => Err(missing()),
            },
            Self::SpringStiffness { name } => match model.spring(name) {
                Some(Spring::Spring { stiffness, .. }) => Ok(*stiffness),
                _ => Err(missing()),
            },
        }
    }

    pub fn quantity(&self, model: &Model) -> Result<Quantity> {
        self.value(model)?;
        Ok(match self {
            Self::Material { constant, .. } => match constant {
                MaterialConstant::Density => Quantity::Density,
                MaterialConstant::Poisson | MaterialConstant::Nu12 | MaterialConstant::Nu13 | MaterialConstant::Nu23 => Quantity::Dimensionless,
                _ => Quantity::Pressure,
            },
            Self::ShellThickness { .. } => Quantity::Length,
            Self::SpringStiffness { name } => match model.spring(name).map(|s| s.direction()) {
                Some(crate::discrete::Direction::Dof(d)) if d > 3 => Quantity::RotationalStiffness,
                _ => Quantity::Stiffness,
            },
        })
    }

    pub fn set(&self, model: &mut Model, value: f64) -> Result<()> {
        if !value.is_finite() {
            return Err(Error::Invalid("design parameter must be finite".into()));
        }
        match self {
            Self::Material { name, constant } => {
                let mut material = model.material(name).cloned().ok_or_else(missing)?;
                *material_slot(&mut material, *constant)? = value;
                model.replace_material(name, material)
            }
            Self::ShellThickness { name } => {
                let mut section = model.section(name).cloned().ok_or_else(missing)?;
                let SectionKind::Shell { thickness, .. } = &mut section.kind else {
                    return Err(missing());
                };
                *thickness = value;
                model.replace_section(name, section)
            }
            Self::SpringStiffness { name } => {
                let mut spring = model.spring(name).cloned().ok_or_else(missing)?;
                let Spring::Spring { stiffness, .. } = &mut spring else {
                    return Err(missing());
                };
                *stiffness = value;
                model.replace_spring(name, spring)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_orthotropic_constants_update_only_the_named_slot() {
        let mut m = Model::default();
        m.add_material(
            "ortho",
            crate::Material {
                elastic: Some(Elastic::EngineeringConstants { young: [100.0, 80.0, 60.0], poisson: [0.1, 0.15, 0.2], shear: [40.0, 30.0, 20.0] }),
                density: Some(2.0),
                ..Default::default()
            },
        )
        .unwrap();
        for c in [
            MaterialConstant::E1,
            MaterialConstant::E2,
            MaterialConstant::E3,
            MaterialConstant::Nu12,
            MaterialConstant::Nu13,
            MaterialConstant::Nu23,
            MaterialConstant::G12,
            MaterialConstant::G13,
            MaterialConstant::G23,
            MaterialConstant::Density,
        ] {
            let p = Parameter::Material { name: "ortho".into(), constant: c };
            let before = m.clone();
            let v = p.value(&m).unwrap();
            p.set(&mut m, v * 1.01).unwrap();
            assert_eq!(p.value(&m).unwrap(), v * 1.01);
            p.set(&mut m, v).unwrap();
            assert_eq!(m, before);
        }
    }
    #[test]
    fn invalid_material_update_does_not_change_model() {
        let mut m = Model::default();
        m.add_material(
            "iso",
            crate::Material { elastic: Some(Elastic::Isotropic { young: 100.0, poisson: 0.3 }), density: Some(1.0), ..Default::default() },
        )
        .unwrap();
        for (c, v) in [(MaterialConstant::Young, -1.0), (MaterialConstant::Poisson, 0.5), (MaterialConstant::Density, f64::NAN)] {
            let before = m.clone();
            assert!(Parameter::Material { name: "iso".into(), constant: c }.set(&mut m, v).is_err());
            assert_eq!(m, before);
        }
    }
}
