//! Sections: what the elements of a set are made of, how thick they are, and how they are
//! integrated.
//!
//! A section is assigned to one set and carries the element formulation for its elements, so one
//! region has one formulation and every solver can express it: the writers map the element kind
//! plus the [`Formulation`] to each format's element names, property cards or formulation
//! numbers.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, Model, Result, context, finite, insert_new, positive, remove_existing, replace_existing, require};

/// Most integration points through a shell's thickness.
pub const MAX_SHELL_POINTS: u32 = 15;

/// A beam's cross-section shape. Dimensions are in the section's local axes: `width` along the
/// beam's orientation vector (local 1), `height` across it (local 2).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum BeamProfile {
    Rectangle {
        width: f64,
        height: f64,
    },
    Circle {
        radius: f64,
    },
    /// A tube: outer radius and wall thickness.
    Pipe {
        radius: f64,
        thickness: f64,
    },
    /// A doubly symmetric I: overall height, flange width, web and flange thicknesses.
    I {
        height: f64,
        width: f64,
        web: f64,
        flange: f64,
    },
    /// A rectangular tube with one wall thickness.
    Box {
        width: f64,
        height: f64,
        thickness: f64,
    },
}

impl BeamProfile {
    pub fn validate(&self) -> Result<()> {
        match *self {
            BeamProfile::Rectangle { width, height } => {
                positive("width", width)?;
                positive("height", height)?;
            }
            BeamProfile::Circle { radius } => {
                positive("radius", radius)?;
            }
            BeamProfile::Pipe { radius, thickness } => {
                positive("radius", radius)?;
                positive("thickness", thickness)?;
                if thickness > radius {
                    return Err(Error::Invalid(format!("a pipe's wall ({thickness}) cannot be thicker than its radius ({radius})")));
                }
            }
            BeamProfile::I { height, width, web, flange } => {
                for (what, v) in [("height", height), ("width", width), ("web", web), ("flange", flange)] {
                    positive(what, v)?;
                }
                if web >= width || 2.0 * flange >= height {
                    return Err(Error::Invalid(
                        "an I profile's web must be narrower than its flanges, and its flanges thinner than half its height".into(),
                    ));
                }
            }
            BeamProfile::Box { width, height, thickness } => {
                for (what, v) in [("width", width), ("height", height), ("thickness", thickness)] {
                    positive(what, v)?;
                }
                if 2.0 * thickness >= width.min(height) {
                    return Err(Error::Invalid("a box's walls must be thinner than half its width and height".into()));
                }
            }
        }
        Ok(())
    }

    /// Cross-sectional area.
    pub fn area(&self) -> f64 {
        match *self {
            BeamProfile::Rectangle { width, height } => width * height,
            BeamProfile::Circle { radius } => std::f64::consts::PI * radius * radius,
            BeamProfile::Pipe { radius, thickness } => std::f64::consts::PI * (radius * radius - (radius - thickness).powi(2)),
            BeamProfile::I { height, width, web, flange } => 2.0 * width * flange + (height - 2.0 * flange) * web,
            BeamProfile::Box { width, height, thickness } => width * height - (width - 2.0 * thickness) * (height - 2.0 * thickness),
        }
    }
}

/// What a section describes.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum SectionKind {
    /// Continuum solid elements.
    Solid {},
    /// Shell elements: thickness, the reference surface's offset from the mid-surface as a
    /// fraction of the thickness (0 mid-surface, 0.5 on the side the normal points to, -0.5 on
    /// the other), and integration points through the thickness.
    Shell {
        thickness: f64,
        #[serde(default)]
        offset: f64,
        #[serde(default = "default_shell_points")]
        points: u32,
    },
    /// Beam elements with a library profile. `orientation` is a vector giving the profile's
    /// local 1 direction (projected normal to each beam's axis); it must not be parallel to a
    /// beam.
    Beam { profile: BeamProfile, orientation: [f64; 3] },
    /// Axial members (rods, trusses) with a cross-sectional area.
    Rod { area: f64 },
    /// 2D continuum in plane stress (σzz = 0: thin plates loaded in their plane), of the given
    /// thickness. The elements lie in the XY plane.
    PlaneStress { thickness: f64 },
    /// 2D continuum in plane strain (εzz = 0: a slice of a long body), of the given thickness
    /// (1 mm by default, so results are per millimetre). The elements lie in the XY plane.
    PlaneStrain {
        #[serde(default = "unit_thickness")]
        thickness: f64,
    },
    /// Axisymmetric continuum: the elements lie in the XY half plane x ≥ 0, which is the r–z
    /// plane of a body of revolution about the Y axis (x = radius, y = axial position; the hoop
    /// direction is the third). Forces, reactions and masses are totals over the full
    /// circumference (2π radians).
    Axisymmetric {},
}

fn default_shell_points() -> u32 {
    5
}

/// A plane-strain section's thickness when none is given: 1 mm, the default display system's
/// unit of length.
const UNIT_THICKNESS: crate::units::Typed = crate::units::Typed::new(1.0, "mm");

fn unit_thickness() -> f64 {
    // `1 mm` always converts; a NaN would be refused as a thickness.
    UNIT_THICKNESS.si(crate::units::Quantity::Length).unwrap_or(f64::NAN)
}

/// How a 2D continuum section extends out of its plane: what turns a length along an element
/// edge into an area, and an element's area into a volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Extent {
    /// Plane stress and plane strain: a uniform thickness.
    Thickness(f64),
    /// Axisymmetric: the full circumference 2πr at radius r = x.
    Revolved,
}

impl Extent {
    /// The out-of-plane measure at a point: the thickness, or 2π x.
    pub fn at(&self, xyz: [f64; 3]) -> f64 {
        match *self {
            Extent::Thickness(t) => t,
            Extent::Revolved => 2.0 * std::f64::consts::PI * xyz[0],
        }
    }
}

impl SectionKind {
    /// The element dimension this section applies to.
    pub fn dimension(&self) -> u8 {
        match self {
            SectionKind::Solid {} => 3,
            SectionKind::Shell { .. } | SectionKind::PlaneStress { .. } | SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric {} => 2,
            SectionKind::Beam { .. } | SectionKind::Rod { .. } => 1,
        }
    }

    /// Whether this is a continuum section (solid or 2D), the kinds with incompatible-mode
    /// formulations.
    pub fn is_continuum(&self) -> bool {
        self.dimension() == 3 || self.is_planar()
    }

    /// Whether this is a 2D continuum section (plane stress, plane strain or axisymmetric).
    pub fn is_planar(&self) -> bool {
        self.extent().is_some()
    }

    /// How a 2D continuum section extends out of its plane (`None` for other kinds).
    pub fn extent(&self) -> Option<Extent> {
        match *self {
            SectionKind::PlaneStress { thickness } | SectionKind::PlaneStrain { thickness } => Some(Extent::Thickness(thickness)),
            SectionKind::Axisymmetric {} => Some(Extent::Revolved),
            _ => None,
        }
    }

    /// `solid`, `shell`, `beam`, `rod`, `plane_stress`, `plane_strain` or `axisymmetric`.
    pub fn name(&self) -> &'static str {
        match self {
            SectionKind::Solid {} => "solid",
            SectionKind::Shell { .. } => "shell",
            SectionKind::Beam { .. } => "beam",
            SectionKind::Rod { .. } => "rod",
            SectionKind::PlaneStress { .. } => "plane_stress",
            SectionKind::PlaneStrain { .. } => "plane_strain",
            SectionKind::Axisymmetric {} => "axisymmetric",
        }
    }

    pub fn validate(&self) -> Result<()> {
        match self {
            SectionKind::Solid {} => {}
            SectionKind::Shell { thickness, offset, points } => {
                positive("thickness", *thickness)?;
                if !(finite("offset", *offset)?.abs() <= 0.5) {
                    return Err(Error::Invalid(format!("a shell offset is a fraction of the thickness in -0.5..=0.5, got {offset}")));
                }
                if !(1..=MAX_SHELL_POINTS).contains(points) {
                    return Err(Error::Invalid(format!(
                        "a shell has 1 to {MAX_SHELL_POINTS} integration points through its thickness, got {points}"
                    )));
                }
            }
            SectionKind::Beam { profile, orientation } => {
                profile.validate()?;
                for c in orientation {
                    finite("orientation", *c)?;
                }
                if orientation.iter().all(|c| *c == 0.0) {
                    return Err(Error::Invalid("a beam's orientation vector must not be zero".into()));
                }
            }
            SectionKind::Rod { area } => {
                positive("area", *area)?;
            }
            SectionKind::PlaneStress { thickness } | SectionKind::PlaneStrain { thickness } => {
                positive("thickness", *thickness)?;
            }
            SectionKind::Axisymmetric {} => {}
        }
        Ok(())
    }
}

/// Number of integration points per element direction.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Integration {
    /// The full rule for the element's order.
    #[default]
    Full,
    /// One order lower (cheaper, no shear locking, needs hourglass control).
    Reduced,
}

/// How the section's elements are formulated, independent of any solver's element names.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Formulation {
    #[serde(default)]
    pub integration: Integration,
    /// Enhanced (incompatible) bending modes: linear hexahedra (solid sections) and linear
    /// quadrilaterals (2D sections) only.
    #[serde(default)]
    pub incompatible_modes: bool,
    /// Mixed displacement–pressure formulation for nearly incompressible materials: solids only.
    #[serde(default)]
    pub hybrid: bool,
}

impl Formulation {
    /// Check the options against each other and against the section kind (the element kinds are
    /// checked when the model is resolved).
    pub fn validate_for(&self, kind: &SectionKind) -> Result<()> {
        if self.incompatible_modes && self.integration == Integration::Reduced {
            return Err(Error::Invalid("incompatible modes use full integration".into()));
        }
        if !kind.is_continuum() && self.incompatible_modes {
            return Err(Error::Invalid(format!("incompatible modes are for solid and 2D sections, not {}", kind.name())));
        }
        let volume_constrained = kind.dimension() == 3 || matches!(kind, SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric {});
        if self.hybrid && !volume_constrained {
            return Err(Error::Invalid(format!("hybrid formulations are for solid, plane strain and axisymmetric sections, not {}", kind.name())));
        }
        Ok(())
    }
}

/// A section: a material and a kind, assigned to the elements of a set.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Section {
    pub set: String,
    pub material: String,
    pub kind: SectionKind,
    #[serde(default)]
    pub formulation: Formulation,
}

impl Section {
    pub fn validate(&self) -> Result<()> {
        self.kind.validate()?;
        self.formulation.validate_for(&self.kind)
    }
}

impl Model {
    pub fn sections(&self) -> &BTreeMap<String, Section> {
        &self.sections
    }

    pub fn section(&self, name: &str) -> Option<&Section> {
        self.sections.get(name)
    }

    fn check_section(&self, name: &str, section: &Section) -> Result<()> {
        section.validate().map_err(|e| context("section", name, e))?;
        require(&self.sets, "set", &section.set).map_err(|e| context("section", name, e))?;
        require(&self.materials, "material", &section.material).map_err(|e| context("section", name, e))
    }

    pub fn add_section(&mut self, name: &str, section: Section) -> Result<()> {
        self.check_section(name, &section)?;
        insert_new(&mut self.sections, "section", name, section)
    }

    pub fn replace_section(&mut self, name: &str, section: Section) -> Result<()> {
        self.check_section(name, &section)?;
        replace_existing(&mut self.sections, "section", name, section)
    }

    pub fn remove_section(&mut self, name: &str) -> Result<()> {
        remove_existing(&mut self.sections, "section", name).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profiles_and_formulations_validate() {
        assert!(BeamProfile::Pipe { radius: 1.0, thickness: 1.5 }.validate().is_err());
        assert!(BeamProfile::I { height: 10.0, width: 5.0, web: 6.0, flange: 1.0 }.validate().is_err());
        assert!(BeamProfile::Box { width: 10.0, height: 4.0, thickness: 2.0 }.validate().is_err());
        assert!(BeamProfile::Rectangle { width: 1.0, height: f64::NAN }.validate().is_err());
        let i = BeamProfile::I { height: 10.0, width: 5.0, web: 1.0, flange: 1.0 };
        i.validate().unwrap();
        assert_eq!(i.area(), 2.0 * 5.0 + 8.0);
        assert!((BeamProfile::Pipe { radius: 2.0, thickness: 2.0 }.area() - BeamProfile::Circle { radius: 2.0 }.area()).abs() < 1e-12);

        let shell = SectionKind::Shell { thickness: 2.0, offset: 0.0, points: 5 };
        assert!(Formulation { hybrid: true, ..Formulation::default() }.validate_for(&shell).is_err());
        assert!(
            Formulation { incompatible_modes: true, integration: Integration::Reduced, hybrid: false }.validate_for(&SectionKind::Solid {}).is_err()
        );
        Formulation { incompatible_modes: true, hybrid: true, ..Formulation::default() }.validate_for(&SectionKind::Solid {}).unwrap();
        for bad in [
            SectionKind::Shell { thickness: 2.0, offset: 0.6, points: 5 },
            SectionKind::Shell { thickness: 2.0, offset: 0.0, points: 0 },
            SectionKind::Shell { thickness: 0.0, offset: 0.0, points: 5 },
            SectionKind::Beam { profile: BeamProfile::Circle { radius: 1.0 }, orientation: [0.0; 3] },
            SectionKind::Rod { area: -1.0 },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let parsed: SectionKind = serde_json::from_value(serde_json::json!({"type": "shell", "thickness": 1.5})).unwrap();
        assert_eq!(parsed, SectionKind::Shell { thickness: 1.5, offset: 0.0, points: 5 });
        assert!(
            serde_json::from_value::<SectionKind>(serde_json::json!({"type": "shell", "thickness": 1.5, "thicknes": 2})).is_err(),
            "typos are errors"
        );
        assert!(serde_json::from_value::<SectionKind>(serde_json::json!({"type": "solid", "x": 1})).is_err());
        assert!(serde_json::from_value::<SectionKind>(serde_json::json!({"type": "membrane"})).is_err());
    }

    #[test]
    fn planar_sections() {
        let strain: SectionKind = serde_json::from_value(serde_json::json!({"type": "plane_strain"})).unwrap();
        // 1 mm by default (results per millimetre, as before SI storage).
        assert_eq!(strain, SectionKind::PlaneStrain { thickness: 1e-3 });
        let stress: SectionKind = serde_json::from_value(serde_json::json!({"type": "plane_stress", "thickness": 2.5})).unwrap();
        assert_eq!(stress.extent(), Some(Extent::Thickness(2.5)));
        assert!(serde_json::from_value::<SectionKind>(serde_json::json!({"type": "plane_stress"})).is_err(), "plane stress needs a thickness");
        let axi: SectionKind = serde_json::from_value(serde_json::json!({"type": "axisymmetric"})).unwrap();
        assert_eq!(axi.extent(), Some(Extent::Revolved));
        assert!((Extent::Revolved.at([2.0, 7.0, 0.0]) - 4.0 * std::f64::consts::PI).abs() < 1e-15);
        assert!(axi.is_planar() && stress.is_planar() && !SectionKind::Solid {}.is_planar());
        assert!(!SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 }.is_planar());
        for bad in [SectionKind::PlaneStress { thickness: 0.0 }, SectionKind::PlaneStrain { thickness: f64::NAN }] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        let im = Formulation { incompatible_modes: true, ..Formulation::default() };
        im.validate_for(&stress).unwrap();
        im.validate_for(&axi).unwrap();
        let hybrid = Formulation { hybrid: true, ..Formulation::default() };
        hybrid.validate_for(&strain).unwrap();
        hybrid.validate_for(&axi).unwrap();
        assert!(hybrid.validate_for(&stress).is_err());
        assert!(im.validate_for(&SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 }).is_err());
    }
}
