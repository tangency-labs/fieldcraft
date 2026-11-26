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

/// A beam's cross-section: a shape, or its section properties given directly (`general`).
/// Dimensions are in the section's local axes: `width` along the beam's orientation vector
/// (local 1), `height` across it (local 2).
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
    /// Any section by its properties in the local axes, as section tables give them: area,
    /// second moments `iy` = ∫ z² dA (about local 1) and `iz` = ∫ y² dA (about local 2), the
    /// Saint-Venant torsion constant and the shear areas for shear along local 1 and local 2
    /// (the meanings of [`BeamProperties`]). The local axes must be the principal axes (no
    /// product of inertia) through the centroid, which is also the shear centre.
    General {
        area: f64,
        iy: f64,
        iz: f64,
        torsion: f64,
        shear_y: f64,
        shear_z: f64,
    },
}

/// Roundings in the largest bound expression of [`general_bounds`]: the disc's A²/2π takes two
/// (A·A, then the division by 2π) and the polar moment I_y + I_z one more (the sum).
const BOUND_ROUNDINGS: f64 = 3.0;

/// The relative slack of [`general_bounds`]: rounding in the arithmetic only. A disc meets two of
/// the bounds with equality (J = I_y + I_z = A²/2π) and a tube one (J = I_y + I_z), and their
/// properties computed in floating point, or converted to SI, may miss the bound by rounding.
/// Each rounding is a relative error of at most `f64::EPSILON` (twice the unit roundoff), so
/// the slack is `BOUND_ROUNDINGS` of them.
const GENERAL_BOUND_SLACK: f64 = BOUND_ROUNDINGS * f64::EPSILON;

/// Bounds every cross-section's properties obey, refused by name when a general profile breaks
/// them:
/// - a shear area is at most the area: with A_s = V² / ∫ τ² dA, the energy-equivalent shear
///   area, and ∫ τ dA = V, Cauchy–Schwarz gives ∫ τ² dA ≥ V² / A, so the shear coefficient
///   A_s / A is at most 1 (Cowper 1966, "The shear coefficient in Timoshenko's beam theory",
///   J. Appl. Mech. 33, finds it below 1 for every shape he treats);
/// - the torsion constant is at most the polar moment I_y + I_z: Saint-Venant's J is I_p less
///   the integral of the squared warping gradient, so equal only for a section that does not
///   warp, a circle or a tube (Timoshenko & Goodier 1970, Theory of Elasticity, ch. 10);
/// - the polar moment I_y + I_z is at least A²/2π, a disc's of the same area: moving any part of
///   the area to a radius nearer the centroid lowers ∫ r² dA, and the disc has all of it nearest
///   (Pólya & Szegő 1951, Isoperimetric Inequalities in Mathematical Physics).
fn general_bounds(area: f64, iy: f64, iz: f64, torsion: f64, shear_y: f64, shear_z: f64) -> Result<()> {
    let above = |v: f64, bound: f64| v > bound * (1.0 + GENERAL_BOUND_SLACK);
    for (what, along, shear) in [("shear_y", 1, shear_y), ("shear_z", 2, shear_z)] {
        if above(shear, area) {
            return Err(Error::Invalid(format!(
                "a general profile's shear area along local {along} ({what} = {shear}) exceeds its area ({area}): no section has a shear coefficient A_s/A above 1"
            )));
        }
    }
    let polar = iy + iz;
    if above(torsion, polar) {
        return Err(Error::Invalid(format!(
            "a general profile's torsion constant ({torsion}) exceeds its polar moment iy + iz ({polar}): no section is stiffer in torsion than a circle or a tube of the same polar moment"
        )));
    }
    let disc = area * area / (2.0 * std::f64::consts::PI);
    if above(disc, polar) {
        return Err(Error::Invalid(format!(
            "a general profile's polar moment iy + iz ({polar}) is below A²/2π ({disc}), a disc's of the same area ({area}), the least any section has"
        )));
    }
    Ok(())
}

impl BeamProfile {
    /// `rectangle`, `circle`, `pipe`, `i`, `box` or `general` (the `type` it is written with).
    pub fn name(&self) -> &'static str {
        match self {
            BeamProfile::Rectangle { .. } => "rectangle",
            BeamProfile::Circle { .. } => "circle",
            BeamProfile::Pipe { .. } => "pipe",
            BeamProfile::I { .. } => "i",
            BeamProfile::Box { .. } => "box",
            BeamProfile::General { .. } => "general",
        }
    }

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
            BeamProfile::General { area, iy, iz, torsion, shear_y, shear_z } => {
                for (what, v) in [("area", area), ("iy", iy), ("iz", iz), ("torsion", torsion), ("shear_y", shear_y), ("shear_z", shear_z)] {
                    positive(what, v)?;
                }
                general_bounds(area, iy, iz, torsion, shear_y, shear_z)?;
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
            BeamProfile::General { area, .. } => area,
        }
    }

    /// The inner and outer radii of a pipe with a bore (a wall thinner than its radius, with a
    /// ring b² − a² the numbers resolve): what a pipe pressure acts on. `None` for every other
    /// profile.
    pub fn bore(&self) -> Option<(f64, f64)> {
        let BeamProfile::Pipe { radius, thickness } = *self else { return None };
        let inner = radius - thickness;
        (inner > 0.0 && radius * radius - inner * inner > 0.0).then_some((inner, radius))
    }

    /// The section properties ([`BeamProperties`]); `poisson` enters the shear coefficients of
    /// the shapes (a `general` profile's are its own).
    pub fn properties(&self, poisson: f64) -> BeamProperties {
        let pi = std::f64::consts::PI;
        let area = self.area();
        match *self {
            BeamProfile::Rectangle { width, height } => {
                let k = 10.0 * (1.0 + poisson) / (12.0 + 11.0 * poisson);
                BeamProperties {
                    area,
                    iy: width * height.powi(3) / 12.0,
                    iz: height * width.powi(3) / 12.0,
                    torsion: rectangle_torsion(width, height),
                    shear_y: k * area,
                    shear_z: k * area,
                }
            }
            BeamProfile::Circle { radius } => {
                let i = pi * radius.powi(4) / 4.0;
                let k = 6.0 * (1.0 + poisson) / (7.0 + 6.0 * poisson);
                BeamProperties { area, iy: i, iz: i, torsion: 2.0 * i, shear_y: k * area, shear_z: k * area }
            }
            BeamProfile::Pipe { radius, thickness } => {
                let inner = radius - thickness;
                let i = pi * (radius.powi(4) - inner.powi(4)) / 4.0;
                let m2 = (inner / radius).powi(2);
                let k = 6.0 * (1.0 + poisson) * (1.0 + m2).powi(2) / ((7.0 + 6.0 * poisson) * (1.0 + m2).powi(2) + (20.0 + 12.0 * poisson) * m2);
                BeamProperties { area, iy: i, iz: i, torsion: 2.0 * i, shear_y: k * area, shear_z: k * area }
            }
            BeamProfile::I { height, width, web, flange } => {
                let inside = height - 2.0 * flange;
                BeamProperties {
                    area,
                    iy: (width * height.powi(3) - (width - web) * inside.powi(3)) / 12.0,
                    iz: (2.0 * flange * width.powi(3) + inside * web.powi(3)) / 12.0,
                    torsion: (2.0 * width * flange.powi(3) + (height - flange) * web.powi(3)) / 3.0,
                    shear_y: 5.0 / 6.0 * 2.0 * width * flange,
                    shear_z: (height - flange) * web,
                }
            }
            BeamProfile::Box { width, height, thickness: t } => {
                let (wi, hi) = (width - 2.0 * t, height - 2.0 * t);
                let (wm, hm) = (width - t, height - t);
                BeamProperties {
                    area,
                    iy: (width * height.powi(3) - wi * hi.powi(3)) / 12.0,
                    iz: (height * width.powi(3) - hi * wi.powi(3)) / 12.0,
                    torsion: 4.0 * (wm * hm).powi(2) * t / (2.0 * (wm + hm)),
                    shear_y: 2.0 * wm * t,
                    shear_z: 2.0 * hm * t,
                }
            }
            BeamProfile::General { area, iy, iz, torsion, shear_y, shear_z } => BeamProperties { area, iy, iz, torsion, shear_y, shear_z },
        }
    }
}

/// A profile's cross-section properties in the beam's local axes: y along the profile's width
/// (local 1, the section's orientation vector made normal to the axis), z along its height
/// (local 2).
///
/// Area, second moments I_y = ∫ z² dA and I_z = ∫ y² dA, the Saint-Venant torsion constant J
/// and the shear areas:
/// - rectangle: J from Saint-Venant's series, J = a b³/3 (1 − 192 b/(π⁵ a) Σ_{n odd}
///   tanh(nπa/2b)/n⁵) for sides a ≥ b (S. P. Timoshenko, J. N. Goodier, *Theory of Elasticity*,
///   3rd ed., McGraw-Hill 1970, §109); shear coefficient 10(1 + ν)/(12 + 11ν);
/// - circle and pipe: J = I_y + I_z; shear coefficients 6(1 + ν)/(7 + 6ν) and, for a tube of
///   radius ratio m, 6(1 + ν)(1 + m²)² / ((7 + 6ν)(1 + m²)² + (20 + 12ν) m²) (G. R. Cowper,
///   *The shear coefficient in Timoshenko's beam theory*, J. Appl. Mech. 33(2):335–340, 1966);
/// - I and box (thin-walled): J = Σ b t³/3 over the midline walls (open I) and Bredt's
///   4 A_m² t / p (closed box); shear areas the walls parallel to the shear, between the midlines
///   of the walls across it (the web for shear along z, 5/6 of the flanges along y; the two
///   webs or the two flanges of a box) (W. D. Pilkey, *Analysis and Design of Elastic Beams*,
///   Wiley 2002, ch. 6).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamProperties {
    pub area: f64,
    /// ∫ z² dA: bending with deflection along z.
    pub iy: f64,
    /// ∫ y² dA: bending with deflection along y.
    pub iz: f64,
    /// Saint-Venant torsion constant.
    pub torsion: f64,
    /// Shear areas for shear along y and along z.
    pub shear_y: f64,
    pub shear_z: f64,
}

/// Saint-Venant's torsion constant of a solid rectangle.
fn rectangle_torsion(width: f64, height: f64) -> f64 {
    let (a, b) = (width.max(height), width.min(height));
    let pi = std::f64::consts::PI;
    let series: f64 = (0..50).map(|k| f64::from(2 * k + 1)).map(|n| libm::tanh(n * pi * a / (2.0 * b)) / n.powi(5)).sum();
    a * b * b * b / 3.0 * (1.0 - 192.0 * b / (pi.powi(5) * a) * series)
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
    /// With `twist` every node also carries the hoop displacement u_θ (degree of freedom 3), so
    /// the body can be twisted about its axis: the shears γ_rθ = ∂u_θ/∂r − u_θ/r and
    /// γ_zθ = ∂u_θ/∂z are stress and strain components xz and yz.
    Axisymmetric {
        #[serde(default, skip_serializing_if = "is_false")]
        twist: bool,
    },
}

fn is_false(flag: &bool) -> bool {
    !*flag
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
            SectionKind::Shell { .. } | SectionKind::PlaneStress { .. } | SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric { .. } => 2,
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

    /// Whether this is an axisymmetric section with twist (a third unknown, the hoop
    /// displacement, at every node).
    pub fn twist(&self) -> bool {
        matches!(self, SectionKind::Axisymmetric { twist: true })
    }

    /// How a 2D continuum section extends out of its plane (`None` for other kinds).
    pub fn extent(&self) -> Option<Extent> {
        match *self {
            SectionKind::PlaneStress { thickness } | SectionKind::PlaneStrain { thickness } => Some(Extent::Thickness(thickness)),
            SectionKind::Axisymmetric { .. } => Some(Extent::Revolved),
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
            SectionKind::Axisymmetric { .. } => "axisymmetric",
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
            SectionKind::Axisymmetric { .. } => {}
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

/// The beam theory of a beam section's elements.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BeamTheory {
    /// Shear-deformable, with rotary inertia.
    #[default]
    Timoshenko,
    /// Shear-rigid, without rotary inertia of the cross-section in bending.
    EulerBernoulli,
}

impl BeamTheory {
    /// `timoshenko` or `euler_bernoulli`.
    pub fn name(&self) -> &'static str {
        match self {
            BeamTheory::Timoshenko => "timoshenko",
            BeamTheory::EulerBernoulli => "euler_bernoulli",
        }
    }

    pub fn is_timoshenko(&self) -> bool {
        *self == BeamTheory::Timoshenko
    }
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
    /// The beam theory: beam sections only.
    #[serde(default, skip_serializing_if = "BeamTheory::is_timoshenko")]
    pub beam_theory: BeamTheory,
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
        let volume_constrained = kind.dimension() == 3 || matches!(kind, SectionKind::PlaneStrain { .. } | SectionKind::Axisymmetric { .. });
        if self.hybrid && !volume_constrained {
            return Err(Error::Invalid(format!("hybrid formulations are for solid, plane strain and axisymmetric sections, not {}", kind.name())));
        }
        if !self.beam_theory.is_timoshenko() && !matches!(kind, SectionKind::Beam { .. }) {
            return Err(Error::Invalid(format!("the {} beam theory is for beam sections, not {}", self.beam_theory.name(), kind.name())));
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
        assert_eq!(BeamProfile::Pipe { radius: 2.0, thickness: 0.5 }.bore(), Some((1.5, 2.0)));
        assert_eq!(BeamProfile::Pipe { radius: 2.0, thickness: 2.0 }.bore(), None, "no bore");
        assert_eq!(BeamProfile::Pipe { radius: 1e20, thickness: 1.0 }.bore(), None, "a wall the numbers cannot resolve");
        assert_eq!(BeamProfile::Circle { radius: 2.0 }.bore(), None);

        let shell = SectionKind::Shell { thickness: 2.0, offset: 0.0, points: 5 };
        assert!(Formulation { hybrid: true, ..Formulation::default() }.validate_for(&shell).is_err());
        assert!(
            Formulation { incompatible_modes: true, integration: Integration::Reduced, ..Formulation::default() }
                .validate_for(&SectionKind::Solid {})
                .is_err()
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
    fn beam_theory_is_a_beam_section_formulation() {
        let beam = SectionKind::Beam { profile: BeamProfile::Circle { radius: 1.0 }, orientation: [0.0, 0.0, 1.0] };
        let shear_rigid = Formulation { beam_theory: BeamTheory::EulerBernoulli, ..Formulation::default() };
        shear_rigid.validate_for(&beam).unwrap();
        for other in [SectionKind::Solid {}, SectionKind::Shell { thickness: 1.0, offset: 0.0, points: 5 }, SectionKind::Rod { area: 1.0 }] {
            let e = shear_rigid.validate_for(&other).unwrap_err().to_string();
            assert!(e.contains("euler_bernoulli beam theory is for beam sections") && e.contains(other.name()), "{e}");
        }
        // Written by name; the default (Timoshenko) is left out, so files without it read and write unchanged.
        assert_eq!(serde_json::to_value(shear_rigid).unwrap()["beam_theory"], "euler_bernoulli");
        assert!(serde_json::to_value(Formulation::default()).unwrap().get("beam_theory").is_none());
        let parsed: Formulation = serde_json::from_value(serde_json::json!({"beam_theory": "euler_bernoulli"})).unwrap();
        assert_eq!(parsed, shear_rigid);
        assert_eq!(serde_json::from_value::<Formulation>(serde_json::json!({})).unwrap().beam_theory, BeamTheory::Timoshenko);
        assert!(serde_json::from_value::<Formulation>(serde_json::json!({"beam_theory": "rayleigh"})).is_err());
    }

    #[test]
    fn profile_properties() {
        // A square: J = 0.1406 a⁴ (Timoshenko & Goodier's table, 0.1406); thin strip J → b t³/3.
        let sq = BeamProfile::Rectangle { width: 2.0, height: 2.0 }.properties(0.3);
        assert!((sq.torsion / 16.0 - 0.1406).abs() < 1e-4, "{}", sq.torsion / 16.0);
        let strip = BeamProfile::Rectangle { width: 1000.0, height: 1.0 }.properties(0.3);
        assert!((strip.torsion / (1000.0 / 3.0) - 1.0).abs() < 1e-3, "{}", strip.torsion);
        // A thin pipe: the shear coefficient tends to 2(1 + ν)/(4 + 3ν), J = 2 π r³ t.
        let (r, t) = (100.0, 0.01);
        let pipe = BeamProfile::Pipe { radius: r, thickness: t }.properties(0.3);
        let k = pipe.shear_y / pipe.area;
        assert!((k - 2.6 / 4.9).abs() < 1e-3, "{k}");
        assert!((pipe.torsion / (2.0 * std::f64::consts::PI * r.powi(3) * t) - 1.0).abs() < 1e-3);
        // A box is a pair of I sections in bending: the same I_y with webs of twice the thickness.
        let b = BeamProfile::Box { width: 10.0, height: 20.0, thickness: 1.0 }.properties(0.3);
        let i = BeamProfile::I { height: 20.0, width: 10.0, web: 2.0, flange: 1.0 }.properties(0.3);
        assert!((b.iy - i.iy).abs() < 1e-9 && (b.area - i.area).abs() < 1e-12);
        let circle = BeamProfile::Circle { radius: 2.0 }.properties(0.3);
        assert!((circle.torsion - std::f64::consts::PI * 8.0).abs() < 1e-12);
    }

    /// A general profile is its section properties: given a shape's, it is that shape to every
    /// consumer of [`BeamProfile::properties`] (the beam element and the deck writers).
    #[test]
    fn a_general_profile_is_the_shape_whose_properties_it_carries() {
        for shape in [
            BeamProfile::Rectangle { width: 20.0, height: 60.0 },
            BeamProfile::Pipe { radius: 20.0, thickness: 2.0 },
            BeamProfile::I { height: 100.0, width: 50.0, web: 5.0, flange: 8.0 },
        ] {
            let p = shape.properties(0.3);
            let general = BeamProfile::General { area: p.area, iy: p.iy, iz: p.iz, torsion: p.torsion, shear_y: p.shear_y, shear_z: p.shear_z };
            general.validate().unwrap();
            assert_eq!(general.properties(0.3), p, "{shape:?}");
            assert_eq!(general.properties(0.0), p, "its own shear areas, whatever the Poisson's ratio");
            assert_eq!(general.area(), shape.area());
        }
        let json = serde_json::json!({"type": "general", "area": 1.0, "iy": 2.0, "iz": 3.0, "torsion": 4.0, "shear_y": 0.5, "shear_z": 0.6});
        let parsed: BeamProfile = serde_json::from_value(json.clone()).unwrap();
        assert_eq!(parsed, BeamProfile::General { area: 1.0, iy: 2.0, iz: 3.0, torsion: 4.0, shear_y: 0.5, shear_z: 0.6 });
        assert_eq!((parsed.name(), serde_json::to_value(parsed).unwrap()), ("general", json));
        assert!(serde_json::from_value::<BeamProfile>(serde_json::json!({"type": "general", "area": 1.0, "iy": 2.0})).is_err(), "every property");
        let good = [1.0, 2.0, 3.0, 4.0, 0.5, 0.6];
        for (i, bad) in [0.0, -1.0, f64::NAN, f64::INFINITY].into_iter().flat_map(|v| (0..6).map(move |i| (i, v))) {
            let mut v = good;
            v[i] = bad;
            let [area, iy, iz, torsion, shear_y, shear_z] = v;
            assert!(BeamProfile::General { area, iy, iz, torsion, shear_y, shear_z }.validate().is_err(), "{v:?}");
        }
    }

    /// The bounds every section obeys: a shear area at most the area, a torsion constant at most
    /// the polar moment, a polar moment at least a disc's. The shapes meet them, a disc and a
    /// tube with equality; numbers that break one are refused, naming it.
    #[test]
    fn a_general_profile_no_section_could_have_is_refused() {
        let shapes = [
            BeamProfile::Rectangle { width: 20.0, height: 60.0 },
            BeamProfile::Rectangle { width: 1.0, height: 1.0 },
            BeamProfile::Circle { radius: 3.0 },
            BeamProfile::Circle { radius: 0.0123 },
            BeamProfile::Pipe { radius: 20.0, thickness: 2.0 },
            BeamProfile::Pipe { radius: 0.05, thickness: 0.05 },
            BeamProfile::I { height: 100.0, width: 50.0, web: 5.0, flange: 8.0 },
            BeamProfile::Box { width: 40.0, height: 80.0, thickness: 3.0 },
        ];
        for shape in shapes {
            for nu in [0.0, 0.3, 0.5] {
                let p = shape.properties(nu);
                let general = BeamProfile::General { area: p.area, iy: p.iy, iz: p.iz, torsion: p.torsion, shear_y: p.shear_y, shear_z: p.shear_z };
                assert!(general.validate().is_ok(), "{shape:?} ν = {nu}: {:?}", general.validate());
            }
        }
        let p = BeamProfile::Rectangle { width: 20.0, height: 60.0 }.properties(0.3);
        let with = |f: &dyn Fn(&mut [f64; 6])| {
            let mut v = [p.area, p.iy, p.iz, p.torsion, p.shear_y, p.shear_z];
            f(&mut v);
            let [area, iy, iz, torsion, shear_y, shear_z] = v;
            BeamProfile::General { area, iy, iz, torsion, shear_y, shear_z }.validate().map_err(|e| e.to_string())
        };
        // The reviewer's case: shear_y = 10 A, which would have written K2 = 10.
        let err = with(&|v| v[4] = 10.0 * v[0]).unwrap_err();
        assert!(err.contains("shear area along local 1") && err.contains("shear_y"), "{err}");
        let err = with(&|v| v[5] = 1.001 * v[0]).unwrap_err();
        assert!(err.contains("shear area along local 2") && err.contains("shear_z"), "{err}");
        assert!(with(&|v| v[4] = v[0]).is_ok(), "a shear coefficient of 1, the bound itself");
        let err = with(&|v| v[3] = 1.001 * (v[1] + v[2])).unwrap_err();
        assert!(err.contains("torsion constant") && err.contains("polar moment"), "{err}");
        let err = with(&|v| {
            v[1] = 0.999 * v[0] * v[0] / (4.0 * std::f64::consts::PI);
            v[2] = v[1];
            v[3] = v[1];
        })
        .unwrap_err();
        assert!(err.contains("below A²/2π"), "{err}");
        // A disc's own numbers sit on two bounds; a hair over either is refused.
        let disc = BeamProfile::Circle { radius: 1.0 }.properties(0.3);
        let general =
            |torsion: f64, i: f64| BeamProfile::General { area: disc.area, iy: i, iz: i, torsion, shear_y: disc.shear_y, shear_z: disc.shear_z };
        assert!(general(disc.torsion, disc.iy).validate().is_ok());
        assert!(general(disc.torsion * 1.000_001, disc.iy).validate().is_err());
        assert!(general(disc.torsion * 0.5, disc.iy * 0.999_999).validate().is_err());
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
