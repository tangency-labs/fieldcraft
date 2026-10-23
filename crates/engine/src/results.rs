//! What the results display shows ([`ResultsView`], set by the `results.*` commands) and the
//! nodal values it draws: displacements, mode shapes, explicit field frames and stress fields,
//! read from a solved step's [`StepResult`]; a fluid step's cell fields are drawn by
//! [`crate::flow`]. An explicit frame's one-point element stresses are
//! averaged to the nodes (unweighted, like the implicit solver's nodal stresses).
//!
//! The view is session state, like the camera: changing it is not an edit, so it is neither
//! journalled nor undone, and it does not change the session revision. A result kept from an
//! earlier revision is still shown, marked stale ([`Session::latest_result`]).

use fieldcraft_explicit::{ExplicitResult, Frame};
use fieldcraft_model::Model;
use fieldcraft_model::units::DisplaySystem;
pub use fieldcraft_model::units::Quantity;
use fieldcraft_solver::{FrequencyResult, HeatResult, StaticResult};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::{EngineError, Result, Session, StepResult};

/// Loads and reactions balance when they differ by less than this fraction of the largest.
pub use fieldcraft_lbm::{ComponentCheck, FluidResult, Peak};
pub use fieldcraft_solver::system::BALANCE_TOLERANCE;

/// The largest deformation scale factor.
pub const MAX_DEFORMATION_SCALE: f64 = 1e12; // unit-factors: a cap on a scale factor

/// A result by its procedure.
pub enum Kind<'a> {
    Static(&'a StaticResult),
    Frequency(&'a FrequencyResult),
    Explicit(&'a ExplicitResult),
    /// Per cell, not per node: drawn by [`crate::flow`].
    Fluid(&'a fieldcraft_lbm::FluidResult),
    Heat(&'a HeatResult),
}

impl StepResult {
    pub fn kind(&self) -> Kind<'_> {
        match self {
            StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) => Kind::Static(r),
            StepResult::Implicit(fieldcraft_solver::StepResult::Frequency(r)) => Kind::Frequency(r),
            StepResult::Implicit(fieldcraft_solver::StepResult::Heat(r)) => Kind::Heat(r),
            StepResult::Explicit(r) => Kind::Explicit(r),
            StepResult::Fluid(r) => Kind::Fluid(r),
        }
    }
}

/// Which part of a result is shown: the mode (from 1, frequency steps) and the field frame
/// (explicit steps; `None` is the last).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Pick {
    pub mode: usize,
    pub frame: Option<usize>,
}

/// The temperature frame of a heat result `pick` asks for.
pub fn heat_frame(r: &HeatResult, pick: Pick) -> Result<(f64, &[f64])> {
    let last = r.frames.len().checked_sub(1).ok_or_else(|| EngineError::Params(format!("step `{}` has no temperatures", r.step)))?;
    let k = pick.frame.unwrap_or(last);
    match (r.times.get(k), r.frames.get(k)) {
        (Some(t), Some(f)) => Ok((*t, f)),
        _ => Err(EngineError::Params(format!("step `{}` has frames 0 to {last}; there is no frame {k}", r.step))),
    }
}

/// The field frame of an explicit result `pick` asks for.
pub fn frame(r: &ExplicitResult, pick: Pick) -> Result<&Frame> {
    let last = r.frames.len().checked_sub(1).ok_or_else(|| EngineError::Params(format!("step `{}` has no field frames", r.step)))?;
    let k = pick.frame.unwrap_or(last);
    r.frames.get(k).ok_or_else(|| EngineError::Params(format!("step `{}` has frames 0 to {last}; there is no frame {k}", r.step)))
}

/// A nodal quantity drawn as contours.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Field {
    /// Displacement magnitude (a mode's shape magnitude for a frequency step).
    Displacement,
    DisplacementX,
    DisplacementY,
    DisplacementZ,
    /// Velocity magnitude (explicit steps; per cell in fluid steps).
    Velocity,
    /// Gauge pressure (per cell, fluid steps).
    Pressure,
    VonMises,
    MaxPrincipal,
    MidPrincipal,
    MinPrincipal,
    StressXx,
    StressYy,
    StressZz,
    StressXy,
    StressXz,
    StressYz,
    /// Temperature (heat transfer and coupled steps).
    Temperature,
    /// Heat flux magnitude (heat transfer steps, at the end of the step).
    HeatFlux,
    /// Signed fluid heat received by an obstacle face per geometric area.
    WallHeatFlux,
    /// Contact pressure (nonlinear static steps with contact: 0 off the contacts' secondary
    /// surfaces).
    ContactPressure,
    /// The normalised contact gap of the secondary surfaces' nodes (negative: penetration).
    ContactGap,
    /// Contact status of the secondary surfaces' nodes: 0 open, 1 closed (frictionless) or
    /// sticking, 2 slipping.
    ContactStatus,
    /// The magnitude of the friction's shear traction (0 off the secondary surfaces).
    ContactShear,
    /// The magnitude of the frictional slip over the last increment (secondary surfaces' nodes).
    ContactSlip,
}

impl Field {
    pub const ALL: [Field; 24] = [
        Field::Displacement,
        Field::DisplacementX,
        Field::DisplacementY,
        Field::DisplacementZ,
        Field::Velocity,
        Field::Pressure,
        Field::VonMises,
        Field::MaxPrincipal,
        Field::MidPrincipal,
        Field::MinPrincipal,
        Field::StressXx,
        Field::StressYy,
        Field::StressZz,
        Field::StressXy,
        Field::StressXz,
        Field::StressYz,
        Field::Temperature,
        Field::HeatFlux,
        Field::WallHeatFlux,
        Field::ContactPressure,
        Field::ContactGap,
        Field::ContactStatus,
        Field::ContactShear,
        Field::ContactSlip,
    ];

    /// The id commands take (`von_mises`…).
    pub fn id(self) -> &'static str {
        match self {
            Field::Displacement => "displacement",
            Field::DisplacementX => "displacement_x",
            Field::DisplacementY => "displacement_y",
            Field::DisplacementZ => "displacement_z",
            Field::Velocity => "velocity",
            Field::Pressure => "pressure",
            Field::VonMises => "von_mises",
            Field::MaxPrincipal => "max_principal",
            Field::MidPrincipal => "mid_principal",
            Field::MinPrincipal => "min_principal",
            Field::StressXx => "stress_xx",
            Field::StressYy => "stress_yy",
            Field::StressZz => "stress_zz",
            Field::StressXy => "stress_xy",
            Field::StressXz => "stress_xz",
            Field::StressYz => "stress_yz",
            Field::Temperature => "temperature",
            Field::HeatFlux => "heat_flux",
            Field::WallHeatFlux => "wall_heat_flux",
            Field::ContactPressure => "contact_pressure",
            Field::ContactGap => "contact_gap",
            Field::ContactStatus => "contact_status",
            Field::ContactShear => "contact_shear",
            Field::ContactSlip => "contact_slip",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Field::Displacement => "Displacement magnitude",
            Field::DisplacementX => "Displacement X",
            Field::DisplacementY => "Displacement Y",
            Field::DisplacementZ => "Displacement Z",
            Field::Velocity => "Velocity magnitude",
            Field::Pressure => "Pressure (gauge)",
            Field::VonMises => "Von Mises stress",
            Field::MaxPrincipal => "Max. principal stress",
            Field::MidPrincipal => "Mid. principal stress",
            Field::MinPrincipal => "Min. principal stress",
            Field::StressXx => "Stress XX",
            Field::StressYy => "Stress YY",
            Field::StressZz => "Stress ZZ",
            Field::StressXy => "Stress XY",
            Field::StressXz => "Stress XZ",
            Field::StressYz => "Stress YZ",
            Field::Temperature => "Temperature",
            Field::HeatFlux => "Heat flux magnitude",
            Field::WallHeatFlux => "Wall heat flux (into wall)",
            Field::ContactPressure => "Contact pressure",
            Field::ContactGap => "Contact gap",
            Field::ContactStatus => "Contact status (0 open, 1 closed or stick, 2 slip)",
            Field::ContactShear => "Contact shear traction",
            Field::ContactSlip => "Contact slip",
        }
    }

    /// What it measures; `None` for a code (a contact status).
    pub fn quantity(self) -> Option<Quantity> {
        Some(match self {
            Field::Displacement | Field::DisplacementX | Field::DisplacementY | Field::DisplacementZ | Field::ContactSlip | Field::ContactGap => {
                Quantity::Length
            }
            Field::Velocity => Quantity::Velocity,
            Field::Pressure => Quantity::Pressure,
            Field::Temperature => Quantity::Temperature,
            Field::HeatFlux | Field::WallHeatFlux => Quantity::HeatFlux,
            Field::ContactStatus => return None,
            Field::ContactPressure
            | Field::ContactShear
            | Field::VonMises
            | Field::MaxPrincipal
            | Field::MidPrincipal
            | Field::MinPrincipal
            | Field::StressXx
            | Field::StressYy
            | Field::StressZz
            | Field::StressXy
            | Field::StressXz
            | Field::StressYz => Quantity::Pressure,
        })
    }

    /// Whether `r` holds it as a mode shape: normalised, with no unit.
    pub fn is_mode_shape(self, r: &StepResult) -> bool {
        matches!(r.kind(), Kind::Frequency(_)) && self.quantity() == Some(Quantity::Length)
    }

    /// Its unit for result `r` in `display`: none for a mode shape (it is normalised).
    pub fn unit(self, r: &StepResult, display: &DisplaySystem) -> Option<String> {
        if self.is_mode_shape(r) { None } else { self.quantity().and_then(|q| unit_of(q, display)) }
    }

    /// Its title over result `r`: the label with its unit (`Von Mises stress [MPa]`), or, for
    /// a mode shape, saying it is normalised.
    pub fn title(self, r: &StepResult, display: &DisplaySystem) -> String {
        if self.is_mode_shape(r) {
            let label = self.label().trim_end_matches(" magnitude");
            return format!("{label} (mode shape, normalised)");
        }
        with_unit(self.label(), self.unit(r, display).as_deref())
    }

    /// A value of this field over `r`, stored in SI, as `display` shows it (a mode shape
    /// normalised as the display system's mass unit makes it, [`mode_scale`]).
    pub fn shown(self, r: &StepResult, si: f64, display: &DisplaySystem) -> f64 {
        if self.is_mode_shape(r) {
            return si * mode_scale(r, display).0;
        }
        match self.quantity() {
            Some(q) => shown(q, si, display),
            None => si,
        }
    }

    /// A value of this field over `r` as `display` shows it, stored in SI: the inverse of
    /// [`Field::shown`].
    pub fn stored(self, r: &StepResult, v: f64, display: &DisplaySystem) -> f64 {
        if self.is_mode_shape(r) {
            let k = mode_scale(r, display).0;
            return if k != 0.0 { v / k } else { v };
        }
        match self.quantity() {
            Some(q) => display.to_si(q, v).unwrap_or(f64::NAN),
            None => v,
        }
    }

    pub fn parse(id: &str) -> Option<Field> {
        Field::ALL.into_iter().find(|f| f.id() == id)
    }

    /// Whether it is a stress.
    pub fn is_stress(self) -> bool {
        !matches!(
            self,
            Field::Displacement
                | Field::DisplacementX
                | Field::DisplacementY
                | Field::DisplacementZ
                | Field::Velocity
                | Field::Pressure
                | Field::Temperature
                | Field::HeatFlux
                | Field::WallHeatFlux
        ) && !self.is_contact()
    }

    /// Whether it is a contact quantity.
    pub fn is_contact(self) -> bool {
        matches!(self, Field::ContactPressure | Field::ContactGap | Field::ContactStatus | Field::ContactShear | Field::ContactSlip)
    }

    /// Whether it is a heat transfer quantity.
    pub fn is_thermal(self) -> bool {
        matches!(self, Field::Temperature | Field::HeatFlux)
    }

    /// The fields a result has: displacements everywhere, stresses for static steps and explicit
    /// frames that kept them (mode shapes have no stress scale), velocity for explicit frames.
    pub fn available(r: &StepResult) -> impl Iterator<Item = Field> + '_ {
        let contact = matches!(r.kind(), Kind::Static(s) if !s.contacts.is_empty());
        let (stress, velocity, temperature, flux) = match r.kind() {
            Kind::Static(s) => (true, false, s.temperature.is_some(), false),
            Kind::Frequency(_) => (false, false, false, false),
            Kind::Explicit(e) => {
                e.frames.last().map_or((false, false, false, false), |f| (!f.stress.is_empty(), !f.velocity.is_empty(), false, false))
            }
            // Cell fields, drawn by `crate::flow`.
            Kind::Fluid(f) => {
                let mut fields = vec![Field::Velocity, Field::Pressure];
                if !f.frames.is_empty() && f.frames.iter().all(|frame| !frame.temperature.is_empty() && frame.temperature.len() == f.grid.len()) {
                    fields.push(Field::Temperature);
                    if f.heat.iter().any(|h| h.body.is_some()) {
                        fields.push(Field::WallHeatFlux);
                    }
                }
                return fields.into_iter();
            }
            Kind::Heat(_) => return Field::ALL.into_iter().filter(|f| f.is_thermal()).collect::<Vec<_>>().into_iter(),
        };
        Field::ALL
            .into_iter()
            .filter(move |f| match f {
                Field::Velocity => velocity,
                Field::Pressure => false,
                Field::WallHeatFlux => false,
                Field::Temperature => temperature,
                Field::HeatFlux => flux,
                f if f.is_contact() => contact,
                f => stress || !f.is_stress(),
            })
            .collect::<Vec<_>>()
            .into_iter()
    }
}

/// The unit of `q` in `display`: `None` for a number without one.
pub fn unit_of(q: Quantity, display: &DisplaySystem) -> Option<String> {
    let symbol = display.symbol(q);
    (symbol != "1").then(|| symbol.to_string())
}

/// An SI value of `q` as `display` shows it (not a number when it cannot be expressed).
pub fn shown(q: Quantity, si: f64, display: &DisplaySystem) -> f64 {
    display.express(q, si).unwrap_or(f64::NAN)
}

/// SI values of `q` as `display` shows them.
pub fn shown_all<const N: usize>(q: Quantity, si: [f64; N], display: &DisplaySystem) -> [f64; N] {
    si.map(|v| shown(q, v, display))
}

/// The factors a natural mode's shape, mass-normalised (φᵀ M φ = 1) in SI, takes to be
/// mass-normalised in `display`: (translations, rotations). The normalisation fixes the scale by
/// the mass unit: φ in a system of mass unit M and length unit L is √M (translations) and L √M
/// (rotations, per length) times φ in SI, where M and L are the units' sizes in kg and m. A
/// buckling mode is scaled to a largest component of 1 and needs neither.
pub fn mode_scale(r: &StepResult, display: &DisplaySystem) -> (f64, f64) {
    match r.kind() {
        Kind::Frequency(f) if f.buckling.is_none() => {
            let mass = display.unit_ref(Quantity::Mass).factor().sqrt();
            (mass, display.unit_ref(Quantity::Length).factor() * mass)
        }
        _ => (1.0, 1.0),
    }
}

/// The displacement in SI that one unit of a result's shape, as `display` shows it, stands for:
/// a deformation scale multiplies the shown shape. A static or explicit displacement is a length
/// (1: the scale is a ratio); a mode shape is drawn as the display units normalise it (a
/// frequency mode by the mass, [`mode_scale`]; a buckling mode to a largest component of one
/// length unit), so a scale of 5 draws it as large in every system as before SI storage.
pub fn deformation_unit(r: &StepResult, display: &DisplaySystem) -> f64 {
    match r.kind() {
        Kind::Frequency(f) if f.buckling.is_none() => mode_scale(r, display).1,
        Kind::Frequency(_) => display.unit_ref(Quantity::Length).factor(),
        Kind::Static(_) | Kind::Explicit(_) | Kind::Fluid(_) | Kind::Heat(_) => 1.0,
    }
}

/// The deformation scale `view` draws `r` with, as it is shown (the view's, or the automatic
/// one, [`auto_scale`], for a model `size` across in metres): the displacement drawn is the SI
/// one times this and [`deformation_unit`].
pub fn shown_scale(view: &ResultsView, r: &StepResult, pick: Pick, size: f64, display: &DisplaySystem) -> f64 {
    view.scale.unwrap_or_else(|| {
        let unit = deformation_unit(r, display);
        if unit > 0.0 { auto_scale(r, pick, size) / unit } else { 0.0 }
    })
}

/// `label` with `unit` in brackets when there is one: `Von Mises stress [MPa]`.
pub fn with_unit(label: &str, unit: Option<&str>) -> String {
    match unit {
        Some(u) => format!("{label} [{u}]"),
        None => label.to_string(),
    }
}

/// Most contour bands.
pub const MAX_BANDS: u32 = 24;
/// Fewest contour bands.
pub const MIN_BANDS: u32 = 2;
/// Contour bands when none are asked for.
pub const DEFAULT_BANDS: u32 = 12;

/// What the results display shows.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResultsView {
    /// The solved step shown.
    pub step: String,
    /// Draw the deformed shape.
    pub deformed: bool,
    /// Deformation scale factor; `None` scales the largest displacement to a tenth of the
    /// model's size.
    pub scale: Option<f64>,
    /// Contours of this field; `None` draws the deformed shape in element colours.
    pub contour: Option<Field>,
    pub bands: u32,
    /// The legend's range, clamped below the field's smallest or above its largest value (each
    /// `None` follows the field): a singular peak no longer squeezes the field into one band.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range_min: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range_max: Option<f64>,
    /// The mode shown, from 1 (frequency steps).
    pub mode: usize,
    /// The field frame shown (explicit steps); `None` is the last.
    #[serde(default)]
    pub frame: Option<usize>,
    /// Animate the deformation (a mode swings through ±, a static step grows from 0).
    pub animate: bool,
    /// A fluid step's slice, surface colours and streamlines.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flow: Option<crate::flow::FlowView>,
}

impl ResultsView {
    /// How a step's result is shown first: a linear static step undeformed with von Mises
    /// contours, so they are read on the true geometry (`results.deformed` shows the deformed
    /// shape, at the automatic scale); a frequency step's first mode deformed with displacement
    /// contours (a mode shape is a deformation); an explicit step at its last frame and a
    /// nonlinear static step at its last increment, deformed to true scale ([`true_scale`];
    /// with von Mises contours when they kept stresses); a heat step undeformed; a fluid step's
    /// velocity magnitude on a slice, the obstacles and streamlines ([`crate::flow::FlowView`]).
    pub fn for_result(r: &StepResult) -> ResultsView {
        let heat = matches!(r.kind(), Kind::Heat(_));
        let flow = match r.kind() {
            Kind::Fluid(f) => Some(crate::flow::FlowView::for_result(f)),
            _ => None,
        };
        let contour = if heat {
            Field::Temperature
        } else if flow.is_some() {
            Field::Velocity
        } else if Field::available(r).any(|f| f == Field::VonMises) {
            Field::VonMises
        } else {
            Field::Displacement
        };
        let true_scale = true_scale(r);
        ResultsView {
            step: r.step().to_string(),
            deformed: true_scale || matches!(r.kind(), Kind::Frequency(_)),
            scale: true_scale.then_some(1.0),
            contour: Some(contour),
            bands: DEFAULT_BANDS,
            range_min: None,
            range_max: None,
            mode: 1,
            frame: None,
            animate: false,
            flow,
        }
    }

    pub fn pick(&self) -> Pick {
        Pick { mode: self.mode, frame: self.frame }
    }

    /// The contour bands for a field ranging from `min` to `max`, clamped to the view's range.
    pub fn bands(&self, min: f64, max: f64) -> fieldcraft_render::contour::Bands {
        fieldcraft_render::contour::Bands::clamped(min, max, self.range_min, self.range_max, self.bands)
    }
}

/// The automatic deformation scale draws the largest displacement this fraction of the model's
/// size.
pub const AUTO_FRACTION: f64 = 0.1;

/// Whether a result's deformation is drawn at true scale even when the scale is automatic: an
/// explicit step's and a nonlinear static step's displacements are already large (finite strain,
/// contact), and magnifying them draws the part crushed through its neighbours. Linear static,
/// frequency and buckling results are magnified.
pub fn true_scale(r: &StepResult) -> bool {
    match r.kind() {
        Kind::Explicit(_) => true,
        Kind::Static(s) => s.nonlinear.is_some(),
        Kind::Frequency(_) | Kind::Fluid(_) | Kind::Heat(_) => false,
    }
}

/// The automatic deformation scale of `r` at `pick` on a model `size` across (its bounding box
/// diagonal): 1 for a [`true_scale`] result, else the factor that draws the largest
/// displacement [`AUTO_FRACTION`] of the size; 0 when nothing moves.
pub fn auto_scale(r: &StepResult, pick: Pick, size: f64) -> f64 {
    if true_scale(r) {
        return 1.0;
    }
    let largest = displacements(r, pick).map_or(0.0, |u| u.iter().map(|v| magnitude(*v)).filter(|m| m.is_finite()).fold(0.0, f64::max));
    if largest > 0.0 && size.is_finite() && size > 0.0 { AUTO_FRACTION * size / largest } else { 0.0 }
}

/// The nodal displacements a result's deformed shape uses: the static displacement, the shape
/// of the picked mode or the picked frame's displacement, in [`StepResult::nodes`] order.
pub fn displacements(r: &StepResult, pick: Pick) -> Result<&[[f64; 3]]> {
    match r.kind() {
        Kind::Static(s) => match (pick.frame, &s.nonlinear) {
            (Some(k), Some(h)) => h.frames.get(k).map(|f| f.displacement.as_slice()).ok_or_else(|| {
                EngineError::Params(format!("step `{}` has frames 0 to {}; there is no frame {k}", s.step, h.frames.len().saturating_sub(1)))
            }),
            _ => Ok(&s.displacement),
        },
        Kind::Frequency(f) => pick
            .mode
            .checked_sub(1)
            .and_then(|k| f.modes.get(k))
            .map(|m| m.shape.as_slice())
            .ok_or_else(|| EngineError::Params(format!("step `{}` has modes 1 to {}; there is no mode {}", f.step, f.modes.len(), pick.mode))),
        Kind::Explicit(e) => {
            let f = frame(e, pick)?;
            if f.displacement.is_empty() {
                return Err(EngineError::Params(format!("step `{}` kept no displacements in its field frames", e.step)));
            }
            Ok(&f.displacement)
        }
        Kind::Fluid(f) => Err(EngineError::Params(format!("step `{}` is a fluid step: its fields are per cell (solve.probe)", f.step))),
        Kind::Heat(h) => Err(EngineError::Params(format!("step `{}` is heat transfer and has no displacements", h.step))),
    }
}

/// A field's value at every node of the result, in [`StepResult::nodes`] order. `model` gives
/// the element connectivity an explicit frame's element stresses are averaged over.
pub fn field_values(model: &Model, r: &StepResult, field: Field, pick: Pick) -> Result<Vec<f64>> {
    let unavailable = || EngineError::Params(format!("{} is not available for step `{}`", field.label(), r.step()));
    let component = |f: fn(&[f64; 3]) -> f64| -> Result<Vec<f64>> { Ok(displacements(r, pick)?.iter().map(f).collect()) };
    let stresses = || -> Result<std::borrow::Cow<'_, [[f64; 6]]>> {
        match r.kind() {
            Kind::Static(s) => Ok(std::borrow::Cow::Borrowed(s.stress.as_slice())),
            Kind::Frequency(_) | Kind::Fluid(_) | Kind::Heat(_) => Err(unavailable()),
            Kind::Explicit(e) => {
                let f = frame(e, pick)?;
                if f.stress.is_empty() {
                    return Err(unavailable());
                }
                Ok(std::borrow::Cow::Owned(nodal_average(model, &e.nodes, &e.elements, &f.stress)))
            }
        }
    };
    let each = |f: &dyn Fn(&[f64; 6]) -> f64| -> Result<Vec<f64>> { Ok(stresses()?.iter().map(f).collect()) };
    let principal = |k: usize| each(&|t| fieldcraft_fem::principal(t).get(k).copied().unwrap_or(f64::NAN));
    match field {
        Field::Displacement => component(|v| magnitude(*v)),
        Field::DisplacementX => component(|v| v[0]),
        Field::DisplacementY => component(|v| v[1]),
        Field::DisplacementZ => component(|v| v[2]),
        Field::Velocity => match r.kind() {
            Kind::Explicit(e) => {
                let f = frame(e, pick)?;
                if f.velocity.is_empty() { Err(unavailable()) } else { Ok(f.velocity.iter().map(|v| magnitude(*v)).collect()) }
            }
            Kind::Static(_) | Kind::Frequency(_) | Kind::Fluid(_) | Kind::Heat(_) => Err(unavailable()),
        },
        // Per cell (crate::flow::cell_values), not per node.
        Field::Pressure | Field::WallHeatFlux => Err(unavailable()),
        Field::Temperature => match r.kind() {
            Kind::Heat(h) => Ok(heat_frame(h, pick)?.1.to_vec()),
            Kind::Static(StaticResult { temperature: Some(t), .. }) => Ok(t.clone()),
            _ => Err(unavailable()),
        },
        Field::HeatFlux => match r.kind() {
            Kind::Heat(h) => Ok(h.flux.iter().map(|v| magnitude(*v)).collect()),
            _ => Err(unavailable()),
        },
        Field::ContactPressure | Field::ContactGap | Field::ContactStatus | Field::ContactShear | Field::ContactSlip => match r.kind() {
            Kind::Static(s) if !s.contacts.is_empty() => Ok(contact_values(s, field)),
            _ => Err(unavailable()),
        },
        Field::VonMises => each(&fieldcraft_fem::von_mises),
        Field::MaxPrincipal => principal(0),
        Field::MidPrincipal => principal(1),
        Field::MinPrincipal => principal(2),
        Field::StressXx => each(&|t| t[0]),
        Field::StressYy => each(&|t| t[1]),
        Field::StressZz => each(&|t| t[2]),
        Field::StressXy => each(&|t| t[3]),
        Field::StressXz => each(&|t| t[4]),
        Field::StressYz => each(&|t| t[5]),
    }
}

/// A contact field at every node of a static result: the largest pressure, shear, slip or status
/// code, or the smallest gap, over the contacts a node is a secondary node of; off them 0
/// pressure and shear, and no gap, status or slip (NaN).
fn contact_values(s: &StaticResult, field: Field) -> Vec<f64> {
    use fieldcraft_solver::contact::Status;
    let zero = matches!(field, Field::ContactPressure | Field::ContactShear);
    let mut out = vec![if zero { 0.0 } else { f64::NAN }; s.nodes.len()];
    let larger = |v: &mut f64, x: f64| *v = if v.is_nan() { x } else { v.max(x) };
    for c in &s.contacts {
        for (k, label) in c.nodes.iter().enumerate() {
            let Some(v) = s.nodes.binary_search(label).ok().and_then(|i| out.get_mut(i)) else { continue };
            match field {
                Field::ContactPressure => larger(v, c.pressure.get(k).copied().unwrap_or(0.0)),
                Field::ContactShear => larger(v, c.shear.get(k).map_or(0.0, |t| magnitude(*t))),
                Field::ContactSlip => {
                    if let Some(g) = c.slip.get(k) {
                        larger(v, magnitude(*g));
                    }
                }
                Field::ContactGap => {
                    if let Some(Some(g)) = c.gap.get(k) {
                        *v = if v.is_nan() { *g } else { v.min(*g) };
                    }
                }
                Field::ContactStatus => {
                    let code = match c.status.get(k) {
                        Some(Status::Closed | Status::Stick) => 1.0,
                        Some(Status::Slip) => 2.0,
                        Some(Status::Open) | None => 0.0,
                    };
                    larger(v, code);
                }
                _ => {}
            }
        }
    }
    out
}

/// Element values (in `elements` order) averaged at the nodes (in `nodes` order) of the
/// elements around them; a node no element reaches gets NaN.
fn nodal_average(model: &Model, nodes: &[u32], elements: &[u32], values: &[[f64; 6]]) -> Vec<[f64; 6]> {
    let mut sum = vec![[0.0; 6]; nodes.len()];
    let mut count = vec![0u32; nodes.len()];
    for (label, v) in elements.iter().zip(values) {
        let Some(e) = model.element(*label) else { continue };
        for n in &e.nodes {
            let Ok(i) = nodes.binary_search(n) else { continue };
            if let (Some(s), Some(c)) = (sum.get_mut(i), count.get_mut(i)) {
                for (a, b) in s.iter_mut().zip(v) {
                    *a += b;
                }
                *c += 1;
            }
        }
    }
    sum.into_iter().zip(count).map(|(s, c)| if c == 0 { [f64::NAN; 6] } else { s.map(|x| x / f64::from(c)) }).collect()
}

/// The smallest and largest finite value, with the nodes they are at (`labels` in the same
/// order as `values`).
pub fn range(labels: &[u32], values: &[f64]) -> Option<((u32, f64), (u32, f64))> {
    let mut out: Option<((u32, f64), (u32, f64))> = None;
    for (&l, &v) in labels.iter().zip(values).filter(|(_, v)| v.is_finite()) {
        out = Some(match out {
            None => ((l, v), (l, v)),
            Some((lo, hi)) => (if v < lo.1 { (l, v) } else { lo }, if v > hi.1 { (l, v) } else { hi }),
        });
    }
    out
}

/// The principal values of a stress (xx, yy, zz, xy, xz, yz), largest first.
pub fn principal(t: &[f64; 6]) -> [f64; 3] {
    fieldcraft_fem::principal(t)
}

pub fn magnitude(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

impl Session {
    /// The latest result of a step and whether it describes the current state (`false`: the
    /// model changed since it was solved; it is stale).
    pub fn latest_result(&self, step: &str) -> Option<(&std::sync::Arc<StepResult>, bool)> {
        self.results.get(step).map(|(revision, r)| (r, *revision == self.revision))
    }

    /// What the results display shows: the view, whether its result is stale, the contour
    /// field's range, a frequency step's modes and an explicit step's frame times (the
    /// `results.*` commands' reply), every number in the display system ([`Session::units`]).
    pub fn results_report(&self) -> Value {
        let display = self.units();
        let Some(v) = &self.results_view else { return json!({"view": null}) };
        let Some((r, current)) = self.latest_result(&v.step) else { return json!({"view": view_json(v, None, display)}) };
        let mut out = json!({"view": view_json(v, Some(r), display), "stale": !current});
        let legend =
            |field: Field| json!({"min": v.range_min.map(|x| field.shown(r, x, display)), "max": v.range_max.map(|x| field.shown(r, x, display))});
        if let Some(f) = v.contour
            && let Ok(values) = field_values(&self.model, r, f, v.pick())
        {
            let ends = range(r.nodes(), &values).map(|((ln, lv), (hn, hv))| {
                json!({"min": {"node": ln, "value": f.shown(r, lv, display)}, "max": {"node": hn, "value": f.shown(r, hv, display)}})
            });
            out["contour"] = json!({"field": f, "label": f.label(), "unit": f.unit(r, display), "range": ends});
            if v.range_min.is_some() || v.range_max.is_some() {
                out["contour"]["legend"] = legend(f);
            }
        }
        let time = |t: f64| shown(Quantity::Time, t, display);
        match r.kind() {
            Kind::Frequency(fr) => {
                out["modes"] = fr
                    .modes
                    .iter()
                    .enumerate()
                    .map(|(i, m)| json!({"mode": i + 1, "frequency": shown(Quantity::Frequency, m.frequency, display)}))
                    .collect();
            }
            Kind::Explicit(e) => {
                out["frames"] = e.frames.iter().map(|f| time(f.time)).collect();
                if let Ok(f) = frame(e, v.pick()) {
                    out["time"] = json!(time(f.time));
                }
            }
            Kind::Fluid(f) => {
                out["frames"] = f.frames.iter().map(|f| time(f.time)).collect();
                if let Ok(frame) = crate::flow::frame_index(f, v.pick()) {
                    out["ports"] = fluid_ports(f, frame, display);
                }
                if let Value::Object(flow) = crate::flow::report(v, f, display) {
                    for (k, x) in flow {
                        out[k.as_str()] = x;
                    }
                }
                // A face held beyond the fluid's range is the extreme, as the legend shows it.
                if v.contour == Some(Field::Temperature) && out["contour"]["range"].is_object() {
                    for h in crate::flow::held_faces(self, f) {
                        let value = shown(Quantity::Temperature, h.value, display);
                        let end = json!({"body": h.body, "face": h.face, "value": value});
                        let range = &mut out["contour"]["range"];
                        if range["max"]["value"].as_f64().is_some_and(|m| value > m) {
                            range["max"] = end;
                        } else if range["min"]["value"].as_f64().is_some_and(|m| value < m) {
                            range["min"] = end;
                        }
                    }
                }
                if v.contour == Some(Field::WallHeatFlux)
                    && let Ok(faces) = crate::flow::wall::faces(self, f, v.pick())
                {
                    let flux = |q: f64| shown(Quantity::HeatFlux, q, display);
                    let values: Vec<f64> = faces.iter().map(|f| f.flux.map_or(f64::NAN, flux)).collect();
                    let end = |(k, value): (usize, f64)| json!({"body": faces.get(k).map(|f| &f.body), "face": faces.get(k).map(|f| &f.face), "value": value});
                    let range = crate::flow::cell_range(&values).map(|(a, b)| json!({"min": end(a), "max": end(b)}));
                    out["contour"] = json!({"field": Field::WallHeatFlux, "label": Field::WallHeatFlux.label(), "unit": unit_of(Quantity::HeatFlux, display), "range": range});
                    if v.range_min.is_some() || v.range_max.is_some() {
                        out["contour"]["legend"] = legend(Field::WallHeatFlux);
                    }
                    let mut faces_json = json!(faces);
                    crate::flow::wall::show(&mut faces_json, display);
                    out["wall_heat"] = json!({"faces": faces_json, "positive": "into wall", "history": "linear interpolation between recorded samples; no extrapolation", "unrecorded": "uncoloured", "units": crate::flow::wall::units(display)});
                }
            }
            Kind::Heat(h) => {
                out["heat_flows"] = heat_flows(h, display);
                out["frames"] = h.times.iter().map(|t| time(*t)).collect();
                if let Ok((t, _)) = heat_frame(h, v.pick()) {
                    out["time"] = json!(time(t));
                }
            }
            Kind::Static(s) => {
                if let Some(h) = &s.nonlinear {
                    out["frames"] = json!(h.load_factors());
                    let shown = v.frame.and_then(|k| h.frames.get(k)).or(h.frames.last());
                    out["load_factor"] = json!(shown.map(|f| f.load_factor));
                }
            }
        }
        out
    }

    /// Every step with a result, current or stale, by name.
    pub fn result_steps(&self) -> impl Iterator<Item = (&str, bool)> {
        self.results.iter().map(|(name, (revision, _))| (name.as_str(), *revision == self.revision))
    }
}

/// Open-face transport, already converted to display units. Both the solve reply and current
/// results view use this, so a reopened result and a selected earlier frame report the same
/// values. No complete energy-balance verdict is inferred from advective transport alone.
pub fn fluid_ports(r: &FluidResult, frame: usize, display: &DisplaySystem) -> Value {
    let ports = match r.ports(frame) {
        Ok(ports) => ports,
        Err(e) => return json!({"error": e.to_string()}),
    };
    let value = |q, v| shown(q, v, display);
    let stream = |s: &fieldcraft_lbm::Stream| {
        json!({
            "volume_flow": value(Quantity::VolumeFlow, s.volume_flow),
            "mass_flow": s.mass_flow.map(|v| value(Quantity::MassFlow, v)),
            "temperature": s.temperature.map(|v| value(Quantity::Temperature, v)),
            "pressure": s.pressure.map(|v| value(Quantity::Pressure, v)),
            "sensible_heat": s.sensible_heat.map(|v| value(Quantity::Power, v)),
        })
    };
    let faces: Vec<Value> = ports
        .iter()
        .map(|p| {
            json!({
                "face": p.face.name(), "area": value(Quantity::Area, p.area),
                "inflow": stream(&p.inflow), "outflow": stream(&p.outflow),
            })
        })
        .collect();
    let net_volume: f64 = ports.iter().map(|p| p.outflow.volume_flow - p.inflow.volume_flow).sum();
    let net_mass = r.port_properties.map(|p| p.density * net_volume);
    let net_heat = ports.iter().try_fold(0.0, |sum, p| Some(sum + p.outflow.sensible_heat? - p.inflow.sensible_heat?));
    let units: serde_json::Map<String, Value> =
        [Quantity::Area, Quantity::VolumeFlow, Quantity::MassFlow, Quantity::Temperature, Quantity::Pressure, Quantity::Power]
            .into_iter()
            .map(|q| (q.name().into(), json!(unit_of(q, display))))
            .collect();
    json!({
        "faces": faces,
        "net_volume_flow": value(Quantity::VolumeFlow, net_volume),
        "net_mass_flow": net_mass.map(|v| value(Quantity::MassFlow, v)),
        "net_sensible_heat": if ports.is_empty() { None } else { net_heat.map(|v| value(Quantity::Power, v)) },
        "positive": "outward",
        "units": units,
        "temperature_reference": value(Quantity::Temperature, r.temperature_reference),
        "mass_convention": "reference density (incompressible)",
        "heat_convention": "advective sensible heat relative to temperature_reference; excludes port conduction, storage, wall heat and sources",
        "source": "saved boundary-layer cell fields; separate inflow and outflow, flow-weighted temperature and static gauge pressure",
    })
}

/// The heat balance of a heat transfer result: each heat load's flow into the body (`loads`:
/// `{name, type, heat}`, negative out of it), the prescribed temperatures' total (`reactions`),
/// and `imbalance`, their sum: about 0 in a steady state, the heat stored per time in a
/// transient one. Heat flows (powers) in `display`.
pub fn heat_flows(h: &HeatResult, display: &DisplaySystem) -> Value {
    let power = |w: f64| shown(Quantity::Power, w, display);
    let loads: Vec<Value> = h.loads.iter().map(|l| json!({"name": l.name, "type": l.kind, "heat": power(l.heat)})).collect();
    let applied: f64 = h.loads.iter().map(|l| l.heat).sum();
    json!({"loads": loads, "reactions": power(h.reaction_total), "imbalance": power(applied + h.reaction_total)})
}

/// A view as a reply shows it: its legend range in its field's unit, its slice position a
/// length, in `display`.
pub fn view_json(v: &ResultsView, r: Option<&StepResult>, display: &DisplaySystem) -> Value {
    let mut out = json!(v);
    if let (Some(field), Some(r)) = (v.contour, r) {
        for (key, x) in [("range_min", v.range_min), ("range_max", v.range_max)] {
            if let Some(x) = x {
                out[key] = json!(field.shown(r, x, display));
            }
        }
    }
    if let Some(position) = v.flow.as_ref().and_then(|f| f.slice.as_ref()).and_then(|s| s.position) {
        out["flow"]["slice"]["position"] = json!(shown(Quantity::Length, position, display));
    }
    if let Some(seeds) = v.flow.as_ref().and_then(|f| f.streamlines.as_ref()).and_then(|s| s.seeds.as_ref()) {
        out["flow"]["streamlines"]["seeds"] = json!(seeds.iter().map(|p| shown_all(Quantity::Length, *p, display)).collect::<Vec<_>>());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmd::solve::tests::cantilever;
    use serde_json::json;

    #[test]
    fn fields_follow_the_result() {
        let mut s = Session::default();
        cantilever(&mut s);
        s.execute("solve.static", &json!({})).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        let (bend, current) = s.latest_result("Bend").unwrap();
        assert!(current);
        let mag = field_values(&s.model, bend, Field::Displacement, Pick { mode: 1, frame: None }).unwrap();
        let uz = field_values(&s.model, bend, Field::DisplacementZ, Pick { mode: 1, frame: None }).unwrap();
        assert_eq!(mag.len(), bend.nodes().len());
        for (m, z) in mag.iter().zip(&uz) {
            assert!(*m >= z.abs() - 1e-12);
        }
        // The tip bends down the most (3.8 mm, in metres); the clamped root does not move.
        let ((lo_node, lo), (_, hi)) = range(bend.nodes(), &uz).unwrap();
        assert!(lo < -3.8e-3 && hi.abs() < 1e-15, "{lo} {hi}");
        assert!([11, 22, 33, 44].contains(&lo_node), "{lo_node}");
        let vm = field_values(&s.model, bend, Field::VonMises, Pick { mode: 1, frame: None }).unwrap();
        let max = field_values(&s.model, bend, Field::MaxPrincipal, Pick { mode: 1, frame: None }).unwrap();
        let min = field_values(&s.model, bend, Field::MinPrincipal, Pick { mode: 1, frame: None }).unwrap();
        for ((v, a), b) in vm.iter().zip(&max).zip(&min) {
            assert!(a >= b && *v >= 0.0, "{v} {a} {b}");
        }
        let (modes, _) = s.latest_result("Modes").unwrap();
        assert_eq!(Field::available(modes).count(), 4);
        assert_eq!(Field::available(bend).count(), Field::ALL.len() - 10, "all but velocity, pressure, the thermal and the contact fields");
        assert!(field_values(&s.model, modes, Field::VonMises, Pick { mode: 1, frame: None }).is_err());
        assert!(field_values(&s.model, modes, Field::Displacement, Pick { mode: 9, frame: None }).is_err());
        assert_eq!(field_values(&s.model, modes, Field::Displacement, Pick { mode: 2, frame: None }).unwrap().len(), modes.nodes().len());
        // An edit leaves the results stale, not gone.
        s.execute("model.material.edit", &json!({"name": "Mat", "elastic": {"type": "isotropic", "young": 2000, "poisson": 0}, "density": 1}))
            .unwrap();
        assert!(s.result("Bend").is_none());
        assert_eq!(s.latest_result("Bend").map(|(_, current)| current), Some(false));
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Bend", false), ("Modes", false)]);
    }

    /// A step's result goes stale only when something it reads changes (#318); undoing the
    /// change makes it current again.
    #[test]
    fn unrelated_edits_keep_results_current() {
        let mut s = Session::default();
        cantilever(&mut s);
        s.execute("solve.static", &json!({})).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        // A set, a boundary condition and a step that neither solved step uses.
        s.execute("model.set.create", &json!({"name": "extra", "nodes": [1]})).unwrap();
        s.execute("model.boundary_condition.create", &json!({"name": "X", "type": "fixed", "target": {"set": "extra"}})).unwrap();
        s.execute("model.step.create", &json!({"name": "Other", "procedure": {"type": "static"}, "boundary_conditions": ["X"]})).unwrap();
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Bend", true), ("Modes", true)]);
        // The load only Bend lists: Bend goes stale, Modes (no loads) stays current.
        s.execute("model.load.edit", &json!({"name": "push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -2]})).unwrap();
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Bend", false), ("Modes", true)]);
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Bend", true), ("Modes", true)]);
        // The clamp's set, which both use.
        s.execute("model.set.edit", &json!({"name": "root", "nodes": [1, 12]})).unwrap();
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Bend", false), ("Modes", false)]);
        s.execute("edit.undo", &json!({})).unwrap();
        // A step's own definition.
        s.execute("model.step.edit", &json!({"name": "Modes", "procedure": {"type": "frequency", "modes": 3}, "boundary_conditions": ["clamp"]}))
            .unwrap();
        assert_eq!(s.result_steps().collect::<Vec<_>>(), vec![("Bend", true), ("Modes", false)]);
    }

    #[test]
    fn units_follow_the_unit_system() {
        let mut s = Session::default();
        cantilever(&mut s);
        s.execute("solve.static", &json!({})).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        let mm = &fieldcraft_model::DisplaySystem::default();
        let si = crate::units::si();
        let (bend, _) = s.latest_result("Bend").unwrap();
        let (modes, _) = s.latest_result("Modes").unwrap();
        assert_eq!(Field::VonMises.title(bend, mm), "Von Mises stress [MPa]");
        assert_eq!(Field::Displacement.title(bend, mm), "Displacement magnitude [mm]");
        assert_eq!(Field::Displacement.title(modes, mm), "Displacement (mode shape, normalised)");
        assert_eq!(Field::DisplacementZ.title(modes, mm), "Displacement Z (mode shape, normalised)");
        assert_eq!(Field::VonMises.title(bend, si), "Von Mises stress [Pa]");
        assert_eq!(Field::ContactStatus.unit(bend, mm), None);
        assert_eq!(Field::ContactGap.unit(bend, mm).as_deref(), Some("mm"));
        assert_eq!(unit_of(Quantity::Power, mm).as_deref(), Some("mW"));
        assert_eq!(unit_of(Quantity::Energy, mm).as_deref(), Some("mJ"));
        assert_eq!(unit_of(Quantity::Moment, mm).as_deref(), Some("N·mm"));
        assert_eq!(unit_of(Quantity::Energy, si).as_deref(), Some("J"));
        assert_eq!(unit_of(Quantity::Temperature, mm).as_deref(), Some("°C"));
        assert_eq!(unit_of(Quantity::Dimensionless, mm), None);
        // A natural mode is mass-normalised: in tonnes it is √1000 times its value in kilograms.
        assert_eq!(mode_scale(modes, mm).0, 1000f64.sqrt());
        assert_eq!(mode_scale(bend, mm), (1.0, 1.0));
    }

    #[test]
    fn field_ids_round_trip() {
        for f in Field::ALL {
            assert_eq!(Field::parse(f.id()), Some(f));
            assert_eq!(serde_json::to_value(f).unwrap(), json!(f.id()));
        }
        assert_eq!(Field::parse("nope"), None);
        assert_eq!(range(&[1, 2], &[f64::NAN, f64::INFINITY]), None);
        assert_eq!(range(&[1, 2, 3], &[2.0, -1.0, 5.0]), Some(((2, -1.0), (3, 5.0))));
    }
}
