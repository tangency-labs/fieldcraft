//! Typed parameter schemas: what a command's JSON parameters are, field by field.
//!
//! A command's [`CommandSpec::form`](crate::CommandSpec::form) lists its [`Field`]s: name, kind,
//! default, whether it is required, its unit dimension, its choices, the model collection it
//! names, nested groups and one-of variants. The schema describes the JSON the command already
//! takes; it never changes it. The UI draws its command forms from it, `fieldcraft-cli commands
//! --json` and the MCP `list_commands` tool give it as JSON Schema ([`json_schema`]).
//!
//! The schema is descriptive: the command's `run` stays the judge of its parameters. What the
//! schema cannot express (a mass given per direction, a table of load-path points) stays reachable
//! through the form's JSON view.

use fieldcraft_model::units::{DisplaySystem, Quantity};
use serde_json::{Map, Value, json};

use crate::Session;

/// Rows a list field holds at most in a form.
pub const MAX_ROWS: usize = 1000;
/// Names a reference list offers at most.
pub const MAX_NAMES: usize = 10_000;

/// One parameter of a command.
#[derive(Debug, Clone, Copy)]
pub struct Field {
    /// Its JSON key. Empty for [`Kind::Variants`] whose tag and fields sit in the parent object
    /// (an internally tagged definition such as a load: `{"name", "type": "force", "vector"}`).
    pub key: &'static str,
    pub label: &'static str,
    pub kind: Kind,
    pub required: bool,
    /// The value the form starts with, as JSON text; `None` leaves the field empty.
    pub default: Option<&'static str>,
    /// What an empty field means ("every body", "1/25 of the diagonal"), or a short note.
    pub hint: &'static str,
    /// What the field means, in a sentence or two, with a typical value: the form's tooltip and
    /// the JSON Schema description. The unit is added from the field's dimension ([`Field::help_text`]).
    pub help: &'static str,
    /// For a [`Kind::Choice`]: a one-line note per option value (what it needs), shown under the
    /// field for the option chosen.
    pub notes: &'static [(&'static str, &'static str)],
}

/// What a field holds.
#[derive(Debug, Clone, Copy)]
pub enum Kind {
    /// A number in the model's unit system (or a fixed unit such as degrees).
    Number {
        dim: Dim,
        min: Option<f64>,
        max: Option<f64>,
    },
    /// A whole number in `min..=max`.
    Integer {
        min: i64,
        max: i64,
    },
    Text,
    /// A file path; `save` when the command writes it.
    Path {
        save: bool,
    },
    Bool,
    /// One of these strings: (value, label).
    Choice(&'static [(&'static str, &'static str)]),
    /// A fixed-length list of numbers, one per label (`["x", "y", "z"]`); `nullable` entries may
    /// be null (a free degree of freedom).
    Vector {
        dim: Dim,
        labels: &'static [&'static str],
        nullable: bool,
    },
    /// The name of an item of the model or the open design.
    Ref(Collection),
    /// A list of such names.
    RefList(Collection),
    /// A choice per item: an object from the names of a collection's items to one of these
    /// strings (value, label); an item left out takes the command's default.
    RefChoices(Collection, &'static [(&'static str, &'static str)]),
    /// A nested object.
    Group(&'static [Field]),
    /// A nested object (or, with an empty key, the parent object) tagged by `tag`: the tag's value
    /// picks the variant whose fields apply.
    Variants {
        tag: &'static str,
        variants: &'static [Variant],
    },
    /// An object with exactly one of these keys (`{"set": …}` or `{"reference_point": …}`), or a
    /// bare string for a [`Kind::Unit`] option.
    Either(&'static [Field]),
    /// An option of [`Kind::Either`] written as its key alone (`"line"`).
    Unit,
    /// A list of objects with these fields, at most `max` long.
    Rows {
        fields: &'static [Field],
        max: usize,
    },
    /// A list of mesh labels, nodes or elements (at most [`MAX_LABELS`]).
    Labels(LabelOf),
    /// A list of points `[x, y, z]`, or of number rows with these `columns` (a load path's
    /// `[fraction, factor]` pairs), at most `max`.
    Points {
        dim: Dim,
        columns: &'static [&'static str],
        max: usize,
    },
    /// Any of these values, as a list: (the value as JSON text, label), e.g. `("1", "UX")` or
    /// `("\"stress\"", "Stress")`.
    Multi(&'static [(&'static str, &'static str)]),
    /// A number at least `min`, or the string `"auto"`.
    NumberOrAuto {
        dim: Dim,
        min: Option<f64>,
    },
    /// One number for every direction, or one per axis (`labels`, e.g. the material axes
    /// 1, 2, 3), each at least `min`.
    PerAxis {
        dim: Dim,
        min: Option<f64>,
        labels: &'static [&'static str],
    },
    /// A list of `min_len` to `labels.len()` positive numbers, one per label in order (a block's
    /// divisions: 3 in 3D, 2 in 2D); `whole` numbers only when set.
    Tuple {
        labels: &'static [&'static str],
        min_len: usize,
        whole: bool,
    },
    /// A name (`name`, a [`Kind::Ref`] or [`Kind::Choice`] whose key is unused) as a bare
    /// string, or an object of `fields` (`object` labels that choice): a mirror plane is a datum
    /// plane's name or a point and a normal.
    NameOr {
        name: &'static Field,
        object: &'static str,
        fields: &'static [Field],
    },
    /// A text per item: an object from the names of a collection's items to a string; an item
    /// left out takes the command's default.
    RefText(Collection),
    /// A list of element faces or edges, `{"element", "face"}` or `{"element", "edge"}` (at most
    /// [`MAX_LABELS`]).
    Sides(Side),
    /// Any JSON: what the schema does not describe.
    Json,
}

/// What a [`Kind::Labels`] list labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelOf {
    Node,
    Element,
}

/// What a [`Kind::Sides`] list holds: element faces or element edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    Face,
    Edge,
}

impl Side {
    /// The key of the local number in each entry, and its word.
    pub fn key(self) -> &'static str {
        match self {
            Side::Face => "face",
            Side::Edge => "edge",
        }
    }
}

impl LabelOf {
    pub fn label(self) -> &'static str {
        match self {
            LabelOf::Node => "node",
            LabelOf::Element => "element",
        }
    }
}

/// Most labels a [`Kind::Labels`] list holds.
pub const MAX_LABELS: usize = 1_000_000;

/// One variant of [`Kind::Variants`].
#[derive(Debug, Clone, Copy)]
pub struct Variant {
    pub tag: &'static str,
    pub label: &'static str,
    pub fields: &'static [Field],
}

/// Where a reference field's names come from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Collection {
    Set,
    Material,
    Section,
    ReferencePoint,
    Coupling,
    Contact,
    Tie,
    Mass,
    Spring,
    BoundaryCondition,
    Load,
    FluidDomain,
    Step,
    /// A solved step.
    SolvedStep,
    /// A body of the open design.
    Body,
    /// A datum plane, axis or point.
    Datum,
    /// A datum plane.
    DatumPlane,
    /// A datum axis.
    DatumAxis,
    /// A datum point.
    DatumPoint,
    Partition,
    /// A face merge (`geometry.faces.merge`).
    FaceMerge,
    /// A stitch (`geometry.stitch`).
    Stitch,
    /// A face split (`geometry.face.split`).
    FaceSplit,
    /// A size control (`mesh.size.set`).
    SizeControl,
    /// A defeature (`geometry.defeature`).
    Defeature,
}

impl Collection {
    pub fn label(self) -> &'static str {
        match self {
            Collection::Set => "set",
            Collection::Material => "material",
            Collection::Section => "section",
            Collection::ReferencePoint => "reference point",
            Collection::Coupling => "coupling",
            Collection::Contact => "contact",
            Collection::Tie => "tie",
            Collection::Mass => "point mass",
            Collection::Spring => "spring",
            Collection::BoundaryCondition => "boundary condition",
            Collection::Load => "load",
            Collection::FluidDomain => "fluid domain",
            Collection::Step => "step",
            Collection::SolvedStep => "solved step",
            Collection::Body => "body",
            Collection::Datum => "datum",
            Collection::DatumPlane => "datum plane",
            Collection::DatumAxis => "datum axis",
            Collection::DatumPoint => "datum point",
            Collection::Partition => "partition",
            Collection::FaceMerge => "face merge",
            Collection::Stitch => "stitch",
            Collection::FaceSplit => "face split",
            Collection::SizeControl => "size control",
            Collection::Defeature => "defeature",
        }
    }
}

impl Session {
    /// The names a reference field may take, in order (at most [`MAX_NAMES`]).
    pub fn names(&self, c: Collection) -> Vec<String> {
        let m = &self.model;
        let names: Box<dyn Iterator<Item = &str>> = match c {
            Collection::Set => Box::new(m.sets().keys().map(String::as_str)),
            Collection::Material => Box::new(m.materials().keys().map(String::as_str)),
            Collection::Section => Box::new(m.sections().keys().map(String::as_str)),
            Collection::ReferencePoint => Box::new(m.reference_points().keys().map(String::as_str)),
            Collection::Coupling => Box::new(m.couplings().keys().map(String::as_str)),
            Collection::Contact => Box::new(m.contacts().keys().map(String::as_str)),
            Collection::Tie => Box::new(m.ties().keys().map(String::as_str)),
            Collection::Mass => Box::new(m.masses().keys().map(String::as_str)),
            Collection::Spring => Box::new(m.springs().keys().map(String::as_str)),
            Collection::BoundaryCondition => Box::new(m.boundary_conditions().keys().map(String::as_str)),
            Collection::Load => Box::new(m.loads().keys().map(String::as_str)),
            Collection::FluidDomain => Box::new(m.fluid_domains().keys().map(String::as_str)),
            Collection::Step => Box::new(m.steps().iter().map(|s| s.name.as_str())),
            Collection::SolvedStep => Box::new(self.solved_steps()),
            Collection::Body => Box::new(self.design.iter().flat_map(|d| d.bodies.iter().map(|b| b.name.as_str()))),
            Collection::Datum => Box::new(m.datums().keys().map(String::as_str)),
            Collection::DatumPlane => Box::new(m.datums().iter().filter(|(_, d)| d.plane().is_some()).map(|(n, _)| n.as_str())),
            Collection::DatumAxis => Box::new(m.datums().iter().filter(|(_, d)| d.axis().is_some()).map(|(n, _)| n.as_str())),
            Collection::DatumPoint => Box::new(m.datums().iter().filter(|(_, d)| d.point().is_some()).map(|(n, _)| n.as_str())),
            Collection::Partition => Box::new(m.partitions().iter().map(|p| p.name.as_str())),
            Collection::FaceMerge => Box::new(m.face_merges().iter().map(|f| f.name.as_str())),
            Collection::Stitch => Box::new(m.stitches().iter().map(|f| f.name.as_str())),
            Collection::FaceSplit => Box::new(m.face_splits().iter().map(|f| f.name.as_str())),
            Collection::SizeControl => Box::new(m.size_controls().keys().map(String::as_str)),
            Collection::Defeature => Box::new(m.defeatures().iter().map(|f| f.name.as_str())),
        };
        names.take(MAX_NAMES).map(str::to_string).collect()
    }
}

/// What a number of a command measures. Most follow the display system (bare numbers are in its
/// unit; [`Dim::quantity`]); [`Dim::Millimetre`] and [`Dim::Degree`] keep the fixed unit their
/// commands always took ([`Dim::fixed`]); a few depend on the value beside them ([`Dim::resolve`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dim {
    /// No unit (a ratio, a count, a factor).
    None,
    Length,
    Area,
    InverseLength,
    Stress,
    Density,
    Force,
    Moment,
    Mass,
    Inertia,
    Time,
    Velocity,
    InverseVelocity,
    Acceleration,
    AngularVelocity,
    StrainRate,
    Stiffness,
    Damping,
    Viscosity,
    ForcePerVolume,
    VolumeFlow,
    Momentum,
    Power,
    HeatGeneration,
    HeatFlux,
    FilmCoefficient,
    Conductivity,
    SpecificHeat,
    Expansion,
    Temperature,
    TemperatureDifference,
    Hertz,
    /// An eigenvalue ω² (`1/s²`).
    Eigenvalue,
    /// Always millimetres, whatever the display system: lengths of the geometry commands, which
    /// SolidCraft geometry (in millimetres) gave them.
    Millimetre,
    /// Always degrees.
    Degree,
    /// A spring's stiffness: torsional (per radian) when the spring's `direction` is a rotation.
    SpringStiffness,
    /// A dashpot's coefficient, torsional likewise.
    SpringDamping,
    /// Prescribed motion `[u1, u2, u3, ur1, ur2, ur3]`: three lengths, then three angles.
    Motion,
    /// A reaction `[f1, f2, f3, m1, m2, m3]`: three forces, then three moments.
    Reaction,
    /// A displacement along the `direction` beside it: a length, or an angle for 4 to 6.
    DofDisplacement,
    /// A fan curve's `[flow, pressure rise]` rows.
    FanCurve,
    /// Any other quantity.
    Of(Quantity),
}

impl Dim {
    /// The quantity, for a dimension that has one of its own (not the dependent ones).
    pub fn quantity(self) -> Option<Quantity> {
        use Quantity as Q;
        Some(match self {
            Dim::None | Dim::SpringStiffness | Dim::SpringDamping | Dim::Motion | Dim::Reaction | Dim::DofDisplacement | Dim::FanCurve => {
                return None;
            }
            Dim::Length | Dim::Millimetre => Q::Length,
            Dim::Area => Q::Area,
            Dim::InverseLength => Q::InverseLength,
            Dim::Stress => Q::Pressure,
            Dim::Density => Q::Density,
            Dim::Force => Q::Force,
            Dim::Moment => Q::Moment,
            Dim::Mass => Q::Mass,
            Dim::Inertia => Q::RotaryInertia,
            Dim::Time => Q::Time,
            Dim::Velocity => Q::Velocity,
            Dim::InverseVelocity => Q::InverseVelocity,
            Dim::Acceleration => Q::Acceleration,
            Dim::AngularVelocity => Q::AngularVelocity,
            Dim::StrainRate => Q::StrainRate,
            Dim::Stiffness => Q::Stiffness,
            Dim::Damping => Q::Damping,
            Dim::Viscosity => Q::DynamicViscosity,
            Dim::ForcePerVolume => Q::ForcePerVolume,
            Dim::VolumeFlow => Q::VolumeFlow,
            Dim::Momentum => Q::Momentum,
            Dim::Power => Q::Power,
            Dim::HeatGeneration => Q::HeatGeneration,
            Dim::HeatFlux => Q::HeatFlux,
            Dim::FilmCoefficient => Q::HeatTransferCoefficient,
            Dim::Conductivity => Q::ThermalConductivity,
            Dim::SpecificHeat => Q::SpecificHeat,
            Dim::Expansion => Q::ThermalExpansion,
            Dim::Temperature => Q::Temperature,
            Dim::TemperatureDifference => Q::TemperatureDifference,
            Dim::Hertz => Q::Frequency,
            Dim::Eigenvalue => Q::Eigenvalue,
            Dim::Degree => Q::Angle,
            Dim::Of(q) => q,
        })
    }

    /// The fixed unit a bare number is in, for the dimensions that never follow the display
    /// system.
    pub fn fixed(self) -> Option<&'static str> {
        match self {
            Dim::Millimetre => Some("mm"),
            Dim::Degree => Some("deg"),
            _ => None,
        }
    }

    /// The quantity of entry `index` of a value of this dimension in `parent` (the object that
    /// holds it, for the dependent dimensions); `None` for a number without a unit.
    pub fn resolve(self, parent: Option<&Value>, index: usize) -> Option<Quantity> {
        use Quantity as Q;
        let dof = |key: &str| parent.and_then(|p| p.get(key)).and_then(|d| d.get("dof").or(Some(d))).and_then(Value::as_u64);
        match self {
            Dim::SpringStiffness => Some(if dof("direction").is_some_and(|d| d >= 4) { Q::RotationalStiffness } else { Q::Stiffness }),
            Dim::SpringDamping => Some(if dof("direction").is_some_and(|d| d >= 4) { Q::RotationalDamping } else { Q::Damping }),
            Dim::Motion => Some(if index < 3 { Q::Length } else { Q::Angle }),
            Dim::Reaction => Some(if index < 3 { Q::Force } else { Q::Moment }),
            Dim::DofDisplacement => {
                Some(if parent.and_then(|p| p.get("direction")).and_then(Value::as_u64).is_some_and(|d| d >= 4) { Q::Angle } else { Q::Length })
            }
            Dim::FanCurve => Some(if index == 0 { Q::VolumeFlow } else { Q::Pressure }),
            d => d.quantity(),
        }
    }

    /// The dimension's name ("pressure", "length"); `None` for [`Dim::None`].
    pub fn name(self) -> Option<&'static str> {
        Some(match self {
            Dim::None => return None,
            Dim::SpringStiffness => "stiffness",
            Dim::SpringDamping => "damping",
            Dim::Motion => "displacement",
            Dim::Reaction => "force",
            Dim::DofDisplacement => "displacement",
            Dim::FanCurve => "volume flow and pressure",
            d => d.quantity().map(Quantity::name)?,
        })
    }

    /// Its unit in `display` (a fixed unit stays), or `None` for a number without one. The
    /// dependent dimensions give the unit of their usual case (a translation, a spring along
    /// a direction).
    pub fn unit(self, display: &DisplaySystem) -> Option<String> {
        if let Some(fixed) = self.fixed() {
            return Some(fixed.to_string());
        }
        let q = self.resolve(None, 0)?;
        let symbol = display.symbol(q);
        (symbol != "1").then(|| symbol.to_string())
    }
}

// ---- Building fields (const, so schemas are static tables) ----

const fn field(key: &'static str, label: &'static str, kind: Kind) -> Field {
    Field { key, label, kind, required: false, default: None, hint: "", help: "", notes: &[] }
}

impl Field {
    /// Required.
    pub const fn req(self) -> Field {
        Field { required: true, ..self }
    }
    /// Starts with this JSON value.
    pub const fn def(self, json: &'static str) -> Field {
        Field { default: Some(json), ..self }
    }
    /// What an empty field means.
    pub const fn hint(self, hint: &'static str) -> Field {
        Field { hint, ..self }
    }
    /// What it means ([`Field::help`]).
    pub const fn help(self, help: &'static str) -> Field {
        Field { help, ..self }
    }

    /// A one-line note per option value ([`Field::notes`]).
    pub const fn notes(self, notes: &'static [(&'static str, &'static str)]) -> Field {
        Field { notes, ..self }
    }
    /// The note for this option value (the default when `None`).
    pub fn note(&self, value: Option<&str>) -> Option<&'static str> {
        let value = value.or_else(|| self.default.map(|d| d.trim_matches('"')))?;
        self.notes.iter().find(|(v, _)| *v == value).map(|(_, n)| *n)
    }

    /// Its physical dimension, for a number or a vector of numbers.
    pub fn dim(&self) -> Option<Dim> {
        match self.kind {
            Kind::Number { dim, .. }
            | Kind::Vector { dim, .. }
            | Kind::Points { dim, .. }
            | Kind::NumberOrAuto { dim, .. }
            | Kind::PerAxis { dim, .. } => Some(dim),
            _ => None,
        }
    }

    /// What its value is measured in, as a form shows it: the unit in `display` ("MPa"); `None`
    /// for a number without a unit.
    pub fn unit_label(&self, display: &DisplaySystem) -> Option<String> {
        self.dim()?.unit(display)
    }

    /// The help as the form's tooltip shows it: [`Field::help`] and the unit.
    pub fn help_text(&self, display: &DisplaySystem) -> String {
        let mut out = self.help.to_string();
        if let Some(unit) = self.unit_label(display) {
            if !out.is_empty() {
                out.push('\n');
            }
            let fixed = self.dim().and_then(Dim::fixed).is_some();
            out.push_str(&if fixed { format!("Unit: {unit} (always).") } else { format!("Unit: {unit}.") });
        }
        out
    }
}

pub const fn num(key: &'static str, label: &'static str, dim: Dim) -> Field {
    field(key, label, Kind::Number { dim, min: None, max: None })
}
pub const fn range(key: &'static str, label: &'static str, dim: Dim, min: f64, max: f64) -> Field {
    field(key, label, Kind::Number { dim, min: Some(min), max: Some(max) })
}
pub const fn positive(key: &'static str, label: &'static str, dim: Dim) -> Field {
    field(key, label, Kind::Number { dim, min: Some(0.0), max: None })
}
pub const fn int(key: &'static str, label: &'static str, min: i64, max: i64) -> Field {
    field(key, label, Kind::Integer { min, max })
}
pub const fn text(key: &'static str, label: &'static str) -> Field {
    field(key, label, Kind::Text)
}
pub const fn path(key: &'static str, label: &'static str, save: bool) -> Field {
    field(key, label, Kind::Path { save })
}
pub const fn flag(key: &'static str, label: &'static str) -> Field {
    field(key, label, Kind::Bool)
}
pub const fn choice(key: &'static str, label: &'static str, options: &'static [(&'static str, &'static str)]) -> Field {
    field(key, label, Kind::Choice(options))
}
pub const fn xyz(key: &'static str, label: &'static str, dim: Dim) -> Field {
    field(key, label, Kind::Vector { dim, labels: &["x", "y", "z"], nullable: false })
}
/// A vector whose entries may be null (left free).
pub const fn free_vector(key: &'static str, label: &'static str, dim: Dim, labels: &'static [&'static str]) -> Field {
    field(key, label, Kind::Vector { dim, labels, nullable: true })
}
pub const fn vector(key: &'static str, label: &'static str, dim: Dim, labels: &'static [&'static str]) -> Field {
    field(key, label, Kind::Vector { dim, labels, nullable: false })
}
pub const fn reference(key: &'static str, label: &'static str, c: Collection) -> Field {
    field(key, label, Kind::Ref(c))
}
pub const fn references(key: &'static str, label: &'static str, c: Collection) -> Field {
    field(key, label, Kind::RefList(c))
}

pub const fn ref_choices(key: &'static str, label: &'static str, c: Collection, options: &'static [(&'static str, &'static str)]) -> Field {
    field(key, label, Kind::RefChoices(c, options))
}
pub const fn group(key: &'static str, label: &'static str, fields: &'static [Field]) -> Field {
    field(key, label, Kind::Group(fields))
}
pub const fn variants(key: &'static str, label: &'static str, variants: &'static [Variant]) -> Field {
    field(key, label, Kind::Variants { tag: "type", variants })
}
pub const fn either(key: &'static str, label: &'static str, options: &'static [Field]) -> Field {
    field(key, label, Kind::Either(options))
}
pub const fn unit(key: &'static str, label: &'static str) -> Field {
    field(key, label, Kind::Unit)
}
pub const fn rows(key: &'static str, label: &'static str, fields: &'static [Field]) -> Field {
    field(key, label, Kind::Rows { fields, max: MAX_ROWS })
}
/// A list of node or element labels.
pub const fn labels(key: &'static str, label: &'static str, of: LabelOf) -> Field {
    field(key, label, Kind::Labels(of))
}
/// A list of element faces or edges.
pub const fn sides(key: &'static str, label: &'static str, side: Side) -> Field {
    field(key, label, Kind::Sides(side))
}
/// The columns of a point, `[x, y, z]`.
pub const XYZ: &[&str] = &["x", "y", "z"];
/// A list of points `[x, y, z]`.
pub const fn points(key: &'static str, label: &'static str, dim: Dim, max: usize) -> Field {
    field(key, label, Kind::Points { dim, columns: XYZ, max })
}
/// A table: a list of number rows with these columns.
pub const fn table(key: &'static str, label: &'static str, dim: Dim, columns: &'static [&'static str], max: usize) -> Field {
    field(key, label, Kind::Points { dim, columns, max })
}
/// Any of `options` (each a value as JSON text, and a label), as a list.
pub const fn multi(key: &'static str, label: &'static str, options: &'static [(&'static str, &'static str)]) -> Field {
    field(key, label, Kind::Multi(options))
}
/// A number at least `min`, or "auto".
pub const fn number_or_auto(key: &'static str, label: &'static str, dim: Dim, min: Option<f64>) -> Field {
    field(key, label, Kind::NumberOrAuto { dim, min })
}
/// One number at least `min`, or one per axis of `labels`.
pub const fn per_axis(key: &'static str, label: &'static str, dim: Dim, min: Option<f64>, labels: &'static [&'static str]) -> Field {
    field(key, label, Kind::PerAxis { dim, min, labels })
}
/// `min_len` to `labels.len()` positive numbers, whole ones when `whole`.
pub const fn tuple(key: &'static str, label: &'static str, labels: &'static [&'static str], min_len: usize, whole: bool) -> Field {
    field(key, label, Kind::Tuple { labels, min_len, whole })
}
/// A name as `name` gives it, or an object of `fields` (labelled `object`).
pub const fn name_or(key: &'static str, label: &'static str, name: &'static Field, object: &'static str, fields: &'static [Field]) -> Field {
    field(key, label, Kind::NameOr { name, object, fields })
}
/// A text per item of `c`.
pub const fn reference_text(key: &'static str, label: &'static str, c: Collection) -> Field {
    field(key, label, Kind::RefText(c))
}
pub const fn json(key: &'static str, label: &'static str) -> Field {
    field(key, label, Kind::Json)
}
pub const fn variant(tag: &'static str, label: &'static str, fields: &'static [Field]) -> Variant {
    Variant { tag, label, fields }
}

/// The `name` of a new or edited definition.
pub const NAME: Field = text("name", "Name").req().help("A name for it, unique among its kind; other definitions refer to it by this name.");

/// A [`Kind::Multi`] option's value: its JSON text parsed (text that is not JSON is a string).
pub fn multi_value(option: &str) -> Value {
    serde_json::from_str(option).unwrap_or_else(|_| Value::String(option.to_string()))
}

impl Kind {
    /// The variant of `tag`, for [`Kind::Variants`].
    pub fn variant(&self, tag: &str) -> Option<&'static Variant> {
        match self {
            Kind::Variants { variants, .. } => variants.iter().find(|v| v.tag == tag),
            _ => None,
        }
    }
}

// ---- Prefilling, switching variants, checking ----

/// `params` with the defaults of the fields it leaves out filled in (required groups created,
/// the active variant's defaults too). What it already has is kept.
pub fn prefill(fields: &[Field], params: &Value) -> Value {
    let mut obj = params.as_object().cloned().unwrap_or_default();
    fill(fields, &mut obj, 0);
    Value::Object(obj)
}

/// Nesting deeper than this is not followed (schemas are shallow; this bounds the recursion).
pub const MAX_DEPTH: usize = 8;

fn fill(fields: &[Field], obj: &mut Map<String, Value>, depth: usize) {
    if depth > MAX_DEPTH {
        return;
    }
    for f in fields {
        match f.kind {
            Kind::Variants { tag, .. } if f.key.is_empty() => {
                if !obj.contains_key(tag)
                    && let Some(d) = parse_default(f)
                {
                    obj.insert(tag.into(), d);
                }
                if let Some(v) = obj.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                    fill(v.fields, obj, depth + 1);
                }
            }
            Kind::Group(_) | Kind::Variants { .. } => {
                if !obj.contains_key(f.key) {
                    match parse_default(f) {
                        Some(d) => {
                            obj.insert(f.key.into(), d);
                        }
                        None if f.required => {
                            obj.insert(f.key.into(), Value::Object(Map::new()));
                        }
                        None => continue,
                    }
                }
                if let Some(Value::Object(child)) = obj.get_mut(f.key) {
                    fill_nested(f, child, depth + 1);
                }
            }
            _ => {
                if !obj.contains_key(f.key)
                    && let Some(d) = parse_default(f)
                {
                    obj.insert(f.key.into(), d);
                }
            }
        }
    }
}

/// Fill a group's or a (keyed) variant field's own object.
fn fill_nested(f: &Field, child: &mut Map<String, Value>, depth: usize) {
    match f.kind {
        Kind::Group(fields) => fill(fields, child, depth),
        Kind::Variants { tag, variants } => {
            if !child.contains_key(tag)
                && let Some(first) = variants.first()
            {
                child.insert(tag.into(), json!(first.tag));
            }
            if let Some(v) = child.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                fill(v.fields, child, depth);
            }
        }
        _ => {}
    }
}

fn parse_default(f: &Field) -> Option<Value> {
    f.default.and_then(|d| serde_json::from_str(d).ok())
}

/// The object for variant `tag` of variants field `f`, from `old`: the tag, the values of `old`
/// the new variant (or `common`, the fields beside the variants) also has, and the new variant's
/// defaults.
pub fn switch_variant(f: &Field, common: &[Field], old: &Map<String, Value>, tag: &str) -> Map<String, Value> {
    let Kind::Variants { tag: tag_key, .. } = f.kind else { return old.clone() };
    let Some(v) = f.kind.variant(tag) else { return old.clone() };
    let keep = |k: &str| v.fields.iter().chain(common).any(|g| g.key == k && g.key != tag_key);
    let mut out: Map<String, Value> = old.iter().filter(|(k, _)| keep(k)).map(|(k, v)| (k.clone(), v.clone())).collect();
    out.insert(tag_key.into(), json!(tag));
    fill(v.fields, &mut out, 0);
    out
}

/// A problem with a field's value: the field, where it is, and what is wrong. [`Problem::for_form`]
/// phrases it by the field's label (the form), [`Problem::for_agents`] by its JSON path (agents,
/// the engine's messages).
#[derive(Debug, Clone)]
pub struct Problem {
    /// The JSON path of the field (`elastic/young`, `rows/0/name`); empty for the parameters as a
    /// whole.
    pub path: String,
    /// The field it is about; `None` for the parameters as a whole. For a variant's tag, the
    /// [`Kind::Variants`] field.
    pub field: Option<&'static Field>,
    pub kind: ProblemKind,
}

/// What is wrong with a value.
#[derive(Debug, Clone, PartialEq)]
pub enum ProblemKind {
    /// Left out (or empty text).
    Required,
    /// Not a finite number.
    NotANumber,
    /// A number below the field's minimum.
    Below(f64),
    /// A number above the field's maximum.
    Above(f64),
    /// Not a whole number in `min..=max`.
    NotAWholeNumber {
        min: i64,
        max: i64,
    },
    NotText,
    NotABool,
    /// Not one of the field's choices (or, for [`Kind::Either`], of its options).
    NotAChoice,
    /// Not a list of `len` numbers.
    NotAVector {
        len: usize,
    },
    NotNames,
    NotAnObject,
    /// A variant tag the field does not have (`None`: no tag at all).
    UnknownVariant(Option<String>),
    /// More rows than `max`, or not a list.
    TooManyRows {
        max: usize,
    },
    /// Not a number greater than 0.
    NotPositive,
    /// Not a number in `min..=max`.
    NotInRange {
        min: f64,
        max: f64,
    },
    /// Not a list of at most `max` of `what` ("[x, y, z]", "labels").
    NotAList {
        max: usize,
        what: &'static str,
    },
    NotAName,
    NotAPath,
    /// Wrong for a reason the schema does not describe: what it must be or do instead
    /// ("names the region of 2D elements; a solid block fills its body").
    Invalid(String),
}

impl PartialEq for Problem {
    /// Problems are equal when they are about the same path and are the same problem (the field
    /// is where the path leads).
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path && self.kind == other.kind
    }
}

impl Problem {
    fn new(path: String, field: &'static Field, kind: ProblemKind) -> Self {
        Problem { path, field: Some(field), kind }
    }

    /// A problem with the value at `path`, its field not looked up yet ([`Problem::in_form`]).
    pub fn at(path: impl Into<String>, kind: ProblemKind) -> Self {
        Problem { path: path.into(), field: None, kind }
    }

    /// The problem with its field looked up in `fields` (when it has none yet and the path is
    /// there).
    pub fn in_form(mut self, fields: &'static [Field]) -> Self {
        if self.field.is_none() && !self.path.is_empty() {
            self.field = field_at(fields, &self.path);
        }
        self
    }

    /// What the value must be, without the field's name ("required", "at least 0").
    pub fn expected(&self) -> String {
        let f = self.field;
        match &self.kind {
            ProblemKind::Required => "required".into(),
            ProblemKind::NotANumber => "a number".into(),
            ProblemKind::Below(lo) => format!("at least {lo}"),
            ProblemKind::Above(hi) => format!("at most {hi}"),
            ProblemKind::NotAWholeNumber { min, max } => format!("a whole number in {min}..={max}"),
            ProblemKind::NotText => "text".into(),
            ProblemKind::NotABool => "true or false".into(),
            ProblemKind::NotAChoice => match f.map(|f| f.kind) {
                Some(Kind::Choice(options) | Kind::RefChoices(_, options)) => {
                    format!("one of {}", options.iter().map(|(o, _)| *o).collect::<Vec<_>>().join(", "))
                }
                Some(Kind::Either(options)) => format!("one of {}", options.iter().map(|o| o.key).collect::<Vec<_>>().join(", ")),
                Some(Kind::Multi(options)) => format!("a list of {}", options.iter().map(|(o, _)| *o).collect::<Vec<_>>().join(", ")),
                _ => "one of the choices".into(),
            },
            ProblemKind::NotAVector { len } => format!("{len} numbers"),
            ProblemKind::NotNames => "a list of names".into(),
            ProblemKind::NotAnObject if f.is_none() => "a JSON object".into(),
            ProblemKind::NotAnObject => "an object".into(),
            ProblemKind::UnknownVariant(Some(t)) => format!("a known {} (not `{t}`)", f.map_or("type".into(), |f| f.label.to_lowercase())),
            ProblemKind::UnknownVariant(None) => match f.map(|f| f.kind) {
                Some(Kind::Variants { tag, .. }) => format!("an object with a known `{tag}`"),
                _ => "an object with a known type".into(),
            },
            ProblemKind::TooManyRows { max } => format!("a list of at most {max} rows"),
            ProblemKind::NotPositive => "a positive number".into(),
            ProblemKind::NotInRange { min, max } => format!("a number in {min}..={max}"),
            ProblemKind::NotAList { max, what } => format!("a list of at most {max} {what}"),
            ProblemKind::NotAName => "a name".into(),
            ProblemKind::NotAPath => "a file path".into(),
            ProblemKind::Invalid(why) => why.clone(),
        }
    }

    /// The problem as agents read it, naming the field by its JSON path:
    /// "`elastic/young` is required", "`bands` must be a whole number in 2..=24".
    pub fn for_agents(&self) -> String {
        self.phrase(Naming::Key)
    }

    /// The problem as a form shows it, naming the field by its label and saying what to do:
    /// "Material: choose a material", "Young's modulus: must be at least 0".
    pub fn for_form(&self) -> String {
        self.phrase(Naming::Label)
    }

    /// The one phrasing of a problem, naming the field by its key or by its label.
    pub fn phrase(&self, naming: Naming) -> String {
        let by_label = naming == Naming::Label;
        let field = self.field.filter(|_| by_label);
        let subject = match field {
            Some(f) => format!("{}:", f.label),
            None if self.path.is_empty() => if by_label { "The parameters" } else { "the parameters" }.to_string(),
            None => format!("`{}`", self.path),
        };
        let required = self.kind == ProblemKind::Required;
        let what = match (field.map(|f| f.kind), &self.kind) {
            (_, ProblemKind::Invalid(why)) => why.clone(),
            (None, _) if required => "is required".into(),
            (None, _) => format!("must be {}", self.expected()),
            (Some(Kind::Ref(c)), _) if required => format!("choose {}", article(c.label())),
            (Some(Kind::RefList(c)), _) if required => format!("choose at least one {}", c.label()),
            (Some(Kind::Either(options)), _) => {
                format!("choose {}", options.iter().map(|o| article(&lower(o.label))).collect::<Vec<_>>().join(" or "))
            }
            (Some(Kind::Choice(options)), _) => format!("choose one of {}", options.iter().map(|(_, l)| *l).collect::<Vec<_>>().join(", ")),
            (Some(Kind::Variants { .. }), _) => format!("choose {}", article(&lower(field.map_or("", |f| f.label)))),
            (Some(Kind::Path { .. }), _) if required => "enter a path".into(),
            (Some(Kind::Text), _) if required && field.is_some_and(|f| f.key == "name") => "enter a name".into(),
            (Some(_), _) if required => "required".into(),
            (Some(_), _) => format!("must be {}", self.expected()),
        };
        format!("{subject} {what}")
    }
}

/// How a problem names its field: by its JSON path (agents) or by its label (a form).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Naming {
    Key,
    Label,
}

/// A command's problem with `params`, made as precise as its form `fields` can: when the
/// command names a group (`plastic`) and the form finds a problem with a field inside it
/// (`plastic/yield_stress` left out), that problem; otherwise the command's own, with its field
/// looked up.
pub fn refine(fields: &'static [Field], params: &Value, problem: Problem) -> Problem {
    let inside = |p: &Problem| p.path.strip_prefix(problem.path.as_str()).is_some_and(|r| r.starts_with('/'));
    if !problem.path.is_empty()
        && let Some(precise) = problems(fields, params).into_iter().find(inside)
    {
        return precise;
    }
    problem.in_form(fields)
}

/// "a material", "an element".
fn article(noun: &str) -> String {
    let an = noun.starts_with(['a', 'e', 'i', 'o', 'u']);
    format!("{} {noun}", if an { "an" } else { "a" })
}

/// A label lower-cased, but not its acronyms ("Reference point" → "reference point", "RBE2" stays).
fn lower(label: &str) -> String {
    let mut chars = label.chars();
    match (chars.next(), chars.next()) {
        (Some(c), Some(d)) if d.is_lowercase() => c.to_lowercase().chain(label.chars().skip(1)).collect(),
        _ => label.to_string(),
    }
}

/// What is missing or out of range in `params`: required fields left out, numbers outside
/// their range, values of the wrong kind. An empty list does not promise the command accepts
/// them (it checks more); a problem means it would not.
pub fn problems(fields: &'static [Field], params: &Value) -> Vec<Problem> {
    match params.as_object() {
        Some(obj) => object_problems(fields, obj),
        None => vec![Problem { path: String::new(), field: None, kind: ProblemKind::NotAnObject }],
    }
}

/// [`problems`] of an object's fields.
pub fn object_problems(fields: &'static [Field], obj: &Map<String, Value>) -> Vec<Problem> {
    let mut out = Vec::new();
    check(fields, obj, "", &mut out, 0);
    out
}

/// The top-level keys `params` may hold for these fields: each field's key, and for a flattened
/// one-of (a load's `type`) its tag and the fields of the variant `params` picks (of every
/// variant when it picks none).
pub fn accepted_keys(fields: &[Field], params: &Value) -> Vec<&'static str> {
    let mut keys = Vec::new();
    for f in fields {
        match f.kind {
            Kind::Variants { tag, variants } if f.key.is_empty() => {
                keys.push(tag);
                let picked = params.get(tag).and_then(Value::as_str).and_then(|t| variants.iter().find(|v| v.tag == t));
                for v in variants.iter().filter(|v| picked.is_none_or(|p| p.tag == v.tag)) {
                    for key in accepted_keys(v.fields, params) {
                        if !keys.contains(&key) {
                            keys.push(key);
                        }
                    }
                }
            }
            _ => keys.push(f.key),
        }
    }
    keys
}

/// The first top-level key of `params` that `fields` do not describe, with the keys they do:
/// a misspelt or misplaced parameter is refused rather than silently ignored.
pub fn unknown_key(fields: &[Field], params: &Value) -> Option<(String, Vec<&'static str>)> {
    let obj = params.as_object()?;
    let accepted = accepted_keys(fields, params);
    obj.keys().find(|k| !accepted.contains(&k.as_str())).map(|k| (k.chars().take(64).collect(), accepted))
}

/// The path of field `key` inside `prefix` (`elastic/young`).
pub fn join(prefix: &str, key: &str) -> String {
    if prefix.is_empty() {
        key.to_string()
    } else if key.is_empty() {
        prefix.to_string()
    } else {
        format!("{prefix}/{key}")
    }
}

/// The field at `path` (`elastic/young`, `rows/0/name`), for a key the engine names. A variant's
/// tag gives its [`Kind::Variants`] field; an [`Kind::Either`] option's key gives the option.
pub fn field_at(fields: &'static [Field], path: &str) -> Option<&'static Field> {
    let mut segments = path.split('/');
    let mut field = find(fields, segments.next()?)?;
    for (depth, segment) in segments.enumerate() {
        if depth > MAX_DEPTH {
            return None;
        }
        field = match field.kind {
            Kind::Rows { .. } if segment.parse::<usize>().is_ok() => field,
            Kind::Group(inner) | Kind::Rows { fields: inner, .. } | Kind::NameOr { fields: inner, .. } => find(inner, segment)?,
            Kind::Variants { tag, .. } if tag == segment => field,
            Kind::Variants { variants, .. } => variants.iter().find_map(|v| find(v.fields, segment))?,
            Kind::Either(options) => options.iter().find(|o| o.key == segment)?,
            _ => return None,
        };
    }
    Some(field)
}

/// The field of `fields` with key `key`, looking into flattened variants.
fn find(fields: &'static [Field], key: &str) -> Option<&'static Field> {
    fields.iter().find_map(|f| match f.kind {
        Kind::Variants { tag, variants } if f.key.is_empty() => {
            if tag == key {
                Some(f)
            } else {
                variants.iter().find_map(|v| find(v.fields, key))
            }
        }
        _ => (f.key == key).then_some(f),
    })
}

fn check(fields: &'static [Field], obj: &Map<String, Value>, prefix: &str, out: &mut Vec<Problem>, depth: usize) {
    if depth > MAX_DEPTH || out.len() > MAX_ROWS {
        return;
    }
    for f in fields {
        if let Kind::Variants { tag, .. } = f.kind
            && f.key.is_empty()
        {
            match obj.get(tag).and_then(Value::as_str) {
                Some(t) => match f.kind.variant(t) {
                    Some(v) => check(v.fields, obj, prefix, out, depth + 1),
                    None => out.push(Problem::new(join(prefix, tag), f, ProblemKind::UnknownVariant(Some(t.to_string())))),
                },
                None if f.required => out.push(Problem::new(join(prefix, tag), f, ProblemKind::Required)),
                None => {}
            }
            continue;
        }
        let path = join(prefix, f.key);
        match obj.get(f.key) {
            None | Some(Value::Null) if f.required => out.push(Problem::new(path, f, ProblemKind::Required)),
            None | Some(Value::Null) => {}
            Some(v) => {
                if let Some(kind) = value_problem(f, v, &path, out, depth) {
                    out.push(Problem::new(path, f, kind));
                }
            }
        }
    }
}

/// What is wrong with `v` as the value of `f`; nested problems go to `out`.
fn value_problem(f: &'static Field, v: &Value, path: &str, out: &mut Vec<Problem>, depth: usize) -> Option<ProblemKind> {
    match f.kind {
        Kind::Number { min, max, .. } => match (v.as_f64().filter(|x| x.is_finite()), min, max) {
            (None, ..) => Some(ProblemKind::NotANumber),
            (Some(x), Some(lo), _) if x < lo => Some(ProblemKind::Below(lo)),
            (Some(x), _, Some(hi)) if x > hi => Some(ProblemKind::Above(hi)),
            _ => None,
        },
        Kind::Integer { min, max } => match v.as_i64() {
            Some(x) if (min..=max).contains(&x) => None,
            _ => Some(ProblemKind::NotAWholeNumber { min, max }),
        },
        Kind::Text | Kind::Path { .. } | Kind::Ref(_) => match v.as_str() {
            Some(s) if f.required && s.trim().is_empty() => Some(ProblemKind::Required),
            Some(_) => None,
            None => Some(ProblemKind::NotText),
        },
        Kind::Bool => (!v.is_boolean()).then_some(ProblemKind::NotABool),
        Kind::Choice(options) => match v.as_str() {
            Some(s) if options.iter().any(|(o, _)| *o == s) => None,
            _ => Some(ProblemKind::NotAChoice),
        },
        Kind::Vector { labels, nullable, .. } => match v.as_array() {
            Some(a) if a.len() == labels.len() && a.iter().all(|x| x.as_f64().is_some_and(f64::is_finite) || (nullable && x.is_null())) => None,
            _ => Some(ProblemKind::NotAVector { len: labels.len() }),
        },
        Kind::RefList(_) => (!v.as_array().is_some_and(|a| a.iter().all(Value::is_string))).then_some(ProblemKind::NotNames),
        Kind::RefChoices(_, options) => match v.as_object() {
            Some(o) if o.values().all(|x| x.as_str().is_some_and(|s| options.iter().any(|(opt, _)| *opt == s))) => None,
            Some(_) => Some(ProblemKind::NotAChoice),
            None => Some(ProblemKind::NotAnObject),
        },
        Kind::Group(fields) => match v.as_object() {
            Some(o) => {
                check(fields, o, path, out, depth + 1);
                None
            }
            None => Some(ProblemKind::NotAnObject),
        },
        Kind::NameOr { name, fields, .. } => match v {
            Value::String(_) => value_problem(name, v, path, out, depth + 1),
            Value::Object(o) => {
                check(fields, o, path, out, depth + 1);
                None
            }
            _ => Some(ProblemKind::Invalid(format!("a {} or an object", name.label.to_lowercase()))),
        },
        Kind::RefText(_) => match v.as_object() {
            Some(o) if o.len() <= MAX_NAMES && o.values().all(Value::is_string) => None,
            Some(_) => Some(ProblemKind::NotText),
            None => Some(ProblemKind::NotAnObject),
        },
        Kind::Variants { tag, .. } => match v.as_object() {
            Some(o) => match o.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                Some(var) => {
                    check(var.fields, o, path, out, depth + 1);
                    None
                }
                None => Some(ProblemKind::UnknownVariant(None)),
            },
            None => Some(ProblemKind::NotAnObject),
        },
        Kind::Either(options) => {
            let ok = match v {
                Value::String(s) => options.iter().any(|o| o.key == s && matches!(o.kind, Kind::Unit)),
                Value::Object(o) if o.len() == 1 => o.iter().all(|(k, inner)| {
                    options.iter().find(|opt| opt.key == k).is_some_and(|opt| value_problem(opt, inner, path, out, depth + 1).is_none())
                }),
                _ => false,
            };
            (!ok).then_some(ProblemKind::NotAChoice)
        }
        Kind::Unit | Kind::Json => None,
        Kind::Sides(side) => match v.as_array() {
            Some(a)
                if a.len() <= MAX_LABELS
                    && a.iter().all(|x| {
                        x.as_object().is_some_and(|o| {
                            o.len() == 2
                                && o.get("element").and_then(Value::as_u64).is_some_and(|n| n >= 1 && n <= u64::from(u32::MAX))
                                && o.get(side.key()).and_then(Value::as_u64).is_some_and(|n| n >= 1 && n <= u64::from(u8::MAX))
                        })
                    }) =>
            {
                None
            }
            _ => Some(ProblemKind::NotAList { max: MAX_LABELS, what: if side == Side::Face { "element faces" } else { "element edges" } }),
        },
        Kind::Labels(_) => match v.as_array() {
            Some(a) if a.len() <= MAX_LABELS && a.iter().all(|x| x.as_u64().is_some_and(|n| n >= 1 && n <= u64::from(u32::MAX))) => None,
            _ => Some(ProblemKind::NotAList { max: MAX_LABELS, what: "labels" }),
        },
        Kind::Points { columns, max, .. } => match v.as_array() {
            Some(a)
                if a.len() <= max
                    && a.iter().all(|p| {
                        p.as_array().is_some_and(|c| c.len() == columns.len() && c.iter().all(|x| x.as_f64().is_some_and(f64::is_finite)))
                    }) =>
            {
                None
            }
            _ => Some(ProblemKind::NotAList { max, what: if columns == XYZ { "points [x, y, z]" } else { "rows of numbers" } }),
        },
        Kind::Multi(options) => match v.as_array() {
            Some(a) if a.iter().all(|x| options.iter().any(|(o, _)| multi_value(o) == *x)) => None,
            _ => Some(ProblemKind::NotAChoice),
        },
        Kind::NumberOrAuto { min, .. } => match (v, v.as_f64().filter(|x| x.is_finite())) {
            (Value::String(s), _) if s == "auto" => None,
            (_, Some(x)) if min.is_none_or(|lo| x >= lo) => None,
            (_, Some(_)) => Some(ProblemKind::Below(min.unwrap_or(0.0))),
            _ => Some(ProblemKind::Invalid("a number or \"auto\"".into())),
        },
        Kind::PerAxis { min, labels, .. } => {
            let values = match v {
                Value::Array(a) if a.len() == labels.len() => a.iter().map(|x| x.as_f64().filter(|x| x.is_finite())).collect::<Option<Vec<_>>>(),
                Value::Array(_) => None,
                _ => v.as_f64().filter(|x| x.is_finite()).map(|x| vec![x]),
            };
            match values {
                Some(xs) if xs.iter().all(|&x| min.is_none_or(|lo| x >= lo)) => None,
                Some(_) => Some(ProblemKind::Below(min.unwrap_or(0.0))),
                None => Some(ProblemKind::Invalid(format!("a number, or {} numbers (one per axis)", labels.len()))),
            }
        }
        Kind::Tuple { labels, min_len, whole } => match v.as_array() {
            Some(a) if (min_len..=labels.len()).contains(&a.len()) => {
                let ok = |x: &Value| x.as_f64().is_some_and(|x| x.is_finite() && x > 0.0 && (!whole || x.fract() == 0.0));
                (!a.iter().all(ok)).then_some(if whole {
                    ProblemKind::Invalid("whole numbers of at least 1".into())
                } else {
                    ProblemKind::NotPositive
                })
            }
            _ => Some(ProblemKind::Invalid(if min_len == labels.len() {
                format!("{min_len} numbers")
            } else {
                format!("{min_len} to {} numbers", labels.len())
            })),
        },
        Kind::Rows { fields, max } => match v.as_array() {
            Some(a) if a.len() <= max => {
                for (i, row) in a.iter().enumerate() {
                    let row_path = format!("{path}/{i}");
                    match row.as_object() {
                        Some(o) => check(fields, o, &row_path, out, depth + 1),
                        None => out.push(Problem::new(row_path, f, ProblemKind::NotAnObject)),
                    }
                }
                None
            }
            _ => Some(ProblemKind::TooManyRows { max }),
        },
    }
}

// ---- JSON Schema ----

/// The fields as a JSON Schema (draft 2020-12) object schema, units and hints in descriptions.
pub fn json_schema(fields: &[Field]) -> Value {
    object_schema(fields, 0)
}

fn object_schema(fields: &[Field], depth: usize) -> Value {
    let mut properties = Map::new();
    let mut required = Vec::new();
    let mut one_of = None;
    for f in fields {
        if let Kind::Variants { tag, variants } = f.kind
            && f.key.is_empty()
        {
            one_of = Some(variant_schemas(tag, variants, depth));
            if f.required {
                required.push(json!(tag));
            }
            continue;
        }
        properties.insert(f.key.into(), field_schema(f, depth));
        if f.required {
            required.push(json!(f.key));
        }
    }
    let mut s = json!({"type": "object", "properties": properties});
    if !required.is_empty() {
        s["required"] = Value::Array(required);
    }
    if let Some(o) = one_of {
        s["oneOf"] = o;
    }
    s
}

fn variant_schemas(tag: &str, variants: &[Variant], depth: usize) -> Value {
    Value::Array(
        variants
            .iter()
            .map(|v| {
                let mut s = if depth > MAX_DEPTH { json!({"type": "object"}) } else { object_schema(v.fields, depth + 1) };
                s["properties"][tag] = json!({"const": v.tag, "title": v.label});
                let mut req = s.get("required").and_then(Value::as_array).cloned().unwrap_or_default();
                req.push(json!(tag));
                s["required"] = Value::Array(req);
                s
            })
            .collect(),
    )
}

fn field_schema(f: &Field, depth: usize) -> Value {
    let mut s = match f.kind {
        Kind::Number { dim, min, max } => {
            let mut s = json!({"type": "number"});
            if let Some(lo) = min {
                s["minimum"] = json!(lo);
            }
            if let Some(hi) = max {
                s["maximum"] = json!(hi);
            }
            if dim != Dim::None {
                s["x-dimension"] = json!(format!("{dim:?}").to_lowercase());
            }
            s
        }
        Kind::Integer { min, max } => json!({"type": "integer", "minimum": min, "maximum": max}),
        Kind::Text | Kind::Path { .. } => json!({"type": "string"}),
        Kind::Bool => json!({"type": "boolean"}),
        Kind::Choice(options) => json!({"enum": options.iter().map(|(o, _)| *o).collect::<Vec<_>>()}),
        Kind::Vector { labels, nullable, .. } => {
            let item = if nullable { json!({"type": ["number", "null"]}) } else { json!({"type": "number"}) };
            json!({"type": "array", "items": item, "minItems": labels.len(), "maxItems": labels.len()})
        }
        Kind::Ref(c) => json!({"type": "string", "x-references": c.label()}),
        Kind::RefList(c) => json!({"type": "array", "items": {"type": "string", "x-references": c.label()}}),
        Kind::RefChoices(c, options) => json!({
            "type": "object",
            "x-references": c.label(),
            "additionalProperties": {"enum": options.iter().map(|(o, _)| *o).collect::<Vec<_>>()},
        }),
        Kind::Group(fields) if depth <= MAX_DEPTH => object_schema(fields, depth + 1),
        Kind::NameOr { name, fields, .. } if depth <= MAX_DEPTH => {
            json!({"oneOf": [field_schema(name, depth + 1), object_schema(fields, depth + 1)]})
        }
        Kind::RefText(c) => json!({"type": "object", "x-references": c.label(), "additionalProperties": {"type": "string"}}),
        Kind::Variants { tag, variants } if depth <= MAX_DEPTH => json!({"oneOf": variant_schemas(tag, variants, depth + 1)}),
        Kind::Either(options) if depth <= MAX_DEPTH => json!({"oneOf": options.iter().map(|o| match o.kind {
            Kind::Unit => json!({"const": o.key}),
            _ => json!({"type": "object", "properties": {o.key: field_schema(o, depth + 1)}, "required": [o.key], "additionalProperties": false}),
        }).collect::<Vec<_>>()}),
        Kind::Rows { fields, max } if depth <= MAX_DEPTH => json!({"type": "array", "items": object_schema(fields, depth + 1), "maxItems": max}),
        Kind::Unit => json!({"const": f.key}),
        Kind::Sides(side) => json!({
            "type": "array",
            "maxItems": MAX_LABELS,
            "items": {
                "type": "object",
                "properties": {"element": {"type": "integer", "minimum": 1, "maximum": u32::MAX}, side.key(): {"type": "integer", "minimum": 1, "maximum": u8::MAX}},
                "required": ["element", side.key()],
                "additionalProperties": false,
            },
        }),
        Kind::Labels(of) => {
            json!({"type": "array", "items": {"type": "integer", "minimum": 1, "maximum": u32::MAX}, "maxItems": MAX_LABELS, "x-labels": of.label()})
        }
        Kind::Points { dim, columns, max } => {
            let n = columns.len();
            let mut item = json!({"type": "array", "items": {"type": "number"}, "minItems": n, "maxItems": n, "x-columns": columns});
            if dim != Dim::None {
                item["x-dimension"] = json!(format!("{dim:?}").to_lowercase());
            }
            json!({"type": "array", "items": item, "maxItems": max})
        }
        Kind::Multi(options) => {
            json!({"type": "array", "items": {"enum": options.iter().map(|(o, _)| multi_value(o)).collect::<Vec<_>>()}, "uniqueItems": true})
        }
        Kind::Tuple { labels, min_len, whole } => {
            let item = if whole { json!({"type": "integer", "minimum": 1}) } else { json!({"type": "number", "exclusiveMinimum": 0}) };
            json!({"type": "array", "items": item, "minItems": min_len, "maxItems": labels.len(), "x-labels": labels})
        }
        Kind::NumberOrAuto { min, .. } => {
            let mut n = json!({"type": "number"});
            if let Some(lo) = min {
                n["minimum"] = json!(lo);
            }
            json!({"oneOf": [n, {"const": "auto"}]})
        }
        Kind::PerAxis { dim, min, labels } => {
            let mut n = json!({"type": "number"});
            if let Some(lo) = min {
                n["minimum"] = json!(lo);
            }
            let mut s =
                json!({"oneOf": [n.clone(), {"type": "array", "items": n, "minItems": labels.len(), "maxItems": labels.len(), "x-labels": labels}]});
            if dim != Dim::None {
                s["x-dimension"] = json!(format!("{dim:?}").to_lowercase());
            }
            s
        }
        _ => json!({}),
    };
    s["title"] = json!(f.label);
    let description = match (f.help, f.hint) {
        ("", "") => None,
        (help, "") => Some(help.to_string()),
        ("", hint) => Some(hint.to_string()),
        (help, hint) => Some(format!("{help} Empty or default: {hint}.")),
    };
    if let Some(d) = description {
        s["description"] = json!(d);
    }
    if let Some(d) = parse_default(f) {
        s["default"] = d;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn units_follow_the_system() {
        use fieldcraft_model::units::Preset;
        let (mm, si, inch) = (DisplaySystem::new(Preset::MmTonneS), DisplaySystem::new(Preset::MKgS), DisplaySystem::new(Preset::InLbfS));
        assert_eq!(Dim::Stress.unit(&mm).as_deref(), Some("MPa"));
        assert_eq!(Dim::Density.unit(&mm).as_deref(), Some("t/mm³"));
        assert_eq!(Dim::Moment.unit(&si).as_deref(), Some("N·m"));
        assert_eq!(Dim::Density.unit(&inch).as_deref(), Some("lbf·s²/in⁴"));
        assert_eq!(Dim::Temperature.unit(&mm).as_deref(), Some("°C"));
        assert_eq!(Dim::Millimetre.unit(&si).as_deref(), Some("mm"));
        assert_eq!(Dim::Degree.unit(&mm).as_deref(), Some("deg"));
        assert_eq!(Dim::None.unit(&si), None);
        // The dependent dimensions follow the value beside them.
        let torsional = json!({"direction": {"dof": 5}});
        assert_eq!(Dim::SpringStiffness.resolve(Some(&torsional), 0), Some(Quantity::RotationalStiffness));
        assert_eq!(Dim::SpringStiffness.resolve(Some(&json!({"direction": "line"})), 0), Some(Quantity::Stiffness));
        assert_eq!(Dim::Motion.resolve(None, 4), Some(Quantity::Angle));
        assert_eq!(Dim::FanCurve.resolve(None, 1), Some(Quantity::Pressure));
    }

    const LOADISH: &[Field] = &[
        NAME,
        variants(
            "",
            "Type",
            &[
                variant(
                    "force",
                    "Force",
                    &[xyz("vector", "Vector", Dim::Force).req(), choice("distribution", "Spread", &[("total", "Total")]).def("\"total\"")],
                ),
                variant("pressure", "Pressure", &[num("magnitude", "Magnitude", Dim::Stress).req(), reference("set", "Set", Collection::Set).req()]),
            ],
        )
        .req()
        .def("\"force\""),
    ];

    /// The typed lists check their values and describe themselves in JSON Schema (#169).
    #[test]
    fn labels_points_multi_and_auto_check() {
        const F: &[Field] = &[
            labels("nodes", "Nodes", LabelOf::Node),
            points("at", "Points", Dim::Length, 2),
            multi("dofs", "DOFs", &[("1", "UX"), ("2", "UY")]),
            multi("output", "Output", &[("\"stress\"", "Stress")]),
            number_or_auto("scale", "Scale", Dim::None, Some(0.0)),
        ];
        assert!(problems(F, &json!({"nodes": [1, 5], "at": [[0, 0, 1e3]], "dofs": [2], "output": ["stress"], "scale": "auto"})).is_empty());
        assert!(problems(F, &json!({"scale": 2.5})).is_empty());
        for bad in [
            json!({"nodes": [0]}),
            json!({"nodes": [1.5]}),
            json!({"nodes": [u64::from(u32::MAX) + 1]}),
            json!({"at": [[0, 0]]}),
            json!({"at": [[0, 0, 0], [0, 0, 0], [0, 0, 0]]}),
            json!({"at": [[0, 0, "x"]]}),
            json!({"dofs": [3]}),
            json!({"dofs": ["1"]}),
            json!({"output": ["strain"]}),
            json!({"scale": -1}),
            json!({"scale": "big"}),
        ] {
            assert_eq!(problems(F, &bad).len(), 1, "{bad}");
        }
        let s = json_schema(F);
        assert_eq!(s["properties"]["dofs"]["items"]["enum"], json!([1, 2]));
        assert_eq!(s["properties"]["output"]["items"]["enum"], json!(["stress"]));
        assert_eq!(s["properties"]["at"]["maxItems"], 2);
        assert_eq!(s["properties"]["scale"]["oneOf"][1]["const"], "auto");
    }

    /// A table is a list of number rows of its column count: a load path's [fraction, factor].
    #[test]
    fn table_rows_have_its_columns() {
        const F: &[Field] = &[table("path", "Path", Dim::None, &["fraction", "factor"], 3)];
        assert!(problems(F, &json!({"path": [[0, 0], [0.5, 1], [1, 1]]})).is_empty());
        for bad in [json!({"path": [[0, 0, 0]]}), json!({"path": [[0]]}), json!({"path": [[0, 0], [0, 0], [0, 0], [1, 1]]})] {
            assert_eq!(problems(F, &bad).len(), 1, "{bad}");
        }
        let s = json_schema(F);
        assert_eq!(s["properties"]["path"]["items"]["minItems"], 2);
        assert_eq!(s["properties"]["path"]["items"]["x-columns"], json!(["fraction", "factor"]));
    }

    /// A per-axis field takes one number or one per axis, each at least its minimum.
    #[test]
    fn per_axis_takes_one_number_or_one_per_axis() {
        const F: &[Field] = &[per_axis("k", "K", Dim::Conductivity, Some(0.0), &["1", "2", "3"])];
        assert!(problems(F, &json!({"k": 50})).is_empty());
        assert!(problems(F, &json!({"k": [50, 20, 20.5]})).is_empty());
        for bad in [json!({"k": -1}), json!({"k": [1, 2]}), json!({"k": [1, 2, -3]}), json!({"k": "x"}), json!({"k": [1, 2, "3"]})] {
            assert_eq!(problems(F, &bad).len(), 1, "{bad}");
        }
        let s = json_schema(F);
        assert_eq!(s["properties"]["k"]["oneOf"][1]["maxItems"], 3);
    }

    /// A name-or-object field takes the name's kind as a string or its fields as an object; a
    /// text per item takes an object of strings.
    #[test]
    fn name_or_object_and_text_per_item() {
        const NAME: Field = reference("", "Datum plane", Collection::DatumPlane);
        const OBJ: &[Field] = &[xyz("point", "Point", Dim::Length).req()];
        const F: &[Field] = &[name_or("plane", "Plane", &NAME, "Point", OBJ), reference_text("as", "Images", Collection::Body)];
        assert!(problems(F, &json!({"plane": "Sym", "as": {"A": "B"}})).is_empty());
        assert!(problems(F, &json!({"plane": {"point": [0, 0, 1]}})).is_empty());
        for bad in [json!({"plane": 3}), json!({"plane": {}}), json!({"plane": {"point": [0, 0]}}), json!({"as": {"A": 1}}), json!({"as": "B"})] {
            assert_eq!(problems(F, &bad).len(), 1, "{bad}");
        }
        assert_eq!(problems(F, &json!({"plane": {}})).first().map(|p| p.path.as_str()), Some("plane/point"));
        let s = json_schema(F);
        assert_eq!(s["properties"]["plane"]["oneOf"][0]["type"], "string");
        assert_eq!(s["properties"]["as"]["additionalProperties"]["type"], "string");
    }

    /// A sides list holds element and local face (or edge) objects, many more than a form's rows.
    #[test]
    fn sides_hold_element_and_local_number() {
        const F: &[Field] = &[sides("faces", "Faces", Side::Face), sides("edges", "Edges", Side::Edge)];
        let many: Vec<Value> = (1..=5000).map(|e| json!({"element": e, "face": 3})).collect();
        assert!(problems(F, &json!({"faces": many, "edges": [{"element": 2, "edge": 12}]})).is_empty());
        for bad in [
            json!({"faces": [{"element": 0, "face": 1}]}),
            json!({"faces": [{"element": 1, "face": 0}]}),
            json!({"faces": [{"element": 1, "edge": 1}]}),
            json!({"faces": [{"element": 1, "face": 1, "x": 1}]}),
            json!({"edges": [3]}),
        ] {
            assert_eq!(problems(F, &bad).len(), 1, "{bad}");
        }
        assert_eq!(json_schema(F)["properties"]["edges"]["items"]["required"], json!(["element", "edge"]));
    }

    /// A tuple takes its minimum to its full length of positive numbers, whole ones when asked.
    #[test]
    fn tuple_takes_two_or_three_positive_numbers() {
        const F: &[Field] = &[tuple("n", "N", &["a", "b", "c"], 2, true), tuple("r", "R", &["a", "b", "c"], 2, false)];
        assert!(problems(F, &json!({"n": [2, 3, 4], "r": [1.5, 1]})).is_empty());
        for bad in [
            json!({"n": [2]}),
            json!({"n": [1, 2, 3, 4]}),
            json!({"n": [1.5, 2]}),
            json!({"n": [0, 2]}),
            json!({"r": [1, -1]}),
            json!({"r": [1, null]}),
        ] {
            assert_eq!(problems(F, &bad).len(), 1, "{bad}");
        }
        let s = json_schema(F);
        assert_eq!(s["properties"]["n"]["minItems"], 2);
        assert_eq!(s["properties"]["n"]["items"]["type"], "integer");
    }

    #[test]
    fn prefill_switch_and_check() {
        let p = prefill(LOADISH, &json!({}));
        assert_eq!(p, json!({"type": "force", "distribution": "total"}));
        let problems = problems(LOADISH, &p);
        let paths: Vec<&str> = problems.iter().map(|p| p.path.as_str()).collect();
        assert_eq!(paths, ["name", "vector"]);
        let root = &LOADISH[1];
        let switched = switch_variant(
            root,
            LOADISH,
            &json!({"name": "a", "type": "force", "vector": [1, 2, 3]}).as_object().cloned().unwrap_or_default(),
            "pressure",
        );
        assert_eq!(Value::Object(switched), json!({"name": "a", "type": "pressure"}));
        assert!(super::problems(LOADISH, &json!({"name": "a", "type": "force", "vector": [1, "x", 3]})).iter().any(|p| p.path == "vector"));
        assert!(super::problems(LOADISH, &json!([1])).len() == 1);
    }

    /// One problem, two phrasings: by label for the form, by path for agents.
    #[test]
    fn problems_carry_their_field_and_kind() {
        let p = problems(LOADISH, &json!({"type": "pressure", "magnitude": "x"}));
        let kinds: Vec<(&str, &ProblemKind)> = p.iter().map(|p| (p.path.as_str(), &p.kind)).collect();
        assert_eq!(kinds, [("name", &ProblemKind::Required), ("magnitude", &ProblemKind::NotANumber), ("set", &ProblemKind::Required)]);
        assert_eq!(p.iter().map(Problem::for_form).collect::<Vec<_>>(), ["Name: enter a name", "Magnitude: must be a number", "Set: choose a set"]);
        assert_eq!(
            p.iter().map(Problem::for_agents).collect::<Vec<_>>(),
            ["`name` is required", "`magnitude` must be a number", "`set` is required"]
        );
        let unknown = problems(LOADISH, &json!({"name": "a", "type": "spin"}));
        assert_eq!(unknown.first().map(Problem::for_form).as_deref(), Some("Type: choose a type"));
        assert_eq!(unknown.first().map(Problem::for_agents).as_deref(), Some("`type` must be a known type (not `spin`)"));
        assert_eq!(problems(LOADISH, &json!([1])).first().map(Problem::for_agents).as_deref(), Some("the parameters must be a JSON object"));
        assert!(field_at(LOADISH, "magnitude").is_some_and(|f| f.label == "Magnitude"));
    }

    #[test]
    fn json_schema_has_variants_and_requirements() {
        let s = json_schema(LOADISH);
        assert_eq!(s["required"], json!(["name", "type"]));
        assert_eq!(s["oneOf"][1]["properties"]["type"]["const"], "pressure");
        assert_eq!(s["oneOf"][1]["required"], json!(["magnitude", "set", "type"]));
    }

    #[test]
    fn unknown_keys_are_found_per_variant() {
        assert_eq!(unknown_key(LOADISH, &json!({"name": "P", "type": "pressure", "magnitude": 1})), None);
        let (key, accepted) = unknown_key(LOADISH, &json!({"name": "P", "type": "pressure", "magnitud": 1})).unwrap();
        assert_eq!(key, "magnitud");
        assert!(accepted.contains(&"magnitude") && accepted.contains(&"name"));
    }

    #[test]
    fn every_command_refuses_an_unknown_parameter() {
        for spec in crate::commands() {
            let mut s = Session::default();
            let e = s.execute(spec.id, &json!({"bogus_param": 1})).unwrap_err().to_string();
            // A command that cannot run yet says so first; otherwise the key is named.
            assert!(e.contains("bogus_param") || s.availability(spec).is_err(), "{}: {e}", spec.id);
        }
        let mut s = Session::default();
        crate::cmd::solve::tests::cantilever(&mut s);
        s.execute("solve.run", &json!({})).unwrap();
        let e = s.execute("results.contour", &json!({"field": "displacement", "deformed": true})).unwrap_err().to_string();
        assert!(e.contains("`deformed`") && e.contains("field"), "{e}");
    }

    /// The top-level keys a params doc names (`"key":` and `"key"?:` directly inside its
    /// leading `{…}`; keys of nested objects and of the reply described after it are not).
    fn doc_keys(doc: &str) -> Vec<&str> {
        let mut keys = Vec::new();
        let mut depth = 0usize;
        let mut chars = doc.char_indices();
        while let Some((i, c)) = chars.next() {
            match c {
                '{' | '[' => depth += 1,
                '}' | ']' => {
                    depth = depth.saturating_sub(1);
                    if depth == 0 {
                        break;
                    }
                }
                '"' => {
                    let rest = &doc[i + 1..];
                    let Some(end) = rest.find('"') else { break };
                    let word = &rest[..end];
                    let tail = &rest[end + 1..];
                    if depth == 1 && (tail.starts_with(':') || tail.starts_with("?:")) && !word.is_empty() {
                        keys.push(word);
                    }
                    for _ in 0..word.chars().count() + 1 {
                        chars.next();
                    }
                }
                _ => {}
            }
        }
        keys
    }

    /// A command's params doc names only parameters its form describes, so a documented
    /// parameter is never refused as unknown (#471: `mesh.node.relax` listed `set`).
    #[test]
    fn every_documented_parameter_is_accepted() {
        let mut wrong = Vec::new();
        for spec in crate::commands() {
            let Some(form) = spec.form else { continue };
            let keys = accepted_keys(form, &Value::Null);
            for key in doc_keys(spec.params) {
                if !keys.contains(&key) {
                    wrong.push(format!("{}: `{key}`", spec.id));
                }
            }
        }
        assert!(wrong.is_empty(), "documented but not in the form:\n{}", wrong.join("\n"));
    }
}
