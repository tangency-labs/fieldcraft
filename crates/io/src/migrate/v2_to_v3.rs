//! Format 2 → 3: every stored quantity in SI. **Released with format 3: never edit.** A mistake
//! is fixed by a new step.
//!
//! Formats 1 and 2 stored numbers as they were typed, in the unit system the model declared in
//! `model.json`'s `units` (`mm-t-s-N-MPa`, `m-kg-s-N-Pa`, `mm-kg-ms-kN-GPa`, `mm-g-ms-N-MPa`,
//! `in-lbf-s-psi`, or `unspecified` when none was chosen). Format 3 stores SI (m, kg, s, N, Pa,
//! K, rad) and shows the model in a display system. This step:
//!
//! - **`model.json`**: scales every dimensional number to SI by the declared system, through its
//!   own field map ([`Walk::model`]) and its own factors ([`System`]), which stand for format 2
//!   as released and never follow later changes to the model's types or the unit table. A
//!   number the map does not know refuses the file rather than pass unscaled. `units` becomes
//!   `display_units`, the declared system, so every number a user sees stays the same.
//!   - `unspecified` is read as mm-t-s-N-MPa: a model's geometry was always in millimetres, and
//!     every default and example was in that system.
//!   - Absolute temperatures were in °C in the millimetre systems, K in m-kg-s-N-Pa and °F in
//!     in-lbf-s-psi (the systems' own temperature units; every solver result depended only on
//!     temperature differences, so nothing computed changes). They go to kelvin with the offset;
//!     temperature differences, expansion coefficients and every other per-degree quantity
//!     scale without it.
//!   - Lengths typed against the geometry (datum points, mesh plans, seeds and size controls)
//!     scale from millimetres in a file built on a design or 3MF meshes (SolidCraft geometry is
//!     in millimetres), from the declared unit in one that was not.
//!   - A mesh made from the geometry in a model of another length unit (metres, inches) holds
//!     millimetres (format 2 never scaled it, and solved them as that unit): its nodes migrate
//!     from millimetres, a block's (and its mirror images' and extrusions') from the declared
//!     unit, by body ([`mesh_units`]), and the file gets a note saying so (`Parts::notes`). A
//!     relaxation field follows its mesh. A node of both kinds refuses the file.
//!   - Defaults the model read for absent dimensional fields (a reference temperature of 0, a
//!     plane strain thickness and a Johnson–Cook reference strain rate of 1) are written out,
//!     converted, since a default means the same number in SI no longer.
//! - **`mesh.bin`**: `MCMESH01` → `MCMESH02`, the same layout with node coordinates scaled from
//!   the declared length unit (millimetres for a mesh of the geometry) to metres ([`mesh`]; its
//!   own reader of the `MCMESH01` layout).
//! - **`journal.json`**: its numbers are not touched. The list of commands gets a header naming
//!   the system its bare numbers were typed in (`units`, the declared system) and the system a
//!   model started in before the first command (`model_units`: `unspecified`, read as
//!   mm-t-s-N-MPa), so replaying it gives the model this step makes:
//!   `{"units": "mm-t-s-N-MPa", "model_units": "mm-t-s-N-MPa", "commands": [...]}`.
//!
//! Scaling by a power of ten works on the decimal a number was written as ([`decimal`]), as the
//! unit table of the release that made format 3 did, so a typed 3.97 mm is the double nearest
//! 0.00397 m.

use std::collections::HashSet;

use serde_json::{Map, Value};

use super::{Parts, set_format};
use crate::{Error, Result, invalid};

/// What a format 2 number is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Q {
    Length,
    /// A length the commands took in millimetres whatever the system.
    Millimetre,
    Area,
    InverseLength,
    Mass,
    RotaryInertia,
    Density,
    Time,
    Frequency,
    AngularVelocity,
    StrainRate,
    Velocity,
    InverseVelocity,
    Acceleration,
    VolumeFlow,
    Force,
    Moment,
    Pressure,
    ForcePerVolume,
    Stiffness,
    RotationalStiffness,
    Damping,
    RotationalDamping,
    DynamicViscosity,
    Temperature,
    ThermalExpansion,
    ThermalConductivity,
    SpecificHeat,
    HeatTransferCoefficient,
    HeatFlux,
    HeatGeneration,
}

/// The exponents of length, mass, time and temperature of a quantity.
fn dimension(q: Q) -> [i32; 4] {
    match q {
        Q::Length | Q::Millimetre => [1, 0, 0, 0],
        Q::Area => [2, 0, 0, 0],
        Q::InverseLength => [-1, 0, 0, 0],
        Q::Mass => [0, 1, 0, 0],
        Q::RotaryInertia => [2, 1, 0, 0],
        Q::Density => [-3, 1, 0, 0],
        Q::Time => [0, 0, 1, 0],
        Q::Frequency | Q::AngularVelocity | Q::StrainRate => [0, 0, -1, 0],
        Q::Velocity => [1, 0, -1, 0],
        Q::InverseVelocity => [-1, 0, 1, 0],
        Q::Acceleration => [1, 0, -2, 0],
        Q::VolumeFlow => [3, 0, -1, 0],
        Q::Force => [1, 1, -2, 0],
        Q::Moment | Q::RotationalStiffness => [2, 1, -2, 0],
        Q::Pressure => [-1, 1, -2, 0],
        Q::ForcePerVolume => [-2, 1, -2, 0],
        Q::Stiffness => [0, 1, -2, 0],
        Q::Damping => [0, 1, -1, 0],
        Q::RotationalDamping => [2, 1, -1, 0],
        Q::DynamicViscosity => [-1, 1, -1, 0],
        Q::Temperature => [0, 0, 0, 1],
        Q::ThermalExpansion => [0, 0, 0, -1],
        Q::ThermalConductivity => [1, 1, -3, -1],
        Q::SpecificHeat => [2, 0, -2, -1],
        Q::HeatTransferCoefficient => [0, 1, -3, -1],
        Q::HeatFlux => [0, 1, -3, 0],
        Q::HeatGeneration => [-1, 1, -3, 0],
    }
}

/// A unit system of format 2, as a conversion to SI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum System {
    MKgS,
    MmTonneS,
    MmKgMs,
    MmGMs,
    InLbfS,
}

/// A factor to SI, `mantissa × 10^exp10`, and the SI value of the unit's zero.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Factor {
    mantissa: f64,
    exp10: i32,
    offset: f64,
}

/// One degree Fahrenheit in kelvin, and Fahrenheit zero.
const FAHRENHEIT_STEP: f64 = 5.0 / 9.0;
const FAHRENHEIT_ZERO: f64 = 459.67 * FAHRENHEIT_STEP;
const CELSIUS_ZERO: f64 = 273.15;

impl System {
    fn of(label: &str) -> Option<System> {
        Some(match label {
            "m-kg-s-N-Pa" => System::MKgS,
            "mm-t-s-N-MPa" | "unspecified" => System::MmTonneS,
            "mm-kg-ms-kN-GPa" => System::MmKgMs,
            "mm-g-ms-N-MPa" => System::MmGMs,
            "in-lbf-s-psi" => System::InLbfS,
            _ => return None,
        })
    }

    /// The label format 3 shows the model in (`unspecified` becomes mm-t-s-N-MPa).
    fn label(self) -> &'static str {
        match self {
            System::MKgS => "m-kg-s-N-Pa",
            System::MmTonneS => "mm-t-s-N-MPa",
            System::MmKgMs => "mm-kg-ms-kN-GPa",
            System::MmGMs => "mm-g-ms-N-MPa",
            System::InLbfS => "in-lbf-s-psi",
        }
    }

    /// The zero of the system's temperature unit in kelvin, when it is not 0 °C (format 3's
    /// default for a fluid's `temperature_zero`): 0 K for m-kg-s-N-Pa, 0 °F for in-lbf-s-psi.
    fn temperature_zero(self) -> Option<f64> {
        match self {
            System::MKgS => Some(0.0),
            System::InLbfS => Some(FAHRENHEIT_ZERO),
            System::MmTonneS | System::MmKgMs | System::MmGMs => None,
        }
    }

    /// The factor of quantity `q` in this system.
    fn factor(self, q: Q) -> Factor {
        // The metric systems' base units are powers of ten of SI's: (length, mass, time).
        let base = match self {
            System::MKgS => [0, 0, 0],
            System::MmTonneS => [-3, 3, 0],
            System::MmKgMs => [-3, 0, -3],
            System::MmGMs => [-3, -3, -3],
            System::InLbfS => return inch_pound(q),
        };
        let [l, m, t, _] = dimension(q);
        let exp10 = if q == Q::Millimetre { -3 } else { l * base[0] + m * base[1] + t * base[2] };
        // Absolute temperatures: kelvin in m-kg-s-N-Pa, °C in the millimetre systems.
        let offset = if q == Q::Temperature && self != System::MKgS { CELSIUS_ZERO } else { 0.0 };
        Factor { mantissa: 1.0, exp10, offset }
    }
}

/// The in-lbf-s-psi factors (inch, lbf·s²/in, second, °F), as the unit table of the release
/// that made format 3 computed them.
fn inch_pound(q: Q) -> Factor {
    let mantissa = match q {
        Q::Millimetre => return Factor { mantissa: 1.0, exp10: -3, offset: 0.0 },
        Q::Length | Q::Velocity | Q::Acceleration => 0.0254,
        Q::Area => 0.00064516,
        Q::InverseLength | Q::InverseVelocity => 39.37007874015748,
        Q::Mass | Q::Stiffness | Q::Damping | Q::HeatFlux => 175.12683524647636,
        Q::RotaryInertia | Q::Moment | Q::RotationalStiffness | Q::RotationalDamping => 0.11298482902761668,
        Q::Density => 10686895.178201316,
        Q::Time | Q::Frequency | Q::AngularVelocity | Q::StrainRate => 1.0,
        Q::VolumeFlow => 1.6387064e-5,
        Q::Force => 4.4482216152605,
        Q::Pressure | Q::DynamicViscosity | Q::HeatGeneration => 6894.757293168361,
        Q::ForcePerVolume => 271447.1375263134,
        Q::Temperature => FAHRENHEIT_STEP,
        Q::ThermalExpansion => 1.7999999999999998,
        Q::ThermalConductivity => 8.006798907468898,
        Q::SpecificHeat => 0.001161288,
        Q::HeatTransferCoefficient => 315.22830344365747,
    };
    let offset = if q == Q::Temperature { FAHRENHEIT_ZERO } else { 0.0 };
    Factor { mantissa, exp10: 0, offset }
}

impl Factor {
    fn to_si(self, v: f64) -> f64 {
        // An absolute °F: (F + 459.67) × 5/9 on the decimal F was written as, one correctly
        // rounded division of exact integers (68 °F is 293.15 K).
        if self.offset == FAHRENHEIT_ZERO {
            return decimal::fahrenheit_to_kelvin(v);
        }
        let scaled = if self.mantissa == 1.0 { v } else { v * self.mantissa };
        let scaled = decimal::shift(scaled, self.exp10);
        if self.offset != 0.0 { decimal::add(scaled, self.offset) } else { scaled }
    }
}

/// Run the step.
pub(super) fn v2_to_v3(parts: &mut Parts) -> Result<()> {
    let Some(model) = parts.model.as_object() else { return invalid("`model.json` is not a JSON object") };
    let label = match model.get("units") {
        None => "unspecified".to_string(),
        Some(Value::String(s)) => s.clone(),
        Some(_) => return invalid("`model.json`: `units` is not a unit system label"),
    };
    let system = System::of(&label)
        .ok_or_else(|| Error::Invalid(format!("`model.json`: unknown unit system `{}`", label.chars().take(64).collect::<String>())))?;
    // Meshes of the geometry are millimetres whatever was declared (read before anything changes).
    let units = mesh_units(system, parts)?;
    // Lengths the commands took against the geometry (datums, plans, seeds, size controls):
    // millimetres in a file built on a design or 3MF meshes, else the declared unit.
    let geometry = parts.manifest.get("geometry");
    let built_on_geometry = geometry.and_then(|g| g.get("design")).is_some_and(|d| !d.is_null())
        || geometry.and_then(|g| g.get("meshes")).and_then(Value::as_array).is_some_and(|m| !m.is_empty());
    let geometry_length = if built_on_geometry { Q::Millimetre } else { Q::Length };
    let Some(model) = parts.model.as_object_mut() else { return invalid("`model.json` is not a JSON object") };
    model.remove("units");
    model.insert("display_units".into(), Value::String(system.label().into()));
    let mut walk = Walk { system, handled: HashSet::new(), geometry_length, units };
    walk.model(&mut parts.model)?;
    walk.check_all_handled(&parts.model)?;
    parts.mesh = mesh(&parts.mesh, system.factor(Q::Length), system.factor(Q::Millimetre), &walk.units.millimetre_nodes)?;
    parts.notes.extend(walk.units.note.take());
    let commands = std::mem::take(&mut parts.journal);
    if !commands.is_array() {
        return invalid("`journal.json` is not a list of commands");
    }
    let mut journal = Map::new();
    journal.insert("units".into(), Value::String(system.label().into()));
    journal.insert("model_units".into(), Value::String(System::MmTonneS.label().into()));
    journal.insert("commands".into(), commands);
    parts.journal = Value::Object(journal);
    set_format(parts, 3)
}

/// Where each mesh body of a format 2 file came from, and so in which unit its node coordinates
/// are: [`Origin::Typed`] in the declared system, [`Origin::Geometry`] in millimetres.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Origin {
    /// Made from numbers typed in the declared system: `mesh.block`, elements tagged onto a body
    /// (`mesh.tag`), and what was mirrored or extruded from those.
    Typed,
    /// Made from the geometry (a SolidCraft design or a 3MF file), always in millimetres.
    Geometry,
}

/// The meshers of format 2 that mesh the geometry, and the parameter naming their bodies.
const GEOMETRY_MESHERS: [(&str, &str); 10] = [
    ("mesh.tet.generate", "body"),
    ("mesh.hex.generate", "body"),
    ("mesh.hex.blocks", "bodies"),
    ("mesh.hex.sweep", "body"),
    ("mesh.hex.mapped", "body"),
    ("mesh.hex.octree", "bodies"),
    ("mesh.hex.project", "body"),
    ("mesh.planar", "body"),
    ("mesh.cells", "bodies"),
    ("mesh.plan.apply", "body"),
];

/// The names a parameter holds: a name or a list of them.
fn names(v: Option<&Value>) -> Vec<String> {
    match v {
        Some(Value::String(s)) => vec![s.clone()],
        Some(Value::Array(a)) => a.iter().filter_map(Value::as_str).map(str::to_string).collect(),
        _ => Vec::new(),
    }
}

/// The origin of each body the journal names, traced through the commands in order: blocks and
/// tagged elements are typed, the geometry meshers' bodies geometry, and a mirror image or an
/// extrusion takes its source's origin. Bounded by the journal's size.
fn origins_from_journal(journal: &Value) -> std::collections::BTreeMap<String, Origin> {
    let mut origin = std::collections::BTreeMap::new();
    for c in journal.as_array().into_iter().flatten() {
        let (Some(id), p) = (c.get("command").and_then(Value::as_str), c.get("params").unwrap_or(&Value::Null)) else { continue };
        match id {
            "mesh.block" => {
                for b in names(p.get("body")) {
                    origin.insert(b, Origin::Typed);
                }
            }
            "mesh.tag" => {
                if let Some(b) = p.get("geometry").and_then(|g| g.get("body")).and_then(Value::as_str) {
                    origin.entry(b.to_string()).or_insert(Origin::Typed);
                }
            }
            "mesh.mirror" => {
                let sources = match p.get("bodies") {
                    Some(_) => names(p.get("bodies")),
                    None => origin.keys().cloned().collect(),
                };
                for s in sources {
                    let Some(o) = origin.get(&s).copied() else { continue };
                    let image = p.get("as").and_then(|a| a.get(&s)).and_then(Value::as_str).map_or_else(|| format!("{s}/mirror"), str::to_string);
                    origin.insert(image, o);
                }
            }
            "mesh.hex.extrude" => {
                let from: Vec<String> = p
                    .get("source")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(|f| f.get("body")?.as_str().map(str::to_string))
                    .collect();
                let Some(out) = p.get("body").and_then(Value::as_str).map(str::to_string).or_else(|| from.first().cloned()) else { continue };
                // From a typed section the layers are typed; from the geometry's, millimetres.
                let o = if from.iter().any(|b| origin.get(b) == Some(&Origin::Typed)) { Origin::Typed } else { Origin::Geometry };
                origin.insert(out, o);
            }
            _ => {
                if let Some((_, key)) = GEOMETRY_MESHERS.iter().find(|(m, _)| *m == id) {
                    for b in names(p.get(*key)) {
                        origin.insert(b, Origin::Geometry);
                    }
                }
            }
        }
    }
    origin
}

/// The bodies a format 2 journal (`[{command, params}]`) made from numbers typed in its own
/// system (blocks, tagged elements, and what was mirrored or extruded from them): with them, a
/// replay can tell the geometry and mesher commands that acted on them.
pub fn format_2_typed_bodies(journal: &Value) -> Vec<String> {
    origins_from_journal(journal).into_iter().filter(|(_, o)| *o == Origin::Typed).map(|(b, _)| b).collect()
}

/// How the step reads a format 2 mesh: the nodes in millimetres, the coordinates (bits) of the
/// nodes of each origin and the box round the millimetre ones (for relaxation fields), and a
/// note for the user.
#[derive(Default)]
struct MeshUnits {
    millimetre_nodes: HashSet<u32>,
    millimetre_points: HashSet<[u64; 3]>,
    typed_points: HashSet<[u64; 3]>,
    millimetre_box: Option<([f64; 3], [f64; 3])>,
    /// The origin of each body.
    origins: std::collections::BTreeMap<String, Origin>,
    note: Option<String>,
}

/// The unit of each mesh body of a file that declared metres or inches and was built on a design
/// or 3MF meshes: format 2 never scaled a mesh of the geometry (its numbers are millimetres,
/// solved as the declared unit, so the part analysed was 1000× or 25.4× the design), while a
/// block's numbers are the declared unit. Each body migrates by its origin
/// ([`origins_from_journal`]; a body the journal does not explain is taken as the geometry's,
/// and the note says so). A node of both kinds (a block joined onto a mesh of the design) has no
/// consistent scale: the file is refused, unchanged.
fn mesh_units(system: System, parts: &Parts) -> Result<MeshUnits> {
    let mut out = MeshUnits::default();
    let (unit, times) = match system {
        System::MKgS => ("metres", "1000"),
        System::InLbfS => ("inches", "25.4"),
        System::MmTonneS | System::MmKgMs | System::MmGMs => return Ok(out),
    };
    let geometry = parts.manifest.get("geometry");
    let design = geometry.and_then(|g| g.get("design")).is_some_and(|d| !d.is_null());
    let meshes = geometry.and_then(|g| g.get("meshes")).and_then(Value::as_array).is_some_and(|m| !m.is_empty());
    if !design && !meshes {
        return Ok(out);
    }
    let traced = origins_from_journal(&parts.journal);
    let layout = MeshLayout::read(&parts.mesh)?;
    // Each element's body (untagged elements were typed: `mesh.element.add`).
    let mut body_of: std::collections::HashMap<u32, &str> = std::collections::HashMap::new();
    for (element, region) in &layout.regions {
        if let Some(b) = layout.bodies.get(*region as usize) {
            body_of.insert(*element, b.as_str());
        }
    }
    let mut inferred: Vec<String> = Vec::new();
    let mut origin_of = |b: &str| -> Origin {
        match traced.get(b) {
            Some(o) => *o,
            None => {
                if !inferred.iter().any(|i| i == b) {
                    inferred.push(b.to_string());
                }
                Origin::Geometry
            }
        }
    };
    // Each node's bodies, by origin: (a millimetre body, a typed body or untagged element).
    let mut node_bodies: std::collections::HashMap<u32, (Option<String>, Option<String>)> = std::collections::HashMap::new();
    let mut geometry_bodies: Vec<String> = Vec::new();
    for (element, nodes) in &layout.elements {
        let (body, origin) = match body_of.get(element) {
            Some(b) => (b.to_string(), origin_of(b)),
            None => ("(elements of no body)".to_string(), Origin::Typed),
        };
        if origin == Origin::Geometry && !geometry_bodies.contains(&body) {
            geometry_bodies.push(body.clone());
        }
        out.origins.insert(body.clone(), origin);
        for n in nodes {
            let slot = node_bodies.entry(*n).or_default();
            match origin {
                Origin::Geometry => slot.0 = Some(body.clone()),
                Origin::Typed => slot.1 = Some(body.clone()),
            }
        }
    }
    if let Some((n, (Some(g), Some(t)))) = node_bodies.iter().find(|(_, (g, t))| g.is_some() && t.is_some()) {
        return invalid(format!(
            "this model declares the unit system {label}, and its node {n} joins the mesh of body `{g}`, made from a SolidCraft design or 3MF file in millimetres, \
             to `{t}`, typed in {unit}: FieldCraft 0.7 and earlier joined the two as if they had one scale, so there is no consistent one to read them in. \
             Nothing was converted and the file is unchanged; it still opens in FieldCraft 0.7.x.",
            label = system.label(),
            g = g.chars().take(64).collect::<String>(),
            t = t.chars().take(64).collect::<String>(),
        ));
    }
    if geometry_bodies.is_empty() {
        return Ok(out);
    }
    for (n, (g, _)) in &node_bodies {
        if g.is_some() {
            out.millimetre_nodes.insert(*n);
        }
    }
    for (label, xyz) in &layout.nodes {
        let bits = xyz.map(f64::to_bits);
        if out.millimetre_nodes.contains(label) {
            out.millimetre_points.insert(bits);
            let (lo, hi) = out.millimetre_box.get_or_insert((*xyz, *xyz));
            for k in 0..3 {
                lo[k] = lo[k].min(xyz[k]);
                hi[k] = hi[k].max(xyz[k]);
            }
        } else {
            out.typed_points.insert(bits);
        }
    }
    let quoted =
        |list: &[String]| list.iter().take(8).map(|b| format!("`{}`", b.chars().take(64).collect::<String>())).collect::<Vec<_>>().join(", ");
    let source = if design { "a SolidCraft design" } else { "a 3MF file" };
    let mut note = format!(
        "This model declared {label}, but the mesh of {which} {names} was made from {source}, which is in millimetres: \
         FieldCraft 0.7 and earlier solved those millimetre coordinates as {unit}, so the part it analysed was {times}× the design. \
         The mesh is now read in millimetres and overlays the design again (block meshes typed in {unit} stay as they were); \
         results will differ from the old ones. The first save keeps the old file beside the new one.",
        label = system.label(),
        which = if geometry_bodies.len() == 1 { "body" } else { "bodies" },
        names = quoted(&geometry_bodies),
    );
    if !inferred.is_empty() {
        note.push_str(&format!(
            " No command in the journal made {} ({}); {} read as made from the {}.",
            if inferred.len() == 1 { "body" } else { "bodies" },
            quoted(&inferred),
            if inferred.len() == 1 { "it is" } else { "they are" },
            if design { "design" } else { "3MF file" },
        ));
    }
    out.note = Some(note);
    Ok(out)
}

/// Most entries `MeshLayout::read` takes from a hostile count before checking it: no more
/// records than the bytes left can hold.
fn records(count: u64, left: usize, smallest: usize) -> Result<usize> {
    usize::try_from(count)
        .ok()
        .filter(|n| n.checked_mul(smallest).is_some_and(|b| b <= left))
        .ok_or_else(|| Error::Invalid("`mesh.bin` lists more records than it holds".into()))
}

/// What the migration needs of an `MCMESH01` part (format 2's layout): node coordinates, each
/// element's nodes, the body name of each entity, and each tagged element's region (an entity).
struct MeshLayout {
    nodes: Vec<(u32, [f64; 3])>,
    elements: Vec<(u32, Vec<u32>)>,
    bodies: Vec<String>,
    regions: Vec<(u32, u32)>,
}

impl MeshLayout {
    fn read(bytes: &[u8]) -> Result<MeshLayout> {
        const HEADER: usize = 8 + 5 * 8;
        const NODE: usize = 4 + 3 * 8;
        let truncated = || Error::Invalid("`mesh.bin` is truncated".into());
        let u32_at = |at: usize| bytes.get(at..at + 4).and_then(|b| <[u8; 4]>::try_from(b).ok()).map(u32::from_le_bytes).ok_or_else(truncated);
        let u64_at = |at: usize| bytes.get(at..at + 8).and_then(|b| <[u8; 8]>::try_from(b).ok()).map(u64::from_le_bytes).ok_or_else(truncated);
        if bytes.get(..8) != Some(b"MCMESH01".as_slice()) {
            return invalid("`mesh.bin` of format 2 must start with MCMESH01");
        }
        let left = |at: usize| bytes.len().saturating_sub(at);
        let node_count = records(u64_at(8)?, left(HEADER), NODE)?;
        let mut nodes = Vec::with_capacity(node_count);
        let mut at = HEADER;
        for _ in 0..node_count {
            let label = u32_at(at)?;
            let c = |k: usize| u64_at(at + 4 + 8 * k).map(f64::from_bits);
            nodes.push((label, [c(0)?, c(1)?, c(2)?]));
            at += NODE;
        }
        // The smallest records: an element of two nodes, a body's entity, an untagged element.
        let element_count = records(u64_at(16)?, left(at), 5 + 8)?;
        let mut elements = Vec::with_capacity(element_count);
        for _ in 0..element_count {
            let label = u32_at(at)?;
            let kind = *bytes.get(at + 4).ok_or_else(truncated)?;
            // Nodes per element kind (format 2's codes).
            let n: usize = match kind {
                1 | 2 => 2,
                3 | 12 => 3,
                4 | 9 => 6,
                5 | 7 => 4,
                6 | 10 => 8,
                8 => 10,
                11 => 20,
                13 => 5,
                _ => return invalid(format!("`mesh.bin`: unknown element kind {kind}")),
            };
            let of = (0..n).map(|k| u32_at(at + 5 + 4 * k)).collect::<Result<Vec<u32>>>()?;
            elements.push((label, of));
            at = at.checked_add(5 + 4 * n).filter(|a| *a <= bytes.len()).ok_or_else(truncated)?;
        }
        let name = |at: usize| -> Result<(String, usize)> {
            let len = bytes.get(at..at + 2).and_then(|b| <[u8; 2]>::try_from(b).ok()).map(u16::from_le_bytes).ok_or_else(truncated)? as usize;
            let text = bytes.get(at + 2..at + 2 + len).ok_or_else(truncated)?;
            Ok((String::from_utf8_lossy(text).into_owned(), at + 2 + len))
        };
        let entity_count = records(u64_at(24)?, left(at), 3)?;
        let mut bodies = Vec::with_capacity(entity_count);
        for _ in 0..entity_count {
            let kind = *bytes.get(at).ok_or_else(truncated)?;
            let (body, next) = name(at + 1)?;
            at = if kind == 0 { next } else { name(next)?.1 };
            bodies.push(body);
        }
        let tagged = records(u64_at(32)?, left(at), 10)?;
        let mut regions = Vec::with_capacity(tagged);
        for _ in 0..tagged {
            let (element, region) = (u32_at(at)?, u32_at(at + 4)?);
            if region != u32::MAX {
                regions.push((element, region));
            }
            at += 8;
            for size in [5usize, 5] {
                let n = *bytes.get(at).ok_or_else(truncated)? as usize;
                at = at.checked_add(1 + n * size).filter(|a| *a <= bytes.len()).ok_or_else(truncated)?;
            }
        }
        Ok(MeshLayout { nodes, elements, bodies, regions })
    }
}

/// The defaults of format 2 that format 3 changed, as a command of a journal written before
/// format 3 meant them when it left them out: `params` of `command` get them written in where
/// they are missing, in the journal's own units (a bare number), so a replay builds the model
/// the step built (which writes the same defaults into `model.json`). Frozen with the step:
/// - Johnson–Cook plasticity's `reference_strain_rate` (under `plastic`): 1 (per the system's
///   time unit);
/// - a material's `thermal.reference_temperature` and a fluid domain's
///   `thermal.initial_temperature`: 0 (in the system's temperature unit);
/// - a transient heat step's `transient.initial_temperature`: 0;
/// - a plane-strain section's `thickness`: 1 (the system's length unit);
/// - a fluid's `thermal.temperature_zero`, which format 3 added: the zero of the system's
///   temperature unit in kelvin (the step writes the same), outside the millimetre systems.
///
/// Each default goes into an object the command gives, never into one it leaves out: that is
/// when format 2 read the default.
pub fn format_2_defaults(command: &str, params: &mut Value, system: &str) {
    let material = command.starts_with("model.material.");
    let fluid = command.starts_with("model.fluid_domain.");
    // A fluid's heat counted from the zero of the system the journal was typed in (in kelvin;
    // a field no form shows, so it is passed on as it is).
    let zero = System::of(system).and_then(System::temperature_zero);
    legacy_defaults(params, "", material, fluid, zero, 0);
}

/// Most nesting [`format_2_defaults`] follows.
const LEGACY_DEPTH: usize = 12;

fn legacy_defaults(v: &mut Value, key: &str, material: bool, fluid: bool, zero: Option<f64>, depth: usize) {
    if depth > LEGACY_DEPTH {
        return;
    }
    match v {
        Value::Object(o) => {
            let default = |o: &mut Map<String, Value>, field: &str, value: f64| {
                if !o.contains_key(field)
                    && let Some(n) = serde_json::Number::from_f64(value)
                {
                    o.insert(field.into(), Value::Number(n));
                }
            };
            match (key, o.get("type").and_then(Value::as_str)) {
                ("plastic", Some("johnson_cook")) => default(o, "reference_strain_rate", 1.0),
                (_, Some("plane_strain")) => default(o, "thickness", 1.0),
                _ => {}
            }
            match key {
                "thermal" if material => default(o, "reference_temperature", 0.0),
                "thermal" if fluid => {
                    default(o, "initial_temperature", 0.0);
                    if let Some(z) = zero {
                        default(o, "temperature_zero", z);
                    }
                }
                "transient" => default(o, "initial_temperature", 0.0),
                _ => {}
            }
            for (k, x) in o.iter_mut() {
                legacy_defaults(x, k, material, fluid, zero, depth + 1);
            }
        }
        Value::Array(items) => items.iter_mut().for_each(|x| legacy_defaults(x, key, material, fluid, zero, depth + 1)),
        _ => {}
    }
}

/// Walks `model.json`, converting each number it knows and recording the JSON pointer of every
/// value it has dealt with (a number, or a whole list or subtree of them, converted or
/// dimensionless): one entry per field, never one per number of a long list.
struct Walk {
    system: System,
    handled: HashSet<String>,
    /// The unit of lengths typed against the geometry (datums, plans, seeds, size controls).
    geometry_length: Q,
    /// The mesh's units by body ([`mesh_units`]).
    units: MeshUnits,
}

/// Whether `v` is a number, a null or a list (nested) of them only.
fn numeric(v: &Value) -> bool {
    match v {
        Value::Number(_) | Value::Null => true,
        Value::Array(items) => items.iter().all(numeric),
        _ => false,
    }
}

/// A JSON pointer token (RFC 6901).
fn token(key: &str) -> String {
    key.replace('~', "~0").replace('/', "~1")
}

impl Walk {
    fn convert(&self, v: f64, q: Q, at: &str) -> Result<Value> {
        let si = self.system.factor(q).to_si(v);
        serde_json::Number::from_f64(si)
            .map(Value::Number)
            .ok_or_else(|| Error::Invalid(format!("`model.json`: the number at {at} is out of range in SI")))
    }

    /// The value at `at`: a number, a list of them (nested) or null, of quantity `q` (`None`:
    /// dimensionless, a count or a label).
    fn value(&mut self, v: &mut Value, at: &str, q: Option<Q>) -> Result<()> {
        match v {
            Value::Number(n) => {
                if let Some(q) = q {
                    let x = n.as_f64().ok_or_else(|| Error::Invalid(format!("`model.json`: the number at {at} is not finite")))?;
                    *v = self.convert(x, q, at)?;
                }
                self.handled.insert(at.to_string());
            }
            // A list of numbers (nested, with nulls) is handled whole: one entry, not one per
            // number. Anything else in it is walked item by item and left for the check.
            Value::Array(_) if numeric(v) => {
                self.convert_all(v, at, q)?;
                self.handled.insert(at.to_string());
            }
            Value::Array(items) => {
                for (i, item) in items.iter_mut().enumerate() {
                    self.value(item, &format!("{at}/{i}"), q)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// Every number under `v`, a number or nested lists of numbers and nulls, converted as `q`.
    fn convert_all(&self, v: &mut Value, at: &str, q: Option<Q>) -> Result<()> {
        match v {
            Value::Number(n) => {
                if let Some(q) = q {
                    let x = n.as_f64().ok_or_else(|| Error::Invalid(format!("`model.json`: the number at {at} is not finite")))?;
                    *v = self.convert(x, q, at)?;
                }
                Ok(())
            }
            Value::Array(items) => items.iter_mut().try_for_each(|item| self.convert_all(item, at, q)),
            _ => Ok(()),
        }
    }

    /// The unit of a relaxation field's cell (of the field made on `bodies`): that of the mesh
    /// it was made on. Its moved corners are nodes of that mesh unless it was remeshed since;
    /// then the field's bodies say, and failing those, whether the cell lies among the
    /// millimetre nodes.
    fn relaxation_unit(&self, cell: &Value, bodies: &[String]) -> Q {
        let u = &self.units;
        if u.millimetre_points.is_empty() {
            return Q::Length;
        }
        let points: Vec<[f64; 3]> = ["moved", "at"]
            .iter()
            .flat_map(|k| cell.get(*k).and_then(Value::as_array).into_iter().flatten())
            .filter_map(|p| <[f64; 3]>::try_from(p.as_array()?.iter().filter_map(Value::as_f64).collect::<Vec<_>>()).ok())
            .collect();
        let bits = |p: &[f64; 3]| p.map(f64::to_bits);
        if points.iter().any(|p| u.millimetre_points.contains(&bits(p))) {
            return Q::Millimetre;
        }
        if points.iter().any(|p| u.typed_points.contains(&bits(p))) {
            return Q::Length;
        }
        let of: Vec<Origin> = bodies.iter().filter_map(|b| u.origins.get(b).copied()).collect();
        if !of.is_empty() {
            return if of.iter().all(|o| *o == Origin::Geometry) { Q::Millimetre } else { Q::Length };
        }
        match u.millimetre_box {
            Some((lo, hi)) => {
                let margin = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max) * 0.01;
                let inside = |p: &[f64; 3]| (0..3).all(|k| p[k] >= lo[k] - margin && p[k] <= hi[k] + margin);
                if !points.is_empty() && points.iter().all(inside) { Q::Millimetre } else { Q::Length }
            }
            None => Q::Length,
        }
    }

    /// Field `key` of `parent` (at `at`), when present.
    fn field(&mut self, parent: &mut Value, at: &str, key: &str, q: Option<Q>) -> Result<()> {
        match parent.get_mut(key) {
            Some(v) => self.value(v, &format!("{at}/{}", token(key)), q),
            None => Ok(()),
        }
    }

    /// Field `key` of quantity `q`, written out as `default` (in the declared system) when
    /// absent: a default that meant this number in the declared system.
    fn field_or(&mut self, parent: &mut Value, at: &str, key: &str, q: Q, default: f64) -> Result<()> {
        if let Some(obj) = parent.as_object_mut()
            && !obj.contains_key(key)
            && let Some(n) = serde_json::Number::from_f64(default)
        {
            obj.insert(key.into(), Value::Number(n));
        }
        self.field(parent, at, key, Some(q))
    }

    /// Every number under `parent`'s `key` is dimensionless (a count, a label, a ratio).
    fn plain(&mut self, parent: &mut Value, at: &str, key: &str) {
        if parent.get(key).is_some() {
            self.handled.insert(format!("{at}/{}", token(key)));
        }
    }

    /// The members of a map of definitions by name, each walked by `f`.
    fn each_named(&mut self, model: &mut Value, key: &str, mut f: impl FnMut(&mut Self, &mut Value, &str) -> Result<()>) -> Result<()> {
        if let Some(Value::Object(items)) = model.get_mut(key) {
            for (name, item) in items.iter_mut() {
                f(self, item, &format!("/{}/{}", token(key), token(name)))?;
            }
        }
        Ok(())
    }

    /// The items of a list, each walked by `f`.
    fn each(&mut self, parent: &mut Value, at: &str, key: &str, mut f: impl FnMut(&mut Self, &mut Value, &str) -> Result<()>) -> Result<()> {
        if let Some(Value::Array(items)) = parent.get_mut(key) {
            for (i, item) in items.iter_mut().enumerate() {
                f(self, item, &format!("{at}/{}/{i}", token(key)))?;
            }
        }
        Ok(())
    }

    fn model(&mut self, m: &mut Value) -> Result<()> {
        self.each_named(m, "sets", |w, set, at| {
            for key in ["nodes", "elements", "faces", "edges"] {
                w.plain(set, at, key);
            }
            Ok(())
        })?;
        self.each_named(m, "materials", Walk::material)?;
        self.each_named(m, "sections", Walk::section)?;
        self.each_named(m, "reference_points", |w, rp, at| {
            w.field(rp, at, "xyz", Some(Q::Length))?;
            w.plain(rp, at, "label");
            Ok(())
        })?;
        self.each_named(m, "couplings", |w, c, at| {
            w.plain(c, at, "dofs");
            Ok(())
        })?;
        self.each_named(m, "contacts", |w, c, at| {
            w.friction(c, at)?;
            w.plain(c, at, "penalty_scale");
            Ok(())
        })?;
        self.each_named(m, "ties", |w, t, at| w.field(t, at, "position_tolerance", Some(Q::Length)))?;
        self.each_named(m, "masses", |w, p, at| {
            w.field(p, at, "mass", Some(Q::Mass))?;
            w.field(p, at, "inertia", Some(Q::RotaryInertia))
        })?;
        self.each_named(m, "springs", Walk::spring)?;
        self.each_named(m, "boundary_conditions", Walk::boundary_condition)?;
        self.each_named(m, "loads", Walk::load)?;
        self.each_named(m, "fluid_domains", Walk::fluid_domain)?;
        self.each(m, "", "steps", Walk::step)?;
        self.each_named(m, "datums", |w, d, at| {
            w.field(d, at, "point", Some(w.geometry_length))?;
            w.plain(d, at, "normal");
            w.plain(d, at, "direction");
            Ok(())
        })?;
        self.each(m, "", "seeds", |w, s, at| {
            w.plain(s, at, "count");
            // Seeds and size controls were read by the meshers, on the geometry: millimetres.
            w.field(s, at, "size", Some(w.geometry_length))?;
            w.plain(s, at, "bias");
            Ok(())
        })?;
        self.each_named(m, "size_controls", |w, c, at| {
            w.field(c, at, "size", Some(w.geometry_length))?;
            // A region is a box, a sphere, or geometry (bodies, faces: names).
            if let Some(region) = c.get_mut("region") {
                let at = format!("{at}/region");
                if let Some(b) = region.get_mut("box") {
                    let at = format!("{at}/box");
                    w.field(b, &at, "min", Some(w.geometry_length))?;
                    w.field(b, &at, "max", Some(w.geometry_length))?;
                }
                if let Some(sphere) = region.get_mut("sphere") {
                    let at = format!("{at}/sphere");
                    w.field(sphere, &at, "center", Some(w.geometry_length))?;
                    w.field(sphere, &at, "radius", Some(w.geometry_length))?;
                }
            }
            Ok(())
        })?;
        self.plain(m, "", "size_growth");
        self.each_named(m, "relaxations", |w, r, at| {
            // A cell is in the coordinates of the mesh it was made on: millimetres for a mesh of
            // the geometry (its moved corners are nodes of it), else the declared length.
            let bodies = names(r.get("bodies"));
            w.each(r, at, "cells", |w, cell, at| {
                let q = w.relaxation_unit(cell, &bodies);
                w.field(cell, at, "at", Some(q))?;
                w.field(cell, at, "moved", Some(q))
            })
        })?;
        self.each_named(m, "mesh_plans", |w, p, at| {
            w.field(p, at, "size", Some(w.geometry_length))?;
            w.field(p, at, "tolerance", Some(w.geometry_length))?;
            for key in ["layers", "direction", "directions", "features", "unmapped"] {
                w.plain(p, at, key);
            }
            w.each(p, at, "cuts", |w, cut, at| {
                w.field(cut, at, "point", Some(w.geometry_length))?;
                w.plain(cut, at, "normal");
                w.plain(cut, at, "feature");
                Ok(())
            })?;
            w.each(p, at, "rings", |w, ring, at| {
                for key in ["center", "radius", "frame"] {
                    w.field(ring, at, key, Some(w.geometry_length))?;
                }
                w.plain(ring, at, "axis");
                w.plain(ring, at, "feature");
                Ok(())
            })
        })?;
        self.each(m, "", "stitches", |w, s, at| w.field(s, at, "tolerance", Some(Q::Length)))?;
        self.each(m, "", "defeatures", |w, d, at| {
            w.plain(d, at, "features");
            Ok(())
        })?;
        Ok(())
    }

    fn material(&mut self, m: &mut Value, at: &str) -> Result<()> {
        if let Some(e) = m.get_mut("elastic") {
            let at = format!("{at}/elastic");
            self.field(e, &at, "young", Some(Q::Pressure))?;
            self.field(e, &at, "shear", Some(Q::Pressure))?;
            self.plain(e, &at, "poisson");
        }
        self.field(m, at, "density", Some(Q::Density))?;
        if let Some(f) = m.get_mut("fluid") {
            self.field(f, &format!("{at}/fluid"), "viscosity", Some(Q::DynamicViscosity))?;
        }
        if let Some(h) = m.get_mut("hyperelastic") {
            let at = format!("{at}/hyperelastic");
            for key in ["c01", "c20", "c30"] {
                self.field(h, &at, key, Some(Q::Pressure))?;
            }
        }
        if let Some(p) = m.get_mut("plastic").filter(|p| p.is_object()) {
            let at = format!("{at}/plastic");
            let johnson_cook = p.get("type").and_then(Value::as_str) == Some("johnson_cook");
            for key in ["yield_stress", "hardening_modulus", "a", "b"] {
                self.field(p, &at, key, Some(Q::Pressure))?;
            }
            self.each(p, &at, "points", |w, point, at| {
                w.plain(point, at, "strain");
                w.field(point, at, "stress", Some(Q::Pressure))
            })?;
            for key in ["n", "c", "m", "taylor_quinney"] {
                self.plain(p, &at, key);
            }
            if johnson_cook {
                self.field_or(p, &at, "reference_strain_rate", Q::StrainRate, 1.0)?;
            }
            self.field(p, &at, "room_temperature", Some(Q::Temperature))?;
            self.field(p, &at, "melt_temperature", Some(Q::Temperature))?;
            self.field(p, &at, "specific_heat", Some(Q::SpecificHeat))?;
        }
        // Damage parameters and erosion limits are dimensionless.
        self.plain(m, at, "damage");
        self.plain(m, at, "erosion");
        if let Some(t) = m.get_mut("thermal").filter(|t| t.is_object()) {
            let at = format!("{at}/thermal");
            self.field(t, &at, "conductivity", Some(Q::ThermalConductivity))?;
            self.field(t, &at, "specific_heat", Some(Q::SpecificHeat))?;
            self.field(t, &at, "expansion", Some(Q::ThermalExpansion))?;
            self.field_or(t, &at, "reference_temperature", Q::Temperature, 0.0)?;
        }
        Ok(())
    }

    fn section(&mut self, s: &mut Value, at: &str) -> Result<()> {
        let Some(kind) = s.get_mut("kind") else { return Ok(()) };
        let at = format!("{at}/kind");
        match kind.get("type").and_then(Value::as_str) {
            Some("plane_strain") => self.field_or(kind, &at, "thickness", Q::Length, 1.0)?,
            _ => self.field(kind, &at, "thickness", Some(Q::Length))?,
        }
        self.plain(kind, &at, "offset");
        self.plain(kind, &at, "points");
        self.plain(kind, &at, "orientation");
        self.field(kind, &at, "area", Some(Q::Area))?;
        if let Some(profile) = kind.get_mut("profile") {
            let at = format!("{at}/profile");
            for key in ["width", "height", "radius", "thickness", "web", "flange"] {
                self.field(profile, &at, key, Some(Q::Length))?;
            }
        }
        Ok(())
    }

    fn friction(&mut self, parent: &mut Value, at: &str) -> Result<()> {
        if let Some(f) = parent.get_mut("friction") {
            let at = format!("{at}/friction");
            self.plain(f, &at, "static");
            self.plain(f, &at, "kinetic");
            self.field(f, &at, "decay", Some(Q::InverseVelocity))?;
        }
        Ok(())
    }

    fn spring(&mut self, s: &mut Value, at: &str) -> Result<()> {
        // A spring along degree of freedom 4 to 6 is torsional.
        let torsional = s.get("direction").and_then(|d| d.get("dof")).and_then(Value::as_u64).is_some_and(|d| (4..=6).contains(&d));
        self.field(s, at, "stiffness", Some(if torsional { Q::RotationalStiffness } else { Q::Stiffness }))?;
        self.field(s, at, "coefficient", Some(if torsional { Q::RotationalDamping } else { Q::Damping }))?;
        if let Some(d) = s.get_mut("direction") {
            let at = format!("{at}/direction");
            self.plain(d, &at, "dof");
            self.plain(d, &at, "vector");
        }
        Ok(())
    }

    fn boundary_condition(&mut self, bc: &mut Value, at: &str) -> Result<()> {
        // Displacements: three lengths, then three rotations in radians.
        if let Some(Value::Array(values)) = bc.get_mut("values") {
            for (i, v) in values.iter_mut().enumerate() {
                self.value(v, &format!("{at}/values/{i}"), (i < 3).then_some(Q::Length))?;
            }
        }
        self.plain(bc, at, "plane");
        self.field(bc, at, "value", Some(Q::Temperature))
    }

    fn load(&mut self, load: &mut Value, at: &str) -> Result<()> {
        let kind = load.get("type").and_then(Value::as_str).unwrap_or_default().to_string();
        let vector = match kind.as_str() {
            "force" => Some(Q::Force),
            "moment" => Some(Q::Moment),
            "traction" => Some(Q::Pressure),
            "body_force" => Some(Q::ForcePerVolume),
            _ => None,
        };
        if let Some(q) = vector {
            self.field(load, at, "vector", Some(q))?;
        }
        self.field(load, at, "magnitude", Some(Q::Pressure))?;
        self.field(load, at, "acceleration", Some(Q::Acceleration))?;
        self.field(load, at, "angular_velocity", Some(Q::AngularVelocity))?;
        self.field(load, at, "point", Some(Q::Length))?;
        self.plain(load, at, "axis");
        self.field(load, at, "flux", Some(Q::HeatFlux))?;
        self.field(load, at, "film", Some(Q::HeatTransferCoefficient))?;
        self.field(load, at, "sink", Some(Q::Temperature))?;
        self.field(load, at, "power", Some(Q::HeatGeneration))?;
        self.field(load, at, "value", Some(Q::Temperature))
    }

    /// A fluid face or surface condition: velocities and a pressure.
    fn face_condition(&mut self, c: &mut Value, at: &str) -> Result<()> {
        self.field(c, at, "velocity", Some(Q::Velocity))?;
        self.field(c, at, "pressure", Some(Q::Pressure))
    }

    /// A thermal condition: a held temperature or a heat flux, both as `value`.
    fn thermal_condition(&mut self, c: &mut Value, at: &str) -> Result<()> {
        let q = match c.get("type").and_then(Value::as_str) {
            Some("flux") => Q::HeatFlux,
            _ => Q::Temperature,
        };
        self.field(c, at, "value", Some(q))
    }

    fn fluid_domain(&mut self, d: &mut Value, at: &str) -> Result<()> {
        for key in ["min", "max", "spacing", "tolerance"] {
            self.field(d, at, key, Some(Q::Length))?;
        }
        self.field(d, at, "initial_velocity", Some(Q::Velocity))?;
        self.field(d, at, "acceleration", Some(Q::Acceleration))?;
        if let Some(Value::Object(faces)) = d.get_mut("faces") {
            for (face, c) in faces.iter_mut() {
                self.face_condition(c, &format!("{at}/faces/{}", token(face)))?;
            }
        }
        self.each(d, at, "surfaces", |w, s, at| match s.get_mut("condition") {
            Some(c) => w.face_condition(c, &format!("{at}/condition")),
            None => Ok(()),
        })?;
        if let Some(t) = d.get_mut("thermal").filter(|t| t.is_object()) {
            let at = format!("{at}/thermal");
            self.field_or(t, &at, "initial_temperature", Q::Temperature, 0.0)?;
            // The lattice ran on the numbers as typed: it keeps counting from the system's zero.
            if let Some(zero) = self.system.temperature_zero()
                && let (Some(obj), Some(n)) = (t.as_object_mut(), serde_json::Number::from_f64(zero))
            {
                obj.insert("temperature_zero".into(), Value::Number(n));
                self.handled.insert(format!("{at}/temperature_zero"));
            }
            self.field(t, &at, "gravity", Some(Q::Acceleration))?;
            self.field(t, &at, "expansion", Some(Q::ThermalExpansion))?;
            if let Some(Value::Object(faces)) = t.get_mut("faces") {
                for (face, c) in faces.iter_mut() {
                    self.thermal_condition(c, &format!("{at}/faces/{}", token(face)))?;
                }
            }
            self.each(t, &at, "surfaces", |w, s, at| match s.get_mut("condition") {
                Some(c) => w.thermal_condition(c, &format!("{at}/condition")),
                None => Ok(()),
            })?;
        }
        self.each(d, at, "porous", |w, z, at| {
            w.field(z, at, "min", Some(Q::Length))?;
            w.field(z, at, "max", Some(Q::Length))?;
            if let Some(r) = z.get_mut("resistance") {
                let at = format!("{at}/resistance");
                w.field(r, &at, "permeability", Some(Q::Area))?;
                w.field(r, &at, "inertial", Some(Q::InverseLength))?;
                w.plain(r, &at, "coefficient");
            }
            Ok(())
        })?;
        self.each(d, at, "fans", |w, f, at| {
            w.field(f, at, "min", Some(Q::Length))?;
            w.field(f, at, "max", Some(Q::Length))?;
            // `[Q, Δp]` points: a volume flow and a pressure rise.
            if let Some(Value::Array(points)) = f.get_mut("curve") {
                for (i, point) in points.iter_mut().enumerate() {
                    if let Value::Array(pair) = point {
                        for (j, v) in pair.iter_mut().enumerate() {
                            let q = if j == 0 { Q::VolumeFlow } else { Q::Pressure };
                            w.value(v, &format!("{at}/curve/{i}/{j}"), Some(q))?;
                        }
                    }
                }
            }
            w.plain(f, at, "relaxation");
            Ok(())
        })
    }

    fn step(&mut self, step: &mut Value, at: &str) -> Result<()> {
        let Some(p) = step.get_mut("procedure") else { return Ok(()) };
        let at = format!("{at}/procedure");
        let kind = p.get("type").and_then(Value::as_str).unwrap_or_default().to_string();
        match kind.as_str() {
            "frequency" => {
                self.plain(p, &at, "modes");
                self.field(p, &at, "min_hz", Some(Q::Frequency))?;
                self.field(p, &at, "max_hz", Some(Q::Frequency))?;
            }
            "buckle" => {
                self.plain(p, &at, "modes");
                self.plain(p, &at, "shift");
            }
            "explicit_dynamic" => {
                for key in ["end_time", "output_interval", "history_interval"] {
                    self.field(p, &at, key, Some(Q::Time))?;
                }
                for key in ["dt_scale", "max_increments", "hourglass", "bulk_viscosity"] {
                    self.plain(p, &at, key);
                }
                self.each(p, &at, "initial_velocities", |w, v, at| w.field(v, at, "velocity", Some(Q::Velocity)))?;
                self.each(p, &at, "rigid_walls", |w, wall, at| {
                    w.field(wall, at, "point", Some(Q::Length))?;
                    w.plain(wall, at, "normal");
                    w.friction(wall, at)
                })?;
                if let Some(ms) = p.get_mut("mass_scaling") {
                    self.field(ms, &format!("{at}/mass_scaling"), "target_increment", Some(Q::Time))?;
                }
            }
            "fluid" => {
                self.field(p, &at, "end_time", Some(Q::Time))?;
                self.field(p, &at, "time_step", Some(Q::Time))?;
                for key in ["steps", "frames", "magic"] {
                    self.plain(p, &at, key);
                }
                if let Some(s) = p.get_mut("steady") {
                    let at = format!("{at}/steady");
                    self.field(s, &at, "window", Some(Q::Time))?;
                    self.field(s, &at, "min_time", Some(Q::Time))?;
                    self.plain(s, &at, "tolerance");
                    self.plain(s, &at, "divergence");
                }
            }
            "heat_transfer" => {
                if let Some(t) = p.get_mut("transient").filter(|t| t.is_object()) {
                    let at = format!("{at}/transient");
                    self.field(t, &at, "end_time", Some(Q::Time))?;
                    self.field(t, &at, "increment", Some(Q::Time))?;
                    self.plain(t, &at, "theta");
                    self.field_or(t, &at, "initial_temperature", Q::Temperature, 0.0)?;
                }
            }
            "nonlinear_static" => {
                // Increments are fractions of the step, tolerances relative, the load path load
                // factors.
                for key in ["initial_increment", "min_increment", "max_increment", "max_increments", "max_iterations", "tolerances", "load_path"] {
                    self.plain(p, &at, key);
                }
                if let Some(imp) = p.get_mut("imperfection") {
                    let at = format!("{at}/imperfection");
                    self.plain(imp, &at, "mode");
                    // The buckling mode, scaled to a largest component of one, times `scale`.
                    self.field(imp, &at, "scale", Some(Q::Length))?;
                }
                if let Some(arc) = p.get_mut("arc_length") {
                    let at = format!("{at}/arc_length");
                    self.plain(arc, &at, "max_load_factor");
                    if let Some(limit) = arc.get_mut("displacement_limit") {
                        let at = format!("{at}/displacement_limit");
                        // Degrees of freedom 4 to 6 are rotations, in radians.
                        let rotation = limit.get("direction").and_then(Value::as_u64).is_some_and(|d| d >= 4);
                        self.plain(limit, &at, "node");
                        self.plain(limit, &at, "direction");
                        if rotation {
                            self.plain(limit, &at, "displacement");
                        } else {
                            self.field(limit, &at, "displacement", Some(Q::Length))?;
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// Refuse a number the map does not know: it would be read as SI unscaled.
    fn check_all_handled(&self, model: &Value) -> Result<()> {
        match self.unhandled(model, "") {
            None => Ok(()),
            Some(p) => invalid(format!(
                "cannot bring `model.json` to format 3: no unit is known for the number at {}",
                p.chars().take(200).collect::<String>()
            )),
        }
    }

    /// The pointer of the first number under `v` (at `at`) outside every handled value; a
    /// handled value is skipped whole.
    fn unhandled(&self, v: &Value, at: &str) -> Option<String> {
        if self.handled.contains(at) {
            return None;
        }
        match v {
            Value::Number(_) => Some(at.to_string()),
            Value::Array(items) => items.iter().enumerate().find_map(|(i, x)| self.unhandled(x, &format!("{at}/{i}"))),
            Value::Object(m) => m.iter().find_map(|(k, x)| self.unhandled(x, &format!("{at}/{}", token(k)))),
            _ => None,
        }
    }
}

/// `MCMESH01` → `MCMESH02`: node coordinates scaled by `length` to metres, every other byte as
/// it was. This reads only the `MCMESH01` header and node table (layout of format 2's
/// docs/formats/fieldcraft.md); the rest is checked when the current reader decodes the result.
fn mesh(bytes: &[u8], length: Factor, millimetre: Factor, in_millimetres: &HashSet<u32>) -> Result<Vec<u8>> {
    const HEADER: usize = 8 + 5 * 8;
    const NODE: usize = 4 + 3 * 8;
    if bytes.get(..8) != Some(b"MCMESH01".as_slice()) {
        return invalid("`mesh.bin` of format 2 must start with MCMESH01");
    }
    let header = bytes.get(..HEADER).ok_or_else(|| Error::Invalid("`mesh.bin` is truncated".into()))?;
    let count = header.get(8..16).and_then(|b| <[u8; 8]>::try_from(b).ok()).map(u64::from_le_bytes).unwrap_or(0);
    let nodes_end = usize::try_from(count)
        .ok()
        .and_then(|n| n.checked_mul(NODE))
        .and_then(|n| n.checked_add(HEADER))
        .filter(|end| *end <= bytes.len())
        .ok_or_else(|| Error::Invalid("`mesh.bin` lists more nodes than it holds".into()))?;
    let mut out = Vec::with_capacity(bytes.len());
    out.extend_from_slice(b"MCMESH02");
    out.extend_from_slice(header.get(8..).unwrap_or_default());
    for record in bytes.get(HEADER..nodes_end).unwrap_or_default().as_chunks::<NODE>().0 {
        let (label, coordinates) = record.split_at(4);
        out.extend_from_slice(label);
        let unit = match <[u8; 4]>::try_from(label).map(u32::from_le_bytes) {
            Ok(l) if in_millimetres.contains(&l) => millimetre,
            _ => length,
        };
        for c in coordinates.as_chunks::<8>().0 {
            let x = f64::from_le_bytes(*c);
            let si = unit.to_si(x);
            if !si.is_finite() {
                return invalid("`mesh.bin` holds a node coordinate that is not finite in metres");
            }
            out.extend_from_slice(&si.to_le_bytes());
        }
    }
    out.extend_from_slice(bytes.get(nodes_end..).unwrap_or_default());
    Ok(out)
}

/// Scaling by powers of ten and adding offsets on the decimal a value was written as: the
/// algorithm of the unit table of the release that made format 3 (fieldcraft-model
/// `units::decimal`), kept here as it was.
mod decimal {
    const EXACT_DIGITS: usize = 15;
    const POW10: [f64; 23] =
        [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22];

    fn short_decimal(v: f64) -> Option<(i128, i32)> {
        if v == 0.0 || !v.is_finite() {
            return None;
        }
        let text = format!("{v:e}");
        let (mantissa, exponent) = text.split_once('e')?;
        let exponent: i32 = exponent.parse().ok()?;
        let (negative, mantissa) = match mantissa.strip_prefix('-') {
            Some(m) => (true, m),
            None => (false, mantissa),
        };
        let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        if whole.len() + fraction.len() > EXACT_DIGITS {
            return None;
        }
        let mut m: i128 = 0;
        for c in whole.chars().chain(fraction.chars()) {
            m = m.checked_mul(10)?.checked_add(i128::from(c.to_digit(10)?))?;
        }
        let fraction_len = i32::try_from(fraction.len()).ok()?;
        Some((if negative { -m } else { m }, exponent.checked_sub(fraction_len)?))
    }

    fn from_decimal(mantissa: i128, exponent: i32) -> f64 {
        format!("{mantissa}e{exponent}").parse().unwrap_or(f64::NAN)
    }

    pub(super) fn shift(v: f64, k: i32) -> f64 {
        if k == 0 || v == 0.0 || !v.is_finite() {
            return v;
        }
        if let Some((m, e)) = short_decimal(v)
            && let Some(e) = e.checked_add(k)
        {
            return from_decimal(m, e);
        }
        let (mut v, mut k) = (v, k);
        let step = POW10.len() as i32 - 1;
        while k != 0 && v.is_finite() && v != 0.0 {
            let n = k.clamp(-step, step);
            let p = POW10.get(n.unsigned_abs() as usize).copied().unwrap_or(1.0);
            v = if n > 0 { v * p } else { v / p };
            k -= n;
        }
        v
    }

    /// An absolute temperature `f` in °F, in kelvin: (F + 459.67) × 5/9 on the decimal F was
    /// written as, as one correctly rounded division of two exact integers.
    pub(super) fn fahrenheit_to_kelvin(f: f64) -> f64 {
        const EXACT_INTEGER: i128 = 1 << 53;
        if f.fract() == 0.0 && f.abs() < 1e12 {
            return ((f as i128 * 100 + 45967) * 5) as f64 / 900.0;
        }
        if let Some((m, e)) = short_decimal(f)
            && let Ok(places) = u32::try_from(-e)
            && places <= 6
        {
            let p = places.max(2);
            let scale = 10i128.pow(p);
            let n = m.checked_mul(10i128.pow(p - places)).and_then(|m| m.checked_add(45967 * scale / 100)).and_then(|n| n.checked_mul(5));
            if let Some(n) = n.filter(|n| n.abs() < EXACT_INTEGER) {
                return n as f64 / (9 * scale) as f64;
            }
        }
        (f + 459.67) * 5.0 / 9.0
    }

    pub(super) fn add(a: f64, b: f64) -> f64 {
        if a == 0.0 || b == 0.0 {
            return a + b;
        }
        let (Some((ma, ea)), Some((mb, eb))) = (short_decimal(a), short_decimal(b)) else { return a + b };
        let e = ea.min(eb);
        let scale = |m: i128, from: i32| -> Option<i128> {
            let n = u32::try_from(from - e).ok()?;
            m.checked_mul(10i128.checked_pow(n)?)
        };
        match (scale(ma, ea), scale(mb, eb)) {
            (Some(x), Some(y)) => match x.checked_add(y) {
                Some(0) => 0.0,
                Some(sum) => from_decimal(sum, e),
                None => a + b,
            },
            _ => a + b,
        }
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    fn parts(model: Value, mesh: Vec<u8>) -> Parts {
        Parts { manifest: json!({"format": "fieldcraft/2", "generator": "FieldCraft 0.7.0"}), model, mesh, journal: json!([]), notes: Vec::new() }
    }

    /// An `MCMESH01` part with `nodes` and no elements.
    fn mesh01(nodes: &[(u32, [f64; 3])]) -> Vec<u8> {
        let mut b = b"MCMESH01".to_vec();
        for n in [nodes.len() as u64, 0, 0, 0, 0] {
            b.extend_from_slice(&n.to_le_bytes());
        }
        for (label, xyz) in nodes {
            b.extend_from_slice(&label.to_le_bytes());
            xyz.iter().for_each(|c| b.extend_from_slice(&c.to_le_bytes()));
        }
        b
    }

    /// About 80 MB of `model.json` (large sets): the step's time and memory are measured with
    /// `/usr/bin/time` (`cargo test -p fieldcraft-io --release migrating_a_large -- --ignored`).
    #[test]
    #[ignore = "a measurement, not a check"]
    fn migrating_a_large_model_takes_little_time_and_memory() {
        let labels: Vec<u32> = (1..=500_000).collect();
        let sets: serde_json::Map<String, Value> = (0..16).map(|k| (format!("set{k}"), json!({"nodes": labels, "elements": labels}))).collect();
        let mut p = parts(json!({"units": "m-kg-s-N-Pa", "sets": sets}), mesh01(&[]));
        println!("model.json: {} MB", serde_json::to_vec(&p.model).unwrap().len() / 1_000_000);
        let start = std::time::Instant::now();
        v2_to_v3(&mut p).unwrap();
        println!("step: {:.2} s", start.elapsed().as_secs_f64());
    }

    /// An `MCMESH01` part with `nodes`, tet4 `elements` (label, nodes), body entities and each
    /// element's region (an index into `bodies`).
    fn mesh01_bodies(nodes: &[(u32, [f64; 3])], elements: &[(u32, [u32; 4], u32)], bodies: &[&str]) -> Vec<u8> {
        let mut b = b"MCMESH01".to_vec();
        for n in [nodes.len(), elements.len(), bodies.len(), elements.len(), 0] {
            b.extend_from_slice(&(n as u64).to_le_bytes());
        }
        for (label, xyz) in nodes {
            b.extend_from_slice(&label.to_le_bytes());
            xyz.iter().for_each(|c| b.extend_from_slice(&c.to_le_bytes()));
        }
        for (label, of, _) in elements {
            b.extend_from_slice(&label.to_le_bytes());
            b.push(7);
            of.iter().for_each(|n| b.extend_from_slice(&n.to_le_bytes()));
        }
        for body in bodies {
            b.push(0);
            b.extend_from_slice(&(body.len() as u16).to_le_bytes());
            b.extend_from_slice(body.as_bytes());
        }
        for (label, _, region) in elements {
            b.extend_from_slice(&label.to_le_bytes());
            b.extend_from_slice(&region.to_le_bytes());
            b.extend_from_slice(&[0, 0]);
        }
        b
    }

    /// A model of metres built on a design, with `journal`, a tet of `Case` (the design's)
    /// on nodes 1–4 and a tet of `body` on nodes 5–8 (or `shared` with Case's node 4).
    fn metres_on_a_design(journal: Value, body: &str, shared: bool) -> Parts {
        let nodes: Vec<(u32, [f64; 3])> = (1..=8).map(|n| (n, [f64::from(n), 0.0, 0.0])).collect();
        let second = if shared { [4, 6, 7, 8] } else { [5, 6, 7, 8] };
        let mut p = parts(json!({"units": "m-kg-s-N-Pa"}), mesh01_bodies(&nodes, &[(1, [1, 2, 3, 4], 0), (2, second, 1)], &["Case", body]));
        p.manifest["geometry"] = json!({"design": {"path": "d.solidcraft", "sha256": "00"}});
        p.journal = journal;
        p
    }

    /// The x of each node after the step, by label.
    fn xs(mesh: &[u8]) -> Vec<f64> {
        (0..8).map(|i| f64::from_le_bytes(mesh[48 + i * 28 + 4..48 + i * 28 + 12].try_into().unwrap())).collect()
    }

    /// A mesh of the design in a model of metres is millimetres; a block, its mirror image and
    /// what was extruded from it stay metres; a body no command made is read as the design's,
    /// and the note says so.
    #[test]
    fn bodies_migrate_by_their_origin() {
        let block = json!({"command": "mesh.block", "params": {"body": "Bar"}});
        let tet = json!({"command": "mesh.tet.generate", "params": {"body": "Case"}});
        let mm = [0.001, 0.002, 0.003, 0.004];
        for (journal, second) in [
            (json!([block, tet]), "Bar"),
            (json!([block, {"command": "mesh.mirror", "params": {"plane": "Sym", "bodies": ["Bar"]}}, tet]), "Bar/mirror"),
            (json!([block, {"command": "mesh.mirror", "params": {"plane": "Sym", "bodies": ["Bar"], "as": {"Bar": "Left"}}}, tet]), "Left"),
            (json!([block, {"command": "mesh.hex.extrude", "params": {"source": [{"body": "Bar", "face": "top"}], "body": "Post"}}, tet]), "Post"),
        ] {
            let mut p = metres_on_a_design(journal, second, false);
            v2_to_v3(&mut p).unwrap();
            assert_eq!(xs(&p.mesh), [mm.as_slice(), &[5.0, 6.0, 7.0, 8.0]].concat(), "{second}");
            assert!(p.notes.len() == 1 && p.notes[0].contains("`Case`") && !p.notes[0].contains(second), "{second}: {:?}", p.notes);
        }
        // Nothing in the journal made `Bar`: read as the design's, saying so.
        let mut p = metres_on_a_design(json!([tet]), "Bar", false);
        v2_to_v3(&mut p).unwrap();
        assert_eq!(xs(&p.mesh), [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008]);
        assert!(p.notes[0].contains("No command in the journal made body (`Bar`)"), "{:?}", p.notes);
    }

    /// A node of both a mesh of the design and a block has no consistent scale: refused,
    /// naming both, with nothing changed.
    #[test]
    fn a_block_joined_to_a_mesh_of_the_design_is_refused() {
        let mut p = metres_on_a_design(json!([{"command": "mesh.block", "params": {"body": "Bar"}}]), "Bar", true);
        let before = p.clone();
        let e = v2_to_v3(&mut p).unwrap_err().to_string();
        for want in ["m-kg-s-N-Pa", "node 4", "`Case`", "`Bar`", "typed in metres", "no consistent", "unchanged", "0.7.x"] {
            assert!(e.contains(want), "`{want}` in {e}");
        }
        assert_eq!(p, before);
    }

    /// A relaxation field whose corners are no longer nodes (made before a remesh) takes its
    /// bodies' unit, else where it lies.
    #[test]
    fn relaxation_fields_follow_their_mesh() {
        let journal = json!([{"command": "mesh.block", "params": {"body": "Bar"}}, {"command": "mesh.tet.generate", "params": {"body": "Case"}}]);
        let cell = |x: f64| json!({"kind": "tet4", "at": [[x, 0.0, 0.0]], "moved": [[x, 0.0, 0.0]]});
        let mut p = metres_on_a_design(journal, "Bar", false);
        p.model["relaxations"] = json!({
            "OnCase": {"bodies": ["Case"], "cells": [cell(9.0)]},
            "OnBar": {"bodies": ["Bar"], "cells": [cell(1.5)]},
            "Inside": {"cells": [cell(2.5)]},
            "Outside": {"cells": [cell(7.5)]},
            "Node": {"cells": [cell(6.0)]},
        });
        v2_to_v3(&mut p).unwrap();
        let x = |f: &str| p.model["relaxations"][f]["cells"][0]["at"][0][0].as_f64().unwrap();
        assert_eq!([x("OnCase"), x("OnBar"), x("Inside"), x("Outside"), x("Node")], [0.009, 1.5, 0.0025, 7.5, 6.0]);
    }

    /// Hostile counts in `mesh.bin` are refused before anything is set aside for them.
    #[test]
    fn hostile_mesh_counts_are_refused_cheaply() {
        for at in [8, 16, 24, 32] {
            let mut b = mesh01_bodies(&[(1, [0.0; 3])], &[], &["Case"]);
            b[at..at + 8].copy_from_slice(&(u64::MAX / 2).to_le_bytes());
            assert!(MeshLayout::read(&b).is_err(), "{at}");
            let mut p = parts(json!({"units": "m-kg-s-N-Pa"}), b);
            p.manifest["geometry"] = json!({"design": {"path": "d.solidcraft", "sha256": "00"}});
            assert!(v2_to_v3(&mut p).is_err(), "{at}");
        }
    }

    #[test]
    fn scales_by_the_declared_system() {
        let mut p = parts(
            json!({"units": "mm-t-s-N-MPa", "materials": {"Steel": {"elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}, "density": 7.85e-9,
                "thermal": {"expansion": 1.2e-5, "conductivity": 0.05, "specific_heat": 4.6e8}}},
                "loads": {"g": {"type": "gravity", "acceleration": [0.0, 0.0, -9810.0]}, "hot": {"type": "temperature", "target": {"set": "all"}, "value": 20.0}},
                "datums": {"Mid": {"type": "plane", "point": [10.0, 0.0, 0.0], "normal": [1.0, 0.0, 0.0]}}}),
            mesh01(&[(1, [40.0, -0.0, 3.97])]),
        );
        v2_to_v3(&mut p).unwrap();
        let m = &p.model;
        assert_eq!(m["display_units"], "mm-t-s-N-MPa");
        assert!(m.get("units").is_none());
        let steel = &m["materials"]["Steel"];
        assert_eq!((steel["elastic"]["young"].as_f64(), steel["elastic"]["poisson"].as_f64()), (Some(2.1e11), Some(0.3)));
        assert_eq!(steel["density"].as_f64(), Some(7850.0));
        assert_eq!(steel["thermal"]["specific_heat"].as_f64(), Some(460.0));
        assert_eq!(steel["thermal"]["conductivity"].as_f64(), Some(0.05));
        assert_eq!(steel["thermal"]["expansion"].as_f64(), Some(1.2e-5));
        // A reference temperature left out meant 0 °C.
        assert_eq!(steel["thermal"]["reference_temperature"].as_f64(), Some(273.15));
        assert_eq!(m["loads"]["g"]["acceleration"], json!([0.0, 0.0, -9.81]));
        assert_eq!(m["loads"]["hot"]["value"].as_f64(), Some(293.15));
        assert_eq!(m["datums"]["Mid"]["point"], json!([0.01, 0.0, 0.0]));
        assert_eq!(&p.mesh[..8], b"MCMESH02");
        let coords: Vec<f64> = p.mesh[48 + 4..48 + 28].chunks(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
        assert_eq!(coords, [0.04, -0.0, 0.00397]);
        assert!(coords[1].is_sign_negative());
        assert_eq!(p.journal, json!({"units": "mm-t-s-N-MPa", "model_units": "mm-t-s-N-MPa", "commands": []}));
        assert_eq!(p.manifest["format"], "fieldcraft/3");
    }

    #[test]
    fn metres_stay_and_datums_follow_what_the_file_was_built_on() {
        let mut p = parts(
            json!({"units": "m-kg-s-N-Pa", "reference_points": {"RP": {"xyz": [0.1, 0.0, 2.0], "label": 7}},
                "boundary_conditions": {"b": {"type": "displacement", "target": {"set": "s"}, "values": [0.001, null, 0.0, 0.5, null, null]},
                    "t": {"type": "temperature", "target": {"set": "s"}, "value": 300.0}},
                "datums": {"P": {"type": "point", "point": [5.0, 0.0, 0.0]}}}),
            mesh01(&[(1, [0.1, 0.05, 0.005])]),
        );
        v2_to_v3(&mut p).unwrap();
        assert_eq!(p.model["reference_points"]["RP"]["xyz"], json!([0.1, 0.0, 2.0]));
        assert_eq!(p.model["reference_points"]["RP"]["label"], 7);
        // Rotations are radians; temperatures were kelvin.
        assert_eq!(p.model["boundary_conditions"]["b"]["values"], json!([0.001, null, 0.0, 0.5, null, null]));
        assert_eq!(p.model["boundary_conditions"]["t"]["value"].as_f64(), Some(300.0));
        // No design: the datum was typed against the block, in metres.
        assert_eq!(p.model["datums"]["P"]["point"], json!([5.0, 0.0, 0.0]));
        assert_eq!(p.mesh[..8], *b"MCMESH02");
        assert_eq!(p.mesh[8..], mesh01(&[(1, [0.1, 0.05, 0.005])])[8..]);
        // Built on a design, a datum was millimetres (the design's).
        let mut p = parts(json!({"units": "m-kg-s-N-Pa", "datums": {"P": {"type": "point", "point": [5.0, 0.0, 0.0]}}}), mesh01(&[]));
        p.manifest["geometry"] = json!({"design": {"path": "d.solidcraft", "sha256": "00"}});
        v2_to_v3(&mut p).unwrap();
        assert_eq!(p.model["datums"]["P"]["point"], json!([0.005, 0.0, 0.0]));
    }

    #[test]
    fn unspecified_reads_as_millimetre_tonne_second() {
        let mut p = parts(json!({"sections": {"s": {"set": "a", "material": "m", "kind": {"type": "plane_strain"}}}}), mesh01(&[]));
        v2_to_v3(&mut p).unwrap();
        assert_eq!(p.model["display_units"], "mm-t-s-N-MPa");
        // The default thickness of 1 was 1 mm.
        assert_eq!(p.model["sections"]["s"]["kind"]["thickness"].as_f64(), Some(0.001));
        assert_eq!(p.journal["units"], "mm-t-s-N-MPa");
    }

    #[test]
    fn inch_pound_and_millisecond_systems() {
        let mut p = parts(
            json!({"units": "in-lbf-s-psi", "loads": {"t": {"type": "temperature", "target": {"set": "s"}, "value": 68.0}, "p": {"type": "pressure", "set": "s", "magnitude": 1.0}}}),
            mesh01(&[(1, [1.0, 0.0, 0.0])]),
        );
        v2_to_v3(&mut p).unwrap();
        assert_eq!(p.model["loads"]["t"]["value"].as_f64(), Some(293.15), "68 °F exactly");
        assert_eq!(p.model["loads"]["p"]["magnitude"].as_f64(), Some(6894.757293168361));
        let x = f64::from_le_bytes(p.mesh[52..60].try_into().unwrap());
        assert_eq!(x, 0.0254);
        let mut p = parts(
            json!({"units": "mm-kg-ms-kN-GPa", "steps": [{"name": "Modes", "procedure": {"type": "frequency", "modes": 3, "max_hz": 2.0}, "output": []}],
                "materials": {"Steel": {"density": 7.85e-6, "elastic": {"type": "isotropic", "young": 210.0, "poisson": 0.3}}}}),
            mesh01(&[]),
        );
        v2_to_v3(&mut p).unwrap();
        // Cycles per millisecond.
        assert_eq!(p.model["steps"][0]["procedure"]["max_hz"].as_f64(), Some(2000.0));
        assert_eq!(p.model["steps"][0]["procedure"]["modes"], 3);
        assert_eq!(p.model["materials"]["Steel"]["density"].as_f64(), Some(7850.0));
        assert_eq!(p.model["materials"]["Steel"]["elastic"]["young"].as_f64(), Some(2.1e11));
    }

    #[test]
    fn refuses_what_it_cannot_scale() {
        for (model, why) in [
            (json!({"units": "furlong-firkin-fortnight"}), "unknown unit system"),
            (json!({"units": 3}), "not a unit system"),
            (json!({"materials": {"m": {"elastic": {"type": "isotropic", "young": 1.0, "poisson": 0.3, "warp": 2.0}}}}), "/materials/m/elastic/warp"),
            (json!({"novel": {"x": 1.0}}), "/novel/x"),
            (json!({"loads": {"l": {"type": "force", "target": {"set": "s"}, "vector": [1e308, 0.0, 0.0]}}}), ""),
            (json!([]), "not a JSON object"),
        ] {
            let mut p = parts(model.clone(), mesh01(&[]));
            let before = p.clone();
            let e = v2_to_v3(&mut p);
            if why.is_empty() {
                // 1e308 N is finite in SI too.
                assert!(e.is_ok(), "{model}");
                continue;
            }
            let e = e.unwrap_err().to_string();
            assert!(e.contains(why), "{model}: {e}");
            let _ = before;
        }
        let mut p = parts(json!({"units": "mm-t-s-N-MPa", "loads": {"l": {"type": "pressure", "set": "s", "magnitude": 1e305}}}), mesh01(&[]));
        assert!(v2_to_v3(&mut p).unwrap_err().to_string().contains("out of range"));
        for bad in [b"MCMESH02".to_vec(), b"MCMESH01".to_vec(), {
            let mut m = mesh01(&[(1, [0.0; 3])]);
            m[8] = 9;
            m
        }] {
            let mut p = parts(json!({}), bad);
            assert!(v2_to_v3(&mut p).is_err());
        }
        let mut p = parts(json!({}), mesh01(&[]));
        p.journal = json!({"not": "a list"});
        assert!(v2_to_v3(&mut p).is_err());
    }

    /// Every released format 1 file (through 1 → 2) takes the step: the map knows every number
    /// in it, and the mesh and the journal come out in format 3's form.
    #[test]
    fn every_format_1_fixture_takes_the_step() {
        use std::io::Read;
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats/v1");
        let mut seen = 0;
        for writer in std::fs::read_dir(&root).unwrap().map(|e| e.unwrap().path()).filter(|p| p.is_dir()) {
            for path in std::fs::read_dir(&writer).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|x| x == "fieldcraft")) {
                let bytes = std::fs::read(&path).unwrap();
                let mut zip = zip::ZipArchive::new(std::io::Cursor::new(&bytes)).unwrap();
                let mut part = |name: &str| {
                    let mut b = Vec::new();
                    zip.by_name(name).unwrap().read_to_end(&mut b).unwrap();
                    b
                };
                let json = |b: Vec<u8>| -> Value { serde_json::from_slice(&b).unwrap() };
                let mut p = Parts {
                    manifest: json(part("fieldcraft.json")),
                    model: json(part("model.json")),
                    mesh: part("mesh.bin"),
                    journal: json(part("journal.json")),
                    notes: Vec::new(),
                };
                super::super::step(&mut p, 1).unwrap();
                let journal = p.journal.clone();
                v2_to_v3(&mut p).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
                // A mesh of the design in a model of metres or inches says what changed.
                let stem = path.file_stem().unwrap().to_string_lossy();
                let noted = ["design-metres", "design-inches", "mixed-metres"].contains(&stem.as_ref());
                assert_eq!(p.notes.iter().any(|n| n.contains("made from a SolidCraft design")), noted, "{stem}: {:?}", p.notes);
                assert_eq!(p.manifest["format"], "fieldcraft/3");
                assert_eq!(p.journal["commands"], journal, "{}: the journal's commands are untouched", path.display());
                assert_eq!(&p.mesh[..8], b"MCMESH02");
                seen += 1;
            }
        }
        assert!(seen >= 30, "{seen}");
    }
}
