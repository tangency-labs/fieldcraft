//! Opening SolidCraft designs (`.solidcraft`): read the design, rebuild its parametric timeline
//! with SolidCraft, and hand back its bodies with persistent face names.
//!
//! The design is rebuilt from its authoritative inputs (sketches, features, parameters); cached
//! display meshes in the file are never used. [`open_design`] takes the file's bytes (reading
//! files is `fieldcraft_io::fs`'s job, which in the browser holds the page's files).

use solidcraft_doc::{Document, Model, naming};
use solidcraft_kernel::{KernelError, guard};

use crate::{Face, Result, Shape};

/// The grid centres are matched on, relative to their spread.
const GRID: f64 = 1e-6; // unit-factors: relative to the spread

/// Largest design file opened (SolidCraft's own limit).
pub const MAX_DESIGN_BYTES: usize = solidcraft_io::MAX_DESIGN_BYTES;

/// A rebuilt design.
#[derive(Debug, Clone)]
pub struct Design {
    pub name: String,
    /// The bodies at the end of the timeline, in the design's order, placed in the world: bodies of
    /// components moved by their occurrences, one body per placement.
    pub bodies: Vec<Shape>,
    /// Features that failed to rebuild. The design still opens with the bodies the rest made,
    /// as it does in SolidCraft.
    pub errors: Vec<FeatureError>,
    /// The partitions applied ([`Design::partitioned`]), in order; empty for the design as opened.
    pub cuts: Vec<crate::partition::CutReport>,
    /// The bodies as opened, before any partition: what [`Design::partitioned`] starts from.
    pub(crate) source: std::sync::Arc<Vec<Shape>>,
    /// The bodies as opened with the face merges on them ([`Design::merged`]): what
    /// [`Design::source_body`] and [`Design::features`] describe.
    pub(crate) opened: std::sync::Arc<Vec<Shape>>,
    /// The timeline the design was rebuilt from ([`Design::features`]).
    pub(crate) history: Option<std::sync::Arc<crate::features::History>>,
    /// The defeatures applied ([`Design::defeatured`]) and the design as opened before them.
    pub(crate) defeat: Option<std::sync::Arc<crate::defeature::Applied>>,
}

/// A feature of the design that failed to rebuild.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureError {
    pub feature: String,
    pub message: String,
    /// The bodies that may differ from the design as SolidCraft builds it: those the feature
    /// names, or every body when it names none of them.
    pub bodies: Vec<String>,
}

impl FeatureError {
    /// One line for the person opening the design, naming the bodies it may change.
    pub fn warning(&self) -> String {
        // The message is SolidCraft's own, in the design's units.
        format!("{} failed to build, so {} may differ from the design (SolidCraft: {})", self.feature, quoted(&self.bodies), self.message)
    }
}

fn quoted(bodies: &[String]) -> String {
    match bodies {
        [] => "the design".into(),
        [b] => format!("body \"{b}\""),
        _ => format!("bodies {}", bodies.iter().map(|b| format!("\"{b}\"")).collect::<Vec<_>>().join(", ")),
    }
}

/// The bodies of `bodies` that `feature` names anywhere in its saved form (targets, tools,
/// selections); all of them when it names none.
fn bodies_of(feature: Option<&solidcraft_doc::Feature>, bodies: &[Shape]) -> Vec<String> {
    fn strings<'a>(v: &'a serde_json::Value, out: &mut Vec<&'a str>) {
        match v {
            serde_json::Value::String(s) => out.push(s),
            serde_json::Value::Array(a) => a.iter().for_each(|v| strings(v, out)),
            serde_json::Value::Object(o) => o.values().for_each(|v| strings(v, out)),
            _ => {}
        }
    }
    let value = feature.and_then(|f| serde_json::to_value(f).ok()).unwrap_or_default();
    let mut named = Vec::new();
    strings(&value, &mut named);
    let hit: Vec<String> = bodies.iter().filter(|b| named.contains(&b.name.as_str())).map(|b| b.name.clone()).collect();
    if hit.is_empty() { bodies.iter().map(|b| b.name.clone()).collect() } else { hit }
}

impl Design {
    pub fn body(&self, name: &str) -> Option<&Shape> {
        self.bodies.iter().find(|b| b.name == name)
    }

    /// The body `name` of the design as opened (before any partition), with its merged faces.
    pub fn source_body(&self, name: &str) -> Option<&Shape> {
        self.opened.iter().find(|b| b.name == name)
    }

    /// The body of the design as opened that `name` is, or that partitions cut it from
    /// (`Base` for `Base/a/b1`).
    pub fn source_of(&self, name: &str) -> Option<&str> {
        self.source
            .iter()
            .map(|b| b.name.as_str())
            .filter(|s| name == *s || name.strip_prefix(s).is_some_and(|rest| rest.starts_with('/')))
            .max_by_key(|s| s.len())
    }
}

/// The newest design format the pinned SolidCraft reads (and the one it writes).
pub const DESIGN_FORMAT: u32 = solidcraft_doc::format::FORMAT;

/// Open a design from its file bytes (a `.solidcraft` container or a JSON design) and rebuild it.
/// A design of a newer format than [`DESIGN_FORMAT`] is [`crate::Error::NewerDesign`].
pub fn open_design(bytes: &[u8]) -> Result<Design> {
    let doc = guard("read design", || solidcraft_io::read_design(bytes).map_err(|e| KernelError::Invalid(format!("not a readable design: {e}"))))
        .map_err(|e| match saved_format(bytes) {
            Some(format) if format > DESIGN_FORMAT => crate::Error::NewerDesign { format, reads: DESIGN_FORMAT },
            _ => e.into(),
        })?;
    rebuild(&doc)
}

/// The format a design file says it has (`"format": "solidcraft/N"`), in a JSON design or in a
/// container's `document.json`; `None` when the file doesn't say. Read only when SolidCraft
/// refused the file, to tell a newer design from a broken one.
fn saved_format(bytes: &[u8]) -> Option<u32> {
    use std::io::Read;
    #[derive(serde::Deserialize)]
    struct Head {
        format: String,
    }
    if bytes.len() > MAX_DESIGN_BYTES {
        return None;
    }
    let json = if solidcraft_io::is_design_container(bytes) {
        let mut zip = zip::ZipArchive::new(std::io::Cursor::new(bytes)).ok()?;
        let mut out = Vec::new();
        zip.by_name("document.json").ok()?.take(MAX_DESIGN_BYTES as u64).read_to_end(&mut out).ok()?;
        std::borrow::Cow::Owned(out)
    } else {
        std::borrow::Cow::Borrowed(bytes)
    };
    let head: Head = serde_json::from_slice(&json).ok()?;
    solidcraft_doc::format::version_of(&head.format)
}

/// Rebuild an in-memory document (the design half of [`open_design`]).
pub(crate) fn rebuild(doc: &Document) -> Result<Design> {
    let model = guard("rebuild design", || {
        let mut m = Model::new();
        m.evaluate(doc);
        Ok(m)
    })?;
    let history = crate::features::History {
        doc: doc.clone(),
        state: model.state(),
        errors: model.results.iter().filter(|r| !r.skipped).filter_map(|r| Some((r.id, r.error.clone()?))).collect(),
    };
    // Placed in the world: a component's bodies through every occurrence of it (an instance is a
    // body of its own, named after its occurrence path), as SolidCraft displays and exports them.
    let state = guard("place components", || Ok(solidcraft_doc::world_state(doc, &model.state())))?;
    let mut bodies = Vec::with_capacity(state.bodies.len());
    for b in &state.bodies {
        let faces = guard("name faces", || {
            let names = naming::face_names(b);
            Ok(names
                .iter()
                .enumerate()
                .map(|(i, name)| Face { name: name.clone(), exact: naming::face_name_exact(b, i), members: Vec::new() })
                .collect())
        })?;
        bodies.push(Shape { name: b.name.clone(), body: b.body.clone(), faces, face_map: Vec::new(), cuts: Vec::new() });
    }
    // A feature that failed leaves a different solid than SolidCraft's: each one names the bodies
    // it would have changed, so opening never hands back a silently different solid.
    let errors: Vec<FeatureError> = model
        .results
        .iter()
        .filter(|r| !r.skipped)
        .filter_map(|r| {
            let message = r.error.clone()?;
            let bodies = bodies_of(doc.features.iter().find(|f| f.id == r.id), &bodies);
            Some(FeatureError { feature: r.name.clone(), message, bodies })
        })
        .collect();
    // The errors name the bodies as SolidCraft does; they follow the pieces' new names.
    let renamed = number_split_bodies(&mut bodies);
    let errors = errors
        .into_iter()
        .map(|mut e| {
            for b in &mut e.bodies {
                if let Some((_, new)) = renamed.iter().find(|(old, _)| old == b) {
                    *b = new.clone();
                }
            }
            e
        })
        .collect();
    Ok(Design {
        name: doc.name.clone(),
        source: std::sync::Arc::new(bodies.clone()),
        opened: std::sync::Arc::new(bodies.clone()),
        bodies,
        errors,
        cuts: Vec::new(),
        history: Some(std::sync::Arc::new(history)),
        defeat: None,
    })
}

/// A body name without the suffixes SolidCraft gives the pieces of a split body: ` (1)` for the
/// new piece and `'` for a piece kept under its old name, any number of times
/// (`Bracket (1) (1)' (1)` → `Bracket`).
pub(crate) fn split_base(name: &str) -> &str {
    let mut s = name;
    loop {
        let t = s.trim_end_matches('\'');
        let t = match t.strip_suffix(')').and_then(|r| r.rsplit_once(" (")) {
            Some((head, n)) if !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) => head,
            _ => t,
        };
        if t.len() == s.len() || t.is_empty() {
            return s;
        }
        s = t;
    }
}

/// Name the pieces of a body split in the design `<body>.1`, `<body>.2`, … by position, instead
/// of the suffixes the splits stacked up (`Bracket (1) (1)' (1)`), which say nothing about where
/// a piece is. Pieces are ordered by centroid, x first, then y, then z, on a grid of a millionth
/// of their spread, so the names follow the geometry, not the order of the design's timeline. A
/// group is left as it is when one of the new names is a body's name already. The pieces keep
/// the group's places in the list of bodies, in the order of their numbers. Returns each
/// renamed body's (old, new) name.
fn number_split_bodies(bodies: &mut [Shape]) -> Vec<(String, String)> {
    let mut renamed = Vec::new();
    let mut groups: std::collections::BTreeMap<String, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, b) in bodies.iter().enumerate() {
        groups.entry(split_base(&b.name).to_string()).or_default().push(i);
    }
    for (base, mut members) in groups {
        if members.len() < 2 || members.iter().all(|&i| bodies.get(i).is_none_or(|b| b.name == base)) {
            continue;
        }
        let names: Vec<String> = (1..=members.len()).map(|k| format!("{base}.{k}")).collect();
        if bodies.iter().any(|b| names.contains(&b.name)) {
            continue;
        }
        let centroid = |i: usize| bodies.get(i).and_then(|b| b.measure().ok()).map_or([0.0; 3], |m| m.centroid);
        let centres: Vec<[f64; 3]> = members.iter().map(|&i| centroid(i)).collect();
        let spread = (0..3)
            .map(|k| {
                let (lo, hi) = centres.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), c| (lo.min(c[k]), hi.max(c[k])));
                hi - lo
            })
            .fold(0.0f64, f64::max);
        let grid = if spread.is_finite() && spread > 0.0 { spread * GRID } else { 1.0 };
        let key = |j: usize| centres.get(j).map_or([0; 3], |c| c.map(|v| (v / grid).round() as i64));
        let mut order: Vec<usize> = (0..members.len()).collect();
        order.sort_by_key(|&j| (key(j), j));
        // The group keeps its places in the list, now in the order of the numbers.
        let ordered: Vec<Shape> = order.into_iter().filter_map(|j| members.get(j).and_then(|&i| bodies.get(i)).cloned()).collect();
        members.sort_unstable();
        for ((slot, mut shape), name) in members.into_iter().zip(ordered).zip(names) {
            renamed.push((std::mem::replace(&mut shape.name, name.clone()), name));
            if let Some(b) = bodies.get_mut(slot) {
                *b = shape;
            }
        }
    }
    renamed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_suffixes_are_stripped() {
        assert_eq!(split_base("Bracket (1) (1)' (1)"), "Bracket");
        assert_eq!(split_base("Bracket (1)'''"), "Bracket");
        assert_eq!(split_base("Bracket"), "Bracket");
        assert_eq!(split_base("Part (a)"), "Part (a)");
        assert_eq!(split_base("(1)"), "(1)");
        assert_eq!(split_base("'''"), "'''");
    }

    #[cfg(feature = "scripts")]
    /// The corpus bracket, split 29 times in the design, opens as `Bracket.1` … `Bracket.30`,
    /// numbered by position.
    #[test]
    fn pieces_of_a_split_body_are_numbered_by_position() {
        let script = include_str!("../../../tests/corpus/filleted-bracket-blocks/recipe.json");
        let d = open_design(&crate::script_design(script).unwrap()).unwrap();
        let names: Vec<&str> = d.bodies.iter().map(|b| b.name.as_str()).collect();
        let expected: Vec<String> = (1..=30).map(|k| format!("Bracket.{k}")).collect();
        assert_eq!(names, expected.iter().map(String::as_str).collect::<Vec<_>>(), "listed in the order of their numbers");
        let x = |n: &str| d.body(n).unwrap().measure().unwrap().centroid[0];
        assert!(x("Bracket.1") <= x("Bracket.30"));
    }
    use crate::sample::{self, Plate};

    fn face_areas(shape: &Shape) -> Vec<f64> {
        let mesh = shape.tessellate(0.01).unwrap();
        let mut areas = vec![0.0; shape.faces().len()];
        for (t, f) in mesh.triangles.iter().zip(&mesh.tri_face) {
            let [a, b, c] = t.map(|i| {
                let [x, y, z] = mesh.positions[i as usize];
                solidcraft_geom::Vec3 { x, y, z }
            });
            areas[*f as usize] += 0.5 * (b - a).cross(c - a).len();
        }
        areas
    }

    #[test]
    fn opens_rebuilds_and_tessellates_a_design() {
        let plate = Plate::default();
        let design = open_design(&sample::plate_design(&plate).unwrap()).unwrap();
        assert!(design.errors.is_empty(), "{:?}", design.errors);
        assert_eq!(design.bodies.len(), 1);
        let body = &design.bodies[0];
        let m = body.measure().unwrap();
        let r = plate.hole_diameter / 2.0;
        let expect = (plate.width * plate.depth - std::f64::consts::PI * r * r) * plate.thickness;
        assert!((m.volume - expect).abs() < 1e-3 * expect, "{} vs {expect}", m.volume);
        // Four sides, top, bottom and the hole (a cylinder is two half faces).
        assert_eq!(m.faces, 8);
        assert_eq!(body.faces().len(), 8);
        let mesh = body.tessellate(0.01).unwrap();
        let mut tags = mesh.tri_face.clone();
        tags.sort_unstable();
        tags.dedup();
        assert_eq!(tags.len(), 8);
        assert!(tags.iter().all(|&t| (t as usize) < body.faces().len()));
        // Every face has its own name.
        let mut names: Vec<&str> = body.faces().iter().map(|f| f.name.as_str()).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), 8, "{names:?}");
    }

    #[test]
    fn face_names_survive_a_parameter_change() {
        let narrow = Plate::default();
        let wide = Plate { width: 70.0, thickness: 15.0, ..Plate::default() };
        let a = open_design(&sample::plate_design(&narrow).unwrap()).unwrap();
        let b = open_design(&sample::plate_design(&wide).unwrap()).unwrap();
        let (a, b) = (&a.bodies[0], &b.bodies[0]);
        let (area_a, area_b) = (face_areas(a), face_areas(b));
        let pi = std::f64::consts::PI;
        // Each half of the hole.
        let hole = |p: &Plate| pi * p.hole_diameter * p.thickness / 2.0;
        let cap = |p: &Plate| p.width * p.depth - pi * p.hole_diameter * p.hole_diameter / 4.0;
        for face in a.faces() {
            let i = a.face_index(&face.name).unwrap() as usize;
            let j = b.face_index(&face.name).unwrap_or_else(|| panic!("{} is gone after the edit", face.name)) as usize;
            // The same face by name: the hole stays the hole, a cap stays a cap, and a side keeps
            // its own dimension (depth sides grow with the thickness only, width sides with both).
            let (x, y) = (area_a[i], area_b[j]);
            let close = |v: f64, e: f64| (v - e).abs() < 0.02 * e;
            let kind = if close(x, hole(&narrow)) {
                assert!(close(y, hole(&wide)), "{}: {y}", face.name);
                "hole"
            } else if close(x, cap(&narrow)) {
                assert!(close(y, cap(&wide)), "{}: {y}", face.name);
                "cap"
            } else if close(x, narrow.width * narrow.thickness) {
                assert!(close(y, wide.width * wide.thickness), "{}: {y}", face.name);
                "width side"
            } else {
                assert!(close(x, narrow.depth * narrow.thickness), "{}: {x}", face.name);
                assert!(close(y, wide.depth * wide.thickness), "{}: {y}", face.name);
                "depth side"
            };
            assert!(face.exact, "{} ({kind}) is named from geometry", face.name);
        }
    }

    #[test]
    fn a_newer_design_asks_for_a_newer_fieldcraft() {
        let newer = DESIGN_FORMAT + 1;
        let container = sample::plate_design(&Plate::default()).unwrap();
        let text = String::from_utf8(solidcraft_io::write_design(&solidcraft_io::read_design(&container).unwrap())).unwrap();
        let json = text.replacen(&format!("\"solidcraft/{DESIGN_FORMAT}\""), &format!("\"solidcraft/{newer}\""), 1);
        assert_ne!(json, text, "the sample design says its format");
        let error = open_design(json.as_bytes()).unwrap_err();
        assert_eq!(error, crate::Error::NewerDesign { format: newer, reads: DESIGN_FORMAT });
        assert_eq!(
            error.to_string(),
            format!(
                "this design was saved by a newer SolidCraft (format {newer}); this FieldCraft reads up to {DESIGN_FORMAT}. Update FieldCraft to open it"
            )
        );
        // The same inside a container, as the desktop app saves it.
        let mut zip = zip::ZipArchive::new(std::io::Cursor::new(&container[..])).unwrap();
        let mut out = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
        for i in 0..zip.len() {
            let mut part = zip.by_index(i).unwrap();
            let name = part.name().unwrap().to_string();
            let mut data = String::new();
            std::io::Read::read_to_string(&mut part, &mut data).unwrap();
            if name == "document.json" {
                data = data.replacen(&format!("\"solidcraft/{DESIGN_FORMAT}\""), &format!("\"solidcraft/{newer}\""), 1);
            }
            out.start_file(name, zip::write::SimpleFileOptions::default()).unwrap();
            std::io::Write::write_all(&mut out, data.as_bytes()).unwrap();
        }
        let container = out.finish().unwrap().into_inner();
        assert_eq!(open_design(&container).unwrap_err(), crate::Error::NewerDesign { format: newer, reads: DESIGN_FORMAT });
        // A broken design of a format this build reads is still SolidCraft's own error.
        let broken = text.replacen("\"features\"", "\"features\": 1, \"x\"", 1);
        assert!(matches!(open_design(broken.as_bytes()), Err(crate::Error::Invalid(_))));
    }

    #[test]
    fn hostile_bytes_are_an_error() {
        for bytes in [&b""[..], b"PK\x03\x04", b"{", b"{}", b"\xff\xfe\x00", b"[1,2,3]", br#"{"format":"solidcraft/9999","name":"x"}"#] {
            assert!(open_design(bytes).is_err(), "{bytes:?}");
        }
        assert!(open_design(&vec![b' '; MAX_DESIGN_BYTES + 1]).is_err());
    }
}
