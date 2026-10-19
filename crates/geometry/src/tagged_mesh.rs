//! Face-tagged meshes in 3MF packages (`docs/formats/tagged-mesh.md`).
//!
//! A 3MF package (3MF Core Specification 1.3) whose triangles may carry `ft:face`, the B-rep face
//! they lie on. SolidCraft defines the format and writes it (its `docs/formats/3mf-face-tags.md`);
//! any other 3MF reads as untagged surfaces (one face each). Reading is plain parsing from bytes (no kernel call, no file
//! system), so it also runs in the browser.
//!
//! Packages are untrusted: the file, its parts, zip entries, XML nesting, vertex and triangle
//! counts and build items are capped, coordinates must be finite and bounded, and every index
//! is range-checked. Every problem is an [`Error::Invalid`].

use std::borrow::Cow;
use std::collections::HashMap;
use std::io::{Cursor, Read};

use quick_xml::NsReader;
use quick_xml::events::{BytesStart, Event};
use quick_xml::name::ResolveResult;

use crate::{Error, Result, SurfaceMesh};

/// Largest package read.
pub const MAX_FILE_BYTES: usize = 512 << 20;
/// Largest uncompressed part read (zip bombs stop here).
const MAX_PART: u64 = 1 << 30;
const MAX_ENTRIES: usize = 10_000;
/// Most vertices, and most triangles, in a file and in all the surfaces read from it.
pub const MAX_ELEMENTS: usize = 20_000_000;
const MAX_BUILD_ITEMS: usize = 10_000;
const MAX_XML_DEPTH: usize = 256;
/// Largest coordinate magnitude accepted, in the file and after units and transforms: the FE
/// model's own limit, so every imported vertex can become a node.
use fieldcraft_model::MAX_COORD;
/// A transform whose determinant is below this fraction of its largest row length cubed
/// flattens the object (relative, so a uniform scale is never refused).
const MIN_RELATIVE_DET: f64 = 1e-12; // unit-factors: relative to the transform's scale

const CORE_NS: &str = "http://schemas.microsoft.com/3dmanufacturing/core/2015/02";
/// SolidCraft's face-tag namespace, version 1.
pub const FACE_TAGS_NS: &str = "https://github.com/bherbruck/solidcraft/formats/3mf-face-tags/2026/10";

/// One build item of a package: a named surface whose `tri_face` holds the face tags.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TaggedSurface {
    pub name: String,
    pub mesh: SurfaceMesh,
    /// Number of faces: `ft:faces` when given, else the largest tag + 1 (1 when untagged).
    pub faces: u32,
    /// The surface as read, when `mesh` is it stitched ([`crate::stitch`]): what stitching
    /// starts from again.
    pub original: Option<std::sync::Arc<SurfaceMesh>>,
}

fn invalid(msg: impl std::fmt::Display) -> Error {
    Error::Invalid(format!("3MF: {msg}"))
}

fn bad<T>(msg: impl std::fmt::Display) -> Result<T> {
    Err(invalid(msg))
}

fn read_part(zip: &mut zip::ZipArchive<Cursor<&[u8]>>, name: &str) -> Result<Option<Vec<u8>>> {
    let want = name.trim_start_matches('/').to_ascii_lowercase();
    let Some(idx) = (0..zip.len())
        .find(|i| zip.name_for_index(*i).is_some_and(|n| n.is_ok_and(|n| n.replace('\\', "/").trim_start_matches('/').to_ascii_lowercase() == want)))
    else {
        return Ok(None);
    };
    let mut f = zip.by_index(idx).map_err(|e| invalid(format!("{name}: {e}")))?;
    // The declared size only sizes the first allocation (it may lie); the read itself is capped.
    let mut out = Vec::with_capacity(usize::try_from(f.size().min(16 << 20)).unwrap_or(0));
    (&mut f).take(MAX_PART + 1).read_to_end(&mut out).map_err(|e| invalid(format!("{name}: {e}")))?;
    if out.len() as u64 > MAX_PART {
        return bad(format!("{name} expands to more than 1 GB"));
    }
    Ok(Some(out))
}

/// Which namespace an element or attribute name is in.
#[derive(Clone, Copy, PartialEq)]
enum Ns {
    /// 3MF Core, or no namespace at all (accepted for lenience).
    Core,
    FaceTags,
    Other,
}

fn ns_of(r: &ResolveResult) -> Ns {
    match r {
        ResolveResult::Bound(n) if n.as_ref() == CORE_NS => Ns::Core,
        ResolveResult::Bound(n) if n.as_ref() == FACE_TAGS_NS => Ns::FaceTags,
        ResolveResult::Unbound => Ns::Core,
        _ => Ns::Other,
    }
}

/// An element's attributes as (namespace, local name, value), borrowed from the event (values
/// are copied only when they hold entity references). Unprefixed attributes are Core.
type Attrs<'a> = Vec<(Ns, &'a str, Cow<'a, str>)>;

fn attrs<'a, R>(rd: &NsReader<R>, e: &'a BytesStart) -> Result<Attrs<'a>> {
    let mut out = Vec::with_capacity(4);
    for a in e.attributes() {
        let a = a.map_err(|err| invalid(format!("model XML: {err}")))?;
        let (ns, local) = rd.resolver().resolve_attribute(a.key);
        let value = a.normalized_value(quick_xml::XmlVersion::Implicit1_0).map_err(|err| invalid(format!("model XML: {err}")))?;
        out.push((ns_of(&ns), local.into_inner(), value));
    }
    Ok(out)
}

fn get<'s>(a: &'s Attrs, ns: Ns, key: &str) -> Option<&'s str> {
    a.iter().find(|(n, k, _)| *n == ns && *k == key).map(|(_, _, v)| v.as_ref())
}

fn int(a: &Attrs, ns: Ns, key: &str, what: &str) -> Result<Option<u32>> {
    match get(a, ns, key) {
        None => Ok(None),
        Some(v) => v.trim().parse::<u32>().map(Some).map_err(|_| invalid(format!("{what} `{key}` is not a whole number: `{v}`"))),
    }
}

fn coord(a: &Attrs, key: &str) -> Result<f64> {
    let v = get(a, Ns::Core, key).ok_or_else(|| invalid(format!("vertex without `{key}`")))?;
    v.trim().parse::<f64>().ok().filter(|x| x.is_finite() && x.abs() <= MAX_COORD).ok_or_else(|| invalid(format!("bad coordinate `{v}`")))
}

/// Row-vector affine transform `[m00 m01 m02 m10 … m32]` (3MF order).
type Xf = [f64; 12];
const IDENTITY: Xf = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];

fn parse_xf(s: Option<&str>) -> Result<Xf> {
    let Some(s) = s else { return Ok(IDENTITY) };
    let v: Vec<f64> = s
        .split_whitespace()
        .take(13)
        .map(|t| t.parse::<f64>().ok().filter(|x| x.is_finite() && x.abs() <= MAX_COORD))
        .collect::<Option<_>>()
        .ok_or_else(|| invalid(format!("bad transform `{s}`")))?;
    <Xf>::try_from(v.as_slice()).map_err(|_| invalid(format!("a transform has 12 numbers, got `{s}`")))
}

fn apply(m: &Xf, p: [f64; 3]) -> [f64; 3] {
    [p[0] * m[0] + p[1] * m[3] + p[2] * m[6] + m[9], p[0] * m[1] + p[1] * m[4] + p[2] * m[7] + m[10], p[0] * m[2] + p[1] * m[5] + p[2] * m[8] + m[11]]
}

fn det(m: &Xf) -> f64 {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6])
}

/// Millimetres per 3MF model unit, from the unit table (`fieldcraft_model::units`).
fn unit_mm(u: &str) -> Option<f64> {
    use fieldcraft_model::units::{Quantity, Unit};
    let symbol = match u {
        "micron" => "µm",
        "millimeter" => "mm",
        "centimeter" => "cm",
        "inch" => "in",
        "foot" => "ft",
        "meter" => "m",
        _ => return None,
    };
    let (from, mm) = (Unit::parse(symbol).ok()?, Unit::parse("mm").ok()?);
    Quantity::Length.express(Quantity::Length.to_si(1.0, &from).ok()?, &mm).ok()
}

#[derive(Clone, Default)]
struct Object {
    name: String,
    positions: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    /// One per triangle, or empty when the object is untagged.
    tags: Vec<u32>,
    faces: Option<u32>,
    components: bool,
}

impl Object {
    /// Checks that need the whole object: tagging is all or nothing, indices exist.
    fn check(&self, id: u32) -> Result<()> {
        if !self.tags.is_empty() && self.tags.len() != self.triangles.len() {
            return bad(format!("object {id}: {} of {} triangles have a face tag", self.tags.len(), self.triangles.len()));
        }
        let n = self.positions.len();
        if let Some(v) = self.triangles.iter().flatten().find(|v| **v as usize >= n) {
            return bad(format!("object {id}: a triangle refers to vertex {v}, but there are {n}"));
        }
        Ok(())
    }
}

#[derive(Default)]
struct Model {
    unit: f64,
    objects: HashMap<u32, Object>,
    build: Vec<(u32, Xf)>,
}

impl Model {
    fn add(&mut self, id: u32, o: Object) -> Result<()> {
        o.check(id)?;
        self.objects.insert(id, o);
        Ok(())
    }
}

fn parse_model(xml: &[u8]) -> Result<Model> {
    let mut rd = NsReader::from_reader(xml);
    rd.config_mut().trim_text(true);
    let mut buf = Vec::new();
    let mut m = Model { unit: 1.0, ..Default::default() };
    let mut depth = 0usize;
    let mut cur: Option<(u32, Object)> = None;
    let (mut vertices, mut triangles) = (0usize, 0usize);
    let mut seen_model = false;
    loop {
        buf.clear();
        let (ns, ev) = rd.read_resolved_event_into(&mut buf).map_err(|e| invalid(format!("model XML: {e}")))?;
        let ns = ns_of(&ns);
        let (e, empty) = match ev {
            Event::Start(e) => (e, false),
            Event::Empty(e) => (e, true),
            Event::End(end) => {
                depth = depth.saturating_sub(1);
                if ns == Ns::Core
                    && end.local_name().as_ref() == "object"
                    && let Some((id, o)) = cur.take()
                {
                    m.add(id, o)?;
                }
                continue;
            }
            Event::Eof => break,
            _ => continue,
        };
        if !empty {
            depth += 1;
            if depth > MAX_XML_DEPTH {
                return bad("XML nested too deeply");
            }
        }
        if ns != Ns::Core {
            continue;
        }
        let a = attrs(&rd, &e)?;
        match e.local_name().as_ref() {
            "model" => {
                seen_model = true;
                if let Some(u) = get(&a, Ns::Core, "unit") {
                    m.unit = unit_mm(u).ok_or_else(|| invalid(format!("unknown unit `{u}`")))?;
                }
            }
            "object" => {
                if cur.is_some() {
                    return bad("an object inside an object");
                }
                let id = int(&a, Ns::Core, "id", "object")?.ok_or_else(|| invalid("object without id"))?;
                if m.objects.contains_key(&id) {
                    return bad(format!("two objects with id {id}"));
                }
                let faces = int(&a, Ns::FaceTags, "faces", "object")?;
                if faces == Some(0) {
                    return bad(format!("object {id} has ft:faces=\"0\""));
                }
                let o = Object { name: get(&a, Ns::Core, "name").unwrap_or_default().trim().to_string(), faces, ..Default::default() };
                if empty {
                    m.add(id, o)?;
                } else {
                    cur = Some((id, o));
                }
            }
            "vertex" => {
                let Some((_, o)) = cur.as_mut() else { return bad("vertex outside an object") };
                vertices += 1;
                if vertices > MAX_ELEMENTS {
                    return bad(format!("more than {MAX_ELEMENTS} vertices"));
                }
                o.positions.push([coord(&a, "x")?, coord(&a, "y")?, coord(&a, "z")?]);
            }
            "triangle" => {
                let Some((id, o)) = cur.as_mut() else { return bad("triangle outside an object") };
                triangles += 1;
                if triangles > MAX_ELEMENTS {
                    return bad(format!("more than {MAX_ELEMENTS} triangles"));
                }
                let mut t = [0u32; 3];
                for (slot, k) in t.iter_mut().zip(["v1", "v2", "v3"]) {
                    *slot = int(&a, Ns::Core, k, "triangle")?.ok_or_else(|| invalid(format!("triangle without `{k}`")))?;
                }
                if let Some(f) = int(&a, Ns::FaceTags, "face", "triangle")? {
                    if let Some(n) = o.faces.filter(|n| f >= *n) {
                        return bad(format!("object {id}: face tag {f} is not below ft:faces=\"{n}\""));
                    }
                    o.tags.push(f);
                }
                o.triangles.push(t);
            }
            "component" => {
                let Some((_, o)) = cur.as_mut() else { return bad("component outside an object") };
                o.components = true;
            }
            "item" => {
                if m.build.len() >= MAX_BUILD_ITEMS {
                    return bad(format!("more than {MAX_BUILD_ITEMS} build items"));
                }
                if get(&a, Ns::Other, "path").is_some() {
                    return bad("not supported yet: build items in other model parts (production extension)");
                }
                let target = int(&a, Ns::Core, "objectid", "build item")?.ok_or_else(|| invalid("build item without objectid"))?;
                m.build.push((target, parse_xf(get(&a, Ns::Core, "transform"))?));
            }
            _ => {}
        }
    }
    if !seen_model {
        return bad("no <model> element");
    }
    if cur.is_some() {
        return bad("the model ends inside an object");
    }
    Ok(m)
}

/// Build item `k` as a surface: object `o` (checked) placed by `xf` and scaled to millimetres.
fn surface(mut o: Object, k: usize, id: u32, xf: &Xf, unit: f64) -> Result<TaggedSurface> {
    let d = det(xf);
    let scale = (0..3).map(|r| (xf[r * 3].powi(2) + xf[r * 3 + 1].powi(2) + xf[r * 3 + 2].powi(2)).sqrt()).fold(0.0, f64::max);
    if !(d.abs() > MIN_RELATIVE_DET * scale.powi(3)) {
        return bad(format!("build item {}: its transform is singular (it flattens object {id} to zero volume)", k + 1));
    }
    if o.components {
        return bad(format!("not supported yet: object {id} is made of components"));
    }
    if o.triangles.is_empty() {
        return bad(format!("object {id} has no triangles"));
    }
    for p in &mut o.positions {
        *p = apply(xf, *p).map(|c| c * unit);
        if !p.iter().all(|c| c.is_finite() && c.abs() <= MAX_COORD) {
            return bad(format!("object {id}: a placed vertex is out of range"));
        }
    }
    if d < 0.0 {
        o.triangles.iter_mut().for_each(|t| t.swap(1, 2));
    }
    if o.tags.is_empty() {
        o.tags = vec![0; o.triangles.len()];
    }
    let faces = o.faces.unwrap_or_else(|| o.tags.iter().max().map_or(1, |f| f.saturating_add(1)));
    let name = if o.name.is_empty() { format!("Mesh{}", k + 1) } else { o.name };
    Ok(TaggedSurface {
        name,
        mesh: SurfaceMesh { positions: o.positions, triangles: o.triangles, tri_face: o.tags, weld_distance: 0.0, ..SurfaceMesh::default() },
        faces,
        original: None,
    })
}

/// The surfaces of a 3MF package, one per build item, in millimetres.
pub fn read_3mf(bytes: &[u8]) -> Result<Vec<TaggedSurface>> {
    if bytes.len() > MAX_FILE_BYTES {
        return bad(format!("file too large ({} MB)", bytes.len() >> 20));
    }
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).map_err(|e| invalid(format!("not a zip package: {e}")))?;
    if zip.len() > MAX_ENTRIES {
        return bad(format!("{} zip entries", zip.len()));
    }
    // The model part named by the package relationships, else the usual name.
    let mut root = "3D/3dmodel.model".to_string();
    if let Some(rels) = read_part(&mut zip, "_rels/.rels")? {
        let mut rd = NsReader::from_reader(rels.as_slice());
        let mut buf = Vec::new();
        loop {
            match rd.read_event_into(&mut buf) {
                Ok(Event::Start(e)) | Ok(Event::Empty(e)) if e.local_name().as_ref() == "Relationship" => {
                    let a = attrs(&rd, &e)?;
                    if get(&a, Ns::Core, "Type").is_some_and(|t| t.ends_with("/3dmodel"))
                        && let Some(t) = get(&a, Ns::Core, "Target")
                    {
                        root = t.trim_start_matches('/').to_string();
                        break;
                    }
                }
                Ok(Event::Eof) | Err(_) => break,
                _ => {}
            }
            buf.clear();
        }
    }
    let xml = read_part(&mut zip, &root)?.ok_or_else(|| invalid(format!("no model part `{root}`")))?;
    let mut m = parse_model(&xml)?;
    drop(xml);
    if m.build.is_empty() {
        return bad("no build items");
    }
    // An object's last placement takes it; earlier ones copy it. Every placement counts against
    // one budget, so placing a big object many times stops at the cap.
    let mut uses: HashMap<u32, usize> = HashMap::new();
    for (id, _) in &m.build {
        *uses.entry(*id).or_default() += 1;
    }
    let mut budget = 0usize;
    let mut out = Vec::with_capacity(m.build.len());
    for (k, (id, xf)) in m.build.iter().enumerate() {
        let left = uses.get_mut(id).map(|n| std::mem::replace(n, n.saturating_sub(1))).unwrap_or(0);
        let o = if left > 1 { m.objects.get(id).cloned() } else { m.objects.remove(id) };
        let o = o.ok_or_else(|| invalid(format!("build item {} refers to missing object {id}", k + 1)))?;
        budget = budget.saturating_add(o.positions.len().max(o.triangles.len()));
        if budget > MAX_ELEMENTS {
            return bad(format!("the build items place more than {MAX_ELEMENTS} vertices or triangles"));
        }
        out.push(surface(o, k, *id, xf, m.unit)?);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;
    use std::io::Write as _;

    use super::*;

    /// The 3MF units, through the unit table, are the exact millimetre factors.
    #[test]
    fn units_of_3mf_are_exact_millimetres() {
        for (u, mm) in [("micron", 1e-3), ("millimeter", 1.0), ("centimeter", 10.0), ("inch", 25.4), ("foot", 304.8), ("meter", 1000.0)] {
            assert_eq!(unit_mm(u), Some(mm), "{u}");
        }
        assert_eq!(unit_mm("furlong"), None);
    }

    const HEAD: &str = r#"<model unit="millimeter" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:ft="https://github.com/bherbruck/solidcraft/formats/3mf-face-tags/2026/10"><resources>"#;

    /// A writer built from docs/formats/tagged-mesh.md alone (SolidCraft's is a separate
    /// program): shortest round-trip coordinates, one object and build item per surface.
    fn write_tagged(surfaces: &[(&str, &SurfaceMesh, u32)]) -> Vec<u8> {
        let mut x = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        let _ = writeln!(
            x,
            "<model unit=\"millimeter\" xml:lang=\"en-US\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\" xmlns:ft=\"{FACE_TAGS_NS}\">"
        );
        x.push_str(" <resources>\n");
        for (k, (name, m, faces)) in surfaces.iter().enumerate() {
            let _ = writeln!(x, "  <object id=\"{}\" type=\"model\" name=\"{name}\" ft:faces=\"{faces}\">\n   <mesh>\n    <vertices>", k + 1);
            for p in &m.positions {
                let _ = writeln!(x, "     <vertex x=\"{}\" y=\"{}\" z=\"{}\"/>", p[0], p[1], p[2]);
            }
            x.push_str("    </vertices>\n    <triangles>\n");
            for (t, f) in m.triangles.iter().zip(&m.tri_face) {
                let _ = writeln!(x, "     <triangle v1=\"{}\" v2=\"{}\" v3=\"{}\" ft:face=\"{f}\"/>", t[0], t[1], t[2]);
            }
            x.push_str("    </triangles>\n   </mesh>\n  </object>\n");
        }
        x.push_str(" </resources>\n <build>\n");
        for k in 0..surfaces.len() {
            let _ = writeln!(x, "  <item objectid=\"{}\"/>", k + 1);
        }
        x.push_str(" </build>\n</model>\n");
        let rels = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/></Relationships>\n";
        let types = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/></Types>\n";
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for (name, data) in [("[Content_Types].xml", types), ("_rels/.rels", rels), ("3D/3dmodel.model", x.as_str())] {
            let o = zip::write::SimpleFileOptions::default()
                .compression_method(zip::CompressionMethod::Deflated)
                .last_modified_time(zip::DateTime::default());
            zip.start_file(name, o).unwrap();
            zip.write_all(data.as_bytes()).unwrap();
        }
        zip.finish().unwrap().into_inner()
    }

    /// A package holding only a model part (no relationships: the usual name is used).
    fn package(model: &str) -> Vec<u8> {
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        zip.start_file("3D/3dmodel.model", zip::write::SimpleFileOptions::default()).unwrap();
        zip.write_all(model.as_bytes()).unwrap();
        zip.finish().unwrap().into_inner()
    }

    fn tet(extra: &str) -> String {
        format!(
            "{HEAD}<object id=\"1\" name=\"Tet\"><mesh><vertices><vertex x=\"0\" y=\"0\" z=\"0\"/><vertex x=\"1\" y=\"0\" z=\"0\"/><vertex x=\"0\" y=\"1\" z=\"0\"/><vertex x=\"0\" y=\"0\" z=\"1\"/></vertices>\
             <triangles><triangle v1=\"0\" v2=\"2\" v3=\"1\"/><triangle v1=\"0\" v2=\"1\" v3=\"3\"/><triangle v1=\"1\" v2=\"2\" v3=\"3\"/><triangle v1=\"0\" v2=\"3\" v3=\"2\"/></triangles></mesh></object>\
             </resources><build><item objectid=\"1\"{extra}/></build></model>"
        )
    }

    /// A SolidCraft block's tessellation goes out and comes back bit for bit, tags included.
    #[test]
    fn tessellation_round_trips_exactly() {
        let block = crate::block([0.1, -0.3, 1.0 / 3.0], [2.0, 3.7, 4.0 + 1e-9]).unwrap();
        let mesh = block.tessellate(0.01).unwrap();
        let faces = u32::try_from(block.measure().unwrap().faces).unwrap();
        let bytes = write_tagged(&[("Block", &mesh, faces), ("Again", &mesh, faces)]);
        assert_eq!(bytes, write_tagged(&[("Block", &mesh, faces), ("Again", &mesh, faces)]), "deterministic bytes");
        let got = read_3mf(&bytes).unwrap();
        assert_eq!(got.len(), 2);
        // The file carries no weld distance and no B-rep edges.
        let mesh = SurfaceMesh { weld_distance: 0.0, edges: Vec::new(), edge_faces: Vec::new(), ..mesh };
        assert_eq!(got[0], TaggedSurface { name: "Block".into(), mesh: mesh.clone(), faces: 6, original: None });
        assert_eq!(got[1].name, "Again");
        assert_eq!(got[1].mesh, mesh);
    }

    /// 3MF from other programs: no tags is one face; units and item transforms apply, and a
    /// mirroring transform keeps the winding outward.
    #[test]
    fn untagged_units_and_transforms() {
        let s = read_3mf(&package(&tet(""))).unwrap();
        // Uniform and uneven (but not flattening) scales are fine.
        for xf in ["0.001 0 0 0 0.001 0 0 0 0.001 0 0 0", "1000 0 0 0 1 0 0 0 0.001 5 5 5"] {
            assert!(read_3mf(&package(&tet(&format!(" transform=\"{xf}\"")))).is_ok(), "{xf}");
        }
        assert_eq!((s[0].faces, s[0].mesh.tri_face.clone()), (1, vec![0; 4]));
        let cm = tet(" transform=\"-1 0 0 0 1 0 0 0 1 0 0 2\"").replace("millimeter", "centimeter");
        let s = read_3mf(&package(&cm)).unwrap();
        let m = &s[0].mesh;
        assert_eq!(m.positions[1], [-10.0, 0.0, 20.0]);
        let p = |i: u32| m.positions[i as usize];
        let six_v: f64 = m
            .triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
                a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
            })
            .sum();
        assert!((six_v / 6.0 - 1000.0 / 6.0).abs() < 1e-9, "{}", six_v / 6.0);
    }

    /// Tags in another namespace, or under another prefix bound to ours, are read by namespace.
    #[test]
    fn tags_are_matched_by_namespace() {
        let other =
            tet("").replace("<triangle v1=\"0\" v2=\"2\" v3=\"1\"/>", "<triangle v1=\"0\" v2=\"2\" v3=\"1\" x:face=\"3\" xmlns:x=\"urn:else\"/>");
        assert_eq!(read_3mf(&package(&other)).unwrap()[0].faces, 1);
        let prefixed = tet("")
            .replace("xmlns:ft=", "xmlns:q=")
            .replace("<object id=\"1\"", "<object q:faces=\"4\" id=\"1\"")
            .replace("\"/><triangle", "\" q:face=\"2\"/><triangle")
            .replace("v3=\"2\"/></triangles>", "v3=\"2\" q:face=\"3\"/></triangles>");
        let s = read_3mf(&package(&prefixed)).unwrap();
        assert_eq!((s[0].faces, s[0].mesh.tri_face.clone()), (4, vec![2, 2, 2, 3]));
    }

    #[test]
    fn hostile_packages_are_errors() {
        let tagged = |a: &str, b: &str| {
            tet("")
                .replace("<triangle v1=\"0\" v2=\"2\" v3=\"1\"/>", &format!("<triangle v1=\"0\" v2=\"2\" v3=\"1\" {a}/>"))
                .replace("name=\"Tet\"", b)
        };
        let cases = [
            // truncated, not XML, nested too deep, no model
            tet("")[..300].to_string(),
            "not xml at all".to_string(),
            format!("{}{}", "<a>".repeat(1000), "</a>".repeat(1000)),
            "<other/>".to_string(),
            // out-of-range index, NaN, infinite and huge coordinates
            tet("").replace("v3=\"3\"/>", "v3=\"4\"/>"),
            tet("").replace("<vertex x=\"1\"", "<vertex x=\"NaN\""),
            tet("").replace("<vertex x=\"1\"", "<vertex x=\"inf\""),
            tet("").replace("<vertex x=\"1\"", "<vertex x=\"1e300\""),
            tet("").replace("v1=\"0\" v2=\"2\"", "v1=\"-1\" v2=\"2\""),
            tet("").replace("v1=\"0\" v2=\"2\"", "v1=\"99999999999\" v2=\"2\""),
            // missing attributes and objects
            tet("").replace(" y=\"0\" z=\"0\"/><vertex x=\"1\"", " z=\"0\"/><vertex x=\"1\""),
            tet("").replace("objectid=\"1\"", "objectid=\"9\""),
            tet("").replace("<item objectid=\"1\"/>", ""),
            // tags: partial, above ft:faces, zero faces, not a number
            tagged("ft:face=\"0\"", "name=\"Tet\""),
            tet("").replace("<triangle ", "<triangle ft:face=\"5\" ").replace("name=\"Tet\"", "name=\"Tet\" ft:faces=\"2\""),
            tet("").replace("name=\"Tet\"", "name=\"Tet\" ft:faces=\"0\""),
            tet("").replace("<triangle ", "<triangle ft:face=\"-2\" "),
            // transforms, units, components, duplicates
            tet(" transform=\"1 2\""),
            tet(" transform=\"1e300 0 0 0 1 0 0 0 1 0 0 0\""),
            tet(" transform=\"1e9 0 0 0 1 0 0 0 1 0 0 0\"").replace("millimeter", "meter"),
            tet(" transform=\"0 0 0 0 0 0 0 0 0 0 0 0\""),
            tet(" transform=\"1 0 0 0 1 0 0 0 1e-300 0 0 0\""),
            tet(" transform=\"1 0 0 1 0 0 0 0 1 0 0 0\""),
            tet("").replace("<vertex x=\"1\"", "<vertex x=\"2e9\""),
            tet("").replace("millimeter", "parsec"),
            format!(
                "{HEAD}<object id=\"1\"><components><component objectid=\"1\"/></components></object></resources><build><item objectid=\"1\"/></build></model>"
            ),
            format!("{HEAD}<object id=\"1\"/><object id=\"1\"/></resources><build><item objectid=\"1\"/></build></model>"),
            format!("{HEAD}<object id=\"1\"><mesh></mesh></object></resources><build><item objectid=\"1\"/></build></model>"),
        ];
        for c in &cases {
            assert!(read_3mf(&package(c)).is_err(), "{c}");
        }
        assert!(read_3mf(b"").is_err());
        assert!(read_3mf(b"PK\x03\x04garbage").is_err());
        let good = write_tagged(&[(
            "T",
            &SurfaceMesh {
                positions: vec![[0.0; 3]; 3],
                triangles: vec![[0, 1, 2]],
                tri_face: vec![0],
                weld_distance: 0.0,
                ..SurfaceMesh::default()
            },
            1,
        )]);
        for cut in [10, good.len() / 2, good.len() - 1] {
            assert!(read_3mf(&good[..cut]).is_err(), "truncated at {cut}");
        }
    }

    /// Placements share one budget, so a big object placed many times stops at the cap.
    #[test]
    fn many_placements_stay_bounded() {
        let verts: String = (0..3_000).map(|i| format!("<vertex x=\"{i}\" y=\"0\" z=\"0\"/>")).collect();
        let model = format!(
            "{HEAD}<object id=\"1\"><mesh><vertices>{verts}</vertices><triangles><triangle v1=\"0\" v2=\"1\" v3=\"2\"/></triangles></mesh></object></resources><build>{}</build></model>",
            "<item objectid=\"1\"/>".repeat(MAX_BUILD_ITEMS)
        );
        let e = read_3mf(&package(&model)).unwrap_err();
        assert!(e.to_string().contains("place more than"), "{e}");
        let too_many = format!("{HEAD}</resources><build>{}</build></model>", "<item objectid=\"1\"/>".repeat(MAX_BUILD_ITEMS + 1));
        assert!(read_3mf(&package(&too_many)).unwrap_err().to_string().contains("build items"));
    }
}
