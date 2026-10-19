//! The binary mesh part of a `.fieldcraft` file (`mesh.bin`): nodes, elements and the geometry
//! tags that say how the mesh lies on the geometry.
//!
//! Little-endian throughout:
//!
//! ```text
//! "MCMESH01"                                   magic, 8 bytes
//! u64 nodes, u64 elements, u64 entities, u64 tagged elements, u64 tagged nodes
//! nodes      × { u32 label, f64 x, f64 y, f64 z }            labels ascending
//! elements   × { u32 label, u8 kind, u32 node × n(kind) }    labels ascending
//! entities   × { u8 kind (0 body, 1 face, 2 edge, 3 vertex),
//!                u16 length, body name (UTF-8), [u16 length, entity name] unless a body }
//! tagged elements × { u32 element, u32 region (0xFFFF_FFFF: none),
//!                     u8 count, count × { u8 face, u32 tag },
//!                     u8 count, count × { u8 edge, u32 tag } }  elements ascending
//! tagged nodes × { u32 node, u32 tag }                        nodes ascending
//! ```
//!
//! Element kinds are numbered by [`KIND_CODES`] (fixed by the format, never by the enum's
//! order). A tag is an index into the entity list. The part's integrity is checked by its SHA-256
//! in the manifest; decoding checks every count against the bytes present before reading, and
//! every value through the model's own checked methods.

use fieldcraft_model::geometry::MAX_ENTITIES;
use fieldcraft_model::{ElementKind, Entity, GeometryRef, MAX_ELEMENTS, MAX_NODES, Model};

use crate::{Error, MAX_FILE_BYTES, Result, invalid, model_error};

/// The magic bytes the part starts with (the `01` is the mesh encoding's version). A new encoding
/// is a new magic, brought in by a new file format whose migration step rewrites older parts.
pub const MAGIC: &[u8; 8] = b"MCMESH01";

/// Largest mesh part read.
pub const MAX_BYTES: u64 = MAX_FILE_BYTES;

const PART: &str = "mesh.bin";

/// No region tag.
const NO_TAG: u32 = u32::MAX;

/// The element kind codes of the format.
pub const KIND_CODES: [(ElementKind, u8); 13] = [
    (ElementKind::Rod2, 1),
    (ElementKind::Beam2, 2),
    (ElementKind::Tri3, 3),
    (ElementKind::Tri6, 4),
    (ElementKind::Quad4, 5),
    (ElementKind::Quad8, 6),
    (ElementKind::Tet4, 7),
    (ElementKind::Tet10, 8),
    (ElementKind::Wedge6, 9),
    (ElementKind::Hex8, 10),
    (ElementKind::Hex20, 11),
    (ElementKind::Rod3, 12),
    (ElementKind::Pyramid5, 13),
];

fn kind_code(kind: ElementKind) -> u8 {
    KIND_CODES.iter().find(|(k, _)| *k == kind).map_or(0, |(_, c)| *c)
}

fn kind_of(code: u8) -> Option<ElementKind> {
    KIND_CODES.iter().find(|(_, c)| *c == code).map(|(k, _)| *k)
}

/// Bytes of the fixed header.
const HEADER: usize = 8 + 5 * 8;
/// Smallest record of each table (for the count checks before reading).
const NODE_BYTES: u64 = 4 + 3 * 8;
const MIN_ELEMENT_BYTES: u64 = 4 + 1 + 2 * 4;
const MIN_ENTITY_BYTES: u64 = 1 + 2 + 1;
const MIN_ELEMENT_TAG_BYTES: u64 = 4 + 4 + 1 + 1;
const NODE_TAG_BYTES: u64 = 4 + 4;

fn put_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, v: usize) {
    out.extend_from_slice(&(v as u64).to_le_bytes());
}

fn put_name(out: &mut Vec<u8>, name: &str) -> Result<()> {
    let len =
        u16::try_from(name.len()).map_err(|_| Error::Limit(format!("geometry name `{}…` is too long", name.chars().take(40).collect::<String>())))?;
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(name.as_bytes());
    Ok(())
}

/// Encode the model's mesh and geometry tags.
pub fn encode(model: &Model) -> Result<Vec<u8>> {
    let c = model.classification();
    let element_bytes: usize = model.elements().values().map(|e| 5 + 4 * e.nodes.len()).sum();
    let mut out = Vec::with_capacity(HEADER + model.nodes().len() * NODE_BYTES as usize + element_bytes + c.nodes().len() * NODE_TAG_BYTES as usize);
    out.extend_from_slice(MAGIC);
    for n in [model.nodes().len(), model.elements().len(), c.entities().len(), c.elements().len(), c.nodes().len()] {
        put_u64(&mut out, n);
    }
    for (&label, node) in model.nodes() {
        put_u32(&mut out, label);
        for x in node.xyz {
            out.extend_from_slice(&x.to_le_bytes());
        }
    }
    for (&label, e) in model.elements() {
        put_u32(&mut out, label);
        out.push(kind_code(e.kind));
        for &n in &e.nodes {
            put_u32(&mut out, n);
        }
    }
    for g in c.entities() {
        out.push(match g.entity {
            Entity::Body => 0,
            Entity::Face(_) => 1,
            Entity::Edge(_) => 2,
            Entity::Vertex(_) => 3,
        });
        put_name(&mut out, &g.body)?;
        if let Some(name) = g.entity.name() {
            put_name(&mut out, name)?;
        }
    }
    for (&element, tags) in c.elements() {
        put_u32(&mut out, element);
        put_u32(&mut out, tags.region.unwrap_or(NO_TAG));
        for map in [&tags.faces, &tags.edges] {
            // An element has at most 12 edges and 6 faces; the model checks the numbers.
            out.push(u8::try_from(map.len()).map_err(|_| Error::Limit(format!("element {element} has too many tags")))?);
            for (&number, &tag) in map {
                out.push(number);
                put_u32(&mut out, tag);
            }
        }
    }
    for (&node, &tag) in c.nodes() {
        put_u32(&mut out, node);
        put_u32(&mut out, tag);
    }
    Ok(out)
}

/// A cursor over the part's bytes; every read fails cleanly at the end.
struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.bytes.len());
        match end.and_then(|e| self.bytes.get(self.pos..e)) {
            Some(s) => {
                self.pos += n;
                Ok(s)
            }
            None => invalid(format!("{PART} ends early (at byte {})", self.pos)),
        }
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N]> {
        let mut a = [0u8; N];
        a.copy_from_slice(self.take(N)?);
        Ok(a)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.array::<1>()?[0])
    }

    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.array()?))
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.array()?))
    }

    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.array()?))
    }

    fn name(&mut self, what: &str) -> Result<String> {
        let len = usize::from(self.u16()?);
        String::from_utf8(self.take(len)?.to_vec()).map_err(|_| Error::Invalid(format!("{PART}: a {what} name is not UTF-8")))
    }

    fn remaining(&self) -> u64 {
        self.bytes.len().saturating_sub(self.pos) as u64
    }
}

/// A count from the header, checked against its cap.
fn count(r: &mut Reader, what: &str, cap: usize) -> Result<usize> {
    let n = r.u64()?;
    match usize::try_from(n) {
        Ok(n) if n <= cap => Ok(n),
        _ => Err(Error::Limit(format!("{PART}: {n} {what} (at most {cap})"))),
    }
}

/// Labels must strictly increase: no duplicates, and one byte order for one mesh.
fn ascending(what: &str, previous: &mut Option<u32>, label: u32) -> Result<()> {
    if previous.is_some_and(|p| label <= p) {
        return invalid(format!("{PART}: {what} {label} is out of order or repeated"));
    }
    *previous = Some(label);
    Ok(())
}

/// Decode a mesh part into `model`, which must have no mesh yet.
pub fn decode_into(model: &mut Model, bytes: &[u8]) -> Result<()> {
    if !model.nodes().is_empty() || !model.elements().is_empty() {
        return invalid("the model already has a mesh");
    }
    if bytes.len() as u64 > MAX_BYTES {
        return Err(Error::Limit(format!("{PART} is larger than {} MB", MAX_BYTES >> 20)));
    }
    let mut r = Reader { bytes, pos: 0 };
    match r.take(MAGIC.len()).ok() {
        Some(m) if m == MAGIC.as_slice() => {}
        // Another encoding of ours: the file's format migrates it before decoding, so one left
        // here is not what the format says.
        Some(m) if m.starts_with(&MAGIC[..6]) => {
            let found = String::from_utf8_lossy(&m[6..]).chars().filter(|c| !c.is_control()).collect::<String>();
            return Err(Error::Version(format!(
                "{PART} is mesh encoding `{found}`; this version decodes `{}`",
                String::from_utf8_lossy(&MAGIC[6..])
            )));
        }
        _ => return invalid(format!("{PART} is not a FieldCraft mesh (no {} header)", String::from_utf8_lossy(MAGIC))),
    }
    let nodes = count(&mut r, "nodes", MAX_NODES)?;
    let elements = count(&mut r, "elements", MAX_ELEMENTS)?;
    let entities = count(&mut r, "geometry entities", MAX_ENTITIES)?;
    let tagged_elements = count(&mut r, "tagged elements", elements)?;
    let tagged_nodes = count(&mut r, "tagged nodes", nodes)?;
    // Every table's smallest size must fit in what is left, before anything is read.
    let least = [
        (nodes, NODE_BYTES),
        (elements, MIN_ELEMENT_BYTES),
        (entities, MIN_ENTITY_BYTES),
        (tagged_elements, MIN_ELEMENT_TAG_BYTES),
        (tagged_nodes, NODE_TAG_BYTES),
    ]
    .iter()
    .try_fold(0u64, |sum, &(n, size)| (n as u64).checked_mul(size).and_then(|b| sum.checked_add(b)));
    if least.is_none_or(|b| b > r.remaining()) {
        return invalid(format!("{PART} is too short for the counts in its header"));
    }

    let mut last = None;
    for _ in 0..nodes {
        let label = r.u32()?;
        ascending("node", &mut last, label)?;
        let xyz = [r.f64()?, r.f64()?, r.f64()?];
        model.set_node(label, xyz).map_err(|e| model_error(PART, e))?;
    }
    let mut last = None;
    for _ in 0..elements {
        let label = r.u32()?;
        ascending("element", &mut last, label)?;
        let code = r.u8()?;
        let kind = kind_of(code).ok_or_else(|| Error::Invalid(format!("{PART}: element {label} has unknown kind {code}")))?;
        let connectivity = (0..kind.node_count()).map(|_| r.u32()).collect::<Result<Vec<u32>>>()?;
        model.set_element(label, kind, connectivity).map_err(|e| model_error(PART, e))?;
    }
    for i in 0..entities {
        let kind = r.u8()?;
        let body = r.name("body")?;
        let entity = match kind {
            0 => Entity::Body,
            1 => Entity::Face(r.name("face")?),
            2 => Entity::Edge(r.name("edge")?),
            3 => Entity::Vertex(r.name("vertex")?),
            k => return invalid(format!("{PART}: geometry entity {i} has unknown kind {k}")),
        };
        let g = GeometryRef::new(body, entity).map_err(|e| model_error(PART, e))?;
        if model.geometry_tag(&g).map_err(|e| model_error(PART, e))? as usize != i {
            return invalid(format!("{PART}: {g} is listed twice among the geometry entities"));
        }
    }
    let mut last = None;
    for _ in 0..tagged_elements {
        let element = r.u32()?;
        ascending("tagged element", &mut last, element)?;
        let region = r.u32()?;
        if region != NO_TAG {
            model.tag_element(element, region).map_err(|e| model_error(PART, e))?;
        }
        for face in [true, false] {
            let mut previous = None;
            for _ in 0..r.u8()? {
                let (number, tag) = (r.u8()?, r.u32()?);
                if previous.is_some_and(|p| number <= p) {
                    return invalid(format!("{PART}: element {element}'s tags are out of order or repeated"));
                }
                previous = Some(number);
                let tagged = if face { model.tag_element_face(element, number, tag) } else { model.tag_element_edge(element, number, tag) };
                tagged.map_err(|e| model_error(PART, e))?;
            }
        }
    }
    let mut last = None;
    for _ in 0..tagged_nodes {
        let node = r.u32()?;
        ascending("tagged node", &mut last, node)?;
        model.tag_node(node, r.u32()?).map_err(|e| model_error(PART, e))?;
    }
    if r.remaining() != 0 {
        return invalid(format!("{PART} has {} bytes after its last table", r.remaining()));
    }
    Ok(())
}
