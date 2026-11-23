//! Triangle files read as imported surfaces ([`TaggedSurface`]): 3MF through FieldCraft's own
//! reader (which keeps SolidCraft's face tags, [`crate::tagged_mesh`]), STL, OBJ, PLY and glTF
//! through SolidCraft's mesh readers (`solidcraft_io::mesh_import_feature`), each mesh one surface
//! of one face. Which files these are comes from SolidCraft's format table
//! ([`crate::formats::Source::Triangles`]).
//!
//! An imported surface is triangles without a B-rep: shown, used for reference and as a fluid
//! obstacle, never volume-meshed (#1049).
//!
//! Units: STL, OBJ and PLY state none, and are read in the unit the caller gives (millimetres by
//! default, as SolidCraft reads them); a 3MF that states no unit too. A 3MF that states its unit,
//! and glTF (metres by its specification), are read in their own, and a unit given for them is an
//! error. Factors come from the unit table (`fieldcraft_model::units`).
//!
//! A binary STL has no signature, so its size is its check: exactly 84 + 50·n bytes for the n
//! triangles its header gives; anything else that is not ASCII STL (SolidCraft's test: it starts
//! with `solid` and has facets) is refused, and so is an ASCII facet without exactly three
//! vertices or a file cut off inside one. Every STL coordinate must be one SolidCraft's reader
//! keeps ([`STL_COORD_BELOW`]). SolidCraft's readers leave such facets, and bad OBJ lines, PLY
//! faces and glTF parts past their caps, out with only a warning; FieldCraft refuses the file
//! instead, naming what was left out: a surface that silently lost triangles is open where the
//! file is damaged.
//!
//! Files are untrusted: the file size and the vertex and triangle counts are capped (SolidCraft's
//! readers cap their own counts as they read, [`MAX_ELEMENTS`] caps the result), coordinates must
//! be finite and within [`MAX_COORD`], indices must exist, and the SolidCraft call runs under
//! `guard`, so a panic inside it is an error.

use solidcraft_kernel::{KernelError, guard};

use crate::tagged_mesh::{MAX_ELEMENTS, MAX_FILE_BYTES, TaggedSurface, read_3mf_in};
use crate::{Error, Result, SurfaceMesh};
use fieldcraft_model::MAX_COORD;
use fieldcraft_model::units::{Quantity, Unit};

/// Bytes of a binary STL's header (80) and triangle count (4).
const STL_HEADER: usize = 84;
/// Bytes of a binary STL triangle: normal and three vertices, 12 `f32`, and a 2-byte attribute.
const STL_TRIANGLE: usize = 50;
/// SolidCraft's STL reader keeps a coordinate only below this magnitude (`solidcraft_io::read_stl`
/// at the pin, its `ok` test) and leaves a facet with another out: FieldCraft refuses it instead.
const STL_COORD_BELOW: f32 = 1e9; // unit-factors: SolidCraft's bound on a coordinate it keeps, not a conversion

/// What a read gives: the surfaces, and what the reader said that left nothing out.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub surfaces: Vec<TaggedSurface>,
    pub warnings: Vec<String>,
}

/// SolidCraft's mesh reader warnings that leave no geometry out (glTF materials' textures, and
/// skins and animations, `solidcraft_io` gltf_in.rs at the pin), by how they start: passed on.
/// Every other warning of its readers says it left geometry out (OBJ lines, PLY faces, glTF parts
/// and placements past its caps), and the read is refused: an import that silently lost triangles
/// would be an open surface. A warning this list does not know refuses too.
const NOTHING_LEFT_OUT: &[&str] = &["textures are not imported", "skins and animations are not imported"];

/// Millimetres per `unit`, a length unit of the unit table.
pub fn mm_per(unit: &Unit) -> Result<f64> {
    let invalid = |e: fieldcraft_model::units::UnitError| Error::Invalid(format!("`unit`: {e}"));
    let mm = Unit::parse("mm").map_err(invalid)?;
    Quantity::Length.express(Quantity::Length.to_si(1.0, unit).map_err(invalid)?, &mm).map_err(invalid)
}

/// The surfaces of the triangle file `file` (its name selects the reader), in millimetres, each
/// carrying the file's name. `unit` is the length unit of a file that states none (`None`:
/// millimetres).
pub fn read_surfaces(file: &str, bytes: &[u8], unit: Option<&Unit>) -> Result<Read> {
    let path = std::path::Path::new(file);
    let name = path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| file.to_string());
    if crate::formats::source_of(file) != Some(crate::formats::Source::Triangles) {
        return Err(Error::Invalid(format!("`{name}` is not a triangle file FieldCraft reads (3MF, STL, OBJ, PLY, glTF)")));
    }
    if bytes.len() > MAX_FILE_BYTES {
        return Err(Error::Invalid(format!("`{name}`: file too large ({} MB)", bytes.len() >> 20)));
    }
    let ext = path.extension().and_then(|e| e.to_str()).map(str::to_ascii_lowercase).unwrap_or_default();
    let unit_mm = unit.map(mm_per).transpose()?;
    let (mut surfaces, warnings) = if ext == "3mf" {
        (read_3mf_in(bytes, unit_mm).map_err(|e| Error::Invalid(format!("`{name}`: {e}")))?, Vec::new())
    } else {
        if unit_mm.is_some() && solidcraft_io::is_gltf_path(file) {
            return Err(Error::Invalid(format!("`{name}` states its own unit: glTF is in metres by its specification; leave `unit` out")));
        }
        if ext == "stl" {
            check_stl(&name, bytes)?;
        }
        solidcraft_meshes(&name, file, bytes, unit_mm)?
    };
    for s in &mut surfaces {
        s.file.clone_from(&name);
    }
    Ok(Read { surfaces, warnings })
}

/// Refuse an STL that is neither a binary STL of exactly the size its header gives nor ASCII STL
/// whose every facet has exactly three vertices, or that has a coordinate SolidCraft's reader
/// would leave out (not finite, or not below [`STL_COORD_BELOW`]).
fn check_stl(name: &str, bytes: &[u8]) -> Result<()> {
    let bad = |what: String| Err(Error::Invalid(format!("`{name}`: {what}")));
    let kept = |v: f32| v.is_finite() && v.abs() < STL_COORD_BELOW;
    let count = bytes.get(80..STL_HEADER).and_then(|b| <[u8; 4]>::try_from(b).ok()).map(|b| u32::from_le_bytes(b) as usize);
    let binary = count.and_then(|n| n.checked_mul(STL_TRIANGLE)).and_then(|b| b.checked_add(STL_HEADER)).is_some_and(|size| size == bytes.len());
    if binary {
        for (i, t) in bytes.get(STL_HEADER..).unwrap_or_default().as_chunks::<STL_TRIANGLE>().0.iter().enumerate() {
            // The normal (bytes 0..12) is not read; the vertices are bytes 12..48.
            let vertices = t.get(12..48).unwrap_or_default();
            if vertices.as_chunks::<4>().0.iter().any(|b| !kept(f32::from_le_bytes(*b))) {
                return bad(format!("triangle {} has a coordinate that is not a finite number in range", i + 1));
            }
        }
        return Ok(());
    }
    // SolidCraft's test for ASCII STL (`solidcraft_io::read_stl`), unchanged.
    let head = String::from_utf8_lossy(bytes.get(..bytes.len().min(512)).unwrap_or_default()).into_owned();
    let text = String::from_utf8_lossy(bytes);
    if !(head.trim_start().starts_with("solid") && text.contains("facet")) {
        return bad(match count {
            Some(n) => format!(
                "not an STL: {} bytes, but a binary STL of the {n} triangles its header gives is {} bytes, and it is not ASCII STL",
                bytes.len(),
                n.saturating_mul(STL_TRIANGLE).saturating_add(STL_HEADER)
            ),
            None => format!("not an STL: {} bytes, shorter than a binary STL's {STL_HEADER}-byte header, and not ASCII STL", bytes.len()),
        });
    }
    // Facet by facet: SolidCraft keeps a facet only with exactly three good vertices.
    let mut open: Option<(usize, usize)> = None;
    let mut facets = 0usize;
    for (i, line) in text.lines().enumerate() {
        let mut words = line.split_whitespace();
        match (words.next(), open.as_mut()) {
            (Some("facet"), None) => open = Some((i + 1, 0)),
            (Some("facet"), Some(_)) => return bad(format!("line {}: a facet starts before the last one ended", i + 1)),
            (Some("vertex"), Some((_, n))) => {
                *n += 1;
                let xyz: Vec<Option<f32>> = words.map(|w| w.parse::<f32>().ok()).collect();
                if xyz.len() != 3 || xyz.iter().any(|v| !v.is_some_and(kept)) {
                    return bad(format!("line {}: a vertex needs three finite numbers in range", i + 1));
                }
            }
            (Some("vertex"), None) => return bad(format!("line {}: a vertex outside a facet", i + 1)),
            (Some("endfacet"), Some((start, n))) => {
                if *n != 3 {
                    return bad(format!("the facet at line {start} has {n} vertices, not 3"));
                }
                facets += 1;
                open = None;
            }
            (Some("endfacet"), None) => return bad(format!("line {}: `endfacet` without its facet", i + 1)),
            _ => {}
        }
    }
    if let Some((start, _)) = open {
        return bad(format!("the file ends inside the facet at line {start}: it is cut off"));
    }
    if facets == 0 {
        return bad("no complete facet".to_string());
    }
    Ok(())
}

/// An STL, OBJ, PLY or glTF file through SolidCraft's reader, an STL, OBJ or PLY file in
/// `unit_mm` millimetres per unit: one surface per mesh, one face each, and the reader's warnings
/// that left nothing out ([`NOTHING_LEFT_OUT`]); one that left geometry out refuses the read.
fn solidcraft_meshes(name: &str, file: &str, bytes: &[u8], unit_mm: Option<f64>) -> Result<(Vec<TaggedSurface>, Vec<String>)> {
    let feature = guard("read mesh file", || {
        let read = match unit_mm {
            Some(mm) => solidcraft_io::mesh_import_feature_in(bytes, file, mm),
            None => solidcraft_io::mesh_import_feature(bytes, file),
        };
        read.map_err(|e| KernelError::Invalid(e.to_string()))
    })?;
    let (kept, left_out): (Vec<String>, Vec<String>) =
        feature.warnings.into_iter().partition(|w| NOTHING_LEFT_OUT.iter().any(|start| w.starts_with(start)));
    if !left_out.is_empty() {
        return Err(Error::Invalid(format!(
            "`{name}`: the reader left geometry out, so the surface would not be the file's: {}",
            left_out.join("; ")
        )));
    }
    let solidcraft_doc::FeatureKind::MeshImport { meshes, .. } = feature.kind else {
        return Err(Error::Failed(format!("`{file}` was not read as a mesh")));
    };
    let mut budget = 0usize;
    let mut out = Vec::with_capacity(meshes.len());
    for (k, m) in meshes.iter().enumerate() {
        let what = if m.name.trim().is_empty() { format!("Mesh{}", k + 1) } else { m.name.trim().to_string() };
        let (positions, []) = m.positions.as_chunks::<3>() else {
            return Err(Error::Invalid(format!("mesh `{what}`: an incomplete vertex")));
        };
        let (triangles, []) = m.triangles.as_chunks::<3>() else {
            return Err(Error::Invalid(format!("mesh `{what}`: an incomplete triangle")));
        };
        budget = budget.saturating_add(positions.len().max(triangles.len()));
        if budget > MAX_ELEMENTS {
            return Err(Error::Invalid(format!("more than {MAX_ELEMENTS} vertices or triangles")));
        }
        let positions = positions.to_vec();
        if positions.iter().flatten().any(|c| !(c.is_finite() && c.abs() <= MAX_COORD)) {
            return Err(Error::Invalid(format!("mesh `{what}`: a vertex is out of range")));
        }
        let triangles = triangles.to_vec();
        if triangles.is_empty() {
            return Err(Error::Invalid(format!("mesh `{what}` has no triangles")));
        }
        if let Some(v) = triangles.iter().flatten().find(|v| **v as usize >= positions.len()) {
            return Err(Error::Invalid(format!("mesh `{what}`: a triangle refers to vertex {v}, but there are {}", positions.len())));
        }
        let tri_face = vec![0; triangles.len()];
        out.push(TaggedSurface {
            name: what,
            file: String::new(),
            mesh: SurfaceMesh { positions, triangles, tri_face, ..SurfaceMesh::default() },
            faces: 1,
            original: None,
        });
    }
    if out.is_empty() {
        return Err(Error::Invalid(format!("`{file}` holds no mesh")));
    }
    Ok((out, kept))
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;

    use super::*;

    /// A unit tetrahedron: vertices, and its four faces wound outward.
    const TETRAHEDRON: [[f32; 3]; 4] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    const TET_TRIANGLES: [[usize; 3]; 4] = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];

    fn tet_triangles() -> [[[f32; 3]; 3]; 4] {
        TET_TRIANGLES.map(|t| t.map(|v| TETRAHEDRON[v]))
    }

    fn ascii_stl() -> String {
        let mut s = String::from("solid tet\n");
        for t in tet_triangles() {
            s.push_str("facet normal 0 0 0\nouter loop\n");
            for [x, y, z] in t {
                writeln!(s, "vertex {x} {y} {z}").unwrap();
            }
            s.push_str("endloop\nendfacet\n");
        }
        s.push_str("endsolid tet\n");
        s
    }

    fn binary_stl(tris: &[[[f32; 3]; 3]]) -> Vec<u8> {
        let mut b = vec![0u8; 80];
        b.extend((tris.len() as u32).to_le_bytes());
        for t in tris {
            b.extend([0u8; 12]);
            for v in t {
                for c in v {
                    b.extend(c.to_le_bytes());
                }
            }
            b.extend([0u8; 2]);
        }
        b
    }

    const TET_OBJ: &str = "o Tet\nv 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 3 2\nf 1 2 4\nf 2 3 4\nf 1 4 3\n";

    const TET_PLY: &str = "ply\nformat ascii 1.0\nelement vertex 4\nproperty float x\nproperty float y\nproperty float z\n\
        element face 4\nproperty list uchar int vertex_indices\nend_header\n\
        0 0 0\n1 0 0\n0 1 0\n0 0 1\n3 0 2 1\n3 0 1 3\n3 1 2 3\n3 0 3 2\n";

    /// The enclosed volume: the sum of a·(b×c)/6 over the triangles.
    fn volume(m: &SurfaceMesh) -> f64 {
        m.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| m.positions[i as usize]);
                let cross = [b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2], b[0] * c[1] - b[1] * c[0]];
                a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]
            })
            .sum::<f64>()
            / 6.0
    }

    #[test]
    fn stl_obj_and_ply_read_as_one_face_surfaces() {
        let cases: [(&str, &str, Vec<u8>, &str); 4] = [
            ("dir/part.stl", "part.stl", ascii_stl().into_bytes(), "part"),
            ("dir/part.STL", "part.STL", binary_stl(&tet_triangles()), "part"),
            ("dir/part.obj", "part.obj", TET_OBJ.as_bytes().to_vec(), "Tet"),
            ("dir/part.ply", "part.ply", TET_PLY.as_bytes().to_vec(), "part"),
        ];
        for (path, file, bytes, name) in cases {
            let surfaces = read_surfaces(path, &bytes, None).unwrap().surfaces;
            assert_eq!(surfaces.len(), 1, "{path}");
            let s = &surfaces[0];
            assert_eq!(s.name, name, "{path}");
            assert_eq!(s.file, file, "{path}");
            assert_eq!(s.faces, 1, "{path}");
            assert_eq!(s.mesh.triangles.len(), 4, "{path}");
            assert_eq!(s.mesh.tri_face, vec![0; 4], "{path}");
            assert!((volume(&s.mesh) - 1.0 / 6.0).abs() < 1e-6, "{path}: volume {}", volume(&s.mesh));
        }
    }

    #[test]
    fn garbage_and_truncated_files_are_errors() {
        let good_stl = ascii_stl().into_bytes();
        let one_triangle = binary_stl(&tet_triangles()[..1]);
        let mut huge_count = vec![0u8; 80];
        huge_count.extend(4_000_000_000u32.to_le_bytes());
        huge_count.extend([0u8; 50]);
        let random: Vec<u8> = (0..4096u32).map(|i| (i.wrapping_mul(2_654_435_761) >> 24) as u8).collect();
        let no_vertices = b"f 1 2 3\n".to_vec();
        let overflow_vertices = b"v 1e400 0 0\nv 0 1 0\nv 0 0 1\nf 1 2 3\n".to_vec();
        let cut_header = b"ply\nformat ascii 1.0\nelement vertex 3\n".to_vec();
        let huge_vertex_count =
            b"ply\nformat ascii 1.0\nelement vertex 4294967295\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n"
                .to_vec();
        let no_faces = b"ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 0\nproperty list uchar int vertex_indices\nend_header\n0 0 0\n1 0 0\n0 1 0\n".to_vec();
        let nan_stl =
            b"solid x\nfacet normal 0 0 0\nouter loop\nvertex nan 0 0\nvertex nan 1 0\nvertex nan 0 1\nendloop\nendfacet\nendsolid x\n".to_vec();
        // A binary STL whose size is not 84 + 50·n for its header's n: one byte short, one over.
        let short = binary_stl(&tet_triangles())[..233].to_vec();
        let mut long = binary_stl(&tet_triangles());
        long.push(0);
        let mut nan_binary = binary_stl(&tet_triangles());
        nan_binary[84 + 12..84 + 16].copy_from_slice(&f32::NAN.to_le_bytes());
        let mut infinite_binary = binary_stl(&tet_triangles());
        infinite_binary[84 + 50 + 20..84 + 50 + 24].copy_from_slice(&f32::INFINITY.to_le_bytes());
        let half_nan_ascii = ascii_stl().replacen("vertex 0 0 0", "vertex 0 nan 0", 1).into_bytes();
        // Cut after its last vertex line, and a facet of two vertices: SolidCraft would drop the facet.
        let good = ascii_stl();
        let cut_ascii = good.as_bytes()[..good.rfind("endloop").unwrap()].to_vec();
        let two_vertices = good.replacen("vertex 1 0 0\n", "", 1).into_bytes();
        // SolidCraft keeps a coordinate only below 1e9 (and FieldCraft's limit is 1e9 too).
        let at_bound = good.replacen("vertex 1 0 0", "vertex 1000000000 0 0", 1).into_bytes();
        // An OBJ cut mid-face, and a PLY face naming a vertex it does not have: the readers
        // would leave the line or face out, and the surface would be open.
        let cut_obj = TET_OBJ.as_bytes()[..TET_OBJ.len() - 3].to_vec();
        let bad_ply_face = TET_PLY.replacen("3 0 3 2", "3 0 3 9", 1).into_bytes();
        let cases: [(&str, Vec<u8>); 27] = [
            ("x.stl", Vec::new()),
            ("x.obj", Vec::new()),
            ("x.ply", Vec::new()),
            ("x.stl", vec![0u8; 84]),
            ("x.stl", one_triangle[..100].to_vec()),
            ("x.stl", huge_count),
            ("x.stl", b"solid x\nfacet\nendsolid".to_vec()),
            ("x.stl", nan_stl),
            ("x.stl", random.clone()),
            ("x.stl", short),
            ("x.stl", long),
            ("x.stl", nan_binary),
            ("x.stl", infinite_binary),
            ("x.stl", half_nan_ascii),
            ("x.stl", cut_ascii),
            ("x.stl", two_vertices),
            ("x.stl", at_bound),
            ("x.obj", cut_obj),
            ("x.ply", bad_ply_face),
            ("x.obj", random.clone()),
            ("x.ply", random),
            ("x.obj", no_vertices),
            ("x.obj", overflow_vertices),
            ("x.ply", cut_header),
            ("x.ply", huge_vertex_count),
            ("x.ply", no_faces),
            ("x.step", good_stl.clone()),
        ];
        let mut accepted = Vec::new();
        for (name, bytes) in cases {
            if let Ok(read) = read_surfaces(name, &bytes, None) {
                accepted.push(format!("{name} ({} bytes) read as {} surfaces", bytes.len(), read.surfaces.len()));
            }
        }
        // A glTF placing its triangle more often than SolidCraft reads (100 000 placements): the
        // placements past the cap are left out, so the file is refused, though the first ones read.
        match read_surfaces("x.gltf", &wide_gltf(100_001), None) {
            Ok(_) => accepted.push("a glTF past the placement cap".to_string()),
            Err(e) => assert!(e.to_string().contains("`x.gltf`: the reader left geometry out") && e.to_string().contains("100000 placements"), "{e}"),
        }
        if read_surfaces("x.json", &good_stl, None).is_ok() {
            accepted.push("x.json (the good ASCII STL) read".to_string());
        }
        assert!(accepted.is_empty(), "hostile input accepted: {accepted:#?}");
    }

    /// A unitless file is read in the unit given; a 3MF without a unit too, and one with its own
    /// unit, or glTF (metres by its specification), refuses another.
    #[test]
    fn units_scale_unitless_files_and_stated_units_win() {
        let inch = Unit::parse("in").unwrap();
        assert_eq!(mm_per(&inch).unwrap(), 25.4);
        assert!(mm_per(&Unit::parse("s").unwrap()).is_err(), "not a length");
        for (path, bytes) in [("p.stl", ascii_stl().into_bytes()), ("p.obj", TET_OBJ.as_bytes().to_vec()), ("p.ply", TET_PLY.as_bytes().to_vec())] {
            let s = read_surfaces(path, &bytes, Some(&inch)).unwrap().surfaces;
            assert!((volume(&s[0].mesh) - 25.4f64.powi(3) / 6.0).abs() < 1e-6, "{path}");
        }
        let e = read_surfaces("p.glb", b"glTF", Some(&inch)).unwrap_err().to_string();
        assert!(e.contains("glTF is in metres by its specification"), "{e}");
        let model = |unit: &str| {
            format!(
                r#"<model {unit} xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02"><resources><object id="1" name="T"><mesh><vertices><vertex x="0" y="0" z="0"/><vertex x="1" y="0" z="0"/><vertex x="0" y="1" z="0"/><vertex x="0" y="0" z="1"/></vertices><triangles><triangle v1="0" v2="2" v3="1"/><triangle v1="0" v2="1" v3="3"/><triangle v1="1" v2="2" v3="3"/><triangle v1="0" v2="3" v3="2"/></triangles></mesh></object></resources><build><item objectid="1"/></build></model>"#
            )
        };
        let package = |xml: String| {
            use std::io::Write as _;
            let mut zip = zip::ZipWriter::new(std::io::Cursor::new(Vec::new()));
            zip.start_file("3D/3dmodel.model", zip::write::SimpleFileOptions::default()).unwrap();
            zip.write_all(xml.as_bytes()).unwrap();
            zip.finish().unwrap().into_inner()
        };
        let unstated = package(model(""));
        let s = read_surfaces("t.3mf", &unstated, Some(&inch)).unwrap().surfaces;
        assert_eq!(s[0].mesh.positions[1], [25.4, 0.0, 0.0]);
        assert_eq!(read_surfaces("t.3mf", &unstated, None).unwrap().surfaces[0].mesh.positions[1], [1.0, 0.0, 0.0]);
        let stated = package(model(r#"unit="centimeter""#));
        assert_eq!(read_surfaces("t.3mf", &stated, None).unwrap().surfaces[0].mesh.positions[1], [10.0, 0.0, 0.0]);
        let e = read_surfaces("t.3mf", &stated, Some(&inch)).unwrap_err().to_string();
        assert!(e.contains("states its own unit (`centimeter`)"), "{e}");
    }

    /// Standard base64 (RFC 4648), for a glTF's embedded buffer.
    fn base64(bytes: &[u8]) -> String {
        const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::new();
        for c in bytes.chunks(3) {
            let n = (u32::from(c[0]) << 16) | (u32::from(*c.get(1).unwrap_or(&0)) << 8) | u32::from(*c.get(2).unwrap_or(&0));
            for k in 0..4 {
                out.push(if k <= c.len() { A[(n >> (18 - 6 * k) & 63) as usize] as char } else { '=' });
            }
        }
        out
    }

    /// A glTF of one triangle placed by `n` root nodes, side by side.
    fn wide_gltf(n: usize) -> Vec<u8> {
        let mut buffer = Vec::new();
        for v in [[0f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]] {
            for c in v {
                buffer.extend(c.to_le_bytes());
            }
        }
        for i in [0u16, 1, 2, 0] {
            buffer.extend(i.to_le_bytes());
        }
        let nodes: Vec<String> = (0..n).map(|i| format!("{{\"name\":\"n{i}\",\"mesh\":0,\"translation\":[{i},0,0]}}")).collect();
        let roots: Vec<String> = (0..n).map(|i| i.to_string()).collect();
        format!(
            r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[{}]}}],"nodes":[{}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1}}]}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3","min":[0,0,0],"max":[1,1,0]}},{{"bufferView":1,"componentType":5123,"count":3,"type":"SCALAR"}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36}},{{"buffer":0,"byteOffset":36,"byteLength":6}}],"buffers":[{{"byteLength":44,"uri":"data:application/octet-stream;base64,{}"}}]}}"#,
            roots.join(","),
            nodes.join(","),
            base64(&buffer)
        )
        .into_bytes()
    }

    /// The same glTF within SolidCraft's caps reads: the count is what refuses it above.
    #[test]
    fn a_shallow_gltf_reads() {
        let read = read_surfaces("x.gltf", &wide_gltf(3), None).unwrap();
        assert!(!read.surfaces.is_empty());
        assert!(read.warnings.is_empty(), "{:?}", read.warnings);
    }
}
