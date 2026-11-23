//! What FieldCraft opens, and what each file becomes, from SolidCraft's format table
//! (`solidcraft_io::formats`): one list, so FieldCraft opens what the pinned SolidCraft reads and
//! never keeps an extension list of its own (#1049).
//!
//! - Rows importing a design open as a [`Source::Design`]: SolidCraft's design files
//!   (`solidcraft_io::DESIGN_EXTS`: `.sc`, `.solidcraft`, `.solvecraft`). The row's plain `.json`
//!   design is left out: a `.json` file in FieldCraft is a command script (File › Run Script).
//! - Rows importing B-rep CAD open as a [`Source::Brep`]: STEP and IGES, and FreeCAD and Rhino
//!   (read through cadmpeg) in native builds ([`NATIVE`]).
//! - Rows importing a mesh, and PLY files (whose faces SolidCraft reads as a mesh too,
//!   `solidcraft_io::read_ply_mesh`), are [`Source::Triangles`]: an imported surface, for
//!   display, reference and fluid obstacles, never volume-meshed.
//!
//! KiCad boards, point clouds (`.xyz`, `.pts`) and 2D drawings are not opened.

use solidcraft_io::formats::{FormatSpec, ImportTarget, importable};

/// Whether this build offers the formats SolidCraft marks native-only (FreeCAD `.fcstd`, Rhino
/// `.3dm`): every build but the browser's.
pub const NATIVE: bool = cfg!(not(target_arch = "wasm32"));

/// What an opened file becomes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// A SolidCraft design, rebuilt from its timeline (`geometry.open`).
    Design,
    /// B-rep CAD read as one import feature (`geometry.open`).
    Brep,
    /// Triangles without a B-rep, read as imported surfaces (`geometry.import.mesh`).
    Triangles,
}

/// One format this build opens: the table's label and lower-case extensions, preferred first.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Format {
    pub label: &'static str,
    pub extensions: Vec<&'static str>,
    pub source: Source,
}

fn source_of_row(row: &FormatSpec) -> Option<Source> {
    match row.import.as_ref()?.target {
        ImportTarget::Design => Some(Source::Design),
        ImportTarget::Brep => Some(Source::Brep),
        ImportTarget::Mesh => Some(Source::Triangles),
        ImportTarget::PointCloud if row.extensions.iter().any(|e| solidcraft_io::is_ply_path(&format!("x.{e}"))) => Some(Source::Triangles),
        ImportTarget::PointCloud | ImportTarget::Kicad | ImportTarget::Sketch => None,
    }
}

/// The extensions of `row` FieldCraft opens as `source`.
fn extensions_of(row: &FormatSpec, source: Source) -> Vec<&'static str> {
    row.extensions
        .iter()
        .copied()
        .filter(|e| match source {
            Source::Design => solidcraft_io::is_design_ext(e),
            Source::Brep => true,
            Source::Triangles => row.import.as_ref().is_some_and(|i| i.target == ImportTarget::Mesh) || solidcraft_io::is_ply_path(&format!("x.{e}")),
        })
        .collect()
}

/// Every format this build opens, in the table's order.
pub fn formats() -> Vec<Format> {
    importable(NATIVE)
        .filter_map(|row| {
            let source = source_of_row(row)?;
            let extensions = extensions_of(row, source);
            (!extensions.is_empty()).then_some(Format { label: row.label, extensions, source })
        })
        .collect()
}

/// What a file becomes, by its extension (any case); `None` for a file FieldCraft does not open
/// as geometry.
pub fn source_of(path: &str) -> Option<Source> {
    let ext = std::path::Path::new(path).extension()?.to_str()?.to_ascii_lowercase();
    formats().into_iter().find(|f| f.extensions.contains(&ext.as_str())).map(|f| f.source)
}

/// The extensions opened as `source`, lower-case, without duplicates, in table order.
pub fn extensions(source: Source) -> Vec<&'static str> {
    let mut out: Vec<&'static str> = Vec::new();
    for e in formats().into_iter().filter(|f| f.source == source).flat_map(|f| f.extensions) {
        if !out.contains(&e) {
            out.push(e);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_cad_and_triangles_follow_the_table() {
        assert_eq!(extensions(Source::Design), vec!["sc", "solidcraft", "solvecraft"]);
        let brep = extensions(Source::Brep);
        for e in ["step", "stp", "igs", "iges"] {
            assert!(brep.contains(&e), "{e}");
        }
        assert_eq!(brep.contains(&"fcstd") && brep.contains(&"3dm"), NATIVE);
        let triangles = extensions(Source::Triangles);
        for e in ["stl", "3mf", "obj", "ply", "glb", "gltf"] {
            assert!(triangles.contains(&e), "{e}");
        }
        for (path, source) in [
            ("a/Part.SC", Some(Source::Design)),
            ("old.SolveCraft", Some(Source::Design)),
            ("p.StP", Some(Source::Brep)),
            ("p.IGES", Some(Source::Brep)),
            ("bin.STL", Some(Source::Triangles)),
            ("scan.ply", Some(Source::Triangles)),
            ("tray.3mf", Some(Source::Triangles)),
            ("script.json", None),
            ("cloud.xyz", None),
            ("board.kicad_pcb", None),
            ("sketch.dxf", None),
            ("model.fieldcraft", None),
            ("noext", None),
        ] {
            assert_eq!(source_of(path), source, "{path}");
        }
        assert_eq!(source_of("part.FCStd").is_some(), NATIVE);
        assert_eq!(source_of("part.3dm").is_some(), NATIVE);
    }
}
