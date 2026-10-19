//! The corpus runner behind `cargo xtask corpus [--part NAME] [--no-write]`: build every part of
//! the meshing test corpus ([`crate::corpus`]), tessellate it, run each [`Stage`] on it and write
//! the results table to docs/meshing/corpus.md.
//!
//! The table is a fixed part columns block (bodies, faces, the boundary surface and the time to
//! build it) followed by one column group per stage. A mesher joins the corpus by implementing
//! [`Stage`] over a [`Loaded`] part (its tessellated boundary surfaces are ready in
//! `Loaded::surfaces`) and adding itself to [`stages`]. A stage that fails on a part shows the
//! error in its first column; the run goes on.

use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

use crate::corpus::{self, CATEGORIES, Part};
use fieldcraft_geometry::{Design, FaceSurface, Measure, SurfaceMesh};
use fieldcraft_mesh::cdt::{CdtOptions, tetrahedralize};
use fieldcraft_mesh::hex::{Cell, HexOptions, Surface, hex_mesh, hex_mesh_cells};
use fieldcraft_mesh::octree::{OctreeOptions, octree_hex_mesh};
use fieldcraft_mesh::size::UniformSize;
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use fieldcraft_mesh::tet::canonical::CanonicalSoup;
use fieldcraft_mesh::tet::{QualityOptions, quality_mesh};

pub const FILE: &str = "docs/meshing/corpus.md";
/// Chordal tolerance of the boundary surface, relative to the part's bounding-box diagonal.
const RELATIVE_TOLERANCE: f64 = 1e-3; // unit-factors: relative to the diagonal

/// The wall-clock milliseconds since `start`, as the tables give times.
fn elapsed_ms(start: Instant) -> f64 {
    const MS_PER_SECOND: f64 = 1e3; // unit-factors: wall-clock time of the report, not a model quantity
    start.elapsed().as_secs_f64() * MS_PER_SECOND
}

/// A built part with its boundary surfaces.
pub struct Loaded<'a> {
    pub part: &'a Part,
    pub design: Design,
    /// SolidCraft's measures, one per body in the design's order.
    pub measures: Vec<Measure>,
    /// One per body, in the design's order.
    pub surfaces: Vec<SurfaceMesh>,
    /// The chordal tolerance the surfaces were made with (mm).
    pub tolerance: f64,
    /// The bounding-box diagonal of all bodies (mm).
    pub diagonal: f64,
}

/// A column group of the results table: a mesher, a check.
pub trait Stage {
    fn headers(&self) -> &'static [&'static str];
    /// One cell per header.
    fn run(&self, part: &Loaded) -> Result<Vec<String>, String>;
}

/// The stages, in column order. Meshers add theirs here.
fn stages() -> Vec<Box<dyn Stage>> {
    vec![Box::new(SurfaceVolume), Box::new(BoundaryCheck), Box::new(TetCdt), Box::new(TetQualityStage), Box::new(StructuredHex), Box::new(OctreeHex)]
}

/// The boundary surface's enclosed volume against SolidCraft's measured volume and the closed
/// form: a surface that is open, inverted or misplaced shows here first.
struct SurfaceVolume;

impl Stage for SurfaceVolume {
    fn headers(&self) -> &'static [&'static str] {
        &["enclosed mm³", "SolidCraft mm³", "rel. diff", "exact mm³"]
    }

    fn run(&self, l: &Loaded) -> Result<Vec<String>, String> {
        let enclosed: f64 = l.surfaces.iter().map(SurfaceMesh::enclosed_volume).sum();
        let measured: f64 = l.measures.iter().map(|m| m.volume).sum();
        let rel = if measured.abs() > 0.0 { (enclosed - measured) / measured } else { f64::NAN };
        Ok(vec![format!("{enclosed:.1}"), format!("{measured:.1}"), format!("{rel:.1e}"), l.part.volume.map_or("–".into(), |v| format!("{v:.1}"))])
    }
}

/// Boundary surface preparation (`fieldcraft_mesh::surface`): welded, closed, manifold, outward
/// and free of self-intersections, or the first defects found, which keep the part from being
/// volume meshed.
struct BoundaryCheck;

impl Stage for BoundaryCheck {
    fn headers(&self) -> &'static [&'static str] {
        &["boundary", "genus", "feature edges", "check ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Vec<String>, String> {
        let start = Instant::now();
        let (mut vertices, mut genus, mut features) = (0, Vec::new(), 0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
            let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).map_err(|e| format!("{}: {e}", body.name))?;
            vertices += s.points().len();
            genus.extend(s.shells().iter().map(|sh| if sh.void { format!("{} (void)", sh.genus) } else { sh.genus.to_string() }));
            features += s.feature_edges().len();
        }
        let ms = elapsed_ms(start);
        Ok(vec![format!("ok, {vertices} vertices"), genus.join(", "), features.to_string(), format!("{ms:.1}")])
    }
}

/// The constrained Delaunay tetrahedralization (`fieldcraft_mesh::cdt`) of each body's boundary:
/// tetrahedra, Steiner points on boundary edges, cavities refilled by gift-wrapping, and the
/// validity gate on the rounded mesh (every tetrahedron positive, volume equal to the boundary's).
struct TetCdt;

impl Stage for TetCdt {
    fn headers(&self) -> &'static [&'static str] {
        &["CDT tets", "Steiner", "cavities", "tet gate", "CDT ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Vec<String>, String> {
        let start = Instant::now();
        let (mut tets, mut steiner, mut cavities, mut inadmissible) = (0, 0, 0, 0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let (s, _) = tet_surface(mesh).map_err(|e| format!("{}: {e}", body.name))?;
            let m = tetrahedralize(&s, &CdtOptions::default()).map_err(|e| format!("{}: {e}", body.name))?;
            let report = m.gate(s.volume()).map_err(|e| format!("{}: rounded mesh: {e}", body.name))?;
            if !report.volume_matches() {
                return Err(format!("{}: rounded mesh: {:?}", body.name, report.into_result().err()));
            }
            inadmissible += report.invalid.len() + report.undetermined.len();
            tets += m.tets.len();
            steiner += m.steiner.len();
            cavities += m.stats.cavities;
        }
        let ms = elapsed_ms(start);
        let gate = if inadmissible == 0 { "pass".to_string() } else { format!("{inadmissible} flat") };
        Ok(vec![tets.to_string(), steiner.to_string(), cavities.to_string(), gate, format!("{ms:.0}")])
    }
}

/// A body's boundary as the tet mesher takes it ([`CanonicalSoup`], as `mesh.tet.generate` does).
pub fn tet_surface(mesh: &SurfaceMesh) -> fieldcraft_mesh::Result<(BoundarySurface, CanonicalSoup)> {
    let soup = CanonicalSoup::new(&mesh.positions, &mesh.triangles, &mesh.tri_face)?;
    let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
    Ok((BoundarySurface::new(&soup.positions, &soup.triangles, &soup.tri_face, &options)?, soup))
}

/// The quality stage (`fieldcraft_mesh::tet`) on each body's CDT, with the default options: the
/// mesh's tetrahedra, its extreme dihedral angles, the slivers (a dihedral angle below 10° or
/// above 170°) of the CDT and of the result, and the validity gate.
struct TetQualityStage;

impl Stage for TetQualityStage {
    fn headers(&self) -> &'static [&'static str] {
        &["quality tets", "min dihedral °", "max dihedral °", "min radius ratio", "slivers CDT → mesh", "quality gate", "quality ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Vec<String>, String> {
        let start = Instant::now();
        let (mut tets, mut lo, mut hi, mut rr, mut before, mut after) = (0, 180.0_f64, 0.0_f64, 1.0_f64, 0, 0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let (s, soup) = tet_surface(mesh).map_err(|e| format!("{}: {e}", body.name))?;
            let m = tetrahedralize(&s, &CdtOptions::default()).map_err(|e| format!("{}: {e}", body.name))?;
            let q = quality_mesh(&s, &m, soup.analytic(&fieldcraft_engine::analytic_surfaces(body)), &QualityOptions::default())
                .map_err(|e| format!("{}: {e}", body.name))?;
            let r = &q.report;
            tets += r.output.tets;
            lo = lo.min(r.output.min_dihedral);
            hi = hi.max(r.output.max_dihedral);
            rr = rr.min(r.output.min_radius_ratio);
            before += r.input.slivers;
            after += r.output.slivers;
        }
        let ms = elapsed_ms(start);
        Ok(vec![
            tets.to_string(),
            format!("{lo:.1}"),
            format!("{hi:.1}"),
            format!("{rr:.3}"),
            format!("{before} → {after}"),
            "pass".into(),
            format!("{ms:.0}"),
        ])
    }
}

/// Element size of the hex stages, relative to the part's bounding-box diagonal.
const HEX_RELATIVE_SIZE: f64 = 1.0 / 25.0;
/// Characters of a "not meshed" reason kept in the table.
const REASON_CHARS: usize = 240;

/// Structured all-hex meshing (`fieldcraft_mesh::hex`): the method that meshed each body, the
/// element count and the smallest scaled Jacobian of the hex8 mesh and of the hex20 mesh (mid-edge
/// nodes on the geometry), or why no structured method fits. A `partitioned` part's bodies are
/// the cells of one conforming mesh (`blocks`).
struct StructuredHex;

impl Stage for StructuredHex {
    fn headers(&self) -> &'static [&'static str] {
        &["hex", "hex elements", "hex8 min SJ", "hex20 min SJ", "hex ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Vec<String>, String> {
        let start = Instant::now();
        let size = UniformSize(HEX_RELATIVE_SIZE * l.diagonal);
        let linear = HexOptions { chordal_tolerance: l.tolerance, ..HexOptions::default() };
        let quadratic = HexOptions { quadratic: true, ..linear.clone() };
        if l.part.categories.iter().any(|c| c == "partitioned") {
            return blocks(l, &size, &linear, &quadratic).map(|mut cells| {
                cells.push(format!("{:.0}", elapsed_ms(start)));
                cells
            });
        }
        let (mut methods, mut elements, mut sj8, mut sj20) = (Vec::new(), 0usize, f64::INFINITY, Ok(f64::INFINITY));
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let surface_options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
            let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &surface_options)
                .map_err(|e| format!("{}: {e}", body.name))?;
            let faces: Vec<Option<Surface>> =
                body.face_surfaces().map_err(|e| format!("{}: {e}", body.name))?.into_iter().map(|f| f.map(analytic)).collect();
            let m_method = match hex_mesh(&s, &faces, &size, &linear) {
                Ok(m) => {
                    methods.push(m.report.method.name());
                    elements += m.element_count();
                    sj8 = sj8.min(m.report.min_scaled_jacobian);
                    m.report.method
                }
                Err(e) => {
                    let reason: String = e.to_string().chars().take(REASON_CHARS).collect();
                    let ms = elapsed_ms(start);
                    return Ok(vec![format!("no: {}: {reason}", body.name), String::new(), String::new(), String::new(), format!("{ms:.0}")]);
                }
            };
            // The hex20 run goes straight to the method that meshed the body.
            let quadratic = HexOptions { method: m_method, ..quadratic.clone() };
            sj20 = match (sj20, hex_mesh(&s, &faces, &size, &quadratic)) {
                (Ok(a), Ok(m)) => Ok(a.min(m.report.min_scaled_jacobian)),
                (Err(e), _) => Err(e),
                (_, Err(e)) => Err(format!("{}: {}", body.name, e.to_string().chars().take(REASON_CHARS).collect::<String>())),
            };
        }
        methods.sort_unstable();
        methods.dedup();
        let ms = elapsed_ms(start);
        let sj20 = sj20.map_or_else(|e| format!("error: {e}"), |v| format!("{v:.3}"));
        Ok(vec![methods.join(", "), elements.to_string(), format!("{sj8:.3}"), sj20, format!("{ms:.0}")])
    }
}

/// The bodies of a partitioned part as the cells of one conforming mesh: method, elements and
/// the hex8 and hex20 smallest scaled Jacobians.
fn blocks(l: &Loaded, size: &UniformSize, linear: &HexOptions, quadratic: &HexOptions) -> Result<Vec<String>, String> {
    let mut surfaces = Vec::with_capacity(l.surfaces.len());
    let mut faces: Vec<Vec<Option<Surface>>> = Vec::with_capacity(l.surfaces.len());
    for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
        let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
        surfaces.push(BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).map_err(|e| format!("{}: {e}", body.name))?);
        faces.push(body.face_surfaces().map_err(|e| format!("{}: {e}", body.name))?.into_iter().map(|f| f.map(analytic)).collect());
    }
    let cells: Vec<Cell> =
        l.design.bodies.iter().zip(&surfaces).zip(&faces).map(|((b, surface), faces)| Cell { name: &b.name, surface, faces }).collect();
    let m = match hex_mesh_cells(&cells, size, linear) {
        Ok(m) => m,
        Err(e) => {
            return Ok(vec![
                format!("no: {}", e.to_string().chars().take(REASON_CHARS).collect::<String>()),
                String::new(),
                String::new(),
                String::new(),
            ]);
        }
    };
    let sj20 =
        hex_mesh_cells(&cells, size, quadratic).map_or_else(|e| format!("error: {e}"), |q| format!("{:.3}", q.mesh.report.min_scaled_jacobian));
    Ok(vec![
        format!("blocks ({} cells, {} shared faces)", cells.len(), m.interfaces),
        m.mesh.element_count().to_string(),
        format!("{:.3}", m.mesh.report.min_scaled_jacobian),
        sj20,
    ])
}

/// The geometry's analytic face surface in the mesher's terms.
pub fn analytic(s: FaceSurface) -> Surface {
    match s {
        FaceSurface::Plane { normal, point } => Surface::Plane { point, normal },
        FaceSurface::Cylinder { origin, axis, radius } => Surface::Cylinder { origin, axis, radius },
        FaceSurface::Cone { apex, axis, half_angle } => Surface::Cone { apex, axis, half_angle },
        FaceSurface::Sphere { center, radius } => Surface::Sphere { center, radius },
        FaceSurface::Torus { center, axis, major, minor } => Surface::Torus { center, axis, major, minor },
    }
}

/// Element cap of the octree hex column: well above every part that meshes, and low enough that a
/// part whose grid explodes (a thin spline wall refined everywhere) is refused in seconds instead
/// of optimised for an hour.
const OCTREE_MAX_ELEMENTS: usize = 500_000;

/// The unstructured all-hex mesher (`fieldcraft_mesh::octree`) at a size of 1/25 of the part's
/// diagonal: elements, the smallest and mean scaled Jacobian, and the mesh volume against
/// SolidCraft's. Every element is proven valid before the mesher returns; a part it cannot mesh
/// gets the reason, as in the structured hex column.
struct OctreeHex;

impl Stage for OctreeHex {
    fn headers(&self) -> &'static [&'static str] {
        &["octree hexes", "min SJ", "mean SJ", "hex vol. diff", "octree ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Vec<String>, String> {
        let start = Instant::now();
        let (mut hexes, mut min_sj, mut sum_sj, mut volume) = (0usize, f64::INFINITY, 0.0, 0.0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
            let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).map_err(|e| format!("{}: {e}", body.name))?;
            let options = OctreeOptions { max_elements: OCTREE_MAX_ELEMENTS, ..OctreeOptions::default() };
            let m = match octree_hex_mesh(&s, &UniformSize(HEX_RELATIVE_SIZE * l.diagonal), &options) {
                Ok(m) => m,
                Err(e) => {
                    let reason: String = e.to_string().chars().take(REASON_CHARS).collect();
                    let ms = elapsed_ms(start);
                    return Ok(vec![format!("no: {}: {reason}", body.name), String::new(), String::new(), String::new(), format!("{ms:.0}")]);
                }
            };
            hexes += m.hexes.len();
            min_sj = min_sj.min(m.report.min_scaled_jacobian);
            sum_sj += m.report.mean_scaled_jacobian * m.hexes.len() as f64;
            volume += m.report.volume;
        }
        let measured: f64 = l.measures.iter().map(|m| m.volume).sum();
        let rel = if measured.abs() > 0.0 { (volume - measured) / measured } else { f64::NAN };
        let ms = elapsed_ms(start);
        Ok(vec![hexes.to_string(), format!("{min_sj:.3}"), format!("{:.3}", sum_sj / hexes.max(1) as f64), format!("{rel:.1e}"), format!("{ms:.0}")])
    }
}

fn load(part: &Part) -> Result<(Loaded<'_>, f64), String> {
    let start = Instant::now();
    let design = corpus::load(part).map_err(|e| e.to_string())?;
    if let Some(e) = design.errors.first() {
        return Err(format!("feature {} failed: {}", e.feature, e.message));
    }
    let measures = design.bodies.iter().map(|b| b.measure().map_err(|e| format!("measure {}: {e}", b.name))).collect::<Result<Vec<_>, _>>()?;
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for m in &measures {
        for i in 0..3 {
            lo[i] = lo[i].min(m.bbox_min[i]);
            hi[i] = hi[i].max(m.bbox_max[i]);
        }
    }
    let diagonal = (0..3).map(|i| (hi[i] - lo[i]).powi(2)).sum::<f64>().sqrt();
    let tolerance = RELATIVE_TOLERANCE * if diagonal.is_finite() && diagonal > 0.0 { diagonal } else { 1.0 };
    let surfaces =
        design.bodies.iter().map(|b| b.tessellate(tolerance).map_err(|e| format!("tessellate {}: {e}", b.name))).collect::<Result<Vec<_>, _>>()?;
    let ms = elapsed_ms(start);
    Ok((Loaded { part, design, measures, surfaces, tolerance, diagonal }, ms))
}

/// A failed part or stage: the error, then blanks to fill `width` cells.
fn error_cells(e: &str, width: usize) -> impl Iterator<Item = String> {
    std::iter::once(format!("error: {e}")).chain(std::iter::repeat_n(String::new(), width.saturating_sub(1)))
}

/// Run the corpus; with `write` (and every part) write [`FILE`] under the workspace `root`.
pub fn run(root: &Path, only: Option<&str>, write: bool) -> Result<(), String> {
    let mut parts = corpus::parts().map_err(|e| e.to_string())?;
    if let Some(name) = only {
        parts.retain(|p| p.name == name);
        if parts.is_empty() {
            return Err(format!("no corpus part named `{name}`"));
        }
    }
    let stages = stages();
    const PART_COLUMNS: &[&str] = &["part", "categories", "bodies", "faces", "triangles", "tol. mm", "build ms"];
    let mut header = PART_COLUMNS.to_vec();
    for s in &stages {
        header.extend_from_slice(s.headers());
    }
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut failures = 0;
    for part in &parts {
        eprintln!("corpus: {}", part.name);
        let mut row = vec![part.name.clone(), part.categories.join(", ")];
        match load(part) {
            Err(e) => {
                failures += 1;
                row.extend(error_cells(&e, header.len() - row.len()));
            }
            Ok((l, ms)) => {
                let faces: usize = l.measures.iter().map(|m| m.faces).sum();
                let triangles: usize = l.surfaces.iter().map(|s| s.triangles.len()).sum();
                row.extend([
                    l.design.bodies.len().to_string(),
                    faces.to_string(),
                    triangles.to_string(),
                    format!("{:.3}", l.tolerance),
                    format!("{ms:.0}"),
                ]);
                for s in &stages {
                    let width = s.headers().len();
                    match s.run(&l) {
                        Ok(cells) if cells.len() == width => row.extend(cells),
                        Ok(cells) => {
                            failures += 1;
                            row.extend(error_cells(&format!("stage gave {} cells for {width} columns", cells.len()), width));
                        }
                        Err(e) => {
                            failures += 1;
                            row.extend(error_cells(&e, width));
                        }
                    }
                }
            }
        }
        rows.push(row);
    }
    let table = markdown_table(&header, &rows);
    println!("{table}");
    if write && only.is_none() {
        let path = root.join(FILE);
        std::fs::write(&path, document(&table, &parts)).map_err(|e| format!("{}: {e}", path.display()))?;
        println!("corpus: wrote {FILE}");
    }
    if failures > 0 { Err(format!("corpus: {failures} part(s) or stage(s) failed")) } else { Ok(()) }
}

fn markdown_table(header: &[&str], rows: &[Vec<String>]) -> String {
    let mut s = format!("| {} |\n|{}\n", header.join(" | "), "---|".repeat(header.len()));
    for r in rows {
        let _ = writeln!(s, "| {} |", r.iter().map(|c| c.replace('|', "\\|")).collect::<Vec<_>>().join(" | "));
    }
    s
}

fn document(table: &str, parts: &[Part]) -> String {
    let mut s = String::from(
        "# Meshing test corpus\n\n\
         Generated by `cargo xtask corpus`; do not edit. The parts are our own: SolidCraft recipes in\n\
         `tests/corpus/<name>/recipe.json` and analytic shapes generated in code\n\
         (`crates/testkit/src/corpus.rs`). Each run rebuilds every part, tessellates it at a chordal\n\
         tolerance of 1e-3 of its bounding-box diagonal and runs every stage on it. Times are wall\n\
         clock on the machine that ran it; every other column is deterministic.\n\n\
         ## Results\n\n",
    );
    s.push_str(table);
    s.push_str("\n## Categories\n\n| category | stresses | parts |\n|---|---|---|\n");
    for (c, what) in CATEGORIES {
        let n = parts.iter().filter(|p| p.categories.iter().any(|x| x == c)).count();
        let _ = writeln!(s, "| {c} | {what} | {n} |");
    }
    s.push_str("\n## Parts\n\n");
    for p in parts {
        let _ = writeln!(s, "- `{}`: {}", p.name, p.description);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The analytic shapes come out of boundary preparation valid, with their genus.
    #[test]
    fn analytic_parts_are_valid_boundaries() {
        let expected =
            [("analytic-sphere", "0"), ("analytic-torus", "1"), ("analytic-cylinder", "0"), ("analytic-pipe", "1"), ("analytic-plate-hole", "1")];
        let parts = corpus::analytic_parts();
        for (name, genus) in expected {
            let part = parts.iter().find(|p| p.name == name).unwrap();
            let (loaded, _) = load(part).unwrap();
            let cells = BoundaryCheck.run(&loaded).unwrap();
            assert_eq!(cells[1], genus, "{name}: {cells:?}");
        }
    }

    /// The partitioned bracket's 30 cells mesh into one conforming hex mesh.
    #[test]
    fn partitioned_bracket_meshes_as_blocks() {
        let part = corpus::read_recipe(&corpus::corpus_dir().join("filleted-bracket-blocks")).unwrap();
        let (loaded, _) = load(&part).unwrap();
        assert_eq!(loaded.design.bodies.len(), 30);
        let cells = StructuredHex.run(&loaded).unwrap();
        assert!(cells[0].starts_with("blocks (30 cells"), "{cells:?}");
        assert!(cells[2].parse::<f64>().unwrap() > 0.5, "{cells:?}");
    }

    /// The constrained Delaunay tetrahedralization passes the gate on parts that exercise each
    /// of its stages: an assembly whose rounding needs a few-ulp repair, cavities refilled by
    /// gift-wrapping, and a holed plate.
    #[test]
    fn cdt_meshes_corpus_parts() {
        let mut parts = corpus::parts().unwrap();
        parts.retain(|p| ["bolted-assembly", "thin-wall-enclosure", "analytic-plate-hole"].contains(&p.name.as_str()));
        assert_eq!(parts.len(), 3);
        for part in &parts {
            let (loaded, _) = load(part).unwrap();
            let cells = TetCdt.run(&loaded).unwrap_or_else(|e| panic!("{}: {e}", part.name));
            assert_eq!(cells[3], "pass", "{}: {cells:?}", part.name);
            assert!(cells[0].parse::<usize>().unwrap() > 0);
        }
    }

    /// The tet mesh of a SolidCraft part, as `mesh.tet.generate` makes it (tet10, default
    /// options): every coordinate's bits, the elements and the boundary faces with their faces' names.
    fn tet_fingerprint(part: &Part) -> Vec<u64> {
        let (l, _) = load(part).unwrap();
        let mut out = Vec::new();
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let (s, soup) = tet_surface(mesh).unwrap();
            let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
            let o = QualityOptions { order: fieldcraft_mesh::tet::Order::Quadratic, ..QualityOptions::default() };
            let mut m = quality_mesh(&s, &cdt, soup.analytic(&fieldcraft_engine::analytic_surfaces(body)), &o).unwrap();
            soup.restore_tags(&mut m);
            out.extend(m.points.iter().flatten().map(|c| c.to_bits()));
            out.extend(m.elements.iter().flatten().map(|&v| u64::from(v)));
            // A face by its persistent name: SolidCraft may number the faces differently each run.
            let faces = body.faces();
            let name = |f: u32| {
                faces
                    .get(f as usize)
                    .map_or(0, |x| x.name.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)))
            };
            out.extend(m.boundary.iter().flat_map(|b| [u64::from(b.tet) << 8 | u64::from(b.opposite), name(b.face)]));
        }
        out
    }

    /// The same part gives the same tet mesh bit for bit on repeated runs and with 1 and 8
    /// threads (#115). SolidCraft lists the rounded box's faces in a different order from run to
    /// run; the mesher's canonical input order ([`CanonicalSoup`]) keeps the mesh the same.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn tet_meshes_do_not_depend_on_the_run_or_the_thread_count() {
        let mut parts = corpus::parts().unwrap();
        parts.retain(|p| ["rounded-box", "filleted-bracket"].contains(&p.name.as_str()));
        assert_eq!(parts.len(), 2);
        for part in &parts {
            let one = fieldcraft_parallel::with_threads(1, || tet_fingerprint(part));
            let many = fieldcraft_parallel::with_threads(8, || tet_fingerprint(part));
            let again = fieldcraft_parallel::with_threads(1, || tet_fingerprint(part));
            assert!(!one.is_empty());
            let first = |a: &[u64], b: &[u64]| a.iter().zip(b).position(|(x, y)| x != y);
            assert!(
                one == many && one == again,
                "{}: the tet mesh changed between runs ({} / {} / {} words, first difference {:?} / {:?})",
                part.name,
                one.len(),
                many.len(),
                again.len(),
                first(&one, &many),
                first(&one, &again)
            );
        }
    }
}
