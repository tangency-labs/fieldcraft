//! The corpus runner behind `cargo xtask corpus`: build every part of the meshing test corpus
//! ([`crate::corpus`]), tessellate it, run each [`Stage`] on it, check every tet mesh it makes
//! against the validity invariants ([`crate::validity`]), compare each part with its fixture
//! ([`crate::fixture`]) and write the results table to docs/meshing/corpus.md.
//!
//! The table is a fixed part columns block (bodies, faces, the boundary surface and the time to
//! build it) followed by one column group per stage. A mesher joins the corpus by implementing
//! [`Stage`] over a [`Loaded`] part (its tessellated boundary surfaces are ready in
//! `Loaded::surfaces`), naming the measures it freezes, and adding itself to [`stages`]. A stage
//! that fails on a part shows the error in its first column, and its fixture keeps the error's
//! kind; the run goes on.

use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

use crate::corpus::{self, CATEGORIES, Part};
use crate::fixture::{self, Fixture, Metric, Results};
use crate::validity;
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

/// What a stage gives for a part.
#[derive(Debug, Clone, Default)]
pub struct Outcome {
    /// One table cell per header.
    pub cells: Vec<String>,
    /// The outcome frozen in the part's fixture: `ok`, `ok: <method>`, or `no: <kind>` when the
    /// mesher refused the part.
    pub status: String,
    /// The deterministic measures frozen in the fixture (never times).
    pub values: Vec<(&'static str, Metric)>,
    /// The validity of the tet meshes the stage made.
    pub validity: Option<validity::Report>,
}

impl Outcome {
    fn ok(cells: Vec<String>, values: Vec<(&'static str, Metric)>) -> Outcome {
        Outcome { cells, status: "ok".into(), values, validity: None }
    }

    /// The mesher refused the part (an expected limit of the method, not a broken stage): the
    /// reason in the first cell, the time in the last.
    fn refused(body: &str, e: &fieldcraft_mesh::Error, width: usize, start: Instant) -> Outcome {
        let reason: String = e.to_string().chars().take(REASON_CHARS).collect();
        let mut cells = vec![format!("no: {body}: {reason}")];
        cells.resize(width.saturating_sub(1), String::new());
        cells.push(format!("{:.0}", elapsed_ms(start)));
        Outcome { cells, status: format!("no: {}", kind_of(e)), values: Vec::new(), validity: None }
    }
}

/// A stage that failed on a part: the error's kind (frozen) and its message (shown).
#[derive(Debug, Clone)]
pub struct Failure {
    pub kind: String,
    pub message: String,
}

impl Failure {
    fn new(kind: &str, message: String) -> Failure {
        Failure { kind: kind.into(), message }
    }

    fn mesh(body: &str, e: &fieldcraft_mesh::Error) -> Failure {
        Failure { kind: kind_of(e), message: format!("{body}: {e}") }
    }

    fn geometry(body: &str, e: &fieldcraft_geometry::Error) -> Failure {
        Failure::new("geometry", format!("{body}: {e}"))
    }
}

/// A mesher error's kind, as its message starts ("limit reached", "not supported yet", …).
fn kind_of(e: &fieldcraft_mesh::Error) -> String {
    e.kind().map_or_else(|| "model error".to_string(), |k| k.to_string())
}

/// A column group of the results table: a mesher, a check.
pub trait Stage {
    /// The stage's key in the fixture.
    fn name(&self) -> &'static str;
    fn headers(&self) -> &'static [&'static str];
    fn run(&self, part: &Loaded) -> Result<Outcome, Failure>;
}

/// The stages, in column order. Meshers add theirs here.
fn stages() -> Vec<Box<dyn Stage>> {
    vec![Box::new(SurfaceVolume), Box::new(BoundaryCheck), Box::new(TetCdt), Box::new(TetQualityStage), Box::new(StructuredHex), Box::new(OctreeHex)]
}

/// The boundary surface's enclosed volume against SolidCraft's measured volume and the closed
/// form: a surface that is open, inverted or misplaced shows here first.
struct SurfaceVolume;

impl Stage for SurfaceVolume {
    fn name(&self) -> &'static str {
        "surface"
    }

    fn headers(&self) -> &'static [&'static str] {
        &["enclosed mm³", "SolidCraft mm³", "rel. diff", "exact mm³"]
    }

    fn run(&self, l: &Loaded) -> Result<Outcome, Failure> {
        let enclosed: f64 = l.surfaces.iter().map(SurfaceMesh::enclosed_volume).sum();
        let measured: f64 = l.measures.iter().map(|m| m.volume).sum();
        let rel = if measured.abs() > 0.0 { (enclosed - measured) / measured } else { f64::NAN };
        let mut values = vec![("enclosed_mm3", Metric::Real(enclosed)), ("solidcraft_mm3", Metric::Real(measured))];
        if let Some(exact) = l.part.volume {
            values.push(("enclosed_vs_exact", Metric::Relative((enclosed - exact) / exact)));
        }
        Ok(Outcome::ok(
            vec![format!("{enclosed:.1}"), format!("{measured:.1}"), format!("{rel:.1e}"), l.part.volume.map_or("–".into(), |v| format!("{v:.1}"))],
            values,
        ))
    }
}

/// Boundary surface preparation (`fieldcraft_mesh::surface`): welded, closed, manifold, outward
/// and free of self-intersections, or the first defects found, which keep the part from being
/// volume meshed.
struct BoundaryCheck;

impl Stage for BoundaryCheck {
    fn name(&self) -> &'static str {
        "boundary"
    }

    fn headers(&self) -> &'static [&'static str] {
        &["boundary", "genus", "feature edges", "check ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Outcome, Failure> {
        let start = Instant::now();
        let (mut vertices, mut genus, mut features, mut total_genus, mut shells) = (0, Vec::new(), 0, 0, 0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
            let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).map_err(|e| Failure::mesh(&body.name, &e))?;
            vertices += s.points().len();
            genus.extend(s.shells().iter().map(|sh| if sh.void { format!("{} (void)", sh.genus) } else { sh.genus.to_string() }));
            total_genus += s.shells().iter().map(|sh| sh.genus as usize).sum::<usize>();
            shells += s.shells().len();
            features += s.feature_edges().len();
        }
        let ms = elapsed_ms(start);
        Ok(Outcome::ok(
            vec![format!("ok, {vertices} vertices"), genus.join(", "), features.to_string(), format!("{ms:.1}")],
            vec![
                ("vertices", Metric::Count(vertices)),
                ("shells", Metric::Count(shells)),
                ("genus", Metric::Count(total_genus)),
                ("feature_edges", Metric::Count(features)),
            ],
        ))
    }
}

/// The constrained Delaunay tetrahedralization (`fieldcraft_mesh::cdt`) of each body's boundary:
/// tetrahedra, Steiner points on boundary edges, cavities refilled by gift-wrapping, and the
/// validity gate on the rounded mesh (every tetrahedron positive, volume equal to the boundary's).
struct TetCdt;

impl Stage for TetCdt {
    fn name(&self) -> &'static str {
        "cdt"
    }

    fn headers(&self) -> &'static [&'static str] {
        &["CDT tets", "Steiner", "cavities", "tet gate", "CDT ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Outcome, Failure> {
        let start = Instant::now();
        let (mut tets, mut steiner, mut cavities, mut inadmissible) = (0, 0, 0, 0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let (s, _) = tet_surface(mesh).map_err(|e| Failure::mesh(&body.name, &e))?;
            let m = tetrahedralize(&s, &CdtOptions::default()).map_err(|e| Failure::mesh(&body.name, &e))?;
            let report = m.gate(s.volume()).map_err(|e| Failure::mesh(&format!("{}: rounded mesh", body.name), &e))?;
            if !report.volume_matches() {
                let message = format!("{}: rounded mesh: volume {} against the boundary's {}", body.name, report.volume, s.volume());
                return Err(Failure::new("meshing failed", message));
            }
            inadmissible += report.invalid.len() + report.undetermined.len();
            tets += m.tets.len();
            steiner += m.steiner.len();
            cavities += m.stats.cavities;
        }
        let ms = elapsed_ms(start);
        let gate = if inadmissible == 0 { "pass".to_string() } else { format!("{inadmissible} flat") };
        Ok(Outcome::ok(
            vec![tets.to_string(), steiner.to_string(), cavities.to_string(), gate, format!("{ms:.0}")],
            vec![
                ("tets", Metric::Count(tets)),
                ("steiner", Metric::Count(steiner)),
                ("cavities", Metric::Count(cavities)),
                ("flat", Metric::Count(inadmissible)),
            ],
        ))
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
/// above 170°) of the CDT and of the result, the validity gate, and the mesh validity invariants
/// checked on the test side ([`validity`]).
struct TetQualityStage;

impl Stage for TetQualityStage {
    fn name(&self) -> &'static str {
        "quality"
    }

    fn headers(&self) -> &'static [&'static str] {
        &["quality tets", "min dihedral °", "max dihedral °", "min radius ratio", "slivers CDT → mesh", "quality gate", "mesh validity", "quality ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Outcome, Failure> {
        let start = Instant::now();
        let (mut tets, mut lo, mut hi, mut rr, mut before, mut after) = (0, 180.0_f64, 0.0_f64, 1.0_f64, 0, 0);
        let mut checked = validity::Report::default();
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let (s, soup) = tet_surface(mesh).map_err(|e| Failure::mesh(&body.name, &e))?;
            let m = tetrahedralize(&s, &CdtOptions::default()).map_err(|e| Failure::mesh(&body.name, &e))?;
            let q = quality_mesh(&s, &m, soup.analytic(&fieldcraft_engine::analytic_surfaces(body)), &QualityOptions::default())
                .map_err(|e| Failure::mesh(&body.name, &e))?;
            let solids = s.shells().iter().filter(|sh| !sh.void).count();
            checked.absorb(&body.name, validity::check(&q.points, &q.elements, &q.boundary, solids));
            let r = &q.report;
            tets += r.output.tets;
            lo = lo.min(r.output.min_dihedral);
            hi = hi.max(r.output.max_dihedral);
            rr = rr.min(r.output.min_radius_ratio);
            before += r.input.slivers;
            after += r.output.slivers;
        }
        let ms = elapsed_ms(start);
        let valid = if checked.is_valid() { "ok".to_string() } else { format!("**{} violations**", checked.total()) };
        let mut outcome = Outcome::ok(
            vec![
                tets.to_string(),
                format!("{lo:.1}"),
                format!("{hi:.1}"),
                format!("{rr:.3}"),
                format!("{before} → {after}"),
                "pass".into(),
                valid,
                format!("{ms:.0}"),
            ],
            vec![
                ("tets", Metric::Count(tets)),
                ("min_dihedral_deg", Metric::Real(lo)),
                ("max_dihedral_deg", Metric::Real(hi)),
                ("min_radius_ratio", Metric::Real(rr)),
                ("slivers_cdt", Metric::Count(before)),
                ("slivers", Metric::Count(after)),
            ],
        );
        outcome.validity = Some(checked);
        Ok(outcome)
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
    fn name(&self) -> &'static str {
        "hex"
    }

    fn headers(&self) -> &'static [&'static str] {
        &["hex", "hex elements", "hex8 min SJ", "hex20 min SJ", "hex ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Outcome, Failure> {
        let start = Instant::now();
        let width = self.headers().len();
        let size = UniformSize(HEX_RELATIVE_SIZE * l.diagonal);
        let linear = HexOptions { chordal_tolerance: l.tolerance, ..HexOptions::default() };
        let quadratic = HexOptions { quadratic: true, ..linear.clone() };
        if l.part.categories.iter().any(|c| c == "partitioned") {
            return blocks(l, &size, &linear, &quadratic, width, start);
        }
        let (mut methods, mut elements, mut sj8, mut sj20) = (Vec::new(), 0usize, f64::INFINITY, Ok(f64::INFINITY));
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let surface_options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
            let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &surface_options)
                .map_err(|e| Failure::mesh(&body.name, &e))?;
            let faces: Vec<Option<Surface>> =
                body.face_surfaces().map_err(|e| Failure::geometry(&body.name, &e))?.into_iter().map(|f| f.map(analytic)).collect();
            let m_method = match hex_mesh(&s, &faces, &size, &linear) {
                Ok(m) => {
                    methods.push(m.report.method.name());
                    elements += m.element_count();
                    sj8 = sj8.min(m.report.min_scaled_jacobian);
                    m.report.method
                }
                Err(e) => return Ok(Outcome::refused(&body.name, &e, width, start)),
            };
            // The hex20 run goes straight to the method that meshed the body.
            let quadratic = HexOptions { method: m_method, ..quadratic.clone() };
            sj20 = match (sj20, hex_mesh(&s, &faces, &size, &quadratic)) {
                (Ok(a), Ok(m)) => Ok(a.min(m.report.min_scaled_jacobian)),
                (Err(e), _) => Err(e),
                (_, Err(e)) => Err((kind_of(&e), format!("{}: {}", body.name, e.to_string().chars().take(REASON_CHARS).collect::<String>()))),
            };
        }
        methods.sort_unstable();
        methods.dedup();
        let ms = elapsed_ms(start);
        let mut status = format!("ok: {}", methods.join(", "));
        let mut values = vec![("elements", Metric::Count(elements)), ("hex8_min_sj", Metric::Real(sj8))];
        let sj20 = match sj20 {
            Ok(v) => {
                values.push(("hex20_min_sj", Metric::Real(v)));
                format!("{v:.3}")
            }
            Err((kind, e)) => {
                let _ = write!(status, "; hex20 error: {kind}");
                format!("error: {e}")
            }
        };
        Ok(Outcome {
            cells: vec![methods.join(", "), elements.to_string(), format!("{sj8:.3}"), sj20, format!("{ms:.0}")],
            status,
            values,
            validity: None,
        })
    }
}

/// The bodies of a partitioned part as the cells of one conforming mesh: method, elements and
/// the hex8 and hex20 smallest scaled Jacobians.
fn blocks(l: &Loaded, size: &UniformSize, linear: &HexOptions, quadratic: &HexOptions, width: usize, start: Instant) -> Result<Outcome, Failure> {
    let mut surfaces = Vec::with_capacity(l.surfaces.len());
    let mut faces: Vec<Vec<Option<Surface>>> = Vec::with_capacity(l.surfaces.len());
    for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
        let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
        surfaces.push(BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).map_err(|e| Failure::mesh(&body.name, &e))?);
        faces.push(body.face_surfaces().map_err(|e| Failure::geometry(&body.name, &e))?.into_iter().map(|f| f.map(analytic)).collect());
    }
    let cells: Vec<Cell> =
        l.design.bodies.iter().zip(&surfaces).zip(&faces).map(|((b, surface), faces)| Cell { name: &b.name, surface, faces }).collect();
    let m = match hex_mesh_cells(&cells, size, linear) {
        Ok(m) => m,
        Err(e) => {
            let mut refused = Outcome::refused("blocks", &e, width, start);
            if let Some(first) = refused.cells.first_mut() {
                *first = format!("no: {}", e.to_string().chars().take(REASON_CHARS).collect::<String>());
            }
            return Ok(refused);
        }
    };
    let mut status = format!("ok: blocks ({} cells)", cells.len());
    let mut values = vec![
        ("elements", Metric::Count(m.mesh.element_count())),
        ("interfaces", Metric::Count(m.interfaces)),
        ("hex8_min_sj", Metric::Real(m.mesh.report.min_scaled_jacobian)),
    ];
    let sj20 = match hex_mesh_cells(&cells, size, quadratic) {
        Ok(q) => {
            values.push(("hex20_min_sj", Metric::Real(q.mesh.report.min_scaled_jacobian)));
            format!("{:.3}", q.mesh.report.min_scaled_jacobian)
        }
        Err(e) => {
            let _ = write!(status, "; hex20 error: {}", kind_of(&e));
            format!("error: {e}")
        }
    };
    let cells = vec![
        format!("blocks ({} cells, {} shared faces)", cells.len(), m.interfaces),
        m.mesh.element_count().to_string(),
        format!("{:.3}", m.mesh.report.min_scaled_jacobian),
        sj20,
        format!("{:.0}", elapsed_ms(start)),
    ];
    Ok(Outcome { cells, status, values, validity: None })
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
    fn name(&self) -> &'static str {
        "octree"
    }

    fn headers(&self) -> &'static [&'static str] {
        &["octree hexes", "min SJ", "mean SJ", "hex vol. diff", "octree ms"]
    }

    fn run(&self, l: &Loaded) -> Result<Outcome, Failure> {
        let start = Instant::now();
        let (mut hexes, mut min_sj, mut sum_sj, mut volume) = (0usize, f64::INFINITY, 0.0, 0.0);
        for (body, mesh) in l.design.bodies.iter().zip(&l.surfaces) {
            let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
            let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).map_err(|e| Failure::mesh(&body.name, &e))?;
            let options = OctreeOptions { max_elements: OCTREE_MAX_ELEMENTS, ..OctreeOptions::default() };
            let m = match octree_hex_mesh(&s, &UniformSize(HEX_RELATIVE_SIZE * l.diagonal), &options) {
                Ok(m) => m,
                Err(e) => return Ok(Outcome::refused(&body.name, &e, self.headers().len(), start)),
            };
            hexes += m.hexes.len();
            min_sj = min_sj.min(m.report.min_scaled_jacobian);
            sum_sj += m.report.mean_scaled_jacobian * m.hexes.len() as f64;
            volume += m.report.volume;
        }
        let measured: f64 = l.measures.iter().map(|m| m.volume).sum();
        let rel = if measured.abs() > 0.0 { (volume - measured) / measured } else { f64::NAN };
        let mean = sum_sj / hexes.max(1) as f64;
        let ms = elapsed_ms(start);
        Ok(Outcome::ok(
            vec![hexes.to_string(), format!("{min_sj:.3}"), format!("{mean:.3}"), format!("{rel:.1e}"), format!("{ms:.0}")],
            vec![
                ("hexes", Metric::Count(hexes)),
                ("min_sj", Metric::Real(min_sj)),
                ("mean_sj", Metric::Real(mean)),
                ("volume_mm3", Metric::Real(volume)),
            ],
        ))
    }
}

fn load(part: &Part) -> Result<(Loaded<'_>, f64), Failure> {
    let start = Instant::now();
    let design = corpus::load(part).map_err(|e| Failure::new("geometry", e.to_string()))?;
    if let Some(e) = design.errors.first() {
        return Err(Failure::new("feature failed", format!("feature {} failed: {}", e.feature, e.message)));
    }
    let measures = design
        .bodies
        .iter()
        .map(|b| b.measure().map_err(|e| Failure::geometry(&format!("measure {}", b.name), &e)))
        .collect::<Result<Vec<_>, _>>()?;
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for m in &measures {
        for i in 0..3 {
            lo[i] = lo[i].min(m.bbox_min[i]);
            hi[i] = hi[i].max(m.bbox_max[i]);
        }
    }
    let diagonal = (0..3).map(|i| (hi[i] - lo[i]).powi(2)).sum::<f64>().sqrt();
    let tolerance = RELATIVE_TOLERANCE * if diagonal.is_finite() && diagonal > 0.0 { diagonal } else { 1.0 };
    let surfaces = design
        .bodies
        .iter()
        .map(|b| b.tessellate(tolerance).map_err(|e| Failure::geometry(&format!("tessellate {}", b.name), &e)))
        .collect::<Result<Vec<_>, _>>()?;
    let ms = elapsed_ms(start);
    Ok((Loaded { part, design, measures, surfaces, tolerance, diagonal }, ms))
}

/// A failed part or stage: the error, then blanks to fill `width` cells.
fn error_cells(e: &str, width: usize) -> impl Iterator<Item = String> {
    std::iter::once(format!("error: {e}")).chain(std::iter::repeat_n(String::new(), width.saturating_sub(1)))
}

/// One part's run: its table row, what its fixture freezes, and the validity of its tet meshes.
pub struct PartRun {
    pub row: Vec<String>,
    pub results: Results,
    pub validity: validity::Report,
}

/// The part columns of the table, before the stages'.
const PART_COLUMNS: &[&str] = &["part", "categories", "bodies", "faces", "triangles", "tol. mm", "build ms"];
/// The fixture's key for building and tessellating the part.
const BUILD: &str = "build";

/// Build, tessellate and mesh one part through every stage.
fn run_part(part: &Part, stages: &[Box<dyn Stage>]) -> PartRun {
    let mut row = vec![part.name.clone(), part.categories.join(", ")];
    let mut results = Results::default();
    let mut checked = validity::Report::default();
    let width: usize = PART_COLUMNS.len() + stages.iter().map(|s| s.headers().len()).sum::<usize>();
    let (l, ms) = match load(part) {
        Ok(loaded) => loaded,
        Err(f) => {
            results.stages.insert(BUILD.into(), format!("error: {}", f.kind));
            results.stages.extend(stages.iter().map(|s| (s.name().to_string(), "not run".to_string())));
            row.extend(error_cells(&f.message, width - row.len()));
            return PartRun { row, results, validity: checked };
        }
    };
    let faces: usize = l.measures.iter().map(|m| m.faces).sum();
    let triangles: usize = l.surfaces.iter().map(|s| s.triangles.len()).sum();
    results.stages.insert(BUILD.into(), "ok".into());
    for (k, n) in [("bodies", l.design.bodies.len()), ("faces", faces), ("triangles", triangles)] {
        results.values.insert(format!("{BUILD}.{k}"), Metric::Count(n));
    }
    row.extend([l.design.bodies.len().to_string(), faces.to_string(), triangles.to_string(), format!("{:.3}", l.tolerance), format!("{ms:.0}")]);
    for s in stages {
        let width = s.headers().len();
        match s.run(&l) {
            Ok(o) if o.cells.len() == width => {
                results.stages.insert(s.name().into(), o.status);
                results.values.extend(o.values.into_iter().map(|(k, m)| (format!("{}.{k}", s.name()), m)));
                if let Some(v) = o.validity {
                    checked.absorb(s.name(), v);
                }
                row.extend(o.cells);
            }
            Ok(o) => {
                let e = format!("stage gave {} cells for {width} columns", o.cells.len());
                results.stages.insert(s.name().into(), format!("error: {e}"));
                row.extend(error_cells(&e, width));
            }
            Err(f) => {
                results.stages.insert(s.name().into(), format!("error: {}", f.kind));
                row.extend(error_cells(&f.message, width));
            }
        }
    }
    PartRun { row, results, validity: checked }
}

/// Options of one corpus run.
#[derive(Debug, Clone, Default)]
pub struct Options {
    /// Run only these parts (all when empty).
    pub parts: Vec<String>,
    /// Rewrite the fixtures of these parts (`all` for every part run).
    pub bless: Vec<String>,
    /// Write [`FILE`] (only after a run of every part).
    pub write: bool,
    /// Parts meshed at once (the meshers are mostly sequential); 0 or 1 runs them one after
    /// another, which keeps the table's times those of a part alone.
    pub jobs: usize,
}

/// Run every part through the stages, `jobs` parts at a time, the results in the parts' order.
/// Each worker thread has its own stages; a part's results do not depend on which thread ran it.
fn run_parts(parts: &[Part], jobs: usize) -> Vec<PartRun> {
    let next = std::sync::atomic::AtomicUsize::new(0);
    let done: std::sync::Mutex<Vec<(usize, PartRun)>> = std::sync::Mutex::new(Vec::with_capacity(parts.len()));
    std::thread::scope(|scope| {
        for _ in 0..jobs.clamp(1, parts.len().max(1)) {
            scope.spawn(|| {
                let stages = stages();
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(part) = parts.get(i) else { break };
                    let start = Instant::now();
                    let run = run_part(part, &stages);
                    eprintln!("corpus: {} ({:.1} s)", part.name, start.elapsed().as_secs_f64());
                    if let Ok(mut d) = done.lock() {
                        d.push((i, run));
                    }
                }
            });
        }
    });
    let mut done = done.into_inner().unwrap_or_default();
    done.sort_by_key(|(i, _)| *i);
    done.into_iter().map(|(_, run)| run).collect()
}

/// Run the corpus under the workspace `root`: print the table, check every part against its
/// fixture (or rewrite the fixtures named by `bless`) and check every tet mesh's validity; with
/// `write` (and every part) write [`FILE`]. Fails when a part's results moved from its fixture,
/// a part has no fixture, a fixture names no part, or a tet mesh breaks a validity invariant (a
/// blessed fixture never hides that).
pub fn run(root: &Path, options: &Options) -> Result<(), String> {
    let corpus_dir = root.join("tests").join("corpus");
    let mut parts = corpus::read_recipes(&corpus_dir).map_err(|e| e.to_string())?;
    parts.extend(corpus::analytic_parts());
    for want in options.parts.iter().chain(options.bless.iter().filter(|b| *b != "all")) {
        if !parts.iter().any(|p| p.name == *want) {
            return Err(format!("no corpus part named `{want}`"));
        }
    }
    let every = options.parts.is_empty();
    if !every {
        parts.retain(|p| options.parts.contains(&p.name));
    }
    let mut header = PART_COLUMNS.to_vec();
    for s in &stages() {
        header.extend_from_slice(s.headers());
    }
    let runs = run_parts(&parts, options.jobs);
    if runs.len() != parts.len() {
        return Err(format!("corpus: {} of {} parts ran", runs.len(), parts.len()));
    }
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut problems = Vec::new();
    let platform = fixture::platform();
    for (part, run) in parts.iter().zip(runs) {
        let path = fixture::path(&corpus_dir, &part.name);
        let old = fixture::read(&path)?;
        if options.bless.iter().any(|b| *b == part.name || b == "all") {
            match Fixture::bless(old.as_ref(), &run.results, &platform) {
                None => eprintln!("corpus: {}: unchanged on {platform}", part.name),
                Some(new) => {
                    if let Some(dir) = path.parent() {
                        std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
                    }
                    std::fs::write(&path, fixture::render(&new)?).map_err(|e| format!("{}: {e}", path.display()))?;
                    let dropped: Vec<&String> = old.iter().flat_map(|o| o.platforms.keys()).filter(|p| !new.platforms.contains_key(*p)).collect();
                    let dropped =
                        if dropped.is_empty() { String::new() } else { format!("; the results of {dropped:?} are stale now: bless it there too") };
                    eprintln!("corpus: {}: blessed on {platform}{dropped}", part.name);
                }
            }
        } else {
            match old {
                None => problems.push(format!("{}: no fixture; bless the part (`cargo xtask corpus --bless {}`)", part.name, part.name)),
                Some(f) => {
                    let moved = fixture::compare(&run.results, &f.expected(&platform));
                    // A fixture without this platform's results may differ by platform alone.
                    let hint = if !moved.is_empty() && !f.platforms.contains_key(&platform) {
                        let on: Vec<&String> = f.platforms.keys().collect();
                        format!(" (blessed on {on:?} only: if {platform} differs by platform alone, bless it here to record that)")
                    } else {
                        String::new()
                    };
                    problems.extend(moved.into_iter().map(|p| format!("{}: {p}{hint}", part.name)));
                }
            }
        }
        if !run.validity.is_valid() {
            problems.push(format!("{}: invalid tet mesh: {}", part.name, run.validity.summary()));
        }
        rows.push(run.row);
    }
    if every {
        problems.extend(stale_fixtures(&corpus_dir, &parts)?);
    }
    let table = markdown_table(&header, &rows);
    println!("{table}");
    if options.write && every {
        let path = root.join(FILE);
        std::fs::write(&path, document(&table, &parts)).map_err(|e| format!("{}: {e}", path.display()))?;
        println!("corpus: wrote {FILE}");
    }
    if problems.is_empty() { Ok(()) } else { Err(format!("corpus: {} problem(s):\n  {}", problems.len(), problems.join("\n  "))) }
}

/// Fixtures under `dir` whose part is gone.
fn stale_fixtures(dir: &Path, parts: &[Part]) -> Result<Vec<String>, String> {
    let rd = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let mut stale: Vec<String> = rd
        .flatten()
        .filter(|e| e.path().join("fixture.json").is_file())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|name| !parts.iter().any(|p| p.name == *name))
        .map(|name| format!("{name}: a fixture for no corpus part; remove tests/corpus/{name}/fixture.json"))
        .collect();
    stale.sort();
    Ok(stale)
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
            let cells = BoundaryCheck.run(&loaded).unwrap().cells;
            assert_eq!(cells[1], genus, "{name}: {cells:?}");
        }
    }

    /// The partitioned bracket's 30 cells mesh into one conforming hex mesh.
    #[test]
    fn partitioned_bracket_meshes_as_blocks() {
        let part = corpus::read_recipe(&corpus::corpus_dir().join("filleted-bracket-blocks")).unwrap();
        let (loaded, _) = load(&part).unwrap();
        assert_eq!(loaded.design.bodies.len(), 30);
        let cells = StructuredHex.run(&loaded).unwrap().cells;
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
            let cells = TetCdt.run(&loaded).unwrap_or_else(|e| panic!("{}: {}", part.name, e.message)).cells;
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
