use crate::schema::ProblemKind;
use std::borrow::Cow;
use std::collections::BTreeMap;

use fieldcraft_geometry::SurfaceMesh;
use fieldcraft_mesh::cdt::{CdtOptions, tetrahedralize};
use fieldcraft_mesh::delaunay::Triangulation;
use fieldcraft_mesh::quality::{self, QualityReport};
use fieldcraft_mesh::surface::{BoundarySurface, Defect, FeatureKind, SurfaceOptions};
use fieldcraft_mesh::tet::canonical::CanonicalSoup;
use fieldcraft_mesh::tet::geometry::AnalyticSurface;
use fieldcraft_mesh::tet::improve::ImproveOptions;
use fieldcraft_mesh::tet::{Order, QualityMesh, QualityOptions, TetQuality, quality_mesh, raw_mesh};
use fieldcraft_model::{ElementKind, GeometryRef, MAX_ELEMENTS};
use serde_json::{Value, json};

use super::{label, labels, listed, optional_number, points, xyz};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const IMPORT_ORPHAN: CommandSpec = CommandSpec {
    id: "mesh.import.orphan",
    label: "Import Orphan Mesh",
    params: r#"{"path": "part.msh", "unmatched": "refuse" | "shell" (default refuse)} — imports ASCII Gmsh MSH 4.1, preserving node/element tags; the file has no unit declaration, so its coordinates use the current input length unit (`coordinate_unit` in the reply). Replaces the whole mesh and its geometry tags while keeping geometry and model definitions; refuses existing geometry/label sets, node-specific step limits, reference-point label collisions and physical-set name collisions rather than silently retargeting them (import into a new model or remove those selections first). Physical volumes become element sets, surfaces matching volume faces become face sets (an internal surface selects both adjacent element faces), curves become edge sets and physical points become node sets. Surface-only meshes keep their surface elements. Unmatched volume-boundary surfaces refuse unless `unmatched: shell` explicitly keeps intentional shells; unmatched curves still refuse. Supports rod2/rod3, tri3/tri6, quad4/quad8, tet4/tet10, wedge6, hex8/hex20. Refuses binary, partitioned, unsupported types and invalid physical names. Replies nodes, elements by kind, sets with membership counts, skipped_sections, replaced mesh counts, and results made stale. One undo step"#,
    form: Some(super::forms::IMPORT_ORPHAN),
    kind: CommandKind::Model,
    requires: &[],
    run: super::orphan::import,
};

pub const NODE_ADD: CommandSpec = CommandSpec {
    id: "mesh.node.add",
    label: "Add Node",
    params: r#"{"label": 1, "xyz": [x, y, z]} — adds the node, or moves it if the label exists"#,
    form: Some(super::forms::NODE_ADD),
    kind: CommandKind::Model,
    requires: &[],
    run: node_add,
};

pub const NODE_REMOVE: CommandSpec = CommandSpec {
    id: "mesh.node.remove",
    label: "Remove Node",
    params: r#"{"label": 1} — only a node no element uses"#,
    form: Some(super::forms::LABEL),
    kind: CommandKind::Model,
    requires: &[Requirement::Nodes],
    run: node_remove,
};

pub const ELEMENT_ADD: CommandSpec = CommandSpec {
    id: "mesh.element.add",
    label: "Add Element",
    params: r#"{"label": 1, "kind": "rod2|beam2|tri3|tri6|quad4|quad8|tet4|tet10|wedge6|pyramid5|hex8|hex20", "nodes": [labels]} — adds or replaces"#,
    form: Some(super::forms::ELEMENT_ADD),
    kind: CommandKind::Model,
    requires: &[Requirement::Nodes],
    run: element_add,
};

pub const ELEMENT_REMOVE: CommandSpec = CommandSpec {
    id: "mesh.element.remove",
    label: "Remove Element",
    params: r#"{"label": 1}"#,
    form: Some(super::forms::LABEL),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: element_remove,
};

pub const QUALITY: CommandSpec = CommandSpec {
    id: "mesh.quality",
    label: "Mesh Quality",
    params: r#"{"elements": [labels] (optional, default all), "bins": 10 (1..=1000)} — certified validity (det J > 0 everywhere), volume, and min/max/mean/worst/histogram per metric (angles in degrees; histograms over each metric's fixed range, out-of-range values in the end bins)"#,
    form: Some(super::forms::QUALITY),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: quality_report,
};

pub const QUALITY_MARK: CommandSpec = CommandSpec {
    id: "mesh.quality.mark",
    label: "Mark Bad Elements",
    params: r#"{"limits"?: {"scaled_jacobian": 0.2, "min_dihedral": 10, "max_dihedral": 165, "aspect_ratio": 10, "min_angle": 15, "max_angle": 165} (the defaults; given limits are merged over them, null drops one; any metric of mesh.quality may be limited: a floor for those where higher is better, a ceiling for max_dihedral, aspect_ratio, edge_ratio, max_angle and skew; angles in degrees), "elements"?: [labels] (default all)} — after meshing, marks the elements that fail a limit (and every solid not proven valid) as the element set `quality_marked`, replaced on every run and removed when none fail; the viewport and `fieldcraft-cli snapshot` draw its elements in the marker colour. Replies `summary`, `measured`, `marked`, `invalid`, the `limits` used and `elements` (worst first: not valid, then lowest scaled Jacobian; at most 1000, `truncated` says when more): each `label`, `kind`, `body`, `at` (its centroid), `valid` and `failed` ({metric, value, limit} per limit it fails). One undo step"#,
    form: Some(super::forms::QUALITY_MARK),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: quality_mark,
};

pub const EXPORT_STL: CommandSpec = CommandSpec {
    id: "mesh.export.stl",
    label: "Export STL",
    params: r#"{"path": "skin.stl", "set": "top", "format": "binary|ascii"} — writes the mesh's outer surface as STL: every solid element face no other solid shares, and the top side (face 1) of every shell; `set` (optional) writes only that set's element faces instead; `format` defaults to binary. Quadrilaterals are split in two, normals point out of the elements, and coordinates are rounded to single precision (STL stores f32). Replies {path, triangles, bytes}"#,
    form: Some(super::forms::EXPORT_STL),
    kind: CommandKind::Session,
    requires: &[Requirement::Mesh],
    run: export_stl,
};

pub const DELAUNAY_POINTS: CommandSpec = CommandSpec {
    id: "mesh.delaunay.points",
    label: "Delaunay of Points",
    params: r#"{"points": [[x, y, z], ...] (at most 20000)} — the exact Delaunay tetrahedralization (a debugging report; the model is not changed): {"dimension", "vertices", "duplicates": [[point, kept point], ...], "tets": [[i, j, k, l], ...]} with 0-based point indices in tet4 order; collinear or coplanar points give no tets"#,
    form: Some(super::forms::DELAUNAY_POINTS),
    kind: CommandKind::Model,
    requires: &[],
    run: delaunay_points,
};

/// The most points `mesh.delaunay.points` takes (its report lists every tetrahedron).
const DELAUNAY_MAX_POINTS: usize = 20_000;

fn delaunay_points(_: &mut Session, p: &Value) -> Result<Value> {
    let points = points(p, "points", DELAUNAY_MAX_POINTS)?;
    let tri = Triangulation::new(&points)?;
    tri.verify()?;
    let duplicates: Vec<[u32; 2]> = (0..points.len() as u32).filter_map(|v| tri.duplicate_of(v).map(|d| [v, d])).collect();
    let tets: Vec<[u32; 4]> = tri.finite_tets().map(|(_, vs)| vs).collect();
    Ok(json!({
        "dimension": tri.dimension(),
        "vertices": points.len() - duplicates.len(),
        "duplicates": duplicates,
        "tets": tets,
    }))
}

pub const SURFACE_CHECK: CommandSpec = CommandSpec {
    id: "mesh.surface.check",
    label: "Check Boundary Surface",
    params: r#"{"body": "Body1", "tolerance": 0.05} or {"surface": 0}, optional "weld": 1e-9 (0 or 1e-15..=1e-3 of the bounding-box diagonal), "weld_distance": length (default: the surface's own, twice SolidCraft's linear tolerance for a tessellated body, 0 for an imported mesh), "feature_angle": 45 (degrees, (0, 180]) — welds the tessellated body or imported surface and checks it can bound a volume mesh: closed, manifold, outward, no self-intersections; reports shells, genus, volume and feature edges, or the located defects"#,
    form: Some(super::forms::SURFACE_CHECK),
    kind: CommandKind::Model,
    requires: &[Requirement::Geometry],
    run: surface_check,
};

pub const TET_GENERATE: CommandSpec = CommandSpec {
    id: "mesh.tet.generate",
    label: "Generate Tet Mesh",
    params: r#"{"body": "Body1" (default: every body of the open design) or "surface": 0 (an imported surface), "tolerance": 0.05 (chordal tolerance of the tessellation, a length; default 1e-3 of the body's bounding-box diagonal), the boundary options of mesh.surface.check ("weld", "weld_distance", "feature_angle"), optional "size": length (target edge length; default a twentieth of the bounding-box diagonal), "order": "linear" | "quadratic" (tet4 or tet10; quadratic mid-edge nodes lie on the B-rep faces), "quality": true (false: the raw constrained Delaunay tetrahedralization), "quality_target": 40 (degrees: tetrahedra whose smallest dihedral-angle sine is at least sin of this are left alone by improvement, 5..=60), "curvature_angle": 22.5 (degrees per element along curved boundary, 0 for none, at most 90), "grading": 0.3 (largest size change per unit length, 0.05..=2; default the model's growth rate, mesh.size.growth), "max_passes": 20 (improvement passes, 0..=200), "max_points": 1000000 (1..=20000000), "max_steiner": 50000 (0..=10000000), "max_tets": 5000000 (1..=50000000) — the constrained Delaunay tetrahedralization of the welded boundary (every boundary triangle a union of tet faces), refined to the size field with its boundary vertices on the B-rep faces, smoothed and improved (Delaunay refinement, ODT smoothing, Klingner–Shewchuk improvement), with a quality report. Replaces the mesh of each body it meshes (its elements and the nodes only they use; other bodies' meshes stay; the reply's `replaced`, `stale_sets` and `results_cleared` say what went), adds the elements and their nodes (labels after the largest existing ones), tags the elements with the body and each boundary element face with its geometry face, and fails, changing nothing, if a surface is invalid or a mesh does not pass the validity gate. Edge seeds (mesh.seed.edge) divide their edges at the seeded nodes, and the size field grades from them; size controls (mesh.size.set) refine it, graded at the model's growth rate (mesh.size.growth) unless `grading` is given. Reports the totals and, under `bodies`, each body's counts, size, quality report and `seeds` (per seeded edge, the mesh edges along it)"#,
    form: Some(super::forms::TET),
    kind: CommandKind::Model,
    requires: &[Requirement::Geometry],
    run: tet_generate,
};

/// The most Steiner points and tetrahedra `mesh.tet.generate` may be asked for.
const TET_MAX_STEINER: u64 = 10_000_000;
const TET_MAX_TETS: u64 = 50_000_000;
const TET_MAX_POINTS: u64 = 20_000_000;
const TET_MAX_PASSES: u64 = 200;

/// An optional whole-number parameter in `range`.
pub(super) fn optional_count(p: &Value, key: &str, default: u64, range: std::ops::RangeInclusive<u64>) -> Result<usize> {
    let v = match p.get(key) {
        None => default,
        Some(v) => v.as_u64().filter(|n| range.contains(n)).ok_or_else(|| EngineError::param(key, whole_in(&range)))?,
    };
    usize::try_from(v).map_err(|_| EngineError::param(key, whole_in(&range)))
}

/// A whole number in `range`, as a problem's kind.
pub(super) fn whole_in(range: &std::ops::RangeInclusive<u64>) -> ProblemKind {
    let int = |n: u64| i64::try_from(n).unwrap_or(i64::MAX);
    ProblemKind::NotAWholeNumber { min: int(*range.start()), max: int(*range.end()) }
}

/// The boundary options of `mesh.surface.check`, `mesh.tet.generate` and the hex meshers.
pub(super) fn surface_options(p: &Value, mesh: &SurfaceMesh) -> Result<SurfaceOptions> {
    let defaults = SurfaceOptions::default();
    Ok(SurfaceOptions {
        weld_tolerance: optional_number(p, "weld", defaults.weld_tolerance)?,
        weld_distance: optional_number(p, "weld_distance", mesh.weld_distance)?,
        feature_angle: optional_number(p, "feature_angle", defaults.feature_angle)?,
    })
}

/// An optional number in `range`.
pub(super) fn number_in(p: &Value, key: &str, default: f64, range: std::ops::RangeInclusive<f64>) -> Result<f64> {
    let v = optional_number(p, key, default)?;
    if range.contains(&v) { Ok(v) } else { Err(EngineError::param(key, ProblemKind::NotInRange { min: *range.start(), max: *range.end() })) }
}

/// The quality options of `mesh.tet.generate`, and whether the quality stage runs at all.
pub(super) fn quality_options(p: &Value, max_tets: u64) -> Result<(QualityOptions, bool)> {
    let d = QualityOptions::default();
    let size = match p.get("size") {
        None => None,
        Some(v) => Some(
            v.as_f64()
                .filter(|h| *h > 0.0 && h.is_finite() && *h <= fieldcraft_model::MAX_COORD)
                .ok_or_else(|| EngineError::param("size", ProblemKind::NotPositive))?,
        ),
    };
    let order = match p.get("order") {
        None => Order::Linear,
        Some(v) => match v.as_str() {
            Some("linear") => Order::Linear,
            Some("quadratic") => Order::Quadratic,
            _ => return Err(EngineError::param("order", ProblemKind::NotAChoice)),
        },
    };
    let target = number_in(p, "quality_target", 40.0, 5.0..=60.0)?;
    let improve_options = ImproveOptions {
        max_passes: optional_count(p, "max_passes", d.improve_options.max_passes as u64, 0..=TET_MAX_PASSES)?,
        good: fieldcraft_mesh::tet::geometry::sin_cos(target.to_radians()).0,
        ..d.improve_options
    };
    let quality = super::opt::<bool>(p, "quality")?.unwrap_or(true);
    if !quality && order == Order::Quadratic {
        return Err(EngineError::Params("quadratic elements need the quality stage (`quality`: true)".into()));
    }
    let options = QualityOptions {
        size,
        curvature_angle: number_in(p, "curvature_angle", d.curvature_angle.to_degrees(), 0.0..=90.0)?.to_radians(),
        grading: number_in(p, "grading", d.grading, 0.05..=2.0)?,
        improve_options,
        order,
        max_tets: optional_count(p, "max_tets", max_tets, 1..=TET_MAX_TETS)?,
        max_points: optional_count(p, "max_points", d.max_points as u64, 1..=TET_MAX_POINTS)?,
        ..d
    };
    Ok((options, quality))
}

/// The exact surfaces of a body's faces, for the quality stage's projections (none when the
/// kernel cannot give them: the stage then keeps to the tessellation).
pub fn analytic_surfaces(shape: &fieldcraft_geometry::Shape) -> Vec<Option<AnalyticSurface>> {
    shape.face_surfaces().unwrap_or_default().into_iter().map(|f| f.map(analytic_surface)).collect()
}

fn analytic_surface(f: fieldcraft_geometry::FaceSurface) -> AnalyticSurface {
    use fieldcraft_geometry::FaceSurface as F;
    match f {
        F::Plane { normal, point } => AnalyticSurface::Plane { normal, point },
        F::Cylinder { origin, axis, radius } => AnalyticSurface::Cylinder { origin, axis, radius },
        F::Cone { apex, axis, half_angle } => AnalyticSurface::Cone { apex, axis, half_angle },
        F::Sphere { center, radius } => AnalyticSurface::Sphere { center, radius },
        F::Torus { center, axis, major, minor } => AnalyticSurface::Torus { center, axis, major, minor },
    }
}

/// The B-rep pieces of `shape`'s merged faces (virtual topology) for a projection target built
/// on its tessellation whose triangles lie on the B-rep faces `tri_piece`
/// ([`fieldcraft_geometry::Shape::tessellate_pieces`]): a merged face on no one surface then
/// projects onto each piece's own; none when no face is merged.
pub(crate) fn pieces(shape: &fieldcraft_geometry::Shape, tri_piece: &[u32]) -> fieldcraft_mesh::tet::geometry::Pieces {
    if !shape.has_merged_faces() {
        return fieldcraft_mesh::tet::geometry::Pieces::default();
    }
    let analytic = shape.brep_face_surfaces().unwrap_or_default().into_iter().map(|f| f.map(analytic_surface)).collect();
    fieldcraft_mesh::tet::geometry::Pieces { tri_piece: tri_piece.to_vec(), analytic }
}

pub(super) fn quality_json(q: &TetQuality) -> Value {
    json!({
        "tets": q.tets,
        "min_dihedral": q.min_dihedral,
        "max_dihedral": q.max_dihedral,
        "min_radius_ratio": q.min_radius_ratio,
        "mean_radius_ratio": q.mean_radius_ratio,
        "min_scaled_jacobian": q.min_scaled_jacobian,
        "slivers": q.slivers,
        "inadmissible": q.inadmissible,
    })
}

/// One body's tetrahedra, ready to add, with what its report needs of the surface and the CDT.
pub(super) struct BodyTets {
    pub(super) name: String,
    pub(super) face_names: Vec<String>,
    pub(super) mesh: QualityMesh,
    steiner: usize,
    cavities: usize,
    surface_volume: f64,
}

impl BodyTets {
    pub(super) fn new(name: String, face_names: Vec<String>, mesh: QualityMesh, steiner: usize, cavities: usize, surface_volume: f64) -> Self {
        BodyTets { name, face_names, mesh, steiner, cavities, surface_volume }
    }
}

/// Mesh one body or imported surface; the model is not touched. With `several` bodies an error
/// names the body.
pub(super) fn mesh_body(
    design: Option<&fieldcraft_geometry::Design>,
    p: &Value,
    input: super::hex::Input,
    cdt: &CdtOptions,
    options: &QualityOptions,
    quality: bool,
    several: bool,
) -> Result<BodyTets> {
    let mesh = || -> Result<BodyTets> {
        // The exact faces of a design body; an imported surface keeps to its triangles.
        let analytic = match p.get("surface") {
            None => design.and_then(|d| d.body(&input.name)).map(analytic_surfaces).unwrap_or_default(),
            Some(_) => Vec::new(),
        };
        fieldcraft_model::progress::stage(&format!("{}: boundary", input.name), None, "");
        // The mesh depends on the geometry only, not on the order the kernel lists it in.
        let soup = CanonicalSoup::new(&input.mesh.positions, &input.mesh.triangles, &input.mesh.tri_face)?;
        let surface = BoundarySurface::new(&soup.positions, &soup.triangles, &soup.tri_face, &surface_options(p, &input.mesh)?)?;
        fieldcraft_model::progress::stage(&format!("{}: Delaunay", input.name), None, "");
        let tets = tetrahedralize(&surface, cdt)?;
        let mut mesh = if quality {
            let frozen_faces = options.frozen_faces.iter().filter_map(|&f| soup.canonical_tag(f)).collect();
            let options = QualityOptions { seeds: input.seeds.clone(), frozen_faces, ..options.clone() };
            quality_mesh(&surface, &tets, soup.analytic(&analytic), &options)?
        } else {
            raw_mesh(&surface, &tets)?
        };
        soup.restore_tags(&mut mesh);
        Ok(BodyTets {
            name: input.name.clone(),
            face_names: input.face_names.clone(),
            mesh,
            steiner: tets.steiner.len(),
            cavities: tets.stats.cavities,
            surface_volume: surface.volume(),
        })
    };
    mesh().map_err(|e| if several { EngineError::context(format!("body `{}`", input.name), e) } else { e })
}

fn tet_generate(s: &mut Session, p: &Value) -> Result<Value> {
    let (options, quality) = quality_options(p, s.limits.max_tets as u64)?;
    let cdt_options =
        CdtOptions { max_steiner: optional_count(p, "max_steiner", s.limits.max_steiner as u64, 0..=TET_MAX_STEINER)?, max_tets: options.max_tets };
    // Every body is meshed before the model changes; with several, an error names its body.
    let inputs = super::hex::inputs(s, p)?;
    let sizing = super::sizes::sizing(s)?;
    let several = inputs.len() > 1;
    let mut meshed = Vec::with_capacity(inputs.len());
    for input in inputs {
        let options = sizing.tet_options(&input.name, p, &options);
        meshed.push(mesh_body(s.design.as_deref(), p, input, &cdt_options, &options, quality, several)?);
    }

    // Replace the bodies' meshes (replace.rs); other bodies and the definitions stay.
    let cleared = super::replace::clear(&mut s.model, meshed.iter().map(|b| b.name.as_str()));
    let mut bodies = Vec::with_capacity(meshed.len());
    for b in &meshed {
        let (first_node, first_element) = add_tets(&mut s.model, &b.name, &b.face_names, &b.mesh)?;
        let r = &b.mesh.report;
        bodies.push(json!({
            "body": b.name,
            "nodes": b.mesh.points.len(),
            "elements": b.mesh.elements.len(),
            "first_node": first_node,
            "first_element": first_element,
            "steiner": b.steiner,
            "cavities": b.cavities,
            "boundary_faces": b.mesh.boundary.len(),
            "volume": r.volume,
            "boundary_volume": r.boundary_volume,
            "surface_volume": b.surface_volume,
            "size": r.size,
            "quality": {
                "cdt": quality_json(&r.input),
                "mesh": quality_json(&r.output),
                "refined": r.refine.iter().map(|x| x.interior + x.boundary).sum::<usize>(),
                "odt_moves": r.odt_moves,
                "improvement": {"passes": r.improve.passes, "smoothed": r.improve.smoothed, "topological": r.improve.topological, "inserted": r.improve.inserted},
                "straightened": r.straightened,
                "stopped": r.stopped,
            },
            "seeds": r.seeds.iter().map(|(name, count, edges)| json!({"edge": name, "count": count, "mesh_edges": edges})).collect::<Vec<_>>(),
        }));
    }
    let sum = |f: fn(&BodyTets) -> f64| meshed.iter().map(f).sum::<f64>();
    let count = |f: fn(&BodyTets) -> usize| meshed.iter().map(f).sum::<usize>();
    let first = bodies.first();
    let reply = json!({
        "kind": meshed.first().map(|b| b.mesh.kind.name()),
        "nodes": count(|b| b.mesh.points.len()),
        "elements": count(|b| b.mesh.elements.len()),
        "first_node": first.map(|b| b["first_node"].clone()),
        "first_element": first.map(|b| b["first_element"].clone()),
        "boundary_faces": count(|b| b.mesh.boundary.len()),
        "volume": sum(|b| b.mesh.report.volume),
        "boundary_volume": sum(|b| b.mesh.report.boundary_volume),
        "surface_volume": sum(|b| b.surface_volume),
        "bodies": bodies,
    });
    Ok(cleared.report(s, reply))
}

/// Add one body's elements after the model's largest labels, tagging them with the body and each
/// boundary element face with its geometry face; returns the first node and element labels.
pub(super) fn add_tets(model: &mut fieldcraft_model::Model, body: &str, face_names: &[String], mesh: &QualityMesh) -> Result<(u32, u32)> {
    let (first_node, first_element) = super::fresh_labels(model, mesh.points.len(), mesh.elements.len())?;
    for (i, xyz) in mesh.points.iter().enumerate() {
        model.set_node(first_node + i as u32, *xyz)?;
    }
    let body_tag = model.geometry_tag(&GeometryRef::body(body)?)?;
    for (i, e) in mesh.elements.iter().enumerate() {
        let label = first_element + i as u32;
        model.set_element(label, mesh.kind, e.iter().map(|&v| first_node + v).collect())?;
        model.tag_element(label, body_tag)?;
    }
    let mut face_tags: Vec<Option<u32>> = vec![None; face_names.len()];
    for b in &mesh.boundary {
        let face = b.face as usize;
        let tag = match face_tags.get(face).copied().flatten() {
            Some(t) => t,
            None => {
                let name = face_names.get(face).ok_or_else(|| EngineError::Params(format!("boundary face {} has no geometry face", b.face)))?;
                let t = model.geometry_tag(&GeometryRef::face(body, name.clone())?)?;
                if let Some(slot) = face_tags.get_mut(face) {
                    *slot = Some(t);
                }
                t
            }
        };
        let number = b.element_face().ok_or_else(|| EngineError::Params(format!("boundary face {} has no element face number", b.opposite)))?;
        model.tag_element_face(first_element + b.tet, number, tag)?;
    }
    Ok((first_node, first_element))
}

pub const PLANAR: CommandSpec = CommandSpec {
    id: "mesh.planar",
    label: "2D Mesh",
    params: r#"{"size": 1.0, "kind": "quad4|quad8|tri3|tri6" (default quad4), and either "polygon": [[x, y], ...] with optional "holes": [[[x, y], ...], ...] and "name": "Plane" (the body name its tags use), or "body": "Body1", "face": "F2:top", "tolerance": 0.05 (a planar face of an open design, parallel to XY); optional "max_elements": 1000000 (1..=4000000)} — meshes the planar region with 2D elements in the XY plane (a face at z = c is meshed at z = 0) for plane stress, plane strain or axisymmetric sections (x the radius, Y the axis): free quadrilaterals (or constrained Delaunay triangles) of about `size`, mid-side nodes on the boundary for quadratic kinds. Replaces the body's mesh (as mesh.tet.generate), then adds nodes and elements after the largest existing labels and tags them: the elements with the face (polygon: face `area` of body `name`), the boundary element edges with edges named after the neighbouring face of the body (polygon: `e1`, `e2`, … side by side, the outer loop then each hole, and vertices `v1`, … at the corners), so sets on those edges carry loads and boundary conditions"#,
    form: Some(super::forms::PLANAR),
    kind: CommandKind::Model,
    requires: &[],
    run: planar,
};

/// Most points of a `mesh.planar` polygon or hole, and most holes.
const PLANAR_MAX_POINTS: usize = 100_000;
const PLANAR_MAX_HOLES: usize = 10_000;
const PLANAR_MAX_ELEMENTS: u64 = 4_000_000;

/// A list of `[x, y]` points.
fn polygon(v: &Value, what: &str) -> Result<Vec<[f64; 2]>> {
    let bad = || EngineError::Params(format!("`{what}` must be a list of 3 to {PLANAR_MAX_POINTS} [x, y] points"));
    let raw = v.as_array().filter(|a| (3..=PLANAR_MAX_POINTS).contains(&a.len())).ok_or_else(bad)?;
    raw.iter()
        .map(|p| match p.as_array().map(Vec::as_slice) {
            Some([x, y]) => x.as_f64().zip(y.as_f64()).map(|(x, y)| [x, y]).filter(|c| c.iter().all(|v| v.is_finite())),
            _ => None,
        })
        .collect::<Option<Vec<_>>>()
        .ok_or_else(bad)
}

fn planar(s: &mut Session, p: &Value) -> Result<Value> {
    use fieldcraft_mesh::planar::{Curve, face_boundary, mesh_region};
    use fieldcraft_model::Entity;

    let size = p
        .get("size")
        .and_then(Value::as_f64)
        .filter(|v| *v > 0.0 && v.is_finite())
        .ok_or_else(|| super::wrong(p, "size", ProblemKind::NotPositive))?;
    let kind = match p.get("kind") {
        None => ElementKind::Quad4,
        Some(k) => k
            .as_str()
            .and_then(ElementKind::from_name)
            .filter(|k| k.dimension() == 2)
            .ok_or_else(|| EngineError::param("kind", ProblemKind::NotAChoice))?,
    };
    let opts = fieldcraft_mesh::quad::QuadOptions {
        max_elements: optional_count(p, "max_elements", 1_000_000, 1..=PLANAR_MAX_ELEMENTS)?,
        ..Default::default()
    };
    // The region's loops of curves, the body name, the face its elements fill, and each curve
    // id's edge name.
    let (loops, body, face, names): (Vec<Vec<Curve>>, String, String, BTreeMap<u32, String>) = match (p.get("polygon"), p.get("body")) {
        (Some(outer), None) => {
            let mut rings = vec![polygon(outer, "polygon")?];
            if let Some(h) = p.get("holes") {
                let list = h
                    .as_array()
                    .filter(|a| a.len() <= PLANAR_MAX_HOLES)
                    .ok_or_else(|| EngineError::Params(format!("`holes` must be a list of at most {PLANAR_MAX_HOLES} polygons")))?;
                for hole in list {
                    rings.push(polygon(hole, "holes")?);
                }
            }
            let name = match p.get("name") {
                None => "Plane".to_string(),
                Some(n) => n.as_str().map(str::to_string).ok_or_else(|| EngineError::param("name", ProblemKind::NotAName))?,
            };
            GeometryRef::body(name.clone())?;
            let mut id = 0u32;
            let mut names = BTreeMap::new();
            let loops = rings
                .into_iter()
                .map(|ring| {
                    ring.into_iter()
                        .map(|pt| {
                            id += 1;
                            names.insert(id, format!("e{id}"));
                            Curve { id: Some(id), points: vec![pt] }
                        })
                        .collect()
                })
                .collect();
            (loops, name, "area".to_string(), names)
        }
        (None, Some(_)) => {
            let face_name = p.get("face").and_then(Value::as_str).ok_or_else(|| EngineError::param("face", ProblemKind::NotAName))?;
            let (shape, mesh) = super::geometry::body_mesh(s, p)?;
            let face = shape.face_index(face_name).ok_or_else(|| EngineError::Params(format!("body `{}` has no face `{face_name}`", shape.name)))?;
            let b = face_boundary(&mesh.positions, &mesh.triangles, &mesh.tri_face, face, mesh.weld_distance)?;
            let names = shape.faces().iter().enumerate().map(|(i, f)| (i as u32, f.name.clone())).collect();
            (b.loops, shape.name.clone(), face_name.to_string(), names)
        }
        _ => return Err(EngineError::Params("give either `polygon` (and `holes`) or `body` with `face` and `tolerance`".into())),
    };
    let m = mesh_region(&loops, kind, size, &opts)?;

    // Replace the body's mesh (replace.rs); other bodies and the definitions stay.
    let cleared = super::replace::clear(&mut s.model, [body.as_str()]);
    let (first_node, first_element) = super::fresh_labels(&s.model, m.nodes.len(), m.elements.len())?;
    for (i, xy) in m.nodes.iter().enumerate() {
        s.model.set_node(first_node + i as u32, [xy[0], xy[1], 0.0])?;
    }
    let region = s.model.geometry_tag(&GeometryRef::face(body.clone(), face.clone())?)?;
    for (i, nodes) in m.elements.iter().enumerate() {
        let label = first_element + i as u32;
        s.model.set_element(label, kind, nodes.iter().map(|v| first_node + v).collect())?;
        s.model.tag_element(label, region)?;
    }
    let mut tags: BTreeMap<u32, u32> = BTreeMap::new();
    for (e, edge, id) in &m.edges {
        let Some(name) = names.get(id) else { continue };
        let tag = match tags.get(id) {
            Some(t) => *t,
            None => {
                let t = s.model.geometry_tag(&GeometryRef::new(body.clone(), Entity::Edge(name.clone()))?)?;
                tags.insert(*id, t);
                t
            }
        };
        s.model.tag_element_edge(first_element + *e as u32, *edge, tag)?;
    }
    let polygon_input = p.get("polygon").is_some();
    if polygon_input {
        for (node, id) in &m.corners {
            let t = s.model.geometry_tag(&GeometryRef::new(body.clone(), Entity::Vertex(format!("v{id}")))?)?;
            s.model.tag_node(first_node + node, t)?;
        }
    }
    let edges: Vec<String> = tags.keys().filter_map(|id| names.get(id)).cloned().collect();
    let reply = json!({
        "body": body,
        "face": face,
        "kind": kind.name(),
        "nodes": m.nodes.len(),
        "elements": m.elements.len(),
        "first_node": first_node,
        "first_element": first_element,
        "edges": edges,
    });
    Ok(cleared.report(s, reply))
}

/// The body `mesh.planar` fills: `body`, or the polygon's `name` (default `Plane`).
pub(super) fn planar_body(_: &Session, p: &Value) -> Result<Vec<String>> {
    let name = match (p.get("polygon"), p.get("body")) {
        (Some(_), None) => p.get("name").map_or(Some("Plane"), Value::as_str),
        (None, Some(b)) => b.as_str(),
        _ => None,
    };
    Ok(vec![name.ok_or_else(|| EngineError::Params("give either `polygon` or `body`".into()))?.to_string()])
}

fn node_add(s: &mut Session, p: &Value) -> Result<Value> {
    let l = label(p, "label")?;
    s.model.set_node(l, xyz(p, "xyz")?)?;
    Ok(json!({"label": l}))
}

fn node_remove(s: &mut Session, p: &Value) -> Result<Value> {
    s.model.remove_node(label(p, "label")?)?;
    Ok(Value::Null)
}

fn element_add(s: &mut Session, p: &Value) -> Result<Value> {
    let l = label(p, "label")?;
    let name = p.get("kind").and_then(Value::as_str).unwrap_or_default();
    let kind = ElementKind::from_name(name).ok_or_else(|| EngineError::Params(format!("unknown element kind `{name}`")))?;
    let nodes = labels(p, "nodes", ElementKind::Hex20.node_count())?;
    if nodes.len() != kind.node_count() {
        return Err(EngineError::Params(format!("{name} needs {} nodes, got {}", kind.node_count(), nodes.len())));
    }
    s.model.set_element(l, kind, nodes)?;
    Ok(json!({"label": l}))
}

fn element_remove(s: &mut Session, p: &Value) -> Result<Value> {
    s.model.remove_element(label(p, "label")?)?;
    Ok(Value::Null)
}

/// A JSON number, or a string for ±∞ (JSON has no infinity).
fn number(v: f64) -> Value {
    if v.is_finite() { json!(v) } else { json!(if v > 0.0 { "inf" } else { "-inf" }) }
}

fn export_stl(s: &mut Session, p: &Value) -> Result<Value> {
    let set: Option<String> = super::opt(p, "set")?;
    let format: Option<String> = super::opt(p, "format")?;
    let ascii = match format.as_deref() {
        None | Some("binary") => false,
        Some("ascii") => true,
        Some(other) => {
            let shown: String = other.chars().take(32).collect();
            return Err(EngineError::Params(format!("unknown STL format `{shown}`; known: binary, ascii")));
        }
    };
    let path = super::file_param(s, p, "path", fieldcraft_io::roots::Access::Write)?;
    let m = &s.model;
    if m.elements().is_empty() {
        return Err(EngineError::Params("the mesh has no elements".into()));
    }
    let faces = match &set {
        Some(name) => {
            let shown: String = name.chars().take(64).collect();
            if !m.sets().contains_key(name) {
                return Err(EngineError::Params(format!("no set `{shown}`")));
            }
            let faces = m.resolve_set(name)?.faces;
            if faces.is_empty() {
                return Err(EngineError::Params(format!("set `{shown}` has no element faces")));
            }
            faces
        }
        None => {
            let faces = m.skin_faces();
            if faces.is_empty() {
                return Err(EngineError::Params("the mesh has only line elements: no surface to write".into()));
            }
            faces
        }
    };
    let bytes = if ascii {
        let name = std::path::Path::new(&path).file_stem().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        fieldcraft_io::stl::write_ascii(m, &faces, &name)?.into_bytes()
    } else {
        fieldcraft_io::stl::write_binary(m, &faces)?
    };
    let triangles = fieldcraft_io::stl::triangles(m, &faces)?.len();
    fieldcraft_io::fs::write_output(std::path::Path::new(&path), &bytes)?;
    // The file has no unit of its own: it is in the display length (the command's frame).
    let length = s.input_units().symbol(fieldcraft_model::units::Quantity::Length).to_string();
    Ok(json!({"path": path, "triangles": triangles, "bytes": bytes.len(), "units": {"length": length}}))
}

fn quality_report(s: &mut Session, p: &Value) -> Result<Value> {
    let bins = match p.get("bins") {
        None => 10,
        Some(v) => v.as_u64().and_then(|b| usize::try_from(b).ok()).ok_or_else(|| EngineError::invalid("bins", "must be a whole number"))?,
    };
    let selection = match p.get("elements") {
        None => None,
        Some(_) => Some(labels(p, "elements", MAX_ELEMENTS)?),
    };
    let r: QualityReport = quality::quality_report(&s.model, selection.as_deref(), bins)?;
    let skipped: serde_json::Map<String, Value> = r.skipped.iter().map(|(k, n)| (k.name().to_string(), json!(n))).collect();
    let metrics: serde_json::Map<String, Value> = r
        .metrics
        .iter()
        .map(|m| {
            let (lo, hi) = m.metric.range();
            let stats = json!({
                "count": m.count,
                "min": number(m.min),
                "max": number(m.max),
                "mean": m.mean.map(number),
                "worst": m.worst,
                "histogram": {"range": [lo, hi], "counts": m.histogram},
            });
            (m.metric.name().to_string(), stats)
        })
        .collect();
    Ok(json!({
        "summary": quality_summary(&r),
        "measured": r.measured,
        "skipped": skipped,
        "valid": r.valid,
        "invalid": listed(&r.invalid),
        "undetermined": listed(&r.undetermined),
        "volume": r.volume,
        "area": r.area,
        "metrics": metrics,
    }))
}

/// The element set `mesh.quality.mark` keeps the marked elements in.
pub const MARKED_SET: &str = "quality_marked";

/// The default limits of `mesh.quality.mark`.
const MARK_LIMITS: [(quality::Metric, f64); 6] = [
    (quality::Metric::ScaledJacobian, 0.2),
    (quality::Metric::MinDihedral, 10.0),
    (quality::Metric::MaxDihedral, 165.0),
    (quality::Metric::AspectRatio, 10.0),
    (quality::Metric::MinAngle, 15.0),
    (quality::Metric::MaxAngle, 165.0),
];

/// Most marked elements `mesh.quality.mark` lists.
const MAX_MARKED_LISTED: usize = 1000;

fn quality_mark(s: &mut Session, p: &Value) -> Result<Value> {
    let mut limits: std::collections::BTreeMap<quality::Metric, f64> = MARK_LIMITS.into_iter().collect();
    match p.get("limits") {
        None | Some(Value::Null) => {}
        Some(Value::Object(given)) => {
            for (name, v) in given {
                let metric = quality::Metric::parse(name).ok_or_else(|| {
                    let known: Vec<&str> = quality::Metric::ALL.iter().map(|m| m.name()).collect();
                    EngineError::Params(format!("no quality metric `{}` (one of {})", name.chars().take(40).collect::<String>(), known.join(", ")))
                })?;
                match v {
                    Value::Null => {
                        limits.remove(&metric);
                    }
                    v => {
                        let limit = v
                            .as_f64()
                            .filter(|x| x.is_finite())
                            .ok_or_else(|| EngineError::Params(format!("the limit of `{name}` must be a number")))?;
                        limits.insert(metric, limit);
                    }
                }
            }
        }
        Some(_) => return Err(EngineError::Params("`limits` must be an object of metric names and numbers".into())),
    }
    let selection = match p.get("elements") {
        None => None,
        Some(_) => Some(labels(p, "elements", MAX_ELEMENTS)?),
    };
    let limits: Vec<(quality::Metric, f64)> = limits.into_iter().collect();
    let (measured, marked) = quality::mark_elements(&s.model, selection.as_deref(), &limits)?;
    let set = fieldcraft_model::Set { elements: marked.iter().map(|m| m.label).collect(), ..Default::default() };
    if set.elements.is_empty() {
        if s.model.set(MARKED_SET).is_some() {
            s.model.remove_set(MARKED_SET)?;
        }
    } else if s.model.set(MARKED_SET).is_some() {
        s.model.replace_set(MARKED_SET, set)?;
    } else {
        s.model.add_set(MARKED_SET, set)?;
    }
    let invalid = marked.iter().filter(|m| m.is_invalid()).count();
    let c = s.model.classification();
    let body = |label: u32| c.elements().get(&label).and_then(|t| t.region).and_then(|r| c.entity(r)).map(|g| g.body.clone());
    let elements: Vec<Value> = marked
        .iter()
        .take(MAX_MARKED_LISTED)
        .map(|m| {
            json!({
                "label": m.label,
                "kind": m.kind.name(),
                "body": body(m.label),
                "at": m.centroid,
                "valid": m.validity.map(|v| v == quality::Validity::Valid),
                "failed": m.failed.iter().map(|(metric, v, limit)| json!({"metric": metric.name(), "value": number(*v), "limit": limit})).collect::<Vec<_>>(),
            })
        })
        .collect();
    let mut summary = format!("{} of {measured} elements marked", marked.len());
    if invalid > 0 {
        summary.push_str(&format!(", {invalid} not valid"));
    }
    if let Some(worst) = marked.first().and_then(|m| m.scaled_jacobian) {
        summary.push_str(&format!("; worst scaled Jacobian {worst:.3}"));
    }
    if !marked.is_empty() {
        summary.push_str(&format!(" (set {MARKED_SET})"));
    }
    Ok(json!({
        "summary": summary,
        "set": (!marked.is_empty()).then_some(MARKED_SET),
        "measured": measured,
        "marked": marked.len(),
        "invalid": invalid,
        "limits": limits.iter().map(|(m, v)| (m.name().to_string(), json!(v))).collect::<serde_json::Map<_, _>>(),
        "elements": elements,
        "truncated": marked.len() > MAX_MARKED_LISTED,
    }))
}

/// A quality report in a line, for the status bar and the CLI: the elements measured, how many
/// are invalid or undetermined, and the scaled Jacobian's minimum and mean.
fn quality_summary(r: &QualityReport) -> String {
    if r.measured == 0 {
        return "no elements to measure".into();
    }
    let mut text = format!("{} elements", r.measured);
    if r.invalid.is_empty() && r.undetermined.is_empty() {
        text.push_str(", all valid");
    }
    if !r.invalid.is_empty() {
        text.push_str(&format!(", {} invalid", r.invalid.len()));
    }
    if !r.undetermined.is_empty() {
        text.push_str(&format!(", {} undetermined", r.undetermined.len()));
    }
    let jacobian = r.metrics.iter().find(|m| m.metric == quality::Metric::ScaledJacobian && m.count > 0);
    if let Some((min, mean)) = jacobian.and_then(|m| m.mean.map(|mean| (m.min, mean))) {
        text.push_str(&format!("; scaled Jacobian min {min:.3}, mean {mean:.3}"));
    }
    text
}

/// The surface named by `body` (tessellated at `tolerance`) or `surface` (an imported index).
fn surface_input<'a>(s: &'a Session, p: &Value) -> Result<Cow<'a, SurfaceMesh>> {
    match (p.get("body"), p.get("surface")) {
        (Some(_), None) => Ok(Cow::Owned(super::geometry::body_mesh(s, p)?.1)),
        (None, Some(i)) => {
            let n = s.surfaces.len();
            let t = i
                .as_u64()
                .and_then(|i| usize::try_from(i).ok())
                .and_then(|i| s.surfaces.get(i))
                .ok_or_else(|| EngineError::Params(format!("`surface` must be the index of an imported surface (there are {n})")))?;
            Ok(Cow::Borrowed(&t.mesh))
        }
        _ => Err(EngineError::Params("give either `body` (with `tolerance`) or `surface`".into())),
    }
}

fn defect_json(d: &Defect, display: &fieldcraft_model::DisplaySystem) -> Value {
    let message = d.message().render(&crate::units::MeshUnits { from: crate::units::geometry(), to: display });
    json!({"kind": d.kind(), "message": message, "at": d.at()})
}

fn surface_check(s: &mut Session, p: &Value) -> Result<Value> {
    let mesh = surface_input(s, p)?;
    let options = surface_options(p, &mesh)?;
    let surface = match BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options) {
        Ok(surface) => surface,
        Err(fieldcraft_mesh::Error::Surface(d)) => {
            return Ok(
                json!({"valid": false, "defects": {"count": d.total, "listed": d.listed.iter().map(|x| defect_json(x, s.units())).collect::<Vec<_>>()}}),
            );
        }
        Err(e) => return Err(e.into()),
    };
    let mut faces: Vec<u32> = surface.tri_face().to_vec();
    faces.sort_unstable();
    faces.dedup();
    let count = |k: FeatureKind| surface.feature_edges().iter().filter(|e| e.kind == k).count();
    let shells: Vec<Value> = surface
        .shells()
        .iter()
        .map(|sh| json!({"triangles": sh.triangles, "vertices": sh.vertices, "genus": sh.genus, "volume": sh.volume, "void": sh.void}))
        .collect();
    Ok(json!({
        "valid": true,
        "vertices": surface.points().len(),
        "triangles": surface.triangles().len(),
        "input": {"vertices": mesh.positions.len(), "triangles": mesh.triangles.len()},
        "faces": faces.len(),
        "volume": surface.volume(),
        "shells": shells,
        "feature_edges": {"face_boundary": count(FeatureKind::FaceBoundary), "sharp": count(FeatureKind::Sharp)},
        "feature_corners": surface.feature_corners().len(),
    }))
}

#[cfg(test)]
mod tests {
    #[test]
    fn export_stl_writes_the_skin_or_a_set() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-stl-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("skin.stl");
        let path_str = path.to_string_lossy().into_owned();
        let mut s = Session::default();
        assert!(s.execute("mesh.export.stl", &json!({"path": path_str})).is_err(), "no mesh");
        s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        let journal = s.journal.len();
        let out = s.execute("mesh.export.stl", &json!({"path": path_str})).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert!(!bytes.starts_with(b"solid"));
        let skin = fieldcraft_io::stl::triangles(&s.model, &s.model.skin_faces()).unwrap().len();
        assert!(skin >= s.model.skin_faces().len() && skin > 0);
        assert_eq!(out["triangles"], skin);
        assert_eq!(bytes.len(), 84 + 50 * skin);
        assert_eq!(out["bytes"], bytes.len());
        assert_eq!(s.journal.len(), journal, "exporting does not change the model");

        // A set with faces writes only those; one without faces, or an unknown one, is refused.
        let faces = s.model.skin_faces();
        let two: Vec<serde_json::Value> = faces.iter().take(2).map(|f| json!({"element": f.element, "face": f.face})).collect();
        let first = *s.model.elements().keys().next().unwrap();
        s.execute("model.set.create", &json!({"name": "stl_faces", "faces": two})).unwrap();
        s.execute("model.set.create", &json!({"name": "stl_elements", "elements": [first]})).unwrap();
        let out = s.execute("mesh.export.stl", &json!({"path": path_str, "set": "stl_faces"})).unwrap();
        assert_eq!(out["triangles"], fieldcraft_io::stl::triangles(&s.model, &faces[..2]).unwrap().len());
        let other = dir.join("other.stl").to_string_lossy().into_owned();
        let e = s.execute("mesh.export.stl", &json!({"path": other, "set": "nope"})).unwrap_err().to_string();
        assert!(e.contains("no set `nope`"), "{e}");
        let e = s.execute("mesh.export.stl", &json!({"path": other, "set": "stl_elements"})).unwrap_err().to_string();
        assert!(e.contains("has no element faces"), "{e}");
        let e = s.execute("mesh.export.stl", &json!({"path": other, "format": "obj"})).unwrap_err().to_string();
        assert!(e.contains("known: binary, ascii"), "{e}");
        assert!(!dir.join("other.stl").exists());

        let out = s.execute("mesh.export.stl", &json!({"path": path_str, "format": "ascii"})).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with("solid skin\n"));
        assert_eq!(text.matches("facet normal").count(), skin);
        assert_eq!(out["triangles"], skin);
        // In the display length: the 100 mm beam reaches x = 100, and 0.1 in SI.
        let reach = |text: &str| {
            text.lines()
                .filter_map(|l| l.trim().strip_prefix("vertex "))
                .filter_map(|v| v.split_whitespace().next()?.parse::<f64>().ok())
                .fold(0.0, f64::max)
        };
        assert_eq!((reach(&text), out["units"]["length"].as_str()), (100.0, Some("mm")));
        s.execute("model.units", &json!({"system": "m-kg-s-N-Pa"})).unwrap();
        let out = s.execute("mesh.export.stl", &json!({"path": path_str, "format": "ascii"})).unwrap();
        assert_eq!((reach(&std::fs::read_to_string(&path).unwrap()), out["units"]["length"].as_str()), (0.1, Some("m")));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn export_stl_refuses_a_mesh_of_lines() {
        let mut s = Session::default();
        s.execute("mesh.node.add", &json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        s.execute("mesh.node.add", &json!({"label": 2, "xyz": [1, 0, 0]})).unwrap();
        s.execute("mesh.element.add", &json!({"label": 1, "kind": "rod2", "nodes": [1, 2]})).unwrap();
        let path = std::env::temp_dir().join(format!("fieldcraft-stl-lines-{}.stl", std::process::id()));
        let e = s.execute("mesh.export.stl", &json!({"path": path.to_string_lossy()})).unwrap_err().to_string();
        assert!(e.contains("only line elements"), "{e}");
        assert!(!path.exists());
    }

    use serde_json::json;

    use crate::Session;

    fn hex(s: &mut Session, label: u32, first_node: u32, mirror: bool) {
        let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        for (i, c) in corners.iter().enumerate() {
            let x = if mirror { -c[0] - 2 } else { c[0] };
            s.execute("mesh.node.add", &json!({"label": first_node + i as u32, "xyz": [x, c[1], c[2]]})).unwrap();
        }
        let nodes: Vec<u32> = (first_node..first_node + 8).collect();
        s.execute("mesh.element.add", &json!({"label": label, "kind": "hex8", "nodes": nodes})).unwrap();
    }

    #[test]
    fn delaunay_of_points_reports_tets_and_duplicates() {
        let mut s = Session::default();
        let cube = json!({"points": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1], [1, 1, 1]]});
        let r = s.execute("mesh.delaunay.points", &cube).unwrap();
        assert_eq!(r["dimension"], 3);
        assert_eq!(r["vertices"], 8);
        assert_eq!(r["duplicates"], json!([[8, 7]]));
        let tets = r["tets"].as_array().unwrap();
        assert!((5..=6).contains(&tets.len()), "{tets:?}");
        assert!(s.journal.is_empty(), "a report does not change the model");

        let flat = s.execute("mesh.delaunay.points", &json!({"points": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]})).unwrap();
        assert_eq!((flat["dimension"].clone(), flat["tets"].clone()), (json!(2), json!([])));
        for bad in [
            json!({}),
            json!({"points": [[0, 0]]}),
            json!({"points": [[0, 0, "x"]]}),
            json!({"points": [[f64::MAX, 0, 0]]}),
            json!({"points": [[1e-300, 0, 0]]}),
            json!({"points": vec![[0, 0, 0]; 20_001]}),
        ] {
            assert!(s.execute("mesh.delaunay.points", &bad).is_err(), "{bad}");
        }
    }

    fn mark_session() -> Session {
        let mut s = Session::default();
        let corners = json!([[0, 0, 0], [2, 0, 0], [2, 1, 0], [0, 1, 0], [0, 0, 1], [2, 0, 1], [2, 1, 1], [0, 1, 1]]);
        s.execute("mesh.block", &json!({"body": "Bar", "kind": "hex8", "corners": corners, "divisions": [2, 1, 1]})).unwrap();
        s
    }

    /// Two good hexes mark nothing; a corner pulled in marks its element, located, as the set
    /// the viewport draws; looser limits empty it; undo brings it back.
    #[test]
    fn marks_bad_elements_as_a_set() {
        let mut s = mark_session();
        let r = s.execute("mesh.quality.mark", &json!({})).unwrap();
        assert_eq!((r["measured"].as_u64(), r["marked"].as_u64()), (Some(2), Some(0)), "{r}");
        assert!(s.model.set(super::MARKED_SET).is_none() && r["set"].is_null());
        assert_eq!(r["limits"]["scaled_jacobian"], 0.2);
        // Pull the far top corner of the second hex nearly onto its middle.
        // Stored in metres: the corner at 2, 1, 1 mm.
        let (&corner, _) = s.model.nodes().iter().find(|(_, n)| n.xyz == [0.002, 0.001, 0.001]).unwrap();
        let (&bad, _) = s.model.elements().iter().find(|(_, e)| e.nodes.contains(&corner)).unwrap();
        s.execute("mesh.node.add", &json!({"label": corner, "xyz": [1.7, 0.7, 0.7]})).unwrap();
        let r = s.execute("mesh.quality.mark", &json!({})).unwrap();
        assert_eq!((r["marked"].as_u64(), r["invalid"].as_u64()), (Some(1), Some(0)), "{r}");
        assert_eq!(r["set"], super::MARKED_SET);
        let e = &r["elements"][0];
        assert_eq!((e["label"].as_u64(), e["kind"].as_str(), e["body"].as_str()), (Some(bad as u64), Some("hex8"), Some("Bar")), "{e}");
        assert!(e["failed"].as_array().unwrap().iter().any(|f| f["metric"] == "scaled_jacobian"), "{e}");
        let at: Vec<f64> = e["at"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        assert!(at[0] > 1.0 && at[0] < 2.0, "{at:?}");
        assert!(r["summary"].as_str().unwrap().starts_with("1 of 2 elements marked"), "{r}");
        assert_eq!(s.model.set(super::MARKED_SET).unwrap().elements, [bad]);
        assert_eq!(crate::view::marked(&s).into_iter().collect::<Vec<_>>(), [bad]);
        // The picture draws its faces in the marker colour, and only its.
        let p = crate::view::Palette::LIGHT;
        let (scene, _, _) = crate::view::scene(&s, &p);
        let red = scene.tris.iter().filter(|t| t.color == p.marked).count();
        assert_eq!(red, 5 * 2, "the bad hex's outer faces (all but the shared one), as triangles");
        // Without limits the (valid) element passes and the set goes.
        let all_off = json!({"limits": {"scaled_jacobian": null, "min_dihedral": null, "max_dihedral": null, "aspect_ratio": null, "min_angle": null, "max_angle": null}});
        let r = s.execute("mesh.quality.mark", &all_off).unwrap();
        assert_eq!(r["marked"], 0, "{r}");
        assert!(s.model.set(super::MARKED_SET).is_none());
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.model.set(super::MARKED_SET).unwrap().elements, [bad], "undo restores the marks");
        // Only the elements asked about.
        let other = s.model.elements().keys().copied().find(|l| *l != bad).unwrap();
        assert_eq!(s.execute("mesh.quality.mark", &json!({"elements": [other]})).unwrap()["marked"], 0);
    }

    #[test]
    fn mark_rejects_hostile_parameters() {
        let mut s = Session::default();
        assert!(s.execute("mesh.quality.mark", &json!({})).is_err(), "needs a mesh");
        let mut s = mark_session();
        for p in [
            json!({"limits": {"nope": 1}}),
            json!({"limits": {"scaled_jacobian": "low"}}),
            json!({"limits": {"scaled_jacobian": true}}),
            json!({"limits": [0.2]}),
            json!({"limits": {"x".repeat(10000): 1}}),
            json!({"elements": [0]}),
            json!({"elements": [99]}),
            json!({"elements": "all"}),
        ] {
            assert!(s.execute("mesh.quality.mark", &p).is_err(), "{p}");
        }
    }

    #[test]
    fn reports_quality_and_inverted_elements() {
        let mut s = Session::default();
        hex(&mut s, 1, 1, false);
        hex(&mut s, 2, 11, true);
        let journal = s.journal.len();
        let r = s.execute("mesh.quality", &json!({"bins": 2})).unwrap();
        assert_eq!(r["measured"], 2);
        assert_eq!(r["valid"], 1);
        assert_eq!(r["invalid"], json!({"count": 1, "labels": [2]}));
        assert_eq!(r["volume"], 0.0);
        let sj = &r["metrics"]["scaled_jacobian"];
        assert_eq!((sj["min"].clone(), sj["max"].clone(), sj["worst"].clone()), (json!(-1.0), json!(1.0), json!(2)));
        assert_eq!(sj["histogram"]["counts"], json!([1, 1]));
        assert_eq!(s.journal.len(), journal, "a report does not change the model");

        let one = s.execute("mesh.quality", &json!({"elements": [1, 1]})).unwrap();
        assert_eq!((one["measured"].clone(), one["volume"].clone()), (json!(1), json!(1.0)));
        assert!(s.execute("mesh.quality", &json!({"elements": [3]})).is_err());
        assert!(s.execute("mesh.quality", &json!({"elements": [0]})).is_err());
        assert!(s.execute("mesh.quality", &json!({"bins": 0})).is_err());
        assert!(s.execute("mesh.quality", &json!({"bins": 1e9})).is_err());
        assert!(s.execute("mesh.quality", &json!({"elements": "all"})).is_err());
    }

    #[test]
    fn a_tessellated_block_is_a_valid_boundary() {
        use fieldcraft_mesh::surface::{BoundarySurface, FeatureKind, SurfaceOptions};
        let mesh = fieldcraft_geometry::block([0.0; 3], [2.0, 3.0, 4.0]).unwrap().tessellate(0.01).unwrap();
        let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
        let s = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options).unwrap();
        assert!(s.points().len() < mesh.positions.len(), "the tessellation is a soup that welding joins");
        assert_eq!((s.shells().len(), s.shells()[0].genus), (1, 0));
        assert!((s.volume() - 24.0).abs() < 1e-9, "{}", s.volume());
        assert!(s.feature_edges().iter().all(|e| e.kind == FeatureKind::FaceBoundary));
        assert_eq!(s.feature_corners().len(), 8);
        let length: f64 = s
            .feature_edges()
            .iter()
            .map(|e| {
                let [a, b] = e.vertices.map(|v| s.points()[v as usize]);
                ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
            })
            .sum();
        assert!((length - 4.0 * (2.0 + 3.0 + 4.0)).abs() < 1e-9, "the twelve edges: {length}");
    }

    #[test]
    fn surface_check_on_an_open_design_and_on_imported_surfaces() {
        use std::sync::Arc;

        use fieldcraft_geometry::sample::{Plate, plate_design};
        use fieldcraft_geometry::tagged_mesh::TaggedSurface;
        use fieldcraft_geometry::{SurfaceMesh, open_design};

        let plate = Plate::default();
        let mut s = Session { design: Some(Arc::new(open_design(&plate_design(&plate).unwrap()).unwrap())), ..Session::default() };
        let r = s.execute("mesh.surface.check", &json!({"body": "Body1", "tolerance": 0.05})).unwrap();
        assert_eq!(r["valid"], true, "{r}");
        assert_eq!(r["faces"], 8);
        assert_eq!(r["shells"].as_array().unwrap().len(), 1);
        assert_eq!(r["shells"][0]["genus"], 1, "the hole is a handle");
        let exact = (plate.width * plate.depth - std::f64::consts::PI * (plate.hole_diameter / 2.0).powi(2)) * plate.thickness;
        let volume = r["volume"].as_f64().unwrap();
        assert!(volume >= exact && volume < exact * 1.01, "a polygonal hole inscribed in the circle: {volume} vs {exact}");
        assert!(r["feature_edges"]["face_boundary"].as_u64().unwrap() > 12);
        assert!(r["feature_corners"].as_u64().unwrap() >= 8);
        assert!(s.journal.is_empty(), "a check changes nothing");

        // A tagged tetrahedron, and the same with a face missing.
        let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let t = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        let tet =
            SurfaceMesh { positions: p.clone(), triangles: t.clone(), tri_face: vec![0, 1, 1, 1], weld_distance: 0.0, ..SurfaceMesh::default() };
        let open = SurfaceMesh { positions: p, triangles: t[..3].to_vec(), tri_face: vec![0, 1, 1], weld_distance: 0.0, ..SurfaceMesh::default() };
        for mesh in [tet, open] {
            s.surfaces.push(Arc::new(TaggedSurface { name: "Tet".into(), mesh, faces: 2, original: None }));
        }
        let r = s.execute("mesh.surface.check", &json!({"surface": 0, "feature_angle": 180})).unwrap();
        assert_eq!((r["valid"].clone(), r["feature_edges"].clone()), (json!(true), json!({"face_boundary": 3, "sharp": 0})));
        let r = s.execute("mesh.surface.check", &json!({"surface": 1})).unwrap();
        assert_eq!(r["valid"], false);
        assert_eq!(r["defects"]["count"], 3);
        assert_eq!(r["defects"]["listed"][0]["kind"], "open_edge");
        assert!(r["defects"]["listed"][0]["at"].is_array());

        for bad in [
            json!({"surface": 2}),
            json!({"surface": -1}),
            json!({"surface": 0, "body": "Body1"}),
            json!({}),
            json!({"surface": 0, "weld": 1.0}),
            json!({"surface": 0, "weld": "fine"}),
            json!({"surface": 0, "weld_distance": -1}),
            json!({"surface": 0, "weld_distance": 1}),
            json!({"surface": 0, "feature_angle": 0}),
            json!({"body": "Body1"}),
            json!({"body": "Nope", "tolerance": 0.05}),
        ] {
            assert!(s.execute("mesh.surface.check", &bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn tet_generate_meshes_a_body_with_tags_and_replaces_its_mesh() {
        use std::sync::Arc;

        use fieldcraft_geometry::open_design;
        use fieldcraft_geometry::sample::{Plate, plate_design};

        let plate = Plate::default();
        let mut s = Session { design: Some(Arc::new(open_design(&plate_design(&plate).unwrap()).unwrap())), ..Session::default() };
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "tolerance": 0.5})).unwrap();
        let elements = r["elements"].as_u64().unwrap() as usize;
        assert!(elements > 0, "{r}");
        assert_eq!(s.model.elements().len(), elements);
        assert!(s.model.elements().values().all(|e| e.kind == fieldcraft_model::ElementKind::Tet4));
        let (volume, boundary) = (r["volume"].as_f64().unwrap(), r["boundary_volume"].as_f64().unwrap());
        assert!((volume - boundary).abs() <= 1e-9 * boundary, "{r}");
        // Every element is on the body; boundary faces carry the 8 face names.
        let c = s.model.classification();
        assert_eq!(c.elements().len(), elements);
        let faces: std::collections::BTreeSet<u32> = c.elements().values().flat_map(|t| t.faces.values().copied()).collect();
        assert_eq!(faces.len(), 8);
        let tagged: usize = c.elements().values().map(|t| t.faces.len()).sum();
        assert_eq!(tagged as u64, r["boundary_faces"].as_u64().unwrap());
        // The certified quality report agrees.
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(q["valid"].as_u64().unwrap() as usize, elements);
        assert!(q["summary"].as_str().unwrap().starts_with(&format!("{elements} elements, all valid; scaled Jacobian min ")), "{q}");
        // A second run replaces the body's mesh (#267).
        let again = s.execute("mesh.tet.generate", &json!({"body": "Body1", "tolerance": 0.5})).unwrap();
        assert_eq!(again["replaced"][0]["elements"].as_u64().unwrap() as usize, elements, "{again}");
        assert_eq!(s.model.elements().len(), elements);
        // The same mesh again, under the same labels: nothing changed, so no undo step.
        assert_eq!(s.journal.len(), 1);

        let journal = s.journal.len();
        let model = s.model.clone();
        for bad in [
            json!({"body": "Body1", "tolerance": 0.5, "max_steiner": -1}),
            json!({"body": "Body1", "tolerance": 0.5, "max_steiner": 1e300}),
            json!({"body": "Body1", "tolerance": 0.5, "max_tets": 0}),
            json!({"body": "Body1", "tolerance": 0.5, "max_tets": 10}),
            json!({"body": "Body1", "tolerance": 0.5, "max_tets": "many"}),
            json!({"body": "Body1", "tolerance": -1}),
            json!({"body": "Nope", "tolerance": 0.5}),
            json!({"body": "Body1", "surface": 0}),
            json!({"surface": 0}),
        ] {
            assert!(s.execute("mesh.tet.generate", &bad).is_err(), "{bad}");
        }
        assert_eq!((s.journal.len(), &s.model), (journal, &model), "failures change nothing");
    }

    /// Without parameters every body of the open design is meshed, each at a tolerance from its
    /// size; an error says which body it is about.
    #[test]
    fn tet_generate_defaults_to_every_body() {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 20, "width": 10, "height": 10, "body_name": "Box"}},
            {"command": "solid.split_body", "params": {"body": "Box", "plane": {"origin": [10, 0, 0], "normal": [1, 0, 0]}}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        assert!(s.execute("mesh.tet.generate", &json!({})).unwrap_err().to_string().contains("Needs geometry"));
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let r = s.execute("mesh.tet.generate", &json!({})).unwrap();
        let bodies = r["bodies"].as_array().unwrap();
        assert_eq!(bodies.len(), 2, "{r}");
        assert_eq!(r["kind"], "tet4");
        assert!((r["volume"].as_f64().unwrap() - 2000.0).abs() < 1e-6, "{r}");
        assert_eq!(s.model.elements().len() as u64, r["elements"].as_u64().unwrap());
        assert_eq!(s.model.classification().elements().len(), s.model.elements().len(), "every element is on its body");
        // The second body's labels follow the first's.
        assert_eq!(bodies[1]["first_element"].as_u64(), Some(bodies[0]["elements"].as_u64().unwrap() + 1));
        assert_eq!(s.journal.len(), 2);
        // A body named alone takes the same default tolerance.
        let one = s.execute("mesh.tet.generate", &json!({"body": bodies[0]["body"]})).unwrap();
        assert_eq!(one["elements"], bodies[0]["elements"]);
        let e = s.execute("mesh.tet.generate", &json!({"max_tets": 1})).unwrap_err().to_string();
        assert!(e.starts_with("body `Box"), "{e}");
    }

    /// Persistent names grow with the design's split history: a tray floor split by joined ribs
    /// names every split, past the old 256-character cap. Such a design meshes, and a set on that
    /// face resolves by its full name, again after a remesh (#374).
    #[test]
    fn long_persistent_face_names_mesh_and_resolve() {
        let mut commands = vec![
            json!({"command": "solid.box", "params": {"length": 40, "width": 30, "height": 20, "corner": [0, 0, 0], "body_name": "Case"}}),
            json!({"command": "solid.shell", "params": {"faces": [[20, 15, 20]], "thickness": 2, "direction": "inside", "body": "Case"}}),
        ];
        // Ribs across the floor split it into strips, each named by every split before it.
        for x in (3..=36).step_by(3) {
            commands.push(json!({"command": "solid.box", "params": {"length": 1, "width": 26, "height": 5, "corner": [x, 2, 2], "operation": "join", "targets": ["Case"]}}));
        }
        let bytes = fieldcraft_geometry::script_design(&json!({"commands": commands}).to_string()).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let long = bodies[0]["faces"].as_array().unwrap().iter().filter_map(|f| f["name"].as_str()).max_by_key(|n| n.len()).unwrap().to_string();
        assert!(long.len() > 256, "the fixture needs a face name past the old cap: {long}");
        s.execute("mesh.tet.generate", &json!({"body": "Case"})).unwrap();
        s.execute("model.set.create", &json!({"name": "floor", "geometry": [{"body": "Case", "face": long}]})).unwrap();
        let first = s.model.resolve_set("floor").unwrap();
        assert!(!first.faces.is_empty());
        s.execute("mesh.tet.generate", &json!({"body": "Case", "tolerance": 0.5})).unwrap();
        let again = s.model.resolve_set("floor").unwrap();
        assert!(!again.faces.is_empty(), "the set follows the remesh by the full name");
    }

    /// The quality report's one-line summary counts the inverted elements.
    #[test]
    fn quality_summary_counts_invalid_elements() {
        let mut s = Session::default();
        assert!(matches!(s.execute("mesh.quality", &json!({})), Err(crate::EngineError::Unavailable(_))), "no elements to measure");
        for (label, xyz) in [(1, [0, 0, 0]), (2, [1, 0, 0]), (3, [0, 1, 0]), (4, [0, 0, 1])] {
            s.execute("mesh.node.add", &json!({"label": label, "xyz": xyz})).unwrap();
        }
        s.execute("mesh.element.add", &json!({"label": 1, "kind": "tet4", "nodes": [1, 2, 3, 4]})).unwrap();
        s.execute("mesh.element.add", &json!({"label": 2, "kind": "tet4", "nodes": [1, 3, 2, 4]})).unwrap();
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        assert!(q["summary"].as_str().unwrap().starts_with("2 elements, 1 invalid; scaled Jacobian min -"), "{q}");
    }

    /// The session's caps are the defaults of `max_steiner` and `max_tets`; a parameter overrides them.
    #[test]
    fn tet_generate_honours_the_session_caps() {
        use std::sync::Arc;

        use fieldcraft_geometry::open_design;
        use fieldcraft_geometry::sample::{Plate, plate_design};

        use crate::{EngineError, Limits};

        let design = Arc::new(open_design(&plate_design(&Plate::default()).unwrap()).unwrap());
        let mut s = Session { design: Some(design), limits: Limits { max_tets: 10, ..Limits::default() }, ..Session::default() };
        let plate = json!({"body": "Body1", "tolerance": 0.5});
        let refused = s.execute("mesh.tet.generate", &plate).unwrap_err();
        assert!(matches!(refused, EngineError::Mesh(ref e) if e.kind() == Some(fieldcraft_mesh::Kind::Limit)), "{refused}");
        assert!(s.model.elements().is_empty());
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "tolerance": 0.5, "max_tets": 1_000_000})).unwrap();
        assert!(r["elements"].as_u64().unwrap() > 0, "{r}");
    }

    /// The quality stage: a graded, improved tet4 mesh with its report, and tet10 with mid-edge
    /// nodes on the bore; the raw tetrahedralization on request.
    #[test]
    fn tet_generate_quality_report_and_quadratic_elements() {
        use std::sync::Arc;

        use fieldcraft_geometry::open_design;
        use fieldcraft_geometry::sample::{Plate, plate_design};

        let plate = Plate::default();
        let design = Arc::new(open_design(&plate_design(&plate).unwrap()).unwrap());
        let mut s = Session { design: Some(design.clone()), ..Session::default() };
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "tolerance": 0.2, "size": 6})).unwrap();
        let q = &r["bodies"][0]["quality"];
        assert_eq!(q["mesh"]["inadmissible"], json!(0), "{r}");
        assert_eq!(q["mesh"]["slivers"], json!(0), "{r}");
        assert!(q["mesh"]["min_dihedral"].as_f64().unwrap() > 15.0, "{r}");
        assert!(q["mesh"]["max_dihedral"].as_f64().unwrap() < 165.0, "{r}");
        assert!(q["cdt"]["slivers"].as_u64().unwrap() > 0, "{r}");
        // The mesh fills what its boundary encloses, which is the plate within the chordal error.
        let (volume, boundary) = (r["volume"].as_f64().unwrap(), r["boundary_volume"].as_f64().unwrap());
        assert!((volume - boundary).abs() <= 1e-9 * boundary, "{r}");
        let exact = plate.width * plate.depth * plate.thickness - std::f64::consts::PI * 0.25 * plate.hole_diameter.powi(2) * plate.thickness;
        // The bore's chordal error at 0.2 mm is at most 2/3 · 0.2 · its area (≈ 42 mm³).
        assert!((volume - exact).abs() < 42.0, "{volume} vs {exact}");

        // Quadratic: mid-edge nodes on the bore's cylinder.
        let mut s = Session { design: Some(design), ..Session::default() };
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "tolerance": 0.2, "size": 6, "order": "quadratic"})).unwrap();
        assert_eq!(r["kind"], json!("tet10"), "{r}");
        assert!(s.model.elements().values().all(|e| e.kind == fieldcraft_model::ElementKind::Tet10));
        let report = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(report["valid"].as_u64().unwrap(), r["elements"].as_u64().unwrap(), "{report}");
        // Nodes on the bore (within the tessellation's reach of its axis) lie on its radius.
        let (cx, cy) = (plate.depth / 2.0, plate.depth / 2.0);
        let radius = plate.hole_diameter / 2.0;
        let mut on_bore = 0;
        for n in s.model.nodes().values() {
            // In millimetres, as the design is.
            let xyz = crate::units::to_geometry_point(n.xyz);
            let d = ((xyz[0] - cx).powi(2) + (xyz[1] - cy).powi(2)).sqrt();
            if (d - radius).abs() < 0.05 {
                on_bore += 1;
                assert!((d - radius).abs() < 1e-9, "a bore node at {xyz:?} is {d} from the axis");
            }
        }
        assert!(on_bore > 20, "{on_bore}");
    }

    #[test]
    fn tet_generate_on_an_imported_surface_and_an_open_one() {
        use std::sync::Arc;

        use fieldcraft_geometry::SurfaceMesh;
        use fieldcraft_geometry::tagged_mesh::TaggedSurface;

        let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let t = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        let tet =
            SurfaceMesh { positions: p.clone(), triangles: t.clone(), tri_face: vec![0, 1, 1, 1], weld_distance: 0.0, ..SurfaceMesh::default() };
        let open = SurfaceMesh { positions: p, triangles: t[..3].to_vec(), tri_face: vec![0, 1, 1], weld_distance: 0.0, ..SurfaceMesh::default() };
        let mut s = Session::default();
        for mesh in [tet, open] {
            s.surfaces.push(Arc::new(TaggedSurface { name: "Tet".into(), mesh, faces: 2, original: None }));
        }
        let r = s.execute("mesh.tet.generate", &json!({"surface": 0, "quality": false})).unwrap();
        assert_eq!((r["elements"].clone(), r["nodes"].clone(), r["boundary_faces"].clone()), (json!(1), json!(4), json!(4)));
        let element = s.model.element(1).unwrap();
        assert_eq!(element.nodes.len(), 4);
        // Refined to a size, the flat-faced tetrahedron keeps its volume exactly.
        let r = s.execute("mesh.tet.generate", &json!({"surface": 0, "size": 0.3})).unwrap();
        assert!(r["elements"].as_u64().unwrap() > 1, "{r}");
        assert!((r["volume"].as_f64().unwrap() - 1.0 / 6.0).abs() < 1e-12, "{r}");
        let total = s.model.elements().len();
        assert!(s.execute("mesh.tet.generate", &json!({"surface": 1})).is_err(), "an open surface bounds no volume");
        assert_eq!(s.model.elements().len(), total);
    }

    #[test]
    fn planar_meshes_a_polygon_that_solves_as_a_2d_model() {
        let mut s = Session::default();
        let hole: Vec<[f64; 2]> = (0..12).map(|k| std::f64::consts::TAU * k as f64 / 12.0).map(|t| [5.0 + t.cos(), 2.0 + t.sin()]).collect();
        let r = s
            .execute(
                "mesh.planar",
                &json!({"polygon": [[0, 0], [10, 0], [10, 4], [0, 4]], "holes": [hole], "size": 0.5, "kind": "quad8", "name": "Strip"}),
            )
            .unwrap();
        assert_eq!(r["kind"], "quad8");
        let elements = r["elements"].as_u64().unwrap() as usize;
        assert!(elements > 20, "{r}");
        assert!(s.model.elements().values().all(|e| e.kind == fieldcraft_model::ElementKind::Quad8));
        assert!(s.model.nodes().values().all(|n| n.xyz[2] == 0.0));
        // Sets on the tagged edges carry the boundary condition and the load.
        for (cmd, params) in [
            ("model.set.create", json!({"name": "all", "geometry": [{"body": "Strip", "face": "area"}]})),
            ("model.set.create", json!({"name": "left", "geometry": [{"body": "Strip", "edge": "e4"}]})),
            ("model.set.create", json!({"name": "right", "geometry": [{"body": "Strip", "edge": "e2"}]})),
            ("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}})),
            ("model.section.create", json!({"name": "plate", "set": "all", "material": "Steel", "kind": {"type": "plane_stress", "thickness": 2.0}})),
            ("model.boundary_condition.create", json!({"name": "hold", "type": "fixed", "target": {"set": "left"}})),
            ("model.load.create", json!({"name": "pull", "type": "traction", "set": "right", "vector": [10.0, 0.0, 0.0]})),
            ("model.step.create", json!({"name": "Pull", "procedure": {"type": "static"}, "boundary_conditions": ["hold"], "loads": ["pull"]})),
        ] {
            s.execute(cmd, &params).unwrap_or_else(|e| panic!("{cmd}: {e}"));
        }
        let r = s.execute("solve.static", &json!({})).unwrap();
        let total: f64 = r["reaction_total"][0].as_f64().unwrap();
        // The traction over the right edge: 10 × 4 × 2.
        assert!((total + 80.0).abs() < 1e-9, "{r}");
        // A second mesh appends; failures change nothing.
        let journal = s.journal.len();
        let model = s.model.clone();
        for bad in [
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 0}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 1e-9}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 0.1, "kind": "hex8"}),
            json!({"polygon": [[0, 0], [1, 0]], "size": 0.1}),
            json!({"polygon": [[0, 0], [1, 0], [1, 1], [0, 1]], "holes": [[[2, 2], [3, 2], [3, 3]]], "size": 0.1}),
            json!({"polygon": [[0, 0], [1, 1], [1, 0], [0, 1]], "size": 0.1}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 0.1, "max_elements": 0}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 0.1, "name": ""}),
            json!({"polygon": [[0, 0], [1, 0], [0, 1]], "size": 0.1, "body": "Body1", "face": "x"}),
            json!({"body": "Body1", "face": "x", "tolerance": 0.1, "size": 1}),
            json!({"size": 1}),
        ] {
            assert!(s.execute("mesh.planar", &bad).is_err(), "{bad}");
        }
        assert_eq!((s.journal.len(), &s.model), (journal, &model), "failures change nothing");
    }

    #[test]
    fn planar_meshes_a_flat_face_of_a_design_with_edges_named_by_their_neighbours() {
        use std::sync::Arc;

        use fieldcraft_geometry::open_design;
        use fieldcraft_geometry::sample::{Plate, plate_design};

        let design = open_design(&plate_design(&Plate::default()).unwrap()).unwrap();
        let names: Vec<String> = design.bodies[0].faces().iter().map(|f| f.name.clone()).collect();
        let mut s = Session { design: Some(Arc::new(design)), ..Session::default() };
        let mut meshed = 0;
        for face in &names {
            let Ok(r) = s.execute("mesh.planar", &json!({"body": "Body1", "face": face, "tolerance": 0.2, "size": 4.0, "kind": "tri6"})) else {
                continue;
            };
            meshed += 1;
            let edges: Vec<&str> = r["edges"].as_array().unwrap().iter().map(|e| e.as_str().unwrap()).collect();
            // The plate's top and bottom faces border its four sides and the hole.
            assert!(edges.len() >= 5 && edges.iter().all(|e| names.iter().any(|n| n == e)), "{face}: {edges:?}");
        }
        assert_eq!(meshed, 2, "the top and bottom faces are flat and parallel to XY; the sides and hole are not");
    }
}
