//! The geometry boundary.
//!
//! Meshing in FieldCraft starts from a triangulated boundary surface. This crate produces it:
//! B-rep bodies come from SolidCraft's kernel (`solidcraft-kernel`, a pinned git dependency;
//! plan/adr/0001-libraries.md) and leave it as a [`SurfaceMesh`] from SolidCraft's tessellator,
//! with each triangle tagged by the B-rep face it lies on (for feature edges and for loads and
//! boundary conditions on faces). No SolidCraft type appears in the public API, so the kernel can
//! be replaced without touching the rest of the workspace. Only this crate may depend on
//! `solidcraft-*` (`cargo xtask layers` enforces it).
//!
//! Geometry enters two ways: shapes built in code (such as [`block`]) and SolidCraft designs
//! opened with [`open_design`], which rebuilds the design's parametric timeline and gives its
//! bodies with persistent face names (a boundary condition on a named face follows that face
//! through an edit of the design).
//!
//! Lengths are millimetres.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod check;
pub mod defeature;
mod design;
mod edges;
pub mod face_split;
pub mod features;
pub mod partition;
pub mod sample;
#[cfg(feature = "scripts")]
mod script;
pub mod stitch;
pub mod virtual_faces;

pub use design::{DESIGN_FORMAT, Design, FeatureError, MAX_DESIGN_BYTES, open_design};
pub use edges::{EDGE_TOLERANCE, Edge, Vertex};
pub use features::{BodyFeatures, Curve, FeatureNode, ProfileGeometry, SketchGeometry};
pub use partition::{Cut, CutReport};
pub use virtual_faces::Merge;
pub mod tagged_mesh;
#[cfg(feature = "scripts")]
pub use script::{MAX_SCRIPT_BYTES, script_design};
use solidcraft_kernel::{Body, KernelError, guard};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    #[error("invalid input: {0}")]
    Invalid(String),
    #[error("the geometry operation failed: {0}")]
    Failed(String),
    /// A design saved by a SolidCraft newer than the one FieldCraft is built with
    /// ([`open_design`]): FieldCraft needs updating, not SolidCraft.
    #[error("this design was saved by a newer SolidCraft (format {format}); this FieldCraft reads up to {reads}. Update FieldCraft to open it")]
    NewerDesign { format: u32, reads: u32 },
    /// From the model layer: a cancelled job ([`fieldcraft_model::progress::check`]).
    #[error(transparent)]
    Model(#[from] fieldcraft_model::Error),
}

pub type Result<T> = std::result::Result<T, Error>;

impl From<KernelError> for Error {
    fn from(e: KernelError) -> Self {
        match e {
            KernelError::Invalid(m) => Error::Invalid(m),
            other => Error::Failed(other.to_string()),
        }
    }
}

/// SolidCraft's linear modelling tolerance (mm): points closer than this are the same point.
const LINEAR_TOLERANCE: f64 = solidcraft_geom::LINEAR_TOL;

/// A solid body.
#[derive(Debug, Clone)]
pub struct Shape {
    pub name: String,
    body: Body,
    /// By face index (the [`SurfaceMesh::tri_face`] tags); empty for shapes built in code. With
    /// merged faces ([`virtual_faces`]) these are the virtual faces, not the B-rep's.
    faces: Vec<Face>,
    /// The virtual face (index into `faces`) of each B-rep face; empty: they are the same.
    face_map: Vec<u32>,
    /// Faces split along a plane ([`face_split`]), applied after the merges, in order.
    cuts: Vec<face_split::Cut>,
}

/// A B-rep face's persistent identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Face {
    /// Stable across rebuilds of the design: named after the feature and sketch curve that made
    /// the face (`F2:side:l4`, `F2:top`, …), not after its position in the body.
    pub name: String,
    /// The name follows the history SolidCraft recorded; `false` when it was worked out from the
    /// geometry, which an upstream edit can fool.
    pub exact: bool,
    /// A merged face's B-rep faces, by persistent name ([`virtual_faces`]); empty for a B-rep face.
    pub members: Vec<String>,
}

/// The analytic surface a B-rep face lies on, when it is one of the elementary kinds. Directions
/// are unit vectors; lengths are millimetres and angles radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FaceSurface {
    /// The outward unit normal and a point on the plane.
    Plane {
        normal: [f64; 3],
        point: [f64; 3],
    },
    Cylinder {
        origin: [f64; 3],
        axis: [f64; 3],
        radius: f64,
    },
    /// The apex, the axis pointing into the cone and the half-angle.
    Cone {
        apex: [f64; 3],
        axis: [f64; 3],
        half_angle: f64,
    },
    Sphere {
        center: [f64; 3],
        radius: f64,
    },
    /// The centre of the ring, its axis, and the ring (major) and tube (minor) radii.
    Torus {
        center: [f64; 3],
        axis: [f64; 3],
        major: f64,
        minor: f64,
    },
}

/// Mass properties and extent of a [`Shape`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Measure {
    pub volume: f64,
    pub area: f64,
    pub centroid: [f64; 3],
    pub bbox_min: [f64; 3],
    pub bbox_max: [f64; 3],
    pub faces: usize,
}

/// A triangulated boundary surface: the input to volume meshing.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SurfaceMesh {
    pub positions: Vec<[f64; 3]>,
    /// Counter-clockwise seen from outside.
    pub triangles: Vec<[u32; 3]>,
    /// The B-rep face each triangle lies on.
    pub tri_face: Vec<u32>,
    /// Copies of one point closer than this are the same point (0: only identical coordinates).
    /// A [`Shape::tessellate`] result is a triangle soup whose faces each carry their own copy of
    /// the boundary points they share, accurate to the kernel's linear tolerance, so two copies
    /// can differ by up to twice it; boundary surface preparation welds at this distance.
    pub weld_distance: f64,
    /// The B-rep edges as polylines, by edge index (empty for a surface without B-rep edges).
    pub edges: Vec<Vec<[f64; 3]>>,
    /// Per edge: the B-rep faces it bounds (one or two), by face index.
    pub edge_faces: Vec<Vec<u32>>,
}

fn xyz(v: solidcraft_geom::Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// A box between two corners.
pub fn block(a: [f64; 3], b: [f64; 3]) -> Result<Shape> {
    let v = |p: [f64; 3]| solidcraft_geom::Vec3 { x: p[0], y: p[1], z: p[2] };
    let body = guard("block", || solidcraft_kernel::box_solid(v(a), v(b)))?;
    Ok(Shape { name: "Block".into(), body, faces: Vec::new(), face_map: Vec::new(), cuts: Vec::new() })
}

impl Shape {
    /// The persistent identities of the faces, by B-rep face index (empty for shapes built in code).
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    /// The index of the face with this persistent name: the [`SurfaceMesh::tri_face`] tag its
    /// triangles carry.
    pub fn face_index(&self, name: &str) -> Option<u32> {
        self.faces.iter().position(|f| f.name == name).and_then(|i| u32::try_from(i).ok())
    }

    /// The area centroid of the face with this persistent name, from a tessellation fine enough
    /// to place it (on a curved face the centroid lies off the surface, inside its curvature).
    pub fn face_centre(&self, name: &str) -> Result<[f64; 3]> {
        let face = self.face_index(name).ok_or_else(|| Error::Invalid(format!("body `{}` has no face `{name}`", self.name)))?;
        let mesh = self.tessellate((self.size()? * edges::EDGE_TOLERANCE).max(1e-4))?;
        let (mut sum, mut area) = ([0.0; 3], 0.0);
        for (t, _) in mesh.triangles.iter().zip(&mesh.tri_face).filter(|(_, f)| **f == face) {
            let corner = |i: u32| mesh.positions.get(i as usize).copied();
            let (Some(a), Some(b), Some(c)) = (corner(t[0]), corner(t[1]), corner(t[2])) else { continue };
            let (u, v) = ([0, 1, 2].map(|k| b[k] - a[k]), [0, 1, 2].map(|k| c[k] - a[k]));
            let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            let w = 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            for (k, s) in sum.iter_mut().enumerate() {
                *s += w * (a[k] + b[k] + c[k]) / 3.0;
            }
            area += w;
        }
        if area > 0.0 && area.is_finite() {
            Ok(sum.map(|s| s / area))
        } else {
            Err(Error::Failed(format!("face `{name}` of body `{}` has no area", self.name)))
        }
    }

    /// The analytic surface of each B-rep face, by face index (the [`SurfaceMesh::tri_face`]
    /// tags): `None` for a free-form face. Meshers project nodes onto these surfaces. A merged
    /// face has the surface its B-rep faces share, `None` when they lie on different ones.
    pub fn face_surfaces(&self) -> Result<Vec<Option<FaceSurface>>> {
        let brep = self.brep_face_surfaces()?;
        self.virtual_surfaces(brep)
    }

    /// [`Shape::face_surfaces`] by B-rep face index (the pieces of [`Shape::tessellate_pieces`]).
    pub fn brep_face_surfaces(&self) -> Result<Vec<Option<FaceSurface>>> {
        guard("face surfaces", || {
            Ok((0..self.body.face_count())
                .map(|f| {
                    solidcraft_kernel::face_surface(&self.body, f).map(|s| match s {
                        solidcraft_kernel::FaceSurface::Plane { normal, point } => FaceSurface::Plane { normal: xyz(normal), point: xyz(point) },
                        solidcraft_kernel::FaceSurface::Cylinder { origin, axis, radius } => {
                            FaceSurface::Cylinder { origin: xyz(origin), axis: xyz(axis), radius }
                        }
                        solidcraft_kernel::FaceSurface::Cone { apex, axis, half_angle } => {
                            FaceSurface::Cone { apex: xyz(apex), axis: xyz(axis), half_angle }
                        }
                        solidcraft_kernel::FaceSurface::Sphere { center, radius } => FaceSurface::Sphere { center: xyz(center), radius },
                        solidcraft_kernel::FaceSurface::Torus { center, axis, major, minor } => {
                            FaceSurface::Torus { center: xyz(center), axis: xyz(axis), major, minor }
                        }
                    })
                })
                .collect())
        })
        .map_err(Error::from)
    }

    /// A characteristic size: the diagonal of the bounding box of the body's B-rep vertices
    /// (cached by the kernel; curved faces may bulge past it). For default tolerances and sizes.
    pub fn size(&self) -> Result<f64> {
        Ok(guard("size", || Ok(self.body.size()))?)
    }

    /// The boundary shells of the body: one for a solid in one piece without voids.
    pub fn shells(&self) -> Result<usize> {
        Ok(guard("shells", || Ok(self.body.shell_count()))?)
    }

    pub fn measure(&self) -> Result<Measure> {
        let m = guard("measure", || solidcraft_kernel::measure(&self.body))?;
        Ok(Measure {
            volume: m.volume,
            area: m.area,
            centroid: xyz(m.centroid),
            bbox_min: xyz(m.bbox.min),
            bbox_max: xyz(m.bbox.max),
            faces: if self.face_map.is_empty() && self.cuts.is_empty() { m.faces } else { self.faces.len() },
        })
    }

    /// Triangulate the boundary within `tolerance` (chordal deviation, mm). Triangles are tagged
    /// by face index ([`Shape::faces`]): merged faces are one face, without the edges between
    /// their B-rep faces.
    pub fn tessellate(&self, tolerance: f64) -> Result<SurfaceMesh> {
        Ok(self.virtual_mesh(self.tessellate_brep(tolerance)?))
    }

    /// [`Shape::tessellate`], and per triangle the B-rep face it lies on: its piece of a merged
    /// face, whose surface is in [`Shape::brep_face_surfaces`] (the same as its tag when no face
    /// is merged).
    pub fn tessellate_pieces(&self, tolerance: f64) -> Result<(SurfaceMesh, Vec<u32>)> {
        let brep = self.tessellate_brep(tolerance)?;
        let pieces = brep.tri_face.clone();
        Ok((self.virtual_mesh(brep), pieces))
    }

    /// The B-rep's own tessellation, tagged by B-rep face.
    fn tessellate_brep(&self, tolerance: f64) -> Result<SurfaceMesh> {
        if !(tolerance > 0.0 && tolerance.is_finite()) {
            return Err(Error::Invalid(format!("tessellation tolerance must be a positive number, got {tolerance}")));
        }
        let m = guard("tessellate", || self.body.tessellate(tolerance))?;
        Ok(SurfaceMesh {
            positions: m.positions.into_iter().map(xyz).collect(),
            triangles: m.triangles,
            tri_face: m.tri_face,
            weld_distance: 2.0 * LINEAR_TOLERANCE,
            edges: m.edges.into_iter().map(|e| e.into_iter().map(xyz).collect()).collect(),
            edge_faces: m.edge_faces,
        })
    }

    /// Whether some face is made of several B-rep faces ([`virtual_faces`]).
    pub fn has_merged_faces(&self) -> bool {
        !self.face_map.is_empty()
    }
}

impl SurfaceMesh {
    /// The volume the triangles enclose (divergence theorem: the sum of the signed tetrahedra each
    /// triangle spans with the origin). Equals the solid's volume for a closed, outward-oriented
    /// surface; a triangle with a vertex index out of range counts as nothing.
    pub fn enclosed_volume(&self) -> f64 {
        let six_v: f64 = self
            .triangles
            .iter()
            .filter_map(|t| {
                let [a, b, c] = [self.positions.get(t[0] as usize)?, self.positions.get(t[1] as usize)?, self.positions.get(t[2] as usize)?];
                Some(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))
            })
            .sum();
        six_v / 6.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The plate's side face at x = 0 (40 × 30 × 10, the hole well inside) has its centre in its
    /// middle; an unknown face is an error.
    #[test]
    fn a_face_centre_is_its_area_centroid() {
        let d = open_design(&sample::plate_design(&sample::Plate::default()).unwrap()).unwrap();
        let body = &d.bodies[0];
        let surfaces = body.face_surfaces().unwrap();
        let side = surfaces
            .iter()
            .position(|s| matches!(s, Some(FaceSurface::Plane { normal, point }) if normal[0] < -0.5 && point[0].abs() < 1e-9))
            .unwrap();
        let c = body.face_centre(&body.faces()[side].name).unwrap();
        assert!((c[0]).abs() < 1e-9 && (c[1] - 15.0).abs() < 1e-6 && (c[2] - 5.0).abs() < 1e-6, "{c:?}");
        assert!(body.face_centre("nope").is_err());
    }

    #[test]
    fn block_measures() {
        let m = block([0., 0., 0.], [10., 20., 30.]).unwrap().measure().unwrap();
        assert!((m.volume - 6000.0).abs() < 1e-6, "{}", m.volume);
        assert!((m.area - 2200.0).abs() < 1e-6, "{}", m.area);
        assert_eq!(m.faces, 6);
        assert_eq!(m.bbox_max, [10., 20., 30.]);
    }

    #[test]
    fn block_faces_are_planes_through_the_box() {
        let surfaces = block([0., 0., 0.], [10., 20., 30.]).unwrap().face_surfaces().unwrap();
        assert_eq!(surfaces.len(), 6);
        for s in surfaces {
            let Some(FaceSurface::Plane { normal, point }) = s else { panic!("{s:?}") };
            // Each plane holds one of the box's faces: the point's coordinate along the normal is
            // 0 or the box's extent there.
            let along: f64 = normal.iter().zip(point).map(|(n, p)| n * p).sum();
            assert!([0.0, 10.0, 20.0, 30.0].iter().any(|v| (along.abs() - v).abs() < 1e-9), "{normal:?} {point:?}");
        }
    }

    #[test]
    fn tessellation_covers_every_face_and_encloses_the_volume() {
        let mesh = block([0., 0., 0.], [2., 3., 4.]).unwrap().tessellate(0.01).unwrap();
        assert_eq!(mesh.triangles.len(), mesh.tri_face.len());
        let mut faces: Vec<u32> = mesh.tri_face.clone();
        faces.sort_unstable();
        faces.dedup();
        assert_eq!(faces.len(), 6);
        // The outward-oriented triangles enclose the block's volume.
        assert!((mesh.enclosed_volume() - 24.0).abs() < 1e-9, "{}", mesh.enclosed_volume());
    }

    #[test]
    fn hostile_input_is_an_error() {
        assert!(block([0.; 3], [0.; 3]).is_err());
        assert!(block([0.; 3], [f64::NAN, 1., 1.]).is_err());
        let b = block([0.; 3], [1.; 3]).unwrap();
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(b.tessellate(t).is_err());
        }
    }
}
