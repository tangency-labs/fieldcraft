//! Curved quadratic MITC8 and MITC6-a: small-displacement layered elasticity and inertia.
//! MITC6-a uses Lee–Bathe (2004), (29)–(38), Figs. 6–7a, on exactly six nodes.
//! Triangles: corners (0,0),(1,0),(0,1), then edges 0–1,1–2,2–0.
//! No centroid DOF; the centroid is only a shear tying location.
//!
//! K. J. Bathe & E. N. Dvorkin, IJNME 22 (1986), 697–722,
//! doi:10.1002/nme.1620220312, Fig. 7 and equations (9)–(22): the in-layer
//! *tensors*, not just their scalar components, are interpolated. Geometry and
//! director kinematics also follow Lee & Bathe, Computers & Structures 82 (2004),
//! 945–962, equations (1)–(8), doi:10.1016/j.compstruc.2004.02.004.
//!
//! Four corners counterclockwise from (-1,-1), then four edge midpoints. No
//! centre node or condensed centre unknown. Global nodal DOFs are ux,uy,uz,rx,ry,rz.
//! The physical director variation is theta cross director (two independent
//! rotations). A separate Hughes–Brezzi drilling penalty, as in `shell`, supplies
//! the sixth rotation; this addition is not part of the original five-DOF MITC8.
//! Standalone nodal directors are normals of the quadratic reference surface;
//! the solver supplies sign-aligned averages at shared nodes for smooth patches.
//! They are interpolated WITHOUT renormalization through the element, preserving
//! exact infinitesimal rigid motions. Thickness is constant, with the same signed
//! reference-surface offset as `ShellSection`. No finite rotations. Mass uses the same director kinematics.

use crate::basis::Basis;
use crate::dense::{Matrix, inverse3};
use crate::frame::{Axes, Vec3, cross, dot, norm, sub, tensor_to_global, unit};
use crate::laminate::{Lamina, Laminate, PlyFields, PlyProperties, section_abscissae};
use crate::quadrature::{Rule, gauss_legendre};
use crate::shell::{DRILLING, ShellSection, ShellStress};
use crate::{Error, Result};
use fieldcraft_model::MaterialOrientation;

pub type Mitc8 = QuadraticShell<8, 48>;
/// Six-node Lee–Bathe MITC6-a (2004), not MITC6-b or a seven-node element.
pub type Mitc6 = QuadraticShell<6, 36>;
type Tensor = [[f64; 3]; 3];
type TensorRows<const D: usize> = [Tensor; D];
/// Engineering strain rows: 11,22,12,13,23 in the integration point's orthonormal frame.
pub type StrainOperator<const D: usize = 48> = [[f64; D]; 5];

#[derive(Debug, Clone)]
pub struct QuadraticShell<const N: usize, const DOFS: usize> {
    xyz: [Vec3; N],
    directors: [Vec3; N],
    section: ShellSection,
    layers: Option<Layers>,
}

#[derive(Debug, Clone)]
struct Layers {
    plies: Vec<PlyProperties>,
    orientation: Option<MaterialOrientation>,
    reference_axis: u8,
    drilling: f64,
}

struct Map {
    g: [Vec3; 3],
    dual: [Vec3; 3],
    axes: Axes,
    det: f64,
}

fn scaled(a: Vec3, b: f64) -> Vec3 {
    a.map(|v| v * b)
}
fn sum(a: Vec3, b: Vec3) -> Vec3 {
    std::array::from_fn(|i| a[i] + b[i])
}
fn outer(a: Vec3, b: Vec3) -> Tensor {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i] * b[j]))
}
fn contract(e: &Tensor, a: Vec3, b: Vec3) -> f64 {
    (0..3).map(|i| (0..3).map(|j| a[i] * e[i][j] * b[j]).sum::<f64>()).sum()
}
fn add_tensor(a: &mut Tensor, b: &Tensor, w: f64) {
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] += w * b[i][j];
        }
    }
}
fn sym_dyad(a: Vec3, b: Vec3) -> Tensor {
    let mut out = outer(a, b);
    add_tensor(&mut out, &outer(b, a), 1.0);
    out
}
fn dual(g: &[Vec3; 3]) -> Result<[Vec3; 3]> {
    let (inverse, _) = inverse3(g).ok_or_else(|| Error::Degenerate("MITC8/MITC6-a singular covariant basis".into()))?;
    Ok(std::array::from_fn(|i| std::array::from_fn(|j| inverse[j][i])))
}
const NATURAL_SLACK: f64 = 1e-12; // unit-factors: dimensionless natural-coordinate tolerance
fn check_point(p: [f64; 3]) -> Result<()> {
    if p.iter().any(|v| !v.is_finite() || v.abs() > 1.0 + NATURAL_SLACK) {
        return Err(Error::Invalid("MITC8/MITC6-a natural coordinates must be finite and in [-1,1]".into()));
    }
    Ok(())
}

impl<const N: usize, const DOFS: usize> QuadraticShell<N, DOFS> {
    fn basis() -> Basis {
        if N == 6 { Basis::Tri6 } else { Basis::Quad8 }
    }
    fn centre() -> [f64; 3] {
        if N == 6 { [1.0 / 3.0, 1.0 / 3.0, 0.0] } else { [0.0; 3] }
    }
    fn surface_points() -> Vec<crate::quadrature::Point> {
        if N == 8 {
            return Rule::Quad(3).points();
        }
        // Symmetrized positive Stroud conical product, exact through degree four.
        // Averaging all barycentric permutations preserves numbering invariance
        // also for curved maps, without adding geometric or displacement nodes.
        let mut out = Vec::new();
        for p in Rule::TriConical(3).points() {
            let l = [1.0 - p.xi[0] - p.xi[1], p.xi[0], p.xi[1]];
            for [i, j] in [[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]] {
                out.push(crate::quadrature::Point { xi: [l[i], l[j], 0.0], weight: p.weight / 6.0 });
            }
        }
        out
    }
    /// Lee–Bathe 2004 Fig.6 and (29)–(34). Components are covariant tensor
    /// strains, not engineering shear. q=(-r+s)/sqrt(2); qq=(rr+ss)/2-rs.
    fn triangle_layer(&self, p: [f64; 3], variation: Option<usize>) -> Result<TensorRows<DOFS>> {
        let root = 3.0_f64.sqrt();
        let a = 0.5 - 0.5 / root;
        let b = 1.0 - a;
        let sample = |r: f64, s: f64, j: usize| -> Result<[f64; DOFS]> {
            let q = [r, s, p[2]];
            let m = self.map(q)?;
            if j < 2 {
                return Ok(self.covariant(q, &m, j, j, variation));
            }
            let rr = self.covariant(q, &m, 0, 0, variation);
            let ss = self.covariant(q, &m, 1, 1, variation);
            let rs = self.covariant(q, &m, 0, 1, variation);
            Ok(std::array::from_fn(|k| 0.5 * (rr[k] + ss[k]) - rs[k]))
        };
        let values = [
            [sample(a, 0.0, 0)?, sample(b, 0.0, 0)?, sample(a, 1.0 / root, 0)?],
            [sample(0.0, a, 1)?, sample(0.0, b, 1)?, sample(1.0 / root, a, 1)?],
            [sample(b, a, 2)?, sample(a, b, 2)?, sample(a, a, 2)?],
        ];
        let m = self.map(p)?;
        Ok(std::array::from_fn(|k| {
            let mut v = [0.0; 3];
            for j in 0..3 {
                let mean = 0.5 * (values[j][0][k] + values[j][1][k]);
                let slope = 0.5 * root * (values[j][1][k] - values[j][0][k]);
                let (aa, bb, x, y) = if j == 2 {
                    (mean + slope, -2.0 * slope, p[0], 1.0 - p[0] - p[1])
                } else if j == 1 {
                    (mean - slope, 2.0 * slope, p[1], p[0])
                } else {
                    (mean - slope, 2.0 * slope, p[0], p[1])
                };
                let cc = root * (values[j][2][k] - aa - bb * a);
                v[j] = aa + bb * x + cc * y;
            }
            let mut e = [[0.0; 3]; 3];
            add_tensor(&mut e, &outer(m.dual[0], m.dual[0]), v[0]);
            add_tensor(&mut e, &outer(m.dual[1], m.dual[1]), v[1]);
            add_tensor(&mut e, &sym_dyad(m.dual[0], m.dual[1]), 0.5 * (v[0] + v[1]) - v[2]);
            e
        }))
    }
    /// Lee–Bathe 2004 Fig.7a, (35)–(38): six edge tying points and
    /// two shear components at the centroid, with quadratic interior terms.
    fn triangle_shear(&self, p: [f64; 3], variation: Option<usize>) -> Result<TensorRows<DOFS>> {
        let root = 3.0_f64.sqrt();
        let a = 0.5 - 0.5 / root;
        let b = 1.0 - a;
        let sample = |r: f64, s: f64, j: usize| -> Result<[f64; DOFS]> {
            let q = [r, s, p[2]];
            let m = self.map(q)?;
            Ok(self.covariant(q, &m, j, 2, variation))
        };
        let r1 = sample(a, 0.0, 0)?;
        let r2 = sample(b, 0.0, 0)?;
        let s1 = sample(0.0, a, 1)?;
        let s2 = sample(0.0, b, 1)?;
        let hr1 = sample(b, a, 0)?;
        let hr2 = sample(a, b, 0)?;
        let hs1 = sample(b, a, 1)?;
        let hs2 = sample(a, b, 1)?;
        let cr = sample(1.0 / 3.0, 1.0 / 3.0, 0)?;
        let cs = sample(1.0 / 3.0, 1.0 / 3.0, 1)?;
        let m = self.map(p)?;
        let r = p[0];
        let s = p[1];
        Ok(std::array::from_fn(|k| {
            let lr = 0.5 * root * (r2[k] - r1[k]);
            let ls = 0.5 * root * (s2[k] - s1[k]);
            let a1 = 0.5 * (r1[k] + r2[k]) - lr;
            let b1 = 2.0 * lr;
            let a2 = 0.5 * (s1[k] + s2[k]) - ls;
            let c2 = 2.0 * ls;
            let mr = 0.5 * (hr1[k] + hr2[k]);
            let ms = 0.5 * (hs1[k] + hs2[k]);
            let lr = 0.5 * root * (hr2[k] - hr1[k]);
            let ls = 0.5 * root * (hs2[k] - hs1[k]);
            let c1 = 6.0 * cr[k] - 3.0 * cs[k] + 2.0 * ms - 2.0 * mr - 4.0 * a1 - b1 + a2;
            let b2 = -3.0 * cr[k] + 6.0 * cs[k] - 2.0 * ms + 2.0 * mr + a1 - 4.0 * a2 - c2;
            let e2 = 3.0 * cr[k] - 6.0 * cs[k] + 3.0 * ms - ls - 3.0 * mr + lr + b1 + 3.0 * a2 + c2;
            let f1 = -6.0 * cr[k] + 3.0 * cs[k] - 3.0 * ms - ls + 3.0 * mr + lr + 3.0 * a1 + b1 + c2;
            let rt = a1 + b1 * r + c1 * s - e2 * r * s + f1 * s * s;
            let st = a2 + b2 * r + c2 * s - f1 * r * s + e2 * r * r;
            let mut e = [[0.0; 3]; 3];
            add_tensor(&mut e, &sym_dyad(m.dual[0], m.dual[2]), rt);
            add_tensor(&mut e, &sym_dyad(m.dual[1], m.dual[2]), st);
            e
        }))
    }
    pub fn new(xyz: &[[f64; 3]], section: ShellSection) -> Result<Self> {
        if !matches!((N, DOFS), (6, 36) | (8, 48)) {
            return Err(Error::Invalid("quadratic shell requires (6,36) or (8,48) nodes/DOFs".into()));
        }
        let xyz: [Vec3; N] = xyz.try_into().map_err(|_| Error::Invalid("quadratic shell has the wrong number of nodes".into()))?;
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        // ShellSection's fields are public; validate them even for direct kernel callers.
        ShellSection::new(
            section.thickness,
            -section.eccentricity / section.thickness,
            &fieldcraft_model::Elastic::Isotropic { young: section.young, poisson: section.poisson },
        )?;
        if !(section.young > 0.0 && section.young.is_finite() && section.poisson > -1.0 && section.poisson < 0.5) {
            return Err(Error::Invalid("MITC8/MITC6-a needs a finite positive modulus and -1 < poisson < 0.5".into()));
        }
        let mut directors = [[0.0; 3]; N];
        for (i, p) in Self::basis().nodes().iter().enumerate() {
            let sh = Self::basis().eval(*p);
            let tangent: [Vec3; 2] = std::array::from_fn(|a| std::array::from_fn(|k| (0..N).map(|j| sh.dn[j][a] * (xyz[j][k] - xyz[0][k])).sum()));
            directors[i] = unit(cross(tangent[0], tangent[1])).ok_or_else(|| Error::Degenerate("MITC8/MITC6-a zero nodal normal".into()))?;
        }
        let centre = Self::basis().eval(Self::centre());
        let tangent: [Vec3; 2] = std::array::from_fn(|a| std::array::from_fn(|k| (0..N).map(|j| centre.dn[j][a] * (xyz[j][k] - xyz[0][k])).sum()));
        let normal = unit(cross(tangent[0], tangent[1])).ok_or_else(|| Error::Degenerate("MITC8/MITC6-a singular centre surface".into()))?;
        if directors.iter().any(|n| dot(*n, normal) <= 1e-8) {
            return Err(Error::Degenerate(
                "MITC8/MITC6-a folded surface or nodal normals turning 90 degrees from the centre; refine the geometry".into(),
            ));
        }
        let out = Self { xyz, directors, section, layers: None };
        // Geometry screen, including both thickness faces and the reference nodes. All
        // quadrature and tying locations are checked again when used. This is a sampled
        // Jacobian check, not a global positivity certificate for an arbitrary curved map.
        for p in Self::basis().nodes().iter().chain(&[Self::centre()]) {
            for t in [-1.0, 0.0, 1.0] {
                out.map([p[0], p[1], t])?;
            }
        }
        for p in Self::surface_points() {
            for t in [-1.0, 0.0, 1.0] {
                out.map([p.xi[0], p.xi[1], t])?;
            }
        }
        Ok(out)
    }

    /// Layerwise plane-stress integration following Nettles, NASA RP-1351
    /// (1994), sections II-IV. Geometry/strain interpolation remains BD86 MITC8.
    pub fn laminated(
        xyz: &[[f64; 3]],
        plies: &[PlyProperties],
        offset: f64,
        orientation: Option<MaterialOrientation>,
        reference_axis: u8,
    ) -> Result<Self> {
        let laminate = Laminate::new(plies, offset)?;
        if orientation.is_none() && reference_axis != 1 {
            return Err(Error::Invalid("MITC8/MITC6-a non-default reference axis needs a material orientation".into()));
        }
        // This section supplies geometry only; the law below always uses the
        // individual laminae and their own thickness intervals.
        let geometry = ShellSection::new(laminate.thickness(), offset, &fieldcraft_model::Elastic::Isotropic { young: 1.0, poisson: 0.0 })?;
        let mut out = Self::new(xyz, geometry)?;
        out.layers = Some(Layers { plies: plies.to_vec(), orientation, reference_axis, drilling: laminate.drilling() });
        Ok(out)
    }

    /// Use mesh-resolved unit directors, signed along this element's normal.
    /// Sharing these vectors and global rotations makes the interpolated director
    /// displacement identical on a common quadratic edge. Geometry is rechecked.
    pub fn with_directors(mut self, directors: [[f64; 3]; N]) -> Result<Self> {
        for (n, old) in directors.iter().zip(&self.directors) {
            if n.iter().any(|v| !v.is_finite()) || (norm(*n) - 1.0).abs() > 1e-10 || dot(*n, *old) < 0.5 {
                return Err(Error::Invalid("MITC8/MITC6-a shared directors must be unit vectors within 60 degrees of local normals".into()));
            }
        }
        self.directors = directors;
        for p in Self::basis().nodes().iter().chain(&[Self::centre()]) {
            for t in [-1.0, 0.0, 1.0] {
                self.map([p[0], p[1], t])?;
            }
        }
        for p in Self::surface_points() {
            for t in [-1.0, 0.0, 1.0] {
                self.map([p.xi[0], p.xi[1], t])?;
            }
        }
        Ok(self)
    }

    pub fn directors(&self) -> &[Vec3; N] {
        &self.directors
    }

    /// Covariant geometry: x=sum N_i (x_i + z n_i), z=e+h*t/2.
    fn map(&self, p: [f64; 3]) -> Result<Map> {
        check_point(p)?;
        if N == 6 && (p[0] < -NATURAL_SLACK || p[1] < -NATURAL_SLACK || p[0] + p[1] > 1.0 + NATURAL_SLACK) {
            return Err(Error::Invalid("MITC6-a natural point outside reference triangle".into()));
        }
        let sh = Self::basis().eval(p);
        let z = self.section.eccentricity + 0.5 * self.section.thickness * p[2];
        let mut g = [[0.0; 3]; 3];
        for i in 0..N {
            let x = sum(sub(self.xyz[i], self.xyz[0]), scaled(self.directors[i], z));
            for a in 0..2 {
                g[a] = sum(g[a], scaled(x, sh.dn[i][a]));
            }
            g[2] = sum(g[2], scaled(self.directors[i], 0.5 * self.section.thickness * sh.n[i]));
        }
        let det = dot(cross(g[0], g[1]), g[2]);
        let scale = norm(g[0]) * norm(g[1]) * norm(g[2]);
        if !(det.is_finite() && det > 1e-10 * scale) {
            return Err(Error::Degenerate("MITC8/MITC6-a inverted, folded or singular thickness Jacobian".into()));
        }
        let n = unit(g[2]).ok_or_else(|| Error::Degenerate("MITC8/MITC6-a zero director".into()))?;
        let e1 = unit(sub(g[0], scaled(n, dot(g[0], n)))).ok_or_else(|| Error::Degenerate("MITC8/MITC6-a singular tangent frame".into()))?;
        Ok(Map { g, dual: dual(&g)?, axes: [e1, cross(n, e1), n], det })
    }

    /// Displacement derivatives for a unit global nodal DOF, using theta cross n.
    fn displacement_gradient(&self, sh: &crate::shape::Shape, z: f64, dof: usize) -> [Vec3; 3] {
        let node = dof / 6;
        let component = dof % 6;
        let mut v = [0.0; 3];
        v[component % 3] = 1.0;
        if component < 3 {
            return [scaled(v, sh.dn[node][0]), scaled(v, sh.dn[node][1]), [0.0; 3]];
        }
        let v = cross(v, self.directors[node]);
        [scaled(v, z * sh.dn[node][0]), scaled(v, z * sh.dn[node][1]), scaled(v, 0.5 * self.section.thickness * sh.n[node])]
    }

    fn covariant(&self, p: [f64; 3], m: &Map, a: usize, b: usize, variation: Option<usize>) -> [f64; DOFS] {
        let sh = Self::basis().eval(p);
        let z = self.section.eccentricity + 0.5 * self.section.thickness * p[2];
        std::array::from_fn(|k| {
            let d = self.displacement_gradient(&sh, z, k);
            if let Some(i) = variation {
                let v = self.displacement_gradient(&sh, z, i);
                let mut second = [[0.0; 3]; 3];
                if i / 6 == k / 6 && i % 6 >= 3 && k % 6 >= 3 {
                    let node = i / 6;
                    let mut ei = [0.0; 3];
                    let mut ek = [0.0; 3];
                    ei[i % 3] = 1.0;
                    ek[k % 3] = 1.0;
                    // Unit-director retraction d(theta)=(d+theta x d)/|d+theta x d|.
                    // Its mixed derivative is -d (ei x d).(ek x d).
                    // The drilling coordinate is a null direction of this retraction.
                    let dd = scaled(self.directors[node], -dot(cross(ei, self.directors[node]), cross(ek, self.directors[node])));
                    second = [scaled(dd, z * sh.dn[node][0]), scaled(dd, z * sh.dn[node][1]), scaled(dd, 0.5 * self.section.thickness * sh.n[node])];
                }
                0.5 * (dot(v[a], d[b]) + dot(v[b], d[a]) + dot(m.g[a], second[b]) + dot(m.g[b], second[a]))
            } else {
                0.5 * (dot(m.g[a], d[b]) + dot(m.g[b], d[a]))
            }
        })
    }

    /// Direct in-layer tensor, BD86 equations (10),(12); membrane and bending
    /// coexist through z. Normal strain is eliminated by the plane-stress law.
    fn in_layer(&self, p: [f64; 3], variation: Option<usize>) -> Result<TensorRows<DOFS>> {
        let m = self.map(p)?;
        let rr = self.covariant(p, &m, 0, 0, variation);
        let ss = self.covariant(p, &m, 1, 1, variation);
        let rs = self.covariant(p, &m, 0, 1, variation);
        Ok(std::array::from_fn(|k| {
            let mut e = [[0.0; 3]; 3];
            add_tensor(&mut e, &outer(m.dual[0], m.dual[0]), rr[k]);
            add_tensor(&mut e, &outer(m.dual[1], m.dual[1]), ss[k]);
            add_tensor(&mut e, &sym_dyad(m.dual[0], m.dual[1]), rs[k]);
            e
        }))
    }

    /// BD86 (11)–(18). Use the repository's corner order throughout. Edge tensors
    /// keep the transverse-to-edge normal strain direct and reconstruct the other
    /// components from adjacent corner tensors in an orthogonalized covariant basis.
    fn assumed_layer(&self, p: [f64; 3], variation: Option<usize>) -> Result<TensorRows<DOFS>> {
        if N == 6 {
            return self.triangle_layer(p, variation);
        }
        let a = 1.0 / 3.0f64.sqrt();
        let points = Self::basis().nodes();
        let mut samples = vec![[[0.0; 3]; 3]; DOFS * 8];
        for i in 0..4 {
            let q = [points[i][0] * a, points[i][1] * a, p[2]];
            let e = self.in_layer(q, variation)?;
            samples[i * DOFS..(i + 1) * DOFS].copy_from_slice(&e);
        }
        for edge in 0..4 {
            let i = edge + 4;
            let next = (edge + 1) % 4;
            let q = [points[i][0] * a, points[i][1] * a, p[2]];
            let m = self.map(q)?;
            // Horizontal edges: preserve ss; vertical edges: preserve rr.
            let keep = if edge % 2 == 0 { 1 } else { 0 };
            let other = 1 - keep;
            let mut g = m.g;
            g[other] = sub(g[other], scaled(g[keep], dot(g[other], g[keep]) / dot(g[keep], g[keep])));
            let d = dual(&g)?;
            let direct = self.covariant(q, &m, keep, keep, variation);
            for k in 0..DOFS {
                let mut mean = samples[edge * DOFS + k];
                for row in &mut mean {
                    for v in row {
                        *v *= 0.5;
                    }
                }
                add_tensor(&mut mean, &samples[next * DOFS + k], 0.5);
                let mut e = [[0.0; 3]; 3];
                add_tensor(&mut e, &outer(d[keep], d[keep]), direct[k]);
                add_tensor(&mut e, &outer(d[other], d[other]), contract(&mean, g[other], g[other]));
                add_tensor(&mut e, &sym_dyad(d[0], d[1]), contract(&mean, g[0], g[1]));
                samples[i * DOFS + k] = e;
            }
        }
        let weights = Self::basis().eval([p[0] / a, p[1] / a, 0.0]);
        let mut out = [[[0.0; 3]; 3]; DOFS];
        for i in 0..8 {
            for k in 0..DOFS {
                add_tensor(&mut out[k], &samples[i * DOFS + k], weights.n[i]);
            }
        }
        Ok(out)
    }

    /// BD86 (19)–(22): four edge samples plus the mean of RA/RB (or SA/SB)
    /// multiplied by the contravariant dyad at the centre, NOT the dyad at RA/RB.
    fn shear(&self, p: [f64; 3], assumed: bool, variation: Option<usize>) -> Result<TensorRows<DOFS>> {
        if N == 6 && assumed {
            return self.triangle_shear(p, variation);
        }
        let mut out = [[[0.0; 3]; 3]; DOFS];
        let a = 1.0 / 3.0f64.sqrt();
        for along in 0..2 {
            if !assumed {
                let m = self.map(p)?;
                let v = self.covariant(p, &m, along, 2, variation);
                for k in 0..DOFS {
                    add_tensor(&mut out[k], &sym_dyad(m.dual[along], m.dual[2]), v[k]);
                }
                continue;
            }
            let across = 1 - along;
            let bubble = (1.0 - (p[along] / a).powi(2)) * (1.0 - p[across].powi(2));
            for signs in [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]] {
                let mut q = [0.0, 0.0, p[2]];
                q[along] = a * signs[0];
                q[across] = signs[1];
                let m = self.map(q)?;
                let v = self.covariant(q, &m, along, 2, variation);
                let w = 0.25 * ((1.0 + signs[0] * p[along] / a) * (1.0 + signs[1] * p[across]) - bubble);
                for k in 0..DOFS {
                    add_tensor(&mut out[k], &sym_dyad(m.dual[along], m.dual[2]), w * v[k]);
                }
            }
            let centre = self.map([0.0, 0.0, p[2]])?;
            for sign in [-1.0, 1.0] {
                let mut q = [0.0, 0.0, p[2]];
                q[along] = sign * a;
                let m = self.map(q)?;
                let v = self.covariant(q, &m, along, 2, variation);
                for k in 0..DOFS {
                    add_tensor(&mut out[k], &sym_dyad(centre.dual[along], centre.dual[2]), 0.5 * bubble * v[k]);
                }
            }
        }
        Ok(out)
    }

    /// Physical engineering-strain operator at a natural point; `assumed=false`
    /// is the direct displacement field (for independent diagnostics/recovery).
    pub fn strain_operator(&self, p: [f64; 3], assumed: bool) -> Result<StrainOperator<DOFS>> {
        self.strain_derivative(p, assumed, None)
    }

    // The MITC projection is linear in covariant strain with reference-geometry
    // coefficients, so the same tying operator applies to both variations.
    fn strain_derivative(&self, p: [f64; 3], assumed: bool, variation: Option<usize>) -> Result<StrainOperator<DOFS>> {
        let m = self.map(p)?;
        let mut e = if assumed { self.assumed_layer(p, variation)? } else { self.in_layer(p, variation)? };
        let shear = self.shear(p, assumed, variation)?;
        let mut b = [[0.0; DOFS]; 5];
        for k in 0..DOFS {
            add_tensor(&mut e[k], &shear[k], 1.0);
            for (row, (i, j, factor)) in [(0, 0, 1.0), (1, 1, 1.0), (0, 1, 2.0), (0, 2, 2.0), (1, 2, 2.0)].into_iter().enumerate() {
                b[row][k] = factor * contract(&e[k], m.axes[i], m.axes[j]);
            }
        }
        if b.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a nonfinite strain operator".into()));
        }
        Ok(b)
    }

    fn intervals(&self) -> Vec<(f64, f64)> {
        let h = self.section.thickness;
        if let Some(layers) = &self.layers {
            let mut lo = -1.0;
            layers
                .plies
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    let hi = if i + 1 == layers.plies.len() { 1.0 } else { lo + 2.0 * p.thickness / h };
                    let span = (lo, hi);
                    lo = hi;
                    span
                })
                .collect()
        } else {
            vec![(-1.0, 1.0)]
        }
    }

    fn thickness_rule(&self) -> Vec<(f64, f64, usize)> {
        self.intervals()
            .iter()
            .enumerate()
            .flat_map(|(i, &(lo, hi))| gauss_legendre(3).iter().map(move |&(x, w)| (0.5 * (lo + hi) + 0.5 * (hi - lo) * x, 0.5 * (hi - lo) * w, i)))
            .collect()
    }

    fn ply_at(&self, t: f64) -> usize {
        self.intervals().iter().position(|&(_, hi)| t <= hi + 1e-14).unwrap_or(0)
    }

    /// Actual fibre angle in the current integration-point strain frame. The
    /// section reference direction is projected at the reference surface, then
    /// held along its director through thickness; it does not rotate with g_r(t).
    fn ply_angle(&self, p: [f64; 3], ply: usize) -> Result<f64> {
        let Some(layers) = &self.layers else { return Ok(0.0) };
        let input = layers.plies.get(ply).ok_or_else(|| Error::Invalid("MITC8/MITC6-a ply index".into()))?;
        let mid = self.map([p[0], p[1], 0.0])?;
        let sh = Self::basis().eval(p);
        let position = std::array::from_fn(|j| (0..N).map(|i| sh.n[i] * self.xyz[i][j]).sum());
        let axes = match &layers.orientation {
            Some(o) => o.shell_axes_at(position, mid.axes[2], layers.reference_axis).map_err(|e| Error::Invalid(e.to_string()))?,
            None => mid.axes,
        };
        let (sn, cs) = libm::sincos(input.angle);
        let fibre = sum(scaled(axes[0], cs), scaled(axes[1], sn));
        let current = self.map(p)?;
        Ok(libm::atan2(dot(fibre, current.axes[1]), dot(fibre, current.axes[0])))
    }

    fn material(&self, p: [f64; 3], ply: usize) -> Result<[[f64; 5]; 5]> {
        let (law, angle) = if let Some(layers) = &self.layers {
            let input = layers.plies.get(ply).ok_or_else(|| Error::Invalid("MITC8/MITC6-a ply index".into()))?;
            (Lamina::new(&input.elastic)?, self.ply_angle(p, ply)?)
        } else {
            (Lamina::new(&fieldcraft_model::Elastic::Isotropic { young: self.section.young, poisson: self.section.poisson })?, 0.0)
        };
        let q = law.rotated(angle)?;
        let shear = law.rotated_shear(angle)?;
        let mut c = [[0.0; 5]; 5];
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] = q[i][j];
            }
        }
        for i in 0..2 {
            for j in 0..2 {
                c[i + 3][j + 3] = shear[i][j];
            }
        }
        Ok(c)
    }

    fn drilling(&self, p: [f64; 3]) -> Result<([f64; DOFS], f64)> {
        let sh = Self::basis().eval(p);
        let mut g = [[0.0; 3]; 2];
        for i in 0..N {
            for a in 0..2 {
                g[a] = sum(g[a], scaled(sub(self.xyz[i], self.xyz[0]), sh.dn[i][a]));
            }
        }
        let area = norm(cross(g[0], g[1]));
        let n = unit(cross(g[0], g[1])).ok_or_else(|| Error::Degenerate("MITC8/MITC6-a drilling surface".into()))?;
        let d = [scaled(cross(g[1], n), 1.0 / area), scaled(cross(n, g[0]), 1.0 / area)];
        let mut row = [0.0; DOFS];
        for i in 0..N {
            let grad = sum(scaled(d[0], sh.dn[i][0]), scaled(d[1], sh.dn[i][1]));
            for k in 0..3 {
                let mut v = [0.0; 3];
                v[k] = 1.0;
                row[6 * i + k] = -0.5 * dot(n, cross(grad, v));
                row[6 * i + 3 + k] = sh.n[i] * n[k];
            }
        }
        Ok((row, area))
    }

    pub fn stiffness(&self) -> Result<Matrix> {
        let mut k = Matrix::zeros(DOFS);
        let rule = self.thickness_rule();
        for q in Self::surface_points() {
            for &(t, tw, ply) in &rule {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let b = self.strain_operator(p, true)?;
                let c = self.material(p, ply)?;
                let cb: StrainOperator<DOFS> = std::array::from_fn(|a| std::array::from_fn(|i| (0..5).map(|j| c[a][j] * b[j][i]).sum()));
                let w = q.weight * tw * m.det;
                for i in 0..DOFS {
                    for j in i..DOFS {
                        let v = w * (0..5).map(|a| b[a][i] * cb[a][j]).sum::<f64>();
                        k.add(i, j, v);
                        if i != j {
                            k.add(j, i, v);
                        }
                    }
                }
            }
            let (b, area) = self.drilling(q.xi)?;
            let penalty = self.layers.as_ref().map_or(DRILLING * self.section.shear_modulus() * self.section.thickness, |l| l.drilling);
            let w = penalty * area * q.weight;
            for i in 0..DOFS {
                for j in i..DOFS {
                    let v = w * b[i] * b[j];
                    k.add(i, j, v);
                    if i != j {
                        k.add(j, i, v);
                    }
                }
            }
        }
        if k.data.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a nonfinite stiffness".into()));
        }
        Ok(k)
    }

    /// Consistent kinetic energy of u=ΣN(u_i+z theta_i×n_i), with physical
    /// translation/rotation coupling integrated through each density layer.
    /// Lee & Bathe (2004), eqs. (1)-(8), inserted in T=∫rho |u_dot|²/2 dV.
    /// As in the flat shell, add ∫rho z² (ΣN theta_i·n_i)²/2: artificial
    /// drilling inertia, not part of the five-DOF shell. Three Gauss points per
    /// direction and per ply; no MITC strain interpolation is used for inertia.
    pub fn mass(&self, density: f64) -> Result<Matrix> {
        let mut out = Matrix::zeros(DOFS);
        for q in Self::surface_points() {
            for (t, tw, ply) in self.thickness_rule() {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let sh = Self::basis().eval(p);
                let z = self.section.eccentricity + 0.5 * self.section.thickness * t;
                let w = self.density(density, ply)? * q.weight * tw * m.det;
                let v: [Vec3; DOFS] = std::array::from_fn(|i| {
                    let mut a = [0.0; 3];
                    a[i % 3] = sh.n[i / 6];
                    if i % 6 < 3 { a } else { scaled(cross(a, self.directors[i / 6]), z) }
                });
                let drill: [f64; DOFS] = std::array::from_fn(|i| if i % 6 < 3 { 0.0 } else { z * sh.n[i / 6] * self.directors[i / 6][i % 3] });
                for i in 0..DOFS {
                    for j in i..DOFS {
                        let a = w * (dot(v[i], v[j]) + drill[i] * drill[j]);
                        out.add(i, j, a);
                        if i != j {
                            out.add(j, i, a);
                        }
                    }
                }
            }
        }
        finite_matrix(out, "mass")
    }

    fn density(&self, homogeneous: f64, ply: usize) -> Result<f64> {
        let rho = match &self.layers {
            Some(l) => l
                .plies
                .get(ply)
                .and_then(|p| p.density)
                .ok_or_else(|| Error::Invalid("every MITC8/MITC6-a ply needs density for inertia or gravity".into()))?,
            None => homogeneous,
        };
        if !(rho > 0.0 && rho.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a density must be positive and finite".into()));
        }
        Ok(rho)
    }

    /// Positive diagonal scaling (Hinton, Rock & Zienkiewicz 1976,
    /// doi:10.1002/eqe.4290040305), separately conserving ∫rho and ∫rho z².
    /// Omits I1 and all inter-node coupling; not row-sum lumping (quad8 row
    /// sums can be negative), and not an exact rigid-rotation energy integral.
    pub fn lumped_mass(&self, density: f64) -> Result<Matrix> {
        let mut diagonal = [[0.0; N]; 2];
        let mut total = [0.0; 2];
        for q in Self::surface_points() {
            for (t, tw, ply) in self.thickness_rule() {
                let p = [q.xi[0], q.xi[1], t];
                let sh = Self::basis().eval(p);
                let z = self.section.eccentricity + 0.5 * self.section.thickness * t;
                let w = self.density(density, ply)? * q.weight * tw * self.map(p)?.det;
                for (a, moment) in [w, w * z * z].into_iter().enumerate() {
                    total[a] += moment;
                    for i in 0..N {
                        diagonal[a][i] += moment * sh.n[i] * sh.n[i];
                    }
                }
            }
        }
        let mut out = Matrix::zeros(DOFS);
        for a in 0..2 {
            let sum: f64 = diagonal[a].iter().sum();
            if !(sum > 0.0 && total[a] > 0.0) {
                return Err(Error::Invalid("MITC8/MITC6-a mass moment underflow".into()));
            }
            for i in 0..N {
                for c in 0..3 {
                    out.add(6 * i + 3 * a + c, 6 * i + 3 * a + c, diagonal[a][i] * total[a] / sum);
                }
            }
        }
        finite_matrix(out, "lumped mass")
    }

    pub fn is_flat(&self) -> Result<bool> {
        let normal = self.plate_normal()?;
        let extent = self.xyz.iter().map(|x| norm(sub(*x, self.xyz[0]))).fold(0.0, f64::max);
        Ok(self.xyz.iter().all(|x| dot(sub(*x, self.xyz[0]), normal).abs() <= 1e-8 * extent)
            && self.directors.iter().all(|n| dot(*n, normal) >= 1.0 - 1e-8))
    }

    /// Reference normal for the flat-plate applicability check. Curvature and
    /// the preload are checked by the solver before assembling plate K_sigma.
    pub fn plate_normal(&self) -> Result<Vec3> {
        Ok(self.map(Self::centre())?.axes[2])
    }

    /// Small-slope flat-plate initial-stress operator, Bathe (2014), §6.8.2:
    /// ∫ grad(w)^T N grad(w) dA; also NASA TP-2528 (1986), eq. (12).
    /// Membrane resultants come from the actual
    /// MITC preload strain and layer laws at each integration point. No curved
    /// shell tangent, director geometric terms, or follower-load stiffness.
    pub fn plate_geometric_stiffness(&self, u: &[f64]) -> Result<Matrix> {
        if u.len() != DOFS || u.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a plate preload has wrong or nonfinite degrees of freedom".into()));
        }
        let normal = self.plate_normal()?;
        let extent = self.xyz.iter().map(|x| norm(sub(*x, self.xyz[0]))).fold(0.0, f64::max);
        if self.xyz.iter().any(|x| dot(sub(*x, self.xyz[0]), normal).abs() > 1e-8 * extent)
            || self.directors.iter().any(|n| dot(*n, normal) < 1.0 - 1e-8)
        {
            return Err(Error::NotSupported("MITC8/MITC6-a curved-shell buckling: only coplanar plates with flat directors are supported".into()));
        }
        let mut out = Matrix::zeros(DOFS);
        for q in Self::surface_points() {
            for (t, tw, ply) in self.thickness_rule() {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let b = self.strain_operator(p, true)?;
                let strain: [f64; 5] = std::array::from_fn(|i| b[i].iter().zip(u).map(|(b, u)| b * u).sum());
                let c = self.material(p, ply)?;
                let stress: [f64; 3] = std::array::from_fn(|i| (0..5).map(|j| c[i][j] * strain[j]).sum());
                let sh = Self::basis().eval(p);
                let g: [[f64; 2]; N] = std::array::from_fn(|i| {
                    let v = sum(scaled(m.dual[0], sh.dn[i][0]), scaled(m.dual[1], sh.dn[i][1]));
                    [dot(v, m.axes[0]), dot(v, m.axes[1])]
                });
                for i in 0..N {
                    for j in 0..N {
                        let v = q.weight
                            * tw
                            * m.det
                            * (g[i][0] * (stress[0] * g[j][0] + stress[2] * g[j][1]) + g[i][1] * (stress[2] * g[j][0] + stress[1] * g[j][1]));
                        for a in 0..3 {
                            for b in 0..3 {
                                out.add(6 * i + a, 6 * j + b, v * normal[a] * normal[b]);
                            }
                        }
                    }
                }
            }
        }
        finite_matrix(out, "plate initial-stress stiffness")
    }

    /// Consistent reference-configuration initial-stress Hessian.
    /// Bathe-Dvorkin (1986), (2)-(4): E_ab=(g_a.g_b-G_a.G_b)/2,
    /// KG_ij=integral S:delta_i delta_j E dV. The same MITC tensor projection
    /// as the elastic operator acts on the second variation, including transverse
    /// shear, thickness geometry and the unit-director second derivative.
    /// This is a linear bifurcation operator, not a nonlinear equilibrium solver.
    /// The external load potential must be dead; no follower-pressure Hessian.
    pub fn geometric_stiffness(&self, u: &[f64]) -> Result<Matrix> {
        if u.len() != DOFS || u.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a geometric preload coordinates".into()));
        }
        let mut out = Matrix::zeros(DOFS);
        for q in Self::surface_points() {
            for (t, tw, ply) in self.thickness_rule() {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let b = self.strain_operator(p, true)?;
                let strain: [f64; 5] = std::array::from_fn(|a| b[a].iter().zip(u).map(|(b, u)| b * u).sum());
                let c = self.material(p, ply)?;
                let stress: [f64; 5] = std::array::from_fn(|a| (0..5).map(|b| c[a][b] * strain[b]).sum());
                let w = q.weight * tw * m.det;
                for i in 0..DOFS {
                    let h = self.strain_derivative(p, true, Some(i))?;
                    for j in i..DOFS {
                        let v = w * (0..5).map(|a| stress[a] * h[a][j]).sum::<f64>();
                        out.add(i, j, v);
                        if i != j {
                            out.add(j, i, v);
                        }
                    }
                }
            }
        }
        finite_matrix(out, "curved initial-stress stiffness")
    }

    pub fn body_force(&self, force: Vec3) -> Result<Vec<f64>> {
        self.volume_load(force, None)
    }

    /// Static weight integrates each ply density, including its first moment;
    /// this is not a mass operator for dynamic analyses.
    pub fn gravity(&self, density: f64, acceleration: Vec3) -> Result<Vec<f64>> {
        if !(density > 0.0 && density.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a gravity needs positive finite density".into()));
        }
        self.volume_load(acceleration, Some(density))
    }

    fn volume_load(&self, force: Vec3, density: Option<f64>) -> Result<Vec<f64>> {
        if force.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a body force must be finite".into()));
        }
        let mut out = vec![0.0; DOFS];
        for q in Self::surface_points() {
            for (t, tw, ply) in self.thickness_rule() {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let sh = Self::basis().eval(p);
                let z = self.section.eccentricity + 0.5 * self.section.thickness * p[2];
                let rho = match (density, &self.layers) {
                    (Some(_), Some(layers)) => layers
                        .plies
                        .get(ply)
                        .and_then(|p| p.density)
                        .ok_or_else(|| Error::Invalid("every MITC8/MITC6-a ply needs density for gravity".into()))?,
                    (Some(rho), None) => rho,
                    (None, _) => 1.0,
                };
                let w = q.weight * tw * m.det * rho;
                for i in 0..N {
                    for j in 0..3 {
                        let mut v = [0.0; 3];
                        v[j] = 1.0;
                        out[6 * i + j] += w * sh.n[i] * force[j];
                        out[6 * i + 3 + j] += w * z * sh.n[i] * dot(cross(v, self.directors[i]), force);
                    }
                }
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a nonfinite body force".into()));
        }
        Ok(out)
    }

    /// Direct nodal stress samples, as in BD86 section 4.4, not extrapolated or
    /// equilibrated interlaminar stresses. Global stresses fit the existing result
    /// layout. Resultants are projected onto the single reported centre frame.
    pub fn stresses(&self, u: &[f64]) -> Result<ShellStress> {
        if u.len() != DOFS || u.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("quadratic shell has wrong or nonfinite displacement/rotation values".into()));
        }
        let axes = self.map(Self::centre())?.axes;
        let mut out = ShellStress { axes, forces: Vec::new(), bottom: Vec::new(), middle: Vec::new(), top: Vec::new(), plies: Vec::new() };
        let stress = |p: [f64; 3], ply: usize| -> Result<Tensor> {
            let c = self.material(p, ply)?;
            let b = self.strain_operator(p, false)?;
            let strain: [f64; 5] = std::array::from_fn(|i| b[i].iter().zip(u).map(|(x, y)| x * y).sum());
            let s: [f64; 5] = std::array::from_fn(|i| (0..5).map(|j| c[i][j] * strain[j]).sum());
            let m = self.map(p)?;
            let global = tensor_to_global(&m.axes, [[s[0], s[2], s[3]], [s[2], s[1], s[4]], [s[3], s[4], 0.0]]);
            Ok([[global[0], global[3], global[4]], [global[3], global[1], global[5]], [global[4], global[5], global[2]]])
        };
        for q in Self::basis().nodes() {
            for (t, values) in [(-1.0, &mut out.bottom), (0.0, &mut out.middle), (1.0, &mut out.top)] {
                let s = stress([q[0], q[1], t], self.ply_at(t))?;
                values.push([s[0][0], s[1][1], s[2][2], s[0][1], s[0][2], s[1][2]]);
            }
            let mut f = [0.0; 8];
            for (t, tw, ply) in self.thickness_rule() {
                let s = stress([q[0], q[1], t], ply)?;
                let w = 0.5 * self.section.thickness * tw;
                let z = self.section.eccentricity + 0.5 * self.section.thickness * t;
                for (i, (a, b)) in [(0, 0), (1, 1), (0, 1)].into_iter().enumerate() {
                    let v = contract(&s, axes[a], axes[b]);
                    f[i] += w * v;
                    f[i + 3] += w * z * v;
                }
                f[6] += w * contract(&s, axes[0], axes[2]);
                f[7] += w * contract(&s, axes[1], axes[2]);
            }
            out.forces.push(f);
            if let Some(layers) = &self.layers {
                let mut fields = Vec::new();
                for (i, (&(lo, hi), input)) in self.intervals().iter().zip(&layers.plies).enumerate() {
                    let at = |t: f64| {
                        let p = [q[0], q[1], t];
                        let b = self.strain_operator(p, false)?;
                        let e = std::array::from_fn(|k| b[k].iter().zip(u).map(|(a, b)| a * b).sum());
                        let z = self.section.eccentricity + 0.5 * self.section.thickness * t;
                        Lamina::new(&input.elastic)?.recover_strain(e, self.ply_angle(p, i)?, z, input.strength)
                    };
                    fields.push(PlyFields {
                        angle: self.ply_angle([q[0], q[1], 0.0], i)?,
                        bottom: at(lo)?,
                        middle: at(0.5 * (lo + hi))?,
                        top: at(hi)?,
                        points: section_abscissae(input.points)?.iter().map(|x| at(0.5 * (lo + hi) + 0.5 * (hi - lo) * x)).collect::<Result<_>>()?,
                    });
                }
                out.plies.push(fields);
            }
        }
        if out.bottom.iter().chain(&out.middle).chain(&out.top).flatten().chain(out.forces.iter().flatten()).any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8/MITC6-a nonfinite stress recovery".into()));
        }
        Ok(out)
    }
}

fn finite_matrix(m: Matrix, name: &str) -> Result<Matrix> {
    if m.data.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid(format!("MITC8/MITC6-a nonfinite {name}")));
    }
    Ok(m)
}

#[cfg(test)]
mod geometric_tests {
    use super::*;
    fn check<const N: usize, const D: usize>() {
        let basis = if N == 8 { Basis::Quad8 } else { Basis::Tri6 };
        let x: Vec<_> = basis.nodes().iter().map(|p| [p[0], p[1], 0.12 * p[0] * p[0] + 0.08 * p[0] * p[1]]).collect();
        let section = ShellSection::new(0.1, 0.2, &fieldcraft_model::Elastic::Isotropic { young: 1200.0, poisson: 0.25 }).unwrap();
        let s = QuadraticShell::<N, D>::new(&x, section).unwrap();
        let p = [0.2, 0.3, 0.4];
        let m = s.map(p).unwrap();
        let v: Vec<_> = (0..D).map(|i| (i as f64 * 0.7).sin() * 0.3).collect();
        let w: Vec<_> = (0..D).map(|i| (i as f64 * 0.4).cos() * 0.2).collect();
        let metric = |a: f64, b: f64| {
            let mut t = s.clone();
            for n in 0..N {
                let u: Vec3 = std::array::from_fn(|j| a * v[6 * n + j] + b * w[6 * n + j]);
                let r: Vec3 = std::array::from_fn(|j| a * v[6 * n + 3 + j] + b * w[6 * n + 3 + j]);
                t.xyz[n] = sum(t.xyz[n], u);
                t.directors[n] = unit(sum(s.directors[n], cross(r, s.directors[n]))).unwrap();
            }
            t.map(p).unwrap().g
        };
        let eps = 1e-4;
        let g = [metric(eps, eps), metric(eps, -eps), metric(-eps, eps), metric(-eps, -eps)];
        for a in 0..3 {
            for b in a..3 {
                let analytic: f64 = (0..D).map(|i| v[i] * s.covariant(p, &m, a, b, Some(i)).iter().zip(&w).map(|(h, w)| h * w).sum::<f64>()).sum();
                let measured = (dot(g[0][a], g[0][b]) - dot(g[1][a], g[1][b]) - dot(g[2][a], g[2][b]) + dot(g[3][a], g[3][b])) / (8.0 * eps * eps);
                assert!((analytic - measured).abs() < 2e-8, "{N}: ({a},{b}) {analytic} != {measured}");
            }
        }
        let preload: Vec<_> = x.iter().flat_map(|p| [0.0, -0.001 * p[1], 0.0, 0.0, 0.0, 0.0]).collect();
        let kg = s.geometric_stiffness(&preload).unwrap();
        for c in 0..3 {
            let rigid: Vec<_> = (0..D).map(|i| if i % 6 == c { 1.0 } else { 0.0 }).collect();
            assert!(kg.mul_vec(&rigid).iter().all(|x| x.abs() < 1e-10));
        }
        for i in 0..D {
            for j in 0..D {
                assert_eq!(kg.data[i * D + j], kg.data[j * D + i]);
            }
        }
        let zero = s.geometric_stiffness(&vec![0.0; D]).unwrap();
        assert!(zero.data.iter().all(|v| *v == 0.0));
    }
    #[test]
    fn curved_green_hessian_includes_unit_director_second_variation() {
        check::<8, 48>();
        check::<6, 36>();
    }
}
