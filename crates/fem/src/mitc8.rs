//! Eight-node curved MITC8, small-displacement homogeneous and layered elastic statics.
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
//! reference-surface offset as `ShellSection`. No finite rotations or mass model.

use crate::basis::Basis;
use crate::dense::{Matrix, inverse3};
use crate::frame::{Axes, Vec3, cross, dot, norm, sub, tensor_to_global, unit};
use crate::laminate::{Lamina, Laminate, PlyFields, PlyProperties, section_abscissae};
use crate::quadrature::{Rule, gauss_legendre};
use crate::shell::{DRILLING, ShellSection, ShellStress};
use crate::{Error, Result};
use fieldcraft_model::MaterialOrientation;

const N: usize = 8;
const DOFS: usize = 48;
type Tensor = [[f64; 3]; 3];
type TensorRows = [Tensor; DOFS];
/// Engineering strain rows: 11,22,12,13,23 in the integration point's orthonormal frame.
pub type StrainOperator = [[f64; DOFS]; 5];

#[derive(Debug, Clone)]
pub struct Mitc8 {
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
    let (inverse, _) = inverse3(g).ok_or_else(|| Error::Degenerate("MITC8 singular covariant basis".into()))?;
    Ok(std::array::from_fn(|i| std::array::from_fn(|j| inverse[j][i])))
}
fn check_point(p: [f64; 3]) -> Result<()> {
    if p.iter().any(|v| !v.is_finite() || v.abs() > 1.0 + 1e-12) {
        return Err(Error::Invalid("MITC8 natural coordinates must be finite and in [-1,1]".into()));
    }
    Ok(())
}

impl Mitc8 {
    pub fn new(xyz: &[[f64; 3]], section: ShellSection) -> Result<Self> {
        let xyz: [Vec3; N] = xyz.try_into().map_err(|_| Error::Invalid("MITC8 requires exactly eight nodes".into()))?;
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        // ShellSection's fields are public; validate them even for direct kernel callers.
        ShellSection::new(
            section.thickness,
            -section.eccentricity / section.thickness,
            &fieldcraft_model::Elastic::Isotropic { young: section.young, poisson: section.poisson },
        )?;
        if !(section.young > 0.0 && section.young.is_finite() && section.poisson > -1.0 && section.poisson < 0.5) {
            return Err(Error::Invalid("MITC8 needs a finite positive modulus and -1 < poisson < 0.5".into()));
        }
        let mut directors = [[0.0; 3]; N];
        for (i, p) in Basis::Quad8.nodes().iter().enumerate() {
            let sh = Basis::Quad8.eval(*p);
            let tangent: [Vec3; 2] = std::array::from_fn(|a| std::array::from_fn(|k| (0..N).map(|j| sh.dn[j][a] * (xyz[j][k] - xyz[0][k])).sum()));
            directors[i] = unit(cross(tangent[0], tangent[1])).ok_or_else(|| Error::Degenerate("MITC8 zero nodal normal".into()))?;
        }
        let centre = Basis::Quad8.eval([0.0; 3]);
        let tangent: [Vec3; 2] = std::array::from_fn(|a| std::array::from_fn(|k| (0..N).map(|j| centre.dn[j][a] * (xyz[j][k] - xyz[0][k])).sum()));
        let normal = unit(cross(tangent[0], tangent[1])).ok_or_else(|| Error::Degenerate("MITC8 singular centre surface".into()))?;
        if directors.iter().any(|n| dot(*n, normal) <= 1e-8) {
            return Err(Error::Degenerate("MITC8 folded surface or nodal normals turning 90 degrees from the centre; refine the geometry".into()));
        }
        let out = Self { xyz, directors, section, layers: None };
        // Geometry screen, including both thickness faces and the reference nodes. All
        // quadrature and tying locations are checked again when used. This is a sampled
        // Jacobian check, not a global positivity certificate for an arbitrary curved map.
        for p in Basis::Quad8.nodes().iter().chain(&[[0.0; 3]]) {
            for t in [-1.0, 0.0, 1.0] {
                out.map([p[0], p[1], t])?;
            }
        }
        for p in Rule::Quad(3).points() {
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
            return Err(Error::Invalid("MITC8 non-default reference axis needs a material orientation".into()));
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
                return Err(Error::Invalid("MITC8 shared directors must be unit vectors within 60 degrees of local normals".into()));
            }
        }
        self.directors = directors;
        for p in Basis::Quad8.nodes().iter().chain(&[[0.0; 3]]) {
            for t in [-1.0, 0.0, 1.0] {
                self.map([p[0], p[1], t])?;
            }
        }
        for p in Rule::Quad(3).points() {
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
        let sh = Basis::Quad8.eval(p);
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
            return Err(Error::Degenerate("MITC8 inverted, folded or singular thickness Jacobian".into()));
        }
        let n = unit(g[2]).ok_or_else(|| Error::Degenerate("MITC8 zero director".into()))?;
        let e1 = unit(sub(g[0], scaled(n, dot(g[0], n)))).ok_or_else(|| Error::Degenerate("MITC8 singular tangent frame".into()))?;
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

    fn covariant(&self, p: [f64; 3], m: &Map, a: usize, b: usize) -> [f64; DOFS] {
        let sh = Basis::Quad8.eval(p);
        let z = self.section.eccentricity + 0.5 * self.section.thickness * p[2];
        std::array::from_fn(|k| {
            let d = self.displacement_gradient(&sh, z, k);
            0.5 * (dot(m.g[a], d[b]) + dot(m.g[b], d[a]))
        })
    }

    /// Direct in-layer tensor, BD86 equations (10),(12); membrane and bending
    /// coexist through z. Normal strain is eliminated by the plane-stress law.
    fn in_layer(&self, p: [f64; 3]) -> Result<TensorRows> {
        let m = self.map(p)?;
        let rr = self.covariant(p, &m, 0, 0);
        let ss = self.covariant(p, &m, 1, 1);
        let rs = self.covariant(p, &m, 0, 1);
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
    fn assumed_layer(&self, p: [f64; 3]) -> Result<TensorRows> {
        let a = 1.0 / 3.0f64.sqrt();
        let points = Basis::Quad8.nodes();
        let mut samples = [[[0.0; 3]; 3]; DOFS * 8];
        for i in 0..4 {
            let q = [points[i][0] * a, points[i][1] * a, p[2]];
            let e = self.in_layer(q)?;
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
            let direct = self.covariant(q, &m, keep, keep);
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
        let weights = Basis::Quad8.eval([p[0] / a, p[1] / a, 0.0]);
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
    fn shear(&self, p: [f64; 3], assumed: bool) -> Result<TensorRows> {
        let mut out = [[[0.0; 3]; 3]; DOFS];
        let a = 1.0 / 3.0f64.sqrt();
        for along in 0..2 {
            if !assumed {
                let m = self.map(p)?;
                let v = self.covariant(p, &m, along, 2);
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
                let v = self.covariant(q, &m, along, 2);
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
                let v = self.covariant(q, &m, along, 2);
                for k in 0..DOFS {
                    add_tensor(&mut out[k], &sym_dyad(centre.dual[along], centre.dual[2]), 0.5 * bubble * v[k]);
                }
            }
        }
        Ok(out)
    }

    /// Physical engineering-strain operator at a natural point; `assumed=false`
    /// is the direct displacement field (for independent diagnostics/recovery).
    pub fn strain_operator(&self, p: [f64; 3], assumed: bool) -> Result<StrainOperator> {
        let m = self.map(p)?;
        let mut e = if assumed { self.assumed_layer(p)? } else { self.in_layer(p)? };
        let shear = self.shear(p, assumed)?;
        let mut b = [[0.0; DOFS]; 5];
        for k in 0..DOFS {
            add_tensor(&mut e[k], &shear[k], 1.0);
            for (row, (i, j, factor)) in [(0, 0, 1.0), (1, 1, 1.0), (0, 1, 2.0), (0, 2, 2.0), (1, 2, 2.0)].into_iter().enumerate() {
                b[row][k] = factor * contract(&e[k], m.axes[i], m.axes[j]);
            }
        }
        if b.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8 nonfinite strain operator".into()));
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
        let input = layers.plies.get(ply).ok_or_else(|| Error::Invalid("MITC8 ply index".into()))?;
        let mid = self.map([p[0], p[1], 0.0])?;
        let sh = Basis::Quad8.eval(p);
        let position = std::array::from_fn(|j| (0..8).map(|i| sh.n[i] * self.xyz[i][j]).sum());
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
            let input = layers.plies.get(ply).ok_or_else(|| Error::Invalid("MITC8 ply index".into()))?;
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
        let sh = Basis::Quad8.eval(p);
        let mut g = [[0.0; 3]; 2];
        for i in 0..N {
            for a in 0..2 {
                g[a] = sum(g[a], scaled(sub(self.xyz[i], self.xyz[0]), sh.dn[i][a]));
            }
        }
        let area = norm(cross(g[0], g[1]));
        let n = unit(cross(g[0], g[1])).ok_or_else(|| Error::Degenerate("MITC8 drilling surface".into()))?;
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
        for q in Rule::Quad(3).points() {
            for &(t, tw, ply) in &rule {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let b = self.strain_operator(p, true)?;
                let c = self.material(p, ply)?;
                let cb: StrainOperator = std::array::from_fn(|a| std::array::from_fn(|i| (0..5).map(|j| c[a][j] * b[j][i]).sum()));
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
            return Err(Error::Invalid("MITC8 nonfinite stiffness".into()));
        }
        Ok(k)
    }

    pub fn body_force(&self, force: Vec3) -> Result<Vec<f64>> {
        self.volume_load(force, None)
    }

    /// Static weight integrates each ply density, including its first moment;
    /// this is not a mass operator for dynamic analyses.
    pub fn gravity(&self, density: f64, acceleration: Vec3) -> Result<Vec<f64>> {
        if !(density > 0.0 && density.is_finite()) {
            return Err(Error::Invalid("MITC8 gravity needs positive finite density".into()));
        }
        self.volume_load(acceleration, Some(density))
    }

    fn volume_load(&self, force: Vec3, density: Option<f64>) -> Result<Vec<f64>> {
        if force.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8 body force must be finite".into()));
        }
        let mut out = vec![0.0; DOFS];
        for q in Rule::Quad(3).points() {
            for (t, tw, ply) in self.thickness_rule() {
                let p = [q.xi[0], q.xi[1], t];
                let m = self.map(p)?;
                let sh = Basis::Quad8.eval(p);
                let z = self.section.eccentricity + 0.5 * self.section.thickness * p[2];
                let rho = match (density, &self.layers) {
                    (Some(_), Some(layers)) => layers
                        .plies
                        .get(ply)
                        .and_then(|p| p.density)
                        .ok_or_else(|| Error::Invalid("every MITC8 ply needs density for gravity".into()))?,
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
            return Err(Error::Invalid("MITC8 nonfinite body force".into()));
        }
        Ok(out)
    }

    /// Direct nodal stress samples, as in BD86 section 4.4, not extrapolated or
    /// equilibrated interlaminar stresses. Global stresses fit the existing result
    /// layout. Resultants are projected onto the single reported centre frame.
    pub fn stresses(&self, u: &[f64]) -> Result<ShellStress> {
        if u.len() != DOFS || u.iter().any(|v| !v.is_finite()) {
            return Err(Error::Invalid("MITC8 needs 48 finite displacement/rotation values".into()));
        }
        let axes = self.map([0.0; 3])?.axes;
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
        for q in Basis::Quad8.nodes() {
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
            return Err(Error::Invalid("MITC8 nonfinite stress recovery".into()));
        }
        Ok(out)
    }
}
