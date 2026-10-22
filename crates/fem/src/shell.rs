//! Flat shell elements: the four-node MITC4 quadrilateral and the three-node MITC3 triangle,
//! with six degrees of freedom per node (three translations, three rotations).
//!
//! Each element is flat: its nodes are projected onto the element's mean plane, where the
//! membrane, the Reissner–Mindlin plate and a drilling stiffness are superposed in local axes and
//! the result rotated to the global ones (the facet approach; Zienkiewicz, Taylor & Fox 2014,
//! *The Finite Element Method for Solid and Structural Mechanics*, 7th ed., ch. 9–10). Curved
//! shells are approximated by facets; quadrilaterals on a cylinder whose edges follow the
//! generators and the circumference are exactly flat.
//!
//! - **Kinematics.** With the reference surface at the nodes and z along the element normal,
//!   u = u₀ + z θ_y, v = v₀ − z θ_x, w = w₀: membrane strains ε = (u₀,x, v₀,y, u₀,y + v₀,x),
//!   curvatures κ = (θ_y,x, −θ_x,y, θ_y,y − θ_x,x) and transverse shear strains
//!   γ = (w,x + θ_y, w,y − θ_x).
//! - **Section.** Linear elastic and homogeneous through the thickness, so the integration
//!   through the thickness is done exactly in closed form (any number of section points gives the
//!   same stiffness): with the mid-surface at z = e from the reference surface (e = −offset · t),
//!   N = A ε + B κ, M = B ε + D κ with A = Q t, B = Q t e, D = Q (t³/12 + t e²) for the plane
//!   stress matrix Q, and Q_s = k G t γ with the shear correction k = 5/6.
//! - **Transverse shear, quadrilateral (MITC4).** The covariant transverse shear strains are
//!   sampled at the mid-sides and interpolated linearly across the element: ẽ_rt from the points
//!   (0, ±1), ẽ_st from (±1, 0), then turned into Cartesian strains with the Jacobian
//!   (E. N. Dvorkin, K.-J. Bathe, *A continuum mechanics based four-node shell element for
//!   general nonlinear analysis*, Engineering Computations 1(1):77–88, 1984; in the flat-plate
//!   form of K.-J. Bathe, E. N. Dvorkin, *A four-node plate bending element based on
//!   Mindlin/Reissner plate theory and a mixed interpolation*, IJNME 21(2):367–383, 1985). This
//!   removes shear locking while passing the patch tests; membrane and bending use the bilinear
//!   fields, all integrated with 2 × 2 Gauss points.
//! - **Transverse shear, triangle (MITC3).** The covariant shear strains are tied at the three
//!   mid-sides: ẽ_rt = e_rt⁽¹⁾ + c s, ẽ_st = e_st⁽²⁾ − c r with c = e_rt⁽³⁾ − e_rt⁽¹⁾ − e_st⁽³⁾ +
//!   e_st⁽²⁾, tying points (½, 0), (0, ½) and (½, ½), which keeps the tangential shear strain of
//!   each edge constant along it (P.-S. Lee, K.-J. Bathe, *Development of MITC isotropic
//!   triangular shell finite elements*, Computers & Structures 82(11–13):945–962, 2004). The
//!   membrane is the constant-strain triangle.
//! - **Drilling rotation.** The rotation about the normal has no stiffness in shell theory; it
//!   is tied to the membrane's in-plane rotation ω = ½ (v,x − u,y) by the penalty
//!   ∫ γ t (θ_z − ω)² dA of the variational formulation of T. J. R. Hughes, F. Brezzi, *On drilling
//!   degrees of freedom*, CMAME 72(1):105–121, 1989, with γ = G/100. Hughes and Brezzi take
//!   γ = G with a drilling-enhanced (Allman-type) membrane; on the plain bilinear membrane the
//!   in-plane rotation of a parallelogram is linear, which θ_z interpolates exactly, so the penalty
//!   does not lock, and the smaller γ keeps distorted elements free of it too. Being a function of
//!   θ_z − ω it leaves rigid motions free. Scaled with G t like the membrane (not with the bending
//!   rigidity), it ties each node's rotation about the surface normal to the membrane as the mesh
//!   is refined, so curved shells modelled with facets converge steadily: the Scordelis–Lo and
//!   pinched-cylinder benchmarks change by less than 0.5 % for γ from G/1000 to G, whereas a
//!   penalty that is a fixed fraction of the bending stiffness lets that rotation go soft on fine
//!   meshes.
//! - **Mass.** Consistent: ρ t on the translations, ρ (t³/12 + t e²) on the rotations and
//!   ρ t e coupling translations and bending rotations when the reference surface is offset. The
//!   drilling rotation gets the same rotary inertia as the bending ones (shell theory gives it
//!   none; it is tiny next to the in-plane inertia of the translations and keeps the mass matrix
//!   positive definite, as the frequency solver's Lanczos iteration needs). Lumped: the diagonal
//!   scaled to the element's mass (Hinton, Rock & Zienkiewicz 1976), the rotary inertia
//!   isotropic about the node so the lumped matrix stays diagonal in global axes.
//! - **Stresses.** Section forces and the stresses at the bottom, middle and top of the section
//!   at the integration points, extrapolated bilinearly to the nodes (quadrilateral) or taken at
//!   the nodes (triangle, constant membrane and bending fields). Transverse shear stresses are
//!   reported at the mid-surface as the peak of the parabolic distribution, 3/2 Q/t.

use fieldcraft_model::{Elastic, ElementKind};

use crate::dense::Matrix;
use crate::dof::{DofSet, ElementDofs};
use crate::frame::{Axes, cross, dot, norm, rotate_matrix, rotate_to_global, rotate_to_local, sub, tensor_to_global, to_local, translation, unit};
use crate::quadrature::Rule;
use crate::{Error, Result};

/// Transverse shear correction factor of a homogeneous section.
pub const SHEAR_CORRECTION: f64 = 5.0 / 6.0;

/// The drilling penalty γ as a fraction of the shear modulus.
pub const DRILLING: f64 = 1e-2;

/// Degrees of freedom per shell node: translations along local (or global) x, y, z, then
/// rotations about them.
pub const NODE_DOFS: usize = 6;

/// A shell section with an isotropic material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellSection {
    pub thickness: f64,
    /// The mid-surface's position along the element normal from the nodes (the reference
    /// surface): −offset · thickness.
    pub eccentricity: f64,
    pub young: f64,
    pub poisson: f64,
}

impl ShellSection {
    /// A section of the model's thickness and offset (a fraction of the thickness, 0.5 puts the
    /// nodes on the side the normal points to).
    pub fn new(thickness: f64, offset: f64, elastic: &Elastic) -> Result<ShellSection> {
        let (young, poisson) = crate::material::isotropic_constants(elastic, "shell")?;
        if !(thickness > 0.0 && thickness.is_finite() && offset.abs() <= 0.5) {
            return Err(Error::Invalid(format!(
                "a shell section needs a positive thickness and an offset in −0.5..=0.5, got {thickness} and {offset}"
            )));
        }
        Ok(ShellSection { thickness, eccentricity: -offset * thickness, young, poisson })
    }

    /// The plane stress matrix Q.
    fn plane_stress(&self) -> [[f64; 3]; 3] {
        let (e, nu) = (self.young, self.poisson);
        let c = e / (1.0 - nu * nu);
        [[c, c * nu, 0.0], [c * nu, c, 0.0], [0.0, 0.0, c * (1.0 - nu) / 2.0]]
    }

    /// The membrane–bending stiffness [[A, B], [B, D]], 6 × 6.
    fn abd(&self) -> [[f64; 6]; 6] {
        let q = self.plane_stress();
        let (t, e) = (self.thickness, self.eccentricity);
        let (a, b, d) = (t, t * e, t * t * t / 12.0 + t * e * e);
        let mut out = [[0.0; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = a * q[i][j];
                out[i][j + 3] = b * q[i][j];
                out[i + 3][j] = b * q[i][j];
                out[i + 3][j + 3] = d * q[i][j];
            }
        }
        out
    }

    pub fn shear_modulus(&self) -> f64 {
        self.young / (2.0 * (1.0 + self.poisson))
    }

    /// The transverse shear stiffness k G t.
    fn shear(&self) -> f64 {
        SHEAR_CORRECTION * self.shear_modulus() * self.thickness
    }
}

/// Section forces and stresses of a shell element, per node in the element's node order.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ShellStress {
    /// The element's local axes: e₁, e₂ in its plane, e₃ its normal (the right-hand normal of
    /// the node order).
    pub axes: Axes,
    /// N_x, N_y, N_xy, M_x, M_y, M_xy, Q_x, Q_y in the local axes, per unit length, about the
    /// reference surface.
    pub forces: Vec<[f64; 8]>,
    /// Stress (global xx, yy, zz, xy, xz, yz) at the bottom (−e₃ side), the middle and the top
    /// (+e₃ side) of the section.
    pub bottom: Vec<[f64; 6]>,
    pub middle: Vec<[f64; 6]>,
    pub top: Vec<[f64; 6]>,
}

/// The shape functions of the element's corners at (r, s) and their derivatives along r and s.
fn shape(kind: ElementKind, r: f64, s: f64) -> ([f64; 4], [[f64; 2]; 4]) {
    if kind == ElementKind::Tri3 {
        ([1.0 - r - s, r, s, 0.0], [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    } else {
        let c = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        let mut n = [0.0; 4];
        let mut d = [[0.0; 2]; 4];
        for (i, ci) in c.iter().enumerate() {
            n[i] = 0.25 * (1.0 + ci[0] * r) * (1.0 + ci[1] * s);
            d[i] = [0.25 * ci[0] * (1.0 + ci[1] * s), 0.25 * ci[1] * (1.0 + ci[0] * r)];
        }
        (n, d)
    }
}

/// The map's Jacobian at a point: J = [[x,r, y,r], [x,s, y,s]], its inverse and determinant.
struct Map {
    n: [f64; 4],
    /// ∂Nᵢ/∂r, ∂Nᵢ/∂s.
    dn: [[f64; 2]; 4],
    /// ∂Nᵢ/∂x, ∂Nᵢ/∂y.
    dx: [[f64; 2]; 4],
    j: [[f64; 2]; 2],
    inv: [[f64; 2]; 2],
    det: f64,
}

/// Most degrees of freedom of a shell element (quad4).
const MAX_DOFS: usize = 4 * NODE_DOFS;

/// A row of a strain operator over the element's local degrees of freedom (zero past them).
type Row = [f64; MAX_DOFS];

/// The covariant shear rows at the element's tying points, which depend on its geometry only:
/// quad4 e_rt at (0, 1), (0, −1) and e_st at (1, 0), (−1, 0); tri3 e_rt⁽¹⁾, e_st⁽²⁾ and the
/// combination c.
struct Tying([Row; 4]);

/// A flat shell element.
#[derive(Debug, Clone)]
pub struct Shell {
    kind: ElementKind,
    nodes: usize,
    axes: Axes,
    /// The nodes in the element plane (local x, y).
    local: [[f64; 2]; 4],
    section: ShellSection,
    area: f64,
}

impl Shell {
    pub fn new(kind: ElementKind, xyz: &[[f64; 3]], section: ShellSection) -> Result<Shell> {
        let nodes = match kind {
            ElementKind::Tri3 => 3,
            ElementKind::Quad4 => 4,
            ElementKind::Tri6 | ElementKind::Quad8 => {
                return Err(Error::NotSupported(format!("{} shells (the solver has quad4 and tri3 shells so far)", kind.name())));
            }
            _ => return Err(Error::Invalid(format!("{} is not a shell element", kind.name()))),
        };
        if xyz.len() != nodes {
            return Err(Error::Invalid(format!("a {} element has {nodes} nodes, got {}", kind.name(), xyz.len())));
        }
        fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
        let degenerate = || Error::Degenerate(format!("a {} shell element is degenerate (zero area, inverted or not convex)", kind.name()));
        let p = |i: usize| xyz.get(i).copied().unwrap_or_default();
        let (normal, along) = if nodes == 3 {
            (cross(sub(p(1), p(0)), sub(p(2), p(0))), sub(p(1), p(0)))
        } else {
            let along = sub([0, 1, 2].map(|k| p(1)[k] + p(2)[k]), [0, 1, 2].map(|k| p(0)[k] + p(3)[k]));
            (cross(sub(p(2), p(0)), sub(p(3), p(1))), along)
        };
        let e3 = unit(normal).ok_or_else(degenerate)?;
        let e1 = unit(sub(along, e3.map(|c| c * dot(along, e3)))).ok_or_else(degenerate)?;
        let e2 = cross(e3, e1);
        let axes = [e1, e2, e3];
        let mut centre = [0.0; 3];
        for x in xyz {
            for k in 0..3 {
                centre[k] += x[k] / nodes as f64;
            }
        }
        let mut local = [[0.0; 2]; 4];
        let mut size = 0.0f64;
        for (l, x) in local.iter_mut().zip(xyz) {
            let d = sub(*x, centre);
            *l = [dot(d, e1), dot(d, e2)];
            size = size.max(norm(d));
        }
        // Corners turn left (counter-clockwise about e₃) by a margin: positive area, convex.
        for i in 0..nodes {
            let (a, b, c) = (local[(i + nodes - 1) % nodes], local[i], local[(i + 1) % nodes]);
            let turn = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
            if !(turn > 1e-10 * size * size) {
                return Err(degenerate());
            }
        }
        // A strongly warped quadrilateral is not a facet.
        let warp = xyz.iter().map(|x| dot(sub(*x, centre), e3).abs()).fold(0.0, f64::max);
        if warp > 0.1 * size {
            return Err(Error::Degenerate(format!(
                "a quad4 shell element is warped by {warp:.3e} out of its plane (more than a tenth of its size); split it into triangles"
            )));
        }
        let mut shell = Shell { kind, nodes, axes, local, section, area: 0.0 };
        shell.area = shell.rule().points().iter().map(|q| shell.map(q.xi[0], q.xi[1]).map_or(0.0, |m| m.det * q.weight)).sum();
        if !(shell.area > 0.0) {
            return Err(degenerate());
        }
        Ok(shell)
    }

    pub fn axes(&self) -> Axes {
        self.axes
    }

    pub fn area(&self) -> f64 {
        self.area
    }

    /// The element's unknowns: three translations and three rotations per node.
    pub fn dofs(&self) -> ElementDofs {
        ElementDofs::new(DofSet::STRUCTURAL)
    }

    /// The size of the element's matrices.
    fn size(&self) -> usize {
        NODE_DOFS * self.nodes
    }

    fn rule(&self) -> Rule {
        if self.kind == ElementKind::Tri3 { Rule::TriConical(2) } else { Rule::Quad(2) }
    }

    fn map(&self, r: f64, s: f64) -> Result<Map> {
        let (n, d) = shape(self.kind, r, s);
        let mut j = [[0.0; 2]; 2];
        for (di, x) in d.iter().zip(&self.local).take(self.nodes) {
            for a in 0..2 {
                for b in 0..2 {
                    j[a][b] += di[a] * x[b];
                }
            }
        }
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if !(det > 0.0 && det.is_finite()) {
            return Err(Error::Degenerate(format!("a {} shell element is inverted or degenerate", self.kind.name())));
        }
        let inv = [[j[1][1] / det, -j[0][1] / det], [-j[1][0] / det, j[0][0] / det]];
        let mut dx = [[0.0; 2]; 4];
        for (g, di) in dx.iter_mut().zip(&d) {
            g[0] = inv[0][0] * di[0] + inv[0][1] * di[1];
            g[1] = inv[1][0] * di[0] + inv[1][1] * di[1];
        }
        Ok(Map { n, dn: d, dx, j, inv, det })
    }

    /// The covariant transverse shear strain e_rt (`along` 0) or e_st (`along` 1) of the
    /// displacement interpolation at (r, s): w,r + θ_y x,r − θ_x y,r (and with s).
    fn covariant(&self, r: f64, s: f64, along: usize) -> Result<Row> {
        let m = self.map(r, s)?;
        let (xr, yr) = (m.j[along][0], m.j[along][1]);
        let mut row = [0.0; MAX_DOFS];
        for i in 0..self.nodes {
            row[6 * i + 2] = m.dn[i][along];
            row[6 * i + 4] = m.n[i] * xr;
            row[6 * i + 3] = -m.n[i] * yr;
        }
        Ok(row)
    }

    fn tying(&self) -> Result<Tying> {
        if self.kind == ElementKind::Tri3 {
            let e1 = self.covariant(0.5, 0.0, 0)?;
            let e2 = self.covariant(0.0, 0.5, 1)?;
            let e3r = self.covariant(0.5, 0.5, 0)?;
            let e3s = self.covariant(0.5, 0.5, 1)?;
            let c: Row = std::array::from_fn(|k| e3r[k] - e1[k] - e3s[k] + e2[k]);
            Ok(Tying([e1, e2, c, [0.0; MAX_DOFS]]))
        } else {
            Ok(Tying([self.covariant(0.0, 1.0, 0)?, self.covariant(0.0, -1.0, 0)?, self.covariant(1.0, 0.0, 1)?, self.covariant(-1.0, 0.0, 1)?]))
        }
    }

    /// The assumed covariant shear strains (ẽ_rt, ẽ_st) at (r, s).
    fn assumed_shear(&self, tying: &Tying, r: f64, s: f64) -> [Row; 2] {
        let lerp = |a: &Row, wa: f64, b: &Row, wb: f64| -> Row { std::array::from_fn(|k| wa * a[k] + wb * b[k]) };
        let [t0, t1, t2, t3] = &tying.0;
        if self.kind == ElementKind::Tri3 {
            [lerp(t0, 1.0, t2, s), lerp(t1, 1.0, t2, -r)]
        } else {
            [lerp(t0, 0.5 * (1.0 + s), t1, 0.5 * (1.0 - s)), lerp(t2, 0.5 * (1.0 + r), t3, 0.5 * (1.0 - r))]
        }
    }

    /// The generalized strain operator at (r, s): membrane strains, curvatures, transverse shear
    /// strains (8 rows) and the drilling row θ_z − ω, with the map.
    fn operator(&self, tying: &Tying, r: f64, s: f64) -> Result<([Row; 9], Map)> {
        let m = self.map(r, s)?;
        let mut b = [[0.0; MAX_DOFS]; 9];
        for i in 0..self.nodes {
            let [gx, gy] = m.dx[i];
            let (u, v, tx, ty, tz) = (6 * i, 6 * i + 1, 6 * i + 3, 6 * i + 4, 6 * i + 5);
            b[0][u] = gx;
            b[1][v] = gy;
            b[2][u] = gy;
            b[2][v] = gx;
            b[3][ty] = gx;
            b[4][tx] = -gy;
            b[5][ty] = gy;
            b[5][tx] = -gx;
            b[8][tz] = m.n[i];
            b[8][u] = 0.5 * gy;
            b[8][v] = -0.5 * gx;
        }
        let [ert, est] = self.assumed_shear(tying, r, s);
        for k in 0..MAX_DOFS {
            b[6][k] = m.inv[0][0] * ert[k] + m.inv[0][1] * est[k];
            b[7][k] = m.inv[1][0] * ert[k] + m.inv[1][1] * est[k];
        }
        Ok((b, m))
    }

    /// The stiffness matrix in global axes, 6n × 6n, degrees of freedom node by node.
    pub fn stiffness(&self) -> Result<Matrix> {
        let n = self.size();
        let abd = self.section.abd();
        let ks = self.section.shear();
        let drill = DRILLING * self.section.shear_modulus() * self.section.thickness;
        let tying = self.tying()?;
        let mut k = Matrix::zeros(n);
        for p in self.rule().points() {
            let (b, m) = self.operator(&tying, p.xi[0], p.xi[1])?;
            let w = m.det * p.weight;
            // C B for the membrane–bending rows.
            let cb: [Row; 6] = std::array::from_fn(|a| std::array::from_fn(|c| (0..6).map(|d| abd[a][d] * b[d][c]).sum()));
            for i in 0..n {
                for j in i..n {
                    let mut v: f64 = (0..6).map(|a| b[a][i] * cb[a][j]).sum();
                    v += ks * (b[6][i] * b[6][j] + b[7][i] * b[7][j]);
                    v += drill * b[8][i] * b[8][j];
                    if v != 0.0 {
                        k.add(i, j, w * v);
                        if i != j {
                            k.add(j, i, w * v);
                        }
                    }
                }
            }
        }
        Ok(rotate_matrix(&k, &self.axes))
    }

    /// The consistent mass matrix in global axes, 6n × 6n.
    pub fn mass(&self, density: f64) -> Result<Matrix> {
        Ok(rotate_matrix(&self.local_mass(density)?, &self.axes))
    }

    /// The consistent mass matrix in the element's axes.
    fn local_mass(&self, density: f64) -> Result<Matrix> {
        let (t, e) = (self.section.thickness, self.section.eccentricity);
        let (m0, m1, m2) = (density * t, density * t * e, density * (t * t * t / 12.0 + t * e * e));
        let mut mm = Matrix::zeros(self.size());
        for p in self.rule().points() {
            let m = self.map(p.xi[0], p.xi[1])?;
            let w = m.det * p.weight;
            for i in 0..self.nodes {
                for j in 0..self.nodes {
                    let nn = w * m.n[i] * m.n[j];
                    for c in 0..3 {
                        mm.add(6 * i + c, 6 * j + c, m0 * nn);
                    }
                    for c in 3..6 {
                        mm.add(6 * i + c, 6 * j + c, m2 * nn);
                    }
                    if m1 != 0.0 {
                        // u θ_y and −v θ_x, both ways.
                        mm.add(6 * i, 6 * j + 4, m1 * nn);
                        mm.add(6 * i + 4, 6 * j, m1 * nn);
                        mm.add(6 * i + 1, 6 * j + 3, -m1 * nn);
                        mm.add(6 * i + 3, 6 * j + 1, -m1 * nn);
                    }
                }
            }
        }
        Ok(mm)
    }

    /// The lumped mass per degree of freedom (the same in local and global axes).
    pub fn lumped_mass(&self, density: f64) -> Result<Vec<f64>> {
        let (t, e) = (self.section.thickness, self.section.eccentricity);
        let mut diag = vec![0.0; self.size()];
        let mut sum = 0.0;
        for p in self.rule().points() {
            let m = self.map(p.xi[0], p.xi[1])?;
            for i in 0..self.nodes {
                let v = m.det * p.weight * m.n[i] * m.n[i];
                diag[6 * i] += v;
                sum += v;
            }
        }
        if !(sum > 0.0) {
            return Err(Error::Degenerate("a shell element has no area".into()));
        }
        let scale = self.area / sum;
        let (m0, m2) = (density * t, density * (t * t * t / 12.0 + t * e * e));
        for i in 0..self.nodes {
            let a = diag[6 * i] * scale;
            for c in 0..3 {
                diag[6 * i + c] = m0 * a;
                diag[6 * i + 3 + c] = m2 * a;
            }
        }
        Ok(diag)
    }

    /// Consistent nodal forces and moments (global) of a uniform force per volume b: the
    /// consistent mass of unit density times the uniform nodal translation b.
    pub fn body_force(&self, b: [f64; 3]) -> Result<Vec<f64>> {
        let r = translation(to_local(&self.axes, b), self.nodes);
        Ok(rotate_to_global(&self.local_mass(1.0)?.mul_vec(&r), &self.axes))
    }

    /// Section forces and surface stresses for the element's nodal displacements and rotations
    /// (global, 6 per node).
    pub fn stresses(&self, u: &[f64]) -> Result<ShellStress> {
        if u.len() != self.size() {
            return Err(Error::Invalid(format!("{} values for a shell element with {} degrees of freedom", u.len(), self.size())));
        }
        let ul = rotate_to_local(u, &self.axes);
        let tying = self.tying()?;
        let strains = |r: f64, s: f64| -> Result<[f64; 8]> {
            let (b, _) = self.operator(&tying, r, s)?;
            Ok(std::array::from_fn(|a| b[a].iter().zip(&ul).map(|(x, y)| x * y).sum()))
        };
        let nodal: Vec<[f64; 8]> = if self.kind == ElementKind::Tri3 {
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]].iter().map(|c| strains(c[0], c[1])).collect::<Result<_>>()?
        } else {
            // 2 × 2 Gauss values, extrapolated with the bilinear functions of the points.
            let g = 1.0 / 3f64.sqrt();
            let corners = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
            let at: Vec<[f64; 8]> = corners.iter().map(|c| strains(c[0] * g, c[1] * g)).collect::<Result<_>>()?;
            corners
                .iter()
                .map(|c| {
                    let (n, _) = shape(ElementKind::Quad4, c[0] / g, c[1] / g);
                    std::array::from_fn(|k| n.iter().zip(&at).map(|(w, v)| w * v[k]).sum())
                })
                .collect()
        };
        let abd = self.section.abd();
        let q = self.section.plane_stress();
        let ks = self.section.shear();
        let (t, e) = (self.section.thickness, self.section.eccentricity);
        let mut out = ShellStress { axes: self.axes, forces: Vec::new(), bottom: Vec::new(), middle: Vec::new(), top: Vec::new() };
        for g in &nodal {
            let mut f = [0.0; 8];
            for a in 0..6 {
                f[a] = (0..6).map(|b| abd[a][b] * g[b]).sum();
            }
            f[6] = ks * g[6];
            f[7] = ks * g[7];
            out.forces.push(f);
            let at = |z: f64, shear: f64| {
                let strain = [g[0] + z * g[3], g[1] + z * g[4], g[2] + z * g[5]];
                let s: [f64; 3] = std::array::from_fn(|a| (0..3).map(|b| q[a][b] * strain[b]).sum());
                let (txz, tyz) = (shear * f[6] / t, shear * f[7] / t);
                tensor_to_global(&self.axes, [[s[0], s[2], txz], [s[2], s[1], tyz], [txz, tyz, 0.0]])
            };
            out.bottom.push(at(e - t / 2.0, 0.0));
            out.middle.push(at(e, 1.5));
            out.top.push(at(e + t / 2.0, 0.0));
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const STEEL: Elastic = Elastic::Isotropic { young: 210e3, poisson: 0.3 };

    fn rigid(xyz: &[[f64; 3]], m: usize) -> Vec<f64> {
        xyz.iter()
            .flat_map(|x| {
                let mut d = [0.0; 6];
                if m < 3 {
                    d[m] = 1.0;
                } else {
                    let mut w = [0.0; 3];
                    w[m - 3] = 1.0;
                    d[..3].copy_from_slice(&cross(w, *x));
                    d[3..].copy_from_slice(&w);
                }
                d
            })
            .collect()
    }

    #[test]
    fn shells_are_free_of_rigid_motion_and_symmetric() {
        let section = ShellSection::new(0.5, 0.25, &STEEL).unwrap();
        // A skew quadrilateral and a triangle in a tilted plane.
        let tilt = |x: f64, y: f64| [x, 0.6 * y, 0.8 * y + 1.0];
        let quad = [tilt(0.0, 0.0), tilt(2.0, 0.2), tilt(2.3, 1.7), tilt(-0.2, 1.4)];
        let tri = [tilt(0.0, 0.0), tilt(2.0, 0.2), tilt(0.4, 1.5)];
        for (kind, xyz) in [(ElementKind::Quad4, &quad[..]), (ElementKind::Tri3, &tri[..])] {
            let shell = Shell::new(kind, xyz, section).unwrap();
            let k = shell.stiffness().unwrap();
            let scale = (0..k.n).map(|i| k.get(i, i)).fold(0.0, f64::max);
            for i in 0..k.n {
                for j in 0..k.n {
                    assert_eq!(k.get(i, j), k.get(j, i));
                }
            }
            for m in 0..6 {
                let f = k.mul_vec(&rigid(xyz, m));
                assert!(f.iter().all(|v| v.abs() < 1e-10 * scale), "{kind:?} rigid {m}: {f:?}");
            }
            // Mass: translations sum to ρ t A per direction (the offset couples rotations only).
            let m = shell.mass(2.0).unwrap();
            let total = 2.0 * 0.5 * shell.area();
            let s: f64 =
                (0..k.n).filter(|i| i % 6 == 0).flat_map(|i| (0..k.n).filter(|j| j % 6 == 0).map(move |j| (i, j))).map(|(i, j)| m.get(i, j)).sum();
            assert!((s - total).abs() < 1e-12 * total, "{s} {total}");
            let lumped = shell.lumped_mass(2.0).unwrap();
            let s: f64 = lumped.iter().step_by(6).sum();
            assert!((s - total).abs() < 1e-12 * total);
            let body = shell.body_force([0.0, 0.0, -1.0]).unwrap();
            let fz: f64 = body.iter().skip(2).step_by(6).sum();
            assert!((fz + 0.5 * shell.area()).abs() < 1e-12, "{fz}");
        }
        let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.5, 0.2, 0.0]];
        assert!(Shell::new(ElementKind::Quad4, &flat, section).is_err(), "not convex");
        let line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(Shell::new(ElementKind::Tri3, &line, section).is_err());
        assert!(Shell::new(ElementKind::Quad8, &flat, section).is_err());
        assert!(ShellSection::new(1.0, 0.0, &Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0; 3] }).is_err());
    }
}
