//! The two-node Timoshenko beam (six degrees of freedom per node), with the section properties of
//! the model's beam profiles. (Axial-only members are the trusses of [`crate::truss`].)
//!
//! **Local axes.** x runs from node 1 to node 2; y is the section's orientation vector (the
//! profile's local 1 direction, its `width`) made normal to x; z = x × y (the profile's local 2
//! direction, its `height`). Degrees of freedom per node: translations along x, y, z, then
//! rotations about them; element matrices are built in local axes and rotated to the global ones
//! ([`crate::frame`]).
//!
//! **Beam.** Axial and torsional stiffness are the exact two-node ones (EA/L, GJ/L). Each section
//! names its bending theory ([`BeamTheory`]). Timoshenko (the default): each bending plane uses
//! the interdependent interpolation of the Timoshenko beam: cubic deflection and quadratic
//! rotation, coupled through Φ = 12 EI / (G A_s L²) so the shear strain is constant and the
//! element is exact for end loads (any number of elements gives the nodal displacements of
//! Timoshenko beam theory; J. S. Przemieniecki, *Theory of Matrix Structural Analysis*,
//! McGraw-Hill 1968, §5.6; Z. Friedman, J. B. Kosmatka, *An improved two-node Timoshenko beam
//! finite element*, Computers & Structures 47(3):473–481, 1993, whose shape functions are used
//! here). Euler–Bernoulli: Φ = 0, the cubic Hermite element of the classical theory, with no
//! shear strain energy and no rotary inertia of the cross-section in bending (Przemieniecki
//! 1968, §5.6; R. D. Cook, D. S. Malkus, M. E. Plesha, R. J. Witt, *Concepts and Applications of
//! Finite Element Analysis*, 4th ed., Wiley 2002, §2.3). With these functions w(ξ) and θ(ξ),
//! ξ = x/L, K = ∫ (EI θ'ᵀ θ' + G A_s (w' − θ)ᵀ (w' − θ)) dx and the consistent mass
//! M = ∫ (ρA N_wᵀ N_w + ρI N_θᵀ N_θ) dx (translational and rotary inertia) are integrated with
//! four Gauss points, exactly; Φ → 0 recovers the Euler–Bernoulli (Hermite) element. The
//! rotation about y is −dw/dx in the x–z plane, so that plane's rotations enter with a minus
//! sign. Lumped mass: half the mass on each node's translations, half the polar inertia on the
//! twist, and the consistent bending rotary diagonal scaled like the translations (Hinton, Rock &
//! Zienkiewicz 1976). End forces are K u minus the consistent loads of the element's body force,
//! reported as the section forces at each end (axial force N, shears V_y, V_z, torque T, moments
//! M_y, M_z on the cut whose outward normal is +x).
//!
//! **Pipes under pressure.** A beam whose profile is a pipe with a bore takes a pressure inside
//! and around it ([`crate::pipe`]): the axial pair (1 − 2ν) A A_w of the closed-end force and
//! the wall's Poisson strain as loads, the wall force N (the closed-end force is an external
//! load, the Poisson strain an initial one) and Lamé's radial and hoop stresses as results.
//!
//! **Section properties** are the profile's ([`fieldcraft_model::BeamProperties`], computed in
//! the model so the deck writers give other solvers the same numbers): area, second moments, the
//! torsion constant and the shear areas; the polar moment I_y + I_z is the rotary inertia of
//! torsion.
//!
//! **2D models.** A beam in the XY plane of a 2D model keeps the in-plane block of these
//! matrices, Ux, Uy and Rz at each node: the plane frame element (Przemieniecki 1968, §5.6, with
//! the shear of the interdependent interpolation). The section's orientation is then along Z, so
//! Z is a principal axis of the profile and the block is uncoupled from Uz, Rx and Ry. The
//! solver takes the block (`fieldcraft_solver::system::structural_kinds`).

use fieldcraft_model::{BeamProfile, BeamProperties, BeamTheory, Elastic};

use crate::dense::Matrix;
use crate::dof::{DofSet, ElementDofs};
use crate::frame::{Axes, Vec3, cross, dot, norm, rotate_matrix, rotate_to_global, rotate_to_local, sub, to_local, translation, unit};
use crate::pipe::{PipePressure, PipeStress, PipeWall};
use crate::quadrature::gauss_legendre;
use crate::{Error, Result};

/// The Young's and shear moduli and Poisson's ratio of an isotropic material.
/// (`what` names the section kind in errors.)
fn moduli(elastic: &Elastic, what: &str) -> Result<(f64, f64, f64)> {
    let (young, poisson) = crate::material::isotropic_constants(elastic, what)?;
    Ok((young, young / (2.0 * (1.0 + poisson)), poisson))
}

/// A beam section: profile properties and an isotropic material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamSection {
    pub properties: BeamProperties,
    pub young: f64,
    pub shear_modulus: f64,
    pub theory: BeamTheory,
    pub poisson: f64,
    /// The wall of a pipe profile with a bore: what a pipe pressure acts on.
    pub wall: Option<PipeWall>,
}

impl BeamSection {
    pub fn new(profile: &BeamProfile, elastic: &Elastic, theory: BeamTheory) -> Result<BeamSection> {
        profile.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        let (young, shear_modulus, poisson) = moduli(elastic, "beam")?;
        Ok(BeamSection { properties: profile.properties(poisson), young, shear_modulus, theory, poisson, wall: PipeWall::of(profile) })
    }
}

/// An axial force pair in local axes: `force` along +x at node 2, against it at node 1.
fn axial_pair(force: f64) -> [f64; 12] {
    let mut f = [0.0; 12];
    f[0] = -force;
    f[6] = force;
    f
}

/// The section forces at the two ends of a beam: N, V_y, V_z, T, M_y, M_z in its local axes.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct BeamForces {
    pub axes: Axes,
    pub ends: [[f64; 6]; 2],
}

/// One bending plane's stiffness and mass (deflection, rotation at node 1, then node 2), with
/// θ = dw/dx in the plane's own convention.
fn bending_plane(ei: f64, shear: Option<f64>, rho_a: f64, rho_i: f64, l: f64) -> ([[f64; 4]; 4], [[f64; 4]; 4]) {
    let (phi, gas) = match shear {
        Some(gas) => (12.0 * ei / (gas * l * l), gas),
        None => (0.0, 0.0),
    };
    let mu = 1.0 / (1.0 + phi);
    let (mut k, mut m) = ([[0.0; 4]; 4], [[0.0; 4]; 4]);
    for &(x, w) in gauss_legendre(4) {
        let xi = 0.5 * (x + 1.0);
        let w = 0.5 * w * l;
        let x2 = xi * xi;
        let x3 = x2 * xi;
        let nw = [
            mu * (2.0 * x3 - 3.0 * x2 - phi * xi + 1.0 + phi),
            mu * l * (x3 - (2.0 + phi / 2.0) * x2 + (1.0 + phi / 2.0) * xi),
            -mu * (2.0 * x3 - 3.0 * x2 - phi * xi),
            mu * l * (x3 - (1.0 - phi / 2.0) * x2 - phi / 2.0 * xi),
        ];
        let nt =
            [mu * 6.0 / l * (x2 - xi), mu * (3.0 * x2 - (4.0 + phi) * xi + 1.0 + phi), -mu * 6.0 / l * (x2 - xi), mu * (3.0 * x2 - (2.0 - phi) * xi)];
        // θ' and γ = w' − θ.
        let dt = [
            mu * 6.0 / (l * l) * (2.0 * xi - 1.0),
            mu / l * (6.0 * xi - 4.0 - phi),
            -mu * 6.0 / (l * l) * (2.0 * xi - 1.0),
            mu / l * (6.0 * xi - 2.0 + phi),
        ];
        let dw = [
            mu / l * (6.0 * x2 - 6.0 * xi - phi),
            mu * (3.0 * x2 - (4.0 + phi) * xi + 1.0 + phi / 2.0),
            -mu / l * (6.0 * x2 - 6.0 * xi - phi),
            mu * (3.0 * x2 - (2.0 - phi) * xi - phi / 2.0),
        ];
        let g: [f64; 4] = std::array::from_fn(|a| dw[a] - nt[a]);
        for a in 0..4 {
            for b in 0..4 {
                k[a][b] += w * (ei * dt[a] * dt[b] + gas * g[a] * g[b]);
                m[a][b] += w * (rho_a * nw[a] * nw[b] + rho_i * nt[a] * nt[b]);
            }
        }
    }
    (k, m)
}

/// A two-node beam element.
#[derive(Debug, Clone, Copy)]
pub struct Beam {
    axes: Axes,
    length: f64,
    section: BeamSection,
}

/// A line element's unit axis (node 1 to node 2) and length.
fn line_axis(xyz: &[[f64; 3]]) -> Result<(Vec3, f64)> {
    let [a, b] = xyz else {
        return Err(Error::Invalid(format!("a two-node line element has 2 nodes, got {}", xyz.len())));
    };
    fieldcraft_model::check_coordinates(xyz.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
    let d = sub(*b, *a);
    let length = norm(d);
    let ex = unit(d).ok_or_else(|| Error::Degenerate("a line element has zero length".into()))?;
    Ok((ex, length))
}

/// An orientation vector within this sine of the beam's axis is parallel to it.
const PARALLEL: f64 = 1e-6; // unit-factors: the sine of an angle

/// The beam's local axes from its end points and orientation vector.
fn beam_axes(xyz: &[[f64; 3]], orientation: Vec3) -> Result<(Axes, f64)> {
    let (ex, length) = line_axis(xyz)?;
    let o = unit(orientation).ok_or_else(|| Error::Invalid("a beam's orientation vector must be nonzero and finite".into()))?;
    let proj = sub(o, ex.map(|c| c * dot(o, ex)));
    let ey = unit(proj)
        .filter(|_| norm(proj) > PARALLEL)
        .ok_or_else(|| Error::Invalid("the beam section's orientation vector is parallel to the beam".into()))?;
    Ok(([ex, ey, cross(ex, ey)], length))
}

impl Beam {
    pub fn new(xyz: &[[f64; 3]], orientation: Vec3, section: BeamSection) -> Result<Beam> {
        let (axes, length) = beam_axes(xyz, orientation)?;
        Ok(Beam { axes, length, section })
    }

    pub fn axes(&self) -> Axes {
        self.axes
    }

    pub fn length(&self) -> f64 {
        self.length
    }

    /// The element's unknowns: three translations and three rotations per node.
    pub fn dofs(&self) -> ElementDofs {
        ElementDofs::new(DofSet::STRUCTURAL)
    }

    /// Local stiffness and consistent mass, 12 × 12.
    fn local(&self, density: f64) -> (Matrix, Matrix) {
        let (p, e, g, l) = (self.section.properties, self.section.young, self.section.shear_modulus, self.length);
        let (mut k, mut m) = (Matrix::zeros(12), Matrix::zeros(12));
        let pair = |k: &mut Matrix, m: &mut Matrix, dofs: [usize; 2], stiff: f64, mass: f64| {
            for (a, &i) in dofs.iter().enumerate() {
                for (b, &j) in dofs.iter().enumerate() {
                    k.add(i, j, if a == b { stiff } else { -stiff });
                    m.add(i, j, mass / 6.0 * if a == b { 2.0 } else { 1.0 });
                }
            }
        };
        pair(&mut k, &mut m, [0, 6], e * p.area / l, density * p.area * l);
        pair(&mut k, &mut m, [3, 9], g * p.torsion / l, density * (p.iy + p.iz) * l);
        // x–y plane: v and θ_z = dv/dx; x–z plane: w and θ_y = −dw/dx.
        let planes = [([1, 5, 7, 11], [1.0; 4], p.iz, p.shear_y), ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0], p.iy, p.shear_z)];
        for (dofs, sign, i, shear) in planes {
            let (shear, rotary) = match self.section.theory {
                BeamTheory::Timoshenko => (Some(g * shear), density * i),
                BeamTheory::EulerBernoulli => (None, 0.0),
            };
            let (kp, mp) = bending_plane(e * i, shear, density * p.area, rotary, l);
            for a in 0..4 {
                for b in 0..4 {
                    k.add(dofs[a], dofs[b], sign[a] * sign[b] * kp[a][b]);
                    m.add(dofs[a], dofs[b], sign[a] * sign[b] * mp[a][b]);
                }
            }
        }
        (k, m)
    }

    pub fn stiffness(&self) -> Matrix {
        rotate_matrix(&self.local(0.0).0, &self.axes)
    }

    pub fn mass(&self, density: f64) -> Matrix {
        rotate_matrix(&self.local(density).1, &self.axes)
    }

    /// The lumped mass, block diagonal (3 × 3 rotation blocks) in global axes.
    pub fn lumped_mass(&self, density: f64) -> Matrix {
        let (_, m) = self.local(density);
        let (p, l) = (self.section.properties, self.length);
        let half = density * p.area * l / 2.0;
        let mut d = Matrix::zeros(12);
        // Bending rotations scaled like their plane's translations (equal in both planes).
        for (rotation, deflection) in [(5, 1), (4, 2)] {
            let scale = 2.0 * half / (m.get(deflection, deflection) + m.get(deflection + 6, deflection + 6));
            d.add(rotation, rotation, scale * m.get(rotation, rotation));
            d.add(rotation + 6, rotation + 6, scale * m.get(rotation + 6, rotation + 6));
        }
        for node in [0, 6] {
            for c in 0..3 {
                d.add(node + c, node + c, half);
            }
            d.add(node + 3, node + 3, density * (p.iy + p.iz) * l / 2.0);
        }
        rotate_matrix(&d, &self.axes)
    }

    /// Consistent nodal loads of a uniform force per volume b, in local axes: the consistent
    /// mass of unit density times the uniform nodal translation.
    fn local_body_force(&self, b: [f64; 3]) -> Vec<f64> {
        self.local(1.0).1.mul_vec(&translation(to_local(&self.axes, b), 2))
    }

    /// Consistent nodal forces and moments (global) of a uniform force per volume b.
    pub fn body_force(&self, b: [f64; 3]) -> Vec<f64> {
        rotate_to_global(&self.local_body_force(b), &self.axes)
    }

    /// The axial force pairs of a pipe pressure ([`crate::pipe`]): the closed-end force
    /// F = A A_w, and the force E A_w ε₀ = −2νA A_w of the wall's Poisson strain.
    fn pressure_pairs(&self, pressure: PipePressure) -> Result<[f64; 2]> {
        let end = self.wall()?.mean_transverse(pressure)? * self.section.properties.area;
        let pairs = [end, -2.0 * self.section.poisson * end];
        if pairs.iter().all(|f| f.is_finite()) { Ok(pairs) } else { Err(Error::Invalid("the loads of a pipe pressure are not finite".into())) }
    }

    fn wall(&self) -> Result<PipeWall> {
        self.section.wall.ok_or_else(|| Error::Invalid("a pipe pressure acts on a beam with a pipe profile that has a bore".into()))
    }

    /// Nodal forces (global, 12) of a pressure inside and around a pipe: the axial pair
    /// (1 − 2ν) A A_w of the closed-end force and the wall's Poisson strain.
    pub fn pressure_load(&self, pressure: PipePressure) -> Result<Vec<f64>> {
        let [end, poisson] = self.pressure_pairs(pressure)?;
        Ok(rotate_to_global(&axial_pair(end + poisson), &self.axes))
    }

    /// The wall stresses of a pipe under `pressure` whose end forces are `forces`.
    pub fn pipe_stress(&self, forces: &BeamForces, pressure: PipePressure) -> Result<PipeStress> {
        let wall = self.wall()?;
        let [inner, outer] = wall.radii();
        let (bore, outside) = (wall.stresses(pressure, inner)?, wall.stresses(pressure, outer)?);
        let area = self.section.properties.area;
        Ok(PipeStress { radial: [bore[0], outside[0]], hoop: [bore[1], outside[1]], axial: forces.ends.map(|end| end[0] / area) })
    }

    /// Section forces at the ends for the nodal displacements and rotations `u` (global, 12)
    /// under a uniform force per volume `body` acting on the element and, in a pipe, the
    /// pressure inside and around it: N is then the force in the wall.
    pub fn forces(&self, u: &[f64], body: [f64; 3], pressure: Option<PipePressure>) -> Result<BeamForces> {
        if u.len() != 12 {
            return Err(Error::Invalid(format!("{} values for a beam element with 12 degrees of freedom", u.len())));
        }
        // The wall's Poisson strain is an initial strain of the element; the closed-end force
        // is an external load, like a force at the nodes.
        let initial = match pressure {
            Some(p) => axial_pair(self.pressure_pairs(p)?[1]),
            None => [0.0; 12],
        };
        let ku = self.local(0.0).0.mul_vec(&rotate_to_local(u, &self.axes));
        let end: Vec<f64> = ku.iter().zip(self.local_body_force(body)).zip(initial).map(|((a, b), c)| a - b - c).collect();
        let at = |i: usize| end.get(i).copied().unwrap_or(0.0);
        Ok(BeamForces { axes: self.axes, ends: [std::array::from_fn(|c| -at(c)), std::array::from_fn(|c| at(6 + c))] })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const STEEL: Elastic = Elastic::Isotropic { young: 210e3, poisson: 0.3 };

    /// Rigid translations and rotations about a point give no force.
    fn rigid_motions(xyz: &[[f64; 3]]) -> Vec<Vec<f64>> {
        (0..6)
            .map(|m| {
                xyz.iter()
                    .flat_map(|x| {
                        let mut d = [0.0; 6];
                        if m < 3 {
                            d[m] = 1.0;
                        } else {
                            let mut w = [0.0; 3];
                            w[m - 3] = 1.0;
                            let t = cross(w, *x);
                            d[..3].copy_from_slice(&t);
                            d[3..].copy_from_slice(&w);
                        }
                        d
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn one_beam_element_is_exact_for_a_tip_load() {
        // Cantilever along a skew axis; tip force along the local z axis: δ = PL³/3EI + PL/(G A_s).
        let profile = BeamProfile::Rectangle { width: 3.0, height: 5.0 };
        let section = BeamSection::new(&profile, &STEEL, BeamTheory::Timoshenko).unwrap();
        let xyz = [[1.0, 2.0, 3.0], [41.0, 32.0, 3.0]];
        let beam = Beam::new(&xyz, [0.0, 0.0, 1.0], section).unwrap();
        let l = beam.length();
        assert!((l - 50.0).abs() < 1e-12);
        let k = beam.stiffness();
        for r in rigid_motions(&xyz) {
            let f = k.mul_vec(&r);
            assert!(f.iter().all(|v| v.abs() < 1e-6), "{f:?}");
        }
        // Clamp node 1: solve the 6 × 6 block of node 2 for a unit force along local z.
        let axes = beam.axes();
        let mut k22 = Matrix::zeros(6);
        for i in 0..6 {
            for j in 0..6 {
                k22.add(i, j, k.get(6 + i, 6 + j));
            }
        }
        let mut f = vec![0.0; 6];
        f[..3].copy_from_slice(&axes[2]);
        crate::dense::Cholesky::new(&k22).unwrap().solve(&mut f);
        let p = section.properties;
        let want = l.powi(3) / (3.0 * section.young * p.iy) + l / (section.shear_modulus * p.shear_z);
        let got = dot([f[0], f[1], f[2]], axes[2]);
        assert!((got - want).abs() < 1e-12 * want, "{got} vs {want}");
        // Section forces at the clamp: V_z = 1, M_y = −L (the moment about y of a +z load at +x).
        let mut u = vec![0.0; 12];
        u[6..].copy_from_slice(&f);
        let forces = beam.forces(&u, [0.0; 3], None).unwrap();
        assert!((forces.ends[0][2] - 1.0).abs() < 1e-9 && (forces.ends[0][4] + l).abs() < 1e-8, "{forces:?}");
        assert!((forces.ends[1][2] - 1.0).abs() < 1e-9 && forces.ends[1][4].abs() < 1e-8, "{forces:?}");
        // Mass: the translational entries sum to the mass along each direction.
        let m = beam.mass(7.8e-9);
        let total = 7.8e-9 * p.area * l;
        for d in 0..3 {
            let s: f64 = [d, 6 + d].iter().flat_map(|&i| [d, 6 + d].map(|j| m.get(i, j))).sum();
            assert!((s - total).abs() < 1e-12 * total, "{s} {total}");
        }
        let lumped = beam.lumped_mass(7.8e-9);
        assert!((lumped.get(0, 0) - total / 2.0).abs() < 1e-12 * total);
        assert!(Beam::new(&xyz, [40.0, 30.0, 0.0], section).is_err(), "orientation along the axis");
        assert!(Beam::new(&[xyz[0], xyz[0]], [0.0, 0.0, 1.0], section).is_err(), "zero length");
    }

    #[test]
    fn shear_rigid_bending_plane_is_the_hermite_element() {
        let (ei, rho_a, l): (f64, f64, f64) = (3.0e5, 2.5e-3, 7.0);
        let (k, m) = bending_plane(ei, None, rho_a, 0.0, l);
        let kc = ei / l.powi(3);
        let want_k = [
            [12.0, 6.0 * l, -12.0, 6.0 * l],
            [6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l],
            [-12.0, -6.0 * l, 12.0, -6.0 * l],
            [6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l],
        ];
        let mc = rho_a * l / 420.0;
        let want_m = [
            [156.0, 22.0 * l, 54.0, -13.0 * l],
            [22.0 * l, 4.0 * l * l, 13.0 * l, -3.0 * l * l],
            [54.0, 13.0 * l, 156.0, -22.0 * l],
            [-13.0 * l, -3.0 * l * l, -22.0 * l, 4.0 * l * l],
        ];
        let (mut ek, mut em) = (0.0f64, 0.0f64);
        for a in 0..4 {
            for b in 0..4 {
                ek = ek.max((k[a][b] - kc * want_k[a][b]).abs() / (kc * want_k[a][b]).abs());
                em = em.max((m[a][b] - mc * want_m[a][b]).abs() / (mc * want_m[a][b]).abs());
            }
        }
        assert!(ek < 1e-12 && em < 1e-12, "{ek} {em}");
        // One element, clamped: the tip deflection under a tip force is L³/(3 E I), no shear term.
        let profile = BeamProfile::Rectangle { width: 3.0, height: 5.0 };
        let section = BeamSection::new(&profile, &STEEL, BeamTheory::EulerBernoulli).unwrap();
        let xyz = [[1.0, 2.0, 3.0], [41.0, 32.0, 3.0]];
        let beam = Beam::new(&xyz, [0.0, 0.0, 1.0], section).unwrap();
        let (kg, axes, len) = (beam.stiffness(), beam.axes(), beam.length());
        let mut k22 = Matrix::zeros(6);
        for i in 0..6 {
            for j in 0..6 {
                k22.add(i, j, kg.get(6 + i, 6 + j));
            }
        }
        let mut f = vec![0.0; 6];
        f[..3].copy_from_slice(&axes[2]);
        crate::dense::Cholesky::new(&k22).unwrap().solve(&mut f);
        let want = len.powi(3) / (3.0 * section.young * section.properties.iy);
        let got = dot([f[0], f[1], f[2]], axes[2]);
        assert!((got - want).abs() < 1e-12 * want, "{got} vs {want}");
        // Lumped rotary inertia of bending: ρ A L³ / 78 per node (the consistent diagonal scaled like the translations).
        let rho = 7.8e-9;
        let lumped = Beam::new(&[[0.0; 3], [len, 0.0, 0.0]], [0.0, 1.0, 0.0], section).unwrap().lumped_mass(rho);
        let want = rho * section.properties.area * len.powi(3) / 78.0;
        assert!((lumped.get(4, 4) - want).abs() < 1e-12 * want && (lumped.get(5, 5) - want).abs() < 1e-12 * want);
    }

    #[test]
    fn a_closed_pipe_under_pressure_stretches_and_carries_the_end_force() {
        // a = 10, b = 20, p_i = 100, p_o = 40: Lamé's A = −20 (Timoshenko & Goodier 1970, §28).
        let profile = BeamProfile::Pipe { radius: 20.0, thickness: 10.0 };
        let section = BeamSection::new(&profile, &STEEL, BeamTheory::Timoshenko).unwrap();
        let xyz = [[1.0, 2.0, 3.0], [41.0, 32.0, 3.0]];
        let beam = Beam::new(&xyz, [0.0, 0.0, 1.0], section).unwrap();
        let pressure = PipePressure { internal: 100.0, external: 40.0 };
        let (mean, area, nu, l) = (-20.0, section.properties.area, 0.3, beam.length());
        let scale = (mean * area).abs();
        // The load: the pair (1 − 2ν) A A_w along the axis, no moments.
        let f = beam.pressure_load(pressure).unwrap();
        let axis = beam.axes()[0];
        let pair = (1.0 - 2.0 * nu) * mean * area;
        for c in 0..3 {
            assert!((f[c] + pair * axis[c]).abs() < 1e-12 * scale && (f[6 + c] - pair * axis[c]).abs() < 1e-12 * scale, "{f:?}");
            assert!(f[3 + c] == 0.0 && f[9 + c] == 0.0, "{f:?}");
        }
        // Clamped at node 1 and free to stretch: ε_x = (1 − 2ν) A/E, N = A A_w, σ_x = A.
        let k = beam.stiffness();
        let mut k22 = Matrix::zeros(6);
        for i in 0..6 {
            for j in 0..6 {
                k22.add(i, j, k.get(6 + i, 6 + j));
            }
        }
        let mut u2 = f[6..].to_vec();
        crate::dense::Cholesky::new(&k22).unwrap().solve(&mut u2);
        let stretch = dot([u2[0], u2[1], u2[2]], axis);
        let want = l * (1.0 - 2.0 * nu) * mean / section.young;
        assert!((stretch - want).abs() < 1e-12 * want.abs(), "{stretch} vs {want}");
        let mut u = vec![0.0; 12];
        u[6..].copy_from_slice(&u2);
        let forces = beam.forces(&u, [0.0; 3], Some(pressure)).unwrap();
        for end in forces.ends {
            assert!((end[0] - mean * area).abs() < 1e-9 * scale, "{forces:?}");
            assert!(end[1..3].iter().all(|v| v.abs() < 1e-9 * scale), "{forces:?}");
        }
        let stress = beam.pipe_stress(&forces, pressure).unwrap();
        for (got, want) in stress.radial.iter().chain(&stress.hoop).chain(&stress.axial).zip([-100.0, -40.0, 60.0, 0.0, mean, mean]) {
            assert!((got - want).abs() < 1e-9 * 100.0, "{stress:?}");
        }
        // Held at both ends: the wall keeps the Poisson force N = 2νA A_w.
        let held = beam.forces(&[0.0; 12], [0.0; 3], Some(pressure)).unwrap();
        assert!(held.ends.iter().all(|end| (end[0] - 2.0 * nu * mean * area).abs() < 1e-12 * scale), "{held:?}");
        // Without a pressure nothing changes; a profile without a bore takes none.
        assert!(beam.forces(&[0.0; 12], [0.0; 3], None).unwrap().ends.iter().flatten().all(|v| *v == 0.0));
        let rod = Beam::new(&xyz, [0.0, 0.0, 1.0], BeamSection::new(&BeamProfile::Circle { radius: 20.0 }, &STEEL, BeamTheory::Timoshenko).unwrap())
            .unwrap();
        assert!(rod.pressure_load(pressure).is_err() && rod.forces(&u, [0.0; 3], Some(pressure)).is_err());
        assert!(beam.pressure_load(PipePressure { internal: 1e308, external: -1e308 }).is_err());
    }
}
