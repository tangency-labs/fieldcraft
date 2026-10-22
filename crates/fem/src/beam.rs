//! The two-node Timoshenko beam (six degrees of freedom per node), with the section properties of
//! the model's beam profiles. (Axial-only members are the trusses of [`crate::truss`].)
//!
//! **Local axes.** x runs from node 1 to node 2; y is the section's orientation vector (the
//! profile's local 1 direction, its `width`) made normal to x; z = x × y (the profile's local 2
//! direction, its `height`). Degrees of freedom per node: translations along x, y, z, then
//! rotations about them; element matrices are built in local axes and rotated to the global ones
//! ([`crate::frame`]).
//!
//! **Beam.** Axial and torsional stiffness are the exact two-node ones (EA/L, GJ/L). Each bending
//! plane uses the interdependent interpolation of the Timoshenko beam: cubic deflection and
//! quadratic rotation, coupled through Φ = 12 EI / (G A_s L²) so the shear strain is constant and
//! the element is exact for end loads (any number of elements gives the nodal displacements of
//! Timoshenko beam theory; J. S. Przemieniecki, *Theory of Matrix Structural Analysis*,
//! McGraw-Hill 1968, §5.6; Z. Friedman, J. B. Kosmatka, *An improved two-node Timoshenko beam
//! finite element*, Computers & Structures 47(3):473–481, 1993, whose shape functions are used
//! here). With these functions w(ξ) and θ(ξ), ξ = x/L,
//! K = ∫ (EI θ'ᵀ θ' + G A_s (w' − θ)ᵀ (w' − θ)) dx and the consistent mass
//! M = ∫ (ρA N_wᵀ N_w + ρI N_θᵀ N_θ) dx (translational and rotary inertia) are integrated with
//! four Gauss points, exactly; Φ → 0 recovers the Euler–Bernoulli (Hermite) element. The
//! rotation about y is −dw/dx in the x–z plane, so that plane's rotations enter with a minus
//! sign. Lumped mass: half the mass on each node's translations, half the polar inertia on the
//! twist, and the consistent bending rotary diagonal scaled like the translations (Hinton, Rock &
//! Zienkiewicz 1976). End forces are K u minus the consistent loads of the element's body force,
//! reported as the section forces at each end (axial force N, shears V_y, V_z, torque T, moments
//! M_y, M_z on the cut whose outward normal is +x).
//!
//! **Section properties.** Area, second moments I_y = ∫ z² dA and I_z = ∫ y² dA, the polar
//! moment I_y + I_z (rotary inertia of torsion), the torsion constant J and the shear areas:
//! - rectangle: J from Saint-Venant's series, J = a b³/3 (1 − 192 b/(π⁵ a) Σ_{n odd}
//!   tanh(nπa/2b)/n⁵) for sides a ≥ b (S. P. Timoshenko, J. N. Goodier, *Theory of Elasticity*,
//!   3rd ed., McGraw-Hill 1970, §109); shear coefficient 10(1 + ν)/(12 + 11ν);
//! - circle and pipe: J = I_y + I_z; shear coefficients 6(1 + ν)/(7 + 6ν) and, for a tube of
//!   radius ratio m, 6(1 + ν)(1 + m²)² / ((7 + 6ν)(1 + m²)² + (20 + 12ν) m²) (G. R. Cowper,
//!   *The shear coefficient in Timoshenko's beam theory*, J. Appl. Mech. 33(2):335–340, 1966);
//! - I and box (thin-walled): J = Σ b t³/3 over the midline walls (open I) and Bredt's
//!   4 A_m² t / p (closed box); shear areas the walls parallel to the shear, between the midlines
//!   of the walls across it (the web for shear along z, 5/6 of the flanges along y; the two
//!   webs or the two flanges of a box) (W. D. Pilkey, *Analysis and Design of Elastic Beams*,
//!   Wiley 2002, ch. 6).

use fieldcraft_model::{BeamProfile, Elastic};

use crate::dense::Matrix;
use crate::dof::{DofSet, ElementDofs};
use crate::frame::{Axes, Vec3, cross, dot, norm, rotate_matrix, rotate_to_global, rotate_to_local, sub, to_local, translation, unit};
use crate::quadrature::gauss_legendre;
use crate::{Error, Result};

/// Cross-section properties in the beam's local axes (y along the profile's width, z along its
/// height).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    pub area: f64,
    /// ∫ z² dA: bending with deflection along z.
    pub iy: f64,
    /// ∫ y² dA: bending with deflection along y.
    pub iz: f64,
    /// Saint-Venant torsion constant.
    pub torsion: f64,
    /// Shear areas for shear along y and along z.
    pub shear_y: f64,
    pub shear_z: f64,
}

/// Saint-Venant's torsion constant of a solid rectangle.
fn rectangle_torsion(width: f64, height: f64) -> f64 {
    let (a, b) = (width.max(height), width.min(height));
    let pi = std::f64::consts::PI;
    let series: f64 = (0..50).map(|k| f64::from(2 * k + 1)).map(|n| (n * pi * a / (2.0 * b)).tanh() / n.powi(5)).sum();
    a * b * b * b / 3.0 * (1.0 - 192.0 * b / (pi.powi(5) * a) * series)
}

/// The section properties of a profile (`poisson` enters the shear coefficients).
pub fn properties(profile: &BeamProfile, poisson: f64) -> Properties {
    let pi = std::f64::consts::PI;
    let area = profile.area();
    match *profile {
        BeamProfile::Rectangle { width, height } => {
            let k = 10.0 * (1.0 + poisson) / (12.0 + 11.0 * poisson);
            Properties {
                area,
                iy: width * height.powi(3) / 12.0,
                iz: height * width.powi(3) / 12.0,
                torsion: rectangle_torsion(width, height),
                shear_y: k * area,
                shear_z: k * area,
            }
        }
        BeamProfile::Circle { radius } => {
            let i = pi * radius.powi(4) / 4.0;
            let k = 6.0 * (1.0 + poisson) / (7.0 + 6.0 * poisson);
            Properties { area, iy: i, iz: i, torsion: 2.0 * i, shear_y: k * area, shear_z: k * area }
        }
        BeamProfile::Pipe { radius, thickness } => {
            let inner = radius - thickness;
            let i = pi * (radius.powi(4) - inner.powi(4)) / 4.0;
            let m2 = (inner / radius).powi(2);
            let k = 6.0 * (1.0 + poisson) * (1.0 + m2).powi(2) / ((7.0 + 6.0 * poisson) * (1.0 + m2).powi(2) + (20.0 + 12.0 * poisson) * m2);
            Properties { area, iy: i, iz: i, torsion: 2.0 * i, shear_y: k * area, shear_z: k * area }
        }
        BeamProfile::I { height, width, web, flange } => {
            let inside = height - 2.0 * flange;
            Properties {
                area,
                iy: (width * height.powi(3) - (width - web) * inside.powi(3)) / 12.0,
                iz: (2.0 * flange * width.powi(3) + inside * web.powi(3)) / 12.0,
                torsion: (2.0 * width * flange.powi(3) + (height - flange) * web.powi(3)) / 3.0,
                shear_y: 5.0 / 6.0 * 2.0 * width * flange,
                shear_z: (height - flange) * web,
            }
        }
        BeamProfile::Box { width, height, thickness: t } => {
            let (wi, hi) = (width - 2.0 * t, height - 2.0 * t);
            let (wm, hm) = (width - t, height - t);
            Properties {
                area,
                iy: (width * height.powi(3) - wi * hi.powi(3)) / 12.0,
                iz: (height * width.powi(3) - hi * wi.powi(3)) / 12.0,
                torsion: 4.0 * (wm * hm).powi(2) * t / (2.0 * (wm + hm)),
                shear_y: 2.0 * wm * t,
                shear_z: 2.0 * hm * t,
            }
        }
    }
}

/// The Young's and shear moduli and Poisson's ratio of an isotropic material.
/// (`what` names the section kind in errors.)
fn moduli(elastic: &Elastic, what: &str) -> Result<(f64, f64, f64)> {
    let (young, poisson) = crate::material::isotropic_constants(elastic, what)?;
    Ok((young, young / (2.0 * (1.0 + poisson)), poisson))
}

/// A beam section: profile properties and an isotropic material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamSection {
    pub properties: Properties,
    pub young: f64,
    pub shear_modulus: f64,
}

impl BeamSection {
    pub fn new(profile: &BeamProfile, elastic: &Elastic) -> Result<BeamSection> {
        profile.validate().map_err(|e| Error::Invalid(e.to_string()))?;
        let (young, shear_modulus, poisson) = moduli(elastic, "beam")?;
        Ok(BeamSection { properties: properties(profile, poisson), young, shear_modulus })
    }
}

/// The section forces at the two ends of a beam: N, V_y, V_z, T, M_y, M_z in its local axes.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct BeamForces {
    pub axes: Axes,
    pub ends: [[f64; 6]; 2],
}

/// One bending plane's stiffness and mass (deflection, rotation at node 1, then node 2), with
/// θ = dw/dx in the plane's own convention.
fn bending_plane(ei: f64, gas: f64, rho_a: f64, rho_i: f64, l: f64) -> ([[f64; 4]; 4], [[f64; 4]; 4]) {
    let phi = 12.0 * ei / (gas * l * l);
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
            let (kp, mp) = bending_plane(e * i, g * shear, density * p.area, density * i, l);
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

    /// Section forces at the ends for the nodal displacements and rotations `u` (global, 12)
    /// under a uniform force per volume `body` acting on the element.
    pub fn forces(&self, u: &[f64], body: [f64; 3]) -> Result<BeamForces> {
        if u.len() != 12 {
            return Err(Error::Invalid(format!("{} values for a beam element with 12 degrees of freedom", u.len())));
        }
        let ku = self.local(0.0).0.mul_vec(&rotate_to_local(u, &self.axes));
        let end: Vec<f64> = ku.iter().zip(self.local_body_force(body)).map(|(a, b)| a - b).collect();
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
        let section = BeamSection::new(&profile, &STEEL).unwrap();
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
        let forces = beam.forces(&u, [0.0; 3]).unwrap();
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
    fn profile_properties() {
        // A square: J = 0.1406 a⁴ (Timoshenko & Goodier's table, 0.1406); thin strip J → b t³/3.
        let sq = properties(&BeamProfile::Rectangle { width: 2.0, height: 2.0 }, 0.3);
        assert!((sq.torsion / 16.0 - 0.1406).abs() < 1e-4, "{}", sq.torsion / 16.0);
        let strip = properties(&BeamProfile::Rectangle { width: 1000.0, height: 1.0 }, 0.3);
        assert!((strip.torsion / (1000.0 / 3.0) - 1.0).abs() < 1e-3, "{}", strip.torsion);
        // A thin pipe: the shear coefficient tends to 2(1 + ν)/(4 + 3ν), J = 2 π r³ t.
        let (r, t) = (100.0, 0.01);
        let pipe = properties(&BeamProfile::Pipe { radius: r, thickness: t }, 0.3);
        let k = pipe.shear_y / pipe.area;
        assert!((k - 2.6 / 4.9).abs() < 1e-3, "{k}");
        assert!((pipe.torsion / (2.0 * std::f64::consts::PI * r.powi(3) * t) - 1.0).abs() < 1e-3);
        // A box is a pair of I sections in bending: the same I_y with webs of twice the thickness.
        let b = properties(&BeamProfile::Box { width: 10.0, height: 20.0, thickness: 1.0 }, 0.3);
        let i = properties(&BeamProfile::I { height: 20.0, width: 10.0, web: 2.0, flange: 1.0 }, 0.3);
        assert!((b.iy - i.iy).abs() < 1e-9 && (b.area - i.area).abs() < 1e-12);
        let circle = properties(&BeamProfile::Circle { radius: 2.0 }, 0.3);
        assert!((circle.torsion - std::f64::consts::PI * 8.0).abs() < 1e-12);
    }
}
