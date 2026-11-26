//! Absorbing layers next to an open face: a sponge of relaxation toward the face's target that
//! lets outgoing waves leave with little reflection (Xu & Sagaut, J. Comput. Phys. 245:14–42, 2013,
//! type II, eq. 24–28, 34; the polynomial profile of M. Israeli, S. A. Orszag, J. Comput. Phys.
//! 41:115–135, 1981). The layer's strength comes from a stated reflection, the profile's shape
//! from Israeli & Orszag's polynomial (their Table 1, n = 4).
//!
//! # Scheme
//! Each layer cell holds a rate λ per time step and a target (ρ_f, u_f). Around the cell's own
//! collision the cell takes two half steps with `h = tanh(λ/2)`:
//!
//! ```text
//! before the collision (implicit):  ρ* = (ρ + h ρ_f)/(1 + h),  j* = (j + h ρ_f u_f)/(1 + h),  u* = j*/ρ*
//!                                   f_q ← f_q + h [f_q^eq(ρ_f, u_f) − f_q^eq(ρ*, u*)]
//! the cell's collision, unchanged
//! after the collision (explicit):   f_q ← f_q + h [f_q^eq(ρ_f, u_f) − f_q^eq(ρ, u)]
//! ```
//!
//! The source is split into the two halves, not added whole. Its weighting is Guo, Zheng and Shi's
//! (Phys. Rev. E 65:046308, 2002), the one every collision here takes its body force with, so no
//! collision changes. Xu & Sagaut write the term for the BGK and MRT collisions with an
//! effective rate that depends on the collision: with two relaxation times the density and the
//! momentum would then be damped at different rates, and the difference couples the outgoing wave
//! to the returning one inside the layer. Here each conserved moment relaxes toward the target by
//! the same factor `(1 − h)/(1 + h) = e^{−λ}` per step, whatever the collision. The exponent is
//! exact, so the factor is monotone for every λ ≥ 0 and there is no strength bound (the paper's
//! χ < 4/s is the bound of its merged form, which this split does not have).
//!
//! # Strength
//! Cell `d` of a layer of `N` cells, `d = 0` at the face and `d = N − 1` at the inner edge, has the
//! centre at `ξ_d = (N − d − ½)/N` from the inner edge over the width. Its rate is
//! `λ_d = Λ s(ξ_d) / Σ_e s(ξ_e)`, with `s(ξ) = (1 − ξ) ξ⁴` and `Λ = c_s ln(1/R) / 2`, so that a
//! wave at rest comes back from the layer with `exp(−2 Σλ / c_s) = R` exactly. The constant of
//! Israeli & Orszag's polynomial cancels in the normalisation. A stream with normal Mach number
//! `M_n` returns `R^{1/(1 − M_n²)}`, since the outgoing wave runs at `c_s + U_n` and the returning
//! one at `c_s − U_n`; the setup reports that value as the stream's reflection.
//!
//! # What the layer itself returns
//! The viscous term sees the layer's rate as a complex frequency, so the wave impedance varies
//! along the layer and sends back a little more: at most `ν_L λ_max / (2 c_s²)`, where `ν_L` is the
//! longitudinal viscosity. The setup reports that bound. The grid's own part is not derived here:
//! its size is measured by [`crate::reflection::layer_reflection`]. Theory gives only the
//! resolution statement that a wave's amplitude must decay by less than a cell's worth per cell,
//! `λ_max ≤ c_s`; the setup warns when a layer falls short of it.
//! Oblique waves return more than `R`: damping every variable equally is a sponge, not a perfectly
//! matched layer.
//!
//! # Where this differs from the paper
//! The source's weighting is Guo's, not the paper's `δt F_i`; its eq. 48 (the sound speed) is not
//! reproduced; the rate uses `tanh(λ/2)` in place of the paper's χ, which keeps the layer's factor
//! exactly `e^{−λ}` for every strength.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.

use fieldcraft_model::{BoxFace, Collision, FaceCondition};

use crate::flow::filled;
use crate::lattice::{CS2, Q, c, equilibrium};
use crate::setup::Setup;
use crate::{Error, Result};

/// The profile's shape: ξ runs from 0 at the layer's inner edge to 1 at the face.
pub fn shape(xi: f64) -> f64 {
    (1.0 - xi) * xi * xi * xi * xi
}

/// The rate per step of each of a layer's `cells` cells, `[0]` at the face. Their sum is
/// `c_s ln(1/R) / 2`, so the layer at rest returns `R`.
pub fn rates(cells: usize, reflection: f64) -> Result<Vec<f64>> {
    if cells == 0 {
        return Err(Error::Invalid("an absorbing layer needs at least one cell".into()));
    }
    if !(reflection.is_finite() && reflection > 0.0 && reflection < 1.0) {
        return Err(Error::Invalid(format!("an absorbing layer's reflection must be above 0 and below 1, not {reflection}")));
    }
    let n = cells as f64;
    let mut shapes = filled(cells, 0.0, "absorbing layer profile")?;
    for (d, s) in shapes.iter_mut().enumerate() {
        *s = shape((n - d as f64 - 0.5) / n);
    }
    let total: f64 = shapes.iter().sum();
    let strength = CS2.sqrt() * libm::log(1.0 / reflection) / 2.0;
    for s in shapes.iter_mut() {
        *s *= strength / total;
    }
    Ok(shapes)
}

/// The longitudinal kinematic viscosity 4ν/3 + ζ a plane sound wave decays with (lattice units).
pub fn longitudinal_viscosity(collision: Collision, nu: f64) -> f64 {
    match collision {
        Collision::Mrt => 4.0 / 3.0 * nu + crate::mrt::bulk_viscosity(),
        _ => 2.0 * nu,
    }
}

/// What one layer is predicted to send back (dimensionless).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LayerReport {
    pub face: BoxFace,
    pub cells: usize,
    /// The stated reflection: the round trip through the layer at rest, exp(−2 Σλ / c_s).
    pub reflection: f64,
    /// The same with the target's stream along the face's normal: reflection^(1 / (1 − M_n²)).
    pub stream_reflection: f64,
    /// At most this much more comes back from the layer itself through the fluid's viscosity: ν_L λ_max / (2 c_s²).
    pub viscous_reflection: f64,
    /// The largest rate per step, λ_max.
    pub peak_rate: f64,
}

/// One layer in lattice units, as the setup holds it.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerSetup {
    pub report: LayerReport,
    /// Each cell's rate per step, `[0]` at the face.
    pub rates: Vec<f64>,
    /// The target in lattice units: the density, and the velocity (a velocity face's: its peak, before the profile).
    pub density: f64,
    pub velocity: [f64; 3],
}

/// The domain's layers in lattice units, and the setup's warnings about them.
pub fn layers(setup: &Setup) -> Result<(Vec<LayerSetup>, Vec<String>)> {
    let cs = CS2.sqrt();
    let mut out = Vec::with_capacity(setup.domain.layers.len());
    let mut warnings = Vec::new();
    for l in &setup.domain.layers {
        let face = l.face.name();
        let a = l.face.axis();
        let cells_f = (l.width / setup.grid.spacing).round();
        let extent = setup.grid.cells.get(a).copied().unwrap_or(0);
        let cells = cells_f as usize;
        if !(cells_f >= 1.0 && cells <= extent) {
            return Err(Error::Invalid(format!("the absorbing layer on {face} is {cells_f} cells wide; it must be 1 to {extent} cells")));
        }
        let (density, velocity) = match setup.domain.faces.get(l.face) {
            FaceCondition::Velocity { velocity, .. } => (1.0, setup.to_lattice(*velocity)),
            FaceCondition::Pressure { pressure, .. } => {
                let v = l.reference_velocity.ok_or_else(|| Error::Invalid(format!("the absorbing layer on {face} needs a reference velocity at its pressure face")))?;
                (setup.lattice_density(*pressure), setup.to_lattice(v))
            }
            _ => return Err(Error::Invalid(format!("the absorbing layer on {face} needs a velocity or pressure face"))),
        };
        if !(density.is_finite() && density > 0.0) {
            return Err(Error::Invalid(format!("the absorbing layer on {face} pulls toward a density that is not positive: {density}")));
        }
        if !velocity.iter().all(|v| v.is_finite()) {
            return Err(Error::Invalid(format!("the absorbing layer on {face} pulls toward a velocity that is not finite")));
        }
        let rates = rates(cells, l.reflection)?;
        let peak = rates.iter().copied().fold(0.0, f64::max);
        let normal = velocity.get(a).copied().unwrap_or(0.0);
        let m2 = normal * normal / CS2;
        if m2 >= 1.0 {
            return Err(Error::Invalid(format!("the absorbing layer on {face} pulls toward a stream at or above the speed of sound")));
        }
        let stream_reflection = libm::exp(libm::log(l.reflection) / (1.0 - m2));
        let nu = (setup.tau_plus - 0.5) * CS2;
        let viscous_reflection = longitudinal_viscosity(setup.collision, nu) * peak / (2.0 * CS2);
        if peak > cs {
            let need = (cells_f * peak / cs).ceil();
            warnings.push(format!(
                "absorbing layer on {face}: {cells} cells do not resolve a reflection of {}: its rate peaks at {peak:.2} per step, above the speed of sound of {cs:.3} cells per step, so a wave falls by more than 1/e within one cell and the layer itself sends part of it back; {need} cells or more resolve it",
                l.reflection
            ));
        }
        if viscous_reflection > l.reflection {
            let need = (cells_f * viscous_reflection / l.reflection).ceil();
            warnings.push(format!(
                "absorbing layer on {face}: through the fluid's viscosity the layer itself sends back up to {viscous_reflection:.1e} of a wave, above the stated {}; about {need} cells bring that under it",
                l.reflection
            ));
        }
        let report = LayerReport { face: l.face, cells, reflection: l.reflection, stream_reflection, viscous_reflection, peak_rate: peak };
        out.push(LayerSetup { report, rates, density, velocity });
    }
    Ok((out, warnings))
}

/// A cell of one or more layers (`Flow`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AbsorbingCell {
    pub cell: u32,
    pub half: f64,
    pub density: f64,
    pub velocity: [f64; 3],
}

/// One half step of a layer cell's populations toward (density, velocity); `implicit` is the half before the
/// collision. Returns false and leaves `fq` alone when its density is not a positive number (the collision
/// then stops the run).
pub(crate) fn relax(fq: &mut [f64; Q], half: f64, density: f64, velocity: [f64; 3], implicit: bool) -> bool {
    let mut rho = 0.0;
    let mut j = [0.0; 3];
    for (q, f) in fq.iter().enumerate() {
        let cq = c(q);
        rho += f;
        for a in 0..3 {
            j[a] += f * cq[a];
        }
    }
    if !(rho > 0.0 && rho.is_finite()) {
        return false;
    }
    let target = [density * velocity[0], density * velocity[1], density * velocity[2]];
    let (r, m) = if implicit {
        let norm = 1.0 + half;
        (
            (rho + half * density) / norm,
            [(j[0] + half * target[0]) / norm, (j[1] + half * target[1]) / norm, (j[2] + half * target[2]) / norm],
        )
    } else {
        (rho, j)
    };
    let u = [m[0] / r, m[1] / r, m[2] / r];
    for (q, f) in fq.iter_mut().enumerate() {
        *f += half * (equilibrium(q, density, velocity) - equilibrium(q, r, u));
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::W;

    fn random_stream() -> impl FnMut() -> f64 {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        move || {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    #[test]
    fn the_rates_give_the_stated_reflection() {
        for cells in [1, 4, 40] {
            for reflection in [0.5, 0.01, 1e-6] {
                let r = rates(cells, reflection).unwrap();
                assert_eq!(r.len(), cells);
                assert!(r.iter().all(|&x| x > 0.0));
                let sum: f64 = r.iter().sum();
                let back = libm::exp(-2.0 * sum / CS2.sqrt()) / reflection;
                assert!((back - 1.0).abs() < 1e-12, "cells {cells}, R {reflection}: {back}");
            }
        }
        let r = rates(40, 0.5).unwrap();
        let peak = (0..40).max_by(|&a, &b| r[a].total_cmp(&r[b])).unwrap();
        assert!(peak == 7 || peak == 8, "peak at {peak}");
        assert!(r[0] < r[7]);
        assert!(rates(0, 0.5).is_err());
        assert!(rates(4, 0.0).is_err());
        assert!(rates(4, 1.0).is_err());
        assert!(rates(4, f64::NAN).is_err());
    }

    #[test]
    fn two_halves_move_the_conserved_moments_by_the_rates_share() {
        let mut random = random_stream();
        for lambda in [1e-3, 0.3, 6.0] {
            let h = libm::tanh(0.5 * lambda);
            let target = (1.02, [0.05, 0.01, -0.02]);
            for _ in 0..50 {
                let mut fq = [0.0; Q];
                for (q, f) in fq.iter_mut().enumerate() {
                    *f = W[q] * (0.9 + 0.2 * random());
                }
                let moments = |f: &[f64; Q]| {
                    let rho: f64 = f.iter().sum();
                    let mut j = [0.0; 3];
                    for (q, v) in f.iter().enumerate() {
                        for a in 0..3 {
                            j[a] += v * c(q)[a];
                        }
                    }
                    (rho, j)
                };
                let (rho0, j0) = moments(&fq);
                assert!(relax(&mut fq, h, target.0, target.1, true));
                assert!(relax(&mut fq, h, target.0, target.1, false));
                let (rho1, j1) = moments(&fq);
                let g = 1.0 - libm::exp(-lambda);
                assert!((rho1 - rho0 - g * (target.0 - rho0)).abs() < 1e-14);
                for a in 0..3 {
                    assert!((j1[a] - j0[a] - g * (target.0 * target.1[a] - j0[a])).abs() < 1e-14);
                }
            }
        }
    }

    #[test]
    fn a_cell_at_its_target_keeps_its_populations() {
        let (density, velocity) = (1.0, [0.03, -0.01, 0.0]);
        let mut fq = [0.0; Q];
        for (q, f) in fq.iter_mut().enumerate() {
            *f = equilibrium(q, density, velocity);
        }
        let before = fq;
        let h = libm::tanh(0.5 * 2.0);
        assert!(relax(&mut fq, h, density, velocity, true));
        assert!(relax(&mut fq, h, density, velocity, false));
        for q in 0..Q {
            assert!((fq[q] - before[q]).abs() < 1e-15, "q {q}");
        }
    }

    #[test]
    fn a_cell_without_a_positive_density_is_left_alone() {
        for fill in [0.0, -0.1, f64::NAN] {
            let mut fq = [fill; Q];
            let before: Vec<u64> = fq.iter().map(|v| v.to_bits()).collect();
            assert!(!relax(&mut fq, 0.3, 1.0, [0.0; 3], true));
            let after: Vec<u64> = fq.iter().map(|v| v.to_bits()).collect();
            assert_eq!(before, after);
        }
    }
}
