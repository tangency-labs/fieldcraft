//! Pressure inside and around a pipe: Lamé's stresses in its wall, and the loads the pressure
//! puts on the pipe's beam ([`crate::beam`]).
//!
//! A pipe is a beam whose profile is a pipe with a bore, of inner radius a and outer radius b,
//! under uniform pressures p_i inside and p_o outside.
//!
//! **Wall stresses.** Lamé's solution of the thick-walled cylinder (S. P. Timoshenko,
//! J. N. Goodier, *Theory of Elasticity*, 3rd ed., McGraw-Hill 1970, §28): σ_r = A − B/r² and
//! σ_θ = A + B/r², with A = (p_i a² − p_o b²)/(b² − a²) and B = (p_i − p_o) a² b²/(b² − a²), so
//! σ_r is −p_i at the bore and −p_o outside. A thin wall is its limit: for mean radius r_m and
//! wall t the hoop stress at the bore is (p r_m/t)(1 + (t/2r_m)²).
//!
//! **Loads on the beam.** σ_r + σ_θ = 2A is the same everywhere in the wall, so by Hooke's law
//! the axial strain ε_x = (σ_x − 2νA)/E is uniform, as a beam's is: the pressure acts on the
//! beam as the initial strain ε₀ = −2νA/E, the axial force pair E A_w ε₀ = −2νA A_w on the
//! wall's area A_w. The fluid also presses on whatever closes the pipe: an element is a free
//! body of wall and enclosed fluid cut at its two ends, pulled outwards at each by the end force
//! F = p_i π a² − p_o π b² = A A_w (the closed-end condition; C. P. Sparks, *The influence of
//! tension, pressure and weight on pipe and riser deformations and stresses*, J. Energy
//! Resources Technology 106(1):46–54, 1984). Between two elements in line the end forces cancel;
//! at a closed end they are the force on the cap, and at a bend the resultant of the pressure on
//! the bend's wall. Together an element carries the axial pair (1 − 2ν) A A_w: a closed pipe
//! free to stretch does so by (1 − 2ν) A/E with the wall force N = A A_w (σ_x = A), and one held
//! at both ends keeps N = 2νA A_w. N is the force in the wall; the effective force N − F is what
//! the beam's equilibrium with its supports carries (Sparks 1984).
//!
//! Not modelled: the ovalisation of bends (elbows), and the stiffening of a pressurised pipe
//! against bending, which is second order.

use fieldcraft_model::BeamProfile;

use crate::{Error, Result};

/// Uniform pressures inside and around a pipe, positive pressing on the wall.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PipePressure {
    pub internal: f64,
    pub external: f64,
}

/// The wall stresses of a pipe under pressure: Lamé's at the bore and at the outside (in that
/// order), and the wall's mean axial stress N/A at the element's two ends.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PipeStress {
    /// Radial stress: −p_i at the bore, −p_o outside.
    pub radial: [f64; 2],
    /// Hoop stress.
    pub hoop: [f64; 2],
    /// Axial stress N/A at end 1 and end 2 (bending adds ± M c/I to it).
    pub axial: [f64; 2],
}

/// A pipe's wall: inner radius a and outer radius b, 0 < a < b.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PipeWall {
    inner: f64,
    outer: f64,
}

impl PipeWall {
    /// The wall of a pipe profile that has a bore; `None` for every other profile.
    pub fn of(profile: &BeamProfile) -> Option<PipeWall> {
        profile.bore().map(|(inner, outer)| PipeWall { inner, outer })
    }

    /// The inner and outer radii.
    pub fn radii(&self) -> [f64; 2] {
        [self.inner, self.outer]
    }

    /// Lamé's constants A and B.
    fn lame(&self, p: PipePressure) -> (f64, f64) {
        let (a2, b2) = (self.inner * self.inner, self.outer * self.outer);
        let ring = b2 - a2;
        ((p.internal * a2 - p.external * b2) / ring, (p.internal - p.external) * a2 * b2 / ring)
    }

    /// Lamé's constant A = (σ_r + σ_θ)/2, the same at every radius: the wall's axial stress
    /// under the closed-end force.
    pub fn mean_transverse(&self, p: PipePressure) -> Result<f64> {
        let (a, _) = self.lame(p);
        if a.is_finite() { Ok(a) } else { Err(self.not_finite(p)) }
    }

    /// The radial and hoop stresses at radius `r` in the wall.
    pub fn stresses(&self, p: PipePressure, r: f64) -> Result<[f64; 2]> {
        if !(r >= self.inner && r <= self.outer) {
            return Err(Error::Invalid(format!("radius {r} is outside the pipe's wall ({} to {})", self.inner, self.outer)));
        }
        let (a, b) = self.lame(p);
        let stresses = [a - b / (r * r), a + b / (r * r)];
        if stresses.iter().all(|s| s.is_finite()) { Ok(stresses) } else { Err(self.not_finite(p)) }
    }

    fn not_finite(&self, p: PipePressure) -> Error {
        Error::Invalid(format!(
            "the wall stresses of a pipe of radii {} and {} under the pressures {} inside and {} outside are not finite",
            self.inner, self.outer, p.internal, p.external
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lame_stresses_meet_the_pressures_and_the_thin_wall_limit() {
        // a = 10, b = 20, p_i = 100, p_o = 40: A = −20, B = 8000.
        let wall = PipeWall::of(&BeamProfile::Pipe { radius: 20.0, thickness: 10.0 }).unwrap();
        let p = PipePressure { internal: 100.0, external: 40.0 };
        assert!((wall.mean_transverse(p).unwrap() + 20.0).abs() < 1e-12);
        let near = |got: [f64; 2], want: [f64; 2]| got.iter().zip(want).all(|(g, w)| (g - w).abs() < 1e-12 * 100.0);
        assert!(near(wall.stresses(p, 10.0).unwrap(), [-100.0, 60.0]));
        assert!(near(wall.stresses(p, 20.0).unwrap(), [-40.0, 0.0]));
        // Equilibrium of half the pipe: ∫ σ_θ dr = p_i a − p_o b (the integral of A + B/r² in closed form).
        let (a, b): (f64, f64) = (-20.0, 8000.0);
        let integral = a * (20.0 - 10.0) + b * (1.0 / 10.0 - 1.0 / 20.0);
        assert!((integral - (100.0 * 10.0 - 40.0 * 20.0)).abs() < 1e-12 * 1000.0, "{integral}");
        let mid = wall.stresses(p, 15.0).unwrap();
        assert!(near(mid, [a - b / 225.0, a + b / 225.0]));
        // Internal pressure alone at the bore: p (a² + b²)/(b² − a²).
        let inside = PipePressure { internal: 100.0, external: 0.0 };
        assert!(near(wall.stresses(inside, 10.0).unwrap(), [-100.0, 100.0 * 500.0 / 300.0]));
        // A thin wall, r_m = 100, t = 5: the bore's hoop stress is (p r_m/t)(1 + (t/2r_m)²).
        let thin = PipeWall::of(&BeamProfile::Pipe { radius: 102.5, thickness: 5.0 }).unwrap();
        let hoop = thin.stresses(PipePressure { internal: 5.0, external: 0.0 }, 97.5).unwrap()[1];
        assert!((hoop - 100.0 * (1.0 + 0.025_f64.powi(2))).abs() < 1e-10 * 100.0, "{hoop}");
        // What is not a wall, a radius outside it, and numbers that overflow.
        assert!(PipeWall::of(&BeamProfile::Circle { radius: 1.0 }).is_none());
        assert!(PipeWall::of(&BeamProfile::Pipe { radius: 1.0, thickness: 1.0 }).is_none());
        assert!(wall.stresses(p, 9.0).is_err() && wall.stresses(p, f64::NAN).is_err());
        let huge = PipePressure { internal: 1e308, external: -1e308 };
        assert!(wall.stresses(huge, 10.0).is_err() && wall.mean_transverse(huge).is_err());
    }
}
