//! Tangent and transpose-adjoint derivatives of linear equilibrium.
//!
//! Differentiate Z(p) u(p) = f(p): Z u' = f' - Z' u. For a fixed
//! complex observation q = w^T u, solve Z^T lambda = w and contract
//! q' = lambda^T (f' - Z' u). There is deliberately no conjugation.
//! This is the discrete direct/adjoint construction of Haug, Choi and
//! Komkov, *Design Sensitivity Analysis of Structural Systems* (1986),
//! applied to the phasor convention in [`crate::complex`]. See also Adelman
//! and Haftka, NASA TM-86333 (1984), equations (4)-(7),
//! <https://ntrs.nasa.gov/citations/19850006877>.
//!
//! Inputs are already reduced to the same fixed free-DOF basis. The caller
//! must differentiate prescribed-displacement elimination, loads, damping,
//! and observation weights when those depend on the parameter. Derivatives
//! carry SI output units per SI parameter unit. These kernels do not infer
//! material derivatives or certify a truncated modal basis.

use crate::complex::{self, Solution, Vector};
use crate::sparse::Csr;
use crate::{Error, Limits, Result};

fn vector(v: &Vector, n: usize) -> Result<()> {
    if v.real.len() != n || v.imaginary.len() != n || v.real.iter().chain(&v.imaginary).any(|x| !x.is_finite()) {
        return Err(Error::Invalid("sensitivity vector must match the system and be finite".into()));
    }
    Ok(())
}

/// Form f' - (A' + i B') u without assuming a positive derivative matrix.
pub fn tangent_load(a: &Csr, b: &Csr, u: &Vector, force_derivative: &Vector, limits: &Limits) -> Result<Vector> {
    complex::validate(a, limits)?;
    complex::validate(b, limits)?;
    if a.n != b.n {
        return Err(Error::Invalid("sensitivity matrices must have the same dimension".into()));
    }
    vector(u, a.n)?;
    vector(force_derivative, a.n)?;
    let mut ar = vec![0.0; a.n];
    let mut ai = vec![0.0; a.n];
    let mut br = vec![0.0; a.n];
    let mut bi = vec![0.0; a.n];
    a.mul(&u.real, &mut ar);
    a.mul(&u.imaginary, &mut ai);
    b.mul(&u.real, &mut br);
    b.mul(&u.imaginary, &mut bi);
    let result = Vector {
        real: force_derivative.real.iter().zip(ar.iter().zip(&bi)).map(|(f, (a, b))| f - a + b).collect(),
        imaginary: force_derivative.imaginary.iter().zip(ai.iter().zip(&br)).map(|(f, (a, b))| f - a - b).collect(),
    };
    vector(&result, a.n)?;
    Ok(result)
}

/// Solve a preassembled tangent load using the original dynamic operator.
/// B=0 also covers linear static equilibrium. The returned residual is for
/// the derivative equation, not the primal equation.
pub fn direct(a: &Csr, b: &Csr, tangent: &Vector, limits: &Limits) -> Result<Solution> {
    complex::solve(a, b, tangent, limits)
}

/// Fixed complex linear observation: transpose (not Hermitian) adjoint.
/// Returns (real, imaginary) derivative. Variable weights need w'^T u added.
pub fn adjoint(a: &Csr, b: &Csr, weights: &Vector, tangent: &Vector, limits: &Limits) -> Result<[f64; 2]> {
    vector(tangent, a.n)?;
    let lambda = complex::solve(a, b, weights, limits)?.displacement;
    let mut out = [0.0, 0.0];
    for ((lr, li), (tr, ti)) in lambda.real.iter().zip(&lambda.imaginary).zip(tangent.real.iter().zip(&tangent.imaginary)) {
        out[0] += lr * tr - li * ti;
        out[1] += lr * ti + li * tr;
    }
    if out.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("overflow in sensitivity adjoint contraction".into()));
    }
    Ok(out)
}

/// Derivatives of magnitude and locally unwrapped phase (radians).
/// At zero response neither derivative is generally defined; refuse it.
/// Scaling by |u| avoids forming |u|², which can overflow near resonance.
pub fn amplitude_phase(value: [f64; 2], derivative: [f64; 2]) -> Result<[f64; 2]> {
    if value.iter().chain(&derivative).any(|x| !x.is_finite()) {
        return Err(Error::Invalid("phasor and derivative must be finite".into()));
    }
    let magnitude = libm::hypot(value[0], value[1]);
    if magnitude == 0.0 || !magnitude.is_finite() {
        return Err(Error::Invalid("amplitude/phase derivative requires a finite nonzero response".into()));
    }
    let r = value[0] / magnitude;
    let i = value[1] / magnitude;
    let out = [r * derivative[0] + i * derivative[1], (r * derivative[1] - i * derivative[0]) / magnitude];
    if out.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("overflow in amplitude/phase derivative".into()));
    }
    Ok(out)
}

mod model;
pub use model::{
    FrequencyDerivative, FrequencySensitivity, ModeShapeDerivative, StaticObservation, StaticSensitivity, frequency_sensitivity,
    frequency_sensitivity_with_vectors, static_sensitivity,
};

mod inverse;
pub use inverse::{FrequencyEstimate, FrequencyFit, inverse_frequency, inverse_frequency_uniform_young};

mod harmonic;
pub use harmonic::{HarmonicDerivative, HarmonicObservation, HarmonicSensitivity, harmonic_sensitivity, inverse_harmonic};

mod nelson;

mod joint;
pub use joint::inverse_joint_uniform_material;
pub use joint::{JointEstimate, JointFit, JointObservation, JointParameter, JointResponse, fit_joint, inverse_joint};

mod pairing;
pub use pairing::{ModePair, ModePairing, pair_modes};
