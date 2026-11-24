//! Stress-life (S-N) fatigue at a material point: the scalar stress measures of a stress tensor,
//! rainflow cycle counting of a load history, the S-N curve's life at a stress amplitude, mean
//! stress corrections, Palmgren–Miner damage and the safety factor against the endurance limit.
//!
//! The loading is proportional: every stress is the reference stress σ_ref of one linear static
//! solution times a load factor λ(t). A signed scalar measure m of σ_ref then gives the stress
//! history m λ(t), and the cycles of that history are the cycles of λ(t) scaled by m: rainflow
//! counting looks only at ranges, so a negative m reverses the history and keeps its ranges.
//! Cycles are counted once on the load factors ([`rainflow`], [`rainflow_repeated`]) and read at
//! every node with its own m ([`assess`]). Non-proportional (multiaxial) loading, strain-life and
//! critical-plane methods are outside this module.
//!
//! Sources, each formula as stated there:
//! - **S-N curve.** Basquin's relation in reversals, σ_a = σ_f′ (2N)^b: O. H. Basquin, *The
//!   exponential law of endurance tests*, Proc. ASTM 10 (1910) 625–630; N. E. Dowling, *Mechanical
//!   Behavior of Materials*, 4th ed., Pearson, 2013, chapter 9. Inverted, N = ½ (σ_a/σ_f′)^(1/b).
//!   A tabulated curve is a straight line between its points on log-log axes (Basquin's relation
//!   piece by piece); beyond its ends the end segments are extended. Below the endurance limit the
//!   life is infinite (the knee of a steel's curve).
//! - **Mean stress.** The equivalent fully reversed amplitude σ_ar a cycle of amplitude σ_a and
//!   mean σ_m > 0 is read at (Dowling, chapter 9; J. A. Bannantine, J. J. Comer, J. L. Handrock,
//!   *Fundamentals of Metal Fatigue Analysis*, Prentice Hall, 1990, chapter 2): Goodman
//!   σ_ar = σ_a / (1 − σ_m/S_u), Gerber σ_ar = σ_a / (1 − (σ_m/S_u)²), Soderberg
//!   σ_ar = σ_a / (1 − σ_m/S_y). A mean at or past the strength is a static failure. A compressive
//!   mean (σ_m ≤ 0) is given no benefit: σ_ar = σ_a, the conservative choice for the lines, which
//!   were drawn through tensile means (R. G. Budynas, J. K. Nisbett, *Shigley's Mechanical
//!   Engineering Design*, 10th ed., McGraw-Hill, 2015, the chapter on fatigue failure from variable loading). The lines are Goodman's (J. Goodman,
//!   *Mechanics Applied to Engineering*, Longmans, Green, London, 1899), Gerber's (W. Gerber, Z.
//!   Bayer. Architekten- u. Ingenieur-Vereins 6 (1874) 101–110) and Soderberg's (C. R. Soderberg,
//!   *Factor of safety and working stress*, Trans. ASME 52 (APM-52-2) (1930) 13–28).
//!   Smith–Watson–Topper (SWT): σ_ar = √(σ_max σ_a) with σ_max = σ_a + σ_m (K. N. Smith, P. Watson,
//!   T. H. Topper, *A stress-strain function for the fatigue of metals*, J. Materials 5(4) (1970)
//!   767–778, in its stress-life form as in Dowling, chapter 9). It needs no strength. It is read
//!   at every mean, a compressive one included (the parameter itself falls with σ_max); a cycle
//!   whose σ_max is not tensile does no damage (σ_ar = 0).
//! - **Safety factor** against the endurance limit S_e: the factor n that scales the cycle (its
//!   amplitude and mean together) onto the failure line (Shigley, the fatigue failure criteria for fluctuating stress, n_f):
//!   no correction n = S_e/σ_a; Goodman 1/n = σ_a/S_e + σ_m/S_u; Soderberg 1/n = σ_a/S_e + σ_m/S_y;
//!   Gerber n = ½ (S_u/σ_m)² (σ_a/S_e) [−1 + √(1 + (2 σ_m S_e/(S_u σ_a))²)]; with a compressive
//!   mean n = S_e/σ_a; SWT n = S_e/√(σ_max σ_a) (the cycle scaled so that √(σ_max σ_a) = S_e).
//!   A history's factor is its worst cycle's.
//! - **Rainflow counting.** T. Matsuishi, T. Endo, *Fatigue of metals subjected to varying
//!   stress*, Japan Soc. Mech. Eng., 1968; the three-point procedure of ASTM E1049-85 (2017),
//!   *Standard Practices for Cycle Counting in Fatigue Analysis* ([`rainflow`]); for a
//!   repeating history, S. D. Downing, D. F. Socie, *Simple rainflow counting algorithms*, Int. J.
//!   Fatigue 4 (1982) 31–40, their method for a repeating history: the history starts at its largest peak or deepest
//!   valley (whichever is larger in magnitude) and runs round to it again, so every cycle closes
//!   ([`rainflow_repeated`]).
//! - **Damage.** The Palmgren–Miner linear rule, D = Σ nᵢ/Nᵢ, failure at D = 1: M. A. Miner,
//!   *Cumulative damage in fatigue*, J. Appl. Mech. 12 (1945) A159–A164.

use fieldcraft_model::{Fatigue, MeanStress, SnCurve, StressMeasure};

use crate::{Error, Result, principal, von_mises};

/// A counted cycle of a load history: its range (largest less smallest value), its mean and how
/// many times it is counted (1 for a full cycle, ½ for a half cycle).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cycle {
    pub range: f64,
    pub mean: f64,
    pub count: f64,
}

impl Cycle {
    fn between(a: f64, b: f64, count: f64) -> Cycle {
        Cycle { range: (a - b).abs(), mean: 0.5 * (a + b), count }
    }
}

/// The relative band around a tie of the two extreme principal stresses (pure shear) in which a
/// state counts as tensile: the sign of a signed measure is that of the principal stress of
/// largest magnitude, and where the two are nearly equal and opposite, rounding noise or a
/// fraction of a per cent in a nodal stress would otherwise set it, and with it whether the mean
/// is tensile. Tension is the conservative side for every correction here.
pub const SIGN_TIE_BAND: f64 = 1e-3; // unit-factors: a relative tolerance on a sign tie, not a conversion

/// The scalar stress `measure` takes of the stress `t` (xx, yy, zz, xy, xz, yz), signed. The sign
/// is that of the principal stress of largest magnitude, except that a state whose largest
/// tensile principal stress is within [`SIGN_TIE_BAND`] (0.1 %) of the largest compressive one in
/// magnitude, or above it, counts as tensile (so pure shear does, and does not flip with noise).
/// A state with no tensile principal stress is compressive.
pub fn measure(measure: StressMeasure, t: &[f64; 6]) -> f64 {
    let [hi, _, lo] = principal(t);
    let tensile = hi.abs() >= lo.abs() || (hi > 0.0 && hi >= lo.abs() * (1.0 - SIGN_TIE_BAND));
    match measure {
        StressMeasure::SignedVonMises => {
            let v = von_mises(t);
            if tensile { v } else { -v }
        }
        StressMeasure::MaxAbsPrincipal => {
            if tensile {
                hi
            } else {
                lo
            }
        }
    }
}

/// The peaks and valleys of `series`, in order: a point between two others on a monotonic run
/// is dropped, and repeated values are kept once. The first and last points are kept.
pub fn reversals(series: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::with_capacity(series.len());
    for &x in series {
        match out.as_mut_slice() {
            [.., last] if *last == x => {}
            [.., a, b] if (*b - *a) * (x - *b) > 0.0 => *b = x,
            _ => out.push(x),
        }
    }
    out
}

/// Ranges counted by the three-point rule: each point of `points` is pushed and, while the most
/// recent range X is at least the one before it, Y, Y is counted and removed. `half_at_start`
/// applies ASTM E1049's starting-point rule (a Y holding the history's first point is counted as
/// half a cycle and only that point is removed); without it every Y is a full cycle. Returns
/// the cycles and the points left over.
fn three_point(points: &[f64], half_at_start: bool) -> (Vec<Cycle>, Vec<f64>) {
    let mut cycles = Vec::new();
    let mut stack: Vec<f64> = Vec::with_capacity(points.len());
    for &p in points {
        stack.push(p);
        while let [.., a, b, c] = stack.as_slice() {
            let (x, y) = ((c - b).abs(), (b - a).abs());
            if x < y {
                break;
            }
            let (a, b, c) = (*a, *b, *c);
            let n = stack.len();
            if half_at_start && n == 3 {
                cycles.push(Cycle::between(a, b, 0.5));
                stack.remove(0);
            } else {
                cycles.push(Cycle::between(a, b, 1.0));
                stack.truncate(n - 3);
                stack.push(c);
            }
        }
    }
    (cycles, stack)
}

/// Rainflow cycles of a history that happens once (the three-point procedure of ASTM E1049-85): full cycles where
/// they close, the ranges left at the end as half cycles. Equal cycles are combined
/// ([`combine`]).
pub fn rainflow(series: &[f64]) -> Vec<Cycle> {
    let (mut cycles, rest) = three_point(&reversals(series), true);
    cycles.extend(rest.windows(2).filter_map(|w| match w {
        [a, b] => Some(Cycle::between(*a, *b, 0.5)),
        _ => None,
    }));
    combine(cycles)
}

/// Rainflow cycles of one block of a history that repeats (Downing & Socie's method for a repeating history): the
/// block is rotated to start at its value of largest magnitude and closed back to it, so every
/// cycle is a full one, the block's join included. Equal cycles are combined.
pub fn rainflow_repeated(series: &[f64]) -> Vec<Cycle> {
    let points = reversals(series);
    let start = points.iter().enumerate().fold(None, |best: Option<(usize, f64)>, (i, v)| match best {
        Some((_, b)) if v.abs() <= b => best,
        _ => Some((i, v.abs())),
    });
    let Some((k, _)) = start else { return Vec::new() };
    let (head, tail) = points.split_at(k);
    let closed: Vec<f64> = tail.iter().chain(head).chain(tail.first()).copied().collect();
    let (mut cycles, rest) = three_point(&reversals(&closed), false);
    // Nothing is left but the starting point; anything else would be counted as a half cycle.
    cycles.extend(rest.windows(2).filter_map(|w| match w {
        [a, b] => Some(Cycle::between(*a, *b, 0.5)),
        _ => None,
    }));
    combine(cycles)
}

/// The cycles with every set of equal ones (the same range and mean) made one with their counts
/// summed, largest range first (then lowest mean): a repeated block reads the curve once per
/// distinct cycle.
pub fn combine(mut cycles: Vec<Cycle>) -> Vec<Cycle> {
    cycles.retain(|c| c.range > 0.0);
    cycles.sort_by(|a, b| b.range.total_cmp(&a.range).then(a.mean.total_cmp(&b.mean)));
    let mut out: Vec<Cycle> = Vec::with_capacity(cycles.len());
    for c in cycles {
        match out.last_mut() {
            Some(last) if last.range == c.range && last.mean == c.mean => last.count += c.count,
            _ => out.push(c),
        }
    }
    out
}

/// The smallest life the curve gives: half a cycle (one reversal), the static end of Basquin's
/// relation (2N = 1 at σ_a = σ_f′). A cycle beyond it is a static failure.
pub const MIN_LIFE: f64 = 0.5;

/// The life the S-N curve gives at an equivalent fully reversed amplitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Life {
    /// Cycles to failure, at least [`MIN_LIFE`]; infinite at or below the endurance limit.
    pub cycles: f64,
    /// The amplitude is at or beyond the curve's static end (or the mean at or past the
    /// strength): failure in the first reversal, not fatigue.
    pub static_failure: bool,
    /// The amplitude is beyond a tabulated curve's points: its end segment was extended.
    pub extrapolated: bool,
}

/// The cycles to failure the curve of `f` gives at the equivalent fully reversed amplitude
/// `amplitude` (≥ 0).
pub fn life(f: &Fatigue, amplitude: f64) -> Life {
    let infinite = Life { cycles: f64::INFINITY, static_failure: false, extrapolated: false };
    if !(amplitude > 0.0) || f.endurance_limit.is_some_and(|e| amplitude <= e) {
        return infinite;
    }
    if amplitude == f64::INFINITY {
        return Life { cycles: MIN_LIFE, static_failure: true, extrapolated: false };
    }
    let (cycles, extrapolated) = match &f.curve {
        SnCurve::Basquin { coefficient, exponent } => (0.5 * (amplitude / coefficient).powf(1.0 / exponent), false),
        SnCurve::Table { points } => {
            // The segment holding the amplitude, else the end segment beyond it.
            let segments: Vec<_> = points.windows(2).filter_map(|w| if let [a, b] = w { Some((*a, *b)) } else { None }).collect();
            let inside = segments.iter().find(|(a, b)| amplitude <= a.amplitude && amplitude >= b.amplitude);
            let pick = inside.or_else(|| if points.first().is_some_and(|p| amplitude > p.amplitude) { segments.first() } else { segments.last() });
            let Some((a, b)) = pick else { return infinite };
            let slope = (b.cycles / a.cycles).ln() / (b.amplitude / a.amplitude).ln();
            (a.cycles * (amplitude / a.amplitude).powf(slope), inside.is_none())
        }
    };
    if !(cycles > MIN_LIFE) {
        return Life { cycles: MIN_LIFE, static_failure: true, extrapolated };
    }
    Life { cycles, static_failure: false, extrapolated }
}

/// Whether `f` has the strengths `correction` needs; the error names what is missing.
pub fn check(f: &Fatigue, correction: MeanStress) -> std::result::Result<(), String> {
    match correction {
        MeanStress::Uncorrected => Ok(()),
        MeanStress::Goodman | MeanStress::Gerber if f.ultimate_strength.is_none() => {
            Err(format!("the {} correction needs the ultimate strength", correction.id()))
        }
        MeanStress::Soderberg if f.yield_strength.is_none() => Err("the Soderberg correction needs the yield strength".into()),
        _ => Ok(()),
    }
}

/// The equivalent fully reversed amplitude σ_ar of a cycle of amplitude `amplitude` (≥ 0) and
/// mean `mean` under `correction`; `None` when the mean reaches the strength the correction is
/// drawn to (a static failure). A strength the correction needs and `f` lacks is an error.
pub fn equivalent_amplitude(f: &Fatigue, correction: MeanStress, amplitude: f64, mean: f64) -> Result<Option<f64>> {
    check(f, correction).map_err(Error::Invalid)?;
    if correction == MeanStress::Swt {
        let max = amplitude + mean;
        return Ok(Some(if max > 0.0 { (max * amplitude).sqrt() } else { 0.0 }));
    }
    if mean <= 0.0 {
        return Ok(Some(amplitude));
    }
    let strength = |s: Option<f64>| s.ok_or_else(|| Error::Invalid("a strength is missing".into()));
    let reduction = match correction {
        MeanStress::Uncorrected | MeanStress::Swt => return Ok(Some(amplitude)),
        MeanStress::Goodman => 1.0 - mean / strength(f.ultimate_strength)?,
        MeanStress::Gerber => 1.0 - (mean / strength(f.ultimate_strength)?).powi(2),
        MeanStress::Soderberg => 1.0 - mean / strength(f.yield_strength)?,
    };
    Ok((reduction > 0.0).then(|| amplitude / reduction))
}

/// The factor that scales a cycle (amplitude and mean together) onto the endurance limit's
/// failure line of `correction`: above 1 the cycle is below the endurance limit. Infinite for a
/// cycle with no amplitude and no tensile mean; `None` without an endurance limit.
pub fn safety_factor(f: &Fatigue, correction: MeanStress, amplitude: f64, mean: f64) -> Result<Option<f64>> {
    check(f, correction).map_err(Error::Invalid)?;
    let Some(se) = f.endurance_limit else { return Ok(None) };
    let plain = if amplitude > 0.0 { se / amplitude } else { f64::INFINITY };
    if correction == MeanStress::Swt {
        let product = (amplitude + mean) * amplitude;
        return Ok(Some(if product > 0.0 { se / product.sqrt() } else { f64::INFINITY }));
    }
    if mean <= 0.0 {
        return Ok(Some(plain));
    }
    let strength = |s: Option<f64>| s.ok_or_else(|| Error::Invalid("a strength is missing".into()));
    Ok(Some(match correction {
        MeanStress::Uncorrected | MeanStress::Swt => plain,
        MeanStress::Goodman => 1.0 / (amplitude / se + mean / strength(f.ultimate_strength)?),
        MeanStress::Soderberg => 1.0 / (amplitude / se + mean / strength(f.yield_strength)?),
        MeanStress::Gerber => {
            let su = strength(f.ultimate_strength)?;
            if amplitude > 0.0 {
                let x = 2.0 * mean * se / (su * amplitude);
                0.5 * (su / mean).powi(2) * (amplitude / se) * ((1.0 + x * x).sqrt() - 1.0)
            } else {
                su / mean
            }
        }
    }))
}

/// The fatigue of one material point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    /// Palmgren–Miner damage of one block (one cycle of a constant amplitude loading).
    pub damage: f64,
    /// Blocks (cycles of a constant amplitude loading) to failure, 1 / damage: infinite with
    /// no damage, and 0 for a static failure (see [`Point::static_failure`]).
    pub life: f64,
    /// The safety factor against the endurance limit (its worst cycle's); `None` without an
    /// endurance limit.
    pub safety_factor: Option<f64>,
    /// The governing cycle (the one with the largest equivalent amplitude): its stress
    /// amplitude (the stress factor applied), mean and equivalent fully reversed amplitude
    /// (infinite for a static failure).
    pub amplitude: f64,
    pub mean: f64,
    pub equivalent: f64,
    /// A static failure: some cycle's equivalent amplitude is beyond the curve's static end or its
    /// mean is past the strength, or the history's peak stress reaches the ultimate strength
    /// (checked whenever the material has one, whatever the correction). The life is then 0,
    /// the damage at least 1/[`MIN_LIFE`] and the safety factor at most 1.
    pub static_failure: bool,
    /// Some cycle was read on an extended end segment of a tabulated curve.
    pub extrapolated: bool,
}

impl Point {
    /// log₁₀ of the life for display: a static failure is shown at the shortest life,
    /// [`MIN_LIFE`], and an infinite life stays infinite.
    pub fn log10_life(&self) -> f64 {
        self.life.max(MIN_LIFE).log10()
    }
}

/// The fatigue at a point whose reference stress measure is `m`, under the load factor
/// `cycles` (from [`rainflow`], [`rainflow_repeated`] or one constant amplitude cycle) of a
/// history whose smallest and largest load factors are `extent` (min, max), with every stress,
/// amplitude and mean alike, multiplied by `stress_factor` (the result is linear in the load).
/// The peak tensile stress `stress_factor · m · λ` over the history is checked against the
/// ultimate strength when the material has one. A measure that is not a finite number is an
/// error.
pub fn assess(f: &Fatigue, correction: MeanStress, stress_factor: f64, m: f64, cycles: &[Cycle], extent: (f64, f64)) -> Result<Point> {
    check(f, correction).map_err(Error::Invalid)?;
    if !m.is_finite() {
        return Err(Error::Invalid(format!("the stress measure {m} is not a finite number")));
    }
    let mut out = Point {
        damage: 0.0,
        life: f64::INFINITY,
        safety_factor: f.endurance_limit.map(|_| f64::INFINITY),
        amplitude: 0.0,
        mean: 0.0,
        equivalent: -1.0,
        static_failure: false,
        extrapolated: false,
    };
    for c in cycles {
        let amplitude = stress_factor * 0.5 * c.range * m.abs();
        let mean = stress_factor * c.mean * m;
        let (equivalent, life) = match equivalent_amplitude(f, correction, amplitude, mean)? {
            Some(e) => (e, life(f, e)),
            None => (f64::INFINITY, Life { cycles: MIN_LIFE, static_failure: true, extrapolated: false }),
        };
        out.damage += c.count / life.cycles;
        out.static_failure |= life.static_failure;
        out.extrapolated |= life.extrapolated;
        if let (Some(n), Some(worst)) = (safety_factor(f, correction, amplitude, mean)?, out.safety_factor.as_mut()) {
            *worst = worst.min(n);
        }
        if equivalent > out.equivalent {
            (out.amplitude, out.mean, out.equivalent) = (amplitude, mean, equivalent);
        }
    }
    if out.equivalent < 0.0 {
        out.equivalent = 0.0;
    }
    // The peak tensile stress of the history against the ultimate strength: a static overload
    // needs no cycle (a constant load has none), and a cycle with a small mean can still peak
    // past the strength.
    let peak = stress_factor * (m * extent.0).max(m * extent.1);
    if let Some(su) = f.ultimate_strength
        && peak >= su
    {
        out.static_failure = true;
        out.damage = out.damage.max(1.0 / MIN_LIFE);
        let n = su / peak;
        out.safety_factor = Some(out.safety_factor.map_or(n, |s| s.min(n)));
    }
    out.life = if out.static_failure {
        0.0
    } else if out.damage > 0.0 {
        1.0 / out.damage
    } else {
        f64::INFINITY
    };
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::SnPoint;

    /// An illustrative steel-like material in MPa (round numbers, not a real material's data).
    fn basquin(endurance: Option<f64>) -> Fatigue {
        Fatigue {
            curve: SnCurve::Basquin { coefficient: 1000.0, exponent: -0.1 },
            endurance_limit: endurance,
            ultimate_strength: Some(800.0),
            yield_strength: Some(600.0),
        }
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1e-300)
    }

    /// Basquin's relation read both ways: the life at σ_f′ (2N)^b is N, and the cycle form
    /// S = A N^b with A = σ_f′ 2^b is the same curve.
    #[test]
    fn the_basquin_curve_is_reproduced_exactly() {
        let f = basquin(None);
        for n in [1.0f64, 10.0, 1e3, 1e5, 1e7, 1e9] {
            let s = 1000.0 * (2.0 * n).powf(-0.1);
            let l = life(&f, s);
            assert!(close(l.cycles, n, 1e-12) && !l.static_failure, "{n}: {l:?}");
            let a = 1000.0 * 2f64.powf(-0.1);
            assert!(close(a * n.powf(-0.1), s, 1e-14));
        }
        // At σ_f′ the life is one reversal, above it a static failure.
        assert_eq!(life(&f, 1000.0).cycles, MIN_LIFE);
        assert!(life(&f, 1500.0).static_failure && life(&f, 1500.0).cycles == MIN_LIFE);
        // No endurance limit: the curve goes on; with one, infinite at and below it.
        assert!(life(&f, 100.0).cycles.is_finite());
        let e = basquin(Some(300.0));
        assert_eq!(life(&e, 300.0).cycles, f64::INFINITY);
        assert_eq!(life(&e, 299.0).cycles, f64::INFINITY);
        assert!(life(&e, 301.0).cycles.is_finite());
        assert_eq!(life(&e, 0.0).cycles, f64::INFINITY);
        assert_eq!(life(&f, f64::NAN).cycles, f64::INFINITY);
    }

    /// A tabulated curve is exact at its points and a straight line between them on log-log
    /// axes: at the geometric mean of two amplitudes the life is the geometric mean of theirs.
    #[test]
    fn a_table_interpolates_on_log_log_axes() {
        let f = Fatigue {
            curve: SnCurve::Table {
                points: vec![
                    SnPoint { cycles: 1e3, amplitude: 500.0 },
                    SnPoint { cycles: 1e5, amplitude: 300.0 },
                    SnPoint { cycles: 1e7, amplitude: 200.0 },
                ],
            },
            endurance_limit: None,
            ultimate_strength: None,
            yield_strength: None,
        };
        for (n, s) in [(1e3, 500.0), (1e5, 300.0), (1e7, 200.0)] {
            assert!(close(life(&f, s).cycles, n, 1e-12), "{s}");
        }
        let mid = life(&f, (300.0f64 * 200.0).sqrt());
        assert!(close(mid.cycles, 1e6, 1e-12) && !mid.extrapolated, "{mid:?}");
        // Beyond the ends: the end segments extended, and said so.
        let low = life(&f, 150.0);
        let slope = (1e7f64 / 1e5).ln() / (200.0f64 / 300.0).ln();
        assert!(close(low.cycles, 1e7 * (150.0f64 / 200.0).powf(slope), 1e-12) && low.extrapolated);
        assert!(life(&f, 600.0).extrapolated);
        let with_limit = Fatigue { endurance_limit: Some(180.0), ..f };
        assert_eq!(life(&with_limit, 170.0).cycles, f64::INFINITY);
    }

    /// The mean stress corrections against their formulas, and the safety factor's point on each
    /// failure line: (n σ_a, n σ_m) satisfies the line's equation.
    #[test]
    fn mean_stress_corrections_follow_their_formulas() {
        let f = basquin(Some(300.0));
        let (su, sy, se) = (800.0, 600.0, 300.0);
        let (a, m) = (150.0, 200.0);
        let eq = |c| equivalent_amplitude(&f, c, a, m).unwrap().unwrap();
        assert_eq!(eq(MeanStress::Uncorrected), a);
        assert!(close(eq(MeanStress::Goodman), a / (1.0 - m / su), 1e-15));
        assert!(close(eq(MeanStress::Gerber), a / (1.0 - (m / su).powi(2)), 1e-15));
        assert!(close(eq(MeanStress::Soderberg), a / (1.0 - m / sy), 1e-15));
        // SWT: σ_max = 350, σ_ar = √(350 · 150) = √52500.
        assert!(close(eq(MeanStress::Swt), 52500.0f64.sqrt(), 1e-15));
        // Hand values: Goodman 150/0.75 = 200, Gerber 150/0.9375 = 160, Soderberg 150/(2/3) = 225.
        assert!(close(eq(MeanStress::Goodman), 200.0, 1e-14));
        assert!(close(eq(MeanStress::Gerber), 160.0, 1e-14));
        assert!(close(eq(MeanStress::Soderberg), 225.0, 1e-14));
        let n = |c| safety_factor(&f, c, a, m).unwrap().unwrap();
        assert!(close(n(MeanStress::Uncorrected), se / a, 1e-15));
        let g = n(MeanStress::Goodman);
        assert!(close(g * a / se + g * m / su, 1.0, 1e-14), "{g}");
        let s = n(MeanStress::Soderberg);
        assert!(close(s * a / se + s * m / sy, 1.0, 1e-14), "{s}");
        let p = n(MeanStress::Gerber);
        assert!(close(p * a / se + (p * m / su).powi(2), 1.0, 1e-13), "{p}");
        // Hand values: Goodman 1/(0.5 + 0.25) = 4/3; Soderberg 1/(0.5 + 1/3) = 1.2.
        assert!(close(g, 4.0 / 3.0, 1e-14) && close(s, 1.2, 1e-14));
        // SWT's safety factor: the scaled cycle reaches √(σ_max σ_a) = S_e at n = S_e/√52500.
        let w = n(MeanStress::Swt);
        assert!(close(w * w * 52500.0, se * se, 1e-14) && close(w, se / 52500.0f64.sqrt(), 1e-15), "{w}");
        // A fully reversed cycle (σ_m = 0) is read at its own amplitude, as without a correction.
        assert!(close(equivalent_amplitude(&f, MeanStress::Swt, a, 0.0).unwrap().unwrap(), a, 1e-15));
        // A cycle whose σ_max is not tensile (σ_m ≤ −σ_a) does no damage under SWT.
        assert_eq!(equivalent_amplitude(&f, MeanStress::Swt, a, -a).unwrap(), Some(0.0));
        assert_eq!(equivalent_amplitude(&f, MeanStress::Swt, a, -2.0 * a).unwrap(), Some(0.0));
        assert_eq!(safety_factor(&f, MeanStress::Swt, a, -a).unwrap(), Some(f64::INFINITY));
        // Needs no strength.
        let bare = Fatigue { ultimate_strength: None, yield_strength: None, ..f.clone() };
        assert!(equivalent_amplitude(&bare, MeanStress::Swt, a, m).is_ok());
        // Gerber's parabola lies outside Goodman's line, Goodman outside Soderberg's.
        assert!(p > g && g > s);
        // A compressive mean is given no benefit.
        for c in [MeanStress::Goodman, MeanStress::Gerber, MeanStress::Soderberg] {
            assert_eq!(equivalent_amplitude(&f, c, a, -m).unwrap(), Some(a));
            assert!(close(safety_factor(&f, c, a, -m).unwrap().unwrap(), se / a, 1e-15));
        }
        // A mean at or past the strength: static failure.
        assert_eq!(equivalent_amplitude(&f, MeanStress::Goodman, a, su).unwrap(), None);
        assert_eq!(equivalent_amplitude(&f, MeanStress::Soderberg, a, 700.0).unwrap(), None);
        // A static mean with no amplitude: Goodman and Gerber reach S_u at n = S_u/σ_m.
        assert!(close(safety_factor(&f, MeanStress::Gerber, 0.0, 400.0).unwrap().unwrap(), 2.0, 1e-15));
        assert!(close(safety_factor(&f, MeanStress::Goodman, 0.0, 400.0).unwrap().unwrap(), 2.0, 1e-15));
        assert_eq!(safety_factor(&f, MeanStress::Goodman, 0.0, 0.0).unwrap(), Some(f64::INFINITY));
        // Missing strengths are named; no endurance limit, no safety factor.
        let bare = Fatigue { ultimate_strength: None, yield_strength: None, ..f.clone() };
        assert!(equivalent_amplitude(&bare, MeanStress::Goodman, a, m).unwrap_err().to_string().contains("ultimate strength"));
        assert!(equivalent_amplitude(&bare, MeanStress::Soderberg, a, m).unwrap_err().to_string().contains("yield strength"));
        assert_eq!(safety_factor(&basquin(None), MeanStress::Goodman, a, m).unwrap(), None);
    }

    /// The example history of ASTM E1049-85 (2017) for rainflow counting, as widely reproduced in
    /// the open literature (the standard itself was not read; these are those counts, each redone
    /// by hand with the three-point procedure): −2, 1, −3, 5,
    /// −1, 3, −4, 4, −2 counts as half a cycle of range 3, one and a half of 4, half of 6, one of 8
    /// and half of 9.
    #[test]
    fn rainflow_counts_the_astm_e1049_example() {
        let history = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0];
        let by_range = |cycles: &[Cycle]| {
            let mut out: Vec<(f64, f64)> = Vec::new();
            for c in cycles {
                match out.iter_mut().find(|(r, _)| *r == c.range) {
                    Some((_, n)) => *n += c.count,
                    None => out.push((c.range, c.count)),
                }
            }
            out.sort_by(|a, b| a.0.total_cmp(&b.0));
            out
        };
        let once = rainflow(&history);
        assert_eq!(by_range(&once), [(3.0, 0.5), (4.0, 1.5), (6.0, 0.5), (8.0, 1.0), (9.0, 0.5)]);
        // Each cycle's mean by hand: 3 from −2 to 1; 4 from 1 to −3 (half) and −1 to 3 (full); 8
        // from −3 to 5 and −4 to 4 (halves); 9 from 5 to −4; 6 from 4 to −2.
        let mut means: Vec<(f64, f64, f64)> = once.iter().map(|c| (c.range, c.mean, c.count)).collect();
        means.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        assert_eq!(means, [(3.0, -0.5, 0.5), (4.0, -1.0, 0.5), (4.0, 1.0, 1.0), (6.0, 1.0, 0.5), (8.0, 0.0, 0.5), (8.0, 1.0, 0.5), (9.0, 0.5, 0.5)]);
        // Eight reversals make four cycles.
        assert_eq!(once.iter().map(|c| c.count).sum::<f64>(), 4.0);
        // Repeated as a block (Downing & Socie's method for a repeating history, counted by hand from the largest
        // value, 5: 5, −1, 3, −4, 4, −2, 1, −3, 5): full cycles of 4 (−1 to 3), 3 (−2 to 1), 7
        // (4 to −3) and 9 (5 to −4).
        let repeated = rainflow_repeated(&history);
        assert!(repeated.iter().all(|c| c.count == 1.0), "{repeated:?}");
        assert_eq!(by_range(&repeated), [(3.0, 1.0), (4.0, 1.0), (7.0, 1.0), (9.0, 1.0)]);
    }

    /// Counting does not depend on the sign or the scale of the history (a node whose measure is
    /// negative sees the reversed history), drops non-reversals and plateaus, and closes a
    /// repeated block wherever it starts.
    #[test]
    fn rainflow_is_invariant_and_counts_every_reversal() {
        let history = [0.0, 2.0, 1.0, 1.0, 3.0, 3.5, -2.0, 0.5, -1.0, 4.0, 0.0];
        assert_eq!(reversals(&history), [0.0, 2.0, 1.0, 3.5, -2.0, 0.5, -1.0, 4.0, 0.0]);
        let once = rainflow(&history);
        let flipped: Vec<f64> = history.iter().map(|x| -2.5 * x).collect();
        let scaled = rainflow(&flipped);
        assert_eq!(once.len(), scaled.len());
        let mut a: Vec<(f64, f64)> = once.iter().map(|c| (2.5 * c.range, c.count)).collect();
        let mut b: Vec<(f64, f64)> = scaled.iter().map(|c| (c.range, c.count)).collect();
        a.sort_by(|x, y| x.0.total_cmp(&y.0));
        b.sort_by(|x, y| x.0.total_cmp(&y.0));
        for ((ra, na), (rb, nb)) in a.iter().zip(&b) {
            assert!(close(*ra, *rb, 1e-14) && na == nb);
        }
        // Once: half a cycle per reversal (9 points, 8 reversals: 4 cycles).
        assert_eq!(once.iter().map(|c| c.count).sum::<f64>(), 4.0);
        // Repeated: every cycle full, the largest range max − min, and a rotation of the block
        // counts the same.
        let r = rainflow_repeated(&history);
        assert!(r.iter().all(|c| c.count.fract() == 0.0));
        assert_eq!(r.first().map(|c| c.range), Some(6.0));
        let rotated: Vec<f64> = history[3..].iter().chain(&history[..3]).copied().collect();
        assert_eq!(rainflow_repeated(&rotated), r);
        // Nothing to count.
        assert!(rainflow(&[1.0, 1.0]).is_empty() && rainflow_repeated(&[2.0]).is_empty() && rainflow(&[]).is_empty());
        assert_eq!(rainflow(&[0.0, 1.0]), [Cycle { range: 1.0, mean: 0.5, count: 0.5 }]);
        assert_eq!(rainflow_repeated(&[0.0, 1.0]), [Cycle { range: 1.0, mean: 0.5, count: 1.0 }]);
    }

    /// Miner's sum of a two-block history against the hand calculation: n₁ cycles at ±1 and n₂ at
    /// ±0.5 times a reference stress, repeated, damage n₁/N₁ + n₂/N₂ per block.
    #[test]
    fn miner_sums_a_two_block_history() {
        let f = basquin(None);
        let (n1, n2) = (3, 5);
        let mut history = Vec::new();
        for _ in 0..n1 {
            history.extend([1.0, -1.0]);
        }
        for _ in 0..n2 {
            history.extend([0.5, -0.5]);
        }
        let cycles = rainflow_repeated(&history);
        assert_eq!(cycles, [Cycle { range: 2.0, mean: 0.0, count: 3.0 }, Cycle { range: 1.0, mean: 0.0, count: 5.0 }]);
        let m = 400.0;
        let p = assess(&f, MeanStress::Goodman, 1.0, m, &cycles, (-1.0, 1.0)).unwrap();
        let big = 0.5 * (400.0f64 / 1000.0).powf(-10.0);
        let small = 0.5 * (200.0f64 / 1000.0).powf(-10.0);
        let hand = 3.0 / big + 5.0 / small;
        assert!(close(p.damage, hand, 1e-12) && close(p.life, 1.0 / hand, 1e-12), "{p:?}");
        assert_eq!((p.amplitude, p.mean, p.equivalent), (400.0, 0.0, 400.0));
        // The stress factor scales every amplitude (not the mean).
        let k = assess(&f, MeanStress::Goodman, 1.5, m, &cycles, (-1.0, 1.0)).unwrap();
        let hand = 3.0 / life(&f, 600.0).cycles + 5.0 / life(&f, 300.0).cycles;
        assert!(close(k.damage, hand, 1e-12));
        // A constant amplitude cycle: the damage of one cycle is 1/N.
        let one = assess(&f, MeanStress::Uncorrected, 1.0, -m, &[Cycle { range: 2.0, mean: 0.0, count: 1.0 }], (-1.0, 1.0)).unwrap();
        assert!(close(one.life, big, 1e-12));
        assert_eq!(one.safety_factor, None);
        // No cycles: no damage.
        assert_eq!(assess(&f, MeanStress::Goodman, 1.0, m, &[], (-1.0, 1.0)).unwrap().life, f64::INFINITY);
    }

    /// The measures of simple stress states: uniaxial tension and compression keep their sign,
    /// pure shear counts as tensile, a biaxial state takes its larger principal stress.
    #[test]
    fn stress_measures_are_signed() {
        let tension = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let compression = tension.map(|x| -x);
        let shear = [0.0, 0.0, 0.0, 50.0, 0.0, 0.0];
        for m in [StressMeasure::SignedVonMises, StressMeasure::MaxAbsPrincipal] {
            assert!(close(measure(m, &tension), 100.0, 1e-14));
            assert!(close(measure(m, &compression), -100.0, 1e-14));
            assert!(measure(m, &shear) > 0.0);
        }
        assert!(close(measure(StressMeasure::SignedVonMises, &shear), 50.0 * 3f64.sqrt(), 1e-12));
        assert!(close(measure(StressMeasure::MaxAbsPrincipal, &shear), 50.0, 1e-12));
        let biaxial = [-120.0, 40.0, 0.0, 0.0, 0.0, 0.0];
        assert!(close(measure(StressMeasure::MaxAbsPrincipal, &biaxial), -120.0, 1e-14));
        assert!(close(measure(StressMeasure::SignedVonMises, &biaxial), -von_mises(&biaxial), 1e-14));
    }

    /// The sign of the signed measures near pure shear does not flip with noise: a shear of 200
    /// with a normal stress of ±1e-6 (or 0) is tensile on both sides of the tie, so the damage
    /// from zero load is the same on both (the life was 11.5 times longer for −1e-6 under
    /// Goodman, and SWT gave none), and a state clearly compressive (σ_xx = −2 against τ = 200,
    /// ratio 0.99) still is.
    #[test]
    fn the_sign_of_a_signed_measure_is_stable_near_pure_shear() {
        let f = basquin(None);
        let stress = |sxx: f64| [sxx, 0.0, 0.0, 200.0, 0.0, 0.0];
        let from_zero = [Cycle { range: 1.0, mean: 0.5, count: 1.0 }];
        let damage = |sxx: f64, m: StressMeasure, c: MeanStress| {
            let v = measure(m, &stress(sxx));
            (v, assess(&f, c, 1.0, v, &from_zero, (0.0, 1.0)).unwrap().damage)
        };
        for m in [StressMeasure::SignedVonMises, StressMeasure::MaxAbsPrincipal] {
            for c in [MeanStress::Goodman, MeanStress::Swt] {
                let reference = damage(0.0, m, c);
                assert!(reference.0 > 0.0 && reference.1 > 0.0, "{m:?} {c:?}: {reference:?}");
                for sxx in [1e-6, -1e-6, -0.05, 0.05] {
                    let got = damage(sxx, m, c);
                    assert!(got.0 > 0.0 && close(got.1, reference.1, 1e-2), "{m:?} {c:?} {sxx}: {got:?} against {reference:?}");
                }
                assert!(damage(-2.0, m, c).0 < 0.0, "{m:?} {c:?}: past the band the state is compressive");
            }
        }
        // An all-compressive state is compressive however close its two largest principal stresses.
        assert!(measure(StressMeasure::SignedVonMises, &[-200.0, -199.9, 0.0, 0.0, 0.0, 0.0]) < 0.0);
    }

    /// A static overload needs no cycle: a constant load whose mean is past the ultimate strength
    /// is a static failure (life 0, safety factor below 1), and so is a cycle that peaks there.
    /// The stress factor scales the mean as well as the amplitude.
    #[test]
    fn a_static_overload_and_a_peak_at_the_strength_are_failures() {
        let f = basquin(Some(300.0));
        // No cycles (max = min): a mean of 1000 against S_u = 800.
        let p = assess(&f, MeanStress::Goodman, 1.0, 1000.0, &[], (1.0, 1.0)).unwrap();
        assert!(p.static_failure && p.life == 0.0 && p.damage >= 2.0, "{p:?}");
        assert!(close(p.safety_factor.unwrap(), 0.8, 1e-15), "{p:?}");
        // Without an endurance limit the static factor is still reported; just under the strength is none.
        let p = assess(&basquin(None), MeanStress::Swt, 1.0, 800.0, &[], (1.0, 1.0)).unwrap();
        assert!(p.static_failure && p.safety_factor == Some(1.0), "{p:?}");
        assert!(!assess(&f, MeanStress::Goodman, 1.0, 799.0, &[], (1.0, 1.0)).unwrap().static_failure);
        // σ_a = 900, σ_m = 100 (peak 1000, past S_u = 800), whatever the correction.
        for corr in [MeanStress::Uncorrected, MeanStress::Goodman, MeanStress::Gerber, MeanStress::Soderberg, MeanStress::Swt] {
            let c = [Cycle { range: 2.0, mean: 1.0 / 9.0, count: 1.0 }];
            let p = assess(&f, corr, 1.0, 900.0, &c, (-1.0 + 1.0 / 9.0, 1.0 + 1.0 / 9.0)).unwrap();
            assert!(p.static_failure && p.life == 0.0, "{corr:?}: {p:?}");
        }
        // Peak below the strength: no static failure, a finite life.
        let c = [Cycle { range: 2.0, mean: 0.0, count: 1.0 }];
        let ok = assess(&f, MeanStress::Goodman, 1.0, 350.0, &c, (-1.0, 1.0)).unwrap();
        assert!(!ok.static_failure && ok.life > 0.0 && ok.life.is_finite(), "{ok:?}");
        // The shortest life is shown at half a cycle, not at minus infinity.
        let failed = assess(&f, MeanStress::Goodman, 1.0, 1000.0, &[], (1.0, 1.0)).unwrap();
        assert!(close(failed.log10_life(), 0.5f64.log10(), 1e-15));
        // The stress factor multiplies amplitude and mean alike: R = 0 at 200, factor 2, Goodman:
        // σ_a = σ_m = 200, σ_ar = 200 / (1 − 200/800) = 266.67 (amplitude only gave 228.6).
        let r0 = [Cycle { range: 1.0, mean: 0.5, count: 1.0 }];
        let k = assess(&f, MeanStress::Goodman, 2.0, 200.0, &r0, (0.0, 1.0)).unwrap();
        assert!((k.mean - 200.0).abs() < 1e-12 && (k.amplitude - 200.0).abs() < 1e-12, "{k:?}");
        assert!(close(k.equivalent, 200.0 / 0.75, 1e-14), "{k:?}");
    }

    /// Static failures and hostile numbers give a defined result, never a panic or NaN damage.
    #[test]
    fn static_failures_and_hostile_values_stay_defined() {
        let f = basquin(Some(300.0));
        let p = assess(&f, MeanStress::Goodman, 1.0, 900.0, &[Cycle { range: 1.0, mean: 1.0, count: 1.0 }], (0.5, 1.5)).unwrap();
        assert!(p.static_failure && p.damage == 2.0 && p.life == 0.0 && p.equivalent == f64::INFINITY, "{p:?}");
        let reversed = [Cycle { range: 2.0, mean: 0.0, count: 1.0 }];
        for m in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
            assert!(assess(&f, MeanStress::Uncorrected, 1.0, m, &reversed, (-1.0, 1.0)).is_err(), "{m}");
        }
        let huge = assess(&f, MeanStress::Uncorrected, 1.0, 1e300, &reversed, (-1.0, 1.0)).unwrap();
        assert!(huge.static_failure && huge.damage == 2.0, "{huge:?}");
        assert!(assess(&basquin(None), MeanStress::Soderberg, 1.0, 1.0, &[], (-1.0, 1.0)).is_ok());
        let bare = Fatigue { yield_strength: None, ..basquin(None) };
        assert!(assess(&bare, MeanStress::Soderberg, 1.0, 1.0, &[], (-1.0, 1.0)).is_err(), "checked even without cycles");
    }
}
