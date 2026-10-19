//! Fans and porous resistances inside a fluid domain: vents, grilles, filters and perforated
//! plates as a porous zone instead of every hole, and fans as a zone that raises the pressure by
//! their pressure–flow curve, so the operating point comes out of the solve.
//!
//! Both are axis-aligned boxes in the domain (the cells whose centres they contain). A
//! [`PorousZone`] resists the flow with the Darcy–Forchheimer law, a force per volume
//! `F = −(μ/K) u − (ρ C₂ / 2) |u| u`, so a plate of thickness `t` crossed at the superficial
//! velocity `u` drops the pressure by `Δp = (μ/K) u t + (ρ C₂/2) u² t`; a loss coefficient `k`
//! (`Δp = k ρ u²/2`) is `C₂ = k / t` over the zone's thickness along its axis. A [`Fan`] pushes
//! along its direction with a uniform force per volume `Δp / t`, `Δp` read from its curve at the
//! volume flow `Q` through it.
//!
//! Units are the model's consistent system, like everything else (permeability a length²,
//! `C₂` per length, `Q` a volume per time, `Δp` a pressure); nothing is converted.
//!
//! Sources: H. Darcy 1856 and P. Forchheimer 1901 through S. Ergun, *Fluid flow through packed
//! columns*, Chem. Eng. Prog. 48:89, 1952 (the quadratic resistance law); I. E. Idelchik,
//! *Handbook of Hydraulic Resistance*, 3rd ed., 1996 (loss coefficients of grids and screens);
//! the fan's operating point is the intersection of its curve with the system's resistance
//! curve (any fan engineering text, e.g. ASHRAE Handbook — Fundamentals, ch. 21).

use serde::{Deserialize, Serialize};

use crate::{Error, Result, finite, positive};

/// Most fans and porous zones together in one domain.
pub const MAX_ZONES: usize = 1000;
/// Most points on one fan curve.
pub const MAX_CURVE_POINTS: usize = 1000;
/// The fan's default under-relaxation of its pressure rise per time step.
pub const DEFAULT_FAN_RELAXATION: f64 = 0.01;

/// A coordinate axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Axis::X => "x",
            Axis::Y => "y",
            Axis::Z => "z",
        }
    }
}

/// How a porous zone resists the flow.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Resistance {
    /// `F = −(μ/K) u − (ρ C₂/2) |u| u`: the permeability `K` (a length²) and the inertial
    /// coefficient `C₂` (per length; 0 for a purely viscous medium).
    DarcyForchheimer {
        permeability: f64,
        #[serde(default)]
        inertial: f64,
    },
    /// `Δp = k ρ u²/2` across the zone along its axis (`C₂ = k / t`, no viscous part).
    LossCoefficient { coefficient: f64 },
}

/// A vent, grille, filter or perforated plate: a box of porous resistance.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PorousZone {
    pub name: String,
    pub min: [f64; 3],
    pub max: [f64; 3],
    /// The direction the flow goes through it: its thickness for a loss coefficient, and the
    /// direction its flow is reported in. The resistance itself acts in every direction.
    pub axis: Axis,
    pub resistance: Resistance,
}

/// The direction a fan blows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FanDirection {
    #[serde(rename = "+x")]
    XPlus,
    #[serde(rename = "-x")]
    XMinus,
    #[serde(rename = "+y")]
    YPlus,
    #[serde(rename = "-y")]
    YMinus,
    #[serde(rename = "+z")]
    ZPlus,
    #[serde(rename = "-z")]
    ZMinus,
}

impl FanDirection {
    /// The axis (0, 1, 2) and the sign along it.
    pub fn axis_sign(self) -> (usize, f64) {
        match self {
            FanDirection::XPlus => (0, 1.0),
            FanDirection::XMinus => (0, -1.0),
            FanDirection::YPlus => (1, 1.0),
            FanDirection::YMinus => (1, -1.0),
            FanDirection::ZPlus => (2, 1.0),
            FanDirection::ZMinus => (2, -1.0),
        }
    }
}

/// A fan: a box (usually one or a few cells thick) that raises the pressure along its direction
/// by its curve's `Δp(Q)`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fan {
    pub name: String,
    pub min: [f64; 3],
    pub max: [f64; 3],
    pub direction: FanDirection,
    /// `[Q, Δp]` points, `Q` (volume per time through the fan) strictly increasing from 0 or
    /// more, `Δp` not increasing; linear between them. Beyond the last point the fan holds its end
    /// value (it never brakes the flow past free delivery), before the first the first value.
    pub curve: Vec<[f64; 2]>,
    /// How far the pressure rise moves toward the curve's value each time step (0 < r ≤ 1);
    /// smaller is steadier.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relaxation: Option<f64>,
}

impl Fan {
    /// The curve's pressure rise at volume flow `q`.
    pub fn pressure_rise(&self, q: f64) -> f64 {
        curve_value(&self.curve, q)
    }

    /// The largest pressure rise (the shut-off pressure, the first point).
    pub fn max_pressure_rise(&self) -> f64 {
        self.curve.iter().map(|p| p[1]).fold(0.0, f64::max)
    }
}

/// `[x, y]` points' piecewise linear value at `x`, the end values beyond them.
pub fn curve_value(points: &[[f64; 2]], x: f64) -> f64 {
    let (Some(first), Some(last)) = (points.first(), points.last()) else { return 0.0 };
    if !(x > first[0]) {
        return first[1];
    }
    if x >= last[0] {
        return last[1];
    }
    for w in points.windows(2) {
        if let [a, b] = w
            && x <= b[0]
        {
            let s = (x - a[0]) / (b[0] - a[0]);
            return a[1] + s * (b[1] - a[1]);
        }
    }
    last[1]
}

fn check_box(min: [f64; 3], max: [f64; 3]) -> Result<()> {
    for a in 0..3 {
        let (lo, hi) = (min.get(a).copied().unwrap_or(0.0), max.get(a).copied().unwrap_or(0.0));
        finite("min", lo)?;
        finite("max", hi)?;
        if !(hi >= lo) {
            return Err(Error::Invalid(format!("the box's {} range {lo}..{hi} is empty", ["x", "y", "z"].get(a).copied().unwrap_or("?"))));
        }
    }
    Ok(())
}

impl PorousZone {
    pub fn validate(&self) -> Result<()> {
        check_box(self.min, self.max)?;
        match self.resistance {
            Resistance::DarcyForchheimer { permeability, inertial } => {
                positive("permeability", permeability)?;
                finite("inertial coefficient", inertial)?;
                if inertial < 0.0 {
                    return Err(Error::Invalid(format!("the inertial coefficient {inertial} is negative")));
                }
            }
            Resistance::LossCoefficient { coefficient } => {
                positive("loss coefficient", coefficient)?;
            }
        }
        Ok(())
    }
}

impl Fan {
    pub fn validate(&self) -> Result<()> {
        check_box(self.min, self.max)?;
        if self.curve.len() < 2 {
            return Err(Error::Invalid("a fan curve needs at least two [Q, Δp] points".into()));
        }
        if self.curve.len() > MAX_CURVE_POINTS {
            return Err(Error::Limit(format!("a fan curve holds at most {MAX_CURVE_POINTS} points")));
        }
        for (i, p) in self.curve.iter().enumerate() {
            finite("fan curve flow", p[0])?;
            finite("fan curve pressure", p[1])?;
            if p[0] < 0.0 {
                return Err(Error::Invalid(format!("fan curve point {i}: the flow {} is negative", p[0])));
            }
            if let Some(prev) = i.checked_sub(1).and_then(|j| self.curve.get(j)) {
                if !(p[0] > prev[0]) {
                    return Err(Error::Invalid(format!("fan curve point {i}: the flows must increase ({} after {})", p[0], prev[0])));
                }
                if p[1] > prev[1] {
                    return Err(Error::Invalid(format!(
                        "fan curve point {i}: the pressure rise must not increase with the flow ({} after {}); a rising branch has no stable operating point",
                        p[1], prev[1]
                    )));
                }
            }
        }
        if !(self.max_pressure_rise() > 0.0) {
            return Err(Error::Invalid("a fan curve needs a positive pressure rise somewhere".into()));
        }
        if let Some(r) = self.relaxation
            && !(r > 0.0 && r <= 1.0)
        {
            return Err(Error::Invalid(format!("the fan relaxation {r} must be in (0, 1]")));
        }
        Ok(())
    }
}

/// Check the zones of a domain: each valid, names unique and not empty, within the caps.
pub(crate) fn validate_zones(porous: &[PorousZone], fans: &[Fan]) -> Result<()> {
    if porous.len() + fans.len() > MAX_ZONES {
        return Err(Error::Limit(format!("a fluid domain holds at most {MAX_ZONES} fans and porous zones")));
    }
    let names: Vec<&str> = porous.iter().map(|z| z.name.as_str()).chain(fans.iter().map(|f| f.name.as_str())).collect();
    for (i, n) in names.iter().enumerate() {
        if n.trim().is_empty() {
            return Err(Error::Invalid("a fan or porous zone needs a name".into()));
        }
        if names.get(..i).is_some_and(|before| before.contains(n)) {
            return Err(Error::Invalid(format!("two fans or porous zones are called `{n}`")));
        }
    }
    for z in porous {
        z.validate().map_err(|e| crate::context("porous zone", &z.name, e))?;
    }
    for f in fans {
        f.validate().map_err(|e| crate::context("fan", &f.name, e))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fan() -> Fan {
        Fan {
            name: "F".into(),
            min: [0.0; 3],
            max: [1.0; 3],
            direction: FanDirection::XPlus,
            curve: vec![[0.0, 10.0], [5.0, 6.0], [10.0, 0.0]],
            relaxation: None,
        }
    }

    #[test]
    fn a_fan_curve_interpolates_and_holds_its_ends() {
        let f = fan();
        assert_eq!(f.pressure_rise(-1.0), 10.0);
        assert_eq!(f.pressure_rise(2.5), 8.0);
        assert_eq!(f.pressure_rise(7.5), 3.0);
        assert_eq!(f.pressure_rise(20.0), 0.0);
        assert_eq!(f.pressure_rise(f64::NAN), 10.0);
        assert!(f.validate().is_ok());
    }

    #[test]
    fn hostile_zones_are_refused() {
        let bad_fans = [
            Fan { curve: vec![[0.0, 1.0]], ..fan() },
            Fan { curve: vec![[0.0, 1.0], [0.0, 0.5]], ..fan() },
            Fan { curve: vec![[0.0, 1.0], [1.0, 2.0]], ..fan() },
            Fan { curve: vec![[-1.0, 1.0], [1.0, 0.0]], ..fan() },
            Fan { curve: vec![[0.0, f64::NAN], [1.0, 0.0]], ..fan() },
            Fan { curve: vec![[0.0, 0.0], [1.0, -1.0]], ..fan() },
            Fan { relaxation: Some(0.0), ..fan() },
            Fan { relaxation: Some(2.0), ..fan() },
            Fan { max: [-1.0, 1.0, 1.0], ..fan() },
            Fan { min: [f64::INFINITY, 0.0, 0.0], ..fan() },
        ];
        for f in &bad_fans {
            assert!(f.validate().is_err(), "{f:?}");
        }
        let zone = |r| PorousZone { name: "V".into(), min: [0.0; 3], max: [1.0; 3], axis: Axis::X, resistance: r };
        for r in [
            Resistance::DarcyForchheimer { permeability: 0.0, inertial: 0.0 },
            Resistance::DarcyForchheimer { permeability: 1.0, inertial: -1.0 },
            Resistance::DarcyForchheimer { permeability: f64::NAN, inertial: 0.0 },
            Resistance::LossCoefficient { coefficient: -2.0 },
        ] {
            assert!(zone(r).validate().is_err(), "{r:?}");
        }
        let ok = zone(Resistance::LossCoefficient { coefficient: 2.0 });
        assert!(validate_zones(std::slice::from_ref(&ok), &[]).is_ok());
        assert!(validate_zones(&[ok.clone(), ok.clone()], &[]).is_err());
        assert!(validate_zones(&[PorousZone { name: "F".into(), ..ok }], &[fan()]).is_err());
    }
}
