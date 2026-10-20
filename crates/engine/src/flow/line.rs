//! A fluid step's fields along a line, in SI: the values at points evenly spaced from one end to
//! the other, interpolated trilinearly between the cell centres around each point (from the
//! fluid cells among them), and their largest and smallest value located between the samples.
//! Within half a cell of a box face there is no cell centre beyond the point, so the value is the
//! nearest centre's: a temperature line run onto a held wall reads the first cell's temperature,
//! not the wall's (`wall_heat_flux` is what reads the wall).
//!
//! The heat flux into a box wall along a line on that wall comes from the wall's thermal
//! condition: a held temperature `T_w` gives `q = k (9 T₁ − T₂ − 8 T_w) / (3 Δx)`, the one-sided
//! second-order difference through the wall (the lattice holds it halfway between the cell
//! layers, on the box face) and the two cell-centre layers beside it, at `Δx/2` and `3Δx/2`; a
//! flux face gives its flux, an adiabatic face zero. Positive into the wall, as the walls' heat
//! histories.
//!
//! An extreme between samples is the vertex of the parabola through the extreme sample and its
//! two neighbours (three-point peak interpolation); at an end of the line, or beside a point
//! without a value, it is the sample itself.

use fieldcraft_lbm::{FluidResult, Frame};
use fieldcraft_model::{BoxFace, Directional, FaceCondition, Model, ThermalCondition};
use fieldcraft_render::Vec3;
use fieldcraft_render::flow::sample_with;
use serde::Deserialize;

use crate::results::magnitude;
use crate::{EngineError, Result};

/// Most points one line samples.
pub const MAX_POINTS: usize = 10_000;

/// A field a line samples.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LineField {
    VelocityX,
    VelocityY,
    VelocityZ,
    Speed,
    Pressure,
    Temperature,
    WallHeatFlux,
}

impl LineField {
    /// The reply's key for the values: the quantity they are (so the display converts them).
    pub fn key(self) -> &'static str {
        match self {
            LineField::VelocityX | LineField::VelocityY | LineField::VelocityZ | LineField::Speed => "velocity",
            LineField::Pressure => "pressure",
            LineField::Temperature => "temperature",
            LineField::WallHeatFlux => "heat_flux",
        }
    }
}

/// An extreme of a profile and where it is.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extreme {
    pub value: f64,
    pub point: [f64; 3],
}

/// A field along a line.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    pub points: Vec<[f64; 3]>,
    /// `None` where no fluid cell is near.
    pub values: Vec<Option<f64>>,
    pub max: Option<Extreme>,
    pub min: Option<Extreme>,
}

/// A box wall of a solved fluid step with heat: its condition and the fluid's conductivity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wall {
    pub face: BoxFace,
    pub condition: ThermalCondition,
    pub conductivity: f64,
}

impl Wall {
    /// Box face `face` of the domain of `r`'s step in `model` (the model `r` was solved from):
    /// an error unless it is a wall of a domain with heat, two cells or more across, in a fluid of
    /// isotropic conductivity.
    pub fn of(model: &Model, r: &FluidResult, face: BoxFace) -> Result<Wall> {
        let (_, domain) = fieldcraft_lbm::fluid_step(model, Some(&r.step))?;
        let thermal = domain.thermal.as_ref().ok_or_else(|| EngineError::Params("`wall_heat_flux` needs a fluid domain with heat".into()))?;
        if !matches!(domain.faces.get(face), FaceCondition::Wall {} | FaceCondition::MovingWall { .. } | FaceCondition::Slip {}) {
            return Err(EngineError::Params(format!("box face {} is not a wall: `wall_heat_flux` is the heat into a wall", face.name())));
        }
        if r.grid.cells.get(face.axis()).is_none_or(|n| *n < 2) {
            return Err(EngineError::Params(format!("box face {}: `wall_heat_flux` needs two cells across the box", face.name())));
        }
        let conductivity = match model.material(&domain.material).and_then(|m| m.thermal.as_ref()).and_then(|t| t.conductivity) {
            Some(Directional::Isotropic(k)) if k.is_finite() && k > 0.0 => k,
            _ => return Err(EngineError::Params("`wall_heat_flux` needs the fluid's isotropic conductivity".into())),
        };
        Ok(Wall { face, condition: *thermal.faces.get(face), conductivity })
    }
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [0, 1, 2].map(|k| a[k] + t * (b[k] - a[k]))
}

/// `count` points (1 to [`MAX_POINTS`]) evenly spaced from `from` to `to` (one point: `from`).
pub fn points(from: [f64; 3], to: [f64; 3], count: usize) -> Result<Vec<[f64; 3]>> {
    if !(1..=MAX_POINTS).contains(&count) {
        return Err(EngineError::Params(format!("a line samples 1 to {MAX_POINTS} points")));
    }
    let last = count.saturating_sub(1).max(1) as f64;
    Ok((0..count).map(|i| lerp(from, to, i as f64 / last)).collect())
}

fn frame(r: &FluidResult, frame: usize) -> Result<&Frame> {
    r.frames.get(frame).ok_or_else(|| EngineError::Params(format!("step `{}` has no frame {frame}", r.step)))
}

/// The frame's temperature (kelvin) in cell `c`, read only where sampled; an error when the frame
/// carries none.
fn temperatures<'a>(r: &'a FluidResult, f: &'a Frame) -> Result<impl Fn(usize) -> Option<f64> + 'a> {
    if f.temperature.len() != r.grid.len() || f.temperature.is_empty() {
        return Err(EngineError::Params(format!("step `{}` carries no heat: give its domain `thermal` settings", r.step)));
    }
    Ok(move |c: usize| f.temperature.get(c).map(|t| fieldcraft_model::units::kelvin_from(f64::from(*t), r.temperature_reference)))
}

/// The profile of `field` (any but [`LineField::WallHeatFlux`]) in frame `frame_index` of `r` at
/// `points`. Only the cells around the points are read.
pub fn profile(r: &FluidResult, frame_index: usize, field: LineField, points: Vec<[f64; 3]>) -> Result<Profile> {
    let f = frame(r, frame_index)?;
    let fluid = |c: usize| !r.solid.get(c).copied().unwrap_or(false);
    let velocity = |c: usize, read: fn([f64; 3]) -> f64| f.velocity.get(c).filter(|_| fluid(c)).map(|u| read(u.map(f64::from)));
    let value: Box<dyn Fn(usize) -> Option<f64> + '_> = match field {
        LineField::VelocityX => Box::new(move |c| velocity(c, |u| u[0])),
        LineField::VelocityY => Box::new(move |c| velocity(c, |u| u[1])),
        LineField::VelocityZ => Box::new(move |c| velocity(c, |u| u[2])),
        LineField::Speed => Box::new(move |c| velocity(c, magnitude)),
        LineField::Pressure => Box::new(move |c| f.pressure.get(c).filter(|_| fluid(c)).map(|p| f64::from(*p))),
        LineField::Temperature => {
            let t = temperatures(r, f)?;
            Box::new(move |c| t(c).filter(|_| fluid(c)))
        }
        LineField::WallHeatFlux => return Err(EngineError::Params("`wall_heat_flux` is read along its wall (`face`)".into())),
    };
    let grid = super::grid(r);
    let values = points.iter().map(|p| Some(sample_with(&grid, &value, Vec3::from(*p))).filter(|v| v.is_finite())).collect();
    Ok(Profile::new(points, values))
}

/// The heat per area into `wall` in frame `frame_index` of `r` at `points`, which must lie on the
/// wall.
pub fn wall_heat_flux(r: &FluidResult, frame_index: usize, wall: &Wall, points: Vec<[f64; 3]>) -> Result<Profile> {
    let h = r.grid.spacing;
    let axis = wall.face.axis();
    let (plane, inward) = if wall.face.is_max() { (r.grid.origin[axis] + h * r.grid.cells[axis] as f64, -1.0) } else { (r.grid.origin[axis], 1.0) };
    if points.iter().any(|p| (p[axis] - plane).abs() > 1e-6 * h) {
        return Err(EngineError::Params(format!("the line must lie on box face {} for `wall_heat_flux`", wall.face.name())));
    }
    let values = match wall.condition {
        ThermalCondition::Temperature { value: held } => {
            let f = frame(r, frame_index)?;
            let t = temperatures(r, f)?;
            let fluid = |c: usize| t(c).filter(|_| !r.solid.get(c).copied().unwrap_or(false));
            let grid = super::grid(r);
            // The temperature `d` cells into the fluid from the wall point `p`.
            let inside = |p: [f64; 3], d: f64| {
                let mut q = p;
                q[axis] = plane + inward * d * h;
                Some(sample_with(&grid, fluid, Vec3::from(q))).filter(|v| v.is_finite())
            };
            points.iter().map(|p| Some(wall.conductivity * (9.0 * inside(*p, 0.5)? - inside(*p, 1.5)? - 8.0 * held) / (3.0 * h))).collect()
        }
        ThermalCondition::Flux { value } => vec![Some(-value); points.len()],
        ThermalCondition::Adiabatic {} => vec![Some(0.0); points.len()],
    };
    Ok(Profile::new(points, values))
}

impl Profile {
    fn new(points: Vec<[f64; 3]>, values: Vec<Option<f64>>) -> Profile {
        // An extreme's point, between the samples at its fractional index.
        let at = |(k, value): (f64, f64)| {
            let i = k.floor() as usize;
            let point = match (points.get(i), points.get(i + 1)) {
                (Some(a), Some(b)) => lerp(*a, *b, k - i as f64),
                (Some(a), None) => *a,
                _ => [f64::NAN; 3],
            };
            Extreme { value, point }
        };
        Profile { max: extreme(&values, 1.0).map(at), min: extreme(&values, -1.0).map(at), points, values }
    }
}

/// The largest (`sign` 1) or smallest (`sign` −1) of `values`, as its fractional sample index
/// and its value: the vertex of the parabola through it and its neighbours when it is a peak
/// between two values, else the sample itself.
pub fn extreme(values: &[Option<f64>], sign: f64) -> Option<(f64, f64)> {
    let (k, best) = values.iter().enumerate().filter_map(|(k, v)| Some((k, (*v)?))).fold(None, |acc: Option<(usize, f64)>, (k, v)| match acc {
        Some((_, b)) if sign * v <= sign * b => acc,
        _ => Some((k, v)),
    })?;
    let neighbours = k.checked_sub(1).and_then(|i| values.get(i).copied().flatten()).zip(values.get(k + 1).copied().flatten());
    let Some((a, c)) = neighbours else { return Some((k as f64, best)) };
    let curvature = a - 2.0 * best + c;
    if sign * curvature >= 0.0 {
        // Flat: the sample itself.
        return Some((k as f64, best));
    }
    let d = (0.5 * (a - c) / curvature).clamp(-0.5, 0.5);
    Some((k as f64 + d, best - 0.25 * (a - c) * d))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A parabola sampled on the line is located exactly, between the samples, at both signs; an
    /// extreme at an end of the line or beside a gap is the sample itself.
    #[test]
    fn extremes_lie_on_the_parabola_through_three_samples() {
        let pts = points([0.0; 3], [1.0, 0.0, 0.0], 11).unwrap();
        assert_eq!((pts.len(), pts[10]), (11, [1.0, 0.0, 0.0]));
        assert_eq!(points([2.0; 3], [5.0; 3], 1).unwrap(), vec![[2.0; 3]]);
        assert!(points([0.0; 3], [1.0; 3], 0).is_err() && points([0.0; 3], [1.0; 3], MAX_POINTS + 1).is_err());
        let parabola = |x: f64| 3.0 - 8.0 * (x - 0.437) * (x - 0.437);
        let values: Vec<Option<f64>> = pts.iter().map(|p| Some(parabola(p[0]))).collect();
        let p = Profile::new(pts.clone(), values.clone());
        let max = p.max.unwrap();
        assert!((max.value - 3.0).abs() < 1e-12 && (max.point[0] - 0.437).abs() < 1e-12, "{max:?}");
        let negated: Vec<Option<f64>> = values.iter().map(|v| v.map(|v| -v)).collect();
        let min = Profile::new(pts.clone(), negated).min.unwrap();
        assert!((min.value + 3.0).abs() < 1e-12 && (min.point[0] - 0.437).abs() < 1e-12, "{min:?}");
        // The smallest of the parabola is at the end: no refinement.
        assert_eq!(p.min.unwrap().point, [1.0, 0.0, 0.0]);
        let mut gap = values.clone();
        gap[3] = None;
        assert_eq!(Profile::new(pts.clone(), gap).max.unwrap().point, pts[4]);
        assert_eq!(extreme(&[None, None], 1.0), None);
    }
}
