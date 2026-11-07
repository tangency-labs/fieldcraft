//! Open-face transport from saved fields, for checking cooling and exchanger runs.
//!
//! Each cell contributes `|u · n| Δx²`. Mass uses the fluid's reference density (the
//! incompressible convention), not a density reconstructed from rounded gauge pressures.
//! Inflow and outflow are integrated separately, so recirculation never cancels a mixed
//! temperature. These are advective sensible-heat rates, not a complete energy balance:
//! conduction through the port, storage, wall heat and sources are separate terms.
//!
//! Sensible heat is `ρ c Q (T − T_in)`, relative to the mixed temperature `T_in` of all the
//! flow entering the box in the frame: the net heat carried out, `Σ_out − Σ_in`, is then the
//! enthalpy the flow gains, whatever the lattice's temperature zero, also while the net volume
//! flow is not zero (a transient). Without inflow there is no reference, and no heat.

use fieldcraft_model::BoxFace;

use crate::{Error, FluidResult, HeatHistory, Result};

/// Properties at solve time, so mass and heat transport never depend on today's (possibly
/// edited) material definition.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PortProperties {
    pub density: f64,
    pub specific_heat: Option<f64>,
}

/// One direction through an open face. Volume/mass rates are nonnegative; sensible heat is
/// signed relative to its reference. No flow has no mixed value.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Stream {
    pub volume_flow: f64,
    pub mass_flow: Option<f64>,
    /// Mass-flow-weighted temperature, kelvin (volume weighted without stored properties).
    pub temperature: Option<f64>,
    /// Flow-weighted gauge static pressure, pascals; not total pressure.
    pub pressure: Option<f64>,
    /// Sensible-heat transport relative to the mixed inflow temperature ([`inlet_temperature`]),
    /// watts.
    pub sensible_heat: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Port {
    pub face: BoxFace,
    /// Fluid area only, excluding obstacle cells, square metres.
    pub area: f64,
    /// Area means over the fluid cells of the saved boundary layer, stagnant ones included;
    /// None for a blocked face.
    pub mean: Option<FaceMean>,
    pub inflow: Stream,
    pub outflow: Stream,
}

/// Area-weighted means over an open face's fluid cells: `Σ value Δx² / Σ Δx²`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceMean {
    /// Gauge static pressure, pascals.
    pub pressure: f64,
    /// Kelvin; None without a thermal field.
    pub temperature: Option<f64>,
    /// Metres per second.
    pub velocity: [f64; 3],
}

#[derive(Default)]
struct Integral {
    flow: f64,
    pressure: f64,
    temperature: f64,
}

impl Integral {
    /// The stream, its temperatures (with heat, `thermal`) above `reference` (kelvin) and its
    /// heat above the mixed inflow temperature `inlet` (also above `reference`), when there is
    /// inflow.
    fn stream(&self, props: Option<PortProperties>, reference: f64, thermal: bool, inlet: Option<f64>) -> Stream {
        let mass_flow = props.map(|p| p.density * self.flow);
        let heat = |p: PortProperties| Some(p.density * p.specific_heat? * (self.temperature - inlet? * self.flow));
        Stream {
            volume_flow: self.flow,
            mass_flow,
            pressure: (self.flow > 0.0).then(|| self.pressure / self.flow),
            temperature: (thermal && self.flow > 0.0).then(|| reference + self.temperature / self.flow),
            sensible_heat: props.and_then(heat),
        }
    }
}

/// The mixed temperature (kelvin) of all the flow entering the box through `ports`: the
/// flow-weighted mean of the inflows' temperatures; `None` without inflow or without heat.
pub fn inlet_temperature(ports: &[Port]) -> Option<f64> {
    let (flow, sum) = ports
        .iter()
        .filter_map(|p| Some((p.inflow.volume_flow, p.inflow.temperature?)))
        .fold((0.0, 0.0), |(q, s), (flow, t)| (q + flow, s + flow * t));
    (flow > 0.0).then(|| sum / flow)
}

impl FluidResult {
    /// Integrate the open faces of one saved frame. Uses fixed face/cell order on every backend
    /// and can be recomputed after reopening without rerunning the flow. Results without
    /// material properties still give volume transport and mixed temperatures.
    pub fn ports(&self, frame: usize) -> Result<Vec<Port>> {
        let f = self.frames.get(frame).ok_or_else(|| Error::Invalid(format!("no fluid frame {frame}")))?;
        let n = self.grid.cells.iter().try_fold(1usize, |n, side| n.checked_mul(*side));
        let Some(n) = n.filter(|n| *n > 0 && !self.grid.cells.contains(&0)) else {
            return Err(Error::Invalid("port integration needs a nonempty grid".into()));
        };
        let thermal = !f.temperature.is_empty();
        if self.solid.len() != n || f.velocity.len() != n || f.pressure.len() != n || (thermal && f.temperature.len() != n) {
            return Err(Error::Invalid("port integration needs complete saved cell fields".into()));
        }
        if !(self.grid.spacing > 0.0 && self.grid.spacing.is_finite() && self.temperature_reference.is_finite())
            || self
                .port_properties
                .is_some_and(|p| !(p.density > 0.0 && p.density.is_finite()) || p.specific_heat.is_some_and(|cp| !(cp > 0.0 && cp.is_finite())))
        {
            return Err(Error::Invalid("invalid saved port properties or cell spacing".into()));
        }
        let cell_area = self.grid.spacing.powi(2);
        let mut faces = Vec::new();
        for (face, open) in BoxFace::ALL.into_iter().zip(self.open_faces) {
            if !open {
                continue;
            }
            let axis = face.axis();
            let layer = if face.is_max() { self.grid.cells[axis] - 1 } else { 0 };
            let sign = if face.is_max() { 1.0 } else { -1.0 };
            let (mut incoming, mut outgoing) = (Integral::default(), Integral::default());
            let mut area = 0.0;
            let (mut pressure_sum, mut temperature_sum, mut velocity_sum) = (0.0, 0.0, [0.0; 3]);
            // Visit only the face, not the volume: this is also used by the results panel.
            let a = (axis + 1) % 3;
            let b = (axis + 2) % 3;
            for j in 0..self.grid.cells[b] {
                for i in 0..self.grid.cells[a] {
                    let mut ijk = [0; 3];
                    ijk[axis] = layer;
                    ijk[a] = i;
                    ijk[b] = j;
                    let c = self.grid.index(ijk);
                    if self.solid.get(c).copied().unwrap_or(true) {
                        continue;
                    }
                    let velocity = f.velocity.get(c).copied().map(|v| v.map(f64::from));
                    let pressure = f.pressure.get(c).copied().map(f64::from);
                    let temperature = if thermal { f.temperature.get(c).copied().map(f64::from) } else { Some(0.0) };
                    let (Some(u), Some(p), Some(t)) = (velocity, pressure, temperature) else {
                        return Err(Error::Invalid("missing saved port cell".into()));
                    };
                    if !(u.iter().all(|v| v.is_finite()) && p.is_finite() && t.is_finite()) {
                        return Err(Error::Invalid(format!("non-finite saved field on {}", face.name())));
                    }
                    area += cell_area;
                    pressure_sum += p * cell_area;
                    temperature_sum += t * cell_area;
                    for (sum, v) in velocity_sum.iter_mut().zip(u) {
                        *sum += v * cell_area;
                    }
                    let q = sign * u[axis] * cell_area;
                    let integral = if q < 0.0 { &mut incoming } else { &mut outgoing };
                    integral.flow += q.abs();
                    integral.pressure += q.abs() * p;
                    integral.temperature += q.abs() * t;
                }
            }
            let mean = (area > 0.0).then(|| FaceMean {
                pressure: pressure_sum / area,
                temperature: thermal.then(|| self.temperature_reference + temperature_sum / area),
                velocity: velocity_sum.map(|v| v / area),
            });
            if mean.is_some_and(|m| m.temperature.into_iter().chain(m.velocity).chain([m.pressure]).any(|v| !v.is_finite())) {
                return Err(Error::Invalid(format!("port average overflow on {}", face.name())));
            }
            faces.push((face, area, mean, incoming, outgoing));
        }
        // The mixed inflow temperature over every open face, above the reference.
        let (q_in, t_in) = faces.iter().fold((0.0, 0.0), |(q, t), (_, _, _, i, _)| (q + i.flow, t + i.temperature));
        let inlet = (thermal && q_in > 0.0).then(|| t_in / q_in);
        let mut ports = Vec::with_capacity(faces.len());
        for (face, area, mean, incoming, outgoing) in faces {
            let inflow = incoming.stream(self.port_properties, self.temperature_reference, thermal, inlet);
            let outflow = outgoing.stream(self.port_properties, self.temperature_reference, thermal, inlet);
            if ![area, inflow.volume_flow, outflow.volume_flow].iter().all(|v| v.is_finite())
                || [inflow, outflow]
                    .iter()
                    .any(|s| [s.mass_flow, s.pressure, s.temperature, s.sensible_heat].into_iter().flatten().any(|v| !v.is_finite()))
            {
                return Err(Error::Invalid(format!("port integration overflow on {}", face.name())));
            }
            ports.push(Port { face, area, mean, inflow, outflow });
        }
        Ok(ports)
    }

    /// The heat `wall` received (watts, positive into it) at `time`, linearly interpolated
    /// between the recorded history samples; None outside them (frame zero predates the
    /// history) or for a malformed history. Never extrapolated.
    pub fn wall_power(&self, wall: &HeatHistory, time: f64) -> Option<f64> {
        interpolate(&self.history_times, &wall.power, time)
    }
}

/// `values` at `time`, linearly between the strictly increasing `times`; None outside them or
/// for malformed input.
fn interpolate(times: &[f64], values: &[f64], time: f64) -> Option<f64> {
    if !time.is_finite()
        || times.len() != values.len()
        || times.is_empty()
        || times.iter().chain(values).any(|v| !v.is_finite())
        || times.windows(2).any(|p| p[0] >= p[1])
    {
        return None;
    }
    let k = times.partition_point(|t| *t < time);
    if times.get(k) == Some(&time) {
        return values.get(k).copied();
    }
    let lo = k.checked_sub(1)?;
    let (a, b) = (*times.get(lo)?, *times.get(k)?);
    let weight = (time - a) / (b - a);
    let value = (1.0 - weight) * values.get(lo)? + weight * values.get(k)?;
    value.is_finite().then_some(value)
}
