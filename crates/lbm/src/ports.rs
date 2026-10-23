//! Open-face transport from saved fields, for checking cooling and exchanger runs.
//!
//! Each cell contributes `|u · n| Δx²`. Mass uses the fluid's reference density (the
//! incompressible convention), not a density reconstructed from rounded gauge pressures.
//! Inflow and outflow are integrated separately, so recirculation never cancels a mixed
//! temperature. These are advective sensible-heat rates, not a complete energy balance:
//! conduction through the port, storage, wall heat and sources are separate terms.

use fieldcraft_model::BoxFace;

use crate::{Error, FluidResult, Result};

/// Properties at solve time. Optional in old saved results, which cannot report mass or heat
/// transport reliably from today's (possibly edited) material definition.
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
    /// Sensible-heat transport relative to the result's temperature reference, watts.
    pub sensible_heat: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Port {
    pub face: BoxFace,
    /// Fluid area only, excluding obstacle cells, square metres.
    pub area: f64,
    pub inflow: Stream,
    pub outflow: Stream,
}

#[derive(Default)]
struct Integral {
    flow: f64,
    pressure: f64,
    temperature: f64,
}

impl Integral {
    fn stream(&self, props: Option<PortProperties>, reference: f64, thermal: bool) -> Stream {
        let mass_flow = props.map(|p| p.density * self.flow);
        Stream {
            volume_flow: self.flow,
            mass_flow,
            pressure: (self.flow > 0.0).then(|| self.pressure / self.flow),
            temperature: (thermal && self.flow > 0.0).then(|| reference + self.temperature / self.flow),
            sensible_heat: if thermal { props.and_then(|p| p.specific_heat.map(|cp| p.density * cp * self.temperature)) } else { None },
        }
    }
}

impl FluidResult {
    /// Integrate the open faces of one saved frame. Uses fixed face/cell order on every backend
    /// and can be recomputed after reopening without rerunning the flow. Old results still give
    /// volume transport and mixed temperatures; absent material properties stay absent.
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
        let mut ports = Vec::new();
        for (face, open) in BoxFace::ALL.into_iter().zip(self.open_faces) {
            if !open {
                continue;
            }
            let axis = face.axis();
            let layer = if face.is_max() { self.grid.cells[axis] - 1 } else { 0 };
            let sign = if face.is_max() { 1.0 } else { -1.0 };
            let (mut incoming, mut outgoing) = (Integral::default(), Integral::default());
            let mut area = 0.0;
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
                    let velocity = f.velocity.get(c).and_then(|v| v.get(axis)).copied().map(f64::from);
                    let pressure = f.pressure.get(c).copied().map(f64::from);
                    let temperature = if thermal { f.temperature.get(c).copied().map(f64::from) } else { Some(0.0) };
                    let (Some(u), Some(p), Some(t)) = (velocity, pressure, temperature) else {
                        return Err(Error::Invalid("missing saved port cell".into()));
                    };
                    if !(u.is_finite() && p.is_finite() && t.is_finite()) {
                        return Err(Error::Invalid(format!("non-finite saved field on {}", face.name())));
                    }
                    area += cell_area;
                    let q = sign * u * cell_area;
                    let integral = if q < 0.0 { &mut incoming } else { &mut outgoing };
                    integral.flow += q.abs();
                    integral.pressure += q.abs() * p;
                    integral.temperature += q.abs() * t;
                }
            }
            let inflow = incoming.stream(self.port_properties, self.temperature_reference, thermal);
            let outflow = outgoing.stream(self.port_properties, self.temperature_reference, thermal);
            if ![area, inflow.volume_flow, outflow.volume_flow].iter().all(|v| v.is_finite())
                || [inflow, outflow]
                    .iter()
                    .any(|s| [s.mass_flow, s.pressure, s.temperature, s.sensible_heat].into_iter().flatten().any(|v| !v.is_finite()))
            {
                return Err(Error::Invalid(format!("port integration overflow on {}", face.name())));
            }
            ports.push(Port { face, area, inflow, outflow });
        }
        Ok(ports)
    }
}
