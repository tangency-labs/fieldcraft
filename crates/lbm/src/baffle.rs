//! Zero-thickness, full-span planar sheets on lattice faces.
//!
//! Halfway no-slip bounce-back uses Ladd (1994); the scalar boundary uses the
//! anti-bounce-back relation of the double-distribution scheme (Krüger et al.,
//! 2017, ch. 8) with Fourier's law through the sheet (Incropera et al. 2011,
//! §3.1); see docs/solvers/literature.md, C3.1 and C9. For the two outgoing
//! populations a and b, the incoming populations are a-H and b+H. With D3Q7
//! weight 1/8 the two surface temperatures are 8a-4H and 8b+4H. Fourier's law
//! H = G (T_low-T_high) therefore gives H = (a-b) G/(G+1/8).
//! This bounded mixing conserves heat pairwise, including at zero conductance,
//! without an explicit thermal time-step restriction or a fictitious solid cell.

use crate::flow::Group;
use crate::setup::Setup;
use crate::voxel::Grid;
use crate::{Error, Result};

/// Final area averages on a baffle's low/high coordinate sides, in SI.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct BaffleResult {
    pub name: String,
    pub area: f64,
    /// Surface temperatures, low side then high side, in kelvin.
    pub temperature: [f64; 2],
    /// Heat flux into the sheet from each fluid side, W/m²; equal and opposite.
    pub heat_flux: [f64; 2],
}

#[derive(Debug, Clone)]
pub(crate) struct Baffle {
    pub axis: usize,
    pub layer: usize,
    pub group: u32,
    mix: f64,
    zero: f64,
    scale: f64,
    pub result: BaffleResult,
}

impl Baffle {
    pub fn prepare(setup: &Setup, groups: &mut Vec<Group>) -> Result<Vec<Self>> {
        let Some(t) = &setup.thermal else { return Ok(Vec::new()) };
        let mut out = Vec::new();
        for b in setup.domain.thermal.iter().flat_map(|t| &t.baffles) {
            let axis = b.axis.index();
            let layer = setup.domain.face_layer(b.axis, b.position).ok().and_then(|l| usize::try_from(l).ok()).unwrap_or(0);
            if layer == 0 || layer >= setup.grid.cells.get(axis).copied().unwrap_or(0) {
                return Err(Error::Invalid(format!("baffle `{}` must lie on an interior lattice face", b.name)));
            }
            // G = conductance * dt/(rho cp dx), using the setup's power scale.
            let g = b.conductance * setup.grid.spacing * setup.grid.spacing / t.heat_scale;
            if !g.is_finite() || (b.conductance > 0.0 && g <= 0.0) {
                return Err(Error::Invalid(format!("baffle `{}` conductance does not fit this lattice", b.name)));
            }
            let group = groups.len() as u32;
            for side in ["low", "high"] {
                groups.push(Group { name: format!("baffle/{}/{side}", b.name), body: None });
            }
            let faces = setup.grid.cells[(axis + 1) % 3] * setup.grid.cells[(axis + 2) % 3];
            out.push(Self {
                axis,
                layer,
                group,
                mix: g / (g + 0.125),
                zero: t.zero,
                scale: t.heat_scale,
                result: BaffleResult {
                    name: b.name.clone(),
                    area: faces as f64 * setup.grid.spacing * setup.grid.spacing,
                    temperature: [t.initial_temperature + t.zero; 2],
                    heat_flux: [0.0; 2],
                },
            });
        }
        Ok(out)
    }

    /// The source coordinate is unwrapped: crossing a tangential periodic seam
    /// does not change which side of this sheet owns a diagonal population.
    pub fn group_for(&self, p: [usize; 3], c: [i32; 3]) -> Option<u32> {
        let dest = p[self.axis] as i64;
        let source = dest - i64::from(c[self.axis]);
        let plane = self.layer as i64;
        ((dest < plane) != (source < plane)).then_some(self.group + u32::from(dest >= plane))
    }

    /// Refuses a sheet with a solid cell on either side of it.
    pub fn validate_cells(&self, grid: &Grid, fluid: &[bool]) -> Result<()> {
        let mut ok = true;
        self.each_face(grid, |low, high| ok &= [low, high].iter().all(|&c| fluid.get(c).copied().unwrap_or(false)));
        if ok {
            Ok(())
        } else {
            Err(Error::Invalid(format!("baffle `{}` touches an obstacle: sheet/solid junctions are not supported yet", self.result.name)))
        }
    }

    /// Calls `f` with the cells either side of each face of the sheet, low side first.
    fn each_face(&self, grid: &Grid, mut f: impl FnMut(usize, usize)) {
        let (u, v) = ((self.axis + 1) % 3, (self.axis + 2) % 3);
        for j in 0..grid.cells[v] {
            for i in 0..grid.cells[u] {
                let mut p = [0; 3];
                p[self.axis] = self.layer - 1;
                p[u] = i;
                p[v] = j;
                let low = grid.index(p);
                p[self.axis] = self.layer;
                f(low, grid.index(p));
            }
        }
    }

    pub fn stream(&mut self, grid: &Grid, stride: usize, g: &[f64], next: &mut [f64], heat: &mut [f64]) {
        let q = 1 + 2 * self.axis;
        let mut sum = [0.0; 2];
        let mut total = 0.0;
        let mix = self.mix;
        self.each_face(grid, |cell, other| {
            let a = g.get(q * stride + cell).copied().unwrap_or(0.0);
            let b = g.get((q + 1) * stride + other).copied().unwrap_or(0.0);
            let h = (a - b) * mix;
            if let Some(v) = next.get_mut((q + 1) * stride + cell) {
                *v = a - h;
            }
            if let Some(v) = next.get_mut(q * stride + other) {
                *v = b + h;
            }
            sum[0] += 8.0 * a - 4.0 * h;
            sum[1] += 8.0 * b + 4.0 * h;
            total += h;
        });
        let faces = (grid.cells[(self.axis + 1) % 3] * grid.cells[(self.axis + 2) % 3]) as f64;
        self.result.temperature = sum.map(|s| s / faces + self.zero);
        let flux = total * self.scale / self.result.area;
        self.result.heat_flux = [flux, -flux];
        for (side, value) in [total, -total].into_iter().enumerate() {
            if let Some(v) = heat.get_mut(self.group as usize + side) {
                *v = value;
            }
        }
    }
}
