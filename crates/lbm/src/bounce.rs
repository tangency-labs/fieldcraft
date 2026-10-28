//! Han, Ooka and Kikumoto's wall-function bounce, on-node, as their eq. 14 states it
//! (`wall_model: "spalding_bounce"`): M. Han, R. Ooka, H. Kikumoto, *A wall function approach in
//! lattice Boltzmann method: algorithm and validation using turbulent channel flow*, Fluid Dyn.
//! Res. 53:045506, 2021 (doi:10.1088/1873-7005/ac1782), §2.1–2.3 and Appendix A.
//!
//! **The scheme.** A boundary node lies on the wall (their node `(i, 0, k)`). Its collision is
//! replaced by eq. 14; the streaming is unchanged. With the inward normal `n` (their `+y`):
//! - every axis population and every population parallel to the wall bounces back,
//!   `f*_q = f_q̄` (their eq. 14c: `f1* = f2`, `f2* = f1`, …, `f9* = f10`, `f15* = f16`, …);
//! - each diagonal pointing into the fluid takes the specular image of the population arriving
//!   from the fluid plus the wall stress (their eqs. 14a–b): with `t = c_q − n` its tangential
//!   part along axis b, `f*_q = f_{q'} + (Δt / 2Δy) τ_{w,b} t_b`, `c_{q'} = t − n`;
//! - each diagonal pointing into the wall takes the value of its opposite (their `f8* = f̃7*`,
//!   `f12* = f̃11*`, `f13* = f̃14*`, `f18* = f̃17*`); it streams out of the box and is dropped;
//! - the rest population is kept.
//!
//! `τ_{w,b} = −(u_b / |u_t|) τ_w` (their eq. 15) opposes the velocity parallel to the wall of the
//! interior node one cell above, `(i, 1, k)`, at distance `Δy`, and `τ_w = ρ u_τ²` with `u_τ`
//! from Spalding's law (their eqs. 16–17, κ 0.41, B 5.5) at that node's velocity and distance, ρ
//! its density. That velocity is the one the collision uses (with the half body force).
//!
//! **Index mapping.** Their Table A-1 numbers D3Q19 differently from [`crate::lattice`]; their
//! populations are mapped by velocity vector (the wall at `y = 0`, the fluid above): 7 (1,1,0) →
//! 7, 8 (−1,−1,0) → 8, 13 (1,−1,0) → 9, 14 (−1,1,0) → 10, 9 (1,0,1) → 11, 10 (−1,0,−1) → 12,
//! 15 (1,0,−1) → 13, 16 (−1,0,1) → 14, 11 (0,1,1) → 15, 12 (0,−1,−1) → 16, 17 (0,1,−1) → 17,
//! 18 (0,−1,1) → 18, 0–6 the same. So eq. 14a–b reads, in our numbering, `f*_7 = f_9 + τ_{w,x}/2`,
//! `f*_10 = f_8 − τ_{w,x}/2`, `f*_15 = f_18 + τ_{w,z}/2`, `f*_17 = f_16 − τ_{w,z}/2` (lattice
//! units, `Δt = Δy = 1`). The code works from the vectors for any of the six faces.
//!
//! **Grid.** The grid is cell-centred, the box faces half a cell outside the outermost centres.
//! At a box face with this wall the outermost cell layer is the boundary nodes: the wall plane
//! passes through their centres, half a cell inside the box face, and the first interior node is
//! one cell further, at `Δy` from the wall. A channel's wall-to-wall distance is therefore the
//! box height less one cell: 2δ on `2N + 1` cells across with `N` cells per δ.
//!
//! **Boundary nodes are not fluid.** They run eq. 14, not the collision, so they have no eddy
//! viscosity, no hybrid blend and no σ, take no part in the finite-difference gradients of their
//! neighbours (which are one-sided there, as at a box face), and are left out of the fluid's
//! cells, mass, mean velocity and fields (reported like solid cells). Their populations that
//! stream out through the box face are dropped. Their populations that would stream in from
//! outside the box are not needed: eq. 14 writes every population of the node from those
//! arriving from the fluid and from the neighbouring boundary nodes. Before eq. 14 each of them
//! is set to the specular image of the one arriving from the fluid (the free-slip ghost); it
//! reaches only the populations that leave the box, and with it the node's 19 populations sum to
//! the same before and after eq. 14, the mass conservation the paper states. The populations
//! parallel to the wall only ever bounce between boundary nodes: they never meet the fluid.
//!
//! **Wall force.** Each step the wall takes the momentum the boundary nodes take from the fluid:
//! `Σ c_q f_q` over the populations arriving from the fluid less `Σ c_q f*_q` over those eq. 14
//! sends back into it. Parallel to the wall that is the stress, `τ_w` per node and step (wall
//! area `Δx²`), to round-off, whatever the normal flow; along the normal it is the bounce-back's
//! `2 Σ f_q` of the populations arriving from the fluid (their pressure and normal momentum flux).
//!
//! **Where it runs.** On plane box faces at rest with the box periodic along the wall (the
//! paper's channel); obstacles, moving walls, a wall that meets another non-periodic face at an
//! edge, and heat are refused ("not supported yet"). The GPU backend runs the same scheme
//! (`crates/lbm-gpu/src/lbm.wgsl`, `node_bounce`).
//!
//! No code of any lattice Boltzmann package was read; this is written from the paper.

use fieldcraft_parallel as parallel;

use crate::lattice::{C, Q, opposite};
use crate::wall::{BOUNCE_B, BOUNCE_KAPPA, spalding_friction};

/// The direction whose velocity is `c_q` with its component along `axis` reversed (the specular
/// image across a plane normal to `axis`), per axis.
const MIRROR: [[u8; Q]; 3] = {
    let mut out = [[0u8; Q]; 3];
    let mut axis = 0;
    while axis < 3 {
        let mut q = 0;
        while q < Q {
            let mut m = C[q];
            m[axis] = -m[axis];
            let mut k = 0;
            while k < Q {
                if C[k][0] == m[0] && C[k][1] == m[1] && C[k][2] == m[2] {
                    out[axis][q] = k as u8;
                }
                k += 1;
            }
            q += 1;
        }
        axis += 1;
    }
    out
};

/// Han et al.'s eq. 14 on one boundary node: its streamed populations `fq` become its
/// post-collision populations. The wall's inward normal is `inward` (±1) along `axis`; `stress` is
/// the wall stress `τ_w` along the velocity parallel to the wall (the momentum the wall takes per
/// step, lattice units: the paper's `τ_{w,b}` is its negative). Returns the momentum the node
/// takes from the fluid (populations arriving from the fluid less those sent back).
pub fn on_node(fq: &mut [f64; Q], axis: usize, inward: i32, stress: [f64; 3]) -> [f64; 3] {
    let axis = axis.min(2);
    let mirror = MIRROR.get(axis).copied().unwrap_or([0; Q]);
    let normal = |q: usize| C.get(q).and_then(|c| c.get(axis)).copied().unwrap_or(0) * inward;
    let get = |f: &[f64; Q], q: usize| f.get(q).copied().unwrap_or(0.0);
    // The free-slip ghost: what would arrive from outside the box is the image of what arrives
    // from the fluid.
    let mut f = *fq;
    for q in 1..Q {
        if normal(q) > 0 {
            let image = usize::from(mirror.get(q).copied().unwrap_or(0));
            f[q] = get(fq, image);
        }
    }
    let mut out = [0.0; Q];
    out[0] = f[0];
    for q in 1..Q {
        let cq = C.get(q).copied().unwrap_or([0; 3]);
        match normal(q) {
            0 => out[q] = get(&f, opposite(q)),
            1 => {
                let image = usize::from(mirror.get(q).copied().unwrap_or(0));
                let mut v = get(&f, image);
                for (b, (cb, s)) in cq.iter().zip(stress).enumerate() {
                    if b != axis && *cb != 0 {
                        v -= 0.5 * f64::from(*cb) * s;
                    }
                }
                out[q] = v;
            }
            _ => {}
        }
    }
    for q in 1..Q {
        if normal(q) < 0 {
            out[q] = get(&out, opposite(q));
        }
    }
    let mut taken = [0.0; 3];
    for q in 1..Q {
        let cq = C.get(q).copied().unwrap_or([0; 3]);
        let v = match normal(q) {
            n if n < 0 => get(&f, q),
            n if n > 0 => -get(&out, q),
            _ => continue,
        };
        for (t, c) in taken.iter_mut().zip(cq) {
            *t += f64::from(c) * v;
        }
    }
    *fq = out;
    taken
}

/// A boundary node of the on-node wall-function bounce.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallNode {
    pub cell: u32,
    /// The interior node one cell into the fluid along the normal, which gives the stress.
    pub interior: u32,
    /// The wall's normal axis and the sign of its inward normal (+1 at a min face).
    pub axis: u8,
    pub inward: i8,
    /// The force group of its box face.
    pub group: u32,
    /// The stress from the last sample (lattice units: the momentum the wall takes per step,
    /// along the interior node's velocity parallel to the wall).
    pub stress: [f64; 3],
    /// The momentum the node took from the fluid in the last step.
    pub taken: [f64; 3],
}

/// The boundary nodes of the on-node wall-function bounce, in cell order.
#[derive(Debug, Clone, PartialEq)]
pub struct WallNodes {
    nodes: Vec<WallNode>,
    /// The fluid's kinematic viscosity (lattice units).
    nu: f64,
}

impl WallNodes {
    /// `nodes` (sorted here by cell) for a fluid of kinematic viscosity `nu` (lattice units).
    pub fn new(mut nodes: Vec<WallNode>, nu: f64) -> WallNodes {
        nodes.sort_by_key(|n| n.cell);
        WallNodes { nodes, nu }
    }

    /// The boundary nodes, in cell order.
    pub fn nodes(&self) -> &[WallNode] {
        &self.nodes
    }

    /// Mutable access, for the collision of the boundary nodes.
    pub(crate) fn nodes_mut(&mut self) -> &mut [WallNode] {
        &mut self.nodes
    }

    /// Before the collision: each node's stress from `moments(cell)`, the density and velocity
    /// (with the half force) of its interior node, Spalding's law (κ 0.41, B 5.5) at distance 1,
    /// in parallel: each node's is its own. Zero at rest or with a number that is not usable.
    pub fn sample(&mut self, moments: impl Fn(usize) -> Option<(f64, [f64; 3])> + Sync) {
        let nu = self.nu;
        parallel::for_each_mut(&mut self.nodes, |_, w| {
            w.stress = [0.0; 3];
            let Some((rho, u)) = moments(w.interior as usize) else { return };
            let axis = usize::from(w.axis).min(2);
            let mut t = u;
            if let Some(v) = t.get_mut(axis) {
                *v = 0.0;
            }
            let speed = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            let ut = spalding_friction(speed, 1.0, nu, BOUNCE_KAPPA, BOUNCE_B);
            let scale = rho * ut * ut / speed;
            if speed > 0.0 && scale.is_finite() {
                w.stress = t.map(|v| scale * v);
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_mirror_reverses_one_component() {
        for (axis, table) in MIRROR.iter().enumerate() {
            for q in 0..Q {
                let m = C[table[q] as usize];
                for a in 0..3 {
                    assert_eq!(m[a], if a == axis { -C[q][a] } else { C[q][a] });
                }
            }
        }
    }
}
