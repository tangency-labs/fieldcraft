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
//! **The law.** With the Spalart–Allmaras model (`wall_model: "spalart_allmaras_bounce"`) `u_τ`
//! comes from the model's own constant-stress law instead ([`crate::sa_law`]), and the flow sets
//! ν̃ = κ u_τ Δy at the interior node, the model's exact solution there; the scheme is otherwise
//! the same.
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
//! **Body force.** The wall plane passes through the node's centre, so the node's cell spans half a
//! cell into the fluid: the fluid between the wall (y = 0) and the first plane the links cross
//! (y = Δx/2) is that half cell. The steady balance of a fully developed channel (S. B. Pope,
//! *Turbulent Flows*, Cambridge University Press 2000, doi:10.1017/CBO9780511840531: the total
//! shear stress falls linearly from `τ_w` at the wall to zero at the centre, `τ(y) = τ_w (1 − y/δ)`,
//! with `τ_w = ρ a δ`) at `y = Δx/2` gives the stress handed to the fluid as `τ_w − ½ ρ a_t Δx`
//! (Han et al. 2021, doi:10.1088/1873-7005/ac1782, leave the body force in their momentum
//! equation's unstated `⋯` term). So the node hands the fluid `stress − body` along the wall, with
//! `body = ½ ρ a` (`ρ` the interior node's density, the one `τ_w = ρ u_τ²` uses; `a` the body
//! acceleration of the node's cell), and the wall receives `body` with the momentum the node takes
//! from the fluid: the stress along the wall, `2 Σ f_q + ½ ρ a_n` along the normal. The node's own
//! populations are no fluid density (two thirds of their sum are its rest and wall-parallel
//! populations, which never meet the fluid), hence the interior density. The split is the momentum
//! balance of the wall's half cell, so it is applied after the wall law's dispatch and holds for any
//! [`NodeLaw`]: Spalding's and the Spalart–Allmaras constant-stress law alike. Each step's share
//! enters that step's wall force, and a recorded wall force is the two-step mean of that per-step
//! total (`flux + ½ρa`), as for every other contribution (`solve`). Summed over the box the
//! walls carry the whole body force, `a (M + ½ Σ ρ)` with `M` the fluid's mass. This is the steady
//! balance: a transient's half-cell `∂(ρu)/∂t` is not represented, while per-step momentum is
//! conserved exactly. Node cells belong to no porous zone or fan, so only the body acceleration
//! enters; heat (buoyancy) is refused. In the laminar limit the profile is the Poiseuille parabola
//! shifted by the slip `a Δx²/(2ν)`: the law's stress `ν u_1 / Δx` equals the wall's `a δ`, so
//! `u_1 = a δ Δx / ν`, where the parabola gives `a (δ − ½Δx) Δx / ν`. The slip alone raises the bulk
//! speed by `6/(M(M+1))` with `M = 2N` the wall-to-wall cells: 8.3 %, 2.2 % and 0.57 % at 4, 8 and 16
//! cells per δ (second order; issue #1101; `the_laminar_on_node_profile_is_the_parabola_shifted_by_the_slip`).
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

/// The wall law the on-node wall takes its stress from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeLaw {
    /// Spalding's (κ 0.41, B 5.5), as Han et al.
    Spalding,
    /// The Spalart–Allmaras model's own constant-stress law ([`crate::sa_law`]).
    SpalartAllmaras,
}

impl NodeLaw {
    /// The friction velocity of a flow of speed `speed` parallel to the wall at distance `y`,
    /// kinematic viscosity `nu` (lattice units); 0 for input that is not usable.
    pub fn friction_velocity(self, speed: f64, y: f64, nu: f64) -> f64 {
        match self {
            NodeLaw::Spalding => spalding_friction(speed, y, nu, BOUNCE_KAPPA, BOUNCE_B),
            NodeLaw::SpalartAllmaras => crate::sa_law::friction_velocity(speed, y, nu),
        }
    }
}

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
/// the momentum along the velocity parallel to the wall the node keeps from the fluid per step (the
/// law's stress less the half cell's body force, see [`WallNode::body`]; lattice units: the paper's
/// `τ_{w,b}` is its negative). Returns the momentum the node takes from the fluid (populations
/// arriving from the fluid less those sent back).
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
    /// The momentum the wall took at the node in the last step: what the node took from the fluid
    /// plus `body`. Along the wall it is `stress` to round-off.
    pub taken: [f64; 3],
    /// The friction velocity u_τ of the last sample (lattice units; 0 when the stress is).
    pub friction: f64,
    /// The body force on the node's half cell from the last sample, `½ ρ a`: ρ the interior node's
    /// density, a the body acceleration of the node's cell. The node hands the fluid `stress − body`
    /// along the wall, and the wall receives `body` with what the node takes from the fluid.
    pub body: [f64; 3],
}

/// The boundary nodes of the on-node wall-function bounce, in cell order.
#[derive(Debug, Clone, PartialEq)]
pub struct WallNodes {
    nodes: Vec<WallNode>,
    /// The fluid's kinematic viscosity (lattice units).
    nu: f64,
    /// The wall law the stress comes from.
    law: NodeLaw,
}

impl WallNodes {
    /// `nodes` (sorted here by cell) for a fluid of kinematic viscosity `nu` (lattice units),
    /// their stress from `law`.
    pub fn new(mut nodes: Vec<WallNode>, nu: f64, law: NodeLaw) -> WallNodes {
        nodes.sort_by_key(|n| n.cell);
        WallNodes { nodes, nu, law }
    }

    /// The wall law the stress comes from.
    pub fn law(&self) -> NodeLaw {
        self.law
    }

    /// The restart state: none beyond the setup's (the law with it: the Spalart–Allmaras law runs
    /// only with the model, whose ν̃ the checkpoint holds). Each node's stress is sampled before
    /// its collision and the momentum it takes is added to the saved forces within the step.
    pub(crate) fn checkpoint(&mut self, a: &mut dyn crate::checkpoint::Archive) -> crate::Result<()> {
        use crate::checkpoint::fixed;
        fixed(a, self.nodes.len() as u64)?;
        for w in &self.nodes {
            let WallNode { cell, interior, axis, inward, group, stress: _, taken: _, friction: _, body: _ } = *w;
            for v in [u64::from(cell), u64::from(interior), u64::from(axis), inward as u64, u64::from(group)] {
                fixed(a, v)?;
            }
        }
        Ok(())
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
    /// (with the half force) of its interior node, and the wall law at distance 1, in parallel:
    /// each node's is its own. Zero at rest or with a number that is not usable. Each node's
    /// `body` is its half cell's body force, `½ ρ force(node cell)` with the same density, after the
    /// wall law's dispatch and whatever the law: zero when that number is not finite. `force(cell)`
    /// is the body acceleration of a cell.
    pub fn sample(&mut self, moments: impl Fn(usize) -> Option<(f64, [f64; 3])> + Sync, force: impl Fn(usize) -> [f64; 3] + Sync) {
        let (nu, law) = (self.nu, self.law);
        parallel::for_each_mut(&mut self.nodes, |_, w| {
            w.stress = [0.0; 3];
            w.friction = 0.0;
            w.body = [0.0; 3];
            let Some((rho, u)) = moments(w.interior as usize) else { return };
            let axis = usize::from(w.axis).min(2);
            let mut t = u;
            if let Some(v) = t.get_mut(axis) {
                *v = 0.0;
            }
            let speed = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            let ut = law.friction_velocity(speed, 1.0, nu);
            let scale = rho * ut * ut / speed;
            if speed > 0.0 && scale.is_finite() {
                w.stress = t.map(|v| scale * v);
                w.friction = ut;
            }
            let body = force(w.cell as usize).map(|a| 0.5 * rho * a);
            if body.iter().all(|v| v.is_finite()) {
                w.body = body;
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
