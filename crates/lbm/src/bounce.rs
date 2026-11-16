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

/// What one boundary node's neighbour along a lattice direction is, for [`roles`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Neighbour {
    /// A cell that collides (bulk fluid, or an inlet or outlet cell): it holds real populations.
    Fluid,
    /// Another boundary node.
    Node,
    /// A solid cell, or outside the box: nothing streams there.
    Wall,
}

/// What a boundary node sends along one direction `q` (the population `f*_q` it streams to its
/// neighbour along `c_q`): agent-hub notes/cfd-wall-obstacles.md §2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    /// Eq. 14a–c: the specular image across the wall of the population arriving from the fluid,
    /// with the stress along the tangential axis of a diagonal (`q` along the inward normal, or a
    /// diagonal of a pair whose two cells and the normal neighbour are all fluid).
    Han,
    /// Halfway bounce-back on a link to the fluid that is not one of Han's (an edge, a corner, the
    /// side of a plate's leading edge, or a pair that is not whole): `f*_q = f_q̄`, its momentum
    /// exchange counted as the node's edge impulse.
    Edge,
    /// Bounce-back on a link to another boundary node (the populations parallel to the wall).
    Between,
    /// Into a solid cell or out of the box: streams nowhere.
    Dropped,
}

/// The roles of a boundary node's directions. `normal` is its wall's normal axis and the sign of
/// its inward normal, when it has exactly one wall direction; `neighbour(q)` is what lies along
/// `c_q`. Han's eq. 14 runs on the normal direction and on each tangential axis `b` whose pair
/// `n ± e_b` reaches the fluid at both ends, provided the normal neighbour is fluid too: the image
/// Han's eq. 14 takes is a fluid population one cell from the wall. Everything else that reaches
/// the fluid bounces back ([`Role::Edge`]).
pub fn roles(normal: Option<(usize, i32)>, neighbour: impl Fn(usize) -> Neighbour) -> [Role; Q] {
    let mut out = [Role::Dropped; Q];
    for (q, r) in out.iter_mut().enumerate().skip(1) {
        *r = match neighbour(q) {
            Neighbour::Fluid => Role::Edge,
            Neighbour::Node => Role::Between,
            Neighbour::Wall => Role::Dropped,
        };
    }
    let Some((axis, inward)) = normal else { return out };
    let axis = axis.min(2);
    let dir = |v: [i32; 3]| C.iter().position(|c| *c == v);
    let mut n = [0; 3];
    n[axis] = inward.signum();
    let Some(qn) = dir(n).filter(|q| neighbour(*q) == Neighbour::Fluid) else { return out };
    if let Some(r) = out.get_mut(qn) {
        *r = Role::Han;
    }
    for b in (0..3).filter(|b| *b != axis) {
        let pair = [1, -1].map(|s| {
            let mut v = n;
            v[b] = s;
            dir(v)
        });
        if let [Some(a), Some(z)] = pair
            && neighbour(a) == Neighbour::Fluid
            && neighbour(z) == Neighbour::Fluid
        {
            out[a] = Role::Han;
            out[z] = Role::Han;
        }
    }
    out
}

/// The roles of a node on a plane wall whose tangential neighbours are all boundary nodes and
/// whose fluid side is whole: exactly what [`on_node`] does.
pub fn flat_roles(axis: usize, inward: i32) -> [Role; Q] {
    let axis = axis.min(2);
    std::array::from_fn(|q| {
        let n = C.get(q).and_then(|c| c.get(axis)).copied().unwrap_or(0) * inward.signum();
        match n {
            _ if q == 0 => Role::Dropped,
            1 => Role::Han,
            0 => Role::Between,
            _ => Role::Dropped,
        }
    })
}

/// Whether the roles carry the stress along tangential axis `b`: its pair is Han's.
pub fn carries(roles: &[Role; Q], axis: usize, inward: i32, b: usize) -> bool {
    if b == axis.min(2) || b > 2 {
        return false;
    }
    let mut v = [0; 3];
    v[axis.min(2)] = inward.signum();
    v[b] = 1;
    C.iter().position(|c| *c == v).is_some_and(|q| roles.get(q) == Some(&Role::Han))
}

/// Eq. 14 by link pairs on one boundary node of any shape (agent-hub notes/cfd-wall-obstacles.md
/// §2): its streamed populations `fq` become its post-collision populations, by `roles`. The wall's
/// normal is along `axis`; `stress` is the wall stress as in [`on_node`], and its
/// component along a tangential axis must be zero unless the roles carry it ([`carries`]).
///
/// Each role is a pair of links whose populations to the fluid sum to those from it, so the
/// fluid's mass is kept to round-off whatever the node's shape; a Han pair along `b` carries
/// `stress[b]` exactly, as on a plane wall. Returns the momentum the node took from the fluid and,
/// of it, the part of its [`Role::Edge`] links.
pub fn on_links(fq: &mut [f64; Q], roles: &[Role; Q], axis: usize, stress: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let axis = axis.min(2);
    let mirror = MIRROR.get(axis).copied().unwrap_or([0; Q]);
    let get = |f: &[f64; Q], q: usize| f.get(q).copied().unwrap_or(0.0);
    let f = *fq;
    let mut out = [0.0; Q];
    out[0] = f[0];
    for q in 1..Q {
        let cq = C.get(q).copied().unwrap_or([0; 3]);
        out[q] = match roles.get(q).copied().unwrap_or(Role::Dropped) {
            Role::Han => {
                let image = usize::from(mirror.get(q).copied().unwrap_or(0));
                let mut v = get(&f, image);
                for (b, (cb, s)) in cq.iter().zip(stress).enumerate() {
                    if b != axis && *cb != 0 {
                        v -= 0.5 * f64::from(*cb) * s;
                    }
                }
                v
            }
            Role::Edge | Role::Between | Role::Dropped => get(&f, opposite(q)),
        };
    }
    let (mut taken, mut edge) = ([0.0; 3], [0.0; 3]);
    for q in 1..Q {
        let role = roles.get(q).copied().unwrap_or(Role::Dropped);
        if !matches!(role, Role::Han | Role::Edge) {
            continue;
        }
        // From the fluid along c_q̄ = −c_q, back to it along c_q.
        let cq = C.get(q).copied().unwrap_or([0; 3]);
        let v = get(&f, opposite(q)) + get(&out, q);
        for a in 0..3 {
            let m = -f64::from(cq[a]) * v;
            taken[a] += m;
            if role == Role::Edge {
                edge[a] += m;
            }
        }
    }
    *fq = out;
    (taken, edge)
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
    /// What it does along each direction ([`roles`]). A node with [`flat_roles`] runs
    /// [`on_node`], one of any other shape [`on_links`].
    pub roles: [Role; Q],
    /// Whether `roles` are [`flat_roles`].
    pub flat: bool,
    /// The stress Spalding's law asked for at the last sample, before the axes the node can't
    /// carry were set to zero in `stress` (equal to it on a flat node).
    pub wanted: [f64; 3],
    /// Of `taken`, the momentum of its [`Role::Edge`] links in the last step (zero on a flat node).
    pub edge: [f64; 3],
}

impl WallNode {
    /// A node at `cell` whose wall's inward normal is `inward` (±1) along `axis`, sampled at
    /// `interior`, its force in `group`, with `roles`.
    pub fn new(cell: u32, interior: u32, axis: u8, inward: i8, group: u32, roles: [Role; Q]) -> WallNode {
        let flat = roles == flat_roles(usize::from(axis), i32::from(inward));
        WallNode { cell, interior, axis, inward, group, stress: [0.0; 3], taken: [0.0; 3], roles, flat, wanted: [0.0; 3], edge: [0.0; 3] }
    }

    /// Eq. 14 on this node's streamed populations `fq`: [`on_node`] when flat, [`on_links`]
    /// otherwise. Sets `taken` and `edge`.
    pub fn collide(&mut self, fq: &mut [f64; Q]) {
        let (axis, inward) = (usize::from(self.axis), i32::from(self.inward));
        if self.flat {
            self.taken = on_node(fq, axis, inward, self.stress);
            self.edge = [0.0; 3];
        } else {
            (self.taken, self.edge) = on_links(fq, &self.roles, axis, self.stress);
        }
    }
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

    /// The restart state: none beyond the setup's. Each node's stress is sampled before its
    /// collision and the momentum it takes is added to the saved forces within the step.
    pub(crate) fn checkpoint(&mut self, a: &mut dyn crate::checkpoint::Archive) -> crate::Result<()> {
        use crate::checkpoint::fixed;
        fixed(a, self.nodes.len() as u64)?;
        for w in &self.nodes {
            let WallNode { cell, interior, axis, inward, group, .. } = *w;
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
    /// (with the half force) of its interior node, Spalding's law (κ 0.41, B 5.5) at distance 1,
    /// in parallel: each node's is its own. Zero at rest or with a number that is not usable.
    pub fn sample(&mut self, moments: impl Fn(usize) -> Option<(f64, [f64; 3])> + Sync) {
        let nu = self.nu;
        parallel::for_each_mut(&mut self.nodes, |_, w| {
            w.stress = [0.0; 3];
            w.wanted = [0.0; 3];
            let axis = usize::from(w.axis).min(2);
            let inward = i32::from(w.inward);
            // A node with no Han link along its normal has no fluid node to sample.
            if !w.flat && !(0..3).any(|b| carries(&w.roles, axis, inward, b)) {
                return;
            }
            let Some((rho, u)) = moments(w.interior as usize) else { return };
            let mut t = u;
            if let Some(v) = t.get_mut(axis) {
                *v = 0.0;
            }
            let speed = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            let ut = spalding_friction(speed, 1.0, nu, BOUNCE_KAPPA, BOUNCE_B);
            let scale = rho * ut * ut / speed;
            if speed > 0.0 && scale.is_finite() {
                w.stress = t.map(|v| scale * v);
                w.wanted = w.stress;
                if !w.flat {
                    for (b, s) in w.stress.iter_mut().enumerate() {
                        if !carries(&w.roles, axis, inward, b) {
                            *s = 0.0;
                        }
                    }
                }
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn populations(seed: u64) -> [f64; Q] {
        let mut x = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
        std::array::from_fn(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            0.01 + 0.1 * (x >> 11) as f64 / (1u64 << 53) as f64
        })
    }

    fn normal_of(q: usize, axis: usize, inward: i32) -> i32 {
        C[q][axis] * inward
    }

    /// A plane wall's node (fluid on its side, nodes beside it, solid behind) gets the flat roles,
    /// and on them `on_links` writes what `on_node` writes on every population that reaches a
    /// neighbour, bit for bit.
    #[test]
    fn flat_roles_are_eq_14_bit_for_bit() {
        for axis in 0..3 {
            for inward in [1, -1] {
                let flat = flat_roles(axis, inward);
                let found = roles(Some((axis, inward)), |q| match normal_of(q, axis, inward) {
                    1 => Neighbour::Fluid,
                    0 => Neighbour::Node,
                    _ => Neighbour::Wall,
                });
                assert_eq!(found, flat, "axis {axis} {inward}");
                let mut s = [2e-4, -3e-4, 5e-5];
                s[axis] = 0.0;
                let f = populations(7 + axis as u64);
                let (mut a, mut b) = (f, f);
                let taken = on_node(&mut a, axis, inward, s);
                let (links, edge) = on_links(&mut b, &flat, axis, s);
                for q in (0..Q).filter(|q| *q == 0 || flat[*q] != Role::Dropped) {
                    assert_eq!(a[q].to_bits(), b[q].to_bits(), "axis {axis} {inward} q {q}");
                }
                assert_eq!(edge, [0.0; 3]);
                for k in 0..3 {
                    assert!((taken[k] - links[k]).abs() < 1e-15, "{taken:?} vs {links:?}");
                }
            }
        }
    }

    /// Whatever the neighbours, the populations a node sends to the fluid sum to those it took
    /// from it (mass to round-off), and each carried axis takes exactly its stress, the rest of
    /// the tangential momentum being the edge links'.
    #[test]
    fn every_shape_keeps_mass_and_carries_exactly_the_stress_of_its_pairs() {
        let mut seed = 1u64;
        for axis in 0..3 {
            for inward in [1, -1] {
                for pattern in 0..400u64 {
                    seed += 1;
                    let kinds = populations(seed ^ 0xabcd);
                    let neighbour = |q: usize| match (kinds[q] * 1e3) as u64 % 5 {
                        0 | 1 if pattern % 4 != 0 || normal_of(q, axis, inward) < 0 => Neighbour::Fluid,
                        2 => Neighbour::Node,
                        _ if normal_of(q, axis, inward) > 0 && pattern % 3 == 0 => Neighbour::Fluid,
                        _ => Neighbour::Wall,
                    };
                    let normal = (pattern % 7 != 0).then_some((axis, inward));
                    let r = roles(normal, neighbour);
                    let mut s = [3e-4, -2e-4, 1e-4];
                    for (b, v) in s.iter_mut().enumerate() {
                        if !carries(&r, axis, inward, b) {
                            *v = 0.0;
                        }
                    }
                    let f = populations(seed);
                    let mut out = f;
                    let (taken, edge) = on_links(&mut out, &r, axis, s);
                    let to_fluid: f64 = (1..Q).filter(|q| matches!(r[*q], Role::Han | Role::Edge)).map(|q| out[q]).sum();
                    let from_fluid: f64 = (1..Q).filter(|q| matches!(r[*q], Role::Han | Role::Edge)).map(|q| f[opposite(q)]).sum();
                    assert!((to_fluid - from_fluid).abs() < 1e-15, "{r:?}: {to_fluid} vs {from_fluid}");
                    for b in (0..3).filter(|b| *b != axis) {
                        let han = taken[b] - edge[b];
                        assert!((han - s[b]).abs() < 1e-15, "{r:?} axis {b}: {han} vs {}", s[b]);
                    }
                    // Roles never send Han's image along a direction that doesn't reach the fluid.
                    for q in 1..Q {
                        if r[q] == Role::Han {
                            assert_eq!(neighbour(q), Neighbour::Fluid);
                            assert!(normal.is_some() && normal_of(q, axis, inward) == 1);
                        }
                    }
                }
            }
        }
    }

    /// A plate's leading-edge node on the floor (z normal): fluid upstream at its own level and
    /// below, nodes downstream. Its fluid side is whole, so it carries both stresses; the
    /// upstream links bounce back as edge links.
    #[test]
    fn a_leading_edge_node_carries_its_stress_and_bounces_upstream() {
        let (axis, inward) = (2, 1);
        let neighbour = |q: usize| {
            let c = C[q];
            match (c[0], c[2]) {
                (_, 1) => Neighbour::Fluid,
                (-1, 0) | (-1, -1) => Neighbour::Fluid,
                (_, 0) => Neighbour::Node,
                _ => Neighbour::Wall,
            }
        };
        let r = roles(Some((axis, inward)), neighbour);
        assert!(carries(&r, axis, inward, 0) && carries(&r, axis, inward, 1));
        let dir = |v: [i32; 3]| C.iter().position(|c| *c == v).unwrap();
        assert_eq!(r[dir([-1, 0, 0])], Role::Edge);
        assert_eq!(r[dir([-1, 0, -1])], Role::Edge);
        assert_eq!(r[dir([1, 0, -1])], Role::Dropped);
        assert_eq!(r[dir([1, 0, 0])], Role::Between);
        // Without the fluid above, nothing is Han's and nothing is carried.
        let blocked = roles(Some((axis, inward)), |q| if C[q][2] == 1 { Neighbour::Node } else { neighbour(q) });
        assert!(blocked.iter().all(|r| *r != Role::Han));
        assert!((0..3).all(|b| !carries(&blocked, axis, inward, b)));
    }

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
