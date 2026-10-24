//! The lattice Boltzmann state and one time step: collide, stream, open boundaries.
//!
//! # Scheme
//! - **Lattice and equilibrium:** D3Q19 with the second-order equilibrium
//!   (Qian, d'Humières & Lallemand 1992; Krüger et al. 2017, ch. 3).
//! - **Collision:** two relaxation times (Ginzburg, Verhaeghe & d'Humières 2008): the symmetric
//!   part of each opposite pair relaxes with `1/τ⁺` (which sets the viscosity), the antisymmetric
//!   part with `1/τ⁻`, related by the parameter `Λ = (τ⁺ − ½)(τ⁻ − ½)`. With `Λ = 3/16`
//!   halfway bounce-back puts the wall exactly halfway between nodes and steady solutions do not
//!   depend on the viscosity. Or, per step, the recursive regularised collision
//!   ([`crate::regularized`]): every population rebuilt from the density, the velocity and the
//!   momentum flux before it relaxes, which stays stable as τ⁺ approaches ½.
//! - **Body force:** Guo, Zheng & Shi 2002 (PRE 65:046308): the velocity is
//!   `(Σ f c + F/2) / ρ` and the forcing term is split into its symmetric part, relaxed with
//!   `1 − 1/(2τ⁺)`, and its antisymmetric part, with `1 − 1/(2τ⁻)` (the TRT form of Guo's
//!   term, Krüger et al. 2017, ch. 10).
//! - **Walls:** halfway bounce-back at every link from a fluid cell to a solid cell or across a
//!   wall face of the box; moving walls add Ladd's term `6 w_q ρ₀ c_q·u_w` (Ladd 1994, JFM
//!   271:285). Slip faces reflect populations specularly.
//! - **Inlets and outlets:** Zou & He 1997 (Phys. Fluids 9:1591) on the first layer of cells: the
//!   unknown incoming populations are the opposite ones plus the equilibrium difference
//!   (bounce-back of the non-equilibrium part), the density (velocity inlet) or the normal velocity
//!   (pressure outlet) follow from mass and normal momentum, and the transverse momentum is
//!   corrected on the diagonal unknowns as M. Hecht and J. Harting do for D3Q19 (*Implementation
//!   of on-site velocity boundary conditions for D3Q19 lattice Boltzmann simulations*, J. Stat.
//!   Mech. 2010, P01018). The prescribed velocity is the one with the half force of Guo's scheme.
//!   With the regularised collision a pressure outlet then holds the equilibrium of its density
//!   and that normal velocity (the equilibrium scheme, Krüger et al. 2017, ch. 5), dropping the
//!   non-equilibrium part the bounce-back gave it: rebuilt from that part, the outlet's momentum
//!   flux feeds an oscillation that flips sign every step along the outlet and grows (measured in
//!   a channel at τ⁺ = 0.52: unstable within 150 steps), which the slowly relaxing antisymmetric
//!   part of the two-relaxation-time scheme damps. The cost is near the outlet: the viscous
//!   stress is not imposed on the outlet layer, so a fully developed channel profile is 3–8 % off
//!   in the cell next to the outlet, more as τ⁺ approaches ½ (3.2 % at τ⁺ = 0.6, 7–7.5 % at 0.52,
//!   9.6 % in the outlet cell itself; the bounce-back outlet, where it stays stable, 0.4 %), and
//!   unchanged at mid-channel. A convective outflow is the better outlet for this collision.
//!   A convective outlet ([`crate::outlet`]) updates the unknown populations by a convective
//!   outflow condition instead, holds each cell at the face's density, and falls back to Zou–He
//!   where the flow comes back in. A cell on the edge of two open faces belongs to one of them.
//! - **Forces:** momentum exchange on every wall link (Ladd 1994; R. Mei, D. Yu, W. Shyy, L.-S.
//!   Luo, *Force evaluation in the lattice Boltzmann method involving curved geometry*, PRE
//!   65:041203, 2002): a population `f*` leaving a fluid cell toward the wall and the one `f'`
//!   coming back transfer `c_out f* − c_in f'` to the wall.
//! - **Fans and porous zones:** a porous zone (a vent, grille, filter or perforated plate) adds
//!   the Darcy–Forchheimer force per volume `−ρ (α u + β |u| u)` (α = ν/K, β = C₂/2) to its
//!   cells through Guo's term; because it depends on the velocity it defines, the velocity comes
//!   from the exact root of that quadratic (Z. Guo, T. S. Zhao, *Lattice Boltzmann model for
//!   incompressible flows through porous media*, PRE 66:036304, 2002). A fan pushes its cells
//!   along its direction with `Δp / t` (t its thickness): before every step its volume flow `Q`
//!   (the mean over its cell layers of `Σ u_axis`) is read, and its pressure rise moves toward
//!   its curve's `Δp(Q)` by its relaxation factor, so the operating point is where the curve
//!   meets the system's resistance. Zones stay off the inlet and outlet layers.
//!
//! # Heat
//! With thermal settings a second lattice carries the temperature: D3Q7 (the rest direction and
//! the six axis directions, the first seven of D3Q19), weights ¼ and ⅛, `c_s² = ¼`, equilibrium
//! `w_q T (1 + 4 c_q·u)` with the velocity the momentum lattice just collided to (Krüger et al.
//! 2017, ch. 8: the passive-scalar double-distribution model). Its collision has two relaxation
//! times (I. Ginzburg, *Equilibrium-type and link-type lattice Boltzmann models for generic
//! advection and anisotropic-dispersion equation*, Adv. Water Resour. 28:1171, 2005): the
//! antisymmetric part of each axis pair, which carries the heat flux, relaxes with
//! `τ_T = 4 α + ½` (that sets the diffusivity), the symmetric part and the rest population with
//! `τ_T⁺ = ½ + Λ_T / (τ_T − ½)`. With `Λ_T` fixed ([`crate::setup::THERMAL_MAGIC`], 1/12) a
//! steady temperature field, its walls included, does not depend on `τ_T`, so neither do the
//! local wall heat fluxes (D. d'Humières, I. Ginzburg, *Viscosity independent numerical errors
//! for Lattice Boltzmann models: From recurrence equations to "magic" collision numbers*, Comput.
//! Math. Appl. 58:823, 2009); with a single relaxation time the anti-bounce-back wall moves with
//! `Λ_T = (τ_T − ½)²`, and the automatic time step's τ_T ≈ 1.4 put the de Vahl Davis cavity's
//! largest local Nusselt number 8 % high (#629). Below `τ_T = ½ + √Λ_T` (0.79) the lattice keeps a
//! single relaxation time, which is as stable as it was there and whose wall error is smaller
//! ([`crate::setup::thermal_magic`]). Its streaming follows the momentum lattice's links: a wall held at a temperature is
//! anti-bounce-back (`g_q = −g*_q̄ + 2 w_q T_w`, which puts the temperature halfway between the
//! cells), an adiabatic wall bounces back, a heat flux through a wall face adds the flux to the
//! bounced population (one axis direction crosses each wall face, so the cell gains exactly that
//! heat), an inlet holds its temperature (equilibrium at the cell's velocity) and an outlet takes
//! its neighbour's. Buoyancy is the Boussinesq force per mass `−β (T − T₀) g`, added to the body
//! force of every cell (Guo forcing as for the constant acceleration). The heat each wall face
//! receives is summed over its links like the momentum exchange, leaving out the links into
//! inlet and outlet cells: the open boundary sets all of their populations after streaming, so
//! that heat never stays in the fluid (counted, a held wall meeting a held inlet reported twice
//! its heat in a heated channel, #675).
//!
//! Obstacles that conduct heat carry their own temperature, conducted on their cells and coupled
//! to this lattice through their faces (conjugate heat transfer, [`crate::conduction`]).
//!
//! # Layout and determinism
//! Two arrays of populations, structure of arrays (`f[q n + cell]`), A–B pattern. Streaming
//! pulls: every population of a cell comes from one source (a neighbour, a periodic image, the
//! cell itself through a wall, or an open boundary). Cells whose 18 neighbours are all fluid
//! cells inside the box stream by one shifted copy per direction; the others follow per-link
//! rules. Collision and streaming run by fixed blocks of cells in parallel, each cell writing
//! only its own populations; the wall forces are summed in link order. There are no parallel
//! reductions, so the result is the same bit for bit on any number of threads and platform.

use fieldcraft_model::{BoxFace, Collision, FaceCondition, Outflow, Profile, ThermalCondition};

use fieldcraft_parallel as parallel;

use crate::conduction::{Interface, Layout, SolidHeat};
use crate::lattice::{C, Q, W, c, equilibrium, moving_wall, opposite};
use crate::outlet::{Action, Outlets, Populations};
use crate::setup::{MAX_LATTICE_SPEED, Obstacle, Setup, ZoneKind, ZoneSetup};
use crate::voxel::{Grid, TriangleTree, voxelise};
use crate::{Error, Result};

/// Where a population of a cell next to a wall or box face comes from (every other population
/// streams from the neighbour at `−c_q`). Public so a device backend can apply the same rules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Rule {
    /// The same direction from another cell (a neighbour or a periodic image).
    Pull { from: u32 },
    /// Bounce-back at a wall of force group `group`; `du` is the moving-wall term.
    Bounce { group: u32, du: f64 },
    /// Specular reflection at a slip face: direction `q` from cell `from`.
    Reflect { from: u32, q: u8, group: u32 },
    /// From outside an inlet or outlet: set by the open boundary after streaming.
    Open,
}

/// A population, direction `q` of `cell`, and where it comes from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Link {
    pub cell: u32,
    pub q: u8,
    pub rule: Rule,
}

/// An inlet or outlet cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenCell {
    pub cell: u32,
    pub axis: usize,
    /// +1 at a min face (the inward normal points up the axis), −1 at a max face.
    pub inward: f64,
    pub target: Target,
}

/// What an open cell holds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Target {
    /// The velocity (lattice units).
    Velocity([f64; 3]),
    /// The density.
    Density(f64),
    /// The density of a convective outlet: each cell is held at it while its velocity follows the
    /// flow leaving it ([`crate::outlet`]).
    Convective(f64),
}

/// A force group: a box face or a face of an obstacle.
#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    /// `x_min` … for box faces; the face name for obstacle faces.
    pub name: String,
    /// The obstacle body, for obstacle faces.
    pub body: Option<String>,
}

/// What a fan or porous zone does to its cells in the current step.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct ZoneAction {
    /// A force per volume (a fan's push, lattice units).
    push: [f64; 3],
    /// Porous resistance: the acceleration `−α u − β |u| u`.
    alpha: f64,
    beta: f64,
}

/// A zone's flow and pressure at one time (lattice units).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ZoneSample {
    /// The volume flow through it along its axis (a fan: along its direction), the mean over its
    /// cell layers.
    pub flow: f64,
    /// The static-pressure difference between the fluid cell layers just outside it (so it
    /// includes the wall friction over its thickness plus one cell): the drop
    /// `p_upstream − p_downstream` for a porous zone (upstream by the flow's sign), the rise
    /// `p_downstream − p_upstream` along a fan's direction.
    pub pressure: f64,
    /// A fan's imposed pressure rise this step (zero for a porous zone).
    pub fan_rise: f64,
}

/// The velocity (with the half force) and the force per volume of a cell from its density `rho`
/// and momentum `j`, with body acceleration `a` and its zone's action. Porous resistance depends
/// on the velocity it defines, `u = (j + F/2)/ρ`; solving that quadratic exactly keeps strong
/// resistances stable (Z. Guo, T. S. Zhao, *Lattice Boltzmann model for incompressible flows
/// through porous media*, PRE 66:036304, 2002): with `v = (j + (ρ a + push)/2)/ρ`,
/// `u = v / (c₀ + √(c₀² + c₁ |v|))`, `c₀ = (1 + α/2)/2`, `c₁ = β/2`.
#[inline(always)]
fn velocity_force(rho: f64, j: [f64; 3], a: [f64; 3], zone: Option<&ZoneAction>) -> ([f64; 3], [f64; 3]) {
    let mut force = [rho * a[0], rho * a[1], rho * a[2]];
    let Some(z) = zone else {
        return ([0, 1, 2].map(|k| (j[k] + 0.5 * force[k]) / rho), force);
    };
    for k in 0..3 {
        force[k] += z.push[k];
    }
    let v = [0, 1, 2].map(|k| (j[k] + 0.5 * force[k]) / rho);
    if z.alpha == 0.0 && z.beta == 0.0 {
        return (v, force);
    }
    let c0 = 0.5 * (1.0 + 0.5 * z.alpha);
    let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let scale = 1.0 / (c0 + (c0 * c0 + 0.5 * z.beta * speed).sqrt());
    let u = v.map(|c| c * scale);
    let us = speed * scale;
    for k in 0..3 {
        force[k] -= rho * (z.alpha + z.beta * us) * u[k];
    }
    (u, force)
}

/// The scalar lattice: directions 0..7 of D3Q19 (rest and the axes).
const QT: usize = 7;
/// Its weights.
const WT: [f64; QT] = [0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125];

/// The scalar equilibrium `w_q T (1 + c_q·u / c_s²)`, `c_s² = ¼`.
#[inline(always)]
fn equilibrium_t(q: usize, t: f64, u: [f64; 3]) -> f64 {
    let cq = c(q);
    WT.get(q).copied().unwrap_or(0.0) * t * (1.0 + 4.0 * (cq[0] * u[0] + cq[1] * u[1] + cq[2] * u[2]))
}

/// The temperature lattice and what it needs (see the module notes).
#[derive(Debug, Clone)]
pub struct ThermalState {
    /// `g[q stride + cell]` and the streaming target, A–B like the momentum lattice.
    g: Vec<f64>,
    g2: Vec<f64>,
    /// The temperature of every cell after the last step (the initial one in solid cells).
    temperature: Vec<f64>,
    /// The symmetric and antisymmetric relaxation rates `1/τ_T⁺` and `1/τ_T⁻`.
    omega_plus: f64,
    omega_minus: f64,
    initial: f64,
    /// `−β g` in lattice units; zero without buoyancy.
    buoyancy: [f64; 3],
    /// Each force group's thermal condition (fluxes in lattice units).
    groups: Vec<ThermalCondition>,
    /// Per open cell: the temperature it holds (an inlet), or `None` to take its neighbour's
    /// (an outlet).
    open: Vec<Option<f64>>,
    /// The open cells, sorted: the open boundary sets every scalar population of each after
    /// streaming, so a wall link into one gives no heat that stays in the fluid.
    reset: Vec<u32>,
    /// The heat each group received in the last step (lattice units: temperature times cell
    /// faces per step).
    heat: Vec<f64>,
    /// The obstacles that conduct heat; `None` without any.
    solid: Option<SolidHeat>,
    /// Per group: whether it is a face of a conducting obstacle, whose links the interface
    /// handles (empty without conducting obstacles).
    conjugate: Vec<bool>,
    /// Whether the solids' temperatures stayed finite in the last step.
    solid_finite: bool,
}

impl ThermalState {
    /// The scalar lattice's symmetric and antisymmetric relaxation rates `1/τ_T⁺` and `1/τ_T⁻`.
    pub fn omegas(&self) -> (f64, f64) {
        (self.omega_plus, self.omega_minus)
    }

    /// The initial temperature (and the reference of the buoyancy).
    pub fn initial(&self) -> f64 {
        self.initial
    }

    /// `−β g` in lattice units; zero without buoyancy.
    pub fn buoyancy(&self) -> [f64; 3] {
        self.buoyancy
    }

    /// Each force group's thermal condition (fluxes in lattice units).
    pub fn groups(&self) -> &[ThermalCondition] {
        &self.groups
    }

    /// Per open cell, the temperature it holds, or `None` to take its neighbour's.
    pub fn open(&self) -> &[Option<f64>] {
        &self.open
    }

    /// The conducting obstacles' cells and state, when there are any.
    pub fn solid(&self) -> Option<&SolidHeat> {
        self.solid.as_ref()
    }

    /// Whether the open boundary sets every scalar population of `cell` after streaming (an
    /// inlet holds its temperature, an outlet copies its neighbour): a wall link into it then
    /// gives the fluid no heat, and the wall heat leaves it out (#675).
    pub fn resets(&self, cell: u32) -> bool {
        self.reset.binary_search(&cell).is_ok()
    }
}

/// The flow state.
#[derive(Debug, Clone)]
pub struct Flow {
    grid: Grid,
    n: usize,
    /// The distance between two directions' arrays in `f` and `g` (at least `n`).
    stride: usize,
    fluid: Vec<bool>,
    f: Vec<f64>,
    g: Vec<f64>,
    /// The links that do not stream by the shifted copies, in cell order.
    links: Vec<Link>,
    open: Vec<OpenCell>,
    omega_plus: f64,
    omega_minus: f64,
    collision: Collision,
    acceleration: [f64; 3],
    /// The momentum transferred to each group in the last step (lattice units).
    forces: Vec<[f64; 3]>,
    groups: Vec<Group>,
    /// Per direction, the index offset of the neighbour a population streams from.
    offset: [isize; Q],
    fluid_cells: usize,
    steps_done: u64,
    /// Obstacles that cover no cell centre (smaller than the spacing, or outside the box).
    empty_obstacles: Vec<String>,
    /// The temperature lattice, with thermal settings.
    thermal: Option<ThermalState>,
    /// With buoyancy: every cell's body acceleration (the constant one plus the buoyant one),
    /// refreshed from the temperature before each step.
    cell_force: Vec<[f64; 3]>,
    /// With heat: every cell's velocity from the last collision, for the scalar equilibrium.
    cell_u: Vec<[f64; 3]>,
    /// The fans and porous zones (porous first), their fluid cells, and per cell its zone
    /// (`NONE` outside them; empty without zones).
    zones: Vec<ZoneSetup>,
    zone_cells: Vec<Vec<u32>>,
    zone_of: Vec<u32>,
    /// Each zone's action this step.
    actions: Vec<ZoneAction>,
    /// Each fan's pressure rise (lattice units; zero for porous zones), relaxed toward its curve.
    fan_rise: Vec<f64>,
    /// What the convective outlet cells keep between steps.
    outlets: Outlets,
    /// Cells per parallel block ([`block_cells`]).
    block: usize,
    /// Each open cell's populations after its boundary, in the step (kept for its capacity).
    open_scratch: Vec<[f64; Q]>,
}

const NONE: u32 = u32::MAX;

/// Cells per parallel block of the collision and the streaming: about a 256th of the grid, so a
/// machine with many threads still has several blocks per thread (4096-cell blocks left 45 for
/// the 184 000 cells of the `flow-past` sample, #643), but at least 512 and at most 4096 cells,
/// and a multiple of 64 so blocks start on cache lines. It depends on the grid alone, and every
/// cell's update reads and writes only its own populations (collision) or writes only its own
/// (streaming), with no sum across cells, so the result does not depend on the number of threads
/// (the largest speed, the one value taken across blocks, is a maximum).
fn block_cells(n: usize) -> usize {
    (n / 256).clamp(512, 4096) / 64 * 64
}

/// Open cells per parallel block of the open boundaries.
const OPEN_BLOCK: usize = 256;

/// The first `n` populations of each direction of an array with direction stride `stride`, cut
/// into blocks of `block` cells: block b holds the Q directions' slices of cells b·block.. .
fn cell_blocks(data: &mut [f64], stride: usize, n: usize, block: usize) -> Vec<Vec<&mut [f64]>> {
    let block = block.max(1);
    let mut blocks: Vec<Vec<&mut [f64]>> = (0..n.div_ceil(block)).map(|_| Vec::with_capacity(Q)).collect();
    for direction in data.chunks_mut(stride.max(1)).take(Q) {
        for (b, cells) in blocks.iter_mut().zip(direction.get_mut(..n).unwrap_or_default().chunks_mut(block)) {
            b.push(cells);
        }
    }
    blocks
}

/// A vector allocated without aborting when memory runs out.
fn filled<T: Clone>(len: usize, value: T, what: &str) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(len).map_err(|_| Error::Limit(format!("not enough memory for the {what} ({len} entries)")))?;
    v.resize(len, value);
    Ok(v)
}

impl Flow {
    /// Build the grid (voxelising the obstacles), the streaming rules and the initial state (at
    /// rest or the domain's initial velocity, density 1).
    pub fn new(setup: &Setup, obstacles: &[Obstacle]) -> Result<Flow> {
        let grid = setup.grid;
        let [nx, ny, nz] = grid.cells;
        let n = grid.len();
        if n == 0 || n >= NONE as usize {
            return Err(Error::Limit(format!("{n} cells cannot be addressed")));
        }
        // Each direction's array starts an odd number of cache lines (8 doubles) after the
        // previous one, so the 19 arrays never fall into the same cache sets (a power-of-two grid
        // would otherwise put them all in one).
        let lines = n.div_ceil(8);
        let stride = 8 * (lines + 1 - lines % 2);
        let q_n = stride.checked_mul(Q).ok_or_else(|| Error::Limit("the grid is too large".into()))?;

        // Solid cells, obstacle by obstacle (the first obstacle that covers a cell owns it).
        let mut fluid = filled(n, true, "cell flags")?;
        let mut owner = filled(n, NONE, "cell owners")?;
        let mut trees = Vec::with_capacity(setup.obstacles.len());
        let mut empty_obstacles = Vec::new();
        // Force groups: the six box faces, then every face of every obstacle.
        let mut groups: Vec<Group> = BoxFace::ALL.iter().map(|f| Group { name: f.name().to_string(), body: None }).collect();
        let mut base = Vec::with_capacity(setup.obstacles.len());
        for (k, faces) in setup.obstacles.iter().enumerate() {
            let o = obstacles.get(faces.index).ok_or_else(|| Error::NotFound(format!("obstacle {}", faces.index)))?;
            let mut inside = filled(n, false, "obstacle cells")?;
            if voxelise(&o.surface, &grid, &mut inside)? == 0 {
                empty_obstacles.push(o.body.clone());
            }
            for (c, &solid) in inside.iter().enumerate() {
                if solid && owner[c] == NONE {
                    owner[c] = k as u32;
                    fluid[c] = false;
                }
            }
            trees.push((TriangleTree::new(&o.surface), o.surface.tri_face()));
            base.push(groups.len() as u32);
            groups.extend(faces.names.iter().map(|name| Group { name: name.clone(), body: Some(o.body.clone()) }));
        }

        // Conducting cells: per grid cell its index among them (empty without conducting solids),
        // and per conducting cell its solid.
        let solids = setup.thermal.as_ref().map_or(&[][..], |t| &t.solids[..]);
        let mut solid_of = vec![None; setup.obstacles.len()];
        for (kind, s) in solids.iter().enumerate() {
            if let Some(slot) = solid_of.get_mut(s.obstacle) {
                *slot = Some(kind);
            }
        }
        let mut conducting = if solids.is_empty() { Vec::new() } else { filled(n, NONE, "conducting cells")? };
        let (mut solid_cells, mut solid_kind) = (Vec::new(), Vec::new());
        for (cell, slot) in conducting.iter_mut().enumerate() {
            if let Some(&Some(kind)) = owner.get(cell).and_then(|k| solid_of.get(*k as usize)) {
                *slot = solid_cells.len() as u32;
                solid_cells.push(cell as u32);
                solid_kind.push(kind);
            }
        }
        // A powered solid that covers no cell would lose its power without a word.
        if let Some(s) = solids.iter().enumerate().find(|(k, s)| s.power != 0.0 && !solid_kind.contains(k)).map(|(_, s)| s) {
            return Err(Error::Invalid(format!(
                "conducting solid `{}` generates heat but covers no cell centre: it is thinner than the spacing ({}) or outside the box; refine the spacing",
                s.body, grid.spacing
            )));
        }
        let conducting_at = |cell: usize| conducting.get(cell).copied().filter(|k| *k != NONE);
        let mut interface: Vec<Interface> = Vec::new();

        let domain = &setup.domain;
        let periodic = [domain.faces.periodic(0), domain.faces.periodic(1), domain.faces.periodic(2)];
        let wall_velocity = |face: BoxFace| match *domain.faces.get(face) {
            FaceCondition::MovingWall { velocity } => Some(setup.to_lattice(velocity)),
            FaceCondition::Wall {} => Some([0.0; 3]),
            _ => None,
        };
        let is_slip = |face: BoxFace| matches!(domain.faces.get(face), FaceCondition::Slip {});
        let dims = [nx as i64, ny as i64, nz as i64];
        let index = |p: [i64; 3]| grid.index([p[0] as usize, p[1] as usize, p[2] as usize]) as u32;
        let mut offset = [0isize; Q];
        for (q, o) in offset.iter_mut().enumerate() {
            let cq = C[q];
            *o = cq[0] as isize + nx as isize * (cq[1] as isize + ny as isize * cq[2] as isize);
        }

        // The group and moving-wall term of a link into fluid cell `cell` (direction q) from the
        // obstacle cell `solid`: the face of the triangle nearest the link's midpoint.
        let obstacle_link = |cell: usize, solid: u32, q: usize| -> Rule {
            let k = owner.get(solid as usize).copied().unwrap_or(NONE) as usize;
            let (Some((tree, tags)), Some(faces), Some(&first)) = (trees.get(k), setup.obstacles.get(k), base.get(k)) else {
                return Rule::Bounce { group: 0, du: 0.0 };
            };
            let (centre, cq) = (grid.centre(cell), c(q));
            let mid = [0, 1, 2].map(|a| centre[a] - 0.5 * cq[a] * grid.spacing);
            let tag = tree.nearest(mid).and_then(|t| tags.get(t as usize).copied()).unwrap_or(0);
            let u = faces.velocities.get(tag as usize).copied().unwrap_or([0.0; 3]);
            Rule::Bounce { group: first + tag, du: moving_wall(q, u) }
        };

        // The links that do not stream from the neighbour the shifted copies read.
        let mut links: Vec<Link> = Vec::new();
        let mut fluid_cells = 0;
        for cell in 0..n {
            if !fluid[cell] {
                continue;
            }
            fluid_cells += 1;
            let ijk = grid.coords(cell);
            let p = [ijk[0] as i64, ijk[1] as i64, ijk[2] as i64];
            for q in 1..Q {
                let cq = C[q];
                let mut s = [p[0] - i64::from(cq[0]), p[1] - i64::from(cq[1]), p[2] - i64::from(cq[2])];
                let mut wrapped = false;
                let mut crossed: [Option<BoxFace>; 3] = [None; 3];
                for a in 0..3 {
                    if s[a] < 0 || s[a] >= dims[a] {
                        if periodic[a] {
                            s[a] = s[a].rem_euclid(dims[a]);
                            wrapped = true;
                        } else {
                            crossed[a] = Some(BoxFace::of(a, s[a] >= dims[a]));
                        }
                    }
                }
                let crossed_faces = || crossed.iter().flatten().copied();
                let rule = match crossed_faces().next() {
                    None => {
                        let from = index(s);
                        match (fluid.get(from as usize).copied().unwrap_or(false), wrapped) {
                            (true, false) => continue,
                            (true, true) => Rule::Pull { from },
                            (false, _) => {
                                let rule = obstacle_link(cell, from, q);
                                if let (true, Some(solid), Rule::Bounce { group, .. }) = (q < QT, conducting_at(from as usize), rule) {
                                    interface.push(Interface { cell: cell as u32, q: q as u8, solid, group, link: links.len() as u32, weight: 0.0 });
                                }
                                rule
                            }
                        }
                    }
                    Some(_) => {
                        if let Some(face) = crossed_faces().find(|f| wall_velocity(*f).is_some()) {
                            Rule::Bounce { group: face.index() as u32, du: moving_wall(q, wall_velocity(face).unwrap_or([0.0; 3])) }
                        } else if let Some(face) = crossed_faces().find(|f| is_slip(*f)) {
                            // The mirror source: the same cell layer, the tangential components shifted.
                            let a = face.axis();
                            let mut m = s;
                            m[a] = p[a];
                            let mut mirrored = C[q];
                            mirrored[a] = -mirrored[a];
                            let q2 = C.iter().position(|v| *v == mirrored).unwrap_or(opposite(q));
                            let single = crossed_faces().count() == 1;
                            let inside = (0..3).all(|b| m[b] >= 0 && m[b] < dims[b]);
                            match inside.then(|| index(m)) {
                                Some(from) if single && fluid.get(from as usize).copied().unwrap_or(false) => {
                                    Rule::Reflect { from, q: q2 as u8, group: face.index() as u32 }
                                }
                                _ => Rule::Bounce { group: face.index() as u32, du: 0.0 },
                            }
                        } else {
                            Rule::Open
                        }
                    }
                };
                links.push(Link { cell: cell as u32, q: q as u8, rule });
            }
        }

        // Inlets and outlets: the fluid cells of each open face's layer. A cell on the edge of two
        // open faces (or the corner of three) belongs to one of them only, so one rule sets it and
        // the GPU's one thread per open cell never races: a velocity face before a pressure face,
        // otherwise the first in `BoxFace::ALL` order (x_min, x_max, y_min, …). Its populations
        // from outside the other faces keep the placeholder, the bounce-back of its own.
        let mut open = Vec::new();
        let mut taken = filled(n, false, "open cells")?;
        let velocity_first = |f: &BoxFace| !matches!(domain.faces.get(*f), FaceCondition::Velocity { .. });
        let mut faces = BoxFace::ALL;
        faces.sort_by_key(velocity_first);
        for face in faces {
            let a = face.axis();
            let layer = if face.is_max() { grid.cells[a] - 1 } else { 0 };
            let (t1, t2) = ((a + 1) % 3, (a + 2) % 3);
            for j in 0..grid.cells[t2] {
                for i in 0..grid.cells[t1] {
                    let mut ijk = [0usize; 3];
                    ijk[a] = layer;
                    ijk[t1] = i;
                    ijk[t2] = j;
                    let cell = grid.index(ijk);
                    if !fluid[cell] || taken.get(cell).copied().unwrap_or(true) {
                        continue;
                    }
                    let target = match *domain.faces.get(face) {
                        FaceCondition::Velocity { velocity, profile } => {
                            let mut scale = 1.0;
                            if profile == Profile::Parabolic {
                                for t in [t1, t2].into_iter().filter(|t| !periodic[*t]) {
                                    let s = (ijk[t] as f64 + 0.5) / grid.cells[t] as f64;
                                    scale *= 4.0 * s * (1.0 - s);
                                }
                            }
                            Target::Velocity(setup.to_lattice(velocity).map(|v| v * scale))
                        }
                        FaceCondition::Pressure { pressure, outflow: Outflow::Uniform } => Target::Density(setup.lattice_density(pressure)),
                        FaceCondition::Pressure { pressure, outflow: Outflow::Convective } => Target::Convective(setup.lattice_density(pressure)),
                        _ => break,
                    };
                    if let Some(t) = taken.get_mut(cell) {
                        *t = true;
                    }
                    open.push(OpenCell { cell: cell as u32, axis: a, inward: if face.is_max() { -1.0 } else { 1.0 }, target });
                }
            }
        }

        // Fans and porous zones: their fluid cells (the first zone that contains a cell owns it),
        // kept off the inlet and outlet layers, whose boundary rule knows only the body force.
        let mut zone_of = if setup.zones.is_empty() { Vec::new() } else { filled(n, NONE, "zone cells")? };
        let mut zone_cells = Vec::with_capacity(setup.zones.len());
        for (z, zone) in setup.zones.iter().enumerate() {
            let mut cells = Vec::new();
            for k in zone.lo[2]..=zone.hi[2] {
                for j in zone.lo[1]..=zone.hi[1] {
                    for i in zone.lo[0]..=zone.hi[0] {
                        let cell = grid.index([i, j, k]);
                        if fluid.get(cell).copied().unwrap_or(false) && zone_of.get(cell) == Some(&NONE) {
                            zone_of[cell] = z as u32;
                            cells.push(cell as u32);
                        }
                    }
                }
            }
            if cells.is_empty() {
                return Err(Error::Invalid(format!("`{}` has no fluid cell: it lies inside obstacles or another zone", zone.name)));
            }
            zone_cells.push(cells);
        }
        if let Some(oc) = open.iter().find(|oc| zone_of.get(oc.cell as usize).is_some_and(|z| *z != NONE)) {
            let name = zone_of.get(oc.cell as usize).and_then(|z| setup.zones.get(*z as usize)).map_or("", |z| z.name.as_str());
            return Err(Error::Invalid(format!(
                "`{name}` reaches the inlet or outlet layer of box face {}: keep fans and porous zones at least one cell inside it",
                BoxFace::of(oc.axis, oc.inward < 0.0).name()
            )));
        }
        let actions: Vec<ZoneAction> = setup
            .zones
            .iter()
            .map(|z| match z.kind {
                ZoneKind::Porous { alpha, beta } => ZoneAction { push: [0.0; 3], alpha, beta },
                ZoneKind::Fan { .. } => ZoneAction::default(),
            })
            .collect();

        // The initial state: equilibrium at density 1 and the initial velocity.
        let mut f = filled(q_n, 0.0, "populations")?;
        // The streaming targets are allocated on the first CPU step: a flow handed to a device
        // backend never needs them, and they are half the CPU state.
        let g = Vec::new();
        for cell in 0..n {
            for q in 0..Q {
                f[q * stride + cell] = if fluid[cell] { equilibrium(q, 1.0, setup.initial_velocity) } else { 0.0 };
            }
        }

        // Heat: the scalar lattice at the initial temperature, each group's condition (the
        // box faces, then every obstacle face), each open cell's temperature.
        let (thermal, cell_force, cell_u) = match &setup.thermal {
            None => (None, Vec::new(), Vec::new()),
            Some(ts) => {
                let mut group_conditions: Vec<ThermalCondition> = ts.faces.to_vec();
                for faces in &setup.obstacles {
                    group_conditions.extend(faces.thermal.iter().copied());
                }
                group_conditions.resize(groups.len(), ThermalCondition::Adiabatic {});
                let open_t: Vec<Option<f64>> = open
                    .iter()
                    .map(|oc| {
                        let face = BoxFace::of(oc.axis, oc.inward < 0.0);
                        match (oc.target, ts.faces.get(face.index()).copied().unwrap_or_default()) {
                            (_, ThermalCondition::Temperature { value }) => Some(value),
                            (Target::Velocity(_), _) => Some(ts.initial_temperature),
                            (Target::Density(_) | Target::Convective(_), _) => None,
                        }
                    })
                    .collect();
                let mut reset: Vec<u32> = open.iter().map(|oc| oc.cell).collect();
                reset.sort_unstable();
                let t_n = stride.checked_mul(QT).ok_or_else(|| Error::Limit("the grid is too large".into()))?;
                let mut tg = filled(t_n, 0.0, "temperature populations")?;
                let tg2 = Vec::new();
                for cell in 0..n {
                    if fluid[cell] {
                        for q in 0..QT {
                            tg[q * stride + cell] = equilibrium_t(q, ts.initial_temperature, setup.initial_velocity);
                        }
                    }
                }
                let solid = if solid_cells.is_empty() {
                    None
                } else {
                    let walls = BoxFace::ALL.map(|face| {
                        let wall = wall_velocity(face).is_some() || is_slip(face);
                        wall.then(|| ts.faces.get(face.index()).copied().unwrap_or_default())
                    });
                    let layout = Layout { grid: &grid, periodic, conducting: &conducting, walls, kf: ts.conductivity, contacts: &ts.contacts };
                    Some(SolidHeat::build(
                        &layout,
                        solids,
                        std::mem::take(&mut solid_cells),
                        &solid_kind,
                        std::mem::take(&mut interface),
                        ts.initial_temperature,
                    )?)
                };
                let mut conjugate = Vec::new();
                if solid.is_some() {
                    conjugate = vec![false; groups.len()];
                    for (k, faces) in setup.obstacles.iter().enumerate() {
                        let (Some(Some(_)), Some(&first)) = (solid_of.get(k), base.get(k)) else { continue };
                        for flag in conjugate.iter_mut().skip(first as usize).take(faces.names.len()) {
                            *flag = true;
                        }
                    }
                }
                let buoyant = ts.buoyancy != [0.0; 3];
                let state = ThermalState {
                    g: tg,
                    g2: tg2,
                    temperature: filled(n, ts.initial_temperature, "temperature")?,
                    omega_plus: 1.0 / ts.tau_plus(),
                    omega_minus: 1.0 / ts.tau,
                    initial: ts.initial_temperature,
                    buoyancy: ts.buoyancy,
                    groups: group_conditions,
                    open: open_t,
                    reset,
                    heat: vec![0.0; groups.len()],
                    solid,
                    conjugate,
                    solid_finite: true,
                };
                (
                    Some(state),
                    if buoyant { filled(n, setup.acceleration, "cell forces")? } else { Vec::new() },
                    filled(n, [0.0; 3], "cell velocities")?,
                )
            }
        };
        // The links and open cells grew by pushes: their spare capacity (up to as much again)
        // goes back before the steps.
        links.shrink_to_fit();
        open.shrink_to_fit();
        let mut flow = Flow {
            grid,
            n,
            stride,
            fluid,
            f,
            g,
            links,
            open,
            omega_plus: 1.0 / setup.tau_plus,
            omega_minus: 1.0 / setup.tau_minus,
            collision: setup.collision,
            acceleration: setup.acceleration,
            forces: vec![[0.0; 3]; groups.len()],
            groups,
            offset,
            fluid_cells,
            steps_done: 0,
            empty_obstacles,
            thermal,
            cell_force,
            cell_u,
            fan_rise: vec![0.0; setup.zones.len()],
            zones: setup.zones.clone(),
            zone_cells,
            zone_of,
            actions,
            outlets: Outlets::default(),
            block: block_cells(n),
            open_scratch: Vec::new(),
        };
        flow.outlets = Outlets::new(&flow.open, &View { flow: &flow, f: &flow.f });
        Ok(flow)
    }

    /// What the convective outlet cells keep between steps.
    pub fn outlets(&self) -> &Outlets {
        &self.outlets
    }

    /// The action of a cell's zone, if it is in one.
    #[inline(always)]
    fn action_at(&self, cell: usize) -> Option<&ZoneAction> {
        self.zone_of.get(cell).and_then(|z| self.actions.get(*z as usize))
    }

    /// The fans and porous zones (porous first).
    pub fn zones(&self) -> &[ZoneSetup] {
        &self.zones
    }

    /// Each zone's fluid cells.
    pub fn zone_cells(&self) -> &[Vec<u32>] {
        &self.zone_cells
    }

    /// The volume flow through zone `z` along its axis, signed along a fan's direction: the
    /// mean over its layers of the flux through each, `Σ u_axis / layers`.
    fn zone_flow(&self, z: usize) -> f64 {
        let (Some(zone), Some(cells)) = (self.zones.get(z), self.zone_cells.get(z)) else { return 0.0 };
        let axis = zone.axis.min(2);
        let sum: f64 = cells.iter().filter_map(|c| self.moments(*c as usize)).map(|(_, u)| u[axis]).sum();
        let sign = if let ZoneKind::Fan { sign, .. } = zone.kind { sign } else { 1.0 };
        sign * sum / zone.layers() as f64
    }

    /// Before a step: every fan's flow, its pressure rise moved toward its curve's value, and
    /// its push (the rise over its thickness, along its direction).
    fn refresh_fans(&mut self) {
        for z in 0..self.zones.len() {
            let Some(ZoneKind::Fan { sign, curve, relaxation }) = self.zones.get(z).map(|zone| &zone.kind) else { continue };
            let (sign, relaxation) = (*sign, *relaxation);
            let target = fieldcraft_model::fluid_zones::curve_value(curve, self.zone_flow(z));
            let (Some(rise), Some(zone)) = (self.fan_rise.get_mut(z), self.zones.get(z)) else { continue };
            *rise += relaxation * (target - *rise);
            let mut push = [0.0; 3];
            push[zone.axis.min(2)] = sign * *rise / zone.layers() as f64;
            if let Some(a) = self.actions.get_mut(z) {
                a.push = push;
            }
        }
    }

    /// Every zone's flow and pressure difference now (lattice units).
    pub fn zone_samples(&self) -> Vec<ZoneSample> {
        (0..self.zones.len()).map(|z| self.zone_sample(z)).collect()
    }

    fn zone_sample(&self, z: usize) -> ZoneSample {
        let Some(zone) = self.zones.get(z) else { return ZoneSample::default() };
        let flow = self.zone_flow(z);
        let axis = zone.axis.min(2);
        // The mean density of the layer `layer` along the axis over the zone's footprint.
        let layer_density = |layer: usize| {
            let (mut sum, mut count) = (0.0, 0usize);
            let (t1, t2) = ((axis + 1) % 3, (axis + 2) % 3);
            for b in zone.lo[t2]..=zone.hi[t2] {
                for a in zone.lo[t1]..=zone.hi[t1] {
                    let mut ijk = [0usize; 3];
                    ijk[axis] = layer;
                    ijk[t1] = a;
                    ijk[t2] = b;
                    if let Some((rho, _)) = self.moments(self.grid.index(ijk)) {
                        sum += rho;
                        count += 1;
                    }
                }
            }
            (count > 0).then(|| sum / count as f64)
        };
        // The layers just outside the zone (its own end layers at the box's ends).
        let before = layer_density(zone.lo[axis].saturating_sub(1)).or_else(|| layer_density(zone.lo[axis])).unwrap_or(1.0);
        let after =
            layer_density((zone.hi[axis] + 1).min(self.grid.cells[axis].saturating_sub(1))).or_else(|| layer_density(zone.hi[axis])).unwrap_or(1.0);
        let up_minus_down = crate::lattice::CS2 * (before - after);
        match zone.kind {
            ZoneKind::Fan { sign, .. } => {
                ZoneSample { flow, pressure: -sign * up_minus_down, fan_rise: self.fan_rise.get(z).copied().unwrap_or(0.0) }
            }
            ZoneKind::Porous { .. } => ZoneSample { flow, pressure: if flow < 0.0 { -up_minus_down } else { up_minus_down }, fan_rise: 0.0 },
        }
    }

    /// The body acceleration of a cell: the constant one, plus the buoyant one with heat.
    #[inline(always)]
    fn force_at(&self, cell: usize) -> [f64; 3] {
        self.cell_force.get(cell).copied().unwrap_or(self.acceleration)
    }

    /// The temperature of a cell (model units), or `None` without heat.
    pub fn temperature(&self, cell: usize) -> Option<f64> {
        self.thermal.as_ref().and_then(|t| t.temperature.get(cell).copied())
    }

    /// Whether the flow carries heat.
    pub fn has_heat(&self) -> bool {
        self.thermal.is_some()
    }

    /// The temperature lattice and its settings, with heat.
    pub fn thermal(&self) -> Option<&ThermalState> {
        self.thermal.as_ref()
    }

    /// Scalar population `q` (0..7) of a cell; zero without heat.
    pub fn thermal_population(&self, q: usize, cell: usize) -> f64 {
        self.thermal.as_ref().and_then(|t| t.g.get(q * self.stride + cell).copied()).unwrap_or(0.0)
    }

    /// The heat each group received in the last step (lattice units), empty without heat.
    pub fn heat(&self) -> &[f64] {
        self.thermal.as_ref().map_or(&[], |t| &t.heat)
    }

    /// Set every fluid cell's temperature (and its populations to the equilibrium at its
    /// velocity). Nothing without heat.
    pub fn set_temperature(&mut self, temperature: impl Fn(usize) -> f64) {
        let n = self.n;
        let s = self.stride;
        let velocities: Vec<[f64; 3]> = (0..n).map(|c| self.moments(c).map_or([0.0; 3], |m| m.1)).collect();
        let fluid = &self.fluid;
        let Some(t) = &mut self.thermal else { return };
        for cell in 0..n {
            if !fluid.get(cell).copied().unwrap_or(false) {
                continue;
            }
            let value = temperature(cell);
            if let Some(slot) = t.temperature.get_mut(cell) {
                *slot = value;
            }
            let u = velocities.get(cell).copied().unwrap_or([0.0; 3]);
            for q in 0..QT {
                if let Some(g) = t.g.get_mut(q * s + cell) {
                    *g = equilibrium_t(q, value, u);
                }
            }
        }
    }

    /// The bytes its per-cell, link and open-cell arrays hold (their capacities): what
    /// [`crate::Memory`] estimates, measured.
    pub fn heap_bytes(&self) -> u64 {
        fn bytes<T>(v: &Vec<T>) -> u64 {
            (v.capacity() * size_of::<T>()) as u64
        }
        let mut total = bytes(&self.f)
            + bytes(&self.g)
            + bytes(&self.fluid)
            + bytes(&self.links)
            + bytes(&self.open)
            + bytes(&self.cell_force)
            + bytes(&self.cell_u)
            + bytes(&self.zone_of)
            + bytes(&self.open_scratch)
            + self.zone_cells.iter().map(bytes).sum::<u64>()
            + self.outlets.heap_bytes();
        if let Some(t) = &self.thermal {
            total += bytes(&t.g) + bytes(&t.g2) + bytes(&t.temperature) + bytes(&t.open) + bytes(&t.reset);
        }
        total
    }

    pub(crate) fn checkpoint_heap_bytes(&mut self) -> Result<u64> {
        self.targets()?;
        let mut bytes = self.heap_bytes();
        bytes = bytes.saturating_add(
            (self.forces.capacity() * size_of::<[f64; 3]>()
                + self.actions.capacity() * size_of::<ZoneAction>()
                + self.fan_rise.capacity() * size_of::<f64>()) as u64,
        );
        if let Some(t) = &self.thermal {
            bytes = bytes.saturating_add((t.heat.capacity() * size_of::<f64>()) as u64);
            if let Some(s) = &t.solid {
                bytes = bytes.saturating_add(s.heap_bytes());
            }
        }
        Ok(bytes)
    }

    pub fn grid(&self) -> &Grid {
        &self.grid
    }

    pub fn fluid_cells(&self) -> usize {
        self.fluid_cells
    }

    pub fn is_fluid(&self, cell: usize) -> bool {
        self.fluid.get(cell).copied().unwrap_or(false)
    }

    pub fn groups(&self) -> &[Group] {
        &self.groups
    }

    /// The momentum each group received in the last step (lattice units: force per step).
    pub fn forces(&self) -> &[[f64; 3]] {
        &self.forces
    }

    /// Obstacles that cover no cell centre.
    pub fn empty_obstacles(&self) -> &[String] {
        &self.empty_obstacles
    }

    pub fn steps_done(&self) -> u64 {
        self.steps_done
    }

    /// The links that do not stream from the neighbour at `−c_q`, in cell order.
    pub fn links(&self) -> &[Link] {
        &self.links
    }

    /// The inlet and outlet cells.
    pub fn open_cells(&self) -> &[OpenCell] {
        &self.open
    }

    /// Per direction, the index offset of the neighbour a population streams from.
    pub fn offsets(&self) -> [isize; Q] {
        self.offset
    }

    /// `1/τ⁺` and `1/τ⁻`.
    pub fn omegas(&self) -> (f64, f64) {
        (self.omega_plus, self.omega_minus)
    }

    /// The collision operator.
    pub fn collision(&self) -> Collision {
        self.collision
    }

    /// The body acceleration in lattice units.
    pub fn acceleration(&self) -> [f64; 3] {
        self.acceleration
    }

    /// Population `q` of `cell` (zero in a solid cell).
    pub fn population(&self, q: usize, cell: usize) -> f64 {
        self.f.get(q * self.stride + cell).copied().unwrap_or(0.0)
    }

    /// Set every fluid cell to the equilibrium of `state(cell) = (density, velocity)` in lattice
    /// units (velocity as defined with the half force).
    pub fn set_equilibrium(&mut self, state: impl Fn(usize) -> (f64, [f64; 3])) {
        let n = self.n;
        let s = self.stride;
        for cell in 0..n {
            if !self.fluid[cell] {
                continue;
            }
            let (rho, u) = state(cell);
            let a = self.force_at(cell);
            let mut force = a.map(|c| rho * c);
            if let Some(z) = self.action_at(cell) {
                let speed = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                for k in 0..3 {
                    force[k] += z.push[k] - rho * (z.alpha + z.beta * speed) * u[k];
                }
            }
            let bare = [0, 1, 2].map(|k| u[k] - 0.5 * force[k] / rho);
            for q in 0..Q {
                self.f[q * s + cell] = equilibrium(q, rho, bare);
            }
        }
        self.outlets = Outlets::new(&self.open, &View { flow: self, f: &self.f });
    }

    /// The density and velocity of a cell (lattice units), or `None` for a solid cell.
    pub fn moments(&self, cell: usize) -> Option<(f64, [f64; 3])> {
        if !self.is_fluid(cell) {
            return None;
        }
        let s = self.stride;
        let (mut rho, mut j) = (0.0, [0.0; 3]);
        for q in 0..Q {
            let fq = self.f.get(q * s + cell).copied().unwrap_or(0.0);
            let cq = c(q);
            rho += fq;
            j[0] += fq * cq[0];
            j[1] += fq * cq[1];
            j[2] += fq * cq[2];
        }
        let (u, _) = velocity_force(rho, j, self.force_at(cell), self.action_at(cell));
        Some((rho, u))
    }

    /// The total mass of the fluid cells (lattice units: the sum of every population).
    pub fn mass(&self) -> f64 {
        let n = self.n;
        let s = self.stride;
        let mut m = 0.0;
        for cell in (0..n).filter(|c| self.fluid[*c]) {
            for q in 0..Q {
                m += self.f[q * s + cell];
            }
        }
        m
    }

    /// One time step. Fails with [`Error::Unstable`] when a density stops being a positive
    /// number or a lattice speed exceeds [`MAX_LATTICE_SPEED`].
    pub fn step(&mut self) -> Result<()> {
        self.targets()?;
        self.refresh_buoyancy();
        self.refresh_fans();
        let max_u2 = self.collide();
        if !(max_u2 <= MAX_LATTICE_SPEED * MAX_LATTICE_SPEED) {
            return Err(Error::Unstable(format!(
                "the flow went unstable at time step {}: the lattice speed reached {:.3} (limit {MAX_LATTICE_SPEED}) or a density stopped being a positive number; lower the time step, refine the spacing or check the boundary conditions",
                self.steps_done + 1,
                max_u2.sqrt()
            )));
        }
        self.stream();
        self.open_boundaries();
        std::mem::swap(&mut self.f, &mut self.g);
        if self.thermal.is_some() {
            self.collide_thermal();
            self.stream_thermal();
            self.open_boundaries_thermal();
            if let Some(t) = &mut self.thermal {
                std::mem::swap(&mut t.g, &mut t.g2);
            }
            if !(self.update_temperature() && self.conduct()) {
                return Err(Error::Unstable(format!(
                    "the heat went unstable at time step {}: a temperature stopped being a number; lower the time step, refine the spacing or check the heat fluxes",
                    self.steps_done + 1
                )));
            }
        }
        self.steps_done += 1;
        Ok(())
    }

    /// The streaming targets (the second arrays of populations), allocated before the first step.
    fn targets(&mut self) -> Result<()> {
        if self.g.len() < self.f.len() {
            self.g = filled(self.f.len(), 0.0, "populations")?;
        }
        if let Some(t) = &mut self.thermal
            && t.g2.len() < t.g.len()
        {
            t.g2 = filled(t.g.len(), 0.0, "temperature populations")?;
        }
        Ok(())
    }

    /// With buoyancy: every cell's acceleration from its temperature.
    fn refresh_buoyancy(&mut self) {
        let Some(t) = &self.thermal else { return };
        if self.cell_force.is_empty() {
            return;
        }
        let (a, b, t0) = (self.acceleration, t.buoyancy, t.initial);
        for (force, temp) in self.cell_force.iter_mut().zip(&t.temperature) {
            let dt = temp - t0;
            *force = [a[0] + b[0] * dt, a[1] + b[1] * dt, a[2] + b[2] * dt];
        }
    }

    /// Two-relaxation-time collision of the scalar lattice, in place, with the velocity of the
    /// momentum collision: each axis pair's symmetric part relaxes with `1/τ_T⁺`, its
    /// antisymmetric part (the heat flux) with `1/τ_T⁻`, the rest population with `1/τ_T⁺`.
    fn collide_thermal(&mut self) {
        let block = self.block;
        let n = self.n;
        let s = self.stride;
        let fluid = &self.fluid;
        let cell_u = &self.cell_u;
        let Some(t) = &mut self.thermal else { return };
        let (wp, wm) = (t.omega_plus, t.omega_minus);
        let mut blocks = cell_blocks(&mut t.g, s, n, block);
        parallel::for_each_mut(&mut blocks, |b, g| {
            let len = g.first().map_or(0, |d| d.len());
            for i in 0..len {
                let cell = b * block + i;
                if !fluid.get(cell).copied().unwrap_or(false) {
                    continue;
                }
                let mut temp = 0.0;
                for q in 0..QT {
                    temp += g.get(q).and_then(|d| d.get(i)).copied().unwrap_or(0.0);
                }
                let u = cell_u.get(cell).copied().unwrap_or([0.0; 3]);
                if let Some(x) = g.get_mut(0).and_then(|d| d.get_mut(i)) {
                    *x -= wp * (*x - equilibrium_t(0, temp, u));
                }
                // the pairs (1, 2), (3, 4), (5, 6)
                for q in [1, 3, 5] {
                    let (a, b) =
                        (g.get(q).and_then(|d| d.get(i)).copied().unwrap_or(0.0), g.get(q + 1).and_then(|d| d.get(i)).copied().unwrap_or(0.0));
                    let (ea, eb) = (equilibrium_t(q, temp, u), equilibrium_t(q + 1, temp, u));
                    let even = wp * 0.5 * (a + b - ea - eb);
                    let odd = wm * 0.5 * (a - b - ea + eb);
                    if let Some(x) = g.get_mut(q).and_then(|d| d.get_mut(i)) {
                        *x = a - even - odd;
                    }
                    if let Some(x) = g.get_mut(q + 1).and_then(|d| d.get_mut(i)) {
                        *x = b - even + odd;
                    }
                }
            }
        });
    }

    /// Pull streaming of the scalar lattice into `g2` with the wall rules of the links (the
    /// first seven directions), then the heat received per group.
    fn stream_thermal(&mut self) {
        let block = self.block;
        let n = self.n;
        let s = self.stride;
        let (links, offset) = (&self.links, &self.offset);
        let Some(t) = &mut self.thermal else { return };
        let (g, groups) = (&t.g, &t.groups);
        let mut blocks = cell_blocks(&mut t.g2, s, n, block);
        parallel::for_each_mut(&mut blocks, |b, out| {
            let start = b * block;
            for (q, dst) in out.iter_mut().enumerate().take(QT) {
                let off = offset.get(q).copied().unwrap_or(0);
                let src = g.get(q * s..q * s + n).unwrap_or_default();
                for (i, x) in dst.iter_mut().enumerate() {
                    if let Some(from) = (start + i).checked_add_signed(-off)
                        && let Some(v) = src.get(from)
                    {
                        *x = *v;
                    }
                }
            }
            let first = links.partition_point(|l| (l.cell as usize) < start);
            for link in links.get(first..).unwrap_or_default().iter().take_while(|l| (l.cell as usize) < start + block) {
                let (cell, q) = (link.cell as usize, link.q as usize);
                if q >= QT {
                    continue;
                }
                let v = thermal_rule(g, s, groups, cell, q, link.rule);
                if let (Some(v), Some(x)) = (v, out.get_mut(q).and_then(|d| d.get_mut(cell - start))) {
                    *x = v;
                }
            }
        });
        t.heat.fill(0.0);
        for link in links {
            let (cell, q) = (link.cell as usize, link.q as usize);
            if q >= QT {
                continue;
            }
            if let (Rule::Bounce { group, .. } | Rule::Reflect { group, .. }, Some(back)) =
                (link.rule, thermal_rule(g, s, groups, cell, q, link.rule))
                && !t.conjugate.get(group as usize).copied().unwrap_or(false)
                && !t.resets(link.cell)
            {
                let out = g.get(opposite(q) * s + cell).copied().unwrap_or(0.0);
                if let Some(h) = t.heat.get_mut(group as usize) {
                    // what left the fluid toward the wall minus what came back
                    *h += out - back;
                }
            }
        }
        // Links into conducting cells: anti-bounce-back at the interface temperature, the heat
        // to the solid cell.
        // The heat a link takes, `out − back` with `back = −out + 2 w T_i`, is affine in the
        // solid's temperature: `a − b T_s`. The solid steps with it implicit, then each link
        // takes its heat at the solid's mean temperature over the step.
        if let Some(sh) = &mut t.solid {
            let link = |i: &Interface| {
                let (cell, q) = (i.cell as usize, i.q as usize);
                let w2 = 2.0 * WT.get(q).copied().unwrap_or(0.0);
                let tf = t.temperature.get(cell).copied().unwrap_or(0.0);
                let out = g.get(opposite(q) * s + cell).copied().unwrap_or(0.0);
                (out, 2.0 * out - w2 * (1.0 - i.weight) * tf, w2 * i.weight)
            };
            for k in 0..sh.interface.len() {
                let Some(i) = sh.interface.get(k).copied() else { continue };
                let (_, a, b) = link(&i);
                sh.link(i.solid, a, b);
            }
            t.solid_finite = sh.step();
            for i in &sh.interface {
                let (out, a, b) = link(i);
                let heat = a - b * sh.effective.get(i.solid as usize).copied().unwrap_or(0.0);
                if let Some(x) = t.g2.get_mut(i.q as usize * s + i.cell as usize) {
                    *x = out - heat;
                }
                if let Some(h) = t.heat.get_mut(i.group as usize) {
                    *h += heat;
                }
            }
            // the box faces' heat from the solids (the box faces are the first six groups)
            for (h, w) in t.heat.iter_mut().zip(sh.wall_heat) {
                *h += w;
            }
        }
    }

    /// Inlets hold their temperature; outlets take their neighbour's populations.
    fn open_boundaries_thermal(&mut self) {
        let s = self.stride;
        let (open, offset, cell_u) = (&self.open, &self.offset, &self.cell_u);
        let Some(t) = &mut self.thermal else { return };
        for (oc, held) in open.iter().zip(&t.open) {
            let cell = oc.cell as usize;
            match held {
                Some(value) => {
                    let u = cell_u.get(cell).copied().unwrap_or([0.0; 3]);
                    for q in 0..QT {
                        if let Some(x) = t.g2.get_mut(q * s + cell) {
                            *x = equilibrium_t(q, *value, u);
                        }
                    }
                }
                None => {
                    // the neighbour one cell inward along the face's axis
                    let axis_dir = 1 + 2 * oc.axis + usize::from(oc.inward < 0.0);
                    let off = offset.get(axis_dir).copied().unwrap_or(0);
                    if let Some(from) = cell.checked_add_signed(off) {
                        for q in 0..QT {
                            let v = t.g2.get(q * s + from).copied();
                            if let (Some(v), Some(x)) = (v, t.g2.get_mut(q * s + cell)) {
                                *x = v;
                            }
                        }
                    }
                }
            }
        }
    }

    /// The solids' temperatures (stepped with the scalar lattice's streaming, through their
    /// links) into the field; false when one is not finite.
    fn conduct(&mut self) -> bool {
        let Some(t) = &mut self.thermal else { return true };
        let Some(sh) = &t.solid else { return true };
        let finite = t.solid_finite;
        for (cell, temp) in sh.cells.iter().zip(&sh.temperature) {
            if let Some(slot) = t.temperature.get_mut(*cell as usize) {
                *slot = *temp;
            }
        }
        finite
    }

    /// The temperature of every fluid cell from its populations; false when one is not finite.
    fn update_temperature(&mut self) -> bool {
        let n = self.n;
        let s = self.stride;
        let fluid = &self.fluid;
        let Some(t) = &mut self.thermal else { return true };
        let mut finite = true;
        for cell in 0..n {
            if !fluid.get(cell).copied().unwrap_or(false) {
                continue;
            }
            let mut temp = 0.0;
            for q in 0..QT {
                temp += t.g.get(q * s + cell).copied().unwrap_or(0.0);
            }
            finite &= temp.is_finite();
            if let Some(slot) = t.temperature.get_mut(cell) {
                *slot = temp;
            }
        }
        finite
    }

    /// The collision with Guo forcing (two relaxation times, or recursive regularised), in place,
    /// by blocks of [`block_cells`] cells in parallel. Returns the largest squared lattice speed, or
    /// NaN when a density is not positive.
    fn collide(&mut self) -> f64 {
        let block = self.block;
        let n = self.n;
        let regularized = self.collision == Collision::Regularized;
        let (wp, wm) = (self.omega_plus, self.omega_minus);
        let (sp, sm) = (1.0 - 0.5 * wp, 1.0 - 0.5 * wm);
        let constant = self.acceleration;
        let cell_force = &self.cell_force;
        let forced = constant != [0.0; 3] || !cell_force.is_empty() || !self.zone_of.is_empty();
        let fluid = &self.fluid;
        let (zone_of, actions) = (&self.zone_of, &self.actions);
        let mut blocks = cell_blocks(&mut self.f, self.stride, n, block);
        let mut maxima = vec![0.0f64; blocks.len()];
        // with heat, each block's velocities go to its slice of cell_u (empty otherwise)
        let mut u_blocks: Vec<&mut [[f64; 3]]> = self.cell_u.chunks_mut(block).collect();
        let mut work: Vec<(&mut Vec<&mut [f64]>, &mut f64, Option<&mut [[f64; 3]]>)> = Vec::with_capacity(blocks.len());
        let mut u_iter = u_blocks.drain(..);
        for (f, m) in blocks.iter_mut().zip(maxima.iter_mut()) {
            work.push((f, m, u_iter.next()));
        }
        parallel::for_each_mut(&mut work, |b, (f, block_max, u_out)| {
            let mut fq = [0.0; Q];
            let mut post = [0.0; Q];
            // The block's largest speed, kept here and stored once: the blocks' maxima share
            // cache lines, so a store per cell bounced those lines between the threads (#643).
            let mut max_u2 = 0.0f64;
            let len = f.first().map_or(0, |d| d.len());
            for i in 0..len {
                if !fluid.get(b * block + i).copied().unwrap_or(false) {
                    continue;
                }
                let a = cell_force.get(b * block + i).copied().unwrap_or(constant);
                let (mut rho, mut j) = (0.0, [0.0; 3]);
                for q in 0..Q {
                    let v = f.get(q).and_then(|d| d.get(i)).copied().unwrap_or(0.0);
                    fq[q] = v;
                    let cq = C[q];
                    rho += v;
                    j[0] += v * f64::from(cq[0]);
                    j[1] += v * f64::from(cq[1]);
                    j[2] += v * f64::from(cq[2]);
                }
                if !(rho > 0.0 && rho.is_finite()) {
                    **block_max = f64::NAN;
                    return;
                }
                let zone = zone_of.get(b * block + i).and_then(|z| actions.get(*z as usize));
                let (u, force) = velocity_force(rho, j, a, zone);
                let uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
                max_u2 = max_u2.max(uu);
                if let Some(slot) = u_out.as_mut().and_then(|u| u.get_mut(i)) {
                    *slot = u;
                }
                let uf = u[0] * force[0] + u[1] * force[1] + u[2] * force[2];
                let mut put = |q: usize, v: f64| {
                    if let Some(x) = f.get_mut(q).and_then(|d| d.get_mut(i)) {
                        *x = v;
                    }
                };
                if regularized {
                    crate::regularized::collide(&fq, rho, u, force, wp, forced, &mut post);
                    for (q, v) in post.iter().enumerate() {
                        put(q, *v);
                    }
                    continue;
                }
                // Rest population: symmetric only.
                let eq0 = W[0] * rho * (1.0 - 1.5 * uu);
                let mut out0 = fq[0] - wp * (fq[0] - eq0);
                if forced {
                    out0 += sp * W[0] * (-3.0 * uf);
                }
                put(0, out0);
                for q in (1..Q).step_by(2) {
                    let o = q + 1;
                    let cq = C[q];
                    let cv = [f64::from(cq[0]), f64::from(cq[1]), f64::from(cq[2])];
                    let cu = cv[0] * u[0] + cv[1] * u[1] + cv[2] * u[2];
                    let w = W[q];
                    let eq_plus = w * rho * (1.0 + 4.5 * cu * cu - 1.5 * uu);
                    let eq_minus = w * rho * 3.0 * cu;
                    let plus = 0.5 * (fq[q] + fq[o]);
                    let minus = 0.5 * (fq[q] - fq[o]);
                    let relax_plus = wp * (plus - eq_plus);
                    let relax_minus = wm * (minus - eq_minus);
                    let (mut s_plus, mut s_minus) = (0.0, 0.0);
                    if forced {
                        let cf = cv[0] * force[0] + cv[1] * force[1] + cv[2] * force[2];
                        s_plus = sp * w * (9.0 * cu * cf - 3.0 * uf);
                        s_minus = sm * w * 3.0 * cf;
                    }
                    put(q, fq[q] - relax_plus - relax_minus + s_plus + s_minus);
                    put(o, fq[o] - relax_plus + relax_minus + s_plus - s_minus);
                }
            }
            **block_max = max_u2;
        });
        // The largest speed of all blocks (a maximum does not depend on the order), or NaN.
        if maxima.iter().any(|m| m.is_nan()) { f64::NAN } else { maxima.iter().fold(0.0, |m, x| m.max(*x)) }
    }

    /// Pull streaming into `g`, with the wall rules (blocks of [`block_cells`] destination cells in
    /// parallel), then the momentum exchange (in link order).
    fn stream(&mut self) {
        let block = self.block;
        let n = self.n;
        let s = self.stride;
        let (f, links, offset) = (&self.f, &self.links, &self.offset);
        let mut blocks = cell_blocks(&mut self.g, s, n, block);
        parallel::for_each_mut(&mut blocks, |b, g| {
            let start = b * block;
            // Every direction as one shifted copy: right for cells whose neighbours are all fluid
            // cells in the box; the others are overwritten by their rules below.
            for (q, dst) in g.iter_mut().enumerate() {
                let off = offset.get(q).copied().unwrap_or(0);
                let src = f.get(q * s..q * s + n).unwrap_or_default();
                for (i, x) in dst.iter_mut().enumerate() {
                    if let Some(from) = (start + i).checked_add_signed(-off)
                        && let Some(v) = src.get(from)
                    {
                        *x = *v;
                    }
                }
            }
            let first = links.partition_point(|l| (l.cell as usize) < start);
            for link in links.get(first..).unwrap_or_default().iter().take_while(|l| (l.cell as usize) < start + block) {
                let (cell, q) = (link.cell as usize, link.q as usize);
                let v = match link.rule {
                    Rule::Pull { from } => f.get(q * s + from as usize).copied(),
                    Rule::Bounce { du, .. } => f.get(opposite(q) * s + cell).map(|out| out + du),
                    Rule::Reflect { from, q: q2, .. } => f.get(q2 as usize * s + from as usize).copied(),
                    // A placeholder the open boundary replaces.
                    Rule::Open => f.get(opposite(q) * s + cell).copied(),
                };
                if let (Some(v), Some(x)) = (v, g.get_mut(q).and_then(|d| d.get_mut(cell - start))) {
                    *x = v;
                }
            }
        });
        self.forces.fill([0.0; 3]);
        for link in &self.links {
            let (cell, q) = (link.cell as usize, link.q as usize);
            let cq = c(q);
            match link.rule {
                Rule::Bounce { group, du } => {
                    let out = self.f[opposite(q) * s + cell];
                    let back = out + du;
                    if let Some(force) = self.forces.get_mut(group as usize) {
                        // c_out f* − c_in f' with c_out = −c_q.
                        for a in 0..3 {
                            force[a] -= cq[a] * (out + back);
                        }
                    }
                }
                Rule::Reflect { from, q: q2, group } => {
                    let q2 = q2 as usize;
                    let out = self.f[q2 * s + from as usize];
                    let c2 = c(q2);
                    if let Some(force) = self.forces.get_mut(group as usize) {
                        for a in 0..3 {
                            force[a] += c2[a] * out - cq[a] * out;
                        }
                    }
                }
                Rule::Pull { .. } | Rule::Open => {}
            }
        }
    }

    /// Zou–He velocity and pressure boundaries on the streamed populations, and convective
    /// outlets ([`crate::outlet`]), whose convected populations are all taken from the streamed
    /// populations before any open cell is changed. Each open cell's boundary reads only its own
    /// populations, so the cells are worked out in parallel (blocks of [`OPEN_BLOCK`]) and
    /// stored afterwards in cell order: done one after another they were a tenth of a step on
    /// one thread, which capped the speed-up of many (#643).
    fn open_boundaries(&mut self) {
        let s = self.stride;
        let regularized = self.collision == Collision::Regularized;
        let mut outlets = std::mem::take(&mut self.outlets);
        let [nx, ny, _] = self.grid.cells;
        outlets.convect(&self.open, &View { flow: self, f: &self.g }, [1, nx, nx * ny]);
        let mut done = std::mem::take(&mut self.open_scratch);
        done.clear();
        done.resize(self.open.len(), [0.0; Q]);
        {
            let (flow, outlets) = (&*self, &outlets);
            parallel::chunks_mut(&mut done, OPEN_BLOCK, |b, chunk| {
                for (i, p) in chunk.iter_mut().enumerate() {
                    let k = b * OPEN_BLOCK + i;
                    let Some(oc) = flow.open.get(k) else { continue };
                    let cell = oc.cell as usize;
                    for (q, x) in p.iter_mut().enumerate() {
                        *x = flow.g.get(q * s + cell).copied().unwrap_or(0.0);
                    }
                    let a_half = flow.force_at(cell).map(|v| 0.5 * v);
                    let (axis, sign) = (oc.axis, oc.inward);
                    // The cell's populations alone: direction stride 1, cell 0.
                    match (oc.target, outlets.action(k)) {
                        (Target::Convective(_), Action::Convect { unknowns, shift }) => {
                            for (q, v) in crate::outlet::unknowns(axis, sign).into_iter().zip(unknowns) {
                                if let Some(x) = p.get_mut(q) {
                                    *x = v;
                                }
                            }
                            for (x, w) in p.iter_mut().zip(W) {
                                *x += w * shift;
                            }
                        }
                        (Target::Convective(rho), Action::Hold) => zou_he(p, 1, 0, axis, sign, a_half, crate::outlet::backflow(rho), regularized),
                        (target, _) => zou_he(p, 1, 0, axis, sign, a_half, target, regularized),
                    }
                }
            });
        }
        for (k, p) in done.iter().enumerate() {
            let Some(oc) = self.open.get(k).copied() else { continue };
            let cell = oc.cell as usize;
            for (q, v) in p.iter().enumerate() {
                if let Some(x) = self.g.get_mut(q * s + cell) {
                    *x = *v;
                }
            }
            if matches!(oc.target, Target::Convective(_)) {
                outlets.record(k, &oc, self.force_at(cell).map(|v| 0.5 * v), |q| p.get(q).copied().unwrap_or(0.0));
            }
        }
        self.open_scratch = done;
        self.outlets = outlets;
    }
}

/// A flow's populations as the outlets read them: `f` is the flow's own array or the one it
/// streamed into.
struct View<'a> {
    flow: &'a Flow,
    f: &'a [f64],
}

impl Populations for View<'_> {
    fn get(&self, q: usize, cell: usize) -> f64 {
        self.f.get(q * self.flow.stride + cell).copied().unwrap_or(0.0)
    }

    fn is_fluid(&self, cell: usize) -> bool {
        self.flow.is_fluid(cell)
    }

    fn half_force(&self, cell: usize) -> [f64; 3] {
        self.flow.force_at(cell).map(|v| 0.5 * v)
    }
}

/// Zou & He 1997 on open cell `cell` (normal `axis`, `sign` +1 at a min face) of the streamed
/// populations `g`: the unknown incoming populations are the opposite ones plus the equilibrium
/// difference (bounce-back of the non-equilibrium part), the density (a velocity inlet) or the
/// normal velocity (a density) follows from mass and normal momentum, and the transverse momentum
/// is corrected on the diagonal unknowns (Hecht & Harting 2010). `a_half` is half the body
/// acceleration. A convective outlet's cell takes the density, with no tangential velocity. With
/// the regularised collision (`regularized`) a density target then holds the equilibrium of its
/// density and that normal velocity (see the module notes): a pressure outlet, and a convective
/// outlet's cell where the flow comes back in ([`crate::outlet::backflow`]).
#[allow(clippy::too_many_arguments)]
fn zou_he(g: &mut [f64], s: usize, cell: usize, axis: usize, sign: f64, a_half: [f64; 3], target: Target, regularized: bool) {
    let ci = |q: usize| f64::from(C[q][axis]);
    let at = |g: &[f64], q: usize| g.get(q * s + cell).copied().unwrap_or(0.0);
    // Tangential (c_n = 0) and outgoing (c_n = −s, known) sums.
    let (mut tangential, mut outgoing) = (0.0, 0.0);
    for q in 0..Q {
        let v = at(g, q);
        let cn = ci(q) * sign;
        if cn == 0.0 {
            tangential += v;
        } else if cn < 0.0 {
            outgoing += v;
        }
    }
    // The bare velocity (without the half force) the populations must carry.
    let (rho, bare) = match target {
        Target::Velocity(u) => {
            let bare = [u[0] - a_half[0], u[1] - a_half[1], u[2] - a_half[2]];
            let un = sign * bare[axis];
            ((tangential + 2.0 * outgoing) / (1.0 - un), bare)
        }
        Target::Density(rho) | Target::Convective(rho) => {
            let un = 1.0 - (tangential + 2.0 * outgoing) / rho;
            let mut bare = [-a_half[0], -a_half[1], -a_half[2]];
            bare[axis] = sign * un;
            (rho, bare)
        }
    };
    // Unknowns: bounce-back of the non-equilibrium part.
    for q in (1..Q).filter(|q| ci(*q) * sign > 0.0) {
        let cq = c(q);
        let cu = cq[0] * bare[0] + cq[1] * bare[1] + cq[2] * bare[2];
        let v = at(g, opposite(q)) + 6.0 * W[q] * rho * cu;
        if let Some(x) = g.get_mut(q * s + cell) {
            *x = v;
        }
    }
    if regularized && matches!(target, Target::Density(_)) {
        hold_equilibrium(g, s, cell, rho, bare);
        return;
    }
    // Transverse momentum: correct the two diagonal unknowns along each tangential axis.
    for t in (0..3).filter(|t| *t != axis) {
        let mut m = -rho * bare[t];
        for q in 0..Q {
            m += at(g, q) * f64::from(C[q][t]);
        }
        for q in (1..Q).filter(|q| ci(*q) * sign > 0.0 && C[*q][t] != 0) {
            if let Some(x) = g.get_mut(q * s + cell) {
                *x -= 0.5 * m * f64::from(C[q][t]);
            }
        }
    }
}

/// Set every streamed population of open cell `cell` to the equilibrium of `rho` and `bare` (the
/// velocity without the half force): what an outflow cell holds under the regularised collision,
/// whose rebuilt non-equilibrium part would otherwise grow there (see the module notes). `g` holds
/// the populations with direction stride `s`.
pub(crate) fn hold_equilibrium(g: &mut [f64], s: usize, cell: usize, rho: f64, bare: [f64; 3]) {
    for q in 0..Q {
        if let Some(x) = g.get_mut(q * s + cell) {
            *x = equilibrium(q, rho, bare);
        }
    }
}

/// What a scalar population of a cell next to a wall or open face becomes (`None`: streamed as
/// usual): a periodic image pulls, an adiabatic wall (or a slip face) bounces back, a held
/// temperature anti-bounces back, a flux adds its heat, an open face is set afterwards.
fn thermal_rule(g: &[f64], s: usize, groups: &[ThermalCondition], cell: usize, q: usize, rule: Rule) -> Option<f64> {
    let back = g.get(opposite(q) * s + cell).copied()?;
    match rule {
        Rule::Pull { from } => g.get(q * s + from as usize).copied(),
        Rule::Bounce { group, .. } | Rule::Reflect { group, .. } => Some(match groups.get(group as usize).copied().unwrap_or_default() {
            ThermalCondition::Adiabatic {} => back,
            ThermalCondition::Temperature { value } => -back + 2.0 * WT.get(q).copied().unwrap_or(0.0) * value,
            ThermalCondition::Flux { value } => back + value,
        }),
        Rule::Open => Some(back),
    }
}

impl crate::backend::Backend for Flow {
    fn checkpoint(&mut self, a: &mut dyn crate::checkpoint::Archive) -> Result<()> {
        use crate::checkpoint::{boolean, doubles, f64, fixed, u64, vectors};
        // Streaming targets are allocated once, not cloned into a separate snapshot.
        self.targets()?;
        fixed(a, self.n as u64)?;
        fixed(a, self.stride as u64)?;
        for v in self.grid.cells {
            fixed(a, v as u64)?;
        }
        for v in self.grid.origin.into_iter().chain([self.grid.spacing, self.omega_plus, self.omega_minus]).chain(self.acceleration) {
            fixed(a, v.to_bits())?;
        }
        for v in &self.fluid {
            fixed(a, u64::from(*v))?;
        }
        fixed(a, self.links.len() as u64)?;
        for l in &self.links {
            fixed(a, l.cell as u64)?;
            fixed(a, l.q as u64)?;
            match l.rule {
                Rule::Pull { from } => {
                    fixed(a, 0)?;
                    fixed(a, from as u64)?;
                }
                Rule::Bounce { group, du } => {
                    fixed(a, 1)?;
                    fixed(a, group as u64)?;
                    fixed(a, du.to_bits())?;
                }
                Rule::Reflect { from, q, group } => {
                    fixed(a, 2)?;
                    fixed(a, from as u64)?;
                    fixed(a, q as u64)?;
                    fixed(a, group as u64)?;
                }
                Rule::Open => fixed(a, 3)?,
            }
        }
        u64(a, &mut self.steps_done)?;
        let n = self.f.len();
        doubles(a, &mut self.f, n)?;
        doubles(a, &mut self.g, n)?;
        vectors(a, &mut self.forces, self.groups.len())?;
        let nf = self.cell_force.len();
        vectors(a, &mut self.cell_force, nf)?;
        let nu = self.cell_u.len();
        vectors(a, &mut self.cell_u, nu)?;
        doubles(a, &mut self.fan_rise, self.zones.len())?;
        for action in &mut self.actions {
            for v in &mut action.push {
                f64(a, v)?;
            }
            f64(a, &mut action.alpha)?;
            f64(a, &mut action.beta)?;
        }
        self.outlets.checkpoint(a)?;
        fixed(a, u64::from(self.thermal.is_some()))?;
        if let Some(t) = &mut self.thermal {
            let nt = t.g.len();
            doubles(a, &mut t.g, nt)?;
            doubles(a, &mut t.g2, nt)?;
            doubles(a, &mut t.temperature, self.n)?;
            doubles(a, &mut t.heat, self.groups.len())?;
            boolean(a, &mut t.solid_finite)?;
            if !t.solid_finite {
                return Err(Error::Invalid("cannot resume a checkpoint with unstable solid temperatures".into()));
            }
            fixed(a, u64::from(t.solid.is_some()))?;
            if let Some(s) = &mut t.solid {
                s.checkpoint(a)?;
            }
        }
        // Scratch open populations are overwritten before use; their capacity is not state.
        Ok(())
    }

    fn name(&self) -> String {
        "cpu".into()
    }

    fn advance(&mut self, steps: u64) -> Result<()> {
        for _ in 0..steps {
            self.step()?;
        }
        Ok(())
    }

    fn steps_done(&self) -> u64 {
        self.steps_done
    }

    fn fields(&mut self) -> Result<crate::backend::Fields> {
        let n = self.n;
        let mut density = filled(n, 0.0, "density field")?;
        let mut velocity = filled(n, [0.0; 3], "velocity field")?;
        for cell in 0..n {
            if let (Some((rho, u)), Some(d), Some(v)) = (self.moments(cell), density.get_mut(cell), velocity.get_mut(cell)) {
                *d = rho;
                *v = u;
            }
        }
        let temperature = self.thermal.as_ref().map(|t| t.temperature.clone()).unwrap_or_default();
        Ok(crate::backend::Fields { density, velocity, temperature })
    }

    fn heat(&mut self) -> Result<Vec<f64>> {
        Ok(Flow::heat(self).to_vec())
    }

    fn zones(&mut self) -> Result<Vec<ZoneSample>> {
        Ok(self.zone_samples())
    }

    fn mass(&mut self) -> Result<f64> {
        Ok(Flow::mass(self))
    }

    fn forces(&mut self) -> Result<Vec<[f64; 3]>> {
        Ok(self.forces.clone())
    }

    fn groups(&self) -> &[Group] {
        &self.groups
    }

    fn fluid_cells(&self) -> usize {
        self.fluid_cells
    }

    fn is_fluid(&self, cell: usize) -> bool {
        Flow::is_fluid(self, cell)
    }

    fn empty_obstacles(&self) -> &[String] {
        &self.empty_obstacles
    }

    fn conducting(&self) -> &[Vec<u32>] {
        self.thermal.as_ref().and_then(|t| t.solid.as_ref()).map_or(&[], |s| &s.by_solid)
    }
}
