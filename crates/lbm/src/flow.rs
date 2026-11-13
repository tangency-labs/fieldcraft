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
//!   271:285). Slip faces reflect populations specularly. Under Han et al.'s on-node
//!   wall-function bounce (`spalding_bounce`, [`crate::bounce`]) the outermost cell layer at a
//!   wall face is the wall instead: those boundary nodes run their eq. 14 in place of the
//!   collision and are not fluid cells.
//! - **Inlets and outlets:** Zou & He 1997 (Phys. Fluids 9:1591) on the first layer of cells: the
//!   unknown incoming populations are the opposite ones plus the equilibrium difference
//!   (bounce-back of the non-equilibrium part), the density (velocity inlet) or the normal velocity
//!   (pressure outlet) follow from mass and normal momentum, and the transverse momentum is
//!   corrected on the diagonal unknowns as M. Hecht and J. Harting do for D3Q19 (*Implementation
//!   of on-site velocity boundary conditions for D3Q19 lattice Boltzmann simulations*, J. Stat.
//!   Mech. 2010, P01018). The prescribed velocity is the one with the half force of Guo's scheme.
//!   Uniform pressure faces instead extrapolate the inner cell's normal velocity and
//!   non-equilibrium populations (Guo, Zheng & Shi 2002, doi:10.1088/1009-1963/11/4/310),
//!   imposing density and zero tangential velocity without discarding viscous stress. The
//!   local Zou–He pressure balance amplifies a longitudinal normal-stress mode even in a
//!   uniform 4×1×1 stream (#819). Every pressure cell reads the pre-boundary streamed state,
//!   including at intersecting open faces. Without a fluid cell immediately inside, the
//!   legacy local closure is used: Zou–He for TRT/MRT, equilibrium for RR/HRR.
//!   A convective outlet ([`crate::outlet`]) updates the unknown populations by a convective
//!   outflow condition instead, holds each cell at the face's density, and falls back to Zou–He
//!   where the flow comes back in. A cell on the edge of two open faces belongs to one of them.
//! - **Forces:** momentum exchange on every wall link (Ladd 1994; R. Mei, D. Yu, W. Shyy, L.-S.
//!   Luo, *Force evaluation in the lattice Boltzmann method involving curved geometry*, PRE
//!   65:041203, 2002): a population `f*` leaving a fluid cell toward the wall and the one `f'`
//!   coming back transfer `c_out f* − c_in f'` to the wall. Under the wall model a wall also
//!   takes the momentum each reconstruction of its wall cells removes ([`crate::wall`]); under
//!   the on-node bounce, the momentum its boundary nodes take from the fluid ([`crate::bounce`]).
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
//! # Melting and freezing
//! With a phase change in the fluid's material the scalar lattice carries the total enthalpy in
//! temperature units, `E = T + (L / c) f_l` (f_l the liquid fraction), instead of the temperature:
//! the rest population's equilibrium becomes `E − T + w₀ T`, the others keep `w_q T (1 + 4 c_q·u)`,
//! so the populations sum to `E`, their flux is `T u` and the diffusion acts on `T`, and the
//! lattice solves `∂E/∂t + ∇·(T u) = ∇·(α ∇T)` with the latent heat released as the enthalpy
//! falls, without a source term or an iteration (R. Huang, H. Wu, P. Cheng, *A new lattice
//! Boltzmann model for solid–liquid phase change*, Int. J. Heat Mass Transfer 59:295–301, 2013;
//! the enthalpy method of V. R. Voller, *Fast implicit finite-difference method for the analysis
//! of phase change problems*, Numer. Heat Transfer B 17:155, 1990). The temperature and liquid
//! fraction follow from `E` (`PhaseSetup::split`): solid below the solidus, liquid above the
//! liquidus plus `L / c`, the liquid fraction linear in between (an isothermal front when solidus
//! and liquidus are equal). Walls, inlets and the conjugate interface work on `T` as before (an
//! inlet holds the enthalpy of its temperature). The flow stops where the fluid has frozen: each
//! cell's velocity carries a Carman–Kozeny resistance `A (1 − f_l)² / (f_l³ + ε)` per step
//! (V. R. Voller, C. Prakash, *A fixed grid numerical modelling methodology for convection–
//! diffusion mushy region phase-change problems*, Int. J. Heat Mass Transfer 30:1709, 1987;
//! [`MUSHY`]), taken implicitly as a porous zone's (Guo and Zhao), so a frozen cell holds still to
//! a few parts in a million of the speed it would take and the mushy zone resists the flow.
//!
//! Obstacles that conduct heat carry their own temperature, conducted on their cells and coupled
//! to this lattice through their faces (conjugate heat transfer, [`crate::conduction`]).
//!
//! # Free surface
//! With a free surface ([`crate::surface`]) the fluid cells are liquid, interface or gas; only the
//! liquid and interface cells collide and stream, the interface cells exchange mass and rebuild
//! the populations the gas does not give them after each streaming, and the cells convert as they
//! fill and empty after each step. Velocity inlets inject liquid.
//!
//! # Layout and determinism
//! Two arrays of populations, structure of arrays (`f[q n + cell]`), A–B pattern. Streaming
//! pulls: every population of a cell comes from one source (a neighbour, a periodic image, the
//! cell itself through a wall, or an open boundary). Cells whose 18 neighbours are all fluid
//! cells inside the box stream by one shifted copy per direction; the others follow per-link
//! rules. Collision and streaming run by fixed blocks of cells in parallel, each cell writing
//! only its own populations; the wall forces are summed in link order. There are no parallel
//! reductions, so the result is the same bit for bit on any number of threads and platform.

use fieldcraft_model::{BoxFace, Collision, FaceCondition, Outflow, Profile, ThermalCondition, WallModel};

use fieldcraft_parallel as parallel;

use crate::bounce::{WallNode, WallNodes};
use crate::conduction::{Interface, Layout, SolidHeat};
use crate::lattice::{C, Q, W, c, equilibrium, moving_wall, opposite};
use crate::outlet::{Action, Outlets, Populations};
use crate::rheology::{Clamp, Clamps, Rheology};
use crate::setup::{MAX_LATTICE_SPEED, Obstacle, PhaseSetup, Setup, ZoneKind, ZoneSetup};
use crate::surface::{Kind, Surface, SurfaceReport};
use crate::voxel::{Grid, TriangleTree, voxelise};
use crate::wall::Walls;
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

/// A fluid cell's state before its collision ([`cell_state`]).
struct CellState {
    rho: f64,
    /// The velocity with the half force, and the force per volume.
    u: [f64; 3],
    force: [f64; 3],
    /// Under the wall model, a rebuilt wall cell's relaxation rate and the momentum the
    /// rebuild took from it.
    wall: Option<(f64, [f64; 3])>,
}

/// The density and momentum of populations `fq` (summed in direction order).
#[inline(always)]
fn density_momentum(fq: &[f64; Q]) -> (f64, [f64; 3]) {
    let (mut rho, mut j) = (0.0, [0.0; 3]);
    for (q, v) in fq.iter().enumerate() {
        let cq = C[q];
        rho += v;
        j[0] += v * f64::from(cq[0]);
        j[1] += v * f64::from(cq[1]);
        j[2] += v * f64::from(cq[2]);
    }
    (rho, j)
}

/// The state cell `cell` collides from, its populations `fq` read: under the wall model (`walls`
/// given only then) a wall cell's populations are first rebuilt from the wall law; then its
/// velocity and force with body acceleration `a` and its zone's action. `None` when the density
/// is not a positive number. The collision and the hybrid collision's velocity pass both call
/// it, so their velocities are the same bit for bit.
#[inline(always)]
fn cell_state(fq: &mut [f64; Q], cell: usize, a: [f64; 3], walls: Option<&Walls>, zone: Option<&ZoneAction>) -> Option<CellState> {
    let (mut rho, mut j) = density_momentum(fq);
    let mut wall = None;
    if rho > 0.0
        && rho.is_finite()
        && let Some(rate) = walls.and_then(|w| w.reconstruct(cell, fq, a))
    {
        let before = j;
        (rho, j) = density_momentum(fq);
        wall = Some((rate, [0, 1, 2].map(|k| before[k] - j[k])));
    }
    if !(rho > 0.0 && rho.is_finite()) {
        return None;
    }
    let (u, force) = velocity_force(rho, j, a, zone);
    Some(CellState { rho, u, force, wall })
}

/// The fluid cell next to the cell at `p` along `axis` (towards `+axis` with `up`): across a
/// periodic face its periodic image, as streaming pulls; `None` when that cell is not a fluid
/// cell or lies outside a face that is not periodic.
#[inline(always)]
fn neighbour(grid: &Grid, periodic: [bool; 3], fluid: &[bool], p: [usize; 3], axis: usize, up: bool) -> Option<usize> {
    let n = *grid.cells.get(axis)?;
    let i = *p.get(axis)?;
    let periodic = periodic.get(axis).copied().unwrap_or(false);
    let next = match (up, i) {
        (true, i) if i + 1 < n => i + 1,
        (true, _) if periodic => 0,
        (false, 0) if periodic => n.checked_sub(1)?,
        (false, 0) | (true, _) => return None,
        (false, i) => i - 1,
    };
    let mut q = p;
    *q.get_mut(axis)? = next;
    let cell = grid.index(q);
    fluid.get(cell).copied().unwrap_or(false).then_some(cell)
}

/// A cell's action: its zone's, with the resistance of its liquid fraction `liquid` (with a phase
/// change) added where it is not wholly liquid.
#[inline(always)]
fn cell_action(zone: Option<&ZoneAction>, liquid: Option<f64>) -> Option<ZoneAction> {
    let r = liquid.map_or(0.0, Phase::resistance);
    if r > 0.0 {
        let mut z = zone.copied().unwrap_or_default();
        z.alpha += r;
        Some(z)
    } else {
        zone.copied()
    }
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

/// The Carman–Kozeny constants `(A, ε)` of the resistance `A (1 − f_l)² / (f_l³ + ε)` per step
/// (lattice units) that stops the flow where the fluid has frozen: a frozen cell's velocity is
/// cut to about `2 ε / A` (2e-6; [`velocity_force`] gives `u = v / (1 + α / 2)` with `α = A / ε`)
/// of what it would be, half-frozen fluid to 1/100, and the flow is barely touched above a liquid
/// fraction of 0.95 (Voller and Prakash 1987).
pub const MUSHY: (f64, f64) = (100.0, 1.0e-4);

/// The fluid's melting and freezing on the temperature lattice (see the module notes).
#[derive(Debug, Clone)]
pub(crate) struct Phase {
    rule: PhaseSetup,
    /// Each cell's liquid fraction after the last step (1 in obstacle cells).
    liquid: Vec<f64>,
}

impl PhaseSetup {
    /// The temperature and the liquid fraction at total enthalpy `e` (temperature units).
    fn split(&self, e: f64) -> (f64, f64) {
        let (low, high) = (self.solidus, self.liquidus + self.latent);
        if e <= low {
            (e, 0.0)
        } else if e >= high {
            (e - self.latent, 1.0)
        } else {
            let f = (e - low) / (high - low);
            (self.solidus + f * (self.liquidus - self.solidus), f)
        }
    }

    /// The liquid fraction at temperature `t`: a fluid at its liquidus counts as liquid.
    fn fraction(&self, t: f64) -> f64 {
        if t >= self.liquidus {
            1.0
        } else if t <= self.solidus {
            0.0
        } else {
            (t - self.solidus) / (self.liquidus - self.solidus)
        }
    }

    /// The total enthalpy (temperature units) and the liquid fraction of the fluid at temperature
    /// `t`.
    fn at(&self, t: f64) -> (f64, f64) {
        let f = self.fraction(t);
        (t + self.latent * f, f)
    }
}

impl Phase {
    /// The flow resistance per step at liquid fraction `f` (none when liquid).
    fn resistance(f: f64) -> f64 {
        if f >= 1.0 {
            return 0.0;
        }
        let (a, eps) = MUSHY;
        let s = 1.0 - f;
        a * s * s / (f * f * f + eps)
    }
}

/// The scalar equilibrium of total enthalpy `e` at temperature `t`: the temperature's, with the
/// latent part `E − T` on the rest population (`e = t` without a phase change).
#[inline(always)]
fn equilibrium_e(q: usize, e: f64, t: f64, u: [f64; 3]) -> f64 {
    equilibrium_t(q, t, u) + if q == 0 { e - t } else { 0.0 }
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
    /// The fluid's melting and freezing, when its material has a phase change.
    phase: Option<Phase>,
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

    /// The total enthalpy at temperature `t`: `t` itself without a phase change.
    fn enthalpy(&self, t: f64) -> f64 {
        self.phase.as_ref().map_or(t, |p| p.rule.at(t).0)
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
    pub(crate) grid: Grid,
    pub(crate) n: usize,
    /// The distance between two directions' arrays in `f` and `g` (at least `n`).
    pub(crate) stride: usize,
    fluid: Vec<bool>,
    pub(crate) f: Vec<f64>,
    pub(crate) g: Vec<f64>,
    /// The links that do not stream by the shifted copies, in cell order.
    links: Vec<Link>,
    open: Vec<OpenCell>,
    omega_plus: f64,
    omega_minus: f64,
    collision: Collision,
    /// Smagorinsky's constant (0: no eddy viscosity) and the turbulent Prandtl number that
    /// turns its eddy viscosity into the temperature lattice's turbulent diffusivity.
    smagorinsky: f64,
    turbulent_prandtl: f64,
    /// Vreman's model constant c = 2.5 C_s² (0: not this model; [`crate::vreman`]).
    vreman: f64,
    /// With heat and an eddy viscosity: each fluid cell's eddy viscosity of the last collision
    /// (lattice units), for the temperature lattice (empty otherwise).
    cell_nu_t: Vec<f64>,
    /// The largest eddy viscosity over the fluid's own of the last collision (0 without one).
    eddy_ratio: f64,
    /// With Van Driest's damping or the wall model: the wall cells (and with damping every
    /// cell's wall distance and the friction velocities).
    walls: Option<Walls>,
    /// Whether the wall cells relax with the wall model's viscosity.
    wall_model: bool,
    /// Whether the wall-function bounce gives the walls their stress (`spalding_bounce_halfway`).
    wall_bounce: bool,
    /// Han et al.'s on-node wall-function bounce (`spalding_bounce`): its boundary nodes.
    nodes: Option<WallNodes>,
    /// With it, the cells that collide: the fluid cells less the boundary nodes (empty
    /// otherwise: every fluid cell collides).
    bulk: Vec<bool>,
    /// A shear-dependent fluid: each cell's relaxation time from its own shear rate.
    rheology: Option<Rheology>,
    /// With one: the cells whose relaxation time was clamped, low and high, in the last
    /// collision and over the run.
    clamps: Clamps,
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
    /// Thin walls: their links are no-slip bounce-back; with heat each passes its conductive
    /// exchange after the thermal streaming.
    baffles: Vec<crate::baffle::Baffle>,
    /// The temperature lattice, with thermal settings.
    thermal: Option<ThermalState>,
    /// With buoyancy: every cell's body acceleration (the constant one plus the buoyant one),
    /// refreshed from the temperature before each step.
    cell_force: Vec<[f64; 3]>,
    /// With heat: every cell's velocity from the last collision, for the scalar equilibrium.
    cell_u: Vec<[f64; 3]>,
    /// The hybrid collision's σ, and with it or with Vreman's model every cell's velocity before
    /// the collision (for the finite-difference gradient; empty otherwise).
    sigma: f64,
    gradient_u: Vec<[f64; 3]>,
    /// Whether σ is dynamic (set per cell from the subgrid model, [`crate::hrr`]), and then its
    /// smallest and mean value over the fluid cells in the last collision.
    dynamic: bool,
    sigma_stats: Option<(f64, f64)>,
    /// Whether each axis is periodic (its faces' neighbours are the periodic images).
    pub(crate) periodic: [bool; 3],
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
    pub(crate) block: usize,
    /// Each open cell's populations after its boundary, in the step (kept for its capacity).
    open_scratch: Vec<[f64; Q]>,
    /// The free surface: which fluid cells hold liquid, with a free surface.
    pub(crate) surface: Option<Surface>,
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

/// The cells that collide: a free surface's liquid and interface cells, else the fluid cells
/// less the boundary nodes of the on-node wall-function bounce (`bulk`, when it is set).
/// Field-level borrows, for passes that hold other fields mutably.
fn populated_mask<'a>(fluid: &'a [bool], bulk: &'a [bool], surface: Option<&'a Surface>) -> &'a [bool] {
    match surface {
        Some(s) => &s.active,
        None if bulk.is_empty() => fluid,
        None => bulk,
    }
}

/// The first `n` populations of each direction of an array with direction stride `stride`, cut
/// into blocks of `block` cells: block b holds the Q directions' slices of cells b·block.. .
pub(crate) fn cell_blocks(data: &mut [f64], stride: usize, n: usize, block: usize) -> Vec<Vec<&mut [f64]>> {
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
pub(crate) fn filled<T: Clone>(len: usize, value: T, what: &str) -> Result<Vec<T>> {
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
        // Each obstacle face group's wall velocity, for Van Driest's damping (lattice units).
        let mut obstacle_walls: Vec<Option<[f64; 3]>> = Vec::new();
        for (k, faces) in setup.obstacles.iter().enumerate() {
            let o = obstacles.get(faces.index).ok_or_else(|| Error::NotFound(format!("obstacle {}", faces.index)))?;
            let mut inside = filled(n, false, "obstacle cells")?;
            voxelise(&o.surface, &grid, &mut inside)?;
            crate::voxel::decide_near_curved(&o.surface, &o.exact, &grid, &mut inside)?;
            if !inside.contains(&true) {
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
            obstacle_walls.extend((0..faces.names.len()).map(|t| Some(setup.to_lattice(faces.velocities.get(t).copied().unwrap_or([0.0; 3])))));
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
        let baffles = crate::baffle::Baffle::prepare(setup, &mut groups)?;
        for b in &baffles {
            b.validate_cells(&grid, &fluid)?;
        }
        let periodic = [domain.faces.periodic(0), domain.faces.periodic(1), domain.faces.periodic(2)];
        let wall_velocity = |face: BoxFace| match *domain.faces.get(face) {
            FaceCondition::MovingWall { velocity } => Some(setup.to_lattice(velocity)),
            FaceCondition::Wall {} => Some([0.0; 3]),
            _ => None,
        };
        let is_slip = |face: BoxFace| matches!(domain.faces.get(face), FaceCondition::Slip {});
        // The wall-function bounce reflects at its walls instead of bouncing back; only at rest.
        let wall_bounce = setup.wall_model == Some(WallModel::SpaldingBounceHalfway);
        if wall_bounce && let Some(face) = BoxFace::ALL.iter().find(|f| matches!(domain.faces.get(**f), FaceCondition::MovingWall { .. })) {
            return Err(Error::Invalid(format!(
                "not supported yet: the wall model `spalding_bounce_halfway` on a moving wall (`{}`); use `musker` or a wall at rest",
                face.name()
            )));
        }
        // Han et al.'s on-node wall-function bounce: the outermost cell layer at each wall face
        // is its boundary nodes ([`crate::bounce`]).
        let on_node = setup.wall_model == Some(WallModel::SpaldingBounce);
        let mut node_list: Vec<WallNode> = Vec::new();
        let mut bulk = Vec::new();
        if on_node {
            let refuse = |what: String| {
                Err(Error::Invalid(format!(
                    "not supported yet: Han et al.'s on-node wall-function bounce (`spalding_bounce`) {what}; choose another supported wall treatment"
                )))
            };
            if let Some(face) = BoxFace::ALL.iter().find(|f| matches!(domain.faces.get(**f), FaceCondition::MovingWall { .. })) {
                return refuse(format!("on a moving wall (`{}`)", face.name()));
            }
            if let Some(o) = setup.obstacles.first() {
                let body = obstacles.get(o.index).map_or("", |o| o.body.as_str());
                return refuse(format!("at the walls of an obstacle (`{body}`): it runs at plane box faces only"));
            }
            if setup.thermal.is_some() {
                return refuse("with heat".to_string());
            }
            for face in BoxFace::ALL.into_iter().filter(|f| wall_velocity(*f).is_some()) {
                let a = face.axis();
                if let Some(other) = BoxFace::ALL.into_iter().find(|f| f.axis() != a && !periodic[f.axis()]) {
                    return refuse(format!(
                        "where the wall `{}` meets the face `{}` that is not periodic (an edge): the box must be periodic along its walls",
                        face.name(),
                        other.name()
                    ));
                }
                if grid.cells[a] < 3 {
                    return refuse(format!("on a box {} cells across its wall `{}`: it needs at least 3", grid.cells[a], face.name()));
                }
                let (t1, t2) = ((a + 1) % 3, (a + 2) % 3);
                let (layer, inward) = if face.is_max() { (grid.cells[a] - 1, -1i32) } else { (0, 1) };
                for j in 0..grid.cells[t2] {
                    for i in 0..grid.cells[t1] {
                        let mut ijk = [0usize; 3];
                        ijk[a] = layer;
                        ijk[t1] = i;
                        ijk[t2] = j;
                        let cell = grid.index(ijk);
                        ijk[a] = if inward > 0 { layer + 1 } else { layer - 1 };
                        node_list.push(WallNode {
                            cell: cell as u32,
                            interior: grid.index(ijk) as u32,
                            axis: a as u8,
                            inward: inward as i8,
                            group: face.index() as u32,
                            stress: [0.0; 3],
                            taken: [0.0; 3],
                        });
                    }
                }
            }
            node_list.sort_by_key(|w| w.cell);
            bulk = fluid.clone();
            for w in &node_list {
                if let Some(b) = bulk.get_mut(w.cell as usize) {
                    *b = false;
                }
            }
        }
        let is_node = |cell: usize| bulk.get(cell).is_some_and(|b| !b) && fluid.get(cell).copied().unwrap_or(false);
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
            if !is_node(cell) {
                fluid_cells += 1;
            }
            let ijk = grid.coords(cell);
            let p = [ijk[0] as i64, ijk[1] as i64, ijk[2] as i64];
            for q in 1..Q {
                let cq = C[q];
                let sheet = baffles.iter().find_map(|b| b.group_for(ijk, cq));
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
                if let Some(group) = sheet {
                    // Through a sheet: bounce-back at it, or at its edge on a box wall, the wall's
                    // (with a moving wall's term, its force to the wall).
                    let rule = match crossed_faces().find_map(|f| wall_velocity(f).map(|u| (f, u))) {
                        Some((face, u)) => Rule::Bounce { group: face.index() as u32, du: moving_wall(q, u) },
                        None => Rule::Bounce { group, du: 0.0 },
                    };
                    links.push(Link { cell: cell as u32, q: q as u8, rule });
                    continue;
                }
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
                        // The specular reflection at `face`: from the same cell layer, the
                        // tangential components shifted; bounce-back at an edge or corner.
                        let reflect = |face: BoxFace| {
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
                        };
                        if let Some(face) = crossed_faces().find(|f| wall_velocity(*f).is_some()) {
                            if is_node(cell) {
                                // From outside the box into a boundary node: not needed, eq. 14
                                // overwrites it (a placeholder, no force).
                                Rule::Open
                            } else if wall_bounce {
                                reflect(face)
                            } else {
                                Rule::Bounce { group: face.index() as u32, du: moving_wall(q, wall_velocity(face).unwrap_or([0.0; 3])) }
                            }
                        } else if let Some(face) = crossed_faces().find(|f| is_slip(*f)) {
                            reflect(face)
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
                        if fluid.get(cell).copied().unwrap_or(false) && !is_node(cell) && zone_of.get(cell) == Some(&NONE) {
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

        // The free surface: liquid in the setup's regions and the inlet cells, gas elsewhere.
        let surface = match &setup.surface {
            None => None,
            Some(ss) => {
                let inlets: Vec<u32> = open.iter().filter(|oc| matches!(oc.target, Target::Velocity(_))).map(|oc| oc.cell).collect();
                let liquid = |cell: usize| {
                    let p = grid.coords(cell);
                    ss.liquid.iter().any(|(lo, hi)| (0..3).all(|a| p[a] >= lo[a] && p[a] <= hi[a]))
                };
                let reflect = BoxFace::ALL.map(|f| is_slip(f) && !periodic[f.axis()]);
                Some(Surface::new(&grid, periodic, reflect, &fluid, &inlets, liquid)?)
            }
        };
        let populated = |cell: usize| surface.as_ref().map_or(fluid[cell], |s| s.active[cell]);
        // The initial state: equilibrium at density 1 and the initial velocity.
        let mut f = filled(q_n, 0.0, "populations")?;
        // The streaming targets are allocated on the first CPU step: a flow handed to a device
        // backend never needs them, and they are half the CPU state.
        let g = Vec::new();
        for cell in 0..n {
            for q in 0..Q {
                f[q * stride + cell] = if populated(cell) { equilibrium(q, 1.0, setup.initial_velocity) } else { 0.0 };
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
                let t0 = ts.initial_temperature;
                let (e0, f0) = ts.phase.map_or((t0, 1.0), |p| p.at(t0));
                let phase = match ts.phase {
                    None => None,
                    Some(rule) => {
                        let mut liquid = filled(n, 1.0, "liquid fractions")?;
                        for (cell, f) in liquid.iter_mut().enumerate() {
                            if fluid[cell] {
                                *f = f0;
                            }
                        }
                        Some(Phase { rule, liquid })
                    }
                };
                for cell in 0..n {
                    if fluid[cell] {
                        for q in 0..QT {
                            tg[q * stride + cell] = equilibrium_e(q, e0, t0, setup.initial_velocity);
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
                    let layout =
                        Layout { grid: &grid, periodic, conducting: &conducting, walls, kf: ts.conductivity, contacts: &ts.contacts, zero: ts.zero };
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
                    phase,
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

        // Van Driest's damping and the wall model: each no-slip group's wall velocity, in group
        // order.
        let nu = (setup.tau_plus - 0.5) / 3.0;
        let walls = if on_node {
            // Van Driest's damping from the interior nodes next to the boundary nodes.
            let cells: Vec<(u32, usize, i32, u32)> =
                node_list.iter().map(|w| (w.interior, usize::from(w.axis), i32::from(w.inward), w.group)).collect();
            setup.van_driest.map(|a| Walls::on_node(&grid, &bulk, &cells, Some(a), nu)).transpose()?
        } else if setup.van_driest.is_some() || setup.wall_model.is_some() {
            let groups: Vec<Option<[f64; 3]>> = BoxFace::ALL.iter().map(|f| wall_velocity(*f)).chain(obstacle_walls).collect();
            Some(Walls::new(&grid, &fluid, &links, &groups, setup.van_driest, nu)?)
        } else {
            None
        };
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
            smagorinsky: setup.smagorinsky,
            turbulent_prandtl: setup.turbulent_prandtl,
            vreman: setup.vreman,
            cell_nu_t: if setup.turbulent() && setup.thermal.is_some() { filled(n, 0.0, "eddy viscosities")? } else { Vec::new() },
            eddy_ratio: 0.0,
            walls,
            wall_model: setup.wall_model == Some(WallModel::Musker),
            wall_bounce: setup.wall_model == Some(WallModel::SpaldingBounceHalfway),
            nodes: on_node.then(|| WallNodes::new(node_list, nu)),
            bulk,
            rheology: setup.rheology,
            clamps: Clamps::default(),
            acceleration: setup.acceleration,
            forces: vec![[0.0; 3]; groups.len()],
            groups,
            offset,
            fluid_cells,
            steps_done: 0,
            empty_obstacles,
            baffles,
            thermal,
            cell_force,
            cell_u,
            sigma: setup.sigma,
            dynamic: setup.dynamic_sigma,
            sigma_stats: None,
            gradient_u: if setup.collision == Collision::Hrr || setup.vreman > 0.0 { filled(n, [0.0; 3], "cell velocities")? } else { Vec::new() },
            periodic,
            fan_rise: vec![0.0; setup.zones.len()],
            zones: setup.zones.clone(),
            zone_cells,
            zone_of,
            actions,
            outlets: Outlets::default(),
            block: block_cells(n),
            open_scratch: Vec::new(),
            surface,
        };
        flow.outlets = Outlets::new(&flow.open, &View { flow: &flow, f: &flow.f });
        if let Some(s) = &mut flow.surface {
            s.initial_mass = s.total_mass(&flow.f, flow.stride);
        }
        Ok(flow)
    }

    /// What the convective outlet cells keep between steps.
    pub fn outlets(&self) -> &Outlets {
        &self.outlets
    }

    /// The free surface, when the flow has one.
    pub fn surface(&self) -> Option<&Surface> {
        self.surface.as_ref()
    }

    /// Whether the flow has a free surface.
    pub fn free_surface(&self) -> bool {
        self.surface.is_some()
    }

    /// With a free surface: how it fared so far ([`SurfaceReport`]).
    pub fn surface_report(&self) -> Option<SurfaceReport> {
        self.surface.as_ref().map(|s| s.report(&self.f, self.stride))
    }

    /// The cells that collide and stream now: with a free surface the liquid and interface
    /// cells, else [`Flow::collides`].
    fn populated(&self) -> &[bool] {
        match &self.surface {
            Some(s) => &s.active,
            None => self.collides(),
        }
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

    /// The liquid fraction of a cell, or `None` without a phase change.
    pub fn liquid_fraction(&self, cell: usize) -> Option<f64> {
        self.phase().and_then(|p| p.liquid.get(cell).copied())
    }

    /// Whether a conducting solid radiates to its surroundings.
    pub fn radiates(&self) -> bool {
        self.thermal.as_ref().and_then(|t| t.solid.as_ref()).is_some_and(|s| !s.radiated.is_empty())
    }

    /// Whether the fluid melts and freezes.
    pub fn phase_change(&self) -> bool {
        self.phase().is_some()
    }

    /// The fluid's phase change, when it has one.
    fn phase(&self) -> Option<&Phase> {
        self.thermal.as_ref().and_then(|t| t.phase.as_ref())
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
            let (e, f) = t.phase.as_ref().map_or((value, 1.0), |p| p.rule.at(value));
            if let Some(slot) = t.phase.as_mut().and_then(|p| p.liquid.get_mut(cell)) {
                *slot = f;
            }
            for q in 0..QT {
                if let Some(g) = t.g.get_mut(q * s + cell) {
                    *g = equilibrium_e(q, e, value, u);
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
            + bytes(&self.cell_nu_t)
            + bytes(&self.gradient_u)
            + bytes(&self.zone_of)
            + bytes(&self.open_scratch)
            + bytes(&self.baffles)
            + self.zone_cells.iter().map(bytes).sum::<u64>()
            + self.outlets.heap_bytes();
        if let Some(t) = &self.thermal {
            total += bytes(&t.g) + bytes(&t.g2) + bytes(&t.temperature) + bytes(&t.open) + bytes(&t.reset);
            total += t.phase.as_ref().map_or(0, |p| bytes(&p.liquid));
        }
        if let Some(s) = &self.surface {
            total += s.heap_bytes();
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

    /// Whether `cell` is a fluid cell (a boundary node of the on-node wall-function bounce is
    /// not: it is the wall).
    pub fn is_fluid(&self, cell: usize) -> bool {
        self.collides().get(cell).copied().unwrap_or(false)
    }

    /// The cells that collide: the fluid cells, less the boundary nodes of the on-node
    /// wall-function bounce.
    fn collides(&self) -> &[bool] {
        if self.bulk.is_empty() { &self.fluid } else { &self.bulk }
    }

    /// Whether Han et al.'s on-node wall-function bounce runs (`spalding_bounce`).
    pub fn on_node_bounce(&self) -> bool {
        self.nodes.is_some()
    }

    /// The boundary nodes of the on-node wall-function bounce (`spalding_bounce`).
    pub fn wall_nodes(&self) -> Option<&WallNodes> {
        self.nodes.as_ref()
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

    /// Internal sheets need a paired thermal update that device backends must support.
    pub fn has_baffles(&self) -> bool {
        !self.baffles.is_empty()
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

    /// Smagorinsky's constant of the eddy viscosity (0: none).
    pub fn smagorinsky(&self) -> f64 {
        self.smagorinsky
    }

    /// Vreman's model constant c = 2.5 C_s² (0: not this model).
    pub fn vreman(&self) -> f64 {
        self.vreman
    }

    /// Whether a subgrid model (Smagorinsky's or Vreman's) gives the cells an eddy viscosity.
    pub fn turbulent(&self) -> bool {
        self.smagorinsky > 0.0 || self.vreman > 0.0
    }

    /// With the hybrid collision's dynamic σ, its smallest and mean value over the fluid cells in
    /// the last step (a diagnostic; `None` without dynamic σ or before the first step).
    pub fn dynamic_sigma(&self) -> Option<(f64, f64)> {
        self.sigma_stats
    }

    /// The hybrid collision's σ, fixed or dynamic (`None` with another collision; [`crate::hrr`]).
    pub fn hybrid_sigma(&self) -> Option<fieldcraft_model::HrrSigma> {
        use fieldcraft_model::HrrSigma;
        (self.collision == Collision::Hrr).then_some(if self.dynamic { HrrSigma::Dynamic } else { HrrSigma::Value(self.sigma) })
    }

    /// Whether each axis is periodic: the finite-difference gradient's neighbours across its faces
    /// are the periodic images ([`crate::hrr`]).
    pub fn periodic(&self) -> [bool; 3] {
        self.periodic
    }

    /// The turbulent Prandtl number of Smagorinsky's model (its turbulent diffusivity ν_t / Pr_t).
    pub fn turbulent_prandtl(&self) -> f64 {
        self.turbulent_prandtl
    }

    /// A shear-dependent fluid's relaxation-time model, when the fluid has one.
    pub fn rheology(&self) -> Option<&Rheology> {
        self.rheology.as_ref()
    }

    /// With a shear-dependent fluid: how many cells' relaxation times were clamped.
    pub fn viscosity_clamps(&self) -> Option<Clamps> {
        self.rheology.is_some().then_some(self.clamps)
    }

    /// The largest eddy viscosity over the fluid's own viscosity in the last step (0 without a
    /// turbulence model).
    pub fn eddy_viscosity_ratio(&self) -> f64 {
        self.eddy_ratio
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
        // every fluid cell, the boundary nodes of the on-node bounce included; not a gas cell
        let populated = |cell: usize| {
            self.surface.as_ref().map_or(self.fluid.get(cell).copied().unwrap_or(false), |s| s.active.get(cell).copied().unwrap_or(false))
        };
        for cell in 0..n {
            if !populated(cell) {
                continue;
            }
            let (rho, u) = state(cell);
            let a = self.force_at(cell);
            let mut force = a.map(|c| rho * c);
            if let Some(z) = cell_action(self.action_at(cell), self.liquid_fraction(cell)) {
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
        if !self.populated().get(cell).copied().unwrap_or(false) {
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
        let action = cell_action(self.action_at(cell), self.liquid_fraction(cell));
        let (u, _) = velocity_force(rho, j, self.force_at(cell), action.as_ref());
        Some((rho, u))
    }

    /// The total mass of the fluid cells (lattice units: the sum of every population; with a free
    /// surface the liquid cells' densities, the interface cells' masses and the reservoir).
    pub fn mass(&self) -> f64 {
        if let Some(surface) = &self.surface {
            return surface.total_mass(&self.f, self.stride);
        }
        let n = self.n;
        let s = self.stride;
        let mut m = 0.0;
        for cell in (0..n).filter(|c| self.is_fluid(*c)) {
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
        self.refresh_damping();
        self.sample_walls();
        self.sample_stress();
        self.sample_nodes();
        self.gradient_velocities();
        let max_u2 = self.collide();
        self.collide_nodes();
        if !(max_u2 <= MAX_LATTICE_SPEED * MAX_LATTICE_SPEED) {
            return Err(Error::Unstable(format!(
                "the flow went unstable at time step {}: the lattice speed reached {:.3} (limit {MAX_LATTICE_SPEED}) or a density stopped being a positive number; lower the time step, refine the spacing or check the boundary conditions",
                self.steps_done + 1,
                max_u2.sqrt()
            )));
        }
        self.stream();
        self.apply_stress();
        self.surface_stream();
        self.open_boundaries();
        std::mem::swap(&mut self.f, &mut self.g);
        self.surface_convert();
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
        // With an eddy viscosity, each cell's turbulent diffusivity α_t = ν_t / Pr_t adds to
        // τ_T = 4 α + ½ (the scalar lattice's c_s² = 1/4).
        let (cell_nu_t, prandtl) = (&self.cell_nu_t, self.turbulent_prandtl);
        let Some(t) = &mut self.thermal else { return };
        let (wp, wm) = (t.omega_plus, t.omega_minus);
        let phase = t.phase.as_ref();
        let mut blocks = cell_blocks(&mut t.g, s, n, block);
        parallel::for_each_mut(&mut blocks, |b, g| {
            let len = g.first().map_or(0, |d| d.len());
            for i in 0..len {
                let cell = b * block + i;
                if !fluid.get(cell).copied().unwrap_or(false) {
                    continue;
                }
                // the populations' sum: the temperature, or with a phase change the total enthalpy
                let mut e = 0.0;
                for q in 0..QT {
                    e += g.get(q).and_then(|d| d.get(i)).copied().unwrap_or(0.0);
                }
                let temp = phase.map_or(e, |p| p.rule.split(e).0);
                let u = cell_u.get(cell).copied().unwrap_or([0.0; 3]);
                // The cell's turbulent diffusivity α_t = ν_t / Pr_t raises τ_T by 4 α_t; τ_T⁺
                // follows with Λ_T kept, and a single relaxation time stays single.
                let (wp, wm) = match cell_nu_t.get(cell) {
                    Some(nu_t) if *nu_t > 0.0 && prandtl > 0.0 => {
                        let tau = 1.0 / wm + 4.0 * nu_t / prandtl;
                        if (wp - wm).abs() <= f64::EPSILON * wm.abs() {
                            (1.0 / tau, 1.0 / tau)
                        } else {
                            let magic = (1.0 / wp - 0.5) * (1.0 / wm - 0.5);
                            (1.0 / (0.5 + magic / (tau - 0.5)), 1.0 / tau)
                        }
                    }
                    _ => (wp, wm),
                };
                if let Some(x) = g.get_mut(0).and_then(|d| d.get_mut(i)) {
                    *x -= wp * (*x - equilibrium_e(0, e, temp, u));
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
        for b in &mut self.baffles {
            b.stream(&self.grid, s, g, &mut t.g2, &mut t.heat);
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
                    let e = t.enthalpy(*value);
                    for q in 0..QT {
                        if let Some(x) = t.g2.get_mut(q * s + cell) {
                            *x = equilibrium_e(q, e, *value, u);
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

    /// The temperature (and with a phase change the liquid fraction) of every fluid cell from its
    /// populations; false when one is not finite.
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
            let mut e = 0.0;
            for q in 0..QT {
                e += t.g.get(q * s + cell).copied().unwrap_or(0.0);
            }
            finite &= e.is_finite();
            let (temp, liquid) = t.phase.as_ref().map_or((e, 1.0), |p| p.rule.split(e));
            if let Some(slot) = t.temperature.get_mut(cell) {
                *slot = temp;
            }
            if let Some(slot) = t.phase.as_mut().and_then(|p| p.liquid.get_mut(cell)) {
                *slot = liquid;
            }
        }
        finite
    }

    /// With Van Driest's damping: each wall cell's friction velocity from the flow before the
    /// collision.
    /// With the wall model: each wall cell's modelled velocity and rate from its second cell's
    /// velocity before the collision.
    fn sample_walls(&mut self) {
        if !self.wall_model {
            return;
        }
        let Some(mut walls) = self.walls.take() else { return };
        walls.sample(|cell| self.moments(cell).map_or([0.0; 3], |m| m.1));
        self.walls = Some(walls);
    }

    /// With the wall-function bounce: each wall cell's stress from its own density and velocity
    /// before the collision.
    fn sample_stress(&mut self) {
        if !self.wall_bounce {
            return;
        }
        let Some(mut walls) = self.walls.take() else { return };
        walls.sample_stress(|cell| self.moments(cell).unwrap_or((0.0, [0.0; 3])));
        self.walls = Some(walls);
    }

    /// With the on-node wall-function bounce: each boundary node's stress from its interior node
    /// before the collision.
    fn sample_nodes(&mut self) {
        let Some(mut nodes) = self.nodes.take() else { return };
        nodes.sample(|cell| self.moments(cell));
        self.nodes = Some(nodes);
    }

    /// With the on-node wall-function bounce, after the collision of the fluid cells: each
    /// boundary node's eq. 14 in place of a collision ([`crate::bounce::on_node`]), one node after
    /// another in cell order (a few per cent of the cells; each reads and writes only its own
    /// populations).
    fn collide_nodes(&mut self) {
        let Some(nodes) = &mut self.nodes else { return };
        let s = self.stride;
        let mut fq = [0.0; Q];
        for w in nodes.nodes_mut() {
            let cell = w.cell as usize;
            for (q, x) in fq.iter_mut().enumerate() {
                *x = self.f.get(q * s + cell).copied().unwrap_or(0.0);
            }
            w.taken = crate::bounce::on_node(&mut fq, usize::from(w.axis), i32::from(w.inward), w.stress);
            for (q, v) in fq.iter().enumerate() {
                if let Some(x) = self.f.get_mut(q * s + cell) {
                    *x = *v;
                }
            }
        }
    }

    /// With the wall-function bounce, after streaming: each wall cell's stress taken from its
    /// streamed populations ([`Walls::stresses`]) and given to its wall's force group, so the
    /// walls carry exactly the stress the law asks for.
    fn apply_stress(&mut self) {
        if !self.wall_bounce {
            return;
        }
        let Some(walls) = &self.walls else { return };
        let s = self.stride;
        for (cell, normal, group, stress) in walls.stresses() {
            for (a, tau) in stress.iter().enumerate().filter(|(_, t)| **t != 0.0) {
                // The two populations into the fluid diagonally along axis a: −τ/2 on the one
                // along +a, +τ/2 on the one along −a.
                for (sign, delta) in [(1.0, -0.5 * tau), (-1.0, 0.5 * tau)] {
                    let mut dir = normal;
                    dir[a] += sign;
                    if let Some(x) = (0..Q).find(|q| c(*q) == dir).and_then(|q| self.g.get_mut(q * s + cell)) {
                        *x += delta;
                    }
                }
            }
            if let Some(force) = self.forces.get_mut(group) {
                for a in 0..3 {
                    force[a] += stress[a];
                }
            }
        }
    }

    fn refresh_damping(&mut self) {
        if !self.van_driest() {
            return;
        }
        let Some(mut walls) = self.walls.take() else { return };
        walls.refresh(|cell| self.moments(cell).map_or([0.0; 3], |m| m.1));
        self.walls = Some(walls);
    }

    /// Whether Van Driest's damping scales the eddy viscosity near walls.
    pub fn van_driest(&self) -> bool {
        self.walls.as_ref().is_some_and(Walls::damped)
    }

    /// Whether the wall cells relax with the wall model's viscosity.
    pub fn wall_model(&self) -> bool {
        self.wall_model
    }

    /// Whether the wall-function bounce gives the walls their stress (`spalding_bounce_halfway`).
    pub fn wall_bounce(&self) -> bool {
        self.wall_bounce
    }

    /// The wall cells, with Van Driest's damping or the wall model.
    pub fn walls(&self) -> Option<&Walls> {
        self.walls.as_ref()
    }

    /// With Van Driest's damping: the largest y⁺ of a wall cell's centre in the last step.
    pub fn wall_y_plus(&self) -> Option<f64> {
        self.walls.as_ref().filter(|w| w.damped()).map(Walls::largest_wall_y_plus)
    }

    /// With Van Driest's damping or the wall model: cell `cell`'s distance to its wall (lattice
    /// units), when a wall reaches it (without damping, the wall cells alone).
    pub fn wall_distance(&self, cell: usize) -> Option<f64> {
        self.walls.as_ref().and_then(|d| d.distance(cell))
    }

    /// With the hybrid collision or Vreman's model, before the collision: every fluid cell's velocity (with the half force),
    /// computed exactly as the collision computes it (the wall model's reconstruction included),
    /// into `gradient_u`, by fixed blocks in parallel (each cell writes only its own slot), for the
    /// finite-difference gradient of its neighbours ([`crate::hrr`], [`crate::vreman`]). Solid cells and cells
    /// whose density is not a positive number get 0 (the collision then stops the run).
    fn gradient_velocities(&mut self) {
        if self.gradient_u.is_empty() {
            return;
        }
        let (f, s, block) = (&self.f, self.stride, self.block);
        let fluid = populated_mask(&self.fluid, &self.bulk, self.surface.as_ref());
        let (cell_force, constant) = (&self.cell_force, self.acceleration);
        let (zone_of, actions) = (&self.zone_of, &self.actions);
        let walls = if self.wall_model { self.walls.as_ref() } else { None };
        // field by field: the pass borrows other fields mutably
        let liquid = self.thermal.as_ref().and_then(|t| t.phase.as_ref()).map(|p| &p.liquid[..]);
        parallel::chunks_mut(&mut self.gradient_u, block, |b, out| {
            let mut fq = [0.0; Q];
            for (i, slot) in out.iter_mut().enumerate() {
                let cell = b * block + i;
                *slot = [0.0; 3];
                if !fluid.get(cell).copied().unwrap_or(false) {
                    continue;
                }
                for (q, x) in fq.iter_mut().enumerate() {
                    *x = f.get(q * s + cell).copied().unwrap_or(0.0);
                }
                let a = cell_force.get(cell).copied().unwrap_or(constant);
                let zone = zone_of.get(cell).and_then(|z| actions.get(*z as usize));
                let action = cell_action(zone, liquid.and_then(|l| l.get(cell)).copied());
                if let Some(state) = cell_state(&mut fq, cell, a, walls, action.as_ref()) {
                    *slot = state.u;
                }
            }
        });
    }

    /// The collision with Guo forcing (two relaxation times, or recursive regularised), in place,
    /// by blocks of [`block_cells`] cells in parallel. Returns the largest squared lattice speed, or
    /// NaN when a density is not positive.
    fn collide(&mut self) -> f64 {
        let block = self.block;
        let n = self.n;
        let hybrid = self.collision == Collision::Hrr;
        // The neighbours' velocities of the pass before the collision, for the finite-difference
        // gradient (the hybrid collision, Vreman's model).
        let gradients = !self.gradient_u.is_empty();
        let regularized = self.collision.is_regularized();
        let mrt = self.collision == Collision::Mrt;
        let (sigma, gradient_u, grid, periodic) = (self.sigma, &self.gradient_u, self.grid, self.periodic);
        // Dynamic σ's cap on the von Kármán length: the box's shortest extent in cells.
        let l_max = self.grid.cells.iter().copied().min().unwrap_or(1) as f64;
        let (wp, wm) = (self.omega_plus, self.omega_minus);
        let (sp, sm) = (1.0 - 0.5 * wp, 1.0 - 0.5 * wm);
        let (smagorinsky, vreman) = (self.smagorinsky, self.vreman);
        let turbulent = self.turbulent();
        // Dynamic σ (set only with the hybrid collision), and the fluid's own viscosity (lattice
        // units) it compares the eddy viscosity with.
        let dynamic = self.dynamic;
        let nu0 = (1.0 / wp - 0.5) / 3.0;
        let walls = self.walls.as_ref();
        let wall_model = self.wall_model;
        let model_walls = walls.filter(|_| wall_model);
        let rheology = self.rheology.as_ref();
        // The lattice temperatures a temperature-dependent viscosity reads (the last step's).
        let temperatures: &[f64] = self.thermal.as_ref().map_or(&[], |t| &t.temperature);
        // Λ, which the two relaxation times keep when the eddy viscosity changes τ⁺.
        let magic = (1.0 / wp - 0.5) * (1.0 / wm - 0.5);
        let constant = self.acceleration;
        let cell_force = &self.cell_force;
        // with a phase change, each cell's liquid fraction after the last step (field by field:
        // the pass borrows other fields mutably)
        let liquid = self.thermal.as_ref().and_then(|t| t.phase.as_ref()).map(|p| &p.liquid[..]);
        let forced = constant != [0.0; 3] || !cell_force.is_empty() || !self.zone_of.is_empty() || liquid.is_some();
        // The boundary nodes of the on-node wall-function bounce run eq. 14 instead; the gas
        // cells of a free surface hold no populations.
        let fluid = populated_mask(&self.fluid, &self.bulk, self.surface.as_ref());
        let (zone_of, actions) = (&self.zone_of, &self.actions);
        let mut blocks = cell_blocks(&mut self.f, self.stride, n, block);
        let mut maxima = vec![0.0f64; blocks.len()];
        // The largest eddy viscosity per block (lattice units).
        let mut eddies = vec![0.0f64; blocks.len()];
        // Each block's clamped relaxation times, low and high.
        let mut clamped = vec![[0u64; 2]; blocks.len()];
        // with heat, each block's velocities go to its slice of cell_u (empty otherwise), and
        // with an eddy viscosity too, its eddy viscosities to its slice of cell_nu_t
        let mut u_blocks: Vec<&mut [[f64; 3]]> = self.cell_u.chunks_mut(block).collect();
        let mut nu_blocks: Vec<&mut [f64]> = self.cell_nu_t.chunks_mut(block).collect();
        // With the wall model, each block's rebuilt wall cells and the momentum their
        // reconstruction took, in cell order.
        let mut removed: Vec<Vec<(u32, [f64; 3])>> = vec![Vec::new(); blocks.len()];
        // With dynamic σ, each block's smallest σ, sum of σ and cell count.
        let mut sigmas = vec![[f64::INFINITY, 0.0, 0.0]; blocks.len()];
        type Work<'a, 'b> = (
            &'a mut Vec<&'b mut [f64]>,
            &'a mut f64,
            Option<&'a mut [[f64; 3]]>,
            Option<&'a mut [f64]>,
            &'a mut f64,
            &'a mut Vec<(u32, [f64; 3])>,
            &'a mut [f64; 3],
            &'a mut [u64; 2],
        );
        let mut work: Vec<Work> = Vec::with_capacity(blocks.len());
        let mut u_iter = u_blocks.drain(..);
        let mut nu_iter = nu_blocks.drain(..);
        for (((((f, m), e), r), sg), c) in
            blocks.iter_mut().zip(maxima.iter_mut()).zip(eddies.iter_mut()).zip(removed.iter_mut()).zip(sigmas.iter_mut()).zip(clamped.iter_mut())
        {
            work.push((f, m, u_iter.next(), nu_iter.next(), e, r, sg, c));
        }
        parallel::for_each_mut(&mut work, |b, (f, block_max, u_out, nu_out, max_nu, removed, sigma_stats, block_clamped)| {
            let mut fq = [0.0; Q];
            let mut post = [0.0; Q];
            // The block's largest speed, kept here and stored once: the blocks' maxima share
            // cache lines, so a store per cell bounced those lines between the threads (#643).
            let mut max_u2 = 0.0f64;
            // The block's clamped relaxation times, likewise kept here and stored once.
            let mut cell_clamps = [0u64; 2];
            let len = f.first().map_or(0, |d| d.len());
            for i in 0..len {
                if !fluid.get(b * block + i).copied().unwrap_or(false) {
                    continue;
                }
                for (q, x) in fq.iter_mut().enumerate() {
                    *x = f.get(q).and_then(|d| d.get(i)).copied().unwrap_or(0.0);
                }
                let a = cell_force.get(b * block + i).copied().unwrap_or(constant);
                let zone = zone_of.get(b * block + i).and_then(|z| actions.get(*z as usize));
                let action = cell_action(zone, liquid.and_then(|l| l.get(b * block + i)).copied());
                let Some(CellState { rho, u, force, wall }) = cell_state(&mut fq, b * block + i, a, model_walls, action.as_ref()) else {
                    **block_max = f64::NAN;
                    return;
                };
                // A wall cell under the wall model relaxes with its own rate.
                let wall_rate = wall.map(|(rate, taken)| {
                    removed.push(((b * block + i) as u32, taken));
                    rate
                });
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
                // The cell's eddy viscosity (from the rate it relaxes with, c_s² = 1/3, but with
                // dynamic σ, where the model's does not enter the rate).
                let from_rate = |w: f64| ((1.0 / w - 1.0 / wp) / 3.0).max(0.0);
                let mut eddy = |nu_t: f64| {
                    **max_nu = max_nu.max(nu_t);
                    if let Some(slot) = nu_out.as_mut().and_then(|v| v.get_mut(i)) {
                        *slot = nu_t;
                    }
                };
                // Smagorinsky's constant of the cell, damped near walls.
                let cell_c = walls.map_or(smagorinsky, |d| smagorinsky * d.factor(b * block + i));
                // The finite-difference velocity gradient, from the neighbours' velocities of the
                // pass before the collision.
                let local = gradients.then(|| {
                    let p = grid.coords(b * block + i);
                    let nb: crate::hrr::Neighbours = [0, 1, 2]
                        .map(|axis| [false, true].map(|up| neighbour(&grid, periodic, fluid, p, axis, up).and_then(|c| gradient_u.get(c).copied())));
                    (crate::hrr::gradient(u, &nb), nb)
                });
                let grad = local.as_ref().map(|(g, _)| g);
                // Vreman's explicit eddy viscosity (not at a wall cell under the wall model, which
                // keeps its own rate), and the rate it gives, τ̄ = 3 (ν + ν_t) + ½; none where
                // ν_t = 0, so the cell keeps the fluid's own rates bit for bit.
                let nu_vreman = match &grad {
                    Some(g) if vreman > 0.0 && wall_rate.is_none() => crate::vreman::eddy_viscosity(g, vreman),
                    _ => 0.0,
                };
                let explicit = (nu_vreman > 0.0).then(|| 1.0 / (1.0 / wp + 3.0 * nu_vreman));
                // Smagorinsky's implicit rate of the cell's own flux (Hou et al.).
                let flux = || crate::turbulence::flux(&fq, rho, u, force, forced);
                let implicit = || crate::turbulence::omega(&flux(), rho, wp, cell_c);
                // A shear-dependent fluid's rate from the cell's flux, counting a clamped one.
                let mut own_rate = |r: &Rheology, p: &[f64; 6]| {
                    let (tau, clamp) = r.relax(p, rho, temperatures.get(b * block + i).copied().unwrap_or(0.0));
                    match clamp {
                        Clamp::None => {}
                        Clamp::Low => cell_clamps[0] += 1,
                        Clamp::High => cell_clamps[1] += 1,
                    }
                    1.0 / tau
                };
                // A wall cell under the wall model relaxes with its viscosity instead.
                if regularized {
                    // With dynamic σ the subgrid model sets σ and stays out of the rate: the
                    // cell relaxes with the fluid's own (or the wall model's).
                    let (rate, c) = match wall_rate {
                        Some(w) => (w, 0.0),
                        None if dynamic => (wp, 0.0),
                        None => (explicit.unwrap_or(wp), cell_c),
                    };
                    // The hybrid collision's σ and finite-difference strain rate. Dynamic: σ from
                    // the model's eddy viscosity (Vreman's, or Hou's implicit Smagorinsky one),
                    // 1 at a wall cell under the wall model.
                    let mut nu_model = None;
                    let blend = local.as_ref().filter(|_| hybrid).map(|(g, nb)| {
                        let cell_sigma = if dynamic {
                            let s = match wall_rate {
                                Some(_) => 1.0,
                                None => {
                                    // Smagorinsky's ν_t itself, 0 below the flux's round-off
                                    // (crate::turbulence, "Round-off").
                                    let nu_t = if smagorinsky > 0.0 {
                                        let p = crate::turbulence::flux(&fq, rho, u, force, forced);
                                        crate::turbulence::eddy_viscosity(&p, rho, wp, cell_c, crate::turbulence::flux_round_off(rho))
                                    } else {
                                        nu_vreman
                                    };
                                    nu_model = Some(nu_t);
                                    crate::hrr::dynamic_sigma(g, &crate::hrr::laplacian(u, nb), nu_t, nu0, l_max)
                                }
                            };
                            sigma_stats[0] = sigma_stats[0].min(s);
                            sigma_stats[1] += s;
                            sigma_stats[2] += 1.0;
                            s
                        } else {
                            sigma
                        };
                        (cell_sigma, crate::hrr::strain(g))
                    });
                    // The cell's rate: a shear-dependent fluid's of its flux (not at a wall cell
                    // under the wall model), Smagorinsky's, or the one above.
                    let w =
                        crate::regularized::rebuild(
                            &fq,
                            rho,
                            u,
                            force,
                            forced,
                            blend.as_ref().map(|(s, st)| (*s, st)),
                            &mut post,
                            |p| match rheology.filter(|_| wall_rate.is_none()) {
                                Some(r) => own_rate(r, p),
                                None if c > 0.0 => crate::turbulence::omega(p, rho, rate, c),
                                None => rate,
                            },
                        );
                    if turbulent {
                        eddy(nu_model.unwrap_or_else(|| from_rate(w)));
                    }
                    for (q, v) in post.iter().enumerate() {
                        put(q, *v);
                    }
                    continue;
                }
                if mrt {
                    // The stresses' rate: the wall model's, Smagorinsky's of the cell's flux, a
                    // shear-dependent fluid's, or the fluid's.
                    let rate = match wall_rate {
                        Some(w) => w,
                        None if smagorinsky > 0.0 => implicit(),
                        None => match rheology {
                            Some(r) => own_rate(r, &flux()),
                            None => explicit.unwrap_or(wp),
                        },
                    };
                    if turbulent {
                        eddy(from_rate(rate));
                    }
                    crate::mrt::collide(&fq, rho, u, force, rate, forced, &mut post);
                    for (q, v) in post.iter().enumerate() {
                        put(q, *v);
                    }
                    continue;
                }
                // The cell's own rates where the eddy viscosity or a shear-dependent fluid sets τ⁺
                // (τ⁻ follows, keeping Λ).
                let cell_rate =
                    wall_rate.or_else(|| (smagorinsky > 0.0).then(implicit)).or(explicit).or_else(|| rheology.map(|r| own_rate(r, &flux())));
                if turbulent {
                    eddy(from_rate(cell_rate.unwrap_or(wp)));
                }
                let (wp, wm, sp, sm) = match cell_rate {
                    Some(wp) => {
                        let wm = 1.0 / (0.5 + magic / (1.0 / wp - 0.5));
                        (wp, wm, 1.0 - 0.5 * wp, 1.0 - 0.5 * wm)
                    }
                    None => (wp, wm, sp, sm),
                };
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
            **block_clamped = cell_clamps;
        });
        drop(work);
        if wall_model && let Some(walls) = &mut self.walls {
            walls.set_removed(&removed.concat());
        }
        // With dynamic σ, its smallest and mean value (the blocks' sums added in block order).
        self.sigma_stats = dynamic.then(|| {
            let (min, sum, count) = sigmas.iter().fold((f64::INFINITY, 0.0, 0.0), |(m, s, c), x| (m.min(x[0]), s + x[1], c + x[2]));
            (min, if count > 0.0 { sum / count } else { 1.0 })
        });
        // The largest eddy viscosity over the fluid's own (a maximum: no order dependence).
        self.eddy_ratio = if turbulent && nu0 > 0.0 { eddies.iter().fold(0.0f64, |m, x| m.max(*x)) / nu0 } else { 0.0 };
        if self.rheology.is_some() {
            // Counts: integer sums, the same in any order.
            let last = clamped.iter().fold([0u64; 2], |s, c| [s[0] + c[0], s[1] + c[1]]);
            self.clamps.add(last);
        }
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
        // The on-node wall-function bounce's share: the momentum each boundary node took from the
        // fluid, in cell order.
        if let Some(nodes) = &self.nodes {
            for w in nodes.nodes() {
                if let Some(force) = self.forces.get_mut(w.group as usize) {
                    for a in 0..3 {
                        force[a] += w.taken[a];
                    }
                }
            }
        }
        // The wall model's share: the momentum each reconstruction took from its wall cell, in
        // cell order.
        if self.wall_model
            && let Some(walls) = &self.walls
        {
            for (group, j) in walls.removed() {
                if let Some(force) = self.forces.get_mut(group) {
                    for a in 0..3 {
                        force[a] += j[a];
                    }
                }
            }
        }
    }

    /// Zou–He velocity, extrapolated pressure and convective boundaries. Every inner-cell read
    /// uses the streamed populations before any open cell is changed, so cells are worked out
    /// in parallel (blocks of [`OPEN_BLOCK`]) and
    /// stored afterwards in cell order: done one after another they were a tenth of a step on
    /// one thread, which capped the speed-up of many (#643).
    fn open_boundaries(&mut self) {
        let s = self.stride;
        let regularized = self.collision.is_regularized();
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
                    // Guo, Zheng & Shi's non-equilibrium extrapolation. Read the streamed
                    // inner cell before any open boundary is written, including at edges.
                    if let Target::Density(rho) = oc.target {
                        let step = [1, nx, nx * ny].get(axis).copied().unwrap_or(0);
                        let inside = if sign > 0.0 { cell.checked_add(step) } else { cell.checked_sub(step) };
                        if let Some(inner) = inside.filter(|&i| flow.is_fluid(i)) {
                            let fi = std::array::from_fn(|q| flow.g.get(q * s + inner).copied().unwrap_or(0.0));
                            let mut inner_half = flow.force_at(inner).map(|v| 0.5 * v);
                            let action = cell_action(flow.action_at(inner), flow.liquid_fraction(inner));
                            if action.is_some() {
                                let inner_rho: f64 = fi.iter().sum();
                                let mut j = [0.0; 3];
                                for (q, f) in fi.iter().enumerate() {
                                    for (ja, ca) in j.iter_mut().zip(c(q)) {
                                        *ja += f * ca;
                                    }
                                }
                                let (_, force) = velocity_force(inner_rho, j, flow.force_at(inner), action.as_ref());
                                inner_half = force.map(|f| 0.5 * f / inner_rho);
                            }
                            *p = crate::outlet::extrapolate_pressure(fi, rho, axis, inner_half, a_half);
                            continue;
                        }
                    }
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
        // With a free surface: the mass the inlets inject, the incoming populations the boundary
        // set less the outgoing ones its placeholders held (which left through the face).
        let mut injected = 0.0;
        for (k, p) in done.iter().enumerate() {
            let Some(oc) = self.open.get(k).copied() else { continue };
            let cell = oc.cell as usize;
            for (q, v) in p.iter().enumerate() {
                if let Some(x) = self.g.get_mut(q * s + cell) {
                    injected += *v - *x;
                    *x = *v;
                }
            }
            if matches!(oc.target, Target::Convective(_)) {
                outlets.record(k, &oc, self.force_at(cell).map(|v| 0.5 * v), |q| p.get(q).copied().unwrap_or(0.0));
            }
        }
        self.open_scratch = done;
        self.outlets = outlets;
        if let Some(surface) = &mut self.surface {
            surface.injected += injected;
        }
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
/// density and that normal velocity: the fallback when a pressure cell has no inner fluid cell,
/// and a convective outlet's cell where the flow comes back in ([`crate::outlet::backflow`]).
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
/// velocity without the half force): the regularised collision's legacy local pressure fallback
/// (also used for convective backflow). `g` holds
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
    fn baffles(&self) -> Vec<crate::baffle::BaffleResult> {
        self.baffles.iter().map(|b| b.result.clone()).collect()
    }

    fn checkpoint(&mut self, a: &mut dyn crate::checkpoint::Archive) -> Result<()> {
        use crate::checkpoint::{boolean, doubles, f64, fixed, u64, vectors};
        // Streaming targets are allocated once, not cloned into a separate snapshot.
        self.targets()?;
        // Every field is named here, without `..`: a new field does not compile until it is
        // saved or said to be fixed by the setup or rebuilt in each step, and a saved one changes
        // the layout tag in `crate::checkpoint`.
        let Flow {
            grid,
            n,
            stride,
            fluid,
            f,
            g,
            links,
            open: _,
            omega_plus,
            omega_minus,
            collision,
            smagorinsky,
            turbulent_prandtl,
            vreman,
            cell_nu_t,
            eddy_ratio,
            walls,
            wall_model,
            wall_bounce,
            nodes,
            bulk: _,
            rheology,
            clamps,
            acceleration,
            forces,
            groups,
            offset: _,
            fluid_cells: _,
            steps_done,
            empty_obstacles: _,
            // `run_checkpointed` refuses them until the sides have a restart layout.
            baffles: _,
            thermal,
            cell_force,
            cell_u,
            sigma,
            // Each cell's velocity before the collision, rebuilt from `f` at the start of a step.
            gradient_u: _,
            dynamic,
            sigma_stats,
            periodic,
            zones,
            zone_cells: _,
            zone_of: _,
            actions,
            fan_rise,
            outlets,
            block: _,
            // Overwritten before use; its capacity is not state.
            open_scratch: _,
            surface,
        } = self;
        let n = *n;
        fixed(a, n as u64)?;
        fixed(a, *stride as u64)?;
        for v in grid.cells {
            fixed(a, v as u64)?;
        }
        for v in grid
            .origin
            .into_iter()
            .chain([grid.spacing, *omega_plus, *omega_minus, *smagorinsky, *turbulent_prandtl, *vreman, *sigma])
            .chain(*acceleration)
        {
            fixed(a, v.to_bits())?;
        }
        let collision = match collision {
            Collision::Trt => 0,
            Collision::Regularized => 1,
            Collision::Mrt => 2,
            Collision::Hrr => 3,
        };
        for flag in [collision, u64::from(*wall_model), u64::from(*wall_bounce), u64::from(*dynamic), u64::from(rheology.is_some())]
            .into_iter()
            .chain(periodic.map(u64::from))
        {
            fixed(a, flag)?;
        }
        for v in fluid.iter() {
            fixed(a, u64::from(*v))?;
        }
        fixed(a, links.len() as u64)?;
        for l in links.iter() {
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
        u64(a, steps_done)?;
        let len = f.len();
        doubles(a, f, len)?;
        doubles(a, g, len)?;
        let nt = cell_nu_t.len();
        doubles(a, cell_nu_t, nt)?;
        f64(a, eddy_ratio)?;
        if *eddy_ratio < 0.0 || cell_nu_t.iter().any(|nu| *nu < 0.0) {
            return Err(Error::Invalid("checkpoint eddy viscosity is negative".into()));
        }
        // The last collision's dynamic σ, reported at the end of the run.
        let mut stats = sigma_stats.is_some();
        boolean(a, &mut stats)?;
        if stats != *dynamic {
            return Err(Error::Invalid("checkpoint dynamic σ statistics differ from this flow".into()));
        }
        if stats {
            let (mut least, mut mean) = sigma_stats.unwrap_or_default();
            f64(a, &mut least)?;
            f64(a, &mut mean)?;
            *sigma_stats = Some((least, mean));
        }
        // The relaxation-time clamps of a shear-dependent fluid over the run.
        for v in clamps.last.iter_mut().chain(clamps.total.iter_mut()).chain([&mut clamps.collisions]) {
            u64(a, v)?;
        }
        if clamps.collisions > *steps_done {
            return Err(Error::Invalid("checkpoint viscosity clamps count more collisions than steps".into()));
        }
        fixed(a, u64::from(walls.is_some()))?;
        if let Some(w) = walls {
            w.checkpoint(a)?;
        }
        fixed(a, u64::from(nodes.is_some()))?;
        if let Some(w) = nodes {
            w.checkpoint(a)?;
        }
        vectors(a, forces, groups.len())?;
        let nf = cell_force.len();
        vectors(a, cell_force, nf)?;
        let nu = cell_u.len();
        vectors(a, cell_u, nu)?;
        doubles(a, fan_rise, zones.len())?;
        for action in actions.iter_mut() {
            for v in &mut action.push {
                f64(a, v)?;
            }
            f64(a, &mut action.alpha)?;
            f64(a, &mut action.beta)?;
        }
        outlets.checkpoint(a)?;
        fixed(a, u64::from(surface.is_some()))?;
        if let Some(s) = surface {
            s.checkpoint(a, n)?;
        }
        fixed(a, u64::from(thermal.is_some()))?;
        if let Some(t) = thermal {
            let ThermalState {
                g,
                g2,
                temperature,
                omega_plus,
                omega_minus,
                initial,
                buoyancy,
                groups: _,
                open: _,
                reset: _,
                heat,
                solid,
                conjugate: _,
                solid_finite,
                phase,
            } = t;
            for v in [*omega_plus, *omega_minus, *initial].into_iter().chain(*buoyancy) {
                fixed(a, v.to_bits())?;
            }
            let nt = g.len();
            doubles(a, g, nt)?;
            doubles(a, g2, nt)?;
            doubles(a, temperature, n)?;
            doubles(a, heat, groups.len())?;
            boolean(a, solid_finite)?;
            if !*solid_finite {
                return Err(Error::Invalid("cannot resume a checkpoint with unstable solid temperatures".into()));
            }
            fixed(a, u64::from(solid.is_some()))?;
            if let Some(s) = solid {
                s.checkpoint(a)?;
            }
            // The liquid fractions of the last step, which the next collision's resistance reads.
            fixed(a, u64::from(phase.is_some()))?;
            if let Some(Phase { rule, liquid }) = phase {
                for v in [rule.latent, rule.solidus, rule.liquidus] {
                    fixed(a, v.to_bits())?;
                }
                doubles(a, liquid, n)?;
                if liquid.iter().any(|f| !(0.0..=1.0).contains(f)) {
                    return Err(Error::Invalid("checkpoint liquid fraction is outside 0 to 1".into()));
                }
            }
        }
        Ok(())
    }

    fn name(&self) -> String {
        "cpu".into()
    }

    fn mean_velocity(&mut self, axis: usize) -> Result<f64> {
        if axis > 2 {
            return Err(Error::Invalid(format!("axis {axis} is not x, y or z")));
        }
        // In cell order, so the sum is the same on every thread count.
        let sum: f64 = (0..self.n).filter_map(|c| self.moments(c)).map(|(_, u)| u.get(axis).copied().unwrap_or(0.0)).sum();
        Ok(sum / self.fluid_cells.max(1) as f64)
    }

    fn set_acceleration(&mut self, a: [f64; 3]) -> Result<()> {
        if !a.iter().all(|v| v.is_finite()) {
            return Err(Error::Invalid(format!("the body acceleration {a:?} is not a number")));
        }
        if !self.cell_force.is_empty() {
            return Err(Error::Invalid("the body acceleration cannot change under buoyancy, which sets each cell's force".into()));
        }
        self.acceleration = a;
        Ok(())
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
        // A cell the wall model rebuilds reports the velocity its reconstruction imposes: its
        // populations between steps are the streamed state the next reconstruction overwrites.
        let modelled =
            |cell: usize| self.walls.as_ref().filter(|_| self.wall_model && self.steps_done > 0).and_then(|w| w.modelled(cell)).map(|m| m.0);
        for cell in 0..n {
            if let (Some((rho, u)), Some(d), Some(v)) = (self.moments(cell), density.get_mut(cell), velocity.get_mut(cell)) {
                *d = rho;
                *v = modelled(cell).unwrap_or(u);
            }
        }
        let temperature = self.thermal.as_ref().map(|t| t.temperature.clone()).unwrap_or_default();
        let mut fill = Vec::new();
        if let Some(s) = &self.surface {
            // a gas cell is at the gas pressure
            for (d, k) in density.iter_mut().zip(&s.kind) {
                if *k == Kind::Gas {
                    *d = crate::surface::GAS_DENSITY;
                }
            }
            fill = s.fill.clone();
        }
        Ok(crate::backend::Fields { density, velocity, temperature, fill })
    }

    fn heat(&mut self) -> Result<Vec<f64>> {
        Ok(Flow::heat(self).to_vec())
    }

    fn radiated(&mut self) -> Result<Vec<f64>> {
        Ok(self.thermal.as_ref().and_then(|t| t.solid.as_ref()).map_or_else(Vec::new, |s| s.radiated.clone()))
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

    fn surface_report(&self) -> Option<SurfaceReport> {
        Flow::surface_report(self)
    }

    fn liquid_extent(&self) -> Option<(Option<([usize; 3], [usize; 3])>, f64)> {
        self.surface.as_ref().map(|s| s.extent(&self.grid))
    }

    fn frozen(&self) -> Option<f64> {
        let phase = self.phase()?;
        let solid: f64 = phase.liquid.iter().zip(&self.fluid).filter(|(_, fluid)| **fluid).map(|(f, _)| 1.0 - f).sum();
        Some(solid / self.fluid_cells.max(1) as f64)
    }

    fn eddy_viscosity_ratio(&self) -> f64 {
        Flow::eddy_viscosity_ratio(self)
    }

    fn wall_y_plus(&self) -> Option<f64> {
        Flow::wall_y_plus(self)
    }

    fn viscosity_clamps(&self) -> Option<Clamps> {
        Flow::viscosity_clamps(self)
    }

    fn empty_obstacles(&self) -> &[String] {
        &self.empty_obstacles
    }

    fn conducting(&self) -> &[Vec<u32>] {
        self.thermal.as_ref().and_then(|t| t.solid.as_ref()).map_or(&[], |s| &s.by_solid)
    }
}
