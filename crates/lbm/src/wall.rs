//! The fluid next to walls under Smagorinsky's model: Van Driest's damping of the length scale,
//! and a wall model for grids that do not resolve the wall.
//!
//! Without damping Smagorinsky's eddy viscosity `(C Δx)² |S|` is largest where the shear is,
//! at the wall, where the turbulence it stands for is suppressed: the friction comes out too
//! high (`docs/verification.md`, "Long studies"). Van Driest's mixing length damps it over the
//! viscous sublayer, `l = C Δx (1 − exp(−y⁺/A⁺))`, with `y⁺ = y u_τ / ν` and A⁺ = 26
//! (E. R. Van Driest, *On turbulent flow near a wall*, J. Aeronaut. Sci. 23(11):1007–1011,
//! 1956), as Smagorinsky large-eddy simulations of channels have used it since
//! (P. Moin, J. Kim, *Numerical investigation of turbulent channel flow*, J. Fluid Mech.
//! 118:341–377, 1982).
//!
//! **Wall distance.** Every fluid cell with a no-slip bounce-back link is a wall cell: its wall
//! lies halfway along the link, and the nearest of those points (an axis link first) is its wall
//! point and the wall's normal. Each other fluid cell takes the nearest wall point of its 26
//! neighbours, spread outward in cell order until no distance shrinks (vector propagation,
//! P.-E. Danielsson, *Euclidean distance mapping*, Comput. Graph. Image Process. 14:227–248, 1980):
//! exact for plane walls, within a fraction of a cell elsewhere. The spreading does not cross
//! periodic faces, so near one a distance may come out longer than across it.
//!
//! **Friction velocity.** Before each collision each wall cell's `u_τ` follows from the speed
//! parallel to its wall (relative to the wall's own), at its distance `y`, through Spalding's law
//! of the wall, `y⁺ = u⁺ + e^{−κB} (e^{κu⁺} − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6)`, κ = 0.4, B = 5.5
//! (D. B. Spalding, *A single formula for the "law of the wall"*, J. Appl. Mech. 28(3):455–458,
//! 1961): the viscous sublayer `u⁺ = y⁺` when the wall cell lies in it, the log law above, so a
//! coarse grid's first cell reads the wall shear too. Each cell is damped with the `u_τ` of the
//! wall cell its wall point belongs to. A still wall has `u_τ = 0`, which damps the eddy
//! viscosity fully.
//!
//! **Wall model.** On a grid whose first cell lies above the viscous sublayer the bounce-back
//! cannot carry the wall stress, and damping cannot help. A wall model takes the near-wall
//! velocity from a law of the wall instead (U. Piomelli, E. Balaras, *Wall-layer models for
//! large-eddy simulations*, Annu. Rev. Fluid Mech. 34:349–374, 2002). This is the lattice
//! Boltzmann wall model of O. Malaspinas, P. Sagaut, *Wall model for large-eddy simulation based
//! on the lattice Boltzmann method*, J. Comput. Phys. 275:25–40, 2014 (doi:10.1016/j.jcp.2014.06.020),
//! §4.2–4.3. Before each collision, for each wall cell x₁ with the wall half a link below it
//! and the second cell x₂ one link further along the normal (distances y₁ and y₂ = 3 y₁):
//! 1. the velocity at x₂, relative to the wall, gives the direction e parallel to the wall and
//!    the speed U₂ there; u_τ solves U₂ = ũ(y₂; u_τ), ũ Musker's profile;
//! 2. x₁'s velocity is set to ũ(y₁; u_τ) e (less half the body force, for the stored moments);
//! 3. x₁'s density is the one its known populations give with that velocity (the Zou–He
//!    relation, ρ = (ρ₀ + 2 ρ₊) / (1 + u·n_out)), and its non-equilibrium stress the one its
//!    known populations carry, the unknown ones (from the wall) taken from their opposites;
//! 4. all of x₁'s populations are rebuilt as the equilibrium at that density and velocity plus
//!    the regularised stress, and x₁ relaxes with ν + ν_t, ν_t the mixing length
//!    (κ y₁ (1 − exp(−y₁⁺/A⁺)))² |∂ũ/∂y| at y₁ (κ = 0.384, A⁺ = 26), not Smagorinsky's.
//!
//! ũ is Musker's profile, as in the paper ([`musker`]). Where this differs from the paper: a
//! moving wall is taken relative to its velocity (the paper's walls are still and flat). A wall
//! cell whose second cell is not fluid inside the box keeps the bounce-back alone.
//!
//! **Wall force.** The force reported on a modelled wall is the momentum the flow hands to it,
//! the drag the simulation carries, not the stress ρ u_τ² the law infers at x₂. It has two
//! parts: the bounce-back's momentum exchange on the wall cell's links, taken from the
//! populations after the rebuilt cell's collision, and the momentum the reconstruction itself
//! takes from the cell, `Σ c f` streamed less `Σ c f` rebuilt. Interior streaming conserves
//! momentum and the collision adds ρ a whole, so the walls' forces close the fluid's momentum
//! budget each step, and at steady state they carry the body force. A wall cell with links into
//! more than one force group (an edge or corner) puts all of its reconstruction's share into
//! the group of its nearest link.
//!
//! No code of any lattice Boltzmann or LES package was read; this is written from the papers.

use fieldcraft_parallel as parallel;

use crate::Result;
use crate::flow::{Link, Rule, filled};
use crate::lattice::{CS2, Q, W, c, equilibrium, opposite};
use crate::voxel::Grid;

/// No wall cell.
const NONE: u32 = u32::MAX;

/// The largest y⁺ of a wall cell's centre for the wall to count as resolved: the viscous
/// sublayer, where Van Driest's damping is meant to act. Above it the run is warned about.
pub const RESOLVED_WALL_Y_PLUS: f64 = 5.0;

/// Kármán's constant and the log law's intercept of Spalding's formula.
const KAPPA: f64 = 0.4;
const B: f64 = 5.5;

/// The most times the spreading may lower each cell's distance on average before it is stopped.
const MAX_UPDATES_PER_CELL: usize = 64;

/// Newton steps (each guarded by bisection) that solve Spalding's law for u⁺: bisection alone
/// narrows the bracket to the last bit in about 60.
const SPALDING_ITERATIONS: usize = 100;

/// Past this y⁺/A⁺, `1 − exp(−y⁺/A⁺)` rounds to 1 in double precision (from about 37.4).
const UNDAMPED: f64 = 40.0;

/// The 26 neighbours of a cell.
const NEIGHBOURS: [[i64; 3]; 26] = {
    let mut out = [[0; 3]; 26];
    let mut k = 0;
    let mut i = 0;
    while i < 27 {
        let o = [i % 3 - 1, i / 3 % 3 - 1, i / 9 - 1];
        if i != 13 {
            out[k] = o;
            k += 1;
        }
        i += 1;
    }
    out
};

/// The von Kármán constant and Van Driest's A⁺ of the wall cell's mixing-length eddy viscosity
/// in the wall model (the values Malaspinas and Sagaut use).
const MIXING_KAPPA: f64 = 0.384;
const MIXING_A_PLUS: f64 = 26.0;

/// A fluid cell next to a no-slip wall.
#[derive(Debug, Clone, Copy, PartialEq)]
struct WallCell {
    cell: u32,
    /// The direction of its nearest wall link (pointing from the wall into the fluid).
    link: u8,
    /// The wall point (lattice coordinates: cell `[i, j, k]` is centred on `[i, j, k]`).
    point: [f64; 3],
    /// The unit normal from the wall point into the fluid, and the cell's distance along it.
    normal: [f64; 3],
    distance: f64,
    /// The wall's velocity (lattice units).
    velocity: [f64; 3],
    /// The force group of its nearest wall link, which takes the momentum its reconstruction
    /// removes.
    group: u32,
    /// The momentum the wall model's reconstruction took from the cell in the last collision
    /// (lattice units).
    removed: [f64; 3],
    /// The friction velocity at the last refresh (lattice units; with damping).
    friction: f64,
    /// The wall model's second cell off the wall, one link further along the normal (`NONE`
    /// when that is not a fluid cell inside the box).
    second: u32,
    /// The wall model's velocity and relaxation rate of this cell, from the last sample.
    boundary: [f64; 3],
    rate: f64,
}

impl WallCell {
    /// The speed of `u` parallel to the wall, relative to the wall's own velocity.
    fn tangential_speed(&self, u: [f64; 3]) -> f64 {
        let rel = [0, 1, 2].map(|a| u[a] - self.velocity[a]);
        let un = rel[0] * self.normal[0] + rel[1] * self.normal[1] + rel[2] * self.normal[2];
        let t = [0, 1, 2].map(|a| rel[a] - un * self.normal[a]);
        (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt()
    }
}

/// The fluid cells next to no-slip walls, for Van Driest's damping (with each fluid cell's
/// distance to the nearest wall and the wall cell that gives its friction velocity) and for the
/// wall model.
#[derive(Debug, Clone, PartialEq)]
pub struct Walls {
    /// A⁺ with damping, and the fluid's kinematic viscosity (lattice units).
    a_plus: Option<f64>,
    nu: f64,
    walls: Vec<WallCell>,
    /// Per cell: its index in `walls` (`NONE`: no wall reached, or not fluid). Without damping
    /// only the wall cells have one.
    wall_of: Vec<u32>,
    /// With damping, per cell: the distance from its centre to its wall point (lattice units);
    /// empty without.
    distance: Vec<f64>,
}

impl Walls {
    /// The wall cells of `grid`'s fluid, for a fluid of kinematic viscosity `nu` (lattice
    /// units), and with Van Driest's `a_plus` every fluid cell's wall distance. `links` are the
    /// flow's links; `walls[g]` is the velocity of force group `g`'s wall when it is no-slip
    /// (`None`: slip, which has no friction).
    pub fn new(grid: &Grid, fluid: &[bool], links: &[Link], walls: &[Option<[f64; 3]>], a_plus: Option<f64>, nu: f64) -> Result<Walls> {
        let n = grid.len();
        let mut found: Vec<WallCell> = Vec::new();
        let mut wall_of = filled(n, NONE, "wall cells")?;
        // The links are in cell order, so each wall cell's links are consecutive.
        for links in links.chunk_by(|a, b| a.cell == b.cell) {
            let Some(cell) = links.first().map(|l| l.cell) else { continue };
            // The nearest wall point: the population of direction q arrives from −c_q, so the
            // wall is halfway there; an axis link first, then the lowest direction.
            let mut best: Option<(f64, usize, [f64; 3], u32)> = None;
            for l in links {
                let Rule::Bounce { group, .. } = l.rule else { continue };
                let Some(Some(velocity)) = walls.get(group as usize) else { continue };
                let q = usize::from(l.q);
                let cq = c(q);
                let len = 0.5 * (cq[0] * cq[0] + cq[1] * cq[1] + cq[2] * cq[2]).sqrt();
                if best.is_none_or(|(d, b, _, _)| len < d || (len == d && q < b)) {
                    best = Some((len, q, *velocity, group));
                }
            }
            let Some((d, q, velocity, group)) = best else { continue };
            let ijk = grid.coords(cell as usize);
            let centre = ijk.map(|v| v as f64);
            let normal = c(q).map(|v| v * 0.5 / d);
            let point = [0, 1, 2].map(|a| centre[a] - normal[a] * d);
            if let Some(w) = wall_of.get_mut(cell as usize) {
                *w = found.len() as u32;
            }
            let step = c(q);
            let next = [0, 1, 2].map(|a| ijk[a] as i64 + step[a] as i64);
            let inside = (0..3).all(|a| next[a] >= 0 && (next[a] as usize) < grid.cells[a]);
            let second = inside.then(|| grid.index(next.map(|v| v as usize))).filter(|s| fluid.get(*s).copied().unwrap_or(false));
            found.push(WallCell {
                cell,
                link: q as u8,
                point,
                normal,
                distance: d,
                velocity,
                group,
                removed: [0.0; 3],
                friction: 0.0,
                second: second.map_or(NONE, |s| s as u32),
                boundary: velocity,
                rate: 1.0 / (3.0 * nu + 0.5),
            });
        }
        let distance = if a_plus.is_some() { spread(grid, fluid, &found, &mut wall_of)? } else { Vec::new() };
        Ok(Walls { a_plus, nu, walls: found, wall_of, distance })
    }

    /// The wall cells, in cell order: each one's cell, the direction of its nearest wall link,
    /// its distance to its wall (lattice units), its wall's velocity (lattice units) and the wall
    /// model's second cell, when there is one.
    pub fn wall_cells(&self) -> impl Iterator<Item = (usize, usize, f64, [f64; 3], Option<usize>)> + '_ {
        self.walls.iter().map(|w| (w.cell as usize, usize::from(w.link), w.distance, w.velocity, (w.second != NONE).then_some(w.second as usize)))
    }

    /// Whether every fluid cell's wall distance is kept, for Van Driest's damping.
    pub fn damped(&self) -> bool {
        self.a_plus.is_some()
    }

    /// Cell `cell`'s distance to its wall (lattice units), when a wall reaches it: every fluid
    /// cell's with damping, the wall cells' alone without.
    pub fn distance(&self, cell: usize) -> Option<f64> {
        let w = self.walls.get(*self.wall_of.get(cell)? as usize)?;
        if self.distance.is_empty() { (w.cell as usize == cell).then_some(w.distance) } else { self.distance.get(cell).copied() }
    }

    /// With damping: refresh every wall cell's friction velocity from `velocity(cell)`, the fluid
    /// velocity of the cell (lattice units), in parallel: each wall cell's is its own.
    pub fn refresh(&mut self, velocity: impl Fn(usize) -> [f64; 3] + Sync) {
        let nu = self.nu;
        parallel::for_each_mut(&mut self.walls, |_, w| {
            w.friction = friction_velocity(w.tangential_speed(velocity(w.cell as usize)), w.distance, nu);
        });
    }

    /// The largest y⁺ of a wall cell's centre at the last refresh (0 without wall cells).
    pub fn largest_wall_y_plus(&self) -> f64 {
        let y = self.walls.iter().map(|w| w.distance * w.friction / self.nu).fold(0.0f64, f64::max);
        if y.is_finite() { y } else { 0.0 }
    }

    /// The factor `1 − exp(−y⁺/A⁺)` of cell `cell`'s length scale: 1 where no wall reaches it,
    /// and without damping.
    #[inline(always)]
    pub fn factor(&self, cell: usize) -> f64 {
        let Some(a_plus) = self.a_plus else { return 1.0 };
        let Some(&w) = self.wall_of.get(cell).filter(|w| **w != NONE) else { return 1.0 };
        let (Some(y), Some(wall)) = (self.distance.get(cell), self.walls.get(w as usize)) else { return 1.0 };
        let x = y * wall.friction / self.nu / a_plus;
        if x.is_finite() && x < UNDAMPED { 1.0 - (-x).exp() } else { 1.0 }
    }

    /// The wall model's sample before a collision: from `velocity(cell)`, the fluid velocity of
    /// each wall cell's second cell (lattice units, with the half force), every wall cell's
    /// velocity and relaxation rate ([`wall_law`]), in parallel: each wall cell's is its own.
    pub fn sample(&mut self, velocity: impl Fn(usize) -> [f64; 3] + Sync) {
        let nu = self.nu;
        parallel::for_each_mut(&mut self.walls, |_, w| {
            if w.second == NONE {
                return;
            }
            let u2 = velocity(w.second as usize);
            let rel = [0, 1, 2].map(|a| u2[a] - w.velocity[a]);
            let un = rel[0] * w.normal[0] + rel[1] * w.normal[1] + rel[2] * w.normal[2];
            let t = [0, 1, 2].map(|a| rel[a] - un * w.normal[a]);
            let speed = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            let (u1, nu_t) = wall_law(speed, w.distance, 3.0 * w.distance, nu);
            let scale = if speed > 0.0 && u1.is_finite() { u1 / speed } else { 0.0 };
            w.boundary = [0, 1, 2].map(|a| w.velocity[a] + scale * t[a]);
            w.rate = 1.0 / (3.0 * (nu + nu_t) + 0.5);
        });
    }

    /// The wall model's velocity and relaxation rate of wall cell `cell` from the last sample
    /// (lattice units), when it is a modelled wall cell.
    pub fn modelled(&self, cell: usize) -> Option<([f64; 3], f64)> {
        let w = self.walls.get(*self.wall_of.get(cell)? as usize).filter(|w| w.cell as usize == cell && w.second != NONE)?;
        Some((w.boundary, w.rate))
    }

    /// Record the momentum the reconstruction took from each wall cell in the last collision:
    /// `removed` holds `(cell, j_before − j_after)` of the cells it rebuilt, in cell order; every
    /// other wall cell took none.
    pub fn set_removed(&mut self, removed: &[(u32, [f64; 3])]) {
        let mut next = removed.iter().peekable();
        for w in &mut self.walls {
            w.removed = match next.next_if(|(cell, _)| *cell == w.cell) {
                Some((_, j)) => *j,
                None => [0.0; 3],
            };
        }
    }

    /// The momentum each wall cell's reconstruction took from the flow in the last collision,
    /// with the force group it goes to, in cell order (none without the wall model).
    pub fn removed(&self) -> impl Iterator<Item = (usize, [f64; 3])> + '_ {
        self.walls.iter().map(|w| (w.group as usize, w.removed))
    }

    /// The wall model's reconstruction of wall cell `cell` (Malaspinas and Sagaut, eqs. 41–46 and
    /// §4.3): its populations `fq` after streaming are replaced by the equilibrium at the modelled
    /// velocity and the density of its known populations, plus the regularised non-equilibrium
    /// stress of its known populations, the unknown ones (from the wall) mirrored from their
    /// opposites. `a` is the cell's body acceleration (the stored velocity is u − a/2). Returns
    /// the rate the cell then relaxes with; `None` when `cell` is not a modelled wall cell.
    pub fn reconstruct(&self, cell: usize, fq: &mut [f64; Q], a: [f64; 3]) -> Option<f64> {
        let w = self.walls.get(*self.wall_of.get(cell)? as usize).filter(|w| w.cell as usize == cell && w.second != NONE)?;
        let link = c(usize::from(w.link));
        // ξ_q · n_out with the outward normal −c_link: known populations have it ≥ 0.
        let out = |q: usize| -dot(c(q), link);
        let u = [0, 1, 2].map(|k| w.boundary[k] - 0.5 * a[k]);
        let (mut rho0, mut rho_plus) = (0.0, 0.0);
        for (q, f) in fq.iter().enumerate() {
            let o = out(q);
            if o == 0.0 {
                rho0 += f;
            } else if o > 0.0 {
                rho_plus += f;
            }
        }
        let u_out = -dot(u, link) / dot(link, link).sqrt();
        let rho = (2.0 * rho_plus + rho0) / (1.0 + u_out);
        if !(rho > 0.0 && rho.is_finite()) {
            return None;
        }
        let eq: [f64; Q] = std::array::from_fn(|q| equilibrium(q, rho, u));
        let known = |q: usize| if out(q) >= 0.0 { q } else { opposite(q) };
        let mut p = [0.0f64; 6];
        for q in 0..Q {
            let k = known(q);
            let f1 = fq.get(k).copied().unwrap_or(0.0) - eq.get(k).copied().unwrap_or(0.0);
            let cq = c(q);
            for (slot, (x, y)) in p.iter_mut().zip(PAIRS) {
                *slot += (cq[x] * cq[y] - if x == y { CS2 } else { 0.0 }) * f1;
            }
        }
        for (q, f) in fq.iter_mut().enumerate() {
            let cq = c(q);
            let qp: f64 =
                PAIRS.iter().zip(&p).map(|(&(x, y), v)| (cq[x] * cq[y] - if x == y { CS2 } else { 0.0 }) * v * if x == y { 1.0 } else { 2.0 }).sum();
            *f = eq.get(q).copied().unwrap_or(0.0) + W.get(q).copied().unwrap_or(0.0) * 4.5 * qp;
        }
        Some(w.rate)
    }
}

/// The six components of a symmetric tensor: xx, yy, zz, xy, xz, yz.
const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The wall law of the wall model (Malaspinas and Sagaut, §4.2–4.3, with Musker's profile ũ, their
/// eq. 59): u_τ from the speed `speed` parallel to the wall at the second cell's distance `y2`,
/// then the first cell's modelled speed ũ(y1) and its mixing-length eddy viscosity
/// `ν_t = (κ y1 (1 − exp(−y1⁺/A⁺)))² |∂ũ/∂y|` (κ = 0.384, A⁺ = 26). Both 0 at rest.
pub fn wall_law(speed: f64, y1: f64, y2: f64, nu: f64) -> (f64, f64) {
    let ut = musker_friction(speed, y2, nu);
    if !(ut > 0.0) {
        return (0.0, 0.0);
    }
    let y_plus = y1 * ut / nu;
    let (u_plus, slope) = musker(y_plus);
    let dudy = ut * ut / nu * slope;
    let l = MIXING_KAPPA * y1 * (1.0 - (-y_plus / MIXING_A_PLUS).exp());
    let nu_t = l * l * dudy.abs();
    if (ut * u_plus).is_finite() && nu_t.is_finite() { (ut * u_plus, nu_t) } else { (0.0, 0.0) }
}

/// Musker's profile u⁺(y⁺) and its slope du⁺/dy⁺ (A. J. Musker, *Explicit expression for the
/// smooth wall velocity distribution in a turbulent boundary layer*, AIAA J. 17(6):655–657, 1979,
/// with the constants Malaspinas and Sagaut print in their eq. 59, 0.434 standing for 1/ln 10):
/// u⁺ = 5.424 atan((2 y⁺ − 8.15) / 16.7) + 0.434 ln((y⁺ + 10.6)^9.6 / (y⁺² − 8.15 y⁺ + 86)²) − 3.50728.
pub fn musker(y_plus: f64) -> (f64, f64) {
    let (a, b) = (2.0 / 16.7, 8.15 / 16.7);
    let t = a * y_plus - b;
    let q = y_plus * y_plus - 8.15 * y_plus + 86.0;
    let u = 5.424 * t.atan() + 0.434 * (9.6 * (y_plus + 10.6).ln() - 2.0 * q.ln()) - 3.507_279_019_362_648_4;
    let du = 5.424 * a / (1.0 + t * t) + 0.434 * (9.6 / (y_plus + 10.6) - 2.0 * (2.0 * y_plus - 8.15) / q);
    (u, du)
}

/// The friction velocity of a flow of speed `speed` parallel to a wall at distance `y` with
/// kinematic viscosity `nu`, from Musker's profile: the root of u_τ u⁺(y u_τ / ν) = U, which rises
/// with u_τ (Newton steps guarded by bisection; the bracket grows from the sublayer's √(ν U / y)
/// until it holds the root). 0 at rest or with a number that is not usable.
pub fn musker_friction(speed: f64, y: f64, nu: f64) -> f64 {
    if !(speed > 0.0 && speed.is_finite() && y > 0.0 && y.is_finite() && nu > 0.0 && nu.is_finite()) {
        return 0.0;
    }
    let g = |ut: f64| {
        let (u, du) = musker(y * ut / nu);
        (ut * u - speed, u + y * ut / nu * du)
    };
    let (mut lo, mut hi) = (0.0, (nu * speed / y).sqrt());
    for _ in 0..SPALDING_ITERATIONS {
        if g(hi).0 > 0.0 {
            break;
        }
        lo = hi;
        hi *= 2.0;
    }
    if !(g(hi).0 > 0.0) {
        return 0.0;
    }
    let mut ut = 0.5 * (lo + hi);
    for _ in 0..SPALDING_ITERATIONS {
        let (value, slope) = g(ut);
        if value == 0.0 {
            break;
        }
        if value > 0.0 {
            hi = ut;
        } else {
            lo = ut;
        }
        if hi - lo <= 1e-14 * hi {
            break;
        }
        let newton = ut - value / slope;
        ut = if newton.is_finite() && newton > lo && newton < hi { newton } else { 0.5 * (lo + hi) };
    }
    ut
}

/// Every fluid cell's distance to the nearest of the wall points of `found`, spread outward over
/// the fluid from the wall cells in a fixed order; `wall_of` takes each cell's wall cell.
fn spread(grid: &Grid, fluid: &[bool], found: &[WallCell], wall_of: &mut [u32]) -> Result<Vec<f64>> {
    let n = grid.len();
    let mut distance = filled(n, f64::INFINITY, "wall distances")?;
    for w in found {
        if let Some(s) = distance.get_mut(w.cell as usize) {
            *s = w.distance;
        }
    }
    let dims = grid.cells.map(|v| v as i64);
    let mut queue: std::collections::VecDeque<u32> = found.iter().map(|w| w.cell).collect();
    let cap = n.saturating_mul(MAX_UPDATES_PER_CELL);
    let mut updates = 0usize;
    while let Some(cell) = queue.pop_front() {
        let Some(&w) = wall_of.get(cell as usize) else { continue };
        let Some(point) = found.get(w as usize).map(|w| w.point) else { continue };
        let p = grid.coords(cell as usize).map(|v| v as i64);
        for o in NEIGHBOURS {
            let s = [p[0] + o[0], p[1] + o[1], p[2] + o[2]];
            if (0..3).any(|a| s[a] < 0 || s[a] >= dims[a]) {
                continue;
            }
            let next = grid.index(s.map(|v| v as usize));
            if !fluid.get(next).copied().unwrap_or(false) {
                continue;
            }
            let d = (0..3).map(|a| (s[a] as f64 - point[a]).powi(2)).sum::<f64>().sqrt();
            let (Some(slot), Some(owner)) = (distance.get_mut(next), wall_of.get_mut(next)) else { continue };
            if d < *slot {
                *slot = d;
                *owner = w;
                queue.push_back(next as u32);
                updates += 1;
                if updates > cap {
                    return Err(crate::Error::Limit(format!("the wall distances did not settle after {cap} updates over {n} cells")));
                }
            }
        }
    }
    Ok(distance)
}

/// The run's warning when the first cells off the walls reached y⁺ `y_plus`, past the
/// sublayer the damping is meant for; none within it.
pub fn unresolved_wall(y_plus: f64) -> Option<String> {
    (y_plus > RESOLVED_WALL_Y_PLUS).then(|| {
        format!(
            "the first cells off the walls reached y⁺ {y_plus:.1}, past the viscous sublayer (y⁺ {RESOLVED_WALL_Y_PLUS}) Van Driest damping is meant for: the grid does not resolve the wall; refine the spacing near the walls, or solve without `van_driest`"
        )
    })
}

/// Spalding's `y⁺(u⁺)` and its slope `dy⁺/du⁺`, with `e = exp(−κB)`.
fn spalding(u_plus: f64, e: f64) -> (f64, f64) {
    let k = KAPPA * u_plus;
    let rest = k.exp() - 1.0 - k - k * k / 2.0;
    (u_plus + e * (rest - k * k * k / 6.0), 1.0 + KAPPA * e * rest)
}

/// The friction velocity of a flow of speed `speed` parallel to a wall at distance `y` with
/// kinematic viscosity `nu`, from Spalding's law: 0 at rest or with a number that is not usable.
pub fn friction_velocity(speed: f64, y: f64, nu: f64) -> f64 {
    let re = speed * y / nu;
    if !(re.is_finite() && re > 0.0) {
        return 0.0;
    }
    // u⁺ y⁺(u⁺) = Re_y rises with u⁺; the sublayer's root √Re_y bounds it from above.
    let e = (-KAPPA * B).exp();
    let (mut lo, mut hi) = (0.0, re.sqrt());
    // The log law at y⁺ = Re_y lies above the root and near it once Re_y is large; Newton walks
    // down an exponential only a step of 1/κ at a time, so it starts there.
    let mut up = (re.ln() / KAPPA + B).clamp(f64::MIN_POSITIVE, hi);
    let mut last = f64::INFINITY;
    for _ in 0..SPALDING_ITERATIONS {
        let (y_plus, slope) = spalding(up, e);
        let value = up * y_plus - re;
        if value == 0.0 {
            break;
        }
        if value > 0.0 {
            hi = up;
        } else {
            lo = up;
        }
        if hi - lo <= 1e-14 * hi {
            break;
        }
        // d/du⁺ (u⁺ y⁺) = y⁺ + u⁺ dy⁺/du⁺. Newton while it lands inside the bracket and at
        // least halves the residual, bisection otherwise.
        let newton = up - value / (y_plus + up * slope);
        let fast = newton.is_finite() && newton > lo && newton < hi && value.abs() <= 0.5 * last;
        last = value.abs();
        up = if fast { newton } else { 0.5 * (lo + hi) };
    }
    if up > 0.0 { speed / up } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::C;

    /// In the viscous sublayer u⁺ = y⁺, so u_τ² = ν U / y; in the log layer u⁺ = ln(y⁺)/κ + B.
    #[test]
    fn spalding_gives_the_sublayer_and_the_log_law() {
        let nu = 1e-3;
        for u_plus in [0.5, 1.0, 2.0] {
            let ut = 0.05;
            let y_plus = spalding(u_plus, (-KAPPA * B).exp()).0;
            assert!((y_plus / u_plus - 1.0).abs() < 2e-3, "the sublayer: {y_plus} at u⁺ {u_plus}");
            let got = friction_velocity(ut * u_plus, y_plus * nu / ut, nu);
            assert!((got / ut - 1.0).abs() < 1e-10, "u⁺ {u_plus}: {got}");
        }
        for y_plus in [100.0f64, 500.0, 2000.0] {
            let (ut, y) = (0.02, y_plus * 1e-3 / 0.02);
            let u = ut * (y_plus.ln() / KAPPA + B);
            let got = friction_velocity(u, y, 1e-3);
            // Spalding blends into the log law: 1 % off it at y⁺ 100, less above.
            assert!((got / ut - 1.0).abs() < 0.02, "y⁺ {y_plus}: {got} vs {ut}");
        }
        // The formula itself, round trip at any y⁺.
        for y_plus in [5.0, 11.0, 30.0, 1e4] {
            let up = {
                let (mut lo, mut hi) = (0.0, y_plus);
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if spalding(mid, (-KAPPA * B).exp()).0 > y_plus { hi = mid } else { lo = mid }
                }
                0.5 * (lo + hi)
            };
            let got = friction_velocity(0.03 * up, y_plus * 1e-4 / 0.03, 1e-4);
            assert!((got / 0.03 - 1.0).abs() < 1e-9, "y⁺ {y_plus}: {got}");
        }
    }

    /// Musker's profile is the sublayer's u⁺ = y⁺ near the wall and the log law above; the wall
    /// law's u_τ makes the second cell sit on it, and sets the first cell on it too.
    #[test]
    fn the_wall_law_follows_musker_s_profile() {
        let (u, du) = musker(0.5);
        assert!((u - 0.5).abs() < 5e-3 && (du - 1.0).abs() < 2e-2, "sublayer: {u} {du}");
        for y_plus in [300.0f64, 1000.0] {
            let log_law = y_plus.ln() / 0.41 + 5.0;
            assert!((musker(y_plus).0 / log_law - 1.0).abs() < 0.02, "log law at {y_plus}: {}", musker(y_plus).0);
        }
        // The slope against a central difference.
        for y_plus in [1.0, 12.0, 80.0] {
            let h = 1e-5;
            let fd = (musker(y_plus + h).0 - musker(y_plus - h).0) / (2.0 * h);
            assert!((musker(y_plus).1 / fd - 1.0).abs() < 1e-8, "{y_plus}");
        }
        let (nu, ut, y1, y2) = (1e-4, 0.05, 0.5, 1.5);
        let u2 = ut * musker(y2 * ut / nu).0;
        assert!((musker_friction(u2, y2, nu) / ut - 1.0).abs() < 1e-12);
        let (u1, nu_t) = wall_law(u2, y1, y2, nu);
        assert!((u1 / (ut * musker(y1 * ut / nu).0) - 1.0).abs() < 1e-12);
        assert!(nu_t > nu, "the mixing length shows at y⁺ {}", y1 * ut / nu);
        assert_eq!(wall_law(0.0, 0.5, 1.5, nu), (0.0, 0.0));
        assert_eq!(wall_law(f64::NAN, 0.5, 1.5, nu), (0.0, 0.0));
        assert_eq!(musker_friction(1.0, 0.0, nu), 0.0);
    }

    /// A wall cell whose populations are already the equilibrium at the modelled velocity plus
    /// a regularised stress comes out of the reconstruction as it went in: the density of its
    /// known populations, the velocity and the stress are kept, the mirrored unknowns agree.
    #[test]
    fn the_reconstruction_keeps_a_consistent_wall_cell() {
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [3, 8, 2] };
        let fluid = vec![true; grid.len()];
        let mut links = Vec::new();
        for cell in 0..grid.len() {
            let [_, j, _] = grid.coords(cell);
            for (q, cq) in C.iter().enumerate().skip(1) {
                let from = j as i64 - i64::from(cq[1]);
                if !(0..8).contains(&from) {
                    links.push(Link { cell: cell as u32, q: q as u8, rule: Rule::Bounce { group: if from < 0 { 2 } else { 3 }, du: 0.0 } });
                }
            }
        }
        let walls = [None, None, Some([0.0; 3]), Some([0.0; 3]), None, None];
        let mut w = Walls::new(&grid, &fluid, &links, &walls, None, 1e-3).unwrap();
        let cell = grid.index([1, 0, 0]);
        let k = w.wall_of[cell] as usize;
        assert_eq!(w.walls[k].second as usize, grid.index([1, 1, 0]));
        let (rho, u, p) = (1.02, [0.03, 0.0, -0.01], [1e-4, -2e-4, 1e-4, 3e-4, 0.0, -1e-4]);
        w.walls[k].boundary = u;
        let mut fq: [f64; Q] = std::array::from_fn(|q| {
            let cq = c(q);
            let qp: f64 =
                PAIRS.iter().zip(&p).map(|(&(x, y), v)| (cq[x] * cq[y] - if x == y { CS2 } else { 0.0 }) * v * if x == y { 1.0 } else { 2.0 }).sum();
            equilibrium(q, rho, u) + W[q] * 4.5 * qp
        });
        let before = fq;
        let rate = w.reconstruct(cell, &mut fq, [0.0; 3]).unwrap();
        assert_eq!(rate, w.walls[k].rate);
        for q in 0..Q {
            assert!((fq[q] - before[q]).abs() < 1e-15, "{q}: {} vs {}", fq[q], before[q]);
        }
        // Not a wall cell: nothing to do.
        assert!(w.reconstruct(grid.index([1, 3, 0]), &mut fq, [0.0; 3]).is_none());
    }

    #[test]
    fn a_still_or_hostile_flow_has_no_friction_velocity() {
        assert_eq!(friction_velocity(0.0, 0.5, 1e-3), 0.0);
        assert_eq!(friction_velocity(f64::NAN, 0.5, 1e-3), 0.0);
        assert_eq!(friction_velocity(1.0, 0.5, 0.0), 0.0);
        assert_eq!(friction_velocity(-1.0, 0.5, 1e-3), 0.0);
        assert!(friction_velocity(1e300, 1e300, 1e-300) == 0.0);
        assert!(friction_velocity(1e10, 0.5, 1e-10).is_finite());
    }

    /// A box with walls at y = −½ and y = ny − ½ (bounce-back halfway): every cell's distance is
    /// to the nearer of the two planes.
    #[test]
    fn plane_walls_give_the_exact_distance() {
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [3, 8, 2] };
        let fluid = vec![true; grid.len()];
        let mut links = Vec::new();
        for cell in 0..grid.len() {
            let [_, j, _] = grid.coords(cell);
            for (q, cq) in C.iter().enumerate().skip(1) {
                let from = j as i64 - i64::from(cq[1]);
                if !(0..8).contains(&from) {
                    let group = if from < 0 { 2 } else { 3 };
                    links.push(Link { cell: cell as u32, q: q as u8, rule: Rule::Bounce { group, du: 0.0 } });
                }
            }
        }
        let walls = [None, None, Some([0.0; 3]), Some([0.0; 3]), None, None];
        let w = Walls::new(&grid, &fluid, &links, &walls, Some(26.0), 1e-3).unwrap();
        assert_eq!(w.walls.len(), 2 * 3 * 2);
        for cell in 0..grid.len() {
            let [_, j, _] = grid.coords(cell);
            let want = (j as f64 + 0.5).min(8.0 - 0.5 - j as f64);
            assert!((w.distance(cell).unwrap() - want).abs() < 1e-12, "{cell}");
        }
        // At rest, full damping; away from the walls no damping once the flow moves.
        assert_eq!(w.factor(0), 0.0);
        let mut w = w;
        w.refresh(|_| [0.1, 0.0, 0.0]);
        // u⁺ ≈ 7 at the wall cell (Re_y 50), so y⁺ ≈ 50 three and a half cells off it.
        assert!(w.factor(grid.index([1, 4, 0])) > 0.8);
        assert!(w.factor(grid.index([1, 0, 0])) < w.factor(grid.index([1, 1, 0])));
        // Slip walls give no wall cells.
        let none = Walls::new(&grid, &fluid, &links, &[None; 6], Some(26.0), 1e-3).unwrap();
        assert_eq!((none.walls.len(), none.factor(0)), (0, 1.0));
        // Far from the walls the factor is exactly 1 (the early return, no rounding change).
        assert_eq!(1.0 - (-UNDAMPED).exp(), 1.0);
    }
}
