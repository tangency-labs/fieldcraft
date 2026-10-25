//! Van Driest's damping of Smagorinsky's length scale near walls.
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
//! No code of any lattice Boltzmann or LES package was read; this is written from the papers.

use fieldcraft_parallel as parallel;

use crate::Result;
use crate::flow::{Link, Rule, filled};
use crate::lattice::c;
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

/// A fluid cell next to a no-slip wall.
#[derive(Debug, Clone, Copy, PartialEq)]
struct WallCell {
    cell: u32,
    /// The wall point (lattice coordinates: cell `[i, j, k]` is centred on `[i, j, k]`).
    point: [f64; 3],
    /// The unit normal from the wall point into the fluid.
    normal: [f64; 3],
    /// The wall's velocity (lattice units).
    velocity: [f64; 3],
    /// The friction velocity at the last refresh (lattice units).
    friction: f64,
}

/// Each fluid cell's distance to the nearest no-slip wall and the wall cell that gives its
/// friction velocity.
#[derive(Debug, Clone, PartialEq)]
pub struct WallDamping {
    /// A⁺, and the fluid's kinematic viscosity (lattice units).
    a_plus: f64,
    nu: f64,
    walls: Vec<WallCell>,
    /// Per cell: its index in `walls` (`NONE`: no wall reached, or not fluid).
    wall_of: Vec<u32>,
    /// Per cell: the distance from its centre to its wall point (lattice units).
    distance: Vec<f64>,
}

impl WallDamping {
    /// The wall distances of `grid`'s fluid cells, for a fluid of kinematic viscosity `nu`
    /// (lattice units). `links` are the flow's links; `walls[g]` is the velocity of force group
    /// `g`'s wall when it is no-slip (`None`: slip, which has no friction to damp by).
    pub fn new(grid: &Grid, fluid: &[bool], links: &[Link], walls: &[Option<[f64; 3]>], a_plus: f64, nu: f64) -> Result<WallDamping> {
        let n = grid.len();
        let mut found: Vec<WallCell> = Vec::new();
        let mut wall_of = filled(n, NONE, "wall cells")?;
        let mut distance = filled(n, f64::INFINITY, "wall distances")?;
        // The links are in cell order, so each wall cell's links are consecutive.
        for links in links.chunk_by(|a, b| a.cell == b.cell) {
            let Some(cell) = links.first().map(|l| l.cell) else { continue };
            // The nearest wall point: the population of direction q arrives from −c_q, so the
            // wall is halfway there; an axis link first, then the lowest direction.
            let mut best: Option<(f64, usize, [f64; 3])> = None;
            for l in links {
                let Rule::Bounce { group, .. } = l.rule else { continue };
                let Some(Some(velocity)) = walls.get(group as usize) else { continue };
                let q = usize::from(l.q);
                let cq = c(q);
                let len = 0.5 * (cq[0] * cq[0] + cq[1] * cq[1] + cq[2] * cq[2]).sqrt();
                if best.is_none_or(|(d, b, _)| len < d || (len == d && q < b)) {
                    best = Some((len, q, *velocity));
                }
            }
            let Some((d, q, velocity)) = best else { continue };
            let centre = grid.coords(cell as usize).map(|v| v as f64);
            let normal = c(q).map(|v| v * 0.5 / d);
            let point = [0, 1, 2].map(|a| centre[a] - normal[a] * d);
            if let (Some(w), Some(s)) = (wall_of.get_mut(cell as usize), distance.get_mut(cell as usize)) {
                *w = found.len() as u32;
                *s = d;
            }
            found.push(WallCell { cell, point, normal, velocity, friction: 0.0 });
        }

        // Spread the nearest wall point outward over the fluid, in a fixed order.
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
        Ok(WallDamping { a_plus, nu, walls: found, wall_of, distance })
    }

    /// Cell `cell`'s distance to its wall (lattice units), when a wall reaches it.
    pub fn distance(&self, cell: usize) -> Option<f64> {
        self.wall_of.get(cell).filter(|w| **w != NONE).and_then(|_| self.distance.get(cell).copied())
    }

    /// Refresh every wall cell's friction velocity from `velocity(cell)`, the fluid velocity of
    /// the cell (lattice units), in parallel: each wall cell's is its own.
    pub fn refresh(&mut self, velocity: impl Fn(usize) -> [f64; 3] + Sync) {
        let (distance, nu) = (&self.distance, self.nu);
        parallel::for_each_mut(&mut self.walls, |_, w| {
            let u = velocity(w.cell as usize);
            let rel = [0, 1, 2].map(|a| u[a] - w.velocity[a]);
            let un = rel[0] * w.normal[0] + rel[1] * w.normal[1] + rel[2] * w.normal[2];
            let t = [0, 1, 2].map(|a| rel[a] - un * w.normal[a]);
            let y = distance.get(w.cell as usize).copied().unwrap_or(0.0);
            w.friction = friction_velocity((t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt(), y, nu);
        });
    }

    /// The largest y⁺ of a wall cell's centre at the last refresh (0 without wall cells).
    pub fn largest_wall_y_plus(&self) -> f64 {
        let y_plus = |w: &WallCell| self.distance.get(w.cell as usize).map_or(0.0, |y| y * w.friction / self.nu);
        let y = self.walls.iter().map(y_plus).fold(0.0f64, f64::max);
        if y.is_finite() { y } else { 0.0 }
    }

    /// The factor `1 − exp(−y⁺/A⁺)` of cell `cell`'s length scale: 1 where no wall reaches it.
    #[inline(always)]
    pub fn factor(&self, cell: usize) -> f64 {
        let Some(&w) = self.wall_of.get(cell).filter(|w| **w != NONE) else { return 1.0 };
        let (Some(y), Some(wall)) = (self.distance.get(cell), self.walls.get(w as usize)) else { return 1.0 };
        let x = y * wall.friction / self.nu / self.a_plus;
        if x.is_finite() && x < UNDAMPED { 1.0 - (-x).exp() } else { 1.0 }
    }
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
        let w = WallDamping::new(&grid, &fluid, &links, &walls, 26.0, 1e-3).unwrap();
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
        let none = WallDamping::new(&grid, &fluid, &links, &[None; 6], 26.0, 1e-3).unwrap();
        assert_eq!((none.walls.len(), none.factor(0)), (0, 1.0));
        // Far from the walls the factor is exactly 1 (the early return, no rounding change).
        assert_eq!(1.0 - (-UNDAMPED).exp(), 1.0);
    }
}
