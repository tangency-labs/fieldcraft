//! Convective pressure outlets: each outlet cell holds the face's pressure, and its velocity is
//! carried out to it by the flow instead of following from a uniform-pressure balance.
//!
//! A Zou–He pressure outlet holds its cell at the face's density and sets the normal velocity
//! from mass and normal momentum, with no tangential velocity. Where a wake that has not
//! recovered reaches it, that balance puts a step into the velocity of the last layer (#590,
//! #599). A convective outlet lets the wake leave.
//!
//! # The condition
//! The populations entering the box from outside an outlet cell N (the unknown ones,
//! `c_q · n < 0` with `n` the outward normal) obey the convective outflow equation
//! `∂f/∂t + U ∂f/∂n = 0`, discretised upwind and implicitly in time (I. Orlanski, *A simple
//! boundary condition for unbounded hyperbolic flows*, J. Comput. Phys. 21(3):251–269, 1976,
//! with the local flow speed as the phase speed instead of his estimate from the interior):
//!
//! `f_q(N, t+1) = (f_q(N, t) + λ f_q(N−1, t+1)) / (1 + λ)`, `λ = U Δt / Δx`,
//!
//! with `f_q(N−1, t+1)` the population of the cell one layer in after streaming and `f_q(N, t)`
//! the outlet cell's own after the last step. With pull streaming, an incoming direction of the
//! cell inside holds what cell N itself sent inward from its last collision, so the update blends
//! the outlet cell's incoming populations before and after its own collision: with λ = 0 they
//! stay as they were, with λ → ∞ they become what its collision sent inward. The update is a
//! convex blend, stable for every `λ ≥ 0`. `U` is the cell's own outward velocity after the last
//! step, clamped to `[0, MAX_LATTICE_SPEED]`: the speed the flow there is carried out at. (On the
//! wake of Schäfer–Turek's cylinder cut 1.5 diameters behind it, the local speed came closer to a
//! long domain's solution at the outlet than the face's mean speed, λ = 0 or λ → ∞.) The known
//! populations stream in as everywhere else.
//!
//! **Pressure.** The cell's density is then set to the face's: every population gains
//! `w_q (ρ_out − ρ)`, which changes the density only, not the momentum. The face holds its
//! pressure cell by cell, as the uniform outlet does, while the velocity, normal and tangential,
//! comes from the flow on both sides of the face.
//!
//! **Backflow guard.** A cell whose convected flow points into the box takes fluid in, and the
//! convective equation carries nothing to it from inside. It is the Zou–He pressure outlet
//! (`flow`): the face's density, no tangential velocity, which stays stable with reversed outlet
//! cells (#598).
//!
//! The rule changes only the streamed populations of the outlet cells, so it does not depend on
//! the collision operator.
//!
//! The convective condition for lattice Boltzmann outlets: Q. Lou, Z. Guo, B. Shi, *Evaluation
//! of outflow boundary conditions for two-phase lattice Boltzmann equation*, PRE 87:063301,
//! 2013; Z. Yang, *Lattice Boltzmann outflow treatments: convective conditions and others*,
//! Comput. Math. Appl. 65(2):160–171, 2013 (both DOIs not checked; no code read).
//!
//! # Determinism
//! Every cell's update reads only populations streamed before any open cell changes, and writes
//! only its own: the order of the cells does not matter.

use crate::flow::{OpenCell, Target};
use crate::lattice::{C, Q};
use crate::setup::MAX_LATTICE_SPEED;

/// The unknown populations of an open cell of D3Q19: one along the normal, four diagonals.
pub const UNKNOWN: usize = 5;

/// The directions entering an open cell from outside (`c_q` along the inward normal of `axis`,
/// `inward` +1 at a min face), in increasing order.
pub fn unknowns(axis: usize, inward: f64) -> [usize; UNKNOWN] {
    let mut out = [0; UNKNOWN];
    let mut k = 0;
    for q in 1..Q {
        if f64::from(C[q][axis.min(2)]) * inward > 0.0
            && let Some(slot) = out.get_mut(k)
        {
            *slot = q;
            k += 1;
        }
    }
    out
}

/// `(f_N + λ f_{N−1}) / (1 + λ)`.
#[inline]
fn convect(own: f64, inside: f64, speed: f64) -> f64 {
    (own + speed * inside) / (1.0 + speed)
}

/// The convection speed `λ` from an outward velocity (lattice units).
fn speed(outward: f64) -> f64 {
    if outward > 0.0 { outward.min(MAX_LATTICE_SPEED) } else { 0.0 }
}

/// The closure of a convective outlet cell where the flow comes back in: the uniform pressure
/// outlet's, so it does whatever a `Target::Density` face does. The one place that decides it
/// (the GPU's `BACKFLOW_KIND` in lbm.wgsl is the same choice).
pub fn backflow(density: f64) -> Target {
    Target::Density(density)
}

/// What a convective outlet cell keeps from one step to the next.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CellState {
    /// Its unknown populations ([`unknowns`]) after the last step.
    pub unknowns: [f64; UNKNOWN],
    /// The convection speed `λ` of the next step: its outward velocity after the last step.
    pub speed: f64,
}

/// What a convective outlet cell does this step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Action {
    /// Its unknown populations are the convected ones, and every population gains `w_q shift`.
    Convect { unknowns: [f64; UNKNOWN], shift: f64 },
    /// The flow comes in: Zou–He at the face's density, no tangential velocity.
    Hold,
}

/// The populations of a flow the outlets read: population `q` of `cell`, whether a cell is
/// fluid, and half the body force per mass of a cell (for the velocity).
pub trait Populations {
    fn get(&self, q: usize, cell: usize) -> f64;
    fn is_fluid(&self, cell: usize) -> bool;
    fn half_force(&self, cell: usize) -> [f64; 3];
}

/// The density and outward velocity (with the half force) of open cell `oc` whose populations
/// are `f(q)`.
fn density_outward(oc: &OpenCell, half_force: [f64; 3], f: impl Fn(usize) -> f64) -> (f64, f64) {
    let axis = oc.axis.min(2);
    let (mut rho, mut jn) = (0.0, 0.0);
    for q in 0..Q {
        let v = f(q);
        rho += v;
        jn += v * f64::from(C[q][axis]);
    }
    let a = half_force.get(axis).copied().unwrap_or(0.0);
    let un = if rho > 0.0 { -oc.inward * (jn + rho * a) / rho } else { 0.0 };
    (rho, un)
}

/// The convective outlet cells of a flow: what each keeps between steps.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Outlets {
    /// Per open cell (convective ones only; default for the others).
    state: Vec<CellState>,
    /// Per open cell: this step's action.
    actions: Vec<Action>,
}

impl Outlets {
    pub(crate) fn checkpoint(&mut self, a: &mut dyn crate::checkpoint::Archive) -> crate::Result<()> {
        use crate::checkpoint::{f64, fixed};
        fixed(a, self.state.len() as u64)?;
        for s in &mut self.state {
            for v in &mut s.unknowns {
                f64(a, v)?;
            }
            f64(a, &mut s.speed)?;
            if !(0.0..=MAX_LATTICE_SPEED).contains(&s.speed) {
                return Err(crate::Error::Invalid("checkpoint outlet speed is outside the existing scheme's bound".into()));
            }
        }
        // Actions are recomputed by convect before the next boundary; only CellState persists.
        Ok(())
    }

    /// The bytes its arrays hold.
    pub fn heap_bytes(&self) -> u64 {
        (self.state.capacity() * size_of::<CellState>() + self.actions.capacity() * size_of::<Action>()) as u64
    }

    /// The state of the convective cells among `open` from the populations `pop`.
    /// Empty without convective cells.
    pub fn new(open: &[OpenCell], pop: &impl Populations) -> Outlets {
        if !open.iter().any(|oc| matches!(oc.target, Target::Convective(_))) {
            return Outlets::default();
        }
        let state = open.iter().map(|oc| state_of(oc, pop.half_force(oc.cell as usize), |q| pop.get(q, oc.cell as usize))).collect();
        Outlets { state, actions: vec![Action::Hold; open.len()] }
    }

    pub fn is_empty(&self) -> bool {
        self.state.is_empty()
    }

    /// Per open cell, its state (meaningful for convective cells only).
    pub fn states(&self) -> &[CellState] {
        &self.state
    }

    /// This step's action of every convective cell from the streamed populations `pop`
    /// (`strides`: the index step along each axis). Reads only, so every cell sees the
    /// populations before any open boundary changes them.
    pub fn convect(&mut self, open: &[OpenCell], pop: &impl Populations, strides: [usize; 3]) {
        for ((oc, own), action) in open.iter().zip(&self.state).zip(&mut self.actions) {
            let Target::Convective(density) = oc.target else { continue };
            let cell = oc.cell as usize;
            let dirs = unknowns(oc.axis, oc.inward);
            let step = strides.get(oc.axis).copied().unwrap_or(0);
            let inside = if oc.inward > 0.0 { cell.checked_add(step) } else { cell.checked_sub(step) }.filter(|c| pop.is_fluid(*c));
            let mut u = own.unknowns;
            for (x, q) in u.iter_mut().zip(dirs) {
                // Without a fluid cell inside, the cell keeps its own.
                if let Some(c) = inside {
                    *x = convect(*x, pop.get(q, c), own.speed);
                }
            }
            let mut f = [0.0; Q];
            for (q, x) in f.iter_mut().enumerate() {
                *x = pop.get(q, cell);
            }
            for (q, x) in dirs.iter().zip(u) {
                if let Some(slot) = f.get_mut(*q) {
                    *slot = x;
                }
            }
            let (rho, un) = density_outward(oc, pop.half_force(cell), |q| f.get(q).copied().unwrap_or(0.0));
            *action = if un > 0.0 { Action::Convect { unknowns: u, shift: density - rho } } else { Action::Hold };
        }
    }

    /// What open cell `k` does this step.
    pub fn action(&self, k: usize) -> Action {
        self.actions.get(k).copied().unwrap_or(Action::Hold)
    }

    /// Record open cell `k`'s state from its populations `f(q)` after the boundary.
    pub fn record(&mut self, k: usize, oc: &OpenCell, half_force: [f64; 3], f: impl Fn(usize) -> f64) {
        if let Some(s) = self.state.get_mut(k) {
            *s = state_of(oc, half_force, f);
        }
    }
}

/// The state of open cell `oc` with populations `f(q)`.
fn state_of(oc: &OpenCell, half_force: [f64; 3], f: impl Fn(usize) -> f64) -> CellState {
    let (_, un) = density_outward(oc, half_force, &f);
    CellState { unknowns: unknowns(oc.axis, oc.inward).map(&f), speed: speed(un) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::{W, equilibrium};

    #[test]
    fn the_unknowns_enter_from_outside_and_the_blend_is_convex() {
        for axis in 0..3 {
            for inward in [1.0, -1.0] {
                let dirs = unknowns(axis, inward);
                assert!(dirs.iter().all(|q| f64::from(C[*q][axis]) * inward > 0.0), "{axis} {inward}: {dirs:?}");
                assert!(dirs.windows(2).all(|w| w[0] < w[1]));
            }
        }
        assert_eq!(convect(1.0, 3.0, 0.0), 1.0, "at rest nothing is carried");
        assert_eq!(convect(1.0, 3.0, 1.0), 2.0);
        assert_eq!(speed(-0.1), 0.0);
        assert_eq!(speed(f64::NAN), 0.0);
        assert_eq!(speed(9.0), MAX_LATTICE_SPEED);
    }

    /// Cells at equilibrium: 0..4 the layer inside an x_max face, 4..8 the outlet layer.
    struct Cells([(f64, [f64; 3]); 8]);

    impl Populations for Cells {
        fn get(&self, q: usize, cell: usize) -> f64 {
            self.0.get(cell).map_or(0.0, |(rho, u)| equilibrium(q, *rho, *u))
        }
        fn is_fluid(&self, cell: usize) -> bool {
            cell < 8
        }
        fn half_force(&self, _: usize) -> [f64; 3] {
            [0.0; 3]
        }
    }

    #[test]
    fn leaving_cells_hold_the_density_and_keep_their_momentum_and_reversed_cells_hold() {
        let mut cells = [(1.0, [0.05, 0.0, 0.0]); 8];
        cells[1] = (1.0, [0.05, 0.01, 0.0]);
        cells[5] = (0.99, [0.05, 0.0, 0.0]);
        cells[6] = (1.0, [-0.02, 0.0, 0.0]);
        cells[7] = (1.03, [0.05, 0.0, 0.0]);
        let pop = Cells(cells);
        let open: Vec<OpenCell> = (4..8).map(|c| OpenCell { cell: c, axis: 0, inward: -1.0, target: Target::Convective(1.0) }).collect();
        let mut o = Outlets::new(&open, &pop);
        assert!((o.states()[0].speed - 0.05).abs() < 1e-12 && o.states()[2].speed == 0.0, "the local outward speeds");
        o.convect(&open, &pop, [4, 1, 1]);
        assert_eq!(o.action(2), Action::Hold, "the reversed cell holds the density");
        let dirs = unknowns(0, -1.0);
        let j = |g: &dyn Fn(usize) -> f64, a: usize| (0..Q).map(|q| g(q) * f64::from(C[q][a])).sum::<f64>();
        for k in [0, 1, 3] {
            let Action::Convect { unknowns: u, shift } = o.action(k) else { panic!("cell {k} convects") };
            let cell = open[k].cell as usize;
            let f = |q: usize| dirs.iter().position(|d| *d == q).map_or(pop.get(q, cell), |i| u[i]);
            let shifted = |q: usize| f(q) + W[q] * shift;
            let (rho, _) = density_outward(&open[k], [0.0; 3], shifted);
            assert!((rho - 1.0).abs() < 1e-14, "cell {k} holds the face's density: {rho}");
            for a in 0..3 {
                assert!((j(&f, a) - j(&shifted, a)).abs() < 1e-15, "the shift keeps the momentum");
            }
            if k == 1 {
                assert!(j(&shifted, 1) > 0.0, "the tangential velocity inside reaches the outlet cell");
            }
        }
    }
}
