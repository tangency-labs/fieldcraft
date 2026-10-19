//! Interval assignment: the number of mesh divisions of every curve, chosen together so that
//! every face mesher downstream can succeed with conforming meshes on shared curves.
//!
//! The problem is the integer program of Tam and Armstrong (*Finite element mesh control by
//! integer programming*, IJNME 36(15):2581–2605, 1993) and Mitchell (*High fidelity interval
//! assignment*, Proc. 6th International Meshing Roundtable, 1997; IJCGA 10(4):399–415, 2000):
//! integer counts `x_i ≥ 1` close to goals `g_i` (the curve length over the target size), subject to
//! - **equalities** `Σ_A x = Σ_B x`: opposite sides of a mapped face or of a sweep's linking
//!   face;
//! - **parity**: `Σ_L x` even on every loop of a face meshed by the free all-quad mesher (a
//!   multiple of 4 around a disk meshed as an O-grid);
//! - **minimum sums** `Σ_L x ≥ m`: a loop needs enough points to be meshed;
//! - **dominance** `Σ_A x ≥ Σ_B x + d`: a side of a triangle primitive shorter than the other two
//!   together.
//!
//! Mitchell's *incremental interval assignment* (Mitchell, *Incremental interval assignment by
//! integer linear algebra*, International Meshing Roundtable 2021; Computer-Aided Design 158:103485, 2023) starts from a
//! rounded solution and repairs violated constraints by unit steps chosen by a priority on the
//! change in relative error. This module uses that scheme with increments only: start from the
//! rounded goals, then repeatedly take the first violated constraint and raise the variable whose
//! increase costs least in squared relative error `((x − g)/g)²`, preferring, for parity, variables that
//! no equality couples to other faces (until a loop keeps coming back: free curves shared only by
//! loops round a point cannot fix their parity together). Increments only make the counts never coarser than asked
//! and make termination provable for feasible systems up to the step cap; an infeasible system
//! (one that only grows) hits the cap and is reported, never looped on.

use crate::{Error, Result};

/// The smallest division goal the error divides by.
const MIN_GOAL: f64 = 1e-3; // unit-factors: a count of divisions

/// The constraint system.
#[derive(Debug, Clone, Default)]
pub(crate) struct Intervals {
    /// Goal count of each variable (curve length over target size).
    pub goals: Vec<f64>,
    /// Lower bound of each variable (at least 1).
    pub mins: Vec<usize>,
    /// `Σ lhs = Σ rhs` (variables may repeat).
    pub equal: Vec<(Vec<usize>, Vec<usize>)>,
    /// `Σ vars` a multiple of the modulus (2: even loops for the free mesher; 4: a disk's
    /// O-grid).
    pub multiple: Vec<(Vec<usize>, usize)>,
    /// `Σ vars ≥ min`.
    pub at_least: Vec<(Vec<usize>, usize)>,
    /// `Σ lhs ≥ Σ rhs + d`.
    pub exceed: Vec<(Vec<usize>, Vec<usize>, usize)>,
    /// Variables held at their goal (an exact edge seed), with the seed in words; empty or
    /// `None`: free.
    pub fixed: Vec<Option<String>>,
}

/// Repairs of one parity constraint by curves no equality couples, before coupled ones may
/// move as well.
const PARITY_PATIENCE: usize = 4;

/// The largest count a variable may reach.
const MAX_COUNT: usize = 1 << 20;

impl Intervals {
    /// Add a variable with its goal and minimum; returns its index.
    pub(crate) fn var(&mut self, goal: f64, min: usize) -> usize {
        self.goals.push(goal);
        self.mins.push(min.max(1));
        self.goals.len() - 1
    }

    /// Apply an edge seed to variable `v`: an exact seed holds it at the seeded count, a minimum
    /// raises its lower bound (and its goal), a free seed sets its goal.
    pub(crate) fn seed(&mut self, v: usize, seed: &crate::seeds::EdgeSeed) {
        use fieldcraft_model::SeedConstraint;
        let n = seed.count.max(1);
        let (Some(goal), Some(min)) = (self.goals.get_mut(v), self.mins.get_mut(v)) else { return };
        match seed.constraint {
            SeedConstraint::Exact => {
                *goal = n as f64;
                *min = n;
                if self.fixed.len() <= v {
                    self.fixed.resize(v + 1, None);
                }
                if let Some(slot) = self.fixed.get_mut(v) {
                    *slot = Some(seed.describe());
                }
            }
            SeedConstraint::Min => {
                *min = (*min).max(n);
                *goal = goal.max(n as f64);
            }
            SeedConstraint::Free => *goal = n as f64,
        }
    }

    fn is_fixed(&self, v: usize) -> bool {
        self.fixed.get(v).is_some_and(Option::is_some)
    }

    /// The conflict error: the exact seeds among `vars` cannot all keep their counts.
    fn conflict(&self, vars: &[usize], what: &str) -> Error {
        let mut seeds: Vec<&str> = vars.iter().filter_map(|&v| self.fixed.get(v).and_then(Option::as_deref)).collect();
        seeds.sort_unstable();
        seeds.dedup();
        Error::Invalid(format!("the exact edge seeds {} conflict: {what}; change a count or make a seed `min` or `free`", seeds.join(", ")))
    }

    pub(crate) fn solve(&self) -> Result<Vec<usize>> {
        let n = self.goals.len();
        if self.mins.len() != n {
            return Err(Error::Failed("interval assignment: goals and minimums differ in length".into()));
        }
        let mut x = Vec::with_capacity(n);
        for (i, &g) in self.goals.iter().enumerate() {
            if !(g >= 0.0 && g.is_finite()) {
                return Err(Error::Invalid(format!("interval goal {g} is not a finite, non-negative number")));
            }
            if g > MAX_COUNT as f64 {
                return Err(Error::Limit(format!("a curve would get {g:.0} divisions; the limit is {MAX_COUNT}")));
            }
            x.push((g.round() as usize).max(self.mins.get(i).copied().unwrap_or(1)));
        }
        let check = |v: &[usize]| v.iter().all(|&i| i < n);
        if !self.equal.iter().all(|(a, b)| check(a) && check(b))
            || !self.multiple.iter().all(|(l, m)| check(l) && *m >= 2)
            || !self.at_least.iter().all(|(l, _)| check(l))
            || !self.exceed.iter().all(|(a, b, _)| check(a) && check(b))
        {
            return Err(Error::Failed("interval assignment: a constraint names an unknown curve".into()));
        }
        let mut coupled = vec![false; n];
        for (a, b) in &self.equal {
            for &i in a.iter().chain(b) {
                if let Some(c) = coupled.get_mut(i) {
                    *c = true;
                }
            }
        }
        let error = |i: usize, c: usize| -> f64 {
            let g = self.goals.get(i).copied().unwrap_or(1.0).max(MIN_GOAL);
            // Squared relative deviation: strictly convex, so repeated increments spread over
            // the candidates instead of piling onto one.
            ((c as f64 - g) / g).powi(2)
        };
        let step_cost = |x: &[usize], i: usize| -> f64 {
            let c = x.get(i).copied().unwrap_or(0);
            error(i, c + 1) - error(i, c)
        };
        // How often each parity constraint has been repaired.
        let mut repaired = vec![0usize; self.multiple.len()];
        let total: usize = x.iter().sum();
        let constraints = self.equal.len() + self.multiple.len() + self.at_least.len() + self.exceed.len();
        let max_steps = total.saturating_mul(8).saturating_add(64 * (n + constraints)).saturating_add(1000);
        for _ in 0..max_steps {
            let sum = |v: &[usize], x: &[usize]| v.iter().map(|&i| x.get(i).copied().unwrap_or(0)).sum::<usize>();
            // The candidates to raise and whether raising them is what the constraint needs.
            let mut raise: Option<Vec<usize>> = None;
            let mut parity = false;
            for (a, b) in &self.equal {
                let (sa, sb) = (sum(a, &x), sum(b, &x));
                if sa == sb {
                    continue;
                }
                // Raise a variable whose net coefficient is on the smaller side.
                let (small, large) = if sa < sb { (a, b) } else { (b, a) };
                let net = |i: usize| small.iter().filter(|&&k| k == i).count() as i64 - large.iter().filter(|&&k| k == i).count() as i64;
                let cands: Vec<usize> = small.iter().copied().filter(|&i| net(i) > 0).collect();
                if cands.is_empty() {
                    return Err(Error::Failed("interval assignment: a side constraint cannot be met (the same curves on both sides)".into()));
                }
                let cands: Vec<usize> = cands.into_iter().filter(|&i| !self.is_fixed(i)).collect();
                if cands.is_empty() {
                    let all: Vec<usize> = a.iter().chain(b).copied().collect();
                    return Err(self.conflict(&all, &format!("opposite sides of a mapped or swept face need equal counts ({sa} against {sb})")));
                }
                raise = Some(cands);
                break;
            }
            if raise.is_none() {
                for (l, m) in &self.at_least {
                    if sum(l, &x) < *m {
                        let free: Vec<usize> = l.iter().copied().filter(|&i| !self.is_fixed(i)).collect();
                        if free.is_empty() {
                            return Err(self.conflict(l, &format!("a face's loop needs at least {m} divisions")));
                        }
                        raise = Some(free);
                        break;
                    }
                }
            }
            if raise.is_none() {
                for (a, b, d) in &self.exceed {
                    let (sa, sb) = (sum(a, &x), sum(b, &x));
                    if sa < sb.saturating_add(*d) {
                        let net = |i: usize| a.iter().filter(|&&k| k == i).count() as i64 - b.iter().filter(|&&k| k == i).count() as i64;
                        let cands: Vec<usize> = a.iter().copied().filter(|&i| net(i) > 0).collect();
                        if cands.is_empty() {
                            return Err(Error::Failed("interval assignment: a dominance constraint cannot be met".into()));
                        }
                        let free: Vec<usize> = cands.into_iter().filter(|&i| !self.is_fixed(i)).collect();
                        if free.is_empty() {
                            let all: Vec<usize> = a.iter().chain(b).copied().collect();
                            return Err(self.conflict(
                                &all,
                                &format!("a side of a three-sided face needs fewer divisions than the other two ({sb} against {sa})"),
                            ));
                        }
                        raise = Some(free);
                        break;
                    }
                }
            }
            if raise.is_none() {
                for (k, (l, m)) in self.multiple.iter().enumerate() {
                    if !sum(l, &x).is_multiple_of(*m) {
                        // Only a variable listed a number of times that is not a multiple of the
                        // modulus moves the remainder.
                        let moving: Vec<usize> = l.iter().copied().filter(|&i| !l.iter().filter(|&&k| k == i).count().is_multiple_of(*m)).collect();
                        if moving.is_empty() {
                            return Err(Error::Failed("interval assignment: a loop's count cannot reach the multiple it needs".into()));
                        }
                        let moving: Vec<usize> = moving.into_iter().filter(|&i| !self.is_fixed(i)).collect();
                        if moving.is_empty() {
                            return Err(self.conflict(l, &format!("a face's loop needs a multiple of {m} divisions")));
                        }
                        raise = Some(moving);
                        // Loops whose free curves only feed each other (the triangles round a
                        // point) cannot all be made even by them: once a loop keeps coming back,
                        // its curves coupled to other faces may move too.
                        let r = repaired.get_mut(k).ok_or_else(|| Error::Failed("parity index".into()))?;
                        *r += 1;
                        parity = *r <= PARITY_PATIENCE;
                        break;
                    }
                }
            }
            let Some(cands) = raise else { return Ok(x) };
            let pick = cands
                .iter()
                .copied()
                .min_by(|&a, &b| {
                    let key = |i: usize| (parity && coupled.get(i).copied().unwrap_or(false), step_cost(&x, i));
                    let (ka, kb) = (key(a), key(b));
                    ka.0.cmp(&kb.0).then(ka.1.total_cmp(&kb.1)).then(a.cmp(&b))
                })
                .ok_or_else(|| Error::Failed("interval assignment: empty constraint".into()))?;
            let v = x.get_mut(pick).ok_or_else(|| Error::Failed("interval index".into()))?;
            *v += 1;
            if *v > MAX_COUNT {
                return Err(Error::Limit(format!("interval assignment drove a curve past {MAX_COUNT} divisions")));
            }
        }
        Err(Error::Failed(format!(
            "interval assignment did not settle in {max_steps} steps: the side and parity constraints of the faces conflict (the region cannot be mapped or swept as decomposed)"
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounds_goals_without_constraints() {
        let mut p = Intervals::default();
        p.var(2.4, 1);
        p.var(0.2, 1);
        p.var(7.6, 3);
        assert_eq!(p.solve().unwrap(), vec![2, 1, 8]);
    }

    #[test]
    fn mapped_sides_become_equal_and_loops_even() {
        // A rectangle: sides 0 and 2 opposite, 1 and 3 opposite; side 0 shared with a free face
        // whose loop is 0 + 4 + 5.
        let mut p = Intervals::default();
        for g in [4.2, 1.1, 3.6, 1.4, 2.0, 2.0] {
            p.var(g, 1);
        }
        p.equal.push((vec![0], vec![2]));
        p.equal.push((vec![1], vec![3]));
        p.multiple.push((vec![0, 4, 5], 2));
        let x = p.solve().unwrap();
        assert_eq!(x[0], x[2]);
        assert_eq!(x[1], x[3]);
        assert_eq!((x[0] + x[4] + x[5]) % 2, 0);
        // The parity fix went to an uncoupled curve.
        assert_eq!(x[0], 4);
    }

    #[test]
    fn sums_and_minimums() {
        let mut p = Intervals::default();
        let a = p.var(1.0, 1);
        let b = p.var(1.0, 1);
        let c = p.var(5.0, 1);
        p.equal.push((vec![a, b], vec![c]));
        p.at_least.push((vec![a, b, c], 12));
        p.multiple.push((vec![a, b, c], 2));
        let x = p.solve().unwrap();
        assert_eq!(x[a] + x[b], x[c]);
        assert!(x[a] + x[b] + x[c] >= 12);
        assert_eq!((x[a] + x[b] + x[c]) % 2, 0);
    }

    #[test]
    fn exact_seeds_hold_and_conflicts_name_them() {
        use crate::seeds::EdgeSeed;
        use fieldcraft_model::SeedConstraint;
        let seed =
            |name: &str, count: usize, constraint| EdgeSeed { name: name.into(), points: Vec::new(), tolerance: 0.0, count, bias: None, constraint };
        // The other side follows an exact seed.
        let mut p = Intervals::default();
        let a = p.var(3.0, 1);
        let b = p.var(3.0, 1);
        p.equal.push((vec![a], vec![b]));
        p.seed(a, &seed("A", 7, SeedConstraint::Exact));
        assert_eq!(p.solve().unwrap(), vec![7, 7]);
        // Two exact seeds on opposite sides that differ.
        p.seed(b, &seed("B", 5, SeedConstraint::Exact));
        let e = p.solve().unwrap_err().to_string();
        assert!(e.contains("A (7 exact)") && e.contains("B (5 exact)"), "{e}");
        // A minimum and a free goal.
        let mut q = Intervals::default();
        let a = q.var(3.0, 1);
        let b = q.var(3.0, 1);
        q.seed(a, &seed("A", 6, SeedConstraint::Min));
        q.seed(b, &seed("B", 9, SeedConstraint::Free));
        assert_eq!(q.solve().unwrap(), vec![6, 9]);
        // An exact odd count on a loop that must be even, with nothing else to move.
        let mut r = Intervals::default();
        let a = r.var(3.0, 1);
        r.multiple.push((vec![a], 2));
        r.seed(a, &seed("A", 3, SeedConstraint::Exact));
        assert!(r.solve().unwrap_err().to_string().contains("multiple of 2"));
    }

    #[test]
    fn dominance_raises_the_other_sides() {
        // A triangle's sides: 6 against 1 + 2 needs the short sides raised to 8 together.
        let mut p = Intervals::default();
        let a = p.var(6.0, 1);
        let b = p.var(1.0, 1);
        let c = p.var(2.0, 1);
        p.exceed.push((vec![b, c], vec![a], 2));
        p.exceed.push((vec![a, c], vec![b], 2));
        p.exceed.push((vec![a, b], vec![c], 2));
        let x = p.solve().unwrap();
        assert_eq!(x[a], 6);
        assert!(x[b] + x[c] >= 8, "{x:?}");
    }

    /// Four triangles round a point (a square cut along both diagonals): each loop must be even.
    /// The half-diagonals are free but every one of them is in two loops, so they cannot fix an odd
    /// perimeter; a side, coupled to the face beyond, has to move.
    #[test]
    fn parity_round_a_point_moves_a_coupled_side() {
        let mut p = Intervals::default();
        let sides: Vec<usize> = [7.0, 7.0, 7.0, 8.0].iter().map(|&g| p.var(g, 1)).collect();
        let halves: Vec<usize> = (0..4).map(|_| p.var(2.0, 1)).collect();
        let beyond: Vec<usize> = (0..4).map(|_| p.var(7.0, 1)).collect();
        for k in 0..4 {
            p.multiple.push((vec![sides[k], halves[k], halves[(k + 1) % 4]], 2));
            p.equal.push((vec![sides[k]], vec![beyond[k]]));
        }
        let x = p.solve().unwrap();
        for (l, m) in &p.multiple {
            assert_eq!(l.iter().map(|&i| x[i]).sum::<usize>() % m, 0, "{x:?}");
        }
    }

    #[test]
    fn infeasible_systems_are_errors_not_loops() {
        let mut p = Intervals::default();
        let a = p.var(3.0, 1);
        let b = p.var(3.0, 1);
        // a = a + b has no solution with b ≥ 1.
        p.equal.push((vec![a], vec![a, b]));
        assert!(p.solve().is_err());
        let mut q = Intervals::default();
        let a = q.var(3.0, 1);
        q.multiple.push((vec![a, a], 2));
        q.multiple.push((vec![a, a, a], 2));
        q.multiple.push((vec![a, a, a], 4));
        assert!(q.solve().is_ok());
        let mut r = Intervals::default();
        r.var(f64::NAN, 1);
        assert!(r.solve().is_err());
        let mut s = Intervals::default();
        s.var(1e30, 1);
        assert!(s.solve().is_err());
    }
}
