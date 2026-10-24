//! Bounded inverse isoparametric maps for sampling solid-element nodal fields.
//!
//! Geometry and the sampled field use the same shape functions (K.-J. Bathe, MIT 2.094,
//! lecture 7, equations 7.1--7.10, https://ocw.mit.edu/courses/2-094-finite-element-analysis-of-solids-and-fluids-ii-spring-2011/resources/mit2_094s11_lec7/).
//! Newton iteration uses normalized local coordinates and a line search. Failure to find a
//! root is **not** evidence that a point is outside: interval subdivision must exclude the
//! whole reference domain, otherwise this module returns an actionable unresolved error.
//!
//! Quadratic Lagrange nodes do not enclose their curved element. Tet10 bounds use quadratic
//! simplex Bernstein controls; hex20 is embedded in tensor Q2 and converted to its 27
//! Bernstein controls. Their nonnegative partition of unity proves the enclosure. See
//! T. W. Sederberg, *Computer Aided Geometric Design*, July 25, 2017, sections 2.5, 12.1,
//! 15.1 (author metadata https://scholarsarchive.byu.edu/facpub/1/; accessible 2017 text:
//! https://flurry.dg.fmph.uniba.sk/webog/SuboryOG/bohdal/ComputerAidedGeometricDesign.pdf).
//! All enclosure arithmetic rounds outward; subdivision is bounded in depth and work.
//! Accepted points have a locally admissible Jacobian (except the pyramid's singular apex
//! chart). This query check is not a proof of global element injectivity; mesh validity
//! remains the model/mesher's responsibility.

use crate::{Basis, Error, Result, dense::positive_inverse, shape::Shape};
use fieldcraft_model::ElementKind;

const RESIDUAL: f64 = 64.0 * f64::EPSILON;
const INSIDE: f64 = 256.0 * f64::EPSILON;
const ITERATIONS: usize = 32;
const LINE_SEARCH: usize = 12;
const EXCLUSION_BOXES: usize = 2048;
const EXCLUSION_DEPTH: usize = 24;

/// Shared across every candidate and station of a query. A unit is one basis evaluation
/// or one interval enclosure (at most 27 nodes); callers also cap index traversal/memory.
#[derive(Debug)]
pub struct Work {
    remaining: usize,
}

impl Work {
    pub fn new(limit: usize) -> Self {
        Self { remaining: limit }
    }

    pub fn remaining(&self) -> usize {
        self.remaining
    }

    fn spend(&mut self) -> Result<()> {
        self.spend_many(1)
    }

    fn spend_many(&mut self, amount: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| Error::Invalid("path location work limit reached; reduce samples or split the path".into()))?;
        Ok(())
    }
}

/// Conservative world-coordinate enclosure, including curved higher-order geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

impl Bounds {
    pub fn contains(self, point: [f64; 3]) -> bool {
        point.iter().enumerate().all(|(k, p)| p.is_finite() && *p >= self.min[k] && *p <= self.max[k])
    }

    /// Candidate queries admit the same scale-relative roundoff band as the inverse map's
    /// reference-domain test. This does not alter the stored conservative enclosure. It
    /// prevents arc-length/units arithmetic at an exterior face from dropping a boundary
    /// station before Newton can validate its residual and natural coordinates.
    pub fn contains_near(self, point: [f64; 3]) -> bool {
        let span = (0..3).map(|k| self.max[k] - self.min[k]).fold(0.0, f64::max);
        let tolerance = if span.is_finite() { INSIDE * span } else { 0.0 };
        point.iter().enumerate().all(|(k, p)| p.is_finite() && *p >= self.min[k] - tolerance && *p <= self.max[k] + tolerance)
    }

    /// Min/max selection introduces no floating-point arithmetic or inward rounding.
    pub fn union(self, other: Self) -> Self {
        Self { min: [0, 1, 2].map(|k| self.min[k].min(other.min[k])), max: [0, 1, 2].map(|k| self.max[k].max(other.max[k])) }
    }
}

/// Natural coordinates and exactly the existing basis used to interpolate nodal scalars.
#[derive(Debug, Clone, Copy)]
pub struct Location {
    pub natural: [f64; 3],
    pub shape: Shape,
}

impl Location {
    /// Undefined contributing nodal data stays undefined. Exact zero weights do not use a
    /// missing value (so a defined value at an exact node/face remains defined). No nearest-
    /// node substitution or renormalization is performed.
    pub fn interpolate(self, values: &[f64]) -> Option<f64> {
        if values.len() != self.shape.count {
            return None;
        }
        let mut value = 0.0;
        for (n, v) in self.shape.values().iter().zip(values) {
            if *n != 0.0 {
                if !v.is_finite() {
                    return None;
                }
                value += n * v;
            }
        }
        value.is_finite().then_some(value)
    }
}

fn validate(kind: ElementKind, xyz: &[[f64; 3]]) -> Result<Basis> {
    let basis = Basis::of(kind);
    if basis.dim() != 3 {
        return Err(Error::NotSupported("path interpolation requires solid elements (tet4/tet10/wedge6/pyramid5/hex8/hex20)".into()));
    }
    if xyz.len() != basis.count() || xyz.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("path element coordinates must be finite and match the element node count".into()));
    }
    Ok(basis)
}

/// Bounds of the mathematical shape-function map over the closed reference element.
/// A standalone call costs at most 27 basis evaluations plus one fixed-size enclosure;
/// index builders must cap their element count. `locate` charges this cost before work.
pub fn bounds(kind: ElementKind, xyz: &[[f64; 3]]) -> Result<Bounds> {
    let basis = validate(kind, xyz)?;
    let origin = xyz.first().copied().ok_or_else(|| Error::Invalid("empty path element".into()))?;
    let mut local = [[Interval::ZERO; 3]; 27];
    for (row, x) in local.iter_mut().zip(xyz) {
        *row = [0, 1, 2].map(|k| Interval::point(x[k]).sub(Interval::point(origin[k])));
    }
    let mut controls = local;
    let count = match kind {
        ElementKind::Tet10 => {
            for edge in kind.edges() {
                if let [a, b, m] = **edge
                    && let (Some(xa), Some(xb), Some(xm)) = (local.get(a), local.get(b), local.get(m))
                    && let Some(control) = controls.get_mut(m)
                {
                    *control = [0, 1, 2].map(|k| xm[k].scale(2.0).sub(xa[k].add(xb[k]).scale(0.5)));
                }
            }
            10
        }
        ElementKind::Hex20 => {
            // At {-1,0,1}^3 every shape coefficient is an exactly represented dyadic.
            // Grid index is x + 3*y + 9*z. Convert each axis from values to controls.
            for z in 0..3 {
                for y in 0..3 {
                    for x in 0..3 {
                        let shape = basis.eval([x as f64 - 1.0, y as f64 - 1.0, z as f64 - 1.0]);
                        controls[x + 3 * y + 9 * z] = interval_position(&local[..20], shape.values());
                    }
                }
            }
            for stride in [1, 3, 9] {
                let previous = controls;
                for i in 0..27 {
                    if i / stride % 3 == 1 {
                        controls[i] =
                            [0, 1, 2].map(|k| previous[i][k].scale(2.0).sub(previous[i - stride][k].add(previous[i + stride][k]).scale(0.5)));
                    }
                }
            }
            27
        }
        _ => basis.count(), // These four families have nonnegative reference weights.
    };
    let mut out = Bounds { min: [f64::INFINITY; 3], max: [f64::NEG_INFINITY; 3] };
    for control in controls.iter().take(count) {
        for k in 0..3 {
            let world = control[k].add(Interval::point(origin[k]));
            if !world.lo.is_finite() || !world.hi.is_finite() {
                return Err(Error::Invalid("path element bounds overflow; use finite coordinates with a representable span".into()));
            }
            out.min[k] = out.min[k].min(world.lo);
            out.max[k] = out.max[k].max(world.hi);
        }
    }
    if out.min.iter().chain(&out.max).any(|x| !x.is_finite()) {
        return Err(Error::Invalid("path element bounds overflow; use finite coordinates with a representable span".into()));
    }
    Ok(out)
}

/// Locate a point, or prove it outside. Unresolved searches and exhausted work are errors,
/// never fabricated gaps. Callers choose deterministic element labels on shared faces.
pub fn locate(kind: ElementKind, xyz: &[[f64; 3]], point: [f64; 3], work: &mut Work) -> Result<Option<Location>> {
    let basis = validate(kind, xyz)?;
    if point.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("path location point must be finite".into()));
    }
    // Debit before any bounds evaluation, including queries outside the enclosure.
    work.spend_many(if kind == ElementKind::Hex20 { 28 } else { 1 })?;
    if !bounds(kind, xyz)?.contains_near(point) {
        return Ok(None);
    }
    let origin = xyz.first().copied().ok_or_else(|| Error::Invalid("empty path element".into()))?;
    let span = xyz.iter().flat_map(|x| (0..3).map(move |k| (x[k] - origin[k]).abs())).fold(0.0, f64::max);
    if !(span > 0.0 && span.is_finite()) {
        return Err(Error::Degenerate("path element has no finite, nonzero extent".into()));
    }
    let mut local = [[0.0; 3]; 27];
    for (row, x) in local.iter_mut().zip(xyz) {
        *row = [0, 1, 2].map(|k| (x[k] - origin[k]) / span);
    }
    let target = [0, 1, 2].map(|k| (point[k] - origin[k]) / span);
    if target.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("path coordinates cannot be normalized".into()));
    }
    let coords = local.get(..basis.count()).ok_or_else(|| Error::Invalid("invalid path basis".into()))?;
    let center = center(kind);
    let (_, j) = evaluate(basis, coords, center, work)?;
    if positive_inverse(&j).is_none() {
        return Err(Error::Degenerate("path element has an invalid center Jacobian".into()));
    }
    // Exact coordinate equality only: this is not a nearest-node fallback. In particular
    // it handles the pyramid apex, where the natural coordinate chart is singular.
    if let Some(i) = xyz.iter().position(|x| *x == point)
        && let Some(natural) = basis.nodes().get(i).copied()
    {
        work.spend()?;
        let shape = basis.eval(natural);
        let (_, j) = geometry(&shape, coords);
        if !(kind == ElementKind::Pyramid5 && natural[2] == 1.0) && positive_inverse(&j).is_none() {
            return Err(Error::Degenerate("path element has an invalid Jacobian at the queried node".into()));
        }
        return Ok(Some(Location { natural, shape }));
    }
    for seed in
        std::iter::once(center).chain(basis.nodes().iter().take(basis.corners()).map(|node| [0, 1, 2].map(|k| 0.75 * center[k] + 0.25 * node[k])))
    {
        if let Some(natural) = newton(basis, coords, target, seed, work)?
            && inside(kind, natural, INSIDE)
        {
            work.spend()?;
            return Ok(Some(Location { natural, shape: basis.eval(natural) }));
        }
    }
    // Interval calculations use world-relative (unscaled) coordinates and target intervals,
    // so normalization roundoff cannot establish a false outside classification.
    if exclude(kind, xyz, point, work)? {
        Ok(None)
    } else {
        Err(Error::Invalid(
            "path point location remains unresolved after bounded inverse/subdivision search; refine the distorted element or move the path".into(),
        ))
    }
}

fn center(kind: ElementKind) -> [f64; 3] {
    match kind {
        ElementKind::Tet4 | ElementKind::Tet10 => [0.25; 3],
        ElementKind::Wedge6 => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        ElementKind::Pyramid5 => [0.0, 0.0, 0.25],
        _ => [0.0; 3],
    }
}

fn inside(kind: ElementKind, xi: [f64; 3], tolerance: f64) -> bool {
    let [r, s, t] = xi;
    match kind {
        ElementKind::Tet4 | ElementKind::Tet10 => r >= -tolerance && s >= -tolerance && t >= -tolerance && r + s + t <= 1.0 + tolerance,
        ElementKind::Wedge6 => r >= -tolerance && s >= -tolerance && r + s <= 1.0 + tolerance && t.abs() <= 1.0 + tolerance,
        ElementKind::Pyramid5 => t >= -tolerance && t <= 1.0 + tolerance && r.abs().max(s.abs()) <= 1.0 - t + tolerance,
        _ => xi.iter().all(|x| x.abs() <= 1.0 + tolerance),
    }
}

fn evaluate(basis: Basis, xyz: &[[f64; 3]], xi: [f64; 3], work: &mut Work) -> Result<([f64; 3], [[f64; 3]; 3])> {
    work.spend()?;
    let shape = basis.eval(xi);
    Ok(geometry(&shape, xyz))
}

fn geometry(shape: &Shape, xyz: &[[f64; 3]]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut position = [0.0; 3];
    let mut jacobian = [[0.0; 3]; 3];
    for ((n, dn), x) in shape.values().iter().zip(shape.derivatives()).zip(xyz) {
        for k in 0..3 {
            position[k] += n * x[k];
            for a in 0..3 {
                jacobian[a][k] += dn[a] * x[k];
            }
        }
    }
    (position, jacobian)
}

fn norm(x: [f64; 3]) -> f64 {
    x.into_iter().map(f64::abs).fold(0.0, f64::max)
}

fn newton(basis: Basis, xyz: &[[f64; 3]], target: [f64; 3], mut xi: [f64; 3], work: &mut Work) -> Result<Option<[f64; 3]>> {
    for _ in 0..ITERATIONS {
        let (position, jacobian) = evaluate(basis, xyz, xi, work)?;
        let residual = [0, 1, 2].map(|k| target[k] - position[k]);
        let Some((inverse, _)) = positive_inverse(&jacobian) else { return Ok(None) };
        let error = norm(residual);
        if error <= RESIDUAL {
            return Ok(Some(xi));
        }
        // The existing convention is J[a][k] = dx_k/dxi_a: invert its transpose.
        let delta = [0, 1, 2].map(|a| (0..3).map(|k| inverse[k][a] * residual[k]).sum::<f64>());
        let mut factor = 1.0;
        let mut accepted = false;
        for _ in 0..LINE_SEARCH {
            let trial = [0, 1, 2].map(|k| xi[k] + factor * delta[k]);
            if trial.iter().all(|x| x.is_finite() && x.abs() <= 4.0) {
                let (mapped, j) = evaluate(basis, xyz, trial, work)?;
                let trial_error = norm([0, 1, 2].map(|k| target[k] - mapped[k]));
                if positive_inverse(&j).is_some() && (trial_error < error || trial_error <= RESIDUAL) {
                    xi = trial;
                    accepted = true;
                    break;
                }
            }
            factor *= 0.5;
        }
        if !accepted {
            return Ok(None);
        }
    }
    Ok(None)
}

#[derive(Debug, Clone, Copy)]
struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    const ZERO: Self = Self { lo: 0.0, hi: 0.0 };
    const ONE: Self = Self { lo: 1.0, hi: 1.0 };

    fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn rounded(lo: f64, hi: f64) -> Self {
        // Overflow/indeterminate arithmetic may weaken exclusion, never prove it.
        if lo.is_nan() || hi.is_nan() {
            return Self { lo: f64::NEG_INFINITY, hi: f64::INFINITY };
        }
        Self { lo: lo.next_down(), hi: hi.next_up() }
    }
    fn add(self, b: Self) -> Self {
        Self::rounded(self.lo + b.lo, self.hi + b.hi)
    }
    fn sub(self, b: Self) -> Self {
        Self::rounded(self.lo - b.hi, self.hi - b.lo)
    }
    fn scale(self, b: f64) -> Self {
        if b == 0.0 {
            return Self::ZERO;
        }
        let (lo, hi) = if b >= 0.0 { (self.lo * b, self.hi * b) } else { (self.hi * b, self.lo * b) };
        Self::rounded(lo, hi)
    }
    fn mul(self, b: Self) -> Self {
        let products = [self.lo * b.lo, self.lo * b.hi, self.hi * b.lo, self.hi * b.hi];
        if products.iter().any(|x| x.is_nan()) {
            return Self { lo: f64::NEG_INFINITY, hi: f64::INFINITY };
        }
        Self::rounded(products.into_iter().fold(f64::INFINITY, f64::min), products.into_iter().fold(f64::NEG_INFINITY, f64::max))
    }
    fn square(self) -> Self {
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 { 0.0 } else { (self.lo * self.lo).min(self.hi * self.hi).next_down() };
        Self { lo, hi: (self.lo * self.lo).max(self.hi * self.hi).next_up() }
    }
}

fn interval_position(xyz: &[[Interval; 3]], weights: &[f64]) -> [Interval; 3] {
    let mut out = [Interval::ZERO; 3];
    for (x, weight) in xyz.iter().zip(weights) {
        for k in 0..3 {
            out[k] = out[k].add(x[k].scale(*weight));
        }
    }
    out
}

// Independent interval extensions of the existing reference functions. This is enclosure
// arithmetic only; every accepted field value uses Basis::eval, never this extension.
fn interval_weights(kind: ElementKind, xi: [Interval; 3]) -> [Interval; 27] {
    let [r, s, t] = xi;
    let one = Interval::ONE;
    let mut n = [Interval::ZERO; 27];
    match kind {
        ElementKind::Tet4 | ElementKind::Tet10 => {
            let l = [one.sub(r).sub(s).sub(t), r, s, t];
            for (out, li) in n.iter_mut().zip(l) {
                *out = if kind == ElementKind::Tet4 { li } else { li.mul(li.scale(2.0).sub(one)) };
            }
            if kind == ElementKind::Tet10 {
                for edge in kind.edges() {
                    if let [a, b, m] = **edge
                        && let (Some(la), Some(lb), Some(out)) = (l.get(a), l.get(b), n.get_mut(m))
                    {
                        *out = la.mul(*lb).scale(4.0);
                    }
                }
            }
        }
        ElementKind::Wedge6 => {
            let l = [one.sub(r).sub(s), r, s];
            for i in 0..3 {
                n[i] = l[i].mul(one.sub(t)).scale(0.5);
                n[i + 3] = l[i].mul(one.add(t)).scale(0.5);
            }
        }
        ElementKind::Pyramid5 => {
            // Collapsed cube coordinates r,s in [-1,1], t in [0,1] avoid division at
            // the apex; the equivalent natural point is [r*(1-t),s*(1-t),t].
            for (i, corner) in crate::shape::HEX_CORNERS.iter().take(4).enumerate() {
                n[i] = one.add(r.scale(corner[0])).mul(one.add(s.scale(corner[1]))).mul(one.sub(t)).scale(0.25);
            }
            n[4] = t;
        }
        ElementKind::Hex8 | ElementKind::Hex20 => {
            for (i, node) in Basis::of(kind).nodes().iter().enumerate() {
                let [a, b, c] = *node;
                let value = if kind == ElementKind::Hex8 {
                    one.add(r.scale(a)).mul(one.add(s.scale(b))).mul(one.add(t.scale(c))).scale(0.125)
                } else if a == 0.0 {
                    one.sub(r.square()).mul(one.add(s.scale(b))).mul(one.add(t.scale(c))).scale(0.25)
                } else if b == 0.0 {
                    one.sub(s.square()).mul(one.add(r.scale(a))).mul(one.add(t.scale(c))).scale(0.25)
                } else if c == 0.0 {
                    one.sub(t.square()).mul(one.add(r.scale(a))).mul(one.add(s.scale(b))).scale(0.25)
                } else {
                    one.add(r.scale(a))
                        .mul(one.add(s.scale(b)))
                        .mul(one.add(t.scale(c)))
                        .mul(r.scale(a).add(s.scale(b)).add(t.scale(c)).sub(Interval::point(2.0)))
                        .scale(0.125)
                };
                if let Some(out) = n.get_mut(i) {
                    *out = value;
                }
            }
        }
        _ => {}
    }
    n
}

fn exclude(kind: ElementKind, xyz: &[[f64; 3]], point: [f64; 3], work: &mut Work) -> Result<bool> {
    let origin = xyz.first().copied().ok_or_else(|| Error::Invalid("empty path element".into()))?;
    let mut coords = [[Interval::ZERO; 3]; 27];
    for (out, x) in coords.iter_mut().zip(xyz) {
        *out = [0, 1, 2].map(|k| Interval::point(x[k]).sub(Interval::point(origin[k])));
    }
    let target = [0, 1, 2].map(|k| Interval::point(point[k]).sub(Interval::point(origin[k])));
    let initial = match kind {
        ElementKind::Tet4 | ElementKind::Tet10 => [Interval { lo: 0.0, hi: 1.0 }; 3],
        ElementKind::Wedge6 => [Interval { lo: 0.0, hi: 1.0 }, Interval { lo: 0.0, hi: 1.0 }, Interval { lo: -1.0, hi: 1.0 }],
        ElementKind::Pyramid5 => [Interval { lo: -1.0, hi: 1.0 }, Interval { lo: -1.0, hi: 1.0 }, Interval { lo: 0.0, hi: 1.0 }],
        _ => [Interval { lo: -1.0, hi: 1.0 }; 3],
    };
    // Depth-first traversal needs at most depth+1 pending siblings, not an unbounded tree.
    let mut stack = Vec::with_capacity(EXCLUSION_DEPTH + 1);
    stack.push((initial, 0));
    for _ in 0..EXCLUSION_BOXES {
        let Some((cell, depth)) = stack.pop() else { return Ok(true) };
        work.spend()?;
        let simplex_min = Interval::point(cell[0].lo).add(Interval::point(cell[1].lo));
        if matches!(kind, ElementKind::Tet4 | ElementKind::Tet10) && simplex_min.add(Interval::point(cell[2].lo)).lo > 1.0
            || kind == ElementKind::Wedge6 && simplex_min.lo > 1.0
        {
            continue;
        }
        let weights = interval_weights(kind, cell);
        let mut mapped = [Interval::ZERO; 3];
        for (x, n) in coords.iter().zip(weights).take(xyz.len()) {
            for k in 0..3 {
                mapped[k] = mapped[k].add(x[k].mul(n));
            }
        }
        if (0..3).any(|k| target[k].lo > mapped[k].hi || target[k].hi < mapped[k].lo) {
            continue;
        }
        if depth >= EXCLUSION_DEPTH {
            return Ok(false);
        }
        let axis = (0..3).max_by(|a, b| (cell[*a].hi - cell[*a].lo).total_cmp(&(cell[*b].hi - cell[*b].lo))).unwrap_or(0);
        let mid = cell[axis].lo + (cell[axis].hi - cell[axis].lo) * 0.5;
        if !(mid > cell[axis].lo && mid < cell[axis].hi) {
            return Ok(false);
        }
        let mut low = cell;
        let mut high = cell;
        low[axis].hi = mid;
        high[axis].lo = mid;
        stack.push((high, depth + 1));
        stack.push((low, depth + 1));
    }
    Ok(stack.is_empty())
}

#[cfg(test)]
mod tests;
