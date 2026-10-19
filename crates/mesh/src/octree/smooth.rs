//! Untangling, fitting and shape optimisation of the hexahedral mesh.
//!
//! # Energy
//! The foldover-free mapping energy of V. Garanzha, I. Kaporin, L. Kudryavtseva, F. Protais,
//! N. Ray, D. Sokolov, *Foldover-free maps in 50 lines of code*, ACM TOG 40(4):102, 2021 (paper
//! only; the authors' code is AGPL and is not read). Each hexahedron is measured at nine points,
//! its eight corners (the trapezoidal rule the paper uses for quads) and its centre, where a
//! twisted element folds first: `J` is the sample's three edge directions against the ideal cube
//! of the element's target edge `s`, `D = det J`, and the sample costs
//!
//! `f_ε + λ g_ε`, `f_ε = |J|³_F / χ(D, ε)`, `g_ε = (D² + 1) / χ(D, ε)`,
//! `χ(D, ε) = (D + √(ε² + D²)) / 2`,
//!
//! weighted by `s³`. `f` is the shape (conformal) term, `g` the volume term. The paper's shape
//! term is `|J|²/D^{2/3}`; its 3/2 power `|J|³/D`, the classical polyconvex distortion in three
//! dimensions (J. M. Ball, 1977), is minimal at the same ideal shapes and needs only square
//! roots, so the optimiser gives the same bits on every platform (no `pow` of the system maths
//! library). `χ` keeps both finite for inverted samples, and `ε` follows the paper's rule (7),
//! `ε = √(10⁻¹² + 0.04 min(0, min D)²)`, so it vanishes once every sample is positive and the
//! energy becomes a barrier against inversion.
//!
//! A vertex on a corner, curve or patch adds `fit · 8 s̄ · dist²` (`s̄` the mean target edge of its
//! elements) to its distance from that feature: with `fit` raised stage by stage the boundary is
//! drawn onto the geometry while the barrier keeps every element valid (Gao, Shen, Panozzo 2019,
//! §3.5). With `fit = 0` constrained vertices are held on their features instead (each step is
//! restricted to the tangent and projected back).
//!
//! # Minimisation
//! Nonlinear Gauss–Seidel: each vertex in turn takes a Newton step on its own energy with the
//! paper's positive approximation of the Hessian (the second derivatives of `χ` and `det` are
//! dropped, §2.3, App. A), accepted only if the energy decreased (halving line search). The paper
//! minimises the same energy globally with L-BFGS; the vertex sweep is the local form of the same
//! descent. Then a compass search raises the worst corner scaled Jacobians, and another repairs
//! elements whose trilinear map the certified check does not prove positive. Sweeps, line-search
//! steps, compass steps and rounds are capped.

use fieldcraft_model::ElementKind;

use crate::quality::{Validity, element_validity};
use crate::vec3::{P, cross, dot, norm, sub};

use super::geom::dist2;

use super::geom::Geometry;

/// The square of the least regularisation of a normalized corner determinant.
const MIN_EPS2: f64 = 1e-12; // unit-factors: the square of a normalized determinant
/// A step this small, relative to the reach, ends the search.
const MIN_STEP: f64 = 1e-9; // unit-factors: relative to the reach
/// Qualities this close are equal.
const SAME_QUALITY: f64 = 1e-12; // unit-factors: a difference of normalized qualities
/// The score of an inverted corner.
const WORST: f64 = -1e9; // unit-factors: a score below every quality

/// How a vertex may move.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Constraint {
    Free,
    /// On a corner of the part (the boundary surface's vertex).
    Corner(u32),
    /// On a feature curve.
    Curve(u32),
    /// On a patch of the boundary.
    Patch(u32),
}

/// Model-order corner `c` of a hex8: its three neighbours along x, y, z and the signs that make
/// the reference corner frame right-handed.
const CORNER_FRAMES: [([usize; 3], [f64; 3]); 8] = [
    ([1, 3, 4], [1.0, 1.0, 1.0]),
    ([0, 2, 5], [-1.0, 1.0, 1.0]),
    ([3, 1, 6], [-1.0, -1.0, 1.0]),
    ([2, 0, 7], [1.0, -1.0, 1.0]),
    ([5, 7, 0], [1.0, 1.0, -1.0]),
    ([4, 6, 1], [-1.0, 1.0, -1.0]),
    ([7, 5, 2], [-1.0, -1.0, -1.0]),
    ([6, 4, 3], [1.0, -1.0, -1.0]),
];

/// Weight of the volume term.
const LAMBDA: f64 = 0.1;

fn chi(d: f64, eps: f64) -> f64 {
    0.5 * (d + (eps * eps + d * d).sqrt())
}

fn chi_d(d: f64, eps: f64) -> f64 {
    0.5 * (1.0 + d / (eps * eps + d * d).sqrt())
}

/// The corner Jacobian (columns) of corner `c` of `x` against an ideal cube of edge `s`.
fn corner_jacobian(x: &[P; 8], c: usize, s: f64) -> [P; 3] {
    let (n, sg) = CORNER_FRAMES[c];
    [0, 1, 2].map(|a| sub(x[n[a]], x[c]).map(|v| v * sg[a] / s))
}

/// The sample points of the energy: the eight corners and the centre. Column `a` of sample `s`'s
/// Jacobian is `Σ_k W[s][a][k] x_k` (before dividing by the target edge). The centre keeps
/// twisted elements (positive corners, folded inside) out of reach, as the nine-point measure of
/// ParaView and of Zhang et al.'s meshers does.
const SAMPLES: usize = 9;

fn sample_weights(s: usize) -> [[f64; 8]; 3] {
    let mut w = [[0.0; 8]; 3];
    if s < 8 {
        let (n, sg) = CORNER_FRAMES[s];
        for a in 0..3 {
            w[a][n[a]] += sg[a];
            w[a][s] -= sg[a];
        }
    } else {
        w = [
            [-0.25, 0.25, 0.25, -0.25, -0.25, 0.25, 0.25, -0.25],
            [-0.25, -0.25, 0.25, 0.25, -0.25, -0.25, 0.25, 0.25],
            [-0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25],
        ];
    }
    w
}

/// The Jacobian (columns) of sample `smp` of `x` against an ideal cube of edge `s`.
fn sample_jacobian(x: &[P; 8], smp: usize, s: f64) -> [P; 3] {
    let w = sample_weights(smp);
    [0, 1, 2].map(|a| {
        let mut c = [0.0; 3];
        for (k, p) in x.iter().enumerate() {
            for r in 0..3 {
                c[r] += w[a][k] * p[r];
            }
        }
        c.map(|v| v / s)
    })
}

fn det3(j: &[P; 3]) -> f64 {
    dot(j[0], cross(j[1], j[2]))
}

/// The cofactor columns: `∂ det / ∂ column a`.
fn cofactors(j: &[P; 3]) -> [P; 3] {
    [cross(j[1], j[2]), cross(j[2], j[0]), cross(j[0], j[1])]
}

/// The energy of one corner.
fn corner_energy(j: &[P; 3], eps: f64) -> f64 {
    let d = det3(j);
    let c = chi(d, eps);
    let fro: f64 = j.iter().map(|v| dot(*v, *v)).sum();
    fro * fro.sqrt() / c + LAMBDA * (d * d + 1.0) / c
}

/// The scaled Jacobian of each corner, smallest (`min over corners of det / product of edge
/// lengths`), the usual hex quality measure.
pub(crate) fn min_scaled_jacobian(x: &[P; 8]) -> f64 {
    (0..8)
        .map(|c| {
            let j = corner_jacobian(x, c, 1.0);
            let l = norm(j[0]) * norm(j[1]) * norm(j[2]);
            if l > 0.0 { det3(&j) / l } else { -1.0 }
        })
        .fold(f64::INFINITY, f64::min)
}

/// The optimiser's state.
pub(crate) struct Optimizer<'a> {
    pub geom: &'a Geometry<'a>,
    pub points: Vec<P>,
    pub hexes: &'a [[u32; 8]],
    /// Target edge per element.
    pub size: &'a [f64],
    pub constraint: Vec<Constraint>,
    /// Elements around each vertex with the vertex's local position.
    around: Vec<Vec<(u32, u8)>>,
    pub eps: f64,
    /// Weight of the fitting term; 0 holds boundary vertices on their features exactly.
    pub fit: f64,
}

impl<'a> Optimizer<'a> {
    pub fn new(geom: &'a Geometry<'a>, points: Vec<P>, hexes: &'a [[u32; 8]], size: &'a [f64], constraint: Vec<Constraint>) -> Optimizer<'a> {
        let mut around: Vec<Vec<(u32, u8)>> = vec![Vec::new(); points.len()];
        for (e, h) in hexes.iter().enumerate() {
            for (k, &v) in h.iter().enumerate() {
                if let Some(a) = around.get_mut(v as usize) {
                    a.push((e as u32, k as u8));
                }
            }
        }
        Optimizer { geom, points, hexes, size, constraint, around, eps: 1.0, fit: 0.0 }
    }

    fn element(&self, e: u32) -> ([P; 8], f64) {
        let h = self.hexes.get(e as usize).copied().unwrap_or([0; 8]);
        (h.map(|v| self.points.get(v as usize).copied().unwrap_or([0.0; 3])), self.size.get(e as usize).copied().unwrap_or(1.0))
    }

    /// The smallest normalised corner determinant over the mesh, and the total energy.
    pub fn measure(&self) -> (f64, f64) {
        let (mut m, mut total) = (f64::INFINITY, 0.0);
        for e in 0..self.hexes.len() as u32 {
            let (x, s) = self.element(e);
            for c in 0..SAMPLES {
                let j = sample_jacobian(&x, c, s);
                m = m.min(det3(&j));
                total += s * s * s * corner_energy(&j, self.eps);
            }
        }
        (m, total)
    }

    /// Where vertex `v`'s feature is closest to `p` (`None` for a free vertex).
    pub fn target(&self, v: u32, p: P) -> Option<P> {
        match self.constraint.get(v as usize).copied().unwrap_or(Constraint::Free) {
            Constraint::Free => None,
            Constraint::Corner(c) => Some(self.geom.xyz(c)),
            Constraint::Curve(c) => self.geom.nearest_on_curve(c, p).map(|x| x.0),
            Constraint::Patch(pt) => self.geom.nearest(p, |t| self.geom.tri_patch.get(t as usize) == Some(&pt)).map(|x| x.1),
        }
    }

    /// The fitting weight of vertex `v`: the fit factor times its elements' size (so a distance
    /// of one element weighs like the shape of its corners), ten times more for a corner.
    fn fit_weight(&self, v: u32) -> f64 {
        let a = self.around.get(v as usize).map_or(&[][..], |a| a.as_slice());
        if a.is_empty() || self.fit == 0.0 {
            return 0.0;
        }
        let s = a.iter().map(|&(e, _)| self.size.get(e as usize).copied().unwrap_or(1.0)).sum::<f64>() / a.len() as f64;
        let corner = if matches!(self.constraint.get(v as usize), Some(Constraint::Corner(_))) { 10.0 } else { 1.0 };
        8.0 * self.fit * corner * s
    }

    /// The energy of the corners that involve vertex `v`, with its fitting term.
    fn local_energy(&self, v: u32) -> f64 {
        let mut total = 0.0;
        if self.fit > 0.0 {
            let p = self.points.get(v as usize).copied().unwrap_or([0.0; 3]);
            if let Some(t) = self.target(v, p) {
                total += self.fit_weight(v) * dist2(p, t);
            }
        }
        for &(e, k) in self.around.get(v as usize).map_or(&[][..], |a| a.as_slice()) {
            let (x, s) = self.element(e);
            let (n, _) = CORNER_FRAMES[k as usize];
            for c in [k as usize, n[0], n[1], n[2], 8] {
                total += s * s * s * corner_energy(&sample_jacobian(&x, c, s), self.eps);
            }
        }
        total
    }

    /// The gradient and approximate Hessian of the local energy at vertex `v`.
    fn derivatives(&self, v: u32) -> (P, [[f64; 3]; 3]) {
        let mut g = [0.0; 3];
        let mut h = [[0.0; 3]; 3];
        for &(e, k) in self.around.get(v as usize).map_or(&[][..], |a| a.as_slice()) {
            let (x, s) = self.element(e);
            let k = k as usize;
            let w = s * s * s;
            let (nk, _) = CORNER_FRAMES[k];
            for c in [k, nk[0], nk[1], nk[2], 8] {
                // ∂ column a / ∂ x_v = coef[a] I.
                let wt = sample_weights(c);
                let coef: [f64; 3] = [0, 1, 2].map(|a| wt[a][k] / s);
                let j = sample_jacobian(&x, c, s);
                let d = det3(&j);
                let ch = chi(d, self.eps);
                let chd = chi_d(d, self.eps);
                let fro: f64 = j.iter().map(|q| dot(*q, *q)).sum();
                let cof = cofactors(&j);
                let root = fro.sqrt();
                // dE/dJ = A J + B cof(J).
                let a_coef = 3.0 * root / ch;
                let f_d = -fro * root * chd / (ch * ch);
                let g_d = 2.0 * d / ch - (d * d + 1.0) * chd / (ch * ch);
                let b_coef = f_d + LAMBDA * g_d;
                let mut dd = [0.0; 3];
                let mut sum_c2 = 0.0;
                for a in 0..3 {
                    if coef[a] == 0.0 {
                        continue;
                    }
                    sum_c2 += coef[a] * coef[a];
                    for r in 0..3 {
                        g[r] += w * coef[a] * (a_coef * j[a][r] + b_coef * cof[a][r]);
                        dd[r] += coef[a] * cof[a][r];
                    }
                }
                // Positive part of the second derivative along D, χ'' and det'' dropped.
                let f_dd = 2.0 * fro * root * chd * chd / (ch * ch * ch);
                let g_dd = 2.0 / ch - 4.0 * d * chd / (ch * ch) + 2.0 * (d * d + 1.0) * chd * chd / (ch * ch * ch);
                let k_dd = (f_dd + LAMBDA * g_dd).max(0.0);
                for r in 0..3 {
                    h[r][r] += w * a_coef * sum_c2;
                    for q in 0..3 {
                        h[r][q] += w * k_dd * dd[r] * dd[q];
                    }
                }
            }
        }
        if self.fit > 0.0 {
            let p = self.points.get(v as usize).copied().unwrap_or([0.0; 3]);
            if let Some(t) = self.target(v, p) {
                let w = self.fit_weight(v);
                for r in 0..3 {
                    g[r] += 2.0 * w * (p[r] - t[r]);
                    h[r][r] += 2.0 * w;
                }
            }
        }
        (g, h)
    }

    /// Project a point onto vertex `v`'s constraint.
    fn project(&self, v: u32, p: P) -> P {
        self.target(v, p).unwrap_or(p)
    }

    /// The direction a constrained vertex may move along: the step restricted to its tangent.
    fn restrict(&self, v: u32, d: P) -> P {
        let p = self.points.get(v as usize).copied().unwrap_or([0.0; 3]);
        match self.constraint.get(v as usize).copied().unwrap_or(Constraint::Free) {
            Constraint::Free => d,
            Constraint::Corner(_) => [0.0; 3],
            Constraint::Curve(c) => {
                let Some((_, i, _)) = self.geom.nearest_on_curve(c, p) else { return [0.0; 3] };
                let Some(cv) = self.geom.curves.get(c as usize) else { return [0.0; 3] };
                let (a, b) = (cv.vertices.get(i).copied().unwrap_or(0), cv.vertices.get(i + 1).copied().unwrap_or(0));
                let t = sub(self.geom.xyz(b), self.geom.xyz(a));
                let l2 = dot(t, t);
                if l2 > 0.0 { t.map(|x| x * dot(d, t) / l2) } else { [0.0; 3] }
            }
            Constraint::Patch(pt) => {
                let Some((t, _, _)) = self.geom.nearest(p, |t| self.geom.tri_patch.get(t as usize) == Some(&pt)) else { return [0.0; 3] };
                let n = self.geom.normals.get(t as usize).copied().unwrap_or([0.0; 3]);
                let dn = dot(d, n);
                [0, 1, 2].map(|k| d[k] - dn * n[k])
            }
        }
    }

    /// One Newton step of vertex `v`; whether it moved.
    fn relax(&mut self, v: u32) -> bool {
        let soft = self.fit > 0.0;
        if (!soft && matches!(self.constraint.get(v as usize), Some(Constraint::Corner(_))))
            || self.around.get(v as usize).is_none_or(|a| a.is_empty())
        {
            return false;
        }
        let (g, h) = self.derivatives(v);
        let Some(step) = solve3(h, g.map(|x| -x)) else { return false };
        let step = if soft { step } else { self.restrict(v, step) };
        if !step.iter().all(|x| x.is_finite()) || norm(step) == 0.0 {
            return false;
        }
        let old = self.points.get(v as usize).copied().unwrap_or([0.0; 3]);
        let e0 = self.local_energy(v);
        let mut tau = 1.0;
        for _ in 0..12 {
            let moved = [0, 1, 2].map(|k| old[k] + tau * step[k]);
            let trial = if soft { moved } else { self.project(v, moved) };
            if let Some(p) = self.points.get_mut(v as usize) {
                *p = trial;
            }
            let e1 = self.local_energy(v);
            if e1 < e0 && e1.is_finite() {
                return true;
            }
            tau *= 0.5;
        }
        if let Some(p) = self.points.get_mut(v as usize) {
            *p = old;
        }
        false
    }

    /// Set the regularisation from the smallest corner determinant (Garanzha et al. 2021: ε
    /// shrinks as the mesh untangles); that determinant and the energy at the new ε.
    fn update_eps(&mut self) -> (f64, f64) {
        let (min_d, energy) = self.measure();
        let old = self.eps;
        self.eps = (MIN_EPS2 + 0.04 * min_d.min(0.0).powi(2)).sqrt();
        (min_d, if self.eps == old { energy } else { self.measure().1 })
    }

    /// Sweep all vertices until every corner is positive and the energy settles, or `max_sweeps`;
    /// the sweeps run.
    pub fn run(&mut self, max_sweeps: usize, tolerance: f64) -> usize {
        let n = self.points.len() as u32;
        let (_, mut energy) = self.update_eps();
        let mut sweeps = 0;
        for sweep in 0..max_sweeps {
            let mut moved = 0usize;
            for v in 0..n {
                if self.relax(v) {
                    moved += 1;
                }
            }
            let old_eps = self.eps;
            let (min_d, e) = self.update_eps();
            sweeps = sweep + 1;
            let settled = old_eps == self.eps && (energy - e).abs() <= tolerance * energy.abs();
            energy = e;
            if moved == 0 || (min_d > 0.0 && settled) {
                break;
            }
        }
        sweeps
    }
}

impl Optimizer<'_> {
    /// The smallest corner scaled Jacobian of the elements around vertex `v`.
    fn local_quality(&self, v: u32) -> f64 {
        self.around
            .get(v as usize)
            .map_or(f64::INFINITY, |a| a.iter().map(|&(e, _)| min_scaled_jacobian(&self.element(e).0)).fold(f64::INFINITY, f64::min))
    }

    /// Move every constrained vertex exactly onto its feature where its elements keep a corner
    /// scaled Jacobian of at least `floor`; the others stay where the fit left them. Returns
    /// how many stayed off and the largest distance of a constrained vertex from its feature.
    pub fn snap(&mut self, floor: f64) -> (usize, f64) {
        let (mut off, mut worst) = (0usize, 0.0f64);
        for v in 0..self.points.len() as u32 {
            let p = self.points.get(v as usize).copied().unwrap_or([0.0; 3]);
            let Some(t) = self.target(v, p) else { continue };
            if p == t {
                continue;
            }
            if let Some(x) = self.points.get_mut(v as usize) {
                *x = t;
            }
            if self.local_quality(v) < floor {
                if let Some(x) = self.points.get_mut(v as usize) {
                    *x = p;
                }
                off += 1;
                worst = worst.max(dist2(p, t).sqrt());
            }
        }
        (off, worst)
    }

    /// Raise the worst elements: every vertex of an element whose corner scaled Jacobian is
    /// below `threshold` moves, within its constraint, by a compass search that maximises the
    /// smallest scaled Jacobian around it (steps shrink by halves; at most `passes` rounds over
    /// the worst elements, 12 step sizes per vertex).
    pub fn improve_worst(&mut self, threshold: f64, passes: usize) {
        for _ in 0..passes {
            let mut verts: Vec<u32> = Vec::new();
            for (e, h) in self.hexes.iter().enumerate() {
                if min_scaled_jacobian(&self.element(e as u32).0) < threshold {
                    verts.extend(h.iter().copied());
                }
            }
            verts.sort_unstable();
            verts.dedup();
            if verts.is_empty() {
                break;
            }
            let mut moved = false;
            for &v in &verts {
                if matches!(self.constraint.get(v as usize), Some(Constraint::Free)) {
                    moved |= self.compass(v, &|o: &Self, v| o.local_quality(v));
                }
            }
            if !moved {
                break;
            }
        }
    }

    /// Compass search on vertex `v` (within its constraint) for the largest `quality`; steps
    /// start at 0.15 of its largest element size and halve when no direction improves, 12 times
    /// at most. Whether it moved.
    fn compass(&mut self, v: u32, quality: &dyn Fn(&Self, u32) -> f64) -> bool {
        let p0 = self.points.get(v as usize).copied().unwrap_or([0.0; 3]);
        let reach = self
            .around
            .get(v as usize)
            .map_or(0.0, |a| a.iter().map(|&(e, _)| self.size.get(e as usize).copied().unwrap_or(0.0)).fold(0.0, f64::max));
        let mut step = 0.15 * reach;
        let mut best = quality(self, v);
        let mut cur = p0;
        let mut moved = false;
        for _ in 0..12 {
            if step <= MIN_STEP * reach.max(1e-300) {
                break;
            }
            let mut improved = false;
            for d in [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]] {
                let dir = self.restrict(v, d.map(|x| x * step));
                if norm(dir) == 0.0 {
                    continue;
                }
                let trial = self.project(v, [0, 1, 2].map(|k| cur[k] + dir[k]));
                if let Some(p) = self.points.get_mut(v as usize) {
                    *p = trial;
                }
                let q = quality(self, v);
                if q > best + SAME_QUALITY {
                    best = q;
                    cur = trial;
                    improved = true;
                    moved = true;
                } else if let Some(p) = self.points.get_mut(v as usize) {
                    *p = cur;
                }
            }
            if !improved {
                step *= 0.5;
            }
        }
        if let Some(p) = self.points.get_mut(v as usize) {
            *p = cur;
        }
        moved
    }

    /// The smallest certified Jacobian ratio of the elements around `v` (negative where the
    /// trilinear map folds somewhere, even with positive corners).
    fn certified_quality(&self, v: u32) -> f64 {
        self.around.get(v as usize).map_or(f64::INFINITY, |a| {
            a.iter()
                .map(|&(e, _)| {
                    let x = self.element(e).0;
                    if min_scaled_jacobian(&x) <= 0.0 {
                        return WORST;
                    }
                    element_validity(ElementKind::Hex8, &x).map_or(-1.0, |c| match c.validity {
                        Validity::Valid => c.ratio.max(1e-300),
                        _ => c.ratio.min(0.0) - SAME_QUALITY,
                    })
                })
                .fold(f64::INFINITY, f64::min)
        })
    }

    /// Repair elements whose trilinear map is not proven positive although their corners are:
    /// the vertices of each such element (free ones, and boundary ones along their features)
    /// move by compass search to raise the smallest certified bound around them; `rounds` rounds
    /// at most. Returns the elements still not proven valid.
    pub fn repair_certified(&mut self, rounds: usize) -> Vec<u32> {
        let check = |o: &Self, e: u32| element_validity(ElementKind::Hex8, &o.element(e).0).is_ok_and(|c| c.validity == Validity::Valid);
        let mut bad: Vec<u32> = (0..self.hexes.len() as u32).filter(|&e| !check(self, e)).collect();
        for _ in 0..rounds {
            if bad.is_empty() {
                break;
            }
            let mut verts: Vec<u32> = bad.iter().flat_map(|&e| self.hexes.get(e as usize).copied().unwrap_or([0; 8])).collect();
            verts.sort_unstable();
            verts.dedup();
            for v in verts {
                if !matches!(self.constraint.get(v as usize), Some(Constraint::Corner(_)) | None) {
                    self.compass(v, &|o: &Self, v| o.certified_quality(v));
                }
            }
            let mut around: Vec<u32> = bad
                .iter()
                .flat_map(|&e| self.hexes.get(e as usize).copied().unwrap_or([0; 8]))
                .flat_map(|v| self.around.get(v as usize).map_or(Vec::new(), |a| a.iter().map(|x| x.0).collect()))
                .collect();
            around.sort_unstable();
            around.dedup();
            bad = around.into_iter().filter(|&e| !check(self, e)).collect();
        }
        bad
    }
}

/// Solve the symmetric 3 × 3 system `h x = b` (Cramer), `None` when singular.
fn solve3(h: [[f64; 3]; 3], b: P) -> Option<P> {
    let d = dot(h[0], cross(h[1], h[2]));
    if !(d.abs() > 1e-300) || !d.is_finite() {
        return None;
    }
    let col = |k: usize| [h[0][k], h[1][k], h[2][k]];
    let (c0, c1, c2) = (col(0), col(1), col(2));
    Some([dot(b, cross(c1, c2)) / d, dot(c0, cross(b, c2)) / d, dot(c0, cross(c1, b)) / d])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::octree::geom::tests::box_surface;

    #[test]
    fn corner_frames_are_right_handed_on_the_unit_cube() {
        let cube: [P; 8] = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
        for c in 0..8 {
            let j = corner_jacobian(&cube, c, 1.0);
            assert!((det3(&j) - 1.0).abs() < 1e-15, "corner {c}");
        }
        assert!((min_scaled_jacobian(&cube) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn a_tangled_interior_vertex_is_untangled() {
        // A 2 × 2 × 2 block of cubes with the centre vertex pushed through a face; the boundary
        // fixed.
        let s = box_surface([0.0; 3], [2.0; 3]);
        let g = Geometry::new(&s, 30.0, 160.0).unwrap();
        let id = |x: u32, y: u32, z: u32| x + 3 * (y + 3 * z);
        let mut p = Vec::new();
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    p.push([x as f64, y as f64, z as f64]);
                }
            }
        }
        let mut h = Vec::new();
        for z in 0..2 {
            for y in 0..2 {
                for x in 0..2 {
                    h.push([
                        id(x, y, z),
                        id(x + 1, y, z),
                        id(x + 1, y + 1, z),
                        id(x, y + 1, z),
                        id(x, y, z + 1),
                        id(x + 1, y, z + 1),
                        id(x + 1, y + 1, z + 1),
                        id(x, y + 1, z + 1),
                    ]);
                }
            }
        }
        let centre = id(1, 1, 1) as usize;
        p[centre] = [1.7, 1.9, 0.2];
        let size = vec![1.0; h.len()];
        // Boundary vertices slide on the box's faces (hard constraints, patch of the nearest
        // face), the centre is free.
        let cons: Vec<Constraint> = (0..p.len())
            .map(|i| {
                if i == centre {
                    return Constraint::Free;
                }
                let t = g.nearest(p[i], |_| true).unwrap().0;
                Constraint::Patch(g.tri_patch[t as usize])
            })
            .collect();
        let mut o = Optimizer::new(&g, p, &h, &size, cons);
        assert!(o.measure().0 < 0.0);
        let st = o.run(200, 1e-9);
        assert!(o.measure().0 > 0.9, "{st:?}");
        assert!(sub(o.points[centre], [1.0; 3]).iter().all(|c| c.abs() < 1e-3), "{:?}", o.points[centre]);
    }
}
