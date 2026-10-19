//! Size fields: the target element edge length at a point, and the edge divisions that follow it.
//!
//! One [`SizeField`] drives every mesher: the quad meshers ([`crate::quad`]) read it for interior
//! density, and the edge-division helpers here turn it into boundary division counts so that
//! faces sharing an edge share its nodes. The tetrahedral refinement reads the same trait.
//!
//! # Divisions
//! The number of divisions of a curve `C` under a size field `h` is the length of `C` measured in
//! units of `h`, `∫_C ds / h(s)` ([`ideal_divisions`]), rounded ([`edge_divisions`]). The nodes
//! of a divided curve split it into pieces of equal measured length ([`divide_polyline`]), so the
//! spacing grades with the field.
//!
//! # Parity
//! An all-quadrilateral mesh of a region needs an even number of boundary edges in total (every
//! quad has four sides, every interior edge is shared by two quads). The free quad mesher
//! ([`crate::quad::free_quad_mesh`]) needs an even count on *every* loop. [`loop_divisions`]
//! rounds the ideal counts of a loop's curves and fixes the parity of the loop with the change
//! that costs the least relative spacing error; [`even_divisions`] makes a single curve even,
//! which also keeps every loop corner at an even position (the free mesher keeps corners exactly
//! when they are there).
//!
//! The integrals use the midpoint rule on [`SAMPLES_PER_SEGMENT`] pieces per polyline segment,
//! with the size taken constant within a piece.
//! Every function is deterministic (no transcendental functions, only `+ − × ÷ √`).

use crate::{Error, Result};

/// The target edge length at a point. Implementations return a finite, positive length; callers
/// treat anything else as invalid input.
pub trait SizeField {
    fn size(&self, p: [f64; 3]) -> f64;

    /// The smallest and largest size anywhere, when the field knows them (a mesher that builds
    /// from a root size, such as the octree, starts from the largest and plans its depth from
    /// the smallest).
    fn range(&self) -> Option<(f64, f64)> {
        None
    }
}

impl<F: Fn([f64; 3]) -> f64> SizeField for F {
    fn size(&self, p: [f64; 3]) -> f64 {
        self(p)
    }
}

/// The same size everywhere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformSize(pub f64);

impl SizeField for UniformSize {
    fn size(&self, _: [f64; 3]) -> f64 {
        self.0
    }

    fn range(&self) -> Option<(f64, f64)> {
        Some((self.0, self.0))
    }
}

/// The largest division count of one curve.
pub const MAX_DIVISIONS: usize = 1 << 20;
/// The most points a polyline handed to the helpers may have.
pub const MAX_POLYLINE_POINTS: usize = 1 << 20;
/// Midpoint-rule samples per polyline segment.
pub const SAMPLES_PER_SEGMENT: usize = 8;

/// The size at `p`, or an error when the field returns something that is not a finite, positive
/// length.
pub fn checked_size(field: &dyn SizeField, p: [f64; 3]) -> Result<f64> {
    let h = field.size(p);
    if h.is_finite() && h > 0.0 {
        Ok(h)
    } else {
        Err(Error::invalid(crate::message![
            "size field returned ",
            crate::message::length(h),
            " at ",
            crate::message::point(p),
            "; sizes must be finite and positive"
        ]))
    }
}

fn check_polyline(polyline: &[[f64; 3]]) -> Result<()> {
    if polyline.len() < 2 {
        return Err(Error::Invalid(format!("a polyline needs at least 2 points, got {}", polyline.len())));
    }
    if polyline.len() > MAX_POLYLINE_POINTS {
        return Err(Error::Limit(format!("polyline has {} points; the limit is {MAX_POLYLINE_POINTS}", polyline.len())));
    }
    fieldcraft_model::check_coordinates(polyline.as_flattened())?;
    Ok(())
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The polyline cut into [`SAMPLES_PER_SEGMENT`] pieces per segment, each with its measured
/// length `∫ ds / h` (midpoint rule): `(start, end, measure)`.
fn pieces(field: &dyn SizeField, polyline: &[[f64; 3]]) -> Result<Vec<([f64; 3], [f64; 3], f64)>> {
    check_polyline(polyline)?;
    let k = SAMPLES_PER_SEGMENT as f64;
    let mut out = Vec::with_capacity(SAMPLES_PER_SEGMENT * polyline.len());
    for w in polyline.windows(2) {
        let [a, b] = w else { continue };
        let piece = dist(*a, *b) / k;
        for s in 0..SAMPLES_PER_SEGMENT {
            let (p, q) = (lerp(*a, *b, s as f64 / k), lerp(*a, *b, (s + 1) as f64 / k));
            let m = piece / checked_size(field, lerp(*a, *b, (s as f64 + 0.5) / k))?;
            if !m.is_finite() {
                return Err(Error::Invalid("size field is too small for this polyline".into()));
            }
            out.push((p, q, m));
        }
    }
    Ok(out)
}

/// The ideal (real) number of divisions of a polyline: its length measured in units of the size
/// field.
pub fn ideal_divisions(field: &dyn SizeField, polyline: &[[f64; 3]]) -> Result<f64> {
    let total: f64 = pieces(field, polyline)?.iter().map(|p| p.2).sum();
    if total.is_finite() { Ok(total) } else { Err(Error::Invalid("size field is too small for this polyline".into())) }
}

fn rounded(ideal: f64) -> Result<usize> {
    if !(ideal >= 0.0) || !ideal.is_finite() {
        return Err(Error::Invalid(format!("ideal division count {ideal} is not a finite, non-negative number")));
    }
    if ideal > MAX_DIVISIONS as f64 {
        return Err(Error::Limit(format!("{ideal:.0} divisions exceed the limit of {MAX_DIVISIONS}")));
    }
    Ok((ideal.round() as usize).max(1))
}

/// The division count of one polyline: its ideal count rounded, at least 1.
pub fn edge_divisions(field: &dyn SizeField, polyline: &[[f64; 3]]) -> Result<usize> {
    rounded(ideal_divisions(field, polyline)?)
}

/// The nearest even count to `ideal`, at least 2 (ties go up).
pub fn even_divisions(ideal: f64) -> Result<usize> {
    let half = rounded(ideal / 2.0)?;
    half.checked_mul(2).ok_or_else(|| Error::Limit("division count overflow".into()))
}

/// The relative spacing error of `n` divisions against `ideal`, as `max(n, x) / min(n, x)`.
fn spacing_error(n: usize, ideal: f64) -> f64 {
    let (n, x) = (n as f64, ideal.max(f64::MIN_POSITIVE));
    if n > x { n / x } else { x / n }
}

/// Division counts for the curves of one loop with an even total: each ideal count is rounded
/// (at least 1); when the sum is odd, the one curve whose count changes by ±1 with the smallest
/// increase of its spacing error is changed (ties: the lower index, then +1 before −1).
pub fn loop_divisions(ideal: &[f64]) -> Result<Vec<usize>> {
    if ideal.is_empty() {
        return Err(Error::Invalid("a loop needs at least one curve".into()));
    }
    let mut counts = ideal.iter().map(|&x| rounded(x)).collect::<Result<Vec<_>>>()?;
    let total = counts.iter().try_fold(0usize, |s, &n| s.checked_add(n)).ok_or_else(|| Error::Limit("division count overflow".into()))?;
    if total % 2 == 0 {
        return Ok(counts);
    }
    let mut best: Option<(f64, usize, bool)> = None;
    for (i, (&n, &x)) in counts.iter().zip(ideal).enumerate() {
        let now = spacing_error(n, x);
        for up in [true, false] {
            if !up && n <= 1 {
                continue;
            }
            let m = if up { n + 1 } else { n - 1 };
            let cost = spacing_error(m, x) - now;
            if best.is_none_or(|(c, _, _)| cost < c) {
                best = Some((cost, i, up));
            }
        }
    }
    if let Some((_, i, up)) = best
        && let Some(n) = counts.get_mut(i)
    {
        *n = if up { *n + 1 } else { *n - 1 };
    }
    Ok(counts)
}

/// `n + 1` points along a polyline (the first and last are its end points, copied exactly) that
/// split it into `n` pieces of equal length measured by the size field.
pub fn divide_polyline(field: &dyn SizeField, polyline: &[[f64; 3]], n: usize) -> Result<Vec<[f64; 3]>> {
    if n == 0 || n > MAX_DIVISIONS {
        return Err(Error::Invalid(format!("division count {n} is outside 1..={MAX_DIVISIONS}")));
    }
    let pieces = pieces(field, polyline)?;
    let total: f64 = pieces.iter().map(|p| p.2).sum();
    let (first, last) = match (polyline.first(), polyline.last()) {
        (Some(a), Some(b)) => (*a, *b),
        _ => return Err(Error::Invalid("empty polyline".into())),
    };
    if !(total > 0.0) || !total.is_finite() {
        return Err(Error::Invalid("polyline has zero length".into()));
    }
    let mut out = Vec::with_capacity(n + 1);
    out.push(first);
    let (mut piece, mut before) = (0usize, 0.0);
    for k in 1..n {
        let target = total * k as f64 / n as f64;
        while let Some(&(_, _, m)) = pieces.get(piece) {
            if before + m >= target || piece + 1 == pieces.len() {
                break;
            }
            before += m;
            piece += 1;
        }
        let (a, b, m) = pieces.get(piece).copied().ok_or_else(|| Error::Failed("polyline division walked off the end".into()))?;
        // The size is taken constant within a piece.
        let t = if m > 0.0 { ((target - before) / m).clamp(0.0, 1.0) } else { 0.0 };
        out.push(lerp(a, b, t));
    }
    out.push(last);
    Ok(out)
}

/// The most triangles the local sizes of one field hold.
pub const MAX_REGION_TRIANGLES: usize = 1 << 20;

/// Where a [`LocalSize`] applies.
#[derive(Debug, Clone, PartialEq)]
pub enum Region {
    /// An axis-aligned box.
    Box { lo: [f64; 3], hi: [f64; 3] },
    /// A ball.
    Sphere { center: [f64; 3], radius: f64 },
    /// A surface, as triangles (a B-rep face's tessellation).
    Triangles(Vec<[[f64; 3]; 3]>),
}

impl Region {
    /// The distance from `p` to the region (0 inside a box or a ball, and on a surface).
    pub fn distance(&self, p: [f64; 3]) -> f64 {
        match self {
            Region::Box { lo, hi } => {
                let d = [0, 1, 2].map(|k| (lo[k] - p[k]).max(p[k] - hi[k]).max(0.0));
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            }
            Region::Sphere { center, radius } => (dist(p, *center) - radius).max(0.0),
            Region::Triangles(tris) => {
                tris.iter().map(|[a, b, c]| dist(p, crate::tet::geometry::closest_on_triangle(p, *a, *b, *c))).fold(f64::INFINITY, f64::min)
            }
        }
    }

    /// The region's bounding box.
    pub fn bounds(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            Region::Box { lo, hi } => (*lo, *hi),
            Region::Sphere { center, radius } => (center.map(|c| c - radius), center.map(|c| c + radius)),
            Region::Triangles(tris) => {
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for p in tris.iter().flatten() {
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k]);
                        hi[k] = hi[k].max(p[k]);
                    }
                }
                (lo, hi)
            }
        }
    }
}

/// A local size: elements of about `size` in `region`.
#[derive(Debug, Clone, PartialEq)]
pub struct LocalSize {
    pub region: Region,
    pub size: f64,
}

/// A global size refined by local sizes and graded: the size at `p` is
/// `min(global, min_i (h_i + g · d_i(p)))` over the local sizes `h_i` with `d_i` the distance to
/// their regions, so it changes by at most `g` per unit length. Exact (no background grid):
/// meant for the structured meshers, which query it along curves; the tetrahedral refinement
/// bakes the same local sizes into its graded grid ([`crate::tet::geometry::GradedSize`]).
#[derive(Debug, Clone, PartialEq)]
pub struct RegionSize {
    global: f64,
    growth: f64,
    /// Each local size with its region's bounds grown by the reach of its refinement.
    locals: Vec<(LocalSize, [f64; 3], [f64; 3])>,
}

impl RegionSize {
    /// A checked field: a finite, positive `global` and local sizes, `growth` finite and
    /// positive, regions with finite coordinates and at most [`MAX_REGION_TRIANGLES`] triangles
    /// in all. Local sizes not below `global` change nothing and are dropped.
    pub fn new(global: f64, growth: f64, locals: Vec<LocalSize>) -> Result<RegionSize> {
        let positive = |what: &str, v: f64| {
            if v > 0.0 && v.is_finite() { Ok(v) } else { Err(Error::Invalid(format!("{what} must be finite and positive, got {v}"))) }
        };
        positive("the global size", global)?;
        positive("the growth rate", growth)?;
        let triangles: usize = locals.iter().map(|l| if let Region::Triangles(t) = &l.region { t.len() } else { 0 }).sum();
        if triangles > MAX_REGION_TRIANGLES {
            return Err(Error::Limit(format!("size regions hold {triangles} triangles; the limit is {MAX_REGION_TRIANGLES}")));
        }
        let mut out = Vec::with_capacity(locals.len());
        for l in locals {
            positive("a local size", l.size)?;
            match &l.region {
                Region::Box { lo, hi } => fieldcraft_model::check_coordinates(&[lo.as_slice(), hi.as_slice()].concat())?,
                Region::Sphere { center, radius } => {
                    fieldcraft_model::check_coordinates(center)?;
                    positive("a sphere's radius", *radius)?;
                }
                Region::Triangles(t) => fieldcraft_model::check_coordinates(t.as_flattened().as_flattened())?,
            }
            if l.size >= global {
                continue;
            }
            let (lo, hi) = l.region.bounds();
            if (0..3).any(|k| !(lo[k] <= hi[k])) {
                continue; // no triangles
            }
            let reach = (global - l.size) / growth;
            out.push((l, lo.map(|c| c - reach), hi.map(|c| c + reach)));
        }
        Ok(RegionSize { global, growth, locals: out })
    }

    pub fn global(&self) -> f64 {
        self.global
    }

    pub fn growth(&self) -> f64 {
        self.growth
    }

    /// Whether any local size refines the global one.
    pub fn is_uniform(&self) -> bool {
        self.locals.is_empty()
    }

    /// The local sizes that refine the global one.
    pub fn locals(&self) -> impl Iterator<Item = &LocalSize> {
        self.locals.iter().map(|l| &l.0)
    }

    /// The smallest size anywhere.
    pub fn finest(&self) -> f64 {
        self.locals().map(|l| l.size).fold(self.global, f64::min)
    }
}

impl SizeField for RegionSize {
    fn size(&self, p: [f64; 3]) -> f64 {
        let mut h = self.global;
        for (l, lo, hi) in &self.locals {
            if (0..3).any(|k| p[k] < lo[k] || p[k] > hi[k]) {
                continue;
            }
            h = h.min(l.size + self.growth * l.region.distance(p));
        }
        h
    }

    fn range(&self) -> Option<(f64, f64)> {
        Some((self.finest(), self.global))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_divisions() {
        let line = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        assert_eq!(edge_divisions(&UniformSize(1.0), &line).unwrap(), 10);
        assert_eq!(edge_divisions(&UniformSize(3.0), &line).unwrap(), 3);
        assert_eq!(edge_divisions(&UniformSize(100.0), &line).unwrap(), 1);
        let pts = divide_polyline(&UniformSize(1.0), &line, 4).unwrap();
        assert_eq!(pts, vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [5.0, 0.0, 0.0], [7.5, 0.0, 0.0], [10.0, 0.0, 0.0]]);
    }

    #[test]
    fn graded_division_is_finer_where_the_size_is_smaller() {
        let field = |p: [f64; 3]| 0.5 + 0.1 * p[0];
        let line = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let pts = divide_polyline(&field, &line, 6).unwrap();
        let gaps: Vec<f64> = pts.windows(2).map(|w| w[1][0] - w[0][0]).collect();
        assert!(gaps.windows(2).all(|g| g[1] > g[0]), "{gaps:?}");
        assert_eq!(pts.last(), Some(&[10.0, 0.0, 0.0]));
    }

    #[test]
    fn parity() {
        assert_eq!(even_divisions(3.2).unwrap(), 4);
        assert_eq!(even_divisions(0.1).unwrap(), 2);
        assert_eq!(even_divisions(5.0).unwrap(), 6);
        let c = loop_divisions(&[3.0, 2.0, 3.0, 2.0]).unwrap();
        assert_eq!(c, vec![3, 2, 3, 2]);
        let c = loop_divisions(&[3.0, 2.0, 3.0, 2.6]).unwrap();
        assert_eq!(c.iter().sum::<usize>() % 2, 0);
        // The curve whose rounding was the worst absorbs the fix.
        assert_eq!(c, vec![3, 2, 3, 2]);
        let c = loop_divisions(&[1.0, 1.0, 1.0]).unwrap();
        assert_eq!(c, vec![2, 1, 1]);
    }

    #[test]
    fn hostile() {
        let line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert!(edge_divisions(&UniformSize(0.0), &line).is_err());
        assert!(edge_divisions(&UniformSize(f64::NAN), &line).is_err());
        assert!(edge_divisions(&UniformSize(-1.0), &line).is_err());
        assert!(edge_divisions(&UniformSize(1e-300), &line).is_err());
        assert!(edge_divisions(&UniformSize(1.0), &[[0.0; 3]]).is_err());
        assert!(edge_divisions(&UniformSize(1.0), &[[0.0; 3], [f64::NAN, 0.0, 0.0]]).is_err());
        assert!(divide_polyline(&UniformSize(1.0), &[[0.0; 3], [0.0; 3]], 2).is_err());
        assert!(divide_polyline(&UniformSize(1.0), &line, 0).is_err());
        assert!(loop_divisions(&[]).is_err());
        assert!(loop_divisions(&[f64::INFINITY]).is_err());
        assert!(even_divisions(f64::NAN).is_err());
    }

    #[test]
    fn region_size_refines_inside_and_grows_away() {
        let boxed = LocalSize { region: Region::Box { lo: [0.0; 3], hi: [2.0; 3] }, size: 0.5 };
        let ball = LocalSize { region: Region::Sphere { center: [10.0, 0.0, 0.0], radius: 1.0 }, size: 0.25 };
        let face = LocalSize { region: Region::Triangles(vec![[[0.0, 10.0, 0.0], [1.0, 10.0, 0.0], [0.0, 11.0, 0.0]]]), size: 1.0 };
        let coarse = LocalSize { region: Region::Box { lo: [0.0; 3], hi: [1.0; 3] }, size: 5.0 };
        let f = RegionSize::new(2.0, 0.5, vec![boxed, ball, face, coarse]).unwrap();
        assert_eq!(f.locals().count(), 3, "the coarser control is dropped");
        assert_eq!(f.size([1.0, 1.0, 1.0]), 0.5);
        assert_eq!(f.size([3.0, 1.0, 1.0]), 1.0, "one unit out: 0.5 + 0.5");
        assert_eq!(f.size([20.0, 20.0, 20.0]), 2.0);
        assert_eq!(f.size([10.5, 0.0, 0.0]), 0.25);
        assert_eq!(f.size([0.2, 10.2, 0.0]), 1.0);
        assert!((f.size([0.2, 10.2, 1.0]) - 1.5).abs() < 1e-12);
        assert_eq!(f.finest(), 0.25);
        // Graded: the size never changes faster than the growth rate.
        let pts: Vec<[f64; 3]> = (0..400).map(|i| [-5.0 + 0.05 * i as f64, 1.0, 1.0]).collect();
        for w in pts.windows(2) {
            assert!((f.size(w[1]) - f.size(w[0])).abs() <= 0.5 * 0.05 + 1e-12);
        }
        // Divisions follow it: a line through the box is cut finer inside.
        let line = [[-4.0, 1.0, 1.0], [6.0, 1.0, 1.0]];
        let n = edge_divisions(&f, &line).unwrap();
        let pts = divide_polyline(&f, &line, n).unwrap();
        let inside = pts.windows(2).filter(|w| w[0][0] >= 0.0 && w[1][0] <= 2.0).map(|w| w[1][0] - w[0][0]).fold(0.0, f64::max);
        assert!(inside < 0.6, "{inside}");
    }

    #[test]
    fn hostile_region_sizes_are_refused() {
        let ok = || LocalSize { region: Region::Sphere { center: [0.0; 3], radius: 1.0 }, size: 0.5 };
        assert!(RegionSize::new(0.0, 0.3, vec![]).is_err());
        assert!(RegionSize::new(1.0, f64::NAN, vec![]).is_err());
        assert!(RegionSize::new(1.0, 0.3, vec![LocalSize { size: -1.0, ..ok() }]).is_err());
        assert!(RegionSize::new(1.0, 0.3, vec![LocalSize { region: Region::Sphere { center: [f64::NAN; 3], radius: 1.0 }, size: 0.5 }]).is_err());
        assert!(RegionSize::new(1.0, 0.3, vec![LocalSize { region: Region::Triangles(vec![[[f64::INFINITY; 3]; 3]]), size: 0.5 }]).is_err());
        let empty = RegionSize::new(1.0, 0.3, vec![LocalSize { region: Region::Triangles(vec![]), size: 0.5 }]).unwrap();
        assert!(empty.is_uniform());
        assert!(RegionSize::new(1.0, 0.3, vec![ok()]).unwrap().size([f64::NAN; 3]) > 0.0);
    }
}
