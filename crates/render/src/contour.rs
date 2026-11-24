//! Banded contour plots: a range of nodal values cut into equal bands, each drawn in one colour.
//!
//! A triangle with values at its corners is split exactly along the band boundaries (the
//! values vary linearly across it), so every piece lies in one band and the band edges are
//! sharp iso-lines on any renderer, with no shader support: the GPU viewport and the CPU
//! rasterizer draw the same flat-coloured pieces.

use crate::{Rgb, SceneTri, Vec3};

/// Equal bands from `min` to `max`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bands {
    pub min: f64,
    pub max: f64,
    pub count: u32,
    /// Values above `max` go to an extra band, [`Bands::count`] (drawn in the "over" colour), not
    /// to the top band: the range was clamped below the field's maximum.
    pub over: bool,
}

impl Bands {
    /// `count` bands over `[min, max]`; a range of one value (or none) is one band around it.
    pub fn new(min: f64, max: f64, count: u32) -> Bands {
        let count = count.max(1);
        if !(min.is_finite() && max.is_finite()) {
            return Bands { min: 0.0, max: 0.0, count, over: false };
        }
        Bands { min: min.min(max), max: max.max(min), count, over: false }
    }

    /// The bands of a field ranging over `[min, max]` drawn over `[lo, hi]` where given (a
    /// clamped legend): values below `lo` go to the bottom band, values above `hi` to the "over"
    /// band. A clamp that leaves no range (lo ≥ hi) is ignored.
    pub fn clamped(min: f64, max: f64, lo: Option<f64>, hi: Option<f64>, count: u32) -> Bands {
        let (a, b) = (lo.filter(|v| v.is_finite()).unwrap_or(min), hi.filter(|v| v.is_finite()).unwrap_or(max));
        if !(a < b) {
            return Bands::new(min, max, count);
        }
        Bands { over: max > b, ..Bands::new(a, b, count) }
    }

    /// Whether band `k` is the "over" band (values above a clamped maximum).
    pub fn is_over(&self, k: u32) -> bool {
        self.over && k >= self.count
    }

    /// The band a value lies in (values outside the range go to the end bands); `None` for a
    /// value that is not a number.
    pub fn band(&self, v: f64) -> Option<u32> {
        if v.is_nan() {
            return None;
        }
        if self.over && v > self.max {
            return Some(self.count);
        }
        if self.max <= self.min {
            return Some(self.count / 2);
        }
        if v <= self.min {
            return Some(0);
        }
        Some(self.crossed_cuts(v, true).min(self.count.saturating_sub(1)))
    }

    /// The first and last bands occupied by the interpolated interior of a line or triangle.
    /// A cut at the highest sample touches the boundary only. A constant field keeps its
    /// point color, including at a clamped maximum.
    pub fn interior_bands(&self, values: &[f64]) -> Option<(u32, u32)> {
        let &first = values.first()?;
        if values.iter().any(|v| v.is_nan()) {
            return None;
        }
        let (low, high) = values.iter().fold((first, first), |(low, high), &v| (low.min(v), high.max(v)));
        if low == high || self.min == self.max {
            let band = self.band(low)?;
            return Some((band, band));
        }
        Some((self.crossed_cuts(low, true), self.crossed_cuts(high, false)))
    }

    fn crossed_cuts(&self, value: f64, inclusive: bool) -> u32 {
        let mut low = 0;
        let mut high = if self.over { self.count } else { self.count.saturating_sub(1) };
        while low < high {
            let k = low + (high - low).div_ceil(2);
            let cut = self.level(k);
            if cut < value || (inclusive && cut == value) {
                low = k;
            } else {
                high = k - 1;
            }
        }
        low
    }

    /// The value between band `k − 1` and band `k` (`level(0)` is `min`, `level(count)` `max`).
    pub fn level(&self, k: u32) -> f64 {
        if k == 0 {
            return self.min;
        }
        if k == self.count {
            return self.max;
        }
        let scaled = (self.max - self.min) * f64::from(k);
        if scaled.is_finite() {
            self.min + scaled / f64::from(self.count)
        } else {
            let fraction = f64::from(k) / f64::from(self.count);
            self.min * (1.0 - fraction) + self.max * fraction
        }
    }
}

/// Split a triangle (corners `p`, values `v`) into pieces that each lie in one band, and hand
/// each piece (a convex polygon of 3 to 5 corners) to `out` with its band. A triangle with a
/// value that is not a number is handed whole, with no band.
pub fn split_triangle(p: [Vec3; 3], v: [f64; 3], bands: &Bands, out: &mut impl FnMut(&[Vec3], Option<u32>)) {
    let Some((lo, hi)) = bands.interior_bands(&v) else {
        out(&p, None);
        return;
    };
    if lo == hi {
        out(&p, Some(lo));
        return;
    }
    let poly: Vec<(Vec3, f64)> = p.into_iter().zip(v).collect();
    for k in lo..=hi {
        // The end bands also hold the values beyond the range, but for the "over" band above a
        // clamped maximum.
        let below = if k == 0 { f64::NEG_INFINITY } else { bands.level(k.min(bands.count)) };
        let above = if k >= bands.count || (k + 1 == bands.count && !bands.over) { f64::INFINITY } else { bands.level(k + 1) };
        if below >= above {
            continue;
        }
        let piece = clip(&clip(&poly, below, true), above, false);
        if piece.len() >= 3 {
            let pts: Vec<Vec3> = piece.iter().map(|(q, _)| *q).collect();
            out(&pts, Some(k));
        }
    }
}

/// [`split_triangle`] with each piece fanned into scene triangles of `normal`, coloured by `color`
/// of its band (`None`: a value that is not a number).
pub fn banded(p: [Vec3; 3], v: [f64; 3], normal: Vec3, bands: &Bands, color: &impl Fn(Option<u32>) -> Rgb, out: &mut Vec<SceneTri>) {
    split_triangle(p, v, bands, &mut |piece, band| {
        let color = color(band);
        for k in 1..piece.len().saturating_sub(1) {
            if let (Some(&a), Some(&b), Some(&c)) = (piece.first(), piece.get(k), piece.get(k + 1)) {
                out.push(SceneTri { p: [a, b, c], normal, color });
            }
        }
    });
}

/// The part of a convex polygon where the value is at least (`keep_above`) or at most `level`
/// (Sutherland–Hodgman against one level).
fn clip(poly: &[(Vec3, f64)], level: f64, keep_above: bool) -> Vec<(Vec3, f64)> {
    if level.is_infinite() {
        return poly.to_vec();
    }
    let inside = |v: f64| if keep_above { v >= level } else { v <= level };
    let mut out = Vec::with_capacity(poly.len() + 1);
    for (i, &(a, va)) in poly.iter().enumerate() {
        let Some(&(b, vb)) = poly.get((i + 1) % poly.len()) else { continue };
        if inside(va) {
            out.push((a, va));
        }
        if va != level && vb != level && inside(va) != inside(vb) {
            let span = vb - va;
            let t = if span.is_finite() { (level - va) / span } else { (level * 0.5 - va * 0.5) / (vb * 0.5 - va * 0.5) };
            out.push((a + (b - a) * t, level));
        }
    }
    out
}

/// The colour of band `k` of `count`: the colour scale `stops` (low to high) sampled at the
/// band's middle.
/// The colour of band `k` of `bands`: [`band_color`], or `over` for the band above a clamped
/// maximum.
pub fn bands_color(stops: &[Rgb], over: Rgb, bands: &Bands, k: u32) -> Rgb {
    if bands.is_over(k) { over } else { band_color(stops, k, bands.count) }
}

pub fn band_color(stops: &[Rgb], k: u32, count: u32) -> Rgb {
    let (Some(first), Some(last)) = (stops.first(), stops.last()) else { return Rgb(128, 128, 128) };
    if stops.len() == 1 {
        return *first;
    }
    let count = count.max(1);
    let t = (f64::from(k.min(count - 1)) + 0.5) / f64::from(count) * (stops.len() - 1) as f64;
    let i = (t.floor() as usize).min(stops.len() - 2);
    let f = t - i as f64;
    let (Some(a), Some(b)) = (stops.get(i), stops.get(i + 1)) else { return *last };
    let mix = |x: u8, y: u8| (f64::from(x) + (f64::from(y) - f64::from(x)) * f).round().clamp(0.0, 255.0) as u8;
    Rgb(mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(p: &[Vec3]) -> f64 {
        let mut n = Vec3::ZERO;
        for i in 1..p.len() - 1 {
            n = n + (p[i] - p[0]).cross(p[i + 1] - p[0]);
        }
        n.len() * 0.5
    }

    #[test]
    fn bands_and_levels() {
        let b = Bands::new(0.0, 10.0, 5);
        assert_eq!(
            (b.band(0.0), b.band(1.99), b.band(2.0), b.band(10.0), b.band(-5.0), b.band(99.0)),
            (Some(0), Some(0), Some(1), Some(4), Some(0), Some(4))
        );
        assert_eq!((b.level(0), b.level(5)), (0.0, 10.0));
        assert_eq!(b.band(f64::NAN), None);
        // One value: one band in the middle; reversed and non-finite ranges are made safe.
        assert_eq!(Bands::new(3.0, 3.0, 6).band(3.0), Some(3));
        assert_eq!(Bands::new(5.0, 1.0, 4), Bands::new(1.0, 5.0, 4));
        assert_eq!(Bands::new(f64::NAN, 1.0, 0).count, 1);
        // Rounded cuts can coincide with the legend endpoints; constant samples keep their endpoint colors.
        let adjacent = Bands::new(1., 1.0_f64.next_up(), 12);
        assert_eq!(adjacent.level(6), adjacent.min);
        assert_eq!(adjacent.interior_bands(&[adjacent.min; 3]), Some((0, 0)));
        assert_eq!(adjacent.interior_bands(&[adjacent.max; 3]), Some((11, 11)));
        assert_eq!(adjacent.interior_bands(&[adjacent.min, adjacent.max]), Some((6, 6)));
    }

    #[test]
    fn a_triangle_splits_along_the_levels_and_keeps_its_area() {
        let p = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), Vec3::new(0.0, 10.0, 0.0)];
        // The value is x: pieces are strips between x = 0, 2, 4, …
        let b = Bands::new(0.0, 10.0, 5);
        let mut pieces = Vec::new();
        split_triangle(p, [0.0, 10.0, 0.0], &b, &mut |q, k| pieces.push((q.to_vec(), k)));
        assert_eq!(pieces.iter().map(|(_, k)| k.unwrap()).collect::<Vec<_>>(), vec![0, 1, 2, 3, 4]);
        let total: f64 = pieces.iter().map(|(q, _)| area(q)).sum();
        assert!((total - 50.0).abs() < 1e-9, "{total}");
        for (q, k) in &pieces {
            let k = k.unwrap() as f64;
            assert!(q.iter().all(|c| c.x >= 2.0 * k - 1e-9 && c.x <= 2.0 * k + 2.0 + 1e-9), "{q:?} in band {k}");
        }
        // A band-uniform triangle stays whole; a NaN value hands it over with no band.
        let mut n = 0;
        split_triangle(p, [0.5, 1.0, 1.5], &b, &mut |q, k| {
            n += 1;
            assert_eq!((q.len(), k), (3, Some(0)));
        });
        split_triangle(p, [0.5, f64::NAN, 1.5], &b, &mut |_, k| {
            n += 1;
            assert_eq!(k, None);
        });
        assert_eq!(n, 2);
        // Fanned into scene triangles, keeping the area and the bands.
        let mut tris = Vec::new();
        banded(p, [0.0, 10.0, 0.0], Vec3::Z, &b, &|k| Rgb(k.map_or(0, |k| k as u8), 0, 0), &mut tris);
        assert!((tris.iter().map(|t| area(&t.p)).sum::<f64>() - 50.0).abs() < 1e-9);
        assert_eq!(tris.iter().map(|t| t.color.0).max(), Some(4));
        assert!(tris.iter().all(|t| t.normal == Vec3::Z));
    }

    #[test]
    fn adjacent_samples_at_a_cut_color_only_the_area_in_each_band() {
        let bands = Bands::new(0., 1., 12);
        let cut = bands.level(5);
        let low = cut.next_down();
        let high = cut.next_up();
        assert_eq!(bands.band(low), Some(4));
        for (values, expected) in [([low, cut, cut], [0.5, 0.]), ([high, cut, cut], [0., 0.5]), ([low, high, high], [0.125, 0.375])] {
            let mut triangles = Vec::new();
            banded([Vec3::ZERO, Vec3::X, Vec3::Y], values, Vec3::Z, &bands, &|k| Rgb(k.unwrap() as u8, 0, 0), &mut triangles);
            let mut areas = [0.; 2];
            for tri in triangles {
                assert!((4..=5).contains(&tri.color.0));
                areas[usize::from(tri.color.0 - 4)] += area(&tri.p);
            }
            assert_eq!(areas, expected);
        }
    }

    /// A clamped range puts the values above it in the "over" band, split exactly at the
    /// clamp; one below it in the bottom band; a clamp leaving no range is ignored.
    #[test]
    fn a_clamped_range_has_an_over_band() {
        let b = Bands::clamped(0.0, 100.0, None, Some(10.0), 5);
        assert_eq!((b.min, b.max, b.over), (0.0, 10.0, true));
        assert_eq!((b.band(9.9), b.band(10.0), b.band(10.1), b.band(1e9)), (Some(4), Some(4), Some(5), Some(5)));
        assert!(b.is_over(5) && !b.is_over(4));
        let p = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(20.0, 0.0, 0.0), Vec3::new(0.0, 20.0, 0.0)];
        let mut pieces = Vec::new();
        split_triangle(p, [0.0, 20.0, 0.0], &b, &mut |q, k| pieces.push((q.to_vec(), k.unwrap())));
        assert_eq!(pieces.iter().map(|(_, k)| *k).collect::<Vec<_>>(), vec![0, 1, 2, 3, 4, 5]);
        let (over, _) = pieces.last().unwrap();
        assert!(over.iter().all(|c| c.x >= 10.0 - 1e-9), "{over:?}");
        assert!((pieces.iter().map(|(q, _)| area(q)).sum::<f64>() - 200.0).abs() < 1e-9);
        let lo = Bands::clamped(0.0, 100.0, Some(50.0), None, 5);
        assert_eq!((lo.min, lo.max, lo.over, lo.band(3.0)), (50.0, 100.0, false, Some(0)));
        assert_eq!(Bands::clamped(0.0, 100.0, Some(50.0), Some(40.0), 5), Bands::new(0.0, 100.0, 5));
        assert!(!Bands::clamped(0.0, 5.0, None, Some(10.0), 5).over, "a clamp above the maximum clamps nothing");
        let stops = [Rgb(0, 0, 255), Rgb(255, 0, 0)];
        assert_eq!(bands_color(&stops, Rgb(1, 2, 3), &b, 5), Rgb(1, 2, 3));
        assert_eq!(bands_color(&stops, Rgb(1, 2, 3), &b, 4), band_color(&stops, 4, 5));
    }

    #[test]
    fn band_colours_run_through_the_stops() {
        let stops = [Rgb(0, 0, 255), Rgb(0, 255, 0), Rgb(255, 0, 0)];
        assert_eq!(band_color(&stops, 0, 2), Rgb(0, 128, 128));
        assert_eq!(band_color(&stops, 1, 2), Rgb(128, 128, 0));
        let first = band_color(&stops, 0, 12);
        let last = band_color(&stops, 11, 12);
        assert!(first.2 > 200 && last.0 > 200, "{first:?} {last:?}");
        assert_eq!(band_color(&[], 0, 3), Rgb(128, 128, 128));
        assert_eq!(band_color(&stops, 99, 3), band_color(&stops, 2, 3));
    }
}
