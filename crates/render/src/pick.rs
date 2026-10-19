//! Picking: what a ray through a pixel hits first, and what lies near a pixel or inside a box on
//! screen.
//!
//! Rays are tested against triangles with the watertight ray–triangle intersection of Woop,
//! Benthin and Wald ("Watertight Ray/Triangle Intersection", Journal of Computer Graphics
//! Techniques 2(1), 2013): the ray is sheared onto the +z axis and the triangle's edges are tested
//! with 2D edge functions. An edge shared by two triangles is evaluated from the same operands in
//! both (once each way round, which IEEE arithmetic negates exactly), so a ray through a shared
//! edge or vertex hits at least one of them and never slips between. Hits are two-sided. The
//! nearest hit wins; equal distances go to the lower index, so the answer is the same on every
//! run and platform.
//!
//! Screen-space picks (nodes, edges, box selection) go through a [`Projector`] (the camera's
//! view-projection, computed once) and, for box selection, a [`DepthMap`] of the drawn faces that
//! tells which points are in front.

use crate::raster::raster_tri;
use crate::{Camera, Mat4, Vec3};

/// The distance along a ray (`origin`, `dir`; `dir` need not be unit) at which it hits triangle
/// `p`, from either side; `None` when it misses, runs in the triangle's plane, or hits behind the
/// origin. Non-finite input misses.
pub fn ray_triangle(origin: Vec3, dir: Vec3, p: &[Vec3; 3]) -> Option<f64> {
    let d = [dir.x, dir.y, dir.z];
    // The dominant axis of the direction becomes z; x and y keep the winding.
    let kz = if d[0].abs() >= d[1].abs() && d[0].abs() >= d[2].abs() {
        0
    } else if d[1].abs() >= d[2].abs() {
        1
    } else {
        2
    };
    let (mut kx, mut ky) = ((kz + 1) % 3, (kz + 2) % 3);
    let dz = *d.get(kz)?;
    if !(dz.abs() > 0.0) || !dz.is_finite() {
        return None;
    }
    if dz < 0.0 {
        std::mem::swap(&mut kx, &mut ky);
    }
    let (sx, sy, sz) = (d.get(kx)? / dz, d.get(ky)? / dz, 1.0 / dz);
    let rel = |q: Vec3| {
        let r = q - origin;
        [r.x, r.y, r.z]
    };
    let (a, b, c) = (rel(p[0]), rel(p[1]), rel(p[2]));
    let shear = |v: &[f64; 3]| -> Option<(f64, f64, f64)> {
        let z = *v.get(kz)?;
        Some((v.get(kx)? - sx * z, v.get(ky)? - sy * z, sz * z))
    };
    let ((ax, ay, az), (bx, by, bz), (cx, cy, cz)) = (shear(&a)?, shear(&b)?, shear(&c)?);
    // Edge functions: e(P, Q) = Px·Qy − Py·Qx for the edges (C, B), (A, C), (B, A).
    let u = cx * by - cy * bx;
    let v = ax * cy - ay * cx;
    let w = bx * ay - by * ax;
    if !(u.is_finite() && v.is_finite() && w.is_finite()) {
        return None;
    }
    if (u < 0.0 || v < 0.0 || w < 0.0) && (u > 0.0 || v > 0.0 || w > 0.0) {
        return None;
    }
    let det = u + v + w;
    if det == 0.0 {
        return None;
    }
    let t = (u * az + v * bz + w * cz) / det;
    (t >= 0.0 && t.is_finite()).then_some(t)
}

/// A ray's first hit among triangles: the triangle's index and the distance along the ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub index: usize,
    pub t: f64,
}

/// The first triangle `tris` (in order) that the ray hits: the nearest, the lowest index among
/// equally near ones.
pub fn first_hit(origin: Vec3, dir: Vec3, tris: impl IntoIterator<Item = [Vec3; 3]>) -> Option<Hit> {
    let mut best: Option<Hit> = None;
    for (index, p) in tris.into_iter().enumerate() {
        if let Some(t) = ray_triangle(origin, dir, &p)
            && best.is_none_or(|b| t < b.t)
        {
            best = Some(Hit { index, t });
        }
    }
    best
}

/// The distance on screen from `p` to the segment `a`–`b` (pixels), and where along it (0 at
/// `a`, 1 at `b`) the nearest point lies.
pub fn segment_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> (f64, f64) {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len2 = dx * dx + dy * dy;
    let k = if len2 > 0.0 { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2).clamp(0.0, 1.0) } else { 0.0 };
    let (qx, qy) = (a[0] + dx * k, a[1] + dy * k);
    (((p[0] - qx).powi(2) + (p[1] - qy).powi(2)).sqrt(), k)
}

/// A camera's projection for one viewport size, computed once: world points to pixels.
#[derive(Clone, Copy, Debug)]
pub struct Projector {
    cam: Camera,
    vp: Mat4,
    w: f64,
    h: f64,
}

impl Projector {
    /// The projection of `cam` onto a `w × h` viewport whose scene has radius `radius` (the
    /// clipping planes', as drawn).
    pub fn new(cam: &Camera, w: f64, h: f64, radius: f64) -> Projector {
        let (w, h) = (if w.is_finite() { w.max(1.0) } else { 1.0 }, if h.is_finite() { h.max(1.0) } else { 1.0 });
        Projector { cam: *cam, vp: cam.view_proj(w / h, radius.max(crate::camera::MIN_SIZE)), w, h }
    }

    /// The pixel (x right, y down) and depth (0 near … 1 far) of a world point in front of the
    /// eye ([`Camera::to_screen`]).
    pub fn screen(&self, p: Vec3) -> Option<[f64; 3]> {
        let c = self.vp.apply(p);
        if !(c[3] > 1e-12) {
            return None;
        }
        let (x, y, z) = (c[0] / c[3], c[1] / c[3], c[2] / c[3]);
        let s = [(x + 1.0) * 0.5 * self.w, (1.0 - y) * 0.5 * self.h, z];
        s.iter().all(|v| v.is_finite()).then_some(s)
    }

    /// The picking ray through pixel (x, y): origin and unit direction ([`Camera::ray`]).
    pub fn ray(&self, x: f64, y: f64) -> (Vec3, Vec3) {
        self.cam.ray(x, y, self.w, self.h)
    }

    pub fn size(&self) -> [f64; 2] {
        [self.w, self.h]
    }
}

/// Most pixels a [`DepthMap`] holds; a larger viewport is sampled more coarsely.
pub const MAX_DEPTH_PIXELS: usize = 2048 * 2048;
/// Depth (of the 0…1 range) a point may lie behind the drawn surface and still count as on it.
const DEPTH_ALLOWANCE: f64 = 2e-3;

/// The depth of the nearest drawn face at every pixel: which points are in front.
#[derive(Clone, Debug)]
pub struct DepthMap {
    proj: Projector,
    /// Map pixels per viewport pixel (1, or less for a large viewport).
    scale: f64,
    w: usize,
    h: usize,
    depth: Vec<f32>,
}

impl DepthMap {
    /// Rasterise `tris` with `proj`. `None` when the buffer cannot be allocated.
    pub fn new(proj: &Projector, tris: impl IntoIterator<Item = [Vec3; 3]>) -> Option<DepthMap> {
        let [vw, vh] = proj.size();
        let scale = ((MAX_DEPTH_PIXELS as f64) / (vw * vh)).sqrt().min(1.0);
        let (w, h) = (((vw * scale).ceil() as usize).max(1), ((vh * scale).ceil() as usize).max(1));
        let n = w.checked_mul(h)?;
        let mut depth = Vec::new();
        depth.try_reserve_exact(n).ok()?;
        depth.resize(n, f32::INFINITY);
        let mut map = DepthMap { proj: *proj, scale, w, h, depth };
        for p in tris {
            map.add(&p);
        }
        Some(map)
    }

    fn add(&mut self, p: &[Vec3; 3]) {
        let s = |q: Vec3| self.proj.screen(q).map(|v| [v[0] * self.scale, v[1] * self.scale, v[2]]);
        let (Some(a), Some(b), Some(c)) = (s(p[0]), s(p[1]), s(p[2])) else { return };
        let (w, depth) = (self.w, &mut self.depth);
        raster_tri(a, b, c, w, self.h, |x, y, z| {
            if let Some(d) = depth.get_mut(y * w + x) {
                *d = d.min(z);
            }
        });
    }

    /// Whether a world point is on screen and not behind a drawn face (a point on a face is in
    /// front of it). A point whose pixel no face covers, on a face's outline, is compared with
    /// the nearest face around it, so what lies behind an outline stays hidden.
    pub fn visible(&self, p: Vec3) -> bool {
        let Some([x, y, z]) = self.proj.screen(p) else { return false };
        let [vw, vh] = self.proj.size();
        if !(0.0..=vw).contains(&x) || !(0.0..=vh).contains(&y) || !(0.0..=1.0).contains(&z) {
            return false;
        }
        let (mx, my) = (((x * self.scale) as usize).min(self.w - 1), ((y * self.scale) as usize).min(self.h - 1));
        let at = |x: usize, y: usize| self.depth.get(y * self.w + x).copied().unwrap_or(f32::INFINITY);
        let mut d = at(mx, my);
        if d == f32::INFINITY {
            for y in my.saturating_sub(1)..=(my + 1).min(self.h - 1) {
                for x in mx.saturating_sub(1)..=(mx + 1).min(self.w - 1) {
                    d = d.min(at(x, y));
                }
            }
        }
        z <= f64::from(d) + DEPTH_ALLOWANCE
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StandardView;

    const SQUARE: [[Vec3; 3]; 2] = [
        [Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 0.0)],
        [Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 0.0), Vec3::new(0.0, 1.0, 0.0)],
    ];

    #[test]
    fn rays_hit_from_either_side_and_report_the_distance() {
        let down = Vec3::new(0.0, 0.0, -1.0);
        assert_eq!(ray_triangle(Vec3::new(0.7, 0.2, 5.0), down, &SQUARE[0]), Some(5.0));
        assert_eq!(ray_triangle(Vec3::new(0.7, 0.2, -2.0), Vec3::Z, &SQUARE[0]), Some(2.0));
        // Behind the origin, beside the triangle, in its plane.
        assert_eq!(ray_triangle(Vec3::new(0.7, 0.2, -2.0), down, &SQUARE[0]), None);
        assert_eq!(ray_triangle(Vec3::new(0.2, 0.7, 5.0), down, &SQUARE[0]), None);
        assert_eq!(ray_triangle(Vec3::new(-1.0, 0.5, 0.0), Vec3::X, &SQUARE[0]), None);
        // A slanted ray, a zero or non-finite one.
        let t = ray_triangle(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.5, 0.25, -1.0), &SQUARE[0]).unwrap();
        assert!((t - 1.0).abs() < 1e-15);
        assert_eq!(ray_triangle(Vec3::new(0.5, 0.2, 1.0), Vec3::ZERO, &SQUARE[0]), None);
        assert_eq!(ray_triangle(Vec3::new(f64::NAN, 0.2, 1.0), down, &SQUARE[0]), None);
        assert_eq!(ray_triangle(Vec3::new(0.5, 0.2, 1.0), Vec3::new(0.0, 0.0, f64::INFINITY), &SQUARE[0]), None);
        let huge = [Vec3::new(1e300, 0.0, 0.0), Vec3::new(-1e300, 1e300, 0.0), Vec3::new(-1e300, -1e300, 0.0)];
        assert_eq!(ray_triangle(Vec3::new(0.0, 0.0, 1.0), down, &huge), None, "overflowing edge functions miss");
    }

    /// Rays through the shared diagonal (inside the square) and corners of two triangles, at
    /// every angle, always hit one of them.
    #[test]
    fn no_ray_slips_through_a_shared_edge_or_vertex() {
        for k in 1..400 {
            let s = k as f64 / 400.0;
            for dir in [Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.3, -0.1, -1.0), Vec3::new(-1e-9, 2e-9, -1.0)] {
                let at = Vec3::new(s, s, 0.0) - dir * 3.0;
                assert!(first_hit(at, dir, SQUARE).is_some(), "{s} {dir:?}");
            }
        }
        for corner in [Vec3::ZERO, Vec3::new(1.0, 1.0, 0.0)] {
            assert!(first_hit(corner + Vec3::Z, Vec3::new(0.0, 0.0, -1.0), SQUARE).is_some());
        }
    }

    #[test]
    fn the_nearest_hit_wins_and_ties_go_to_the_lower_index() {
        let lift = |t: [Vec3; 3], z: f64| t.map(|p| p + Vec3::new(0.0, 0.0, z));
        let tris = [lift(SQUARE[0], 0.0), lift(SQUARE[0], 2.0), lift(SQUARE[0], 1.0)];
        let hit = first_hit(Vec3::new(0.7, 0.2, 10.0), Vec3::new(0.0, 0.0, -1.0), tris).unwrap();
        assert_eq!(hit, Hit { index: 1, t: 8.0 });
        let same = [SQUARE[0], SQUARE[0]];
        assert_eq!(first_hit(Vec3::new(0.7, 0.2, 1.0), Vec3::new(0.0, 0.0, -1.0), same).map(|h| h.index), Some(0));
        assert_eq!(first_hit(Vec3::new(0.7, 0.2, 1.0), Vec3::new(0.0, 0.0, -1.0), []), None);
    }

    #[test]
    fn screen_distances_to_segments() {
        assert_eq!(segment_distance([1.0, 1.0], [0.0, 0.0], [2.0, 0.0]), (1.0, 0.5));
        assert_eq!(segment_distance([3.0, 0.0], [0.0, 0.0], [2.0, 0.0]), (1.0, 1.0));
        assert_eq!(segment_distance([0.0, 2.0], [0.0, 0.0], [0.0, 0.0]), (2.0, 0.0));
    }

    #[test]
    fn the_projector_agrees_with_the_camera_and_the_depth_map_hides_what_is_behind() {
        let mut cam = Camera { distance: 10.0, ..Camera::default() };
        cam.set_view(StandardView::Top);
        cam.target = Vec3::new(0.5, 0.5, 0.0);
        let p = Projector::new(&cam, 400.0, 300.0, 2.0);
        let q = Vec3::new(0.25, 0.75, 0.0);
        let (x, y, z) = cam.to_screen(q, 400.0, 300.0, 2.0).unwrap();
        assert_eq!(p.screen(q), Some([x, y, z]));
        // The ray through a point's pixel passes through the point.
        let (o, d) = p.ray(x, y);
        assert!((q - o).cross(d).len() < 1e-9);
        let map = DepthMap::new(&p, SQUARE).unwrap();
        assert!(map.visible(q), "a point on the face");
        assert!(map.visible(q + Vec3::Z * 0.5), "a point above it");
        assert!(!map.visible(q - Vec3::Z * 0.5), "a point under it");
        assert!(map.visible(Vec3::new(1.5, 0.5, -0.5)), "beside it, nothing covers the point");
        assert!(!map.visible(Vec3::new(500.0, 0.5, 0.0)), "off screen");
        assert!(!map.visible(Vec3::new(f64::NAN, 0.5, 0.0)));
        // A huge viewport is sampled coarsely.
        let big = DepthMap::new(&Projector::new(&cam, 8192.0, 8192.0, 2.0), SQUARE).unwrap();
        assert!(big.depth.len() <= MAX_DEPTH_PIXELS + 2 * 2048 + 1);
        let hostile = Projector::new(&cam, f64::NAN, -1.0, f64::INFINITY);
        assert!(DepthMap::new(&hostile, SQUARE).is_some());
    }
}
