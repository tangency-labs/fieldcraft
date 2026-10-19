//! An orbit camera about a target point, Z up, with perspective or orthographic projection.
//!
//! Adapted from SolidCraft's `crates/render/src/camera.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.

use serde::{Deserialize, Serialize};

/// The smallest size the camera frames or clips at: a micrometre (model lengths are metres).
pub const MIN_SIZE: f64 = 1e-6;

use crate::pick::Projector;
use crate::{Aabb3, Vec3};

/// Column-major 4×4 matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mat4(pub [[f64; 4]; 4]);

impl Mat4 {
    pub fn mul(&self, o: &Mat4) -> Mat4 {
        let mut r = [[0.0; 4]; 4];
        for (c, col) in r.iter_mut().enumerate() {
            for (rr, cell) in col.iter_mut().enumerate() {
                *cell = (0..4).map(|k| self.0[k][rr] * o.0[c][k]).sum();
            }
        }
        Mat4(r)
    }
    /// Transform a point; returns (x, y, z, w).
    pub fn apply(&self, p: Vec3) -> [f64; 4] {
        let v = [p.x, p.y, p.z, 1.0];
        let mut out = [0.0; 4];
        for (r, o) in out.iter_mut().enumerate() {
            *o = (0..4).map(|c| self.0[c][r] * v[c]).sum();
        }
        out
    }
    pub fn to_f32(&self) -> [[f32; 4]; 4] {
        self.0.map(|c| c.map(|x| x as f32))
    }
}

/// Named view orientations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StandardView {
    Front,
    Back,
    Top,
    Bottom,
    Left,
    Right,
    Iso,
}

impl StandardView {
    pub const ALL: [StandardView; 7] = [
        StandardView::Front,
        StandardView::Back,
        StandardView::Left,
        StandardView::Right,
        StandardView::Top,
        StandardView::Bottom,
        StandardView::Iso,
    ];

    pub fn parse(s: &str) -> Option<StandardView> {
        Some(match s.to_ascii_lowercase().as_str() {
            "front" => StandardView::Front,
            "back" => StandardView::Back,
            "top" => StandardView::Top,
            "bottom" => StandardView::Bottom,
            "left" => StandardView::Left,
            "right" => StandardView::Right,
            "iso" | "isometric" => StandardView::Iso,
            _ => return None,
        })
    }
    pub fn label(self) -> &'static str {
        match self {
            StandardView::Front => "Front",
            StandardView::Back => "Back",
            StandardView::Top => "Top",
            StandardView::Bottom => "Bottom",
            StandardView::Left => "Left",
            StandardView::Right => "Right",
            StandardView::Iso => "Iso",
        }
    }
    /// (yaw, pitch) in radians. Yaw 0 looks along +Y (front view); the iso view looks from the
    /// front-right-top corner. Pitch is the elevation above the XY plane.
    pub fn angles(self) -> (f64, f64) {
        use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
        match self {
            StandardView::Front => (0.0, 0.0),
            StandardView::Back => (PI, 0.0),
            StandardView::Right => (-FRAC_PI_2, 0.0),
            StandardView::Left => (FRAC_PI_2, 0.0),
            StandardView::Top => (0.0, FRAC_PI_2),
            StandardView::Bottom => (0.0, -FRAC_PI_2),
            // atan(1/√2): the elevation of the cube diagonal.
            StandardView::Iso => (-FRAC_PI_4, 0.6154797086703874),
        }
    }
}

/// Where a camera looks from, in degrees as people and agents give it: `yaw` turns about Z (0
/// looks along +Y, the front view; −45 is the iso view's), `pitch` is the elevation above the XY
/// plane and `distance`, when given, the distance from the target (a zoom in orthographic views).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Orientation {
    pub yaw: f64,
    pub pitch: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub distance: Option<f64>,
}

impl Orientation {
    /// Closest and farthest camera distance (metres: a micrometre to a million kilometres).
    pub const DISTANCE: (f64, f64) = (MIN_SIZE, 1e9);

    /// The orientation of a standard view.
    pub fn of(view: StandardView) -> Self {
        let (yaw, pitch) = view.angles();
        Orientation { yaw: yaw.to_degrees(), pitch: pitch.to_degrees(), distance: None }
    }

    /// A checked orientation: `None` when an angle or the distance is not finite or the distance
    /// is not above 0. The yaw is wrapped into [−180, 180), the pitch clamped to ±90 and the
    /// distance to [`Orientation::DISTANCE`].
    pub fn checked(yaw: f64, pitch: f64, distance: Option<f64>) -> Option<Self> {
        if !(yaw.is_finite() && pitch.is_finite()) || distance.is_some_and(|d| !(d.is_finite() && d > 0.0)) {
            return None;
        }
        let yaw = (yaw + 180.0).rem_euclid(360.0) - 180.0;
        let (near, far) = Self::DISTANCE;
        Some(Orientation { yaw, pitch: pitch.clamp(-90.0, 90.0), distance: distance.map(|d| d.clamp(near, far)) })
    }
}

/// The whole view, as the `free` part of [`Camera::fit_in`].
pub const FULL_VIEW: [f64; 4] = [0.0, 0.0, 1.0, 1.0];

/// The share of the free part's width and height a fit leaves empty on each side.
pub const FIT_MARGIN: f64 = 0.06;

/// `free` checked: finite, within the view, at least a tenth of it wide and high (else the whole
/// view).
fn checked_free(free: [f64; 4]) -> [f64; 4] {
    let [x0, y0, x1, y1] = free.map(|v| if v.is_finite() { v.clamp(0.0, 1.0) } else { f64::NAN });
    if x1 - x0 >= 0.1 && y1 - y0 >= 0.1 { [x0, y0, x1, y1] } else { FULL_VIEW }
}

/// Orbit camera around `target`, Z up.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Camera {
    pub target: Vec3,
    /// Rotation about Z; 0 = looking along +Y (from the front).
    pub yaw: f64,
    /// Elevation above the XY plane, in [−π/2, π/2].
    pub pitch: f64,
    /// Distance from target to eye.
    pub distance: f64,
    /// Vertical field of view (radians); 0 = orthographic.
    pub fov: f64,
}

impl Default for Camera {
    fn default() -> Self {
        let (yaw, pitch) = StandardView::Iso.angles();
        Camera { target: Vec3::ZERO, yaw, pitch, distance: 200.0, fov: 0.0 }
    }
}

impl Camera {
    /// Direction from the target toward the eye.
    pub fn back(&self) -> Vec3 {
        let (cp, sp) = (self.pitch.cos(), self.pitch.sin());
        Vec3::new(-self.yaw.sin() * cp, -self.yaw.cos() * cp, sp).normalized().unwrap_or(Vec3::Z)
    }
    pub fn eye(&self) -> Vec3 {
        self.target + self.back() * self.distance
    }
    /// Camera basis: right, up, back (orthonormal).
    pub fn basis(&self) -> (Vec3, Vec3, Vec3) {
        let b = self.back();
        let right = Vec3::new(self.yaw.cos(), -self.yaw.sin(), 0.0);
        let up = b.cross(right).normalized().unwrap_or(Vec3::Z);
        (right, up, b)
    }
    pub fn view(&self) -> Mat4 {
        let (r, u, b) = self.basis();
        let e = self.eye();
        Mat4([[r.x, u.x, b.x, 0.0], [r.y, u.y, b.y, 0.0], [r.z, u.z, b.z, 0.0], [-r.dot(e), -u.dot(e), -b.dot(e), 1.0]])
    }
    /// Half the visible height at the target distance.
    pub fn half_height(&self) -> f64 {
        if self.fov > 0.0 { self.distance * (self.fov * 0.5).tan() } else { self.distance * 0.4 }
    }
    /// Model length one pixel spans at the target, in a view `height_px` pixels high (the scale
    /// symbols of a fixed screen size are drawn at); `None` for a size or camera that has none.
    pub fn length_per_pixel(&self, height_px: f64) -> Option<f64> {
        let l = 2.0 * self.half_height() / height_px;
        (height_px > 0.0 && l.is_finite() && l > 0.0).then_some(l)
    }
    /// Near/far planes that keep a scene of radius `radius` around the target visible.
    pub fn clip(&self, radius: f64) -> (f64, f64) {
        let r = radius.max(MIN_SIZE);
        if self.fov > 0.0 {
            ((self.distance - r * 2.0).max(self.distance * 1e-3).max(MIN_SIZE), self.distance + r * 2.0)
        } else {
            // A margin of a fiftieth of the radius beyond the scene: relative, so the depth
            // range (and what a depth difference means) is the same at every scale.
            (self.distance - r * 4.02, self.distance + r * 4.02)
        }
    }
    /// Projection to clip space (wgpu style, depth 0..1).
    pub fn proj(&self, aspect: f64, near: f64, far: f64) -> Mat4 {
        let aspect = if aspect.is_finite() && aspect > 1e-6 { aspect } else { 1.0 };
        let (n, f) = (near, far.max(near + 1e-6));
        if self.fov > 0.0 {
            let t = 1.0 / (self.fov * 0.5).tan();
            Mat4([[t / aspect, 0.0, 0.0, 0.0], [0.0, t, 0.0, 0.0], [0.0, 0.0, f / (n - f), -1.0], [0.0, 0.0, n * f / (n - f), 0.0]])
        } else {
            let h = self.half_height();
            let w = h * aspect;
            Mat4([[1.0 / w, 0.0, 0.0, 0.0], [0.0, 1.0 / h, 0.0, 0.0], [0.0, 0.0, 1.0 / (n - f), 0.0], [0.0, 0.0, n / (n - f), 1.0]])
        }
    }
    pub fn view_proj(&self, aspect: f64, radius: f64) -> Mat4 {
        let (n, f) = self.clip(radius);
        self.proj(aspect, n, f).mul(&self.view())
    }
    /// World point → pixel coordinates (x right, y down) and depth (0 near … 1 far).
    /// For many points, a [`Projector`] computes the projection once.
    pub fn to_screen(&self, p: Vec3, w: f64, h: f64, radius: f64) -> Option<(f64, f64, f64)> {
        Projector::new(self, w, h, radius).screen(p).map(|[x, y, z]| (x, y, z))
    }
    /// Picking ray through pixel (x, y): origin and unit direction.
    pub fn ray(&self, x: f64, y: f64, w: f64, h: f64) -> (Vec3, Vec3) {
        let (r, u, b) = self.basis();
        let hh = self.half_height();
        let aspect = w / h.max(1.0);
        let nx = (x / w.max(1.0)) * 2.0 - 1.0;
        let ny = 1.0 - (y / h.max(1.0)) * 2.0;
        if self.fov > 0.0 {
            let d = (r * (nx * hh * aspect) + u * (ny * hh) - b * self.distance).normalized().unwrap_or(-b);
            (self.eye(), d)
        } else {
            let o = self.target + r * (nx * hh * aspect) + u * (ny * hh) + b * (self.distance * 4.0);
            (o, -b)
        }
    }
    /// The point under pixel (x, y) on the plane through the target facing the eye: the anchor
    /// that stays put when zooming at the cursor.
    pub fn target_plane_point(&self, x: f64, y: f64, w: f64, h: f64) -> Option<Vec3> {
        let (o, d) = self.ray(x, y, w, h);
        let b = self.back();
        let dn = d.dot(b);
        if dn.abs() < 1e-12 {
            return None;
        }
        let p = o + d * ((self.target - o).dot(b) / dn);
        p.is_finite().then_some(p)
    }
    /// Zoom by `factor` (< 1 zooms in) toward a world point that stays put on screen.
    pub fn zoom_at(&mut self, factor: f64, anchor: Option<Vec3>) {
        let f = if factor.is_finite() { factor.clamp(0.05, 20.0) } else { 1.0 };
        let nd = (self.distance * f).clamp(Orientation::DISTANCE.0, Orientation::DISTANCE.1);
        if let Some(a) = anchor.filter(|a| a.is_finite()) {
            self.target = a + (self.target - a) * (nd / self.distance);
        }
        self.distance = nd;
    }
    /// Pan by a screen delta in pixels for a viewport of height `h`.
    pub fn pan(&mut self, dx: f64, dy: f64, h: f64) {
        if !(dx.is_finite() && dy.is_finite()) {
            return;
        }
        let (r, u, _) = self.basis();
        let k = 2.0 * self.half_height() / h.max(1.0);
        self.target = self.target - r * (dx * k) + u * (dy * k);
    }
    /// Orbit by pixel deltas. The model follows the cursor: dragging right turns the side
    /// facing the viewer to the right.
    pub fn orbit(&mut self, dx: f64, dy: f64) {
        if !(dx.is_finite() && dy.is_finite()) {
            return;
        }
        let lim = std::f64::consts::FRAC_PI_2 - 1e-4;
        self.yaw = (self.yaw + dx * 0.008) % std::f64::consts::TAU;
        self.pitch = (self.pitch + dy * 0.008).clamp(-lim, lim);
    }
    /// Orbit by pixel deltas about `pivot`, which keeps its place on screen: the eye turns
    /// about the pivot and the target follows, so the pivot's view coordinates (and so its
    /// screen position, in perspective and orthographic views) are the same after the turn.
    pub fn orbit_about(&mut self, pivot: Vec3, dx: f64, dy: f64) {
        if !pivot.is_finite() {
            self.orbit(dx, dy);
            return;
        }
        let (r, u, b) = self.basis();
        let v = pivot - self.eye();
        let (vx, vy, vz) = (v.dot(r), v.dot(u), v.dot(b));
        self.orbit(dx, dy);
        let (r, u, b) = self.basis();
        let target = pivot - (r * vx + u * vy + b * vz) - b * self.distance;
        if target.is_finite() {
            self.target = target;
        }
    }
    pub fn set_view(&mut self, v: StandardView) {
        let lim = std::f64::consts::FRAC_PI_2 - 1e-6;
        let (y, p) = v.angles();
        self.yaw = y;
        self.pitch = p.clamp(-lim, lim);
    }
    /// Turn to `o` (and move to its distance, when it has one); the pitch stays just short of the
    /// poles, as an orbit's does.
    pub fn orient(&mut self, o: &Orientation) {
        let lim = std::f64::consts::FRAC_PI_2 - 1e-6;
        if o.yaw.is_finite() && o.pitch.is_finite() {
            self.yaw = o.yaw.to_radians();
            self.pitch = o.pitch.to_radians().clamp(-lim, lim);
        }
        if let Some(d) = o.distance.filter(|d| d.is_finite() && *d > 0.0) {
            self.distance = d.clamp(Orientation::DISTANCE.0, Orientation::DISTANCE.1);
        }
    }
    /// Frame a bounding box in the whole view of a square viewport (see [`Camera::fit_in`]).
    pub fn fit(&mut self, b: &Aabb3) {
        self.fit_in(b, 1.0, FULL_VIEW);
    }

    /// Frame a bounding box in part of the view, keeping the direction of view: the smallest
    /// distance (zoom, in orthographic views) at which all eight corners of `b` project inside
    /// `free` shrunk by [`FIT_MARGIN`] on each side, and the target moved across the view so the
    /// box sits centred in it. `free` is `[x0, y0, x1, y1]` as fractions of the view's width and
    /// height, y down (the part overlays leave free); `aspect` is the view's width over height.
    ///
    /// The corners' projections are exact, for orthographic and perspective cameras alike, so
    /// a long beam seen side on fills the width instead of a sphere around it the height. With
    /// `s(z) = α·d − β·z` the half height of the view at depth `z` toward the eye (orthographic
    /// α = 0.4, β = 0; perspective α = β = tan(fov/2)), a corner at `(x, z)` in view
    /// coordinates (relative to the box centre) shows at `(x − t)/(aspect·s(z))` in normalized
    /// device coordinates for a target offset `t`. Keeping every corner between the free part's
    /// edges `a0 … a1` is linear in `(t, d)`: a pair of corners `i`, `j` needs
    /// `d ≥ ((x_i − x_j)/aspect + β·(a1·z_i − a0·z_j)) / (α·(a1 − a0))`; the largest over all
    /// pairs and both axes is the distance, and `t` is the middle of what is then allowed.
    pub fn fit_in(&mut self, b: &Aabb3, aspect: f64, free: [f64; 4]) {
        if b.is_empty() || !(b.min.is_finite() && b.max.is_finite()) {
            return;
        }
        let aspect = if aspect.is_finite() && aspect > 1e-6 { aspect.min(1e6) } else { 1.0 };
        let [x0, y0, x1, y1] = checked_free(free);
        // The free part's edges in normalized device coordinates (y up), inside the margin.
        let (mx, my) = ((x1 - x0) * FIT_MARGIN, (y1 - y0) * FIT_MARGIN);
        let ax = (2.0 * (x0 + mx) - 1.0, 2.0 * (x1 - mx) - 1.0);
        let ay = (1.0 - 2.0 * (y1 - my), 1.0 - 2.0 * (y0 + my));
        let (r, u, back) = self.basis();
        let c = b.center();
        let corners = b.corners().map(|p| {
            let p = p - c;
            [p.dot(r), p.dot(u), p.dot(back)]
        });
        // The half height per unit of distance (`half_height`), and its change with depth.
        let alpha = Camera { distance: 1.0, ..*self }.half_height();
        let beta = if self.fov > 0.0 { alpha } else { 0.0 };
        let size = (b.diagonal() * 0.5).max(MIN_SIZE);
        // Nearest distance: a tiny box still has a size, and in perspective every corner is in
        // front of the eye (an orthographic distance is only a zoom).
        let mut d = size * 1e-3 / alpha;
        if beta > 0.0 {
            d = corners.iter().map(|q| q[2] + size * 1e-3).fold(d, f64::max);
        }
        for (axis, (a0, a1), k) in [(0, ax, aspect), (1, ay, 1.0)] {
            for i in &corners {
                for j in &corners {
                    let need = ((i[axis] - j[axis]) / k + beta * (a1 * i[2] - a0 * j[2])) / (alpha * (a1 - a0));
                    d = d.max(need);
                }
            }
        }
        let d = d.clamp(Orientation::DISTANCE.0, Orientation::DISTANCE.1);
        // The target offset along each axis: the middle of the offsets that keep every corner in.
        let offset = |axis: usize, (a0, a1): (f64, f64), k: f64| {
            let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
            for q in &corners {
                let s = k * (alpha * d - beta * q[2]);
                lo = lo.max(q[axis] - a1 * s);
                hi = hi.min(q[axis] - a0 * s);
            }
            0.5 * (lo + hi)
        };
        let (tx, ty) = (offset(0, ax, aspect), offset(1, ay, 1.0));
        let target = c + r * tx + u * ty;
        if target.is_finite() && d.is_finite() {
            self.target = target;
            self.distance = d;
        }
    }

    /// The distance (zoom, in orthographic views) at which all eight corners of `b` project
    /// inside the view shrunk by [`FIT_MARGIN`] on each side, keeping the target where it is (at
    /// the view's centre): [`Camera::fit_in`] for a camera aimed at a chosen point, which the box
    /// need not be centred on. With `s(z) = α·d − β·z` as there, a corner at `(x, y, z)` relative
    /// to the target stays in when `|x| ≤ a·aspect·s(z)` and `|y| ≤ a·s(z)` for the margin's
    /// edge `a`, so `d ≥ (|x|/(a·aspect) + β·z)/α`, and likewise for `y`.
    pub fn fit_about(&mut self, b: &Aabb3, aspect: f64) {
        if b.is_empty() || !(b.min.is_finite() && b.max.is_finite() && self.target.is_finite()) {
            return;
        }
        let aspect = if aspect.is_finite() && aspect > 1e-6 { aspect.min(1e6) } else { 1.0 };
        let a = 1.0 - 2.0 * FIT_MARGIN;
        let (r, u, back) = self.basis();
        let alpha = Camera { distance: 1.0, ..*self }.half_height();
        let beta = if self.fov > 0.0 { alpha } else { 0.0 };
        let size = (b.diagonal() * 0.5).max(MIN_SIZE);
        let mut d = size * 1e-3 / alpha;
        for p in b.corners() {
            let p = p - self.target;
            let (x, y, z) = (p.dot(r), p.dot(u), p.dot(back));
            d = d.max((x.abs() / (a * aspect) + beta * z) / alpha).max((y.abs() / a + beta * z) / alpha);
            if beta > 0.0 {
                d = d.max(z + size * 1e-3);
            }
        }
        if d.is_finite() {
            self.distance = d.clamp(Orientation::DISTANCE.0, Orientation::DISTANCE.1);
        }
    }

    pub fn is_valid(&self) -> bool {
        self.target.is_finite() && self.yaw.is_finite() && self.pitch.is_finite() && self.distance.is_finite() && self.distance > 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A box off to one side of the target fits with the target still at the centre and the
    /// box's farthest corner at the margin.
    #[test]
    fn fit_about_keeps_the_target_centred() {
        let b = Aabb3 { min: Vec3::new(10.0, -5.0, 0.0), max: Vec3::new(30.0, 5.0, 8.0) };
        for fov in [0.0, 0.6] {
            let mut c = Camera { target: Vec3::new(12.0, 0.0, 2.0), fov, ..Default::default() };
            c.fit_about(&b, 1.5);
            let (w, h) = (300.0, 200.0);
            let (x, y, _) = c.to_screen(c.target, w, h, 50.0).unwrap();
            assert!((x - 150.0).abs() < 1e-6 && (y - 100.0).abs() < 1e-6);
            let worst = b
                .corners()
                .iter()
                .map(|p| c.to_screen(*p, w, h, 50.0).unwrap())
                .map(|(x, y, _)| ((x - 150.0).abs() / 150.0).max((y - 100.0).abs() / 100.0))
                .fold(0.0, f64::max);
            assert!((worst - (1.0 - 2.0 * FIT_MARGIN)).abs() < 1e-6, "fov {fov}: {worst}");
        }
    }

    #[test]
    fn front_view_projects_target_to_centre() {
        let mut c = Camera { target: Vec3::new(1.0, 2.0, 3.0), distance: 50.0, ..Default::default() };
        c.set_view(StandardView::Front);
        let (x, y, z) = c.to_screen(c.target, 800.0, 600.0, 10.0).unwrap();
        assert!((x - 400.0).abs() < 1e-6 && (y - 300.0).abs() < 1e-6 && (0.0..=1.0).contains(&z));
        // Front view looks along +Y: +X is to the right, +Z up.
        let (xr, _, _) = c.to_screen(c.target + Vec3::X, 800.0, 600.0, 10.0).unwrap();
        let (_, yu, _) = c.to_screen(c.target + Vec3::Z, 800.0, 600.0, 10.0).unwrap();
        assert!(xr > 400.0 && yu < 300.0);
        let (o, d) = c.ray(400.0, 300.0, 800.0, 600.0);
        assert!(d.dist(Vec3::Y) < 1e-9, "{d:?}");
        assert!((o - c.target).cross(d).len() < 1e-6);
    }

    #[test]
    fn horizontal_orbit_follows_the_cursor() {
        let mut c = Camera { distance: 50.0, ..Default::default() };
        c.set_view(StandardView::Front);
        let near = Vec3::new(0.0, -1.0, 0.0);
        let x0 = c.to_screen(near, 800.0, 600.0, 10.0).unwrap().0;
        c.orbit(20.0, 0.0);
        let x1 = c.to_screen(near, 800.0, 600.0, 10.0).unwrap().0;
        assert!(x1 > x0, "dragging right must move the near side right: {x0} -> {x1}");
    }

    #[test]
    fn top_view_and_perspective() {
        let mut c = Camera { distance: 100.0, fov: 0.8, ..Default::default() };
        c.set_view(StandardView::Top);
        let (xr, _, _) = c.to_screen(Vec3::X, 800.0, 600.0, 10.0).unwrap();
        let (_, yu, _) = c.to_screen(Vec3::Y, 800.0, 600.0, 10.0).unwrap();
        assert!(xr > 400.0 && yu < 300.0, "{xr} {yu}");
        let near = c.to_screen(Vec3::new(0.0, 0.0, 10.0), 800.0, 600.0, 20.0).unwrap().2;
        let far = c.to_screen(Vec3::new(0.0, 0.0, -10.0), 800.0, 600.0, 20.0).unwrap().2;
        assert!(near < far);
    }

    #[test]
    fn zoom_at_the_cursor_keeps_the_point_under_it() {
        let mut c = Camera { distance: 100.0, ..Default::default() };
        c.set_view(StandardView::Front);
        let a = c.target_plane_point(600.0, 200.0, 800.0, 600.0).unwrap();
        c.zoom_at(0.5, Some(a));
        let (x, y, _) = c.to_screen(a, 800.0, 600.0, 10.0).unwrap();
        assert!((x - 600.0).abs() < 1e-6 && (y - 200.0).abs() < 1e-6, "{x} {y}");
    }

    /// Orbiting about a pivot off the target keeps the pivot where it was on screen, in both
    /// projections and from any view.
    #[test]
    fn orbit_about_keeps_the_pivot_on_screen() {
        let pivot = Vec3::new(37.0, -12.0, 8.0);
        for fov in [0.0, 0.7] {
            for v in StandardView::ALL {
                let mut c = Camera { target: Vec3::new(5.0, 3.0, -2.0), distance: 150.0, fov, ..Default::default() };
                c.set_view(v);
                let before = c.to_screen(pivot, 800.0, 600.0, 100.0).unwrap();
                for (dx, dy) in [(40.0, 0.0), (0.0, -25.0), (-60.0, 33.0)] {
                    c.orbit_about(pivot, dx, dy);
                    let after = c.to_screen(pivot, 800.0, 600.0, 100.0).unwrap();
                    assert!(
                        (after.0 - before.0).abs() < 1e-6 && (after.1 - before.1).abs() < 1e-6,
                        "{} fov {fov}: {before:?} -> {after:?}",
                        v.label()
                    );
                }
                assert!(c.is_valid());
            }
        }
        let mut c = Camera::default();
        c.orbit_about(Vec3::new(f64::NAN, 0.0, 0.0), 10.0, 10.0);
        assert!(c.is_valid() && c.target == Vec3::ZERO);
    }

    /// A wheel notch zooms by the same factor in orthographic and perspective views: the visible
    /// height scales by it, and the point under the cursor stays put (#273).
    #[test]
    fn zoom_scales_the_view_alike_in_both_projections() {
        for fov in [0.0, 0.7] {
            let mut c = Camera { distance: 300.0, fov, ..Default::default() };
            c.set_view(StandardView::Front);
            let h0 = c.half_height();
            let a = c.target_plane_point(650.0, 120.0, 800.0, 600.0).unwrap();
            c.zoom_at(0.8, Some(a));
            assert!((c.half_height() / h0 - 0.8).abs() < 1e-12, "fov {fov}");
            let (x, y, _) = c.to_screen(a, 800.0, 600.0, 10.0).unwrap();
            assert!((x - 650.0).abs() < 1e-6 && (y - 120.0).abs() < 1e-6, "fov {fov}: {x} {y}");
        }
    }

    #[test]
    fn fit_and_hostile() {
        let mut c = Camera::default();
        let b = Aabb3 { min: Vec3::ZERO, max: Vec3::new(40.0, 30.0, 20.0) };
        c.fit(&b);
        assert!(c.target.dist(Vec3::new(20.0, 15.0, 10.0)) < 1e-9);
        c.zoom_at(f64::NAN, Some(Vec3::new(f64::NAN, 0.0, 0.0)));
        c.orbit(1e9, -1e9);
        c.orbit(f64::NAN, 0.0);
        c.pan(f64::INFINITY, 0.0, 100.0);
        assert!(c.is_valid());
        c.fit(&Aabb3::EMPTY);
        assert!(c.is_valid());
        for v in StandardView::ALL {
            assert_eq!(StandardView::parse(v.label()), Some(v));
        }
    }

    /// Every corner of the box projects inside `free` (in pixels of a `w` × `h` view), and the
    /// fit is tight: the corners reach the margin on one axis.
    fn assert_framed(c: &Camera, b: &Aabb3, w: f64, h: f64, free: [f64; 4], what: &str) {
        let [x0, y0, x1, y1] = free;
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        let proj = crate::pick::Projector::new(c, w, h, b.diagonal());
        for (k, p) in b.corners().into_iter().enumerate() {
            let [x, y, z] = proj.screen(p).unwrap();
            assert!((0.0..=1.0).contains(&z), "{what}: corner {k} clipped (depth {z})");
            let eps = 1e-6 * (w + h);
            assert!(
                x >= x0 * w - eps && x <= x1 * w + eps && y >= y0 * h - eps && y <= y1 * h + eps,
                "{what}: corner {k} at ({x}, {y}) outside {free:?} of {w}×{h}"
            );
            lo = [lo[0].min(x), lo[1].min(y)];
            hi = [hi[0].max(x), hi[1].max(y)];
        }
        // Tight: the box spans the free width or height less the margins.
        let span = ((hi[0] - lo[0]) / ((x1 - x0) * w)).max((hi[1] - lo[1]) / ((y1 - y0) * h));
        assert!((span - (1.0 - 2.0 * FIT_MARGIN)).abs() < 1e-6, "{what}: spans {span} of the free part");
    }

    #[test]
    fn fit_in_frames_every_corner_in_the_free_part() {
        let boxes = [
            ("long beam", Aabb3 { min: Vec3::new(0.0, -5.0, -10.0), max: Vec3::new(200.0, 5.0, 10.0) }),
            ("cube off the origin", Aabb3 { min: Vec3::new(100.0, 50.0, -20.0), max: Vec3::new(140.0, 90.0, 20.0) }),
            ("tall post", Aabb3 { min: Vec3::new(-2.0, -2.0, 0.0), max: Vec3::new(2.0, 2.0, 300.0) }),
            ("flat plate", Aabb3 { min: Vec3::ZERO, max: Vec3::new(120.0, 80.0, 0.0) }),
        ];
        // Whole view; right of a legend and above the navigation bar; a narrow middle strip.
        let frees = [FULL_VIEW, [0.22, 0.0, 0.85, 0.93], [0.3, 0.1, 0.7, 0.9]];
        let sizes = [(1280.0, 800.0), (400.0, 1000.0), (800.0, 800.0), (2400.0, 600.0)];
        for (name, b) in &boxes {
            for v in StandardView::ALL {
                for fov in [0.0, 0.7] {
                    for (w, h) in sizes {
                        for free in frees {
                            let mut c = Camera { fov, ..Default::default() };
                            c.set_view(v);
                            c.fit_in(b, w / h, free);
                            assert!(c.is_valid());
                            assert_framed(&c, b, w, h, free, &format!("{name}, {} view, fov {fov}, {w}×{h}, {free:?}", v.label()));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn fit_in_refuses_hostile_input() {
        let b = Aabb3 { min: Vec3::ZERO, max: Vec3::new(10.0, 10.0, 10.0) };
        let mut c = Camera::default();
        c.fit_in(&b, f64::NAN, [f64::NAN, 0.0, 2.0, -1.0]);
        assert!(c.is_valid());
        let before = c;
        c.fit_in(&Aabb3 { min: Vec3::new(f64::NEG_INFINITY, 0.0, 0.0), max: Vec3::ZERO }, 1.0, FULL_VIEW);
        assert_eq!(c, before);
        // A point: framed at a tiny distance, still valid.
        c.fit_in(&Aabb3 { min: Vec3::ZERO, max: Vec3::ZERO }, 1.0, FULL_VIEW);
        assert!(c.is_valid());
    }

    #[test]
    fn orientation_in_degrees() {
        // The iso view's angles, given in degrees, turn the camera as the iso view does.
        let (mut a, mut b) = (Camera::default(), Camera::default());
        a.set_view(StandardView::Front);
        a.orient(&Orientation::of(StandardView::Iso));
        b.set_view(StandardView::Iso);
        assert!((a.yaw - b.yaw).abs() < 1e-12 && (a.pitch - b.pitch).abs() < 1e-12);
        let o = Orientation::checked(-45.0, 35.264, Some(300.0)).unwrap();
        a.orient(&o);
        assert_eq!(a.distance, 300.0);
        // Wrapped, clamped, refused.
        let o = Orientation::checked(540.0, 1e6, Some(1e300)).unwrap();
        assert_eq!((o.yaw, o.pitch, o.distance), (-180.0, 90.0, Some(Orientation::DISTANCE.1)));
        a.orient(&o);
        assert!(a.is_valid() && a.pitch < std::f64::consts::FRAC_PI_2, "{a:?}");
        for (y, p, d) in
            [(f64::NAN, 0.0, None), (0.0, f64::INFINITY, None), (0.0, 0.0, Some(0.0)), (0.0, 0.0, Some(-1.0)), (0.0, 0.0, Some(f64::NAN))]
        {
            assert_eq!(Orientation::checked(y, p, d), None);
        }
        // A hostile orientation built by hand leaves the camera as it was.
        let before = a;
        a.orient(&Orientation { yaw: f64::NAN, pitch: 0.0, distance: Some(f64::INFINITY) });
        assert_eq!(a, before);
    }
}
