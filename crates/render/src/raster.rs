//! A drawable scene (coloured triangles, lines and points) and a CPU rasterizer for headless
//! snapshots: z-buffered, lit triangles and depth-tested lines, rendered at 2× and downsampled
//! for anti-aliasing; a scene's [`Scene::on_top`] lines are drawn over everything. The GPU
//! viewport draws the same [`Scene`].
//!
//! Adapted from SolidCraft's `crates/render/src/raster.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.

use crate::pick::Projector;
use crate::{Camera, Vec3, font};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rgb(pub u8, pub u8, pub u8);

impl Rgb {
    fn f(self) -> [f64; 3] {
        [self.0 as f64 / 255.0, self.1 as f64 / 255.0, self.2 as f64 / 255.0]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SceneTri {
    pub p: [Vec3; 3],
    /// Unit normal; lighting is two-sided.
    pub normal: Vec3,
    pub color: Rgb,
}

/// A translucent triangle (a datum plane's fill): drawn after the opaque faces, depth-tested
/// against them without hiding what is behind, blended at `alpha` (0–1).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SceneFill {
    pub tri: SceneTri,
    pub alpha: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SceneLine {
    pub a: Vec3,
    pub b: Vec3,
    pub color: Rgb,
    /// Width in output pixels.
    pub width: f64,
}

/// A node marker: a square `size` output pixels across.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScenePoint {
    pub p: Vec3,
    pub color: Rgb,
    pub size: f64,
}

/// Text beside a point (a reference point's name), drawn over everything else. The GPU viewport
/// draws labels with the interface's font; the CPU rasterizer with [`crate::font`].
#[derive(Clone, Debug, PartialEq)]
pub struct SceneLabel {
    pub p: Vec3,
    pub text: String,
    pub color: Rgb,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Scene {
    pub tris: Vec<SceneTri>,
    /// Translucent triangles ([`SceneFill`]).
    pub fills: Vec<SceneFill>,
    pub lines: Vec<SceneLine>,
    pub points: Vec<ScenePoint>,
    pub labels: Vec<SceneLabel>,
    /// Lines drawn over everything, after the depth-tested layers (a concentrated load at a
    /// point inside the part stays visible from any side). A zero-length one is a round dot.
    pub on_top: Vec<SceneLine>,
    /// Background gradient top → bottom.
    pub background: Option<(Rgb, Rgb)>,
    /// Scene radius around the camera target, for the clipping planes.
    pub radius: f64,
}

impl Scene {
    /// Add `other`'s primitives (not its background or radius) to this scene.
    pub fn append(&mut self, other: Scene) {
        self.tris.extend(other.tris);
        self.fills.extend(other.fills);
        self.lines.extend(other.lines);
        self.points.extend(other.points);
        self.labels.extend(other.labels);
        self.on_top.extend(other.on_top);
    }

    /// Line widths and point sizes (given in points) in pixels at `pixels_per_point`, as the GPU
    /// viewport draws them: a CPU render at device resolution then matches the window.
    pub fn scale_marks(&mut self, pixels_per_point: f64) {
        if !(pixels_per_point.is_finite() && pixels_per_point > 0.0) {
            return;
        }
        for l in self.lines.iter_mut().chain(self.on_top.iter_mut()) {
            l.width *= pixels_per_point;
        }
        for p in &mut self.points {
            p.size *= pixels_per_point;
        }
    }
}

/// The lit colour of a surface with normal `n` seen along `view_dir` (toward the eye). The GPU
/// shader computes the same.
pub fn shade(n: Vec3, view_dir: Vec3, base: [f64; 3]) -> [f64; 3] {
    let n = if n.dot(view_dir) < 0.0 { -n } else { n };
    let key = (view_dir * 0.8 + Vec3::new(-0.3, 0.2, 0.6)).normalized().unwrap_or(view_dir);
    let diff = n.dot(key).max(0.0);
    let fill = n.dot(view_dir).max(0.0);
    let h = (key + view_dir).normalized().unwrap_or(view_dir);
    let spec = n.dot(h).max(0.0).powi(40) * 0.25;
    let k = 0.42 + 0.38 * diff + 0.25 * fill;
    [(base[0] * k + spec).min(1.0), (base[1] * k + spec).min(1.0), (base[2] * k + spec).min(1.0)]
}

/// RGBA pixels plus depth.
pub struct Canvas {
    pub w: usize,
    pub h: usize,
    pub rgba: Vec<u8>,
    depth: Vec<f32>,
}

const SS: usize = 2;
/// Largest output side in pixels.
pub const MAX_DIM: usize = 8192;
/// Most output pixels (`w × h`): 16 MP, drawn at 2× into 64 MP of colour and depth (~512 MB).
pub const MAX_PIXELS: usize = 4096 * 4096;
/// Most primitives (triangles, lines and points together) one render draws.
pub const MAX_PRIMITIVES: usize = 20_000_000;

/// Why a render was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RenderError {
    /// The image size is zero or over [`MAX_DIM`] per side or [`MAX_PIXELS`] in all.
    Size { w: usize, h: usize },
    /// The scene has more than [`MAX_PRIMITIVES`] primitives.
    SceneTooLarge { primitives: usize },
    /// The image buffers could not be allocated.
    OutOfMemory,
    /// The PNG encoder failed.
    Encode,
}

impl std::fmt::Display for RenderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RenderError::Size { w, h } => {
                write!(f, "image size {w} x {h} is out of range (1..={MAX_DIM} per side, at most {MAX_PIXELS} pixels)")
            }
            RenderError::SceneTooLarge { primitives } => {
                write!(f, "the scene has {primitives} primitives; the CPU renderer draws at most {MAX_PRIMITIVES}")
            }
            RenderError::OutOfMemory => f.write_str("not enough memory for the image"),
            RenderError::Encode => f.write_str("PNG encoding failed"),
        }
    }
}

impl std::error::Error for RenderError {}

/// Whether `w × h` is an image size [`render`] accepts.
fn size_ok(w: usize, h: usize) -> bool {
    (1..=MAX_DIM).contains(&w) && (1..=MAX_DIM).contains(&h) && w.checked_mul(h).is_some_and(|n| n <= MAX_PIXELS)
}

/// The largest size with the aspect of `w × h` that [`render`] accepts (for drawing a window of
/// any size: the image is scaled up on screen).
pub fn fit_size(w: usize, h: usize) -> (usize, usize) {
    let (w, h) = (w.max(1), h.max(1));
    if size_ok(w, h) {
        return (w, h);
    }
    let by_side = MAX_DIM as f64 / w.max(h) as f64;
    let by_area = (MAX_PIXELS as f64 / (w as f64 * h as f64)).sqrt();
    let k = by_side.min(by_area);
    let fit = |v: usize| ((v as f64 * k).floor() as usize).clamp(1, MAX_DIM);
    (fit(w), fit(h))
}
/// Largest half-width of a line or point, in supersampled pixels.
const MAX_MARK: f64 = 64.0;

/// A triangle of less (twice signed) area than this, in pixels², is degenerate.
const MIN_AREA: f64 = 1e-12; // unit-factors: an area in pixels
/// A pixel centre this far outside an edge (in barycentric weight) is still inside: shared edges leave no gap.
const EDGE_ROUND_OFF: f64 = 1e-9; // unit-factors: a barycentric weight
/// Most labels drawn, and most characters drawn of one.
const MAX_LABELS: usize = 200;
const MAX_LABEL_CHARS: usize = 48;
/// Output pixels per font pixel, and a label's offset from its point (output pixels).
const LABEL_SCALE: usize = 2;
pub const LABEL_OFFSET: f64 = 8.0;
/// Most samples drawn along one line segment (a segment from far off screen is bounded).
const MAX_LINE_STEPS: usize = 20_000;

impl Canvas {
    /// A white canvas; the buffers are reserved fallibly, so a size the machine cannot hold is
    /// an error rather than an abort.
    fn new(w: usize, h: usize) -> Result<Canvas, RenderError> {
        let n = w.checked_mul(h).ok_or(RenderError::Size { w, h })?;
        let bytes = n.checked_mul(4).ok_or(RenderError::Size { w, h })?;
        let mut rgba = Vec::new();
        rgba.try_reserve_exact(bytes).map_err(|_| RenderError::OutOfMemory)?;
        rgba.resize(bytes, 255);
        let mut depth = Vec::new();
        depth.try_reserve_exact(n).map_err(|_| RenderError::OutOfMemory)?;
        depth.resize(n, f32::INFINITY);
        Ok(Canvas { w, h, rgba, depth })
    }
    fn put(&mut self, x: usize, y: usize, c: [f64; 3], a: f64) {
        let i = (y * self.w + x) * 4;
        if let Some(px) = self.rgba.get_mut(i..i + 3) {
            for (k, v) in px.iter_mut().enumerate() {
                let old = *v as f64 / 255.0;
                *v = ((old * (1.0 - a) + c[k] * a).clamp(0.0, 1.0) * 255.0).round() as u8;
            }
        }
    }
}

/// Render a scene with the camera into an RGBA canvas of `w × h` pixels. Sizes outside
/// 1..=[`MAX_DIM`] per side or over [`MAX_PIXELS`], and scenes over [`MAX_PRIMITIVES`], are
/// refused ([`fit_size`] gives an accepted size).
pub fn render(scene: &Scene, cam: &Camera, w: usize, h: usize) -> Result<Canvas, RenderError> {
    render_layers(scene, &[], cam, w, h)
}

/// Depth (of the 0…1 range) an overlay face may lie behind the scene's and still be drawn: it
/// covers the faces it repeats.
const OVERLAY_DEPTH: f32 = 2e-4;
/// The depth bias of lines and points over faces: of the 0…1 range in orthographic views, and
/// times eight of the distance to the eye in perspective (the GPU shader's `screen.z`).
const DEPTH_BIAS: f64 = 2e-4;

/// [`render`] with overlays drawn over the scene (the setup's symbols, a selection highlight):
/// their faces over the scene's faces they coincide with (a highlight repeats the faces it
/// marks), then every layer's translucent faces blended over them, then every layer's lines and points depth-tested, then every layer's
/// [`Scene::on_top`] lines and labels over everything. The overlays' background and radius are
/// not used.
pub fn render_layers(scene: &Scene, overlays: &[&Scene], cam: &Camera, w: usize, h: usize) -> Result<Canvas, RenderError> {
    if !size_ok(w, h) {
        return Err(RenderError::Size { w, h });
    }
    let count = |s: &Scene| {
        s.tris
            .len()
            .saturating_add(s.fills.len())
            .saturating_add(s.lines.len())
            .saturating_add(s.points.len())
            .saturating_add(s.labels.len())
            .saturating_add(s.on_top.len())
    };
    let primitives = overlays.iter().fold(count(scene), |n, o| n.saturating_add(count(o)));
    if primitives > MAX_PRIMITIVES {
        return Err(RenderError::SceneTooLarge { primitives });
    }
    let (sw, sh) = (w * SS, h * SS);
    let mut c = Canvas::new(sw, sh)?;
    let (top, bot) = scene.background.unwrap_or((Rgb(236, 240, 245), Rgb(205, 212, 222)));
    let (a, b) = (top.f(), bot.f());
    for y in 0..sh {
        let t = y as f64 / sh.max(1) as f64;
        let col = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
        for x in 0..sw {
            c.put(x, y, col, 1.0);
        }
    }
    let proj = Projector::new(cam, sw as f64, sh as f64, scene.radius);
    let project = |p: Vec3| proj.screen(p);
    // Edges, nodes and overlay faces win the depth test against the faces they lie on by a bias
    // toward the eye, `s` biases of it. Orthographic depth is linear, so a fixed depth offset
    // is the same everywhere; perspective depth is not (a fixed offset spans a large distance far
    // from the near plane, which let hidden edges show through, #282), so the point is moved
    // toward the eye by a fraction of its distance instead, as the GPU shader does.
    let eye = cam.eye();
    let biased = |p: Vec3, s: f64| -> Option<[f64; 3]> {
        let mut q = project(p)?;
        if cam.fov > 0.0 {
            q[2] = project(p + (eye - p) * (DEPTH_BIAS * 8.0 * s))?[2] + DEPTH_BIAS;
        }
        Some(q)
    };
    let view_dir = cam.back();
    for (t, overlay) in scene.tris.iter().map(|t| (t, false)).chain(overlays.iter().flat_map(|o| &o.tris).map(|t| (t, true))) {
        let corner = |p: Vec3| if overlay { biased(p, 0.5) } else { project(p) };
        let (Some(a), Some(b), Some(cc)) = (corner(t.p[0]), corner(t.p[1]), corner(t.p[2])) else { continue };
        let allowance = if overlay { OVERLAY_DEPTH } else { 0.0 };
        let col = shade(t.normal, view_dir, t.color.f());
        raster_tri(a, b, cc, sw, sh, |x, y, z| {
            if let Some(d) = c.depth.get_mut(y * sw + x)
                && z < *d + allowance
            {
                *d = d.min(z);
                c.put(x, y, col, 1.0);
            }
        });
    }
    let layers = || std::iter::once(scene).chain(overlays.iter().copied());
    // Translucent faces over the opaque ones, without writing depth: what is behind shows.
    for f in layers().flat_map(|s| &s.fills) {
        let t = &f.tri;
        let (Some(a), Some(b), Some(cc)) = (project(t.p[0]), project(t.p[1]), project(t.p[2])) else { continue };
        let alpha = if f.alpha.is_finite() { f.alpha.clamp(0.0, 1.0) } else { 0.0 };
        let col = shade(t.normal, view_dir, t.color.f());
        raster_tri(a, b, cc, sw, sh, |x, y, z| {
            if c.depth.get(y * sw + x).is_some_and(|d| z < *d) {
                c.put(x, y, col, alpha);
            }
        });
    }
    let line = |c: &mut Canvas, l: &SceneLine, depth_test: bool| {
        let (Some(a), Some(b)) = (biased(l.a, 1.0), biased(l.b, 1.0)) else { return };
        let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
        let steps = if len.is_finite() { (len.ceil() as usize).clamp(1, MAX_LINE_STEPS) } else { return };
        let col = l.color.f();
        for k in 0..=steps {
            let t = k as f64 / steps as f64;
            let p = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
            dot(c, p, l.width * SS as f64, col, true, depth_test);
        }
    };
    for l in layers().flat_map(|s| &s.lines) {
        line(&mut c, l, true);
    }
    for p in layers().flat_map(|s| &s.points) {
        if let Some(q) = biased(p.p, 1.0) {
            dot(&mut c, q, p.size * SS as f64, p.color.f(), false, true);
        }
    }
    for l in layers().flat_map(|s| &s.on_top) {
        line(&mut c, l, false);
    }
    for l in layers().flat_map(|s| &s.labels).take(MAX_LABELS) {
        if let Some(q) = project(l.p) {
            label(&mut c, q, &l.text, l.color.f());
        }
    }
    downsample(&c, w, h)
}

/// `layers`' labels drawn over finished `w × h` RGBA pixels (a GPU picture), as [`render_layers`]
/// draws them: the pixels are drawn at the supersampled size, so pixels away from a label come
/// back unchanged. `radius` is the scene's clipping radius.
pub fn draw_labels(rgba: &[u8], w: usize, h: usize, layers: &[&Scene], cam: &Camera, radius: f64) -> Result<Vec<u8>, RenderError> {
    if !size_ok(w, h) || rgba.len() != w * h * 4 {
        return Err(RenderError::Size { w, h });
    }
    let labels: Vec<&SceneLabel> = layers.iter().flat_map(|s| &s.labels).take(MAX_LABELS).collect();
    if labels.is_empty() {
        return Ok(rgba.to_vec());
    }
    let (sw, sh) = (w * SS, h * SS);
    let mut c = Canvas::new(sw, sh)?;
    for y in 0..sh {
        for x in 0..sw {
            let (i, o) = (((y / SS) * w + x / SS) * 4, (y * sw + x) * 4);
            if let (Some(src), Some(dst)) = (rgba.get(i..i + 4), c.rgba.get_mut(o..o + 4)) {
                dst.copy_from_slice(src);
            }
        }
    }
    let proj = Projector::new(cam, sw as f64, sh as f64, radius);
    for l in labels {
        if let Some(q) = proj.screen(l.p) {
            label(&mut c, q, &l.text, l.color.f());
        }
    }
    Ok(downsample(&c, w, h)?.rgba)
}

/// Cover a triangle with screen corners `a`, `b`, `c` (x, y, depth) on a `w × h` pixel grid:
/// `pixel(x, y, depth)` for every pixel whose centre is inside it (its edges included), the depth
/// interpolated. Degenerate and off-grid triangles cover nothing. The renderer and the picking
/// depth map ([`crate::pick::DepthMap`]) share it, so a box selects what is drawn.
pub(crate) fn raster_tri(a: [f64; 3], b: [f64; 3], c: [f64; 3], w: usize, h: usize, mut pixel: impl FnMut(usize, usize, f32)) {
    let (x0, x1) = (a[0].min(b[0]).min(c[0]).floor().max(0.0), a[0].max(b[0]).max(c[0]).ceil().min(w as f64 - 1.0));
    let (y0, y1) = (a[1].min(b[1]).min(c[1]).floor().max(0.0), a[1].max(b[1]).max(c[1]).ceil().min(h as f64 - 1.0));
    if !(x0 <= x1 && y0 <= y1) {
        return;
    }
    let area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if !(area.abs() >= MIN_AREA) {
        return;
    }
    for y in y0 as usize..=y1 as usize {
        for x in x0 as usize..=x1 as usize {
            let (px, py) = (x as f64 + 0.5, y as f64 + 0.5);
            let w0 = ((b[0] - px) * (c[1] - py) - (b[1] - py) * (c[0] - px)) / area;
            let w1 = ((c[0] - px) * (a[1] - py) - (c[1] - py) * (a[0] - px)) / area;
            let w2 = 1.0 - w0 - w1;
            if w0 < -EDGE_ROUND_OFF || w1 < -EDGE_ROUND_OFF || w2 < -EDGE_ROUND_OFF {
                continue;
            }
            pixel(x, y, (w0 * a[2] + w1 * b[2] + w2 * c[2]) as f32);
        }
    }
}

/// Text above and right of screen point `p`, over everything, with a halo of the opposite
/// lightness so it reads on any background.
fn label(c: &mut Canvas, p: [f64; 3], text: &str, col: [f64; 3]) {
    if !(p[0].is_finite() && p[1].is_finite()) {
        return;
    }
    let px = LABEL_SCALE * SS;
    let (x0, y0) = (p[0] + LABEL_OFFSET * SS as f64, p[1] - LABEL_OFFSET * SS as f64 - (font::H * px) as f64);
    let (sw, sh) = (c.w as f64, c.h as f64);
    if x0 > sw || y0 > sh || x0 < -sw || y0 < -sh {
        return;
    }
    let luma = 0.3 * col[0] + 0.6 * col[1] + 0.1 * col[2];
    let halo = if luma > 0.5 { [0.0; 3] } else { [1.0; 3] };
    let (x0, y0) = (x0.floor() as i64, y0.floor() as i64);
    let chars: Vec<char> = text.chars().take(MAX_LABEL_CHARS).collect();
    // The halo, spread a supersampled pixel pair around each font pixel, then the text.
    for (grow, color, a) in [(SS as i64, halo, 0.7), (0, col, 1.0)] {
        for (i, ch) in chars.iter().enumerate() {
            for fy in 0..font::H {
                for fx in (0..font::W).filter(|&fx| font::on(*ch, fx, fy)) {
                    let (bx, by) = (x0 + ((i * font::ADVANCE + fx) * px) as i64, y0 + (fy * px) as i64);
                    for y in by - grow..by + px as i64 + grow {
                        for x in bx - grow..bx + px as i64 + grow {
                            if x >= 0 && y >= 0 && (x as usize) < c.w && (y as usize) < c.h {
                                c.put(x as usize, y as usize, color, a);
                            }
                        }
                    }
                }
            }
        }
    }
}

/// A disc (`round`) or square of diameter `size` at screen point `p` (x, y, depth), depth-tested
/// unless `depth_test` is off.
fn dot(c: &mut Canvas, p: [f64; 3], size: f64, col: [f64; 3], round: bool, depth_test: bool) {
    // Widths are capped so a hostile width cannot make the loop below unbounded.
    let half = if size.is_finite() { (size * 0.5).clamp(0.5, MAX_MARK) } else { 0.5 };
    let r = half.ceil() as i64;
    let (sw, sh) = (c.w as i64, c.h as i64);
    if !(p[0].is_finite() && p[1].is_finite()) || p[0] < -half || p[1] < -half || p[0] > sw as f64 + half || p[1] > sh as f64 + half {
        return;
    }
    for dy in -r..=r {
        for dx in -r..=r {
            let (ix, iy) = (p[0].floor() as i64 + dx, p[1].floor() as i64 + dy);
            if ix < 0 || iy < 0 || ix >= sw || iy >= sh {
                continue;
            }
            let (ox, oy) = (ix as f64 + 0.5 - p[0], iy as f64 + 0.5 - p[1]);
            let d = if round { (ox * ox + oy * oy).sqrt() } else { ox.abs().max(oy.abs()) };
            let cov = (half + 0.5 - d).clamp(0.0, 1.0);
            let i = iy as usize * c.w + ix as usize;
            // Lines and points sit on faces: a small depth allowance keeps them in front.
            if cov > 0.0 && (!depth_test || c.depth.get(i).is_some_and(|dd| (p[2] as f32) <= *dd + 2e-4)) {
                c.put(ix as usize, iy as usize, col, cov);
            }
        }
    }
}

fn downsample(c: &Canvas, w: usize, h: usize) -> Result<Canvas, RenderError> {
    let mut out = Canvas::new(w, h)?;
    for y in 0..h {
        for x in 0..w {
            let mut acc = [0u32; 4];
            for sy in 0..SS {
                for sx in 0..SS {
                    let i = ((y * SS + sy) * c.w + x * SS + sx) * 4;
                    for (k, a) in acc.iter_mut().enumerate() {
                        *a += c.rgba.get(i + k).copied().unwrap_or(0) as u32;
                    }
                }
            }
            let o = (y * w + x) * 4;
            for (k, a) in acc.iter().enumerate() {
                if let Some(v) = out.rgba.get_mut(o + k) {
                    *v = (*a / (SS * SS) as u32) as u8;
                }
            }
        }
    }
    Ok(out)
}

/// Render to PNG bytes.
pub fn render_png(scene: &Scene, cam: &Camera, w: usize, h: usize) -> Result<Vec<u8>, RenderError> {
    render_layers_png(scene, &[], cam, w, h)
}

/// [`render_layers`] as PNG bytes.
pub fn render_layers_png(scene: &Scene, overlays: &[&Scene], cam: &Camera, w: usize, h: usize) -> Result<Vec<u8>, RenderError> {
    let c = render_layers(scene, overlays, cam, w, h)?;
    png(c.w, c.h, &c.rgba).ok_or(RenderError::Encode)
}

/// RGBA pixels (`w × h`, 8 bits per channel) to PNG bytes; `None` when the sizes disagree.
pub fn png(w: usize, h: usize, rgba: &[u8]) -> Option<Vec<u8>> {
    if w.checked_mul(h)?.checked_mul(4)? != rgba.len() || w == 0 || h == 0 {
        return None;
    }
    let mut buf = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut buf, u32::try_from(w).ok()?, u32::try_from(h).ok()?);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header().ok()?;
        writer.write_image_data(rgba).ok()?;
        writer.finish().ok()?;
    }
    Some(buf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StandardView;

    fn red_triangle() -> Scene {
        Scene {
            tris: vec![SceneTri {
                p: [Vec3::new(-10.0, 0.0, -10.0), Vec3::new(10.0, 0.0, -10.0), Vec3::new(0.0, 0.0, 10.0)],
                normal: Vec3::new(0.0, -1.0, 0.0),
                color: Rgb(200, 30, 30),
            }],
            radius: 20.0,
            ..Default::default()
        }
    }

    #[test]
    fn renders_a_triangle_lines_and_points() {
        let mut scene = red_triangle();
        let mut cam = Camera { distance: 50.0, ..Default::default() };
        cam.set_view(StandardView::Front);
        let c = render(&scene, &cam, 64, 64).unwrap();
        let px = |c: &Canvas, x: usize, y: usize| c.rgba[(y * c.w + x) * 4..(y * c.w + x) * 4 + 4].to_vec();
        let centre = px(&c, 32, 32);
        assert!(centre[0] > centre[1] + 50, "centre pixel is red: {centre:?}");
        // A black line across the middle, in front of the triangle, and a point.
        scene.lines.push(SceneLine { a: Vec3::new(-10.0, -0.1, 0.0), b: Vec3::new(10.0, -0.1, 0.0), color: Rgb(0, 0, 0), width: 3.0 });
        scene.points.push(ScenePoint { p: Vec3::new(0.0, -0.1, 5.0), color: Rgb(0, 0, 255), size: 5.0 });
        let c = render(&scene, &cam, 64, 64).unwrap();
        assert!(px(&c, 32, 32)[0] < 40, "{:?}", px(&c, 32, 32));
        let (x, y, _) = cam.to_screen(Vec3::new(0.0, 0.0, 5.0), 64.0, 64.0, 20.0).unwrap();
        let p = px(&c, x as usize, y as usize);
        assert!(p[2] > 200 && p[0] < 60, "{p:?}");
    }

    /// A translucent face tints what is behind it and hides nothing; one behind an opaque face
    /// does not show (#302).
    #[test]
    fn translucent_faces_blend_over_what_is_in_front_only() {
        let mut cam = Camera { distance: 50.0, ..Default::default() };
        cam.set_view(StandardView::Front);
        let px = |c: &Canvas| c.rgba[(32 * c.w + 32) * 4..(32 * c.w + 32) * 4 + 3].to_vec();
        let red = red_triangle();
        let plain = px(&render(&red, &cam, 64, 64).unwrap());
        let fill = |y: f64| SceneFill {
            tri: SceneTri { p: red.tris[0].p.map(|p| Vec3::new(p.x, y, p.z)), color: Rgb(0, 0, 255), ..red.tris[0] },
            alpha: 0.25,
        };
        let mut front = red.clone();
        front.fills.push(fill(-1.0));
        let tinted = px(&render(&front, &cam, 64, 64).unwrap());
        assert!(tinted[2] > plain[2] + 20 && tinted[0] > 120, "red still shows, tinted blue: {tinted:?} over {plain:?}");
        let mut back = red.clone();
        back.fills.push(fill(1.0));
        assert_eq!(px(&render(&back, &cam, 64, 64).unwrap()), plain);
    }

    /// In perspective, close to the part (the near plane at its smallest), a line a little
    /// behind a face stays hidden and one on it shows (#282).
    #[test]
    fn perspective_hides_lines_behind_faces() {
        let mut scene = red_triangle();
        scene.radius = 40.0;
        let mut cam = Camera { distance: 30.0, fov: 0.6, ..Default::default() };
        cam.set_view(StandardView::Front);
        let px = |c: &Canvas| c.rgba[(32 * c.w + 32) * 4..(32 * c.w + 32) * 4 + 3].to_vec();
        // Front view looks along +y: y = 1 is behind the face, 3 % of the distance.
        scene.lines.push(SceneLine { a: Vec3::new(-10.0, 1.0, 0.0), b: Vec3::new(10.0, 1.0, 0.0), color: Rgb(0, 0, 0), width: 3.0 });
        let behind = render(&scene, &cam, 64, 64).unwrap();
        assert!(px(&behind)[0] > 150, "hidden: {:?}", px(&behind));
        scene.lines[0].a.y = 0.0;
        scene.lines[0].b.y = 0.0;
        let on = render(&scene, &cam, 64, 64).unwrap();
        assert!(px(&on)[0] < 60, "on the face: {:?}", px(&on));
    }

    /// An overlay face repeating a scene face covers it; one behind another face does not show.
    #[test]
    fn an_overlay_covers_the_faces_it_repeats() {
        let scene = red_triangle();
        let mut cam = Camera { distance: 50.0, ..Default::default() };
        cam.set_view(StandardView::Front);
        let px = |c: &Canvas| c.rgba[(32 * c.w + 32) * 4..(32 * c.w + 32) * 4 + 3].to_vec();
        let green = |s: &Scene| Scene { tris: s.tris.iter().map(|t| SceneTri { color: Rgb(0, 255, 0), ..*t }).collect(), ..Scene::default() };
        let over = render_layers(&scene, &[&green(&scene)], &cam, 64, 64).unwrap();
        assert!(px(&over)[1] > px(&over)[0] + 50, "{:?}", px(&over));
        // The same face moved behind the red one (front view looks along +y) stays hidden.
        let mut behind = green(&scene);
        for t in &mut behind.tris {
            t.p = t.p.map(|p| p + Vec3::Y * 5.0);
        }
        let hidden = render_layers(&scene, &[&behind], &cam, 64, 64).unwrap();
        assert_eq!(px(&hidden), px(&render(&scene, &cam, 64, 64).unwrap()));
    }

    /// A line behind a face is hidden; the same line on top is drawn over it.
    #[test]
    fn on_top_lines_draw_over_the_faces_in_front_of_them() {
        let scene = red_triangle();
        let mut cam = Camera { distance: 50.0, ..Default::default() };
        cam.set_view(StandardView::Front);
        let px = |c: &Canvas| c.rgba[(32 * c.w + 32) * 4..(32 * c.w + 32) * 4 + 3].to_vec();
        // Behind the triangle (the front view looks along +y).
        let line = SceneLine { a: Vec3::new(-10.0, 5.0, 0.0), b: Vec3::new(10.0, 5.0, 0.0), color: Rgb(0, 0, 255), width: 4.0 };
        let behind = Scene { lines: vec![line], ..Scene::default() };
        let hidden = render_layers(&scene, &[&behind], &cam, 64, 64).unwrap();
        assert_eq!(px(&hidden), px(&render(&scene, &cam, 64, 64).unwrap()));
        let over = Scene { on_top: vec![line], ..Scene::default() };
        let shown = render_layers(&scene, &[&over], &cam, 64, 64).unwrap();
        assert!(px(&shown)[2] > 200 && px(&shown)[0] < 60, "{:?}", px(&shown));
    }

    #[test]
    fn png_bytes_and_hostile_sizes() {
        let scene = red_triangle();
        let cam = Camera::default();
        assert!(render_png(&scene, &cam, 16, 16).unwrap().starts_with(&[0x89, b'P', b'N', b'G']));
        // Out-of-range sizes are refused before anything is allocated.
        for (w, h) in [(0, 16), (16, 0), (0, usize::MAX), (MAX_DIM + 1, 16), (100_000, 100_000), (usize::MAX, usize::MAX), (8192, 8192)] {
            assert_eq!(render(&scene, &cam, w, h).err(), Some(RenderError::Size { w, h }));
        }
        assert!(render(&scene, &cam, MAX_DIM, MAX_PIXELS / MAX_DIM).is_ok_and(|c| c.w == MAX_DIM));
        assert!(png(2, 2, &[0; 3]).is_none());
        assert!(png(usize::MAX, 2, &[]).is_none());
        // Far-away and non-finite geometry is skipped or bounded, never a panic.
        let wild = Scene {
            lines: vec![SceneLine { a: Vec3::new(-1e12, 0.0, 0.0), b: Vec3::new(1e12, 0.0, 0.0), color: Rgb(0, 0, 0), width: 1e6 }],
            points: vec![ScenePoint { p: Vec3::new(f64::NAN, 0.0, 0.0), color: Rgb(0, 0, 0), size: 4.0 }],
            radius: f64::NAN,
            ..scene
        };
        let _ = render(&wild, &cam, 32, 32);
    }

    #[test]
    fn labels_draw_over_finished_pixels() {
        let (w, h) = (64, 32);
        let px: Vec<u8> = (0..w * h).flat_map(|i| [(i % 251) as u8, 40, 200, 255]).collect();
        let cam = Camera::default();
        let mut scene = Scene::default();
        assert_eq!(draw_labels(&px, w, h, &[&scene], &cam, 10.0).unwrap(), px);
        scene.labels.push(SceneLabel { p: Vec3::ZERO, text: "RP".into(), color: Rgb(0, 0, 0) });
        let out = draw_labels(&px, w, h, &[&scene], &cam, 10.0).unwrap();
        let changed = out.chunks(4).zip(px.chunks(4)).filter(|(a, b)| a != b).count();
        assert!(changed > 10 && changed < w * h / 2, "{changed}");
        assert!(draw_labels(&px, w + 1, h, &[&scene], &cam, 10.0).is_err());
    }

    #[test]
    fn labels_draw_over_everything_beside_their_point() {
        let mut scene = red_triangle();
        let mut cam = Camera { distance: 50.0, ..Default::default() };
        cam.set_view(StandardView::Front);
        let plain = render(&scene, &cam, 64, 64).unwrap();
        scene.labels.push(SceneLabel { p: Vec3::new(0.0, 0.0, 0.0), text: "RP".into(), color: Rgb(0, 0, 0) });
        let c = render(&scene, &cam, 64, 64).unwrap();
        // Up and right of the centre, some pixels turned dark; left of it nothing changed.
        let dark = |c: &Canvas, x0: usize, x1: usize| {
            (x0..x1).flat_map(|x| (0..32).map(move |y| (x, y))).filter(|&(x, y)| c.rgba[(y * 64 + x) * 4] < 40).count()
        };
        assert!(dark(&c, 32, 64) > 10 && dark(&plain, 32, 64) == 0, "{}", dark(&c, 32, 64));
        assert_eq!(dark(&c, 0, 32), 0);
        // Labels far off screen, empty or very long, or at a non-finite point, are bounded.
        scene.labels = vec![
            SceneLabel { p: Vec3::new(1e9, 0.0, 0.0), text: "far".into(), color: Rgb(0, 0, 0) },
            SceneLabel { p: Vec3::new(f64::NAN, 0.0, 0.0), text: "nan".into(), color: Rgb(0, 0, 0) },
            SceneLabel { p: Vec3::ZERO, text: String::new(), color: Rgb(255, 255, 255) },
            SceneLabel { p: Vec3::ZERO, text: "x".repeat(100_000), color: Rgb(255, 255, 255) },
        ];
        assert!(render(&scene, &cam, 64, 64).is_ok());
    }

    #[test]
    fn huge_scenes_are_refused() {
        let line = SceneLine { a: Vec3::ZERO, b: Vec3::X, color: Rgb(0, 0, 0), width: 1.0 };
        let scene = Scene { lines: vec![line; MAX_PRIMITIVES + 1], radius: 1.0, ..Default::default() };
        let r = render_png(&scene, &Camera::default(), 16, 16);
        assert_eq!(r, Err(RenderError::SceneTooLarge { primitives: MAX_PRIMITIVES + 1 }));
    }

    #[test]
    fn fit_size_keeps_the_aspect_within_the_limits() {
        assert_eq!(fit_size(1280, 800), (1280, 800));
        assert_eq!(fit_size(0, 0), (1, 1));
        for (w, h) in [(7680, 4320), (100_000, 100_000), (usize::MAX, 3), (20_000, 10)] {
            let (fw, fh) = fit_size(w, h);
            assert!(size_ok(fw, fh), "{w}x{h} -> {fw}x{fh}");
        }
        let (fw, fh) = fit_size(7680, 4320);
        assert!(((fw as f64 / fh as f64) - 16.0 / 9.0).abs() < 1e-2, "{fw}x{fh}");
    }
}
