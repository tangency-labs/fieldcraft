//! The whole window drawn without a window: the app's UI pass run in an offscreen egui context,
//! tessellated, and rasterised on the CPU (panels, toolbar, dialogs, and the viewport through the
//! CPU model renderer).
//!
//! `ui.screenshot` needs a presented frame, and a hidden window or a locked screen (macOS marks
//! the window occluded) presents none: eframe then skips the UI pass. `ui.capture` works in any
//! of those states, and headless tests use it to look at the real UI.
//!
//! The offscreen context keeps its own textures (font atlas, icons, the viewport image) and
//! memory between captures. Its rasteriser follows egui's painting contract: premultiplied sRGBA
//! vertex colours multiplied by the texel, blended "over" in gamma space, triangles clipped to
//! their clip rectangle; egui's feathering supplies the anti-aliasing.

use std::collections::HashMap;

use egui::{Color32, ColorImage, Pos2, Rect, TextureId, ViewportId, ViewportInfo};

use crate::icons::IconCache;
use crate::theme::Tokens;
use crate::{ModelApp, Viewport};

/// Largest capture side in physical pixels.
pub const MAX_SIDE: usize = 8192;
/// Most physical pixels in one capture.
pub const MAX_PIXELS: usize = 4096 * 4096;
/// Pixels per point a capture may ask for.
pub const SCALE_RANGE: std::ops::RangeInclusive<f32> = 0.5..=4.0;
/// Frames run before the one rasterised: egui sizes some panels and windows in a first pass.
const WARMUP_FRAMES: usize = 2;
/// A triangle of less (twice signed) area than this, in pixels², is degenerate.
const MIN_AREA: f32 = 1e-6; // unit-factors: an area in pixels

/// A texture the offscreen context uploaded.
struct Texture {
    size: [usize; 2],
    pixels: Vec<Color32>,
    linear: bool,
}

/// The offscreen context and what belongs to it (textures, icon textures, the viewport image),
/// kept between captures.
#[derive(Default)]
pub struct Offscreen {
    ctx: egui::Context,
    textures: HashMap<TextureId, Texture>,
    icons: IconCache,
    viewport: Viewport,
}

/// A capture's size: points and pixels per point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaptureSize {
    pub width: f32,
    pub height: f32,
    pub scale: f32,
}

impl CaptureSize {
    /// The image size in physical pixels, or why it is refused.
    pub fn pixels(&self) -> Result<[usize; 2], String> {
        if !SCALE_RANGE.contains(&self.scale) {
            return Err(format!("`scale` must be from {} to {}", SCALE_RANGE.start(), SCALE_RANGE.end()));
        }
        let px = |v: f32| (v * self.scale).round();
        let (w, h) = (px(self.width), px(self.height));
        let ok = |v: f32| v.is_finite() && (1.0..=MAX_SIDE as f32).contains(&v);
        if !(ok(w) && ok(h)) || (w as usize).saturating_mul(h as usize) > MAX_PIXELS {
            return Err(format!(
                "capture size {} x {} points at scale {} is out of range (at most {MAX_SIDE} pixels a side, {MAX_PIXELS} in all)",
                self.width, self.height, self.scale
            ));
        }
        Ok([w as usize, h as usize])
    }
}

impl ModelApp {
    /// The whole window as it would be drawn at `size`, rendered offscreen.
    pub fn capture(&mut self, size: CaptureSize) -> Result<ColorImage, String> {
        let [w, h] = size.pixels()?;
        let mut off = self.offscreen.take().unwrap_or_default();
        // The UI pass draws with the offscreen context's own icon textures and a viewport with no
        // GPU target, which renders the model on the CPU into an offscreen texture.
        std::mem::swap(&mut self.icons, &mut off.icons);
        std::mem::swap(&mut self.viewport, &mut off.viewport);
        let toolbar = std::mem::take(&mut self.toolbar_rects);
        let layout = std::mem::take(&mut self.layout);
        let panels = std::mem::take(&mut self.panel_rects);
        self.capturing = true;
        crate::theme::apply(&off.ctx, self.ui.dark);
        // A capture is one still frame: windows and menus appear without fading in.
        off.ctx.global_style_mut(|s| s.animation_time = 0.0);
        let mut shapes = Vec::new();
        let mut ppp = size.scale;
        // Textures freed by the frame painted, released after it is painted.
        let mut free = Vec::new();
        for _ in 0..WARMUP_FRAMES {
            let mut raw = egui::RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, egui::vec2(size.width, size.height))),
                time: Some(self.now),
                ..Default::default()
            };
            raw.viewports.insert(ViewportId::ROOT, ViewportInfo { native_pixels_per_point: Some(size.scale), ..Default::default() });
            let mut out = off.ctx.run_ui(raw, |ui| self.ui(ui));
            off.free(free.drain(..));
            let mut delta = std::mem::take(&mut out.textures_delta);
            for (id, deltas) in &delta.set {
                for d in deltas {
                    off.set_texture(*id, d);
                }
            }
            free.extend(delta.free.iter().copied());
            // Every delta is applied (egui checks that none is dropped unhandled).
            delta.clear();
            shapes = out.shapes;
            ppp = out.pixels_per_point;
        }
        std::mem::swap(&mut self.icons, &mut off.icons);
        std::mem::swap(&mut self.viewport, &mut off.viewport);
        self.capturing = false;
        self.toolbar_rects = toolbar;
        self.layout = layout;
        self.panel_rects = panels;
        let primitives = off.ctx.tessellate(shapes, ppp);
        let mut canvas = Canvas::new(w, h, Tokens::get(self.ui.dark).panel)?;
        for p in &primitives {
            if let egui::epaint::Primitive::Mesh(mesh) = &p.primitive {
                canvas.mesh(mesh, p.clip_rect, ppp, off.textures.get(&mesh.texture_id));
            }
        }
        off.free(free.drain(..));
        self.offscreen = Some(off);
        Ok(canvas.into_image())
    }
}

impl Offscreen {
    fn set_texture(&mut self, id: TextureId, d: &egui::epaint::ImageDelta) {
        let egui::ImageData::Color(image) = &d.image;
        let linear = d.options.magnification == egui::TextureFilter::Linear;
        match d.pos {
            None => {
                self.textures.insert(id, Texture { size: image.size, pixels: image.pixels.clone(), linear });
            }
            Some([x0, y0]) => {
                let Some(t) = self.textures.get_mut(&id) else { return };
                let [pw, ph] = image.size;
                for y in 0..ph {
                    for x in 0..pw {
                        let (tx, ty) = (x0 + x, y0 + y);
                        if tx < t.size[0]
                            && ty < t.size[1]
                            && let (Some(dst), Some(src)) = (t.pixels.get_mut(ty * t.size[0] + tx), image.pixels.get(y * pw + x))
                        {
                            *dst = *src;
                        }
                    }
                }
            }
        }
    }

    fn free(&mut self, ids: impl IntoIterator<Item = TextureId>) {
        for id in ids {
            self.textures.remove(&id);
        }
    }
}

/// Premultiplied RGBA in gamma space, 0–255 per channel.
type Rgba = [f32; 4];

fn rgba(c: Color32) -> Rgba {
    [f32::from(c.r()), f32::from(c.g()), f32::from(c.b()), f32::from(c.a())]
}

impl Texture {
    /// The texel at `uv` (clamped to the edge), bilinear for a linear texture.
    fn sample(&self, u: f32, v: f32) -> Rgba {
        let [w, h] = self.size;
        let at = |x: usize, y: usize| self.pixels.get(y.min(h.saturating_sub(1)) * w + x.min(w.saturating_sub(1))).copied().map_or([0.0; 4], rgba);
        let (x, y) = (u * w as f32 - 0.5, v * h as f32 - 0.5);
        if !self.linear || !(x.is_finite() && y.is_finite()) {
            let (x, y) = ((u * w as f32).max(0.0) as usize, (v * h as f32).max(0.0) as usize);
            return at(x, y);
        }
        let (x0, y0) = (x.floor(), y.floor());
        let (fx, fy) = (x - x0, y - y0);
        let (xa, ya) = (x0.max(0.0) as usize, y0.max(0.0) as usize);
        let (xb, yb) = ((x0 + 1.0).max(0.0) as usize, (y0 + 1.0).max(0.0) as usize);
        let (a, b, c, d) = (at(xa, ya), at(xb, ya), at(xa, yb), at(xb, yb));
        std::array::from_fn(|k| {
            let top = a[k] + (b[k] - a[k]) * fx;
            let bottom = c[k] + (d[k] - c[k]) * fx;
            top + (bottom - top) * fy
        })
    }
}

/// The capture's pixels.
struct Canvas {
    w: usize,
    h: usize,
    px: Vec<Rgba>,
}

impl Canvas {
    fn new(w: usize, h: usize, background: Color32) -> Result<Canvas, String> {
        let n = w.checked_mul(h).ok_or("capture too large")?;
        let mut px = Vec::new();
        px.try_reserve_exact(n).map_err(|_| "not enough memory for the capture".to_string())?;
        px.resize(n, rgba(background));
        Ok(Canvas { w, h, px })
    }

    /// Rasterise a mesh (positions in points) inside `clip` (points).
    fn mesh(&mut self, mesh: &egui::Mesh, clip: Rect, ppp: f32, texture: Option<&Texture>) {
        // Clip rectangles snap to whole pixels, as a GPU scissor rectangle does.
        let cx0 = (clip.min.x * ppp).round().clamp(0.0, self.w as f32);
        let cy0 = (clip.min.y * ppp).round().clamp(0.0, self.h as f32);
        let cx1 = (clip.max.x * ppp).round().clamp(0.0, self.w as f32);
        let cy1 = (clip.max.y * ppp).round().clamp(0.0, self.h as f32);
        if !(cx0 < cx1 && cy0 < cy1) {
            return;
        }
        for &[i, j, k] in mesh.indices.as_chunks::<3>().0 {
            let v = |i: u32| mesh.vertices.get(i as usize);
            let (Some(a), Some(b), Some(c)) = (v(i), v(j), v(k)) else { continue };
            self.triangle([a, b, c], [cx0, cy0, cx1, cy1], ppp, texture);
        }
    }

    fn triangle(&mut self, v: [&egui::epaint::Vertex; 3], clip: [f32; 4], ppp: f32, texture: Option<&Texture>) {
        let p = v.map(|v| [v.pos.x * ppp, v.pos.y * ppp]);
        let [a, b, c] = p;
        let area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if !(area.abs() > MIN_AREA) {
            return;
        }
        let x0 = a[0].min(b[0]).min(c[0]).floor().max(clip[0]);
        let x1 = a[0].max(b[0]).max(c[0]).ceil().min(clip[2]);
        let y0 = a[1].min(b[1]).min(c[1]).floor().max(clip[1]);
        let y1 = a[1].max(b[1]).max(c[1]).ceil().min(clip[3]);
        if !(x0 < x1 && y0 < y1) {
            return;
        }
        // Edge functions oriented so inside is positive; the top-left rule gives a pixel centre
        // on an edge shared by two triangles to exactly one of them.
        let s = area.signum();
        let edge = |p: [f32; 2], q: [f32; 2], x: f32, y: f32| s * ((q[0] - p[0]) * (y - p[1]) - (q[1] - p[1]) * (x - p[0]));
        let top_left = |p: [f32; 2], q: [f32; 2]| {
            let (dx, dy) = (s * (q[0] - p[0]), s * (q[1] - p[1]));
            (dy == 0.0 && dx < 0.0) || dy > 0.0
        };
        let bias = [top_left(b, c), top_left(c, a), top_left(a, b)];
        let col = v.map(|v| rgba(v.color));
        let uv = v.map(|v| [v.uv.x, v.uv.y]);
        let inv = 1.0 / area.abs();
        for y in y0 as usize..y1 as usize {
            for x in x0 as usize..x1 as usize {
                let (px, py) = (x as f32 + 0.5, y as f32 + 0.5);
                let w = [edge(b, c, px, py), edge(c, a, px, py), edge(a, b, px, py)];
                if w.iter().zip(bias).any(|(w, tl)| *w < 0.0 || (*w == 0.0 && !tl)) {
                    continue;
                }
                let l = w.map(|w| w * inv);
                let mut src: Rgba = std::array::from_fn(|k| l[0] * col[0][k] + l[1] * col[1][k] + l[2] * col[2][k]);
                if let Some(t) = texture {
                    let (u, vv) = (l[0] * uv[0][0] + l[1] * uv[1][0] + l[2] * uv[2][0], l[0] * uv[0][1] + l[1] * uv[1][1] + l[2] * uv[2][1]);
                    let texel = t.sample(u, vv);
                    for k in 0..4 {
                        src[k] *= texel[k] / 255.0;
                    }
                }
                if let Some(d) = self.px.get_mut(y * self.w + x) {
                    let keep = 1.0 - src[3] / 255.0;
                    for k in 0..4 {
                        d[k] = src[k] + d[k] * keep;
                    }
                }
            }
        }
    }

    fn into_image(self) -> ColorImage {
        let pixels = self
            .px
            .iter()
            .map(|p| {
                Color32::from_rgb(p[0].round().clamp(0.0, 255.0) as u8, p[1].round().clamp(0.0, 255.0) as u8, p[2].round().clamp(0.0, 255.0) as u8)
            })
            .collect();
        ColorImage { size: [self.w, self.h], source_size: egui::vec2(self.w as f32, self.h as f32), pixels }
    }
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;

    use super::*;

    fn size(width: f32, height: f32) -> CaptureSize {
        CaptureSize { width, height, scale: 1.0 }
    }

    fn pixel(image: &ColorImage, x: usize, y: usize) -> Color32 {
        image.pixels[y * image.size[0] + x]
    }

    #[test]
    fn a_textured_mesh_is_clipped_and_blended() {
        let mut c = Canvas::new(8, 8, Color32::BLACK).unwrap();
        let mut mesh = egui::Mesh::default();
        mesh.add_colored_rect(Rect::from_min_max(Pos2::ZERO, egui::pos2(8.0, 8.0)), Color32::from_rgba_premultiplied(100, 0, 0, 128));
        c.mesh(&mesh, Rect::from_min_max(Pos2::ZERO, egui::pos2(4.0, 8.0)), 1.0, None);
        let image = c.into_image();
        // Half-transparent red over black inside the clip rectangle, black outside.
        assert_eq!(pixel(&image, 1, 1), Color32::from_rgb(100, 0, 0));
        assert_eq!(pixel(&image, 6, 1), Color32::BLACK);
        // The two triangles of the rectangle meet on its diagonal without blending twice.
        assert_eq!(pixel(&image, 2, 2), Color32::from_rgb(100, 0, 0));
        let mut c = Canvas::new(4, 4, Color32::BLACK).unwrap();
        let t = Texture { size: [1, 1], pixels: vec![Color32::from_rgb(0, 200, 0)], linear: true };
        let mut mesh = egui::Mesh::with_texture(TextureId::default());
        mesh.add_rect_with_uv(
            Rect::from_min_max(Pos2::ZERO, egui::pos2(4.0, 4.0)),
            Rect::from_min_max(Pos2::ZERO, egui::pos2(1.0, 1.0)),
            Color32::WHITE,
        );
        c.mesh(&mesh, Rect::EVERYTHING, 1.0, Some(&t));
        assert_eq!(pixel(&c.into_image(), 3, 3), Color32::from_rgb(0, 200, 0));
    }

    #[test]
    fn the_whole_window_is_captured_with_panels_and_viewport() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        let image = app.capture(size(900.0, 600.0)).unwrap();
        assert_eq!(image.size, [900, 600]);
        let t = Tokens::dark();
        // The title bar along the top, the timeline along the bottom (their own colours).
        assert_eq!(pixel(&image, 450, 595), t.timeline, "timeline");
        assert_eq!(pixel(&image, 899, 1), t.app_bar, "title bar");
        // The viewport shows the model: the middle is neither panel nor background.
        let middle = pixel(&image, 450, 330);
        assert!(middle != t.panel && middle != t.app_bar, "{middle:?}");
        // Text was drawn: the title bar has pixels brighter than its fill.
        let bright = (0..300).any(|x| (2..30).any(|y| pixel(&image, x, y).r() > 150));
        assert!(bright, "title text");
        // The live window's state is untouched.
        assert!(app.viewport.rect.is_none() && app.toolbar_rects.is_empty() && app.layout.is_empty());
        // A second capture reuses the offscreen context; at scale 2 it doubles the pixels.
        let twice = app.capture(CaptureSize { width: 450.0, height: 300.0, scale: 2.0 }).unwrap();
        assert_eq!(twice.size, [900, 600]);
    }

    #[test]
    fn hostile_capture_sizes_are_refused() {
        let mut app = ModelApp::new(Session::default());
        for s in [
            size(0.0, 10.0),
            size(f32::NAN, 10.0),
            size(f32::INFINITY, 10.0),
            size(-5.0, 10.0),
            size(10_000.0, 10.0),
            size(8000.0, 8000.0),
            CaptureSize { width: 10.0, height: 10.0, scale: 0.0 },
            CaptureSize { width: 10.0, height: 10.0, scale: f32::NAN },
            CaptureSize { width: 10.0, height: 10.0, scale: 100.0 },
        ] {
            assert!(app.capture(s).is_err(), "{s:?}");
        }
    }
}
