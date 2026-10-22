//! Offscreen pictures with the viewport's pipelines ([`crate::view`]): the scene, the setup's
//! symbols and the overlay drawn with the window's MSAA and depth test into a texture, copied back
//! and laid over the viewport's background gradient. Every GPU call sits inside error scopes,
//! so a device error (out of memory, validation) returns as an error, never a panic.

use std::sync::{Arc, OnceLock};

use fieldcraft_render::{Camera, Rgb, Scene};

use crate::select::{self, Policy};
use crate::view::{Frame, GpuScene, GpuTarget, Layer, Pipelines};
use crate::{GpuError, Result};
use crate::{Renderer, Tier};

/// Largest picture side and pixel count drawn (the CPU rasterizer's limits).
pub const MAX_SIDE: u32 = fieldcraft_render::MAX_DIM as u32;
pub const MAX_PIXELS: u64 = fieldcraft_render::MAX_PIXELS as u64;
/// The viewport's MSAA sample count, when the device offers it.
pub const SAMPLES: u32 = 4;
const FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
const DEPTH: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
/// The background when a scene gives none (the CPU rasterizer's).
const BACKGROUND: (Rgb, Rgb) = (Rgb(236, 240, 245), Rgb(205, 212, 222));

/// A device ready to draw pictures.
pub struct Capturer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipelines: Pipelines,
    target: GpuTarget,
    /// The adapter's name and tier, for reports.
    pub adapter: String,
    pub tier: Tier,
}

/// The error a popped error scope caught: a native device answers once it has been polled.
fn popped(f: impl std::future::Future<Output = Option<wgpu::Error>>) -> Option<wgpu::Error> {
    let mut f = std::pin::pin!(f);
    match f.as_mut().poll(&mut std::task::Context::from_waker(std::task::Waker::noop())) {
        std::task::Poll::Ready(e) => e,
        std::task::Poll::Pending => None,
    }
}

impl Capturer {
    /// The capturer `policy` picks ([`select::adapter`]); `Ok(None)` for the CPU rasterizer.
    pub fn open(policy: Policy) -> Result<Option<Capturer>> {
        match select::adapter(policy)? {
            Some((adapter, tier)) => Capturer::new(&adapter, tier).map(Some),
            None => Ok(None),
        }
    }

    /// A device on `adapter` with the viewport's pipelines.
    pub fn new(adapter: &wgpu::Adapter, tier: Tier) -> Result<Capturer> {
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("fieldcraft capture"),
            required_limits: wgpu::Limits::downlevel_defaults().using_resolution(adapter.limits()),
            ..Default::default()
        }))
        .map_err(|e| GpuError::Device(e.to_string()))?;
        // Errors outside a scope are logged, replacing wgpu's default handler, which panics.
        device.on_uncaptured_error(Arc::new(|e: wgpu::Error| log::error!("graphics device error: {e}")));
        Capturer::on_device(adapter, device, queue, tier)
    }

    /// A capturer on a device already open on `adapter` (the desktop app's), which keeps its
    /// error handler: pictures share the window's device instead of opening another.
    pub fn on_device(adapter: &wgpu::Adapter, device: wgpu::Device, queue: wgpu::Queue, tier: Tier) -> Result<Capturer> {
        let samples = if adapter.get_texture_format_features(FORMAT).flags.sample_count_supported(SAMPLES)
            && adapter.get_texture_format_features(DEPTH).flags.sample_count_supported(SAMPLES)
        {
            SAMPLES
        } else {
            1
        };
        let target = GpuTarget { format: FORMAT, depth: Some(DEPTH), samples };
        let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
        let pipelines = Pipelines::new(&device, target);
        if let Some(e) = popped(scope.pop()) {
            return Err(GpuError::Device(e.to_string()));
        }
        let info = adapter.get_info();
        Ok(Capturer { device, queue, pipelines, target, adapter: format!("{} ({:?})", info.name, info.backend), tier })
    }

    /// The MSAA sample count pictures are drawn with.
    pub fn samples(&self) -> u32 {
        self.target.samples
    }

    /// The device pictures are drawn on, its queue and the target its passes have, for a caller
    /// that draws through its own [`Pipelines`] in [`Capturer::offscreen`].
    pub fn device(&self) -> (&wgpu::Device, &wgpu::Queue, GpuTarget) {
        (&self.device, &self.queue, self.target)
    }

    /// RGBA pixels (`w × h`, rows top to bottom) of `scene` with `layers` (the symbols, then
    /// the overlay) seen by `cam`, over the scene's background. Widths are taken as pixels.
    pub fn rgba(&self, scene: &Scene, layers: &[&Scene], cam: &Camera, w: u32, h: u32) -> Result<Vec<u8>> {
        let (d, q) = (&self.device, &self.queue);
        let drawn = self.offscreen(w, h, |pass| {
            let mut uploaded = [Layer::default(), Layer::default(), Layer::default()];
            for (layer, s) in uploaded.iter_mut().zip(std::iter::once(scene).chain(layers.iter().copied())) {
                layer.upload(d, q, &GpuScene::from_scene(s, 1.0));
            }
            self.pipelines.write_frame(q, &Frame::of(cam, w as f32, h as f32, scene.radius));
            let [s, y, o] = &uploaded;
            self.pipelines.paint(pass, s, y, o);
        })?;
        finish(&drawn, scene, layers, cam, w, h)
    }

    /// The premultiplied RGBA pixels (`w × h`, rows top to bottom, transparent where nothing
    /// is drawn) of one render pass with the viewport's MSAA and depth, whose contents `paint`
    /// records (uploads included). Sizes are checked first; device errors in `paint` and the
    /// read-back are caught and returned.
    pub fn offscreen(&self, w: u32, h: u32, paint: impl FnOnce(&mut wgpu::RenderPass<'static>)) -> Result<Vec<u8>> {
        let max = self.device.limits().max_texture_dimension_2d.min(MAX_SIDE);
        if w == 0 || h == 0 || w > max || h > max || u64::from(w) * u64::from(h) > MAX_PIXELS {
            return Err(GpuError::Size(format!("picture size {w} x {h} is out of range (1..={max} per side, at most {MAX_PIXELS} pixels)")));
        }
        let oom = self.device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
        let validation = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let drawn = self.draw(w, h, paint);
        let caught = popped(validation.pop()).or_else(|| popped(oom.pop()));
        if let Some(e) = caught {
            return Err(GpuError::Failed(e.to_string()));
        }
        drawn
    }

    fn draw(&self, w: u32, h: u32, paint: impl FnOnce(&mut wgpu::RenderPass<'static>)) -> Result<Vec<u8>> {
        let d = &self.device;
        let size = wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 };
        let texture = |label: &str, format: wgpu::TextureFormat, samples: u32, usage: wgpu::TextureUsages| {
            d.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size,
                mip_level_count: 1,
                sample_count: samples,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats: &[],
            })
        };
        let samples = self.target.samples;
        let out = texture("fc_capture", FORMAT, 1, wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC);
        let out_view = out.create_view(&Default::default());
        let msaa = (samples > 1).then(|| texture("fc_capture_msaa", FORMAT, samples, wgpu::TextureUsages::RENDER_ATTACHMENT));
        let msaa_view = msaa.as_ref().map(|t| t.create_view(&Default::default()));
        let depth = texture("fc_capture_depth", DEPTH, samples, wgpu::TextureUsages::RENDER_ATTACHMENT);
        let depth_view = depth.create_view(&Default::default());
        let (view, resolve) = match &msaa_view {
            Some(m) => (m, Some(&out_view)),
            None => (&out_view, None),
        };
        let row = padded_row(w);
        let buffer = d.create_buffer(&wgpu::BufferDescriptor {
            label: Some("fc_capture_read"),
            size: u64::from(row) * u64::from(h),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut encoder = d.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("fc_capture") });
        {
            let mut pass = encoder
                .begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("fc_capture"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view,
                        depth_slice: None,
                        resolve_target: resolve,
                        ops: wgpu::Operations { load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT), store: wgpu::StoreOp::Store },
                    })],
                    depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                        view: &depth_view,
                        depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Discard }),
                        stencil_ops: None,
                    }),
                    timestamp_writes: None,
                    occlusion_query_set: None,
                    multiview_mask: None,
                })
                .forget_lifetime();
            paint(&mut pass);
        }
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo { texture: &out, mip_level: 0, origin: wgpu::Origin3d::ZERO, aspect: wgpu::TextureAspect::All },
            wgpu::TexelCopyBufferInfo {
                buffer: &buffer,
                layout: wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(row), rows_per_image: Some(h) },
            },
            size,
        );
        self.queue.submit([encoder.finish()]);
        let slice = buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        d.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| GpuError::Failed(e.to_string()))?;
        rx.recv().map_err(|e| GpuError::Failed(e.to_string()))?.map_err(|e| GpuError::Failed(e.to_string()))?;
        let mapped = slice.get_mapped_range().map_err(|e| GpuError::Failed(e.to_string()))?;
        let pixels = unpad(&mapped, w as usize, h as usize, row as usize);
        drop(mapped);
        buffer.unmap();
        pixels.ok_or_else(|| GpuError::Failed("the picture read back has the wrong size".into()))
    }
}

/// The process's capturer, opened on first use by the `FIELDCRAFT_GPU` policy: `Ok(None)` for
/// the CPU rasterizer. Its error is kept, so a missing device is not looked for again.
pub fn shared() -> &'static Result<Option<Capturer>> {
    static SHARED: OnceLock<Result<Option<Capturer>>> = OnceLock::new();
    SHARED.get_or_init(|| Policy::from_env().and_then(Capturer::open))
}

/// A PNG of `scene` with `layers` seen by `cam`, `w` × `h` pixels, drawn by `capturer` (a
/// [`shared`] result): on its GPU, or by the CPU rasterizer when it is `Ok(None)`. Under the
/// `auto` policy a GPU failure falls back to the CPU and the reason is reported; under a strict
/// one it is the error.
pub fn png(
    capturer: &Result<Option<Capturer>>,
    policy: Policy,
    scene: &Scene,
    layers: &[&Scene],
    cam: &Camera,
    w: u32,
    h: u32,
) -> Result<(Vec<u8>, Renderer)> {
    let (px, renderer) = rgba(capturer, policy, scene, layers, cam, w, h)?;
    let png = fieldcraft_render::png(w as usize, h as usize, &px).ok_or_else(|| GpuError::Failed("the picture could not be encoded".into()))?;
    Ok((png, renderer))
}

/// [`png`]'s picture as RGBA pixels, for a caller that draws over it (a clip's caption bar).
pub fn rgba(
    capturer: &Result<Option<Capturer>>,
    policy: Policy,
    scene: &Scene,
    layers: &[&Scene],
    cam: &Camera,
    w: u32,
    h: u32,
) -> Result<(Vec<u8>, Renderer)> {
    let gpu = match capturer {
        Ok(Some(c)) => c.rgba(scene, layers, cam, w, h).map(|px| Some((px, c))),
        Ok(None) => Ok(None),
        Err(e) => Err(e.clone()),
    };
    let why = match gpu {
        Ok(Some((px, c))) => return Ok((px, Renderer::Gpu { adapter: c.adapter.clone(), tier: c.tier, samples: c.samples() })),
        Ok(None) => (policy == Policy::Auto).then(|| "no GPU adapter".to_string()),
        Err(e @ GpuError::Size(_)) => return Err(e),
        Err(e) if policy == Policy::Auto => Some(e.to_string()),
        Err(e) => return Err(e),
    };
    let c = fieldcraft_render::render_layers(scene, layers, cam, w as usize, h as usize).map_err(|e| GpuError::Size(e.to_string()))?;
    Ok((c.rgba, Renderer::Cpu { why }))
}

/// [`rgba`] on the process's [`shared`] capturer under the `FIELDCRAFT_GPU` policy, with what drew
/// it as JSON text: what a host registers as the engine's picture drawer
/// (`fieldcraft_engine::view::register_drawer`).
pub fn rgba_shared(scene: &Scene, layers: &[&Scene], cam: &Camera, w: u32, h: u32) -> std::result::Result<(Vec<u8>, String), String> {
    let policy = Policy::from_env().map_err(|e| e.to_string())?;
    let (px, r) = rgba(shared(), policy, scene, layers, cam, w, h).map_err(|e| e.to_string())?;
    Ok((px, r.report().to_string()))
}

/// A picture from [`Capturer::offscreen`]'s pixels of `scene` and `layers`: laid over the
/// scene's background gradient, with the name labels drawn in the CPU rasterizer's font (the
/// window draws them with the interface's).
pub fn finish(drawn: &[u8], scene: &Scene, layers: &[&Scene], cam: &Camera, w: u32, h: u32) -> Result<Vec<u8>> {
    let pixels = over_background(drawn, w as usize, h as usize, scene.background.unwrap_or(BACKGROUND));
    let layers: Vec<&Scene> = std::iter::once(scene).chain(layers.iter().copied()).collect();
    fieldcraft_render::draw_labels(&pixels, w as usize, h as usize, &layers, cam, scene.radius).map_err(|e| GpuError::Failed(e.to_string()))
}

/// Bytes per row of a `w` pixel RGBA copy, padded to wgpu's 256-byte alignment.
pub fn padded_row(w: u32) -> u32 {
    let a = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    (w * 4).div_ceil(a) * a
}

/// The `w × h` RGBA pixels of rows `row` bytes apart; `None` when `bytes` is too short.
fn unpad(bytes: &[u8], w: usize, h: usize, row: usize) -> Option<Vec<u8>> {
    let mut out = Vec::new();
    out.try_reserve_exact(w.checked_mul(h)?.checked_mul(4)?).ok()?;
    for y in 0..h {
        out.extend_from_slice(bytes.get(y * row..y * row + w * 4)?);
    }
    Some(out)
}

/// Premultiplied pixels over the vertical gradient `bg` (top, bottom), as the viewport draws
/// its background behind the model. Opaque pixels are unchanged.
fn over_background(px: &[u8], w: usize, h: usize, bg: (Rgb, Rgb)) -> Vec<u8> {
    let f = |c: Rgb| [c.0, c.1, c.2].map(f64::from);
    let (top, bottom) = (f(bg.0), f(bg.1));
    let mut out = px.to_vec();
    for (i, [r, g, b, a]) in out.as_chunks_mut::<4>().0.iter_mut().enumerate() {
        let t = (i / w.max(1)) as f64 / h.max(1) as f64;
        let k = 1.0 - f64::from(*a) / 255.0;
        for (c, (top, bottom)) in [r, g, b].into_iter().zip(top.iter().zip(&bottom)) {
            let back = top + (bottom - top) * t;
            *c = (f64::from(*c) + back * k).round().clamp(0.0, 255.0) as u8;
        }
        *a = 255;
    }
    out
}

#[cfg(test)]
mod tests {
    use fieldcraft_render::{Aabb3, SceneLine, SceneTri, Vec3};

    use super::*;

    #[test]
    fn rows_are_padded_and_unpadded() {
        assert_eq!(padded_row(1), 256);
        assert_eq!(padded_row(64), 256);
        assert_eq!(padded_row(65), 512);
        let bytes: Vec<u8> = (0..512u32).map(|i| i as u8).collect();
        assert_eq!(unpad(&bytes, 2, 2, 256), Some(vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7]));
        assert_eq!(unpad(&bytes, 2, 3, 256), None);
    }

    #[test]
    fn transparent_pixels_take_the_background() {
        let bg = (Rgb(0, 0, 0), Rgb(200, 100, 0));
        let px = [0, 0, 0, 0, 10, 20, 30, 255, 0, 0, 0, 0, 0, 0, 0, 0];
        let out = over_background(&px, 2, 2, bg);
        assert_eq!(&out[0..4], &[0, 0, 0, 255]);
        assert_eq!(&out[4..8], &[10, 20, 30, 255]);
        assert_eq!(&out[8..12], &[100, 50, 0, 255]);
    }

    /// A red triangle with a black edge, framed.
    fn triangle() -> (Scene, Camera) {
        let p = [Vec3::new(-1.0, 0.0, -1.0), Vec3::new(1.0, 0.0, -1.0), Vec3::new(0.0, 0.0, 1.0)];
        let scene = Scene {
            tris: vec![SceneTri { p, normal: Vec3::new(0.0, -1.0, 0.0), color: Rgb(220, 30, 30) }],
            lines: vec![SceneLine { a: p[0], b: p[2], color: Rgb(0, 0, 0), width: 3.0 }],
            background: Some((Rgb(255, 255, 255), Rgb(255, 255, 255))),
            radius: 4.0,
            ..Default::default()
        };
        let mut cam = Camera::default();
        cam.set_view(fieldcraft_render::StandardView::Front);
        cam.fit_in(&Aabb3 { min: Vec3::new(-1.0, 0.0, -1.0), max: Vec3::new(1.0, 0.0, 1.0) }, 1.0, fieldcraft_render::FULL_VIEW);
        (scene, cam)
    }

    /// A picture on whatever adapter the machine has: the right size, the triangle's red and
    /// the edge's black on the white background. A machine with no adapter is reported, not
    /// failed (`auto` then draws on the CPU), unless `FIELDCRAFT_REQUIRE_GPU=1` (CI).
    #[test]
    fn draws_offscreen_when_an_adapter_exists() {
        let cap = match Capturer::open(Policy::from_env().unwrap()) {
            Ok(Some(c)) => c,
            Ok(None) if crate::select::gpu_required() => panic!("FIELDCRAFT_REQUIRE_GPU=1 but no GPU adapter"),
            Err(e) if crate::select::gpu_required() => panic!("FIELDCRAFT_REQUIRE_GPU=1 but no GPU device: {e}"),
            Ok(None) => return eprintln!("no GPU adapter: skipped"),
            Err(e) => return eprintln!("no GPU device ({e}): skipped"),
        };
        let (scene, cam) = triangle();
        let px = cap.rgba(&scene, &[], &cam, 96, 64).unwrap();
        assert_eq!(px.len(), 96 * 64 * 4);
        let count = |f: &dyn Fn(&[u8]) -> bool| px.chunks(4).filter(|p| f(p)).count();
        assert!(count(&|p| p[0] > 150 && p[1] < 80 && p[2] < 80) > 200, "red faces");
        assert!(count(&|p| p[0] < 40 && p[1] < 40 && p[2] < 40) > 20, "black edge");
        assert!(count(&|p| p == [255, 255, 255, 255]) > 1000, "background");
        // Hostile sizes are refused before anything is allocated.
        for (w, h) in [(0, 10), (10, 0), (MAX_SIDE + 1, 10), (u32::MAX, u32::MAX), (8192, 8192)] {
            assert!(matches!(cap.rgba(&scene, &[], &cam, w, h), Err(GpuError::Size(_))), "{w} x {h}");
        }
    }
}
