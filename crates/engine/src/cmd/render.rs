//! Renders (`results.render.save`, `results.render.delete`) and `results.export.animation`: how
//! a saved view becomes a clip, kept in the file's view state ([`fieldcraft_io::views::Render`]),
//! and the clip written as numbered PNG pictures with the ffmpeg command line that encodes them
//! (FieldCraft bundles no video encoder).
//!
//! Each picture is drawn offscreen as the MCP `screenshot` tool draws one ([`crate::view`]: the
//! host's GPU pipeline when it registered one, [`crate::view::register_drawer`], else the CPU
//! rasterizer), with the view's camera and results display at one frame, and a caption bar under
//! it: the caption, the contour legend with its range and unit, and the frame's time.

use fieldcraft_io::roots::Access;
use fieldcraft_io::views::{Frames, Hold, Render, RenderFormat};
use fieldcraft_render::Rgb;
use fieldcraft_render::caption::{Bar, Colors, Legend};
use serde_json::{Value, json};

use super::{file_param, req};
use crate::requirement::Requirement;
use crate::results::{Kind, Quantity, unit_of};
use crate::view::{Aim, Palette};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const SAVE: CommandSpec = CommandSpec {
    id: "results.render.save",
    label: "Save Render",
    params: r#"{"name": "Airflow clip", "view": "Cutaway", "step"?: "Flow", "frames"?: "all" | [from, to], "fps"?: 12, "hold"?: {"each": 2, "last": 16}, "size"?: [1280, 720], "caption"?: "AIO radiator · fan air", "legend"?: true, "time_label"?: true, "format"?: "png-sequence"} — keeps how a saved view becomes a clip, saved with the model (a render of that name is replaced): the step whose frames it plays (the view's by default), which frames (from 0, inclusive), the frame rate, how many clip frames each picture and the last one show, the picture size in pixels (the caption bar adds its height under it), a caption, the contour legend and the frame's time under the picture; `format` is png-sequence, numbered PNG pictures with the ffmpeg command line that encodes them (FieldCraft bundles no video encoder). View state, outside the definition"#,
    form: Some(super::forms::RENDER_SAVE),
    kind: CommandKind::Session,
    requires: &[],
    run: save,
};

pub const DELETE: CommandSpec = CommandSpec {
    id: "results.render.delete",
    label: "Delete Render",
    params: r#"{"name": "Airflow clip"} — deletes a saved render"#,
    form: Some(super::forms::RENDER_NAME),
    kind: CommandKind::Session,
    requires: &[],
    run: delete,
};

pub const EXPORT: CommandSpec = CommandSpec {
    id: "results.export.animation",
    label: "Export Animation",
    params: r#"{"render": "Airflow clip", "path": "clip.mp4"} — draws a saved render offscreen, picture by picture: its view's camera and results display at each frame, with the caption, the legend and the time under it, and writes them as numbered PNG pictures beside `path` (`clip_00000.png`, `clip_00001.png`, …; each picture repeated for its hold, so the sequence plays at the render's frame rate). FieldCraft bundles no video encoder: the reply's `ffmpeg` is the command line that encodes the sequence as `path` (an .mp4 or .gif; `<stem>.mp4` otherwise). The session's own display and camera are left as they were. Replies {render, pictures, files, first, size, fps, renderer, ffmpeg}"#,
    form: Some(super::forms::EXPORT_ANIMATION),
    kind: CommandKind::Session,
    requires: &[Requirement::Result],
    run: export,
};

/// Most pictures one clip draws, and most files it writes (holds count).
pub const MAX_PICTURES: usize = 100_000;
pub const MAX_FILES: usize = 1_000_000;

fn save(s: &mut Session, p: &Value) -> Result<Value> {
    let mut params = p.clone();
    // The defaults fill in what is left out (`Render`'s own): the file reads the same.
    if let Some(o) = params.as_object_mut() {
        o.retain(|_, v| !v.is_null());
    }
    let render: Render = serde_json::from_value(params).map_err(|e| EngineError::Params(format!("render: {e}")))?;
    render.check().map_err(|e| EngineError::Params(e.to_string()))?;
    if s.view.saved(&render.view).is_none() {
        return Err(EngineError::Params(format!("no saved view `{}` (view.save keeps one)", render.view)));
    }
    let replaced = match s.view.renders.iter_mut().find(|r| r.name == render.name) {
        Some(r) => {
            *r = render.clone();
            true
        }
        None => {
            if s.view.renders.len() >= fieldcraft_io::views::MAX_VIEWS {
                return Err(EngineError::Params(format!("a model keeps at most {} renders", fieldcraft_io::views::MAX_VIEWS)));
            }
            s.view.renders.push(render.clone());
            false
        }
    };
    Ok(json!({"render": render, "replaced": replaced}))
}

fn delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    let before = s.view.renders.len();
    s.view.renders.retain(|r| r.name != name);
    if s.view.renders.len() == before {
        return Err(EngineError::Params(format!("no render `{}`", name.chars().take(80).collect::<String>())));
    }
    Ok(json!({"renders": s.view.renders.iter().map(|r| &r.name).collect::<Vec<_>>()}))
}

/// The bar's colours: the picture's dark end with light text.
const BAR: Colors = Colors { background: Rgb(0x1e, 0x22, 0x28), text: Rgb(0xf2, 0xf4, 0xf7), dim: Rgb(0xb8, 0xbf, 0xc9) };

/// A number of a legend or label: four significant figures, no trailing zeros.
fn number(v: f64) -> String {
    if !v.is_finite() {
        return "–".into();
    }
    let digits = if v == 0.0 { 0 } else { (3 - v.abs().log10().floor() as i32).clamp(0, 12) as usize };
    let text = format!("{v:.digits$}");
    let text = if text.contains('.') { text.trim_end_matches('0').trim_end_matches('.').to_string() } else { text };
    if text == "-0" { "0".into() } else { text }
}

/// The caption bar of a picture showing `report` (the results display's reply).
fn bar(s: &Session, render: &Render, report: &Value) -> Bar {
    let display = s.units();
    let legend = report.get("contour").filter(|_| render.legend).and_then(|c| {
        let unit = c["unit"].as_str().map(|u| format!(" {u}")).unwrap_or_default();
        // The legend's clamp where it has one, else the field's range.
        let end = |k: &str| c["legend"][k].as_f64().or_else(|| c["range"][k]["value"].as_f64());
        let (low, high) = (end("min")?, end("max")?);
        let bands = s.results_view.as_ref().map_or(crate::results::DEFAULT_BANDS, |v| v.bands);
        let palette = Palette::LIGHT;
        let stops = palette.contour;
        let colors = (0..bands).map(|k| fieldcraft_render::contour::band_color(&stops, k, bands)).collect();
        let title = c["label"].as_str().map(|l| format!("{l}  ")).unwrap_or_default();
        Some(Legend { low: format!("{title}{}{unit}", number(low)), high: format!("{}{unit}", number(high)), bands: colors })
    });
    let detail = if render.time_label {
        match (report.get("time").and_then(Value::as_f64), report.get("load_factor").and_then(Value::as_f64)) {
            (Some(t), _) => Some(format!("t = {}{}", number(t), unit_of(Quantity::Time, display).map(|u| format!(" {u}")).unwrap_or_default())),
            (None, Some(f)) => Some(format!("load factor {}", number(f))),
            _ => None,
        }
    } else {
        None
    };
    Bar { caption: render.caption.clone(), legend, detail }
}

/// `path` with `_NNNNN` before its extension, and its folder and stem: the pictures' names.
fn picture_path(path: &str, n: usize) -> String {
    let p = std::path::Path::new(path);
    let stem = p.file_stem().map_or_else(|| "clip".into(), |s| s.to_string_lossy().into_owned());
    p.with_file_name(format!("{stem}_{n:05}.png")).to_string_lossy().into_owned()
}

/// A path quoted for a POSIX shell.
fn quoted(path: &str) -> String {
    format!("'{}'", path.replace('\'', r"'\''"))
}

/// The ffmpeg command line that encodes the pictures beside `path` at `fps` as `path` (an .mp4 or
/// .gif), or as `<stem>.mp4`.
fn ffmpeg(path: &str, fps: u32) -> String {
    let p = std::path::Path::new(path);
    let pattern = picture_path(path, 0).replace("_00000.png", "_%05d.png");
    let out = match p.extension().and_then(|e| e.to_str()) {
        Some("mp4" | "gif" | "webm" | "mov" | "mkv") => path.to_string(),
        _ => p.with_extension("mp4").to_string_lossy().into_owned(),
    };
    let codec = if out.ends_with(".gif") {
        String::new()
    } else {
        // Even sides, as yuv420p needs.
        " -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -c:v libx264 -pix_fmt yuv420p -movflags +faststart".to_string()
    };
    format!("ffmpeg -framerate {fps} -i {}{codec} {}", quoted(&pattern), quoted(&out))
}

fn export(s: &mut Session, p: &Value) -> Result<Value> {
    let name: String = req(p, "render")?;
    let render = s.view.render(&name).cloned().ok_or_else(|| {
        let names: Vec<&str> = s.view.renders.iter().map(|r| r.name.as_str()).collect();
        EngineError::Params(format!(
            "no render `{}` (renders: {})",
            name.chars().take(80).collect::<String>(),
            if names.is_empty() { "none: results.render.save keeps one".into() } else { names.join(", ") }
        ))
    })?;
    let path = file_param(s, p, "path", Access::Write)?;
    // What the session shows is put back afterwards: exporting is not looking.
    let kept = (s.results_view.clone(), s.camera.clone(), s.view.current.clone(), s.pending_view.clone(), s.camera_requests);
    let out = draw_clip(s, &render, &path);
    (s.results_view, s.camera, s.view.current, s.pending_view, s.camera_requests) = kept;
    out
}

fn draw_clip(s: &mut Session, render: &Render, path: &str) -> Result<Value> {
    let shown_view = super::views::show(s, &render.view)?;
    if shown_view["pending"] == true {
        return Err(EngineError::Params(format!("view `{}` shows a step that has no result yet: solve it first", render.view)));
    }
    let step = match (&render.step, s.view.saved(&render.view).and_then(|v| v.results.as_ref())) {
        (Some(step), _) => step.clone(),
        (None, Some(r)) => r.step.clone(),
        (None, None) => s
            .results_view
            .as_ref()
            .map(|v| v.step.clone())
            .ok_or_else(|| EngineError::Params(format!("view `{}` shows no results: give the render a `step`", render.view)))?,
    };
    // The step's frames (one still when it has none).
    let count = match s.latest_result(&step).map(|(r, _)| r.kind()) {
        None => return Err(EngineError::Params(format!("step `{step}` has no result: solve it first"))),
        Some(Kind::Explicit(e)) => e.frames.len(),
        Some(Kind::Fluid(f)) => f.frames.len(),
        Some(Kind::Static(r)) => r.nonlinear.as_ref().map_or(0, |h| h.frames.len()),
        Some(Kind::Frequency(_) | Kind::Heat(_)) => 0,
    };
    let units = s.units().clone();
    let frames: Vec<Option<usize>> = if count == 0 {
        vec![None]
    } else {
        let (from, to) = match render.frames {
            Frames::All => (0, count - 1),
            Frames::Range([a, b]) => {
                let (a, b) = (usize::try_from(a).unwrap_or(usize::MAX), usize::try_from(b).unwrap_or(usize::MAX));
                if b >= count {
                    return Err(EngineError::Params(format!("render `{}`: step `{step}` has frames 0 to {}", render.name, count - 1)));
                }
                (a, b)
            }
        };
        (from..=to).map(Some).collect()
    };
    let Hold { each, last } = render.hold;
    let files = frames.len().saturating_sub(1).saturating_mul(each as usize).saturating_add(last as usize);
    if frames.len() > MAX_PICTURES || files > MAX_FILES {
        return Err(EngineError::Params(format!(
            "render `{}` would draw {} pictures into {files} files (at most {MAX_PICTURES} and {MAX_FILES})",
            render.name,
            frames.len()
        )));
    }
    let RenderFormat::PngSequence = render.format;
    let aim = s
        .camera
        .as_ref()
        .map_or_else(|| Aim::new(fieldcraft_render::Orientation::of(fieldcraft_render::StandardView::Iso)).targeted(), Aim::targeted);
    let [w, h] = render.size;
    let mut written = 0usize;
    let mut size = [w, h];
    let mut renderer = String::new();
    for (i, k) in frames.iter().enumerate() {
        let report = match k {
            Some(k) => super::views::run_in(s, "results.frame", &json!({"step": step, "frame": k}), &units)?,
            None => s.results_report(),
        };
        let framed = crate::view::frame_at(s, &aim, w as usize, h as usize, &Palette::LIGHT, None)?;
        let (rgba, drawn_by) = crate::view::draw_rgba(&framed, w, h)?;
        renderer = drawn_by;
        let (pixels, total) = fieldcraft_render::caption::under(rgba, w as usize, h as usize, &bar(s, render, &report), &BAR);
        size = [w, u32::try_from(total).unwrap_or(u32::MAX)];
        let png = fieldcraft_render::png(w as usize, total, &pixels)
            .ok_or_else(|| EngineError::Params("render: the picture could not be encoded".into()))?;
        let repeat = if i + 1 == frames.len() { last } else { each };
        for _ in 0..repeat {
            fieldcraft_io::fs::write_output(std::path::Path::new(&picture_path(path, written)), &png)?;
            written += 1;
        }
    }
    Ok(json!({
        "render": render.name,
        "pictures": frames.len(),
        "files": written,
        "first": picture_path(path, 0),
        "size": size,
        "fps": render.fps,
        "renderer": renderer,
        "ffmpeg": ffmpeg(path, render.fps),
    }))
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("fieldcraft-render-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// The size of a PNG.
    fn png_size(bytes: &[u8]) -> (u32, u32) {
        let d = png::Decoder::new(std::io::Cursor::new(bytes)).read_info().unwrap();
        (d.info().width, d.info().height)
    }

    /// #646: a saved render draws its view at every frame, with the caption bar, as numbered
    /// pictures held as it says, and the ffmpeg line that encodes them; the session's own display
    /// stays; the render is saved with the model.
    #[test]
    fn a_render_writes_its_clip_as_a_png_sequence() {
        let d = dir("clip");
        let mut s = crate::cmd::views::tests::cutaway();
        s.execute("view.save", &json!({"name": "Cutaway", "camera": {"yaw": -70, "pitch": 18, "fit": true}})).unwrap();
        let r = s
            .execute(
                "results.render.save",
                &json!({"name": "Clip", "view": "Cutaway", "hold": {"each": 2, "last": 3}, "size": [96, 64], "caption": "Couette gap · heated wall"}),
            )
            .unwrap();
        assert_eq!(r["render"]["fps"], 12);
        s.execute("results.frame", &json!({"frame": 0})).unwrap();
        let before = (s.results_view.clone(), s.camera.clone(), s.camera_requests());
        let path = d.join("clip.mp4");
        let r = s.execute("results.export.animation", &json!({"render": "Clip", "path": path})).unwrap();
        assert_eq!((r["pictures"].as_u64(), r["files"].as_u64()), (Some(5), Some(4 * 2 + 3)), "{r}");
        let height = 64 + fieldcraft_render::caption::height(64) as u32;
        assert_eq!(r["size"], json!([96, height]));
        let ffmpeg = r["ffmpeg"].as_str().unwrap();
        assert!(ffmpeg.starts_with("ffmpeg -framerate 12 -i ") && ffmpeg.contains("clip_%05d.png") && ffmpeg.ends_with("clip.mp4'"), "{ffmpeg}");
        let file = |n: usize| std::fs::read(d.join(format!("clip_{n:05}.png"))).unwrap();
        assert_eq!(png_size(&file(0)), (96, height));
        assert_eq!(file(0), file(1), "each picture held for two frames");
        assert_ne!(file(1), file(2), "the next frame");
        assert_eq!(file(8), file(10), "the last held for three");
        assert!(!d.join("clip_00011.png").exists());
        assert_eq!((s.results_view.clone(), s.camera.clone(), s.camera_requests()), before, "exporting is not looking");

        // Saved with the model, it draws the same from the file.
        let model = d.join("couette.fieldcraft");
        s.execute("file.save_as", &json!({"path": model})).unwrap();
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": model})).unwrap();
        assert_eq!(t.view.renders, s.view.renders);
        let again = d.join("again.mp4");
        t.execute("results.export.animation", &json!({"render": "Clip", "path": again})).unwrap();
        assert_eq!(std::fs::read(d.join("again_00010.png")).unwrap(), file(10));
        // A view a render starts from is not deleted.
        assert!(t.execute("view.delete", &json!({"name": "Cutaway"})).unwrap_err().to_string().contains("render `Clip`"));
        t.execute("results.render.delete", &json!({"name": "Clip"})).unwrap();
        t.execute("view.delete", &json!({"name": "Cutaway"})).unwrap();
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn hostile_render_params_change_nothing() {
        let d = dir("hostile");
        let mut s = crate::cmd::views::tests::cutaway();
        s.execute("view.save", &json!({"name": "V"})).unwrap();
        for p in [
            json!({}),
            json!({"name": "R"}),
            json!({"name": "R", "view": "Nope"}),
            json!({"name": "R", "view": "V", "fps": 0}),
            json!({"name": "R", "view": "V", "fps": 1e300}),
            json!({"name": "R", "view": "V", "size": [1, 1]}),
            json!({"name": "R", "view": "V", "size": [100000, 64]}),
            json!({"name": "R", "view": "V", "hold": {"each": 0, "last": 1}}),
            json!({"name": "R", "view": "V", "frames": [3, 1]}),
            json!({"name": "R", "view": "V", "frames": "some"}),
            json!({"name": "R", "view": "V", "format": "mp4"}),
            json!({"name": "R", "view": "V", "caption": "a\nb"}),
            json!({"name": "R", "view": "V", "step": "1x"}),
        ] {
            assert!(s.execute("results.render.save", &p).is_err(), "{p}");
        }
        assert!(s.view.renders.is_empty());
        s.execute("results.render.save", &json!({"name": "Far", "view": "V", "frames": [0, 99], "size": [32, 32]})).unwrap();
        let path = d.join("x.mp4");
        for p in [
            json!({"render": "Nope", "path": path}),
            json!({"render": "Far", "path": path}),
            json!({"render": "Far"}),
            json!({"render": "Far", "path": ""}),
        ] {
            assert!(s.execute("results.export.animation", &p).is_err(), "{p}");
        }
        assert_eq!(std::fs::read_dir(&d).unwrap().count(), 0, "nothing written");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn labels_and_command_lines_read_well() {
        assert_eq!(number(0.25), "0.25");
        assert_eq!(number(45.0), "45");
        assert_eq!(number(12345.678), "12346");
        assert_eq!(number(-0.0001234), "-0.0001234");
        assert_eq!(number(f64::NAN), "–");
        assert_eq!(picture_path("/a/b/clip.mp4", 7), "/a/b/clip_00007.png");
        assert_eq!(ffmpeg("/a/it's.gif", 10), "ffmpeg -framerate 10 -i '/a/it'\\''s_%05d.png' '/a/it'\\''s.gif'");
        assert!(ffmpeg("/a/clip.png", 12).ends_with("'/a/clip.mp4'"));
    }
}
