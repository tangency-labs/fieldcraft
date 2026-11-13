//! `view.image.save`: a picture of the viewport written as a PNG, in any workspace and with no
//! solved step (the window's equivalent of the MCP `screenshot` tool).
//!
//! The picture is drawn offscreen as a render's pictures are ([`super::render`]): the host's GPU
//! pipeline when it registered one ([`crate::view::register_drawer`]), else the CPU rasterizer,
//! with the session's camera (the app's viewport, or the last MCP picture's; a `camera` given
//! here instead) and what the results display shows, in the desktop app's light palette.

use fieldcraft_io::roots::Access;
use fieldcraft_io::views::{MAX_SIDE, MIN_SIDE};
use serde_json::{Value, json};

use super::{file_param, opt};
use crate::view::{Aim, Palette};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const SAVE: CommandSpec = CommandSpec {
    id: "view.image.save",
    label: "Save Image",
    params: r#"{"path": "view.png", "size"?: [1600, 1000], "camera"?: {"yaw": -70, "pitch": 18, "distance"?: 300, "target"?: [x, y, z] | {"center_of": "mass" | "bbox", "bodies"?}, "fit"?: true, "fov"?: 0}} — writes a PNG picture of the model as the viewport shows it: the geometry or mesh, the analysis symbols and the results display, with the session's camera (the app's viewport, or the last MCP picture's; `camera` takes ui.camera's parameters instead, in the display units). `size` is [width, height] in pixels, each 16 to 4096. Needs no solved step. Replies {path, size, renderer}"#,
    form: Some(super::forms::IMAGE_SAVE),
    kind: CommandKind::Session,
    requires: &[],
    run: save,
};

/// The picture's size when none is given, in pixels.
pub const DEFAULT_SIZE: [u32; 2] = [1600, 1000];

fn save(s: &mut Session, p: &Value) -> Result<Value> {
    let [w, h] = opt::<[u32; 2]>(p, "size")?.unwrap_or(DEFAULT_SIZE);
    if [w, h].iter().any(|side| !(MIN_SIDE..=MAX_SIDE).contains(side)) {
        return Err(EngineError::invalid("size", format!("[width, height], each {MIN_SIDE} to {MAX_SIDE} pixels")));
    }
    let aim = match p.get("camera").filter(|v| !v.is_null()) {
        Some(c) => crate::view::camera(c)?,
        None => s
            .camera
            .as_ref()
            .map_or_else(|| Aim::new(fieldcraft_render::Orientation::of(fieldcraft_render::StandardView::Iso)).targeted(), Aim::targeted),
    };
    let path = file_param(s, p, "path", Access::Write)?;
    let framed = crate::view::frame_at(s, &aim, w as usize, h as usize, &Palette::LIGHT, None)?;
    let (rgba, renderer) = crate::view::draw_rgba(&framed, w, h)?;
    let png = fieldcraft_render::png(w as usize, h as usize, &rgba)
        .ok_or_else(|| EngineError::Params("view.image.save: the picture could not be encoded".into()))?;
    fieldcraft_io::fs::write_output(std::path::Path::new(&path), &png)?;
    Ok(json!({"path": path, "size": [w, h], "renderer": renderer}))
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use crate::Session;

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("fieldcraft-image-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn png_size(bytes: &[u8]) -> (u32, u32) {
        let d = png::Decoder::new(std::io::Cursor::new(bytes)).read_info().unwrap();
        (d.info().width, d.info().height)
    }

    /// A1-13: a picture of an unsolved model, at the size asked, from the session's camera or
    /// the one given; the session's camera is left as it was.
    #[test]
    fn saves_a_picture_of_an_unsolved_model() {
        let d = dir("model");
        let mut s = Session::default();
        s.run_script(&crate::sample::cantilever()).unwrap();
        let before = (s.camera.clone(), s.camera_requests());
        let path = d.join("view.png");
        let r = s.execute("view.image.save", &json!({"path": path, "size": [96, 64]})).unwrap();
        assert_eq!(r["size"], json!([96, 64]), "{r}");
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(png_size(&bytes), (96, 64));
        let top = d.join("top.png");
        s.execute("view.image.save", &json!({"path": top, "size": [96, 64], "camera": {"yaw": 0, "pitch": 89, "fit": true}})).unwrap();
        assert_ne!(std::fs::read(&top).unwrap(), bytes, "another camera, another picture");
        assert_eq!((s.camera.clone(), s.camera_requests()), before, "saving a picture moves no camera");
        // The default size.
        let full = d.join("full.png");
        s.execute("view.image.save", &json!({"path": full})).unwrap();
        assert_eq!(png_size(&std::fs::read(&full).unwrap()), (super::DEFAULT_SIZE[0], super::DEFAULT_SIZE[1]));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn hostile_image_params_write_nothing() {
        let d = dir("hostile");
        let mut s = Session::default();
        let path = d.join("x.png");
        for p in [
            json!({}),
            json!({"path": ""}),
            json!({"path": path, "size": [1, 1]}),
            json!({"path": path, "size": [100_000, 64]}),
            json!({"path": path, "size": [64]}),
            json!({"path": path, "size": [-5, 64]}),
            json!({"path": path, "size": "big"}),
            json!({"path": path, "camera": {"yaw": "x"}}),
            json!({"path": path, "camera": {"zoom": 2}}),
            json!({"path": path, "camera": {"yaw": 1e300, "pitch": 0, "fov": 1e9}}),
        ] {
            assert!(s.execute("view.image.save", &p).is_err(), "{p}");
        }
        assert_eq!(std::fs::read_dir(&d).unwrap().count(), 0, "nothing written");
        let _ = std::fs::remove_dir_all(&d);
    }
}
