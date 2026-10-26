//! `fieldcraft-cli snapshot`: a script in, a PNG out, with no window and no display. The picture
//! is the viewport's: the engine's scene ([`fieldcraft_engine::view::frame`]) drawn with the
//! window's GPU pipelines ([`fieldcraft_gpu::capture`]) on a hardware GPU or a software one
//! (Vulkan: lavapipe or SwiftShader), and by the CPU rasterizer as the last tier (`FIELDCRAFT_GPU`,
//! [`fieldcraft_gpu::select`]).

use fieldcraft_engine::view::{self, Aim, Palette, Section, Target};
use fieldcraft_gpu::Renderer;
use fieldcraft_gpu::capture::{self, Capturer};
use fieldcraft_gpu::select::Policy;
use fieldcraft_render::{Orientation, StandardView};

/// The picture size when `--size` is not given.
const DEFAULT_SIZE: (u32, u32) = (1280, 800);

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub script: String,
    pub out: String,
    pub view: StandardView,
    /// A camera aimed as the MCP `screenshot`'s (`--camera JSON`; not with `--view`).
    pub camera: Option<Aim>,
    pub size: (u32, u32),
    /// A section through the mesh (`--section x|y|z[=VALUE]`).
    pub section: Option<Section>,
    /// A nodal field contoured on the last solved step (`--field von_mises`).
    pub field: Option<String>,
}

/// A view name: a standard view, or the axis looked along (`x`: from the right, `y`: from the
/// front, `z`: from the top).
fn view_of(s: &str) -> Option<StandardView> {
    match s.to_ascii_lowercase().as_str() {
        "x" => Some(StandardView::Right),
        "y" => Some(StandardView::Front),
        "z" => Some(StandardView::Top),
        other => StandardView::parse(other),
    }
}

/// `WxH`, each side 1..=the renderer's limit.
fn size_of(s: &str) -> Option<(u32, u32)> {
    let (w, h) = s.split_once(['x', 'X'])?;
    let (w, h) = (w.trim().parse::<u32>().ok()?, h.trim().parse::<u32>().ok()?);
    let max = fieldcraft_gpu::capture::MAX_SIDE;
    ((1..=max).contains(&w) && (1..=max).contains(&h) && u64::from(w) * u64::from(h) <= fieldcraft_gpu::capture::MAX_PIXELS).then_some((w, h))
}

/// The options: the script, then each flag once, as `--flag value` or `--flag=value`.
pub fn options(args: &[String]) -> Result<Options, String> {
    let (mut script, mut out, mut view, mut size, mut section, mut field, mut camera) = (None, None, None, None, None, None, None);
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        if !arg.starts_with("--") {
            if script.replace(arg.clone()).is_some() {
                return Err(format!("snapshot takes one script, not `{arg}` too"));
            }
            continue;
        }
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) => (f, Some(v.to_string())),
            None => (arg.as_str(), None),
        };
        let value = inline.or_else(|| it.next().cloned()).filter(|v| !v.is_empty()).ok_or_else(|| format!("{flag} needs a value"))?;
        let twice = match flag {
            "--out" => out.replace(value).is_some(),
            "--view" => view
                .replace(
                    view_of(&value).ok_or_else(|| format!("--view: unknown view `{value}` (iso, x, y, z, front, back, left, right, top, bottom)"))?,
                )
                .is_some(),
            "--size" => size
                .replace(size_of(&value).ok_or_else(|| {
                    format!(
                        "--size: `{value}` is not WxH with each side 1..={} and at most {} pixels",
                        fieldcraft_gpu::capture::MAX_SIDE,
                        fieldcraft_gpu::capture::MAX_PIXELS
                    )
                })?)
                .is_some(),
            "--section" => section
                .replace(Section::parse(&value).ok_or_else(|| format!("--section: `{value}` is not x, y, z or AXIS=VALUE (x=12.5)"))?)
                .is_some(),
            "--field" => field.replace(value).is_some(),
            "--camera" => {
                let v: serde_json::Value = serde_json::from_str(&value).map_err(|e| format!("--camera: not JSON: {e}"))?;
                camera.replace(view::camera(&v).map_err(|e| format!("--camera: {e}"))?).is_some()
            }
            other => return Err(format!("snapshot: unknown option `{other}`")),
        };
        if twice {
            return Err(format!("{flag} is given twice"));
        }
    }
    if view.is_some() && camera.is_some() {
        return Err("give --view or --camera, not both".into());
    }
    Ok(Options {
        script: script.ok_or("snapshot needs a script")?,
        out: out.ok_or("snapshot needs --out FILE.png")?,
        view: view.unwrap_or(StandardView::Iso),
        camera,
        size: size.unwrap_or(DEFAULT_SIZE),
        section,
        field,
    })
}

/// The PNG of the session `o` describes, drawn as `policy` says ([`capture::png`]).
pub fn picture(o: &Options, policy: Policy) -> Result<(Vec<u8>, Renderer, view::Shown), String> {
    let mut session = crate::session_of(&o.script)?;
    if let Some(field) = &o.field {
        session.execute("results.contour", &serde_json::json!({"field": field})).map_err(|e| format!("--field: {e}"))?;
    }
    let (w, h) = o.size;
    let f = match &o.camera {
        Some(aim) => view::frame_at(&session, aim, w as usize, h as usize, &Palette::LIGHT, o.section).map_err(|e| e.to_string())?,
        None => view::frame(&session, &Orientation::of(o.view), w as usize, h as usize, &Palette::LIGHT, o.section),
    };
    let capturer = Capturer::open(policy);
    let (png, renderer) = capture::png(&capturer, policy, &f.scene, &[&f.symbols], &f.cam, w, h).map_err(|e| e.to_string())?;
    Ok((png, renderer, f.shown))
}

pub fn run(args: &[String]) -> Result<(), String> {
    let o = options(args)?;
    let policy = Policy::from_env().map_err(|e| e.to_string())?;
    let (png, renderer, shown) = picture(&o, policy)?;
    fieldcraft_engine::write_output(std::path::Path::new(&o.out), &png).map_err(|e| e.to_string())?;
    let mut report = crate::picture_report(&renderer, &shown);
    if let Some(r) = report.as_object_mut() {
        r.insert("out".into(), o.out.clone().into());
        r.insert("width".into(), o.size.0.into());
        r.insert("height".into(), o.size.1.into());
        match &o.camera {
            Some(aim) => r.insert("camera".into(), Aim { target: Some(aim.target.clone().unwrap_or(Target::MODEL)), ..aim.clone() }.json()),
            None => r.insert("view".into(), o.view.label().into()),
        };
    }
    println!("{report}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(args: &[&str]) -> Result<Options, String> {
        options(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
    }

    #[test]
    fn options_are_checked() {
        let o = opts(&["s.json", "--out", "a.png"]).unwrap();
        assert_eq!(
            o,
            Options {
                script: "s.json".into(),
                out: "a.png".into(),
                view: StandardView::Iso,
                camera: None,
                size: DEFAULT_SIZE,
                section: None,
                field: None
            }
        );
        let o = opts(&["s", "--out", "o", "--section", "y=2.5", "--field", "displacement"]).unwrap();
        assert_eq!((o.section, o.field.as_deref()), (Some(Section { axis: 1, at: Some(2.5) }), Some("displacement")));
        assert_eq!(opts(&["s", "--out", "o", "--section=Z"]).unwrap().section, Some(Section { axis: 2, at: None }));
        let o = opts(&["--view=z", "s.json", "--size", "640x480", "--out=b.png"]).unwrap();
        assert_eq!((o.view, o.size), (StandardView::Top, (640, 480)));
        assert_eq!(opts(&["s", "--out", "o", "--view", "X"]).unwrap().view, StandardView::Right);
        assert_eq!(opts(&["s", "--out", "o", "--view", "front"]).unwrap().view, StandardView::Front);
        let o = opts(&["s", "--out", "o", "--camera", r#"{"yaw": 30, "pitch": 10, "target": [1, 2, 3], "fit": true}"#]).unwrap();
        assert_eq!(o.camera.and_then(|a| a.target), Some(Target::Point(fieldcraft_render::Vec3::new(1.0, 2.0, 3.0))));
        for bad in [
            vec![],
            vec!["s.json"],
            vec!["--out", "a.png"],
            vec!["s", "t", "--out", "o"],
            vec!["s", "--out"],
            vec!["s", "--out", "o", "--out", "p"],
            vec!["s", "--out", "o", "--view", "diagonal"],
            vec!["s", "--out", "o", "--size", "0x10"],
            vec!["s", "--out", "o", "--size", "10"],
            vec!["s", "--out", "o", "--size", "99999x10"],
            vec!["s", "--out", "o", "--size", "8192x8192"],
            vec!["s", "--out", "o", "--size", "-5x10"],
            vec!["s", "--out", "o", "--nope", "1"],
            vec!["s", "--out", "o", "--section", "w"],
            vec!["s", "--out", "o", "--section", "x=inf"],
            vec!["s", "--out", "o", "--section", "x=", "--field", "f"],
            vec!["s", "--out", "o", "--camera", "{yaw: 1}"],
            vec!["s", "--out", "o", "--camera", r#"{"yaw": 0}"#],
            vec!["s", "--out", "o", "--camera", r#"{"yaw": 0, "pitch": 0}"#, "--view", "x"],
        ] {
            assert!(opts(&bad).is_err(), "{bad:?}");
        }
    }

    /// PNG width, height and pixels (RGBA).
    fn decode(png: &[u8]) -> (u32, u32, Vec<u8>) {
        let mut r = png::Decoder::new(std::io::Cursor::new(png)).read_info().unwrap();
        let mut buf = vec![0; r.output_buffer_size().unwrap()];
        let info = r.next_frame(&mut buf).unwrap();
        buf.truncate(info.buffer_size());
        (info.width, info.height, buf)
    }

    fn script(dir: &std::path::Path, name: &str, solve: bool) -> String {
        let mut steps = fieldcraft_engine::sample::cantilever();
        if solve && let Some(a) = steps.as_array_mut() {
            a.push(serde_json::json!({"command": "solve.static", "params": {}}));
        }
        let path = dir.join(name);
        std::fs::write(&path, steps.to_string()).unwrap();
        path.to_string_lossy().into()
    }

    fn opts_for(script: String) -> Options {
        Options { script, out: String::new(), view: StandardView::Iso, camera: None, size: (320, 200), section: None, field: None }
    }

    /// The cantilever drawn headless: a PNG of the size asked for, with the model in it (pixels
    /// of the hex colour) and its element edges (dark pixels). `auto` draws on the GPU when the
    /// machine has one and on the CPU otherwise (it must be the GPU under `FIELDCRAFT_REQUIRE_GPU=1`); `cpu` always
    /// on the CPU, and both pictures show the same model. A section and a field change it.
    #[test]
    fn snapshots_the_model_and_its_mesh_headless() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-snapshot-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let o = opts_for(script(&dir, "cantilever.json", false));
        for policy in [Policy::Auto, Policy::Cpu] {
            let (png, drawn, _) = picture(&o, policy).unwrap();
            let (w, h, px) = decode(&png);
            assert_eq!((w, h), (320, 200), "{drawn:?}");
            let hex = Palette::LIGHT.hex;
            let near = |p: &[u8], c: fieldcraft_render::Rgb, d: u8| p[0].abs_diff(c.0) <= d && p[1].abs_diff(c.1) <= d && p[2].abs_diff(c.2) <= d;
            // Shaded faces vary in brightness around the hex colour; the edges are near black.
            let faces = px.chunks(4).filter(|p| near(p, hex, 60) && p[2] > p[0]).count();
            let edges = px.chunks(4).filter(|p| p[0] < 90 && p[1] < 90 && p[2] < 90).count();
            assert!(faces > 500, "{drawn:?}: {faces} face pixels");
            assert!(edges > 100, "{drawn:?}: {edges} edge pixels");
            if policy == Policy::Cpu {
                assert_eq!(drawn, Renderer::Cpu { why: None });
            } else if fieldcraft_gpu::select::gpu_required() {
                assert!(matches!(drawn, Renderer::Gpu { .. }), "FIELDCRAFT_REQUIRE_GPU=1 but drawn by {drawn:?}");
            }
        }
        // A section along the beam draws half of it.
        let whole = picture(&o, Policy::Cpu).unwrap().0;
        let cut = picture(&Options { section: Some(Section { axis: 0, at: None }), ..o.clone() }, Policy::Cpu).unwrap().0;
        assert_ne!(whole, cut);
        // A field needs a solved step; on a solved one it is the contour reported.
        assert!(picture(&Options { field: Some("von_mises".into()), ..o.clone() }, Policy::Cpu).is_err());
        let solved = Options { field: Some("displacement".into()), ..opts_for(script(&dir, "solved.json", true)) };
        let (_, _, shown) = picture(&solved, Policy::Auto).unwrap();
        assert_eq!(shown.contour.map(|c| c.0).as_deref(), Some("displacement"));
        std::fs::remove_dir_all(&dir).ok();
    }
}
