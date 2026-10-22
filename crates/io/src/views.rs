//! Saved result views and renders, kept in the manifest's `view` ([`crate::View`]): view state,
//! outside the model's definition (and its hash, [`crate::results::definition_hash`]).
//!
//! A saved view ([`SavedView`]) is a camera and what the results display shows, each in the
//! parameter shape of the command that sets it: `camera` is `ui.camera`'s (and the MCP
//! `screenshot`'s), `results.contour` is `results.contour`'s and so on, so a view is those
//! commands' parameters, applied in turn. `units` names the unit system their bare numbers are
//! in (as a journal's header does). A render ([`Render`]) says how a saved view becomes a clip:
//! which step and frames, the frame rate and holds, the picture size, a caption, the legend and
//! the time label. The engine checks every parameter when it applies them, as the commands do;
//! this module checks what the file itself must hold: counts, names, references and ranges.

use std::collections::BTreeSet;

use fieldcraft_model::DisplaySystem;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

use crate::{Error, Result, invalid};

/// Most saved views, and most renders.
pub const MAX_VIEWS: usize = 1000;

/// Longest name of a view or render, and longest caption.
pub const MAX_VIEW_NAME: usize = 64;
pub const MAX_CAPTION: usize = 200;

/// The largest frame rate, hold and side of a render.
pub const MAX_FPS: u32 = 120;
pub const MAX_HOLD: u32 = 10_000;
pub const MIN_SIDE: u32 = 16;
pub const MAX_SIDE: u32 = 4096;

/// The parameters of a command, as an object.
pub type Params = Map<String, Value>;

/// What the results display shows in a saved view: the step, then the parameters of each
/// `results.*` command that sets part of it (without their `step`), applied in this order.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResultsSetup {
    pub step: String,
    /// `results.deformed`'s parameters (`show`, `scale`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub deformed: Option<Params>,
    /// `results.modes`' `mode` (frequency and buckling steps).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mode: Option<u64>,
    /// `results.frame`'s `frame` (steps with frames).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frame: Option<u64>,
    /// `results.contour`'s parameters (`field`, `bands`, `range`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contour: Option<Params>,
    /// `results.slice`'s parameters (fluid steps: `field`, `axis`, `position`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub slice: Option<Params>,
    /// `results.surface`'s parameters (fluid steps: `show`, `field`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface: Option<Params>,
    /// `results.streamlines`' parameters (fluid steps: `show`, `count`, `max_steps`, `seeds`, `line`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub streamlines: Option<Params>,
}

/// A named view: a camera and what the results display shows.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SavedView {
    pub name: String,
    /// The unit system the bare numbers of `camera` and `results` are in.
    pub units: DisplaySystem,
    /// `ui.camera`'s parameters (`yaw`, `pitch`, `distance`, `target`, `fit`, `fov`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub camera: Option<Params>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub results: Option<ResultsSetup>,
}

/// Which frames a render draws: every one, or those from `from` to `to` (inclusive, from 0).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Frames {
    #[default]
    #[serde(with = "all")]
    All,
    Range([u64; 2]),
}

/// `"all"`, the only name of [`Frames::All`].
mod all {
    use serde::{Deserialize, Deserializer, Serializer};

    pub fn serialize<S: Serializer>(s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str("all")
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<(), D::Error> {
        let s = String::deserialize(d)?;
        if s == "all" { Ok(()) } else { Err(serde::de::Error::custom("frames must be \"all\" or [from, to]")) }
    }
}

/// How long each picture of a render shows, in frames of the clip.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hold {
    /// Every picture but the last.
    pub each: u32,
    /// The last picture (the settled state).
    pub last: u32,
}

impl Default for Hold {
    fn default() -> Self {
        Hold { each: 1, last: 1 }
    }
}

/// What a render writes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum RenderFormat {
    /// Numbered PNG pictures in a folder, with the ffmpeg command line that encodes them.
    #[default]
    #[serde(rename = "png-sequence")]
    PngSequence,
}

/// How a saved view becomes a clip (`results.export.animation`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Render {
    pub name: String,
    /// The saved view it starts from.
    pub view: String,
    /// The step whose frames it plays: the view's step when left out.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub step: Option<String>,
    pub frames: Frames,
    /// Frames per second of the clip.
    pub fps: u32,
    pub hold: Hold,
    /// Width and height in pixels.
    pub size: [u32; 2],
    /// A line of text under the picture.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub caption: Option<String>,
    /// The contour legend under the picture.
    pub legend: bool,
    /// The frame's time (`t = 0.25 s`) under the picture.
    pub time_label: bool,
    pub format: RenderFormat,
}

impl Default for Render {
    fn default() -> Self {
        Render {
            name: String::new(),
            view: String::new(),
            step: None,
            frames: Frames::All,
            fps: 12,
            hold: Hold::default(),
            size: [1280, 720],
            caption: None,
            legend: true,
            time_label: true,
            format: RenderFormat::PngSequence,
        }
    }
}

/// A view or render name: 1 to [`MAX_VIEW_NAME`] characters, no control characters, not only
/// spaces.
pub fn check_view_name(what: &str, name: &str) -> Result<()> {
    if name.trim().is_empty() || name.chars().count() > MAX_VIEW_NAME || name.chars().any(char::is_control) {
        return invalid(format!(
            "a {what} name is 1 to {MAX_VIEW_NAME} characters with no control characters (`{}`)",
            name.chars().filter(|c| !c.is_control()).take(MAX_VIEW_NAME).collect::<String>()
        ));
    }
    Ok(())
}

impl Render {
    /// Whether its numbers are in range and its texts are names and captions.
    pub fn check(&self) -> Result<()> {
        check_view_name("render", &self.name)?;
        check_view_name("view", &self.view)?;
        if let Some(step) = &self.step {
            fieldcraft_model::check_name("step", step).map_err(|e| Error::Invalid(e.to_string()))?;
        }
        if let Frames::Range([from, to]) = self.frames
            && from > to
        {
            return invalid(format!("render `{}`: frames [{from}, {to}] run backwards", self.name));
        }
        if !(1..=MAX_FPS).contains(&self.fps) {
            return invalid(format!("render `{}`: fps is 1 to {MAX_FPS}", self.name));
        }
        if !(1..=MAX_HOLD).contains(&self.hold.each) || !(1..=MAX_HOLD).contains(&self.hold.last) {
            return invalid(format!("render `{}`: hold.each and hold.last are 1 to {MAX_HOLD} frames", self.name));
        }
        if self.size.iter().any(|s| !(MIN_SIDE..=MAX_SIDE).contains(s)) {
            return invalid(format!("render `{}`: size is [width, height], each {MIN_SIDE} to {MAX_SIDE} pixels", self.name));
        }
        if let Some(c) = &self.caption
            && (c.chars().count() > MAX_CAPTION || c.chars().any(char::is_control))
        {
            return invalid(format!("render `{}`: a caption is at most {MAX_CAPTION} characters on one line", self.name));
        }
        Ok(())
    }
}

/// The saved views, the current one and the renders checked together: names unique, the current
/// view and every render's view among the views.
pub(crate) fn check(views: &[SavedView], current: Option<&str>, renders: &[Render]) -> Result<()> {
    if views.len() > MAX_VIEWS || renders.len() > MAX_VIEWS {
        return Err(Error::Limit(format!("a model keeps at most {MAX_VIEWS} views and {MAX_VIEWS} renders")));
    }
    let mut names = BTreeSet::new();
    for v in views {
        check_view_name("view", &v.name)?;
        if !names.insert(v.name.as_str()) {
            return invalid(format!("view `{}` is listed twice", v.name));
        }
        if let Some(r) = &v.results {
            fieldcraft_model::check_name("step", &r.step).map_err(|e| Error::Invalid(format!("view `{}`: {e}", v.name)))?;
        }
    }
    if let Some(c) = current
        && !names.contains(c)
    {
        return invalid(format!("the current view `{}` is not among the saved views", c.chars().take(MAX_VIEW_NAME).collect::<String>()));
    }
    let mut render_names = BTreeSet::new();
    for r in renders {
        r.check()?;
        if !render_names.insert(r.name.as_str()) {
            return invalid(format!("render `{}` is listed twice", r.name));
        }
        if !names.contains(r.view.as_str()) {
            return invalid(format!("render `{}` starts from view `{}`, which is not saved", r.name, r.view));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    #[test]
    fn a_render_reads_with_its_defaults_and_writes_back() {
        let r: Render = serde_json::from_value(json!({"name": "Airflow clip", "view": "Cutaway"})).unwrap();
        assert_eq!((r.fps, r.hold, r.size, r.frames, r.legend), (12, Hold { each: 1, last: 1 }, [1280, 720], Frames::All, true));
        let v = serde_json::to_value(&r).unwrap();
        assert_eq!(v["frames"], "all");
        assert_eq!(v["format"], "png-sequence");
        assert_eq!(serde_json::from_value::<Render>(v).unwrap(), r);
        let r: Render = serde_json::from_value(json!({"name": "R", "view": "V", "frames": [2, 9], "hold": {"each": 2, "last": 16}})).unwrap();
        assert_eq!(r.frames, Frames::Range([2, 9]));
        for bad in [
            json!({"name": "R", "view": "V", "frames": "some"}),
            json!({"name": "R", "view": "V", "format": "mp4"}),
            json!({"name": "R", "view": "V", "speed": 2}),
            json!({"name": "R", "view": "V", "hold": {"each": 1}}),
        ] {
            assert!(serde_json::from_value::<Render>(bad.clone()).is_err(), "{bad}");
        }
        for bad in [
            Render { fps: 0, ..r.clone() },
            Render { fps: 1000, ..r.clone() },
            Render { size: [8, 720], ..r.clone() },
            Render { hold: Hold { each: 0, last: 1 }, ..r.clone() },
            Render { frames: Frames::Range([5, 1]), ..r.clone() },
            Render { caption: Some("a\nb".into()), ..r.clone() },
            Render { name: " ".into(), ..r.clone() },
            Render { step: Some("1 bad".into()), ..r.clone() },
        ] {
            assert!(bad.check().is_err(), "{bad:?}");
        }
    }

    #[test]
    fn views_and_renders_refer_to_each_other() {
        let view = |name: &str| SavedView { name: name.into(), units: DisplaySystem::default(), camera: None, results: None };
        let render = |name: &str, view: &str| Render { name: name.into(), view: view.into(), ..Render::default() };
        assert!(check(&[view("A"), view("B")], Some("B"), &[render("R", "A")]).is_ok());
        assert!(check(&[view("A"), view("A")], None, &[]).is_err());
        assert!(check(&[view("A")], Some("C"), &[]).is_err());
        assert!(check(&[view("A")], None, &[render("R", "C")]).is_err());
        assert!(check(&[view("A")], None, &[render("R", "A"), render("R", "A")]).is_err());
        assert!(check(&[view(" ")], None, &[]).is_err());
        assert!(check(&[view("a\nb")], None, &[]).is_err());
        let many: Vec<SavedView> = (0..=MAX_VIEWS).map(|i| view(&format!("V{i}"))).collect();
        assert!(matches!(check(&many, None, &[]), Err(Error::Limit(_))));
    }
}
