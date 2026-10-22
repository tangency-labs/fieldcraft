//! Saved result views, kept in the manifest's `view` ([`crate::View`]): view state,
//! outside the model's definition (and its hash, [`crate::results::definition_hash`]).
//!
//! A saved view ([`SavedView`]) is a camera and what the results display shows, each in the
//! parameter shape of the command that sets it: `camera` is `ui.camera`'s (and the MCP
//! `screenshot`'s), `results.contour` is `results.contour`'s and so on, so a view is those
//! commands' parameters, applied in turn. `units` names the unit system their bare numbers are
//! in (as a journal's header does). The engine checks every parameter when it applies them, as the commands do;
//! this module checks what the file itself must hold: counts, names, references and ranges.

use std::collections::BTreeSet;

use fieldcraft_model::DisplaySystem;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

use crate::{Error, Result, invalid};

/// Most saved views.
pub const MAX_VIEWS: usize = 1000;

/// Longest name of a view.
pub const MAX_VIEW_NAME: usize = 64;

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

/// The saved views and the current one checked together: names unique, the current view among
/// the views.
pub(crate) fn check(views: &[SavedView], current: Option<&str>) -> Result<()> {
    if views.len() > MAX_VIEWS {
        return Err(Error::Limit(format!("a model keeps at most {MAX_VIEWS} views")));
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
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn views_are_named_once_and_the_current_one_is_saved() {
        let view = |name: &str| SavedView { name: name.into(), units: DisplaySystem::default(), camera: None, results: None };
        assert!(check(&[view("A"), view("B")], Some("B")).is_ok());
        assert!(check(&[view("A"), view("A")], None).is_err());
        assert!(check(&[view("A")], Some("C")).is_err());
        assert!(check(&[view(" ")], None).is_err());
        assert!(check(&[view("a\nb")], None).is_err());
        let many: Vec<SavedView> = (0..=MAX_VIEWS).map(|i| view(&format!("V{i}"))).collect();
        assert!(matches!(check(&many, None), Err(Error::Limit(_))));
    }
}
