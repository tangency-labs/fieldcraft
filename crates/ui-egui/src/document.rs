//! Another document replaces the open one: nothing of the old model's interface survives.
//!
//! [`ModelApp::document_replaced`] is the one place that says, for every field of [`ModelApp`],
//! whether it belongs to the document (reset), to the window or the user (kept), or is a cache
//! keyed by what it was built for (kept: a new document never matches its key). It takes the
//! struct apart without `..`, so a field added later does not compile until it is classified
//! here. The trigger is the engine's document identity ([`fieldcraft_engine::Session::opened`]): every path that
//! replaces the document (File › New, Open, the start page, a sample, a dropped file, Read Deck,
//! Recover, and `file.*` from scripts, MCP and the control channel) changes it.

use fieldcraft_engine::catalogue::{self, Tab};
use fieldcraft_engine::progress::{StageState, Status};

use crate::files::Prompt;
use crate::{ModelApp, PickState, UiState, Viewport};

/// The browser sections unfolded on a fresh document.
pub(crate) const TREE_OPEN: [&str; 2] = ["mesh", "results"];

/// The tab where the document's work is: that of the first stage of its setup not yet done
/// (an empty model's geometry, a meshed one's material, a model ready to solve its Solve), or
/// Results once every stage is. `None` for a stage whose command is on the Tools tab.
pub(crate) fn stage_tab(progress: &[StageState]) -> Option<Tab> {
    let at = progress.iter().find(|s| !matches!(s.status, Status::Done | Status::Optional)).or(progress.last())?;
    catalogue::stage_tab(at.command)
}

impl ModelApp {
    /// The open document was replaced ([`fieldcraft_engine::Session::opened`] changed): reset everything of the
    /// old one's interface and show the workspace the new one's stage needs.
    pub(crate) fn document_replaced(&mut self) {
        let tab = stage_tab(&self.session.progress());
        let results_shown = self.session.results_view.is_some();
        let ModelApp {
            // The document itself, replaced by the command.
            session: _,
            ui,
            // The camera is framed home on the new document ([`ModelApp::frame_new_document`]).
            cam: _,
            // The status line says what replaced the document.
            status: _,
            dialog,
            preview,
            merge_preview,
            relax_preview,
            fluid_note,
            node_move,
            node_drag,
            // Numbers every form opening: never reused.
            forms_opened: _,
            viewport,
            // Window: last frame's layout, devices, counters.
            shell: _,
            toolbar_rects: _,
            // Whether this frame's input came from the control channel.
            driven: _,
            // The title bar's offer span as drawn last frame.
            offer_drawn: _,
            gpu_info: _,
            #[cfg(not(target_arch = "wasm32"))]
                capturer: _,
            revision: _,
            // Caches keyed by the model revision or design they were built for.
            display: _,
            glyphs: _,
            plot: _,
            probe,
            // The units the probe was read in, followed by [`ModelApp::follow_units`].
            shown_units: _,
            reactions,
            xy,
            path_editor,
            reply,
            features_open,
            results_shown: shown,
            shape_shown,
            // The framing the new document asks for ([`ModelApp::frame_new_document`]).
            reframe: _,
            // The hosts count camera requests across documents.
            camera_seen: _,
            // The document identity the view follows.
            framed_document: _,
            // Caches keyed by the design, revision or step they were built for.
            bodies: _,
            prebuilt: _,
            prepared: _,
            scene: _,
            symbols: _,
            step_glyphs: _,
            // Window and host: theme, control channel, input, files, services, jobs.
            styled: _,
            synthetic: _,
            synthetic_mods: _,
            synthetic_drops: _,
            control_rx: _,
            shots: _,
            offscreen: _,
            now: _,
            quit_requested: _,
            services: _,
            picker: _,
            recent: _,
            start_page: _,
            prefs: _,
            // The Preferences window: the user's settings, not the document's.
            settings: _,
            usage: _,
            consent_asked: _,
            updates: _,
            reload: _,
            setup: _,
            loading: _,
            refused,
            jobs: _,
            queued: _,
            prompt,
            report: _,
            last_error: _,
            held,
            // Keyed by the document identity it marks.
            imported_document: _,
            custom_titlebar: _,
            integrated_titlebar: _,
            palette_text: _,
            tree_open,
            timeline_open,
            selected,
            selected_seed,
            selected_mesh,
            pick,
            capturing: _,
            // The navigation bar's mode, the user's.
            nav: _,
            layout: _,
            panel_rects: _,
            timeline,
            // Keyed by the model revision they were read at.
            element_counts: _,
            properties: _,
            mesh_made: _,
            search,
        } = self;
        ui.document_replaced(tab);
        *dialog = None;
        *preview = Default::default();
        *merge_preview = Default::default();
        *relax_preview = Default::default();
        *fluid_note = Default::default();
        *node_move = Default::default();
        *node_drag = None;
        viewport.document_replaced();
        *probe = None;
        *reactions = None;
        *xy = None;
        *path_editor = None;
        *reply = None;
        features_open.clear();
        // The new document's stage chose the workspace: results it opens with bring no tab.
        *shown = results_shown;
        *shape_shown = None;
        *refused = None;
        // A question about the old document (its unsaved changes, what a deletion or a remesh
        // of its items would do, a design over its model) has nothing left to answer for.
        if matches!(prompt, Some(Prompt::Unsaved(_) | Prompt::ReplaceMesh(_) | Prompt::Delete(_) | Prompt::Design(_))) {
            *prompt = None;
            *held = None;
        }
        *tree_open = TREE_OPEN.into_iter().map(str::to_string).collect();
        timeline_open.clear();
        *timeline = None;
        *selected = None;
        *selected_seed = None;
        *selected_mesh = false;
        pick.document_replaced();
        // Hits may name the old document's items.
        *search = None;
    }
}

impl UiState {
    /// What of the window shows the old document: the workspace goes to the new one's stage
    /// (`tab`) and an open menu closes. The display choices and the panels are the user's.
    fn document_replaced(&mut self, tab: Option<Tab>) {
        let UiState {
            style: _,
            show_nodes: _,
            show_geometry: _,
            show_mesh: _,
            show_loads: _,
            show_interactions: _,
            result_symbols: _,
            result_extremes: _,
            symbol_size: _,
            perspective: _,
            dark: _,
            show_tree: _,
            browser_width: _,
            browser_collapsed: _,
            browser_float: _,
            browser_eye_left: _,
            show_timeline: _,
            history_dock: _,
            status_bar: _,
            nav_preset: _,
            reverse_zoom: _,
            workspace: _,
            tab: _,
            // The command search is the window's.
            palette: _,
            menu,
        } = self;
        *menu = None;
        if let Some(tab) = tab {
            self.select_tab(tab);
        }
    }
}

impl Viewport {
    /// The context menu's item and the orbit's pivot were the old document's; the renderer and
    /// its uploads are the window's (keyed by the scene).
    fn document_replaced(&mut self) {
        let Viewport {
            rect: _,
            gpu: _,
            slot: _,
            uploaded: _,
            overlay_slot: _,
            overlay_uploaded: _,
            symbols_slot: _,
            symbols_uploaded: _,
            cpu: _,
            cpu_failed: _,
            gpu_errors: _,
            overlays: _,
            recheck: _,
            // A drag in progress keeps its gesture; what it turns about goes.
            drag: _,
            orbit_pivot,
            menu_target,
        } = self;
        *orbit_pivot = None;
        *menu_target = None;
    }
}

impl PickState {
    /// Picks and hovers are by name, and the names of another document are other things. The
    /// filter and the box mode are the user's; the filter an open form set is restored when the
    /// form closes ([`ModelApp::follow_form_filter`]).
    fn document_replaced(&mut self) {
        let PickState {
            filter: _,
            selection,
            hover,
            feature_hover,
            box_mode: _,
            dragging,
            hover_at: _,
            form_hover,
            generation,
            hover_generation,
            // Keyed by the model revision they were built at.
            index: _,
            names: _,
            selected_overlay: _,
            form_filter: _,
            overlay: _,
            hover_key,
            checked: _,
        } = self;
        if !selection.is_empty() {
            selection.clear();
            *generation = generation.wrapping_add(1);
        }
        if hover.is_some() || !feature_hover.is_empty() {
            *hover = None;
            feature_hover.clear();
            *hover_generation = hover_generation.wrapping_add(1);
        }
        *dragging = None;
        *form_hover = false;
        *hover_key = None;
    }
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;
    use fieldcraft_engine::catalogue::{Tab, Workspace};
    use serde_json::json;

    use crate::ModelApp;
    use crate::browser::Def;

    /// The cantilever solved, its load selected and hovered, a form open, von Mises shown with
    /// a fixed legend, a probe, a History run open, a browser section folded, the Results
    /// workspace up and a refusal listed.
    fn worked_on() -> ModelApp {
        let mut app = ModelApp::new(Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.select_definition(Some((Def::Load, "Push".into())));
        app.run("results.contour", json!({"step": "Load", "field": "von_mises", "range": {"min": 0, "max": 80}})).unwrap();
        app.ui.select_workspace(Workspace::Results);
        app.open_form("model.material.create", json!({}));
        app.pick.hover = Some(crate::selection::Item::Node(1));
        app.timeline_open.insert(0);
        app.tree_open.remove("mesh");
        app.features_open.insert(("Body".into(), 1));
        app.refused = fieldcraft_engine::find_command("results.xy").and_then(|spec| Session::default().availability(spec).err());
        app.ui.menu = Some("view".into());
        assert!(app.dialog.is_some() && app.selected.is_some() && app.session.results_view.is_some() && app.refused.is_some());
        app
    }

    fn assert_fresh(app: &ModelApp, opened: &str) {
        assert!(app.selected.is_none(), "{opened}: the old load is still selected");
        assert!(app.pick.selection.is_empty() && app.pick.hover.is_none(), "{opened}: the old pick or hover survives");
        assert!(app.dialog.is_none(), "{opened}: the old form is still open");
        assert!(app.refused.is_none(), "{opened}: the old refusal is still listed");
        assert!(app.reply.is_none() && app.probe.is_none() && app.xy.is_none(), "{opened}: an old reply, probe or plot survives");
        assert!(app.timeline_open.is_empty(), "{opened}: the old History run is still open");
        assert!(app.features_open.is_empty(), "{opened}: old Feature Graph rows are still open");
        assert!(app.tree_open.contains("mesh"), "{opened}: the browser keeps the old document's folding");
        assert!(app.ui.menu.is_none(), "{opened}: a menu stays open over the new document");
    }

    /// File › New after work on the cantilever: an empty model, on the Geometry tab where its
    /// first stage is, with nothing of the cantilever's interface left.
    #[test]
    fn a_new_model_leaves_nothing_of_the_old_one() {
        let mut app = worked_on();
        app.run("file.new", json!({})).unwrap();
        assert!(app.session.results_view.is_none());
        assert_fresh(&app, "file.new");
        assert_eq!((app.ui.workspace, app.ui.tab), (Workspace::Prepare, Tab::Geometry), "an empty model starts at its geometry");
    }

    /// File › New on unsaved work asks first and resets nothing until it is answered; once the
    /// work is let go, the new model starts afresh.
    #[test]
    fn new_on_unsaved_work_asks_before_it_resets() {
        let mut app = worked_on();
        app.run("model.reference_point.create", json!({"name": "RP", "xyz": [0, 0, 0]})).unwrap();
        assert!(app.has_unsaved_work(), "the tab's dot");
        let journal = app.session.journal.len();
        app.new_model();
        assert_eq!(app.prompt, Some(crate::files::Prompt::Unsaved(crate::files::Replace::New)), "it asks to save first");
        assert_eq!(app.session.journal.len(), journal, "the model is still there");
        assert!(app.selected.is_some() && app.session.results_view.is_some(), "nothing is reset before the answer");
        app.answer_unsaved(crate::files::Answer::Cancel);
        assert_eq!(app.session.journal.len(), journal, "Cancel keeps the model");
        assert!(app.selected.is_some(), "and its interface");
        app.new_model();
        app.answer_unsaved(crate::files::Answer::Discard);
        assert_fresh(&app, "New after Don't Save");
    }

    /// Another sample opened over the worked-on cantilever: the same reset, and the workspace of
    /// the stage it is at (the impact sample is ready to solve).
    #[test]
    fn another_sample_leaves_nothing_of_the_old_one() {
        let mut app = worked_on();
        app.run("file.sample", json!({"name": "impact"})).unwrap();
        assert_fresh(&app, "file.sample");
        assert_eq!((app.ui.workspace, app.ui.tab), (Workspace::Analysis, Tab::Solve), "a sample ready to solve opens at Solve");
    }
}
