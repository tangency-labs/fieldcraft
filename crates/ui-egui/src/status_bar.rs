//! The status bar (Settings › Interface › Status bar; off by default), after SolidCraft's
//! `crates/ui-egui/src/status_bar.rs` (MIT OR Apache-2.0, see ATTRIBUTION.md): a thin bar along
//! the window's bottom edge, below the History card, with the model's size, the projection and
//! the graphics device, right-aligned, and an update on offer before them. Status messages are
//! the view's toast ([`crate::viewport`]), not the bar's.
//!
//! The bar is window chrome like the application bar: full width and flat. It is allocated before
//! the viewport, so the view and everything anchored to its bottom (History, the navigation bar,
//! the toast) sit above it.

use egui::RichText;
use solidcraft_ui_kit as kit;

use crate::ModelApp;
use crate::theme::Tokens;

/// The bar's height (points), SolidCraft's.
pub const HEIGHT: f32 = 24.0;

/// The bar is drawn: the setting is on.
pub fn shown(app: &ModelApp) -> bool {
    app.ui.status_bar
}

/// What the bar says, right to left: the graphics device, the projection and the model's size.
pub fn items(app: &ModelApp) -> [String; 3] {
    let (nodes, elements) = (app.session.model.nodes().len(), app.session.model.elements().len());
    [
        app.gpu_info.clone().unwrap_or_else(|| "CPU renderer".into()),
        (if app.ui.perspective { "Perspective" } else { "Orthographic" }).into(),
        format!("{nodes} nodes  {elements} elements"),
    ]
}

/// Draw the bar along the window's bottom edge. Call it before the shell and the viewport take
/// the rest of the window: the bar then spans the whole width and sits below History.
pub fn show(app: &mut ModelApp, ui: &mut egui::Ui) {
    if !shown(app) {
        return;
    }
    let t = Tokens::get(app.ui.dark);
    let k = solidcraft_theme::Tokens::get();
    let host = kit::Container::new(ui.ctx(), kit::Surface::ApplicationBar);
    let panel = host.edge("fc_status_bar", kit::Edge::Bottom).exact_size(HEIGHT).show_separator_line(false);
    let frame = host.edge_frame(t.toolbar);
    host.show_edge(ui, panel, frame, |ui| {
        let bar = ui.max_rect();
        app.place("status_bar", bar);
        ui.painter().hline(bar.x_range(), bar.top() + 0.5, k.primary_border(k.border));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_space(k.space_6);
            let item = |text: &str| RichText::new(text).size(k.font_body).color(t.text_dim);
            for (i, text) in items(app).iter().enumerate() {
                if i > 0 {
                    rule(ui);
                }
                ui.label(item(text));
            }
            crate::updates::notice(app, ui);
        });
    });
}

/// A rule between two of the bar's items, as SolidCraft's status bar divides its own.
pub(crate) fn rule(ui: &mut egui::Ui) {
    let k = solidcraft_theme::Tokens::get();
    let (r, _) = ui.allocate_exact_size(egui::vec2(1.0, ui.available_height()), egui::Sense::hover());
    k.secondary_vline(ui.painter(), r.center().x, (r.top() + k.space_4)..=(r.bottom() - k.space_4), k.border);
}

#[cfg(test)]
mod tests {
    use egui::{Pos2, Rect, vec2};
    use fieldcraft_engine::Session;

    use super::*;

    fn draw(app: &mut ModelApp, ctx: &egui::Context) {
        let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, vec2(1400.0, 900.0))), ..Default::default() };
        let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
        out.textures_delta.clear();
    }

    /// Off by default; on, it runs along the window's bottom edge, under the History card, and
    /// says the model's size, the projection and the device.
    #[test]
    fn the_status_bar_is_off_by_default_and_sits_below_history() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        for _ in 0..2 {
            draw(&mut app, &ctx);
        }
        assert!(!app.ui.status_bar && !app.layout.contains_key("status_bar"));
        app.ui.status_bar = true;
        for _ in 0..2 {
            draw(&mut app, &ctx);
        }
        let bar = app.layout.get("status_bar").copied().unwrap_or(Rect::NOTHING);
        let history = app.layout.get("timeline").copied().unwrap_or(Rect::NOTHING);
        assert!((bar.bottom() - 900.0).abs() < 0.5 && (bar.height() - HEIGHT).abs() < 0.5 && bar.width() > 1399.0, "{bar:?}");
        assert!(history.bottom() < bar.top(), "{history:?} {bar:?}");
        assert_eq!(items(&app), ["CPU renderer".to_string(), "Orthographic".into(), "0 nodes  0 elements".into()]);
    }
}
