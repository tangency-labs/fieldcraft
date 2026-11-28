//! The window's shell around the viewport, as SolidCraft's (`crates/ui-egui/src/shell.rs`): the
//! viewport fills the window below the application bar, and the workspace and ribbon cards, the
//! browser card and the History card are kit `Container` cards drawn over it. Each card records
//! where it was drawn this frame here, so the viewport keeps its overlays and its framing clear of
//! them; egui's layers keep pointer input on the card over the view.

use egui::{Rect, pos2};
use solidcraft_panel_layout::Zone;
use solidcraft_theme::Tokens;
use solidcraft_ui_kit as kit;

/// The History card's docks: along the bottom of the viewport (the default) or under the top
/// cards. The grip in its header drags the bar to the other dock (`timeline` / `panels`);
/// persisted across sessions.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum HistoryDock {
    Top,
    #[default]
    Bottom,
}

/// Zone-keyed exclusion geometry over the viewport this frame. Window state, never saved.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct State {
    /// The cards along the top (the workspace switcher and the ribbon), when they float over the
    /// viewport.
    pub cards: Vec<Rect>,
    /// The dock History stands at ([`UiState::history_dock`](crate::UiState::history_dock)).
    pub history_dock: HistoryDock,
    /// Rightmost edge of docked `Zone::Left` panels, if any.
    pub left_dock_right: Option<f32>,
    /// Leftmost edge of docked `Zone::Right` panels, if any.
    pub right_dock_left: Option<f32>,
    /// Bottommost edge of top-docked chrome (such as top-docked History), if any.
    pub top_dock_bottom: Option<f32>,
    /// Topmost edge of docked `Zone::Bottom` panels, if any.
    pub bottom_dock_top: Option<f32>,
    /// Docked persistent panels drawn inside the main window this frame.
    pub docked: Vec<Rect>,
    /// Floating persistent panels drawn inside the main window this frame.
    pub floating: Vec<Rect>,
    /// Transient right-side cards (`app.dialog`, `app.refused`) drawn this frame.
    pub transient: Vec<Rect>,
}

impl State {
    /// Forget last frame's geometry; each pass records its zones again. History stays at `dock`.
    pub fn clear(&mut self, dock: HistoryDock) {
        *self = State { history_dock: dock, ..State::default() };
    }

    /// Record one persistent panel rectangle returned by `PanelHost` in the given zone.
    pub fn set_panel_geometry(&mut self, zone: Zone, rect: Rect) {
        if !rect.is_finite() || !rect.is_positive() {
            return;
        }
        match zone {
            Zone::Left => {
                self.left_dock_right = Some(self.left_dock_right.map_or(rect.right(), |v| v.max(rect.right())));
                self.docked.push(rect);
            }
            Zone::Right => {
                self.right_dock_left = Some(self.right_dock_left.map_or(rect.left(), |v| v.min(rect.left())));
                self.docked.push(rect);
            }
            Zone::Bottom => {
                self.bottom_dock_top = Some(self.bottom_dock_top.map_or(rect.top(), |v| v.min(rect.top())));
                self.docked.push(rect);
            }
            Zone::Floating => {
                self.floating.push(rect);
            }
            Zone::Window => {}
        }
    }

    /// Record a top-docked strip (such as `HistoryDock::Top`) above the side columns.
    pub fn record_top_dock(&mut self, rect: Rect) {
        if !rect.is_finite() || !rect.is_positive() {
            return;
        }
        self.top_dock_bottom = Some(self.top_dock_bottom.map_or(rect.bottom(), |v| v.max(rect.bottom())));
        self.docked.push(rect);
    }

    /// Record a transient card (`app.dialog` or `app.refused`) drawn above the viewport.
    pub fn add_transient(&mut self, rect: Rect) {
        if rect.is_finite() && rect.is_positive() {
            self.transient.push(rect);
        }
    }

    /// Every surface over the viewport this frame: the viewport's overlays and framing keep clear.
    pub fn covered(&self) -> impl Iterator<Item = Rect> + '_ {
        self.cards.iter().copied().chain(self.docked.iter().copied()).chain(self.floating.iter().copied()).chain(self.transient.iter().copied())
    }

    /// The vertical slice of `available` below top cards / top-docked history and above bottom-docked history.
    pub fn panel_space(&self, mut available: Rect) -> Rect {
        for card in &self.cards {
            available.min.y = kit::within(card.bottom(), available.top(), available.bottom());
        }
        if let Some(top) = self.top_dock_bottom {
            available.min.y = kit::within(top, available.top(), available.bottom());
        }
        if let Some(bottom) = self.bottom_dock_top {
            available.max.y = kit::within(bottom, available.top(), available.bottom());
        }
        available
    }

    /// The horizontal room `(left_edge, right_edge)` left between left and right docked panels in `available`.
    pub fn side_bounds(&self, available: Rect) -> (f32, f32) {
        let left = self.left_dock_right.map_or(available.left(), |v| kit::within(v, available.left(), available.right()));
        let right = self.right_dock_left.map_or(available.right(), |v| kit::within(v, left, available.right()));
        (left, right)
    }
}

/// The viewport's top corners are rounded concentrically with the cards: the card margin plus
/// the card radius (SolidCraft `shell::corner_radius`).
fn corner_radius() -> f32 {
    Tokens::outset_radius(f32::from(Tokens::get().radius_lg), kit::SHELL_GAP)
}

/// Fill the two concave corners where the application bar flows into the viewport's top edge, so
/// the viewport reads as a rounded surface under the bar (SolidCraft `shell::paint_top_corners`).
pub fn paint_top_corners(painter: &egui::Painter, viewport: Rect, color: egui::Color32) {
    let r = corner_radius().min(viewport.width() / 2.0).min(viewport.height() / 2.0);
    if !(r.is_finite() && r > 0.5) {
        return;
    }
    const STEPS: usize = 32;
    for left in [true, false] {
        let (corner, centre, dir) = if left {
            (viewport.left_top(), pos2(viewport.left() + r, viewport.top() + r), -1.0)
        } else {
            (viewport.right_top(), pos2(viewport.right() - r, viewport.top() + r), 1.0)
        };
        let mut mesh = egui::Mesh::default();
        mesh.colored_vertex(corner, color);
        for i in 0..=STEPS {
            // From the top tangent point round to the side tangent point.
            let a = std::f32::consts::FRAC_PI_2 * i as f32 / STEPS as f32;
            mesh.colored_vertex(pos2(centre.x + dir * r * a.sin(), centre.y - r * a.cos()), color);
        }
        for i in 0..STEPS as u32 {
            mesh.add_triangle(0, i + 1, i + 2);
        }
        painter.add(egui::Shape::mesh(mesh));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui::vec2;

    #[test]
    fn the_panel_space_is_below_the_top_cards_and_above_history() {
        let full = Rect::from_min_size(pos2(0.0, 40.0), vec2(1400.0, 860.0));
        let mut s = State::default();
        assert_eq!(s.panel_space(full), full);
        s.cards = vec![Rect::from_min_size(pos2(7.0, 47.0), vec2(90.0, 110.0)), Rect::from_min_size(pos2(104.0, 47.0), vec2(1289.0, 110.0))];
        s.set_panel_geometry(Zone::Bottom, Rect::from_min_max(pos2(7.0, 853.0), pos2(1393.0, 893.0)));
        let space = s.panel_space(full);
        assert_eq!((space.top(), space.bottom()), (157.0, 853.0));
        assert_eq!(s.covered().count(), 3);
        // A card squeezed past the window gives an empty space, never a panic.
        s.cards = vec![Rect::from_min_size(pos2(0.0, 2000.0), vec2(10.0, 10.0))];
        assert_eq!(s.panel_space(full).height(), 0.0);
    }

    #[test]
    fn with_history_at_the_top_the_panel_space_is_below_it() {
        let full = Rect::from_min_size(pos2(0.0, 40.0), vec2(1400.0, 860.0));
        let mut s = State::default();
        s.clear(HistoryDock::Top);
        s.cards = vec![Rect::from_min_size(pos2(104.0, 47.0), vec2(1289.0, 110.0))];
        s.record_top_dock(Rect::from_min_max(pos2(7.0, 164.0), pos2(1393.0, 204.0)));
        let space = s.panel_space(full);
        assert_eq!((space.top(), space.bottom()), (204.0, 900.0));
    }
}
