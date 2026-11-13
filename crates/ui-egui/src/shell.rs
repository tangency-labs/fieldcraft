//! The window's shell around the viewport, as SolidCraft's (`crates/ui-egui/src/shell.rs`): the
//! viewport fills the window below the application bar, and the workspace and ribbon cards, the
//! browser card and the History card are kit `Container` cards drawn over it. Each card records
//! where it was drawn this frame here, so the viewport keeps its overlays and its framing clear of
//! them; egui's layers keep pointer input on the card over the view.

use egui::{Rect, pos2};
use solidcraft_theme::Tokens;
use solidcraft_ui_kit as kit;

/// Where the cards over the viewport were drawn this frame. Window state, never saved.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct State {
    /// The cards along the top (the workspace switcher and the ribbon), when they float over the
    /// viewport.
    pub cards: Vec<Rect>,
    /// The History card along the bottom, when it floats over the viewport.
    pub history: Option<Rect>,
    /// The browser card, while it is shown.
    pub browser: Option<Rect>,
}

impl State {
    /// Forget last frame's cards; each draws itself again and records where.
    pub fn clear(&mut self) {
        *self = State::default();
    }

    /// Every card over the viewport this frame: the viewport's overlays and framing keep clear.
    pub fn covered(&self) -> impl Iterator<Item = Rect> + '_ {
        self.cards.iter().copied().chain(self.history).chain(self.browser)
    }

    /// The part of `available` the browser card may use: below the top cards, above History.
    pub fn panel_space(&self, mut available: Rect) -> Rect {
        for card in &self.cards {
            available.min.y = kit::within(card.bottom(), available.top(), available.bottom());
        }
        if let Some(history) = self.history.filter(|h| h.intersects(available)) {
            available.max.y = kit::within(history.top(), available.top(), available.bottom());
        }
        available
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
        s.history = Some(Rect::from_min_max(pos2(7.0, 853.0), pos2(1393.0, 893.0)));
        let space = s.panel_space(full);
        assert_eq!((space.top(), space.bottom()), (157.0, 853.0));
        assert_eq!(s.covered().count(), 3);
        // A card squeezed past the window gives an empty space, never a panic.
        s.cards = vec![Rect::from_min_size(pos2(0.0, 2000.0), vec2(10.0, 10.0))];
        assert_eq!(s.panel_space(full).height(), 0.0);
    }
}
