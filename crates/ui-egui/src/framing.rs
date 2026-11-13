//! Framing: when the view frames the model, and the part of the viewport it frames it in.
//!
//! One rule for every fit: the model as drawn (the deformed shape when results show one, and the
//! setup's symbols) is framed with [`fieldcraft_render::Camera::fit_in`] in the part of the
//! viewport the overlays (results legend, selection bar, orientation gizmo, navigation bar) leave
//! free, with a margin, keeping the direction of view. The view frames itself
//! ([`Reframe`]) when another model, sample or design is opened ([`fieldcraft_engine::Session::opened`];
//! the same design opened again after an edit keeps the camera) and when the first model content
//! appears (home view), and when results come up, another step or mode is shown or the viewport
//! is resized with part of the model out of the free part (a resize only when the model was
//! framed before it).
//! Fit, home and the standard views frame on request.

use egui::Rect;

/// A framing the viewport does when it next draws, once it knows its size and its overlays.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Reframe {
    /// Frame the model only if part of it is out of the free part (results coming up).
    IfOffscreen,
    /// Frame the model, keeping the direction of view.
    Fit,
    /// The home view: isometric, framed (a model opened).
    Home,
}

impl Reframe {
    /// Two requests make the stronger one.
    pub fn merge(a: Option<Reframe>, b: Reframe) -> Option<Reframe> {
        Some(a.map_or(b, |a| a.max(b)))
    }
}

/// After a resize or another step or mode shown, the frames in which the view checks again that
/// the model is not cut off: the overlays (legend, selection bar) take their new size a frame or
/// two later.
pub const RECHECK_FRAMES: u8 = 3;

/// The free part of `view`: the largest rectangle left after cutting each overlay off (largest
/// overlays first, each cut along the side that keeps the most area), as fractions
/// `[x0, y0, x1, y1]` of the view's width and height. Overlays outside the view are ignored.
pub fn free_part(view: Rect, overlays: &[Rect]) -> [f64; 4] {
    if !(view.width() > 1.0 && view.height() > 1.0) {
        return fieldcraft_render::FULL_VIEW;
    }
    let mut sorted: Vec<Rect> = overlays.iter().map(|o| o.intersect(view)).filter(|o| o.is_positive()).collect();
    sorted.sort_by(|a, b| b.area().total_cmp(&a.area()));
    let mut free = view;
    for o in sorted {
        let o = o.intersect(free);
        if !o.is_positive() {
            continue;
        }
        let cuts = [
            Rect::from_min_max(free.min, egui::pos2(o.left(), free.bottom())),
            Rect::from_min_max(egui::pos2(o.right(), free.top()), free.max),
            Rect::from_min_max(free.min, egui::pos2(free.right(), o.top())),
            Rect::from_min_max(egui::pos2(free.left(), o.bottom()), free.max),
        ];
        if let Some(best) = cuts.into_iter().filter(|c| c.is_positive()).max_by(|a, b| a.area().total_cmp(&b.area())) {
            free = best;
        }
    }
    let (w, h) = (f64::from(view.width()), f64::from(view.height()));
    [
        f64::from(free.left() - view.left()) / w,
        f64::from(free.top() - view.top()) / h,
        f64::from(free.right() - view.left()) / w,
        f64::from(free.bottom() - view.top()) / h,
    ]
}

#[cfg(test)]
mod tests {
    use egui::{pos2, vec2};

    use super::*;

    fn close(a: [f64; 4], b: [f64; 4]) -> bool {
        a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-6)
    }

    #[test]
    fn overlays_are_cut_off_along_their_cheapest_side() {
        let view = Rect::from_min_size(pos2(250.0, 126.0), vec2(1000.0, 600.0));
        assert!(close(free_part(view, &[]), fieldcraft_render::FULL_VIEW));
        // A legend column on the left: the free part is right of it.
        let legend = Rect::from_min_size(pos2(262.0, 138.0), vec2(170.0, 295.0));
        let f = free_part(view, &[legend]);
        assert!(close(f, [0.182, 0.0, 1.0, 1.0]), "{f:?}");
        // With the navigation bar (bottom centre) and the gizmo (top right) too.
        let nav = Rect::from_center_size(pos2(750.0, 704.0), vec2(220.0, 30.0));
        let gizmo = Rect::from_center_size(pos2(1182.0, 194.0), vec2(116.0, 116.0));
        let f = free_part(view, &[legend, nav, gizmo]);
        assert!(f[0] >= 0.182 - 1e-9 && f[2] <= 0.874 + 1e-6 && f[3] <= (689.0 - 126.0) / 600.0 + 1e-6, "{f:?}");
        // An overlay outside the view, or an empty view, changes nothing.
        assert!(close(free_part(view, &[Rect::from_min_size(pos2(0.0, 0.0), vec2(100.0, 100.0))]), fieldcraft_render::FULL_VIEW));
        assert!(close(free_part(Rect::NOTHING, &[legend]), fieldcraft_render::FULL_VIEW));
    }

    #[test]
    fn the_stronger_request_wins() {
        assert_eq!(Reframe::merge(None, Reframe::IfOffscreen), Some(Reframe::IfOffscreen));
        assert_eq!(Reframe::merge(Some(Reframe::Home), Reframe::IfOffscreen), Some(Reframe::Home));
        assert_eq!(Reframe::merge(Some(Reframe::Fit), Reframe::Home), Some(Reframe::Home));
    }
}
