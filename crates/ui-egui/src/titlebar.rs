//! The window's own title bar: the kit's window chrome (`solidcraft_ui_kit::window_chrome`) with
//! FieldCraft's answers. On Windows and Linux the window has no OS decorations
//! ([`ModelApp::custom_titlebar`]): empty space in the application bar drags the window and a
//! double-click there maximizes or restores it, caption buttons at the far right minimize,
//! maximize/restore and close, and thin invisible zones along the window's edges resize it.
//! macOS keeps its traffic lights over the bar ([`ModelApp::integrated_titlebar`]).

pub use solidcraft_ui_kit::window_chrome::{WindowAction, release_press_taken_by_os, resize_zones};

use crate::ModelApp;

/// Do what a caption button asked. Close quits the way the window's own close does: the app
/// asks about unsaved work before the window goes.
pub fn run(app: &mut ModelApp, ctx: &egui::Context, action: WindowAction) {
    if action == WindowAction::Close {
        app.quit_requested = true;
    }
    action.send(ctx);
}

#[cfg(test)]
mod tests {
    use egui::{Event, PointerButton, Pos2, RawInput, Rect, ResizeDirection, ViewportCommand, ViewportId, pos2, vec2};
    use fieldcraft_engine::Session;
    use solidcraft_ui_kit::window_chrome::{Caption, caption_id};

    use super::*;

    /// A headless window `width` pt wide showing the application bar and the resize zones.
    struct Win {
        ctx: egui::Context,
        app: ModelApp,
        width: f32,
        time: f64,
    }

    impl Win {
        fn new(width: f32) -> Self {
            let mut app = ModelApp::new(Session::default());
            app.custom_titlebar = true;
            Self { ctx: egui::Context::default(), app, width, time: 0.0 }
        }

        fn frame(&mut self, events: Vec<Event>) -> egui::FullOutput {
            self.time += 0.05;
            let input = RawInput {
                screen_rect: Some(Rect::from_min_size(Pos2::ZERO, vec2(self.width, 700.0))),
                time: Some(self.time),
                events,
                ..Default::default()
            };
            let app = &mut self.app;
            let mut out = self.ctx.run_ui(input, |ui| {
                crate::appbar::app_bar(app, ui);
                resize_zones(ui);
            });
            out.textures_delta.clear();
            out
        }

        fn gesture(&mut self, path: &[Pos2]) -> Vec<ViewportCommand> {
            let button = |p: Pos2, pressed| Event::PointerButton { pos: p, button: PointerButton::Primary, pressed, modifiers: Default::default() };
            let mut events = vec![vec![Event::PointerMoved(path[0])], vec![button(path[0], true)]];
            events.extend(path[1..].iter().map(|p| vec![Event::PointerMoved(*p)]));
            events.push(vec![button(*path.last().unwrap_or(&path[0]), false)]);
            events.into_iter().flat_map(|e| self.frame(e).viewport_output.remove(&ViewportId::ROOT).map(|v| v.commands).unwrap_or_default()).collect()
        }

        fn caption(&self, c: Caption) -> Rect {
            self.ctx.read_response(caption_id(c)).map(|r| r.rect).unwrap_or(Rect::NOTHING)
        }
    }

    #[test]
    fn caption_buttons_sit_in_the_corner_and_close_quits() {
        let mut w = Win::new(1200.0);
        w.frame(vec![]);
        w.frame(vec![]);
        let close = w.caption(Caption::Close);
        assert_eq!((close.right(), close.top()), (1200.0, 0.0));
        assert!(w.gesture(&[w.caption(Caption::Minimize).center()]).contains(&ViewportCommand::Minimized(true)));
        assert!(w.gesture(&[w.caption(Caption::Maximize).center()]).contains(&ViewportCommand::Maximized(true)));
        assert!(w.gesture(&[close.center()]).contains(&ViewportCommand::Close));
        assert!(w.app.quit_requested);
    }

    #[test]
    fn empty_bar_space_drags_and_edges_resize() {
        let mut w = Win::new(1600.0);
        w.frame(vec![]);
        // Between the quick-access button and the document tab.
        let empty = pos2(500.0, 17.0);
        let cmds = w.gesture(&[empty, empty + vec2(12.0, 4.0), empty + vec2(30.0, 8.0)]);
        assert_eq!(cmds.iter().filter(|c| matches!(c, ViewportCommand::StartDrag)).count(), 1, "{cmds:?}");
        let cmds = w.gesture(&[pos2(1.0, 300.0), pos2(8.0, 300.0)]);
        assert_eq!(cmds, vec![ViewportCommand::BeginResize(ResizeDirection::West)]);
    }
}
