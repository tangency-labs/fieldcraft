//! Select the stored frame, array, and optional displacement of a VTK result.
use fieldcraft_engine::external_view::{Association, Component, visible};
use serde_json::json;
use solidcraft_ui_kit as kit;

use crate::{ModelApp, results::number};

pub(crate) fn controls(app: &mut ModelApp, ui: &mut egui::Ui) {
    let Some((frame, view)) = visible(&app.session) else { return };
    let frames = app.session.external_results().map_or(0, |s| s.frames.len());
    let mut change = None;
    let mut playback = None;
    if let Some(note) = &app.external_display_error {
        kit::status_line(ui, note, kit::StatusTone::Error);
    }
    kit::scope(ui, |ui| {
        if frames > 1 {
            kit::row(ui, "Frame", |ui| {
                let mut k = view.frame;
                if kit::slider(
                    ui,
                    &mut k,
                    0..=frames.saturating_sub(1),
                    kit::Slider { integer: true, format: |v| format!("{v:.0}"), ..Default::default() },
                )
                .changed()
                {
                    change = Some(json!({"frame":k}));
                }
            });
            kit::row(ui, "Playback", |ui| {
                if kit::button(ui, if view.animate { "Pause" } else { "Play" }, kit::Button::Normal).clicked() {
                    playback = Some(!view.animate);
                }
            });
            kit::row(ui, "Time", |ui| kit::text(ui, &number(frame.time), kit::Text::Cell));
        }
        let selected = view.field.as_deref().map_or_else(|| "None".into(), |name| format!("{name} ({})", view.association.label()));
        kit::row(ui, "Field", |ui| {
            kit::bounded_dropdown_menu(ui, "fc_external_field", &selected, ui.available_width(), |ui| {
                if ui.selectable_label(view.field.is_none(), "None").clicked() {
                    change = Some(json!({"field":""}));
                }
                for a in [Association::Point, Association::Cell] {
                    for (name, f) in a.fields(frame) {
                        if let Some(&component) = Component::for_components(f.components).first()
                            && ui
                                .selectable_label(view.association == a && view.field.as_ref() == Some(name), format!("{name} ({})", a.label()))
                                .clicked()
                        {
                            change = Some(json!({"association":a,"field":name,"component":component}));
                        }
                    }
                }
            })
        });
        if let Some(f) = view.field.as_ref().and_then(|n| view.association.fields(frame).get(n)).filter(|f| f.components > 1) {
            kit::row(ui, "Component", |ui| {
                kit::bounded_dropdown_menu(ui, "fc_external_component", view.component.label(), ui.available_width(), |ui| {
                    for &c in Component::for_components(f.components) {
                        if ui.selectable_label(view.component == c, c.label()).clicked() {
                            change = Some(json!({"component":c}));
                        }
                    }
                })
            });
        }
        if let Some(range) = app.session.external_bands() {
            kit::row(ui, "Range", |ui| kit::text(ui, &format!("{} to {}", number(range.min), number(range.max)), kit::Text::Cell));
            kit::row(ui, "", |ui| {
                let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), 8.0), egui::Sense::hover());
                let palette = crate::theme::Tokens::get(app.ui.dark).contour_stops();
                let bands = fieldcraft_engine::results::DEFAULT_BANDS;
                for k in 0..bands {
                    let c = fieldcraft_render::contour::band_color(&palette, k, bands);
                    let x = |i: u32| rect.left() + rect.width() * i as f32 / bands as f32;
                    let swatch = egui::Rect::from_min_max(egui::pos2(x(k), rect.top()), egui::pos2(x(k + 1), rect.bottom()));
                    ui.painter().rect_filled(swatch, 0.0, egui::Color32::from_rgb(c.0, c.1, c.2));
                }
            });
        }
        kit::row(ui, "Displacement", |ui| {
            kit::bounded_dropdown_menu(ui, "fc_external_displacement", view.displacement.as_deref().unwrap_or("None"), ui.available_width(), |ui| {
                if ui.selectable_label(view.displacement.is_none(), "None").clicked() {
                    change = Some(json!({"displacement":""}));
                }
                for (name, _) in frame.point_data.iter().filter(|(_, f)| matches!(f.components, 2 | 3)) {
                    if ui.selectable_label(view.displacement.as_ref() == Some(name), name).clicked() {
                        change = Some(json!({"displacement":name}));
                    }
                }
            })
        });
        if view.displacement.is_some() {
            kit::row(ui, "Scale", |ui| {
                let mut scale = view.scale;
                if kit::number_field(ui, &mut scale, 0.0..=fieldcraft_engine::results::MAX_DEFORMATION_SCALE, 0.1, "×").changed() {
                    change = Some(json!({"scale":scale}));
                }
            });
        }
    });
    if let Some(on) = playback {
        let _ = app.run("results.animate", json!({"on":on}));
    }
    if let Some(p) = change {
        let _ = app.run("results.external.view", p);
    }
}

impl ModelApp {
    /// Advance the read-only display using the PVD's physical sample times. The same active
    /// frame is reported to scripts and saved with the project. Seeking manually pauses it.
    pub(crate) fn advance_external_playback(&mut self, ctx: &egui::Context) {
        let Some((_, v)) = visible(&self.session).filter(|(_, v)| v.animate) else {
            self.external_clock = None;
            return;
        };
        let generation = self.session.external_generation();
        let clock = self.external_clock.get_or_insert((generation, self.now, v.frame));
        if clock.0 != generation {
            *clock = (generation, self.now, v.frame);
        }
        let Some(series) = self.session.external_results() else { return };
        let Some(frame) =
            fieldcraft_engine::external_view::ExternalView::playback_frame(series, clock.2, (self.now - clock.1) / crate::results::PERIOD)
        else {
            return;
        };
        let is_last = frame + 1 == series.frames.len();
        if let Some(v) = self.session.external_view.as_mut() {
            v.frame = frame;
            if is_last {
                v.animate = false;
                self.external_clock = None;
            }
        }
        if !is_last {
            ctx.request_repaint_after(std::time::Duration::from_secs_f64(crate::results::PERIOD / crate::results::FRAMES as f64));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_engine::external_view::{Cell, Field, Frame, Grid, Series};
    use std::collections::BTreeMap;

    fn app_with_frames() -> ModelApp {
        let frames = [0., 1., 4.].map(|time| Frame {
            displacement_unit: None,
            time,
            grid: Grid { points: vec![[time, 0., 0.]], cells: vec![Cell { vtk_type: 1, points: vec![0] }] },
            point_data: BTreeMap::from([("T".into(), Field { components: 1, values: vec![time] })]),
            cell_data: BTreeMap::new(),
        });
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.session.set_external_results(Series { frames: Vec::from(frames) }).unwrap();
        app
    }

    #[test]
    fn scripts_restart_playback_only_after_a_successful_playback_command() {
        let mut app = app_with_frames();
        let ctx = egui::Context::default();
        let seek_and_play = json!([{"command":"results.external.view","params":{"frame":0}},{"command":"results.animate","params":{"on":true}}]);
        let mut with_header = json!(app.session.journal_header);
        with_header["commands"] = seek_and_play.clone();
        for script in [seek_and_play.clone(), with_header] {
            app.run("results.animate", json!({"on":true})).unwrap();
            app.advance_external_playback(&ctx);
            app.now += crate::results::PERIOD * 0.5;
            app.advance_external_playback(&ctx);
            assert_eq!(app.session.external_view.as_ref().unwrap().frame, 1);
            let clock = app.external_clock;
            app.run_script(&json!([{"command":"model.info"}])).unwrap();
            assert_eq!(app.external_clock, clock);
            assert!(app.run_script(&json!([{"command":"results.animate","params":{"on":"bad"}}])).is_err());
            assert_eq!(app.external_clock, clock);
            app.run_script(&script).unwrap();
            app.advance_external_playback(&ctx);
            assert_eq!(app.session.external_view.as_ref().unwrap().frame, 0);
        }
        let mut partial = seek_and_play;
        partial.as_array_mut().unwrap().push(json!({"command":"no.such.command"}));
        assert!(app.run_script(&partial).is_err());
        assert!(app.external_clock.is_none());
    }

    #[test]
    fn play_reports_the_visible_frame_and_manual_seeking_pauses_it() {
        let mut app = app_with_frames();
        let ctx = egui::Context::default();
        app.run("results.animate", json!({"on":true})).unwrap();
        app.now = 50.0;
        app.advance_external_playback(&ctx);
        app.now += crate::results::PERIOD * 0.25;
        app.advance_external_playback(&ctx);
        assert_eq!(app.run("results.external.view", json!({})).unwrap()["time"], 1.0);
        assert_eq!(app.scene().1.points[0].p.x, 1.0);
        app.run("results.external.view", json!({"frame":0})).unwrap();
        app.now += crate::results::PERIOD;
        app.advance_external_playback(&ctx);
        assert_eq!(app.scene().1.points[0].p.x, 0.0);
        assert!(!app.session.external_view.as_ref().unwrap().animate);
        app.run("results.clear", json!({})).unwrap();
        assert!(app.scene().1.points.is_empty());
    }
}
