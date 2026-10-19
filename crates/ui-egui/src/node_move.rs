//! The Move Nodes tool (`mesh.node.move`): with its form open, a drag on the mesh grabs the node
//! under the pointer (or the nodes the form lists, when the press is on one of them) and moves it
//! in the screen plane; the nodes round it follow as the form says (radius, Alt + scroll over
//! the view; falloff and strength, or harmonic). While dragging, the view draws the original
//! element edges as a ghost under the moved ones, those coloured by how each element's scaled
//! Jacobian changed (the error colour worse, the ok colour better), and the sphere of influence;
//! the form shows the region's numbers before and now. Release runs the command (one undo
//! step) and keeps the form open for the next drag; Esc during a drag cancels it.
//!
//! The tool only previews: what is drawn comes from [`fieldcraft_engine::cmd::NodeMovePrepared`],
//! the command's own move worked out without doing it. The part that does not depend on the
//! displacement (who follows, the design's faces) is worked out once per form and model
//! revision, so a drag frame only moves the region's nodes and measures its elements; a region
//! of more than [`PREVIEW_NODES`] nodes is not previewed (the move applies on release).

use std::sync::Arc;

use fieldcraft_engine::cmd::{NodeMoveCapped, NodeMovePrepared, NodeMoved};
use fieldcraft_render::{Scene, SceneLine, Vec3};
use serde_json::{Value, json};

use crate::theme::{self, Tokens};
use crate::{ModelApp, scene};

/// The command the tool runs.
pub const COMMAND: &str = "mesh.node.move";

/// Most nodes a live preview moves: a larger region is moved on release only. On the crank arm's
/// hex20 mesh (8264 elements) a region of 492 nodes prepares in 20 ms (falloff) to 27 ms
/// (harmonic) and draws a frame in 12 to 52 ms; one of 4051 nodes takes 0.8 s to prepare
/// harmonically and 0.3 s a frame, too slow to follow the pointer.
pub const PREVIEW_NODES: usize = 1_000;

/// How long the form must stay the same (seconds) before a changed move is prepared again: Alt +
/// scroll changes the radius every frame, and preparing (the design's faces, the harmonic
/// solve) takes longer than a frame on a big region.
pub const DEBOUNCE: f64 = 0.15;

/// The prepared part of the move.
#[derive(Clone)]
enum Prep {
    Ready(Arc<NodeMovePrepared>),
    /// More than [`PREVIEW_NODES`] would follow: not prepared, applied on release only.
    TooLarge(usize),
    Failed(String),
}

/// A drag of the Move Nodes tool in progress.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Drag {
    /// The point under the press, on the plane the drag moves in.
    pub start: Vec3,
    /// That plane's normal: the view direction when the drag started.
    pub normal: Vec3,
}

/// What the open form's move would do, and what it was worked out for.
#[derive(Default)]
pub struct Cache {
    /// The move up to its displacement, by model revision and form text without `by` and `to`.
    prepared: Option<((u64, String), Prep)>,
    /// The form (without the displacement) last seen and when (app time): a change waits
    /// [`DEBOUNCE`] before it is prepared again.
    seen: Option<((u64, String), f64)>,
    /// The move itself, by model revision and form text.
    preview: Option<((u64, String), Option<Arc<Preview>>)>,
}

/// The open form's move, without doing it.
pub struct Preview {
    pub sphere: Option<([f64; 3], f64)>,
    /// Nodes that take a share of the move.
    pub nodes: usize,
    /// The move; `None` when the region is too large to preview or the form is not complete.
    pub moved: Option<NodeMoved>,
    /// Why there is no move to show.
    pub note: Option<String>,
}

/// Where the ray `origin + t·dir` meets the plane through `point` with normal `normal`; `None`
/// when it runs along the plane.
pub fn on_plane(origin: Vec3, dir: Vec3, point: Vec3, normal: Vec3) -> Option<Vec3> {
    let denom = dir.dot(normal);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (point - origin).dot(normal) / denom;
    let p = origin + dir * t;
    p.is_finite().then_some(p)
}

/// The form text without the displacement: the key of the prepared move.
fn without_displacement(params: &Value) -> String {
    let mut p = params.clone();
    if let Some(o) = p.as_object_mut() {
        o.remove("by");
        o.remove("to");
    }
    p.to_string()
}

impl ModelApp {
    /// Whether the Move Nodes form is open.
    pub fn moving_nodes(&self) -> bool {
        self.dialog.as_ref().is_some_and(|d| d.command == COMMAND)
    }

    fn move_params(&self) -> Option<Value> {
        let d = self.dialog.as_ref().filter(|d| d.command == COMMAND)?;
        serde_json::from_str::<Value>(&d.text).ok().filter(Value::is_object)
    }

    fn set_move_params(&mut self, params: Value) {
        if let Some(d) = self.dialog.as_mut().filter(|d| d.command == COMMAND) {
            d.text = params.to_string();
            d.error = None;
            d.refused = None;
        }
    }

    /// The open form's move up to its displacement, cached by form and revision; `None` while a
    /// changed form waits out [`DEBOUNCE`] (the last one is kept meanwhile, never prepared on
    /// every frame of an Alt + scroll).
    fn move_prepared(&self, params: &Value) -> Option<Prep> {
        let key = (self.revision, without_displacement(params));
        let mut cache = self.node_move.borrow_mut();
        if let Some((k, p)) = &cache.prepared
            && *k == key
        {
            return Some(p.clone());
        }
        if cache.prepared.is_some() {
            match &cache.seen {
                Some((k, at)) if *k == key => {
                    if self.now - at < DEBOUNCE {
                        return None;
                    }
                }
                _ => {
                    cache.seen = Some((key, self.now));
                    return None;
                }
            }
        }
        let p = match fieldcraft_engine::cmd::node_move_prepare_capped(&self.session, params, Some(PREVIEW_NODES)) {
            Ok(NodeMoveCapped::Ready(p)) => Prep::Ready(Arc::from(p)),
            Ok(NodeMoveCapped::TooLarge(n)) => Prep::TooLarge(n),
            Err(e) => Prep::Failed(e.to_string()),
        };
        cache.prepared = Some((key, p.clone()));
        cache.seen = None;
        Some(p)
    }

    /// Whether a changed Move Nodes form waits to be prepared (the view draws again after
    /// [`DEBOUNCE`]).
    pub fn node_move_pending(&self) -> bool {
        self.node_move.borrow().seen.is_some()
    }

    /// The sphere the form describes, without preparing the move: its centre, or the listed
    /// nodes' centroid, and its radius.
    fn quick_sphere(&self, params: &Value) -> Option<([f64; 3], f64)> {
        let r = params.get("radius").and_then(Value::as_f64).filter(|r| *r > 0.0 && r.is_finite())?;
        let xyz = |v: &Value| -> Option<[f64; 3]> {
            let a = v.as_array()?;
            Some([a.first()?.as_f64()?, a.get(1)?.as_f64()?, a.get(2)?.as_f64()?])
        };
        if let Some(c) = params.get("center").and_then(xyz) {
            return Some((c, r));
        }
        let pts: Vec<[f64; 3]> = params
            .get("nodes")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|n| n.as_u64().and_then(|n| u32::try_from(n).ok()))
            .filter_map(|n| self.session.model.node(n).map(|p| p.xyz))
            .collect();
        if pts.is_empty() {
            return None;
        }
        let k = 1.0 / pts.len() as f64;
        Some((pts.iter().fold([0.0; 3], |a, p| [a[0] + p[0] * k, a[1] + p[1] * k, a[2] + p[2] * k]), r))
    }

    /// What the open Move Nodes form would do, worked out once per form text and model revision.
    pub fn node_move_preview(&self) -> Option<Arc<Preview>> {
        let d = self.dialog.as_ref().filter(|d| d.command == COMMAND)?;
        let key = (self.revision, d.text.clone());
        if let Some((k, p)) = &self.node_move.borrow().preview
            && *k == key
        {
            return p.clone();
        }
        let params = serde_json::from_str::<Value>(&d.text).ok()?;
        let Some(prep) = self.move_prepared(&params) else {
            // Waiting out the debounce: the sphere as the form has it, nothing cached.
            return Some(Arc::new(Preview { sphere: self.quick_sphere(&params), nodes: 0, moved: None, note: Some("Preview: updating…".into()) }));
        };
        let preview = match prep {
            Prep::Failed(e) => Some(Preview { sphere: self.quick_sphere(&params), nodes: 0, moved: None, note: Some(e) }),
            Prep::TooLarge(n) => Some(Preview {
                sphere: self.quick_sphere(&params),
                nodes: n,
                moved: None,
                note: Some(format!("{n} nodes follow: too many for a live preview (over {PREVIEW_NODES}); the move applies on release")),
            }),
            Prep::Ready(prep) => {
                let mut out = Preview { sphere: prep.sphere, nodes: prep.nodes(), moved: None, note: None };
                match prep.displacement(&params).and_then(|d| prep.preview_live(&self.session, d)) {
                    Ok(m) => out.moved = Some(m),
                    Err(e) => out.note = Some(e.to_string()),
                }
                Some(out)
            }
        }
        .map(Arc::new);
        self.node_move.borrow_mut().preview = Some((key, preview.clone()));
        preview
    }

    /// The open Move Nodes form's preview in a line: the region's quality before and now.
    pub fn node_move_note(&self) -> Option<String> {
        let p = self.node_move_preview()?;
        let Some(m) = &p.moved else { return p.note.clone() };
        let mut line = format!(
            "Preview: {} nodes move (largest {:.3}); scaled Jacobian min {:.3} to {:.3}, mean {:.3} to {:.3}; edge ratio max {:.2} to {:.2}; {} worse, {} better",
            m.moved,
            m.max_displacement,
            m.before.min_quality,
            m.after.min_quality,
            m.mean_before,
            m.mean_after,
            m.before.max_edge_ratio,
            m.after.max_edge_ratio,
            m.worse,
            m.better
        );
        if m.scale < 1.0 {
            line.push_str(&format!("; scaled to {:.0}% so that no element inverts", m.scale * 100.0));
        }
        Some(line)
    }

    /// Start a drag at window point `at` on the mesh: grab the listed nodes when the node under
    /// the press is one of them, else that node. Whether a drag started.
    pub fn node_move_grab(&mut self, at: egui::Pos2) -> bool {
        let Some(mut params) = self.move_params() else { return false };
        let (Some(rect), Some(hit)) = (self.viewport.rect, self.surface_point(at)) else { return false };
        let Some(node) = fieldcraft_engine::cmd::node_move_nearest(&self.session, [hit.x, hit.y, hit.z]) else { return false };
        let listed: Vec<u64> =
            params.get("nodes").and_then(Value::as_array).map(|a| a.iter().filter_map(Value::as_u64).collect()).unwrap_or_default();
        if !listed.contains(&u64::from(node)) {
            params["nodes"] = json!([node]);
        }
        if let Some(o) = params.as_object_mut() {
            for key in ["center", "faces", "to"] {
                o.remove(key);
            }
        }
        params["by"] = json!([0.0, 0.0, 0.0]);
        let (o, d) = self.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
        let normal = self.cam.back();
        let Some(start) = on_plane(o, d, hit, normal) else { return false };
        self.set_move_params(params);
        self.node_drag = Some(Drag { start, normal });
        true
    }

    /// The pointer of the drag is at window point `at`: the move follows it in the screen plane.
    pub fn node_move_drag(&mut self, at: egui::Pos2) {
        let (Some(drag), Some(rect), Some(mut params)) = (self.node_drag, self.viewport.rect, self.move_params()) else { return };
        let (o, d) = self.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
        let Some(p) = on_plane(o, d, drag.start, drag.normal) else { return };
        // The form is in the display units; the scene in metres.
        params["by"] = json!(self.shown_point(p - drag.start));
        self.set_move_params(params);
    }

    /// The drag ends: run the move (one undo step) and keep the form for the next one.
    pub fn node_move_release(&mut self) -> Result<(), String> {
        if self.node_drag.take().is_none() {
            return Ok(());
        }
        let Some(mut params) = self.move_params() else { return Ok(()) };
        let by = params.get("by").and_then(Value::as_array).map(|a| a.iter().filter_map(Value::as_f64).map(|x| x * x).sum::<f64>()).unwrap_or(0.0);
        if by > 0.0 {
            self.run(COMMAND, params.clone())?;
        }
        params["by"] = json!([0.0, 0.0, 0.0]);
        self.set_move_params(params);
        Ok(())
    }

    /// Esc during a drag: nothing moves.
    pub fn node_move_cancel(&mut self) {
        if self.node_drag.take().is_some()
            && let Some(mut params) = self.move_params()
        {
            params["by"] = json!([0.0, 0.0, 0.0]);
            self.set_move_params(params);
        }
    }

    /// Scale the open Move Nodes form's radius by `factor` (Alt + scroll over the view); a form
    /// without one starts at a twentieth of the mesh's size.
    pub fn node_move_scale_radius(&mut self, factor: f64) -> bool {
        let Some(mut params) = self.move_params() else { return false };
        let radius = match params.get("radius").and_then(Value::as_f64).filter(|r| *r > 0.0 && r.is_finite()) {
            Some(r) => r * factor,
            None => {
                let mut bounds = fieldcraft_render::Aabb3::EMPTY;
                for n in self.session.model.nodes().values() {
                    bounds.add(Vec3::from(n.xyz));
                }
                // A twentieth of the mesh's size (metres), in the form's display units.
                self.shown_length(bounds.diagonal() / 20.0)
            }
        };
        if !(radius > 0.0 && radius.is_finite()) {
            return false;
        }
        params["radius"] = json!(radius);
        self.set_move_params(params);
        true
    }

    /// The open Move Nodes form drawn over the scene: the sphere of influence, the original edges
    /// of the elements it changes as a ghost, and their moved edges coloured by how their scaled
    /// Jacobian changed.
    pub(crate) fn node_move_overlay(&self, s: &mut Scene, t: &Tokens) {
        let Some(p) = self.node_move_preview() else { return };
        if let Some((c, r)) = p.sphere {
            scene::sphere(s, Vec3::from(c), r, theme::rgb(t.pick_hover));
        }
        let Some(m) = &p.moved else { return };
        let model = &self.session.model;
        let ghost = theme::rgb(t.text_dim);
        let (worse, better, grabbed) = (theme::rgb(t.error), theme::rgb(t.ok), theme::rgb(t.pick_selected));
        let moved = |n: u32| m.positions.get(&n).copied().or_else(|| model.node(n).map(|x| x.xyz)).map(Vec3::from);
        let original = |n: u32| model.node(n).map(|x| Vec3::from(x.xyz));
        for c in m.elements.iter().filter(|c| c.worse() || c.better()) {
            let Some(e) = model.element(c.label) else { continue };
            let color = if c.worse() { worse } else { better };
            for edge in e.kind.edges() {
                let (Some(&a), Some(&b)) = (edge.first().and_then(|&i| e.nodes.get(i)), edge.get(1).and_then(|&i| e.nodes.get(i))) else { continue };
                if let (Some(pa), Some(pb)) = (original(a), original(b)) {
                    s.lines.push(SceneLine { a: pa, b: pb, color: ghost, width: 1.5 });
                }
                if let (Some(pa), Some(pb)) = (moved(a), moved(b)) {
                    s.lines.push(SceneLine { a: pa, b: pb, color, width: 3.0 });
                }
            }
        }
        // The grabbed nodes where they go, over everything.
        if let Some(params) = self.move_params() {
            for n in params.get("nodes").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_u64) {
                if let Some(x) = u32::try_from(n).ok().and_then(moved) {
                    s.on_top.push(SceneLine { a: x, b: x, color: grabbed, width: 7.0 });
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ray meets the drag plane where it should; one along the plane does not.
    #[test]
    fn the_drag_plane_meets_the_pointer_ray() {
        let p = on_plane(Vec3::new(1.0, 2.0, 10.0), Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, 0.0, 3.0), Vec3::Z).unwrap();
        assert!(p.dist(Vec3::new(1.0, 2.0, 3.0)) < 1e-12);
        let slanted = on_plane(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, -1.0), Vec3::new(5.0, 5.0, 0.0), Vec3::Z).unwrap();
        assert!(slanted.dist(Vec3::new(1.0, 0.0, 0.0)) < 1e-12);
        assert!(on_plane(Vec3::new(0.0, 0.0, 1.0), Vec3::X, Vec3::new(0.0, 0.0, 0.0), Vec3::Z).is_none());
    }

    #[test]
    fn the_prepared_key_leaves_out_the_displacement() {
        let a = json!({"nodes": [3], "by": [1.0, 0.0, 0.0], "radius": 2.0});
        let b = json!({"nodes": [3], "by": [0.0, 5.0, 0.0], "radius": 2.0});
        let c = json!({"nodes": [3], "to": [0.0, 5.0, 0.0], "radius": 3.0});
        assert_eq!(without_displacement(&a), without_displacement(&b));
        assert_ne!(without_displacement(&a), without_displacement(&c));
    }
}
