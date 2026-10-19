//! What the viewport draws: the model's display data in the chosen style and the theme's
//! colours. The GPU viewport and the headless render (`ui.render`) draw the same [`Scene`].
//!
//! Geometry (the open design, imported surfaces) draws as faces with their boundaries while
//! the model has no element faces; once it is meshed, as its boundaries drawn along the mesh
//! (the element edges on geometry edges), so the mesh is not hidden behind coincident faces and
//! every line lies on the element faces. The analysis setup's symbols ([`symbols`]) are a layer of
//! their own over the model, a constant size on screen, muted over a results plot.

use fieldcraft_render::glyphs::{Definition, Drawn, GlyphShow, GlyphView, Glyphs};
use fieldcraft_render::{Aabb3, Camera, DisplayMesh, Rgb, Scene, SceneFill, SceneLabel, SceneLine, ScenePoint, SceneTri, Vec3};

use crate::theme::{Tokens, rgb};

pub use fieldcraft_engine::view::{Show, Style};

/// Draw the analysis setup's symbols into `s` at the nodes' positions `pos` and symbol unit
/// `unit` (model length, [`fieldcraft_render::glyphs::symbol_unit`]), muted over a results plot,
/// the `highlight`ed definition's in the selection colour.
#[allow(clippy::too_many_arguments)]
pub fn symbols(
    s: &mut Scene,
    glyphs: &Glyphs,
    pos: impl Fn(u32) -> Option<Vec3>,
    show: Show,
    results: bool,
    unit: f64,
    highlight: Option<(Definition, &str)>,
    t: &Tokens,
) -> Drawn {
    let colors = t.glyph_colors();
    let colors = if results { colors.over_results(rgb(t.viewport_bottom)) } else { colors };
    let view = GlyphView { show: GlyphShow { loads: show.loads, interactions: show.interactions }, colors: &colors, unit, highlight };
    glyphs.draw(pos, &view, s)
}

/// Draw the datums the session's view shows at the model (`bounds`) into `s`, and `selected` (the
/// datum selected in the browser), shown or hidden, highlighted; returns the box around them.
pub fn datums(s: &mut Scene, session: &fieldcraft_engine::Session, bounds: &Aabb3, selected: Option<&str>, t: &Tokens) -> Aabb3 {
    let shown = |name: &str| session.shows(fieldcraft_engine::ItemKind::Datum, name);
    fieldcraft_render::glyphs::datums(&session.model, bounds, &shown, selected, &t.glyph_colors(), s)
}

/// A datum being defined in a form, drawn highlighted (through the part) where OK would put it.
pub fn datum_preview(s: &mut Scene, datum: fieldcraft_model::Datum, bounds: &Aabb3, t: &Tokens) -> Aabb3 {
    const NAME: &str = "preview";
    let mut model = fieldcraft_model::Model::default();
    if model.add_datum(NAME, datum).is_err() {
        return Aabb3::EMPTY;
    }
    fieldcraft_render::glyphs::datums(&model, bounds, &|_| true, Some(NAME), &t.glyph_colors(), s)
}

/// A sphere drawn translucent (a tool's region of influence, such as Relax Nodes'), with its
/// three great circles as lines.
pub fn sphere(s: &mut Scene, c: Vec3, r: f64, color: Rgb) {
    if !(r > 0.0 && r.is_finite()) {
        return;
    }
    const RINGS: usize = 12;
    const SEGMENTS: usize = 24;
    let at = |i: usize, j: usize| {
        let theta = std::f64::consts::PI * i as f64 / RINGS as f64;
        let phi = std::f64::consts::TAU * j as f64 / SEGMENTS as f64;
        Vec3 { x: theta.sin() * phi.cos(), y: theta.sin() * phi.sin(), z: theta.cos() }
    };
    for i in 0..RINGS {
        for j in 0..SEGMENTS {
            let (a, b, d, e) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            for tri in [[a, b, d], [a, d, e]] {
                let n = (tri[0] + tri[1] + tri[2]).normalized().unwrap_or(tri[0]);
                s.fills.push(SceneFill { tri: SceneTri { p: tri.map(|v| c + v * r), normal: n, color }, alpha: 0.18 });
            }
        }
    }
    for k in 0..SEGMENTS {
        let (p, q) = (std::f64::consts::TAU * k as f64 / SEGMENTS as f64, std::f64::consts::TAU * (k + 1) as f64 / SEGMENTS as f64);
        for (u, v) in [
            (Vec3 { x: p.cos(), y: p.sin(), z: 0.0 }, Vec3 { x: q.cos(), y: q.sin(), z: 0.0 }),
            (Vec3 { x: p.cos(), y: 0.0, z: p.sin() }, Vec3 { x: q.cos(), y: 0.0, z: q.sin() }),
            (Vec3 { x: 0.0, y: p.cos(), z: p.sin() }, Vec3 { x: 0.0, y: q.cos(), z: q.sin() }),
        ] {
            s.lines.push(SceneLine { a: c + u * r, b: c + v * r, color, width: 1.0 });
        }
    }
}

/// The scene for `display` drawn in `style` in the theme's colours: the engine's builder, which
/// the headless pictures use too.
pub fn build(display: &DisplayMesh, style: Style, show: Show, t: &Tokens) -> Scene {
    fieldcraft_engine::view::build(display, style, show, &t.palette())
}

/// Size of a seed's tick (points).
const SEED_TICK: f64 = 6.0;

/// Draw the model's edge seeds as ticks at the nodes they put along their edges (the drawn
/// bodies' B-rep edges, [`DisplayMesh::curves`], named by `names`): the seeded count, spaced by
/// the bias, so a grading shows as ticks crowding toward an end. Where the mesh drawn has
/// another number of element edges along the edge ([`meshed_divisions`]) the ticks are in the
/// error colour, with a label such as `6 seeded, 20 meshed` (#342).
pub fn seeds(s: &mut Scene, model: &fieldcraft_model::Model, names: &[crate::selection::BodyNames], display: &DisplayMesh, t: &Tokens) {
    if model.seeds().is_empty() {
        return;
    }
    let color = rgb(t.warm);
    for run in display.curves.chunk_by(|a, b| a.body == b.body && a.edge == b.edge) {
        let Some(first) = run.first() else { continue };
        let Some(edge) = names.get(first.body as usize).and_then(|b| b.edge(first.edge)) else { continue };
        let Some(seed) = model.seeds().iter().find(|seed| seed.edge.selects(&edge)) else { continue };
        let lengths: Vec<f64> = run.iter().map(|c| c.a.dist(c.b)).collect();
        let total: f64 = lengths.iter().sum();
        if !(total > 0.0 && total.is_finite()) {
            continue;
        }
        let count = seed.count(total);
        let divisions = meshed_divisions(display, run, total);
        let color = if divisions > 0 && divisions != count as usize {
            let mid = run.get(run.len() / 2).map_or(first.a, |c| (c.a + c.b) * 0.5);
            s.labels.push(SceneLabel { p: mid, text: format!("{count} seeded, {divisions} meshed"), color: rgb(t.error) });
            rgb(t.error)
        } else {
            color
        };
        // Walk the polyline to each node's arc length.
        let (mut k, mut walked) = (0, 0.0);
        for f in seed.spacing(count) {
            let target = f * total;
            while k + 1 < run.len() && walked + lengths.get(k).copied().unwrap_or(0.0) < target {
                walked += lengths.get(k).copied().unwrap_or(0.0);
                k += 1;
            }
            let (Some(c), Some(l)) = (run.get(k), lengths.get(k)) else { continue };
            let u = if *l > 0.0 { ((target - walked) / l).clamp(0.0, 1.0) } else { 0.0 };
            s.points.push(ScenePoint { p: c.a + (c.b - c.a) * u, color, size: SEED_TICK });
        }
    }
}

/// Draw the model's size controls (mesh.size.set) where they apply: the faces and bodies they
/// refine tinted on the drawn geometry, a box as its twelve edges, a sphere as
/// [`sphere`] draws one, each labelled with its name and size.
pub fn sizes(s: &mut Scene, model: &fieldcraft_model::Model, names: &[crate::selection::BodyNames], display: &DisplayMesh, t: &Tokens) {
    use fieldcraft_model::{Entity, GeometryRef, SizeRegion};
    if model.size_controls().is_empty() {
        return;
    }
    let color = rgb(t.warm);
    for (name, c) in model.size_controls() {
        let text = format!("{name} {}", c.size);
        let v = |p: [f64; 3]| Vec3 { x: p[0], y: p[1], z: p[2] };
        match &c.region {
            SizeRegion::Box { min, max } => {
                let corner = |k: usize| v([0, 1, 2].map(|a| if (k >> a) & 1 == 1 { max[a] } else { min[a] }));
                for k in 0..8usize {
                    for a in 0..3 {
                        if (k >> a) & 1 == 0 {
                            s.lines.push(SceneLine { a: corner(k), b: corner(k | (1 << a)), color, width: 1.5 });
                        }
                    }
                }
                s.labels.push(SceneLabel { p: (corner(0) + corner(7)) * 0.5, text, color });
            }
            SizeRegion::Sphere { center, radius } => {
                sphere(s, v(*center), *radius, color);
                s.labels.push(SceneLabel { p: v(*center), text, color });
            }
            SizeRegion::Geometry(refs) => {
                let (mut sum, mut n) = (Vec3 { x: 0.0, y: 0.0, z: 0.0 }, 0usize);
                for tri in &display.surface_tris {
                    let Some(body) = names.get(tri.body as usize) else { continue };
                    let face = body.faces.get(tri.face as usize).and_then(|f| GeometryRef::face(body.name.clone(), f.clone()).ok());
                    let hit = refs.iter().any(|r| match &r.entity {
                        Entity::Body => r.body == body.name,
                        _ => face.as_ref().is_some_and(|f| r.selects(f)),
                    });
                    if hit {
                        s.fills.push(SceneFill { tri: SceneTri { p: tri.p, normal: tri.normal, color }, alpha: 0.3 });
                        sum = sum + (tri.p[0] + tri.p[1] + tri.p[2]) * (1.0 / 3.0);
                        n += 1;
                    }
                }
                if n > 0 {
                    s.labels.push(SceneLabel { p: sum * (1.0 / n as f64), text, color });
                }
            }
        }
    }
}

/// How many element edges of the drawn mesh lie along the geometry edge `run` (its polyline,
/// `length` long): those whose ends and middle are all on it, counted once per element edge
/// (a quadratic edge is drawn in pieces). 0 when the mesh has none there (not meshed).
fn meshed_divisions(display: &DisplayMesh, run: &[fieldcraft_render::SurfaceCurve], length: f64) -> usize {
    // Within a small part of the edge: a curved edge's polyline is a chord of the true curve.
    let tol = 0.01 * length;
    let near = |p: Vec3| run.iter().any(|c| distance_to_segment(p, c.a, c.b) <= tol);
    let mut on: std::collections::BTreeSet<[u32; 2]> = std::collections::BTreeSet::new();
    for e in &display.edges {
        if near(e.p[0]) && near(e.p[1]) && near((e.p[0] + e.p[1]) * 0.5) {
            on.insert(e.corners);
        }
    }
    on.len()
}

fn distance_to_segment(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let ab = b - a;
    let len2 = ab.dot(ab);
    let k = if len2 > 0.0 { ((p - a).dot(ab) / len2).clamp(0.0, 1.0) } else { 0.0 };
    p.dist(a + ab * k)
}

/// The radius around the camera target that holds the bounds (for the clipping planes).
pub fn radius(b: &Aabb3, cam: &Camera) -> f64 {
    if b.is_empty() { 1.0 } else { (cam.target.dist(b.center()) + b.diagonal() * 0.5).max(fieldcraft_render::MIN_SIZE) }
}

#[cfg(test)]
mod tests {
    use fieldcraft_model::{ElementKind, Model};

    use super::*;

    #[test]
    fn seeds_draw_ticks_at_their_nodes() {
        use fieldcraft_model::{Bias, BiasToward, Divisions, Entity, GeometryRef, Seed, SeedConstraint};
        let mut d = DisplayMesh::default();
        // Two edges of 10 mm; the first in two segments.
        d.add_curves(&[vec![[0.0; 3], [4.0, 0.0, 0.0], [10.0, 0.0, 0.0]], vec![[0.0; 3], [0.0, 10.0, 0.0]]]);
        let names = [crate::selection::BodyNames { name: "B".into(), faces: vec![], edges: vec!["a|b".into(), "a|c".into()], vertices: Vec::new() }];
        let mut m = Model::default();
        let edge = GeometryRef::new("B", Entity::Edge("a|b".into())).unwrap();
        m.set_seed(Seed::new(edge, Divisions::Count(4), Some(Bias { ratio: 8.0, toward: BiasToward::Start }), SeedConstraint::Exact).unwrap())
            .unwrap();
        let mut s = Scene::default();
        seeds(&mut s, &m, &names, &d, &Tokens::dark());
        let xs: Vec<f64> = s.points.iter().map(|p| p.p.x).collect();
        let want = [0.0, 10.0 / 15.0, 30.0 / 15.0, 70.0 / 15.0, 10.0];
        assert_eq!(xs.len(), 5);
        assert!(xs.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-9), "{xs:?}");
        assert!(s.points.iter().all(|p| p.p.y == 0.0), "only the seeded edge has ticks");
        assert!(s.labels.is_empty());
        // A mesh with other divisions on the edge: the ticks turn to the error colour, labelled;
        // with as many, they stay as they are.
        let mesh = |n: u32| {
            let mut d = d.clone();
            d.edges = (0..n)
                .map(|k| {
                    let x = |k: u32| Vec3::new(10.0 * f64::from(k) / f64::from(n), 0.0, 0.0);
                    fieldcraft_render::DisplayEdge { p: [x(k), x(k + 1)], nodes: [k, k + 1], corners: [k, k + 1], feature: true }
                })
                .collect();
            // An element edge across the seeded edge is not along it.
            d.edges.push(fieldcraft_render::DisplayEdge {
                p: [Vec3::ZERO, Vec3 { x: 0.0, y: 10.0, z: 0.0 }],
                nodes: [0, 99],
                corners: [0, 99],
                feature: true,
            });
            d
        };
        let mut s = Scene::default();
        seeds(&mut s, &m, &names, &mesh(20), &Tokens::dark());
        let error = rgb(Tokens::dark().error);
        assert!(s.points.iter().all(|p| p.color == error));
        assert_eq!(s.labels.iter().map(|l| l.text.as_str()).collect::<Vec<_>>(), ["4 seeded, 20 meshed"]);
        let mut s = Scene::default();
        seeds(&mut s, &m, &names, &mesh(4), &Tokens::dark());
        assert!(s.labels.is_empty(), "the mesh follows the seed");
    }

    /// Over the mesh the geometry's outline is drawn along the element edges on geometry edges,
    /// never the design's tessellation, which lies off the element faces (#260).
    #[test]
    fn geometry_edges_over_the_mesh_follow_its_element_edges() {
        let mut m = Model::default();
        for (i, p) in
            [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]].into_iter().enumerate()
        {
            m.set_node(i as u32 + 1, p).unwrap();
        }
        m.set_element(1, ElementKind::Hex8, (1..=8).collect()).unwrap();
        let edge = fieldcraft_model::GeometryRef::new("B", fieldcraft_model::Entity::Edge("e".into())).unwrap();
        let tag = m.geometry_tag(&edge).unwrap();
        m.tag_element_edge(1, 1, tag).unwrap();
        let mut d = DisplayMesh::from_model(&m);
        // The design's own edge, a little off the element edge (its chord error).
        d.add_surface(&[[0.0, -0.01, 0.0], [1.0, -0.01, 0.0], [0.5, -0.01, 1.0]], &[[0, 1, 2]], &[0]);
        let t = Tokens::dark();
        let geometry = |s: &Scene| s.lines.iter().filter(|l| l.color == rgb(t.geometry_edge)).map(|l| [l.a, l.b]).collect::<Vec<_>>();
        let s = build(&d, Style::ShadedEdges, Show::DEFAULT, &t);
        assert_eq!(geometry(&s), vec![[Vec3 { x: 0.0, y: 0.0, z: 0.0 }, Vec3 { x: 1.0, y: 0.0, z: 0.0 }]]);
        // Without the mesh the design's edges are drawn.
        let s = build(&d, Style::ShadedEdges, Show { mesh: false, ..Show::DEFAULT }, &t);
        assert_eq!(geometry(&s).len(), 3);
    }

    #[test]
    fn styles_choose_faces_edges_and_nodes() {
        let mut m = Model::default();
        for (i, p) in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into_iter().enumerate() {
            m.set_node(i as u32 + 1, p).unwrap();
        }
        m.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 4]).unwrap();
        let d = DisplayMesh::from_model(&m);
        let t = Tokens::dark();
        let s = build(&d, Style::ShadedEdges, Show::DEFAULT, &t);
        assert_eq!((s.tris.len(), s.lines.len(), s.points.len()), (4, 6, 0));
        let s = build(&d, Style::Shaded, Show { nodes: true, ..Show::DEFAULT }, &t);
        assert_eq!((s.tris.len(), s.lines.len(), s.points.len()), (4, 0, 4));
        let s = build(&d, Style::ShadedEdges, Show { mesh: false, ..Show::DEFAULT }, &t);
        assert_eq!((s.tris.len(), s.lines.len()), (0, 0), "a hidden mesh draws nothing");
        let s = build(&d, Style::Wireframe, Show::DEFAULT, &t);
        assert_eq!((s.tris.len(), s.lines.len()), (0, 6));
        assert_eq!(s.tris.first().map(|t| t.color), None);
        assert!(radius(&d.bounds, &Camera::default()) > 0.5);
    }

    #[test]
    fn size_controls_draw_where_they_apply() {
        use fieldcraft_model::{SizeControl, SizeRegion};
        let mut m = Model::default();
        m.set_size_control("Zone", SizeControl { region: SizeRegion::Box { min: [0.0; 3], max: [1.0; 3] }, size: 0.5 }).unwrap();
        m.set_size_control("Ball", SizeControl { region: SizeRegion::Sphere { center: [0.0; 3], radius: 2.0 }, size: 0.5 }).unwrap();
        let mut s = Scene::default();
        sizes(&mut s, &m, &[], &DisplayMesh::default(), &Tokens::dark());
        assert!(s.lines.len() >= 12 && !s.fills.is_empty(), "the box's edges and the sphere");
        let mut labels: Vec<&str> = s.labels.iter().map(|l| l.text.as_str()).collect();
        labels.sort();
        assert_eq!(labels, ["Ball 0.5", "Zone 0.5"]);
        let mut empty = Scene::default();
        sizes(&mut empty, &Model::default(), &[], &DisplayMesh::default(), &Tokens::dark());
        assert!(empty.lines.is_empty() && empty.labels.is_empty());
    }
}
