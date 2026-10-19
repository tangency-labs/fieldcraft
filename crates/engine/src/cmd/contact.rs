//! `interaction.contact.*`: contact definitions (pairs of surfaces, general contact) and the
//! search for surfaces that touch; `interaction.tie*`: ties between non-matching meshes;
//! `interaction.rigid_wall*`: rigid walls, which belong to an explicit dynamics step (its
//! `rigid_walls`, also edited in the step's form).
//!
//! A contact's surfaces are sets of element faces, or of element edges in a 2D model
//! ([`fieldcraft_model::contact`]). The commands
//! take the friction, penalty scale and initial-penetration choice of the model's form; `edit`
//! replaces a definition with the model's JSON form, `delete` removes it.
//!
//! `find` looks for the bodies of the mesh that touch: the solid elements are split into
//! connected parts (elements sharing a node belong to the same part), and a node of one part's
//! exterior faces within `tolerance` of an exterior face of another part makes the two parts
//! touch there. Each touching pair is reported with the faces of each side that are close to the
//! other; with `create`, each becomes two face sets and a contact pair.

use crate::schema::ProblemKind;
use std::collections::BTreeMap;

use fieldcraft_fem::contact::{Grid, bounds, closest_point};
use fieldcraft_model::{Contact, ContactSurfaces, ElementFace, Friction, InitialPenetration, Model, Procedure, RigidWall, Set, Step, Tie};
use serde_json::{Value, json};

use super::{name, opt};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// The options every contact takes, from the parameters.
fn contact(p: &Value, surfaces: ContactSurfaces) -> Result<Contact> {
    Ok(Contact {
        surfaces,
        friction: opt::<Friction>(p, "friction")?.unwrap_or(Friction::NONE),
        penalty_scale: opt::<f64>(p, "penalty_scale")?.unwrap_or(1.0),
        initial_penetration: opt::<InitialPenetration>(p, "initial_penetration")?.unwrap_or_default(),
    })
}

/// A required name-valued parameter.
fn text<'a>(p: &'a Value, key: &str) -> Result<&'a str> {
    p.get(key).and_then(Value::as_str).ok_or_else(|| EngineError::param(key, ProblemKind::NotAName))
}

pub const RIGID_WALL: CommandSpec = CommandSpec {
    id: "interaction.rigid_wall",
    label: "Rigid Wall",
    params: r#"{"step": "Impact" (an explicit dynamics step; default: the model's only one), "set": "front" (the nodes the wall stops), "point": [x, y, z] (on the wall), "normal": [nx, ny, nz] (toward the side the nodes stay on), "friction": {"static": 0.3, "kinetic": 0.2, "decay": 10} (Coulomb; default frictionless)} — a fixed rigid plane the set's nodes cannot pass, added to the step's rigid walls (the same list as model.step.create's `rigid_walls`). Replies the step and the wall's index in it"#,
    form: Some(super::forms::RIGID_WALL),
    kind: CommandKind::Model,
    requires: &[Requirement::Step, Requirement::Set],
    run: rigid_wall,
};

pub const RIGID_WALL_DELETE: CommandSpec = CommandSpec {
    id: "interaction.rigid_wall.delete",
    label: "Delete Rigid Wall",
    params: r#"{"step": "Impact", "index": 0} — removes a rigid wall from an explicit dynamics step (its index in the step's `rigid_walls`, from 0)"#,
    form: Some(super::forms::RIGID_WALL_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: rigid_wall_delete,
};

/// The explicit dynamics step `step` names, or the model's only one when it names none.
fn explicit_step(m: &Model, p: &Value) -> Result<Step> {
    let explicit = |s: &&Step| matches!(s.procedure, Procedure::ExplicitDynamic(_));
    let step = match p.get("step") {
        None | Some(Value::Null) => {
            let mut all = m.steps().iter().filter(explicit);
            match (all.next(), all.next()) {
                (Some(s), None) => s,
                (None, _) => {
                    return Err(EngineError::Params(
                        "the model has no explicit dynamics step, where rigid walls act: create one (Steps › Step, Dynamic, explicit)".into(),
                    ));
                }
                (Some(_), Some(_)) => {
                    return Err(EngineError::param("step", ProblemKind::Invalid("name the explicit dynamics step: the model has several".into())));
                }
            }
        }
        Some(_) => {
            let n = text(p, "step")?;
            let s = m.steps().iter().find(|s| s.name == n).ok_or_else(|| EngineError::Params(format!("no step `{n}`")))?;
            if !explicit(&s) {
                return Err(EngineError::Params(format!("step `{n}` is not an explicit dynamics step: rigid walls act only there")));
            }
            s
        }
    };
    Ok(step.clone())
}

fn rigid_wall(s: &mut Session, p: &Value) -> Result<Value> {
    let mut step = explicit_step(&s.model, p)?;
    let set = text(p, "set")?.to_string();
    if s.model.set(&set).is_none() {
        return Err(EngineError::Params(format!("no set `{set}`")));
    }
    let wall = RigidWall {
        set,
        point: super::req(p, "point")?,
        normal: super::req(p, "normal")?,
        friction: opt::<Friction>(p, "friction")?.unwrap_or(Friction::NONE),
    };
    let Procedure::ExplicitDynamic(e) = &mut step.procedure else { return Err(EngineError::Params("not an explicit dynamics step".into())) };
    e.rigid_walls.push(wall);
    let index = e.rigid_walls.len() - 1;
    let name = step.name.clone();
    s.model.replace_step(step)?;
    Ok(json!({"step": name, "index": index}))
}

fn rigid_wall_delete(s: &mut Session, p: &Value) -> Result<Value> {
    let mut step = explicit_step(&s.model, p)?;
    let index: usize = super::req(p, "index")?;
    let Procedure::ExplicitDynamic(e) = &mut step.procedure else { return Err(EngineError::Params("not an explicit dynamics step".into())) };
    if index >= e.rigid_walls.len() {
        return Err(EngineError::param("index", ProblemKind::Invalid(format!("step `{}` has {} rigid walls", step.name, e.rigid_walls.len()))));
    }
    e.rigid_walls.remove(index);
    s.model.replace_step(step)?;
    Ok(Value::Null)
}

pub const PAIR: CommandSpec = CommandSpec {
    id: "interaction.contact.pair",
    label: "Contact Pair",
    params: concat!(
        r#"{"name": "impact", "main": "plate" (a set of element faces), "secondary": "ball" (another), "#,
        r#""friction": {"static": 0.3, "kinetic": 0.2, "decay": 10} (Coulomb, μ = μk + (μs − μk) e^(−decay × slip rate); default frictionless), "penalty_scale": 1 (0 < scale ≤ 1: a factor on the penalty: explicit, the stiffest that leaves the time increment alone; nonlinear static, the augmented Lagrangian's), "initial_penetration": "adjust" (nodes that start inside: explicit moves them onto the surface, strain-free; nonlinear static takes the overlap out of their gap; both report them) | "refuse" (refuse the analysis, naming the first)"#,
        r#"} — contact between two surfaces: in explicit dynamics each one's nodes are kept out of the other's faces; in nonlinear static steps `secondary` is the mortar (non-mortar) side whose pressure field is enforced against `main` by an augmented Lagrangian, with Coulomb friction by a return map (static coefficient while a node sticks, kinetic once it slips) (fails if the name is taken)"#
    ),
    form: Some(super::forms::CONTACT_PAIR),
    kind: CommandKind::Model,
    requires: &[Requirement::Set],
    run: |s: &mut Session, p: &Value| {
        let n = name(p)?;
        let surfaces = ContactSurfaces::Pair { main: text(p, "main")?.to_string(), secondary: text(p, "secondary")?.to_string() };
        s.model.add_contact(n, contact(p, surfaces)?)?;
        Ok(json!({"name": n}))
    },
};

pub const GENERAL: CommandSpec = CommandSpec {
    id: "interaction.contact.general",
    label: "General Contact",
    params: concat!(
        r#"{"name": "all", "surface": "skin" (optional: a set of element faces, each in contact with every other: self contact; default every exterior face of the solids), "#,
        r#""friction": {"static": 0.3, "kinetic": 0.2, "decay": 10} (Coulomb, μ = μk + (μs − μk) e^(−decay × slip rate); default frictionless), "penalty_scale": 1 (0 < scale ≤ 1: a factor on the penalty: explicit, the stiffest that leaves the time increment alone; nonlinear static, the augmented Lagrangian's), "initial_penetration": "adjust" (nodes that start inside: explicit moves them onto the surface, strain-free; nonlinear static takes the overlap out of their gap; both report them) | "refuse" (refuse the analysis, naming the first)"#,
        r#"} — single-surface contact of everything that may touch (explicit dynamics, and mortar self-contact in nonlinear static steps; fails if the name is taken)"#
    ),
    form: Some(super::forms::CONTACT_GENERAL),
    kind: CommandKind::Model,
    requires: &[],
    run: |s: &mut Session, p: &Value| {
        let n = name(p)?;
        let surface = match p.get("surface") {
            None | Some(Value::Null) => None,
            Some(_) => Some(text(p, "surface")?.to_string()),
        };
        s.model.add_contact(n, contact(p, ContactSurfaces::General { surface })?)?;
        Ok(json!({"name": n}))
    },
};

pub const EDIT: CommandSpec = CommandSpec {
    id: "interaction.contact.edit",
    label: "Edit Contact",
    params: r#"{"name": "impact", "surfaces": {"type": "pair", "main": "plate", "secondary": "ball"} | {"type": "general", "surface": "skin" (optional)}, "friction": {…}, "penalty_scale": 1, "initial_penetration": "adjust|refuse"} — replaces the contact of that name (the model's form; see interaction.contact.pair)"#,
    form: Some(super::forms::CONTACT_EDIT),
    kind: CommandKind::Model,
    requires: &[Requirement::Contact],
    run: |s: &mut Session, p: &Value| {
        let n = name(p)?;
        let mut obj = p.as_object().cloned().ok_or_else(|| EngineError::Params("the parameters of a contact are a JSON object".into()))?;
        obj.remove("name");
        let c: Contact = serde_json::from_value(Value::Object(obj)).map_err(|e| EngineError::Params(format!("contact: {e}")))?;
        s.model.replace_contact(n, c)?;
        Ok(json!({"name": n}))
    },
};

pub const DELETE: CommandSpec = CommandSpec {
    id: "interaction.contact.delete",
    label: "Delete Contact",
    params: r#"{"name": "impact"}"#,
    form: Some(super::forms::CONTACT_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Contact],
    run: |s: &mut Session, p: &Value| {
        s.model.remove_contact(name(p)?)?;
        Ok(Value::Null)
    },
};

/// A tie from the parameters.
fn tie(p: &Value) -> Result<Tie> {
    Ok(Tie {
        main: text(p, "main")?.to_string(),
        secondary: text(p, "secondary")?.to_string(),
        position_tolerance: opt::<f64>(p, "position_tolerance")?,
    })
}

pub const TIE: CommandSpec = CommandSpec {
    id: "interaction.tie",
    label: "Tie",
    params: r#"{"name": "glue", "main": "plate" (a set of element faces), "secondary": "rib" (another), "position_tolerance": 0.1 (optional: how far a secondary node may lie from the main surface and still be tied; default half the main element's thickness)} — bonds the secondary surface to the main one across non-matching meshes: each secondary node follows the main surface under it (dual mortar constraints in the implicit analyses; a constant stress crosses the tie exactly); fails if the name is taken"#,
    form: Some(super::forms::TIE),
    kind: CommandKind::Model,
    requires: &[Requirement::Set],
    run: |s: &mut Session, p: &Value| {
        let n = name(p)?;
        s.model.add_tie(n, tie(p)?)?;
        Ok(json!({"name": n}))
    },
};

pub const TIE_EDIT: CommandSpec = CommandSpec {
    id: "interaction.tie.edit",
    label: "Edit Tie",
    params: r#"{"name": "glue", "main": "plate", "secondary": "rib", "position_tolerance": 0.1} — replaces the tie of that name (see interaction.tie)"#,
    form: Some(super::forms::TIE_EDIT),
    kind: CommandKind::Model,
    requires: &[Requirement::Tie],
    run: |s: &mut Session, p: &Value| {
        let n = name(p)?;
        s.model.replace_tie(n, tie(p)?)?;
        Ok(json!({"name": n}))
    },
};

pub const TIE_DELETE: CommandSpec = CommandSpec {
    id: "interaction.tie.delete",
    label: "Delete Tie",
    params: r#"{"name": "glue"}"#,
    form: Some(super::forms::TIE_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Tie],
    run: |s: &mut Session, p: &Value| {
        s.model.remove_tie(name(p)?)?;
        Ok(Value::Null)
    },
};

pub const FIND: CommandSpec = CommandSpec {
    id: "interaction.contact.find",
    label: "Find Contacts",
    params: concat!(
        r#"{"tolerance": 0.1 (distance within which surfaces touch; default a tenth of the mean exterior face size), "create": false (true: add two face sets and a contact pair per touching pair of parts, named find-N-main, find-N-secondary and find-N), "#,
        r#""friction", "penalty_scale", "initial_penetration" (of the created pairs, as interaction.contact.pair)} — the parts of the mesh (solid elements connected through their nodes) whose exterior faces come within the tolerance of each other"#
    ),
    form: Some(super::forms::CONTACT_FIND),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: find,
};

/// Most touching part pairs `find` reports.
pub const MAX_FOUND: usize = fieldcraft_model::contact::MAX_CONTACTS;

/// Most point–face checks `find` makes.
pub const MAX_CHECKS: usize = 200_000_000;

/// The connected parts of the solid elements: per element label, its part (numbered from 0 in
/// the order of their lowest element label).
fn parts(model: &Model) -> BTreeMap<u32, usize> {
    // Union–find over the node labels of the solid elements.
    let solids: Vec<(u32, &[u32])> =
        model.elements().iter().filter(|(_, e)| e.kind.dimension() == 3).map(|(l, e)| (*l, e.nodes.as_slice())).collect();
    let mut index: BTreeMap<u32, usize> = BTreeMap::new();
    for (_, nodes) in &solids {
        for n in *nodes {
            let next = index.len();
            index.entry(*n).or_insert(next);
        }
    }
    let mut parent: Vec<usize> = (0..index.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while let Some(&p) = parent.get(i) {
            if p == i {
                break;
            }
            let grand = parent.get(p).copied().unwrap_or(p);
            if let Some(slot) = parent.get_mut(i) {
                *slot = grand;
            }
            i = grand;
        }
        i
    }
    for (_, nodes) in &solids {
        let mut it = nodes.iter().filter_map(|n| index.get(n).copied());
        let Some(first) = it.next() else { continue };
        for other in it {
            let (a, b) = (root(&mut parent, first), root(&mut parent, other));
            if a != b {
                let (lo, hi) = (a.min(b), a.max(b));
                if let Some(slot) = parent.get_mut(hi) {
                    *slot = lo;
                }
            }
        }
    }
    // Number the parts in element order.
    let mut numbers: BTreeMap<usize, usize> = BTreeMap::new();
    let mut out = BTreeMap::new();
    for (label, nodes) in &solids {
        let Some(r) = nodes.first().and_then(|n| index.get(n).copied()).map(|i| root(&mut parent, i)) else { continue };
        let next = numbers.len();
        out.insert(*label, *numbers.entry(r).or_insert(next));
    }
    out
}

fn find(s: &mut Session, p: &Value) -> Result<Value> {
    let model = &s.model;
    let faces = model.exterior_faces();
    let part_of = parts(model);
    // Each exterior face's corner coordinates and part.
    let mut face_xyz: Vec<Vec<[f64; 3]>> = Vec::with_capacity(faces.len());
    for f in &faces {
        let labels = model.face_nodes(*f).unwrap_or_default();
        let corners = if labels.len() == 3 || labels.len() == 6 { 3 } else { 4 };
        face_xyz.push(labels.iter().take(corners).filter_map(|l| model.node(*l)).map(|n| n.xyz).collect());
    }
    let size = {
        let sizes: Vec<f64> = face_xyz
            .iter()
            .filter_map(|x| {
                let b = bounds(x, 0.0);
                let d = (0..3).map(|k| b[1][k] - b[0][k]).fold(0.0, f64::max);
                d.is_finite().then_some(d)
            })
            .collect();
        if sizes.is_empty() { 0.0 } else { sizes.iter().sum::<f64>() / sizes.len() as f64 }
    };
    let tolerance = opt::<f64>(p, "tolerance")?.unwrap_or(0.1 * size);
    if !(tolerance >= 0.0 && tolerance.is_finite()) {
        return Err(EngineError::Params(format!("`tolerance` must be a finite distance ≥ 0, got {tolerance}")));
    }
    let create = opt::<bool>(p, "create")?.unwrap_or(false);
    let options = contact(p, ContactSurfaces::General { surface: None })?;
    options.validate()?;

    let boxes: Vec<_> = face_xyz.iter().map(|x| bounds(x, tolerance)).collect();
    let grid = Grid::new(&boxes, 4 * faces.len() + 64);
    let part = |c: usize| faces.get(c).and_then(|f| part_of.get(&f.element).copied());
    // Each face's centroid and outward normal there.
    let centroid = |x: &[[f64; 3]]| -> [f64; 3] { [0, 1, 2].map(|k| x.iter().map(|p| p[k]).sum::<f64>() / x.len().max(1) as f64) };
    let normals: Vec<Option<[f64; 3]>> = face_xyz.iter().map(|x| closest_point(x, centroid(x)).map(|q| q.normal)).collect();
    // (lower part, higher part) → the faces of each side that touch the other: a face whose
    // centroid or a corner is within the tolerance of a face of the other part facing it
    // (normals more than 60° apart from parallel).
    let mut touching: BTreeMap<(usize, usize), [std::collections::BTreeSet<usize>; 2]> = BTreeMap::new();
    let mut checks = 0usize;
    for (c, xyz) in face_xyz.iter().enumerate() {
        let (Some(k), Some(normal)) = (part(c), normals.get(c).copied().flatten()) else { continue };
        for point in xyz.iter().copied().chain(std::iter::once(centroid(xyz))) {
            let candidates = grid.candidates(point);
            checks += candidates.len();
            if checks > MAX_CHECKS {
                return Err(EngineError::Params(format!(
                    "the tolerance {tolerance:e} puts too many faces near each other (more than {MAX_CHECKS} checks); lower it"
                )));
            }
            for &d in candidates {
                let d = d as usize;
                let Some(other) = part(d).filter(|o| *o != k) else { continue };
                let facing = normals.get(d).copied().flatten().is_some_and(|n| (0..3).map(|i| n[i] * normal[i]).sum::<f64>() < -0.5);
                let near = facing && face_xyz.get(d).and_then(|x| closest_point(x, point)).is_some_and(|q| q.distance <= tolerance);
                if near {
                    let key = (k.min(other), k.max(other));
                    touching.entry(key).or_default()[usize::from(k != key.0)].insert(c);
                    if touching.len() > MAX_FOUND {
                        return Err(EngineError::Params(format!("more than {MAX_FOUND} pairs of parts touch; lower the tolerance")));
                    }
                }
            }
        }
    }
    let face_list = |set: &std::collections::BTreeSet<usize>| -> Vec<ElementFace> { set.iter().filter_map(|c| faces.get(*c).copied()).collect() };
    let mut found = Vec::new();
    let mut created = Vec::new();
    for (n, ((a, b), sides)) in touching.iter().enumerate() {
        let (main, secondary) = (face_list(&sides[0]), face_list(&sides[1]));
        let mut entry = json!({"parts": [a, b], "faces": [main.len(), secondary.len()]});
        if create {
            let names = [format!("find-{}-main", n + 1), format!("find-{}-secondary", n + 1), format!("find-{}", n + 1)];
            created.push((names, main, secondary));
            entry["contact"] = json!(format!("find-{}", n + 1));
        }
        found.push(entry);
    }
    let part_count = part_of.values().max().map_or(0, |m| m + 1);
    for ([main_name, secondary_name, contact_name], main, secondary) in created {
        s.model.add_set(&main_name, Set { faces: main, ..Set::default() })?;
        s.model.add_set(&secondary_name, Set { faces: secondary, ..Set::default() })?;
        let surfaces = ContactSurfaces::Pair { main: main_name, secondary: secondary_name };
        s.model.add_contact(&contact_name, Contact { surfaces, ..options.clone() })?;
    }
    Ok(json!({"parts": part_count, "exterior_faces": faces.len(), "tolerance": tolerance, "touching": found}))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rigid walls go into the explicit step (the model's only one, or the one named) and come
    /// out by index; a model without an explicit step, or a step of another kind, is refused.
    #[test]
    fn rigid_walls_belong_to_the_explicit_step() {
        let mut s = Session::default();
        super::super::define::tests::defined_bar(&mut s);
        let wall = json!({"set": "end", "point": [3, 0, 0], "normal": [-1, 0, 0], "friction": {"static": 0.2}});
        assert!(s.execute("interaction.rigid_wall", &wall).is_err(), "no explicit step yet");
        let mut on_bend = wall.clone();
        on_bend["step"] = json!("Bend");
        assert!(s.execute("interaction.rigid_wall", &on_bend).is_err(), "a static step");
        s.execute("model.step.create", &json!({"name": "Hit", "procedure": {"type": "explicit_dynamic", "end_time": 1e-4}})).unwrap();
        let r = s.execute("interaction.rigid_wall", &wall).unwrap();
        assert_eq!(r, json!({"step": "Hit", "index": 0}));
        let walls = |s: &Session| match &s.model.step("Hit").unwrap().procedure {
            Procedure::ExplicitDynamic(e) => e.rigid_walls.clone(),
            _ => panic!(),
        };
        assert_eq!(walls(&s).len(), 1);
        assert_eq!(walls(&s)[0].set, "end");
        for bad in [json!({"set": "nope", "point": [0, 0, 0], "normal": [1, 0, 0]}), json!({"set": "end", "point": [0, 0, 0], "normal": [0, 0, 0]})] {
            assert!(s.execute("interaction.rigid_wall", &bad).is_err(), "{bad}");
        }
        assert!(s.execute("interaction.rigid_wall.delete", &json!({"step": "Hit", "index": 1})).is_err());
        s.execute("interaction.rigid_wall.delete", &json!({"step": "Hit", "index": 0})).unwrap();
        assert!(walls(&s).is_empty());
        s.execute("model.step.create", &json!({"name": "Hit2", "procedure": {"type": "explicit_dynamic", "end_time": 1e-4}})).unwrap();
        assert!(s.execute("interaction.rigid_wall", &wall).is_err(), "two explicit steps: name one");
    }

    /// Two unit cubes of 2×2×2 hex8 each, the second resting on the first (separate nodes), and
    /// a third cube far away.
    fn stacked(s: &mut Session) {
        let mut next_node = 1u32;
        let mut next_element = 1u32;
        for (z0, x0) in [(0.0, 0.0), (1.0, 0.0), (0.0, 10.0)] {
            let base = next_node;
            for k in 0..3 {
                for j in 0..3 {
                    for i in 0..3 {
                        let xyz = [x0 + 0.5 * i as f64, 0.5 * j as f64, z0 + 0.5 * k as f64];
                        s.execute("mesh.node.add", &json!({"label": next_node, "xyz": xyz})).unwrap();
                        next_node += 1;
                    }
                }
            }
            let n = |i: u32, j: u32, k: u32| base + i + 3 * j + 9 * k;
            for k in 0..2 {
                for j in 0..2 {
                    for i in 0..2 {
                        let nodes = [
                            n(i, j, k),
                            n(i + 1, j, k),
                            n(i + 1, j + 1, k),
                            n(i, j + 1, k),
                            n(i, j, k + 1),
                            n(i + 1, j, k + 1),
                            n(i + 1, j + 1, k + 1),
                            n(i, j + 1, k + 1),
                        ];
                        s.execute("mesh.element.add", &json!({"label": next_element, "kind": "hex8", "nodes": nodes})).unwrap();
                        next_element += 1;
                    }
                }
            }
        }
    }

    #[test]
    fn contacts_are_created_edited_and_deleted() {
        let mut s = Session::default();
        stacked(&mut s);
        s.execute("model.set.create", &json!({"name": "a", "faces": [{"element": 5, "face": 2}]})).unwrap();
        s.execute("model.set.create", &json!({"name": "b", "faces": [{"element": 9, "face": 1}]})).unwrap();
        s.execute("interaction.contact.pair", &json!({"name": "ab", "main": "a", "secondary": "b", "friction": {"static": 0.2}})).unwrap();
        let c = s.model.contact("ab").unwrap();
        assert_eq!(c.friction, Friction::coulomb(0.2));
        assert!(s.execute("interaction.contact.pair", &json!({"name": "ab", "main": "a", "secondary": "b"})).is_err(), "taken");
        assert!(s.execute("interaction.contact.pair", &json!({"name": "aa", "main": "a", "secondary": "a"})).is_err(), "self pair");
        assert!(s.execute("interaction.contact.pair", &json!({"name": "x", "main": "a", "secondary": "nope"})).is_err());
        assert!(s.execute("interaction.contact.pair", &json!({"name": "x", "main": "a", "secondary": "b", "penalty_scale": 2})).is_err());
        assert!(s.execute("interaction.contact.pair", &json!({"name": "x", "main": "a", "secondary": "b", "friction": {"static": -1}})).is_err());
        s.execute("interaction.contact.general", &json!({"name": "all"})).unwrap();
        s.execute("interaction.contact.general", &json!({"name": "self", "surface": "a", "initial_penetration": "refuse"})).unwrap();
        assert_eq!(s.model.contact("self").unwrap().initial_penetration, InitialPenetration::Refuse);
        s.execute("interaction.contact.edit", &json!({"name": "all", "surfaces": {"type": "general", "surface": "b"}, "penalty_scale": 0.5}))
            .unwrap();
        assert_eq!(s.model.contact("all").unwrap().surfaces, ContactSurfaces::General { surface: Some("b".into()) });
        assert!(s.execute("interaction.contact.edit", &json!({"name": "nope", "surfaces": {"type": "general"}})).is_err());
        assert!(s.execute("model.set.delete", &json!({"name": "a"})).unwrap_err().to_string().contains("contact `ab`"));
        s.execute("model.set.create", &json!({"name": "n", "nodes": [1, 2]})).unwrap();
        s.execute("interaction.contact.general", &json!({"name": "bad", "surface": "n"})).unwrap();
        let check = s.execute("model.check", &json!({})).unwrap();
        assert!(check.to_string().contains("set `n` has no element faces"), "{check}");
        for n in ["ab", "all", "self", "bad"] {
            s.execute("interaction.contact.delete", &json!({"name": n})).unwrap();
        }
        assert!(s.model.contacts().is_empty());
    }

    #[test]
    fn ties_are_created_edited_and_deleted() {
        let mut s = Session::default();
        stacked(&mut s);
        s.execute("model.set.create", &json!({"name": "a", "faces": [{"element": 5, "face": 2}]})).unwrap();
        s.execute("model.set.create", &json!({"name": "b", "faces": [{"element": 9, "face": 1}]})).unwrap();
        s.execute("interaction.tie", &json!({"name": "glue", "main": "a", "secondary": "b"})).unwrap();
        assert_eq!(s.model.tie("glue"), Some(&Tie::new("a", "b")));
        for bad in [
            json!({"name": "glue", "main": "a", "secondary": "b"}),
            json!({"name": "x", "main": "a", "secondary": "a"}),
            json!({"name": "x", "main": "a", "secondary": "nope"}),
            json!({"name": "x", "main": "a", "secondary": "b", "position_tolerance": -1}),
            json!({"name": "x", "main": "a", "secondary": "b", "position_tolerance": "far"}),
            json!({"name": "x", "main": 3, "secondary": "b"}),
        ] {
            assert!(s.execute("interaction.tie", &bad).is_err(), "{bad}");
        }
        s.execute("interaction.tie.edit", &json!({"name": "glue", "main": "b", "secondary": "a", "position_tolerance": 0.01})).unwrap();
        assert_eq!(s.model.tie("glue").and_then(|t| t.position_tolerance), Some(1e-5), "0.01 mm, stored in metres");
        assert!(s.execute("interaction.tie.edit", &json!({"name": "nope", "main": "a", "secondary": "b"})).is_err());
        assert!(s.execute("model.set.delete", &json!({"name": "a"})).unwrap_err().to_string().contains("tie `glue`"));
        s.execute("interaction.tie.delete", &json!({"name": "glue"})).unwrap();
        assert!(s.model.ties().is_empty() && s.execute("interaction.tie.delete", &json!({"name": "glue"})).is_err());
    }

    #[test]
    fn find_reports_and_creates_the_touching_parts() {
        let mut s = Session::default();
        stacked(&mut s);
        let r = s.execute("interaction.contact.find", &json!({})).unwrap();
        assert_eq!(r["parts"], 3, "{r}");
        assert_eq!(r["exterior_faces"], 3 * 24, "{r}");
        // The two stacked cubes touch over four faces each; the far one touches nothing.
        assert_eq!(r["touching"], json!([{"parts": [0, 1], "faces": [4, 4]}]), "{r}");
        assert!(s.model.contacts().is_empty(), "nothing created without `create`");
        let r = s.execute("interaction.contact.find", &json!({"create": true, "friction": {"static": 0.1}})).unwrap();
        assert_eq!(r["touching"][0]["contact"], "find-1");
        let c = s.model.contact("find-1").unwrap();
        assert_eq!(c.friction, Friction::coulomb(0.1));
        let main = s.model.resolve_set("find-1-main").unwrap();
        assert!(main.faces.iter().all(|f| (5..=8).contains(&f.element) && f.face == 2) && main.faces.len() == 4, "{main:?}");
        let secondary = s.model.resolve_set("find-1-secondary").unwrap();
        assert!(secondary.faces.iter().all(|f| (9..=12).contains(&f.element) && f.face == 1) && secondary.faces.len() == 4, "{secondary:?}");
        // A tight tolerance on a gap finds nothing; a second create reuses no names.
        assert!(s.execute("interaction.contact.find", &json!({"create": true})).is_err(), "the names are taken");
        for bad in [json!({"tolerance": -1}), json!({"tolerance": "x"}), json!({"tolerance": f64::MAX, "create": "yes"})] {
            assert!(s.execute("interaction.contact.find", &bad).is_err(), "{bad}");
        }
        let empty = s.execute("interaction.contact.find", &json!({"tolerance": 0})).unwrap();
        assert_eq!(empty["touching"].as_array().map(Vec::len), Some(1), "touching faces are at distance 0: {empty}");
    }
}
