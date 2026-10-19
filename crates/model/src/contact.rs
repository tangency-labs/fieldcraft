//! Contact: which surfaces may touch, and the friction between them.
//!
//! A contact is a model definition, active in every step that can take it (explicit dynamics,
//! and nonlinear static steps; the linear analyses refuse contact by name). Its
//! surfaces are sets with element
//! faces (in a 2D model, element edges), so they survive remeshing like a pressure load does:
//! - a **pair**: the faces of `main` against those of `secondary`, both ways (each surface's
//!   nodes are kept out of the other's faces, so neither surface is privileged);
//! - **general** contact: every face against every other face of one surface (single-surface or
//!   self contact, for folding and buckling), or of every exterior face of the model's solid
//!   elements (in a 2D model every exterior edge of its 2D elements) when no surface is named.
//!
//! **Friction** is Coulomb's law with a coefficient that falls from its static to its kinetic
//! value with the slip rate γ̇: μ = μ_k + (μ_s − μ_k) e^{−d γ̇} (the exponential decay of
//! E. Rabinowicz's observations, as LS-DYNA's and Abaqus's exponential decay friction write it).
//! A static step has no slip rate: its nodes take the static coefficient while they stick and the
//! kinetic one once they slip.
//! The same friction applies to the planar rigid walls of an explicit step
//! ([`crate::RigidWall`]).
//!
//! **Penalty.** Explicit contact pushes penetrating nodes back with a penalty spring whose
//! stiffness the solver chooses from the masses and the stable time increment, the stiffest
//! that leaves the increment alone; `penalty_scale` (at most 1) softens it. Nonlinear static
//! steps enforce a pair's mortar pressure by an augmented Lagrangian whose penalty the solver
//! chooses from the elements under the secondary surface; `penalty_scale` scales it too.
//!
//! **Initial penetration.** Nodes of a surface that start inside the other surface are either
//! moved onto it before the analysis, strain-free (`adjust`, the default, each one reported;
//! nonlinear static steps take the overlap out of their gap instead), or
//! the analysis is refused naming the first node, its depth and the face it is behind
//! (`refuse`).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, Model, Result, context, insert_new, remove_existing, replace_existing, require};

/// Most contact definitions in a model.
pub const MAX_CONTACTS: usize = 1000;

/// Largest friction coefficient (a coefficient of 10 already locks every practical interface).
pub const MAX_FRICTION: f64 = 10.0;

/// Largest exponential decay coefficient, in s/m (time per length: 1/d is the slip rate at which
/// the coefficient has gone 63% of the way from static to kinetic): 10¹⁵ s/m, a slip rate of
/// a femtometre per second, is a limit no interface reaches.
pub const MAX_DECAY: f64 = 1e15;

/// Penalty scale factors outside this range are refused (above 1 the contact springs would
/// shorten the stable time increment; see the module docs).
pub const PENALTY_SCALE_RANGE: (f64, f64) = (1e-6, 1.0); // unit-factors: a factor on the penalty stiffness

/// Coulomb friction with an exponential decay from the static to the kinetic coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Friction {
    /// The static coefficient μ_s (0: frictionless).
    #[serde(rename = "static", default)]
    pub static_coefficient: f64,
    /// The kinetic coefficient μ_k (by default the static one: no decay).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kinetic: Option<f64>,
    /// The decay coefficient d of μ = μ_k + (μ_s − μ_k) e^{−d γ̇}, per slip rate (time per length).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub decay: f64,
}

fn is_zero(v: &f64) -> bool {
    *v == 0.0
}

impl Friction {
    /// No friction.
    pub const NONE: Friction = Friction { static_coefficient: 0.0, kinetic: None, decay: 0.0 };

    /// A constant coefficient.
    pub fn coulomb(mu: f64) -> Friction {
        Friction { static_coefficient: mu, ..Friction::NONE }
    }

    pub fn is_frictionless(&self) -> bool {
        self.static_coefficient == 0.0 && self.kinetic.is_none_or(|k| k == 0.0)
    }

    /// The kinetic coefficient (the static one when none is given).
    pub fn kinetic_coefficient(&self) -> f64 {
        self.kinetic.unwrap_or(self.static_coefficient)
    }

    /// The coefficient at slip rate `rate` (≥ 0).
    pub fn coefficient(&self, rate: f64) -> f64 {
        let (s, k) = (self.static_coefficient, self.kinetic_coefficient());
        if self.decay == 0.0 { s } else { k + (s - k) * (-self.decay * rate.abs()).exp() }
    }

    pub fn validate(&self) -> Result<()> {
        for (what, v) in [("static friction coefficient", self.static_coefficient), ("kinetic friction coefficient", self.kinetic_coefficient())] {
            if !(0.0..=MAX_FRICTION).contains(&v) {
                return Err(Error::Invalid(format!("the {what} is in 0..={MAX_FRICTION}, got {v}")));
            }
        }
        if !(0.0..=MAX_DECAY).contains(&self.decay) {
            return Err(Error::Invalid(format!("the friction decay coefficient is in 0..={MAX_DECAY:e}, got {}", self.decay)));
        }
        if self.decay == 0.0 && self.kinetic_coefficient() != self.static_coefficient {
            return Err(Error::Invalid("a kinetic friction coefficient other than the static one needs a decay coefficient".into()));
        }
        Ok(())
    }
}

/// Which faces a contact puts against which.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum ContactSurfaces {
    /// The faces of two sets, each surface's nodes against the other's faces.
    Pair { main: String, secondary: String },
    /// Every face against every other face of a set's faces (single-surface or self contact), or
    /// of every exterior face of the model's solid elements when no set is named.
    General {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        surface: Option<String>,
    },
}

impl ContactSurfaces {
    /// The sets it names.
    pub fn sets(&self) -> Vec<&str> {
        match self {
            ContactSurfaces::Pair { main, secondary } => vec![main.as_str(), secondary.as_str()],
            ContactSurfaces::General { surface } => surface.iter().map(String::as_str).collect(),
        }
    }
}

/// What happens to nodes that start inside the surface they are in contact with.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InitialPenetration {
    /// Move them onto the surface before the analysis (strain-free), and report them.
    #[default]
    Adjust,
    /// Refuse the analysis, naming the first one.
    Refuse,
}

/// A contact definition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Contact {
    pub surfaces: ContactSurfaces,
    #[serde(default, skip_serializing_if = "Friction::is_frictionless")]
    pub friction: Friction,
    /// A factor on the solver's penalty stiffness.
    #[serde(default = "one")]
    pub penalty_scale: f64,
    #[serde(default)]
    pub initial_penetration: InitialPenetration,
}

fn one() -> f64 {
    1.0
}

impl Contact {
    /// A frictionless contact with the default penalty and initial adjustment.
    pub fn new(surfaces: ContactSurfaces) -> Contact {
        Contact { surfaces, friction: Friction::NONE, penalty_scale: 1.0, initial_penetration: InitialPenetration::Adjust }
    }

    pub fn validate(&self) -> Result<()> {
        if let ContactSurfaces::Pair { main, secondary } = &self.surfaces
            && main == secondary
        {
            return Err(Error::Invalid(format!(
                "a contact pair needs two different surfaces, got `{main}` twice (self contact is general contact on that surface)"
            )));
        }
        self.friction.validate()?;
        let (lo, hi) = PENALTY_SCALE_RANGE;
        if !(lo..=hi).contains(&self.penalty_scale) {
            return Err(Error::Invalid(format!("the penalty scale is in {lo:e}..={hi}, got {}", self.penalty_scale)));
        }
        Ok(())
    }
}

impl Model {
    pub fn contacts(&self) -> &BTreeMap<String, Contact> {
        &self.contacts
    }

    pub fn contact(&self, name: &str) -> Option<&Contact> {
        self.contacts.get(name)
    }

    fn check_contact(&self, name: &str, c: &Contact) -> Result<()> {
        c.validate().map_err(|e| context("contact", name, e))?;
        for set in c.surfaces.sets() {
            require(&self.sets, "set", set).map_err(|e| context("contact", name, e))?;
        }
        Ok(())
    }

    pub fn add_contact(&mut self, name: &str, contact: Contact) -> Result<()> {
        self.check_contact(name, &contact)?;
        if self.contacts.len() >= MAX_CONTACTS {
            return Err(Error::Limit(format!("a model holds at most {MAX_CONTACTS} contacts")));
        }
        insert_new(&mut self.contacts, "contact", name, contact)
    }

    pub fn replace_contact(&mut self, name: &str, contact: Contact) -> Result<()> {
        self.check_contact(name, &contact)?;
        replace_existing(&mut self.contacts, "contact", name, contact)
    }

    pub fn remove_contact(&mut self, name: &str) -> Result<()> {
        remove_existing(&mut self.contacts, "contact", name).map(|_| ())
    }

    /// The exterior edges of the 2D elements (the edges no other 2D element shares), in
    /// element order: the boundary of a 2D model, which general contact takes by default.
    pub fn exterior_edges(&self) -> Vec<crate::ElementEdge> {
        self.unshared(2, |kind| kind.edges(), 2, |element, i| crate::ElementEdge { element, edge: i })
    }

    /// The exterior faces of the solid elements (the faces no other solid element shares), in
    /// element order.
    pub fn exterior_faces(&self) -> Vec<crate::ElementFace> {
        // A face is known by its sorted corner nodes (three or four).
        self.unshared(3, |kind| kind.faces(), 4, |element, face| crate::ElementFace { element, face })
    }

    /// The mesh's outer skin, in element order: the [`exterior_faces`](Model::exterior_faces)
    /// of the solids and face 1 (the side its normal points to) of every 2D element. Each face's
    /// corners run counter-clockwise seen from outside ([`crate::ElementKind::faces`]).
    pub fn skin_faces(&self) -> Vec<crate::ElementFace> {
        let shells = self.elements.iter().filter(|(_, e)| e.kind.dimension() == 2).map(|(&element, _)| crate::ElementFace { element, face: 1 });
        let mut out: Vec<crate::ElementFace> = self.exterior_faces().into_iter().chain(shells).collect();
        out.sort_unstable();
        out
    }

    /// The sides (`sides` of each element kind) of the elements of `dimension` that no other
    /// such element shares, each known by its sorted corners (at most `corners`; triangles of
    /// quadratic kinds by their first three), made by `make` from the element and the 1-based
    /// side number, sorted.
    fn unshared<T: Ord>(
        &self,
        dimension: u8,
        sides: impl Fn(crate::ElementKind) -> &'static [&'static [usize]],
        corners: usize,
        make: impl Fn(u32, u8) -> T,
    ) -> Vec<T> {
        let mut keyed: Vec<([u32; 4], u32, u8)> = Vec::new();
        for (label, e) in self.elements.iter().filter(|(_, e)| e.kind.dimension() == dimension) {
            for (i, local) in sides(e.kind).iter().enumerate() {
                let count = if local.len() == 6 || local.len() == 3 { 3 } else { 4 }.min(corners);
                let mut key = [u32::MAX; 4];
                for (k, n) in key.iter_mut().zip(local.iter().take(count)) {
                    *k = e.nodes.get(*n).copied().unwrap_or(u32::MAX);
                }
                key.sort_unstable();
                keyed.push((key, *label, i as u8 + 1));
            }
        }
        keyed.sort_unstable();
        let mut out: Vec<T> = keyed
            .chunk_by(|a, b| a.0 == b.0)
            .filter_map(|run| match run {
                [(_, element, side)] => Some(make(*element, *side)),
                _ => None,
            })
            .collect();
        out.sort_unstable();
        out
    }

    /// The contacts that name a set, as ``contact `name` `` strings.
    pub(crate) fn contact_users<'a>(&'a self, set: &'a str) -> impl Iterator<Item = String> + 'a {
        self.contacts.iter().filter(move |(_, c)| c.surfaces.sets().contains(&set)).map(|(n, _)| format!("contact `{n}`"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Set;

    fn model() -> Model {
        let mut m = Model::default();
        for s in ["top", "bottom"] {
            m.add_set(s, Set { nodes: vec![1], ..Set::default() }).unwrap();
        }
        m
    }

    #[test]
    fn contacts_validate_refer_to_sets_and_round_trip() {
        let mut m = model();
        let json = serde_json::json!({
            "surfaces": {"type": "pair", "main": "top", "secondary": "bottom"},
            "friction": {"static": 0.3, "kinetic": 0.2, "decay": 5.0}, "penalty_scale": 0.5
        });
        let c: Contact = serde_json::from_value(json).unwrap();
        assert_eq!((c.friction.static_coefficient, c.friction.kinetic_coefficient(), c.initial_penetration), (0.3, 0.2, InitialPenetration::Adjust));
        assert_eq!(serde_json::from_value::<Contact>(serde_json::to_value(&c).unwrap()).unwrap(), c);
        m.add_contact("impact", c.clone()).unwrap();
        assert!(matches!(m.remove_set("top"), Err(Error::InUse(_))), "the contact refers to the set");
        assert!(m.add_contact("impact", c.clone()).is_err(), "names stay unique");
        let general: Contact = serde_json::from_value(serde_json::json!({"surfaces": {"type": "general"}})).unwrap();
        assert_eq!(general.surfaces, ContactSurfaces::General { surface: None });
        assert!(general.friction.is_frictionless() && general.penalty_scale == 1.0);
        m.add_contact("everything", general).unwrap();
        m.remove_contact("impact").unwrap();
        m.remove_set("top").unwrap();
        assert!(m.remove_contact("impact").is_err());
        assert!(serde_json::from_value::<Contact>(serde_json::json!({"surfaces": {"type": "general"}, "mu": 1})).is_err(), "typos are errors");

        let base = Contact::new(ContactSurfaces::Pair { main: "bottom".into(), secondary: "top".into() });
        let bad = [
            Contact { friction: Friction::coulomb(-0.1), ..base.clone() },
            Contact { friction: Friction::coulomb(f64::NAN), ..base.clone() },
            Contact { friction: Friction { kinetic: Some(f64::INFINITY), ..Friction::coulomb(0.1) }, ..base.clone() },
            Contact { friction: Friction { decay: -1.0, ..Friction::coulomb(0.1) }, ..base.clone() },
            Contact { penalty_scale: 0.0, ..base.clone() },
            Contact { penalty_scale: f64::NAN, ..base.clone() },
            Contact { penalty_scale: 1.5, ..base.clone() },
            Contact::new(ContactSurfaces::Pair { main: "bottom".into(), secondary: "bottom".into() }),
        ];
        for b in bad {
            assert!(b.validate().is_err(), "{b:?}");
        }
        let missing = Contact::new(ContactSurfaces::General { surface: Some("nope".into()) });
        assert!(matches!(m.add_contact("x", missing), Err(Error::NotFound(_))));
    }

    #[test]
    fn friction_decays_from_static_to_kinetic() {
        let f = Friction { static_coefficient: 0.5, kinetic: Some(0.2), decay: 2.0 };
        assert_eq!(f.coefficient(0.0), 0.5);
        assert!((f.coefficient(1.0) - (0.2 + 0.3 * (-2.0f64).exp())).abs() < 1e-15);
        assert!((f.coefficient(1e6) - 0.2).abs() < 1e-15);
        assert_eq!(Friction::coulomb(0.3).coefficient(10.0), 0.3);
        // A kinetic coefficient acts through the decay.
        assert!(Friction { kinetic: Some(0.1), ..Friction::coulomb(0.3) }.validate().is_err());
        assert!(Friction { kinetic: Some(0.3), ..Friction::coulomb(0.3) }.validate().is_ok());
        assert!(Friction::NONE.is_frictionless() && !Friction::coulomb(0.1).is_frictionless());
    }

    #[test]
    fn the_check_wants_solid_faces_on_contact_surfaces() {
        let mut m = Model::default();
        let corners = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
        for (i, x) in corners.iter().enumerate() {
            m.set_node(i as u32 + 1, *x).unwrap();
        }
        m.set_element(1, crate::ElementKind::Hex8, (1..=8).collect()).unwrap();
        m.add_set("nodes", Set { nodes: vec![1, 2], ..Set::default() }).unwrap();
        m.add_set("top", Set { faces: vec![crate::ElementFace { element: 1, face: 2 }], ..Set::default() }).unwrap();
        m.add_contact("c", Contact::new(ContactSurfaces::Pair { main: "top".into(), secondary: "nodes".into() })).unwrap();
        let errors: Vec<String> = m.check().iter().filter(|d| d.severity == crate::Severity::Error).map(ToString::to_string).collect();
        assert!(errors.iter().any(|e| e.contains("contact `c`") && e.contains("set `nodes` has no element faces")), "{errors:?}");
        assert!(!errors.iter().any(|e| e.contains("set `top`")), "{errors:?}");
    }

    #[test]
    fn exterior_faces_are_those_of_one_solid() {
        // Two hexahedra sharing a face, and a tetrahedron on the top of the second.
        let mut m = Model::default();
        for k in 0..3 {
            for (i, xy) in [[0., 0.], [1., 0.], [1., 1.], [0., 1.]].iter().enumerate() {
                m.set_node(4 * k + i as u32 + 1, [xy[0], xy[1], k as f64]).unwrap();
            }
        }
        m.set_node(13, [0.5, 0.5, 3.0]).unwrap();
        m.set_element(1, crate::ElementKind::Hex8, (1..=8).collect()).unwrap();
        m.set_element(2, crate::ElementKind::Hex8, (5..=12).collect()).unwrap();
        m.set_element(3, crate::ElementKind::Tet4, vec![9, 10, 11, 13]).unwrap();
        let faces = m.exterior_faces();
        let of = |e: u32| faces.iter().filter(|f| f.element == e).map(|f| f.face).collect::<Vec<_>>();
        assert_eq!(of(1), [1, 3, 4, 5, 6], "all but the shared top");
        assert_eq!(of(2), [2, 3, 4, 5, 6], "all but the shared bottom (its top meets the tet over half its area)");
        assert_eq!(of(3), [1, 2, 3, 4], "a triangle never matches a quad");
        assert!(Model::default().exterior_faces().is_empty());
    }

    #[test]
    fn skin_faces_count_solid_and_shell_sides() {
        let cube = |m: &mut Model, base: u32, z: f64| {
            for (i, xy) in [[0., 0.], [1., 0.], [1., 1.], [0., 1.]].iter().enumerate() {
                m.set_node(base + i as u32, [xy[0], xy[1], z]).unwrap();
            }
        };
        let mut m = Model::default();
        cube(&mut m, 1, 0.0);
        cube(&mut m, 5, 1.0);
        m.set_element(1, crate::ElementKind::Hex8, (1..=8).collect()).unwrap();
        assert_eq!(m.skin_faces().len(), 6, "one hex8");
        cube(&mut m, 9, 2.0);
        m.set_element(2, crate::ElementKind::Hex8, (5..=12).collect()).unwrap();
        assert_eq!(m.skin_faces().len(), 10, "two hex8 sharing a face");

        let mut t = Model::default();
        for (i, p) in [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [0.5, 0., 0.],
            [0.5, 0.5, 0.],
            [0., 0.5, 0.],
            [0., 0., 0.5],
            [0.5, 0., 0.5],
            [0., 0.5, 0.5],
        ]
        .iter()
        .enumerate()
        {
            t.set_node(i as u32 + 1, *p).unwrap();
        }
        t.set_element(1, crate::ElementKind::Tet10, (1..=10).collect()).unwrap();
        assert_eq!(t.skin_faces().len(), 4, "a tet10");

        let mut s = Model::default();
        cube(&mut s, 1, 0.0);
        s.set_element(7, crate::ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
        assert_eq!(s.skin_faces(), [crate::ElementFace { element: 7, face: 1 }], "a quad4 shell, its top side");
    }

    #[test]
    fn exterior_edges_are_those_of_one_2d_element() {
        // Two quads sharing an edge and a triangle on the second's top edge.
        let mut m = Model::default();
        for (i, xy) in [[0., 0.], [1., 0.], [2., 0.], [0., 1.], [1., 1.], [2., 1.], [1.5, 2.]].iter().enumerate() {
            m.set_node(i as u32 + 1, [xy[0], xy[1], 0.0]).unwrap();
        }
        m.set_element(1, crate::ElementKind::Quad4, vec![1, 2, 5, 4]).unwrap();
        m.set_element(2, crate::ElementKind::Quad4, vec![2, 3, 6, 5]).unwrap();
        m.set_element(3, crate::ElementKind::Tri3, vec![5, 6, 7]).unwrap();
        let edges = m.exterior_edges();
        let of = |e: u32| edges.iter().filter(|f| f.element == e).map(|f| f.edge).collect::<Vec<_>>();
        assert_eq!(of(1), [1, 3, 4]);
        assert_eq!(of(2), [1, 2]);
        assert_eq!(of(3), [2, 3]);
    }

    #[test]
    fn a_model_holds_a_bounded_number_of_contacts() {
        let mut m = model();
        for i in 0..MAX_CONTACTS {
            m.add_contact(&format!("c{i}"), Contact::new(ContactSurfaces::General { surface: None })).unwrap();
        }
        assert!(matches!(m.add_contact("one-more", Contact::new(ContactSurfaces::General { surface: None })), Err(Error::Limit(_))));
    }
}
