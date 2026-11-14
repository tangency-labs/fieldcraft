//! Point masses, rotary inertia, springs and dashpots in the implicit system.
//!
//! The model resolves each definition to its nodes ([`fieldcraft_model::Resolved::masses`] and
//! `springs`); here each copy is turned into global degrees of freedom and its matrix:
//!
//! - a point mass adds m_x, m_y, m_z on the node's translations and its inertia tensor J (3 × 3,
//!   off-diagonal products included) on its rotations, the same in the consistent and the lumped
//!   mass matrix;
//! - a spring with stiffness k and unit direction coefficients c on each end's six degrees of
//!   freedom adds k c cᵀ to each end and −k c cᵀ between them (Cook, Malkus, Plesha, Witt 2002,
//!   §2.2): the elongation is cᵀ(u_second − u_first);
//! - a dashpot has the same form on velocities; static and frequency steps have no velocities, so
//!   they ignore dashpots with a warning (damping enters dynamic steps).
//!
//! What each copy needs at its nodes ([`needs`]) goes into the mesh's [`fieldcraft_fem::DofLayout`]
//! like an element's declared kinds; [`dof_of`] then finds its unknowns there.

use fieldcraft_fem::{DofKind, DofSet};
use fieldcraft_model::{Resolved, ResolvedSpring};

use crate::sparse::{Csr, add_entry};
use crate::system::Mesh;
use crate::{Error, Result};

/// A point mass on global degrees of freedom.
#[derive(Debug, Clone, PartialEq)]
pub struct MassAt {
    pub name: String,
    /// Translations along X, Y, Z: (degree of freedom, mass).
    pub translations: Vec<(usize, f64)>,
    /// The rotations' degrees of freedom and the inertia tensor on them (empty without inertia).
    pub rotations: Vec<usize>,
    pub inertia: [[f64; 3]; 3],
}

/// A spring or dashpot on global degrees of freedom: the elongation is Σ c u over `terms`, whose
/// coefficients are the direction's on the second end and their negatives on the first.
#[derive(Debug, Clone, PartialEq)]
pub struct SpringAt {
    pub name: String,
    pub dashpot: bool,
    /// Connected to the ground (one end only).
    pub grounded: bool,
    pub terms: Vec<(usize, f64)>,
    pub value: f64,
}

/// Every discrete element of the model on the mesh's degrees of freedom.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Discrete {
    pub masses: Vec<MassAt>,
    pub springs: Vec<SpringAt>,
}

/// The kinds each discrete copy needs at the labels it acts on: translations for a mass,
/// rotations for an inertia, translations or rotations for a spring or dashpot.
/// `translations` are the model's (X and Y in a 2D model).
pub fn needs(resolved: &Resolved, translations: DofSet) -> impl Iterator<Item = (u32, DofSet)> + '_ {
    let masses = resolved.masses.iter().flat_map(move |m| {
        let mass = m.mass.iter().any(|c| *c != 0.0).then_some((m.node, translations));
        let inertia = m.inertia.iter().any(|c| *c != 0.0).then_some((m.node, DofSet::ROTATIONS));
        mass.into_iter().chain(inertia)
    });
    // A node of a 3D model moves in all three directions (a spring along X does not take Y and
    // Z away), so a spring asks for every translation, or every rotation.
    let springs = resolved.springs.iter().flat_map(move |s| {
        let kinds = if s.direction.is_rotational() { DofSet::ROTATIONS } else { translations };
        std::iter::once((s.first, kinds)).chain(s.second.map(|b| (b, kinds)))
    });
    masses.chain(springs)
}

/// The global degree of freedom of direction `d` (0–2 translations, 3–5 rotations) of a node or
/// reference point label, or `None` when the node has no such direction.
pub fn dof_of(mesh: &Mesh, label: u32, d: usize) -> Option<usize> {
    mesh.layout.dof(mesh.layout.slot(label)?, DofKind::from_direction(d)?)
}

fn missing(what: &str, name: &str, label: u32, d: usize) -> Error {
    let direction = DofKind::from_direction(d).map_or("?", DofKind::describe);
    Error::Invalid(format!("{what} `{name}`: node {label} has no degree of freedom in {direction}"))
}

/// The spring's terms at one end, with the direction's coefficients times `sign`.
fn ends(mesh: &Mesh, s: &ResolvedSpring, label: u32, sign: f64) -> Result<Vec<(usize, f64)>> {
    s.direction
        .coefficients()
        .iter()
        .enumerate()
        .filter(|(_, c)| **c != 0.0)
        .map(|(d, c)| dof_of(mesh, label, d).map(|dof| (dof, sign * c)).ok_or_else(|| missing(s.kind(), &s.name, label, d)))
        .collect()
}

impl Discrete {
    /// The resolved discrete elements on the mesh's degrees of freedom.
    pub fn new(mesh: &Mesh, resolved: &Resolved) -> Result<Discrete> {
        let mut out = Discrete::default();
        for m in &resolved.masses {
            // The model's translations only (a 2D model's masses act in its plane).
            let planar = mesh.translation_kinds();
            let translations: Vec<(usize, f64)> = (0..3)
                .filter(|d| m.mass[*d] != 0.0 && DofKind::translation(*d).is_some_and(|k| planar.contains(k)))
                .map(|d| dof_of(mesh, m.node, d).map(|dof| (dof, m.mass[d])).ok_or_else(|| missing("point mass", &m.name, m.node, d)))
                .collect::<Result<_>>()?;
            let inertia = fieldcraft_model::discrete::inertia_matrix(m.inertia);
            let rotations = if m.inertia.iter().any(|c| *c != 0.0) {
                (3..6).map(|d| dof_of(mesh, m.node, d).ok_or_else(|| missing("point mass", &m.name, m.node, d))).collect::<Result<_>>()?
            } else {
                Vec::new()
            };
            out.masses.push(MassAt { name: m.name.clone(), translations, rotations, inertia });
        }
        for s in &resolved.springs {
            let mut terms = ends(mesh, s, s.first, -1.0)?;
            if let Some(b) = s.second {
                terms.extend(ends(mesh, s, b, 1.0)?);
            }
            out.springs.push(SpringAt { name: s.name.clone(), dashpot: s.dashpot, grounded: s.second.is_none(), terms, value: s.value });
        }
        Ok(out)
    }

    /// The groups of degrees of freedom each element couples, for the matrix pattern.
    pub fn couplings(&self) -> impl Iterator<Item = Vec<usize>> + '_ {
        let masses = self.masses.iter().map(|m| m.translations.iter().map(|t| t.0).chain(m.rotations.iter().copied()).collect());
        masses.chain(self.springs.iter().map(|s| s.terms.iter().map(|t| t.0).collect()))
    }

    /// Every degree of freedom a discrete element acts on.
    pub fn dofs(&self) -> impl Iterator<Item = usize> + '_ {
        self.couplings().flatten()
    }

    /// Names of the dashpots, which static and frequency steps ignore.
    pub fn dashpots(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.springs.iter().filter(|s| s.dashpot).map(|s| s.name.as_str()).collect();
        names.dedup();
        names
    }

    /// Add the springs' stiffness.
    pub fn add_stiffness(&self, k: &mut Csr) -> Result<()> {
        for s in self.springs.iter().filter(|s| !s.dashpot) {
            for (i, ci) in &s.terms {
                for (j, cj) in &s.terms {
                    add_entry(k, *i, *j, s.value * ci * cj)?;
                }
            }
        }
        Ok(())
    }

    /// The sum of the forces springs to the ground exert on the model, along X, Y and Z (the
    /// loads they carry to the ground, besides the boundary conditions' reactions).
    pub fn ground_force(&self, u: &[f64], mesh: &Mesh) -> [f64; 3] {
        let mut total = [0.0; 3];
        for s in self.springs.iter().filter(|s| !s.dashpot && s.grounded) {
            // The spring's force on the model: −k (cᵀu) c over its terms.
            let stretch: f64 = s.terms.iter().map(|(d, c)| c * u.get(*d).copied().unwrap_or(0.0)).sum();
            for (d, c) in &s.terms {
                if let Some(t) = mesh.kind(*d).and_then(DofKind::direction).filter(|d| *d < 3).and_then(|d| total.get_mut(d)) {
                    *t -= s.value * stretch * c;
                }
            }
        }
        total
    }

    /// Gravity on the point masses: m_d g_d on each translation.
    pub fn add_gravity(&self, f: &mut [f64], g: [f64; 3], mesh: &Mesh) {
        for p in &self.masses {
            for (dof, m) in &p.translations {
                if let (Some(x), Some(a)) = (f.get_mut(*dof), mesh.kind(*dof).and_then(DofKind::direction).and_then(|d| g.get(d))) {
                    *x += m * a;
                }
            }
        }
    }
}
