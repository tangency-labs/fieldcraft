use super::lex::{canonical, error};
use crate::Result;
use fieldcraft_model::{ElementKind as K, Formulation, Integration, SectionKind};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Family {
    Solid,
    Shell,
    Stress,
    Strain,
    Axis,
    Beam,
    Rod,
    PlanarRod,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Type {
    pub kind: K,
    pub family: Family,
    pub formulation: Formulation,
}

impl Type {
    pub fn parse(name: &str, line: usize) -> Result<Self> {
        let normalized = canonical(name);
        let (family, tail) = if let Some(t) = normalized.strip_prefix("C3D") {
            (Family::Solid, t)
        } else if let Some(t) = normalized.strip_prefix("CPS") {
            (Family::Stress, t)
        } else if let Some(t) = normalized.strip_prefix("CPE") {
            (Family::Strain, t)
        } else if let Some(t) = normalized.strip_prefix("CAX") {
            (Family::Axis, t)
        } else if let Some(t) = normalized.strip_prefix("T3D") {
            (Family::Rod, t)
        } else if let Some(t) = normalized.strip_prefix("T2D") {
            (Family::PlanarRod, t)
        } else if normalized == "B31" {
            (Family::Beam, "2")
        } else if let Some(t) = normalized.strip_prefix('S') {
            (Family::Shell, t)
        } else {
            return Err(error(line, format!("unsupported element type `{name}`")));
        };
        let hybrid = tail.ends_with('H');
        let tail = if hybrid { tail.strip_suffix('H').unwrap_or(tail) } else { tail };
        let (count, reduced, incompatible) = if let Some(t) = tail.strip_suffix('R') {
            (t, true, false)
        } else if let Some(t) = tail.strip_suffix('I') {
            (t, false, true)
        } else {
            (tail, false, false)
        };
        let kind = match (family, count) {
            (Family::Solid, "4") => K::Tet4,
            (Family::Solid, "10") => K::Tet10,
            (Family::Solid, "6") => K::Wedge6,
            (Family::Solid, "5") => K::Pyramid5,
            (Family::Solid, "8") => K::Hex8,
            (Family::Solid, "20") => K::Hex20,
            (Family::Shell | Family::Stress | Family::Strain | Family::Axis, "3") => K::Tri3,
            (Family::Stress | Family::Strain | Family::Axis, "6") => K::Tri6,
            (Family::Shell | Family::Stress | Family::Strain | Family::Axis, "4") => K::Quad4,
            (Family::Shell | Family::Stress | Family::Strain | Family::Axis, "8") => K::Quad8,
            (Family::Rod | Family::PlanarRod, "2") => K::Rod2,
            (Family::Rod | Family::PlanarRod, "3") => K::Rod3,
            (Family::Beam, "2") => K::Beam2,
            _ => return Err(error(line, format!("unsupported element type `{name}`"))),
        };
        let out = Self {
            kind,
            family,
            formulation: Formulation {
                integration: if reduced { Integration::Reduced } else { Integration::Full },
                incompatible_modes: incompatible,
                hybrid,
            },
        };
        let section =
            fieldcraft_model::Section { set: "Set".into(), material: "Material".into(), kind: out.section(1.0), formulation: out.formulation };
        let written = super::super::element_type(kind, &section).map_err(|_| error(line, format!("unsupported element type `{name}`")))?;
        let written = if family == Family::PlanarRod { written.replacen("T3D", "T2D", 1) } else { written };
        if written != normalized {
            return Err(error(line, format!("unsupported element type `{name}`")));
        }
        Ok(out)
    }
    pub fn section(self, thickness: f64) -> SectionKind {
        match self.family {
            Family::Solid => SectionKind::Solid {},
            Family::Shell => SectionKind::Shell { thickness, offset: 0.0, points: 5 },
            Family::Stress => SectionKind::PlaneStress { thickness },
            Family::Strain => SectionKind::PlaneStrain { thickness },
            Family::Axis => SectionKind::Axisymmetric {},
            Family::Rod | Family::PlanarRod => SectionKind::Rod { area: thickness },
            Family::Beam => SectionKind::Beam { profile: fieldcraft_model::BeamProfile::Circle { radius: thickness }, orientation: [0., 0., -1.] },
        }
    }
}
