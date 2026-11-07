//! Linear laminate export from the public LS-DYNA R14 manuals, Vol. I PART_COMPOSITE
//! (IRPL=103), Vol. II MAT_021 and its material-direction rules. No vendor code or decks.
//! https://lsdyna.ansys.com/manuals-download/
use super::{Deck, Error, FORMAT, Field, Part, Plan, Result, WIDTH, broken};
use fieldcraft_model::units::{Quantity, Unit};
use fieldcraft_model::{Elastic, Layup, Material, MaterialOrientation, Model, Resolved, SectionKind};
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn check(model: &Model) -> Result<()> {
    let mut bad = BTreeSet::new();
    for (name, m) in model.materials() {
        if m.ply_strength.is_some() {
            bad.insert(format!(
                "ply strengths of material `{name}`: MAT_021's optional failure polynomial does not preserve all four model criteria"
            ));
        }
    }
    for (name, s) in model.sections() {
        if matches!(s.kind, SectionKind::Solid {}) && s.orientation.is_some_and(|o| matches!(o.frame, MaterialOrientation::Cylindrical { .. })) {
            bad.insert(format!("material orientation of section `{name}`: cylindrical solid directions evaluated at each integration point cannot be replaced by an element-centre frame"));
        }
        if s.orientation.is_some() || matches!(s.kind, SectionKind::CompositeShell { .. }) {
            for material in s.material_names() {
                if model
                    .material(material)
                    .is_some_and(|m| m.plastic.is_some() || m.hyperelastic.is_some() || m.damage.is_some() || m.erosion.is_some())
                {
                    bad.insert(format!("section `{name}` material `{material}`: oriented MAT_021 output supports linear elasticity only"));
                }
            }
        }
        if let SectionKind::CompositeShell { layup } = &s.kind {
            if layup.plies.iter().any(|p| p.points != 3) {
                bad.insert(format!("composite shell section `{name}`: PART_COMPOSITE IRPL=103 requires three section points per ply; its default midpoint rule loses each ply's bending moment"));
            }
            let densities: Vec<_> = layup.plies.iter().map(|p| model.material(&p.material).and_then(|m| m.density)).collect();
            if densities.windows(2).any(|p| p.first() != p.get(1)) {
                bad.insert(format!("composite shell section `{name}`: unequal ply densities need independently verified I1/I2 mass moments; thickness-averaged density is not a lossless substitute"));
            }
        }
    }
    crate::unsupported(FORMAT, bad)
}

/// MAT orientation belongs to the material in this format, but to the section in the
/// model. Give each oriented section its own material definitions. Leave load sets intact.
pub(super) fn prepare(model: &Model, r: &Resolved) -> Result<Option<Model>> {
    let explicit = crate::composite::explicit_shell_frames(model, r, FORMAT)?;
    let source = explicit.as_ref().unwrap_or(model);
    let selected: Vec<_> = source.sections().iter().filter(|(_, s)| s.orientation.is_some()).map(|(n, s)| (n.clone(), s.clone())).collect();
    if selected.is_empty() {
        return Ok(explicit);
    }
    let additional = selected
        .iter()
        .try_fold(0usize, |n, (_, s)| n.checked_add(s.material_names().count()).ok_or_else(|| Error::Limit("too many material frames".into())))?;
    if source.materials().len().saturating_add(additional) > fieldcraft_model::MAX_DEFINITIONS {
        return Err(Error::Limit("section material frames exceed the material definition cap".into()));
    }
    let mut out = source.clone();
    let mut names: BTreeSet<_> = source.materials().keys().map(|n| n.to_ascii_uppercase()).collect();
    let mut next = 0usize;
    for (name, mut section) in selected {
        fieldcraft_model::progress::check()?;
        let mut replacements = BTreeMap::new();
        for old in section.material_names() {
            if replacements.contains_key(old) {
                continue;
            }
            let generated = loop {
                next += 1;
                if next > 2 * fieldcraft_model::MAX_DEFINITIONS {
                    return Err(Error::Limit("cannot name section materials".into()));
                }
                let candidate = format!("FC_LS_M_{next}");
                if names.insert(candidate.clone()) {
                    break candidate;
                }
            };
            let material = source.material(old).ok_or_else(|| broken(format!("material `{old}` is missing")))?;
            out.add_material(&generated, material.clone())?;
            replacements.insert(old.to_owned(), generated);
        }
        match &mut section.kind {
            SectionKind::CompositeShell { layup } => {
                for p in &mut layup.plies {
                    p.material = replacements.get(&p.material).cloned().ok_or_else(|| broken("missing ply material replacement".into()))?;
                }
            }
            _ => {
                section.material =
                    replacements.get(&section.material).cloned().ok_or_else(|| broken("missing section material replacement".into()))?
            }
        }
        out.replace_section(&name, section)?;
    }
    Ok(Some(out))
}

pub(super) fn frames(model: &Model) -> Result<BTreeMap<&str, [[f64; 3]; 3]>> {
    let mut out = BTreeMap::new();
    for s in model.sections().values() {
        let Some(o) = s.orientation else { continue };
        let axes = o.frame.axes_at([0.0; 3])?;
        let i = if matches!(s.kind, SectionKind::Shell { .. } | SectionKind::CompositeShell { .. }) {
            usize::from(o.shell_reference_axis - 1)
        } else {
            0
        };
        let first = axes.get(i).copied().ok_or_else(|| broken("invalid shell reference axis".into()))?;
        let rotated = [first, axes[(i + 1) % 3], axes[(i + 2) % 3]];
        for name in s.material_names() {
            out.insert(name, rotated);
        }
    }
    Ok(out)
}

pub(super) fn material(d: &mut Deck, name: &str, mid: u32, m: &Material, axes: [[f64; 3]; 3]) -> Result<()> {
    let (e, nu, g) = match m.elastic {
        Some(Elastic::EngineeringConstants { young, poisson, shear }) => (young, poisson, shear),
        Some(Elastic::Isotropic { young, poisson }) => ([young; 3], [poisson; 3], [young / (2.0 * (1.0 + poisson)); 3]),
        None => return Err(broken(format!("material `{name}` has no elastic constants"))),
    };
    d.line("*MAT_ORTHOTROPIC_THERMAL_TITLE");
    d.line(name);
    d.header(&["mid", "ro", "ea", "eb", "ec", "prba", "prca", "prcb"]);
    d.card(
        [Field::from(mid), Field::from(m.density.unwrap_or_default())]
            .into_iter()
            .chain([e[0], e[1], e[2], nu[0] * e[1] / e[0], nu[1] * e[2] / e[0], nu[2] * e[2] / e[1]].map(Field::from)),
    );
    // Plain MAT_021, zero expansion: strictly linear elasticity, no damage or curing.
    d.header(&["gab", "gbc", "gca", "aa", "ab", "ac", "aopt", "macf"]);
    d.card([g[0], g[2], g[1], 0.0, 0.0, 0.0, 2.0].map(Field::from).into_iter().chain([Field::from(1)]));
    d.header(&["xp", "yp", "zp", "a1", "a2", "a3"]);
    d.card([const { Field::Blank }; 3].into_iter().chain(axes[0].map(Field::from)));
    d.header(&["v1", "v2", "v3", "d1", "d2", "d3", "beta", "ref"]);
    d.card([const { Field::Blank }; 3].into_iter().chain(axes[1].map(Field::from)));
    Ok(())
}

pub(super) fn part(d: &mut Deck, p: &Part<'_>, layup: &Layup, plan: &Plan<'_>) -> Result<()> {
    d.line("*PART_COMPOSITE_LONG");
    d.line(p.name);
    d.line(&format!("{:<WIDTH$}{:>WIDTH$}", "OPTCARD", 103));
    d.header(&["pid", "elform", "shrf", "nloc", "marea", "hgid"]);
    d.card([
        Field::from(p.id),
        Field::from(p.elform),
        Field::from(5.0 / 6.0),
        Field::from(2.0 * layup.offset),
        Field::Blank,
        p.hourglass.map_or(Field::Blank, Field::from),
    ]);
    d.header(&["mid", "thick", "angle", "tmid", "plyid"]);
    for (i, ply) in layup.plies.iter().enumerate() {
        let mid = plan.materials.get(ply.material.as_str()).copied().ok_or_else(|| broken("missing ply MID".into()))?;
        let angle = Unit::parse("deg")
            .and_then(|u| Quantity::Angle.express(ply.angle, &u))
            .map_err(|e| Error::Model(fieldcraft_model::Error::Invalid(e.to_string())))?;
        let ply_id = super::id(i)?;
        d.card([Field::from(mid), Field::from(ply.thickness), Field::from(angle), Field::Blank, Field::from(ply_id)]);
    }
    Ok(())
}
