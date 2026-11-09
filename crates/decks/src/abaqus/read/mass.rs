//! Scope-aware mass requirements. These describe evidence, never silently select a solver option.
use super::{MAX_RECORDS, canonical, lex};
use crate::Result;
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MassRequirement {
    pub element_type: String,
    /// ELSET on the source ELEMENT card, not a resolved section assignment.
    pub elset: Option<String>,
    /// The existing consistent/lumped contract, if the formula family is identified.
    /// This is NOT an applied setting or a statement that native lumping uses that formula.
    pub choice: Option<&'static str>,
    pub rule: &'static str,
    pub conditions: &'static str,
    pub native_gap: &'static str,
}

pub(super) fn inspect(text: &str) -> Result<Vec<MassRequirement>> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for b in lex::blocks(text, MAX_RECORDS)?.iter().filter(|b| b.keyword == "ELEMENT") {
        let element_type = canonical(b.required("TYPE")?);
        let elset = b.value("ELSET").map(|s| lex::unquote(s, b.line).map(|s| s.to_ascii_uppercase())).transpose()?;
        if !seen.insert((element_type.clone(), elset.clone())) {
            continue;
        }
        let (choice, rule, conditions, native_gap) = classify(&element_type);
        out.push(MassRequirement { element_type, elset, choice, rule, conditions, native_gap });
    }
    Ok(out)
}
fn classify(t: &str) -> (Option<&'static str>, &'static str, &'static str, &'static str) {
    match t {
        "C3D5" | "C3D5H" => (None, "unknown_pyramid_distribution", "distorted pyramid; uniform density", "U-PYR: neither row-sum nor HRZ identified"),
        "S3" | "S3R" | "S4" | "S4R" => (
            None,
            "unknown_shell_projection",
            "flat constant-thickness shell; zero offset",
            "U-SHELL: full translational/rotational distribution unknown",
        ),
        "B31" => (
            Some("lumped"),
            "conditional_row_sum_hrz_equivalence",
            "straight uniform prismatic beam; no offset; constant orientation",
            "verify rotational inertia and section conditions; no automatic selection",
        ),
        "C3D4" | "C3D4H" | "CPS3" | "CPE3" | "T2D2" | "T3D2" => (
            Some("lumped"),
            "conditional_row_sum_hrz_equivalence",
            "affine simplex or straight constant-area rod; uniform density; plane constant thickness",
            "conditions not automatically verified; no automatic selection",
        ),
        "C3D6" | "C3D6H" | "C3D8" | "C3D8H" | "C3D8I" | "C3D8IH" | "C3D8R" | "C3D8RH" | "CAX3" | "CAX4" | "CAX4I" | "CAX4R" | "CPE4" | "CPE4I"
        | "CPE4R" | "CPS4" | "CPS4I" | "CPS4R" => (
            Some("lumped"),
            "row_sum",
            "uniform density; tested geometry and section; axisymmetric full ring",
            "native default is HRZ; explicit section row_sum is available for mass=lumped, but source applicability is not automatically verified",
        ),
        "CAX6" => (
            Some("consistent"),
            "triangle_six_point_degree_four",
            "straight edges with midpoint nodes; uniform density; full axisymmetric ring",
            "native TriConical(4) is not the identified six-point quadrature",
        ),
        "C3D20" | "C3D20H" | "C3D20R" | "C3D20RH" | "CAX8" | "CAX8R" => (
            Some("consistent"),
            "consistent_gauss_three",
            "uniform density; tested geometry; axisymmetric full ring",
            "quadrature family available; applicability and mixed-type selection not verified",
        ),
        "C3D10" | "C3D10H" | "CPS6" | "CPE6" | "CPS8" | "CPS8R" | "CPE8" | "CPE8R" | "T2D3" | "T3D3" => (
            Some("consistent"),
            "consistent_interpolation",
            "uniform density; tested geometry and section",
            "formula family available; applicability and mixed-type selection not verified",
        ),
        "S8R" => (
            Some("consistent"),
            "consistent_translation_only",
            "flat constant-thickness shell; zero offset",
            "full rotational/drilling inertia rule not established",
        ),
        "MASS" => (None, "specified_point_mass", "explicit mass property", "not a distributed mass selection"),
        "SPRING1" => (
            None,
            "zero_incremental_mass",
            "grounded spring; difference from identical point-mass control",
            "not a direct standalone native mass export",
        ),
        _ => (None, "unknown", "outside established mass evidence", "no mass-default mapping"),
    }
}

pub(super) const POLICY_MARKER: &str = "** FIELDCRAFT LUMPED MASS,";

/// Native-only policy metadata. Abaqus treats this as a comment, never a vendor keyword.
pub(super) fn restore(text: &str, model: &mut fieldcraft_model::Model) -> Result<()> {
    let mut assigned = BTreeSet::new();
    for (index, line) in text.lines().enumerate() {
        let Some(params) = line.trim().strip_prefix(POLICY_MARKER) else { continue };
        let synthetic = format!("*FIELDCRAFTLUMPEDMASS,{params}\n");
        let blocks = lex::blocks(&synthetic, MAX_RECORDS)?;
        let b = blocks.first().ok_or_else(|| lex::error(index + 1, "empty mass policy"))?;
        b.only(&["ELSET", "STRATEGY"]).map_err(|e| lex::error(index + 1, e.to_string()))?;
        let set = lex::unquote(b.required("ELSET").map_err(|e| lex::error(index + 1, e.to_string()))?, index + 1)?;
        let strategy = match canonical(b.required("STRATEGY").map_err(|e| lex::error(index + 1, e.to_string()))?).as_str() {
            "HRZ" => fieldcraft_model::LumpedMass::Hrz,
            "ROWSUM" | "ROW_SUM" => fieldcraft_model::LumpedMass::RowSum,
            "CONSISTENT" => fieldcraft_model::LumpedMass::Consistent,
            _ => return Err(lex::error(index + 1, "unknown FIELDCRAFT LUMPED MASS STRATEGY")),
        };
        let names: Vec<_> = model.sections().iter().filter(|(_, s)| s.set.eq_ignore_ascii_case(&set)).map(|(n, _)| n.clone()).collect();
        if names.is_empty() {
            return Err(lex::error(index + 1, format!("mass policy ELSET {set} has no section")));
        }
        for name in names {
            if !assigned.insert(name.clone()) {
                return Err(lex::error(index + 1, format!("duplicate mass policy for section {name}")));
            }
            let mut section = model.sections().get(&name).cloned().ok_or_else(|| lex::error(index + 1, "missing section"))?;
            section.formulation.lumped_mass = strategy;
            model.replace_section(&name, section).map_err(|e| lex::error(index + 1, e.to_string()))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fifty_one_types_have_scoped_rules_and_six_named_unknowns() {
        let types = "B31 C3D10 C3D10H C3D20 C3D20H C3D20R C3D20RH C3D4 C3D4H C3D5 C3D5H C3D6 C3D6H C3D8 C3D8H C3D8I C3D8IH C3D8R C3D8RH CAX3 CAX4 CAX4I CAX4R CAX6 CAX8 CAX8R CPE3 CPE4 CPE4I CPE4R CPE6 CPE8 CPE8R CPS3 CPS4 CPS4I CPS4R CPS6 CPS8 CPS8R MASS S3 S3R S4 S4R S8R SPRING1 T2D2 T2D3 T3D2 T3D3";
        let deck = types.split_whitespace().map(|t| format!("*ELEMENT,TYPE={t},ELSET={t}\n")).collect::<String>();
        let report = inspect(&deck).unwrap();
        assert_eq!(report.len(), 51);
        assert_eq!(report.iter().filter(|r| r.rule.starts_with("unknown")).count(), 6);
        assert_eq!(report.iter().filter(|r| r.rule == "conditional_row_sum_hrz_equivalence").count(), 7);
        assert!(report.iter().all(|r| r.choice.is_none_or(|c| matches!(c, "consistent" | "lumped"))));
        assert_eq!(classify("NOT_AN_ELEMENT").1, "unknown");
    }
    #[test]
    fn mixed_types_need_different_contract_choices_even_in_one_set() {
        let report = inspect("*ELEMENT,TYPE=C3D8,ELSET=Mixed\n*ELEMENT,TYPE=C3D10,ELSET=Mixed\n").unwrap();
        let choices: BTreeSet<_> = report.iter().filter_map(|r| r.choice).collect();
        assert_eq!(choices, BTreeSet::from(["consistent", "lumped"]));
        assert_eq!(report[0].elset, report[1].elset);
    }
}
