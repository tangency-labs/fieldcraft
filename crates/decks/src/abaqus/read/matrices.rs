//! Public MATRIX GENERATE / MATRIX OUTPUT syntax; only unloaded linear global operators.
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-matrixgenerate.htm
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-matrixoutput.htm
use super::history::{allow, unsupported, value};
use super::{Block, canonical, error};
use crate::Result;
use fieldcraft_model::{Procedure, steps::MatrixKind};
use std::collections::BTreeSet;

pub(super) fn procedure(generate: &Block<'_>, body: &[Block<'_>]) -> Result<Procedure> {
    let flags = [
        ("STIFFNESS", MatrixKind::Stiffness),
        ("MASS", MatrixKind::Mass),
        ("VISCOUSDAMPING", MatrixKind::ViscousDamping),
        ("STRUCTURALDAMPING", MatrixKind::StructuralDamping),
    ];
    allow(generate, &["STIFFNESS", "MASS", "VISCOUSDAMPING", "STRUCTURALDAMPING", "MPC"])?;
    if canonical(value(generate, "MPC", "YES")?) != "YES" || !generate.rows.is_empty() {
        return Err(unsupported(generate.line, "matrix generation requires default constraint elimination and no data rows"));
    }
    let mut matrices = BTreeSet::new();
    for (key, kind) in flags {
        if let Some(v) = generate.params.get(key) {
            if v.is_some() {
                return Err(error(generate.line, "matrix operator selectors are flags"));
            }
            matrices.insert(kind);
        }
    }
    if matrices.is_empty() {
        return Err(error(generate.line, "no matrix requested"));
    }
    let outputs: Vec<_> = body.iter().filter(|b| b.keyword == "MATRIXOUTPUT").collect();
    if outputs.len() > 1 {
        return Err(error(generate.line, "multiple matrix output requests"));
    }
    if let Some(output) = outputs.first() {
        allow(output, &["STIFFNESS", "MASS", "VISCOUSDAMPING", "STRUCTURALDAMPING", "FORMAT"])?;
        if canonical(value(output, "FORMAT", "MATRIXINPUT")?) != "MATRIXINPUT" || !output.rows.is_empty() {
            return Err(unsupported(output.line, "only MATRIX INPUT output is mapped"));
        }
        for (key, kind) in flags {
            if let Some(v) = output.params.get(key)
                && (v.is_some() || !matrices.contains(&kind))
            {
                return Err(error(output.line, "output operator must be a generated flag"));
            }
        }
    }
    Ok(Procedure::MatrixGeneration { matrices })
}
