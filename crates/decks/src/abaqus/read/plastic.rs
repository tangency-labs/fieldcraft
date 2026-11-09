//! Constant isotropic hardening from the public *PLASTIC keyword reference:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-plastic.htm
//! Stress/plastic-strain rows map to the existing native piecewise-linear J2 table.
use super::{Block, Result, canonical, error, number, quantity, row_fields};
use fieldcraft_model::units::Quantity;
use fieldcraft_model::{DisplaySystem, ElementKind, Model, Plastic};

pub(super) fn material(b: &Block<'_>, units: &DisplaySystem) -> Result<Plastic> {
    b.only(&["HARDENING"])?;
    let kind = if b.params.contains_key("HARDENING") { b.required("HARDENING")? } else { "ISOTROPIC" };
    if canonical(kind) != "ISOTROPIC" {
        return Err(error(b.line, format!("*PLASTIC HARDENING={kind} is not imported; only constant ISOTROPIC tables are mapped")));
    }
    if b.rows.is_empty() || b.rows.len() > fieldcraft_model::materials::MAX_TABLE_POINTS {
        return Err(error(b.line, "*PLASTIC needs a bounded nonempty hardening table"));
    }
    let mut points = Vec::with_capacity(b.rows.len());
    for row in &b.rows {
        let f = row_fields(row.text, row.line)?;
        if f.len() != 2 {
            return Err(error(row.line, "*PLASTIC requires yield stress and plastic strain only; temperature/field dependence is not imported"));
        }
        points.push(fieldcraft_model::materials::HardeningPoint {
            stress: quantity(units, Quantity::Pressure, number(f.first().map(String::as_str).unwrap_or_default(), row.line)?, row.line)?,
            strain: number(f.get(1).map(String::as_str).unwrap_or_default(), row.line)?,
        });
    }
    if points.windows(2).any(|p| p[1].stress < p[0].stress) {
        return Err(error(b.line, "*PLASTIC softening table is not imported"));
    }
    let plastic = Plastic::Table { points };
    plastic.validate().map_err(|e| error(b.line, e.to_string()))?;
    Ok(plastic)
}

pub(super) fn scope(step: &Block<'_>, body: &[Block<'_>], model: &Model, nonlinear_geometry: bool) -> Result<()> {
    use super::history::{unsupported, value};
    if nonlinear_geometry || step.params.contains_key("PERTURBATION") || !body.iter().any(|b| b.keyword == "STATIC") {
        return Err(unsupported(step.line, "PLASTIC requires one small-strain STATIC step with NLGEOM=NO, without PERTURBATION"));
    }
    if canonical(value(step, "AMPLITUDE", "RAMP")?) != "RAMP" {
        return Err(unsupported(step.line, "PLASTIC requires STEP AMPLITUDE=RAMP"));
    }
    if model.elements().is_empty()
        || model.sections().is_empty()
        || model.elements().values().any(|e| e.kind != ElementKind::Tet4)
        || model.sections().values().any(|s| s.orientation.is_some())
    {
        return Err(unsupported(step.line, "PLASTIC import requires isotropic C3D4 without ORIENTATION"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn unassigned_plastic_material_does_not_create_a_runnable_empty_equilibrium() {
        let text = "*NODE\n1,0,0,0\n*MATERIAL,NAME=M\n*ELASTIC\n1000,.3\n*PLASTIC\n10,0\n*STEP\n*STATIC\n*END STEP\n";
        let read = super::super::read_analysis_limited(text, fieldcraft_model::units::Preset::MKgS, super::super::Limits::default()).unwrap();
        assert!(!read.analysis_complete);
        assert!(read.model.steps().is_empty());
        assert!(read.warnings.iter().any(|w| w.contains("PLASTIC") && w.contains("C3D4")));
    }
    #[test]
    fn table_units_and_unmapped_options_are_not_silently_defaulted() {
        let text = "*MATERIAL,NAME=M\n*ELASTIC\n.001,.25\n*PLASTIC\n.000001,0\n.000002,.01\n";
        let got = super::super::read_in(text, fieldcraft_model::units::Preset::MmTonneS).unwrap();
        let plastic = got.model.material("M").unwrap().plastic.as_ref().unwrap();
        assert_eq!(plastic.initial_yield(), 1.);
        for (from, to) in [
            ("*PLASTIC", "*PLASTIC,HARDENING"),
            ("*PLASTIC", "*PLASTIC,HARDENING=KINEMATIC"),
            ("*PLASTIC", "*PLASTIC,RATE=1"),
            (".000001,0", ".000001,0,20"),
            (".000001,0", ".000001,.001"),
            (".000002,.01", ".0000005,.01"),
        ] {
            assert!(super::super::read_in(&text.replace(from, to), fieldcraft_model::units::Preset::MmTonneS).is_err());
        }
    }
}
