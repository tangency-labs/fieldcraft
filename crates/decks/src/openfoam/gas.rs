//! A compressible case's gas and dictionaries (rhoCentralFoam).
use std::fmt::Write as _;

use fieldcraft_model::units::{Quantity, Typed};
use fieldcraft_model::{BoxFace, Directional, FaceCondition, Fluid, FluidDomain, Model, SurfaceCondition};

use super::{Surface, dictionary, field_file, num, unsupported, vector};
use crate::{Error, Result};

/// The density-based compressible solver (UG §A.1; DOC "rhoCentralFoam": fields p [Pa], U, T [K]).
pub const SOLVER: &str = "rhoCentralFoam";
/// The molar gas constant R = N_A k, J/(mol K): exact since the 2019 SI (BIPM SI Brochure, 9th ed.).
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;
/// Default Courant limit, 0.2: the CFL number Greenshields et al. (2010, IJNMF 63, §4.1 and §4.2)
/// ran their central-scheme cases at (shock tube and forward step; their figures 2 and 3 use 0.2).
pub const COURANT: f64 = 0.2;
/// Largest Courant limit, ½: the stability limit of the central schemes that Greenshields et al.
/// (2010, §4.1) state ("below the stability limit of 0.5 of the central schemes"). The paper states
/// it and does not derive it; it is not proved here (docs/manual/fluid/openfoam.md).
pub const MAX_COURANT: f64 = 0.5;
/// Largest ratio of specific heats of an ideal gas, 5/3: c_v = c_p − R_s is at least the three
/// translational quadratic terms, 3R_s/2 (equipartition; Kittel & Kroemer, *Thermal Physics*, 2nd
/// ed., ch. 3), so γ = 1 + R_s/c_v ≤ 5/3, with equality for a monatomic gas. A larger γ is noted on
/// the case, not refused: rounded data for a monatomic gas can sit just above it.
pub const MAX_IDEAL_GAMMA: f64 = 5.0 / 3.0;
/// Most initial-pressure boxes: a cap on a hostile input (each is a region of setFieldsDict).
pub const MAX_PRESSURE_BOXES: usize = 1024;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PressureBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
    /// Gauge pressure, Pa.
    pub pressure: f64,
}

/// What a compressible case needs beyond the model (SI): none has a default.
#[derive(Debug, Clone, PartialEq)]
pub struct GasConditions {
    /// Absolute pressure the model's gauge pressures count from, Pa.
    pub ambient: f64,
    /// Initial, inflow and reservoir (total) temperature, K.
    pub temperature: f64,
    /// Specific gas constant R_s, J/(kg K).
    pub gas_constant: f64,
    /// The far-field relaxation length the non-reflecting outlets write as `lInf` (m): the distance
    /// over which an outlet's pressure relaxes towards the ambient.
    pub far_field: f64,
    pub initial_pressure: Vec<PressureBox>,
}

/// The ideal gas of a case: the conditions and the material's properties (SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gas {
    pub ambient: f64,
    pub temperature: f64,
    pub gas_constant: f64,
    pub specific_heat: f64,
    pub viscosity: f64,
    pub conductivity: f64,
    pub far_field: f64,
}

/// Writes `value` as a positive number, or refuses it with `message`.
fn check_positive(value: f64, message: &str) -> Result<()> {
    if value.is_finite() && value > 0.0 { Ok(()) } else { Err(unsupported(message)) }
}

impl Gas {
    /// The gas of a fluid domain's material and the conditions. Every property is read from the
    /// material or the conditions; none is defaulted. The molar-mass-free gas law and the
    /// constant-property transport of OpenFOAM's `hConst`/`const` package (DOC "Constant") need
    /// c_p, μ, k and R_s; a temperature-dependent table of any of them is refused, as the case
    /// holds them constant.
    pub fn of(model: &Model, domain: &FluidDomain, c: &GasConditions) -> Result<Gas> {
        let name = &domain.material;
        let material = model.materials().get(name).ok_or_else(|| unsupported(format!("the fluid domain names no material `{name}`")))?;
        let viscosity = match material.fluid {
            Some(Fluid::Newtonian { viscosity }) if viscosity.is_finite() && viscosity > 0.0 => viscosity,
            Some(Fluid::Newtonian { .. }) => return Err(unsupported(format!("material `{name}`: a viscosity that is not positive"))),
            Some(other) => return Err(unsupported(format!("material `{name}`: {} viscosity (only Newtonian fluids so far)", other.name()))),
            None => return Err(unsupported(format!("material `{name}` has no fluid behaviour"))),
        };
        let thermal = material.thermal.as_ref();
        let specific_heat = thermal
            .and_then(|t| t.specific_heat)
            .ok_or_else(|| unsupported(format!("material `{name}` has no specific heat (a compressible case needs it)")))?;
        let conductivity = match thermal.and_then(|t| t.conductivity.as_ref()) {
            None => return Err(unsupported(format!("material `{name}` has no conductivity (a compressible case needs it)"))),
            Some(Directional::Isotropic(k)) => *k,
            Some(Directional::Orthotropic(_)) => {
                return Err(unsupported(format!("material `{name}`: a conductivity per axis (a gas conducts the same in every direction)")));
            }
        };
        if let Some(t) = &material.temperature_dependent
            && let Some((what, _)) = t.tables().find(|(what, _)| matches!(*what, "conductivity" | "specific heat" | "viscosity"))
        {
            return Err(unsupported(format!("material `{name}`: a {what} that follows the temperature (the case holds it constant)")));
        }
        check_positive(c.ambient, "an ambient pressure that is not positive")?;
        check_positive(c.temperature, "a temperature that is not positive")?;
        check_positive(c.gas_constant, "a gas constant that is not positive")?;
        check_positive(specific_heat, &format!("material `{name}`: a specific heat that is not positive"))?;
        check_positive(conductivity, &format!("material `{name}`: a conductivity that is not positive"))?;
        check_positive(c.far_field, "a far-field length that is not a positive length")?;
        if specific_heat <= c.gas_constant {
            return Err(unsupported("a specific heat at constant pressure that is not above the gas constant (cv = cp - R must be positive)"));
        }
        let gas = Gas {
            ambient: c.ambient,
            temperature: c.temperature,
            gas_constant: c.gas_constant,
            specific_heat,
            viscosity,
            conductivity,
            far_field: c.far_field,
        };
        // Each is a quotient of positive finite numbers, which a tiny divisor still overflows.
        for (what, value) in [
            ("ratio of specific heats", gas.gamma()),
            ("Prandtl number", gas.prandtl()),
            ("speed of sound", gas.sound_speed()),
            ("reference density", gas.reference_density()),
        ] {
            check_positive(
                value,
                &format!(
                    "a {what} that is not finite and positive (the gas constant, specific heat, viscosity, conductivity or temperature is out of range)"
                ),
            )?;
        }
        Ok(gas)
    }

    /// The ratio of specific heats γ = c_p / c_v = c_p / (c_p − R_s), c_v by Mayer's relation.
    pub fn gamma(&self) -> f64 {
        self.specific_heat / (self.specific_heat - self.gas_constant)
    }

    /// The Prandtl number Pr = μ c_p / k, the one number OpenFOAM's `const` transport takes beside μ.
    pub fn prandtl(&self) -> f64 {
        self.viscosity * self.specific_heat / self.conductivity
    }

    /// The speed of sound c = √(γ R_s T) at the case's temperature (m/s).
    pub fn sound_speed(&self) -> f64 {
        (self.gamma() * self.gas_constant * self.temperature).sqrt()
    }

    /// The density of the gas at the ambient pressure and the case's temperature, ρ = p / (R_s T).
    pub fn reference_density(&self) -> f64 {
        self.ambient / (self.gas_constant * self.temperature)
    }
}

/// The heat of formation of the gas, J/kg: no reaction in the case (DOC "Constant (enthalpy)": Hf [J/kg]).
const HEAT_OF_FORMATION: f64 = 0.0;

/// `constant/thermophysicalProperties`: the ideal gas with constant c_p, μ and Pr, an internal-energy
/// (`sensibleInternalEnergy`) form of the `hePsiThermo` package. The keywords are the documented
/// ones of the thermophysical models (DOC "Thermophysical", "Specie", "Perfect gas", "Constant
/// (enthalpy)", "Constant" transport, "Pure Mixture"; UG §5.2 listing). The molar mass is R_u / R_s
/// in g/mol: the units table has no amount of substance, so the gram is converted through it.
pub(super) fn thermophysical(gas: &Gas) -> Result<String> {
    let grams_per_kilogram = Typed::new(1.0, "g").si(Quantity::Mass).map_err(|e| unsupported(e.to_string()))?;
    let molar_weight = MOLAR_GAS_CONSTANT / gas.gas_constant / grams_per_kilogram;
    check_positive(molar_weight, "a molar weight (R_u / R_s) that is not finite and positive: the gas constant is too small")?;
    let body = format!(
        "thermoType\n{{\n    type hePsiThermo;\n    mixture pureMixture;\n    transport const;\n    thermo hConst;\n    \
         equationOfState perfectGas;\n    specie specie;\n    energy sensibleInternalEnergy;\n}}\n\nmixture\n{{\n    \
         specie {{ molWeight {}; }}\n    equationOfState {{ }}\n    thermodynamics {{ Cp {}; Hf {}; }}\n    \
         transport {{ mu {}; Pr {}; }}\n}}\n",
        num(molar_weight),
        num(gas.specific_heat),
        num(HEAT_OF_FORMATION),
        num(gas.viscosity),
        num(gas.prandtl()),
    );
    Ok(dictionary("thermophysicalProperties", &body))
}

/// `system/fvSchemes` for rhoCentralFoam (Kurganov–Noelle–Petrova's central-upwind flux, the
/// van Leer limiter on ρ, U and T, the non-orthogonal correction limited to `psi`).
/// The flux is the central-upwind scheme of KNP (KNP §2–3; G10 §3.1, eqs. (7)–(10)); the
/// reconstruction is the van Leer limiter (vL74; G10 §3.1); `vanLeerV` is the vector variant
/// (G10 eq. (17); UG §6.2.1). Time, gradient and interpolation schemes follow G10 §3 and UG Tables
/// 6.3–6.8. The keywords `fluxScheme` and `reconstruct(...)` are keywords of the format in no
/// documentation page read: they are spelled as the public material names them and are not
/// verified against a case run here.
pub(super) fn fv_schemes(psi: f64) -> String {
    let body = format!(
        "fluxScheme Kurganov;\n\nddtSchemes\n{{\n    default Euler;\n}}\n\ngradSchemes\n{{\n    default Gauss linear;\n}}\n\n\
         divSchemes\n{{\n    default Gauss linear;\n}}\n\nlaplacianSchemes\n{{\n    default Gauss linear limited corrected {psi};\n}}\n\n\
         interpolationSchemes\n{{\n    default linear;\n    reconstruct(rho) vanLeer;\n    reconstruct(U) vanLeerV;\n    \
         reconstruct(T) vanLeer;\n}}\n\nsnGradSchemes\n{{\n    default limited corrected {psi};\n}}\n",
        psi = num(psi),
    );
    dictionary("fvSchemes", &body)
}

/// `system/fvSolution` for rhoCentralFoam: the explicit conserved variables take a diagonal solver
/// (UG Table 6.9; G10 §4: ρ, ρU and ρE are solved explicitly), and the implicit diffusion
/// corrections of U and e a symmetric Gauss–Seidel smoother to `tolerance` (G10 eqs. (19), (21);
/// UG §6.3.2). Field names follow G10's variables; they are not verified against a case run here.
pub(super) fn fv_solution(tolerance: f64) -> String {
    let body = format!(
        "solvers\n{{\n    \"(rho|rhoU|rhoE)\" {{ solver diagonal; }}\n    \"(U|e)\" {{ solver smoothSolver; smoother symGaussSeidel; \
         nSweeps 1; tolerance {}; relTol 0; }}\n}}\n",
        num(tolerance),
    );
    dictionary("fvSolution", &body)
}

/// The field files `0/U`, `0/p` and `0/T`. Pressures are absolute (`[1 -1 -2 0 0 0 0]`); a face
/// with a gauge pressure above zero is a reservoir (`totalPressure`, `totalTemperature`), any other
/// pressure face a non-reflecting outlet (`waveTransmissive`, relaxing to its pressure over
/// `far_field`). The wave condition is the one OpenFOAM's `advective`/`acousticWaveTransmissive`
/// entries describe (DOC v2212; UG §5.1 listing, §A.4 Table A.10): its `lInf` is the far-field
/// length the conditions give. Velocity faces are refused when they are supersonic (|u·n| ≥ c).
pub(super) fn fields(d: &FluidDomain, gas: &Gas, surfaces: &[Surface], notes: &mut Vec<String>) -> Result<[(&'static str, String); 3]> {
    let g = num(gas.gamma());
    let t = num(gas.temperature);
    let c = gas.sound_speed();
    let (mut u, mut p, mut temp) = (String::new(), String::new(), String::new());
    for face in BoxFace::ALL {
        let n = face.name();
        let (fu, fp, ft) = match *d.faces.get(face) {
            FaceCondition::Velocity { velocity, .. } => {
                let [vx, vy, vz] = velocity;
                let normal = match face.axis() {
                    0 => vx,
                    1 => vy,
                    _ => vz,
                };
                if normal.abs() >= c {
                    return Err(unsupported(format!(
                        "{n}: a supersonic velocity face (Mach {:.3}): its pressure would have to be given too (not mapped yet)",
                        normal.abs() / c
                    )));
                }
                (
                    format!("type fixedValue; value uniform {};", vector(velocity)),
                    "type zeroGradient;".to_string(),
                    format!("type inletOutlet; inletValue uniform {t}; value uniform {t};"),
                )
            }
            FaceCondition::Pressure { pressure, .. } => {
                let total = gas.ambient + pressure;
                if !(total.is_finite() && total > 0.0) {
                    return Err(unsupported(format!("{n}: a pressure at or below absolute zero (ambient + gauge)")));
                }
                let big_p = num(total);
                let inlet_outlet_u = "type pressureInletOutletVelocity; tangentialVelocity uniform (0 0 0); value uniform (0 0 0);".to_string();
                if pressure > 0.0 {
                    notes.push(format!("{n}: its pressure is held as the total pressure of a reservoir at the case's temperature"));
                    (
                        inlet_outlet_u,
                        format!("type totalPressure; psi thermo:psi; gamma {g}; p0 uniform {big_p}; value uniform {big_p};"),
                        format!("type totalTemperature; psi thermo:psi; gamma {g}; T0 uniform {t}; value uniform {t};"),
                    )
                } else {
                    notes.push(format!("{n}: a non-reflecting outlet relaxing to its pressure over the far-field length"));
                    (
                        inlet_outlet_u,
                        format!(
                            "type waveTransmissive; field p; psi thermo:psi; gamma {g}; fieldInf {big_p}; lInf {}; value uniform {big_p};",
                            num(gas.far_field)
                        ),
                        format!("type inletOutlet; inletValue uniform {t}; value uniform {t};"),
                    )
                }
            }
            FaceCondition::Slip {} => ("type slip;".into(), "type zeroGradient;".into(), "type zeroGradient;".into()),
            FaceCondition::Symmetry {} => ("type symmetryPlane;".into(), "type symmetryPlane;".into(), "type symmetryPlane;".into()),
            FaceCondition::Wall {} => ("type noSlip;".into(), "type zeroGradient;".into(), "type zeroGradient;".into()),
            FaceCondition::MovingWall { velocity } => {
                (format!("type fixedValue; value uniform {};", vector(velocity)), "type zeroGradient;".into(), "type zeroGradient;".into())
            }
            FaceCondition::Periodic {} => return Err(unsupported(format!("a periodic face ({n}; not mapped yet)"))),
        };
        let _ = writeln!(u, "    {n} {{ {fu} }}");
        let _ = writeln!(p, "    {n} {{ {fp} }}");
        let _ = writeln!(temp, "    {n} {{ {ft} }}");
    }
    for surface in surfaces {
        for r in &surface.regions {
            let fu = match r.condition {
                SurfaceCondition::Wall {} => "type noSlip;".to_string(),
                SurfaceCondition::MovingWall { velocity } => format!("type fixedValue; value uniform {};", vector(velocity)),
            };
            let _ = writeln!(u, "    {} {{ {fu} }}", r.patch);
            let _ = writeln!(p, "    {} {{ type zeroGradient; }}", r.patch);
            let _ = writeln!(temp, "    {} {{ type zeroGradient; }}", r.patch);
        }
    }
    let file = |class: &str, object: &str, dim: &str, internal: String, b: &str| {
        field_file(class, object, &format!("dimensions {dim};\n\ninternalField uniform {internal};\n\nboundaryField\n{{\n{b}}}\n"))
    };
    Ok([
        ("U", file("volVectorField", "U", "[0 1 -1 0 0 0 0]", vector(d.initial_velocity), &u)),
        ("p", file("volScalarField", "p", "[1 -1 -2 0 0 0 0]", num(gas.ambient), &p)),
        ("T", file("volScalarField", "T", "[0 0 0 1 0 0 0]", t.clone(), &temp)),
    ])
}

/// `system/setFieldsDict`: the absolute initial pressure in each box, and the case's temperature
/// everywhere else (TG §4.1.4: defaultFieldValues, regions, boxToCell).
pub(super) fn set_fields(gas: &Gas, boxes: &[PressureBox]) -> Result<String> {
    if boxes.len() > MAX_PRESSURE_BOXES {
        return Err(Error::Limit(format!("{} initial pressure boxes (at most {MAX_PRESSURE_BOXES})", boxes.len())));
    }
    let mut regions = String::new();
    for b in boxes {
        let finite = b.min.iter().chain(b.max.iter()).all(|x| x.is_finite()) && b.pressure.is_finite();
        if !finite {
            return Err(unsupported("an initial pressure box with a number that is not finite"));
        }
        if b.min.iter().zip(b.max.iter()).any(|(lo, hi)| lo >= hi) {
            return Err(unsupported("an initial pressure box whose min is not below its max"));
        }
        let total = gas.ambient + b.pressure;
        if !(total.is_finite() && total > 0.0) {
            return Err(unsupported("an initial pressure box at or below absolute zero"));
        }
        let _ =
            writeln!(regions, "    boxToCell {{ box {} {}; fieldValues ( volScalarFieldValue p {} ); }}", vector(b.min), vector(b.max), num(total));
    }
    let body = format!(
        "defaultFieldValues ( volScalarFieldValue p {} volScalarFieldValue T {} );\n\nregions\n(\n{regions});\n",
        num(gas.ambient),
        num(gas.temperature),
    );
    Ok(dictionary("setFieldsDict", &body))
}

/// The axes whose wave speeds sum into OpenFOAM's Courant number (one per axis of a cell).
const AXES: f64 = 3.0;

/// The speed that sizes the first time step (m/s): the sound speed on every axis, plus the largest
/// imposed speed (the L1 norm of a velocity, as `speed_scale` in mod.rs counts it), and for a
/// reservoir the isentropic speed from its total pressure to the lowest pressure of the domain,
/// u = √(2 c_p T (1 − (p_min/p0)^((γ−1)/γ))) (the energy equation c_p T0 = c_p T + u²/2 with the
/// isentrope). The lowest pressure is the ambient when only inlets exist, or the domain's own
/// lowest outlet pressure when it is below the ambient.
pub(super) fn courant_speed(d: &FluidDomain, gas: &Gas) -> f64 {
    let norm = |v: [f64; 3]| v[0].abs() + v[1].abs() + v[2].abs();
    let mut speed = norm(d.initial_velocity);
    let (mut p_lo, mut gauge_hi) = (gas.ambient, f64::NEG_INFINITY);
    for face in BoxFace::ALL {
        match *d.faces.get(face) {
            FaceCondition::Velocity { velocity, .. } | FaceCondition::MovingWall { velocity } => speed = speed.max(norm(velocity)),
            FaceCondition::Pressure { pressure, .. } => {
                p_lo = p_lo.min(gas.ambient + pressure);
                gauge_hi = gauge_hi.max(pressure);
            }
            _ => {}
        }
    }
    for s in &d.surfaces {
        if let SurfaceCondition::MovingWall { velocity } = s.condition {
            speed = speed.max(norm(velocity));
        }
    }
    if gauge_hi > 0.0 {
        let p0 = gas.ambient + gauge_hi;
        let ratio = libm::pow(p_lo / p0, (gas.gamma() - 1.0) / gas.gamma());
        speed = speed.max((2.0 * gas.specific_heat * gas.temperature * (1.0 - ratio)).sqrt());
    }
    AXES * gas.sound_speed() + speed
}

#[cfg(test)]
mod tests {
    use super::super::{grammar, num};
    use super::*;
    use serde_json::json;

    /// A 10 × 4 × 4 mm channel (SI) with the given material and faces, no surfaces.
    fn channel(material: serde_json::Value, faces: serde_json::Value) -> Model {
        serde_json::from_value(json!({
            "materials": {"Air": material},
            "fluid_domains": {"Flow": {"material": "Air", "min": [0.0, 0.0, 0.0], "max": [0.01, 0.004, 0.004], "spacing": 0.0005, "faces": faces}},
            "steps": [{"name": "Jet", "procedure": {"type": "fluid", "domain": "Flow", "end_time": 0.01, "frames": 5}, "output": []}]
        }))
        .unwrap()
    }

    fn air() -> serde_json::Value {
        json!({"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}, "thermal": {"conductivity": 0.0257, "specific_heat": 1005.0}})
    }

    fn conditions() -> GasConditions {
        GasConditions { ambient: 101_325.0, temperature: 293.15, gas_constant: 287.05, far_field: 0.05, initial_pressure: vec![] }
    }

    fn gas_of(material: serde_json::Value, c: &GasConditions) -> Result<Gas> {
        let model = channel(material, json!({}));
        Gas::of(&model, model.fluid_domains().get("Flow").unwrap(), c)
    }

    fn gas() -> Gas {
        gas_of(air(), &conditions()).unwrap()
    }

    #[test]
    fn a_derived_number_that_overflows_is_refused_by_name() {
        // A subnormal gas constant: every input is finite and positive, the density is not.
        let tiny = GasConditions { gas_constant: 1e-310, ..conditions() };
        let e = gas_of(air(), &tiny).unwrap_err().to_string();
        assert!(e.contains("reference density"), "{e}");
        // The molar weight R_u / R_s overflows on its own, where the density is not asked.
        let g = Gas { gas_constant: 1e-310, ..gas() };
        let e = thermophysical(&g).unwrap_err().to_string();
        assert!(e.contains("molar weight"), "{e}");
        // A conductivity this small overflows the Prandtl number.
        let k = json!({"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}, "thermal": {"conductivity": 1e-310, "specific_heat": 1005.0}});
        let e = gas_of(k, &conditions()).unwrap_err().to_string();
        assert!(e.contains("Prandtl number"), "{e}");
    }

    #[test]
    fn the_gas_comes_from_the_material_and_nothing_is_defaulted() {
        let g = gas();
        assert!((g.gamma() - 1005.0 / (1005.0 - 287.05)).abs() < 1e-12);
        assert!((g.prandtl() - 1.8e-5 * 1005.0 / 0.0257).abs() < 1e-12);
        assert!((g.reference_density() - 101_325.0 / (287.05 * 293.15)).abs() < 1e-12);

        let no_thermal = json!({"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}});
        let no_cp = json!({"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}, "thermal": {"conductivity": 0.0257}});
        let no_k = json!({"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}, "thermal": {"specific_heat": 1005.0}});
        let k_axes = json!({"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}, "thermal": {"conductivity": [0.02, 0.02, 0.03], "specific_heat": 1005.0}});
        let refused = [
            (gas_of(no_thermal, &conditions()), "specific heat"),
            (gas_of(no_cp, &conditions()), "specific heat"),
            (gas_of(no_k, &conditions()), "conductivity"),
            (gas_of(k_axes, &conditions()), "per axis"),
        ];
        for (result, words) in refused {
            match result {
                Err(e) => assert!(e.to_string().contains(words), "{e} does not name {words}"),
                Ok(g) => panic!("refused input gave {g:?}"),
            }
        }

        let with = |edit: fn(&mut GasConditions)| {
            let mut c = conditions();
            edit(&mut c);
            c
        };
        let bad = [
            (with(|c| c.gas_constant = 1005.0), "cv"),
            (with(|c| c.gas_constant = 2000.0), "cv"),
            (with(|c| c.ambient = 0.0), "ambient"),
            (with(|c| c.temperature = -1.0), "temperature"),
            (with(|c| c.gas_constant = f64::NAN), "gas constant"),
            (with(|c| c.far_field = 0.0), "far-field"),
            (with(|c| c.far_field = f64::NAN), "far-field"),
        ];
        for (c, words) in bad {
            match gas_of(air(), &c) {
                Err(e) => assert!(e.to_string().contains(words), "{e} does not name {words}"),
                Ok(g) => panic!("refused conditions gave {g:?}"),
            }
        }

        let argon =
            json!({"density": 1.6, "fluid": {"type": "newtonian", "viscosity": 2.2e-5}, "thermal": {"conductivity": 0.0177, "specific_heat": 520.3}});
        let argon_gas = gas_of(argon, &GasConditions { gas_constant: 208.13, ..conditions() });
        assert!(argon_gas.is_ok(), "argon is refused: {argon_gas:?}");
    }

    #[test]
    fn thermophysical_properties_are_the_documented_package() {
        let text = thermophysical(&gas()).unwrap();
        let dict = grammar::check_file(&text).unwrap();
        assert_eq!(dict.text("thermoType/type").as_deref(), Some("hePsiThermo"));
        assert_eq!(dict.text("thermoType/energy").as_deref(), Some("sensibleInternalEnergy"));
        assert_eq!(dict.text("thermoType/equationOfState").as_deref(), Some("perfectGas"));
        assert!((dict.number("mixture/specie/molWeight").unwrap() - 28.9652).abs() < 1e-4);
        assert_eq!(dict.number("mixture/thermodynamics/Cp"), Some(1005.0));
        assert_eq!(dict.number("mixture/thermodynamics/Hf"), Some(0.0));
        assert_eq!(dict.number("mixture/transport/mu"), Some(1.8e-5));
        assert!((dict.number("mixture/transport/Pr").unwrap() - gas().prandtl()).abs() < 1e-12);
        assert!(matches!(dict.get("mixture/equationOfState"), Some(grammar::Entry::Dict(_))));
    }

    #[test]
    fn schemes_and_solvers_are_the_cited_ones() {
        let psi = 0.333;
        let text = fv_schemes(psi);
        let dict = grammar::check_file(&text).unwrap();
        assert_eq!(dict.text("fluxScheme").as_deref(), Some("Kurganov"));
        assert_eq!(dict.text("ddtSchemes/default").as_deref(), Some("Euler"));
        assert_eq!(dict.text("divSchemes/default").as_deref(), Some("Gauss linear"));
        assert!(text.contains(&format!("default Gauss linear limited corrected {};", num(psi))));
        assert_eq!(dict.text("interpolationSchemes/reconstruct(U)").as_deref(), Some("vanLeerV"));
        assert_eq!(dict.text("interpolationSchemes/reconstruct(rho)").as_deref(), Some("vanLeer"));
        assert_eq!(dict.text("interpolationSchemes/reconstruct(T)").as_deref(), Some("vanLeer"));
        assert!(!text.contains("tauMC"));

        let text = fv_solution(1e-7);
        grammar::check_file(&text).unwrap();
        assert!(text.contains("\"(rho|rhoU|rhoE)\" { solver diagonal; }"));
        assert!(text.contains("\"(U|e)\""));
        assert!(text.contains("relTol 0;"));
        assert!(!text.contains("1e-9"));
    }

    /// The field files of `faces` on the channel.
    fn written(faces: serde_json::Value, c: &GasConditions) -> ([(&'static str, String); 3], Vec<String>) {
        let model = channel(air(), faces);
        let mut notes = Vec::new();
        let g = Gas::of(&model, model.fluid_domains().get("Flow").unwrap(), c).unwrap();
        let files = fields(model.fluid_domains().get("Flow").unwrap(), &g, &[], &mut notes).unwrap();
        (files, notes)
    }

    #[test]
    fn faces_become_reservoirs_non_reflecting_outlets_and_adiabatic_walls() {
        let faces = json!({"x_min": {"type": "pressure", "pressure": 2.0e5}, "x_max": {"type": "pressure", "pressure": 0.0}});
        let c = GasConditions { far_field: 0.05, ..conditions() };
        let ([(_, u), (_, p), (_, t)], notes) = written(faces.clone(), &c);
        assert!(p.contains("x_min { type totalPressure; psi thermo:psi; gamma"));
        assert!(p.contains("p0 uniform 3.01325e5;"));
        assert!(p.contains("x_max { type waveTransmissive;"));
        assert!(p.contains("fieldInf 1.01325e5;"));
        assert!(p.contains("lInf 5e-2;"));
        assert!(p.contains("internalField uniform 1.01325e5"));
        assert!(t.contains("x_min { type totalTemperature;"));
        assert!(t.contains("y_min { type zeroGradient; }"));
        assert!(t.contains("x_max { type inletOutlet; inletValue uniform 2.9315e2;"));
        assert!(u.contains("x_max { type pressureInletOutletVelocity; tangentialVelocity uniform (0 0 0);"));
        assert_eq!(notes.len(), 2);
        assert!(notes.iter().any(|n| n.contains("x_min") && n.contains("total pressure")));
        assert!(notes.iter().any(|n| n.contains("x_max")));
        for file in [&u, &p, &t] {
            grammar::check_file(file).unwrap();
        }

        let characteristic = json!({"x_min": {"type": "pressure", "pressure": 2.0e5}, "x_max": {"type": "pressure", "pressure": 0.0, "outflow": "characteristic", "sigma": 0.5}});
        let ([_, (_, p), _], _) = written(characteristic, &c);
        assert!(p.contains("lInf 5e-2;"));

        let supersonic = json!({"x_min": {"type": "velocity", "velocity": [400.0, 0.0, 0.0]}, "x_max": {"type": "pressure", "pressure": 0.0}});
        let model = channel(air(), supersonic);
        let mut notes = Vec::new();
        let g = Gas::of(&model, model.fluid_domains().get("Flow").unwrap(), &c).unwrap();
        match fields(model.fluid_domains().get("Flow").unwrap(), &g, &[], &mut notes) {
            Err(e) => assert!(e.to_string().contains("supersonic"), "{e}"),
            Ok(_) => panic!("a supersonic inlet was written"),
        }
        let subsonic = json!({"x_min": {"type": "velocity", "velocity": [300.0, 0.0, 0.0]}, "x_max": {"type": "pressure", "pressure": 0.0}});
        let ([_, _, (_, t)], _) = written(subsonic, &c);
        assert!(t.contains("x_min { type inletOutlet;"));

        let below_zero = json!({"x_max": {"type": "pressure", "pressure": -2.0e5}});
        let model = channel(air(), below_zero);
        let mut notes = Vec::new();
        let g = Gas::of(&model, model.fluid_domains().get("Flow").unwrap(), &c).unwrap();
        match fields(model.fluid_domains().get("Flow").unwrap(), &g, &[], &mut notes) {
            Err(e) => assert!(e.to_string().contains("absolute zero"), "{e}"),
            Ok(_) => panic!("a pressure below absolute zero was written"),
        }
    }

    #[test]
    fn initial_pressure_boxes_are_set_fields_regions() {
        let g = gas();
        let boxes = [PressureBox { min: [-1.0; 3], max: [0.003, 1.0, 1.0], pressure: 2.99e7 }];
        let text = set_fields(&g, &boxes).unwrap();
        grammar::check_file(&text).unwrap();
        assert!(text.contains("box (-1e0 -1e0 -1e0) (3e-3 1e0 1e0); fieldValues ( volScalarFieldValue p 3.0001325e7 );"));
        assert!(text.contains("volScalarFieldValue T 2.9315e2"));

        let degenerate = [PressureBox { min: [0.0; 3], max: [0.0, 1.0, 1.0], pressure: 1.0 }];
        assert!(set_fields(&g, &degenerate).is_err());
        let below_zero = [PressureBox { min: [0.0; 3], max: [1.0; 3], pressure: -2.0e5 }];
        assert!(set_fields(&g, &below_zero).is_err());
        let not_a_number = [PressureBox { min: [0.0; 3], max: [1.0; 3], pressure: f64::NAN }];
        assert!(set_fields(&g, &not_a_number).is_err());
        let too_many = vec![PressureBox { min: [0.0; 3], max: [1.0; 3], pressure: 1.0 }; MAX_PRESSURE_BOXES + 1];
        assert!(matches!(set_fields(&g, &too_many), Err(Error::Limit(_))));
    }

    #[test]
    fn the_first_step_counts_the_sound_speed_on_every_axis() {
        let g = gas();
        let c = g.sound_speed();
        assert!((c - (g.gamma() * 287.05 * 293.15).sqrt()).abs() < 1e-9);

        let model = channel(air(), json!({}));
        let d = model.fluid_domains().get("Flow").unwrap();
        assert!((courant_speed(d, &g) - 3.0 * c).abs() < 1e-9);

        let moving = channel(air(), json!({"x_min": {"type": "velocity", "velocity": [2.0, 1.0, 0.0]}}));
        assert!((courant_speed(moving.fluid_domains().get("Flow").unwrap(), &g) - (3.0 * c + 3.0)).abs() < 1e-9);

        let driven = channel(air(), json!({"x_min": {"type": "pressure", "pressure": 2.0e5}, "x_max": {"type": "pressure", "pressure": 0.0}}));
        let expected = 3.0 * c + (2.0 * 1005.0 * 293.15 * (1.0 - (101_325.0_f64 / 301_325.0).powf((g.gamma() - 1.0) / g.gamma()))).sqrt();
        assert!((courant_speed(driven.fluid_domains().get("Flow").unwrap(), &g) - expected).abs() < 1e-9);
    }
}
