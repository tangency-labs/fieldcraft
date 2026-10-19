//! Materials: linear elasticity, mass density, fluid behaviour, thermal properties (conductivity,
//! specific heat, expansion), and the inelastic behaviours: von Mises (J2) plasticity with a
//! bilinear, tabulated or Johnson–Cook flow stress (explicit dynamics and nonlinear static
//! steps), Johnson–Cook damage and element erosion criteria (explicit dynamics).
//!
//! The records are solver-neutral and named for what they describe, never after a product's
//! keyword or material number; the deck writers map them to each format.
//!
//! - Johnson–Cook flow stress: G. R. Johnson, W. H. Cook, *A constitutive model and data for
//!   metals subjected to large strains, high strain rates and high temperatures*, Proc. 7th Int.
//!   Symp. on Ballistics, The Hague, 541–547, 1983:
//!   σ_y = (A + B ε_pⁿ)(1 + C ln ε̇*)(1 − T*ᵐ), ε̇* = ε̇_p / ε̇₀, T* = (T − T_room)/(T_melt − T_room).
//! - Johnson–Cook fracture: G. R. Johnson, W. H. Cook, *Fracture characteristics of three metals
//!   subjected to various strains, strain rates, temperatures and pressures*, Eng. Fract. Mech.
//!   21(1):31–48, 1985: ε_f = (D1 + D2 e^{D3 σ*})(1 + D4 ln ε̇*)(1 + D5 T*), D = Σ Δε_p / ε_f.

use serde::{Deserialize, Serialize};

use crate::{Error, Model, Result, context, finite, insert_new, not_in_use, positive, remove_existing, replace_existing};

/// Names of the engineering constants' shear moduli, in order.
pub(crate) const SHEAR_MODULI: [&str; 3] = ["G12", "G13", "G23"];

/// Linear elastic behaviour.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Elastic {
    /// Young's modulus and Poisson's ratio.
    Isotropic { young: f64, poisson: f64 },
    /// Orthotropic, by engineering constants in the material axes 1, 2, 3 (the global X, Y, Z
    /// until material orientations land): Young's moduli `[E1, E2, E3]`, Poisson's ratios
    /// `[ν12, ν13, ν23]` (ν_ij is the contraction in j under stress in i) and shear moduli
    /// `[G12, G13, G23]`.
    EngineeringConstants { young: [f64; 3], poisson: [f64; 3], shear: [f64; 3] },
}

impl Elastic {
    /// Check the constants describe a stable material (positive-definite stiffness).
    pub fn validate(&self) -> Result<()> {
        match *self {
            Elastic::Isotropic { young, poisson } => {
                positive("Young's modulus", young)?;
                if !(poisson > -1.0 && poisson < 0.5) {
                    return Err(Error::Invalid(format!("Poisson's ratio must be above -1 and below 0.5, got {poisson}")));
                }
            }
            Elastic::EngineeringConstants { young, poisson, shear } => {
                for (i, e) in young.iter().enumerate() {
                    positive(&format!("E{}", i + 1), *e)?;
                }
                for (what, g) in SHEAR_MODULI.iter().zip(shear) {
                    positive(what, g)?;
                }
                let [e1, e2, e3] = young;
                let [n12, n13, n23] = poisson;
                if poisson.iter().any(|n| !n.is_finite()) {
                    return Err(Error::Invalid("Poisson's ratios must be finite numbers".into()));
                }
                // Reciprocal ratios ν_ji = ν_ij E_j / E_i; stability needs |ν_ij| < sqrt(E_i/E_j)
                // and a positive determinant of the compliance's normal block.
                let (n21, n31, n32) = (n12 * e2 / e1, n13 * e3 / e1, n23 * e3 / e2);
                for (name, nij, ei, ej) in [("ν12", n12, e1, e2), ("ν13", n13, e1, e3), ("ν23", n23, e2, e3)] {
                    if !(nij.abs() < (ei / ej).sqrt()) {
                        return Err(Error::Invalid(format!("{name} = {nij} is too large for the moduli (must be below {})", (ei / ej).sqrt())));
                    }
                }
                let det = 1.0 - n12 * n21 - n23 * n32 - n31 * n13 - 2.0 * n21 * n32 * n13;
                if !(det > 0.0) {
                    return Err(Error::Invalid("the Poisson's ratios make the material unstable (stiffness not positive definite)".into()));
                }
            }
        }
        Ok(())
    }
}

/// Hyperelastic behaviour of an isotropic material at large strain, in nonlinear steps.
///
/// The constants are the material's isotropic elastic ones: the shear modulus μ = E / (2(1 + ν))
/// and the bulk modulus κ = E / (3(1 − 2ν)), so at small strain the material is that linear
/// elastic one (linear steps use it so). Without this option an elastic material in a nonlinear
/// step is St Venant–Kirchhoff (the second Piola–Kirchhoff stress linear in the Green–Lagrange
/// strain), which is right for large rotations with small strains.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Hyperelastic {
    /// Compressible neo-Hookean with the isochoric–volumetric split: the strain energy per
    /// reference volume W = μ/2 (Ī₁ − 3) + κ/2 (J − 1)², with Ī₁ = J^(−2/3) tr C (C₁₀ = μ/2,
    /// D₁ = 2/κ in the usual polynomial notation).
    NeoHookean {},
    /// Mooney–Rivlin (M. Mooney, J. Appl. Phys. 11:582–592, 1940; R. S. Rivlin, Phil. Trans. R.
    /// Soc. A 241:379–397, 1948), decoupled: W = C₁₀ (Ī₁ − 3) + C₀₁ (Ī₂ − 3) + κ/2 (J − 1)². The
    /// small-strain shear modulus 2 (C₁₀ + C₀₁) is the elastic μ, so `c01` sets the second
    /// invariant's share and C₁₀ = μ/2 − C₀₁ (0 ≤ C₀₁ ≤ μ/2, both constants non-negative).
    MooneyRivlin { c01: f64 },
    /// Yeoh (O. H. Yeoh, Rubber Chem. Technol. 66(5):754–771, 1993), decoupled:
    /// W = C₁₀ (Ī₁ − 3) + C₂₀ (Ī₁ − 3)² + C₃₀ (Ī₁ − 3)³ + κ/2 (J − 1)², with C₁₀ = μ/2 from the
    /// elastic constants (the small-strain shear modulus); `c20` (often negative) and `c30`
    /// shape the large-strain response.
    Yeoh { c20: f64, c30: f64 },
}

impl Hyperelastic {
    /// The shear and bulk moduli μ = E / (2(1 + ν)) and κ = E / (3(1 − 2ν)) the model takes from
    /// the isotropic constants.
    pub fn moduli(young: f64, poisson: f64) -> (f64, f64) {
        (young / (2.0 * (1.0 + poisson)), young / (3.0 * (1.0 - 2.0 * poisson)))
    }

    /// The name of the model, for messages.
    pub fn name(&self) -> &'static str {
        match self {
            Hyperelastic::NeoHookean {} => "neo-Hookean",
            Hyperelastic::MooneyRivlin { .. } => "Mooney–Rivlin",
            Hyperelastic::Yeoh { .. } => "Yeoh",
        }
    }

    /// The polynomial constants (C₁₀, C₀₁, C₂₀, C₃₀) for the shear modulus μ.
    pub fn constants(&self, shear: f64) -> [f64; 4] {
        match *self {
            Hyperelastic::NeoHookean {} => [shear / 2.0, 0.0, 0.0, 0.0],
            Hyperelastic::MooneyRivlin { c01 } => [shear / 2.0 - c01, c01, 0.0, 0.0],
            Hyperelastic::Yeoh { c20, c30 } => [shear / 2.0, 0.0, c20, c30],
        }
    }

    /// Check the option's constants against the elastic shear modulus μ.
    pub fn validate(&self, shear: f64) -> Result<()> {
        match *self {
            Hyperelastic::NeoHookean {} => {}
            Hyperelastic::MooneyRivlin { c01 } => {
                if !(c01 >= 0.0 && c01 <= shear / 2.0) {
                    return Err(Error::Invalid(format!(
                        "the Mooney–Rivlin C01 = {c01} must be within [0, μ/2 = {}] (C10 = μ/2 − C01 and C01 both non-negative)",
                        shear / 2.0
                    )));
                }
            }
            Hyperelastic::Yeoh { c20, c30 } => {
                finite("Yeoh C20", c20)?;
                finite("Yeoh C30", c30)?;
            }
        }
        Ok(())
    }
}

/// Fluid behaviour.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Fluid {
    /// A Newtonian fluid of this dynamic viscosity μ (shear stress per shear rate; for example
    /// 1.0e-9 MPa·s for water in mm, t, s). Its kinematic viscosity is μ over the density.
    Newtonian { viscosity: f64 },
}

impl Fluid {
    /// The dynamic viscosity μ.
    pub fn viscosity(&self) -> f64 {
        match *self {
            Fluid::Newtonian { viscosity } => viscosity,
        }
    }
}

/// Most points of a plastic hardening table.
pub const MAX_TABLE_POINTS: usize = 1000;

/// How the yield surface moves with plastic flow.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Hardening {
    /// The yield surface grows (the yield stress rises with the equivalent plastic strain).
    #[default]
    Isotropic,
    /// The yield surface translates (linear kinematic hardening of the back stress, Prager), so
    /// reversed loading yields early (Bauschinger effect).
    Kinematic,
}

/// One point of a hardening table: the yield stress at an equivalent plastic strain.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HardeningPoint {
    /// Equivalent plastic strain.
    pub strain: f64,
    /// Yield stress (von Mises equivalent).
    pub stress: f64,
}

/// Rate-independent or Johnson–Cook von Mises plasticity (needs isotropic elasticity).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Plastic {
    /// Linear hardening: the yield stress `yield_stress + hardening_modulus · ε_p` (isotropic), or
    /// a yield surface of radius `yield_stress` whose centre moves with modulus
    /// `hardening_modulus` (kinematic). The hardening modulus is the slope of the uniaxial stress
    /// against the plastic strain (0: perfectly plastic).
    Bilinear {
        yield_stress: f64,
        #[serde(default)]
        hardening_modulus: f64,
        #[serde(default)]
        hardening: Hardening,
    },
    /// Isotropic hardening by a table of yield stress against equivalent plastic strain, linear
    /// between the points and constant after the last. The first point is at plastic strain 0;
    /// the strains increase.
    Table { points: Vec<HardeningPoint> },
    /// Johnson–Cook flow stress (isotropic hardening): `a`, `b`, `n` the strain hardening, `c`
    /// the strain-rate sensitivity relative to `reference_strain_rate` (rates below it count as
    /// it), `m` the thermal softening between `room_temperature` and `melt_temperature`. With a
    /// `specific_heat`, plastic work heats the material adiabatically (the fraction
    /// `taylor_quinney` of it); without one the temperature stays at room temperature.
    JohnsonCook {
        a: f64,
        b: f64,
        n: f64,
        #[serde(default)]
        c: f64,
        #[serde(default = "one")]
        m: f64,
        #[serde(default = "one")]
        reference_strain_rate: f64,
        room_temperature: f64,
        melt_temperature: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        specific_heat: Option<f64>,
        #[serde(default = "default_taylor_quinney")]
        taylor_quinney: f64,
    },
}

fn one() -> f64 {
    1.0
}

fn default_taylor_quinney() -> f64 {
    0.9
}

impl Plastic {
    /// The name of the kind, for messages.
    pub fn name(&self) -> &'static str {
        match self {
            Plastic::Bilinear { .. } => "bilinear plasticity",
            Plastic::Table { .. } => "tabulated plasticity",
            Plastic::JohnsonCook { .. } => "Johnson–Cook plasticity",
        }
    }

    /// The initial yield stress (at zero plastic strain, the reference rate and room
    /// temperature).
    pub fn initial_yield(&self) -> f64 {
        match self {
            Plastic::Bilinear { yield_stress, .. } => *yield_stress,
            Plastic::Table { points } => points.first().map_or(0.0, |p| p.stress),
            Plastic::JohnsonCook { a, .. } => *a,
        }
    }

    pub fn validate(&self) -> Result<()> {
        let nonnegative = |what: &str, v: f64| {
            if v >= 0.0 && v.is_finite() { Ok(v) } else { Err(Error::Invalid(format!("{what} must be a finite number ≥ 0, got {v}"))) }
        };
        match self {
            Plastic::Bilinear { yield_stress, hardening_modulus, .. } => {
                positive("yield stress", *yield_stress)?;
                nonnegative("hardening modulus", *hardening_modulus)?;
            }
            Plastic::Table { points } => {
                if points.is_empty() || points.len() > MAX_TABLE_POINTS {
                    return Err(Error::Invalid(format!("a hardening table has 1 to {MAX_TABLE_POINTS} points, got {}", points.len())));
                }
                if points.first().is_some_and(|p| p.strain != 0.0) {
                    return Err(Error::Invalid("a hardening table starts at plastic strain 0".into()));
                }
                for p in points {
                    nonnegative("plastic strain", p.strain)?;
                    positive("yield stress", p.stress)?;
                }
                if points.windows(2).any(|w| w.get(1).map(|b| b.strain) <= w.first().map(|a| a.strain)) {
                    return Err(Error::Invalid("the plastic strains of a hardening table must increase".into()));
                }
            }
            Plastic::JohnsonCook { a, b, n, c, m, reference_strain_rate, room_temperature, melt_temperature, specific_heat, taylor_quinney } => {
                positive("A", *a)?;
                nonnegative("B", *b)?;
                nonnegative("n", *n)?;
                nonnegative("C", *c)?;
                positive("m", *m)?;
                positive("reference strain rate", *reference_strain_rate)?;
                finite("room temperature", *room_temperature)?;
                finite("melt temperature", *melt_temperature)?;
                if !(melt_temperature > room_temperature) {
                    return Err(Error::Invalid(format!(
                        "the melt temperature {melt_temperature} must be above the room temperature {room_temperature}"
                    )));
                }
                if let Some(cp) = specific_heat {
                    positive("specific heat", *cp)?;
                }
                if !(0.0..=1.0).contains(taylor_quinney) {
                    return Err(Error::Invalid(format!("the Taylor–Quinney fraction is in [0, 1], got {taylor_quinney}")));
                }
            }
        }
        Ok(())
    }
}

/// A property given per material axis: one value for every direction, or three along the
/// material axes 1, 2, 3 (the global X, Y, Z until material orientations land).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Directional {
    Isotropic(f64),
    Orthotropic([f64; 3]),
}

impl Directional {
    /// The values along the three axes.
    pub fn axes(self) -> [f64; 3] {
        match self {
            Directional::Isotropic(v) => [v; 3],
            Directional::Orthotropic(v) => v,
        }
    }
}

/// Heat conduction and thermal expansion.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Thermal {
    /// Conductivity k (power per length per temperature): q = −k ∇T.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<Directional>,
    /// Specific heat c (energy per mass per temperature); with the density, the heat capacity
    /// of transient heat transfer.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<f64>,
    /// Linear expansion coefficient α (strain per temperature): ε_th = α (T − T_ref).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expansion: Option<Directional>,
    /// The temperature at which thermal strain is zero (kelvin; 0 °C when left out).
    #[serde(default = "celsius_zero")]
    pub reference_temperature: f64,
}

fn celsius_zero() -> f64 {
    crate::units::CELSIUS_ZERO
}

impl Default for Thermal {
    fn default() -> Self {
        Thermal { conductivity: None, specific_heat: None, expansion: None, reference_temperature: celsius_zero() }
    }
}

impl Thermal {
    pub fn validate(&self) -> Result<()> {
        if let Some(k) = self.conductivity {
            for v in k.axes() {
                positive("conductivity", v)?;
            }
        }
        if let Some(c) = self.specific_heat {
            positive("specific heat", c)?;
        }
        if let Some(a) = self.expansion {
            for v in a.axes() {
                finite("expansion coefficient", v)?;
            }
        }
        finite("reference temperature", self.reference_temperature)?;
        if self.conductivity.is_none() && self.specific_heat.is_none() && self.expansion.is_none() {
            return Err(Error::Invalid("thermal properties need a conductivity, a specific heat or an expansion coefficient".into()));
        }
        Ok(())
    }
}

/// Damage that accumulates with plastic strain; a material point fails when it reaches 1.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Damage {
    /// Johnson–Cook fracture strain ε_f = (d1 + d2 e^{d3 σ*})(1 + d4 ln ε̇*)(1 + d5 T*), with σ*
    /// the mean stress over the von Mises stress; ε̇* and T* as in Johnson–Cook plasticity (`d4`
    /// and `d5` need it).
    JohnsonCook {
        d1: f64,
        d2: f64,
        d3: f64,
        #[serde(default)]
        d4: f64,
        #[serde(default)]
        d5: f64,
    },
}

/// When elements of the material are removed from the analysis (eroded), besides failing by
/// [`Damage`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Erosion {
    /// The equivalent plastic strain at which an element erodes.
    pub plastic_strain: f64,
}

/// A material.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elastic: Option<Elastic>,
    /// Large-strain behaviour in nonlinear steps, built on the isotropic elastic constants.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hyperelastic: Option<Hyperelastic>,
    /// Mass per volume, needed by frequency steps, gravity and transient heat transfer.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density: Option<f64>,
    /// How it flows, for fluid domains (needs a density).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fluid: Option<Fluid>,
    /// Von Mises plasticity (explicit dynamics and nonlinear static steps; the linear solvers use
    /// the elastic constants alone).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plastic: Option<Plastic>,
    /// Damage to failure (explicit dynamics: a failed element erodes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub damage: Option<Damage>,
    /// Erosion criteria (explicit dynamics).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub erosion: Option<Erosion>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thermal: Option<Thermal>,
}

impl Material {
    pub fn validate(&self) -> Result<()> {
        if self.elastic.is_none() && self.density.is_none() && self.fluid.is_none() && self.thermal.is_none() {
            return Err(Error::Invalid("a material needs elastic constants, a density, fluid behaviour or thermal properties".into()));
        }
        if let Some(f) = &self.fluid {
            positive("viscosity", f.viscosity())?;
            if self.density.is_none() {
                return Err(Error::Invalid("a fluid needs a density".into()));
            }
        }
        if let Some(e) = &self.elastic {
            e.validate()?;
        }
        if let Some(h) = &self.hyperelastic {
            if !matches!(self.elastic, Some(Elastic::Isotropic { .. })) {
                return Err(Error::Invalid(format!("{} hyperelasticity needs isotropic elastic constants", h.name())));
            }
            if self.plastic.is_some() {
                return Err(Error::Invalid(format!("plasticity of a {} hyperelastic material is not supported yet", h.name())));
            }
            if let Some(Elastic::Isotropic { young, poisson }) = self.elastic {
                h.validate(Hyperelastic::moduli(young, poisson).0)?;
            }
        }
        if let Some(d) = self.density {
            positive("density", d)?;
        }
        if let Some(p) = &self.plastic {
            if !matches!(self.elastic, Some(Elastic::Isotropic { .. })) {
                return Err(Error::Invalid(format!("{} needs isotropic elastic constants", p.name())));
            }
            p.validate()?;
        }
        if let Some(Damage::JohnsonCook { d1, d2, d3, d4, d5 }) = self.damage {
            for (what, v) in [("D1", d1), ("D2", d2), ("D3", d3), ("D4", d4), ("D5", d5)] {
                finite(what, v)?;
            }
            match &self.plastic {
                None => return Err(Error::Invalid("damage needs plasticity (it accumulates with the plastic strain)".into())),
                Some(Plastic::JohnsonCook { .. }) => {}
                Some(_) if d4 != 0.0 || d5 != 0.0 => {
                    return Err(Error::Invalid(
                        "the damage terms D4 and D5 need Johnson–Cook plasticity (its reference strain rate and temperatures)".into(),
                    ));
                }
                Some(_) => {}
            }
        }
        if let Some(e) = self.erosion {
            positive("erosion plastic strain", e.plastic_strain)?;
            if self.plastic.is_none() {
                return Err(Error::Invalid("a plastic-strain erosion criterion needs plasticity".into()));
            }
        }
        if let Some(t) = &self.thermal {
            t.validate()?;
        }
        Ok(())
    }
}

impl Model {
    pub fn materials(&self) -> &std::collections::BTreeMap<String, Material> {
        &self.materials
    }

    pub fn material(&self, name: &str) -> Option<&Material> {
        self.materials.get(name)
    }

    /// The hyperelastic materials some section uses, by name.
    pub fn hyperelastic_in_use(&self) -> std::collections::BTreeSet<&str> {
        self.sections.values().map(|s| s.material.as_str()).filter(|m| self.material(m).is_some_and(|m| m.hyperelastic.is_some())).collect()
    }

    pub fn add_material(&mut self, name: &str, material: Material) -> Result<()> {
        material.validate().map_err(|e| context("material", name, e))?;
        insert_new(&mut self.materials, "material", name, material)
    }

    pub fn replace_material(&mut self, name: &str, material: Material) -> Result<()> {
        material.validate().map_err(|e| context("material", name, e))?;
        self.check_material_fits_fluid_domains(name, &material)?;
        replace_existing(&mut self.materials, "material", name, material)
    }

    /// The definitions that refer to a material, as `kind \`name\`` strings.
    pub fn material_users(&self, name: &str) -> Vec<String> {
        let mut users: Vec<String> = self.sections.iter().filter(|(_, s)| s.material == name).map(|(n, _)| format!("section `{n}`")).collect();
        users.extend(self.fluid_domains_of_material(name).map(|n| format!("fluid domain `{n}`")));
        users
    }

    pub fn remove_material(&mut self, name: &str) -> Result<()> {
        not_in_use("material", name, self.material_users(name))?;
        remove_existing(&mut self.materials, "material", name).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iso(young: f64, poisson: f64) -> Material {
        Material { elastic: Some(Elastic::Isotropic { young, poisson }), ..Material::default() }
    }

    #[test]
    fn materials_validate_and_name_the_bad_constant() {
        let mut m = Model::default();
        m.add_material("Steel", Material { density: Some(7.85e-9), ..iso(2.1e5, 0.3) }).unwrap();
        assert!(m.add_material("Steel", iso(1.0, 0.3)).is_err(), "duplicate");
        for bad in [iso(0.0, 0.3), iso(-1.0, 0.3), iso(f64::NAN, 0.3), iso(1.0, 0.5), iso(1.0, -1.0), iso(1.0, f64::INFINITY), Material::default()] {
            let e = m.add_material("Bad", bad).unwrap_err();
            assert!(e.to_string().contains("material `Bad`"), "{e}");
        }
        assert!(m.add_material("Light", Material { density: Some(-1.0), ..Material::default() }).is_err());
        let fluid = |viscosity, density| Material { density, fluid: Some(Fluid::Newtonian { viscosity }), ..Material::default() };
        m.add_material("Water", fluid(1.0e-9, Some(1.0e-9))).unwrap();
        for bad in [fluid(0.0, Some(1.0)), fluid(-1.0, Some(1.0)), fluid(f64::NAN, Some(1.0)), fluid(1.0, None)] {
            assert!(m.add_material("Bad", bad).is_err());
        }
        assert!(m.add_material("1st", iso(1.0, 0.3)).is_err(), "names start with a letter");
        assert!(m.add_material("has space", iso(1.0, 0.3)).is_err());
        m.replace_material("Steel", iso(2.0e5, 0.29)).unwrap();
        assert!(m.replace_material("Nope", iso(2.0e5, 0.29)).is_err());
        assert_eq!(m.materials().len(), 2);
    }

    #[test]
    fn orthotropic_stability() {
        let ortho = |poisson| Elastic::EngineeringConstants { young: [140e3, 10e3, 10e3], poisson, shear: [5e3, 5e3, 3.8e3] };
        ortho([0.3, 0.3, 0.4]).validate().unwrap();
        // ν21 = ν12 E2/E1 is small; ν12 itself may exceed 0.5 in a stiff direction.
        ortho([0.9, 0.3, 0.4]).validate().unwrap();
        assert!(ortho([4.0, 0.3, 0.4]).validate().is_err(), "|ν12| < sqrt(E1/E2) = 3.74");
        assert!(ortho([0.3, 0.3, 1.0]).validate().is_err());
        assert!(ortho([3.5, 3.5, 0.9]).validate().is_err(), "determinant");
        assert!(ortho([0.3, f64::NAN, 0.4]).validate().is_err());
        let zero_shear = Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0, 0.0, 1.0] };
        assert!(zero_shear.validate().unwrap_err().to_string().contains("G13"));
    }

    fn steel() -> Material {
        Material { density: Some(7.85e-9), ..iso(2.1e5, 0.3) }
    }

    fn jc() -> Plastic {
        Plastic::JohnsonCook {
            a: 792.0,
            b: 510.0,
            n: 0.26,
            c: 0.014,
            m: 1.03,
            reference_strain_rate: 1.0,
            room_temperature: 294.0,
            melt_temperature: 1793.0,
            specific_heat: Some(4.77e8),
            taylor_quinney: 0.9,
        }
    }

    #[test]
    fn plasticity_damage_and_erosion_validate() {
        let ok = |p: Plastic| Material { plastic: Some(p), ..steel() }.validate();
        ok(Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 1000.0, hardening: Hardening::Kinematic }).unwrap();
        ok(Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: 250.0 }, HardeningPoint { strain: 0.1, stress: 300.0 }] }).unwrap();
        ok(jc()).unwrap();
        for bad in [
            Plastic::Bilinear { yield_stress: 0.0, hardening_modulus: 0.0, hardening: Hardening::Isotropic },
            Plastic::Bilinear { yield_stress: 1.0, hardening_modulus: -1.0, hardening: Hardening::Isotropic },
            Plastic::Table { points: vec![] },
            Plastic::Table { points: vec![HardeningPoint { strain: 0.1, stress: 250.0 }] },
            Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: 250.0 }, HardeningPoint { strain: 0.0, stress: 300.0 }] },
            Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: f64::NAN }] },
            Plastic::Table { points: vec![HardeningPoint { strain: 0.0, stress: 1.0 }; MAX_TABLE_POINTS + 1] },
        ] {
            assert!(ok(bad.clone()).is_err(), "{bad:?}");
        }
        let mut cold = jc();
        if let Plastic::JohnsonCook { melt_temperature, room_temperature, .. } = &mut cold {
            *melt_temperature = *room_temperature;
        }
        assert!(ok(cold).is_err(), "melt above room");
        // Plasticity needs isotropic elasticity; damage and erosion need plasticity.
        let ortho = Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0; 3] };
        assert!(Material { elastic: Some(ortho), plastic: Some(jc()), ..steel() }.validate().is_err());
        let damage = Damage::JohnsonCook { d1: 0.05, d2: 3.44, d3: -2.12, d4: 0.002, d5: 0.61 };
        Material { plastic: Some(jc()), damage: Some(damage), erosion: Some(Erosion { plastic_strain: 2.0 }), ..steel() }.validate().unwrap();
        assert!(Material { damage: Some(damage), ..steel() }.validate().is_err());
        let bilinear = Plastic::Bilinear { yield_stress: 250.0, hardening_modulus: 0.0, hardening: Hardening::Isotropic };
        assert!(Material { plastic: Some(bilinear.clone()), damage: Some(damage), ..steel() }.validate().is_err(), "D4, D5 need Johnson–Cook");
        let plain = Damage::JohnsonCook { d1: 0.5, d2: 0.0, d3: 0.0, d4: 0.0, d5: 0.0 };
        Material { plastic: Some(bilinear), damage: Some(plain), ..steel() }.validate().unwrap();
        assert!(Material { erosion: Some(Erosion { plastic_strain: 1.0 }), ..steel() }.validate().is_err());
        // JSON: defaults fill in, unknown fields are refused.
        let m: Material = serde_json::from_str(
            r#"{"elastic":{"type":"isotropic","young":1,"poisson":0.3},"plastic":{"type":"johnson_cook","a":1,"b":2,"n":0.5,"room_temperature":0,"melt_temperature":100}}"#,
        )
        .unwrap();
        assert!(
            matches!(m.plastic, Some(Plastic::JohnsonCook { m, reference_strain_rate, taylor_quinney, .. }) if m == 1.0 && reference_strain_rate == 1.0 && taylor_quinney == 0.9)
        );
        assert!(serde_json::from_str::<Material>(r#"{"density":1,"plastic":{"type":"bilinear","yield_stress":1,"oops":2}}"#).is_err());
        let text = serde_json::to_string(&Material { plastic: Some(jc()), damage: Some(damage), ..steel() }).unwrap();
        assert_eq!(serde_json::from_str::<Material>(&text).unwrap(), Material { plastic: Some(jc()), damage: Some(damage), ..steel() });
    }

    #[test]
    fn hyperelasticity_builds_on_isotropic_constants() {
        let m: Material = serde_json::from_value(serde_json::json!({
            "elastic": {"type": "isotropic", "young": 3.0, "poisson": 0.49}, "hyperelastic": {"type": "neo_hookean"}
        }))
        .unwrap();
        m.validate().unwrap();
        assert_eq!(m.hyperelastic, Some(Hyperelastic::NeoHookean {}));
        let back: Material = serde_json::from_value(serde_json::to_value(&m).unwrap()).unwrap();
        assert_eq!(back, m);
        assert!(serde_json::from_value::<Material>(serde_json::json!({"hyperelastic": {"type": "neo_hookean", "c10": 1}})).is_err());
        let alone = Material { hyperelastic: Some(Hyperelastic::NeoHookean {}), density: Some(1.0), ..Material::default() };
        assert!(alone.validate().is_err(), "needs isotropic constants");
        let ortho = Elastic::EngineeringConstants { young: [1.0; 3], poisson: [0.0; 3], shear: [1.0; 3] };
        assert!(Material { elastic: Some(ortho), ..alone }.validate().is_err());
        let plastic = Plastic::Bilinear { yield_stress: 1.0, hardening_modulus: 0.0, hardening: Hardening::Isotropic };
        assert!(Material { plastic: Some(plastic), ..m.clone() }.validate().is_err());
        // Mooney–Rivlin and Yeoh: μ = E/(2(1 + ν)) = 1.00671 here, so C01 is within [0, μ/2].
        let mr: Material = serde_json::from_value(serde_json::json!({
            "elastic": {"type": "isotropic", "young": 3.0, "poisson": 0.49}, "hyperelastic": {"type": "mooney_rivlin", "c01": 0.1}
        }))
        .unwrap();
        mr.validate().unwrap();
        let mu = Hyperelastic::moduli(3.0, 0.49).0;
        let c = Hyperelastic::MooneyRivlin { c01: 0.1 }.constants(mu);
        assert!((2.0 * (c[0] + c[1]) - mu).abs() < 1e-15 && c[1] == 0.1);
        let yeoh = Hyperelastic::Yeoh { c20: -0.01, c30: 0.001 };
        assert_eq!(yeoh.constants(mu), [mu / 2.0, 0.0, -0.01, 0.001]);
        assert!(Material { hyperelastic: Some(yeoh), ..m.clone() }.validate().is_ok());
        for bad in [
            Hyperelastic::MooneyRivlin { c01: -0.1 },
            Hyperelastic::MooneyRivlin { c01: 0.6 },
            Hyperelastic::MooneyRivlin { c01: f64::NAN },
            Hyperelastic::Yeoh { c20: f64::INFINITY, c30: 0.0 },
            Hyperelastic::Yeoh { c20: 0.0, c30: f64::NAN },
        ] {
            assert!(Material { hyperelastic: Some(bad), ..m.clone() }.validate().is_err(), "{bad:?}");
        }
    }

    #[test]
    fn thermal_properties() {
        let m: Material = serde_json::from_value(serde_json::json!({
            "thermal": {"conductivity": [1, 2, 3], "specific_heat": 460, "expansion": 1.2e-5, "reference_temperature": 20}
        }))
        .unwrap();
        let t = m.thermal.unwrap();
        assert_eq!(t.conductivity.unwrap().axes(), [1.0, 2.0, 3.0]);
        assert_eq!(t.expansion.unwrap().axes(), [1.2e-5; 3]);
        assert!(m.validate().is_ok(), "a heat-only material needs no elasticity");
        let back: Material = serde_json::from_value(serde_json::to_value(&m).unwrap()).unwrap();
        assert_eq!(back, m);
        for bad in [
            Thermal::default(),
            Thermal { conductivity: Some(Directional::Isotropic(0.0)), ..Thermal::default() },
            Thermal { conductivity: Some(Directional::Orthotropic([1.0, -1.0, 1.0])), ..Thermal::default() },
            Thermal { specific_heat: Some(f64::NAN), ..Thermal::default() },
            Thermal { expansion: Some(Directional::Isotropic(f64::INFINITY)), ..Thermal::default() },
            Thermal { expansion: Some(Directional::Isotropic(1e-5)), reference_temperature: f64::NAN, ..Thermal::default() },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
        assert!(serde_json::from_value::<Material>(serde_json::json!({"thermal": {"conductivity": [1, 2]}})).is_err());
        assert!(serde_json::from_value::<Material>(serde_json::json!({"thermal": {"heat": 1}})).is_err());
        // Negative expansion (invar-like behaviour of some composites) is allowed.
        assert!(Thermal { expansion: Some(Directional::Isotropic(-1e-6)), ..Thermal::default() }.validate().is_ok());
    }
}
