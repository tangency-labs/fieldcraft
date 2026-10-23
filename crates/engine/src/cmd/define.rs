//! Model definition commands: create, edit and delete sets, materials, sections, reference
//! points, couplings, point masses, springs and dashpots, boundary conditions, loads, fluid
//! domains and steps, and the step activity table.
//!
//! Each definition's parameters are its JSON form in the model (the same form a `.fieldcraft`
//! file stores), plus its `name`. `create` fails when the name is taken, `edit` replaces an
//! existing definition, `delete` refuses while another definition refers to it.

use fieldcraft_model::{
    BoundaryCondition, Coupling, FaceCondition, FluidBoundary, FluidDomain, Load, Material, Model, PointMass, ReferencePoint, Section, Set, Spring,
    Step, ThermalCondition,
};
use serde::de::DeserializeOwned;
use serde_json::{Value, json};

use fieldcraft_io::{Item, ItemKind};
use fieldcraft_model::rename::Kind as RenameKind;

use super::{name, opt, req};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// The parameters other than `name`, read as a definition's JSON form.
fn fields<T: DeserializeOwned>(p: &Value, what: &str) -> Result<T> {
    let mut obj = p.as_object().cloned().ok_or_else(|| EngineError::Params(format!("the parameters of a {what} are a JSON object")))?;
    obj.remove("name");
    serde_json::from_value(Value::Object(obj)).map_err(|e| EngineError::Params(format!("{what}: {e}")))
}

/// Create, edit and delete commands for one kind of definition.
macro_rules! definition_commands {
    ($what:literal, $ty:ty, $id:literal, $label:literal, $doc:literal, $form:expr, $delete_form:expr, $add:ident, $replace:ident, $remove:ident,
     $CREATE:ident, $EDIT:ident, $DELETE:ident, requires: $requires:expr, exists: $exists:ident) => {
        pub const $CREATE: CommandSpec = CommandSpec {
            id: concat!("model.", $id, ".create"),
            label: concat!("New ", $label),
            params: concat!($doc, " — fails if the name is taken"),
            form: Some($form),
            kind: CommandKind::Model,
            requires: $requires,
            run: |s: &mut Session, p: &Value| {
                let n = name(p)?;
                s.model.$add(n, fields::<$ty>(p, $what)?)?;
        Ok(json!({"name": n}))
            },
        };
        pub const $EDIT: CommandSpec = CommandSpec {
            id: concat!("model.", $id, ".edit"),
            label: concat!("Edit ", $label),
            params: concat!($doc, " — replaces the existing definition of that name"),
            form: Some($form),
            kind: CommandKind::Model,
            requires: &[Requirement::$exists],
            run: |s: &mut Session, p: &Value| {
                let n = name(p)?;
                s.model.$replace(n, fields::<$ty>(p, $what)?)?;
        Ok(json!({"name": n}))
            },
        };
        pub const $DELETE: CommandSpec = CommandSpec {
            id: concat!("model.", $id, ".delete"),
            label: concat!("Delete ", $label),
            params: r#"{"name": "…"} — refused while another definition refers to it"#,
            form: Some($delete_form),
            kind: CommandKind::Model,
            requires: &[Requirement::$exists],
            run: |s: &mut Session, p: &Value| {
                s.model.$remove(name(p)?)?;
                Ok(Value::Null)
            },
        };
    };
}

definition_commands!(
    "set",
    Set,
    "set",
    "Set",
    r#"{"name": "root", "geometry": [{"body": "Body1", "face": "F2:top"}, {"body": "Body1", "edge": "…"}, {"body": "Body1", "vertex": "…"}, {"body": "Body1"}], "nodes": [labels], "elements": [labels], "faces": [{"element": 1, "face": 3}], "edges": [{"element": 1, "edge": 2}]} — any of the member lists (element edges: where a 2D model's surface loads act); geometry members resolve against the mesh's geometry tags after every remesh (a face name also selects its split pieces `name#0`, `name#1`, …)"#,
    super::forms::SET,
    super::forms::SET_DELETE,
    add_set,
    replace_set,
    remove_set,
    SET_CREATE,
    SET_EDIT,
    SET_DELETE,
    requires: &[],
    exists: Set
);

definition_commands!(
    "material",
    Material,
    "material",
    "Material",
    r#"{"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3} | {"type": "engineering_constants", "young": [E1, E2, E3], "poisson": [ν12, ν13, ν23], "shear": [G12, G13, G23]}, "density": 7.85e-9, "fluid": {"type": "newtonian", "viscosity": 1.0e-9} (dynamic viscosity; a fluid needs a density), "plastic": {…} (see model.material.plasticity), "damage": {…}, "erosion": {…} (see model.material.damage), "hyperelastic": {"type": "neo_hookean"} | {"type": "mooney_rivlin", "c01": 0.1} (C10 = μ/2 − C01) | {"type": "yeoh", "c20": -0.01, "c30": 0.001} (C10 = μ/2) (large strain in nonlinear steps, from the isotropic constants: μ = E/(2(1 + ν)), κ = E/(3(1 − 2ν)); without it an elastic material is St Venant–Kirchhoff there), "thermal": {"conductivity": 0.05 | [kx, ky, kz], "specific_heat": 4.6e8, "expansion": 1.2e-5 | [αx, αy, αz], "reference_temperature": 20 (zero thermal strain; default 0 °C)} (heat transfer and thermal strain; each key optional)} — numbers in the display units (model.units), or with their unit ("210 GPa")"#,
    super::forms::MATERIAL,
    super::forms::MATERIAL_DELETE,
    add_material,
    replace_material,
    remove_material,
    MATERIAL_CREATE,
    MATERIAL_EDIT,
    MATERIAL_DELETE,
    requires: &[],
    exists: Material
);

/// A copy of material `name` with `change` applied, put back through the model's validation.
fn change_material(s: &mut Session, p: &Value, change: impl FnOnce(&mut Material, &Value) -> Result<()>) -> Result<Value> {
    let n = name(p)?;
    let mut m = s.model.material(n).cloned().ok_or_else(|| EngineError::Params(format!("material `{n}` does not exist")))?;
    change(&mut m, p)?;
    s.model.replace_material(n, m)?;
    Ok(json!({"name": n}))
}

pub const MATERIAL_PLASTICITY: CommandSpec = CommandSpec {
    id: "model.material.plasticity",
    label: "Plasticity",
    params: r#"{"name": "Steel", "plastic": {"type": "bilinear", "yield_stress": 250, "hardening_modulus": 1000, "hardening": "isotropic|kinematic"} | {"type": "table", "points": [{"strain": 0, "stress": 250}, {"strain": 0.1, "stress": 300}]} | {"type": "johnson_cook", "a": 792, "b": 510, "n": 0.26, "c": 0.014, "m": 1.03, "reference_strain_rate": 1, "room_temperature": 294, "melt_temperature": 1793, "specific_heat": 4.77e8, "taylor_quinney": 0.9} | null} — von Mises plasticity of an isotropic material (null removes it); explicit dynamics uses it, the linear solver the elastic constants alone"#,
    form: Some(super::forms::MATERIAL_PLASTICITY),
    kind: CommandKind::Model,
    requires: &[Requirement::Material],
    run: |s: &mut Session, p: &Value| {
        change_material(s, p, |m, p| {
            m.plastic = req(p, "plastic")?;
            Ok(())
        })
    },
};

pub const MATERIAL_DAMAGE: CommandSpec = CommandSpec {
    id: "model.material.damage",
    label: "Damage",
    params: r#"{"name": "Steel", "damage": {"type": "johnson_cook", "d1": 0.05, "d2": 3.44, "d3": -2.12, "d4": 0.002, "d5": 0.61} | null, "erosion": {"plastic_strain": 1.0} | null} — damage to failure and element erosion of a plastic material (each key optional; null removes it); a failed element is eroded in explicit dynamics"#,
    form: Some(super::forms::MATERIAL_DAMAGE),
    kind: CommandKind::Model,
    requires: &[Requirement::Material],
    run: |s: &mut Session, p: &Value| {
        change_material(s, p, |m, p| {
            // Absent leaves a field as it is; null removes it.
            let (damage, erosion) = (opt(p, "damage")?, opt(p, "erosion")?);
            if damage.is_none() && erosion.is_none() {
                return Err(EngineError::Params("give `damage`, `erosion` or both".into()));
            }
            if let Some(d) = damage {
                m.damage = d;
            }
            if let Some(e) = erosion {
                m.erosion = e;
            }
            Ok(())
        })
    },
};

definition_commands!(
    "section",
    Section,
    "section",
    "Section",
    r#"{"name": "web", "set": "all", "material": "Steel", "kind": {"type": "solid"} | {"type": "shell", "thickness": 2, "offset": 0 (-0.5..=0.5 of the thickness), "points": 5} | {"type": "beam", "profile": {"type": "rectangle", "width", "height"} | {"type": "circle", "radius"} | {"type": "pipe", "radius", "thickness"} | {"type": "i", "height", "width", "web", "flange"} | {"type": "box", "width", "height", "thickness"}, "orientation": [x, y, z]} | {"type": "rod", "area": 10} | {"type": "plane_stress", "thickness": 1} | {"type": "plane_strain", "thickness": 1 (default)} | {"type": "axisymmetric"} (2D: elements in the XY plane; axisymmetric about Y with x the radius, forces totals over the circumference), "formulation": {"integration": "full|reduced", "incompatible_modes": false, "hybrid": false}}"#,
    super::forms::SECTION,
    super::forms::SECTION_DELETE,
    add_section,
    replace_section,
    remove_section,
    SECTION_CREATE,
    SECTION_EDIT,
    SECTION_DELETE,
    requires: &[Requirement::Material, Requirement::Set],
    exists: Section
);

definition_commands!(
    "reference point",
    ReferencePoint,
    "reference_point",
    "Reference Point",
    r#"{"name": "RP", "xyz": [x, y, z], "label": 0} — label: its node label in decks (default: the highest free label, kept through remeshes)"#,
    super::forms::REFERENCE_POINT,
    super::forms::REFERENCE_POINT_DELETE,
    add_reference_point,
    replace_reference_point,
    remove_reference_point,
    REFERENCE_POINT_CREATE,
    REFERENCE_POINT_EDIT,
    REFERENCE_POINT_DELETE,
    requires: &[],
    exists: ReferencePoint
);

definition_commands!(
    "coupling",
    Coupling,
    "coupling",
    "Coupling",
    r#"{"name": "pad", "reference_point": "RP", "set": "hole", "kind": "kinematic|distributing", "dofs": [1, 2, 3, 4, 5, 6] (default all)}"#,
    super::forms::COUPLING,
    super::forms::COUPLING_DELETE,
    add_coupling,
    replace_coupling,
    remove_coupling,
    COUPLING_CREATE,
    COUPLING_EDIT,
    COUPLING_DELETE,
    requires: &[Requirement::ReferencePoint, Requirement::Set],
    exists: Coupling
);

definition_commands!(
    "point mass",
    PointMass,
    "mass",
    "Point Mass",
    r#"{"name": "motor", "target": {"set": "…"} | {"reference_point": "…"}, "mass": 2.5 | [mx, my, mz] (per global direction), "inertia": [I11, I22, I33, I12, I13, I23] (tensor about the node in global axes; needs rotations: a reference point)} — at every node of the target, each with the full value"#,
    super::forms::MASS,
    super::forms::MASS_DELETE,
    add_mass,
    replace_mass,
    remove_mass,
    MASS_CREATE,
    MASS_EDIT,
    MASS_DELETE,
    requires: &[Requirement::Target],
    exists: Mass
);

definition_commands!(
    "spring",
    Spring,
    "spring",
    "Spring or Dashpot",
    r#"{"name": "mount", "type": "spring", "stiffness": 1000 | "type": "dashpot", "coefficient": 5, "first": {"set": "…"} | {"reference_point": "…"}, "second": (optional; without it each node of `first` connects to the ground; with it each end is one node), "direction": {"dof": 1..6} (4–6: torsional) | {"vector": [x, y, z]} | "line" (along the two ends, two-point only)} — dashpots act in dynamic steps; static and frequency steps ignore them"#,
    super::forms::SPRING,
    super::forms::SPRING_DELETE,
    add_spring,
    replace_spring,
    remove_spring,
    SPRING_CREATE,
    SPRING_EDIT,
    SPRING_DELETE,
    requires: &[Requirement::Target],
    exists: Spring
);

definition_commands!(
    "boundary condition",
    BoundaryCondition,
    "boundary_condition",
    "Boundary Condition",
    r#"{"name": "clamp", "type": "fixed", "target": {"set": "root"} | {"reference_point": "RP"}} | {"name", "type": "displacement", "target", "values": [u1, u2, u3, ur1, ur2, ur3] (null leaves a degree of freedom free)} | {"name", "type": "symmetry", "target", "plane": "x|y|z" | {"normal": [nx, ny, nz]} (the plane's normal; an inclined plane holds n · u = 0 and the rotations about its in-plane axes)} — active in the steps that list it (model.step.activate)"#,
    super::forms::BOUNDARY_CONDITION,
    super::forms::BOUNDARY_CONDITION_DELETE,
    add_boundary_condition,
    replace_boundary_condition,
    remove_boundary_condition,
    BOUNDARY_CONDITION_CREATE,
    BOUNDARY_CONDITION_EDIT,
    BOUNDARY_CONDITION_DELETE,
    requires: &[Requirement::Target],
    exists: BoundaryCondition
);

definition_commands!(
    "load",
    Load,
    "load",
    "Load",
    r#"{"name": "push", "type": "force|moment", "target": {"set": "…"} | {"reference_point": "…"}, "vector": [x, y, z], "distribution": "total" (default: the vector is the total, split by area over the set's faces, or equally over its nodes) | "per_node" (at every node: depends on the mesh)} | {"name", "type": "traction", "set", "vector" (per area, on the set's faces)} | {"name", "type": "pressure", "set", "magnitude" (positive pushes in)} | {"name", "type": "gravity", "acceleration": [x, y, z], "set" (optional, default the whole model)} | {"name", "type": "body_force", "set", "vector" (per volume)} | {"name", "type": "centrifugal", "angular_velocity" (rad per time unit), "point": [x, y, z] (on the axis), "axis": [x, y, z], "set" (optional)} | {"name", "type": "surface_flux", "set", "flux" (heat per area into the set's faces)} | {"name", "type": "convection", "set", "film" (film coefficient h), "sink" (surrounding temperature): heat flux h (sink − T) into the set's faces} | {"name", "type": "heat_source", "set", "power" (heat per volume in the set's elements)} | {"name", "type": "temperature", "target": {"set": "…"} | {"reference_point": "…"}, "value" (a temperature that strains the material in stress steps)} — heat loads act in heat transfer and coupled steps; 2D models: in-plane vectors, surface loads on element edges; active in the steps that list it (model.step.activate)"#,
    super::forms::LOAD,
    super::forms::LOAD_DELETE,
    add_load,
    replace_load,
    remove_load,
    LOAD_CREATE,
    LOAD_EDIT,
    LOAD_DELETE,
    requires: &[],
    exists: Load
);

definition_commands!(
    "fluid domain",
    FluidDomain,
    "fluid_domain",
    "Fluid Domain",
    r#"{"name": "Flow", "material": "Water" (a density and fluid behaviour), "min": [x, y, z], "max": [x, y, z], "spacing": 1.0 (cubic cells; every extent a whole number of spacings), "faces": {"x_min": {"type": "velocity", "velocity": [u, v, w], "profile": "uniform|parabolic"}, "x_max": {"type": "pressure", "pressure": 0, "outflow": "uniform|convective" (default uniform: the face at the pressure, the flow crossing it normally; convective: each cell at the pressure, the flow leaving as it arrives, with reversed cells held normal: for outlets a wake or vortices reach)}, "y_min": {"type": "wall"} | {"type": "moving_wall", "velocity": [u, v, w] (in the face's plane)} | {"type": "slip"} | {"type": "periodic"} (in opposite pairs), …} (default walls), "obstacles": ["Body1"] (solid bodies of the design or imported surfaces, by name), "surfaces": [{"face": {"body": "Body1", "face": "F2:top"}, "condition": {"type": "moving_wall", "velocity": [u, v, w]}}] (obstacle faces are walls unless listed), "initial_velocity": [u, v, w], "acceleration": [x, y, z] (gravity or a driving force per mass), "tolerance": 0.1 (obstacle triangulation; default a tenth of the spacing), "thermal": {"initial_temperature": 20, "gravity": [0, -9810, 0], "expansion": 2.1e-4 (β per temperature; buoyancy needs both, Boussinesq), "faces": {"y_min": {"type": "temperature", "value": 80} | {"type": "flux", "value": 1000 (heat per area into the fluid)} | {"type": "adiabatic"}, …} (default adiabatic; a velocity inlet carries its held temperature in, a pressure outlet lets heat leave), "surfaces": [{"face": {"body": "Body1", "face": "F2:top"}, "condition": {"type": "temperature", "value": 80}}] (obstacle faces are adiabatic unless listed), "solids": [{"body": "Sink", "material": "Aluminium" (a density, a conductivity: isotropic or [kx, ky, kz] along the box's axes, and a specific heat), "power": 5 (heat generated, spread over its volume; default 0), "limit": 95 (optional: the highest temperature it may reach, a junction or case limit; the fluid solve reports its peak against it)}] (obstacles that conduct heat, coupled to the fluid across their faces, which take no thermal condition; a step run to steady state scales their heat capacity down, which the steady temperatures do not depend on), "contacts": [{"bodies": ["Die", "Lid"] (two touching conducting solids), "resistance": 1e-5 (per area: temperature drop over heat flux)} | {"bodies", "thickness": 0.05, "conductivity": 5 (a paste's or pad's layer: resistance thickness / conductivity)}] (thermal interface layers, not resolved in the lattice: each face between the two bodies' cells adds the resistance in series; solids not listed touch perfectly)} (heat in the fluid, carried by the flow and diffused with the material's conductivity and specific heat; left out, the flow is isothermal), "porous": [{"name": "Vent", "min": [x, y, z], "max": [x, y, z] (the cells whose centres it contains), "axis": "x|y|z" (the flow's direction through it), "resistance": {"type": "darcy_forchheimer", "permeability": K (length²), "inertial": C₂ (per length; default 0)} | {"type": "loss_coefficient", "coefficient": k (Δp = k ρ u²/2 across the zone)}}] (vents, grilles, filters, perforated plates: a force per volume −(μ/K) u − (ρ C₂/2)|u| u), "fans": [{"name": "Fan", "min": [x, y, z], "max": [x, y, z], "direction": "+x|-x|+y|-y|+z|-z", "curve": [[Q, Δp], …] (volume flow increasing from 0, pressure rise not increasing; linear between points; beyond the last point the fan holds its end value, so it never brakes the flow past free delivery), "relaxation": 0.01 (how far the rise moves toward the curve per time step)}] (the pressure rise pushes the fan's cells; the solve reports each fan's operating point and each vent's flow); fans and zones stay off the inlet and outlet layers; numbers in the display units (model.units), or with their unit ("2 m/s")}"#,
    super::forms::FLUID_DOMAIN,
    super::forms::FLUID_DOMAIN_DELETE,
    add_fluid_domain,
    replace_fluid_domain,
    remove_fluid_domain,
    FLUID_DOMAIN_CREATE,
    FLUID_DOMAIN_EDIT,
    FLUID_DOMAIN_DELETE,
    requires: &[Requirement::Material],
    exists: FluidDomain
);

pub const FLUID_BOUNDARY: CommandSpec = CommandSpec {
    id: "model.fluid_domain.boundary",
    label: "Inlet, Outlet or Wall",
    params: r#"{"domain": "Flow", "at": "x_min" (a box face) | {"body": "Body1", "face": "F2:top"} (an obstacle face), "condition": {"type": "velocity", "velocity": [u, v, w], "profile": "uniform|parabolic"} | {"type": "pressure", "pressure": 0, "outflow": "uniform|convective"} | {"type": "wall"} | {"type": "moving_wall", "velocity": [u, v, w]} | {"type": "slip"} | {"type": "periodic"}} — sets one boundary's condition (obstacle faces take walls and moving walls)"#,
    form: Some(super::forms::FLUID_BOUNDARY),
    kind: CommandKind::Model,
    requires: &[Requirement::FluidDomain],
    run: fluid_boundary,
};

fn fluid_boundary(s: &mut Session, p: &Value) -> Result<Value> {
    let domain: &str = req(p, "domain")?;
    let at: FluidBoundary = req(p, "at")?;
    let condition: FaceCondition = req(p, "condition")?;
    s.model.set_fluid_boundary(domain, at, condition)?;
    Ok(json!({"domain": domain}))
}

pub const FLUID_THERMAL_BOUNDARY: CommandSpec = CommandSpec {
    id: "model.fluid_domain.thermal_boundary",
    label: "Temperature or Heat Flux",
    params: r#"{"domain": "Flow", "at": "x_min" (a box face) | {"body": "Body1", "face": "F2:top"} (an obstacle face), "condition": {"type": "temperature", "value": 80} | {"type": "flux", "value": 1000 (heat per area into the fluid)} | {"type": "adiabatic"}} — sets one boundary's thermal condition; the domain needs `thermal` settings (model.fluid_domain.edit)"#,
    form: Some(super::forms::FLUID_THERMAL_BOUNDARY),
    kind: CommandKind::Model,
    requires: &[Requirement::FluidDomain],
    run: fluid_thermal_boundary,
};

fn fluid_thermal_boundary(s: &mut Session, p: &Value) -> Result<Value> {
    let domain: &str = req(p, "domain")?;
    let at: FluidBoundary = req(p, "at")?;
    let condition: ThermalCondition = req(p, "condition")?;
    s.model.set_fluid_thermal_boundary(domain, at, condition)?;
    Ok(json!({"domain": domain}))
}

pub const FLUID_COMPONENT_LIMIT: CommandSpec = CommandSpec {
    id: "model.fluid_domain.limit",
    label: "Temperature Limit",
    params: r#"{"domain": "Flow", "body": "Die" (a conducting solid of the domain), "limit": 95 (the highest temperature it may reach, a junction or case limit; left out, none)} — sets one conducting solid's temperature limit; the fluid solve reports each solid's peak, its margin to the limit and whether it passes, and whether every component does (`components_ok`)"#,
    form: Some(super::forms::FLUID_COMPONENT_LIMIT),
    kind: CommandKind::Model,
    requires: &[Requirement::FluidDomain],
    run: fluid_component_limit,
};

fn fluid_component_limit(s: &mut Session, p: &Value) -> Result<Value> {
    let domain: &str = req(p, "domain")?;
    let body: &str = req(p, "body")?;
    let limit: Option<f64> = opt(p, "limit")?;
    s.model.set_component_limit(domain, body, limit)?;
    Ok(json!({"domain": domain, "body": body}))
}

pub const STEP_CREATE: CommandSpec = CommandSpec {
    id: "model.step.create",
    label: "New Step",
    params: r#"{"name": "Bend", "procedure": {"type": "static"} | {"type": "nonlinear_static", "nonlinear_geometry": true, "initial_increment": 0.1, "min_increment": 1e-5, "max_increment": 1.0 (fractions of the step), "max_increments": 100, "max_iterations": 16 (per increment), "tolerances": {"force": 1e-6, "displacement": 1e-5, "energy": 1e-10}, "line_search": false, "load_path": [[0, 0], [0.5, 1], [1, 0]] (optional: the load factor against the step fraction, piecewise linear from [0, 0] to fraction 1; increments end at its points)} (loads and prescribed displacements are the load factor times their values; it ramps from 0 to 1 without a load path) | {"type": "frequency", "modes": 10, "min_hz": 0, "max_hz": 500} | {"type": "buckle", "modes": 5, "shift"?: 0} | {"type": "heat_transfer", "transient"?: {"end_time": 10, "increment": 0.5, "theta": 1, "initial_temperature": 20}} (steady without `transient`) | {"type": "coupled_temperature_displacement"} (steady heat and the static equilibrium of its thermal strains) | {"type": "fluid", "domain": "Flow", "end_time": 1.0 | "steps": 1000, "time_step": (default: chosen for Mach 0.087 and τ ≤ 1), "frames": 10, "collision": "trt" (default) | "regularized" (recursive regularised: stable with τ close to ½), "magic": 0.1875 (the two-relaxation-time Λ), "turbulence": {"type": "smagorinsky", "constant": 0.1, "prandtl": 0.9} (optional: Smagorinsky's subgrid eddy viscosity (C Δx)² |S|, C above 0 and at most 0.5; with heat the turbulent diffusivity ν_t / Pr_t, Pr_t 0.1 to 10, default 0.9; it lets the fluid's own τ go down to 0.5001; the reply's eddy_viscosity_ratio is the largest ν_t / ν; CPU backend only), "steady": {"tolerance": 1e-3, "window": (flow time the change is measured over, sampled 8 times; default a twentieth of the step), "min_time": 0, "divergence": (optional relative change that stops the run as diverged), "monitors": [{"type": "mass_flow", "face": "x_max"} | {"type": "force", "faces": ["x_min", {"body": "Body1", "face": "F3"}]} | {"type": "heat", "faces": […]} | {"type": "mean_temperature"}]} (optional: stop once the velocity residual, the temperature residual with heat, and every monitor vary by less than the tolerance (relative change per window) for two windows in a row; end_time / steps is then the most it runs)} | {"type": "explicit_dynamic", "end_time": 0.001, "dt_scale": 0.9, "max_increments": 10000000, "output_interval", "history_interval", "history_nodes": "set", "initial_velocities": [{"set", "velocity": [vx, vy, vz]}], "rigid_walls": [{"set", "point", "normal", "friction": {"static": 0.3, "kinetic": 0.2, "decay": 10} (Coulomb; default frictionless)}], "hourglass": {"type": "viscous|stiffness", "coefficient": 0.1}, "bulk_viscosity": {"linear": 0.06, "quadratic": 1.5}, "mass_scaling": {"target_increment": 1e-6, "set": "…" (optional, default every element)} (opt-in: elements below the target increment get the density that reaches it; the added mass is reported)}, "output": ["displacement", "reaction_force", "stress", "strain"] (default by procedure), "boundary_conditions": [names], "loads": [names] (active in this step; nothing carries over between steps), "before": "step" (default: last)}"#,
    form: Some(super::forms::STEP_CREATE),
    kind: CommandKind::Model,
    requires: &[],
    run: step_create,
};

pub const STEP_EDIT: CommandSpec = CommandSpec {
    id: "model.step.edit",
    label: "Edit Step",
    params: r#"{"name": "Bend", "procedure", "output", "boundary_conditions", "loads"} — as model.step.create; what is left out keeps its value"#,
    form: Some(super::forms::STEP_EDIT),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: step_edit,
};

pub const STEP_DELETE: CommandSpec = CommandSpec {
    id: "model.step.delete",
    label: "Delete Step",
    params: r#"{"name": "Bend"}"#,
    form: Some(super::forms::STEP_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: |s, p| s.model.remove_step(name(p)?).map(|_| Value::Null).map_err(Into::into),
};

pub const STEP_ACTIVATE: CommandSpec = CommandSpec {
    id: "model.step.activate",
    label: "Activate in Step",
    params: r#"{"step": "Bend", "boundary_conditions": [names], "loads": [names], "active": true} — switches the listed items on (or off with false) in that step only"#,
    form: Some(super::forms::STEP_ACTIVATE),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: step_activate,
};

pub const RENAME: CommandSpec = CommandSpec {
    id: "model.rename",
    label: "Rename",
    params: r#"{"kind": "set" | "material" | "section" | "step" | "boundary_condition" | "load" | "reference_point" | "coupling" | "contact" | "tie" | "point_mass" | "spring" | "fluid_domain" | "datum", "from": "end", "to": "Tip"} — renames the definition and every reference to it (sets in sections, couplings, contacts, ties, discrete elements, boundary conditions, loads and explicit steps; materials in sections and fluid domains; reference points; boundary conditions and loads in the steps they are active in; steps in imperfections; fluid domains in fluid steps; datum planes in partitions). A name already taken is refused; the viewport's hidden items follow the new name"#,
    form: Some(super::forms::RENAME),
    kind: CommandKind::Model,
    requires: &[],
    run: rename,
};

/// `model.rename`'s kinds, as the command takes them.
pub const RENAME_KINDS: &[(&str, RenameKind)] = &[
    ("set", RenameKind::Set),
    ("material", RenameKind::Material),
    ("section", RenameKind::Section),
    ("step", RenameKind::Step),
    ("boundary_condition", RenameKind::BoundaryCondition),
    ("load", RenameKind::Load),
    ("reference_point", RenameKind::ReferencePoint),
    ("coupling", RenameKind::Coupling),
    ("contact", RenameKind::Contact),
    ("tie", RenameKind::Tie),
    ("point_mass", RenameKind::PointMass),
    ("spring", RenameKind::Spring),
    ("fluid_domain", RenameKind::FluidDomain),
    ("datum", RenameKind::Datum),
];

fn rename(s: &mut Session, p: &Value) -> Result<Value> {
    let kind: &str = req(p, "kind")?;
    let from: &str = req(p, "from")?;
    let to: &str = req(p, "to")?;
    let Some((_, k)) = RENAME_KINDS.iter().find(|(id, _)| *id == kind) else {
        let kinds: Vec<&str> = RENAME_KINDS.iter().map(|(id, _)| *id).collect();
        return Err(EngineError::Params(format!("`kind` is one of {}", kinds.join(", "))));
    };
    s.model.rename(*k, from, to)?;
    // What the viewport hides follows the new name.
    let item = match k {
        RenameKind::Datum => Some(ItemKind::Datum),
        RenameKind::BoundaryCondition => Some(ItemKind::BoundaryCondition),
        RenameKind::Load => Some(ItemKind::Load),
        RenameKind::ReferencePoint => Some(ItemKind::ReferencePoint),
        RenameKind::Coupling => Some(ItemKind::Coupling),
        _ => None,
    };
    if let Some(kind) = item
        && s.view.hidden.remove(&Item { kind, name: from.to_string() })
    {
        s.view.hidden.insert(Item { kind, name: to.to_string() });
    }
    Ok(json!({"kind": kind, "from": from, "name": to}))
}

pub const ACTIVITY: CommandSpec = CommandSpec {
    id: "model.activity",
    label: "Step Table",
    params: "{}",
    form: Some(&[]),
    kind: CommandKind::Model,
    requires: &[],
    run: |s, _| Ok(activity(&s.model)),
};

fn step_create(s: &mut Session, p: &Value) -> Result<Value> {
    let n = name(p)?;
    let mut step = Step::new(n, req(p, "procedure")?);
    if let Some(o) = opt(p, "output")? {
        step.output = o;
    }
    step.boundary_conditions = opt(p, "boundary_conditions")?.unwrap_or_default();
    step.loads = opt(p, "loads")?.unwrap_or_default();
    s.model.add_step(step, opt(p, "before")?)?;
    Ok(json!({"name": n}))
}

fn step_edit(s: &mut Session, p: &Value) -> Result<Value> {
    let n = name(p)?;
    let mut step = s.model.step(n).cloned().ok_or_else(|| EngineError::Params(format!("no step `{n}`")))?;
    if let Some(procedure) = opt(p, "procedure")? {
        step.set_procedure(procedure);
    }
    if let Some(o) = opt(p, "output")? {
        step.output = o;
    }
    if let Some(b) = opt(p, "boundary_conditions")? {
        step.boundary_conditions = b;
    }
    if let Some(l) = opt(p, "loads")? {
        step.loads = l;
    }
    s.model.replace_step(step)?;
    Ok(json!({"name": n}))
}

fn step_activate(s: &mut Session, p: &Value) -> Result<Value> {
    let step: &str = req(p, "step")?;
    let active = opt(p, "active")?.unwrap_or(true);
    let bcs: Vec<&str> = opt(p, "boundary_conditions")?.unwrap_or_default();
    let loads: Vec<&str> = opt(p, "loads")?.unwrap_or_default();
    if bcs.is_empty() && loads.is_empty() {
        return Err(EngineError::Params("list `boundary_conditions` or `loads` to switch".into()));
    }
    for b in &bcs {
        s.model.set_boundary_condition_active(step, b, active)?;
    }
    for l in &loads {
        s.model.set_load_active(step, l, active)?;
    }
    Ok(activity(&s.model))
}

/// The step × boundary-condition/load table: the steps in order, and for each item the steps it
/// is active in.
fn activity(m: &Model) -> Value {
    let row =
        |active: &dyn Fn(&Step) -> bool| -> Value { m.steps().iter().filter(|s| active(s)).map(|s| s.name.as_str()).collect::<Vec<_>>().into() };
    let bcs: serde_json::Map<String, Value> =
        m.boundary_conditions().keys().map(|n| (n.clone(), row(&|s| s.boundary_conditions.contains(n)))).collect();
    let loads: serde_json::Map<String, Value> = m.loads().keys().map(|n| (n.clone(), row(&|s| s.loads.contains(n)))).collect();
    json!({
        "steps": m.steps().iter().map(|s| json!({"name": s.name, "procedure": s.procedure.name()})).collect::<Vec<_>>(),
        "boundary_conditions": bcs,
        "loads": loads,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use serde_json::json;

    use crate::Session;

    /// A two-element bar tagged with body `Bar` and faces `left` and `right`.
    pub(crate) fn tagged_bar(s: &mut Session) {
        for (i, x) in [0.0, 1.0, 2.0].into_iter().enumerate() {
            for (j, (y, z)) in [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)].into_iter().enumerate() {
                s.execute("mesh.node.add", &json!({"label": 1 + i * 4 + j, "xyz": [x, y, z]})).unwrap();
            }
        }
        // Node (i, y, z) = 1 + 4i + 2y + z.
        let n = |i: u32, y: u32, z: u32| 1 + 4 * i + 2 * y + z;
        for e in 0..2 {
            let nodes = [n(e, 0, 0), n(e + 1, 0, 0), n(e + 1, 1, 0), n(e, 1, 0), n(e, 0, 1), n(e + 1, 0, 1), n(e + 1, 1, 1), n(e, 1, 1)];
            s.execute("mesh.element.add", &json!({"label": e + 1, "kind": "hex8", "nodes": nodes})).unwrap();
        }
        s.execute("mesh.tag", &json!({"geometry": {"body": "Bar"}, "elements": [1, 2]})).unwrap();
        s.execute("mesh.tag", &json!({"geometry": {"body": "Bar", "face": "left"}, "faces": [{"element": 1, "face": 6}]})).unwrap();
        s.execute("mesh.tag", &json!({"geometry": {"body": "Bar", "face": "right"}, "faces": [{"element": 2, "face": 4}]})).unwrap();
    }

    pub(crate) fn defined_bar(s: &mut Session) {
        tagged_bar(s);
        let script = json!([
            {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
            {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Bar"}]}},
            {"command": "model.set.create", "params": {"name": "root", "geometry": [{"body": "Bar", "face": "left"}]}},
            {"command": "model.set.create", "params": {"name": "end", "geometry": [{"body": "Bar", "face": "right"}]}},
            {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}, "density": 7.85e-9}},
            {"command": "model.section.create", "params": {"name": "bar", "set": "all", "material": "Steel", "kind": {"type": "solid"}, "formulation": {"incompatible_modes": true}}},
            {"command": "model.reference_point.create", "params": {"name": "RP", "xyz": [3, 0.5, 0.5]}},
            {"command": "model.coupling.create", "params": {"name": "pad", "reference_point": "RP", "set": "end", "kind": "kinematic"}},
            {"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}},
            {"command": "model.load.create", "params": {"name": "push", "type": "force", "target": {"reference_point": "RP"}, "vector": [0, 0, -100]}},
            {"command": "model.load.create", "params": {"name": "press", "type": "pressure", "set": "end", "magnitude": 2}},
            {"command": "model.step.create", "params": {"name": "Bend", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["push"]}},
            {"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 6}, "boundary_conditions": ["clamp"]}},
        ]);
        s.run_script(&script).unwrap();
    }

    /// A rename rewrites the references, checks clean, is undone in one step, follows the hidden
    /// state and refuses a name in use (#290).
    #[test]
    fn rename_keeps_the_model_whole_and_undoes() {
        let mut s = Session::default();
        defined_bar(&mut s);
        s.execute("view.hide", &json!({"items": [{"kind": "load", "name": "push"}]})).unwrap();
        s.execute("model.rename", &json!({"kind": "set", "from": "end", "to": "Tip"})).unwrap();
        s.execute("model.rename", &json!({"kind": "load", "from": "push", "to": "Push"})).unwrap();
        s.execute("model.rename", &json!({"kind": "step", "from": "Bend", "to": "Load"})).unwrap();
        assert!(s.model.sets().get("end").is_none() && !s.model.set_users("Tip").is_empty());
        assert!(s.model.step("Load").is_some_and(|st| st.loads.contains("Push")));
        assert!(s.view.is_hidden(fieldcraft_io::ItemKind::Load, "Push") && !s.view.is_hidden(fieldcraft_io::ItemKind::Load, "push"));
        assert_eq!(s.execute("model.check", &json!({})).unwrap()["ready"], true);
        let e = s.execute("model.rename", &json!({"kind": "set", "from": "Tip", "to": "root"})).unwrap_err().to_string();
        assert!(e.contains("already exists"), "{e}");
        assert!(s.execute("model.rename", &json!({"kind": "partition", "from": "a", "to": "b"})).is_err());
        s.execute("edit.undo", &json!({})).unwrap();
        assert!(s.model.step("Bend").is_some() && s.model.step("Load").is_none());
    }

    #[test]
    fn a_script_defines_a_model_that_checks_clean() {
        let mut s = Session::default();
        defined_bar(&mut s);
        let check = s.execute("model.check", &json!({})).unwrap();
        assert_eq!(check["ready"], true, "{check}");
        // `press` is defined but in no step: a warning, not an error.
        assert_eq!(check["warnings"], 1, "{check}");
        assert!(check["diagnostics"][0]["subject"] == "load `press`", "{check}");

        let r = s.execute("model.set.resolve", &json!({"name": "end"})).unwrap();
        assert_eq!(r["nodes"], json!({"count": 4, "labels": [9, 10, 11, 12]}));
        assert_eq!(r["faces"]["labels"], json!([{"element": 2, "face": 4}]));

        let table = s.execute("model.step.activate", &json!({"step": "Bend", "loads": ["press"]})).unwrap();
        assert_eq!(table["loads"]["press"], json!(["Bend"]));
        assert_eq!(table["boundary_conditions"]["clamp"], json!(["Bend", "Modes"]));
        assert_eq!(s.execute("model.check", &json!({})).unwrap()["warnings"], 0);
        let off = s.execute("model.step.activate", &json!({"step": "Bend", "loads": ["press", "push"], "active": false})).unwrap();
        assert_eq!(off["loads"], json!({"press": [], "push": []}));

        let info = s.execute("model.info", &json!({})).unwrap();
        assert_eq!((info["display_units"].as_str(), info["own_units"].as_bool()), (Some("mm-t-s-N-MPa"), Some(true)));
        assert_eq!(info["definitions"]["sets"], 3);
        assert_eq!(info["definitions"]["steps"], 2);
    }

    #[test]
    fn edits_and_deletes_keep_the_model_consistent() {
        let mut s = Session::default();
        defined_bar(&mut s);
        let journal = s.journal.len();

        // Refused: a taken name, a missing reference, a referenced item, a typo in a field.
        assert!(s.execute("model.set.create", &json!({"name": "all", "nodes": [1]})).is_err());
        assert!(s.execute("model.section.create", &json!({"name": "x", "set": "nope", "material": "Steel", "kind": {"type": "solid"}})).is_err());
        let e = s.execute("model.set.delete", &json!({"name": "end"})).unwrap_err().to_string();
        assert!(e.contains("coupling `pad`") && e.contains("load `press`"), "{e}");
        assert!(
            s.execute("model.material.create", &json!({"name": "Al", "elastic": {"type": "isotropic", "young": 70000, "poison": 0.33}})).is_err()
        );
        assert!(s.execute("model.step.edit", &json!({"name": "Bend", "loads": ["nope"]})).is_err());
        assert!(s.execute("model.step.create", &json!({"name": "Late", "procedure": {"type": "static"}, "before": "nope"})).is_err());
        assert_eq!(s.journal.len(), journal, "refusals change nothing");

        // Edits replace; the step keeps its table entries unless told otherwise.
        s.execute(
            "model.material.edit",
            &json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 200000, "poisson": 0.29}, "density": 7.8e-9}),
        )
        .unwrap();
        s.execute("model.step.edit", &json!({"name": "Modes", "procedure": {"type": "frequency", "modes": 12, "max_hz": 1000}})).unwrap();
        assert!(s.model.step("Modes").unwrap().boundary_conditions.contains("clamp"));
        s.execute(
            "model.step.create",
            &json!({"name": "Pre", "procedure": {"type": "static"}, "before": "Bend", "boundary_conditions": ["clamp"], "loads": ["press"]}),
        )
        .unwrap();
        assert_eq!(s.model.steps().iter().map(|st| st.name.as_str()).collect::<Vec<_>>(), ["Pre", "Bend", "Modes"]);

        // A load listed by a step is not deleted until no step lists it.
        assert!(s.execute("model.load.delete", &json!({"name": "press"})).unwrap_err().to_string().contains("step `Pre`"));
        s.execute("model.step.activate", &json!({"step": "Pre", "loads": ["press"], "active": false})).unwrap();
        s.execute("model.load.delete", &json!({"name": "press"})).unwrap();
        s.execute("model.coupling.delete", &json!({"name": "pad"})).unwrap();
        s.execute("model.set.delete", &json!({"name": "end"})).unwrap();
        s.execute("model.step.delete", &json!({"name": "Pre"})).unwrap();
        assert!(s.execute("model.step.delete", &json!({"name": "Pre"})).is_err());
    }

    #[test]
    fn the_check_names_the_bad_items_and_sets_follow_a_retag() {
        let mut s = Session::default();
        defined_bar(&mut s);
        // Point the section at a face set: it has no elements, and the elements lose their section.
        s.execute("model.section.edit", &json!({"name": "bar", "set": "root", "material": "Steel", "kind": {"type": "solid"}})).unwrap();
        let check = s.execute("model.check", &json!({})).unwrap();
        assert_eq!(check["ready"], false);
        let messages: Vec<String> = check["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .map(|d| format!("{}: {}", d["subject"].as_str().unwrap(), d["message"].as_str().unwrap()))
            .collect();
        assert!(messages.contains(&"section `bar`: set `root` has no elements".to_string()), "{messages:#?}");
        assert!(messages.contains(&"model: elements 1, 2 have no section".to_string()), "{messages:#?}");

        // Clearing the tags leaves geometry sets unresolved, and says so.
        s.execute("mesh.tag.clear", &json!({})).unwrap();
        let e = s.execute("model.set.resolve", &json!({"name": "root"})).unwrap_err().to_string();
        assert!(e.contains("face left of body Bar is not on the mesh (the mesh carries no geometry tags)"), "{e}");
    }

    #[test]
    fn plasticity_and_damage_edit_a_material() {
        let mut s = Session::default();
        defined_bar(&mut s);
        let bilinear = json!({"type": "bilinear", "yield_stress": 250, "hardening_modulus": 1000});
        s.execute("model.material.plasticity", &json!({"name": "Steel", "plastic": bilinear})).unwrap();
        assert!(s.model.material("Steel").unwrap().plastic.is_some());
        // Damage needs plasticity; D4 and D5 need Johnson–Cook: refused, the model unchanged.
        let before = s.model.clone();
        let jc_damage = json!({"type": "johnson_cook", "d1": 0.05, "d2": 3.44, "d3": -2.12, "d4": 0.002, "d5": 0.61});
        assert!(s.execute("model.material.damage", &json!({"name": "Steel", "damage": jc_damage})).is_err());
        assert_eq!(s.model, before);
        s.execute(
            "model.material.damage",
            &json!({"name": "Steel", "damage": {"type": "johnson_cook", "d1": 0.3, "d2": 0, "d3": 0}, "erosion": {"plastic_strain": 0.5}}),
        )
        .unwrap();
        let m = s.model.material("Steel").unwrap();
        assert!(m.damage.is_some() && m.erosion.is_some());
        s.execute("model.material.damage", &json!({"name": "Steel", "erosion": null})).unwrap();
        assert!(s.model.material("Steel").unwrap().erosion.is_none() && s.model.material("Steel").unwrap().damage.is_some());
        // Removing plasticity under damage is refused.
        assert!(s.execute("model.material.plasticity", &json!({"name": "Steel", "plastic": null})).is_err());
        for bad in
            [json!({"name": "Steel"}), json!({"name": "Nope", "plastic": null}), json!({"name": "Steel", "plastic": {"type": "table", "points": []}})]
        {
            assert!(s.execute("model.material.plasticity", &bad).is_err(), "{bad}");
        }
        assert!(s.execute("model.material.damage", &json!({"name": "Steel"})).is_err());
    }
}
