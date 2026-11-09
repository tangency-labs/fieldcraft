//! Original StVK rod user-element export from the public Abaqus UEL/*DYNAMIC contracts.
//! Standard T3D2 finite-strain elasticity is not the same constitutive equation.
use crate::{Error, Result};
use fieldcraft_model::dynamics::{Integration, Method};
use fieldcraft_model::{BoundaryCondition, Elastic, ElementKind, Load, MassKind, Model, Procedure, SectionKind};
use std::fmt::Write;
mod uel;
pub use uel::SOURCE as USER_ELEMENT_SOURCE;

fn refuse(option: impl std::fmt::Display) -> Error {
    Error::Unsupported { format: "Abaqus nonlinear rod dynamics", what: option.to_string() }
}

/// SI input plus an original companion user element; never a native results payload.
pub fn write(model: &Model, step: Option<&str>, mass: MassKind, representation: &str) -> Result<String> {
    if representation != "uel" {
        return Err(refuse(format!(
            "representation `{representation}`: T3D2 ELASTIC is not the StVK Green-strain/second-Piola rod; use the explicit uel companion"
        )));
    }
    if model.nodes().len() > 10_000 || model.elements().len() > 10_000 {
        return Err(refuse("nonlinear companion export limit: 10000 nodes/elements"));
    }
    let r = model.resolve()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::Transient(_)), "transient")?;
    let Procedure::Transient(t) = &step.procedure else {
        return Err(refuse("transient step required"));
    };
    t.validate()?;
    for n in t.initial_displacement.keys().chain(t.initial_velocity.keys()) {
        if !model.nodes().contains_key(n) {
            return Err(refuse(format!("initial state node {n} missing")));
        }
    }
    if !matches!(t.method, Method::Direct) || t.adaptive.is_some() || t.base_acceleration.is_some() || t.base_displacement.is_some() {
        return Err(refuse("modal/adaptive/base-motion options"));
    }
    if t.damping != Default::default() || !t.load_amplitudes.is_empty() {
        return Err(refuse("physical damping or load_amplitudes: this user-element export supports undamped constant dead forces only"));
    }
    if t.initial_velocity.values().any(|v| v.iter().any(|x| *x != 0.)) {
        return Err(refuse("nonzero initial_velocity: the displacement preparation step has zero final velocity"));
    }
    if !model.contacts().is_empty()
        || !model.couplings().is_empty()
        || !model.ties().is_empty()
        || !model.foundations().is_empty()
        || !model.initial_fields().is_empty()
        || !model.pretensions().is_empty()
        || !model.reference_points().is_empty()
        || !model.springs().is_empty()
        || !model.masses().is_empty()
    {
        return Err(refuse("contact/coupling/tie/foundation/initial-field/pretension/reference-point/spring/point-mass definitions"));
    }
    if step.output.iter().any(|v| !matches!(v, fieldcraft_model::OutputVariable::Displacement)) {
        return Err(refuse("requested stress/other field: UEL has no native Cauchy stress field; recover Green/second-Piola measures explicitly"));
    }
    let connected: std::collections::BTreeSet<_> = model.elements().values().flat_map(|e| e.nodes.iter().copied()).collect();
    for n in model.nodes().keys() {
        if !connected.contains(n) {
            return Err(refuse(format!("unconnected node {n}")));
        }
    }
    let (alpha, beta, gamma) = match t.integration {
        Integration::AverageAcceleration => (0., 0.25, 0.5),
        Integration::Hht { alpha } => (alpha, (1. - alpha).powi(2) / 4., 0.5 - alpha),
        _ => return Err(refuse("linear_acceleration integration")),
    };
    let mut out = String::from(
        "*HEADING\nFieldCraft nonlinear StVK rod dynamics; SI m kg s N Pa\n** Requires the generated UEL companion; not a T3D2 material mapping.\n** Run with output_precision=full; history output and timestamps have separate precision limits.\n** Compare step FC_MOTION only; FC_INITIALIZE imposes u0 and the same dead loads.\n*NODE\n",
    );
    for (n, x) in model.nodes() {
        writeln!(out, "{n}, {:.17e}, {:.17e}, {:.17e}", x.xyz[0], x.xyz[1], x.xyz[2]).ok();
    }
    out.push_str("*NSET, NSET=FC_NODES\n");
    for n in model.nodes().keys() {
        writeln!(out, "{n}").ok();
    }
    out.push_str("*USER ELEMENT, TYPE=U1, NODES=2, COORDINATES=3, PROPERTIES=4, VARIABLES=6\n1, 2, 3\n");
    for (id, e) in model.elements() {
        if e.kind != ElementKind::Rod2 || e.nodes.len() != 2 {
            return Err(refuse(format!("element {id} `{}`", e.kind.name())));
        }
        let section = model
            .sections()
            .values()
            .find(|s| r.sets.get(&s.set).is_some_and(|set| set.elements.contains(id)))
            .ok_or_else(|| refuse(format!("element {id} section")))?;
        let SectionKind::Rod { area } = section.kind else {
            return Err(refuse(format!("element {id} non-rod section")));
        };
        let material = model.material(&section.material).ok_or_else(|| refuse("missing material"))?;
        let Some(Elastic::Isotropic { young, .. }) = material.elastic else {
            return Err(refuse("non-isotropic axial material"));
        };
        if material.plastic.is_some() || material.hyperelastic.is_some() || material.damping.is_some_and(|d| d != Default::default()) {
            return Err(refuse(format!("material `{}` plastic/hyperelastic/damping", section.material)));
        }
        let rho = material.density.ok_or_else(|| refuse("rod density"))?;
        writeln!(
            out,
            "*ELEMENT, TYPE=U1, ELSET=FC_E{id}\n{id}, {}, {}\n*UEL PROPERTY, ELSET=FC_E{id}\n{young:.17e}, {area:.17e}, {rho:.17e}, {}",
            e.nodes[0],
            e.nodes[1],
            if mass == MassKind::Lumped { 1 } else { 0 }
        )
        .ok();
    }
    out.push_str("*ELSET, ELSET=FC_ELEMENTS\n");
    for id in model.elements().keys() {
        writeln!(out, "{id}").ok();
    }
    let mut fixed = Vec::new();
    for name in &step.boundary_conditions {
        let bc = model.boundary_condition(name).ok_or_else(|| refuse(format!("boundary `{name}`")))?;
        if matches!(bc, BoundaryCondition::Temperature { .. }) || (1..=6).any(|d| bc.value(d).is_some_and(|v| v != 0.)) {
            return Err(refuse(format!("boundary `{name}`: only homogeneous structural supports")));
        }
        for n in r.target_nodes(bc.target()) {
            for d in 1..=3 {
                if bc.dofs().contains(d) {
                    if t.initial_displacement.get(n).is_some_and(|u| u[usize::from(d - 1)] != 0.) {
                        return Err(refuse(format!("initial_displacement conflicts with boundary `{name}` at node {n}")));
                    }
                    fixed.push((*n, d));
                }
            }
        }
    }
    for name in &step.loads {
        if !matches!(model.load(name), Some(Load::Force { .. })) {
            return Err(refuse(format!("load `{name}`: only dead nodal force")));
        }
    }
    out.push_str("*FILE FORMAT, ZERO INCREMENT\n*STEP, NAME=FC_INITIALIZE, NLGEOM=YES\n*STATIC\n1., 1.\n*BOUNDARY\n");
    for n in model.nodes().keys() {
        let u = t.initial_displacement.get(n).copied().unwrap_or([0.; 3]);
        for (i, v) in u.iter().enumerate() {
            writeln!(out, "{n}, {}, {}, {v:.17e}", i + 1, i + 1).ok();
        }
    }
    if !step.loads.is_empty() {
        out.push_str("*CLOAD\n");
    }
    for name in &step.loads {
        for (n, f) in r.nodal_loads.get(name).ok_or_else(|| refuse(format!("load `{name}` not resolved")))? {
            for (i, v) in f.iter().enumerate() {
                if *v != 0. {
                    writeln!(out, "{n}, {}, {v:.17e}", i + 1).ok();
                }
            }
        }
    }
    out.push_str("*OUTPUT, FIELD, FREQUENCY=1\n*NODE OUTPUT, NSET=FC_NODES\nU, RF\n*END STEP\n");
    writeln!(out,"*STEP, NAME=FC_MOTION, NLGEOM=YES, AMPLITUDE=STEP, INC={}\n*DYNAMIC, DIRECT, NOHAF, APPLICATION=TRANSIENT FIDELITY, ALPHA={alpha:.17e}, BETA={beta:.17e}, GAMMA={gamma:.17e}\n{:.17e}, {:.17e}",t.increments,t.time_step,t.time_step*f64::from(t.increments)).ok();
    out.push_str("** Initial acceleration is recalculated; INITIAL=NO is deliberately absent.\n*BOUNDARY, OP=NEW\n");
    for (n, d) in fixed {
        writeln!(out, "{n}, {d}, {d}, 0.").ok();
    }
    out.push_str("** Constant CLOADs propagate unchanged from FC_INITIALIZE.\n");
    writeln!(out,"*OUTPUT, FIELD, FREQUENCY={}\n*NODE OUTPUT, NSET=FC_NODES\nU, V, A, RF\n*OUTPUT, HISTORY, FREQUENCY={}\n*ENERGY OUTPUT\nALLKE, ALLSE, ALLWK, ALLVD, ETOTAL\n*END STEP",t.output_every,t.output_every).ok();
    // SDV1..6 are element nodal resisting forces; user elements have no ODB stress field.
    let end = out.rfind("*END STEP").ok_or_else(|| refuse("missing dynamic end"))?;
    out.insert_str(end, &format!("*EL FILE, ELSET=FC_ELEMENTS, FREQUENCY={}\nSDV\n*ENERGY FILE, FREQUENCY={}\n", t.output_every, t.output_every));
    Ok(out)
}
