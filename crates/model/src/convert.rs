//! Every dimensional number of the model with its quantity: the model expressed in another unit
//! system ([`Model::converted`]: a deck in a chosen system, the millimetres the meshers work in)
//! and back ([`Model::from_converted`]).
//!
//! The visitors destructure every type field by field, so a new field does not compile until it
//! is given its quantity here (or said to be dimensionless).

use crate::units::Quantity as Q;
use crate::*;

/// Calls `f` with every dimensional number of a value and its quantity.
pub(crate) trait Quantities {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64));
}

fn each(values: &mut [f64], q: Q, f: &mut dyn FnMut(Q, &mut f64)) {
    for v in values {
        f(q, v);
    }
}

impl Quantities for Material {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let Material { elastic, hyperelastic, density, fluid, plastic, damage, erosion, thermal } = self;
        match elastic {
            Some(Elastic::Isotropic { young, poisson: _ }) => f(Q::Pressure, young),
            Some(Elastic::EngineeringConstants { young, poisson: _, shear }) => {
                each(young, Q::Pressure, f);
                each(shear, Q::Pressure, f);
            }
            None => {}
        }
        match hyperelastic {
            Some(Hyperelastic::NeoHookean {}) | None => {}
            Some(Hyperelastic::MooneyRivlin { c01 }) => f(Q::Pressure, c01),
            Some(Hyperelastic::Yeoh { c20, c30 }) => {
                f(Q::Pressure, c20);
                f(Q::Pressure, c30);
            }
        }
        if let Some(rho) = density {
            f(Q::Density, rho);
        }
        if let Some(Fluid::Newtonian { viscosity }) = fluid {
            f(Q::DynamicViscosity, viscosity);
        }
        match plastic {
            Some(Plastic::Bilinear { yield_stress, hardening_modulus, hardening: _ }) => {
                f(Q::Pressure, yield_stress);
                f(Q::Pressure, hardening_modulus);
            }
            Some(Plastic::Table { points }) => {
                for HardeningPoint { strain: _, stress } in points {
                    f(Q::Pressure, stress);
                }
            }
            Some(Plastic::JohnsonCook {
                a,
                b,
                n: _,
                c: _,
                m: _,
                reference_strain_rate,
                room_temperature,
                melt_temperature,
                specific_heat,
                taylor_quinney: _,
            }) => {
                f(Q::Pressure, a);
                f(Q::Pressure, b);
                f(Q::StrainRate, reference_strain_rate);
                f(Q::Temperature, room_temperature);
                f(Q::Temperature, melt_temperature);
                if let Some(c) = specific_heat {
                    f(Q::SpecificHeat, c);
                }
            }
            None => {}
        }
        // Johnson–Cook damage constants and the erosion strain are dimensionless.
        let _: (&mut Option<Damage>, &mut Option<Erosion>) = (damage, erosion);
        if let Some(Thermal { conductivity, specific_heat, expansion, reference_temperature }) = thermal {
            directional(conductivity, Q::ThermalConductivity, f);
            if let Some(c) = specific_heat {
                f(Q::SpecificHeat, c);
            }
            directional(expansion, Q::ThermalExpansion, f);
            f(Q::Temperature, reference_temperature);
        }
    }
}

fn directional(d: &mut Option<Directional>, q: Q, f: &mut dyn FnMut(Q, &mut f64)) {
    match d {
        Some(Directional::Isotropic(v)) => f(q, v),
        Some(Directional::Orthotropic(v)) => each(v, q, f),
        None => {}
    }
}

impl Quantities for Section {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let Section { set: _, material: _, kind, formulation: _ } = self;
        match kind {
            SectionKind::Solid {} | SectionKind::Axisymmetric {} => {}
            SectionKind::Shell { thickness, offset: _, points: _ } => f(Q::Length, thickness),
            SectionKind::Beam { profile, orientation: _ } => match profile {
                BeamProfile::Rectangle { width, height } => {
                    f(Q::Length, width);
                    f(Q::Length, height);
                }
                BeamProfile::Circle { radius } => f(Q::Length, radius),
                BeamProfile::Pipe { radius, thickness } => {
                    f(Q::Length, radius);
                    f(Q::Length, thickness);
                }
                BeamProfile::I { height, width, web, flange } => {
                    for v in [height, width, web, flange] {
                        f(Q::Length, v);
                    }
                }
                BeamProfile::Box { width, height, thickness } => {
                    for v in [width, height, thickness] {
                        f(Q::Length, v);
                    }
                }
            },
            SectionKind::Rod { area } => f(Q::Area, area),
            SectionKind::PlaneStress { thickness } | SectionKind::PlaneStrain { thickness } => f(Q::Length, thickness),
        }
    }
}

impl Quantities for Friction {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let Friction { static_coefficient: _, kinetic: _, decay } = self;
        f(Q::InverseVelocity, decay);
    }
}

impl Quantities for PointMass {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let PointMass { target: _, mass, inertia } = self;
        match mass {
            Some(Mass::Isotropic(m)) => f(Q::Mass, m),
            Some(Mass::Directional(m)) => each(m, Q::Mass, f),
            None => {}
        }
        if let Some(i) = inertia {
            each(i, Q::RotaryInertia, f);
        }
    }
}

impl Quantities for Spring {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        // Along degree of freedom 4 to 6 a spring is torsional.
        let torsional = |d: &Direction| matches!(d, Direction::Dof(4..=6));
        match self {
            Spring::Spring { first: _, second: _, direction, stiffness } => {
                f(if torsional(direction) { Q::RotationalStiffness } else { Q::Stiffness }, stiffness)
            }
            Spring::Dashpot { first: _, second: _, direction, coefficient } => {
                f(if torsional(direction) { Q::RotationalDamping } else { Q::Damping }, coefficient)
            }
        }
    }
}

impl Quantities for BoundaryCondition {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        match self {
            BoundaryCondition::Fixed { target: _ } | BoundaryCondition::Symmetry { target: _, plane: _ } => {}
            // Three displacements, then three rotations.
            BoundaryCondition::Displacement { target: _, values } => {
                for (i, v) in values.iter_mut().enumerate() {
                    if let Some(v) = v {
                        f(if i < 3 { Q::Length } else { Q::Angle }, v);
                    }
                }
            }
            BoundaryCondition::Temperature { target: _, value } => f(Q::Temperature, value),
        }
    }
}

impl Quantities for Load {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        match self {
            Load::Force { target: _, vector, distribution: _ } => each(vector, Q::Force, f),
            Load::Moment { target: _, vector, distribution: _ } => each(vector, Q::Moment, f),
            Load::Traction { set: _, vector } => each(vector, Q::Pressure, f),
            Load::Pressure { set: _, magnitude } => f(Q::Pressure, magnitude),
            Load::Gravity { acceleration, set: _ } => each(acceleration, Q::Acceleration, f),
            Load::BodyForce { set: _, vector } => each(vector, Q::ForcePerVolume, f),
            Load::Centrifugal { angular_velocity, point, axis: _, set: _ } => {
                f(Q::AngularVelocity, angular_velocity);
                each(point, Q::Length, f);
            }
            Load::SurfaceFlux { set: _, flux } => f(Q::HeatFlux, flux),
            Load::Convection { set: _, film, sink } => {
                f(Q::HeatTransferCoefficient, film);
                f(Q::Temperature, sink);
            }
            Load::HeatSource { set: _, power } => f(Q::HeatGeneration, power),
            Load::Temperature { target: _, value } => f(Q::Temperature, value),
        }
    }
}

fn face_condition(c: &mut FaceCondition, f: &mut dyn FnMut(Q, &mut f64)) {
    match c {
        FaceCondition::Wall {} | FaceCondition::Slip {} | FaceCondition::Periodic {} => {}
        FaceCondition::MovingWall { velocity } | FaceCondition::Velocity { velocity, profile: _ } => each(velocity, Q::Velocity, f),
        FaceCondition::Pressure { pressure, .. } => f(Q::Pressure, pressure),
    }
}

fn thermal_condition(c: &mut ThermalCondition, f: &mut dyn FnMut(Q, &mut f64)) {
    match c {
        ThermalCondition::Adiabatic {} => {}
        ThermalCondition::Temperature { value } => f(Q::Temperature, value),
        ThermalCondition::Flux { value } => f(Q::HeatFlux, value),
    }
}

impl Quantities for FluidDomain {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let FluidDomain {
            material: _,
            min,
            max,
            spacing,
            faces,
            obstacles: _,
            surfaces,
            initial_velocity,
            acceleration,
            tolerance,
            thermal,
            porous,
            fans,
        } = self;
        each(min, Q::Length, f);
        each(max, Q::Length, f);
        f(Q::Length, spacing);
        let BoxFaces { x_min, x_max, y_min, y_max, z_min, z_max } = faces;
        for c in [x_min, x_max, y_min, y_max, z_min, z_max] {
            face_condition(c, f);
        }
        for SurfaceBoundary { face: _, condition } in surfaces {
            match condition {
                SurfaceCondition::Wall {} => {}
                SurfaceCondition::MovingWall { velocity } => each(velocity, Q::Velocity, f),
            }
        }
        each(initial_velocity, Q::Velocity, f);
        each(acceleration, Q::Acceleration, f);
        if let Some(t) = tolerance {
            f(Q::Length, t);
        }
        if let Some(FluidThermal { initial_temperature, gravity, expansion, faces, surfaces, solids, contacts, temperature_zero }) = thermal {
            f(Q::Temperature, initial_temperature);
            f(Q::Temperature, temperature_zero);
            each(gravity, Q::Acceleration, f);
            if let Some(b) = expansion {
                f(Q::ThermalExpansion, b);
            }
            let BoxThermalFaces { x_min, x_max, y_min, y_max, z_min, z_max } = faces;
            for c in [x_min, x_max, y_min, y_max, z_min, z_max] {
                thermal_condition(c, f);
            }
            for SurfaceThermalBoundary { face: _, condition } in surfaces {
                thermal_condition(condition, f);
            }
            for ConductingSolid { body: _, material: _, power, limit } in solids {
                f(Q::Power, power);
                if let Some(t) = limit {
                    f(Q::Temperature, t);
                }
            }
            for ThermalContact { bodies: _, resistance, thickness, conductivity } in contacts {
                if let Some(r) = resistance {
                    f(Q::ThermalContactResistance, r);
                }
                if let Some(t) = thickness {
                    f(Q::Length, t);
                }
                if let Some(k) = conductivity {
                    f(Q::ThermalConductivity, k);
                }
            }
        }
        for PorousZone { name: _, min, max, axis: _, resistance } in porous {
            each(min, Q::Length, f);
            each(max, Q::Length, f);
            match resistance {
                Resistance::DarcyForchheimer { permeability, inertial } => {
                    f(Q::Area, permeability);
                    f(Q::InverseLength, inertial);
                }
                Resistance::LossCoefficient { coefficient: _ } => {}
            }
        }
        for Fan { name: _, min, max, direction: _, curve, relaxation: _ } in fans {
            each(min, Q::Length, f);
            each(max, Q::Length, f);
            for [flow, rise] in curve {
                f(Q::VolumeFlow, flow);
                f(Q::Pressure, rise);
            }
        }
    }
}

impl Quantities for Step {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let Step { name: _, procedure, output: _, boundary_conditions: _, loads: _ } = self;
        match procedure {
            Procedure::Static {} | Procedure::CoupledTemperatureDisplacement {} => {}
            Procedure::Frequency { modes: _, min_hz, max_hz } => {
                for v in [min_hz, max_hz].into_iter().flatten() {
                    f(Q::Frequency, v);
                }
            }
            // The shift is a load factor.
            Procedure::Buckle { modes: _, shift: _ } => {}
            Procedure::ExplicitDynamic(e) => {
                let ExplicitDynamic {
                    end_time,
                    dt_scale: _,
                    max_increments: _,
                    output_interval,
                    history_interval,
                    history_nodes: _,
                    initial_velocities,
                    rigid_walls,
                    hourglass: _,
                    bulk_viscosity: _,
                    nonlinear_geometry: _,
                    stress_rate: _,
                    mass_scaling,
                } = e;
                f(Q::Time, end_time);
                for v in [output_interval, history_interval].into_iter().flatten() {
                    f(Q::Time, v);
                }
                for InitialVelocity { set: _, velocity } in initial_velocities {
                    each(velocity, Q::Velocity, f);
                }
                for RigidWall { set: _, point, normal: _, friction } in rigid_walls {
                    each(point, Q::Length, f);
                    friction.quantities(f);
                }
                if let Some(MassScaling { target_increment, set: _ }) = mass_scaling {
                    f(Q::Time, target_increment);
                }
            }
            Procedure::Fluid { domain: _, end_time, steps: _, time_step, frames: _, magic: _, collision: _, turbulence: _, steady } => {
                for v in [end_time, time_step].into_iter().flatten() {
                    f(Q::Time, v);
                }
                if let Some(Steady { monitors: _, tolerance: _, window, min_time, divergence: _, averaging }) = steady {
                    for v in [window, min_time].into_iter().flatten() {
                        f(Q::Time, v);
                    }
                    if let Some(a) = averaging {
                        f(Q::Time, &mut a.start);
                    }
                }
            }
            Procedure::HeatTransfer { transient } => {
                if let Some(Transient { end_time, increment, theta: _, initial_temperature }) = transient {
                    f(Q::Time, end_time);
                    f(Q::Time, increment);
                    f(Q::Temperature, initial_temperature);
                }
            }
            Procedure::NonlinearStatic(n) => {
                // Increments are fractions of the step, tolerances relative, the load path load
                // factors.
                let NonlinearStatic {
                    nonlinear_geometry: _,
                    initial_increment: _,
                    min_increment: _,
                    max_increment: _,
                    max_increments: _,
                    max_iterations: _,
                    tolerances: _,
                    line_search: _,
                    load_path: _,
                    imperfection,
                    arc_length,
                } = n;
                // The buckling mode, scaled to a largest component of one, times `scale`.
                if let Some(Imperfection { step: _, mode: _, scale }) = imperfection {
                    f(Q::Length, scale);
                }
                if let Some(ArcLength {
                    constraint: _,
                    max_load_factor: _,
                    displacement_limit: Some(DisplacementLimit { node: _, direction, displacement }),
                }) = arc_length
                {
                    f(if *direction >= 4 { Q::Angle } else { Q::Length }, displacement);
                }
            }
        }
    }
}

impl Quantities for Datum {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        match self {
            Datum::Plane { point, normal: _ } | Datum::Axis { point, direction: _ } | Datum::Point { point } => each(point, Q::Length, f),
        }
    }
}

impl Quantities for Seed {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let Seed { edge: _, divisions, bias: _, constraint: _ } = self;
        match divisions {
            Divisions::Count(_) => {}
            Divisions::Size(h) => f(Q::Length, h),
        }
    }
}

impl Quantities for SizeControl {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let SizeControl { region, size } = self;
        f(Q::Length, size);
        match region {
            SizeRegion::Geometry(_) => {}
            SizeRegion::Box { min, max } => {
                each(min, Q::Length, f);
                each(max, Q::Length, f);
            }
            SizeRegion::Sphere { center, radius } => {
                each(center, Q::Length, f);
                f(Q::Length, radius);
            }
        }
    }
}

impl Quantities for Relaxation {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let Relaxation { bodies: _, cells, reapply: _ } = self;
        for RelaxationCell { kind: _, at, moved } in cells {
            for p in at.iter_mut().chain(moved.iter_mut()) {
                each(p, Q::Length, f);
            }
        }
    }
}

impl Quantities for MeshPlan {
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        let MeshPlan {
            body: _,
            size,
            tolerance,
            layers: _,
            mesher: _,
            direction: _,
            directions: _,
            features: _,
            cuts,
            rings,
            copies: _,
            unmapped: _,
            notes: _,
        } = self;
        f(Q::Length, size);
        if let Some(t) = tolerance {
            f(Q::Length, t);
        }
        for PlanCut { point, normal: _, role: _, feature: _, reason: _ } in cuts {
            each(point, Q::Length, f);
        }
        for PlanRing { feature: _, kind: _, center, axis: _, radius, frame, faces: _ } in rings {
            each(center, Q::Length, f);
            f(Q::Length, radius);
            f(Q::Length, frame);
        }
    }
}

impl Model {
    /// Every dimensional number of the model, mesh included, with its quantity.
    fn quantities(&mut self, f: &mut dyn FnMut(Q, &mut f64)) {
        // Destructured so that a new collection does not compile until it is listed here.
        let Model {
            nodes,
            elements: _,
            classification: _,
            display_units: _,
            sets: _,
            materials,
            sections,
            reference_points,
            couplings: _,
            contacts,
            ties,
            masses,
            springs,
            boundary_conditions,
            loads,
            fluid_domains,
            steps,
            datums,
            partitions: _,
            suppressed_datums: _,
            seeds,
            size_controls,
            size_growth: _,
            relaxations,
            mesh_plans,
            face_merges: _,
            stitches,
            face_splits: _,
            defeatures: _,
        } = self;
        for n in nodes.values_mut() {
            each(&mut n.xyz, Q::Length, f);
        }
        materials.values_mut().for_each(|x| x.quantities(f));
        sections.values_mut().for_each(|x| x.quantities(f));
        for rp in reference_points.values_mut() {
            each(&mut rp.xyz, Q::Length, f);
        }
        contacts.values_mut().for_each(|c| c.friction.quantities(f));
        for t in ties.values_mut() {
            if let Some(p) = &mut t.position_tolerance {
                f(Q::Length, p);
            }
        }
        masses.values_mut().for_each(|x| x.quantities(f));
        springs.values_mut().for_each(|x| x.quantities(f));
        boundary_conditions.values_mut().for_each(|x| x.quantities(f));
        loads.values_mut().for_each(|x| x.quantities(f));
        fluid_domains.values_mut().for_each(|x| x.quantities(f));
        steps.iter_mut().for_each(|x| x.quantities(f));
        datums.values_mut().for_each(|x| x.quantities(f));
        seeds.iter_mut().for_each(|x| x.quantities(f));
        size_controls.values_mut().for_each(|x| x.quantities(f));
        relaxations.values_mut().for_each(|x| x.quantities(f));
        mesh_plans.values_mut().for_each(|x| x.quantities(f));
        for s in stitches {
            f(Q::Length, &mut s.tolerance);
        }
    }

    /// The model with every dimensional number `v` of quantity `q` replaced by `convert(q, v)`:
    /// stored SI expressed in another system, or such a copy brought back to SI. Fails when a
    /// converted number is not finite.
    pub fn converted(&self, convert: &dyn Fn(Q, f64) -> Option<f64>) -> Result<Model> {
        let mut out = self.clone();
        let mut bad = None;
        out.quantities(&mut |q, v| match convert(q, *v).filter(|x| x.is_finite()) {
            Some(x) => *v = x,
            None => {
                bad.get_or_insert((q, *v));
            }
        });
        match bad {
            None => Ok(out),
            Some((q, v)) => Err(Error::Invalid(format!("the {} {v} cannot be expressed in that unit system", q.name().replace('_', " ")))),
        }
    }

    /// Bring `self`, a copy of `base` expressed in another system ([`Model::converted`] of `base`
    /// as `converted_base`) and then changed, back with `to_base`: what is unchanged since the
    /// conversion keeps `base`'s exact value (a round trip through a unit can move a value by an
    /// ulp), the rest is converted.
    pub fn from_converted(&self, converted_base: &Model, base: &Model, to_base: &dyn Fn(Q, f64) -> Option<f64>) -> Result<Model> {
        let mut out = self.converted(to_base)?;
        // Nodes, one by one.
        for (label, node) in &self.nodes {
            if converted_base.nodes.get(label) == Some(node)
                && let (Some(kept), Some(slot)) = (base.nodes.get(label), out.nodes.get_mut(label))
            {
                *slot = *kept;
            }
        }
        // Definitions, whole: one unchanged keeps its value.
        fn keep_map<T: Clone + PartialEq>(
            now: &BTreeMap<String, T>,
            was: &BTreeMap<String, T>,
            base: &BTreeMap<String, T>,
            out: &mut BTreeMap<String, T>,
        ) {
            for (name, item) in now {
                if was.get(name) == Some(item)
                    && let (Some(kept), Some(slot)) = (base.get(name), out.get_mut(name))
                {
                    *slot = kept.clone();
                }
            }
        }
        fn keep_list<T: Clone + PartialEq>(now: &[T], was: &[T], base: &[T], out: &mut [T]) {
            for (i, item) in now.iter().enumerate() {
                if was.get(i) == Some(item)
                    && let (Some(kept), Some(slot)) = (base.get(i), out.get_mut(i))
                {
                    *slot = kept.clone();
                }
            }
        }
        let (c, b) = (converted_base, base);
        keep_map(&self.materials, &c.materials, &b.materials, &mut out.materials);
        keep_map(&self.sections, &c.sections, &b.sections, &mut out.sections);
        keep_map(&self.reference_points, &c.reference_points, &b.reference_points, &mut out.reference_points);
        keep_map(&self.contacts, &c.contacts, &b.contacts, &mut out.contacts);
        keep_map(&self.ties, &c.ties, &b.ties, &mut out.ties);
        keep_map(&self.masses, &c.masses, &b.masses, &mut out.masses);
        keep_map(&self.springs, &c.springs, &b.springs, &mut out.springs);
        keep_map(&self.boundary_conditions, &c.boundary_conditions, &b.boundary_conditions, &mut out.boundary_conditions);
        keep_map(&self.loads, &c.loads, &b.loads, &mut out.loads);
        keep_map(&self.fluid_domains, &c.fluid_domains, &b.fluid_domains, &mut out.fluid_domains);
        keep_list(&self.steps, &c.steps, &b.steps, &mut out.steps);
        keep_map(&self.datums, &c.datums, &b.datums, &mut out.datums);
        keep_list(&self.seeds, &c.seeds, &b.seeds, &mut out.seeds);
        keep_map(&self.size_controls, &c.size_controls, &b.size_controls, &mut out.size_controls);
        keep_map(&self.relaxations, &c.relaxations, &b.relaxations, &mut out.relaxations);
        keep_map(&self.mesh_plans, &c.mesh_plans, &b.mesh_plans, &mut out.mesh_plans);
        keep_list(&self.stitches, &c.stitches, &b.stitches, &mut out.stitches);
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_keep_unchanged_values_exactly() {
        let mut m = Model::default();
        // 0.1 + 0.2 has 17 significant digits: through millimetres and back it may move by an ulp.
        for (label, x) in [(1, 0.1 + 0.2), (2, 3.97e-3), (3, 1.0 / 3.0)] {
            m.set_node(label, [x, -x, 0.0]).unwrap();
        }
        let to_mm = |q: Q, v: f64| (q == Q::Length).then_some(v * 1e3).or(Some(v));
        let to_m = |q: Q, v: f64| (q == Q::Length).then(|| v / 1e3).or(Some(v));
        let mm = m.converted(&to_mm).unwrap();
        assert!((mm.node(2).unwrap().xyz[0] - 3.97).abs() < 1e-14);
        let mut moved = mm.clone();
        moved.set_node(3, [5.0, 0.0, 0.0]).unwrap();
        let back = moved.from_converted(&mm, &m, &to_m).unwrap();
        assert_eq!(back.node(1), m.node(1));
        assert_eq!(back.node(2), m.node(2));
        assert_eq!(back.node(3).unwrap().xyz, [0.005, 0.0, 0.0]);
        assert!(m.converted(&|_, v| Some(v * f64::MAX * 4.0)).is_err());
    }
}
