//! What an explicit dynamic run keeps: histories sampled at the history interval and field
//! frames at the output interval.

use serde::Serialize;

use crate::elements::Failure;

/// The energy account at one time (Belytschko, Liu, Moran, Elkhodary 2014, §6.2.3). The works of
/// the internal and external forces are integrated with the trapezoidal rule over each
/// increment's displacement; the constraints' with their impulses (see [`crate::integrate`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, serde::Deserialize)]
pub struct Energies {
    /// ½ vᵀ M v.
    pub kinetic: f64,
    /// Work of the material stresses: the strain energy for an elastic material.
    pub internal: f64,
    /// Work of the hourglass forces (stored in the stiffness form, dissipated in the viscous one).
    pub hourglass: f64,
    /// Work of the bulk viscosity (dissipated).
    pub viscous: f64,
    /// Work of the loads.
    pub external: f64,
    /// Work of the boundary conditions and rigid walls on the model.
    pub constraint: f64,
    /// Energy taken by contact: the work the model does against the contact forces (the penalty
    /// energy stored in the contact springs, plus what friction dissipates there).
    pub contact: f64,
    /// The kinetic energy of the mass added by mass scaling (part of `kinetic`).
    pub added_kinetic: f64,
    /// Friction dissipation, in contact and on rigid walls (part of `contact` and of
    /// `constraint`).
    pub friction: f64,
    /// Plastic work: the part of `internal` dissipated by plastic flow.
    pub plastic: f64,
    /// Internal energy of the eroded elements when they eroded (part of `internal`).
    pub eroded: f64,
    /// kinetic + internal + hourglass + viscous + contact − external − constraint − the initial
    /// kinetic energy: zero when energy is conserved.
    pub balance: f64,
}

impl Energies {
    /// The balance relative to the largest energy in the account.
    pub fn relative_error(&self, initial_kinetic: f64) -> f64 {
        let scale = [
            initial_kinetic,
            self.kinetic,
            self.internal + self.hourglass + self.viscous,
            self.external.abs(),
            self.constraint.abs(),
            self.contact.abs(),
        ]
        .into_iter()
        .fold(0.0, f64::max);
        if scale > 0.0 { self.balance.abs() / scale } else { 0.0 }
    }
}

/// The fields at one time. Vectors that the step's output does not request are empty.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Frame {
    pub time: f64,
    /// Per node in [`ExplicitResult::nodes`] order.
    pub displacement: Vec<[f64; 3]>,
    pub velocity: Vec<[f64; 3]>,
    /// Per element in [`ExplicitResult::elements`] order (one-point elements: one value each):
    /// xx, yy, zz, xy, xz, yz, engineering shear strains.
    pub stress: Vec<[f64; 6]>,
    /// The accumulated strain increments (small strain: the strain; large deformation: the
    /// integral of the rate of deformation, rotated with the material).
    pub strain: Vec<[f64; 6]>,
    /// Per element, with stress output in a model with plasticity: equivalent plastic strain,
    /// damage and temperature (NaN for an element whose material has no temperature).
    pub plastic_strain: Vec<f64>,
    pub damage: Vec<f64>,
    pub temperature: Vec<f64>,
    /// Per element, with stress output: whether it has eroded.
    pub failed: Vec<bool>,
}

/// An element removed from the analysis.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
pub struct ErodedElement {
    pub label: u32,
    pub time: f64,
    pub failure: Failure,
    /// Its internal energy (stress work) when it eroded.
    pub internal_energy: f64,
}

/// Quantities sampled through the run.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct History {
    pub time: Vec<f64>,
    /// The time increment taken at each sample.
    pub increment: Vec<f64>,
    pub energies: Vec<Energies>,
    /// Total momentum M v.
    pub momentum: Vec<[f64; 3]>,
    /// Sum of the forces the boundary conditions and rigid walls exert on the model.
    pub constraint_force: Vec<[f64; 3]>,
    /// Per contact (in [`ExplicitResult::contacts`] order): the normal force across it (the sum
    /// of its constraints' normal forces).
    pub contact_force: Vec<Vec<f64>>,
    /// The nodes with histories, and per sample their displacement and velocity (in this order).
    pub nodes: Vec<u32>,
    pub displacement: Vec<Vec<[f64; 3]>>,
    pub velocity: Vec<Vec<[f64; 3]>>,
}

/// What mass scaling added to a run.
#[derive(Debug, Clone, Default, PartialEq, Serialize, serde::Deserialize)]
pub struct MassScalingSummary {
    pub target_increment: f64,
    /// Elements whose density was raised.
    pub elements: usize,
    pub added_mass: f64,
    /// The added mass over the model's mass (with it).
    pub fraction: f64,
    /// The added mass per section, by name.
    pub sections: Vec<(String, f64)>,
}

/// What a contact did over the run.
#[derive(Debug, Clone, Default, PartialEq, Serialize, serde::Deserialize)]
pub struct ContactSummary {
    pub name: String,
    /// Nodes moved onto the other surface before the analysis.
    pub adjusted: Vec<crate::contact::Adjustment>,
    /// The largest normal force across it, and when.
    pub peak_force: f64,
    pub peak_time: f64,
    /// The first and the last time it carried a force (none when it never did).
    pub first_contact: Option<f64>,
    pub last_contact: Option<f64>,
    /// The deepest penetration of a node into a face.
    pub max_penetration: f64,
    /// Friction dissipated in it.
    pub friction: f64,
}

/// An element that limits the stable time increment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, serde::Deserialize)]
pub struct ControllingElement {
    pub label: u32,
    /// Its critical increment at rest (no bulk viscosity), before the scale factor.
    pub increment: f64,
}

/// The result of an explicit dynamic step.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExplicitResult {
    pub step: String,
    /// Node labels with degrees of freedom, sorted.
    pub nodes: Vec<u32>,
    /// Element labels, sorted.
    pub elements: Vec<u32>,
    pub frames: Vec<Frame>,
    pub history: History,
    /// The elements with the smallest critical increments (at most ten), smallest first.
    pub controlling: Vec<ControllingElement>,
    /// The first time increment taken (the scale times the smallest critical increment).
    pub initial_increment: f64,
    pub increments: u64,
    /// The time reached: the step's end time unless the run stopped at the increment cap.
    pub time: f64,
    pub completed: bool,
    pub total_mass: f64,
    /// The energy account at the end.
    pub energy: Energies,
    pub initial_kinetic: f64,
    /// Elements eroded during the run, in the order they eroded.
    pub eroded: Vec<ErodedElement>,
    /// Per contact of the model, in name order.
    pub contacts: Vec<ContactSummary>,
    /// What mass scaling added, when the step asks for it.
    pub mass_scaling: Option<MassScalingSummary>,
    pub warnings: Vec<String>,
}
