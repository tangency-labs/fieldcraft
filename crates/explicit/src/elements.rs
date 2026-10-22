//! Element internal forces, the material update and critical time increments.
//!
//! **Kinematics.** With nonlinear geometry the elements are updated Lagrangian (T. Belytschko,
//! W. K. Liu, B. Moran, K. Elkhodary, *Nonlinear Finite Elements for Continua and Structures*,
//! 2nd ed., Wiley 2014, ch. 4 and 6): the strain and spin increments of an increment come from
//! the displacement increment Δu and the gradients of the **midpoint** configuration
//! x_{n−½} = X + u_n − ½Δu, ΔG = (1/V) Σ_I Δu_I ⊗ b_I(x_{n−½}), which makes the update
//! incrementally objective (T. J. R. Hughes, J. Winget, IJNME 15(12):1862–1867, 1980); the nodal
//! forces of the Cauchy stress use the gradients of the current configuration,
//! f_iI = σ_ij b_jI(x_n). Without it, both are the gradients of the undeformed mesh and the
//! update is small-strain elasticity (and plasticity).
//!
//! **Stress update** in rate form: the material state (stress, back stress, the accumulated
//! strain and the stiffness-form hourglass forces) is rotated with the increment's rotation —
//! the Cayley transform of the spin increment (Jaumann rate, Hughes & Winget 1980) or the
//! increment of the polar rotation (Green–Naghdi rate, D. P. Flanagan, L. M. Taylor, CMAME
//! 62(3):305–320, 1987), see [`fieldcraft_fem::rotation`] — then the constitutive increment is
//! added: hypoelastic, or J2 radial return with Johnson–Cook damage
//! ([`fieldcraft_fem::plasticity`]). Orthotropic elasticity is updated in the polar frame (its
//! axes turn with the material) whichever rate the step names.
//!
//! **One-point hexahedron** (D. P. Flanagan, T. Belytschko, IJNME 17(5):679–706, 1981): the
//! hourglass modes are resisted through the shape vectors γ_α of the current configuration:
//! with the hourglass rates q̇_iα = Σ_I v_iI γ_αI, the generalized forces are Q_iα = C q̇_iα
//! (viscous) or accumulate as Q_iα += k Σ_I Δu_iI γ_αI (stiffness, rotated with the element),
//! and give the nodal forces f_iI = Σ_α Q_iα γ_αI. Because the γ_α are orthogonal to linear
//! fields, these forces do no work in a constant strain state. The coefficients are scaled from
//! the dimensionless hourglass coefficient κ with the element's current size
//! h = V / √(2 λ_max(Σ b bᵀ)) (the edge, for a cube), density ρ and wave speed c:
//! C = κ ρ c h² / 16 damps a cube's hourglass modes at a fraction κ of critical at the frequency
//! 2c/h, and k = κ ρ c² h / 16 puts their frequency at √κ · 2c/h. (These normalizations are ours;
//! the paper leaves the scale of the coefficient to the user.)
//!
//! **Average nodal pressure tetrahedron** (J. Bonet, A. J. Burton, Commun. Numer. Methods Eng.
//! 14(5):437–449, 1998), in rate form: the volumetric strain increment of each node is the
//! volume-weighted mean of its tetrahedra's, and each tetrahedron takes the mean of its nodes' as
//! its own, so its pressure increment is K times that mean while the deviator is its own. It does
//! not lock when the material is nearly incompressible or flows plastically.
//!
//! **Bulk viscosity** (J. von Neumann, R. D. Richtmyer, J. Appl. Phys. 21(3):232–237, 1950; M. L.
//! Wilkins, J. Comput. Phys. 36(3):281–303, 1980): in compression (volumetric strain rate
//! ε̇ᵥ < 0) a pressure q = ρ L (C_q L ε̇ᵥ² − C_l c ε̇ᵥ) spreads shocks over a few elements and
//! damps the ringing behind wave fronts. It lowers the stable increment to
//! Δt = L / (Q + √(Q² + c²)), Q = C_l c + C_q L |ε̇ᵥ| (Wilkins 1980).
//!
//! **Critical increment.** Each element's stable length L = 2c/ω_max comes from the bound on its
//! highest frequency in its current shape ([`fieldcraft_fem::uniform`], with the element's lumped
//! masses, plus the hourglass stiffness when it is used), so without bulk viscosity its critical
//! increment is 2/ω_max (Belytschko et al. 2014, §6.3). The elastic moduli bound the tangent
//! moduli of plastic flow, so the bound holds for the inelastic materials too.
//!
//! **Erosion** (G. R. Johnson, R. A. Stryk, Int. J. Impact Eng. 5:411–421, 1987; T. Belytschko,
//! J. I. Lin, ibid. 5:111–127, 1987): an element whose material point fails (damage 1) or
//! reaches its material's erosion plastic strain is removed: it carries no stress from then on,
//! its nodes keep their mass and momentum, and its internal energy at that moment is reported.

use fieldcraft_fem::constitutive::MaterialModel;
use fieldcraft_fem::material::stress;
use fieldcraft_fem::plasticity::{Increment, PointState};
use fieldcraft_fem::rotation::{self, IDENTITY, Mat3, Polar};
use fieldcraft_fem::uniform::{Hex8, Tet4, UniformStrain};
use fieldcraft_model::{BulkViscosity, HourglassKind, StressRate};

use crate::system::{Material, NODE_DOFS, System};
use crate::{Error, Result};

/// Scratch vectors of one evaluation.
#[derive(Debug, Clone)]
pub struct Workspace {
    /// Internal forces per degree of freedom, by origin.
    pub stress_force: Vec<f64>,
    pub hourglass_force: Vec<f64>,
    pub viscous_force: Vec<f64>,
    /// Contact forces on the model per degree of freedom (empty without contact).
    pub contact_force: Vec<f64>,
    /// Per node: the volume-weighted volumetric strain increments of its tetrahedra, and their
    /// volume.
    theta: Vec<f64>,
    tet_volume: Vec<f64>,
    /// Per node, for hyperelastic tetrahedra: the current and the rest volume of its
    /// tetrahedra, then their ratio.
    tet_ratio: Vec<f64>,
    tet_rest: Vec<f64>,
    /// Per tetrahedron (in [`System::tets`] order): its increment between the two passes.
    tet_increment: Vec<Option<TetIncrement>>,
}

impl Workspace {
    pub fn new(sys: &System) -> Workspace {
        let n = sys.labels.len();
        let dofs = sys.dofs();
        let hyper_tets = sys.tets.iter().any(|t| sys.materials.get(t.material).is_some_and(|m| m.hyperelastic.is_some()));
        Workspace {
            stress_force: vec![0.0; dofs],
            hourglass_force: vec![0.0; dofs],
            viscous_force: vec![0.0; dofs],
            contact_force: if sys.contacts.is_empty() { Vec::new() } else { vec![0.0; dofs] },
            theta: vec![0.0; n],
            tet_volume: vec![0.0; n],
            tet_ratio: if hyper_tets { vec![1.0; n] } else { Vec::new() },
            tet_rest: if hyper_tets { vec![0.0; n] } else { Vec::new() },
            tet_increment: vec![None; sys.tets.len()],
        }
    }
}

/// Why an element was eroded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Failure {
    /// Its damage reached 1.
    Damage,
    /// Its equivalent plastic strain reached the material's erosion strain.
    PlasticStrain,
}

/// The state of one element (its single material point).
#[derive(Debug, Clone)]
pub struct ElementState {
    pub point: PointState,
    /// The accumulated (rotated) strain increments: engineering shears.
    pub strain: [f64; 6],
    /// The polar rotation and stretch (Green–Naghdi rate and orthotropic materials).
    pub polar: Polar,
    /// Stiffness-form hourglass generalized forces Q_α (one-point hexahedra).
    pub hourglass: [[f64; 3]; 4],
    /// Work of the stress on the element (its internal energy, hourglass apart).
    pub energy: f64,
    /// Plastic work on the element.
    pub plastic_work: f64,
    pub failed: Option<Failure>,
    /// Whether its material follows a temperature (Johnson–Cook): the others have none to
    /// report (their point's temperature is a placeholder the flow law never reads).
    pub has_temperature: bool,
}

/// The element states and the erosions of a run.
#[derive(Debug, Clone)]
pub struct State {
    /// Per element, in [`System::elements`] order.
    pub elements: Vec<ElementState>,
    /// Elements eroded during the latest evaluation: (index, failure).
    pub eroded: Vec<(usize, Failure)>,
}

impl State {
    pub fn new(sys: &System) -> State {
        let mut material = vec![usize::MAX; sys.elements.len()];
        for (index, m) in sys.hexes.iter().map(|h| (h.index, h.material)).chain(sys.tets.iter().map(|t| (t.index, t.material))) {
            if let Some(slot) = material.get_mut(index) {
                *slot = m;
            }
        }
        let elements = material
            .into_iter()
            .map(|m| {
                let j2 = sys.materials.get(m).and_then(|m| m.j2.as_ref());
                let temperature = j2.map_or(0.0, |j| j.initial_temperature());
                ElementState {
                    has_temperature: j2.is_some_and(|j| j.is_temperature_dependent()),
                    point: PointState::new(temperature),
                    strain: [0.0; 6],
                    polar: Polar::default(),
                    hourglass: [[0.0; 3]; 4],
                    energy: 0.0,
                    plastic_work: 0.0,
                    failed: None,
                }
            })
            .collect();
        State { elements, eroded: Vec::new() }
    }
}

pub(crate) fn gather<const N: usize>(x: &[f64], nodes: &[usize; N]) -> [[f64; 3]; N] {
    nodes.map(|n| [0, 1, 2].map(|k| x.get(NODE_DOFS * n + k).copied().unwrap_or(0.0)))
}

/// Add a node's vector into a per-degree-of-freedom vector.
pub(crate) fn add_at(f: &mut [f64], node: usize, v: &[f64; 3]) {
    for (k, c) in v.iter().enumerate() {
        if let Some(x) = f.get_mut(NODE_DOFS * node + k) {
            *x += c;
        }
    }
}

fn scatter<const N: usize>(f: &mut [f64], nodes: &[usize; N], v: &[[f64; 3]; N]) {
    for (n, vi) in nodes.iter().zip(v) {
        add_at(f, *n, vi);
    }
}

/// The bulk viscosity pressure (positive in compression) for a volumetric strain rate.
pub fn bulk_pressure(rate: f64, density: f64, length: f64, wave_speed: f64, b: &BulkViscosity) -> f64 {
    if rate < 0.0 { density * length * (b.quadratic * length * rate * rate - b.linear * wave_speed * rate) } else { 0.0 }
}

/// The critical increment of an element of stable length `length` at a volumetric strain rate
/// (Wilkins 1980).
pub fn critical_increment(rate: f64, length: f64, wave_speed: f64, b: &BulkViscosity) -> f64 {
    let q = if rate < 0.0 { b.linear * wave_speed + b.quadratic * length * rate.abs() } else { 0.0 };
    length / (q + (q * q + wave_speed * wave_speed).sqrt())
}

/// The largest absolute row sum of the hourglass vectors' Gram matrix Σ_I γ_αI γ_βI / m_I (an
/// upper bound of its largest eigenvalue, exact for parallelepipeds where the γ_α are
/// orthogonal).
pub fn hourglass_gram_bound(gamma: &[[f64; 8]; 4], mass: &[f64; 8]) -> f64 {
    let mut worst: f64 = 0.0;
    for a in gamma {
        let row: f64 = gamma.iter().map(|b| a.iter().zip(b).zip(mass).map(|((x, y), m)| x * y / m).sum::<f64>().abs()).sum();
        worst = worst.max(row);
    }
    worst
}

/// What one evaluation of the internal forces is given.
pub struct Motion<'a> {
    /// Displacements u^n.
    pub u: &'a [f64],
    /// The displacement increment that led to `u`, Δt^{n−½} v^{n−½}.
    pub du: &'a [f64],
    /// The velocity v^{n−½}.
    pub v: &'a [f64],
    /// The increment's length Δt^{n−½} (0 at the start).
    pub dt: f64,
    pub time: f64,
    pub bulk: &'a BulkViscosity,
}

/// The element's nodal values: rest positions, displacements, displacement increments and
/// velocities, gathered once per evaluation.
struct Nodal<const N: usize> {
    x: [[f64; 3]; N],
    u: [[f64; 3]; N],
    du: [[f64; 3]; N],
    v: [[f64; 3]; N],
}

impl<const N: usize> Nodal<N> {
    fn gather(sys: &System, nodes: &[usize; N], m: &Motion) -> Nodal<N> {
        Nodal { x: nodes.map(|n| sys.xyz.get(n).copied().unwrap_or_default()), u: gather(m.u, nodes), du: gather(m.du, nodes), v: gather(m.v, nodes) }
    }

    /// The positions X + u − s·Δu.
    fn positions(&self, s: f64) -> [[f64; 3]; N] {
        std::array::from_fn(|i| [0, 1, 2].map(|k| self.x[i][k] + self.u[i][k] - s * self.du[i][k]))
    }
}

fn inverted(label: u32, time: f64) -> Error {
    Error::Unstable(format!("element {label} turned inside out at time {time:.6e} (its volume is no longer positive)"))
}

/// The rotation of an increment for an element's state, advancing its polar decomposition when
/// it is tracked.
fn rotation_of(sys: &System, mat: &Material, es: &mut ElementState, g: &Mat3, label: u32, time: f64) -> Result<Mat3> {
    if !sys.nonlinear {
        return Ok(IDENTITY);
    }
    if sys.stress_rate == StressRate::GreenNaghdi || mat.j2.is_none() {
        es.polar.advance(g).ok_or_else(|| inverted(label, time))
    } else {
        let (_, w) = rotation::split(g);
        rotation::cayley(&w).ok_or_else(|| inverted(label, time))
    }
}

/// Rotate an element's state with `q`, then add the constitutive increment for the strain
/// increment `d` (engineering shears) with volumetric part `vol` (the tetrahedron's averaged
/// one; the trace of `d` for the hexahedron), track its energy at the midpoint volume `volume`,
/// and check for erosion.
#[allow(clippy::too_many_arguments)]
fn update(mat: &Material, es: &mut ElementState, q: &Mat3, d: &[f64; 6], vol: f64, dt: f64, density: f64, volume: f64) {
    let p = &mut es.point;
    if *q != IDENTITY {
        // Only the parts of the state that are in use (a back stress with kinematic hardening,
        // hourglass forces with stiffness-form control).
        p.stress = rotation::rotate(q, &p.stress);
        if p.back != [0.0; 6] {
            p.back = rotation::rotate(q, &p.back);
        }
        es.strain = rotation::rotate_strain(q, &es.strain);
        if es.hourglass != [[0.0; 3]; 4] {
            for h in &mut es.hourglass {
                *h = rotation::apply(q, h);
            }
        }
    }
    let before = p.stress;
    let inc = match &mat.j2 {
        Some(j2) => j2.update(p, d, vol, dt, density),
        None => {
            // Orthotropic: in the polar frame, σ̂ += D d̂.
            let r = es.polar.rotation;
            if r == IDENTITY {
                for (s, x) in p.stress.iter_mut().zip(stress(&mat.elasticity, d)) {
                    *s += x;
                }
            } else {
                let local = rotation::rotate_strain(&rotation::transpose(&r), d);
                for (s, x) in p.stress.iter_mut().zip(rotation::rotate(&r, &stress(&mat.elasticity, &local))) {
                    *s += x;
                }
            }
            Increment::default()
        }
    };
    // Stress work over the increment, ½(σ_n + σ_{n+1}) : Δε, with the volumetric strain the
    // stress took.
    let shift = (d[0] + d[1] + d[2] - vol) / 3.0;
    let mut work = 0.0;
    for k in 0..6 {
        let dk = if k < 3 { d[k] - shift } else { d[k] };
        work += 0.5 * (before[k] + p.stress[k]) * dk;
        es.strain[k] += dk;
    }
    es.energy += work * volume;
    es.plastic_work += inc.plastic_work * volume;
    if es.failed.is_none() {
        if p.damage >= 1.0 {
            es.failed = Some(Failure::Damage);
        } else if mat.erosion.is_some_and(|e| p.plastic_strain >= e) {
            es.failed = Some(Failure::PlasticStrain);
        }
    }
}

/// The hyperelastic update of an element at deformation gradient `f` (the hexahedron's mean,
/// or the tetrahedron's with its averaged volume ratio): the material point's Kirchhoff stress τ
/// gives the stress of the nodal forces τ/`j` (`j` the element's own volume ratio, so the forces
/// are ∫ τ ∇ₓN dV₀), the logarithmic strain is kept for output, and the stress work over the
/// increment `d` (on the midpoint volume `volume`) is added to the element's energy. Returns the
/// ratio of the tangent's largest normal stiffness (with the stress) to the material's
/// small-strain one, at least 1, which scales the element's frequency bound.
#[allow(clippy::too_many_arguments)]
fn hyperelastic(
    model: &MaterialModel,
    mat: &Material,
    es: &mut ElementState,
    f: &Mat3,
    j: f64,
    pressure: Option<f64>,
    d: &[f64; 6],
    volume: f64,
) -> Result<f64> {
    let unstable = |e: fieldcraft_fem::Error| Error::Unstable(e.to_string());
    if !(j > 0.0 && j.is_finite()) {
        return Err(Error::Unstable(format!("the material is inverted (volume ratio {j:.3e})")));
    }
    let r = model.response(f, None, &[], &mut []).map_err(unstable)?;
    let before = es.point.stress;
    es.point.stress = match pressure {
        // The isochoric part of the element's own F with the averaged pressure.
        Some(p) => {
            let iso = model.isochoric(f, None, &[], &mut []).map_err(unstable)?;
            let mut s = iso.stress.map(|t| t / j);
            for v in s.iter_mut().take(3) {
                *v += p;
            }
            s
        }
        None => r.stress.map(|t| t / j),
    };
    let work: f64 = (0..6).map(|k| 0.5 * (before[k] + es.point.stress[k]) * d[k]).sum();
    es.energy += work * volume;
    if let Some(e) = fieldcraft_fem::constitutive::logarithmic_strain(f) {
        es.strain = e;
    }
    let normal = (0..3).map(|k| r.tangent[k][k]).fold(0.0, f64::max) / j;
    let stress = es.point.stress.iter().fold(0.0f64, |m, s| m.max(s.abs()));
    Ok(((normal + stress) / mat.stiffest).max(1.0))
}

/// The mechanical part of a strain increment: less the thermal strain increment `thermal` in
/// every normal direction (isotropic expansion, the same in any frame).
fn mechanical(mut d: [f64; 6], thermal: f64) -> [f64; 6] {
    if thermal != 0.0 {
        for x in d.iter_mut().take(3) {
            *x -= thermal;
        }
    }
    d
}

/// The density and dilatational wave speed of an element of rest volume `rest`, now of volume
/// `volume`.
fn current(mat: &Material, rest: f64, volume: f64) -> (f64, f64) {
    let rho = mat.density * rest / volume;
    (rho, (mat.stiffest / rho).sqrt())
}

/// The bound ω² ≤ k λ_max(Σ b bᵀ) / (V m_min) on an element's highest frequency
/// ([`fieldcraft_fem::uniform`]), with the lumped masses of the undeformed element and
/// `gradient_norm` = λ_max(Σ b bᵀ).
pub fn frequency_bound<const N: usize>(k: &UniformStrain<N>, gradient_norm: f64, modulus: f64, mass: &[f64; N]) -> f64 {
    let lightest = mass.iter().copied().fold(f64::INFINITY, f64::min);
    modulus * gradient_norm / (k.volume * lightest)
}

/// A hexahedron's size h, hourglass coefficient (damping or stiffness) and frequency bound ω²
/// (with the hourglass stiffness), at density `density`.
#[derive(Debug, Clone, Copy)]
pub struct HexScales {
    pub size: f64,
    pub coefficient: f64,
    pub omega2: f64,
}

/// The [`HexScales`] of a hexahedron in configuration `k` with shape vectors `gamma`.
pub fn hex_scales(
    k: &Hex8,
    gamma: &[[f64; 8]; 4],
    mass: &[f64; 8],
    mat: &Material,
    density: f64,
    hourglass: fieldcraft_model::Hourglass,
) -> HexScales {
    let norm = k.gradient_norm();
    // h: the edge of the cube with the element's volume and gradient norm.
    let size = k.volume / (2.0 * norm).sqrt();
    let wave2 = mat.stiffest / density;
    let mut omega2 = frequency_bound(k, norm, mat.modulus, mass);
    let kappa = hourglass.coefficient;
    let coefficient = match hourglass.kind {
        HourglassKind::Viscous => kappa * density * wave2.sqrt() * size * size / 16.0,
        HourglassKind::Stiffness => {
            let stiffness = kappa * density * wave2 * size / 16.0;
            omega2 += stiffness * hourglass_gram_bound(gamma, mass);
            stiffness
        }
    };
    HexScales { size, coefficient, omega2 }
}

/// Elements per parallel task.
const BLOCK: usize = 128;

/// Tasks whose element forces are computed before they are added (bounds the memory of the
/// element forces held at once).
const ROUND: usize = 256;

/// What one element gives an evaluation.
enum Outcome<const N: usize> {
    /// Eroded earlier: nothing.
    Skipped,
    /// Eroded now.
    Eroded(Failure),
    Forces(Forces<N>),
}

/// An element's nodal forces by origin and its critical increment.
struct Forces<const N: usize> {
    hourglass: Option<[[f64; 3]; N]>,
    stress: [[f64; 3]; N],
    viscous: Option<[[f64; 3]; N]>,
    critical: f64,
}

/// A block of consecutive elements of one kind (positions `items` in their list) with the
/// slice of element states that holds theirs, from state index `base` on.
struct Piece<'s, T> {
    items: std::ops::Range<usize>,
    base: usize,
    states: &'s mut [ElementState],
    out: Vec<Result<T>>,
}

/// Cut the element states into one slice per block of [`BLOCK`] elements of a kind whose state
/// indices (`index(i)` for item i) increase, as [`System`] builds them.
fn pieces<'s, T>(states: &'s mut [ElementState], count: usize, index: &impl Fn(usize) -> usize) -> Result<Vec<Piece<'s, T>>> {
    let ranges: Vec<std::ops::Range<usize>> = (0..count.div_ceil(BLOCK)).map(|k| k * BLOCK..((k + 1) * BLOCK).min(count)).collect();
    let mut out = Vec::with_capacity(ranges.len());
    let mut rest = states;
    let mut at = 0;
    for (k, items) in ranges.iter().enumerate() {
        let base = index(items.start);
        let end = ranges.get(k + 1).map_or(at + rest.len(), |r| index(r.start));
        if base < at || end < base || end - at > rest.len() {
            return Err(Error::Internal("element states out of order".into()));
        }
        let (_, tail) = std::mem::take(&mut rest).split_at_mut(base - at);
        let (mine, tail) = tail.split_at_mut(end - base);
        out.push(Piece { items: items.clone(), base, states: mine, out: Vec::new() });
        rest = tail;
        at = end;
    }
    Ok(out)
}

/// `f(i, state)` for every element i of a kind, in parallel by blocks and a round of blocks at a
/// time; the outcomes go to `apply` in element order, so whatever it sums is summed as a
/// sequential loop would. The first error in element order stops the evaluation.
fn run<T: Send>(
    states: &mut [ElementState],
    count: usize,
    index: impl Fn(usize) -> usize + Sync,
    f: impl Fn(usize, &mut ElementState) -> Result<T> + Sync + Send,
    mut apply: impl FnMut(usize, T) -> Result<()>,
) -> Result<()> {
    let mut all = pieces::<T>(states, count, &index)?;
    for round in all.chunks_mut(ROUND) {
        fieldcraft_parallel::for_each_mut(round, |_, piece| {
            let Piece { items, base, states, out } = piece;
            *out = items
                .clone()
                .map(|i| match states.get_mut(index(i).wrapping_sub(*base)) {
                    Some(es) => f(i, es),
                    None => Err(Error::Internal("element state".into())),
                })
                .collect();
        });
        for piece in round.iter_mut() {
            for (i, outcome) in piece.items.clone().zip(std::mem::take(&mut piece.out)) {
                apply(i, outcome?)?;
            }
        }
    }
    Ok(())
}

/// Add an element's forces into the workspace and return its critical increment.
fn add<const N: usize>(ws: &mut Workspace, nodes: &[usize; N], f: &Forces<N>) -> f64 {
    if let Some(h) = &f.hourglass {
        scatter(&mut ws.hourglass_force, nodes, h);
    }
    scatter(&mut ws.stress_force, nodes, &f.stress);
    if let Some(q) = &f.viscous {
        scatter(&mut ws.viscous_force, nodes, q);
    }
    f.critical
}

/// One hexahedron's update and forces.
fn hexahedron(sys: &System, h: &crate::system::Hex, m: &Motion, es: &mut ElementState) -> Result<Outcome<8>> {
    let Some(mat) = sys.materials.get(h.material) else { return Ok(Outcome::Skipped) };
    if es.failed.is_some() {
        return Ok(Outcome::Skipped);
    }
    let hourglass = fieldcraft_model::Hourglass { kind: sys.hourglass, coefficient: sys.hourglass_coefficient };
    let e = Nodal::gather(sys, &h.nodes, m);
    let x_now = e.positions(0.0);
    let (mid, now) = if sys.nonlinear {
        let at = |x: &[[f64; 3]; 8]| Hex8::hex8_deformed(x).map_err(|_| inverted(h.label, m.time));
        (at(&e.positions(0.5))?, at(&x_now)?)
    } else {
        (h.kinematics, h.kinematics)
    };
    let g = mid.gradient(&e.du);
    let q = rotation_of(sys, mat, es, &g, h.label, m.time)?;
    let (d, _) = rotation::split(&g);
    let d = mechanical(d, h.thermal_rate * m.dt);
    let (rho, c) = current(mat, h.kinematics.volume, now.volume);
    // Mass scaling raises the density the dynamics see (not the one adiabatic heating uses).
    let (rho_s, c_s) = (rho * h.mass_factor, c / h.mass_factor.sqrt());
    let stiffer = match &mat.hyperelastic {
        Some(model) => {
            let mut f = h.kinematics.gradient(&e.u);
            for k in 0..3 {
                f[k][k] += 1.0;
            }
            let j = fieldcraft_fem::dense::det3(&f);
            hyperelastic(model, mat, es, &f, j, None, &d, mid.volume).map_err(|_| inverted(h.label, m.time))?
        }
        None => {
            update(mat, es, &q, &d, d[0] + d[1] + d[2], m.dt, rho, mid.volume);
            1.0
        }
    };
    if let Some(failure) = es.failed {
        return Ok(Outcome::Eroded(failure));
    }
    // Hourglass forces with the current shape vectors and size (fixed without nonlinear
    // geometry).
    let (gamma, scales) = if sys.nonlinear {
        let gamma = now.hourglass_vectors(&x_now);
        let scales = hex_scales(&now, &gamma, &h.mass, mat, rho_s, hourglass);
        (gamma, scales)
    } else {
        (h.gamma, h.rest)
    };
    let mut f = [[0.0; 3]; 8];
    let modes = if sys.hourglass == HourglassKind::Viscous { &e.v } else { &e.du };
    for (k, g) in gamma.iter().enumerate() {
        let rate: [f64; 3] = [0, 1, 2].map(|i| modes.iter().zip(g).map(|(x, gk)| x[i] * gk).sum());
        // Viscous: Q = C q̇; stiffness: Q accumulates k Δq.
        let generalized = match (sys.hourglass, es.hourglass.get_mut(k)) {
            (HourglassKind::Stiffness, Some(qa)) => {
                for i in 0..3 {
                    qa[i] += scales.coefficient * rate[i];
                }
                *qa
            }
            _ => rate.map(|r| scales.coefficient * r),
        };
        for (fk, gk) in f.iter_mut().zip(g) {
            for i in 0..3 {
                fk[i] += generalized[i] * gk;
            }
        }
    }
    let element = Element { now: &now, mid: &mid, density: rho_s, wave_speed: c_s, length: 2.0 * c_s / (stiffer * scales.omega2).sqrt() };
    let mut forces = element.finish(&es.point.stress, &e.v, m.bulk);
    forces.hourglass = Some(f);
    Ok(Outcome::Forces(forces))
}

/// A tetrahedron's first pass: its increment (none when it has eroded).
fn tetrahedron_increment(sys: &System, t: &crate::system::Tet, m: &Motion, state: &State) -> Result<Option<TetIncrement>> {
    if state.elements.get(t.index).is_none_or(|e| e.failed.is_some()) {
        return Ok(None);
    }
    let e = Nodal::gather(sys, &t.nodes, m);
    let (mid, now) = if sys.nonlinear {
        let at = |x: &[[f64; 3]; 4]| Tet4::tet4(x).map_err(|_| inverted(t.label, m.time));
        (at(&e.positions(0.5))?, at(&e.positions(0.0))?)
    } else {
        (t.kinematics, t.kinematics)
    };
    let g = mid.gradient(&e.du);
    Ok(Some(TetIncrement { g, mid, now, v: e.v }))
}

/// A tetrahedron's second pass: its update with the mean of its nodes' volumetric increments.
fn tetrahedron(
    sys: &System,
    t: &crate::system::Tet,
    inc: Option<TetIncrement>,
    theta: &[f64],
    ratio: &[f64],
    m: &Motion,
    es: &mut ElementState,
) -> Result<Outcome<4>> {
    let (Some(mat), Some(inc)) = (sys.materials.get(t.material), inc) else { return Ok(Outcome::Skipped) };
    let averaged = t.nodes.iter().map(|n| theta.get(*n).copied().unwrap_or(0.0)).sum::<f64>() / 4.0;
    let q = rotation_of(sys, mat, es, &inc.g, t.label, m.time)?;
    let (d, _) = rotation::split(&inc.g);
    let thermal = t.thermal_rate * m.dt;
    let d = mechanical(d, thermal);
    let (rho, c) = current(mat, t.kinematics.volume, inc.now.volume);
    let (rho_s, c_s) = (rho * t.mass_factor, c / t.mass_factor.sqrt());
    let stiffer = match &mat.hyperelastic {
        // The average nodal pressure of a decoupled material (Bonet & Burton 1998, at finite
        // strain): each node's pressure U′(J_a) from its tetrahedra's volume ratio, the element's
        // the mean of its nodes', with the isochoric stress of its own F.
        Some(model) => {
            let u = gather(m.u, &t.nodes);
            let mut f = t.kinematics.gradient(&u);
            for k in 0..3 {
                f[k][k] += 1.0;
            }
            let j = fieldcraft_fem::dense::det3(&f);
            let pressure = if model.is_decoupled() {
                let p: Option<f64> = t.nodes.iter().map(|n| model.volumetric(ratio.get(*n).copied().unwrap_or(1.0), false).map(|(p, _)| p)).sum();
                p.map(|p| p / 4.0)
            } else {
                None
            };
            hyperelastic(model, mat, es, &f, j, pressure, &d, inc.mid.volume).map_err(|_| inverted(t.label, m.time))?
        }
        None => {
            update(mat, es, &q, &d, averaged - 3.0 * thermal, m.dt, rho, inc.mid.volume);
            1.0
        }
    };
    if let Some(failure) = es.failed {
        return Ok(Outcome::Eroded(failure));
    }
    let omega2 = stiffer * if sys.nonlinear { frequency_bound(&inc.now, inc.now.gradient_norm(), mat.modulus, &t.mass) } else { t.omega2 };
    let element = Element { now: &inc.now, mid: &inc.mid, density: rho_s, wave_speed: c_s, length: 2.0 * c_s / omega2.sqrt() };
    Ok(Outcome::Forces(element.finish(&es.point.stress, &inc.v, m.bulk)))
}

/// Evaluate the internal forces at the motion `m` into `ws`, updating the element states, and
/// return the smallest critical increment of the elements still in the analysis.
///
/// The elements are updated in parallel (each owns its state); their nodal forces are added in
/// element order, hexahedra first, so every sum is the one of a sequential loop and the result
/// does not depend on the number of threads.
pub fn internal_forces(sys: &System, m: &Motion, state: &mut State, ws: &mut Workspace) -> Result<f64> {
    for f in [&mut ws.stress_force, &mut ws.hourglass_force, &mut ws.viscous_force] {
        f.fill(0.0);
    }
    state.eroded.clear();
    let mut critical = f64::INFINITY;
    let State { elements, eroded } = state;
    let hexes = &sys.hexes;
    run(
        elements,
        hexes.len(),
        |i| hexes.get(i).map_or(usize::MAX, |h| h.index),
        |i, es| match hexes.get(i) {
            Some(h) => hexahedron(sys, h, m, es),
            None => Ok(Outcome::Skipped),
        },
        |i, outcome| {
            let Some(h) = hexes.get(i) else { return Ok(()) };
            match outcome {
                Outcome::Skipped => {}
                Outcome::Eroded(failure) => eroded.push((h.index, failure)),
                Outcome::Forces(f) => critical = critical.min(add(ws, &h.nodes, &f)),
            }
            Ok(())
        },
    )?;
    if !sys.tets.is_empty() {
        let tets = &sys.tets;
        // Pass 1: each tetrahedron's increment, and each node's mean volumetric increment.
        let increments = {
            let st: &State = state;
            fieldcraft_parallel::map_grained(tets.len(), BLOCK, |k| match tets.get(k) {
                Some(t) => tetrahedron_increment(sys, t, m, st),
                None => Ok(None),
            })
        };
        ws.theta.fill(0.0);
        ws.tet_volume.fill(0.0);
        for ((t, slot), inc) in tets.iter().zip(ws.tet_increment.iter_mut()).zip(increments) {
            *slot = inc?;
            let Some(inc) = slot else { continue };
            for node in t.nodes {
                if let (Some(x), Some(w)) = (ws.theta.get_mut(node), ws.tet_volume.get_mut(node)) {
                    *x += inc.mid.volume * (inc.g[0][0] + inc.g[1][1] + inc.g[2][2]);
                    *w += inc.mid.volume;
                }
            }
        }
        for (x, w) in ws.theta.iter_mut().zip(&ws.tet_volume) {
            if *w > 0.0 {
                *x /= w;
            }
        }
        if !ws.tet_ratio.is_empty() {
            ws.tet_ratio.fill(0.0);
            ws.tet_rest.fill(0.0);
            for (t, inc) in tets.iter().zip(&ws.tet_increment) {
                let Some(inc) = inc else { continue };
                for node in t.nodes {
                    if let (Some(x), Some(w)) = (ws.tet_ratio.get_mut(node), ws.tet_rest.get_mut(node)) {
                        *x += inc.now.volume;
                        *w += t.kinematics.volume;
                    }
                }
            }
            for (x, w) in ws.tet_ratio.iter_mut().zip(&ws.tet_rest) {
                *x = if *w > 0.0 { *x / w } else { 1.0 };
            }
        }
        // Pass 2: each tetrahedron's update.
        let State { elements, eroded } = state;
        let (theta, ratio, increments) = (&ws.theta, &ws.tet_ratio, &ws.tet_increment);
        let mut forces: Vec<(usize, Forces<4>)> = Vec::new();
        run(
            elements,
            tets.len(),
            |k| tets.get(k).map_or(usize::MAX, |t| t.index),
            |k, es| match tets.get(k) {
                Some(t) => tetrahedron(sys, t, increments.get(k).copied().flatten(), theta, ratio, m, es),
                None => Ok(Outcome::Skipped),
            },
            |k, outcome| {
                match (tets.get(k), outcome) {
                    (Some(t), Outcome::Eroded(failure)) => eroded.push((t.index, failure)),
                    (Some(_), Outcome::Forces(f)) => forces.push((k, f)),
                    _ => {}
                }
                Ok(())
            },
        )?;
        for (k, f) in &forces {
            if let Some(t) = tets.get(*k) {
                critical = critical.min(add(ws, &t.nodes, f));
            }
        }
    }
    for r in &sys.rods {
        let (sigma, strain) = r.forces(&sys.xyz, m.u, sys.nonlinear, &mut ws.stress_force);
        if let Some(es) = state.elements.get_mut(r.index) {
            es.point.stress = sigma;
            es.strain = strain;
        }
        if let Some(mat) = sys.materials.get(r.material) {
            critical = critical.min(r.length / mat.wave_speed);
        }
    }
    for l in &sys.discrete.links {
        l.forces(&sys.xyz, m.u, m.v, sys.nonlinear, &mut ws.stress_force, &mut ws.viscous_force);
    }
    if let Some((dt, _)) = &sys.discrete_increment {
        critical = critical.min(*dt);
    }
    Ok(critical)
}

/// A tetrahedron's increment between the two passes: the gradient of its displacement increment,
/// its midpoint and current kinematics, and its nodal velocities.
#[derive(Debug, Clone, Copy)]
struct TetIncrement {
    g: Mat3,
    mid: Tet4,
    now: Tet4,
    v: [[f64; 3]; 4],
}

/// What the hexahedron and the tetrahedron share once their stress is known.
struct Element<'a, const N: usize> {
    now: &'a UniformStrain<N>,
    mid: &'a UniformStrain<N>,
    density: f64,
    wave_speed: f64,
    length: f64,
}

impl<const N: usize> Element<'_, N> {
    /// The stress and bulk viscosity forces and the critical increment.
    fn finish(&self, sigma: &[f64; 6], v: &[[f64; 3]; N], bulk: &BulkViscosity) -> Forces<N> {
        let rate = self.mid.divergence(v);
        let q = bulk_pressure(rate, self.density, self.length, self.wave_speed, bulk);
        Forces {
            hourglass: None,
            stress: self.now.forces(sigma),
            viscous: (q != 0.0).then(|| self.now.pressure_forces(q)),
            critical: critical_increment(rate, self.length, self.wave_speed, bulk),
        }
    }
}
