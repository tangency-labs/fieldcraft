//! Nédélec edge elements of the first kind on tetrahedra (J.-C. Nédélec, *Mixed finite elements
//! in ℝ³*, Numer. Math. 35:315–341, 1980), for the electric field of time-harmonic Maxwell
//! problems (docs/rfc/0001-electromagnetics.md §3.1), with the hierarchical basis of J. P. Webb,
//! *Hierarchal vector basis functions of arbitrary order for triangular and tetrahedral finite
//! elements*, IEEE Trans. Antennas Propag. 47(8):1244–1253, 1999:
//!
//! - **order 1** (6 functions): the Whitney edge functions w_ij = λ_i∇λ_j − λ_j∇λ_i (H. Whitney,
//!   *Geometric Integration Theory*, Princeton 1957), one per edge, whose tangential component
//!   integrates to 1 along their own edge and to 0 along the others;
//! - **order 2** (20 functions, Webb's "order 1.5", complete to first order): the Whitney
//!   functions, the gradients ∇(λ_iλ_j) of the edges, and two functions per face i < j < k,
//!   λ_i w_jk and λ_j w_ki (the third, λ_k w_ij, is dependent on them).
//!
//! **Orientation** (RFC §8.1). The element is built on its corners **sorted by global label**
//! ([`EdgeTet::of_labels`]): every local edge then runs from its lower to its higher label and
//! every local face is its sorted triple, the same from every tetrahedron that shares them, so
//! the element's functions are the global ones and assembly needs no sign and no 2 × 2
//! combination; [`entity_keys`] names the global edges and faces of the local ones. A sign per
//! entity would be right for the Whitney functions but not for the face functions: a face seen
//! in another vertex order holds a pair that is a combination of the global pair. The sorted
//! corners may make a left-handed tetrahedron: the gradients take the signed Jacobian, the
//! volume its magnitude.
//!
//! **Local unknowns** ([`Order::role`]), entity by entity as the degree-of-freedom layout numbers
//! them (an entity's unknowns contiguous, RFC §8.2): at order 1 function k is the Whitney
//! function of edge [`EDGES`]`[k]`; at order 2 each edge of [`EDGES`] has its Whitney function
//! then its gradient, and each face of [`FACES`] then its two functions. The Whitney functions
//! span a + b × x and the gradients of the linear functions; order 2 adds every linear field and
//! the gradients of the quadratics.
//!
//! **Matrices** for constant coefficients, diagonal in the global axes (an isotropic material is
//! the same value three times; a perfectly matched layer has three different ones): the
//! curl–curl matrix ∫ curl Nᵢ · diag(ν) curl Nⱼ dV and the mass matrix ∫ Nᵢ · diag(ε) Nⱼ dV, the
//! two of the weak form of curl ν curl E − k₀² ε E = 0. They are real; a complex coefficient (a
//! lossy material, a stretched coordinate) is the sum of the real matrices of its parts, formed
//! where the system is assembled. Every function is a polynomial in the barycentric coordinates
//! times barycentric gradients, so the products integrate exactly by
//!
//! ∫_T λ₀^a λ₁^b λ₂^c λ₃^d dV = 3! V a! b! c! d! / (a + b + c + d + 3)!
//!
//! (O. C. Zienkiewicz, R. L. Taylor, *The Finite Element Method*, vol. 1, 5th ed. 2000, §8.8;
//! Z. J. Cendes, *Vector finite elements for electromagnetic field computation*, IEEE Trans.
//! Magn. 27(5):3958–3966, 1991), with no quadrature. Curls come from
//! curl(m(λ) ∇λ_k) = Σ_p ∂m/∂λ_p ∇λ_p × ∇λ_k.
//!
//! A tetrahedron is accepted when FieldCraft's one validity criterion holds for its corner
//! Jacobian in either handedness ([`fieldcraft_model::MIN_NORMALIZED_JACOBIAN`]).
//!
//! Adapted from EMCraft's `crates/fem/src/nedelec.rs` and `poly.rs` (MIT OR Apache-2.0, by this
//! module's author). Curved (quadratic-geometry) tetrahedra through the covariant Piola map are
//! not here yet.

use crate::dense::{Matrix, inverse3};
use crate::dof::{DofKind, DofLayout, DofSet, Entity, LayoutBuilder};
use crate::frame::{Vec3, add, cross, dot, scale};
use crate::{Error, Result};

/// Barycentric coordinates λ₀..λ₃ of a point of a tetrahedron (they sum to 1).
pub type Barycentric = [f64; 4];

/// A mesh edge by its nodes' labels, the lower first ([`Entity::Edge`]).
pub type EdgeKey = [u32; 2];

/// A mesh face by its nodes' labels, ascending ([`Entity::Face`]).
pub type FaceKey = [u32; 3];

/// The local edges: edge k joins corners `EDGES[k]`, lower first.
pub const EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

/// The local faces: face f has corners `FACES[f]`, ascending.
pub const FACES: [[usize; 3]; 4] = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];

/// The order of an edge element.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Order {
    One,
    Two,
}

/// What a local function is: which entity it belongs to, and which of that entity's unknowns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Role {
    /// The Whitney function of local edge `edge`.
    Whitney { edge: usize },
    /// The gradient ∇(λ_iλ_j) of local edge `edge` (order 2).
    Gradient { edge: usize },
    /// Function `which` (0: λ_i w_jk, 1: λ_j w_ki) of local face `face` (order 2).
    Face { face: usize, which: usize },
}

impl Order {
    /// Unknowns per edge: the Whitney function, and at order 2 the edge's gradient.
    pub fn per_edge(self) -> usize {
        match self {
            Order::One => 1,
            Order::Two => 2,
        }
    }

    /// Unknowns per face.
    pub fn per_face(self) -> usize {
        match self {
            Order::One => 0,
            Order::Two => 2,
        }
    }

    /// The unknowns each of its edges and faces carries.
    pub fn kinds(self) -> DofSet {
        match self {
            Order::One => DofSet::FIELD,
            Order::Two => DofSet::FIELDS,
        }
    }

    /// Basis functions per tetrahedron.
    pub fn functions(self) -> usize {
        EDGES.len() * self.per_edge() + FACES.len() * self.per_face()
    }

    /// What local function `k` is; `None` past the last one.
    pub fn role(self, k: usize) -> Option<Role> {
        let (per_edge, edge_unknowns) = (self.per_edge(), EDGES.len() * self.per_edge());
        match k {
            _ if k >= self.functions() => None,
            _ if k < edge_unknowns && k.is_multiple_of(per_edge) => Some(Role::Whitney { edge: k / per_edge }),
            _ if k < edge_unknowns => Some(Role::Gradient { edge: k / per_edge }),
            _ => Some(Role::Face { face: (k - edge_unknowns) / 2, which: (k - edge_unknowns) % 2 }),
        }
    }
}

/// The powers of λ₀..λ₃ in a monomial.
type Exponents = [u8; 4];

/// One term of a basis function before the element's shape is known:
/// `coef · Π λᵢ^exp[i] · ∇λ_grad`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Term {
    coef: f64,
    exp: Exponents,
    grad: usize,
}

const fn power(i: usize, j: Option<usize>) -> Exponents {
    let mut e = [0u8; 4];
    e[i] += 1;
    if let Some(j) = j {
        e[j] += 1;
    }
    e
}

const fn term(coef: f64, i: usize, j: Option<usize>, grad: usize) -> Term {
    Term { coef, exp: power(i, j), grad }
}

/// Whitney: λ_i∇λ_j − λ_j∇λ_i.
const fn whitney(i: usize, j: usize) -> [Term; 2] {
    [term(1.0, i, None, j), term(-1.0, j, None, i)]
}

/// The order-1 functions, one per edge.
const FIRST: [[Term; 2]; 6] = {
    let mut out = [[term(0.0, 0, None, 0); 2]; 6];
    let mut e = 0;
    while e < 6 {
        let [i, j] = EDGES[e];
        out[e] = whitney(i, j);
        e += 1;
    }
    out
};

/// The order-2 functions as their two terms, entity by entity ([`Order::role`]).
const SECOND: [[Term; 2]; 20] = {
    let mut out = [[term(0.0, 0, None, 0); 2]; 20];
    let mut e = 0;
    while e < 6 {
        let [i, j] = EDGES[e];
        // The edge's gradient ∇(λ_iλ_j) = λ_i∇λ_j + λ_j∇λ_i follows its Whitney function.
        out[2 * e] = whitney(i, j);
        out[2 * e + 1] = [term(1.0, i, None, j), term(1.0, j, None, i)];
        e += 1;
    }
    let mut f = 0;
    while f < 4 {
        let [i, j, k] = FACES[f];
        // λ_i w_jk = λ_iλ_j∇λ_k − λ_iλ_k∇λ_j and λ_j w_ki = λ_jλ_k∇λ_i − λ_jλ_i∇λ_k.
        out[12 + 2 * f] = [term(1.0, i, Some(j), k), term(-1.0, i, Some(k), j)];
        out[13 + 2 * f] = [term(1.0, j, Some(k), i), term(-1.0, j, Some(i), k)];
        f += 1;
    }
    out
};

/// The basis functions of an order as terms (the same on every tetrahedron).
fn terms(order: Order) -> &'static [[Term; 2]] {
    match order {
        Order::One => &FIRST,
        Order::Two => &SECOND,
    }
}

/// A vector polynomial in the barycentric coordinates: Σ c · Π λᵢ^exp[i], one term per exponent.
#[derive(Debug, Clone, Default, PartialEq)]
struct Poly {
    terms: Vec<(Exponents, Vec3)>,
}

/// n! for the small n the integrals meet (a product of two order-2 functions has degree 4).
fn factorial(n: u32) -> f64 {
    (1..=n).map(f64::from).product()
}

/// ∫_T Π λᵢ^αᵢ dV over a tetrahedron of volume `volume`.
fn monomial_integral(alpha: Exponents, volume: f64) -> f64 {
    let total: u32 = alpha.iter().map(|&a| u32::from(a)).sum();
    6.0 * volume * alpha.iter().map(|&a| factorial(u32::from(a))).product::<f64>() / factorial(total + 3)
}

/// Π λᵢ^exp[i], by multiplications (no `powi`, whose rounding is the platform's).
fn monomial(exp: Exponents, l: Barycentric) -> f64 {
    exp.iter().zip(l).fold(1.0, |m, (&e, x)| (0..e).fold(m, |m, _| m * x))
}

/// Σ cᵢ vᵢ.
fn combine(pairs: impl IntoIterator<Item = (f64, Vec3)>) -> Vec3 {
    pairs.into_iter().fold([0.0; 3], |sum, (c, v)| add(sum, scale(v, c)))
}

fn sum_exponents(a: Exponents, b: Exponents) -> Exponents {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
}

/// The polynomial of some terms, those with the same exponents merged (in the order they come).
impl FromIterator<(Exponents, Vec3)> for Poly {
    fn from_iter<I: IntoIterator<Item = (Exponents, Vec3)>>(terms: I) -> Poly {
        let terms = terms.into_iter().fold(Vec::<(Exponents, Vec3)>::new(), |mut merged, (exp, coef)| {
            match merged.iter_mut().find(|t| t.0 == exp) {
                Some(t) => t.1 = add(t.1, coef),
                None => merged.push((exp, coef)),
            }
            merged
        });
        Poly { terms }
    }
}

impl Poly {
    fn eval(&self, l: Barycentric) -> Vec3 {
        combine(self.terms.iter().map(|(exp, c)| (monomial(*exp, l), *c)))
    }

    /// ∫_T self · diag(w) other dV. The volume scales one factor before the other multiplies
    /// it: a curl's coefficients are of order 1/h², and their product alone would overflow on a
    /// tetrahedron of edges below about 1e-77 whose integral is finite.
    fn inner(&self, other: &Poly, w: Vec3, volume: f64) -> f64 {
        self.terms
            .iter()
            .map(|(ea, ca)| {
                let ca = scale(*ca, volume);
                other
                    .terms
                    .iter()
                    .map(|(eb, cb)| {
                        let c = w[0] * ca[0] * cb[0] + w[1] * ca[1] * cb[1] + w[2] * ca[2] * cb[2];
                        c * monomial_integral(sum_exponents(*ea, *eb), 1.0)
                    })
                    .sum::<f64>()
            })
            .sum()
    }
}

/// The function of `terms` on a tetrahedron with barycentric gradients `g`.
fn realise(terms: &[Term], g: &[Vec3; 4]) -> Poly {
    terms.iter().filter_map(|t| g.get(t.grad).map(|gk| (t.exp, scale(*gk, t.coef)))).collect()
}

/// Its curl: curl(m(λ) ∇λ_k) = Σ_p ∂m/∂λ_p ∇λ_p × ∇λ_k.
fn realise_curl(terms: &[Term], g: &[Vec3; 4]) -> Poly {
    terms
        .iter()
        .filter_map(|t| g.get(t.grad).map(|gk| (t, *gk)))
        .flat_map(|(t, gk)| {
            g.iter().zip(t.exp).enumerate().filter(move |&(q, (_, power))| power > 0 && q != t.grad).map(move |(q, (gq, power))| {
                let lowered: Exponents = std::array::from_fn(|i| t.exp.get(i).copied().unwrap_or(0) - u8::from(i == q));
                (lowered, scale(cross(*gq, gk), t.coef * f64::from(power)))
            })
        })
        .collect()
}

/// The order in which a tetrahedron's corners are listed for [`EdgeTet::new`]: by ascending
/// global label, so `sorted[k]` is the position in `labels` of the k-th lowest. Refused when two
/// labels are the same.
pub fn sort_corners(labels: [u32; 4]) -> Result<[usize; 4]> {
    let mut order = [0, 1, 2, 3];
    order.sort_by_key(|&k| labels.get(k).copied().unwrap_or(u32::MAX));
    let ascending = order.windows(2).all(|w| labels.get(w[0]) < labels.get(w[1]));
    if !ascending {
        return Err(Error::Invalid(format!("a tetrahedron's corners must be four different nodes, got {labels:?}")));
    }
    Ok(order)
}

/// The global edges and faces of a tetrahedron's local ones ([`EDGES`], [`FACES`] of its sorted
/// corners), as sorted pairs and triples of its corners' labels: what a numbering of the mesh's
/// edges and faces is looked up by.
pub fn entity_keys(labels: [u32; 4]) -> Result<([EdgeKey; 6], [FaceKey; 4])> {
    let sorted = sort_corners(labels)?.map(|k| labels.get(k).copied().unwrap_or_default());
    let at = |i: usize| sorted.get(i).copied().unwrap_or_default();
    Ok((EDGES.map(|e| e.map(at)), FACES.map(|f| f.map(at))))
}

/// Declares what a tetrahedron of nodes `labels` needs of a layout: its edges, and at order 2
/// its faces, each with [`Order::kinds`]. A mesh's elements are all of one order: [`slots`]
/// refuses one whose edge or face an element of the other order shares.
pub fn declare(layout: &mut LayoutBuilder, labels: [u32; 4], order: Order) -> Result<()> {
    let (edges, faces) = entity_keys(labels)?;
    for e in edges {
        layout.edge(e, order.kinds());
    }
    if order.per_face() > 0 {
        for f in faces {
            layout.face(f, order.kinds());
        }
    }
    Ok(())
}

/// The layout slots of a tetrahedron's edges, then (order 2) its faces, in the element's local
/// order ([`Order::role`]): with [`Order::kinds`], what `add_block` assembles its matrices by.
/// `None` when the layout lacks one of them, or one carries field unknowns other than this
/// order's: an element of the other order shares it, and the field would not be continuous
/// across it (an order-1 element has no second unknown on the edge to match). Unknowns of
/// other kinds on a slot are no concern of the element.
pub fn slots(layout: &DofLayout, labels: [u32; 4], order: Order) -> Option<Vec<usize>> {
    let (edges, faces) = entity_keys(labels).ok()?;
    let faces = faces.iter().take(if order.per_face() > 0 { FACES.len() } else { 0 }).map(|f| layout.face(*f));
    let slots: Vec<usize> = edges.iter().map(|e| layout.edge(*e)).chain(faces).collect::<Option<_>>()?;
    let fields = |slot: &usize| layout.kinds(*slot).iter().filter(|k| DofSet::FIELDS.contains(*k)).collect::<DofSet>();
    slots.iter().all(|s| fields(s) == order.kinds()).then_some(slots)
}

/// What a column of the discrete gradient ([`gradient`]) is: a nodal potential's function.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Potential {
    /// The linear function of a node, by label.
    Node(u32),
    /// The quadratic function λ_aλ_b of an edge, by its layout slot (order 2).
    Edge(usize),
}

/// The discrete gradient G of a layout's edge unknowns (RFC §8.2, §8.4): the edge-element field
/// of ∇φ, φ a sum of nodal potentials, is G φ, with
/// - the Whitney unknown (E1) of edge a → b taking φ_b − φ_a: w_ab's integral along its edge is
///   1, so ∫_ab ∇φ · t = φ(b) − φ(a) for a linear φ;
/// - its gradient unknown (E2) taking the coefficient of its own edge function: ∇(λ_aλ_b) is
///   that function at order 2.
///
/// As (row: the unknown, column, coefficient) triples, by row. The curl–curl matrix K has
/// K G = 0, the cheapest strong test of a numbering (a wrong orientation or incidence breaks it),
/// and G's columns span the gradients in K's kernel, which a cavity's eigen step projects out.
/// They are all of the kernel unless the domain has a handle (a ring) and walls left without a
/// conductor, where fields that are no gradient also lie in it.
pub fn gradient(layout: &DofLayout) -> Vec<(usize, Potential, f64)> {
    (0..layout.slots())
        .filter_map(|slot| match layout.entity(slot) {
            Some(Entity::Edge([a, b])) => Some((slot, a, b)),
            _ => None,
        })
        .flat_map(|(slot, a, b)| {
            let whitney = layout.dof(slot, DofKind::E1).into_iter().flat_map(move |r| [(r, Potential::Node(a), -1.0), (r, Potential::Node(b), 1.0)]);
            let own = layout.dof(slot, DofKind::E2).map(|r| (r, Potential::Edge(slot), 1.0));
            whitney.chain(own)
        })
        .collect()
}

fn flat() -> Error {
    Error::Degenerate("a tetrahedron is flat".into())
}

/// The unsigned volume of a tetrahedron and its barycentric gradients ∇λᵢ, its corners in the
/// order given.
fn geometry(corners: &[Vec3; 4]) -> Result<(f64, [Vec3; 4])> {
    // The rows are the edges from corner 0, so ∇λ₁, ∇λ₂, ∇λ₃ are the columns of J⁻¹.
    let j = fieldcraft_model::tet4_jacobian(corners);
    let (inv, det) = inverse3(&j).ok_or_else(flat)?;
    let column = |a: usize| inv.map(|row| row.get(a).copied().unwrap_or_default());
    let (g1, g2, g3) = (column(0), column(1), column(2));
    Ok((det.abs() / 6.0, [scale(add(add(g1, g2), g3), -1.0), g1, g2, g3]))
}

/// Coordinates that are finite numbers in range, and a tetrahedron that is not flat by
/// FieldCraft's validity criterion with its corners in the order given (the criterion is taken
/// at corner 0).
fn check(corners: &[Vec3; 4]) -> Result<()> {
    fieldcraft_model::check_coordinates(corners.as_flattened()).map_err(|e| Error::Invalid(e.to_string()))?;
    let j = fieldcraft_model::tet4_jacobian(corners);
    if fieldcraft_model::normalized_jacobian(&j).abs() > fieldcraft_model::MIN_NORMALIZED_JACOBIAN { Ok(()) } else { Err(flat()) }
}

/// A straight tetrahedron with its edge-element basis and the curls of its functions.
#[derive(Debug, Clone)]
pub struct EdgeTet {
    order: Order,
    corners: [Vec3; 4],
    volume: f64,
    grads: [Vec3; 4],
    basis: Vec<Poly>,
    curls: Vec<Poly>,
}

impl EdgeTet {
    /// The element on `corners`, which must be sorted by global label (see the module notes; or
    /// use [`EdgeTet::of_labels`]). Refused when a coordinate is not a finite number in range,
    /// when the tetrahedron, its corners in this order, is flat by FieldCraft's validity
    /// criterion, or when its basis is not finite.
    pub fn new(order: Order, corners: [Vec3; 4]) -> Result<EdgeTet> {
        check(&corners)?;
        let (volume, grads) = geometry(&corners)?;
        EdgeTet::on(order, corners, volume, grads)
    }

    /// The element of a mesh's tetrahedron, its corners at `xyz` with global node labels
    /// `labels`, sorted here; the order of its corners is returned with it ([`sort_corners`]).
    /// The validity criterion is applied to the corners in mesh order, as the meshers' gate
    /// applies it (`fieldcraft_model::Mesh`), and the volume and the gradients ∇λᵢ are computed
    /// in that order too, then listed in the sorted one (each belongs to its corner), so whether
    /// an element is accepted, and its matrices, do not depend on its nodes' labels.
    pub fn of_labels(order: Order, xyz: [Vec3; 4], labels: [u32; 4]) -> Result<(EdgeTet, [usize; 4])> {
        let sorted = sort_corners(labels)?;
        check(&xyz)?;
        let (volume, grads) = geometry(&xyz)?;
        let corners = sorted.map(|k| xyz.get(k).copied().unwrap_or_default());
        let grads = sorted.map(|k| grads.get(k).copied().unwrap_or_default());
        Ok((EdgeTet::on(order, corners, volume, grads)?, sorted))
    }

    /// The element on `corners` with their gradients ∇λᵢ, already checked. Refused when its
    /// basis is not finite: edges of 1e-155 m pass the validity criterion, which does not depend
    /// on size, but their gradients' cross products overflow.
    fn on(order: Order, corners: [Vec3; 4], volume: f64, grads: [Vec3; 4]) -> Result<EdgeTet> {
        let terms = terms(order);
        let basis: Vec<Poly> = terms.iter().map(|t| realise(t, &grads)).collect();
        let curls: Vec<Poly> = terms.iter().map(|t| realise_curl(t, &grads)).collect();
        let finite = |p: &Poly| p.terms.iter().flat_map(|(_, c)| c).all(|c| c.is_finite());
        if !(grads.as_flattened().iter().all(|g| g.is_finite()) && basis.iter().chain(&curls).all(finite)) {
            return Err(Error::Degenerate("a tetrahedron is too small or too thin for its matrices to be represented".into()));
        }
        Ok(EdgeTet { order, corners, volume, grads, basis, curls })
    }

    pub fn order(&self) -> Order {
        self.order
    }

    /// The corners, sorted by global label.
    pub fn corners(&self) -> [Vec3; 4] {
        self.corners
    }

    /// The unsigned volume.
    pub fn volume(&self) -> f64 {
        self.volume
    }

    /// ∇λᵢ.
    pub fn grads(&self) -> [Vec3; 4] {
        self.grads
    }

    /// ∫ curl Nᵢ · diag(ν) curl Nⱼ dV, ν (1/μᵣ) constant and diagonal in the global axes.
    /// Refused when an entry is not finite (a tetrahedron too small or too thin, or ν too large).
    pub fn curl_curl(&self, nu: Vec3) -> Result<Matrix> {
        gram(&self.curls, nu, self.volume)
    }

    /// ∫ Nᵢ · diag(ε) Nⱼ dV, ε (εᵣ) constant and diagonal in the global axes. Refused when an
    /// entry is not finite.
    pub fn mass(&self, eps: Vec3) -> Result<Matrix> {
        gram(&self.basis, eps, self.volume)
    }

    /// The basis functions at barycentric coordinates `l`.
    pub fn values(&self, l: Barycentric) -> Vec<Vec3> {
        self.basis.iter().map(|p| p.eval(l)).collect()
    }

    /// Their curls at `l`.
    pub fn curl_values(&self, l: Barycentric) -> Vec<Vec3> {
        self.curls.iter().map(|p| p.eval(l)).collect()
    }

    /// The barycentric coordinates of a point.
    pub fn barycentric(&self, p: Vec3) -> Barycentric {
        let [c0, ..] = self.corners;
        let [_, g1, g2, g3] = self.grads;
        let d = crate::frame::sub(p, c0);
        let (l1, l2, l3) = (dot(g1, d), dot(g2, d), dot(g3, d));
        [1.0 - l1 - l2 - l3, l1, l2, l3]
    }

    /// The point at barycentric coordinates `l`.
    pub fn point(&self, l: Barycentric) -> Vec3 {
        combine(l.into_iter().zip(self.corners))
    }
}

/// [∫ pᵢ · diag(w) pⱼ dV]: the upper triangle computed once, row by row, and mirrored. Refused
/// when an entry is not finite.
fn gram(p: &[Poly], w: Vec3, volume: f64) -> Result<Matrix> {
    let n = p.len();
    let upper: Vec<f64> = p.iter().enumerate().flat_map(|(i, pi)| p.iter().skip(i).map(move |pj| pi.inner(pj, w, volume))).collect();
    if !upper.iter().all(|v| v.is_finite()) {
        return Err(Error::Degenerate(
            "an element's matrix is not finite: its tetrahedron is too small or too thin, or its material too large, for it to be represented".into(),
        ));
    }
    // Row i of the upper triangle starts after the n + (n − 1) + … + (n − i + 1) entries above it.
    let at = |i: usize, j: usize| {
        let (i, j) = (i.min(j), i.max(j));
        i * n - i * i.saturating_sub(1) / 2 + (j - i)
    };
    Ok(Matrix { n, data: (0..n * n).map(|k| upper.get(at(k / n, k % n)).copied().unwrap_or(0.0)).collect() })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::{Cholesky, det3};
    use crate::frame::sub;
    use crate::quadrature::{Rule, gauss_legendre};
    use crate::shape::tet4;

    const ONE: Vec3 = [1.0; 3];

    fn tet() -> [Vec3; 4] {
        [[0.1, 0.0, 0.2], [1.3, 0.1, -0.1], [0.2, 0.9, 0.3], [0.4, 0.3, 1.1]]
    }

    /// Barycentric coordinates and weights summing to the volume, `n` points per direction.
    fn rule(e: &EdgeTet, n: u8) -> Vec<(Barycentric, f64)> {
        Rule::TetConical(n)
            .points()
            .into_iter()
            .map(|p| {
                let s = tet4(p.xi);
                ([s.n[0], s.n[1], s.n[2], s.n[3]], 6.0 * e.volume() * p.weight)
            })
            .collect()
    }

    fn close(a: Vec3, b: Vec3, tol: f64) -> bool {
        a.iter().zip(&b).all(|(x, y)| (x - y).abs() <= tol)
    }

    /// The coefficients of a field's projection: Σ xₖ Nₖ at `l`.
    fn field(x: &[f64], values: Vec<Vec3>) -> Vec3 {
        combine(x.iter().copied().zip(values))
    }

    #[test]
    fn roles_name_every_local_function_entity_by_entity() {
        assert_eq!((Order::One.functions(), Order::Two.functions()), (6, 20));
        assert_eq!(Order::One.role(5), Some(Role::Whitney { edge: 5 }));
        assert_eq!(Order::One.role(6), None);
        assert_eq!(Order::Two.role(2), Some(Role::Whitney { edge: 1 }));
        assert_eq!(Order::Two.role(3), Some(Role::Gradient { edge: 1 }));
        assert_eq!(Order::Two.role(12), Some(Role::Face { face: 0, which: 0 }));
        assert_eq!(Order::Two.role(19), Some(Role::Face { face: 3, which: 1 }));
        assert_eq!(Order::Two.role(20), None);
        // Every role once, and an entity's unknowns next to each other.
        let roles: Vec<Role> = (0..20).filter_map(|k| Order::Two.role(k)).collect();
        let entity = |r: &Role| match *r {
            Role::Whitney { edge } | Role::Gradient { edge } => edge,
            Role::Face { face, .. } => 6 + face,
        };
        assert!(roles.windows(2).all(|w| w[0] != w[1] && entity(&w[0]) <= entity(&w[1])));
        assert_eq!(roles.iter().map(entity).max(), Some(9));
    }

    #[test]
    fn a_tiny_tetrahedron_has_its_scaled_matrices() {
        // Edges of 1e-90: the curl–curl matrix scales as 1/h and the mass matrix as h, with no
        // overflow on the way.
        let h = 1e-90;
        for order in [Order::One, Order::Two] {
            let unit = EdgeTet::new(order, tet()).unwrap();
            let tiny = EdgeTet::new(order, tet().map(|p| scale(p, h))).unwrap();
            let (k, m) = (unit.curl_curl([1.0; 3]).unwrap(), unit.mass([1.0; 3]).unwrap());
            let (kt, mt) = (tiny.curl_curl([1.0; 3]).unwrap(), tiny.mass([1.0; 3]).unwrap());
            for i in 0..order.functions() {
                let close = |a: f64, b: f64| a.is_finite() && (a - b).abs() <= 1e-9 * b.abs();
                assert!(close(kt.get(i, i) * h, k.get(i, i)), "{order:?} {i}: {} against {}", kt.get(i, i) * h, k.get(i, i));
                assert!(close(mt.get(i, i) / h, m.get(i, i)), "{order:?} {i}: {} against {}", mt.get(i, i) / h, m.get(i, i));
            }
        }
    }

    #[test]
    fn whitney_functions_are_dual_to_the_edges() {
        let e = EdgeTet::new(Order::One, tet()).unwrap();
        let c = e.corners();
        for f in 0..EDGES.len() {
            for (k, [p, q]) in EDGES.iter().enumerate() {
                // ∫ N_f · t ds along edge k, t = c_q − c_p over the parameter s ∈ [0, 1].
                let t = sub(c[*q], c[*p]);
                let s: f64 = gauss_legendre(4)
                    .iter()
                    .map(|&(x, w)| {
                        let s = (x + 1.0) / 2.0;
                        let mut l = [0.0; 4];
                        l[*p] = 1.0 - s;
                        l[*q] = s;
                        w / 2.0 * dot(e.values(l)[f], t)
                    })
                    .sum();
                let expect = if f == k { 1.0 } else { 0.0 };
                assert!((s - expect).abs() < 1e-13, "function {f} on edge {k}: {s}");
            }
        }
    }

    #[test]
    fn order_one_matrices_match_the_closed_form() {
        // curl w_ij = 2 ∇λ_i × ∇λ_j (Whitney 1957) and ∫ λ_a λ_b = V (1 + δ_ab)/20
        // (Zienkiewicz & Taylor, vol. 1, §8.8).
        let e = EdgeTet::new(Order::One, tet()).unwrap();
        let (s, m) = (e.curl_curl(ONE).unwrap(), e.mass(ONE).unwrap());
        let (g, v) = (e.grads(), e.volume());
        let ll = |a: usize, b: usize| v * if a == b { 2.0 } else { 1.0 } / 20.0;
        for (p, [i, j]) in EDGES.iter().enumerate() {
            for (q, [k, l]) in EDGES.iter().enumerate() {
                let cs = 4.0 * v * dot(cross(g[*i], g[*j]), cross(g[*k], g[*l]));
                let cm =
                    ll(*i, *k) * dot(g[*j], g[*l]) - ll(*i, *l) * dot(g[*j], g[*k]) - ll(*j, *k) * dot(g[*i], g[*l]) + ll(*j, *l) * dot(g[*i], g[*k]);
                assert!((s.get(p, q) - cs).abs() <= 1e-12 * cs.abs().max(1.0), "S {p} {q}");
                assert!((m.get(p, q) - cm).abs() <= 1e-14, "M {p} {q}");
            }
        }
    }

    #[test]
    fn diagonal_coefficients_add_up_by_axis() {
        // diag(a) is Σ_d a_d of the one-axis matrices, and an isotropic a scales the matrix.
        let e = EdgeTet::new(Order::Two, tet()).unwrap();
        let (a, b) = ([2.0, 0.5, 3.0], [1.5, 4.0, 0.25]);
        let axis = |d: usize| {
            let mut w = [0.0; 3];
            w[d] = 1.0;
            w
        };
        for matrix in [EdgeTet::curl_curl, EdgeTet::mass] {
            let build = |e: &EdgeTet, w: Vec3| matrix(e, w).unwrap();
            let (whole, parts) = (build(&e, a), [0, 1, 2].map(|d| build(&e, axis(d))));
            let (iso, unit) = (build(&e, [2.5; 3]), build(&e, ONE));
            for k in 0..whole.data.len() {
                let sum: f64 = (0..3).map(|d| a[d] * parts[d].data[k]).sum();
                assert!((whole.data[k] - sum).abs() <= 1e-13 * sum.abs().max(1.0), "{k}");
                assert!((iso.data[k] - 2.5 * unit.data[k]).abs() <= 1e-13 * unit.data[k].abs().max(1.0), "{k}");
            }
            // An anisotropic material is not a multiple of the isotropic matrix.
            let (m, u) = (build(&e, b), unit);
            assert!(m.data.iter().zip(&u.data).any(|(x, y)| (x - b[0] * y).abs() > 1e-6 * y.abs().max(1e-3)));
        }
    }

    #[test]
    fn order_two_matrices_agree_with_quadrature() {
        // The mass integrand has degree 4: four points per direction integrate degree 5.
        let e = EdgeTet::new(Order::Two, tet()).unwrap();
        let (s, m) = (e.curl_curl(ONE).unwrap(), e.mass(ONE).unwrap());
        let n = Order::Two.functions();
        let (mut sq, mut mq) = (vec![0.0; n * n], vec![0.0; n * n]);
        for (l, w) in rule(&e, 4) {
            let (v, c) = (e.values(l), e.curl_values(l));
            for i in 0..n {
                for j in 0..n {
                    mq[i * n + j] += w * dot(v[i], v[j]);
                    sq[i * n + j] += w * dot(c[i], c[j]);
                }
            }
        }
        for (k, ((s, sq), (m, mq))) in s.data.iter().zip(&sq).zip(m.data.iter().zip(&mq)).enumerate() {
            assert!((s - sq).abs() < 1e-11, "S[{k}] {s} {sq}");
            assert!((m - mq).abs() < 1e-13, "M[{k}] {m} {mq}");
        }
        assert!((0..n).all(|i| (0..n).all(|j| s.get(i, j) == s.get(j, i) && m.get(i, j) == m.get(j, i))));
    }

    #[test]
    fn curls_match_finite_differences() {
        let e = EdgeTet::new(Order::Two, tet()).unwrap();
        let p = e.point([0.2, 0.3, 0.1, 0.4]);
        let h = 1e-6;
        let shifted = |axis: usize, by: f64| {
            let mut q = p;
            q[axis] += by;
            e.values(e.barycentric(q))
        };
        let (plus, minus): (Vec<_>, Vec<_>) = (0..3).map(|a| (shifted(a, h), shifted(a, -h))).unzip();
        // ∂_axis of component `comp` of function f.
        let d = |f: usize, axis: usize, comp: usize| (plus[axis][f][comp] - minus[axis][f][comp]) / (2.0 * h);
        for (f, c) in e.curl_values(e.barycentric(p)).iter().enumerate() {
            let fd = [d(f, 1, 2) - d(f, 2, 1), d(f, 2, 0) - d(f, 0, 2), d(f, 0, 1) - d(f, 1, 0)];
            assert!(close(fd, *c, 1e-6), "function {f}: {fd:?} vs {c:?}");
        }
    }

    /// The L² projection of a field onto the element's span; the largest error at sample points.
    fn projection_error(order: Order, f: impl Fn(Vec3) -> Vec3) -> f64 {
        let e = EdgeTet::new(order, tet()).unwrap();
        let mut x = vec![0.0; order.functions()];
        for (l, w) in rule(&e, 4) {
            let fl = f(e.point(l));
            for (r, v) in x.iter_mut().zip(e.values(l)) {
                *r += w * dot(v, fl);
            }
        }
        Cholesky::new(&e.mass(ONE).unwrap()).unwrap().solve(&mut x);
        [[0.25; 4], [0.7, 0.1, 0.1, 0.1], [0.0, 0.5, 0.5, 0.0], [0.1, 0.2, 0.3, 0.4]]
            .into_iter()
            .map(|l| sub(field(&x, e.values(l)), f(e.point(l))).iter().fold(0.0f64, |a, b| a.max(b.abs())))
            .fold(0.0, f64::max)
    }

    #[test]
    fn order_one_spans_rotations_and_order_two_linear_fields() {
        let constant = |_: Vec3| [0.3, -1.2, 0.7];
        assert!(projection_error(Order::One, constant) < 1e-12);
        // The Whitney functions span a + b × x.
        let rotation = |x: Vec3| add([0.3, -1.2, 0.7], cross([0.5, 0.1, -0.4], x));
        assert!(projection_error(Order::One, rotation) < 1e-12);
        let linear = |x: Vec3| [1.0 + 2.0 * x[0] - x[2], 0.5 * x[1] + 3.0 * x[2], x[0] - x[1] + 0.2];
        assert!(projection_error(Order::One, linear) > 1e-2);
        assert!(projection_error(Order::Two, linear) < 1e-11);
    }

    #[test]
    fn a_constant_curl_is_reproduced_exactly() {
        // E = ½ b × x has curl b: its edge integrals as order-1 coefficients give that curl at
        // every point (RFC §3.1).
        let b = [0.4, -1.1, 0.8];
        let e = EdgeTet::new(Order::One, tet()).unwrap();
        let c = e.corners();
        let f = |x: Vec3| scale(cross(b, x), 0.5);
        // ∫ E · t ds along an edge of a linear field: the field at the midpoint times t.
        let x: Vec<f64> = EDGES.iter().map(|&[p, q]| dot(f(scale(add(c[p], c[q]), 0.5)), sub(c[q], c[p]))).collect();
        for l in [[0.25; 4], [0.6, 0.2, 0.1, 0.1], [0.0, 0.0, 0.3, 0.7]] {
            let curl = field(&x, e.curl_values(l));
            assert!(close(curl, b, 1e-13), "{curl:?}");
        }
    }

    #[test]
    fn gradients_are_in_the_kernel_of_the_curl() {
        // Order 2 holds the gradients of quadratics: of the linear ones through the Whitney
        // functions (coefficients φ_j − φ_i), of λ_iλ_j through the edge gradients.
        let (phi, psi) = ([0.3, -1.0, 2.0, 0.5], [0.7, -0.2, 0.1, 1.5, -0.9, 0.4]);
        for order in [Order::One, Order::Two] {
            let e = EdgeTet::new(order, tet()).unwrap();
            let s = e.curl_curl(ONE).unwrap();
            let x: Vec<f64> = (0..order.functions())
                .map(|k| match order.role(k) {
                    Some(Role::Whitney { edge }) => EDGES.get(edge).map_or(0.0, |&[i, j]| phi[j] - phi[i]),
                    Some(Role::Gradient { edge }) => psi[edge],
                    _ => 0.0,
                })
                .collect();
            let norm = s.data.iter().fold(0.0f64, |a, b| a.max(b.abs()));
            for (i, r) in s.mul_vec(&x).iter().enumerate() {
                assert!(r.abs() < 1e-13 * norm, "{order:?} row {i}: {r}");
            }
        }
    }

    #[test]
    fn tangential_components_agree_across_a_face_seen_in_other_orders() {
        // Two tetrahedra share the face of nodes 12, 31, 40. The mesh lists their corners in
        // orders of its own, the face in a different vertex order from each side: sorting by
        // label makes the face's functions the same global ones on both (RFC §8.1), at order 2
        // where a sign per entity would not.
        let xyz = |n: u32| match n {
            5 => [0.0, 0.0, 0.0],
            12 => [1.0, 0.1, 0.0],
            31 => [0.1, 1.0, 0.2],
            40 => [0.2, 0.3, 1.0],
            _ => [1.0, 1.0, 1.0],
        };
        let (la, lb) = ([40, 5, 31, 12], [31, 77, 12, 40]);
        let (a, _) = EdgeTet::of_labels(Order::Two, la.map(xyz), la).unwrap();
        let (b, _) = EdgeTet::of_labels(Order::Two, lb.map(xyz), lb).unwrap();
        // The shared unknowns are those whose entities have the same keys in both.
        let ((ea, fa), (eb, fb)) = (entity_keys(la).unwrap(), entity_keys(lb).unwrap());
        let key = |edges: &[EdgeKey; 6], faces: &[FaceKey; 4], k: usize| match Order::Two.role(k) {
            Some(Role::Whitney { edge }) => (0, edges[edge].to_vec(), 0),
            Some(Role::Gradient { edge }) => (1, edges[edge].to_vec(), 0),
            Some(Role::Face { face, which }) => (2, faces[face].to_vec(), which),
            None => (3, vec![], 0),
        };
        let pairs: Vec<(usize, usize)> =
            (0..20).flat_map(|i| (0..20).map(move |j| (i, j))).filter(|&(i, j)| key(&ea, &fa, i) == key(&eb, &fb, j)).collect();
        assert_eq!(pairs.len(), 3 * 2 + 2, "three edges of two unknowns and the face's two");
        assert_eq!((fa[3], fb[0]), ([12, 31, 40], [12, 31, 40]));
        let normal = cross(sub(xyz(31), xyz(12)), sub(xyz(40), xyz(12)));
        for (u, v, w) in [(0.2, 0.3, 0.5), (0.6, 0.2, 0.2), (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)] {
            let x = combine([(u, xyz(12)), (v, xyz(31)), (w, xyz(40))]);
            let (va, vb) = (a.values(a.barycentric(x)), b.values(b.barycentric(x)));
            for &(i, j) in &pairs {
                let (ta, tb) = (cross(normal, va[i]), cross(normal, vb[j]));
                assert!(close(ta, tb, 1e-12), "{i}/{j}: {ta:?} {tb:?}");
            }
            // The functions of b not on the face have no tangential component on it.
            for (k, f) in vb.iter().enumerate() {
                if !pairs.iter().any(|(_, j)| *j == k) {
                    assert!(close(cross(normal, *f), [0.0; 3], 1e-12), "b {k}");
                }
            }
        }
    }

    #[test]
    fn a_left_handed_sorted_tetrahedron_has_the_same_matrices() {
        // Swapping two corners makes the tetrahedron left-handed: the volume and the matrices are
        // those of the same functions with edges renumbered and signs flipped where an edge
        // reverses. Here corners 1 and 2 swap: edges (0,1) ↔ (0,2), (1,3) ↔ (2,3), and (1,2)
        // reverses.
        let c = tet();
        let (r, l) = (EdgeTet::new(Order::One, c).unwrap(), EdgeTet::new(Order::One, [c[0], c[2], c[1], c[3]]).unwrap());
        assert!((r.volume() - l.volume()).abs() < 1e-15);
        let map = [(0, 1, 1.0), (1, 0, 1.0), (2, 2, 1.0), (3, 3, -1.0), (4, 5, 1.0), (5, 4, 1.0)];
        let (sr, sl, mr, ml) = (r.curl_curl(ONE).unwrap(), l.curl_curl(ONE).unwrap(), r.mass(ONE).unwrap(), l.mass(ONE).unwrap());
        for (p, pl, sp) in map {
            for (q, ql, sq) in map {
                assert!((sr.get(p, q) - sp * sq * sl.get(pl, ql)).abs() < 1e-12, "S {p} {q}");
                assert!((mr.get(p, q) - sp * sq * ml.get(pl, ql)).abs() < 1e-14, "M {p} {q}");
            }
        }
    }

    #[test]
    fn bad_tetrahedra_and_labels_are_refused() {
        let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        assert!(matches!(EdgeTet::new(Order::One, flat), Err(Error::Degenerate(_))));
        let point = [[0.5; 3]; 4];
        assert!(matches!(EdgeTet::new(Order::Two, point), Err(Error::Degenerate(_))));
        for bad in [f64::NAN, f64::INFINITY, 1e300] {
            let mut c = tet();
            c[2][1] = bad;
            assert!(matches!(EdgeTet::new(Order::One, c), Err(Error::Invalid(_))), "{bad}");
        }
        assert!(matches!(sort_corners([3, 9, 3, 1]), Err(Error::Invalid(_))));
        assert_eq!(sort_corners([40, 5, 31, 12]).unwrap(), [1, 3, 2, 0]);
        assert!(EdgeTet::of_labels(Order::One, tet(), [1, 2, 2, 4]).is_err());
        assert!(entity_keys([7, 7, 1, 2]).is_err());
        assert_eq!(entity_keys([40, 5, 31, 12]).unwrap().0[0], [5, 12]);
        // A tiny but sound tetrahedron is accepted: the criterion does not depend on size.
        let tiny = tet().map(|p| scale(p, 1e-9));
        assert!(EdgeTet::new(Order::Two, tiny).is_ok());
        // Edges of 1e-155, 1e-155 and 1e-5 pass it too, but ∇λ₁ × ∇λ₂ overflows: refused, not a
        // matrix of infinities.
        let thin = [[0.0; 3], [1e-155, 0.0, 0.0], [0.0, 1e-155, 0.0], [0.0, 0.0, 1e-5]];
        for order in [Order::One, Order::Two] {
            assert!(matches!(EdgeTet::new(order, thin), Err(Error::Degenerate(_))), "{order:?}");
            assert!(matches!(EdgeTet::of_labels(order, thin, [1, 2, 3, 4]), Err(Error::Degenerate(_))), "{order:?}");
            // A sound element whose material makes its matrix overflow (each curl coefficient is
            // of order 1 here, so ν · c · V · c and their sum pass f64::MAX): refused too.
            let e = EdgeTet::new(order, tet()).unwrap();
            assert!(matches!(e.curl_curl([f64::MAX; 3]), Err(Error::Degenerate(_))), "{order:?}");
        }
    }

    #[test]
    fn whether_a_tetrahedron_is_accepted_does_not_depend_on_its_labels() {
        // Sound at o (normalized Jacobian 1), nearly flat at a (8.3e-13 < 1e-12): the meshers'
        // gate takes corner 0 in mesh order, so the element must too, whichever node is lowest.
        let (o, a, b, c) = ([0.0; 3], [1.1, 0.0, 0.0], [0.0, 1e-6, 0.0], [0.0, 0.0, 1e-6]);
        let j = |p: [Vec3; 4]| fieldcraft_model::normalized_jacobian(&fieldcraft_model::tet4_jacobian(&p)).abs();
        assert!(j([o, a, b, c]) > fieldcraft_model::MIN_NORMALIZED_JACOBIAN);
        assert!(j([a, o, b, c]) < fieldcraft_model::MIN_NORMALIZED_JACOBIAN);
        let volumes: Vec<f64> = [[2, 1, 3, 4], [1, 2, 3, 4], [4, 3, 2, 1]]
            .into_iter()
            .map(|labels| EdgeTet::of_labels(Order::Two, [o, a, b, c], labels).unwrap().0.volume())
            .collect();
        assert!(volumes.iter().all(|v| (v - 1.1e-12 / 6.0).abs() < 1e-24), "{volumes:?}");
        // Refused in mesh order whatever the labels: flat at its own corner 0.
        assert!(matches!(EdgeTet::of_labels(Order::One, [a, o, b, c], [1, 2, 3, 4]), Err(Error::Degenerate(_))));
        // Sound in mesh order (determinant 3d², exact), but listed from (1, 1, 1) its Jacobian's
        // determinant rounds to exactly 0: the volume and the gradients are taken in mesh order.
        let d = 2f64.powi(-28);
        let xyz = [[0.0; 3], [1.0; 3], [-d, d, 0.0], [0.0, -d, d]];
        assert_eq!(det3(&fieldcraft_model::tet4_jacobian(&[xyz[1], xyz[0], xyz[2], xyz[3]])), 0.0);
        for order in [Order::One, Order::Two] {
            let (e, _) = EdgeTet::of_labels(order, xyz, [1, 2, 3, 4]).unwrap();
            let (f, sorted) = EdgeTet::of_labels(order, xyz, [2, 1, 3, 4]).unwrap();
            assert_eq!((e.volume(), f.volume()), (d * d / 2.0, d * d / 2.0));
            // Each gradient is its corner's: Σ ∇λᵢ = 0 and ∇λᵢ · (xⱼ − x₀) = δᵢⱼ − δᵢ₀ in either order.
            assert_eq!(sorted, [1, 0, 2, 3]);
            for (g, corners) in [(e.grads(), e.corners()), (f.grads(), f.corners())] {
                for (i, gi) in g.into_iter().enumerate() {
                    for (j, xj) in corners.into_iter().enumerate() {
                        let expected = f64::from(u8::from(i == j)) - f64::from(u8::from(i == 0));
                        assert!((dot(gi, sub(xj, corners[0])) - expected).abs() < 1e-6, "{order:?} {i} {j}");
                    }
                }
            }
            assert!(f.curl_curl([1.0; 3]).is_ok() && f.mass([1.0; 3]).is_ok());
        }
    }

    #[test]
    fn a_layout_numbers_the_elements_unknowns_and_its_gradient_is_in_the_kernel() {
        // Two tetrahedra sharing a face, declared to a layout: its slots give each element's
        // unknowns in local order, and the assembled curl–curl matrix annihilates the discrete
        // gradient of every potential (RFC §8.2), which a wrong orientation or slot would break.
        let xyz = |n: u32| match n {
            2 => [0.0, 0.0, 0.0],
            8 => [1.0, 0.1, 0.0],
            3 => [0.1, 1.0, 0.2],
            5 => [0.2, 0.3, 1.0],
            _ => [1.0, 1.0, 1.1],
        };
        let tets = [[5, 2, 8, 3], [3, 9, 8, 5]];
        for order in [Order::One, Order::Two] {
            let mut b = DofLayout::builder();
            for t in tets {
                declare(&mut b, t, order).unwrap();
            }
            let layout = b.build();
            assert_eq!((layout.edges(), layout.faces()), (9, if order == Order::Two { 7 } else { 0 }));
            let n = layout.dofs();
            let mut k = vec![0.0; n * n];
            for t in tets {
                let s = slots(&layout, t, order).unwrap();
                assert_eq!(s.len() * order.kinds().len(), order.functions());
                let dofs = layout.element_dofs(&s, order.kinds()).unwrap();
                let (e, _) = EdgeTet::of_labels(order, t.map(xyz), t).unwrap();
                let m = e.curl_curl(ONE).unwrap();
                for (i, r) in dofs.iter().enumerate() {
                    for (j, c) in dofs.iter().enumerate() {
                        k[r * n + c] += m.get(i, j);
                    }
                }
            }
            let g = gradient(&layout);
            let mut columns: Vec<Potential> = g.iter().map(|t| t.1).collect();
            columns.sort();
            columns.dedup();
            // Order 1: the five nodes; order 2: and the nine edges.
            assert_eq!(columns.len(), if order == Order::Two { 14 } else { 5 });
            let scale = k.iter().fold(0.0f64, |a, b| a.max(b.abs()));
            for column in columns {
                let mut x = vec![0.0; n];
                for (r, _, v) in g.iter().filter(|t| t.1 == column) {
                    x[*r] += v;
                }
                for row in 0..n {
                    let kx: f64 = (0..n).map(|c| k[row * n + c] * x[c]).sum();
                    assert!(kx.abs() < 1e-13 * scale, "{order:?} {column:?} row {row}: {kx}");
                }
            }
        }
        // A layout that lacks an entity gives no slots.
        let mut b = DofLayout::builder();
        declare(&mut b, tets[0], Order::One).unwrap();
        assert_eq!(slots(&b.build(), tets[1], Order::One), None);
        assert!(declare(&mut DofLayout::builder(), [1, 1, 2, 3], Order::One).is_err());
        // Orders 1 and 2 sharing a face: the order-1 element is refused (its shared edges carry
        // a second unknown it has no function for); the order-2 one finds its own.
        let mut b = DofLayout::builder();
        declare(&mut b, tets[0], Order::Two).unwrap();
        declare(&mut b, tets[1], Order::One).unwrap();
        b.edge([2, 8], DofSet::TEMPERATURE);
        let mixed = b.build();
        assert_eq!(slots(&mixed, tets[1], Order::One), None);
        assert_eq!(slots(&mixed, tets[0], Order::Two).map(|s| s.len()), Some(10));
        // Another kind on a slot (the temperature on edge 2–8) is no concern of the element.
        let mut b = DofLayout::builder();
        declare(&mut b, tets[0], Order::One).unwrap();
        b.edge([2, 8], DofSet::TEMPERATURE);
        assert_eq!(slots(&b.build(), tets[0], Order::One).map(|s| s.len()), Some(6));
    }
}
