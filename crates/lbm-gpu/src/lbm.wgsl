// The lattice Boltzmann step of fieldcraft_lbm::flow on the GPU: the same D3Q19 two-relaxation-
// time collision with Guo forcing, pull streaming with the per-link rules the CPU prepared
// (bounce-back with Ladd's moving-wall term, specular reflection, open placeholders), the Zou–He
// open boundaries, and the momentum exchange per force group. Populations are stored direction-
// major. The default AA pattern uses one population buffer (Bailey et al., ICPP 2009;
// Wittmann et al., Comput. Math. Appl. 65 (2013), 924–935). Even steps read a cell and
// write its opposite directions; odd steps read opposite directions at the upstream
// neighbours and write each outgoing direction at the downstream neighbour. Each thread
// writes exactly the addresses it read, so distinct cells never race. Physical coordinates
// wrap; the existing per-link rules still supply walls, reflection and open placeholders.
// A gather pass saves those boundary values before a scatter pass overrides their targets:
// overwriting in the gather itself could destroy another link's post-collision source.
// The optional A–B diagnostic path retains the original two-buffer collide/pull scheme.
//
// Single precision, with every population stored as its deviation from the weight, h_q = f_q − w_q
// (Krüger et al. 2017, §13.3.4: the deviations are small, so f32 keeps their digits instead of
// the weight's): ρ = 1 + Σ h, momentum Σ h c (Σ w c = 0), and every rule below is written for h.
// The CPU scheme (double precision, fixed sums) is the reference; see the crate's tests for the
// tolerances.
//
// With heat, a D3Q7 scalar lattice (the rest direction and the six axes of D3Q19, weights ¼ and
// ⅛) rides on the flow in the buffers t and t2, stored as deviations from w_q T₀ (T₀ the initial
// temperature) so that f32 keeps the digits of the differences: θ = T − T₀ = Σ t. It collides
// (two relaxation times, as the CPU) with the velocity the momentum collision left in `moments`, streams with the links' rules
// read through the group's thermal condition (anti-bounce-back for a held temperature, bounce-
// back for an adiabatic face, a flux added to the bounced population), takes its inlet
// temperatures and outlet neighbours, and gives each group the heat its links exchanged (in
// contrib.w and forces.w). The temperature feeds back as the buoyant acceleration −β g θ.
//
// Conducting solids (fieldcraft_lbm::conduction) step as the CPU's SolidHeat does, in f32 with
// temperatures as θ = T − T₀: `solid_ab` sums each conducting cell's fluid links (their heat is
// a − b θ_s), the sub-steps alternate between the two halves of the temperatures (`solid_sub0`,
// `solid_sub1`; a sub-step needs the whole previous one, so each is its own dispatch) with the
// links implicit and the faces explicit, and `solid_scatter` sets each link's population at the
// solid's mean temperature over the step, so the heat the fluid loses is what the solid gains.
//
// Smagorinsky's eddy viscosity (fieldcraft_lbm::turbulence) sets each cell's own relaxation rate
// from its non-equilibrium flux, τ = ½(τ₀ + √(τ₀² + 18√2 C² |P| / ρ)): in the regularised
// collision the rate the flux relaxes with, in the two-relaxation-time one τ⁺ (τ⁻ follows,
// keeping Λ). It is local, so it needs no neighbour and no buffer. With heat the scalar lattice
// needs the cell's ν_t too: the momentum collision leaves it in `moments`' w (the density there
// is only read after `cell_moments` has rewritten it), so it costs no memory and no binding.
// The chunk's largest ν_t goes to stats[3].

override AA_MODE: bool = true;
// The recursive regularised collision instead of the two-relaxation-time one
// (fieldcraft_lbm::regularized).
override REGULARIZED: bool = false;
// Smagorinsky's eddy viscosity (its constant is in `params.smag`). An override rather than a
// uniform test: with it off the kernels compile to what they were without the model, so their
// results stay bit-for-bit (a branch on a uniform changed the driver's roundings of the rest).
override SMAGORINSKY: bool = false;
// The largest lattice speed (fieldcraft_lbm::setup::MAX_LATTICE_SPEED): the cap of a convective
// outlet's convection speed.
override MAX_SPEED: f32 = 0.5;

const Q: u32 = 19u;
const QT: u32 = 7u;
const WG: u32 = 256u;
const KIND_PULL: u32 = 0u;
const KIND_BOUNCE: u32 = 1u;
const KIND_REFLECT: u32 = 2u;
const KIND_OPEN: u32 = 3u;
// Set in a link's first word when its cell is an inlet or outlet cell, whose scalar populations
// the open boundary sets after streaming: the heat the link gives never stays in the fluid, and
// the wall heat leaves it out (fieldcraft_lbm::flow::ThermalState::resets).
const RESET: u32 = 0x1000000u;

struct Params {
    n: u32,
    nx: u32,
    ny: u32,
    nz: u32,
    // The row width of every dispatch wider than one row (MAX_GROUPS_X in native.rs).
    groups_x: u32,
    links: u32,
    open: u32,
    groups: u32,
    omega_plus: f32,
    omega_minus: f32,
    s_plus: f32,
    s_minus: f32,
    // The body acceleration (xyz) and whether there is one (w).
    acc: vec4<f32>,
    // The buoyancy −β g in lattice units (xyz) and whether there is any (w).
    buoy: vec4<f32>,
    // The scalar lattice: its antisymmetric relaxation rate 1/τ_T⁻ (x), T₀ (y), whether there
    // is heat (z), its symmetric relaxation rate 1/τ_T⁺ (w).
    therm: vec4<f32>,
    // Current AA orientation (x: 0 normal, 1 opposite directions at upstream neighbours).
    phase: vec4<u32>,
    // Conducting solids: cells m (x), interface links k (y), padding (z), and 1 / sub-steps (w,
    // as f32 bits).
    solid: vec4<u32>,
    // Smagorinsky's model: C² (x; 0 without it), 1 / Pr_t (y), the two relaxation times' Λ (z).
    smag: vec4<f32>,
};

// Each kernel's bind group holds only the bindings it reaches (native.rs reads them from this
// shader), so no kernel binds more storage buffers than WebGPU's default of 8 per stage; all of
// them in one layout would need 22, more than GL allows (16, #640).
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> fluid: array<u32>;
@group(0) @binding(2) var<storage, read_write> f: array<f32>;
@group(0) @binding(3) var<storage, read_write> g: array<f32>;
// Per link: (q | kind << 8 | q2 << 16 | RESET, from, group, cell) and the moving-wall term.
@group(0) @binding(5) var<storage, read> link_u: array<vec4<u32>>;
// Per convective outlet cell (fieldcraft_lbm::outlet), OUTLET floats at its slot: its five unknown
// populations after the last step (h, 0..5) and its convection speed (5), kept between steps;
// this step's convected unknowns (6..11), the density shift (11) and whether it convects (12),
// from `outlet_convect` for `open_boundaries`.
@group(0) @binding(4) var<storage, read_write> outlet: array<f32>;
@group(0) @binding(6) var<storage, read> link_f: array<f32>;
// Per open cell: (cell, axis, kind: 0 velocity / 1 density / 2 convective outlet, its slot in `outlet`) and (target xyz or density in x, inward sign in w).
@group(0) @binding(7) var<storage, read> open_u: array<vec4<u32>>;
@group(0) @binding(8) var<storage, read> open_f: array<vec4<f32>>;
// Per cell: velocity (with the half force) and density, from the moments kernel.
@group(0) @binding(9) var<storage, read_write> moments: array<vec4<f32>>;
// [0]: the largest squared lattice speed of the chunk (as f32 bits); [1]: a density was not positive;
// [2]: a temperature was not a number; [3]: the chunk's largest eddy viscosity (f32 bits).
@group(0) @binding(10) var<storage, read_write> stats: array<atomic<u32>>;
// Per link, the momentum it transferred to its group in the last step.
@group(0) @binding(11) var<storage, read_write> contrib: array<vec4<f32>>;
// Link indices ordered by group, and each group's range in that order (groups + 1 entries).
@group(0) @binding(12) var<storage, read> by_group: array<u32>;
@group(0) @binding(13) var<storage, read> group_start: array<u32>;
@group(0) @binding(14) var<storage, read_write> forces: array<vec4<f32>>;
// Per workgroup, the mass of its cells.
@group(0) @binding(15) var<storage, read_write> partials: array<f32>;
// The scalar populations remain A–B (t[q n + cell]; see the notes above).
@group(0) @binding(16) var<storage, read_write> t: array<f32>;
@group(0) @binding(17) var<storage, read_write> t2: array<f32>;
// Per force group: its thermal condition (x: 0 adiabatic, 1 temperature, 2 flux) and the value
// (y: the held temperature − T₀, or the flux as a temperature per step).
@group(0) @binding(18) var<storage, read> group_t: array<vec4<f32>>;
// Per open cell: whether it holds a temperature (x) and that temperature − T₀ (y).
@group(0) @binding(19) var<storage, read> open_t: array<vec4<f32>>;
// Conducting solids, m cells and k interface links, in three buffers whose layouts the `s_*`
// functions below name: u32 data (cells, each cell's first link, the links, the faces), f32
// values (capacities, sources, the faces' conductances and values, the links' weights) and the
// state (θ, the links' sums, the mean θ, the heat into each box face, the links' terms).
@group(0) @binding(20) var<storage, read> solid_u: array<u32>;
@group(0) @binding(21) var<storage, read> solid_f: array<f32>;
@group(0) @binding(22) var<storage, read_write> solid_s: array<f32>;

var<workgroup> red: array<vec4<f32>, 256>;
// The largest squared speed of a workgroup's cells (collide): one atomic per workgroup, not
// per cell, which serialise on a single address.
var<workgroup> speed: array<f32, 256>;
// The largest eddy viscosity of a workgroup's cells, likewise.
var<workgroup> eddy: array<f32, 256>;

fn c_of(q: u32) -> vec3<i32> {
    switch q {
        case 0u: { return vec3<i32>(0, 0, 0); }
        case 1u: { return vec3<i32>(1, 0, 0); }
        case 2u: { return vec3<i32>(-1, 0, 0); }
        case 3u: { return vec3<i32>(0, 1, 0); }
        case 4u: { return vec3<i32>(0, -1, 0); }
        case 5u: { return vec3<i32>(0, 0, 1); }
        case 6u: { return vec3<i32>(0, 0, -1); }
        case 7u: { return vec3<i32>(1, 1, 0); }
        case 8u: { return vec3<i32>(-1, -1, 0); }
        case 9u: { return vec3<i32>(1, -1, 0); }
        case 10u: { return vec3<i32>(-1, 1, 0); }
        case 11u: { return vec3<i32>(1, 0, 1); }
        case 12u: { return vec3<i32>(-1, 0, -1); }
        case 13u: { return vec3<i32>(1, 0, -1); }
        case 14u: { return vec3<i32>(-1, 0, 1); }
        case 15u: { return vec3<i32>(0, 1, 1); }
        case 16u: { return vec3<i32>(0, -1, -1); }
        case 17u: { return vec3<i32>(0, 1, -1); }
        default: { return vec3<i32>(0, -1, 1); }
    }
}

fn opp(q: u32) -> u32 {
    if q == 0u { return 0u; }
    if (q & 1u) == 1u { return q + 1u; }
    return q - 1u;
}

fn w_of(q: u32) -> f32 {
    if q == 0u { return 1.0 / 3.0; }
    if q < 7u { return 1.0 / 18.0; }
    return 1.0 / 36.0;
}

fn wt_of(q: u32) -> f32 {
    if q == 0u { return 0.25; }
    return 0.125;
}

fn has_heat() -> bool {
    return params.therm.z > 0.5;
}

// The temperature of a cell above T₀: the sum of its scalar populations.
fn theta_of(cell: u32) -> f32 {
    var s = 0.0;
    for (var q = 0u; q < QT; q = q + 1u) {
        s = s + t[q * params.n + cell];
    }
    return s;
}

// The body acceleration of a cell: the constant one plus the buoyant one.
fn acc_of(cell: u32) -> vec3<f32> {
    var a = params.acc.xyz;
    if params.buoy.w > 0.5 { a = a + params.buoy.xyz * theta_of(cell); }
    return a;
}

fn forced() -> bool {
    return params.acc.w > 0.5 || params.buoy.w > 0.5;
}

// The scalar equilibrium w_q T (1 + 4 c_q·u), shifted by w_q T₀ like the populations.
fn eq_t(q: u32, theta: f32, u: vec3<f32>) -> f32 {
    let cu = 4.0 * dot(vec3<f32>(c_of(q)), u);
    return wt_of(q) * (theta * (1.0 + cu) + params.therm.y * cu);
}

fn cell_of(lid: u32, wg: vec3<u32>) -> u32 {
    return (wg.y * params.groups_x + wg.x) * WG + lid;
}

fn is_fluid(cell: u32) -> bool {
    return fluid[cell] != 0u;
}

// A cell's grid coordinates (integer divisions once per cell, not once per direction).
fn coords(cell: u32) -> vec3<i32> {
    let plane = params.nx * params.ny;
    return vec3<i32>(i32(cell % params.nx), i32((cell % plane) / params.nx), i32(cell / plane));
}

// The cell at `p`, wrapped periodically. A lattice step moves at most one cell per axis, so a
// compare and an add per axis replace the remainder.
fn wrapped(p: vec3<i32>) -> u32 {
    let size = vec3<i32>(i32(params.nx), i32(params.ny), i32(params.nz));
    var w = select(p, p + size, p < vec3<i32>(0));
    w = select(w, w - size, w >= size);
    return u32((w.z * size.y + w.y) * size.x + w.x);
}

// Where population q of the cell at `xyz` is stored in the given AA orientation: in place when
// even, reversed at the upstream neighbour when odd.
fn population_index(q: u32, cell: u32, xyz: vec3<i32>, phase: u32) -> u32 {
    if phase == 0u { return q * params.n + cell; }
    return opp(q) * params.n + wrapped(xyz - c_of(q));
}

fn population(q: u32, cell: u32, xyz: vec3<i32>) -> f32 {
    if AA_MODE { return f[population_index(q, cell, xyz, params.phase.x)]; }
    return f[q * params.n + cell];
}

// Collision's output at its source cell, before the boundary rules.
fn post_index(q: u32, cell: u32, xyz: vec3<i32>) -> u32 {
    if !AA_MODE { return q * params.n + cell; }
    if params.phase.x == 0u { return opp(q) * params.n + cell; }
    return q * params.n + wrapped(xyz + c_of(q));
}

fn post(q: u32, cell: u32) -> f32 {
    return f[post_index(q, cell, aa_coords(cell))];
}

fn streamed(q: u32, cell: u32, xyz: vec3<i32>) -> f32 {
    if AA_MODE { return f[population_index(q, cell, xyz, params.phase.x ^ 1u)]; }
    return g[q * params.n + cell];
}

fn store_streamed(q: u32, cell: u32, xyz: vec3<i32>, v: f32) {
    if AA_MODE { f[population_index(q, cell, xyz, params.phase.x ^ 1u)] = v; }
    else { g[q * params.n + cell] = v; }
}

// The coordinates the AA addressing needs; none in the two-buffer layout.
fn aa_coords(cell: u32) -> vec3<i32> {
    if AA_MODE { return coords(cell); }
    return vec3<i32>(0);
}

// Two-relaxation-time collision with Guo forcing, in place in f (Ginzburg et al. 2008; Guo et
// al. 2002; Krüger et al. 2017 ch. 10), and the chunk's largest speed.
@compute @workgroup_size(256)
fn collide(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let cell = cell_of(lid, wg);
    speed[lid] = 0.0;
    eddy[lid] = 0.0;
    if cell < params.n && is_fluid(cell) {
        collide_cell(cell, lid);
    }
    workgroupBarrier();
    for (var k = 128u; k > 0u; k = k >> 1u) {
        if lid < k {
            speed[lid] = max(speed[lid], speed[lid + k]);
            eddy[lid] = max(eddy[lid], eddy[lid + k]);
        }
        workgroupBarrier();
    }
    if lid == 0u {
        atomicMax(&stats[0], bitcast<u32>(speed[0]));
        if turbulent() { atomicMax(&stats[3], bitcast<u32>(eddy[0])); }
    }
}

fn turbulent() -> bool {
    return SMAGORINSKY;
}

// The non-equilibrium flux P = Σ c c h − δρ c_s² I − ρ u u + ½ (u F + F u): its diagonal (xx, yy,
// zz) and off-diagonal (xy, xz, yz) parts.
struct Flux {
    d: vec3<f32>,
    o: vec3<f32>,
};

fn flux(fq: ptr<function, array<f32, 19>>, drho: f32, rho: f32, u: vec3<f32>, force: vec3<f32>, forced: bool) -> Flux {
    // (Opposite directions share c c, so by pairs.)
    var d = vec3<f32>(0.0);
    var o = vec3<f32>(0.0);
    for (var q = 1u; q < Q; q = q + 2u) {
        let c = vec3<f32>(c_of(q));
        let v = (*fq)[q] + (*fq)[q + 1u];
        d = d + c * c * v;
        o = o + vec3<f32>(c.x * c.y, c.x * c.z, c.y * c.z) * v;
    }
    let cs2 = 1.0 / 3.0;
    var pd = d - vec3<f32>(drho * cs2) - rho * u * u;
    var po = o - rho * vec3<f32>(u.x * u.y, u.x * u.z, u.y * u.z);
    if forced {
        pd = pd + u * force;
        po = po + 0.5 * vec3<f32>(u.x * force.y + force.x * u.y, u.x * force.z + force.x * u.z, u.y * force.z + force.y * u.z);
    }
    return Flux(pd, po);
}

// A cell's relaxation rate under Smagorinsky's model (fieldcraft_lbm::turbulence::omega): never
// faster than the fluid's `omega0`, and `omega0` where the flux is not a usable number.
fn smagorinsky_omega(p: Flux, rho: f32, omega0: f32) -> f32 {
    let tau0 = 1.0 / omega0;
    let pp = dot(p.d, p.d) + 2.0 * dot(p.o, p.o);
    let tau = 0.5 * (tau0 + sqrt(tau0 * tau0 + 18.0 * 1.4142135 * params.smag.x * sqrt(pp) / rho));
    if tau >= tau0 && tau < 1e30 { return 1.0 / tau; }
    return omega0;
}

// The eddy viscosity of a cell relaxing at `w` (c_s² = ⅓), kept for the reply and, with heat, for
// the scalar lattice.
fn record_eddy(cell: u32, lid: u32, w: f32, u: vec3<f32>) {
    let nu_t = max((1.0 / w - 1.0 / params.omega_plus) / 3.0, 0.0);
    eddy[lid] = nu_t;
    if has_heat() { moments[cell] = vec4<f32>(u, nu_t); }
}

fn collide_cell(cell: u32, lid: u32) {
    let n = params.n;
    // Only the odd AA orientation addresses neighbours.
    var xyz = vec3<i32>(0);
    if AA_MODE && params.phase.x == 1u { xyz = coords(cell); }
    var fq: array<f32, 19>;
    var drho = 0.0;
    var j = vec3<f32>(0.0);
    for (var q = 0u; q < Q; q = q + 1u) {
        let v = population(q, cell, xyz);
        fq[q] = v;
        drho = drho + v;
        j = j + vec3<f32>(c_of(q)) * v;
    }
    let rho = 1.0 + drho;
    if !(rho > 0.0) || rho != rho || abs(rho) > 1e30 {
        atomicStore(&stats[1], 1u);
        return;
    }
    let a = acc_of(cell);
    let force = rho * a;
    let u = (j + 0.5 * force) / rho;
    let uu = dot(u, u);
    speed[lid] = uu;
    if has_heat() { moments[cell] = vec4<f32>(u, rho); }
    let uf = dot(u, force);
    let forced = forced();
    if REGULARIZED {
        collide_regularized(cell, lid, xyz, &fq, drho, rho, u, force, forced);
        return;
    }
    if turbulent() {
        collide_trt_smagorinsky(cell, lid, xyz, &fq, drho, rho, u, force, forced);
        return;
    }
    let wp = params.omega_plus;
    let wm = params.omega_minus;
    let sp = params.s_plus;
    let sm = params.s_minus;
    // Rest population: symmetric only. Equilibria are shifted by the weight like the populations.
    let eq0 = w_of(0u) * (drho - 1.5 * rho * uu);
    var out0 = fq[0] - wp * (fq[0] - eq0);
    if forced { out0 = out0 + sp * w_of(0u) * (-3.0 * uf); }
    f[post_index(0u, cell, xyz)] = out0;
    for (var q = 1u; q < Q; q = q + 2u) {
        let o = q + 1u;
        let cv = vec3<f32>(c_of(q));
        let cu = dot(cv, u);
        let w = w_of(q);
        let eq_plus = w * (drho + rho * (4.5 * cu * cu - 1.5 * uu));
        let eq_minus = w * rho * 3.0 * cu;
        let plus = 0.5 * (fq[q] + fq[o]);
        let minus = 0.5 * (fq[q] - fq[o]);
        let relax_plus = wp * (plus - eq_plus);
        let relax_minus = wm * (minus - eq_minus);
        var s_plus = 0.0;
        var s_minus = 0.0;
        if forced {
            let cf = dot(cv, force);
            s_plus = sp * w * (9.0 * cu * cf - 3.0 * uf);
            s_minus = sm * w * 3.0 * cf;
        }
        f[post_index(q, cell, xyz)] = fq[q] - relax_plus - relax_minus + s_plus + s_minus;
        f[post_index(o, cell, xyz)] = fq[o] - relax_plus + relax_minus + s_plus - s_minus;
    }
}

// The two-relaxation-time collision with each cell's own rates under Smagorinsky's model: τ⁺
// from the eddy viscosity, τ⁻ following it so that Λ stays the fluid's.
fn collide_trt_smagorinsky(cell: u32, lid: u32, xyz: vec3<i32>, fq: ptr<function, array<f32, 19>>, drho: f32, rho: f32, u: vec3<f32>, force: vec3<f32>, forced: bool) {
    let wp = smagorinsky_omega(flux(fq, drho, rho, u, force, forced), rho, params.omega_plus);
    record_eddy(cell, lid, wp, u);
    let wm = 1.0 / (0.5 + params.smag.z / (1.0 / wp - 0.5));
    let sp = 1.0 - 0.5 * wp;
    let sm = 1.0 - 0.5 * wm;
    let uu = dot(u, u);
    let uf = dot(u, force);
    let eq0 = w_of(0u) * (drho - 1.5 * rho * uu);
    var out0 = (*fq)[0] - wp * ((*fq)[0] - eq0);
    if forced { out0 = out0 + sp * w_of(0u) * (-3.0 * uf); }
    f[post_index(0u, cell, xyz)] = out0;
    for (var q = 1u; q < Q; q = q + 2u) {
        let o = q + 1u;
        let cv = vec3<f32>(c_of(q));
        let cu = dot(cv, u);
        let w = w_of(q);
        let eq_plus = w * (drho + rho * (4.5 * cu * cu - 1.5 * uu));
        let eq_minus = w * rho * 3.0 * cu;
        let plus = 0.5 * ((*fq)[q] + (*fq)[o]);
        let minus = 0.5 * ((*fq)[q] - (*fq)[o]);
        let relax_plus = wp * (plus - eq_plus);
        let relax_minus = wm * (minus - eq_minus);
        var s_plus = 0.0;
        var s_minus = 0.0;
        if forced {
            let cf = dot(cv, force);
            s_plus = sp * w * (9.0 * cu * cf - 3.0 * uf);
            s_minus = sm * w * 3.0 * cf;
        }
        f[post_index(q, cell, xyz)] = (*fq)[q] - relax_plus - relax_minus + s_plus + s_minus;
        f[post_index(o, cell, xyz)] = (*fq)[o] - relax_plus + relax_minus + s_plus - s_minus;
    }
}

// The recursive regularised collision with Guo forcing (fieldcraft_lbm::regularized, whose notes
// derive it): every population rebuilt from ρ, u, the non-equilibrium flux
// P = Π − ρ (c_s² I + u u) + ½ (u F + F u) and the third-order moments D3Q19 holds (equilibrium
// ρ u_a u_b u_c, non-equilibrium u_a P_bc + u_b P_ac + u_c P_ab, as the orthogonal sums and
// differences of each pair), then f* = f^eq + (1 − ω) f^(1)[P] + ½ F_q. In deviations from the
// weights, Π − ρ c_s² I = Σ c c h − δρ c_s² I.
fn collide_regularized(cell: u32, lid: u32, xyz: vec3<i32>, fq: ptr<function, array<f32, 19>>, drho: f32, rho: f32, u: vec3<f32>, force: vec3<f32>, forced: bool) {
    let p = flux(fq, drho, rho, u, force, forced);
    let pd = p.d;
    let po = p.o;
    let trace = pd.x + pd.y + pd.z;
    var omega = params.omega_plus;
    if turbulent() {
        omega = smagorinsky_omega(p, rho, omega);
        record_eddy(cell, lid, omega, u);
    }
    let relax = 1.0 - omega;
    // Third-order coefficients (equilibrium plus (1 − ω) times the recursive part).
    let xyy = rho * u.x * u.y * u.y + relax * (2.0 * u.y * po.x + u.x * pd.y);
    let xzz = rho * u.x * u.z * u.z + relax * (2.0 * u.z * po.y + u.x * pd.z);
    let xxy = rho * u.x * u.x * u.y + relax * (2.0 * u.x * po.x + u.y * pd.x);
    let yzz = rho * u.y * u.z * u.z + relax * (2.0 * u.z * po.z + u.y * pd.z);
    let xxz = rho * u.x * u.x * u.z + relax * (2.0 * u.x * po.y + u.z * pd.x);
    let yyz = rho * u.y * u.y * u.z + relax * (2.0 * u.y * po.z + u.z * pd.y);
    let sums = 13.5 * vec3<f32>(xyy + xzz, xxy + yzz, xxz + yyz);
    let diffs = 4.5 * vec3<f32>(xyy - xzz, xxy - yzz, xxz - yyz);
    let uu = dot(u, u);
    let uf = dot(u, force);
    let third = 1.0 / 3.0;
    // The rest population, then each pair of opposite directions: the even terms are shared, the
    // odd ones (orders 1 and 3, the force's first) change sign.
    var rest = w_of(0u) * (drho - 1.5 * rho * uu - relax * 1.5 * trace);
    if forced { rest = rest - 1.5 * w_of(0u) * uf; }
    f[post_index(0u, cell, xyz)] = rest;
    for (var q = 1u; q < Q; q = q + 2u) {
        let c = vec3<f32>(c_of(q));
        let w = w_of(q);
        let cu = dot(c, u);
        // H_abb = c_a c_b² − c_a/3, paired as in fieldcraft_lbm::regularized.
        let hxyy = c.x * c.y * c.y - c.x * third;
        let hxzz = c.x * c.z * c.z - c.x * third;
        let hxxy = c.x * c.x * c.y - c.y * third;
        let hyzz = c.y * c.z * c.z - c.y * third;
        let hxxz = c.x * c.x * c.z - c.z * third;
        let hyyz = c.y * c.y * c.z - c.z * third;
        let h3 = dot(vec3<f32>(hxyy + hxzz, hxxy + hyzz, hxxz + hyyz), sums) + dot(vec3<f32>(hxyy - hxzz, hxxy - hyzz, hxxz - hyyz), diffs);
        let cpc = dot(c * c, pd) + 2.0 * dot(vec3<f32>(c.x * c.y, c.x * c.z, c.y * c.z), po);
        let second = 4.5 * (cpc - trace * third);
        var even = w * (drho + rho * (4.5 * cu * cu - 1.5 * uu) + relax * second);
        var odd = w * (rho * 3.0 * cu + h3);
        if forced {
            let cf = dot(c, force);
            even = even + 0.5 * w * (9.0 * cu * cf - 3.0 * uf);
            odd = odd + 1.5 * w * cf;
        }
        f[post_index(q, cell, xyz)] = even + odd;
        f[post_index(q + 1u, cell, xyz)] = even - odd;
    }
}

// Pull streaming of f into g, first pass: every population from the neighbour at −c_q
// (wrapped; a wall or an open face overrides it in the second pass).
@compute @workgroup_size(256)
fn stream(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let cell = cell_of(lid, wg);
    if cell >= params.n { return; }
    if !is_fluid(cell) { return; }
    let n = params.n;
    let nx = i32(params.nx);
    let ny = i32(params.ny);
    let nz = i32(params.nz);
    let i = i32(cell % params.nx);
    let jj = i32((cell / params.nx) % params.ny);
    let k = i32(cell / (params.nx * params.ny));
    for (var q = 0u; q < Q; q = q + 1u) {
        let c = c_of(q);
        let si = (i - c.x + nx) % nx;
        let sj = (jj - c.y + ny) % ny;
        let sk = (k - c.z + nz) % nz;
        let src = u32((sk * ny + sj) * nx + si);
        g[q * n + cell] = f[q * n + src];
    }
}

// Second pass, one thread per link: the populations with a rule (a periodic image, a wall, a
// slip face, an open face). Every (cell, direction) has at most one link, so no two threads
// write the same population.
@compute @workgroup_size(256)
fn stream_links(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    gather_link(cell_of(lid, wg), false);
}

// The same on the step whose wall forces are read back (a chunk's last): the AA layout
// overwrites the post-collision sources in the scatter, so the momentum exchange is taken here.
@compute @workgroup_size(256)
fn stream_links_exchange(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    gather_link(cell_of(lid, wg), true);
}

fn gather_link(l: u32, exchange: bool) {
    if l >= params.links { return; }
    let n = params.n;
    let lu = link_u[l];
    let q = lu.x & 0xffu;
    let kind = (lu.x >> 8u) & 0xffu;
    let q2 = (lu.x >> 16u) & 0xffu;
    let cell = lu.w;
    var v = 0.0;
    if kind == KIND_PULL {
        v = post(q, lu.y);
    } else if kind == KIND_BOUNCE {
        v = post(opp(q), cell) + link_f[l];
    } else if kind == KIND_REFLECT {
        // The mirrored direction has the same weight, so h reflects like f.
        v = post(q2, lu.y);
    } else {
        v = post(opp(q), cell);
    }
    if AA_MODE {
        g[l] = v;
        if exchange { contrib[l] = vec4<f32>(momentum_exchange(l, lu, q, kind, q2), 0.0); }
    } else {
        g[q * n + cell] = v;
    }
}

// The gather above must complete before any AA boundary population is overwritten.
@compute @workgroup_size(256)
fn scatter_links(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let l = cell_of(lid, wg);
    if l >= params.links { return; }
    let lu = link_u[l];
    store_streamed(lu.x & 0xffu, lu.w, aa_coords(lu.w), g[l]);
}

// Every streamed population of an open cell set to the equilibrium of `rho` and `bare` (the
// velocity without the half force), as deviations from the weights: what an outflow cell holds
// under the regularised collision (fieldcraft_lbm::flow::hold_equilibrium). Only where `apply`
// (see OPEN_CELLS_WALK_TOGETHER).
fn hold_equilibrium(cell: u32, xyz: vec3<i32>, rho: f32, bare: vec3<f32>, apply: bool) {
    let uu = dot(bare, bare);
    for (var q = 0u; q < Q; q = q + 1u) {
        let cu = dot(vec3<f32>(c_of(q)), bare);
        if apply {
            store_streamed(q, cell, xyz, w_of(q) * ((rho - 1.0) + rho * (3.0 * cu + 4.5 * cu * cu - 1.5 * uu)));
        }
    }
}

// OPEN_CELLS_WALK_TOGETHER. Every thread of the open-cell kernels below walks the same loops
// over the directions, whatever its cell's kind: the kind decides which values a thread reads and
// stores, never whether a loop that touches a buffer runs. lavapipe before Mesa 26.2 crashes on
// the other shape. It runs threads in SIMD groups; a value set inside a loop that only some
// threads of a group entered, the group's first not among them, is read for that first thread
// wherever the group's mask is empty, and when the value is a buffer's descriptor address it is
// unset there, so the access dereferences null (fixed by Mesa commit 63b89856, "gallivm: don't
// deref a NULL buffer descriptor with an empty exec mask"). Zou–He run only by the convective
// cells whose flow came back in crashed every CI run on a software GPU (Ubuntu 24.04, Mesa
// 25.2.8). Threads past the last open cell may return early: they are the last of their group.

const OUTLET: u32 = 16u;
// The open-cell kind whose closure a convective outlet cell takes where the flow comes back in:
// the uniform pressure outlet's (fieldcraft_lbm::outlet::backflow).
const BACKFLOW_KIND: u32 = 1u;

// Whether direction q enters an open cell of normal `axis` (inward sign `sign`) from outside.
fn is_unknown(q: u32, axis: u32, sign: f32) -> bool {
    return f32(c_of(q)[axis]) * sign > 0.0;
}

// The index step of `axis` and the cell one layer in from an open cell.
fn inside_of(cell: u32, axis: u32, sign: f32) -> u32 {
    var step = 1u;
    if axis == 1u { step = params.nx; }
    if axis == 2u { step = params.nx * params.ny; }
    if sign > 0.0 { return cell + step; }
    return cell - step;
}

// The outward velocity (with the half force) of a cell of an open face from its density and
// normal momentum.
fn outward(rho: f32, jn: f32, cell: u32, axis: u32, sign: f32) -> f32 {
    if !(rho > 0.0) { return 0.0; }
    return -sign * (jn + 0.5 * rho * acc_of(cell)[axis]) / rho;
}

// The convective outlet (fieldcraft_lbm::outlet), first pass, one thread per open cell: the
// unknown populations convected from the cell's own after the last step and the streamed ones of
// the cell inside, (own + λ inside) / (1 + λ); the density shift to the face's density, and
// whether the flow leaves. Reads only populations, so every cell sees them before any open
// boundary changes them.
// The other open cells walk along (OPEN_CELLS_WALK_TOGETHER) and store nothing.
@compute @workgroup_size(256)
fn outlet_convect(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let o = cell_of(lid, wg);
    if o >= params.open { return; }
    let ou = open_u[o];
    let convective = ou.z == 2u;
    let cell = ou.x;
    let axis = ou.y;
    let sign = open_f[o].w;
    let xyz = aa_coords(cell);
    let inside = inside_of(cell, axis, sign);
    let has_inside = is_fluid(inside);
    var inside_xyz = vec3<i32>(0);
    if has_inside { inside_xyz = aa_coords(inside); }
    let base = ou.w * OUTLET;
    var lambda = 0.0;
    if convective { lambda = outlet[base + 5u]; }
    var k = 0u;
    var drho = 0.0;
    var jn = 0.0;
    for (var q = 0u; q < Q; q = q + 1u) {
        var v = 0.0;
        if !is_unknown(q, axis, sign) {
            v = streamed(q, cell, xyz);
        } else {
            if convective {
                v = outlet[base + k];
                if has_inside { v = (v + lambda * streamed(q, inside, inside_xyz)) / (1.0 + lambda); }
                outlet[base + 6u + k] = v;
            }
            k = k + 1u;
        }
        drho = drho + v;
        jn = jn + f32(c_of(q)[axis]) * v;
    }
    if convective {
        outlet[base + 11u] = (open_f[o].x - 1.0) - drho;
        outlet[base + 12u] = select(0.0, 1.0, outward(1.0 + drho, jn, cell, axis, sign) > 0.0);
    }
}

// Zou–He velocity and pressure boundaries on the streamed populations in g (Zou & He 1997;
// Hecht & Harting 2010 for the D3Q19 transverse correction), and the second pass of the
// convective outlet: where the flow leaves, the convected unknowns and every population shifted
// by w_q times the density shift; where it comes back in, the Zou–He pressure outlet. A
// convective cell then keeps its unknowns and its outward velocity for the next step. Every
// cell walks through both (OPEN_CELLS_WALK_TOGETHER): Zou–He stores nothing where the outlet
// convects, and Zou–He leaves the known populations the convection reads as they are (the
// regularised hold, which sets them all, is not applied there either).
@compute @workgroup_size(256)
fn open_boundaries(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let o = cell_of(lid, wg);
    if o >= params.open { return; }
    let ou = open_u[o];
    let cell = ou.x;
    let axis = ou.y;
    let sign = open_f[o].w;
    let xyz = aa_coords(cell);
    let base = ou.w * OUTLET;
    let convective = ou.z == 2u;
    var convects = false;
    var shift = 0.0;
    if convective {
        convects = outlet[base + 12u] > 0.5;
        shift = outlet[base + 11u];
    }
    zou_he(o, select(ou.z, BACKFLOW_KIND, convective), cell, axis, sign, xyz, !convects);
    var k = 0u;
    var drho = 0.0;
    var jn = 0.0;
    for (var q = 0u; q < Q; q = q + 1u) {
        let unknown = is_unknown(q, axis, sign);
        var v = streamed(q, cell, xyz);
        if convects {
            if unknown { v = outlet[base + 6u + k]; }
            v = v + w_of(q) * shift;
            store_streamed(q, cell, xyz, v);
        }
        if unknown {
            if convective { outlet[base + k] = v; }
            k = k + 1u;
        }
        drho = drho + v;
        jn = jn + f32(c_of(q)[axis]) * v;
    }
    if convective {
        outlet[base + 5u] = clamp(outward(1.0 + drho, jn, cell, axis, sign), 0.0, MAX_SPEED);
    }
}

// Zou–He on open cell `o`: the velocity (kind 0) or the density (kind 1) of its target, stored
// only where `apply` (OPEN_CELLS_WALK_TOGETHER).
fn zou_he(o: u32, kind: u32, cell: u32, axis: u32, sign: f32, xyz: vec3<i32>, apply: bool) {
    let ot = open_f[o];
    let a_half = 0.5 * acc_of(cell);
    // The sums of h over the tangential (c_n = 0) and outgoing (c_n < 0) directions. The sums of
    // f add the weights' (2/3 and 1/6 for D3Q19, so tangential + 2 outgoing gains exactly 1): the
    // weight is added after the deviations are summed, not before, which in single precision would
    // round away the deviations' last digits at low speed.
    var tangential = 0.0;
    var outgoing = 0.0;
    for (var q = 0u; q < Q; q = q + 1u) {
        let v = streamed(q, cell, xyz);
        let cn = f32(c_of(q)[axis]) * sign;
        if cn == 0.0 {
            tangential = tangential + v;
        } else if cn < 0.0 {
            outgoing = outgoing + v;
        }
    }
    var rho = 0.0;
    var bare = vec3<f32>(0.0);
    if kind == 0u {
        bare = ot.xyz - a_half;
        let un = sign * bare[axis];
        rho = (1.0 + (tangential + 2.0 * outgoing)) / (1.0 - un);
    } else {
        rho = ot.x;
        let un = ((rho - 1.0) - (tangential + 2.0 * outgoing)) / rho;
        // Not `bare[axis] = …`: DX12's shader compiler (FXC) refuses a write through a
        // dynamic index into a local vector (X3500).
        bare = select(-a_half, vec3<f32>(sign * un), vec3<u32>(0u, 1u, 2u) == vec3<u32>(axis));
    }
    // With the regularised collision a density target holds the equilibrium of its density and
    // the normal velocity just found (fieldcraft_lbm::flow::zou_he): a pressure outlet, and a
    // convective outlet's cell where the flow comes back in (BACKFLOW_KIND).
    let hold = REGULARIZED && kind == 1u;
    hold_equilibrium(cell, xyz, rho, bare, apply && hold);
    let closure = apply && !hold;
    for (var q = 1u; q < Q; q = q + 1u) {
        let c = c_of(q);
        if closure && f32(c[axis]) * sign > 0.0 {
            let cu = dot(vec3<f32>(c), bare);
            store_streamed(q, cell, xyz, streamed(opp(q), cell, xyz) + 6.0 * w_of(q) * rho * cu);
        }
    }
    // The tangential momentum along each axis but the normal one (the normal's pass stores
    // nothing).
    for (var t = 0u; t < 3u; t = t + 1u) {
        var m = -rho * bare[t];
        for (var q = 0u; q < Q; q = q + 1u) {
            m = m + streamed(q, cell, xyz) * f32(c_of(q)[t]);
        }
        for (var q = 1u; q < Q; q = q + 1u) {
            let c = c_of(q);
            if closure && t != axis && f32(c[axis]) * sign > 0.0 && c[t] != 0 {
                store_streamed(q, cell, xyz, streamed(q, cell, xyz) - 0.5 * m * f32(c[t]));
            }
        }
    }
}

fn momentum_exchange(l: u32, lu: vec4<u32>, q: u32, kind: u32, q2: u32) -> vec3<f32> {
    var v = vec3<f32>(0.0);
    if kind == KIND_BOUNCE {
        let out = post(opp(q), lu.w);
        let back = out + link_f[l];
        v = -vec3<f32>(c_of(q)) * (out + back);
    } else if kind == KIND_REFLECT {
        let out = post(q2, lu.y);
        v = vec3<f32>(c_of(q2)) * out - vec3<f32>(c_of(q)) * out;
    }
    return v;
}

// The momentum each link transferred to its group in the last step, from the post-collision
// populations in f (Ladd 1994; Mei et al. 2002): c_out f* − c_in f'.
@compute @workgroup_size(256)
fn link_forces(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let l = cell_of(lid, wg);
    if l >= params.links { return; }
    let n = params.n;
    let lu = link_u[l];
    let q = lu.x & 0xffu;
    let kind = (lu.x >> 8u) & 0xffu;
    let q2 = (lu.x >> 16u) & 0xffu;
    let cell = lu.w;
    // The deviations only: the weights' share of the exchange (2 w_q c_q per bounced link, which
    // cancels between a wall's links in exact arithmetic but is 10⁴ times the shear they carry)
    // is added on the host in double precision, so the sum keeps the deviations' digits.
    var v = contrib[l].xyz;
    if !AA_MODE { v = momentum_exchange(l, lu, q, kind, q2); }
    var heat = 0.0;
    // The heat the link gave its wall: what left the fluid minus what came back.
    if has_heat() && q < QT && (kind == KIND_BOUNCE || kind == KIND_REFLECT) && (lu.x & RESET) == 0u {
        heat = t[opp(q) * n + cell] - thermal_rule(lu, q, kind, cell);
    }
    contrib[l] = vec4<f32>(v, heat);
}

// ---- The scalar lattice ----

// Two-relaxation-time collision of the scalar lattice, in place in t, with the velocity of the
// momentum collision: each axis pair's symmetric part relaxes with 1/τ_T⁺, its antisymmetric
// part with 1/τ_T⁻, the rest population with 1/τ_T⁺ (fieldcraft_lbm::flow).
@compute @workgroup_size(256)
fn collide_thermal(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let cell = cell_of(lid, wg);
    if cell >= params.n { return; }
    if !is_fluid(cell) { return; }
    let n = params.n;
    var tq: array<f32, 7>;
    var theta = 0.0;
    for (var q = 0u; q < QT; q = q + 1u) {
        let v = t[q * n + cell];
        tq[q] = v;
        theta = theta + v;
    }
    let u = moments[cell].xyz;
    var wm = params.therm.x;
    var wp = params.therm.w;
    // With Smagorinsky's model the turbulent diffusivity ν_t / Pr_t raises τ_T by 4 ν_t / Pr_t
    // (fieldcraft_lbm::flow::collide_thermal): a single rate stays single, two keep their Λ.
    let nu_t = moments[cell].w;
    if turbulent() && nu_t > 0.0 && params.smag.y > 0.0 {
        let tau = 1.0 / wm + 4.0 * nu_t * params.smag.y;
        if abs(wp - wm) <= 1.1920929e-7 * abs(wm) {
            wp = 1.0 / tau;
        } else {
            let magic = (1.0 / wp - 0.5) * (1.0 / wm - 0.5);
            wp = 1.0 / (0.5 + magic / (tau - 0.5));
        }
        wm = 1.0 / tau;
    }
    t[cell] = tq[0] - wp * (tq[0] - eq_t(0u, theta, u));
    // the pairs (1, 2), (3, 4), (5, 6)
    for (var q = 1u; q < QT; q = q + 2u) {
        let a = tq[q];
        let b = tq[q + 1u];
        let ea = eq_t(q, theta, u);
        let eb = eq_t(q + 1u, theta, u);
        let even = wp * 0.5 * (a + b - ea - eb);
        let odd = wm * 0.5 * (a - b - ea + eb);
        t[q * n + cell] = a - even - odd;
        t[(q + 1u) * n + cell] = b - even + odd;
    }
}

// Pull streaming of t into t2, first pass (as `stream`).
@compute @workgroup_size(256)
fn stream_thermal(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let cell = cell_of(lid, wg);
    if cell >= params.n { return; }
    if !is_fluid(cell) { return; }
    let n = params.n;
    let nx = i32(params.nx);
    let ny = i32(params.ny);
    let nz = i32(params.nz);
    let i = i32(cell % params.nx);
    let jj = i32((cell / params.nx) % params.ny);
    let k = i32(cell / (params.nx * params.ny));
    for (var q = 0u; q < QT; q = q + 1u) {
        let c = c_of(q);
        let si = (i - c.x + nx) % nx;
        let sj = (jj - c.y + ny) % ny;
        let sk = (k - c.z + nz) % nz;
        let src = u32((sk * ny + sj) * nx + si);
        t2[q * n + cell] = t[q * n + src];
    }
}

// The scalar population a link's rule gives, from t: a pull as the flow's; at a wall or a slip
// face the group's thermal condition on the bounced population; an open placeholder.
fn thermal_rule(lu: vec4<u32>, q: u32, kind: u32, cell: u32) -> f32 {
    let n = params.n;
    if kind == KIND_PULL { return t[q * n + lu.y]; }
    let back = t[opp(q) * n + cell];
    if kind == KIND_OPEN { return back; }
    let cond = group_t[lu.z];
    let k = u32(cond.x);
    if k == 1u { return -back + 2.0 * wt_of(q) * cond.y; }
    if k == 2u { return back + cond.y; }
    return back;
}

// Second pass, one thread per link of the first seven directions.
@compute @workgroup_size(256)
fn stream_links_thermal(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let l = cell_of(lid, wg);
    if l >= params.links { return; }
    let lu = link_u[l];
    let q = lu.x & 0xffu;
    if q >= QT { return; }
    let kind = (lu.x >> 8u) & 0xffu;
    let cell = lu.w;
    t2[q * params.n + cell] = thermal_rule(lu, q, kind, cell);
}

// Inlets hold their temperature (the equilibrium at the cell's velocity); outlets take the
// streamed populations of the neighbour one cell inward.
@compute @workgroup_size(256)
fn open_thermal(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let o = cell_of(lid, wg);
    if o >= params.open { return; }
    let n = params.n;
    let ou = open_u[o];
    let cell = ou.x;
    let held = open_t[o];
    if held.x > 0.5 {
        let u = moments[cell].xyz;
        for (var q = 0u; q < QT; q = q + 1u) {
            t2[q * n + cell] = eq_t(q, held.y, u);
        }
    } else {
        var step = 1;
        if ou.y == 1u { step = i32(params.nx); }
        if ou.y == 2u { step = i32(params.nx * params.ny); }
        let src = u32(i32(cell) + i32(open_f[o].w) * step);
        for (var q = 0u; q < QT; q = q + 1u) {
            t2[q * n + cell] = t2[q * n + src];
        }
    }
}

// Flags a temperature that stopped being a number.
@compute @workgroup_size(256)
fn thermal_check(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let cell = cell_of(lid, wg);
    if cell >= params.n { return; }
    if !is_fluid(cell) { return; }
    let theta = theta_of(cell);
    if theta != theta || abs(theta) > 1e30 { atomicStore(&stats[2], 1u); }
}

// One workgroup per group: the sum of its links' contributions.
@compute @workgroup_size(256)
fn group_forces(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let group = wg.y * params.groups_x + wg.x;
    if group >= params.groups { return; }
    let lo = group_start[group];
    let hi = group_start[group + 1u];
    var s = vec4<f32>(0.0);
    for (var i = lo + lid; i < hi; i = i + WG) {
        s = s + contrib[by_group[i]];
    }
    red[lid] = s;
    workgroupBarrier();
    for (var k = 128u; k > 0u; k = k >> 1u) {
        if lid < k { red[lid] = red[lid] + red[lid + k]; }
        workgroupBarrier();
    }
    if lid == 0u { forces[group] = red[0]; }
}

// The density and velocity (with the half force) of every cell from f; zero in solid cells.
@compute @workgroup_size(256)
fn cell_moments(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let cell = cell_of(lid, wg);
    if cell >= params.n { return; }
    if !is_fluid(cell) {
        moments[cell] = vec4<f32>(0.0);
        return;
    }
    let n = params.n;
    var drho = 0.0;
    var j = vec3<f32>(0.0);
    let xyz = aa_coords(cell);
    for (var q = 0u; q < Q; q = q + 1u) {
        let v = population(q, cell, xyz);
        drho = drho + v;
        j = j + vec3<f32>(c_of(q)) * v;
    }
    let rho = 1.0 + drho;
    let u = (j + 0.5 * rho * acc_of(cell)) / rho;
    moments[cell] = vec4<f32>(u, rho);
}

// Per workgroup, the sum of its cells' densities (the mass, since every cell has unit volume).
@compute @workgroup_size(256)
fn mass(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let group = wg.y * params.groups_x + wg.x;
    let cell = group * WG + lid;
    var v = 0.0;
    if cell < params.n && is_fluid(cell) {
        v = moments[cell].w;
    }
    red[lid] = vec4<f32>(v, 0.0, 0.0, 0.0);
    workgroupBarrier();
    for (var k = 128u; k > 0u; k = k >> 1u) {
        if lid < k { red[lid] = red[lid] + red[lid + k]; }
        workgroupBarrier();
    }
    if lid == 0u { partials[group] = red[0].x; }
}

// ---- Conducting solids ----

const FACE_NONE: u32 = 0u;
const FACE_SOLID: u32 = 1u;
const FACE_HELD: u32 = 2u;
const FACE_FLUX: u32 = 3u;

fn s_m() -> u32 { return params.solid.x; }
fn s_k() -> u32 { return params.solid.y; }
fn s_share() -> f32 { return bitcast<f32>(params.solid.w); }

// solid_u: the cells' grid indices [0, m), each cell's first link [m, 2m + 1), the links ordered
// by cell, four each (fluid cell, q, conducting cell, link index), then per cell its six faces
// (+x, −x, +y, −y, +z, −z), two each: the kind and the other cell (a wall's: its box face).
fn s_first(c: u32) -> u32 { return solid_u[s_m() + c]; }
fn s_link(i: u32) -> vec4<u32> {
    let o = 2u * s_m() + 1u + 4u * i;
    return vec4<u32>(solid_u[o], solid_u[o + 1u], solid_u[o + 2u], solid_u[o + 3u]);
}
fn s_face(c: u32, k: u32) -> vec2<u32> {
    let o = 2u * s_m() + 1u + 4u * s_k() + 12u * c + 2u * k;
    return vec2<u32>(solid_u[o], solid_u[o + 1u]);
}
// solid_f: capacity [0, m), source per step [m, 2m), per cell its faces' (g, value) [2m, 14m),
// the links' weights k_s / (k_s + k_f) [14m, 14m + k).
fn s_capacity(c: u32) -> f32 { return solid_f[c]; }
fn s_source(c: u32) -> f32 { return solid_f[s_m() + c]; }
fn s_face_value(c: u32, k: u32) -> vec2<f32> {
    let o = 2u * s_m() + 12u * c + 2u * k;
    return vec2<f32>(solid_f[o], solid_f[o + 1u]);
}
fn s_weight(i: u32) -> f32 { return solid_f[14u * s_m() + i]; }
// solid_s: θ in two halves [0, 2m), the links' sums a [2m, 3m) and b [3m, 4m), the mean θ over
// the step [4m, 5m), the heat into box face f per cell [5m + 6c + f], then per link its heat
// [11m, 11m + k), its population leaving the fluid and its a [11m + k, 11m + 3k).
fn s_theta(half: u32, c: u32) -> u32 { return half * s_m() + c; }
fn s_a(c: u32) -> u32 { return 2u * s_m() + c; }
fn s_b(c: u32) -> u32 { return 3u * s_m() + c; }
fn s_mean(c: u32) -> u32 { return 4u * s_m() + c; }
fn s_wall(c: u32, f: u32) -> u32 { return 5u * s_m() + 6u * c + f; }
fn s_heat(i: u32) -> u32 { return 11u * s_m() + i; }
fn s_out(i: u32) -> u32 { return 11u * s_m() + s_k() + 2u * i; }

// Per conducting cell: the sums a and b of its fluid links (a link's heat is
// `2 out − 2 w (1 − weight) θ_f − 2 w weight θ_s`, kept per link for `solid_scatter`), and the
// step's accumulators zeroed.
@compute @workgroup_size(256)
fn solid_ab(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let c = cell_of(lid, wg);
    if c >= s_m() { return; }
    var a = 0.0;
    var b = 0.0;
    for (var i = s_first(c); i < s_first(c + 1u); i = i + 1u) {
        let l = s_link(i);
        let out = t[opp(l.y) * params.n + l.x];
        let w2 = 2.0 * wt_of(l.y);
        let weight = s_weight(i);
        let ai = 2.0 * out - w2 * (1.0 - weight) * theta_of(l.x);
        solid_s[s_out(i)] = out;
        solid_s[s_out(i) + 1u] = ai;
        a = a + ai;
        b = b + w2 * weight;
    }
    solid_s[s_a(c)] = a;
    solid_s[s_b(c)] = b;
    solid_s[s_mean(c)] = 0.0;
    for (var f = 0u; f < 6u; f = f + 1u) { solid_s[s_wall(c, f)] = 0.0; }
}

// One sub-step from half `src` into half `dst`: C (θ' − θ) = share (a + source + Σ faces − b θ'),
// the faces explicit at θ; a wall's heat goes to its box face.
fn solid_sub(c: u32, src: u32, dst: u32) {
    let share = s_share();
    let theta = solid_s[s_theta(src, c)];
    let cap = s_capacity(c);
    var heat = solid_s[s_a(c)] + s_source(c);
    for (var k = 0u; k < 6u; k = k + 1u) {
        let face = s_face(c, k);
        let gv = s_face_value(c, k);
        var into = 0.0;
        if face.x == FACE_SOLID {
            heat = heat + gv.x * (solid_s[s_theta(src, face.y)] - theta);
        } else if face.x == FACE_HELD {
            into = gv.x * (gv.y - theta);
        } else if face.x == FACE_FLUX {
            into = gv.y;
        }
        if face.x == FACE_HELD || face.x == FACE_FLUX {
            solid_s[s_wall(c, face.y)] = solid_s[s_wall(c, face.y)] - into * share;
            heat = heat + into;
        }
    }
    let next = (theta + heat * share / cap) / (1.0 + solid_s[s_b(c)] * share / cap);
    solid_s[s_theta(dst, c)] = next;
    solid_s[s_mean(c)] = solid_s[s_mean(c)] + next * share;
    if next != next || abs(next) > 1e30 { atomicStore(&stats[2], 1u); }
}

@compute @workgroup_size(256)
fn solid_sub0(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let c = cell_of(lid, wg);
    if c >= s_m() { return; }
    solid_sub(c, 0u, 1u);
}

@compute @workgroup_size(256)
fn solid_sub1(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let c = cell_of(lid, wg);
    if c >= s_m() { return; }
    solid_sub(c, 1u, 0u);
}

// After an odd number of sub-steps: the temperatures back into the first half.
@compute @workgroup_size(256)
fn solid_copy(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let c = cell_of(lid, wg);
    if c >= s_m() { return; }
    solid_s[s_theta(0u, c)] = solid_s[s_theta(1u, c)];
}

// Per interface link: its population at the solid's mean temperature, and the heat it took.
@compute @workgroup_size(256)
fn solid_scatter(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let i = cell_of(lid, wg);
    if i >= s_k() { return; }
    let l = s_link(i);
    let heat = solid_s[s_out(i) + 1u] - 2.0 * wt_of(l.y) * s_weight(i) * solid_s[s_mean(l.z)];
    t2[l.y * params.n + l.x] = solid_s[s_out(i)] - heat;
    solid_s[s_heat(i)] = heat;
}

// The last step's link heats into their groups' exchange (after `link_forces`).
@compute @workgroup_size(256)
fn solid_contrib(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let i = cell_of(lid, wg);
    if i >= s_k() { return; }
    let l = s_link(i).w;
    contrib[l] = vec4<f32>(contrib[l].xyz, solid_s[s_heat(i)]);
}

// One workgroup per box face (groups 0..6): the heat the solids gave it, added to its exchange
// (after `group_forces`).
@compute @workgroup_size(256)
fn solid_wall_heat(@builtin(local_invocation_index) lid: u32, @builtin(workgroup_id) wg: vec3<u32>) {
    let face = wg.x;
    if face >= 6u { return; }
    var s = 0.0;
    for (var c = lid; c < s_m(); c = c + WG) { s = s + solid_s[s_wall(c, face)]; }
    red[lid] = vec4<f32>(s, 0.0, 0.0, 0.0);
    workgroupBarrier();
    for (var k = 128u; k > 0u; k = k >> 1u) {
        if lid < k { red[lid] = red[lid] + red[lid + k]; }
        workgroupBarrier();
    }
    if lid == 0u && face < params.groups { forces[face] = forces[face] + vec4<f32>(0.0, 0.0, 0.0, red[0].x); }
}
