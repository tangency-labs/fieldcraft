//! The free surface: liquid moving through a gas that carries nothing, by the volume-of-fluid
//! lattice Boltzmann scheme of Körner et al. 2005 and Thürey 2007.
//!
//! # Scheme
//! Every fluid cell (a cell that is not an obstacle) is **gas**, **interface** or **liquid**.
//! Liquid cells run the ordinary collision and streaming; their mass is their density. Gas cells
//! hold no populations and do nothing. Interface cells, the one-cell layer between the two, run
//! the collision too and carry a mass `m` apart from their density `ρ`; their fill level
//! `ε = m / ρ` says how much of the cell is liquid.
//! - **Mass exchange** (Körner et al. 2005, eq. 9): after streaming, an interface cell's mass
//!   changes by `Δm = Σ_i (f_ī(x + c_i) − f_i(x))` over its liquid neighbours, and by the same
//!   weighted with `(ε(x) + ε(x + c_i)) / 2` over its interface neighbours; nothing passes to a
//!   gas cell or a wall. Each exchange is antisymmetric between its two cells, so the total mass,
//!   the liquid cells' densities plus the interface cells' masses, is conserved to round-off.
//! - **Reconstruction** (Körner et al. 2005, eq. 14; Thürey 2007, §4.2): the populations an
//!   interface cell would receive from a gas cell at `x + c_i` do not exist and are set to
//!   `f_ī = f_i^eq(ρ_A, u) + f_ī^eq(ρ_A, u) − f_i(x)`: the equilibria at the gas density `ρ_A`
//!   (the gauge pressure 0, so `ρ_A` = 1) and the cell's velocity, less the population leaving
//!   toward the gas. That is a free-slip, constant-pressure condition at the surface. The
//!   populations arriving from the gas side of the surface normal, `n · c_i > 0` with `n` pointing
//!   from the liquid into the gas, are reconstructed too (Thürey 2007), so the balance at the
//!   surface does not depend on which cells happen to be gas.
//! - **Normal**: `n_a = ½ (ε(x − e_a) − ε(x + e_a))`, the central difference of the fill levels
//!   (a neighbour outside the box or inside an obstacle takes the cell's own fill).
//! - **Conversion** (Thürey 2007, §4.3): an interface cell whose mass exceeds `(1 + κ) ρ`, or
//!   that has no gas neighbour left, becomes liquid; one whose mass falls below `−κ ρ`, or that
//!   has no liquid neighbour left, becomes gas (`κ` is [`FILL_OFFSET`], which keeps a cell from
//!   converting back in the next step). A cell that fills turns its gas neighbours into interface
//!   cells (mass 0; populations at the equilibrium of the mean density and velocity of their
//!   liquid and interface neighbours); a cell that empties turns its liquid neighbours into
//!   interface cells (mass: their density) and drops its own populations. A cell next to one that
//!   fills does not empty in the same step. The mass a conversion leaves over (`m − ρ` of a filled
//!   cell, `m` of an emptied one) goes to the interface neighbours, weighted by `n · c_i` toward
//!   the gas for a filled cell and away from it for an emptied one, evenly when no weight is
//!   positive, and into a reservoir shared out over every interface cell when there is none.
//! - **Slip faces**: a population leaving a cell toward a slip face is reflected specularly and
//!   lands in the cell shifted by its tangential velocity, so that mirror cell is the partner of
//!   the exchange along that direction (and the source of the reconstruction when it is gas); a
//!   wall returns the cell's own population and exchanges nothing.
//! - **Inlets**: a velocity inlet's cells are liquid for good and inject liquid through the
//!   Zou–He condition as usual; a cell next to an inlet never empties, so an inlet cell always
//!   streams from populated cells. The mass that enters through the face in a step, the incoming
//!   populations the boundary sets less the outgoing ones that left, is summed in cell order
//!   ([`Surface::injected`]).
//!
//! Surface tension is not modelled (`ρ_A` carries no curvature term): the scheme serves flows
//! gravity and inertia drive, such as a collapsing column, a filling cavity or a wave, where
//! the Bond and Weber numbers are large.
//!
//! # Sources
//! - C. Körner, M. Thies, T. Hofmann, N. Thürey, U. Rüde, *Lattice Boltzmann model for free
//!   surface flow for modeling foaming*, J. Stat. Phys. 121(1–2):179–196, 2005.
//! - N. Thürey, *Physically based animation of free surface flows with the lattice Boltzmann
//!   method*, PhD thesis, Universität Erlangen-Nürnberg, 2007, ch. 4.
//! - N. Thürey, U. Rüde, *Stable free surface flows with the lattice Boltzmann method on
//!   adaptively coarsened grids*, Comput. Visual. Sci. 12(5):247–263, 2009 (the conversion rules
//!   and the excess mass).
//!
//! No code of any free-surface lattice Boltzmann implementation was read.
//!
//! # Determinism
//! The mass exchange and the reconstruction run by fixed blocks of cells in parallel, each cell
//! reading the fill levels and populations of the last step and writing only its own mass and
//! populations; the conversions run on one thread in cell order.

use fieldcraft_parallel as parallel;

use crate::flow::{Flow, cell_blocks, filled};
use crate::lattice::{C, Q, c, equilibrium, opposite};
use crate::voxel::Grid;
use crate::{Error, Result};

/// `κ`: the share of a cell's density its mass must overshoot (1 + κ) or undershoot (−κ) before
/// the cell converts (Thürey 2007).
pub const FILL_OFFSET: f64 = 1.0e-3; // unit-factors: a share of a cell, not a conversion

/// The gas density: the gauge pressure 0.
pub const GAS_DENSITY: f64 = 1.0;

/// What a cell holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Kind {
    /// Inside an obstacle (not a fluid cell).
    Obstacle = 0,
    /// Gas: no populations.
    Gas = 1,
    /// The layer between liquid and gas: populations and a mass of its own.
    Interface = 2,
    /// Liquid: populations, whose sum is its mass.
    Liquid = 3,
}

impl Kind {
    fn from_u64(v: u64) -> Option<Kind> {
        Some(match v {
            0 => Kind::Obstacle,
            1 => Kind::Gas,
            2 => Kind::Interface,
            3 => Kind::Liquid,
            _ => return None,
        })
    }
}

/// A conversion marker.
const NONE: u8 = 0;
const FILLS: u8 = 1;
const EMPTIES: u8 = 2;

/// The free surface's state.
#[derive(Debug, Clone)]
pub struct Surface {
    pub(crate) kind: Vec<Kind>,
    /// An interface cell's mass (0 in every other cell: a liquid cell's mass is its density).
    pub(crate) mass: Vec<f64>,
    /// Each cell's fill level after the last step: 1 liquid, 0 gas, `m / ρ` clamped to 0..1 in an
    /// interface cell.
    pub(crate) fill: Vec<f64>,
    /// Liquid or interface: the cells that collide and stream.
    pub(crate) active: Vec<bool>,
    /// Whether a cell is an inlet cell or next to one: it never empties.
    near_inlet: Vec<bool>,
    /// Per box face (in `BoxFace::ALL` order): whether it reflects specularly (a slip face that
    /// is not periodic).
    reflect: [bool; 6],
    /// A conversion marker per cell (scratch, cleared after use).
    marks: Vec<u8>,
    /// Mass left over by conversions with no interface neighbour, shared out over the interface
    /// cells when there are any.
    pub(crate) reservoir: f64,
    /// The mass the inlets injected so far (lattice units).
    pub(crate) injected: f64,
    /// The total mass at the start.
    pub(crate) initial_mass: f64,
    /// Conversions over the run: cells filled, cells emptied.
    pub(crate) filled: u64,
    pub(crate) emptied: u64,
}

/// How the free surface fared, for the run statistics (lattice units of mass: cell volumes at
/// density 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceReport {
    pub initial_mass: f64,
    pub mass: f64,
    pub injected: f64,
    /// `(mass − initial_mass − injected) / max(initial_mass, injected)`: what the scheme lost or
    /// made, which should be round-off.
    pub mass_error: f64,
    /// The share of the fluid cells' volume that is liquid: the mean fill level.
    pub fill: f64,
    pub liquid_cells: usize,
    pub interface_cells: usize,
    pub gas_cells: usize,
    pub filled: u64,
    pub emptied: u64,
}

/// The cell at `p + c_q` (across a periodic face its image), or `None` outside the box.
#[inline(always)]
pub(crate) fn neighbour(grid: &Grid, periodic: [bool; 3], p: [usize; 3], q: usize) -> Option<usize> {
    let cq = C.get(q)?;
    let mut s = [0usize; 3];
    for a in 0..3 {
        let n = *grid.cells.get(a)? as i64;
        let mut v = *p.get(a)? as i64 + i64::from(cq[a]);
        if v < 0 || v >= n {
            if !periodic.get(a).copied().unwrap_or(false) {
                return None;
            }
            v = v.rem_euclid(n);
        }
        s[a] = v as usize;
    }
    Some(grid.index(s))
}

/// The cell whose population arrives at `p` from direction `q` after streaming, and which the
/// population leaving `p` along `c_q` reaches: the neighbour at `p + c_q` (across a periodic face
/// its image); across one specularly reflecting face the mirror cell, `p` shifted by the
/// tangential part of `c_q`, when that is a fluid cell (as the streaming rules reflect); `None` at
/// a wall, an edge of two faces or an obstacle, where the population bounces back.
#[inline(always)]
fn partner(grid: &Grid, periodic: [bool; 3], reflect: [bool; 6], kind: &[Kind], p: [usize; 3], q: usize) -> Option<usize> {
    let cq = C.get(q)?;
    let mut s = [0i64; 3];
    let mut crossed: Option<usize> = None;
    for a in 0..3 {
        let n = *grid.cells.get(a)? as i64;
        let mut v = *p.get(a)? as i64 + i64::from(cq[a]);
        if v < 0 || v >= n {
            if periodic.get(a).copied().unwrap_or(false) {
                v = v.rem_euclid(n);
            } else {
                // the face index in `BoxFace::ALL` order: 2 a for the min face, 2 a + 1 for the max
                let face = 2 * a + usize::from(v >= n);
                if crossed.is_some() || !reflect.get(face).copied().unwrap_or(false) {
                    return None;
                }
                crossed = Some(a);
                v = *p.get(a)? as i64;
            }
        }
        s[a] = v;
    }
    let cell = grid.index([s[0] as usize, s[1] as usize, s[2] as usize]);
    match crossed {
        Some(_) if kind.get(cell) == Some(&Kind::Obstacle) => None,
        _ => Some(cell),
    }
}

/// The surface normal at `p`, from the liquid into the gas: the central difference of the fill
/// levels, a missing neighbour taking the cell's own.
#[inline(always)]
fn normal(grid: &Grid, periodic: [bool; 3], fill: &[f64], p: [usize; 3], own: f64) -> [f64; 3] {
    let at = |q: usize| neighbour(grid, periodic, p, q).and_then(|c| fill.get(c)).copied().unwrap_or(own);
    // directions 1, 3, 5 are +x, +y, +z; 2, 4, 6 their opposites
    [0.5 * (at(2) - at(1)), 0.5 * (at(4) - at(3)), 0.5 * (at(6) - at(5))]
}

/// The density and momentum of the populations of `cell` in `f` (direction stride `s`).
#[inline(always)]
fn moments(f: &[f64], s: usize, cell: usize) -> (f64, [f64; 3]) {
    let (mut rho, mut j) = (0.0, [0.0; 3]);
    for q in 0..Q {
        let v = f.get(q * s + cell).copied().unwrap_or(0.0);
        let cq = c(q);
        rho += v;
        j[0] += v * cq[0];
        j[1] += v * cq[1];
        j[2] += v * cq[2];
    }
    (rho, j)
}

impl Surface {
    /// The surface of a grid whose `fluid` cells are not obstacles: liquid where `liquid(cell)`
    /// or in an `inlet` cell, interface in every gas cell next to a liquid cell, gas elsewhere.
    pub(crate) fn new(
        grid: &Grid,
        periodic: [bool; 3],
        reflect: [bool; 6],
        fluid: &[bool],
        inlets: &[u32],
        liquid: impl Fn(usize) -> bool,
    ) -> Result<Surface> {
        let n = grid.len();
        let mut kind = filled(n, Kind::Gas, "cell kinds")?;
        for (cell, k) in kind.iter_mut().enumerate() {
            if !fluid.get(cell).copied().unwrap_or(false) {
                *k = Kind::Obstacle;
            } else if liquid(cell) {
                *k = Kind::Liquid;
            }
        }
        let mut near_inlet = filled(n, false, "inlet marks")?;
        for &cell in inlets {
            let cell = cell as usize;
            if let Some(k) = kind.get_mut(cell).filter(|k| **k != Kind::Obstacle) {
                *k = Kind::Liquid;
            }
            if let Some(m) = near_inlet.get_mut(cell) {
                *m = true;
            }
            let p = grid.coords(cell);
            for q in 1..Q {
                if let Some(m) = neighbour(grid, periodic, p, q).and_then(|c| near_inlet.get_mut(c)) {
                    *m = true;
                }
            }
        }
        // the interface: every gas cell with a liquid neighbour
        let mut interface = Vec::new();
        for cell in 0..n {
            if kind.get(cell) != Some(&Kind::Liquid) {
                continue;
            }
            let p = grid.coords(cell);
            for q in 1..Q {
                if let Some(nb) = neighbour(grid, periodic, p, q)
                    && kind.get(nb) == Some(&Kind::Gas)
                {
                    interface.push(nb);
                }
            }
        }
        for nb in interface {
            if let Some(k) = kind.get_mut(nb) {
                *k = Kind::Interface;
            }
        }
        let mut fill = filled(n, 0.0, "fill levels")?;
        let mut active = filled(n, false, "active cells")?;
        for cell in 0..n {
            match kind.get(cell) {
                Some(Kind::Liquid) => {
                    fill[cell] = 1.0;
                    active[cell] = true;
                }
                Some(Kind::Interface) => active[cell] = true,
                _ => {}
            }
        }
        Ok(Surface {
            kind,
            mass: filled(n, 0.0, "cell masses")?,
            fill,
            active,
            near_inlet,
            reflect,
            marks: filled(n, NONE, "conversion marks")?,
            reservoir: 0.0,
            injected: 0.0,
            initial_mass: 0.0,
            filled: 0,
            emptied: 0,
        })
    }

    pub fn kind(&self, cell: usize) -> Option<Kind> {
        self.kind.get(cell).copied()
    }

    /// The fill level of a cell: 1 liquid, 0 gas (and obstacles), between in an interface cell.
    pub fn fill(&self, cell: usize) -> f64 {
        self.fill.get(cell).copied().unwrap_or(0.0)
    }

    /// Every cell's fill level.
    pub fn fills(&self) -> &[f64] {
        &self.fill
    }

    /// The mass the inlets injected so far (lattice units).
    pub fn injected(&self) -> f64 {
        self.injected
    }

    /// The total mass: the liquid cells' densities (from their populations `f`, direction stride
    /// `s`), the interface cells' masses and the reservoir, summed in cell order.
    pub(crate) fn total_mass(&self, f: &[f64], s: usize) -> f64 {
        let mut total = self.reservoir;
        for (cell, k) in self.kind.iter().enumerate() {
            match k {
                Kind::Liquid => total += moments(f, s, cell).0,
                Kind::Interface => total += self.mass.get(cell).copied().unwrap_or(0.0),
                Kind::Gas | Kind::Obstacle => {}
            }
        }
        total
    }

    /// The report of the state now, the populations `f` read as in [`Surface::total_mass`].
    pub(crate) fn report(&self, f: &[f64], s: usize) -> SurfaceReport {
        let mass = self.total_mass(f, s);
        let (mut liquid, mut interface, mut gas) = (0usize, 0usize, 0usize);
        let mut fill = 0.0;
        for (k, e) in self.kind.iter().zip(&self.fill) {
            match k {
                Kind::Liquid => liquid += 1,
                Kind::Interface => interface += 1,
                Kind::Gas => gas += 1,
                Kind::Obstacle => continue,
            }
            fill += e;
        }
        let cells = (liquid + interface + gas).max(1) as f64;
        let scale = self.initial_mass.abs().max(self.injected.abs());
        SurfaceReport {
            initial_mass: self.initial_mass,
            mass,
            injected: self.injected,
            mass_error: if scale > 0.0 { (mass - self.initial_mass - self.injected) / scale } else { 0.0 },
            fill: fill / cells,
            liquid_cells: liquid,
            interface_cells: interface,
            gas_cells: gas,
            filled: self.filled,
            emptied: self.emptied,
        }
    }

    /// The liquid's extent now: the lowest and highest cell index along each axis of the cells at
    /// least half full (`None` when no cell is), and the liquid's volume in cells (the fill
    /// levels summed).
    pub(crate) fn extent(&self, grid: &Grid) -> (Option<([usize; 3], [usize; 3])>, f64) {
        let (mut lo, mut hi) = ([usize::MAX; 3], [0usize; 3]);
        let mut volume = 0.0;
        let mut any = false;
        for (cell, e) in self.fill.iter().enumerate() {
            volume += e;
            if *e >= 0.5 {
                any = true;
                let p = grid.coords(cell);
                for a in 0..3 {
                    lo[a] = lo[a].min(p[a]);
                    hi[a] = hi[a].max(p[a]);
                }
            }
        }
        (any.then_some((lo, hi)), volume)
    }

    pub(crate) fn heap_bytes(&self) -> u64 {
        ((self.kind.capacity() + self.active.capacity() + self.marks.capacity() + self.near_inlet.capacity()) as u64)
            + ((self.mass.capacity() + self.fill.capacity()) as u64) * 8
    }

    pub(crate) fn checkpoint(&mut self, a: &mut dyn crate::checkpoint::Archive, n: usize) -> Result<()> {
        use crate::checkpoint::{doubles, f64, fixed, u64};
        fixed(a, n as u64)?;
        for k in self.kind.iter_mut() {
            let mut v = *k as u64;
            u64(a, &mut v)?;
            *k = Kind::from_u64(v).ok_or_else(|| Error::Invalid("checkpoint cell kind is invalid".into()))?;
        }
        for flag in self.near_inlet.iter().chain(self.reflect.iter()) {
            fixed(a, u64::from(*flag))?;
        }
        doubles(a, &mut self.mass, n)?;
        doubles(a, &mut self.fill, n)?;
        if self.fill.iter().any(|e| !(0.0..=1.0).contains(e)) {
            return Err(Error::Invalid("checkpoint fill level is outside 0 to 1".into()));
        }
        f64(a, &mut self.reservoir)?;
        f64(a, &mut self.injected)?;
        f64(a, &mut self.initial_mass)?;
        u64(a, &mut self.filled)?;
        u64(a, &mut self.emptied)?;
        if a.loading() {
            for (active, k) in self.active.iter_mut().zip(&self.kind) {
                *active = matches!(k, Kind::Liquid | Kind::Interface);
            }
            if self.kind.iter().zip(&self.mass).any(|(k, m)| !matches!(k, Kind::Interface) && *m != 0.0) {
                return Err(Error::Invalid("checkpoint has mass outside the interface cells".into()));
            }
        }
        Ok(())
    }
}

impl Flow {
    /// After streaming (`f` the populations before it, `g` after): every interface cell's mass
    /// exchange with its liquid and interface neighbours, the reconstruction of its populations
    /// from the gas, and the gas cells' populations set to zero; by blocks in parallel, each cell
    /// writing its own mass and populations.
    pub(crate) fn surface_stream(&mut self) {
        let Some(surface) = &mut self.surface else { return };
        let (n, s, block) = (self.n, self.stride, self.block);
        let (grid, periodic) = (self.grid, self.periodic);
        let f = &self.f;
        let (kind, fill, reflect) = (&surface.kind, &surface.fill, surface.reflect);
        let mut blocks = cell_blocks(&mut self.g, s, n, block);
        let mut masses: Vec<&mut [f64]> = surface.mass.chunks_mut(block).collect();
        let mut work: Vec<(Vec<&mut [f64]>, &mut [f64])> = Vec::with_capacity(blocks.len());
        for (g, m) in blocks.drain(..).zip(masses.drain(..)) {
            work.push((g, m));
        }
        parallel::for_each_mut(&mut work, |b, (g, mass)| {
            let len = g.first().map_or(0, |d| d.len());
            let mut fq = [0.0; Q];
            for i in 0..len {
                let cell = b * block + i;
                match kind.get(cell) {
                    Some(Kind::Gas) => {
                        for d in g.iter_mut() {
                            if let Some(x) = d.get_mut(i) {
                                *x = 0.0;
                            }
                        }
                        continue;
                    }
                    Some(Kind::Interface) => {}
                    _ => continue,
                }
                for (q, x) in fq.iter_mut().enumerate() {
                    *x = f.get(q * s + cell).copied().unwrap_or(0.0);
                }
                let (rho, j) = {
                    let (mut rho, mut j) = (0.0, [0.0; 3]);
                    for (q, v) in fq.iter().enumerate() {
                        let cq = c(q);
                        rho += v;
                        j[0] += v * cq[0];
                        j[1] += v * cq[1];
                        j[2] += v * cq[2];
                    }
                    (rho, j)
                };
                let u = if rho > 0.0 { j.map(|v| v / rho) } else { [0.0; 3] };
                let own = fill.get(cell).copied().unwrap_or(0.0);
                let p = grid.coords(cell);
                // The mass exchange, from the streamed populations: `g[ī]` came from the
                // partner along `+c_q`.
                let mut dm = 0.0;
                for q in 1..Q {
                    let Some(nb) = partner(&grid, periodic, reflect, kind, p, q) else { continue };
                    let weight = match kind.get(nb) {
                        Some(Kind::Liquid) => 1.0,
                        Some(Kind::Interface) => 0.5 * (own + fill.get(nb).copied().unwrap_or(0.0)),
                        _ => continue,
                    };
                    let incoming = g.get(opposite(q)).and_then(|d| d.get(i)).copied().unwrap_or(0.0);
                    dm += weight * (incoming - fq[q]);
                }
                if let Some(m) = mass.get_mut(i) {
                    *m += dm;
                }
                // The reconstruction: from gas neighbours and from the gas side of the normal.
                let nrm = normal(&grid, periodic, fill, p, own);
                for q in 1..Q {
                    let cq = c(q);
                    let from_gas = partner(&grid, periodic, reflect, kind, p, q).is_some_and(|nb| kind.get(nb) == Some(&Kind::Gas));
                    let toward_gas = nrm[0] * cq[0] + nrm[1] * cq[1] + nrm[2] * cq[2] > 0.0;
                    if !(from_gas || toward_gas) {
                        continue;
                    }
                    let o = opposite(q);
                    let v = equilibrium(q, GAS_DENSITY, u) + equilibrium(o, GAS_DENSITY, u) - fq[q];
                    if let Some(x) = g.get_mut(o).and_then(|d| d.get_mut(i)) {
                        *x = v;
                    }
                }
            }
        });
    }

    /// After the step's populations are in place: the cells that fill or empty, their
    /// neighbours' new kinds, the excess mass shared out, and every cell's fill level; on one
    /// thread in cell order.
    pub(crate) fn surface_convert(&mut self) {
        let Some(surface) = &mut self.surface else { return };
        let (n, s) = (self.n, self.stride);
        let (grid, periodic) = (self.grid, self.periodic);
        let f = &mut self.f;
        let density = |f: &[f64], cell: usize| moments(f, s, cell).0;
        let mut fills: Vec<u32> = Vec::new();
        let mut empties: Vec<u32> = Vec::new();
        for cell in 0..n {
            if surface.kind.get(cell) != Some(&Kind::Interface) {
                continue;
            }
            let rho = density(f, cell);
            let m = surface.mass.get(cell).copied().unwrap_or(0.0);
            let p = grid.coords(cell);
            let (mut gas, mut liquid) = (false, false);
            for q in 1..Q {
                match neighbour(&grid, periodic, p, q).and_then(|nb| surface.kind.get(nb)) {
                    Some(Kind::Gas) => gas = true,
                    Some(Kind::Liquid) => liquid = true,
                    _ => {}
                }
            }
            if m > (1.0 + FILL_OFFSET) * rho || !gas {
                fills.push(cell as u32);
            } else if (m < -FILL_OFFSET * rho || !liquid) && !surface.near_inlet.get(cell).copied().unwrap_or(false) {
                empties.push(cell as u32);
            }
        }
        if fills.is_empty() && empties.is_empty() && surface.reservoir == 0.0 {
            // The fill levels of the interface cells, whose masses moved.
            for cell in 0..n {
                if surface.kind.get(cell) == Some(&Kind::Interface) {
                    let rho = density(f, cell);
                    let m = surface.mass.get(cell).copied().unwrap_or(0.0);
                    if let Some(e) = surface.fill.get_mut(cell) {
                        *e = if rho > 0.0 { (m / rho).clamp(0.0, 1.0) } else { 0.0 };
                    }
                }
            }
            return;
        }
        for &cell in &fills {
            if let Some(m) = surface.marks.get_mut(cell as usize) {
                *m = FILLS;
            }
        }
        for &cell in &empties {
            if let Some(m) = surface.marks.get_mut(cell as usize) {
                *m = EMPTIES;
            }
        }
        // A cell next to one that fills does not empty in the same step.
        for &cell in &fills {
            let p = grid.coords(cell as usize);
            for q in 1..Q {
                if let Some(m) = neighbour(&grid, periodic, p, q).and_then(|nb| surface.marks.get_mut(nb))
                    && *m == EMPTIES
                {
                    *m = NONE;
                }
            }
        }
        // (cell, excess, whether it filled)
        let mut leftover: Vec<(u32, f64, bool)> = Vec::with_capacity(fills.len() + empties.len());
        for &cell in &fills {
            let cell = cell as usize;
            let rho = density(f, cell);
            let m = surface.mass.get(cell).copied().unwrap_or(0.0);
            leftover.push((cell as u32, m - rho, true));
            surface.kind[cell] = Kind::Liquid;
            surface.active[cell] = true;
            surface.mass[cell] = 0.0;
            surface.fill[cell] = 1.0;
            surface.filled += 1;
            let p = grid.coords(cell);
            for q in 1..Q {
                let Some(nb) = neighbour(&grid, periodic, p, q) else { continue };
                if surface.kind.get(nb) != Some(&Kind::Gas) {
                    continue;
                }
                // A new interface cell: the equilibrium of its populated neighbours' mean state.
                let (mut rho_sum, mut u_sum, mut count) = (0.0, [0.0; 3], 0.0);
                let pn = grid.coords(nb);
                for r in 1..Q {
                    if let Some(c2) = neighbour(&grid, periodic, pn, r)
                        && matches!(surface.kind.get(c2), Some(Kind::Liquid | Kind::Interface))
                    {
                        let (rho2, j2) = moments(f, s, c2);
                        if rho2 > 0.0 {
                            rho_sum += rho2;
                            for a in 0..3 {
                                u_sum[a] += j2[a] / rho2;
                            }
                            count += 1.0;
                        }
                    }
                }
                let (rho_bar, u_bar) = if count > 0.0 { (rho_sum / count, u_sum.map(|v| v / count)) } else { (GAS_DENSITY, [0.0; 3]) };
                for q2 in 0..Q {
                    if let Some(x) = f.get_mut(q2 * s + nb) {
                        *x = equilibrium(q2, rho_bar, u_bar);
                    }
                }
                surface.kind[nb] = Kind::Interface;
                surface.active[nb] = true;
                surface.mass[nb] = 0.0;
                surface.fill[nb] = 0.0;
            }
        }
        for &cell in &empties {
            let cell = cell as usize;
            if surface.marks.get(cell) != Some(&EMPTIES) {
                continue;
            }
            let m = surface.mass.get(cell).copied().unwrap_or(0.0);
            leftover.push((cell as u32, m, false));
            surface.kind[cell] = Kind::Gas;
            surface.active[cell] = false;
            surface.mass[cell] = 0.0;
            surface.fill[cell] = 0.0;
            surface.emptied += 1;
            for q in 0..Q {
                if let Some(x) = f.get_mut(q * s + cell) {
                    *x = 0.0;
                }
            }
            let p = grid.coords(cell);
            for q in 1..Q {
                let Some(nb) = neighbour(&grid, periodic, p, q) else { continue };
                if surface.kind.get(nb) == Some(&Kind::Liquid) {
                    surface.kind[nb] = Kind::Interface;
                    surface.mass[nb] = density(f, nb);
                    surface.fill[nb] = 1.0;
                }
            }
        }
        for &cell in fills.iter().chain(&empties) {
            if let Some(m) = surface.marks.get_mut(cell as usize) {
                *m = NONE;
            }
        }
        // The excess mass, to the interface neighbours by the normal, else evenly, else the
        // reservoir.
        for (cell, excess, filled) in leftover {
            let cell = cell as usize;
            let p = grid.coords(cell);
            let own = surface.fill.get(cell).copied().unwrap_or(0.0);
            let nrm = normal(&grid, periodic, &surface.fill, p, own);
            let sign = if filled { 1.0 } else { -1.0 };
            let mut weights = [0.0; Q];
            let (mut total, mut count) = (0.0, 0usize);
            for q in 1..Q {
                let Some(nb) = neighbour(&grid, periodic, p, q) else { continue };
                if surface.kind.get(nb) != Some(&Kind::Interface) {
                    continue;
                }
                let cq = c(q);
                weights[q] = (sign * (nrm[0] * cq[0] + nrm[1] * cq[1] + nrm[2] * cq[2])).max(0.0);
                total += weights[q];
                count += 1;
            }
            if count == 0 {
                surface.reservoir += excess;
                continue;
            }
            for q in 1..Q {
                let Some(nb) = neighbour(&grid, periodic, p, q) else { continue };
                if surface.kind.get(nb) != Some(&Kind::Interface) {
                    continue;
                }
                let share = if total > 0.0 { weights[q] / total } else { 1.0 / count as f64 };
                if let Some(m) = surface.mass.get_mut(nb) {
                    *m += excess * share;
                }
            }
        }
        // The reservoir, shared out evenly over the interface cells when there are any.
        if surface.reservoir != 0.0 {
            let count = surface.kind.iter().filter(|k| **k == Kind::Interface).count();
            if count > 0 {
                let share = surface.reservoir / count as f64;
                for (k, m) in surface.kind.iter().zip(surface.mass.iter_mut()) {
                    if *k == Kind::Interface {
                        *m += share;
                    }
                }
                surface.reservoir = 0.0;
            }
        }
        for cell in 0..n {
            if surface.kind.get(cell) == Some(&Kind::Interface) {
                let rho = density(f, cell);
                let m = surface.mass.get(cell).copied().unwrap_or(0.0);
                if let Some(e) = surface.fill.get_mut(cell) {
                    *e = if rho > 0.0 { (m / rho).clamp(0.0, 1.0) } else { 0.0 };
                }
            }
        }
    }
}
