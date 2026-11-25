//! `results.openfoam`: the results of a case `deck.openfoam` wrote and OpenFOAM ran, read back
//! onto the step's lattice as its fluid results, for viewing with FieldCraft's slices, probes,
//! lines and measurements. The files are read in OpenFOAM's documented ASCII format
//! ([`fieldcraft_decks::openfoam::fields`]); what OpenFOAM computed is shown, never used to judge
//! FieldCraft's own results. A compressible case (rhoCentralFoam) also gives its temperature and
//! each cell's density and Mach number, from the gas it was written with; a pressure or temperature
//! that is not positive is refused, naming its cell.
//!
//! Every count a file gives is bounded before it is used: the cell centres by the memory the
//! fluid results may take, each field by the centres' count, the processor directories by the
//! ranks a case is written for, and the frames by the frame memory.

use std::path::{Path, PathBuf};

use fieldcraft_decks::openfoam::fields::{Count, FieldValues, read_field};
use fieldcraft_decks::openfoam::{MAX_PROCESSES, sample_reach};
use fieldcraft_io::roots::Access;
use fieldcraft_lbm::solve::{FluidResult, FluidStats, Frame};
use fieldcraft_lbm::voxel::Grid;
use fieldcraft_lbm::{Limits, PortProperties};
use fieldcraft_model::units::Quantity;
use fieldcraft_model::{BoxFace, FaceCondition};
use serde_json::{Value, json};

use super::{MARKER, Marker, MarkerGas, nearest};
use crate::cmd::opt;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};

pub const READ: CommandSpec = CommandSpec {
    id: "results.openfoam",
    label: "Read OpenFOAM Results",
    params: r#"{"path": "case_dir", "spacing"?: 0.4, "max_cells"?: 16000000, "every"?: 1, "from"?: 0, "to"?: 1} — reads the results of a case `deck.openfoam` wrote and OpenFOAM ran (ASCII; reconstructed, or else straight from its processor directories) as the fluid results of its step, for viewing: the written times from `from` to `to` (input units, inclusive), every `every`-th of them and always the last, each a frame sampled at the centres of a lattice over the step's domain box (its spacing, or `spacing` in the input units) from the nearest OpenFOAM cell; the obstacles' cells are solid. Velocity, and the gauge pressure (OpenFOAM's kinematic pressure times the density the case was written with). A compressible case (rhoCentralFoam): the gauge pressure (the absolute pressure less the ambient the case was written with), the temperature above the case's temperature, each cell's density p/(R T) from the gas the case was written with, and `max_mach`, the largest Mach number of the frames read. A pressure or temperature that is not positive is an error naming its cell. A lattice cell no OpenFOAM cell is near is an error, never a zero. Replaces the step's results; the model is unchanged"#,
    form: Some(crate::cmd::forms::RESULTS_OPENFOAM),
    kind: CommandKind::Session,
    requires: &[Requirement::FluidDomain, Requirement::Step],
    run: read,
};

/// The most bytes one field file may have: the memory the fluid results' frames may take
/// ([`Limits::MAX_FRAME_BYTES`]); its text is read whole.
const MAX_FIELD_BYTES: u64 = Limits::MAX_FRAME_BYTES;

/// The most OpenFOAM cells a case may have: as many centres as fit the frame memory.
const MAX_OPENFOAM_CELLS: usize = (Limits::MAX_FRAME_BYTES / CENTRE_BYTES) as usize;

/// The bytes an OpenFOAM cell's centre takes (3 × f64).
const CENTRE_BYTES: u64 = std::mem::size_of::<[f64; 3]>() as u64;

/// The bytes a lattice cell takes beside its frames: its solid flag, its nearest OpenFOAM cell
/// (`Option<u32>`) and its bucket offset (`u32`).
const LATTICE_BYTES_PER_CELL: u64 = (std::mem::size_of::<bool>() + std::mem::size_of::<Option<u32>>() + std::mem::size_of::<u32>()) as u64;

/// The bytes an OpenFOAM cell takes once its centre is bucketed: the centre and its bucket id (`u32`).
const BUCKETED_CELL_BYTES: u64 = CENTRE_BYTES + std::mem::size_of::<u32>() as u64;

/// The bytes an OpenFOAM cell takes while a field is joined: the part's values and their copy in
/// the joined vector, each a 3 × f64 (a scalar field takes less).
const JOINING_CELL_BYTES: u64 = 2 * CENTRE_BYTES;

/// What a read may take at its peak: the session's memory bound, or no bound when the machine's
/// memory is unknown (the fixed caps then hold, the defaults `solve.fluid` uses).
#[derive(Clone, Copy)]
struct Budget {
    bound: u64,
    /// The bytes held from the start: the lattice and every frame.
    resident: u64,
}

impl Budget {
    /// The bytes left under the bound once `held` more are taken; `None` without a bound.
    fn room(&self, held: u64) -> Option<u64> {
        (self.bound != 0).then(|| self.bound.saturating_sub(self.resident).saturating_sub(held))
    }

    /// Refuses `held` bytes beyond the resident ones when together they pass the bound.
    fn check(&self, held: u64, what: &str) -> Result<()> {
        let peak = self.resident.saturating_add(held);
        if self.bound != 0 && peak > self.bound {
            return Err(EngineError::invalid(
                "path",
                format!(
                    "{what} would take {} at its peak, above this machine's bound for a run ({}, half its memory): coarsen `spacing`, read fewer frames, or use a smaller case",
                    fieldcraft_lbm::memory::megabytes(peak),
                    fieldcraft_lbm::memory::megabytes(self.bound)
                ),
            ));
        }
        Ok(())
    }

    /// The most bytes a file read now may have, given `held` other bytes.
    fn file_cap(&self, held: u64) -> u64 {
        self.room(held).map_or(MAX_FIELD_BYTES, |r| r.min(MAX_FIELD_BYTES))
    }
}

/// The bytes a frame takes per lattice cell: a velocity (3 × f32) and a pressure (f32); a
/// compressible frame also a temperature and a density (f32 each).
pub(super) fn frame_bytes_per_cell(gas: bool) -> u64 {
    (if gas { 6 } else { 4 }) * std::mem::size_of::<f32>() as u64
}

/// What a compressible frame keeps of one fluid cell, in f64 (the frame stores f32).
#[derive(Clone, Copy)]
struct GasCell {
    pressure: f64,
    temperature: f64,
    density: f64,
    mach: f64,
}

impl GasCell {
    /// A solid cell holds no gas: every value is zero.
    const SOLID: GasCell = GasCell { pressure: 0.0, temperature: 0.0, density: 0.0, mach: 0.0 };
}

/// A fluid cell's values from its absolute pressure `p` (Pa), absolute temperature `t` (K) and
/// speed (m/s): the gauge pressure, the temperature above the case's, the density of the ideal
/// gas ρ = p/(R_s T) (doc.openfoam.com, "Perfect gas") and the Mach number |u|/c, with
/// c = √(γ R_s T). `p` and `t` are positive: [`refuse_unless_positive`] checks every cell first.
fn gas_cell(gas: &MarkerGas, gamma: f64, p: f64, t: f64, speed: f64) -> GasCell {
    GasCell {
        pressure: p - gas.ambient,
        temperature: t - gas.temperature,
        density: p / (gas.gas_constant * t),
        mach: speed / (gamma * gas.gas_constant * t).sqrt(),
    }
}

/// Refuses the first cell of `values`, a field joined over `parts` (of `sizes` cells each), whose
/// `what` is not finite and positive, naming the part's file (`rel` under it) and the cell's
/// index in that file; `need` says why it must be. A run that failed leaves such a cell, and the
/// density and sound speed are undefined there. Every cell of the case is checked, before any
/// is sampled.
fn refuse_unless_positive(
    units: &fieldcraft_model::DisplaySystem,
    (parts, sizes): (&[PathBuf], &[usize]),
    rel: &Path,
    values: &[f64],
    what: &str,
    q: Quantity,
    need: &str,
) -> Result<()> {
    let Some((global, v)) = values.iter().enumerate().find(|(_, v)| !(v.is_finite() && **v > 0.0)) else { return Ok(()) };
    // The part holding the cell, and its index there.
    let (mut part, mut local) = (parts.first(), global);
    for (p, &n) in parts.iter().zip(sizes) {
        if local < n {
            part = Some(p);
            break;
        }
        local -= n;
    }
    let file = part.map_or_else(|| rel.to_path_buf(), |p| p.join(rel));
    Err(EngineError::invalid(
        "path",
        format!("{}: cell {local} has a {what} of {}: {need}", file.display(), crate::units::shown(units, crate::units::si(), q, *v)),
    ))
}

fn field(path: &Path, count: Count, max_bytes: u64) -> Result<FieldValues> {
    let bytes = super::super::file::read_file(&path.to_string_lossy(), max_bytes)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::invalid("path", format!("{} is not text", path.display())))?;
    read_field(text, count).map_err(|e| EngineError::context(format!("{}", path.display()), e))
}

/// The written times of a case part, in order (`0`, the initial state, left out).
fn times(root: &Path) -> Result<Vec<(f64, PathBuf)>> {
    let mut out: Vec<(f64, PathBuf)> = fieldcraft_io::fs::list_dirs(root)?
        .into_iter()
        .filter_map(|d| {
            let t: f64 = d.file_name()?.to_str()?.parse().ok()?;
            (t.is_finite() && t > 0.0).then_some((t, d))
        })
        .collect();
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    if let Some(w) = out.windows(2).find(|w| w[0].0 == w[1].0) {
        return Err(EngineError::invalid("path", format!("two time directories for the same time: {} and {}", w[0].1.display(), w[1].1.display())));
    }
    Ok(out)
}

/// Where a case's results are: the case itself (reconstructed) or its processor directories in
/// rank order (decomposed), each with its own `0/C` and time directories.
fn parts(root: &Path) -> Result<(Vec<PathBuf>, bool)> {
    if !times(root)?.is_empty() {
        return Ok((vec![root.to_path_buf()], false));
    }
    let mut ranks: Vec<(usize, PathBuf)> = fieldcraft_io::fs::list_dirs(root)?
        .into_iter()
        .filter_map(|d| Some((d.file_name()?.to_str()?.strip_prefix("processor")?.parse().ok()?, d)))
        .collect();
    if ranks.is_empty() {
        return Err(EngineError::invalid("path", "the case has no written times yet: run it (sh Allrun)"));
    }
    if ranks.len() > MAX_PROCESSES as usize {
        return Err(EngineError::invalid("path", format!("{} processor directories (a case is written for at most {MAX_PROCESSES})", ranks.len())));
    }
    ranks.sort_by_key(|r| r.0);
    if ranks.iter().enumerate().any(|(i, r)| r.0 != i) {
        return Err(EngineError::invalid("path", "the processor directories are not numbered 0 to N − 1"));
    }
    Ok((ranks.into_iter().map(|r| r.1).collect(), true))
}

/// A field over every part, joined in part order; `sizes` are the parts' cell counts.
fn joined(parts: &[PathBuf], sizes: &[usize], rel: &Path, vector: bool, budget: Budget) -> Result<(Vec<f64>, Vec<[f64; 3]>)> {
    let cells = sizes.iter().map(|&n| n as u64).sum::<u64>();
    let held = cells.saturating_mul(BUCKETED_CELL_BYTES + JOINING_CELL_BYTES);
    budget.check(held, "joining the fields")?;
    // Sized once to the cells, as the budget counts them (a growing vector would hold up to
    // twice that while it doubles).
    let all = usize::try_from(cells).map_err(|_| EngineError::invalid("path", "more OpenFOAM cells than this machine addresses"))?;
    let (mut scalars, mut vectors) = if vector { (Vec::new(), Vec::with_capacity(all)) } else { (Vec::with_capacity(all), Vec::new()) };
    for (part, &n) in parts.iter().zip(sizes) {
        let path = part.join(rel);
        let f = field(&path, Count::Exactly(n), budget.file_cap(held))?;
        let missing = || EngineError::invalid("path", format!("{}: not a {} field", path.display(), if vector { "vector" } else { "scalar" }));
        for i in 0..n {
            if vector {
                vectors.push(f.vector(i).ok_or_else(missing)?);
            } else {
                scalars.push(f.scalar(i).ok_or_else(missing)?);
            }
        }
    }
    Ok((scalars, vectors))
}

/// The times read: those in `[from, to]`, every `every`-th of them and always the last. Time
/// directory names carry `digits` significant digits, so a bound matches a name within half a
/// unit of its last digit, 0.5·10^(1 − digits) relative.
fn select(written: &[(f64, PathBuf)], every: usize, from: f64, to: f64, digits: u32) -> Vec<usize> {
    let rounding = 0.5 * 10f64.powi(1 - digits as i32);
    let tol = |t: f64| rounding * t.abs();
    let inside: Vec<usize> = written.iter().enumerate().filter(|(_, w)| w.0 >= from - tol(from) && w.0 <= to + tol(to)).map(|(i, _)| i).collect();
    let mut out: Vec<usize> = inside.iter().copied().step_by(every.max(1)).collect();
    if let Some(&last) = inside.last()
        && out.last() != Some(&last)
    {
        out.push(last);
    }
    out
}

fn read(s: &mut Session, p: &Value) -> Result<Value> {
    let dir = super::super::file_param(s, p, "path", Access::Read)?;
    let root = PathBuf::from(&dir);
    let marker_path = root.join(MARKER);
    let marker_bytes = super::super::file::read_file(&marker_path.to_string_lossy(), Marker::MAX_BYTES)
        .map_err(|e| EngineError::context(format!("{} is not a case deck.openfoam wrote (no {MARKER})", root.display()), e))?;
    let marker: Marker = serde_json::from_slice(&marker_bytes).map_err(|e| EngineError::invalid("path", format!("{MARKER}: {e}")))?;
    marker.check()?;
    let (step, domain) = fieldcraft_lbm::fluid_step(&s.model, Some(&marker.step))?;
    let step_name = step.name.clone();
    let every: usize = opt(p, "every")?.unwrap_or(1);
    if every == 0 {
        return Err(EngineError::invalid("every", "must be 1 or more"));
    }
    let from = opt::<f64>(p, "from")?.unwrap_or(f64::NEG_INFINITY);
    let to = opt::<f64>(p, "to")?.unwrap_or(f64::INFINITY);
    if from.is_nan() || to.is_nan() || from > to {
        return Err(EngineError::invalid("from", "must not be after `to`"));
    }
    let (parts, decomposed) = parts(&root)?;
    let written = times(parts.first().map_or(&root, |p| p.as_path()))?;
    if written.is_empty() {
        return Err(EngineError::invalid("path", "the case has no written times yet: run it (sh Allrun)"));
    }
    let chosen = select(&written, every, from, to, marker.time_precision);
    if chosen.is_empty() {
        let shown = |t: Option<&(f64, PathBuf)>| crate::units::shown(s.units(), crate::units::si(), Quantity::Time, t.map_or(0.0, |w| w.0));
        return Err(EngineError::invalid(
            "from",
            format!("no written time between `from` and `to`: the case has {} to {}", shown(written.first()), shown(written.last())),
        ));
    }

    // The lattice over the domain box.
    let spacing = opt::<f64>(p, "spacing")?.unwrap_or(domain.spacing);
    if !(spacing.is_finite() && spacing > 0.0) {
        return Err(EngineError::invalid("spacing", "must be a positive length"));
    }
    let bound = s.limits.memory;
    let ceilings = Limits { max_memory: bound, ..Limits::default() };
    let max_cells = super::super::solve::raisable(p, "max_cells", Limits::MAX_CELLS, ceilings.max_cells_ceiling(), bound)?;
    // The model's own rule: a whole number of spacings along each axis, within the cells per side.
    let mut sized = domain.clone();
    sized.spacing = spacing;
    let cells = sized.cells().map_err(|_| {
        EngineError::invalid(
            "spacing",
            format!(
                "{} does not divide the domain's box into a whole number of spacings along each axis, of at most {} cells",
                crate::units::shown(s.units(), crate::units::si(), Quantity::Length, spacing),
                fieldcraft_model::fluids::MAX_CELLS_PER_SIDE
            ),
        )
    })?;
    let total = cells.iter().try_fold(1u64, |acc, &n| acc.checked_mul(n)).unwrap_or(u64::MAX);
    if total > max_cells {
        return Err(EngineError::invalid(
            "spacing",
            format!("gives {total} lattice cells, above max_cells {max_cells}: coarsen it or raise max_cells"),
        ));
    }
    // At most `max_cells` (and `MAX_CELLS_PER_SIDE` along an axis), which an index addresses.
    let cells = cells.map(|n| n as usize);
    let grid = Grid { origin: domain.min, spacing, cells };
    let frame_bytes = total.saturating_mul(frame_bytes_per_cell(marker.gas.is_some())).saturating_mul(chosen.len() as u64);
    if frame_bytes > ceilings.max_frame_ceiling() {
        return Err(EngineError::invalid(
            "spacing",
            format!(
                "{} frames of {total} cells take more than {} bytes: coarsen it, or read fewer frames (every, from, to)",
                chosen.len(),
                ceilings.max_frame_ceiling()
            ),
        ));
    }
    let budget = Budget { bound, resident: frame_bytes.saturating_add(total.saturating_mul(LATTICE_BYTES_PER_CELL)) };
    budget.check(0, "the lattice and its frames")?;
    let mut solid = vec![false; grid.len()];
    for o in super::super::solve::obstacles(s, domain)? {
        fieldcraft_lbm::voxel::voxelise(&o.surface, &grid, &mut solid)?;
    }

    // The OpenFOAM cells' centres (every part's, in part order), and the nearest one to each
    // lattice cell.
    let mut centres = Vec::new();
    let mut sizes = Vec::with_capacity(parts.len());
    for part in &parts {
        let held = (centres.len() as u64).saturating_mul(BUCKETED_CELL_BYTES);
        let fits = budget.room(held).map_or(MAX_OPENFOAM_CELLS, |r| usize::try_from(r / BUCKETED_CELL_BYTES).unwrap_or(usize::MAX));
        let left = MAX_OPENFOAM_CELLS.min(fits).saturating_sub(centres.len());
        match field(&part.join("0/C"), Count::AtMost(left), budget.file_cap(held))? {
            FieldValues::Vectors(v) => {
                sizes.push(v.len());
                centres.extend(v);
            }
            _ => {
                return Err(EngineError::invalid(
                    "path",
                    format!(
                        "{} holds no cell centres (Allrun writes them with `postProcess -func writeCellCentres -time 0`, before decomposePar)",
                        part.join("0/C").display()
                    ),
                ));
            }
        }
    }
    budget.check((centres.len() as u64).saturating_mul(BUCKETED_CELL_BYTES), "the case's cell centres")?;
    let reach = sample_reach(marker.background, marker.spacing);
    let nearest = nearest::nearest_cells(&grid, &solid, &centres, reach).map_err(|f| match f {
        nearest::Failure::Refused(e) => e,
        nearest::Failure::Uncovered { position: [x, y, z], reach } => {
            let len = |m: f64| crate::units::shown(s.units(), crate::units::si(), Quantity::Length, m);
            EngineError::invalid(
                "path",
                format!(
                    "no OpenFOAM cell centre within {} of the lattice cell at ({}, {}, {}): the case does not cover the step's domain there",
                    len(reach),
                    len(x),
                    len(y),
                    len(z)
                ),
            )
        }
    })?;

    let gas = marker.gas.as_ref().map(|g| (g, g.gamma()));
    let mut frames = Vec::with_capacity(chosen.len());
    let mut max_speed = 0.0f64;
    let mut max_mach = 0.0f64;
    for &k in &chosen {
        let Some((time, path)) = written.get(k) else { continue };
        let name = path.file_name().map(PathBuf::from).unwrap_or_default();
        let (_, u) = joined(&parts, &sizes, &name.join("U"), true, budget)?;
        let (pf, _) = joined(&parts, &sizes, &name.join("p"), false, budget)?;
        let tf = match gas {
            Some(_) => Some(joined(&parts, &sizes, &name.join("T"), false, budget)?.0),
            None => None,
        };
        if let Some(tf) = &tf {
            // Every cell of the case, before any is sampled: a failed run leaves one behind.
            let absolute = |what: &str| format!("an absolute {what} must be finite and positive (the run failed there)");
            refuse_unless_positive(s.units(), (&parts, &sizes), &name.join("p"), &pf, "pressure", Quantity::Pressure, &absolute("pressure"))?;
            refuse_unless_positive(s.units(), (&parts, &sizes), &name.join("T"), tf, "temperature", Quantity::Temperature, &absolute("temperature"))?;
            // The density the frame stores (f32): p/(R_s T) of two positive numbers still overflows.
            if let Some((g, _)) = gas {
                let rho: Vec<f64> = pf.iter().zip(tf).map(|(p, t)| f64::from((p / (g.gas_constant * t)) as f32)).collect();
                let need = "the density p/(R T) of the gas must be finite and positive as a single-precision number";
                refuse_unless_positive(s.units(), (&parts, &sizes), &name.join("T"), &rho, "density", Quantity::Density, need)?;
            }
        }
        let mut frame = Frame { time: *time, step: k as u64, ..Frame::default() };
        frame.velocity = Vec::with_capacity(grid.len());
        frame.pressure = Vec::with_capacity(grid.len());
        if tf.is_some() {
            frame.temperature = Vec::with_capacity(grid.len());
            frame.density = Vec::with_capacity(grid.len());
        }
        let mut fluid = 0usize;
        let mut gas_mass = 0.0f64;
        for slot in &nearest {
            let (v, p, t) = match slot {
                Some(id) => {
                    let id = *id as usize;
                    // `joined` read exactly one value per centre, so every id has one.
                    let v = u.get(id).copied().ok_or_else(|| EngineError::invalid("path", "a cell without a velocity"))?;
                    let p = pf.get(id).copied().ok_or_else(|| EngineError::invalid("path", "a cell without a pressure"))?;
                    let t = match &tf {
                        Some(tf) => Some(tf.get(id).copied().ok_or_else(|| EngineError::invalid("path", "a cell without a temperature"))?),
                        None => None,
                    };
                    fluid += 1;
                    (v, p, t)
                }
                // A solid cell: no flow in it, and no gas.
                None => ([0.0; 3], 0.0, tf.as_ref().map(|_| 0.0)),
            };
            let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            max_speed = max_speed.max(speed);
            frame.velocity.push(v.map(|c| c as f32));
            match (gas, t) {
                (Some((g, gamma)), Some(t)) => {
                    let cell = if slot.is_some() { gas_cell(g, gamma, p, t, speed) } else { GasCell::SOLID };
                    max_mach = max_mach.max(cell.mach);
                    gas_mass += cell.density * spacing.powi(3);
                    frame.pressure.push(cell.pressure as f32);
                    frame.temperature.push(cell.temperature as f32);
                    frame.density.push(cell.density as f32);
                }
                _ => frame.pressure.push((p * marker.density) as f32),
            }
        }
        frame.mass = if gas.is_some() { gas_mass } else { marker.density * spacing.powi(3) * fluid as f64 };
        frames.push(frame);
    }

    let fluid_cells = nearest.iter().filter(|c| c.is_some()).count();
    let open = |f: BoxFace| matches!(domain.faces.get(f), FaceCondition::Velocity { .. } | FaceCondition::Pressure { .. });
    let end_time = frames.last().map_or(0.0, |f| f.time);
    let stats = FluidStats {
        validity: None,
        cells: grid.len(),
        fluid_cells,
        steps: 0,
        time_step: 0.0,
        end_time,
        tau_plus: 0.0,
        tau_minus: 0.0,
        magic: 0.0,
        collision: fieldcraft_model::Collision::default(),
        shear_dependent_tau: None,
        velocity: max_speed,
        lattice_velocity: 0.0,
        mach: max_mach,
        cell_updates: 0,
        backend: format!("openfoam ({})", marker.solver),
        tau_thermal: None,
        eddy_viscosity_ratio: None,
        frozen: None,
        free_surface: None,
        wall_links: None,
    };
    let history_times = frames.iter().map(|f| f.time).collect();
    let n_frames = frames.len();
    let result = FluidResult {
        port_properties: Some(PortProperties { density: marker.density, specific_heat: gas.map(|(g, _)| g.specific_heat) }),
        step: step_name.clone(),
        grid,
        solid,
        open_faces: BoxFace::ALL.map(open),
        frames,
        history_times,
        forces: Vec::new(),
        heat: Vec::new(),
        baffles: Box::default(),
        zones: Vec::new(),
        stats,
        warnings: vec![format!("read from the OpenFOAM case {} ({}), not solved on the lattice", root.display(), marker.solver)],
        steady: None,
        temperature_reference: gas.map_or(fieldcraft_model::units::CELSIUS_ZERO, |(g, _)| g.temperature),
        solids: Vec::new(),
        surface: None,
    };
    s.store_result(StepResult::Fluid(Box::new(result)));
    let source = if decomposed { format!("decomposed ({} processors)", parts.len()) } else { "reconstructed".to_string() };
    let mut reply = json!({"step": step_name, "path": dir, "solver": marker.solver, "source": source, "frames": n_frames,
        "written_times": written.len(), "end_time": end_time, "cells": cells, "spacing": spacing, "fluid_cells": fluid_cells,
        "openfoam_cells": centres.len(), "max_speed": max_speed});
    if gas.is_some() {
        reply["max_mach"] = json!(max_mach);
    }
    Ok(reply)
}
