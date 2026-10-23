//! The solved results part of a `.fieldcraft` file (`results.bin`, optional): per solved step
//! what its solver reported, with the hash of the definition it was solved from, so a model
//! opens with its results and shows them without solving again.
//!
//! Little-endian:
//!
//! ```text
//! "MCRSLT01"                     magic, 8 bytes (the results encoding's version is 01)
//! u32 steps
//! steps × { u64 length, step name (UTF-8),
//!           32 bytes: SHA-256 of the definition it was solved from ([`definition_hash`]),
//!           u8 kind (0 static, 1 frequency, 2 heat transfer, 3 explicit dynamic,
//!                    4 legacy fluid, 5 fluid with port properties),
//!           u64 length, the result in the [`crate::codec`] encoding }   names ascending
//! ```
//!
//! The result is what the solver reported, every number in SI: a static step's displacements,
//! rotations, reactions, stresses (nodal, element and integration point, shells and beams),
//! reference points, contacts and nonlinear increments with their frames; a frequency or
//! buckling step's modes; a heat transfer step's temperature frames, fluxes and heat flows; an
//! explicit step's field frames, histories, energies, contacts and erosion; a fluid step's grid
//! (origin, spacing, cells), cell fields per frame (velocity, pressure, temperature), wall
//! forces, heat and zone histories, statistics and steady-state monitors.
//!
//! The payload follows the field and variant order of the solvers' result types, so the
//! encoding's version changes with any change to them ([`MAGIC`]). A file whose results are of
//! another version opens without them ([`Decoded::OtherVersion`]): results are what solving
//! gives again, unlike the definition.
//! Kind 5 adds an envelope around the unchanged kind-4 fluid payload. New kinds do not alter
//! existing payloads; older readers refuse an unknown kind and keep the model definition.
//!
//! Reading treats the part as hostile, as `mesh.bin`: the number of steps is capped, every length
//! is checked against the bytes present before anything is read for it, nothing may be left
//! over, and every result is checked ([`check`]) before it is used: names, ascending labels, and
//! per node, element and cell arrays the lengths of what they follow.

use fieldcraft_explicit::ExplicitResult;
use fieldcraft_lbm::FluidResult;
use fieldcraft_model::Model;
use sha2::{Digest, Sha256};

use crate::{Error, Geometry, MAX_FILE_BYTES, Result, codec, invalid};

/// The magic bytes the part starts with (the `01` is the results encoding's version).
pub const MAGIC: &[u8; 8] = b"MCRSLT01";

/// The part's name.
pub const PART: &str = "results.bin";

/// Largest results part read or written.
pub const MAX_BYTES: u64 = MAX_FILE_BYTES;

/// Most steps a results part holds.
pub const MAX_STEPS: usize = 10_000;

/// Most cells of a stored fluid grid.
pub const MAX_CELLS: usize = u32::MAX as usize;

/// A step's result as a solver reported it.
#[derive(Debug, Clone, PartialEq)]
pub enum Solved {
    /// A static, frequency, buckling or heat transfer step.
    Implicit(fieldcraft_solver::StepResult),
    Explicit(ExplicitResult),
    Fluid(FluidResult),
}

/// A step's result to write, borrowed.
#[derive(Debug, Clone, Copy)]
pub enum SolvedRef<'a> {
    Implicit(&'a fieldcraft_solver::StepResult),
    Explicit(&'a ExplicitResult),
    Fluid(&'a FluidResult),
}

// Keep the original kind-4 field order intact. Kind 5 adds immutable solve-time properties
// outside that payload, so old packages still open and stale results never use edited materials.
#[derive(serde::Serialize)]
struct FluidPortsRef<'a> {
    flow: &'a FluidResult,
    properties: Option<fieldcraft_lbm::PortProperties>,
}

#[derive(serde::Deserialize)]
struct FluidPorts {
    flow: FluidResult,
    properties: Option<fieldcraft_lbm::PortProperties>,
}

impl Solved {
    pub fn as_ref(&self) -> SolvedRef<'_> {
        match self {
            Solved::Implicit(r) => SolvedRef::Implicit(r),
            Solved::Explicit(r) => SolvedRef::Explicit(r),
            Solved::Fluid(r) => SolvedRef::Fluid(r),
        }
    }

    pub fn step(&self) -> &str {
        self.as_ref().step()
    }
}

impl<'a> SolvedRef<'a> {
    pub fn step(&self) -> &'a str {
        match self {
            SolvedRef::Implicit(r) => r.step(),
            SolvedRef::Explicit(r) => &r.step,
            SolvedRef::Fluid(r) => &r.step,
        }
    }

    /// The kind code of the format.
    fn kind(&self) -> u8 {
        match self {
            SolvedRef::Implicit(fieldcraft_solver::StepResult::Static(_)) => 0,
            SolvedRef::Implicit(fieldcraft_solver::StepResult::Frequency(_)) => 1,
            SolvedRef::Implicit(fieldcraft_solver::StepResult::Heat(_)) => 2,
            SolvedRef::Explicit(_) => 3,
            SolvedRef::Fluid(_) => 5,
        }
    }

    fn encode_into(&self, out: &mut Vec<u8>) -> std::result::Result<(), codec::CodecError> {
        match self {
            SolvedRef::Implicit(fieldcraft_solver::StepResult::Static(r)) => codec::encode_into(r, out),
            SolvedRef::Implicit(fieldcraft_solver::StepResult::Frequency(r)) => codec::encode_into(r, out),
            SolvedRef::Implicit(fieldcraft_solver::StepResult::Heat(r)) => codec::encode_into(r, out),
            SolvedRef::Explicit(r) => codec::encode_into(r, out),
            SolvedRef::Fluid(r) => codec::encode_into(&FluidPortsRef { flow: r, properties: r.port_properties }, out),
        }
    }
}

/// A step's result as a file holds it: the definition it was solved from, by its hash.
#[derive(Debug, Clone, PartialEq)]
pub struct Stored {
    /// [`definition_hash`] when the step was solved.
    pub definition: [u8; 32],
    pub result: Solved,
}

/// A step's result to write.
#[derive(Debug, Clone, Copy)]
pub struct StoredRef<'a> {
    pub definition: [u8; 32],
    pub result: SolvedRef<'a>,
}

/// What a results part gave.
#[derive(Debug, Clone, PartialEq)]
pub enum Decoded {
    /// Every step's result, by name.
    Results(Vec<Stored>),
    /// The results were written in another results encoding (its magic, printable): the file
    /// opens without them.
    OtherVersion(String),
}

/// The SHA-256 of what solving step `step` of `model` reads ([`Model::step_inputs`]: the mesh,
/// the definitions the step uses and the step itself) and of the geometry files the model is
/// built on (by their content's hash; a fluid step's obstacles come from them), or `None` when
/// there is no such step. A result whose step's hash is the same describes the model as it is.
pub fn definition_hash(model: &Model, step: &str, geometry: &Geometry) -> Option<[u8; 32]> {
    /// Hashes what is written to it.
    struct HashWriter(Sha256);
    impl std::io::Write for HashWriter {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.update(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let inputs = model.step_inputs(step)?;
    let mut h = HashWriter(Sha256::new());
    serde_json::to_writer(&mut h, &inputs).ok()?;
    for source in geometry.design.iter().chain(&geometry.meshes) {
        h.0.update(b"\0");
        h.0.update(source.sha256().as_bytes());
    }
    Some(h.0.finalize().into())
}

/// The results part holding `results` (in step name order; names must be unique).
pub fn encode(results: &[(&str, StoredRef<'_>)]) -> Result<Vec<u8>> {
    if results.len() > MAX_STEPS {
        return Err(Error::Limit(format!("{PART}: at most {MAX_STEPS} steps' results")));
    }
    let mut sorted: Vec<&(&str, StoredRef<'_>)> = results.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    if sorted.windows(2).any(|w| w.first().map(|a| a.0) == w.get(1).map(|b| b.0)) {
        return invalid(format!("{PART}: a step's result is listed twice"));
    }
    let mut out = MAGIC.to_vec();
    out.extend_from_slice(&u32::try_from(sorted.len()).unwrap_or(u32::MAX).to_le_bytes());
    for (name, stored) in sorted {
        if stored.result.step() != *name {
            return invalid(format!("{PART}: the result of step `{name}` is step `{}`'s", stored.result.step()));
        }
        out.extend_from_slice(&(name.len() as u64).to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&stored.definition);
        out.push(stored.result.kind());
        // The payload is encoded in place, its length written before it once known.
        let at = out.len();
        out.extend_from_slice(&[0; 8]);
        stored.result.encode_into(&mut out).map_err(|e| Error::Invalid(format!("{PART}: step `{name}`: {e}")))?;
        let len = (out.len() - at - 8) as u64;
        if let Some(slot) = out.get_mut(at..at + 8) {
            slot.copy_from_slice(&len.to_le_bytes());
        }
        if out.len() as u64 > MAX_BYTES {
            return Err(Error::Limit(format!(
                "the results are more than {} MB (at most that in a file): save without them (`results: false`)",
                MAX_BYTES >> 20
            )));
        }
    }
    Ok(out)
}

/// Reads the part's fields in order, each length checked against the bytes left.
struct Reader<'a> {
    bytes: &'a [u8],
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: u64, what: &str) -> Result<&'a [u8]> {
        let n = usize::try_from(n)
            .ok()
            .filter(|n| *n <= self.bytes.len())
            .ok_or_else(|| Error::Invalid(format!("{PART} is damaged: {what} of {n} bytes where {} are left", self.bytes.len())))?;
        let (head, rest) = self.bytes.split_at(n);
        self.bytes = rest;
        Ok(head)
    }

    fn array<const N: usize>(&mut self, what: &str) -> Result<[u8; N]> {
        let mut a = [0u8; N];
        a.copy_from_slice(self.take(N as u64, what)?);
        Ok(a)
    }
}

/// The results part `bytes` read and checked.
pub fn decode(bytes: &[u8]) -> Result<Decoded> {
    if bytes.len() as u64 > MAX_BYTES {
        return Err(Error::Limit(format!("{PART} is {} MB (at most {} MB)", bytes.len() >> 20, MAX_BYTES >> 20)));
    }
    let mut r = Reader { bytes };
    let magic: [u8; 8] = r.array("the magic")?;
    if &magic != MAGIC {
        if magic.starts_with(b"MCRSLT") {
            return Ok(Decoded::OtherVersion(magic.iter().map(|b| if b.is_ascii_graphic() { *b as char } else { '?' }).collect()));
        }
        return invalid(format!("{PART} is not a results part (no {} magic)", String::from_utf8_lossy(MAGIC)));
    }
    let count = u32::from_le_bytes(r.array("the step count")?) as usize;
    if count > MAX_STEPS {
        return Err(Error::Limit(format!("{PART} holds {count} steps' results (at most {MAX_STEPS})")));
    }
    let mut out: Vec<Stored> = Vec::new();
    let mut last: Option<&str> = None;
    for _ in 0..count {
        let len = u64::from_le_bytes(r.array("a step name's length")?);
        let name = std::str::from_utf8(r.take(len, "a step name")?).map_err(|_| Error::Invalid(format!("{PART}: a step name is not UTF-8")))?;
        fieldcraft_model::check_name("step", name).map_err(|e| Error::Invalid(format!("{PART}: {e}")))?;
        if last.is_some_and(|l| l >= name) {
            return invalid(format!("{PART}: steps must be in name order, each once (`{name}`)"));
        }
        last = Some(name);
        let definition: [u8; 32] = r.array("a definition hash")?;
        let [kind] = r.array::<1>("a result kind")?;
        let len = u64::from_le_bytes(r.array("a result's length")?);
        let payload = r.take(len, "a result")?;
        let damaged = |e: codec::CodecError| Error::Invalid(format!("{PART}: the result of step `{name}` is damaged: {e}"));
        let result = match kind {
            0 => Solved::Implicit(fieldcraft_solver::StepResult::Static(codec::decode(payload).map_err(damaged)?)),
            1 => Solved::Implicit(fieldcraft_solver::StepResult::Frequency(codec::decode(payload).map_err(damaged)?)),
            2 => Solved::Implicit(fieldcraft_solver::StepResult::Heat(codec::decode(payload).map_err(damaged)?)),
            3 => Solved::Explicit(codec::decode(payload).map_err(damaged)?),
            4 => Solved::Fluid(codec::decode(payload).map_err(damaged)?),
            5 => {
                let mut stored: FluidPorts = codec::decode(payload).map_err(damaged)?;
                stored.flow.port_properties = stored.properties;
                Solved::Fluid(stored.flow)
            }
            k => return invalid(format!("{PART}: step `{name}` has a result of unknown kind {k}")),
        };
        if result.step() != name {
            return invalid(format!(
                "{PART}: the result listed for step `{name}` is step `{}`'s",
                result.step().chars().take(80).collect::<String>()
            ));
        }
        check(&result).map_err(|why| Error::Invalid(format!("{PART}: the result of step `{name}`: {why}")))?;
        out.push(Stored { definition, result });
    }
    if !r.bytes.is_empty() {
        return invalid(format!("{PART} has {} trailing bytes", r.bytes.len()));
    }
    Ok(Decoded::Results(out))
}

/// Whether labels strictly ascend (results look them up by binary search).
fn ascending(labels: &[u32]) -> bool {
    labels.windows(2).all(|w| matches!(w, [a, b] if a < b))
}

/// `Err(what)` unless every array `lens` names is `n` long (or, for those marked optional,
/// empty).
fn lengths(n: usize, lens: &[(&str, usize, bool)]) -> std::result::Result<(), String> {
    match lens.iter().find(|(_, len, optional)| *len != n && !(*optional && *len == 0)) {
        Some((what, len, _)) => Err(format!("{what} has {len} values for {n}")),
        None => Ok(()),
    }
}

/// A stored result checked for what reading it relies on: ascending labels, arrays of the
/// lengths of what they follow, and a fluid grid that is a grid.
pub fn check(r: &Solved) -> std::result::Result<(), String> {
    use fieldcraft_solver::StepResult as S;
    match r {
        Solved::Implicit(S::Static(s)) => {
            if !ascending(&s.nodes) {
                return Err("its nodes do not ascend".into());
            }
            let n = s.nodes.len();
            lengths(
                n,
                &[
                    ("the displacement", s.displacement.len(), false),
                    ("the rotation", s.rotation.len(), true),
                    ("the reaction", s.reaction.len(), true),
                    ("the reaction moment", s.reaction_moment.len(), true),
                    ("the stress", s.stress.len(), true),
                    ("the von Mises stress", s.von_mises.len(), true),
                    ("the temperature", s.temperature.as_ref().map_or(0, Vec::len), true),
                ],
            )?;
            if let Some(h) = &s.nonlinear {
                for f in &h.frames {
                    lengths(n, &[("a frame's displacement", f.displacement.len(), false)])?;
                    if !f.reaction.windows(2).all(|w| matches!(w, [a, b] if a.0 < b.0)) || f.reaction.iter().any(|(i, _)| *i >= n) {
                        return Err("a frame's reactions are not at ascending nodes of the result".into());
                    }
                }
            }
            for c in &s.contacts {
                if !ascending(&c.nodes) {
                    return Err(format!("contact `{}`'s nodes do not ascend", c.name));
                }
                let m = c.nodes.len();
                lengths(
                    m,
                    &[
                        ("a contact's pressure", c.pressure.len(), false),
                        ("a contact's gap", c.gap.len(), false),
                        ("a contact's status", c.status.len(), false),
                        ("a contact's multipliers", c.multiplier.len(), true),
                        ("a contact's shear", c.shear.len(), true),
                        ("a contact's slip", c.slip.len(), true),
                    ],
                )?;
            }
            Ok(())
        }
        Solved::Implicit(S::Frequency(f)) => {
            if !ascending(&f.nodes) {
                return Err("its nodes do not ascend".into());
            }
            for m in &f.modes {
                lengths(f.nodes.len(), &[("a mode's shape", m.shape.len(), false), ("a mode's rotation", m.rotation.len(), true)])?;
            }
            Ok(())
        }
        Solved::Implicit(S::Heat(h)) => {
            if !ascending(&h.nodes) {
                return Err("its nodes do not ascend".into());
            }
            if h.times.len() != h.frames.len() {
                return Err(format!("{} frame times for {} frames", h.times.len(), h.frames.len()));
            }
            for f in &h.frames {
                lengths(h.nodes.len(), &[("a temperature frame", f.len(), false)])?;
            }
            lengths(h.nodes.len(), &[("the heat flux", h.flux.len(), true), ("the reaction", h.reaction.len(), true)])
        }
        Solved::Explicit(e) => {
            if !ascending(&e.nodes) || !ascending(&e.elements) {
                return Err("its node or element labels do not ascend".into());
            }
            for f in &e.frames {
                lengths(e.nodes.len(), &[("a frame's displacement", f.displacement.len(), true), ("a frame's velocity", f.velocity.len(), true)])?;
                lengths(
                    e.elements.len(),
                    &[
                        ("a frame's stress", f.stress.len(), true),
                        ("a frame's strain", f.strain.len(), true),
                        ("a frame's plastic strain", f.plastic_strain.len(), true),
                        ("a frame's damage", f.damage.len(), true),
                        ("a frame's temperature", f.temperature.len(), true),
                        ("a frame's erosion", f.failed.len(), true),
                    ],
                )?;
            }
            Ok(())
        }
        Solved::Fluid(f) => check_fluid(f),
    }
}

fn check_fluid(f: &FluidResult) -> std::result::Result<(), String> {
    if f.port_properties
        .is_some_and(|p| !(p.density > 0.0 && p.density.is_finite()) || p.specific_heat.is_some_and(|cp| !(cp > 0.0 && cp.is_finite())))
    {
        return Err("its port properties need a finite positive density and specific heat".into());
    }
    let g = &f.grid;
    if !(g.spacing.is_finite() && g.spacing > 0.0) || g.origin.iter().any(|x| !x.is_finite()) {
        return Err("its grid has no finite origin and positive spacing".into());
    }
    let cells = g.cells.iter().try_fold(1usize, |n, c| n.checked_mul(*c).filter(|_| *c > 0)).filter(|n| *n <= MAX_CELLS);
    let Some(n) = cells else { return Err(format!("its grid of {:?} cells is not 1 to {MAX_CELLS} cells", g.cells)) };
    lengths(n, &[("the solid cells", f.solid.len(), false)])?;
    for frame in &f.frames {
        lengths(
            n,
            &[
                ("a frame's velocity", frame.velocity.len(), false),
                ("a frame's pressure", frame.pressure.len(), false),
                ("a frame's temperature", frame.temperature.len(), true),
            ],
        )?;
    }
    let samples = f.history_times.len();
    let histories = f.forces.iter().map(|h| ("a force history", h.force.len())).chain(f.heat.iter().map(|h| ("a heat history", h.power.len())));
    for (what, len) in histories {
        lengths(samples, &[(what, len, false)])?;
    }
    for z in &f.zones {
        lengths(
            samples,
            &[("a zone's flow", z.flow.len(), false), ("a zone's pressure", z.pressure.len(), false), ("a fan's rise", z.fan_rise.len(), true)],
        )?;
    }
    if f.solids.iter().any(|s| s.cells.iter().any(|c| *c as usize >= n)) {
        return Err("a conducting solid has a cell outside the grid".into());
    }
    Ok(())
}
