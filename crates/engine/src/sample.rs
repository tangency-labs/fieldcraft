//! Sample models as scripts of ordinary commands, so they build the same way a user's script
//! does.
//!
//! [`script`], the model the desktop app opens with `--sample`: two parts side by side, a
//! 40 × 20 × 10 block of 5 mm hex8 cells, and a 15 mm cube split into 162 tet4 elements (each
//! 5 mm cell into six tetrahedra around its main diagonal, the Kuhn subdivision, which is
//! conforming between neighbouring cells). Every element is positively oriented.
//!
//! [`cantilever`], a model ready to solve: a steel beam clamped at one end with a load on the
//! other, a static and a frequency step. [`impact`]: the same beam flying onto a rigid wall, an
//! explicit dynamic step.
//!
//! The examples gallery (docs/examples.md), each ready to solve with `solve.run`:
//! [`ibeam`], [`drop_test`], [`rubber_seal`], [`flow_past`], [`heated_bracket`], [`bracket_hex`] and [`crank_arm`]. Their
//! geometry is a SolidCraft script run when the sample opens ([`Sample::design`]): a corpus
//! recipe (`tests/corpus`) or [`UPRIGHT_BRACKET`].

use fieldcraft_geometry::Shape;
use fieldcraft_model::units::{DisplaySystem, Preset, Quantity, Typed};
use serde_json::{Value, json};

mod cooling;
mod gpu_card;
mod heated_flow;

/// A sample the app offers (Help › Samples, the empty viewport, `file.sample`, `--sample NAME`).
pub struct Sample {
    pub name: &'static str,
    pub label: &'static str,
    /// One line: what it shows.
    pub about: &'static str,
    /// The SolidCraft script (`{"commands": [...]}`, the corpus recipe form) that builds the
    /// sample's geometry, opened as its design before [`Sample::script`] runs.
    pub design: Option<&'static str>,
    /// Commands run on the opened design before [`Sample::script`] (datums and partitions), so
    /// the script is given the bodies as they cut it.
    pub prepare: Option<fn() -> Value>,
    /// The model's commands, given the design's bodies (none without a design).
    pub script: fn(&[Shape]) -> Value,
}

/// Every sample, the default first.
pub const SAMPLES: &[Sample] = &[
    Sample {
        name: "parts",
        label: "Hex and Tet Parts",
        about: "A hex8 block beside a tet4 cube: meshes to look at, nothing to solve",
        design: None,
        prepare: None,
        script: |_| script(),
    },
    Sample {
        name: "cantilever",
        label: "Cantilever",
        about: "A clamped steel beam with an end load: a static and a frequency step, ready to solve",
        design: None,
        prepare: None,
        script: |_| cantilever(),
    },
    Sample {
        name: "impact",
        label: "Bar Impact",
        about: "The steel beam flying onto a rigid wall: an explicit dynamic step, ready to solve",
        design: None,
        prepare: None,
        script: |_| impact(),
    },
    Sample {
        name: "ibeam",
        label: "I-Beam Cantilever",
        about: "A 200 mm steel I-beam swept all-hex, clamped at one end with 1 kN on the other through a coupling: static, ready to solve",
        design: Some(CANTILEVER_I_BEAM),
        prepare: None,
        script: |_| ibeam(),
    },
    Sample {
        name: "drop-test",
        label: "Drop Test",
        about: "An aluminium L-bracket dropped 2 m onto a tilted rigid floor: explicit dynamics, the strike and its stress waves",
        design: Some(UPRIGHT_BRACKET),
        prepare: None,
        script: |_| drop_test(),
    },
    Sample {
        name: "rubber-seal",
        label: "Rubber Seal",
        about: "A rubber gasket ring squeezed 25 % between bonded plates: nonlinear static, Mooney–Rivlin, large strain",
        design: None,
        prepare: None,
        script: |_| rubber_seal(),
    },
    Sample {
        name: "flow-past",
        label: "Flow Past a Bracket",
        about: "Air flowing past an L-bracket on the floor: lattice Boltzmann flow and the force on each face",
        design: Some(UPRIGHT_BRACKET),
        prepare: None,
        script: |_| flow_past(),
    },
    Sample {
        name: "heated-flow",
        label: "Heated Flow Past a Bracket",
        about: "Cool air passing an 80 °C bracket: temperature slices and signed wall heat, a small transient demonstration",
        design: Some(UPRIGHT_BRACKET),
        prepare: None,
        script: heated_flow::script,
    },
    Sample {
        name: "cooling-heat-sink",
        label: "Heat Sink in Forced Air",
        about: "A six-fin heat sink held at 50 °C in a duct of 25 °C air at 0.1 m/s, buoyancy on: thermal lattice Boltzmann, the heat it sheds and the air's rise (minutes on a CPU)",
        design: Some(cooling::HEAT_SINK_DESIGN),
        prepare: None,
        script: cooling::heat_sink,
    },
    Sample {
        name: "cooling-enclosure",
        label: "Sealed Enclosure Cooling",
        about: "A 40 °C chip on a PCB in a sealed box with 25 °C walls, no fan: the buoyant plume and the heat natural convection carries, to steady state (minutes on a CPU)",
        design: Some(cooling::ENCLOSURE_DESIGN),
        prepare: None,
        script: cooling::enclosure,
    },
    Sample {
        name: "cooling-aio",
        label: "AIO Cooler Radiator",
        about: "Fan air through the half of an AIO cooler's radiator held at a 45 °C coolant temperature: the heat the air takes and the exchanger effectiveness (5 M cells: solve on a GPU)",
        design: Some(cooling::AIO_HALF_DESIGN),
        prepare: None,
        script: cooling::aio,
    },
    Sample {
        name: "cooling-gpu-card",
        label: "Graphics Card Cooling",
        about: "A blower-cooled two-slot graphics card: 275 W conducting from the die through a copper base, heat pipes and fins into the air, with memory and VRMs (conjugate heat transfer; hours on a CPU)",
        design: Some(gpu_card::DESIGN),
        prepare: None,
        script: gpu_card::card,
    },
    Sample {
        name: "heated-bracket",
        label: "Heated Bracket",
        about: "A steel bracket clamped at its holes and heated at one end: steady heat conduction and the thermal stress",
        design: Some(FILLETED_BRACKET),
        prepare: None,
        script: |_| heated_bracket(),
    },
    Sample {
        name: "bracket-hex",
        label: "Hex-Meshed Bracket",
        about: "The filleted bracket split into 30 cells and meshed all-hex (hex20): clamped, loaded, von Mises",
        design: Some(PARTITIONED_BRACKET),
        prepare: None,
        script: bracket_hex,
    },
    Sample {
        name: "crank-mirror",
        label: "Mirrored Crank Arm",
        about: "The crank arm of crank-arm under the name of its mirror workflow: a partitioned half meshed structured all-hex (hex20) and mirrored into one conforming mesh",
        design: Some(CRANK_ARM),
        prepare: Some(crank_prepare),
        script: crank_arm,
    },
    Sample {
        name: "crank-arm",
        label: "Bicycle Crank Arm",
        about: "An aluminium bicycle crank meshed structured all-hex (hex20) from a partitioned half and its mirror image: square bore fixed, a rider's weight on the pedal, von Mises",
        design: Some(CRANK_ARM),
        prepare: Some(crank_prepare),
        script: crank_arm,
    },
];

/// The sample called `name`.
pub fn find(name: &str) -> Option<&'static Sample> {
    SAMPLES.iter().find(|s| s.name == name)
}

/// Hex block: cells along x, y, z and the cell size (mm).
const HEX_CELLS: [u32; 3] = [8, 4, 2];
/// Tet cube: cells per side and its offset from the origin (mm).
const TET_CELLS: u32 = 3;
const TET_ORIGIN: [f64; 3] = [50.0, 2.5, 0.0];
const CELL: f64 = 5.0;

/// The script: `mesh.node.add` and `mesh.element.add` steps.
pub fn script() -> Value {
    let mut steps = Vec::new();
    let mut node = 0u32;
    let mut element = 0u32;
    // Hex block.
    let [nx, ny, nz] = HEX_CELLS;
    let first = node + 1;
    let hid = |i: u32, j: u32, k: u32| first + i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                node += 1;
                steps.push(json!({"command": "mesh.node.add", "params": {"label": node, "xyz": [f64::from(i) * CELL, f64::from(j) * CELL, f64::from(k) * CELL]}}));
            }
        }
    }
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                element += 1;
                let nodes = [
                    hid(i, j, k),
                    hid(i + 1, j, k),
                    hid(i + 1, j + 1, k),
                    hid(i, j + 1, k),
                    hid(i, j, k + 1),
                    hid(i + 1, j, k + 1),
                    hid(i + 1, j + 1, k + 1),
                    hid(i, j + 1, k + 1),
                ];
                steps.push(json!({"command": "mesh.element.add", "params": {"label": element, "kind": "hex8", "nodes": nodes}}));
            }
        }
    }
    // Tet cube.
    let n = TET_CELLS;
    let first = node + 1;
    let tid = |i: u32, j: u32, k: u32| first + i + j * (n + 1) + k * (n + 1) * (n + 1);
    let at = |i: u32, j: u32, k: u32| [TET_ORIGIN[0] + f64::from(i) * CELL, TET_ORIGIN[1] + f64::from(j) * CELL, TET_ORIGIN[2] + f64::from(k) * CELL];
    for k in 0..=n {
        for j in 0..=n {
            for i in 0..=n {
                node += 1;
                steps.push(json!({"command": "mesh.node.add", "params": {"label": node, "xyz": at(i, j, k)}}));
            }
        }
    }
    // The six axis orders: each walks from the cell's corner (0,0,0) to (1,1,1) one axis at a time.
    const ORDERS: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                for order in ORDERS {
                    let mut c = [0u32; 3];
                    let mut corners = vec![c];
                    for axis in order {
                        if let Some(x) = c.get_mut(axis) {
                            *x = 1;
                        }
                        corners.push(c);
                    }
                    let mut nodes: Vec<u32> = corners.iter().map(|c| tid(i + c[0], j + c[1], k + c[2])).collect();
                    let p: Vec<[f64; 3]> = corners.iter().map(|c| at(i + c[0], j + c[1], k + c[2])).collect();
                    if signed_volume(&p) < 0.0 {
                        nodes.swap(1, 2);
                    }
                    element += 1;
                    steps.push(json!({"command": "mesh.element.add", "params": {"label": element, "kind": "tet4", "nodes": nodes}}));
                }
            }
        }
    }
    Value::Array(steps)
}

// The samples type their geometry in mm-t-s-N-MPa (each script starts with `model.units`) and
// their physical constants with their own units, converted through the table like any input.

/// The samples' system, for a constant a script computes with ([`sample_value`]).
fn sample_units() -> DisplaySystem {
    DisplaySystem::new(Preset::MmTonneS)
}

/// `t`, a value of quantity `q`, in the samples' system (mm-t-s-N-MPa). A constant of a known
/// unit always converts; a NaN would be refused by the command it went into.
pub fn sample_value(q: Quantity, t: Typed) -> f64 {
    t.express(q, &sample_units()).unwrap_or(f64::NAN)
}

/// Structural steel: E = 210 GPa (ν = 0.3), ρ = 7850 kg/m³.
const STEEL_YOUNG: Typed = Typed::new(210.0, "GPa");
const STEEL_DENSITY: Typed = Typed::new(7850.0, "kg/m^3");

/// The cantilever's length (x), width (y) and height (z) in mm, and its hex8 cells along each.
pub const CANTILEVER_SIZE: [f64; 3] = [100.0, 10.0, 20.0];
const CANTILEVER_CELLS: [u32; 3] = [20, 2, 4];

/// The beam's mesh in steel ([`STEEL_YOUNG`], [`STEEL_DENSITY`]): nodes,
/// hex8 elements, sets `beam` (elements), `root` (nodes at x = 0) and `tip` (at x = L), the
/// material `Steel` and a solid section `Beam` with the given formulation.
fn beam(formulation: Value) -> Vec<Value> {
    let [nx, ny, nz] = CANTILEVER_CELLS;
    let [lx, ly, lz] = CANTILEVER_SIZE;
    let id = |i: u32, j: u32, k: u32| 1 + i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
    let mut steps = vec![json!({"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}})];
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                let xyz = [lx * f64::from(i) / f64::from(nx), ly * f64::from(j) / f64::from(ny), lz * f64::from(k) / f64::from(nz)];
                steps.push(json!({"command": "mesh.node.add", "params": {"label": id(i, j, k), "xyz": xyz}}));
            }
        }
    }
    let mut element = 0u32;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                element += 1;
                let nodes = [
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i + 1, j + 1, k + 1),
                    id(i, j + 1, k + 1),
                ];
                steps.push(json!({"command": "mesh.element.add", "params": {"label": element, "kind": "hex8", "nodes": nodes}}));
            }
        }
    }
    let end = |i: u32| -> Vec<u32> { (0..=nz).flat_map(|k| (0..=ny).map(move |j| id(i, j, k))).collect() };
    steps.extend([
        json!({"command": "model.set.create", "params": {"name": "beam", "elements": (1..=element).collect::<Vec<_>>()}}),
        json!({"command": "model.set.create", "params": {"name": "root", "nodes": end(0)}}),
        json!({"command": "model.set.create", "params": {"name": "tip", "nodes": end(nx)}}),
        json!({"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": STEEL_YOUNG, "poisson": 0.3}, "density": STEEL_DENSITY}}),
        json!({"command": "model.section.create", "params": {"name": "Beam", "set": "beam", "material": "Steel", "kind": {"type": "solid"}, "formulation": formulation}}),
    ]);
    steps
}

/// A steel cantilever of hex8 elements with incompatible modes ([`CANTILEVER_SIZE`]), clamped at
/// x = 0 (set `root`), with a total force of 1000 N down on the free end (set `tip`); step
/// `Load` (static) and step `Modes` (six modes). Euler–Bernoulli puts the tip
/// 4 P L³ / (E w h³) = 0.238 mm down.
pub fn cantilever() -> Value {
    let mut steps = beam(json!({"incompatible_modes": true}));
    steps.extend([
        json!({"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}}),
        json!({"command": "model.load.create", "params": {"name": "Push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -1000]}}),
        json!({"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Push"]}}),
        json!({"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 6}, "boundary_conditions": ["Clamp"]}}),
    ]);
    Value::Array(steps)
}

/// The same steel beam in one-point hex8 elements, flying at 10 m/s onto a rigid wall at x = 0:
/// step `Impact` (explicit dynamics, 40 µs, ten field frames, energies and the tip's history).
/// The compression wave runs the 100 mm bar in about 19 µs (c = √(E/ρ) ≈ 5.2 mm/µs).
pub fn impact() -> Value {
    let mut steps = beam(json!({"integration": "reduced"}));
    let step = json!({
        "type": "explicit_dynamic", "end_time": 4e-5, "output_interval": 4e-6, "history_interval": 2e-7, "history_nodes": "tip",
        "initial_velocities": [{"set": "all_nodes", "velocity": [-10000, 0, 0]}],
        "rigid_walls": [{"set": "root", "point": [0, 0, 0], "normal": [1, 0, 0]}]
    });
    let [nx, ny, nz] = CANTILEVER_CELLS;
    let nodes = (nx + 1) * (ny + 1) * (nz + 1);
    steps.extend([
        json!({"command": "model.set.create", "params": {"name": "all_nodes", "nodes": (1..=nodes).collect::<Vec<_>>()}}),
        json!({"command": "model.step.create", "params": {"name": "Impact", "procedure": step}}),
    ]);
    Value::Array(steps)
}

/// The 200 mm I-beam lying along x, web upright (`tests/corpus/cantilever-i-beam`): 40 × 40 mm,
/// 3 mm flanges, a 2 mm web. Body `Beam`; its end faces are `F2:start` (x = 0) and `F2:top`
/// (x = 200).
pub const CANTILEVER_I_BEAM: &str = include_str!("../../../tests/corpus/cantilever-i-beam/recipe.json");

/// The I-beam's tip load.
pub const IBEAM_LOAD: Typed = Typed::new(1.0, "kN");

/// [`CANTILEVER_I_BEAM`] swept all-hex (hex8 of 4 mm, incompatible modes) in steel, clamped at
/// x = 0 (set `root`), [`IBEAM_LOAD`] down (−z) on a reference point on the beam's axis at the
/// free end, kinematically coupled to the end face (set `tip`) so the section stays plane. Step
/// `Load` (linear static). The README demo's model (tools/demo-video).
pub fn ibeam() -> Value {
    let face = |f: &str| json!([{"body": "Beam", "face": f}]);
    let load = sample_value(Quantity::Force, IBEAM_LOAD);
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "mesh.hex.sweep", "params": {"body": "Beam", "size": 4}},
        {"command": "model.set.create", "params": {"name": "beam", "geometry": [{"body": "Beam"}]}},
        {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": STEEL_YOUNG, "poisson": 0.3}, "density": STEEL_DENSITY}},
        {"command": "model.section.create", "params": {"name": "Beam", "set": "beam", "material": "Steel", "kind": {"type": "solid"}, "formulation": {"incompatible_modes": true}}},
        {"command": "model.set.create", "params": {"name": "root", "geometry": face("F2:start")}},
        {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}},
        {"command": "model.set.create", "params": {"name": "tip", "geometry": face("F2:top")}},
        {"command": "model.reference_point.create", "params": {"name": "End", "xyz": [200, 0, 20]}},
        {"command": "model.coupling.create", "params": {"name": "EndPlate", "reference_point": "End", "set": "tip", "kind": "kinematic"}},
        {"command": "model.load.create", "params": {"name": "Push", "type": "force", "target": {"reference_point": "End"}, "vector": [0, 0, -load]}},
        {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Push"]}},
    ])
}

/// An L-bracket standing on its foot, Z up, as a SolidCraft script: a 60 × 40 mm foot and a 50 mm
/// upright, both 8 mm thick, with an R6 fillet inside the bend (the corpus `filleted-bracket`'s
/// profile without its holes, sketched on XZ so it sweeps along y). Body `Bracket`; its faces
/// `F2:side:l1` (the foot's underside, z = 0), `l2` (the foot's end, x = 60), `l3` (the foot's top),
/// `l4` (the upright's inside, x = 8), `l5` (the upright's top, z = 50) and `l6` (its back, x = 0).
pub const UPRIGHT_BRACKET: &str = r#"{"commands": [
    {"command": "sketch.create", "params": {"plane": "XZ", "name": "Profile"}},
    {"command": "sketch.line", "params": {"points": [[0, 0], [60, 0], [60, -8], [8, -8], [8, -50], [0, -50]], "closed": true}},
    {"command": "sketch.finish", "params": {}},
    {"command": "solid.extrude", "params": {"distance": 40, "body_name": "Bracket"}},
    {"command": "solid.fillet", "params": {"edges": [[8, 20, 8]], "radius": 6}}
]}"#;

/// Aluminium (E = 70 GPa, ν = 0.33, ρ = 2700 kg/m³, yield 240 MPa with 700 MPa linear hardening).
fn aluminium(name: &str) -> Value {
    json!({"command": "model.material.create", "params": {"name": name, "elastic": {"type": "isotropic", "young": "70 GPa", "poisson": 0.33}, "density": ALUMINIUM_DENSITY,
        "plastic": {"type": "bilinear", "yield_stress": 240, "hardening_modulus": 700, "hardening": "isotropic"}}})
}

/// Aluminium's density.
pub const ALUMINIUM_DENSITY: Typed = Typed::new(2700.0, "kg/m^3");

/// Gravity at the Earth's surface, as the samples take it.
pub const GRAVITY: Typed = Typed::new(9.81, "m/s^2"); // unit-factors: the samples' gravity, written in its unit
/// The drop test's height of fall, mm: the bracket meets the floor at √(2 g h) = 6.26 m/s.
pub const DROP_HEIGHT: f64 = 2000.0;
/// The drop test's floor tilt about y, degrees: the foot's end strikes first.
const DROP_TILT: f64 = 15.0;

/// [`UPRIGHT_BRACKET`] in aluminium, hex8 elements swept 3 mm (one-point, stiffness hourglass
/// control), falling at √(2 g h) onto a rigid floor tilted 15° so the foot's end lands first,
/// 0.1 mm below that corner; gravity on. Step `Drop` (explicit dynamics, 1 ms: the strike, the
/// stress waves through the bend and the bracket pivoting onto its foot; 20 field frames, the
/// upright's top (set `top`) in the histories).
pub fn drop_test() -> Value {
    let gravity = sample_value(Quantity::Acceleration, GRAVITY);
    let speed = (2.0 * gravity * DROP_HEIGHT).sqrt();
    let tilt = DROP_TILT.to_radians();
    let step = json!({
        "type": "explicit_dynamic", "end_time": "1 ms", "output_interval": "50 µs", "history_interval": "5 µs", "history_nodes": "top",
        "initial_velocities": [{"set": "bracket", "velocity": [0, 0, -speed]}],
        "rigid_walls": [{"set": "bracket", "point": [60, 0, -0.1], "normal": [-tilt.sin(), 0, tilt.cos()]}],
        "hourglass": {"type": "stiffness", "coefficient": 0.05}
    });
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "mesh.hex.generate", "params": {"body": "Bracket", "size": 3, "method": "extrude"}},
        {"command": "model.set.create", "params": {"name": "bracket", "geometry": [{"body": "Bracket"}]}},
        {"command": "model.set.create", "params": {"name": "top", "geometry": [{"body": "Bracket", "face": "F2:side:l5"}]}},
        aluminium("Aluminium"),
        {"command": "model.section.create", "params": {"name": "Bracket", "set": "bracket", "material": "Aluminium", "kind": {"type": "solid"}, "formulation": {"integration": "reduced"}}},
        {"command": "model.load.create", "params": {"name": "Gravity", "type": "gravity", "acceleration": [0, 0, -gravity]}},
        {"command": "model.step.create", "params": {"name": "Drop", "procedure": step, "loads": ["Gravity"]}},
    ])
}

/// The rubber seal: inner and outer radius and height of the ring (mm), and how far it is squeezed.
pub const SEAL: [f64; 3] = [20.0, 30.0, 10.0];
pub const SEAL_SQUEEZE: f64 = 2.5;

/// A rubber gasket ring ([`SEAL`]: 20 mm inside radius, 30 mm outside, 10 mm high) of hybrid hex8
/// elements, four mapped quarter blocks with arc edges (body `Seal`, block faces `bottom` and
/// `top`), squeezed by [`SEAL_SQUEEZE`] (25 %) between plates bonded to it: `bottom` fixed, `top`
/// moved down with its in-plane motion held. Generic rubber: E = 6 MPa, ν = 0.49, Mooney–Rivlin
/// with C01 = 0.2 MPa (C10 = μ/2 − C01 = 0.81 MPa). Step `Squeeze` (nonlinear static).
pub fn rubber_seal() -> Value {
    let [inner, outer, height] = SEAL;
    let mut steps = vec![json!({"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}})];
    for quarter in 0..4u8 {
        let (a0, a1) = (f64::from(quarter) * std::f64::consts::FRAC_PI_2, f64::from(quarter + 1) * std::f64::consts::FRAC_PI_2);
        let at = |r: f64, a: f64, z: f64| [r * a.cos(), r * a.sin(), z];
        let corners = [
            at(inner, a0, 0.0),
            at(outer, a0, 0.0),
            at(outer, a1, 0.0),
            at(inner, a1, 0.0),
            at(inner, a0, height),
            at(outer, a0, height),
            at(outer, a1, height),
            at(inner, a1, height),
        ];
        let arc = |a: u8, b: u8, z: f64| json!({"corners": [a, b], "arc": {"center": [0, 0, z]}});
        steps.push(json!({"command": "mesh.block", "params": {
            "body": "Seal", "kind": "hex8", "corners": corners, "divisions": [5, 10, 5],
            "edges": [arc(1, 4, 0.0), arc(2, 3, 0.0), arc(5, 8, height), arc(6, 7, height)],
            "faces": {"1": "bottom", "2": "top"}
        }}));
    }
    steps.extend([
        json!({"command": "model.set.create", "params": {"name": "seal", "geometry": [{"body": "Seal"}]}}),
        json!({"command": "model.set.create", "params": {"name": "bottom", "geometry": [{"body": "Seal", "face": "bottom"}]}}),
        json!({"command": "model.set.create", "params": {"name": "top", "geometry": [{"body": "Seal", "face": "top"}]}}),
        json!({"command": "model.material.create", "params": {"name": "Rubber", "elastic": {"type": "isotropic", "young": "6 MPa", "poisson": 0.49}, "density": "1100 kg/m^3",
            "hyperelastic": {"type": "mooney_rivlin", "c01": 0.2}}}),
        json!({"command": "model.section.create", "params": {"name": "Seal", "set": "seal", "material": "Rubber", "kind": {"type": "solid"}, "formulation": {"hybrid": true}}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "Base", "type": "fixed", "target": {"set": "bottom"}}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "Press", "type": "displacement", "target": {"set": "top"}, "values": [0, 0, -SEAL_SQUEEZE, null, null, null]}}),
        json!({"command": "model.step.create", "params": {"name": "Squeeze", "procedure": {"type": "nonlinear_static", "initial_increment": 0.1}, "boundary_conditions": ["Base", "Press"]}}),
    ]);
    Value::Array(steps)
}

/// The flow sample's free-stream speed, mm/s.
pub const FLOW_SPEED: f64 = 25.0;
/// Air at 20 °C: density and dynamic viscosity.
const AIR: [Typed; 2] = [Typed::new(1.2, "kg/m^3"), Typed::new(1.8e-5, "Pa s")];

/// Air at [`FLOW_SPEED`] past [`UPRIGHT_BRACKET`] standing on the floor, the upright's back face
/// to the wind: a 240 × 120 × 100 mm box of 2.5 mm cells (184 320), a uniform inlet at x_min, an
/// outlet at zero pressure at x_max, a no-slip floor and slip walls at the sides and top. Reynolds
/// number ρ U h / μ = 83 on the upright's 50 mm height: laminar, a steady wake. Step `Flow` (16 s,
/// about 6.7 flow-throughs; ten frames) reports the force on each named face of the bracket.
pub fn flow_past() -> Value {
    let [density, viscosity] = AIR;
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "model.material.create", "params": {"name": "Air", "density": density, "fluid": {"type": "newtonian", "viscosity": viscosity}}},
        {"command": "model.fluid_domain.create", "params": {
            "name": "Wind", "material": "Air", "min": [-60, -40, 0], "max": [180, 80, 100], "spacing": 2.5,
            "faces": {
                "x_min": {"type": "velocity", "velocity": [FLOW_SPEED, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0},
                "y_min": {"type": "slip"}, "y_max": {"type": "slip"}, "z_min": {"type": "wall"}, "z_max": {"type": "slip"}
            },
            "obstacles": ["Bracket"], "initial_velocity": [FLOW_SPEED, 0, 0]
        }},
        {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Wind", "end_time": 16, "frames": 10}}},
    ])
}

/// The corpus filleted bracket (`tests/corpus/filleted-bracket`): an 8 mm L-bracket, foot along
/// x (60 mm) and upright along y (50 mm), extruded 40 mm along z, an R6 fillet, an 8 mm hole
/// through each leg (`F4` through the foot at x = 40, `F5` through the upright at y = 32).
pub const FILLETED_BRACKET: &str = include_str!("../../../tests/corpus/filleted-bracket/recipe.json");

/// The heated bracket: the heater's and the surrounding air's temperatures (°C, also the stress
/// free temperature), and steel's expansion coefficient (1/K).
pub const HEATER: f64 = 120.0;
pub const AMBIENT: f64 = 20.0;
pub const STEEL_EXPANSION: f64 = 1.2e-5;

/// [`FILLETED_BRACKET`] in steel (k = 50 W/(m K), α = 12 × 10⁻⁶/K, stress free at 20 °C), tet10
/// elements of 4 mm: its upright's end (`F2:side:l5`, set `hot`) held at [`HEATER`], the rest of
/// its outside (set `skin`) cooled by still air at [`AMBIENT`] (h = 10 W/(m² K)), both hole bores
/// clamped (set `holes`). Step `Heat` (coupled temperature–displacement, steady): the
/// temperature field and the stress of the expansion the clamps hold back.
pub fn heated_bracket() -> Value {
    let face = |f: &str| json!({"body": "Bracket", "face": f});
    let skin: Vec<Value> =
        ["F2:start", "F2:top", "F2:side:l1", "F2:side:l2", "F2:side:l3", "F2:side:l4", "F2:side:l6", "F3:blend:F2:side:l3|F2:side:l4"]
            .into_iter()
            .map(face)
            .collect();
    let holes: Vec<Value> = ["F4:bore:side:0", "F4:bore:side:1", "F5:bore:side:0", "F5:bore:side:1"].into_iter().map(face).collect();
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "mesh.tet.generate", "params": {"body": "Bracket", "size": 4, "order": "quadratic", "tolerance": 0.05}},
        {"command": "model.set.create", "params": {"name": "bracket", "geometry": [{"body": "Bracket"}]}},
        {"command": "model.set.create", "params": {"name": "holes", "geometry": holes}},
        {"command": "model.set.create", "params": {"name": "hot", "geometry": [face("F2:side:l5")]}},
        {"command": "model.set.create", "params": {"name": "skin", "geometry": skin}},
        {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": STEEL_YOUNG, "poisson": 0.3}, "density": STEEL_DENSITY,
            "thermal": {"conductivity": "50 W/(m K)", "specific_heat": "460 J/(kg K)", "expansion": STEEL_EXPANSION, "reference_temperature": AMBIENT}}},
        {"command": "model.section.create", "params": {"name": "Bracket", "set": "bracket", "material": "Steel", "kind": {"type": "solid"}}},
        {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "holes"}}},
        {"command": "model.boundary_condition.create", "params": {"name": "Heater", "type": "temperature", "target": {"set": "hot"}, "value": HEATER}},
        {"command": "model.load.create", "params": {"name": "Air", "type": "convection", "set": "skin", "film": "10 W/(m^2 K)", "sink": AMBIENT}},
        {"command": "model.step.create", "params": {"name": "Heat", "procedure": {"type": "coupled_temperature_displacement"}, "boundary_conditions": ["Clamp", "Heater"], "loads": ["Air"]}},
    ])
}

/// The corpus filleted bracket split into 30 cells (`tests/corpus/filleted-bracket-blocks`),
/// each a mapped block: the part `mesh.hex.blocks` meshes all-hex.
pub const PARTITIONED_BRACKET: &str = include_str!("../../../tests/corpus/filleted-bracket-blocks/recipe.json");

/// [`PARTITIONED_BRACKET`] meshed all-hex as one conforming mesh of hex20 elements of 3 mm, in
/// steel: the upright's end (`F2:side:l5` in every cell that has it, set `top`) clamped, 1 kN down
/// (−y) on the foot's end (`F2:side:l2`, set `tip`). Step `Load` (linear static). `bodies` are the
/// design's cells.
pub fn bracket_hex(bodies: &[Shape]) -> Value {
    let faces = |name: &str| -> Vec<Value> {
        bodies.iter().filter(|b| b.faces().iter().any(|f| f.name == name)).map(|b| json!({"body": b.name, "face": name})).collect()
    };
    let cells: Vec<Value> = bodies.iter().map(|b| json!({"body": b.name})).collect();
    json!([
        {"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}},
        {"command": "mesh.hex.blocks", "params": {"size": 3, "order": 2}},
        {"command": "model.set.create", "params": {"name": "bracket", "geometry": cells}},
        {"command": "model.set.create", "params": {"name": "top", "geometry": faces("F2:side:l5")}},
        {"command": "model.set.create", "params": {"name": "tip", "geometry": faces("F2:side:l2")}},
        {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": STEEL_YOUNG, "poisson": 0.3}, "density": STEEL_DENSITY}},
        {"command": "model.section.create", "params": {"name": "Bracket", "set": "bracket", "material": "Steel", "kind": {"type": "solid"}}},
        {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "top"}}},
        {"command": "model.load.create", "params": {"name": "Push", "type": "force", "target": {"set": "tip"}, "vector": [0, -1000, 0]}},
        {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Push"]}},
    ])
}

/// The bicycle crank arm (`crates/engine/samples/crank-arm.json`): a 170 mm aluminium arm, 16 mm
/// thick and tapering from 38 to 24 mm, joined by tangent root arcs to a Ø46 × 22 spindle boss
/// with a 12.7 mm square bore (R1.5 corners) and a Ø30 × 20 pedal boss with a Ø14 hole, and an
/// 82 × 14 × 8 pocket in its back (R5 corners). Made to be swept all-hex by `mesh.hex.generate`
/// as a stack of extrusions along z: no pocket floor fillet or boss chamfers, and the root arcs
/// meet the bosses about 20° off tangent (tangent to a circle 0.5 mm inside each boss), since a
/// tangent meeting leaves a zero-angle cusp on the arm's top face.
pub const CRANK_ARM: &str = include_str!("../samples/crank-arm.json");

/// The crank arm's pedal load (N): a standing rider's weight on one pedal.
pub const PEDAL_LOAD: f64 = 1500.0;

/// The half of the crank [`crank_prepare`] cuts off on the symmetry plane (y < 0), the one
/// meshed; its mirror image is `CrankArm/b`.
pub const CRANK_HALF: &str = "CrankArm/a";

/// [`CRANK_ARM`] meshed structured all-hex (hex20 of 3 mm) from its half and mirrored, in
/// 6061-T6 aluminium: the square bore's faces and corner fillets (set `spindle`) fixed,
/// [`PEDAL_LOAD`] on the pedal hole (set `pedal`), in the crank's plane and across the arm (−y).
/// Step `Pedal` (linear static). The faces are picked by their surfaces, not their names: the
/// square bore's are the vertical planes 6.35 mm from the spindle axis and the R1.5 cylinders
/// beside it, the pedal hole is the Ø14 cylinder on the pedal axis.
///
/// The mesh: [`crank_prepare`] cuts the crank on its symmetry plane only; `mesh.hex.generate`
/// extrudes the half [`CRANK_HALF`], decomposing its cross-section into blocks itself (half
/// rings round both bores, the pocket's corners as triangle primitives, four-sided blocks along
/// the arm), with no seeds; `mesh.mirror` reflects it onto `CrankArm/b`, joined on the plane.
pub fn crank_arm(bodies: &[Shape]) -> Value {
    let mesh = vec![
        json!({"command": "mesh.hex.generate", "params": {"body": CRANK_HALF, "method": "extrude", "size": 3, "order": 2}}),
        json!({"command": "mesh.mirror", "params": {"plane": "Sym", "bodies": [CRANK_HALF], "as": {CRANK_HALF: "CrankArm/b"}}}),
    ];
    crank_model(bodies, mesh)
}

/// The crank arm's symmetry plane (`Sym`, y = 0) and the partition by it (`Half`): the one datum
/// the structured mesh needs.
fn crank_prepare() -> Value {
    json!([
        {"command": "geometry.datum", "params": {"name": "Sym", "by": {"type": "offset", "plane": "zx", "offset": 0}}},
        {"command": "geometry.partition", "params": {"name": "Half", "datum": "Sym"}},
    ])
}

/// The crank's model after the `mesh` commands: material, the spindle fixed, the rider on the
/// pedal, step `Pedal`.
fn crank_model(bodies: &[Shape], mesh: Vec<Value>) -> Value {
    let [spindle, pedal] = crank_faces(bodies);
    let mut script = vec![json!({"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}})];
    script.extend(mesh);
    script.extend(json_array(json!([
        {"command": "model.set.create", "params": {"name": "crank", "geometry": bodies.iter().map(|b| json!({"body": b.name})).collect::<Vec<_>>()}},
        {"command": "model.set.create", "params": {"name": "spindle", "geometry": spindle}},
        {"command": "model.set.create", "params": {"name": "pedal", "geometry": pedal}},
        {"command": "model.material.create", "params": {"name": "Al6061-T6", "elastic": {"type": "isotropic", "young": "68.9 GPa", "poisson": 0.33}, "density": ALUMINIUM_DENSITY}},
        {"command": "model.section.create", "params": {"name": "Crank", "set": "crank", "material": "Al6061-T6", "kind": {"type": "solid"}}},
        {"command": "model.boundary_condition.create", "params": {"name": "Spindle", "type": "fixed", "target": {"set": "spindle"}}},
        {"command": "model.load.create", "params": {"name": "Rider", "type": "force", "target": {"set": "pedal"}, "vector": [0, -PEDAL_LOAD, 0]}},
        {"command": "model.step.create", "params": {"name": "Pedal", "procedure": {"type": "static"}, "boundary_conditions": ["Spindle"], "loads": ["Rider"]}},
    ])));
    Value::Array(script)
}

fn json_array(v: Value) -> Vec<Value> {
    match v {
        Value::Array(a) => a,
        other => vec![other],
    }
}

/// The crank's square-bore faces and pedal-hole faces as set members, found by their surfaces.
fn crank_faces(bodies: &[Shape]) -> [Vec<Value>; 2] {
    use fieldcraft_geometry::FaceSurface;
    // Directions are unit vectors; the geometry is in millimetres.
    const SAME_DIRECTION: f64 = 1e-6; // unit-factors: a component of a unit vector
    let on = Typed::new(1.0, "nm").express(Quantity::Length, crate::units::geometry()).unwrap_or(f64::NAN);
    let vertical = |axis: [f64; 3]| axis[2].abs() > 1.0 - SAME_DIRECTION;
    let (mut spindle, mut pedal) = (Vec::new(), Vec::new());
    for body in bodies {
        let surfaces = body.face_surfaces().unwrap_or_default();
        for (face, surface) in body.faces().iter().zip(surfaces) {
            let member = json!({"body": body.name, "face": face.name});
            match surface {
                Some(FaceSurface::Plane { normal, point })
                    if normal[2].abs() < SAME_DIRECTION && ((normal[0] * point[0] + normal[1] * point[1]).abs() - 6.35).abs() < on =>
                {
                    spindle.push(member)
                }
                Some(FaceSurface::Cylinder { origin, axis, radius })
                    if vertical(axis) && (radius - 1.5).abs() < on && origin[0].hypot(origin[1]) < 8.0 =>
                {
                    spindle.push(member)
                }
                Some(FaceSurface::Cylinder { origin, axis, radius })
                    if vertical(axis) && (radius - 7.0).abs() < on && (origin[0] - 170.0).hypot(origin[1]) < on =>
                {
                    pedal.push(member)
                }
                _ => {}
            }
        }
    }
    [spindle, pedal]
}

/// Six times the signed volume of the tetrahedron `p[0..4]` (positive when p3 lies on the side
/// p0 → p1 → p2 turns counter-clockwise around).
fn signed_volume(p: &[[f64; 3]]) -> f64 {
    let (Some(a), Some(b), Some(c), Some(d)) = (p.first(), p.get(1), p.get(2), p.get(3)) else { return 0.0 };
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let w = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Session;

    /// The `file.sample` form offers every sample, by name.
    #[test]
    fn the_form_lists_every_sample() {
        let Some(crate::schema::Kind::Choice(options)) = crate::cmd::forms::SAMPLE.first().map(|f| f.kind) else { panic!("not a choice") };
        let mut offered: Vec<&str> = options.iter().map(|(n, _)| *n).collect();
        let mut names: Vec<&str> = SAMPLES.iter().map(|s| s.name).collect();
        offered.sort_unstable();
        names.sort_unstable();
        assert_eq!(offered, names);
    }

    #[test]
    fn the_sample_builds_positive_elements_that_fill_both_parts() {
        let mut s = Session::default();
        s.run_script(&script()).unwrap();
        let info = s.execute("model.info", &json!({})).unwrap();
        assert_eq!(info["elements"]["hex8"], 64);
        assert_eq!(info["elements"]["tet4"], 162);
        assert_eq!(info["nodes"], 9 * 5 * 3 + 4 * 4 * 4);
        // Tet volumes are positive and add up to the cube's.
        let mut volume = 0.0;
        for e in s.model.elements().values().filter(|e| e.kind == fieldcraft_model::ElementKind::Tet4) {
            // In millimetres, as the sample is drawn.
            let p: Vec<[f64; 3]> = e.nodes.iter().map(|n| crate::units::to_geometry_point(s.model.node(*n).unwrap().xyz)).collect();
            let v = signed_volume(&p) / 6.0;
            assert!(v > 0.0);
            volume += v;
        }
        assert!((volume - 15.0f64.powi(3)).abs() < 1e-9, "{volume}");
    }

    #[test]
    fn the_cantilever_solves_close_to_beam_theory() {
        let mut s = Session::default();
        s.run_script(&cantilever()).unwrap();
        let r = s.execute("solve.static", &json!({})).unwrap();
        let [l, w, h] = CANTILEVER_SIZE;
        let beam = 4.0 * 1000.0 * l.powi(3) / (210000.0 * w * h.powi(3));
        let tip = r["max_displacement"]["value"].as_f64().unwrap();
        // Shear deformation adds a little (L/h = 5).
        assert!(tip > beam && tip < 1.1 * beam, "{tip} vs {beam}");
        let f = s.execute("solve.frequency", &json!({})).unwrap();
        assert_eq!(f["modes"].as_array().map(Vec::len), Some(6));
    }

    /// The I-beam opens swept all-hex and ready; solved, its tip deflection is the beam theory's
    /// bending plus shear (P L³ / 3 E I + P L / G A_web: 0.149 + 0.035 mm) within 5 %.
    #[test]
    fn the_ibeam_bends_as_beam_theory() {
        let mut s = open("ibeam");
        assert!(s.model.elements().values().all(|e| e.kind == fieldcraft_model::ElementKind::Hex8));
        s.execute("solve.static", &json!({})).unwrap();
        let p = s.execute("solve.probe", &json!({"step": "Load", "at": [[200.0, 0.0, 20.0]]})).unwrap();
        let tip = -p["nodes"][0]["displacement"][2].as_f64().unwrap();
        let (e, g) = (210000.0, 210000.0 / 2.6);
        let i = (40.0 * 40.0f64.powi(3) - 38.0 * 34.0f64.powi(3)) / 12.0;
        let load = sample_value(Quantity::Force, IBEAM_LOAD);
        let theory = load * 200.0f64.powi(3) / (3.0 * e * i) + load * 200.0 / (g * 2.0 * 37.0);
        assert!((tip - theory).abs() < 0.05 * theory, "{tip} against {theory}");
    }

    #[test]
    fn the_impact_runs_and_keeps_frames_with_stresses() {
        let mut s = Session::default();
        s.run_script(&impact()).unwrap();
        let r = s.execute("solve.explicit", &json!({})).unwrap();
        assert_eq!((r["completed"].as_bool(), r["frames"].as_u64()), (Some(true), Some(11)), "{r}");
        assert!(r["energy_balance"].as_f64().unwrap() < 0.05, "{r}");
    }

    /// Hex20 extrudes the crank with its root arcs tangent to the spindle boss (#197), at the
    /// sizes where it once failed the validity gate with inverted elements and hex8 did not: the
    /// sample moves the root arcs off tangent, so this builds the tangent variant (root arcs of
    /// R10 whose centres lie 33 mm from the spindle axis, R23 + R10). The elements at the cusp,
    /// whose two curved sides meet tangentially, stay valid.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn hex20_extrudes_the_crank_with_tangent_root_arcs() {
        let mut script = CRANK_ARM.to_string();
        for (from, to) in [
            ("15.883151316213851, 28.354461805265395", "16.9566899202338, 28.310257274511628"),
            ("15.883151316213851, -28.354461805265395", "16.9566899202338, -28.310257274511628"),
            ("15.471735240951828, 18.362928548901984", "16.545273844971774, 18.31872401814822"),
            ("15.471735240951828, -18.362928548901984", "16.545273844971774, -18.31872401814822"),
            ("13.455063072993756, 18.65372021077189", "11.818299035314467, 19.731391433750527"),
            ("13.455063072993756, -18.65372021077189", "11.818299035314467, -19.731391433750527"),
        ] {
            assert!(script.contains(from), "{from}");
            script = script.replace(from, to);
        }
        let bytes = fieldcraft_geometry::script_design(&script).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        for size in [2.5, 3.0, 5.0] {
            let r = s.execute("mesh.hex.generate", &json!({"body": "CrankArm", "size": size, "order": 2, "method": "extrude"})).unwrap();
            let body = &r["bodies"][0];
            assert_eq!(body["kind"], "hex20", "{r}");
            assert!(body["quality"]["min_scaled_jacobian"].as_f64().unwrap() > 0.0, "{size}: {r}");
        }
    }

    /// The crank's half y > 0 divides at x = 110 … 148 and by the diagonal plane through
    /// (125, 2) as its mirror half does (#283): the oversized cutting boxes once made the kernel's
    /// boolean fail there, reported as a plane that does not divide. Each cut leaves two cells
    /// whose volumes add up to the half's (within the measure's tessellation).
    #[test]
    fn the_crank_half_divides_where_its_mirror_does() {
        let mut s = open("crank-mirror");
        let volume = |s: &mut Session, body: &str| -> f64 {
            let r = s.execute("geometry.bodies", &json!({})).unwrap();
            let b = r.as_array().unwrap().iter().find(|b| b["name"] == body).unwrap_or_else(|| panic!("{body}: {r}"));
            b["volume"].as_f64().unwrap_or_else(|| panic!("{b}"))
        };
        let half = volume(&mut s, "CrankArm/b");
        let (sin, cos) = (32f64.to_radians().sin(), 32f64.to_radians().cos());
        for (name, by) in [
            ("X110", json!({"type": "offset", "plane": "yz", "offset": 110})),
            ("X148", json!({"type": "offset", "plane": "yz", "offset": 148})),
            ("X148R", json!({"type": "point_normal", "point": [148, 0, 0], "normal": [-1, 0, 0]})),
            ("Diag", json!({"type": "point_normal", "point": [125, 2, 0], "normal": [sin, -cos, 0]})),
        ] {
            s.execute("geometry.datum", &json!({"name": name, "by": by})).unwrap();
            let r = s.execute("geometry.partition", &json!({"name": name, "datum": name, "bodies": ["CrankArm/b"]})).unwrap();
            assert_eq!(r["cells"], json!(["CrankArm/b/a", "CrankArm/b/b"]), "{name}: {r}");
            let cells = volume(&mut s, "CrankArm/b/a") + volume(&mut s, "CrankArm/b/b");
            assert!((cells - half).abs() < 1e-5 * half, "{name}: {cells} against {half}");
            s.execute("geometry.partition.delete", &json!({"name": name})).unwrap();
        }
    }

    /// A cut at x = 130 runs along the pocket's end wall (#427), and the kernel only makes it
    /// with its cutting box nudged 1.3 µm and 2 µrad off the plane: on the half y > 0 that left a
    /// cell bounded by slivers next to the wall, which `mesh.cells` refused with an open boundary
    /// at x = 130. The partition refuses such a cut, saying why, on the whole crank and on the
    /// half, and leaves the body whole; the plane moved past the pocket (x = 136) cuts the half
    /// into cells that mesh.
    #[test]
    fn a_cut_the_kernel_can_only_nudge_off_its_plane_is_refused() {
        let datum = |name: &str, x: f64| json!({"name": name, "by": {"type": "offset", "plane": "yz", "offset": x}});
        let refused = |s: &mut Session, body: &str| {
            let e = s.execute("geometry.partition", &json!({"name": "X130", "datum": "X130", "bodies": [body]})).unwrap_err().to_string();
            assert!(e.contains(&format!("could not cut `{body}`")) && e.contains("only made off the plane"), "{e}");
        };
        let names = |s: &mut Session| -> Vec<String> {
            let r = s.execute("geometry.bodies", &json!({})).unwrap();
            r.as_array().unwrap().iter().filter_map(|b| b["name"].as_str().map(String::from)).collect()
        };
        let bytes = fieldcraft_geometry::script_design(CRANK_ARM).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut whole = Session::default();
        whole.execute("geometry.open", &json!({"path": path})).unwrap();
        whole.execute("geometry.datum", &datum("X130", 130.0)).unwrap();
        refused(&mut whole, "CrankArm");
        assert_eq!(names(&mut whole), ["CrankArm"]);

        let mut s = open("crank-mirror");
        s.execute("geometry.datum", &datum("X130", 130.0)).unwrap();
        refused(&mut s, "CrankArm/b");
        assert_eq!(names(&mut s), ["CrankArm/a", "CrankArm/b"]);
        s.execute("geometry.datum", &datum("X136", 136.0)).unwrap();
        let r = s.execute("geometry.partition", &json!({"name": "X136", "datum": "X136", "bodies": ["CrankArm/b"]})).unwrap();
        assert_eq!(r["cells"], json!(["CrankArm/b/a", "CrankArm/b/b"]), "{r}");
        s.execute("mesh.cells", &json!({"bodies": ["CrankArm/b/a", "CrankArm/b/b"], "default": "extrude", "size": 3})).unwrap();
    }

    /// The whole crank extruded all-hex (#198): its back pocket's floor is a rounded rectangle,
    /// straight sides running tangent into corner arcs, so the region is a disk whose O-grid's
    /// middle block follows its 6:1 elongation; a square middle block squashed its elements to a
    /// scaled Jacobian of 0.16 (now 0.49 at 2 mm, the worst elements elsewhere at 3 and 5 mm).
    #[test]
    fn the_crank_pocket_floor_keeps_its_elements_square() {
        let bytes = fieldcraft_geometry::script_design(CRANK_ARM).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        for size in [2.0, 3.0, 5.0] {
            let r = s.execute("mesh.hex.generate", &json!({"body": "CrankArm", "size": size, "method": "extrude"})).unwrap();
            let min = r["bodies"][0]["quality"]["min_scaled_jacobian"].as_f64().unwrap();
            assert!(min > 0.35, "{size}: {r}");
        }
    }

    /// The one-to-one sweep and the crank (#444). The arm's profile extruded alone is a
    /// one-to-one sweep (one source, one target) and sweeps all-hex as well as extrude does. The
    /// crank itself ends at several levels (pocket floor at 8 mm, arm 16, pedal boss 20, spindle
    /// boss 22), so no one-to-one sweep takes it: the sweep says so, naming a cap only part of a
    /// ring reaches, where it reported a "linking face with no top side"; extrude meshes it.
    #[test]
    fn the_sweep_takes_the_arm_profile_and_names_the_crank_levels() {
        let session = |script: &str| {
            let bytes = fieldcraft_geometry::script_design(script).unwrap();
            let (dir, path) = crate::tests::design_file(&bytes);
            let mut s = Session::default();
            s.execute("geometry.open", &json!({"path": path})).unwrap();
            (dir, s)
        };
        let quality = |s: &mut Session, method: &str| -> crate::Result<f64> {
            let r = s.execute("mesh.hex.generate", &json!({"body": "CrankArm", "size": 3, "method": method}))?;
            assert_eq!(r["bodies"][0]["method"], method, "{r}");
            Ok(r["bodies"][0]["quality"]["min_scaled_jacobian"].as_f64().unwrap())
        };
        // The script up to the arm's extrusion: the profile alone.
        let end = CRANK_ARM.find("{\"command\": \"solid.cylinder\"").unwrap();
        let profile = format!("{}]}}", CRANK_ARM[..end].trim_end().trim_end_matches(','));
        let (_dir, mut s) = session(&profile);
        let (sweep, extrude) = (quality(&mut s, "sweep").unwrap(), quality(&mut s, "extrude").unwrap());
        assert!(sweep > 0.3 && sweep > 0.95 * extrude, "sweep {sweep}, extrude {extrude}");

        let (_dir, mut s) = session(CRANK_ARM);
        let e = quality(&mut s, "sweep").unwrap_err().to_string();
        assert!(e.contains("is parallel to the source and only part of the ring reaches it") && e.contains("method `extrude`"), "{e}");
        assert!(quality(&mut s, "extrude").unwrap() > 0.35);
    }

    /// `file.sample NAME`, checked clean: ready, no errors or warnings.
    fn open(name: &str) -> Session {
        let mut s = Session::default();
        s.execute("file.sample", &json!({"name": name})).unwrap();
        let check = s.execute("model.check", &json!({})).unwrap();
        assert_eq!((check["ready"].as_bool(), check["errors"].as_u64(), check["warnings"].as_u64()), (Some(true), Some(0), Some(0)), "{check}");
        s
    }

    /// The squeeze force of the ring, against the closed form for an incompressible Mooney–Rivlin
    /// solid in uniaxial compression, P = 2 (λ − λ⁻²)(C10 + C01/λ) per unit of original area: the
    /// same ring between frictionless plates (each face held only along z, the rigid motions
    /// held at three nodes on the symmetry planes) compresses evenly and matches it within 3 %
    /// (ν = 0.49 is nearly incompressible); bonded, the bulge needs more force (shape factor
    /// ¼ (b − a)/h = 0.5 by the loaded over the free area), within 2.5 times.
    #[test]
    fn the_rubber_seal_squeezes_as_mooney_rivlin_rubber() {
        let [inner, outer, height] = SEAL;
        let area = std::f64::consts::PI * (outer * outer - inner * inner);
        let lambda: f64 = 1.0 - SEAL_SQUEEZE / height;
        let (c01, c10) = (0.2, 6.0 / (2.0 * 1.49) / 2.0 - 0.2);
        let closed = 2.0 * (lambda - lambda.powi(-2)) * (c10 + c01 / lambda) * area;
        let squeeze = |s: &mut Session| -> f64 {
            s.execute("solve.nonlinear", &json!({})).unwrap();
            let r = s.execute("results.reactions", &json!({"step": "Squeeze", "sets": ["top"]})).unwrap();
            r["sets"][0]["force"][2].as_f64().unwrap()
        };
        let mut s = open("rubber-seal");
        let bonded = squeeze(&mut s);
        let node = |s: &Session, xyz: [f64; 3]| -> u32 {
            let mm = crate::units::to_geometry_point;
            *s.model.nodes().iter().find(|(_, n)| (0..3).all(|k| (mm(n.xyz)[k] - xyz[k]).abs() < 1e-9)).unwrap().0
        };
        let held = [
            ("hold_y", [inner, 0.0, 0.0], [None, Some(0.0), None]),
            ("hold_x", [0.0, inner, 0.0], [Some(0.0), None, None]),
            ("hold_y2", [-inner, 0.0, 0.0], [None, Some(0.0), None]),
        ];
        for (name, at, [ux, uy, _]) in held {
            let n = node(&s, at);
            s.execute("model.set.create", &json!({"name": name, "nodes": [n]})).unwrap();
            s.execute(
                "model.boundary_condition.create",
                &json!({"name": name, "type": "displacement", "target": {"set": name}, "values": [ux, uy, null, null, null, null]}),
            )
            .unwrap();
        }
        s.execute(
            "model.boundary_condition.edit",
            &json!({"name": "Base", "type": "displacement", "target": {"set": "bottom"}, "values": [null, null, 0, null, null, null]}),
        )
        .unwrap();
        s.execute(
            "model.boundary_condition.edit",
            &json!({"name": "Press", "type": "displacement", "target": {"set": "top"}, "values": [null, null, -SEAL_SQUEEZE, null, null, null]}),
        )
        .unwrap();
        s.execute("model.step.edit", &json!({"name": "Squeeze", "boundary_conditions": ["Base", "Press", "hold_y", "hold_x", "hold_y2"]})).unwrap();
        let free = squeeze(&mut s);
        assert!((free - closed).abs() < 0.03 * closed.abs(), "frictionless {free} N against {closed} N");
        let stiffer = bonded / free;
        assert!(stiffer > 1.05 && stiffer < 2.5, "bonded {bonded} N, frictionless {free} N");
    }

    /// The seal's force–travel curve (#204): every frame keeps the reaction on `top`. Squeezed
    /// by a prescribed travel, the top's mean travel is the load factor times the squeeze and its
    /// reaction grows with it at every frame (a smooth, stiffening curve), the base reacting
    /// equal and opposite; results.reactions of a frame and results.xy read the same numbers.
    /// Pressed by a force instead, the base reacts with the load factor times that force at
    /// every frame.
    #[test]
    fn the_rubber_seal_keeps_a_force_travel_curve() {
        let mut s = open("rubber-seal");
        s.execute("mesh.node.add", &json!({"label": 99_999, "xyz": [0, 0, 50]})).unwrap();
        s.execute("model.set.create", &json!({"name": "loose", "nodes": [99_999]})).unwrap();
        s.execute("solve.nonlinear", &json!({})).unwrap();
        let h = s.execute("solve.history", &json!({"step": "Squeeze", "sets": ["top", "bottom"]})).unwrap();
        let floats = |v: &Value| -> Vec<f64> { v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect() };
        let lambda = floats(&h["load_factor"]);
        assert!(lambda.len() >= 6 && lambda[0] == 0.0 && lambda.last() == Some(&1.0), "{lambda:?}");
        let (top, bottom) = (&h["sets"][0], &h["sets"][1]);
        let force: Vec<f64> = top["reaction"].as_array().unwrap().iter().map(|f| f[2].as_f64().unwrap()).collect();
        let travel: Vec<f64> = top["displacement"].as_array().unwrap().iter().map(|u| u[2].as_f64().unwrap()).collect();
        let scale = force.last().unwrap().abs();
        for (k, l) in lambda.iter().enumerate() {
            assert!((travel[k] + l * SEAL_SQUEEZE).abs() < 1e-9, "frame {k}: travel {} at λ {l}", travel[k]);
            let base = bottom["reaction"][k][2].as_f64().unwrap();
            assert!((force[k] + base).abs() < 1e-6 * scale, "frame {k}: top {} base {base}", force[k]);
        }
        // The press pushes the top down (the reaction on it is negative), harder at every frame,
        // and the rubber stiffens in compression: the secant stiffness never falls.
        for k in 1..force.len() {
            assert!(force[k] < force[k - 1], "{force:?}");
            if k > 1 {
                let (a, b) = (force[k - 1] / travel[k - 1], force[k] / travel[k]);
                assert!(b > a * (1.0 - 1e-3), "secant stiffness fell from {a} to {b}: {force:?}");
            }
        }
        let last = s.execute("results.reactions", &json!({"step": "Squeeze", "sets": ["top"]})).unwrap();
        assert!((last["sets"][0]["force"][2].as_f64().unwrap() + scale).abs() < 1e-9 * scale, "{last}");
        let k = lambda.len() / 2;
        let at = s.execute("results.reactions", &json!({"step": "Squeeze", "sets": ["top"], "frame": k})).unwrap();
        assert_eq!((at["frame"].as_u64(), at["load_factor"].as_f64()), (Some(k as u64), Some(lambda[k])), "{at}");
        assert!((at["sets"][0]["force"][2].as_f64().unwrap() - force[k]).abs() < 1e-9 * scale, "{at}");
        let xy = s
            .execute(
                "results.xy",
                &json!({"step": "Squeeze", "y": {"reaction": "top", "component": "z"}, "x": {"displacement": "top", "component": "z"}}),
            )
            .unwrap();
        assert_eq!(xy["x"]["axis"], "displacement");
        assert_eq!(floats(&xy["x"]["values"]), travel);
        assert!(floats(&xy["y"]["values"]).iter().zip(&force).all(|(a, b)| (a - b).abs() < 1e-9 * scale), "{xy}");
        // Hostile parameters are errors; so is a set without result nodes.
        for (id, p) in [
            ("results.reactions", json!({"step": "Squeeze", "frame": 100000})),
            ("results.reactions", json!({"step": "Squeeze", "frame": -1})),
            ("results.reactions", json!({"step": "Squeeze", "sets": ["nope"], "frame": 1})),
            ("solve.history", json!({"step": "Squeeze", "sets": ["nope"]})),
            ("solve.history", json!({"step": "Squeeze", "node": 0})),
            ("solve.history", json!({"step": "Squeeze", "every": 0})),
            ("results.xy", json!({"step": "Squeeze", "y": {"reaction": "nope"}})),
            ("results.xy", json!({"step": "Squeeze", "y": {"reaction": "top", "component": "w"}})),
            ("results.xy", json!({"step": "Squeeze", "y": {"reaction": "top"}, "x": {"displacement": 3}})),
            ("results.xy", json!({"step": "Squeeze", "y": {"reaction": "top"}, "x": {}})),
        ] {
            assert!(s.execute(id, &p).is_err(), "{id} {p}");
        }
        // A set whose nodes have no results (a node outside the mesh) is an error too.
        let e = s.execute("solve.history", &json!({"step": "Squeeze", "sets": ["loose"]})).unwrap_err().to_string();
        assert!(e.contains("set `loose` has no nodes with results"), "{e}");
        // Load control: the top held laterally and pressed by a force.
        let press = -0.5 * scale;
        s.execute(
            "model.boundary_condition.edit",
            &json!({"name": "Press", "type": "displacement", "target": {"set": "top"}, "values": [0, 0, null, null, null, null]}),
        )
        .unwrap();
        s.execute("model.load.create", &json!({"name": "Push", "type": "force", "target": {"set": "top"}, "vector": [0, 0, press]})).unwrap();
        s.execute("model.step.edit", &json!({"name": "Squeeze", "loads": ["Push"]})).unwrap();
        s.execute("solve.nonlinear", &json!({})).unwrap();
        let h = s.execute("solve.history", &json!({"step": "Squeeze", "sets": ["bottom"]})).unwrap();
        let lambda = floats(&h["load_factor"]);
        assert!(lambda.len() >= 3, "{h}");
        for (k, l) in lambda.iter().enumerate() {
            let base = h["sets"][0]["reaction"][k][2].as_f64().unwrap();
            assert!((base + l * press).abs() < 1e-6 * press.abs(), "frame {k}: base {base} at λ {l}, applied {}", l * press);
        }
    }

    #[test]
    fn the_flow_sample_opens_ready_to_solve() {
        let s = open("flow-past");
        assert_eq!(s.model.fluid_domains().len(), 1);
    }

    /// The drag on the bracket, as a coefficient on the upright's frontal area (40 × 50 mm), is
    /// that of a bluff plate standing on a wall at a Reynolds number of 83, raised by the 17 %
    /// blockage of the box: between 1 and 3.5 (a flat plate broadside to the flow is about 1.2 at
    /// high Reynolds number and rises as it falls). The lattice conserves mass.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_flow_past_the_bracket_has_a_bluff_body_drag() {
        let mut s = open("flow-past");
        let r = s.execute("solve.run", &json!({})).unwrap();
        let r = &r["steps"][0]["summary"];
        let drag: f64 = r["forces"].as_array().unwrap().iter().filter(|f| f["body"] == "Bracket").map(|f| f["mean"][0].as_f64().unwrap()).sum();
        let density = sample_value(Quantity::Density, AIR[0]);
        let cd = drag / (0.5 * density * FLOW_SPEED * FLOW_SPEED * 40.0 * 50.0);
        assert!(cd > 1.0 && cd < 3.5, "Cd {cd}: {r}");
        assert!(r["mass"]["relative_change"].as_f64().unwrap().abs() < 0.01, "{r}");
    }

    /// The heated bracket's 4 mm tet mesh stays near its 23 000 elements. #303: refinement piled
    /// face vertices against a feature curve, and the improvement then grew the mesh to 29 257
    /// elements in about two minutes. The element count is the check: wall time depends on the
    /// machine and the build profile (about 5 s in a release build, far longer on CI).
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_heated_bracket_opens_with_its_usual_mesh() {
        let s = open("heated-bracket");
        let elements = s.model.elements().len();
        assert!((20_000..=26_000).contains(&elements), "{elements} elements");
    }

    /// The heated bracket's temperatures lie between the air's and the heater's, falling along
    /// the bracket away from the heater. Its thermal expansion: with the air at the heater's
    /// temperature the bracket warms evenly by ΔT, and held at the foot's hole alone (at x = 40)
    /// it grows freely, so a point moves α ΔT times its distance from the hole's axis along x:
    /// the upright's top by −0.0432 mm (within 3 %), the foot's end by +0.024 mm (within 6 %).
    /// The clamped bore cannot grow, which holds both points back a little; how much depends on
    /// the tet mesh, which differs between platforms (#115): 2.1 % at the top on the Linux CI.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_heated_bracket_warms_and_expands_as_steel() {
        let mut s = open("heated-bracket");
        let r = s.execute("solve.run", &json!({})).unwrap();
        let r = &r["steps"][0]["summary"];
        assert!((r["max_temperature"]["value"].as_f64().unwrap() - HEATER).abs() < 1e-9, "{r}");
        let probe = |s: &mut Session| s.execute("solve.probe", &json!({"step": "Heat", "at": [[4.0, 50.0, 20.0], [60.0, 4.0, 20.0]]})).unwrap();
        let p = probe(&mut s);
        let (upright, foot) = (p["nodes"][0]["temperature"].as_f64().unwrap(), p["nodes"][1]["temperature"].as_f64().unwrap());
        assert!(upright > foot && foot > AMBIENT && upright <= HEATER, "{upright} {foot}");
        let foot_hole = json!([{"body": "Bracket", "face": "F4:bore:side:0"}, {"body": "Bracket", "face": "F4:bore:side:1"}]);
        s.execute("model.set.create", &json!({"name": "foot_hole", "geometry": foot_hole})).unwrap();
        s.execute("model.boundary_condition.edit", &json!({"name": "Clamp", "type": "fixed", "target": {"set": "foot_hole"}})).unwrap();
        s.execute("model.load.edit", &json!({"name": "Air", "type": "convection", "set": "skin", "film": 0.01, "sink": HEATER})).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        let p = probe(&mut s);
        let grow = STEEL_EXPANSION * (HEATER - AMBIENT);
        // The bore itself cannot grow, which holds back the foot's end, 16 mm from its edge, a
        // little more than the upright's top, 32 mm away.
        for (i, tolerance) in [(0, 0.03), (1, 0.06)] {
            let node = &p["nodes"][i];
            let (ux, expected) = (node["displacement"][0].as_f64().unwrap(), grow * (node["xyz"][0].as_f64().unwrap() - 40.0));
            assert!((ux - expected).abs() < tolerance * expected.abs(), "{node}: {ux} against {expected}");
        }
    }

    /// The hex20 bracket's tip deflection agrees with a fine tet10 solution of the unpartitioned
    /// part (0.5093 mm, crates/testkit/tests/bracket_blocks.rs) within 1 %, its mesh conforms
    /// across all 30 cells, and the load reaches the clamp.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_hex_bracket_bends_like_the_fine_tet_solution() {
        let mut s = open("bracket-hex");
        let elements = s.model.elements().values().filter(|e| e.kind == fieldcraft_model::ElementKind::Hex20).count();
        assert_eq!(elements, s.model.elements().len());
        let r = s.execute("solve.run", &json!({})).unwrap();
        let r = &r["steps"][0]["summary"];
        let probe = s.execute("solve.probe", &json!({"step": "Load", "at": [[60.0, 4.0, 20.0]]})).unwrap();
        let tip = probe["nodes"][0]["displacement"][1].as_f64().unwrap();
        assert!((tip + 0.5093).abs() < 0.01 * 0.5093, "{tip}");
        let reaction = s.execute("results.reactions", &json!({"step": "Load", "sets": ["top"]})).unwrap();
        assert!((reaction["sets"][0]["force"][1].as_f64().unwrap() - 1000.0).abs() < 1e-6, "{reaction}");
        assert!(r["max_von_mises"]["value"].as_f64().unwrap() > 200.0, "{r}");
    }

    /// The crank opens as one conforming structured hex20 mesh: the half meshed in blocks
    /// without partitions or seeds, no ties, the half and its mirror
    /// image alike and joined on the symmetry plane, the spindle and pedal faces found on both
    /// halves (the image's faces carry the names of the half they are the image of), and the
    /// quality floor: scaled Jacobian at least 0.37, edge ratio at most 10.
    #[test]
    fn the_crank_opens_as_one_structured_mesh_of_two_alike_halves() {
        let (arm, mirror) = (find("crank-arm").unwrap(), find("crank-mirror").unwrap());
        let prepared = |x: &Sample| x.prepare.map(|f| f());
        assert!(
            arm.design == mirror.design && prepared(arm) == prepared(mirror) && (arm.script)(&[]) == (mirror.script)(&[]),
            "one model under two names"
        );
        {
            let name = "crank-arm";
            let mut s = open(name);
            let bodies: Vec<String> = s.design.as_ref().unwrap().bodies.iter().map(|b| b.name.clone()).collect();
            assert_eq!(bodies, [CRANK_HALF, "CrankArm/b"]);
            // One datum, the symmetry plane, hidden.
            assert!(s.model.datums().len() == 1 && s.model.datums().keys().all(|d| !s.shows(crate::ItemKind::Datum, d)), "datums shown");
            assert!(!s.is_dirty(), "a sample just opened has nothing to save");
            assert!(s.model.ties().is_empty(), "no ties");
            let elements = s.model.elements().values().filter(|e| e.kind == fieldcraft_model::ElementKind::Hex20).count();
            assert_eq!(elements, s.model.elements().len());
            let mut set = |name: &str, bodies: &[&str]| {
                s.execute("model.set.create", &json!({"name": name, "geometry": bodies.iter().map(|b| json!({"body": b})).collect::<Vec<_>>()}))
                    .unwrap();
                s.model.resolve_set(name).unwrap()
            };
            let (a, b) = (set("half_a", &[CRANK_HALF]), set("half_b", &["CrankArm/b"]));
            assert_eq!(a.elements.len(), b.elements.len());
            assert_eq!(a.elements.len() + b.elements.len(), elements);
            assert_eq!(a.faces.len(), b.faces.len(), "the halves' skins are alike");
            let shared = a.nodes.iter().filter(|n| b.nodes.binary_search(n).is_ok()).count();
            assert!(shared > 100, "{shared}");
            assert_eq!(a.nodes.len() + b.nodes.len() - shared, s.model.nodes().len());
            let q = s.execute("mesh.quality", &json!({})).unwrap();
            assert_eq!(q["valid"].as_u64().unwrap() as usize, elements, "{q}");
            let metric = |m: &str, k: &str| q["metrics"][m][k].as_f64().unwrap();
            assert!(metric("scaled_jacobian", "min") >= 0.37, "{q}");
            assert!(metric("edge_ratio", "max") <= 10.0, "{q}");
            let members = |set: &str| s.model.sets().get(set).map_or(Vec::new(), |x| x.geometry.iter().map(|g| g.body.clone()).collect::<Vec<_>>());
            for set in ["spindle", "pedal"] {
                let bodies = members(set);
                assert!(bodies.iter().any(|b| b == CRANK_HALF) && bodies.iter().any(|b| b == "CrankArm/b"), "{set}: {bodies:?}");
                assert!(!s.model.resolve_set(set).unwrap().faces.is_empty(), "{set}");
            }
            let faces = |set: &str| s.model.sets().get(set).map_or(0, |x| x.geometry.len());
            assert_eq!((faces("spindle"), faces("pedal")), (10, 4), "{name}");
        }
    }

    /// Solved, the crank's root bends as a beam: σxx on the arm's sides at x = 29 mm (35.6 mm
    /// wide, 16 thick, 141 mm from the pedal) within 3 % of M/Z = 1500 · 141 / (16 · 35.6² / 6) =
    /// 62.6 MPa, equal and opposite on the two mirror-image sides, and the spindle carries the
    /// whole load.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_crank_arm_root_bends_as_the_beam_formula_says() {
        let mut s = open("crank-arm");
        s.execute("solve.run", &json!({})).unwrap();
        let half_width: f64 = 19.0 - 7.0 * 29.0 / 170.0;
        let probe = s.execute("solve.probe", &json!({"step": "Pedal", "at": [[29.0, -half_width, 8.0], [29.0, half_width, 8.0]]})).unwrap();
        let sxx = |i: usize| probe["nodes"][i]["stress"][0].as_f64().unwrap();
        let hand = PEDAL_LOAD * (170.0 - 29.0) / (16.0 * (2.0 * half_width).powi(2) / 6.0);
        for (i, sign) in [(0, -1.0), (1, 1.0)] {
            assert!((sign * sxx(i) - hand).abs() < 0.03 * hand, "{} against {hand} MPa: {probe}", sxx(i));
        }
        assert!((sxx(0) + sxx(1)).abs() < 0.01 * hand, "{probe}");
        let reaction = s.execute("results.reactions", &json!({"step": "Pedal", "sets": ["spindle"]})).unwrap();
        assert!((reaction["sets"][0]["force"][1].as_f64().unwrap() - PEDAL_LOAD).abs() < 1e-6, "{reaction}");
    }

    /// The conformity of `s`'s mesh ([`fieldcraft_mesh::quality::conformity`]).
    fn conformity(s: &Session) -> fieldcraft_mesh::quality::Conformity {
        let nodes: Vec<(u32, [f64; 3])> = s.model.nodes().iter().map(|(l, n)| (*l, n.xyz)).collect();
        let elements = s.model.elements().iter().map(|(l, e)| (*l, e.kind, e.nodes.as_slice()));
        fieldcraft_mesh::quality::conformity(&nodes, elements, 1e-6).unwrap()
    }

    /// The crank arm's hex20 mesh has no crack (#260): its free faces close up across the cells
    /// and the mirror plane, no two nodes coincide and every element edge has one mid-side node.
    #[test]
    fn the_crank_arm_hex20_mesh_is_conforming() {
        let s = open("crank-arm");
        let c = conformity(&s);
        assert!(c.is_conforming(), "{:?} {:?} {:?}", c.open_edges.first(), c.coincident.first(), c.midside_mismatches.first());
    }

    #[test]
    fn the_drop_test_opens_ready_to_solve() {
        let s = open("drop-test");
        assert!(s.model.elements().len() > 1000);
    }

    /// The bracket's mass is ρ V (V: foot, upright and the fillet's corner, 36 − 9π mm² by 40 mm),
    /// it meets the floor at √(2 g h) with ½ m v² of kinetic energy and conserves energy; by the
    /// end the floor has taken half its downward momentum and turned kinetic energy into
    /// vibration.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn the_drop_test_lands_at_the_free_fall_speed_and_bounces() {
        let mut s = open("drop-test");
        let r = s.execute("solve.run", &json!({})).unwrap();
        let r = &r["steps"][0]["summary"];
        let volume = 60.0 * 8.0 * 40.0 + 8.0 * 42.0 * 40.0 + (36.0 - 9.0 * std::f64::consts::PI) * 40.0;
        let mass = sample_value(Quantity::Density, ALUMINIUM_DENSITY) * volume;
        let close = |a: f64, b: f64, rel: f64| (a - b).abs() <= rel * b.abs();
        assert!(close(r["total_mass"].as_f64().unwrap(), mass, 0.005), "{r}");
        assert!(r["completed"].as_bool() == Some(true) && r["energy_balance"].as_f64().unwrap() < 0.01, "{r}");
        let h = s.execute("solve.history", &json!({"step": "Drop"})).unwrap();
        let speed = (2.0 * sample_value(Quantity::Acceleration, GRAVITY) * DROP_HEIGHT).sqrt();
        let pz = |i: usize| h["momentum"][i][2].as_f64().unwrap();
        let kinetic = |i: usize| h["energies"][i]["kinetic"].as_f64().unwrap();
        assert!(close(pz(0), -mass * speed, 0.005), "{}", pz(0));
        assert!(close(kinetic(0), 0.5 * mass * speed * speed, 0.01), "{}", kinetic(0));
        let last = h["time"].as_array().unwrap().len() - 1;
        let kept = pz(last) / pz(0);
        assert!(kept > 0.4 && kept < 0.6 && kinetic(last) < kinetic(0), "{} {}", pz(last), kinetic(last));
    }
}
