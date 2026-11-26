//! `results.openfoam` of a compressible case (rhoCentralFoam): the pressure, temperature, density
//! and Mach number it reads back, the refusal of a pressure or temperature that is not positive,
//! and the marker that decides what is read.

use std::path::{Path, PathBuf};

use fieldcraft_model::units::Quantity;
use serde_json::{Value, json};

use super::MARKER;
use crate::Session;

/// The lattice spacing of the test, in the session's display units.
const SPACING: f64 = 5.0;

/// The ambient pressure (Pa), the gas's specific heat at constant pressure and its specific gas
/// constant (J/(kg·K)), and the initial temperature (K) the case is written with.
const AMBIENT: f64 = 101_325.0;
const CP: f64 = 1005.0;
const R_S: f64 = 287.05;
const T_INITIAL: f64 = 293.15;

/// The inlet side's and outlet side's absolute pressure (Pa) and temperature (K) and the speed (m/s)
/// of the run written by hand into the case.
const P_INLET: f64 = 2.0e5;
const P_OUTLET: f64 = 1.0e5;
const T_INLET: f64 = 300.0;
const T_OUTLET: f64 = 250.0;
const SPEED: f64 = 100.0;

/// A value of `q` given in SI, in the session's display units, as a command parameter takes it.
fn show(s: &Session, q: Quantity, si: f64) -> f64 {
    s.units().express(q, si).unwrap()
}

/// A value of `q` given in the session's display units, in SI.
fn si(s: &Session, q: Quantity, shown: f64) -> f64 {
    s.units().to_si(q, shown).unwrap()
}

fn close(got: f64, want: f64) -> bool {
    (got / want - 1.0).abs() < 1e-6
}

/// The plate's compressible case, written by `deck.openfoam`, with a run written by hand into it:
/// the velocity everywhere, pressure and temperature by the side of the box's middle.
struct Case {
    _guard: crate::tests::TempDir,
    s: Session,
    case: PathBuf,
    min: Vec<f64>,
    max: Vec<f64>,
}

/// The sample plate in a box of air: a pressure face of 2 bar gauge on the inlet side, the ambient
/// on the outlet side, the plate as the obstacle; written as a compressible case.
fn written_case() -> Case {
    let (guard, design) = crate::tests::sample_design_file(&fieldcraft_geometry::sample::Plate::default());
    let dir = Path::new(&design).parent().unwrap().to_path_buf();
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": design})).unwrap();
    let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
    let [lo, hi] = [0, 1].map(|k| bodies[0]["bbox"][k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>());
    let min: Vec<f64> = lo.iter().map(|v| (v / SPACING).floor() * SPACING - 20.0).collect();
    let max: Vec<f64> = hi.iter().map(|v| (v / SPACING).ceil() * SPACING + 20.0).collect();
    let script = json!([
        {"command": "model.material.create", "params": {"name": "Air", "density": show(&s, Quantity::Density, 1.2),
            "fluid": {"type": "newtonian", "viscosity": show(&s, Quantity::DynamicViscosity, 1.8e-5)},
            "thermal": {"conductivity": show(&s, Quantity::ThermalConductivity, 0.0257), "specific_heat": show(&s, Quantity::SpecificHeat, CP)}}},
        {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Air", "min": min, "max": max, "spacing": SPACING, "obstacles": ["Body1"],
            "faces": {"x_min": {"type": "pressure", "pressure": show(&s, Quantity::Pressure, P_INLET - AMBIENT)}, "x_max": {"type": "pressure", "pressure": 0}}}},
        {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10}}},
    ]);
    s.run_script(&script).unwrap();
    let case = dir.join("case");
    let gas = json!({
        "ambient": show(&s, Quantity::Pressure, AMBIENT),
        "temperature": show(&s, Quantity::Temperature, T_INITIAL),
        "gas_constant": show(&s, Quantity::SpecificHeat, R_S),
        "far_field": show(&s, Quantity::Length, 0.05),
    });
    s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "compressible": gas})).unwrap();

    // OpenFOAM cells at the centres of the lattice's cells (metres), the box's middle in metres.
    let n: Vec<usize> = (0..3).map(|a| ((max[a] - min[a]) / SPACING).round() as usize).collect();
    let mut centres: Vec<[f64; 3]> = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let c = [i, j, k];
                centres.push(std::array::from_fn(|a| si(&s, Quantity::Length, min[a] + (c[a] as f64 + 0.5) * SPACING)));
            }
        }
    }
    let middle = si(&s, Quantity::Length, (min[0] + max[0]) / 2.0);
    let point = |c: &[f64; 3]| format!("({} {} {})", c[0], c[1], c[2]);
    std::fs::create_dir_all(case.join("0")).unwrap();
    std::fs::write(case.join("0/C"), field_file("volVectorField", "C", "vector", &centres.iter().map(point).collect::<Vec<_>>())).unwrap();
    for time in ["0.25", "0.5"] {
        write_time(&case, time, &centres, middle);
    }
    Case { _guard: guard, s, case, min, max }
}

/// A run's fields at one time: the velocity everywhere, the pressure and temperature of the cells
/// on the inlet side of `middle` (metres) at the inlet values, the others at the outlet values.
fn write_time(case: &Path, time: &str, centres: &[[f64; 3]], middle: f64) {
    let dir = case.join(time);
    std::fs::create_dir_all(&dir).unwrap();
    let inlet = |c: &[f64; 3]| c[0] < middle;
    let u: Vec<String> = centres.iter().map(|_| format!("({SPEED} 0 0)")).collect();
    std::fs::write(dir.join("U"), field_file("volVectorField", "U", "vector", &u)).unwrap();
    let p: Vec<String> = centres.iter().map(|c| (if inlet(c) { P_INLET } else { P_OUTLET }).to_string()).collect();
    std::fs::write(dir.join("p"), field_file("volScalarField", "p", "scalar", &p)).unwrap();
    let t: Vec<String> = centres.iter().map(|c| (if inlet(c) { T_INLET } else { T_OUTLET }).to_string()).collect();
    std::fs::write(dir.join("T"), field_file("volScalarField", "T", "scalar", &t)).unwrap();
}

/// A field file in OpenFOAM's ASCII format, its `internalField` the items, one per line.
fn field_file(class: &str, object: &str, kind: &str, items: &[String]) -> String {
    format!(
        "FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {class};\n    object {object};\n}}\ninternalField nonuniform List<{kind}> {}\n(\n{}\n)\n;\nboundaryField\n{{\n}}\n",
        items.len(),
        items.join("\n")
    )
}

/// The probe at `x` (display units), on the lattice's middle line, in the last frame.
fn probe_at(s: &mut Session, min: &[f64], x: f64) -> Value {
    let point = [x, min[1] + SPACING / 2.0, min[2] + SPACING / 2.0];
    s.execute("solve.probe", &json!({"step": "Flow", "points": [point], "frame": 1})).unwrap()["points"][0].clone()
}

fn read(s: &mut Session, case: &Path) -> Result<Value, String> {
    s.execute("results.openfoam", &json!({"path": case.to_str().unwrap()})).map_err(|e| e.to_string())
}

#[test]
fn a_compressible_case_reads_back_pressure_temperature_and_density() {
    let mut c = written_case();
    let r = read(&mut c.s, &c.case).unwrap();
    assert_eq!(r["frames"], 2, "{r}");
    let gamma = CP / (CP - R_S);
    let mach = SPEED / (gamma * R_S * T_OUTLET).sqrt();
    assert!(close(r["max_mach"].as_f64().unwrap(), mach), "the largest Mach number is on the cooler outlet side: {r}");
    assert_eq!(r["mach_basis"], "physical U/c", "the Mach number of an OpenFOAM case is physical, not lattice: {r}");

    let (min, max) = (c.min.clone(), c.max.clone());
    let s = &mut c.s;
    let inlet = probe_at(s, &min, min[0] + SPACING / 2.0);
    let p = si(s, Quantity::Pressure, inlet["pressure"].as_f64().unwrap());
    assert!(close(p, P_INLET - AMBIENT), "gauge pressure on the inlet side: {p}");
    let t = si(s, Quantity::Temperature, inlet["temperature"].as_f64().unwrap());
    assert!(close(t, T_INLET), "temperature on the inlet side: {t}");
    let rho = si(s, Quantity::Density, inlet["density"].as_f64().unwrap());
    assert!(close(rho, P_INLET / (R_S * T_INLET)), "density p/(R T) on the inlet side: {rho}");

    let outlet = probe_at(s, &min, max[0] - SPACING / 2.0);
    let p = si(s, Quantity::Pressure, outlet["pressure"].as_f64().unwrap());
    assert!(close(p, P_OUTLET - AMBIENT), "a negative gauge pressure on the outlet side: {p}");
    let t = si(s, Quantity::Temperature, outlet["temperature"].as_f64().unwrap());
    assert!(close(t, T_OUTLET), "temperature on the outlet side: {t}");
    let rho = si(s, Quantity::Density, outlet["density"].as_f64().unwrap());
    assert!(close(rho, P_OUTLET / (R_S * T_OUTLET)), "density on the outlet side: {rho}");
}

#[test]
fn a_pressure_or_temperature_that_is_not_positive_is_refused() {
    let mut c = written_case();
    read(&mut c.s, &c.case).unwrap();
    let min = c.min.clone();
    let before = probe_at(&mut c.s, &min, min[0] + SPACING / 2.0);
    let time = c.case.join("0.5");
    let (p, t) = (std::fs::read_to_string(time.join("p")).unwrap(), std::fs::read_to_string(time.join("T")).unwrap());

    // The first cell's temperature is zero: the cell and the value are named, and no result is kept.
    let first_t = format!("\n{T_INLET}\n");
    std::fs::write(time.join("T"), t.replacen(&first_t, "\n0\n", 1)).unwrap();
    let e = read(&mut c.s, &c.case).unwrap_err();
    assert!(e.contains("temperature") && e.contains("cell"), "{e}");
    assert_eq!(probe_at(&mut c.s, &min, min[0] + SPACING / 2.0), before, "a refused read keeps no result");
    std::fs::write(time.join("T"), &t).unwrap();

    // The first cell's pressure is negative.
    let first_p = format!("\n{P_INLET}\n");
    std::fs::write(time.join("p"), p.replacen(&first_p, "\n-1\n", 1)).unwrap();
    let e = read(&mut c.s, &c.case).unwrap_err();
    assert!(e.contains("pressure") && e.contains("cell"), "{e}");
    assert_eq!(probe_at(&mut c.s, &min, min[0] + SPACING / 2.0), before);
    std::fs::write(time.join("p"), &p).unwrap();

    // A time without its temperature is refused.
    std::fs::remove_file(time.join("T")).unwrap();
    assert!(read(&mut c.s, &c.case).is_err(), "a time directory without T");
    assert_eq!(probe_at(&mut c.s, &min, min[0] + SPACING / 2.0), before);
}

#[test]
fn a_density_that_overflows_single_precision_is_refused() {
    let mut c = written_case();
    let time = c.case.join("0.5");
    let t = std::fs::read_to_string(time.join("T")).unwrap();
    // A temperature near zero is positive, and p/(R T) is infinite as an f32.
    std::fs::write(time.join("T"), t.replacen(&format!("\n{T_INLET}\n"), "\n1e-40\n", 1)).unwrap();
    let e = read(&mut c.s, &c.case).unwrap_err();
    assert!(e.contains("density") && e.contains("cell 0") && e.contains("0.5") && e.contains("T"), "{e}");
}

/// The case's cells split over two processors by the parity of their index, each processor with
/// its own centres and fields; the marker is copied.
fn decomposed(c: &Case) -> (PathBuf, [Vec<[f64; 3]>; 2], f64) {
    let s = &c.s;
    let n: Vec<usize> = (0..3).map(|a| ((c.max[a] - c.min[a]) / SPACING).round() as usize).collect();
    let mut cells: [Vec<[f64; 3]>; 2] = [Vec::new(), Vec::new()];
    let mut index = 0usize;
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let at = [i, j, k];
                cells[index % 2].push(std::array::from_fn(|a| si(s, Quantity::Length, c.min[a] + (at[a] as f64 + 0.5) * SPACING)));
                index += 1;
            }
        }
    }
    let middle = si(s, Quantity::Length, (c.min[0] + c.max[0]) / 2.0);
    let split = c.case.parent().unwrap().join("split");
    std::fs::create_dir_all(&split).unwrap();
    std::fs::copy(c.case.join(MARKER), split.join(MARKER)).unwrap();
    for (rank, centres) in cells.iter().enumerate() {
        let part = split.join(format!("processor{rank}"));
        std::fs::create_dir_all(part.join("0")).unwrap();
        let point = |c: &[f64; 3]| format!("({} {} {})", c[0], c[1], c[2]);
        std::fs::write(part.join("0/C"), field_file("volVectorField", "C", "vector", &centres.iter().map(point).collect::<Vec<_>>())).unwrap();
        for time in ["0.25", "0.5"] {
            write_time(&part, time, centres, middle);
        }
    }
    (split, cells, middle)
}

#[test]
fn a_bad_cell_of_a_decomposed_case_is_named_by_its_processor_file_and_local_index() {
    let mut c = written_case();
    let (split, _, _) = decomposed(&c);
    let r = read(&mut c.s, &split).unwrap();
    assert_eq!(r["source"], "decomposed (2 processors)", "{r}");

    // Processor 1's cell 1: the first value of the inlet side is cell 0 of it, the second is 1.
    let file = |rank: usize, name: &str| split.join(format!("processor{rank}/0.5/{name}"));
    let edit = |rank: usize, name: &str, local: usize, value: &str| {
        let text = std::fs::read_to_string(file(rank, name)).unwrap();
        let (head, body) = text.split_once("(\n").unwrap();
        let (values, tail) = body.split_once("\n)\n").unwrap();
        let mut lines: Vec<&str> = values.split('\n').collect();
        lines[local] = value;
        std::fs::write(file(rank, name), format!("{head}(\n{}\n)\n{tail}", lines.join("\n"))).unwrap();
        text
    };
    for (rank, name, value, what) in [(1, "T", "0", "temperature"), (0, "p", "-1", "pressure"), (1, "T", "1e-40", "density")] {
        let original = edit(rank, name, 1, value);
        let e = read(&mut c.s, &split).unwrap_err();
        assert!(e.contains(&format!("processor{rank}")) && e.contains(&format!("0.5/{name}")), "the processor's file: {e}");
        assert!(e.contains("cell 1") && e.contains(what), "the local index and the quantity: {e}");
        std::fs::write(file(rank, name), original).unwrap();
    }
    assert!(read(&mut c.s, &split).is_ok());
}

#[test]
fn the_marker_decides_what_is_read() {
    let mut c = written_case();
    let original = std::fs::read_to_string(c.case.join(MARKER)).unwrap();
    let rewrite = |case: &Path, solver: &str, keep_gas: bool| {
        let mut m: Value = serde_json::from_str(&original).unwrap();
        m["solver"] = json!(solver);
        if !keep_gas {
            m.as_object_mut().unwrap().remove("gas");
        }
        std::fs::write(case.join(MARKER), m.to_string()).unwrap();
    };
    rewrite(&c.case, "rhoCentralFoam", false);
    let e = read(&mut c.s, &c.case).unwrap_err();
    assert!(e.contains("fieldcraft.json"), "a compressible solver without its gas: {e}");
    rewrite(&c.case, "pimpleFoam", true);
    let e = read(&mut c.s, &c.case).unwrap_err();
    assert!(e.contains("fieldcraft.json"), "an incompressible solver with a gas: {e}");
    std::fs::write(c.case.join(MARKER), &original).unwrap();
    assert!(read(&mut c.s, &c.case).is_ok(), "the case as written reads");
}

#[test]
fn a_compressible_frame_is_24_bytes_a_cell() {
    assert_eq!(super::read::frame_bytes_per_cell(true), 24);
    assert_eq!(super::read::frame_bytes_per_cell(false), 16);
}
