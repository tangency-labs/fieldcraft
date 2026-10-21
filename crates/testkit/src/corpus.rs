//! The meshing test corpus.
//!
//! Two kinds of part, both regenerated on every load (no geometry is stored):
//! - **Recipes**: `tests/corpus/<name>/recipe.json`, a SolidCraft command script (the file
//!   `solidcraft-cli run` takes) with a `corpus` object beside its `commands` giving the
//!   categories it stresses and, when it has one, its closed-form volume; `README.md` beside it is
//!   one line on what it stresses. The script runs in SolidCraft's engine, the design is saved and
//!   opened again through [`fieldcraft_geometry::open_design`], the path a person's design takes.
//! - **Analytic shapes** ([`Analytic`]): sphere, torus, cylinder, thick-walled pipe and a plate
//!   with a hole, built in code with their exact volume, the classic solver verification
//!   geometries.
//!
//! Every part's accepted mesher results are frozen in `tests/corpus/<name>/fixture.json`
//! ([`crate::fixture`]); an analytic shape's directory holds only that file.

use std::path::{Path, PathBuf};

use fieldcraft_geometry::{Design, open_design, sample, script_design};
use serde_json::{Value, json};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0}")]
    Io(String),
    #[error("{part}: {message}")]
    Part { part: String, message: String },
    #[error("{part}: {source}")]
    Geometry { part: String, source: fieldcraft_geometry::Error },
}

pub type Result<T> = std::result::Result<T, Error>;

/// What a part stresses in a mesher. A recipe names one or more.
pub const CATEGORIES: &[(&str, &str)] = &[
    ("thin-wall", "walls thin against the part's size"),
    ("small-feature", "features much smaller than the part"),
    ("fillet", "blends between faces"),
    ("genus", "through holes (genus > 0)"),
    ("sharp", "sharp edges and acute dihedral angles"),
    ("near-tangent", "faces meeting at a shallow angle"),
    ("boolean", "intersection curves from booleans"),
    ("sweepable", "extrusions and sweeps (one-to-one sweep meshing)"),
    ("blocky", "axis-aligned blocks (polycube-friendly)"),
    ("organic", "free-form, lofted or revolved spline surfaces"),
    ("assembly", "several bodies placed by component transforms"),
    ("partitioned", "split into cells that each mesh as a mapped block (conforming multi-block hex)"),
    ("analytic", "closed-form geometry for solver verification"),
];

/// Largest recipe or README read.
const MAX_FILE_BYTES: u64 = 4 << 20;
/// Most parts in the corpus directory.
const MAX_PARTS: usize = 10_000;

/// A part of the corpus.
#[derive(Debug, Clone, PartialEq)]
pub struct Part {
    /// The recipe's directory name, or the analytic shape's name (`analytic-…`).
    pub name: String,
    /// One line on what it stresses.
    pub description: String,
    /// Entries of [`CATEGORIES`].
    pub categories: Vec<String>,
    /// Closed-form volume (mm³), when the part has one.
    pub volume: Option<f64>,
    pub source: Source,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Source {
    /// The script, read from `recipe.json`.
    Recipe(String),
    Analytic(Analytic),
}

/// Closed-form shapes (mm).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Analytic {
    Sphere {
        radius: f64,
    },
    Torus {
        major: f64,
        minor: f64,
    },
    Cylinder {
        radius: f64,
        height: f64,
    },
    /// A thick-walled cylinder (Lamé).
    Pipe {
        inner: f64,
        outer: f64,
        length: f64,
    },
    /// A plate with a round through hole (Kirsch), built as a parametric design in code.
    PlateWithHole(sample::Plate),
}

impl Analytic {
    /// The exact volume.
    pub fn volume(&self) -> f64 {
        use std::f64::consts::PI;
        match *self {
            Analytic::Sphere { radius } => 4.0 / 3.0 * PI * radius.powi(3),
            Analytic::Torus { major, minor } => 2.0 * PI * PI * major * minor * minor,
            Analytic::Cylinder { radius, height } => PI * radius * radius * height,
            Analytic::Pipe { inner, outer, length } => PI * (outer * outer - inner * inner) * length,
            Analytic::PlateWithHole(p) => (p.width * p.depth - PI * p.hole_diameter * p.hole_diameter / 4.0) * p.thickness,
        }
    }

    /// The design as `.solidcraft` bytes.
    pub fn design(&self) -> fieldcraft_geometry::Result<Vec<u8>> {
        let script = match *self {
            Analytic::Sphere { radius } => json!([{"command": "solid.sphere", "params": {"radius": radius, "body_name": "Sphere"}}]),
            Analytic::Torus { major, minor } => json!([{"command": "solid.torus", "params": {"major": major, "minor": minor, "body_name": "Torus"}}]),
            Analytic::Cylinder { radius, height } => {
                json!([{"command": "solid.cylinder", "params": {"radius": radius, "height": height, "body_name": "Cylinder"}}])
            }
            Analytic::Pipe { inner, outer, length } => json!([
                {"command": "solid.cylinder", "params": {"radius": outer, "height": length, "body_name": "Pipe"}},
                {"command": "solid.cylinder", "params": {"radius": inner, "height": length, "operation": "cut", "targets": ["Pipe"]}}
            ]),
            Analytic::PlateWithHole(p) => return sample::plate_design(&p),
        };
        script_design(&script.to_string())
    }
}

/// The analytic shapes of the corpus.
pub fn analytic_parts() -> Vec<Part> {
    let plate = sample::Plate { width: 200.0, depth: 100.0, thickness: 5.0, hole_diameter: 20.0 };
    let shapes: [(&str, &str, &[&str], Analytic); 5] = [
        ("analytic-sphere", "Sphere r 20: curvature everywhere, no edges.", &[], Analytic::Sphere { radius: 20.0 }),
        ("analytic-torus", "Torus 30/8: genus 1 with doubly curved faces and no edges.", &["genus"], Analytic::Torus { major: 30.0, minor: 8.0 }),
        (
            "analytic-cylinder",
            "Cylinder r 10, h 40: the round bar of beam and torsion checks.",
            &["sweepable"],
            Analytic::Cylinder { radius: 10.0, height: 40.0 },
        ),
        (
            "analytic-pipe",
            "Thick-walled pipe 20/30, length 60: the Lamé cylinder.",
            &["genus", "sweepable"],
            Analytic::Pipe { inner: 20.0, outer: 30.0, length: 60.0 },
        ),
        ("analytic-plate-hole", "Plate 200 x 100 x 5 with a 20 mm hole: the Kirsch plate.", &["genus", "sweepable"], Analytic::PlateWithHole(plate)),
    ];
    shapes
        .into_iter()
        .map(|(name, description, more, shape)| Part {
            name: name.into(),
            description: description.into(),
            categories: std::iter::once("analytic").chain(more.iter().copied()).map(str::to_string).collect(),
            volume: Some(shape.volume()),
            source: Source::Analytic(shape),
        })
        .collect()
}

/// The workspace root (this crate is `crates/testkit` in it).
pub fn workspace_root() -> PathBuf {
    let here = Path::new(env!("CARGO_MANIFEST_DIR"));
    here.ancestors().nth(2).unwrap_or(here).to_path_buf()
}

/// The workspace's corpus directory, `tests/corpus`.
pub fn corpus_dir() -> PathBuf {
    workspace_root().join("tests").join("corpus")
}

pub(crate) fn read(path: &Path) -> Result<String> {
    let io = |e: std::io::Error| Error::Io(format!("{}: {e}", path.display()));
    let len = std::fs::metadata(path).map_err(io)?.len();
    if len > MAX_FILE_BYTES {
        return Err(Error::Io(format!("{}: larger than {MAX_FILE_BYTES} bytes", path.display())));
    }
    std::fs::read_to_string(path).map_err(io)
}

/// The recipe in `dir` (named after the directory).
pub fn read_recipe(dir: &Path) -> Result<Part> {
    let name = dir.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let bad = |message: String| Error::Part { part: name.clone(), message };
    let script = read(&dir.join("recipe.json"))?;
    let value: Value = serde_json::from_str(&script).map_err(|e| bad(format!("recipe.json is not JSON: {e}")))?;
    let meta = value.get("corpus").ok_or_else(|| bad("recipe.json has no `corpus` object".into()))?;
    let categories: Vec<String> = meta
        .get("categories")
        .and_then(Value::as_array)
        .map(|a| a.iter().filter_map(Value::as_str).map(str::to_string).collect())
        .unwrap_or_default();
    if categories.is_empty() {
        return Err(bad("`corpus.categories` names none".into()));
    }
    if let Some(c) = categories.iter().find(|c| !CATEGORIES.iter().any(|(k, _)| k == c)) {
        return Err(bad(format!("unknown category `{c}`")));
    }
    let volume = meta
        .get("volume_mm3")
        .map(|v| v.as_f64().filter(|v| v.is_finite() && *v > 0.0).ok_or_else(|| bad("`corpus.volume_mm3` must be a positive number".into())))
        .transpose()?;
    let readme = read(&dir.join("README.md"))?;
    let description = readme.lines().map(|l| l.trim_start_matches('#').trim()).find(|l| !l.is_empty()).unwrap_or_default().to_string();
    if description.is_empty() {
        return Err(bad("README.md is empty".into()));
    }
    Ok(Part { name, description, categories, volume, source: Source::Recipe(script) })
}

/// Every recipe under `dir`, by name.
pub fn read_recipes(dir: &Path) -> Result<Vec<Part>> {
    let rd = std::fs::read_dir(dir).map_err(|e| Error::Io(format!("{}: {e}", dir.display())))?;
    let mut dirs: Vec<PathBuf> = rd.flatten().map(|e| e.path()).filter(|p| p.join("recipe.json").is_file()).take(MAX_PARTS + 1).collect();
    if dirs.len() > MAX_PARTS {
        return Err(Error::Io(format!("{}: more than {MAX_PARTS} parts", dir.display())));
    }
    dirs.sort();
    dirs.iter().map(|d| read_recipe(d)).collect()
}

/// The whole corpus: the recipes in [`corpus_dir`], then the analytic shapes.
pub fn parts() -> Result<Vec<Part>> {
    let mut parts = read_recipes(&corpus_dir())?;
    parts.extend(analytic_parts());
    Ok(parts)
}

/// A part's design as `.solidcraft` bytes: its recipe run (or its shape generated).
pub fn design_bytes(part: &Part) -> Result<Vec<u8>> {
    match &part.source {
        Source::Recipe(script) => script_design(script),
        Source::Analytic(shape) => shape.design(),
    }
    .map_err(|source| Error::Geometry { part: part.name.clone(), source })
}

/// Build a part: run its recipe (or generate the shape) and open the design.
pub fn load(part: &Part) -> Result<Design> {
    open_design(&design_bytes(part)?).map_err(|source| Error::Geometry { part: part.name.clone(), source })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analytic_shapes_have_their_closed_form_volume() {
        for part in analytic_parts() {
            let design = load(&part).unwrap();
            assert!(design.errors.is_empty(), "{}: {:?}", part.name, design.errors);
            let expect = part.volume.unwrap();
            let got: f64 = design.bodies.iter().map(|b| b.measure().unwrap().volume).sum();
            assert!((got - expect).abs() < 1e-3 * expect, "{}: {got} vs {expect}", part.name);
        }
    }

    /// Every committed recipe builds without a failed feature, and its volume matches the closed
    /// form where it states one.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn every_recipe_builds() {
        let parts = read_recipes(&corpus_dir()).unwrap();
        assert!(parts.len() >= 15, "{}", parts.len());
        for part in &parts {
            let design = load(part).unwrap_or_else(|e| panic!("{e}"));
            assert!(design.errors.is_empty(), "{}: {:?}", part.name, design.errors);
            assert!(!design.bodies.is_empty(), "{}: no bodies", part.name);
            if let Some(expect) = part.volume {
                let got: f64 = design.bodies.iter().map(|b| b.measure().unwrap().volume).sum();
                assert!((got - expect).abs() < 1e-3 * expect, "{}: {got} vs {expect}", part.name);
            }
        }
        // Every category is covered.
        for (c, _) in CATEGORIES {
            assert!(parts.iter().chain(&analytic_parts()).any(|p| p.categories.iter().any(|x| x == c)), "no part covers {c}");
        }
    }

    /// The assembly's clamps sit where their occurrences put them: one moved, one instance moved
    /// and turned 90 degrees about Z.
    #[test]
    fn assembly_bodies_are_placed() {
        let part = read_recipe(&corpus_dir().join("bolted-assembly")).unwrap();
        let design = load(&part).unwrap();
        let mut boxes: Vec<([f64; 3], [f64; 3])> = design
            .bodies
            .iter()
            .map(|b| {
                let m = b.measure().unwrap();
                (m.bbox_min, m.bbox_max)
            })
            .collect();
        boxes.sort_by(|a, b| a.0[0].total_cmp(&b.0[0]).then(a.1[0].total_cmp(&b.1[0])));
        let expect = [([0.0, 0.0, 0.0], [100.0, 60.0, 5.0]), ([10.0, 10.0, 5.0], [30.0, 20.0, 13.0]), ([65.0, 5.0, 5.0], [75.0, 25.0, 13.0])];
        assert_eq!(boxes.len(), expect.len());
        for (got, want) in boxes.iter().zip(&expect) {
            let close = |a: &[f64; 3], b: &[f64; 3]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-6);
            assert!(close(&got.0, &want.0) && close(&got.1, &want.1), "{got:?} vs {want:?}");
        }
    }

    #[test]
    fn hostile_recipes_are_an_error() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-corpus-test-{}", std::process::id()));
        let part = dir.join("bad");
        std::fs::create_dir_all(&part).unwrap();
        std::fs::write(part.join("README.md"), "A bad part.\n").unwrap();
        for recipe in [
            "",
            "{",
            r#"{"commands": []}"#,
            r#"{"corpus": {"categories": []}, "commands": []}"#,
            r#"{"corpus": {"categories": ["nope"]}, "commands": []}"#,
            r#"{"corpus": {"categories": ["blocky"], "volume_mm3": -1}, "commands": []}"#,
        ] {
            std::fs::write(part.join("recipe.json"), recipe).unwrap();
            assert!(read_recipe(&part).is_err(), "{recipe}");
        }
        std::fs::write(part.join("recipe.json"), r#"{"corpus": {"categories": ["blocky"]}, "commands": [{"command": "no.such"}]}"#).unwrap();
        assert!(load(&read_recipe(&part).unwrap()).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
