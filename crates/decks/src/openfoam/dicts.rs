//! The fixed dictionaries of a case: meshing controls, discretisation schemes, solution
//! controls, run controls and the turbulence model. Every value is a named constant that states
//! the public page documenting it, or how it is derived from the case; the few keyword names the
//! pages read do not show are marked where they are written.
//!
//! Sources, all public documentation of ESI OpenFOAM (nothing of its source or tutorials):
//! - UG: the User Guide, `openfoam.com/documentation/user-guide` (§3.2 parallel runs, §4.2
//!   boundaries, §4.4 snappyHexMesh, §5.3 turbulence, §6.1 controlDict, §6.2 fvSchemes, §6.3
//!   fvSolution);
//! - TG: the Tutorial Guide, `openfoam.com/documentation/tutorial-guide`, §2.1 (the cavity);
//! - DOC: `doc.openfoam.com`, the pages of snappyHexMesh (castellation, snapping, mesh quality,
//!   global settings), the linear solvers and residuals, the surface-normal gradient schemes, the
//!   function objects, the wall-distance methods, and the Smagorinsky model, its length scales and
//!   Van Driest's damping.
//!
//! Read on: UG and TG from their Internet Archive snapshots (2026-11; the user guide carries no
//! version of its own); DOC from the v2212 and v2312 pages. No v2406 page was readable.
//!
//! A value the guides give only in an "Example" column is called an example here: the guides
//! state no default for it.

use std::fmt::Write as _;

use fieldcraft_model::Turbulence;

use super::{RefineBox, Surface, dictionary, num, unsupported, vector};
use crate::Result;

/// The incompressible transient solver the case runs (UG §A.1; DOC "pimpleFoam": fields `p`,
/// the kinematic pressure, and `U`).
pub const SOLVER: &str = "pimpleFoam";

/// The default Courant limit of the adaptive time step, and the largest a case takes. TG
/// §2.1.1.4: "To achieve temporal accuracy and numerical stability when running icoFoam, a
/// Courant number of less than 1 is required"; the time step is chosen for "a Courant number less
/// than or equal to 1 throughout the domain". With one outer corrector `pimpleFoam` runs the same
/// predictor–corrector sequence as `icoFoam` (Issa's PISO).
pub const MAX_COURANT: f64 = 1.0;

/// `surfaceFeatureExtract`: edges whose two faces meet at less than this angle (degrees) are
/// feature edges. UG §4.4.3 (its example, `-includedAngle 150`).
const INCLUDED_ANGLE: f64 = 150.0;

// castellatedMeshControls (UG §4.4.3 Table 4.5; DOC snappyHexMesh "Castellation").

/// The largest cell count an OpenFOAM build with 32-bit labels can hold; the cell limits are
/// capped at it.
pub(super) const MAX_LABEL: u64 = 2_147_483_647;
/// `minRefinementCells`: refinement stops when no more than this many cells are left to refine
/// (DOC "Castellation"). 0: every requested level is reached.
const MIN_REFINEMENT_CELLS: u64 = 0;
/// `maxLoadUnbalance`: the cell imbalance between processors that triggers rebalancing in a
/// parallel meshing run (UG Table 4.5, its example). `Allrun` meshes in serial: no effect.
const MAX_LOAD_UNBALANCE: f64 = 0.1;
/// `nCellsBetweenLevels`: buffer cells between two refinement levels. DOC "Castellation": "A good
/// compromise between limiting the cell size jump and number of cells is with
/// nCellsBetweenLevels 3."
const CELLS_BETWEEN_LEVELS: u32 = 3;
/// `resolveFeatureAngle`: where surface normals differ by more than this (degrees), the surface's
/// second level applies (UG Table 4.5). The same definition of a feature as [`INCLUDED_ANGLE`]:
/// faces meeting at 150° have normals 30° apart. The two levels written are equal, so it changes
/// no refinement.
const RESOLVE_FEATURE_ANGLE: f64 = 180.0 - INCLUDED_ANGLE;

// snapControls (UG §4.4.6 Table 4.6; DOC snappyHexMesh "Snapping", whose block of values is what
// it says to start from: "The good practice is to start with the default values").

/// `nSmoothPatch`: patch smoothing iterations before the surface is matched.
const SMOOTH_PATCH: u32 = 3;
/// `tolerance`: points within this many local edge lengths of the surface are attracted to it
/// (DOC "Snapping": "by default set to 2").
const SNAP_TOLERANCE: f64 = 2.0;
/// `nSolveIter`: mesh displacement relaxation iterations.
const SNAP_SOLVE_ITER: u32 = 30;
/// `nRelaxIter`: the most snapping relaxation iterations.
const SNAP_RELAX_ITER: u32 = 5;
/// `nFeatureSnapIter`: feature-edge snapping iterations (none written disables feature snapping).
const FEATURE_SNAP_ITER: u32 = 10;

/// `meshQualityControls`, written in full so the limits do not come from the installation: UG
/// §4.4.8 Table 4.8, each entry with the table's value, except `minVol` and `minFaceFlatness`.
const QUALITY_LIMITS: [(&str, f64); 13] = [
    // Maximum face non-orthogonality (degrees); 180 disables.
    ("maxNonOrtho", 65.0),
    // Maximum boundary face skewness; 0 disables.
    ("maxBoundarySkewness", 20.0),
    // Maximum internal face skewness; 0 disables.
    ("maxInternalSkewness", 4.0),
    // Maximum concaveness (degrees); 180 disables.
    ("maxConcave", 80.0),
    // Minimum pyramid volume (m³). The table's 1e-13 is an absolute volume, larger than a whole
    // pyramid of any cell under 84 µm; 0 keeps the limit's purpose at every scale: no inverted
    // pyramid. Flat cells are bounded by the dimensionless limits below.
    ("minVol", 0.0),
    // Minimum quality of the tetrahedra of the cell decomposition; a small positive value.
    ("minTetQuality", 1e-13),
    // Minimum face area; negative: not checked.
    ("minArea", -1.0),
    // Minimum face twist; -1 disables.
    ("minTwist", 0.05),
    // Minimum normalised cell determinant; 1 is a hexahedron, 0 an illegal cell.
    ("minDeterminant", MIN_DETERMINANT),
    // Minimum face interpolation weight, 0 to 0.5.
    ("minFaceWeight", 0.05),
    // Minimum ratio of the volumes of neighbouring cells, 0 to 1.
    ("minVolRatio", 0.01),
    // Minimum twist of a face's triangles; -1 disables.
    ("minTriangleTwist", -1.0),
    // Minimum face flatness; -1: not checked. DOC "Mesh quality" (optional, shown as -1); written
    // so the installation's default does not decide.
    ("minFaceFlatness", -1.0),
];
/// `minDeterminant`: the least normalised cell determinant a cell may have (UG Table 4.8; 1 is a
/// hexahedron, 0 an illegal cell).
const MIN_DETERMINANT: f64 = 0.001; // unit-factors: dimensionless, a normalised cell determinant

/// `nSmoothScale`: iterations that spread a displacement error to neighbouring points (Table 4.8).
const SMOOTH_SCALE: u32 = 4;
/// `errorReduction`: the factor a displacement is scaled back by at a point that breaks a limit
/// (Table 4.8).
const ERROR_REDUCTION: f64 = 0.75;

/// `mergeTolerance`: points closer than this fraction of the initial mesh's bounding box are one
/// point (UG §4.4.1 Table 4.4, its example). DOC "Global settings": "the write tolerance needs to
/// be higher than this", which [`WRITE_PRECISION`] digits give (10⁻⁸ against 10⁻⁶).
pub(super) const MERGE_TOLERANCE: f64 = 1e-6; // unit-factors: relative to the initial mesh's bounding box

/// `Ce` of OpenFOAM's Smagorinsky model: the documented default (DOC "Smagorinsky", 1.048, with
/// Ck 0.094 there). The value is on DOC v2312's "Smagorinsky" page; the citation, Yoshizawa (1986,
/// Phys. Fluids 29:2152), is on DOC v2312's "k equation" page; that paper was not read
/// (paywalled). The model is ν_t = Ck Δ √k with the subgrid energy k in local balance,
/// Ce k^(3/2)/Δ = 2 ν_t D:D; for incompressible flow that is ν_t = Ck √(Ck/Ce) Δ² |S| with
/// |S| = √(2 D:D), so only Cs² = Ck √(Ck/Ce) reaches the flow and the value of Ce is free. It is
/// written as the documented default; the step's constant Cs gives Ck = (Cs² √Ce)^(2/3).
const CE: f64 = 1.048;

/// `writePrecision`: significant digits of written values. Results are read back into `f32`
/// frames, and 9 significant decimal digits round-trip every `f32` (IEEE 754-2008 §5.12.2).
const WRITE_PRECISION: u32 = 9;

/// ψ of the `limited` surface-normal gradient scheme. UG §6.2.2 eq. (6.1): 0.333 holds the
/// explicit non-orthogonal correction to half the orthogonal part, so each
/// non-orthogonal corrector at least halves the lag of the correction on the worst face (the mesh
/// allows faces up to `maxNonOrtho`, where the unlimited correction is twice the orthogonal part).
const NON_ORTHOGONAL_LIMIT: f64 = 0.333;

/// Convection under LES: central differencing, second order, adds no numerical dissipation to the
/// subgrid model's (Mittal & Moin 1997, AIAA J. 35:1415–1417; UG §6.2.5 Table 6.7, `linear`).
const LES_CONVECTION: &str = "Gauss linear";

/// Convection of a laminar step: linear upwind differencing, bounded second order (UG Table 6.7,
/// "First/second order, bounded"; Warming & Beam, AIAA J. 14:1241, 1976), with the velocity
/// gradient it corrects by.
const LAMINAR_CONVECTION: &str = "Gauss linearUpwind grad(U)";

/// von Kármán's constant of the van Driest length scale (DOC "Van Driest": optional entry
/// `kappa`, 0.41), written out so the installation's default does not decide.
const KAPPA: f64 = 0.41;

/// `nCorrectors`: Issa's PISO is an implicit momentum predictor and two pressure correctors per
/// time step (R. I. Issa, J. Comput. Phys. 62:40–65, 1986); UG §6.3.3: "PISO requires more than 1".
const PISO_CORRECTORS: u32 = 2;

/// `nNonOrthogonalCorrectors` of a mesh with snapped or refined cells: the explicit correction is
/// evaluated once more from the pressure of the same corrector. UG §6.3.3: "0 for an orthogonal
/// mesh and increasing with the degree of non-orthogonality".
const NON_ORTHOGONAL_CORRECTORS: u32 = 1;

/// The tolerance of the linear solvers, on OpenFOAM's normalised residual (UG §6.3.1.1; DOC
/// "Residuals": Σ|b − Ax| over Σ(|Ax − Ax̄| + |b − Ax̄|), dimensionless). Results are read back
/// into `f32` frames: a residual below an `f32`'s relative spacing changes nothing they keep.
fn solver_tolerance() -> f64 {
    f64::from(f32::EPSILON)
}

/// `system/surfaceFeatureExtractDict`: one sub-dictionary per obstacle surface, named after its STL
/// (API guide v2112, "surfaceFeatureExtract.C File Reference", Detailed Description, the page's text
/// only, read as public documentation on the owner's ruling of 2026-11-24): `extractionMethod
/// extractFromSurface` extracts features from the surface's geometry; `includedAngle` marks edges
/// whose adjacent normals meet at less than it; `geometricTestOnly no` also marks region edges, which
/// here are the borders between the model's named faces (one STL solid per face), so every face's
/// border is a feature the snap follows (the page states no default, so it is written); `writeObj`
/// is an output switch for post-processing, off.
pub(super) fn feature_extract(surfaces: &[Surface]) -> String {
    let mut s = String::new();
    for surface in surfaces {
        let _ = writeln!(
            s,
            "{}.stl\n{{\n    extractionMethod extractFromSurface;\n    includedAngle {};\n    geometricTestOnly no;\n    writeObj no;\n}}\n",
            surface.name,
            num(INCLUDED_ANGLE)
        );
    }
    dictionary("surfaceFeatureExtractDict", &s)
}

/// The most cells refinement can make: every background cell at the finest level, capped at the
/// largest label. `maxGlobalCells` stops refinement where it binds (DOC "Castellation"), which
/// would give a mesh that is not the one asked for; at this value it cannot bind.
fn cell_limit(background_cells: [u64; 3], finest_level: u32) -> u64 {
    let background = background_cells.iter().fold(1u64, |n, c| n.saturating_mul(*c));
    let per_cell = 8u64.checked_pow(finest_level).unwrap_or(u64::MAX);
    background.saturating_mul(per_cell).min(MAX_LABEL)
}

/// `system/snappyHexMeshDict`. `level` is the obstacle surfaces' refinement level, `boxes` the
/// refinement boxes with theirs, `background` the background cell (m) and `background_cells` the
/// background block's cells along x, y and z.
pub(super) fn snappy(
    surfaces: &[Surface],
    level: u32,
    boxes: &[(RefineBox, u32)],
    location: [f64; 3],
    background: f64,
    background_cells: [u64; 3],
) -> String {
    let mut geometry = String::new();
    let mut features = String::new();
    let mut refinement = String::new();
    for surface in surfaces {
        let name = &surface.name;
        let _ = writeln!(geometry, "    {name}.stl {{ type triSurfaceMesh; name {name}; }}");
        let _ = writeln!(features, "        {{ file \"{name}.eMesh\"; level {level}; }}");
        let _ = writeln!(refinement, "        {name} {{ level ({level} {level}); patchInfo {{ type wall; }} }}");
    }
    let mut regions = String::new();
    for (i, (b, box_level)) in boxes.iter().enumerate() {
        let _ = writeln!(geometry, "    refine{i} {{ type searchableBox; min {}; max {}; }}", vector(b.min), vector(b.max));
        // An `inside` region's distance "is not required so is ignored (but it must be specified)"
        // (UG §4.4.5): the cell size the level gives.
        let cell = background / f64::from(1u32 << *box_level);
        let _ = writeln!(regions, "        refine{i} {{ mode inside; levels (({} {box_level})); }}", num(cell));
    }
    let finest_level = boxes.iter().map(|(_, l)| *l).fold(level, u32::max);
    // `maxLocalCells` changes the balancing strategy of a parallel meshing run above it (DOC
    // "Castellation"); at the global limit it never does.
    let limit = cell_limit(background_cells, finest_level);
    let mut quality = String::new();
    for (entry, value) in QUALITY_LIMITS {
        let _ = writeln!(quality, "    {entry} {};", num(value));
    }
    // `singleRegionName false`: a one-face body's patch is `<body>_<face>` like every other (DOC
    // "Global settings"). No `addLayersControls` while `addLayers` is false, and no `relaxed`
    // limits, which only layer addition reads (UG Table 4.8).
    let body = format!(
        "castellatedMesh true;\nsnap true;\naddLayers false;\nsingleRegionName false;\n\n\
geometry\n{{\n{geometry}}}\n\n\
castellatedMeshControls\n{{\n    maxLocalCells {limit};\n    maxGlobalCells {limit};\n    minRefinementCells {MIN_REFINEMENT_CELLS};\n    maxLoadUnbalance {};\n    nCellsBetweenLevels {CELLS_BETWEEN_LEVELS};\n\
    features\n    (\n{features}    );\n    refinementSurfaces\n    {{\n{refinement}    }}\n    resolveFeatureAngle {};\n\
    refinementRegions\n    {{\n{regions}    }}\n    locationInMesh {};\n    allowFreeStandingZoneFaces true;\n}}\n\n\
snapControls\n{{\n    nSmoothPatch {SMOOTH_PATCH};\n    tolerance {};\n    nSolveIter {SNAP_SOLVE_ITER};\n    nRelaxIter {SNAP_RELAX_ITER};\n    nFeatureSnapIter {FEATURE_SNAP_ITER};\n    implicitFeatureSnap false;\n    explicitFeatureSnap true;\n    multiRegionFeatureSnap false;\n}}\n\n\
meshQualityControls\n{{\n{quality}    nSmoothScale {SMOOTH_SCALE};\n    errorReduction {};\n}}\n\nmergeTolerance {};\n",
        num(MAX_LOAD_UNBALANCE),
        num(RESOLVE_FEATURE_ANGLE),
        vector(location),
        num(SNAP_TOLERANCE),
        num(ERROR_REDUCTION),
        num(MERGE_TOLERANCE),
    );
    dictionary("snappyHexMeshDict", &body)
}

/// `constant/transportProperties`: a Newtonian fluid of kinematic viscosity `nu` (m²/s; TG
/// §2.1.1.3: "The keyword for kinematic viscosity is nu"). `transportModel Newtonian` is a keyword
/// of the format; not in the documentation pages read (U3).
pub(super) fn transport(nu: f64) -> String {
    dictionary("transportProperties", &format!("transportModel Newtonian;\nnu {};\n", num(nu)))
}

/// `constant/turbulenceProperties` (UG §5.3, Tables 5.3 and 5.4): no model, or OpenFOAM's
/// Smagorinsky model with `Ck` from the step's constant ([`CE`]) and the cube root of the cell
/// volume as its length scale (DOC "Cube-root volume": Δ = c V^(1/3); on a hexahedral cell of side
/// h that is h, the lattice's Δx, with `deltaCoeff` c = 1). Van Driest's damping (DOC "Van Driest")
/// replaces that length scale near walls, with the step's constant as `Cdelta`.
pub(super) fn turbulence_properties(t: Option<&Turbulence>, notes: &mut Vec<String>) -> Result<String> {
    let body = match t {
        None => "simulationType laminar;\n".to_string(),
        Some(Turbulence::Smagorinsky { constant, van_driest, wall_model, .. }) => {
            if !(constant.is_finite() && *constant > 0.0) {
                return Err(unsupported("a Smagorinsky constant that is not positive"));
            }
            if wall_model.is_some() {
                notes.push("the wall model is not mapped: walls are resolved no-slip".into());
            }
            let ck = libm::pow(constant * constant * CE.sqrt(), 2.0 / 3.0);
            // DOC "Van Driest": Δ = min(κ y D / Cdelta, Δ_g), D = 1 − exp(−y⁺/A⁺). With Cdelta the
            // step's constant C, C Δ = min(κ y D, C Δ_g): Van Driest's damped mixing length at the
            // wall (J. Aeronaut. Sci. 23:1007, 1956) under the Smagorinsky length.
            let delta = match van_driest {
                Some(a) => {
                    if !(a.is_finite() && *a > 0.0) {
                        return Err(unsupported("a Van Driest constant A+ that is not positive"));
                    }
                    notes.push("Van Driest's damping is OpenFOAM's law, min(κ y D, C Δ), not the lattice's C Δx D".into());
                    format!(
                        "delta vanDriest;\n    vanDriestCoeffs\n    {{\n        delta cubeRootVol;\n        cubeRootVolCoeffs {{ deltaCoeff 1; }}\n        kappa {};\n        Aplus {};\n        Cdelta {};\n    }}",
                        num(KAPPA),
                        num(*a),
                        num(*constant)
                    )
                }
                None => "delta cubeRootVol;\n    cubeRootVolCoeffs { deltaCoeff 1; }".to_string(),
            };
            format!(
                "simulationType LES;\n\nLES\n{{\n    LESModel Smagorinsky;\n    turbulence on;\n    printCoeffs on;\n    {delta}\n    SmagorinskyCoeffs {{ Ck {}; Ce {}; }}\n}}\n",
                num(ck),
                num(CE)
            )
        }
        Some(Turbulence::Vreman { .. }) => return Err(unsupported("Vreman's subgrid model (not mapped yet)")),
        Some(Turbulence::SpalartAllmaras { .. }) => return Err(unsupported("the Spalart–Allmaras model (not mapped yet)")),
    };
    Ok(dictionary("turbulenceProperties", &body))
}

/// `system/fvSchemes` (UG §6.2): second order in time (`backward`, Table 6.8) and space (`Gauss
/// linear`, Tables 6.3 and 6.5). Convection: central under LES ([`LES_CONVECTION`], owner-approved
/// fork F) and linear upwind when laminar ([`LAMINAR_CONVECTION`]). The surface-normal gradient
/// with the limited non-orthogonal correction ([`NON_ORTHOGONAL_LIMIT`]; syntax
/// `limited corrected <psi>`: DOC "Limited surface-normal gradient scheme"). No limiter on
/// `grad(U)`: it is the strain rate the subgrid model reads. `wall_distance`: the wall distance
/// that Van Driest's damping reads, by the mesh-wave method (DOC "Wall distance calculation
/// methods").
pub(super) fn fv_schemes(les: bool, wall_distance: bool) -> String {
    let psi = num(NON_ORTHOGONAL_LIMIT);
    let convection = if les { LES_CONVECTION } else { LAMINAR_CONVECTION };
    let wall = if wall_distance { "wallDist { method meshWave; }\n" } else { "" };
    dictionary(
        "fvSchemes",
        &format!(
            "ddtSchemes {{ default backward; }}\n\
gradSchemes {{ default Gauss linear; }}\n\
divSchemes {{ default Gauss linear; div(phi,U) {convection}; }}\n\
laplacianSchemes {{ default Gauss linear limited corrected {psi}; }}\n\
interpolationSchemes {{ default linear; }}\n\
snGradSchemes {{ default limited corrected {psi}; }}\n{wall}"
        ),
    )
}

/// `system/fvSolution` (UG §6.3). Pressure: multigrid with the Gauss–Seidel smoother (§6.3.1.3:
/// "Generally GaussSeidel is the most reliable option"; DOC "GAMG Solver": `solver`, `smoother`,
/// `tolerance` and `relTol` are its mandatory entries). Velocity: the smooth solver with symmetric
/// Gauss–Seidel sweeps, the residual checked after each (`nSweeps 1`); with Co ≤ 1 the momentum
/// matrix is diagonally dominant. `relTol 0` (§6.3.1.1: "In transient simulations, it is usual to
/// set the solver relative tolerance to 0 to force the solution to converge to the solver
/// tolerance in each time step"). `PIMPLE` with one outer corrector is Issa's sequence
/// ([`PISO_CORRECTORS`]). `pRefCell`/`pRefValue`: the pressure level of a domain no face fixes
/// the pressure of (§6.3.3.1); the initial pressure is 0. `orthogonal`: the mesh has no snapped
/// and no refined cell. Keywords of the format that are not in the documentation pages read
/// (U3): `momentumPredictor` (the predictor is Issa's scheme's first step) and the `UFinal` of the
/// `"(U|UFinal)"` solver key (it follows `pFinal`).
pub(super) fn fv_solution(orthogonal: bool) -> String {
    let tolerance = num(solver_tolerance());
    let non_orthogonal = if orthogonal { 0 } else { NON_ORTHOGONAL_CORRECTORS };
    dictionary(
        "fvSolution",
        &format!(
            "solvers\n{{\n    p {{ solver GAMG; smoother GaussSeidel; tolerance {tolerance}; relTol 0; }}\n    pFinal {{ $p; }}\n\
    \"(U|UFinal)\" {{ solver smoothSolver; smoother symGaussSeidel; nSweeps 1; tolerance {tolerance}; relTol 0; }}\n}}\n\n\
PIMPLE\n{{\n    momentumPredictor yes;\n    nOuterCorrectors 1;\n    nCorrectors {PISO_CORRECTORS};\n    nNonOrthogonalCorrectors {non_orthogonal};\n    pRefCell 0;\n    pRefValue 0;\n}}\n"
        ),
    )
}

/// What `system/controlDict` is written from.
pub(super) struct Control {
    pub solver: &'static str,
    /// The run's length (s).
    pub duration: f64,
    /// Results written over the run, evenly (at least 1).
    pub frames: u32,
    pub courant: f64,
    /// The finest cell (m).
    pub finest: f64,
    /// The largest speed the boundaries impose (m/s; zero for a fluid at rest).
    pub speed: f64,
    /// Where the step's running means start (s), if it averages.
    pub average_from: Option<f64>,
}

/// `system/controlDict` (UG §6.1 Table 6.1). The run starts from the step's initial state every
/// time (`startFrom startTime`: `Allrun` re-meshes, so a restart from a later time would put old
/// fields on a new mesh). The time step adapts to the Courant limit (`adjustTimeStep`, `maxCo`),
/// is never longer than the time between two frames (`maxDeltaT`; a fluid at rest has no Courant
/// limit) and lands on the frame times (`adjustableRunTime`). Time names keep enough digits to
/// tell two neighbouring frames apart: they differ by a relative 1/frames at the end, and d digits
/// resolve 10^(1−d). A running mean (DOC "fieldAverage"; `writeTime`: "alongside standard field
/// output", DOC "Function objects") starts where the step's own averaging does: the velocity's
/// mean and fluctuation, as the step reports them, and the pressure's mean.
pub(super) fn control(c: &Control) -> String {
    let interval = c.duration / f64::from(c.frames);
    // The first step at the Courant limit on the finest cell for the boundaries' speed; a fluid
    // at rest has no limit: the write interval.
    let first = if c.speed > 0.0 { (c.courant * c.finest / c.speed).min(interval) } else { interval };
    let time_precision = WRITE_PRECISION.max(c.frames.max(1).ilog10() + 2);
    let functions = match c.average_from {
        Some(start) => format!(
            "\nfunctions\n{{\n    fieldAverage\n    {{\n        type fieldAverage;\n        libs (fieldFunctionObjects);\n        timeStart {};\n\
        executeControl timeStep;\n        executeInterval 1;\n        writeControl writeTime;\n        writeInterval 1;\n\
        fields ( U {{ mean on; prime2Mean on; base time; }} p {{ mean on; prime2Mean off; base time; }} );\n    }}\n}}\n",
            num(start)
        ),
        None => String::new(),
    };
    dictionary(
        "controlDict",
        &format!(
            "application {};\nstartFrom startTime;\nstartTime 0;\nstopAt endTime;\nendTime {};\ndeltaT {};\n\
adjustTimeStep yes;\nmaxCo {};\nmaxDeltaT {};\n\
writeControl adjustableRunTime;\nwriteInterval {};\npurgeWrite 0;\nwriteFormat ascii;\nwritePrecision {WRITE_PRECISION};\nwriteCompression off;\ntimeFormat general;\ntimePrecision {time_precision};\n\
runTimeModifiable true;\n{functions}",
            c.solver,
            num(c.duration),
            num(first),
            num(c.courant),
            num(interval),
            num(interval),
        ),
    )
}

/// `system/decomposeParDict` (UG §3.2.1 Table 3.1): `scotch` "requires no geometric input from the
/// user and attempts to minimise the number of processor boundaries".
pub(super) fn decompose(processes: u32) -> String {
    dictionary("decomposeParDict", &format!("numberOfSubdomains {processes};\nmethod scotch;\n"))
}

#[cfg(test)]
mod tests {
    use super::super::grammar;
    use super::*;

    #[test]
    fn the_cell_limit_is_every_cell_at_the_finest_level_within_a_label() {
        assert_eq!(cell_limit([20, 8, 8], 0), 1280);
        assert_eq!(cell_limit([5, 2, 2], 2), 1280);
        assert_eq!(cell_limit([1000, 1000, 1000], 10), MAX_LABEL);
        assert_eq!(cell_limit([u64::MAX, 2, 2], 40), MAX_LABEL);
    }

    #[test]
    fn a_running_mean_starts_where_the_step_averages() {
        let input = |average_from: Option<f64>, frames: u32| Control {
            solver: SOLVER,
            duration: 2.0,
            frames,
            courant: 1.0,
            finest: 0.01,
            speed: 0.0,
            average_from,
        };
        let d = grammar::check_file(&control(&input(Some(0.5), 4))).unwrap();
        assert_eq!(d.text("functions/fieldAverage/type").as_deref(), Some("fieldAverage"));
        assert_eq!(d.number("functions/fieldAverage/timeStart"), Some(0.5));
        assert_eq!(d.text("functions/fieldAverage/writeControl").as_deref(), Some("writeTime"));
        assert!(d.get("functions/fieldAverage/fields").is_some());
        // A fluid at rest: the first step is the write interval, and so is the longest step.
        assert_eq!(d.number("deltaT"), Some(0.5));
        assert_eq!(d.number("maxDeltaT"), Some(0.5));
        // 2 000 000 frames have 7 digits: 8 would tell neighbours apart, 9 are written anyway.
        let many = grammar::check_file(&control(&input(None, 2_000_000))).unwrap();
        assert_eq!(many.number("timePrecision"), Some(9.0));
        let more = grammar::check_file(&control(&input(None, u32::MAX))).unwrap();
        assert_eq!(more.number("timePrecision"), Some(11.0));
    }
}
