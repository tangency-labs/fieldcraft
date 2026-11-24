use super::*;

/// A 10 × 4 × 4 mm channel (SI) with a 2 mm cube in it, air, a velocity inlet at x_min and a
/// pressure outlet at x_max.
fn channel(extra: serde_json::Value) -> Model {
    let mut domain = serde_json::json!({
        "material": "Air", "min": [0.0, 0.0, 0.0], "max": [0.01, 0.004, 0.004], "spacing": 0.0005,
        "faces": {"x_min": {"type": "velocity", "velocity": [2.0, 0.0, 0.0]},
                  "x_max": {"type": "pressure", "pressure": 0.0},
                  "y_min": {"type": "symmetry"}},
        "obstacles": ["Cube"]
    });
    if let (Some(d), Some(e)) = (domain.as_object_mut(), extra.as_object()) {
        for (k, v) in e {
            d.insert(k.clone(), v.clone());
        }
    }
    let model = serde_json::json!({
        "materials": {"Air": {"density": 1.2, "fluid": {"type": "newtonian", "viscosity": 1.8e-5}}},
        "fluid_domains": {"Flow": domain},
        "steps": [{"name": "Jet", "procedure": {"type": "fluid", "domain": "Flow", "end_time": 0.01, "frames": 5,
                    "turbulence": {"type": "smagorinsky", "constant": 0.16}}, "output": []}]
    });
    serde_json::from_value(model).unwrap()
}

/// A cube from 4 to 6 mm in x and 1 to 3 mm in y and z: 12 triangles, a face per side.
fn cube() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<u32>, Vec<String>) {
    let (x0, x1, y0, y1) = (0.004, 0.006, 0.001, 0.003);
    let p = vec![[x0, y0, y0], [x1, y0, y0], [x1, y1, y0], [x0, y1, y0], [x0, y0, y1], [x1, y0, y1], [x1, y1, y1], [x0, y1, y1]];
    let quads = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]];
    let mut t = Vec::new();
    let mut f = Vec::new();
    for (i, q) in quads.iter().enumerate() {
        t.push([q[0], q[1], q[2]]);
        t.push([q[0], q[2], q[3]]);
        f.extend([i as u32, i as u32]);
    }
    let names = ["F1:-z", "F1:+z", "F1:-y", "F1:+y", "F1:-x", "F1:+x"].map(String::from).to_vec();
    (p, t, f, names)
}

/// The default options around a point in the fluid by the outlet, at a finest cell's centre.
fn options() -> CaseOptions {
    CaseOptions::new(MeshOptions::default(), [0.00925, 0.00225, 0.00225])
}

fn write(model: &Model, options: &CaseOptions) -> Result<Case> {
    let (p, t, f, n) = cube();
    let o = ObstacleSurface { body: "Cube", positions: &p, triangles: &t, tri_face: &f, faces: &n };
    write_case(model, None, &[o], options)
}

#[test]
fn writes_a_complete_case_with_a_patch_per_face() {
    let case = write(&channel(serde_json::json!({})), &options()).unwrap();
    let names: Vec<&str> = case.files.keys().map(String::as_str).collect();
    for f in [
        "0/U",
        "0/p",
        "0/nut",
        "Allrun",
        "case.foam",
        "constant/transportProperties",
        "constant/triSurface/Cube.stl",
        "constant/turbulenceProperties",
        "system/blockMeshDict",
        "system/controlDict",
        "system/decomposeParDict",
        "system/fvSchemes",
        "system/fvSolution",
        "system/snappyHexMeshDict",
        "system/surfaceFeatureExtractDict",
    ] {
        assert!(names.contains(&f), "{f} missing: {names:?}");
    }
    assert_eq!(case.background_cells, [20, 8, 8], "a 0.5 mm background (cell factor 1) over 10 × 4 × 4 mm");
    let block = &case.files["system/blockMeshDict"];
    assert!(block.contains("hex (0 1 2 3 4 5 6 7) (20 8 8)"), "{block}");
    assert!(block.contains("y_min\n    {\n        type symmetryPlane;"), "{block}");
    assert!(block.contains("y_max\n    {\n        type wall;"), "unset faces are walls: {block}");
    // One STL solid and one patch per face, with names OpenFOAM accepts.
    let stl = &case.files["constant/triSurface/Cube.stl"];
    assert_eq!(stl.matches("endsolid").count(), 6);
    assert_eq!(stl.matches("endfacet").count(), 12);
    // Every solid's border (a model face's) is a feature edge the snap follows.
    assert!(case.files["system/surfaceFeatureExtractDict"].contains("geometricTestOnly no;"));
    assert!(stl.contains("solid F1_mx\n") && stl.contains("solid F1_px\n"), "{stl}");
    assert_eq!(case.patches.get("Cube_F1_mz").map(String::as_str), Some("Cube/F1:-z"));
    assert_eq!(case.patches.len(), 12);
    // Kinematic viscosity and pressure, the inlet velocity, the obstacle walls.
    assert!(case.files["constant/transportProperties"].contains("nu 1.5e-5;"));
    assert!(case.files["0/U"].contains("x_min { type fixedValue; value uniform (2e0 0 0); }"), "{}", case.files["0/U"]);
    assert!(case.files["0/U"].contains("Cube_F1_mz { type noSlip; }"));
    assert!(case.files["0/p"].contains("x_max { type fixedValue; value uniform 0; }"));
    // Smagorinsky 0.16 → Ck ≈ 0.0882.
    let les = &case.files["constant/turbulenceProperties"];
    assert!(les.contains("LESModel Smagorinsky") && les.contains("Ck 8.82"), "{les}");
    // The step's length and frames.
    let control = &case.files["system/controlDict"];
    assert!(control.contains("endTime 1e-2;") && control.contains("writeInterval 2e-3;"), "{control}");
    let snappy = &case.files["system/snappyHexMeshDict"];
    assert!(snappy.contains("level (0 0)") && snappy.contains("locationInMesh (9.25e-3 2.25e-3 2.25e-3);"), "{snappy}");
}

#[test]
fn the_same_model_gives_the_same_case() {
    let m = channel(serde_json::json!({}));
    assert_eq!(write(&m, &options()).unwrap(), write(&m, &options()).unwrap());
}

#[test]
fn named_moving_faces_and_refinement_boxes() {
    let m = channel(serde_json::json!({"surfaces": [{"face": {"body": "Cube", "face": "F1:+z"},
        "condition": {"type": "moving_wall", "velocity": [0.5, 0.0, 0.0]}}]}));
    let options = CaseOptions {
        mesh: MeshOptions { cell_factor: 4, refine_boxes: vec![RefineBox { min: [0.0, 0.0, 0.0], max: [0.01, 0.002, 0.002], spacing: None }] },
        processes: 8,
        ..options()
    };
    let case = write(&m, &options).unwrap();
    assert!(case.files["0/U"].contains("Cube_F1_pz { type fixedValue; value uniform (5e-1 0 0); }"), "{}", case.files["0/U"]);
    let snappy = &case.files["system/snappyHexMeshDict"];
    assert!(snappy.contains("refine0 { type searchableBox;") && snappy.contains("levels ((5e-4 2))"), "{snappy}");
    assert!(case.files["system/decomposeParDict"].contains("numberOfSubdomains 8;"));
    // The ranks are decomposeParDict's, written into the script: nothing can run more or fewer.
    let allrun = &case.files["Allrun"];
    assert!(allrun.contains("mpirun -np 8 pimpleFoam -parallel") && !allrun.contains("NP"), "{allrun}");
    let serial = write(&m, &CaseOptions { processes: 1, ..options.clone() }).unwrap();
    assert!(serial.files["Allrun"].contains("\nrun pimpleFoam\n") && !serial.files["Allrun"].contains("mpirun"));
}

#[test]
fn refuses_what_it_cannot_map() {
    let cases = [
        (serde_json::json!({"faces": {"x_min": {"type": "periodic"}, "x_max": {"type": "periodic"}}}), "periodic"),
        (serde_json::json!({"faces": {"x_min": {"type": "velocity", "velocity": [1.0, 0.0, 0.0], "profile": "parabolic"}}}), "parabolic"),
        (serde_json::json!({"acceleration": [0.0, 0.0, -9.81]}), "acceleration"),
    ];
    for (extra, what) in cases {
        let e = write(&channel(extra), &options()).unwrap_err().to_string();
        assert!(e.contains("not supported by OpenFOAM") && e.contains(what), "{e}");
    }
    let e = write(&channel(serde_json::json!({})), &CaseOptions { mesh: MeshOptions { cell_factor: 3, ..MeshOptions::default() }, ..options() })
        .unwrap_err();
    assert!(e.to_string().contains("cell factor"), "{e}");
    // An obstacle with no surface.
    let e = write_case(&channel(serde_json::json!({})), None, &[], &options()).unwrap_err();
    assert!(e.to_string().contains("obstacle `Cube` has no surface"), "{e}");
}

#[test]
fn hostile_surfaces_never_panic() {
    let m = channel(serde_json::json!({}));
    let (p, t, f, n) = cube();
    let bad_index = [[0u32, 1, 99]];
    let short_faces = [0u32; 3];
    let nan = [[f64::NAN; 3]; 8];
    for (pos, tris, faces) in [(&p[..], &bad_index[..], &f[..1]), (&p[..], &t[..], &short_faces[..]), (&nan[..], &t[..], &f[..])] {
        let o = ObstacleSurface { body: "Cube", positions: pos, triangles: tris, tri_face: faces, faces: &n };
        assert!(write_case(&m, None, &[o], &options()).is_err());
    }
}

#[test]
fn every_written_dictionary_follows_the_grammar() {
    let case = write(&channel(serde_json::json!({})), &options()).unwrap();
    for (path, text) in &case.files {
        // Dictionaries and fields only: not the STL, the script or the empty case.foam marker.
        if path.ends_with(".stl") || path == "Allrun" || path == "case.foam" {
            continue;
        }
        let tree = super::grammar::check_file(text).unwrap_or_else(|e| panic!("{path}: {e}\n{text}"));
        assert!(tree.all_directives().is_empty(), "{path}: every setting is written into the case, no #include: {:?}", tree.all_directives());
    }
}

#[test]
fn snappy_writes_every_quality_limit_and_no_layer_block() {
    let case = write(&channel(serde_json::json!({})), &options()).unwrap();
    let d = super::grammar::check_file(&case.files["system/snappyHexMeshDict"]).unwrap();
    // User guide Table 4.8, and minVol 0 (scale-free).
    for (entry, value) in [
        ("maxNonOrtho", 65.0),
        ("maxBoundarySkewness", 20.0),
        ("maxInternalSkewness", 4.0),
        ("maxConcave", 80.0),
        ("minVol", 0.0),
        ("minTetQuality", 1e-13),
        ("minArea", -1.0),
        ("minTwist", 0.05),
        ("minDeterminant", 0.001),
        ("minFaceWeight", 0.05),
        ("minVolRatio", 0.01),
        ("minTriangleTwist", -1.0),
        ("minFaceFlatness", -1.0),
        ("nSmoothScale", 4.0),
        ("errorReduction", 0.75),
    ] {
        assert_eq!(d.number(&format!("meshQualityControls/{entry}")), Some(value), "{entry}");
    }
    assert!(d.get("meshQualityControls/relaxed").is_none(), "relaxed limits are for layers only");
    assert!(d.get("addLayersControls").is_none(), "no layer block while addLayers is false");
    assert_eq!(d.text("addLayers").as_deref(), Some("false"));
    assert_eq!(d.text("singleRegionName").as_deref(), Some("false"), "a one-face body keeps its <body>_<face> patch");
    // 20 x 8 x 8 background cells, nothing refined: the cell limits cannot bind.
    assert_eq!(d.number("castellatedMeshControls/maxGlobalCells"), Some(1280.0));
    assert_eq!(d.number("castellatedMeshControls/maxLocalCells"), Some(1280.0));
    assert_eq!(d.number("castellatedMeshControls/minRefinementCells"), Some(0.0));
    assert_eq!(d.number("castellatedMeshControls/nCellsBetweenLevels"), Some(3.0));
    assert_eq!(d.number("castellatedMeshControls/resolveFeatureAngle"), Some(30.0));
    assert_eq!(d.number("snapControls/tolerance"), Some(2.0));
    assert_eq!(d.number("snapControls/nSolveIter"), Some(30.0));
    assert_eq!(d.number("mergeTolerance"), Some(1e-6));
    let features = super::grammar::check_file(&case.files["system/surfaceFeatureExtractDict"]).unwrap();
    assert_eq!(features.number("Cube.stl/includedAngle"), Some(150.0));
    // The written digits resolve the merge tolerance (10^(1 - writePrecision) below it).
    let control = super::grammar::check_file(&case.files["system/controlDict"]).unwrap();
    let digits = control.number("writePrecision").unwrap();
    assert!(10f64.powf(1.0 - digits) < d.number("mergeTolerance").unwrap());
}

#[test]
fn schemes_and_solution_are_the_cited_ones() {
    let case = write(&channel(serde_json::json!({})), &options()).unwrap();
    let s = super::grammar::check_file(&case.files["system/fvSchemes"]).unwrap();
    assert_eq!(s.text("ddtSchemes/default").as_deref(), Some("backward"));
    assert_eq!(s.text("gradSchemes/default").as_deref(), Some("Gauss linear"));
    assert!(s.get("gradSchemes/grad(U)").is_none(), "no limiter on the strain rate the subgrid model reads");
    assert_eq!(s.text("divSchemes/default").as_deref(), Some("Gauss linear"));
    // The channel is an LES: central convection, adding no dissipation to the subgrid model's.
    assert_eq!(s.text("divSchemes/div(phi,U)").as_deref(), Some("Gauss linear"));
    assert_eq!(s.text("laplacianSchemes/default").as_deref(), Some("Gauss linear limited corrected 0.333"));
    assert_eq!(s.text("snGradSchemes/default").as_deref(), Some("limited corrected 0.333"));
    assert_eq!(s.text("interpolationSchemes/default").as_deref(), Some("linear"));
    assert!(s.get("wallDist").is_none(), "no wall-distance model is written");

    // A laminar step convects by linear upwind, the bounded second-order scheme.
    let mut laminar = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
    laminar["steps"][0]["procedure"]["turbulence"] = serde_json::Value::Null;
    let laminar: Model = serde_json::from_value(laminar).unwrap();
    let l = super::grammar::check_file(&write(&laminar, &options()).unwrap().files["system/fvSchemes"]).unwrap();
    assert_eq!(l.text("divSchemes/div(phi,U)").as_deref(), Some("Gauss linearUpwind grad(U)"));

    let f = super::grammar::check_file(&case.files["system/fvSolution"]).unwrap();
    let tolerance = Some(f64::from(f32::EPSILON));
    assert_eq!(f.text("solvers/p/solver").as_deref(), Some("GAMG"));
    assert_eq!(f.text("solvers/p/smoother").as_deref(), Some("GaussSeidel"));
    assert_eq!(f.number("solvers/p/tolerance"), tolerance);
    assert_eq!(f.number("solvers/p/relTol"), Some(0.0));
    assert!(f.get("solvers/pFinal").is_some());
    assert_eq!(f.text("solvers/(U|UFinal)/smoother").as_deref(), Some("symGaussSeidel"));
    assert_eq!(f.number("solvers/(U|UFinal)/nSweeps"), Some(1.0));
    assert_eq!(f.number("solvers/(U|UFinal)/tolerance"), tolerance);
    assert_eq!(f.number("solvers/(U|UFinal)/relTol"), Some(0.0));
    assert_eq!(f.number("PIMPLE/nOuterCorrectors"), Some(1.0));
    assert_eq!(f.number("PIMPLE/nCorrectors"), Some(2.0));
    assert_eq!(f.number("PIMPLE/nNonOrthogonalCorrectors"), Some(1.0), "a snapped mesh");
    assert_eq!(f.number("PIMPLE/pRefCell"), Some(0.0));
    assert_eq!(f.number("PIMPLE/pRefValue"), Some(0.0));

    // No obstacle and no refinement: an orthogonal mesh needs no non-orthogonal corrector.
    let empty = write_case(&channel(serde_json::json!({"obstacles": []})), None, &[], &options()).unwrap();
    let f = super::grammar::check_file(&empty.files["system/fvSolution"]).unwrap();
    assert_eq!(f.number("PIMPLE/nNonOrthogonalCorrectors"), Some(0.0));
}

#[test]
fn run_controls_are_derived_from_the_step() {
    let case = write(&channel(serde_json::json!({})), &options()).unwrap();
    let c = super::grammar::check_file(&case.files["system/controlDict"]).unwrap();
    assert_eq!(c.text("application").as_deref(), Some("pimpleFoam"));
    assert_eq!(c.text("startFrom").as_deref(), Some("startTime"));
    assert_eq!(c.number("startTime"), Some(0.0));
    assert_eq!(c.text("writeControl").as_deref(), Some("adjustableRunTime"));
    assert_eq!(c.text("adjustTimeStep").as_deref(), Some("yes"));
    assert_eq!(c.number("maxCo"), Some(1.0));
    assert_eq!(c.number("writeInterval"), Some(0.002));
    assert_eq!(c.number("maxDeltaT"), Some(0.002), "no step longer than the time between two frames");
    // Courant 1 on the 0.5 mm cell at the 2 m/s inlet.
    assert_eq!(c.number("deltaT"), Some(0.00025));
    assert_eq!(c.number("writePrecision"), Some(9.0), "nine digits round-trip an f32");
    assert_eq!(c.number("timePrecision"), Some(9.0));
    assert_eq!(c.text("writeFormat").as_deref(), Some("ascii"));
    assert!(c.get("functions").is_none(), "the step does not average");
}

#[test]
fn smagorinsky_and_van_driest_are_mapped() {
    let case = write(&channel(serde_json::json!({})), &options()).unwrap();
    let t = super::grammar::check_file(&case.files["constant/turbulenceProperties"]).unwrap();
    assert_eq!(t.text("simulationType").as_deref(), Some("LES"));
    assert_eq!(t.text("LES/LESModel").as_deref(), Some("Smagorinsky"));
    assert_eq!(t.text("LES/delta").as_deref(), Some("cubeRootVol"));
    assert_eq!(t.number("LES/cubeRootVolCoeffs/deltaCoeff"), Some(1.0));
    let (ck, ce) = (t.number("LES/SmagorinskyCoeffs/Ck").unwrap(), t.number("LES/SmagorinskyCoeffs/Ce").unwrap());
    assert_eq!(ce, 1.048);
    // Cs² = Ck √(Ck/Ce) gives back the step's constant, 0.16.
    assert!(((ck * (ck / ce).sqrt()).sqrt() - 0.16).abs() < 1e-12, "Ck {ck}");
    assert!((ck - 0.08823).abs() < 1e-5, "Ck {ck}");

    // Van Driest's damping: the step's A+ and its constant as Cdelta, with OpenFOAM's law.
    let mut model = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
    model["steps"][0]["procedure"]["turbulence"]["van_driest"] = serde_json::json!(26.0);
    let model: Model = serde_json::from_value(model).unwrap();
    let damped = write(&model, &options()).unwrap();
    let v = super::grammar::check_file(&damped.files["constant/turbulenceProperties"]).unwrap();
    assert_eq!(v.text("LES/delta").as_deref(), Some("vanDriest"));
    assert_eq!(v.number("LES/vanDriestCoeffs/Aplus"), Some(26.0));
    assert_eq!(v.number("LES/vanDriestCoeffs/kappa"), Some(0.41));
    assert_eq!(v.number("LES/vanDriestCoeffs/Cdelta"), Some(0.16));
    let f = super::grammar::check_file(&damped.files["system/fvSchemes"]).unwrap();
    assert_eq!(f.text("wallDist/method").as_deref(), Some("meshWave"));
    assert!(damped.notes.iter().any(|n| n.contains("Van Driest")), "{:?}", damped.notes);
}

/// A one-triangle surface of a body with the given face name.
fn tri(face: &str) -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<u32>, Vec<String>) {
    (vec![[0.004, 0.001, 0.001], [0.006, 0.001, 0.001], [0.004, 0.003, 0.001]], vec![[0, 1, 2]], vec![0], vec![face.to_string()])
}

/// A model whose fluid domain lists these obstacles, written from one-triangle surfaces named by
/// `(body, face)`.
fn write_named(bodies: &[(&'static str, &str)], boxes: Vec<RefineBox>) -> Case {
    let mut model = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
    model["fluid_domains"]["Flow"]["obstacles"] = serde_json::json!(bodies.iter().map(|(b, _)| *b).collect::<Vec<_>>());
    let model: Model = serde_json::from_value(model).unwrap();
    let surfaces: Vec<_> = bodies.iter().map(|(b, f)| (*b, tri(f))).collect();
    let obstacles: Vec<ObstacleSurface<'_>> =
        surfaces.iter().map(|(b, (p, t, f, n))| ObstacleSurface { body: b, positions: p, triangles: t, tri_face: f, faces: n }).collect();
    let options = CaseOptions { mesh: MeshOptions { cell_factor: 1, refine_boxes: boxes }, ..options() };
    write_case(&model, None, &obstacles, &options).unwrap()
}

/// Every patch of `0/U`'s `boundaryField` exactly once, the case's patches as many, and the
/// snappy geometry names all different.
fn assert_names_are_unique(case: &Case, patches: usize) {
    assert_eq!(case.patches.len(), patches, "{:?}", case.patches);
    let u = super::grammar::check_file(&case.files["0/U"]).unwrap();
    let Some(super::grammar::Entry::Dict(field)) = u.get("boundaryField") else { panic!("no boundaryField") };
    let mut names: Vec<&str> = field.entries.iter().map(|(k, _)| k.as_str()).collect();
    names.sort_unstable();
    let listed = names.len();
    names.dedup();
    assert_eq!(listed, names.len(), "a patch listed twice in 0/U: {names:?}");
    assert_eq!(names, case.patches.keys().map(String::as_str).collect::<Vec<_>>());
    let snappy = super::grammar::check_file(&case.files["system/snappyHexMeshDict"]).unwrap();
    let Some(super::grammar::Entry::Dict(geometry)) = snappy.get("geometry") else { panic!("no geometry") };
    let mut keys: Vec<&str> = geometry.entries.iter().map(|(k, _)| k.as_str()).collect();
    let all = keys.len();
    keys.sort_unstable();
    keys.dedup();
    assert_eq!(all, keys.len(), "a geometry listed twice: {keys:?}");
}

#[test]
fn a_body_patch_never_takes_a_box_patch_or_another_bodys_name() {
    // `x` with a face `min` would be `x_min`, a box patch.
    assert_names_are_unique(&write_named(&[("x", "min")], Vec::new()), 7);
    // `a_b` with `c` and `a` with `b_c` both give `a_b_c`.
    assert_names_are_unique(&write_named(&[("a_b", "c"), ("a", "b_c")], Vec::new()), 8);
    // A body named `refine0` against the first refinement box's geometry.
    let boxes = vec![RefineBox { min: [0.0; 3], max: [0.002; 3], spacing: None }];
    assert_names_are_unique(&write_named(&[("refine0", "f")], boxes), 7);
}

#[test]
fn the_finest_centre_follows_the_blocks_cells() {
    // An x extent n h (1 + 5e-7), within what `cells` takes as whole.
    let mut model = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
    model["fluid_domains"]["Flow"]["max"] = serde_json::json!([0.01 * (1.0 + 5e-7), 0.004, 0.004]);
    let model: Model = serde_json::from_value(model).unwrap();
    let domain = model.fluid_domains().get("Flow").unwrap();
    let mesh = MeshOptions { cell_factor: 1, refine_boxes: vec![RefineBox { min: [0.0; 3], max: [0.002; 3], spacing: Some(0.0005 / 8.0) }] };
    let k = [1_000_000_000u64, 0, 0];
    let c = finest_centre(domain, &mesh, k).unwrap();
    // 20 background cells, each halved three times.
    let want = domain.min[0] + (k[0] as f64 + 0.5) * (domain.max[0] - domain.min[0]) / (20.0 * 8.0);
    assert!(((c[0] - want) / want).abs() < 1e-12, "{} against {want}", c[0]);
}

#[test]
fn a_facet_of_text_is_within_the_size_the_cap_assumes() {
    let (p, t, f, n) = tri("F");
    let o = ObstacleSurface { body: "Cube", positions: &p, triangles: &t, tri_face: &f, faces: &n };
    let (text, _) = stl(&o).unwrap();
    let solid_lines = "solid F\nendsolid F\n".len();
    assert!(text.len() - solid_lines <= FACET_TEXT_BYTES, "{} bytes", text.len() - solid_lines);
}

#[test]
fn refuses_each_case_it_cannot_write() {
    let patch_material = |patch: serde_json::Value| {
        let mut model = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
        if let (Some(m), Some(p)) = (model["materials"]["Air"].as_object_mut(), patch.as_object()) {
            for (k, v) in p {
                if v.is_null() {
                    m.remove(k);
                } else {
                    m.insert(k.clone(), v.clone());
                }
            }
        }
        serde_json::from_value::<Model>(model).unwrap()
    };
    let with_turbulence = |t: serde_json::Value| {
        let mut model = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
        model["steps"][0]["procedure"]["turbulence"] = t;
        serde_json::from_value::<Model>(model).unwrap()
    };
    let zone = serde_json::json!({"name": "Z", "min": [0.001, 0.001, 0.001], "max": [0.002, 0.002, 0.002]});
    let table = serde_json::json!({"viscosity": [[250.0, 1.7e-5], [350.0, 2.1e-5]]});
    let mut fine_box = options();
    fine_box.mesh.refine_boxes = vec![RefineBox { min: [0.0; 3], max: [0.01, 0.004, 0.004], spacing: Some(0.0005 / 1024.0) }];
    let mut long = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
    long["fluid_domains"]["Flow"]["max"] = serde_json::json!([1.0, 0.004, 0.004]);
    let long: Model = serde_json::from_value(long).unwrap();
    let mut heat = serde_json::to_value(channel(serde_json::json!({}))).unwrap();
    heat["fluid_domains"]["Flow"]["thermal"] = serde_json::json!({});
    heat["materials"]["Air"]["thermal"] = serde_json::json!({"conductivity": 0.026, "specific_heat": 1005.0});
    let heat: Model = serde_json::from_value(heat).unwrap();
    let cases: Vec<(&str, Model, CaseOptions, &str)> = vec![
        ("merge tolerance", long, fine_box, "merge tolerance"),
        ("Courant above the limit", channel(serde_json::json!({})), CaseOptions { max_courant: 1.5, ..options() }, "Courant"),
        ("Courant zero", channel(serde_json::json!({})), CaseOptions { max_courant: 0.0, ..options() }, "Courant"),
        ("no processes", channel(serde_json::json!({})), CaseOptions { processes: 0, ..options() }, "0 processes"),
        ("location outside", channel(serde_json::json!({})), CaseOptions { location: [1.0, 0.002, 0.002], ..options() }, "location"),
        ("heat", heat, options(), "heat"),
        (
            "porous zone",
            channel(serde_json::json!({"porous": [{"name": "Z", "min": [0.001, 0.001, 0.001], "max": [0.002, 0.002, 0.002],
                "axis": "x", "resistance": {"type": "loss_coefficient", "coefficient": 1.0}}]})),
            options(),
            "porous",
        ),
        (
            "fan",
            channel(serde_json::json!({"fans": [{"name": "F", "min": zone["min"], "max": zone["max"], "direction": "+x",
                "curve": [[0.0, 10.0], [1.0, 0.0]]}]})),
            options(),
            "fans",
        ),
        (
            "baffle",
            channel(
                serde_json::json!({"faces": {"y_min": {"type": "symmetry"}}, "baffles": [{"name": "B", "axis": "x", "position": 0.005, "conductance": 1.0}]}),
            ),
            options(),
            "baffles",
        ),
        (
            "free surface",
            channel(serde_json::json!({"free_surface": {"liquid": [{"min": [0.0, 0.0, 0.0], "max": [0.005, 0.004, 0.004]}]}})),
            options(),
            "free surface",
        ),
        ("Vreman", with_turbulence(serde_json::json!({"type": "vreman", "constant": 0.17})), options(), "Vreman"),
        ("Spalart-Allmaras", with_turbulence(serde_json::json!({"type": "spalart_allmaras"})), options(), "Spalart"),
        (
            "power law",
            patch_material(serde_json::json!({"fluid": {"type": "power_law", "consistency": 1e-3, "index": 0.8}})),
            options(),
            "power-law viscosity",
        ),
        (
            "temperature-dependent viscosity",
            patch_material(serde_json::json!({"temperature_dependent": table})),
            options(),
            "viscosity that follows the temperature",
        ),
        (
            "temperature-dependent density",
            patch_material(serde_json::json!({"temperature_dependent": {"density": [[250.0, 1.3], [350.0, 1.0]]}})),
            options(),
            "density that follows the temperature",
        ),
    ];
    for (what, model, options, needle) in cases {
        let e = write(&model, &options).expect_err(what).to_string();
        assert!(e.contains(needle), "{what}: {e}");
    }
}
