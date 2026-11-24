//! #785 / #1104: a subgrid model on halfway obstacle walls warns where a connected oblique region of the surface
//! rises at least one lattice cell along its axis. A curved body covers ±45° about its axis, so it rises 0.29 R
//! along it: a smooth sphere warns from R ≥ Δx / (1 − cos 45°) ≈ 3.4 Δx, about 7 cells across; a faceted one a
//! little earlier, since its regions reach past 45°.
mod common;
use common::*;
use fieldcraft_lbm::{Limits, Setup};
use fieldcraft_mesh::surface::BoundarySurface;
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile, Turbulence, WallModel, WallTreatment};

/// The solve's setup for a channel with one obstacle named `Body`, with the given turbulence and walls.
fn warnings_of(body: BoundarySurface, turbulence: Option<Turbulence>, walls: WallTreatment) -> Vec<String> {
    let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], 0.05);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.obstacles = vec!["Body".into()];
    let mut p = collide(1.0, 0.005, None, Collision::Regularized);
    if let Procedure::Fluid { turbulence: t, walls: w, .. } = &mut p {
        *t = turbulence;
        *w = walls;
    }
    let m = model(1.0, 1e-3, d, p);
    let obstacles = [obstacle("Body", body, &["wall"])];
    let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
    setup.warnings
}

fn smagorinsky(wall_model: Option<WallModel>) -> Option<Turbulence> {
    Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model })
}

/// A prism along x from `x0` to `x1` over the convex cross-section `section`, its (y, z) corners in order: the side
/// walls and two caps, each cap a fan from its first corner.
fn prism(x0: f64, x1: f64, section: &[[f64; 2]]) -> BoundarySurface {
    let n = section.len() as u32;
    let positions: Vec<[f64; 3]> = [x0, x1].iter().flat_map(|&x| section.iter().map(move |p| [x, p[0], p[1]])).collect();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        triangles.push([i, j, n + j]);
        triangles.push([i, n + j, n + i]);
    }
    for k in 1..n - 1 {
        triangles.push([0, k, k + 1]);
        triangles.push([n, n + k + 1, n + k]);
    }
    let count = section.len() as f64;
    let centre = [0.5 * (x0 + x1), section.iter().map(|p| p[0]).sum::<f64>() / count, section.iter().map(|p| p[1]).sum::<f64>() / count];
    let tags = vec![0; triangles.len()];
    surface(positions, triangles, tags, move |c, _| [c[0] - centre[0], c[1] - centre[1], c[2] - centre[2]])
}

#[test]
fn a_subgrid_model_on_a_staircase_warns_once_and_only_there() {
    // Δx is 0.05: radii of 3.6 Δx (warns) and 2 Δx (does not), clear of the smooth threshold 3.4 Δx.
    let sphere_body = || sphere([0.7, 0.5, 0.5], 0.18, 3);
    let small_sphere = || sphere([0.7, 0.5, 0.5], 0.10, 3);
    let block_body = || cuboid([0.5, 0.3, 0.3], [0.9, 0.7, 0.7]);
    let vreman = Some(Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: None });
    // Oblique, with edges 0.010 and 0.011 long: a feature shorter than one cell (spacing 0.05) rises by less than a cell.
    let sub_cell = || parallelepiped([0.5, 0.3, 0.3], [[0.01, 0.0, 0.0], [0.004, 0.01, 0.0], [0.0, 0.004, 0.01]]);
    // A large box, x from 0.3 to 1.5, with one edge along x chamfered 0.012 across (sub-cell).
    let chamfer = || prism(0.3, 1.5, &[[0.2, 0.2], [0.8, 0.2], [0.8, 0.788], [0.788, 0.8], [0.2, 0.8]]);
    // The same box with two such chamfers at opposite edges, 0.6 apart: each is far more than a cell from the other.
    let two_chamfers = || prism(0.3, 1.5, &[[0.212, 0.2], [0.8, 0.2], [0.8, 0.788], [0.788, 0.8], [0.2, 0.8], [0.2, 0.212]]);
    // A box 0.8 long along x and tilted 10 degrees about y: its bottom face rises 0.8 sin 10° = 0.14 along z, 2.8 cells.
    let tilted = || {
        let (c, s) = (10f64.to_radians().cos(), 10f64.to_radians().sin());
        parallelepiped([0.6, 0.3, 0.3], [[0.8 * c, 0.0, 0.8 * s], [-0.4 * s, 0.0, 0.4 * c], [0.0, 0.4, 0.0]])
    };
    let cases: [(&str, BoundarySurface, Option<Turbulence>, WallTreatment, bool); 11] = [
        ("sphere, Smagorinsky, halfway", sphere_body(), smagorinsky(None), WallTreatment::Halfway, true),
        ("sphere of radius 2 Δx, Smagorinsky, halfway", small_sphere(), smagorinsky(None), WallTreatment::Halfway, false),
        ("sphere, Vreman, halfway", sphere_body(), vreman, WallTreatment::Halfway, true),
        ("sphere, Smagorinsky, interpolated", sphere_body(), smagorinsky(None), WallTreatment::Interpolated, false),
        ("sphere, no model, halfway", sphere_body(), None, WallTreatment::Halfway, false),
        ("sphere, Smagorinsky with a wall model, halfway", sphere_body(), smagorinsky(Some(WallModel::Musker)), WallTreatment::Halfway, false),
        ("cuboid, Smagorinsky, halfway", block_body(), smagorinsky(None), WallTreatment::Halfway, false),
        ("sub-cell oblique faces, Smagorinsky, halfway", sub_cell(), smagorinsky(None), WallTreatment::Halfway, false),
        ("large box with a sub-cell chamfer, Smagorinsky, halfway", chamfer(), smagorinsky(None), WallTreatment::Halfway, false),
        ("two sub-cell chamfers, Smagorinsky, halfway", two_chamfers(), smagorinsky(None), WallTreatment::Halfway, false),
        ("box tilted 10 degrees, Smagorinsky, halfway", tilted(), smagorinsky(None), WallTreatment::Halfway, true),
    ];
    for (name, body, turbulence, walls, hit) in cases {
        let warnings = warnings_of(body, turbulence, walls);
        let staircase: Vec<&String> = warnings.iter().filter(|w| w.contains("staircase")).collect();
        if hit {
            assert_eq!(staircase.len(), 1, "{name}: {warnings:?}");
            assert!(staircase[0].contains("`Body`"), "{name}: {}", staircase[0]);
        } else {
            assert!(staircase.is_empty(), "{name}: {warnings:?}");
        }
    }
}
