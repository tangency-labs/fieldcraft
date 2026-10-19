use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::*;
use crate::quality::Validity;
use crate::size::{UniformSize, divide_polyline, even_divisions, ideal_divisions, loop_divisions};

fn opts() -> QuadOptions {
    QuadOptions::default()
}

/// A `w × h` rectangle at `(x0, y0)`, counterclockwise, `nx` divisions along x and `ny` along y.
fn rect(x0: f64, y0: f64, w: f64, h: f64, nx: usize, ny: usize) -> Vec<[f64; 2]> {
    let mut p = Vec::new();
    for i in 0..nx {
        p.push([x0 + w * i as f64 / nx as f64, y0]);
    }
    for j in 0..ny {
        p.push([x0 + w, y0 + h * j as f64 / ny as f64]);
    }
    for i in 0..nx {
        p.push([x0 + w - w * i as f64 / nx as f64, y0 + h]);
    }
    for j in 0..ny {
        p.push([x0, y0 + h - h * j as f64 / ny as f64]);
    }
    p
}

fn circle(cx: f64, cy: f64, r: f64, n: usize) -> Vec<[f64; 2]> {
    (0..n).map(|k| 2.0 * PI * k as f64 / n as f64).map(|t| [cx + r * t.cos(), cy + r * t.sin()]).collect()
}

/// A closed polygon with every side divided evenly into `per_side` pieces.
fn subdivide(corners: &[[f64; 2]], per_side: usize) -> Vec<[f64; 2]> {
    let mut p = Vec::new();
    for (i, a) in corners.iter().enumerate() {
        let b = corners[(i + 1) % corners.len()];
        for k in 0..per_side {
            let t = k as f64 / per_side as f64;
            p.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]);
        }
    }
    p
}

fn signed_area(p: &[[f64; 2]]) -> f64 {
    (0..p.len()).map(|i| p[i][0] * p[(i + 1) % p.len()][1] - p[(i + 1) % p.len()][0] * p[i][1]).sum::<f64>() / 2.0
}

/// Every mesh invariant: boundary nodes unchanged and first, every element valid by the exact
/// test and counterclockwise, conformity (each interior edge used once in each direction, each
/// boundary edge once), the boundary edges exactly the input loops, and the area.
fn check(mesh: &QuadMesh, loops: &[Vec<[f64; 2]>]) {
    let input: Vec<[f64; 2]> = loops.iter().flatten().copied().collect();
    assert_eq!(mesh.boundary_nodes, input.len());
    assert_eq!(&mesh.nodes[..input.len()], input.as_slice());
    let mut area = 0.0;
    let mut directed: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for q in &mesh.quads {
        let c = q.map(|v| mesh.nodes[v as usize]);
        assert_eq!(crate::quality::planar_quad_validity(&c).unwrap(), Validity::Valid, "{c:?}");
        area += signed_area(&c);
        for k in 0..4 {
            *directed.entry((q[k], q[(k + 1) % 4])).or_default() += 1;
        }
    }
    assert!(directed.values().all(|&n| n == 1), "an edge is used twice in the same direction");
    let boundary: Vec<(u32, u32)> = directed.keys().filter(|(a, b)| !directed.contains_key(&(*b, *a))).copied().collect();
    let mut expected = Vec::new();
    let mut base = 0u32;
    for l in loops {
        let n = l.len() as u32;
        let ccw = signed_area(l) > 0.0;
        let outer = base == 0;
        for i in 0..n {
            let (a, b) = (base + i, base + (i + 1) % n);
            // The outer loop runs counterclockwise around the mesh, holes clockwise.
            expected.push(if ccw == outer { (a, b) } else { (b, a) });
        }
        base += n;
    }
    expected.sort_unstable();
    assert_eq!(boundary, expected, "the mesh boundary is not the input loops");
    let want: f64 = loops.iter().enumerate().map(|(k, l)| if k == 0 { signed_area(l).abs() } else { -signed_area(l).abs() }).sum();
    assert!((area - want).abs() <= 1e-9 * want.abs().max(1.0), "area {area} vs {want}");
}

#[test]
fn rectangle_free() {
    let loops = vec![rect(0.0, 0.0, 10.0, 5.0, 20, 10)];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
    assert!(m.stats.min_scaled_jacobian > 0.3, "{:?}", m.stats);
    // About (10 × 5) / 0.5² = 200 quads.
    assert!((100..400).contains(&m.quads.len()), "{}", m.quads.len());
}

#[test]
fn single_quad_and_a_six_point_loop() {
    let quad = vec![vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]];
    let m = free_quad_mesh(&quad, &UniformSize(1.0), &opts()).unwrap();
    assert_eq!(m.quads.len(), 1);
    check(&m, &quad);
    let six = vec![rect(0.0, 0.0, 2.0, 1.0, 2, 1)];
    let m = free_quad_mesh(&six, &UniformSize(1.0), &opts()).unwrap();
    check(&m, &six);
}

#[test]
fn clockwise_input_is_accepted() {
    let mut l = rect(0.0, 0.0, 4.0, 2.0, 8, 4);
    l.reverse();
    let loops = vec![l];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
}

#[test]
fn l_shape() {
    let corners: [[f64; 2]; 6] = [[0.0, 0.0], [10.0, 0.0], [10.0, 4.0], [4.0, 4.0], [4.0, 10.0], [0.0, 10.0]];
    // Even divisions per side keep the corners at even positions.
    let mut l = Vec::new();
    for (i, a) in corners.iter().enumerate() {
        let b = corners[(i + 1) % corners.len()];
        let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
        let n = even_divisions(len / 0.5).unwrap();
        for k in 0..n {
            let t = k as f64 / n as f64;
            l.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]);
        }
    }
    let loops = vec![l];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
    assert!(m.stats.min_scaled_jacobian > 0.2, "{:?}", m.stats);
}

#[test]
fn plate_with_a_hole() {
    let loops = vec![rect(0.0, 0.0, 20.0, 10.0, 40, 20), circle(10.0, 5.0, 2.5, 32)];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
    assert!(m.stats.min_scaled_jacobian > 0.2, "{:?}", m.stats);
}

#[test]
fn plate_with_many_holes() {
    let mut loops = vec![rect(0.0, 0.0, 30.0, 30.0, 30, 30)];
    for i in 0..3 {
        for j in 0..3 {
            loops.push(circle(5.0 + 10.0 * i as f64, 5.0 + 10.0 * j as f64, 2.0 + 0.5 * ((i + j) % 2) as f64, 16 + 4 * i));
        }
    }
    // A square hole too, given clockwise.
    loops[5] = rect(13.0, 13.0, 4.0, 4.0, 4, 4).into_iter().rev().collect();
    let m = free_quad_mesh(&loops, &UniformSize(1.0), &opts()).unwrap();
    check(&m, &loops);
}

#[test]
fn thin_sliver() {
    let loops = vec![rect(0.0, 0.0, 100.0, 1.0, 200, 2)];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
    assert_eq!(m.quads.len(), 400);
    // A sliver one division across: the split puts a mid node on each short side.
    let loops = vec![rect(0.0, 0.0, 50.0, 0.25, 100, 1)];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
}

fn star(r_out: f64, r_in: f64) -> Vec<[f64; 2]> {
    (0..10).map(|k| (k as f64 * PI / 5.0 + PI / 2.0, if k % 2 == 0 { r_out } else { r_in })).map(|(t, r)| [r * t.cos(), r * t.sin()]).collect()
}

#[test]
fn concave_star() {
    let loops = vec![subdivide(&star(10.0, 4.0), 4)];
    let m = free_quad_mesh(&loops, &UniformSize(1.0), &opts()).unwrap();
    check(&m, &loops);
    // A star with a star-shaped hole.
    let loops = vec![subdivide(&star(10.0, 4.0), 4), subdivide(&star(2.5, 1.0), 2)];
    let m = free_quad_mesh(&loops, &UniformSize(1.0), &opts()).unwrap();
    check(&m, &loops);
}

/// A polygon with every side divided evenly at about `h` (even counts, so corners stay at even
/// positions).
fn divided(corners: &[[f64; 2]], h: f64) -> Vec<[f64; 2]> {
    let mut l = Vec::new();
    for (i, a) in corners.iter().enumerate() {
        let b = corners[(i + 1) % corners.len()];
        let n = even_divisions(((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt() / h).unwrap();
        l.extend((0..n).map(|k| k as f64 / n as f64).map(|t| [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]));
    }
    l
}

#[test]
fn i_beam_end_face() {
    // Flanges 10 × 1, web 1 thick, 12 deep overall: thin walls meeting at reflex corners.
    let c = [
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 1.0],
        [5.5, 1.0],
        [5.5, 11.0],
        [10.0, 11.0],
        [10.0, 12.0],
        [0.0, 12.0],
        [0.0, 11.0],
        [4.5, 11.0],
        [4.5, 1.0],
        [0.0, 1.0],
    ];
    for h in [0.5, 0.25] {
        let loops = vec![divided(&c, h)];
        let m = free_quad_mesh(&loops, &UniformSize(h), &opts()).unwrap();
        check(&m, &loops);
        assert!(m.stats.min_scaled_jacobian > 0.2, "{h}: {:?}", m.stats);
    }
}

#[test]
fn perforated_plate() {
    // A 40 × 20 plate with a 7 × 3 grid of holes of radius 1 on a 5 mm pitch.
    let h = 0.5;
    let mut loops = vec![divided(&[[0.0, 0.0], [40.0, 0.0], [40.0, 20.0], [0.0, 20.0]], h)];
    for i in 0..7 {
        for j in 0..3 {
            loops.push(circle(5.0 + 5.0 * i as f64, 5.0 + 5.0 * j as f64, 1.0, 2 * (PI / h).ceil() as usize));
        }
    }
    let m = free_quad_mesh(&loops, &UniformSize(h), &opts()).unwrap();
    check(&m, &loops);
    assert!(m.stats.min_scaled_jacobian > 0.2, "{:?}", m.stats);
}

#[test]
fn comb_with_deep_slots() {
    // Teeth of width 1 between slots of width 1, 4 deep: strongly concave.
    let mut corners = vec![[0.0, 0.0], [9.0, 0.0]];
    for k in (0..5).rev() {
        let x = 2.0 * k as f64;
        corners.extend([[x + 1.0, 6.0], [x, 6.0]]);
        if k > 0 {
            corners.extend([[x, 2.0], [x - 1.0, 2.0]]);
        }
    }
    let loops = vec![subdivide(&corners, 4)];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &opts()).unwrap();
    check(&m, &loops);
}

#[test]
fn graded_size_field() {
    let loops = vec![rect(0.0, 0.0, 20.0, 20.0, 20, 20)];
    let near = |p: [f64; 3]| 0.25 + 0.15 * (p[0] * p[0] + p[1] * p[1]).sqrt();
    let graded = free_quad_mesh(&loops, &near, &opts()).unwrap();
    check(&graded, &loops);
    let uniform = free_quad_mesh(&loops, &UniformSize(1.0), &opts()).unwrap();
    // Denser than the uniform mesh near the origin, coarser far from it.
    let within = |m: &QuadMesh, r: f64| m.nodes.iter().filter(|p| p[0] * p[0] + p[1] * p[1] < r * r).count();
    assert!(within(&graded, 8.0) > within(&uniform, 8.0), "{} {}", within(&graded, 8.0), within(&uniform, 8.0));
    assert!(graded.nodes.len() - within(&graded, 16.0) < uniform.nodes.len() - within(&uniform, 16.0));
}

#[test]
fn boundary_from_the_size_helpers() {
    // An arc-and-line boundary whose divisions come from the size field, loop parity fixed.
    let field = |p: [f64; 3]| 0.3 + 0.05 * p[0].abs();
    let arc: Vec<[f64; 3]> = (0..=32).map(|k| PI * k as f64 / 32.0).map(|t| [5.0 * t.cos(), 5.0 * t.sin(), 0.0]).collect();
    let base = vec![[-5.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
    let curves = [arc, base];
    let ideal: Vec<f64> = curves.iter().map(|c| ideal_divisions(&field, c).unwrap()).collect();
    let counts = loop_divisions(&ideal).unwrap();
    let mut l = Vec::new();
    for (c, n) in curves.iter().zip(counts) {
        let pts = divide_polyline(&field, c, n).unwrap();
        l.extend(pts[..pts.len() - 1].iter().map(|p| [p[0], p[1]]));
    }
    assert_eq!(l.len() % 2, 0);
    let loops = vec![l];
    let m = free_quad_mesh(&loops, &field, &opts()).unwrap();
    check(&m, &loops);
}

#[test]
fn deterministic() {
    let loops = vec![rect(0.0, 0.0, 20.0, 10.0, 40, 20), circle(7.0, 5.0, 2.0, 24), circle(14.0, 5.0, 3.0, 30)];
    let a = free_quad_mesh(&loops, &UniformSize(0.6), &opts()).unwrap();
    let b = free_quad_mesh(&loops, &UniformSize(0.6), &opts()).unwrap();
    assert_eq!(a, b);
}

#[test]
fn no_smoothing_is_still_valid() {
    let loops = vec![rect(0.0, 0.0, 20.0, 10.0, 40, 20), circle(10.0, 5.0, 2.5, 32)];
    let m = free_quad_mesh(&loops, &UniformSize(0.5), &QuadOptions { smoothing_sweeps: 0, ..opts() }).unwrap();
    check(&m, &loops);
    assert_eq!(m.stats.smoothing_sweeps, 0);
}

#[test]
fn hostile_free_input() {
    let s = UniformSize(1.0);
    let ok = rect(0.0, 0.0, 4.0, 4.0, 4, 4);
    let invalid = |loops: Vec<Vec<[f64; 2]>>| matches!(free_quad_mesh(&loops, &s, &opts()), Err(ref e) if e.kind() == Some(crate::Kind::Invalid));
    assert!(invalid(vec![]));
    assert!(invalid(vec![vec![[0.0, 0.0], [1.0, 0.0]]]));
    // Odd count: an error, never a silent repair.
    let mut odd = ok.clone();
    odd.push([0.0, 0.5]);
    assert!(invalid(vec![odd]));
    let mut nan = ok.clone();
    nan[3] = [f64::NAN, 1.0];
    assert!(invalid(vec![nan]));
    let mut inf = ok.clone();
    inf[3] = [f64::INFINITY, 1.0];
    assert!(invalid(vec![inf]));
    let mut far = ok.clone();
    far[3] = [1e300, 1.0];
    assert!(invalid(vec![far]));
    // Zero-length edge and repeated point.
    let mut zero = ok.clone();
    zero[2] = zero[1];
    assert!(invalid(vec![zero]));
    let mut repeated = ok.clone();
    repeated[10] = repeated[2];
    assert!(invalid(vec![repeated]));
    // Self-intersecting (bow tie) loop.
    assert!(invalid(vec![vec![[0.0, 0.0], [2.0, 2.0], [2.0, 1.0], [2.0, 0.0], [0.0, 2.0], [0.0, 1.0]]]));
    // A point touching another edge.
    assert!(invalid(vec![vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [2.0, 0.0], [0.0, 4.0], [0.0, 2.0]]]));
    // Collinear fold-back.
    assert!(invalid(vec![vec![[0.0, 0.0], [4.0, 0.0], [2.0, 0.0], [2.0, 2.0], [1.0, 2.0], [0.0, 2.0]]]));
    // All collinear.
    assert!(invalid(vec![vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0], [5.0, 0.0]]]));
    // Holes: outside, crossing the outer loop, nested, too few points.
    assert!(invalid(vec![ok.clone(), circle(10.0, 10.0, 1.0, 8)]));
    assert!(invalid(vec![ok.clone(), circle(4.0, 2.0, 1.0, 8)]));
    assert!(invalid(vec![rect(0.0, 0.0, 10.0, 10.0, 10, 10), circle(5.0, 5.0, 3.0, 16), circle(5.0, 5.0, 1.0, 8)]));
    assert!(invalid(vec![ok.clone(), rect(1.0, 1.0, 2.0, 2.0, 1, 1)]));
    // Two holes crossing each other.
    assert!(invalid(vec![rect(0.0, 0.0, 10.0, 10.0, 10, 10), circle(4.0, 5.0, 2.0, 16), circle(6.0, 5.0, 2.0, 16)]));
    // Bad size fields.
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(free_quad_mesh(std::slice::from_ref(&ok), &UniformSize(h), &opts()), Err(ref e) if e.kind() == Some(crate::Kind::Invalid)),
            "{h}"
        );
    }
    // Limits.
    assert!(matches!(free_quad_mesh(std::slice::from_ref(&ok), &UniformSize(1e-4), &opts()), Err(ref e) if e.kind() == Some(crate::Kind::Limit)));
    assert!(
        matches!(free_quad_mesh(&[ok], &UniformSize(0.1), &QuadOptions { max_elements: 10, ..opts() }), Err(ref e) if e.kind() == Some(crate::Kind::Limit))
    );
}

#[test]
fn random_boundaries_never_panic() {
    // Pseudo-random jagged loops (a fixed linear congruential sequence), every point a corner:
    // every result is a valid mesh or an input error.
    let mut seed = 12345u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) as f64) / (1u64 << 31) as f64
    };
    let (mut meshed, mut rejected) = (0, 0);
    for case in 0..200 {
        let n = 6 + 2 * (case % 10);
        let mut loops = vec![
            (0..n)
                .map(|k| {
                    let t = 2.0 * PI * k as f64 / n as f64;
                    let r = 1.0 + 3.0 * next();
                    [r * t.cos(), r * t.sin()]
                })
                .collect::<Vec<_>>(),
        ];
        if case % 3 == 0 {
            loops.push(circle(0.0, 0.0, 0.5 * next(), 6));
        }
        if case % 7 == 0 {
            loops[0] = (0..n).map(|_| [10.0 * next() - 5.0, 10.0 * next() - 5.0]).collect();
        }
        match free_quad_mesh(&loops, &UniformSize(0.2 + next()), &opts()) {
            Ok(m) => {
                check(&m, &loops);
                meshed += 1;
            }
            Err(ref e) if e.kind() == Some(crate::Kind::Invalid) => rejected += 1,
            Err(e) => panic!("case {case}: {e}"),
        }
    }
    // Every valid region is meshed; the rest are self-intersecting loops.
    assert_eq!(meshed + rejected, 200);
    assert!(meshed > 150, "{meshed}");
}

#[test]
fn mapped_rectangle_is_the_grid() {
    let l = rect(0.0, 0.0, 4.0, 2.0, 4, 2);
    let m = mapped_quad_mesh(&l, [0, 4, 6, 10], &opts()).unwrap();
    assert_eq!((m.nu, m.nv), (4, 2));
    assert_eq!(m.mesh.quads.len(), 8);
    check(&m.mesh, std::slice::from_ref(&l));
    for j in 0..=2 {
        for i in 0..=4 {
            assert_eq!(m.mesh.nodes[m.grid[j * 5 + i] as usize], [i as f64, j as f64]);
        }
    }
    assert!((m.mesh.stats.min_scaled_jacobian - 1.0).abs() < 1e-12);
    assert!(m.mesh.stats.max_skew < 1e-12);
}

#[test]
fn mapped_curved_and_graded() {
    // A quarter annulus, radii 1 and 2, graded radially.
    let (na, nr) = (8, 5);
    let r = |j: usize| 1.0 + (j as f64 / nr as f64).powi(2);
    let mut l = Vec::new();
    for j in 0..nr {
        l.push([r(j), 0.0]);
    }
    for i in 0..na {
        let t = 0.5 * PI * i as f64 / na as f64;
        l.push([2.0 * t.cos(), 2.0 * t.sin()]);
    }
    for j in 0..nr {
        l.push([0.0, r(nr - j)]);
    }
    for i in 0..na {
        let t = 0.5 * PI * (na - i) as f64 / na as f64;
        l.push([t.cos(), t.sin()]);
    }
    let corners = [0, nr, nr + na, 2 * nr + na];
    let m = mapped_quad_mesh(&l, corners, &opts()).unwrap();
    check(&m.mesh, std::slice::from_ref(&l));
    // The same region given clockwise.
    let mut cw = l.clone();
    cw.reverse();
    let n = cw.len();
    let m = mapped_quad_mesh(&cw, [n - 1 - corners[3], n - 1 - corners[2], n - 1 - corners[1], n - 1 - corners[0]], &opts()).unwrap();
    check(&m.mesh, &[cw]);
}

#[test]
fn mapped_hostile() {
    let l = rect(0.0, 0.0, 4.0, 2.0, 4, 2);
    assert!(mapped_quad_mesh(&l, [0, 4, 6, 11], &opts()).is_err()); // unequal opposite sides
    assert!(mapped_quad_mesh(&l, [0, 6, 4, 10], &opts()).is_err()); // out of order
    assert!(mapped_quad_mesh(&l, [0, 4, 6, 99], &opts()).is_err());
    assert!(mapped_quad_mesh(&l, [0, 0, 6, 10], &opts()).is_err());
    let mut nan = l.clone();
    nan[5] = [f64::NAN, 0.0];
    assert!(mapped_quad_mesh(&nan, [0, 4, 6, 10], &opts()).is_err());
    // A region whose TFI folds: the top side dips almost to the bottom between its corners.
    let fold = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0], [4.0, 4.0], [3.0, 4.0], [2.2, 0.2], [1.8, 0.2], [1.0, 4.0]];
    let r = mapped_quad_mesh(&fold, [0, 4, 5, 9], &opts());
    assert!(matches!(r, Err(ref e) if matches!(e.kind(), Some(crate::Kind::Failed | crate::Kind::Invalid))), "{r:?}");
}

#[test]
fn transfinite_in_space() {
    let bottom = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let top = [[0.0, 1.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]];
    let left = [[0.0, 0.0, 0.0], [0.0, 1.0, 1.0]];
    let right = [[2.0, 0.0, 0.0], [2.0, 1.0, 1.0]];
    let g = transfinite_grid(&bottom, &right, &top, &left).unwrap();
    assert_eq!(g.len(), 6);
    assert_eq!(g[4], [1.0, 1.0, 1.0]);
    let side = |a: [f64; 2], b: [f64; 2]| [a, [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0], b];
    let g =
        transfinite_grid(&side([0.0, 0.0], [2.0, 0.0]), &side([2.0, 0.0], [2.0, 2.0]), &side([0.0, 2.0], [2.0, 2.0]), &side([0.0, 0.0], [0.0, 2.0]))
            .unwrap();
    assert_eq!(g[4], [1.0, 1.0]);
    // Corners that do not match, a zero-length side, empty sides.
    assert!(transfinite_grid(&bottom, &right, &top, &[[0.0, 0.0, 0.1], [0.0, 1.0, 1.0]]).is_err());
    assert!(transfinite_grid(&[[0.0; 2], [0.0; 2]], &[[0.0; 2], [0.0, 1.0]], &[[0.0, 1.0], [0.0, 1.0]], &[[0.0; 2], [0.0, 1.0]]).is_err());
    assert!(transfinite_grid::<2>(&[], &[], &[], &[]).is_err());
}

#[test]
fn planar_quad_metrics() {
    use crate::quality::{Metric, planar_quad_quality, planar_quad_validity};
    let rect = [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]];
    let q = planar_quad_quality(&rect).unwrap();
    assert_eq!(q.validity(), Some(Validity::Valid));
    assert_eq!(q.size, 2.0);
    assert_eq!(q.get(Metric::AspectRatio), Some(2.0));
    assert_eq!(q.get(Metric::Skew), Some(0.0));
    assert_eq!(q.get(Metric::ScaledJacobian), Some(1.0));
    assert_eq!(q.get(Metric::JacobianRatio), Some(1.0));
    // A rhombus with 60° corners: its principal axes are its edge directions, so skew cos 60°
    // and aspect 1; the scaled Jacobian is sin 60°.
    let h = 3.0_f64.sqrt() / 2.0;
    let rh = [[0.0, 0.0], [1.0, 0.0], [1.5, h], [0.5, h]];
    let q = planar_quad_quality(&rh).unwrap();
    assert!((q.get(Metric::ScaledJacobian).unwrap() - h).abs() < 1e-12);
    assert!((q.get(Metric::Skew).unwrap() - 0.5).abs() < 1e-12);
    assert!((q.get(Metric::AspectRatio).unwrap() - 1.0).abs() < 1e-12);
    // Clockwise, bow tie, concave (dart), degenerate (a straight corner): invalid by the exact test.
    let cw = [rect[0], rect[3], rect[2], rect[1]];
    assert_eq!(planar_quad_validity(&cw).unwrap(), Validity::Invalid);
    assert_eq!(planar_quad_validity(&[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]]).unwrap(), Validity::Invalid);
    assert_eq!(planar_quad_validity(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.7, 0.3]]).unwrap(), Validity::Invalid);
    assert_eq!(planar_quad_validity(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0]]).unwrap(), Validity::Invalid);
    // A corner off by one ulp from straight is decided exactly.
    let x = 1.0 + f64::EPSILON;
    assert_eq!(planar_quad_validity(&[[0.0, 0.0], [1.0, 0.0], [2.0, x - 1.0], [0.0, 1.0]]).unwrap(), Validity::Valid);
    assert!(planar_quad_validity(&[[0.0, 0.0], [1.0, f64::NAN], [1.0, 1.0], [0.0, 1.0]]).is_err());
    let q = planar_quad_quality(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0]]).unwrap();
    assert_eq!(q.get(Metric::JacobianRatio), Some(0.0));
}

#[test]
#[ignore = "timing and quality report: cargo test --release -p fieldcraft-mesh quad_report -- --ignored --nocapture"]
fn quad_report() {
    let cases: Vec<(&str, Vec<Vec<[f64; 2]>>, f64)> = vec![
        ("rectangle", vec![rect(0.0, 0.0, 10.0, 5.0, 20, 10)], 0.5),
        ("plate with hole", vec![rect(0.0, 0.0, 20.0, 10.0, 40, 20), circle(10.0, 5.0, 2.5, 32)], 0.5),
        ("star", vec![subdivide(&star(10.0, 4.0), 4)], 1.0),
        ("large plate", vec![rect(0.0, 0.0, 200.0, 100.0, 1000, 500), circle(50.0, 50.0, 20.0, 400), circle(150.0, 50.0, 20.0, 400)], 0.2),
    ];
    for (name, loops, h) in cases {
        let t = std::time::Instant::now();
        let m = free_quad_mesh(&loops, &UniformSize(h), &opts()).unwrap();
        println!("{name}: {} quads in {:?}, {:?}", m.quads.len(), t.elapsed(), m.stats);
        let t = std::time::Instant::now();
        let _ = free_quad_mesh(&loops, &UniformSize(h), &QuadOptions { smoothing_sweeps: 0, ..opts() }).unwrap();
        println!("  without smoothing {:?}", t.elapsed());
    }
}

#[test]
fn unsmoothed_meshes_are_valid() {
    // Regression: recombination once accepted nearly straight corners that the split turned
    // degenerate after rounding the mid nodes (this plate had one before the quality floor).
    let loops = vec![rect(0.0, 0.0, 60.0, 30.0, 300, 150), circle(15.0, 15.0, 6.0, 240)];
    let m = free_quad_mesh(&loops, &UniformSize(0.2), &QuadOptions { smoothing_sweeps: 0, ..opts() }).unwrap();
    check(&m, &loops);
}

#[test]
fn free_triangle_mesh_of_a_plate_with_a_hole() {
    // Odd point counts are fine for triangles.
    let outer = rect(0.0, 0.0, 10.0, 6.0, 9, 5);
    let hole: Vec<[f64; 2]> = circle(4.0, 3.0, 1.5, 17).into_iter().rev().collect();
    let loops = vec![outer.clone(), hole.clone()];
    let m = free_tri_mesh(&loops, &UniformSize(0.7), &opts()).unwrap();
    assert_eq!(m.boundary_nodes, outer.len() + hole.len());
    let boundary: Vec<[f64; 2]> = outer.iter().chain(&hole).copied().collect();
    assert_eq!(&m.nodes[..m.boundary_nodes], boundary.as_slice(), "boundary nodes are the input points");
    let area: f64 = m
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| m.nodes[v as usize]);
            let s = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
            assert!(s > 0.0, "counterclockwise");
            s
        })
        .sum();
    let want = signed_area(&outer).abs() - signed_area(&hole).abs();
    assert!((area - want).abs() < 1e-9 * want, "{area} vs {want}");
    assert!(m.triangles.len() > 100, "refined to the size: {}", m.triangles.len());
    // Deterministic.
    assert_eq!(free_tri_mesh(&loops, &UniformSize(0.7), &opts()).unwrap(), m);
    // Hostile.
    assert!(free_tri_mesh(&[vec![[0.0, 0.0], [1.0, 0.0]]], &UniformSize(0.5), &opts()).is_err());
    assert!(free_tri_mesh(std::slice::from_ref(&outer), &UniformSize(0.001), &QuadOptions { max_elements: 100, ..opts() }).is_err());
    assert!(free_tri_mesh(&[vec![[0.0, 0.0], [f64::NAN, 0.0], [0.0, 1.0]]], &UniformSize(0.5), &opts()).is_err());
}
