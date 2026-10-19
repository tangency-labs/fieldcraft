use fieldcraft_model::ElementKind;

use super::*;
use crate::quality::{GateOptions, enclosed_volume, validity_gate_elements};

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(1);
        order::mix(self.0)
    }
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            v.swap(i, (self.next_u64() % (i as u64 + 1)) as usize);
        }
    }
}

fn random(seed: u64, n: usize) -> Vec<[f64; 3]> {
    let mut r = Rng(seed);
    (0..n).map(|_| [r.unit(), r.unit(), r.unit()]).collect()
}

fn grid(nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                out.push([x as f64, y as f64, z as f64]);
            }
        }
    }
    out
}

/// The finite tetrahedra as sorted coordinate quadruples: comparable across vertex numberings.
fn tet_set(tri: &Triangulation) -> Vec<[[u64; 3]; 4]> {
    let mut out: Vec<[[u64; 3]; 4]> = tri
        .finite_tets()
        .map(|(_, vs)| {
            let mut q = vs.map(|v| tri.points()[v as usize].coords().map(f64::to_bits));
            q.sort_unstable();
            q
        })
        .collect();
    out.sort_unstable();
    out
}

/// The validity gate on the finite tetrahedra, with the hull volume expected; returns the volume.
fn gate(tri: &Triangulation, expected: f64) -> f64 {
    let elements = tri.finite_tets().map(|(t, vs)| (t + 1, ElementKind::Tet4, vs.map(|v| tri.points()[v as usize].coords())));
    let report = validity_gate_elements(elements, &GateOptions { expected_volume: Some(expected), volume_tolerance: 1e-9 }).unwrap();
    assert_eq!(report.checked, tri.finite_tet_count());
    report.into_result().unwrap().volume
}

/// The volume enclosed by the hull faces, outward oriented, through the quality module.
fn hull_volume(tri: &Triangulation) -> f64 {
    let points: Vec<[f64; 3]> = tri.points().iter().map(Point3::coords).collect();
    let mut triangles = Vec::new();
    for t in 0..tri.slots() as u32 {
        let Some(vs) = tri.tet(t) else { continue };
        let Some(g) = vs.iter().position(|&v| v == GHOST) else { continue };
        // Moving the ghost's slot to the end takes 3 - g transpositions.
        let mut f: Vec<usize> = vs.iter().filter(|&&v| v != GHOST).map(|&v| v as usize).collect();
        if (3 - g) % 2 == 1 {
            f.swap(0, 1);
        }
        triangles.push([f[0], f[1], f[2]]);
    }
    enclosed_volume(&points, &triangles).unwrap()
}

/// Full check: invariants, brute-force empty spheres, the tet4 validity gate.
fn check(tri: &Triangulation, hull_volume: f64) {
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
    gate(tri, hull_volume);
}

/// Delaunay tetrahedra by definition: every positively oriented 4-subset whose perturbed
/// circumsphere excludes every other point.
fn brute_force(points: &[[f64; 3]]) -> Vec<[[u64; 3]; 4]> {
    let p: Vec<Point3> = points.iter().map(|&c| Point3::try_from(c).unwrap()).collect();
    let n = p.len();
    let mut out = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                for l in k + 1..n {
                    let (a, b, c, d) = match orient3d(&p[i], &p[j], &p[k], &p[l]) {
                        Sign::Zero => continue,
                        Sign::Positive => (&p[i], &p[j], &p[k], &p[l]),
                        Sign::Negative => (&p[j], &p[i], &p[k], &p[l]),
                    };
                    let empty = (0..n).filter(|m| ![i, j, k, l].contains(m)).all(|m| insphere_sos(a, b, c, d, &p[m]) == Sign::Negative);
                    if empty {
                        let mut q = [i, j, k, l].map(|v| points[v].map(f64::to_bits));
                        q.sort_unstable();
                        out.push(q);
                    }
                }
            }
        }
    }
    out.sort_unstable();
    out
}

#[test]
fn one_tetrahedron_and_its_ghosts() {
    let tri = Triangulation::new(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).unwrap();
    assert_eq!(tri.dimension(), 3);
    assert_eq!(tri.finite_tet_count(), 1);
    assert_eq!((0..tri.slots() as u32).filter(|&t| tri.is_ghost(t)).count(), 4);
    check(&tri, 1.0 / 6.0);
}

#[test]
fn random_points_are_delaunay_and_valid() {
    let points = random(1, 1500);
    let tri = Triangulation::new(&points).unwrap();
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
    // The tetrahedra fill the volume the hull faces (the ghosts' finite faces) enclose.
    let vol = hull_volume(&tri);
    assert!(vol > 0.9 && vol < 1.0, "{vol}");
    gate(&tri, vol);
}

#[test]
fn matches_the_definition_on_small_sets() {
    let sets: Vec<Vec<[f64; 3]>> = vec![
        random(7, 9),
        random(8, 14),
        grid(3, 3, 2),
        grid(3, 3, 3),
        // A cube's corners plus its centre: every face cospherical.
        {
            let mut g = grid(2, 2, 2);
            g.push([0.5, 0.5, 0.5]);
            g
        },
        // Coplanar on two planes.
        vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [2.0, 2.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 1.0], [2.0, 2.0, 1.0]],
    ];
    for points in sets {
        let tri = Triangulation::new(&points).unwrap();
        tri.verify().unwrap();
        assert_eq!(tet_set(&tri), brute_force(&points), "{points:?}");
    }
}

#[test]
fn the_result_does_not_depend_on_the_input_order() {
    for points in [random(2, 400), grid(5, 4, 6)] {
        let reference = tet_set(&Triangulation::new(&points).unwrap());
        let mut r = Rng(99);
        for _ in 0..3 {
            let mut shuffled = points.clone();
            r.shuffle(&mut shuffled);
            let tri = Triangulation::new(&shuffled).unwrap();
            tri.verify().unwrap();
            assert_eq!(tet_set(&tri), reference);
        }
    }
}

#[test]
fn deterministic_bit_for_bit() {
    let points = random(3, 3000);
    let a = Triangulation::new(&points).unwrap();
    let b = Triangulation::new(&points).unwrap();
    assert_eq!(a.tets, b.tets);
    assert_eq!(a.adj, b.adj);
}

#[test]
fn grids_are_cospherical_everywhere() {
    let tri = Triangulation::new(&grid(7, 6, 5)).unwrap();
    check(&tri, 6.0 * 5.0 * 4.0);
    assert_eq!(hull_volume(&tri), 120.0);
    // A rotated-in-place grid with non-integer spacing.
    let skew: Vec<[f64; 3]> = grid(6, 6, 6).into_iter().map(|[x, y, z]| [0.1 * x, 0.3 * y + 1e3, 0.7 * z - 5.0]).collect();
    let tri = Triangulation::new(&skew).unwrap();
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
}

#[test]
fn cospherical_and_cocircular_points() {
    // Every integer point on the sphere x² + y² + z² = 9 (all cospherical).
    let mut sphere = Vec::new();
    for x in -3i32..=3 {
        for y in -3i32..=3 {
            for z in -3i32..=3 {
                if x * x + y * y + z * z == 9 {
                    sphere.push([f64::from(x), f64::from(y), f64::from(z)]);
                }
            }
        }
    }
    assert_eq!(sphere.len(), 30);
    let tri = Triangulation::new(&sphere).unwrap();
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
    assert_eq!(tet_set(&tri), brute_force(&sphere));

    // Two regular 12-gons on parallel planes (a CAD cylinder's rings).
    let mut cylinder = Vec::new();
    for z in [0.0, 4.0] {
        for i in 0..12 {
            let a = f64::from(i) * std::f64::consts::TAU / 12.0;
            cylinder.push([3.0 * a.cos(), 3.0 * a.sin(), z]);
        }
    }
    let tri = Triangulation::new(&cylinder).unwrap();
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
    assert_eq!(tet_set(&tri), brute_force(&cylinder));
}

#[test]
fn lower_dimensional_inputs_have_no_tetrahedra() {
    assert_eq!(Triangulation::new(&[]).unwrap().dimension(), 0);
    let one = Triangulation::new(&[[1.0, 2.0, 3.0]]).unwrap();
    assert_eq!((one.dimension(), one.finite_tet_count()), (0, 0));
    let line: Vec<[f64; 3]> = (0..50).map(|i| [f64::from(i), 2.0 * f64::from(i), -f64::from(i)]).collect();
    let tri = Triangulation::new(&line).unwrap();
    assert_eq!((tri.dimension(), tri.finite_tet_count()), (1, 0));
    tri.verify().unwrap();

    let plane: Vec<[f64; 3]> = grid(10, 10, 1).into_iter().map(|[x, y, _]| [x, y, x + 2.0 * y]).collect();
    let mut tri = Triangulation::new(&plane).unwrap();
    assert_eq!((tri.dimension(), tri.finite_tet_count()), (2, 0));
    assert_eq!(tri.locate([0.0, 0.0, 0.0]).unwrap(), None);
    // Duplicates of the planar points are still found.
    assert_eq!(tri.insert([3.0, 4.0, 11.0]).unwrap(), Insertion::Duplicate(34));
    // One point off the plane makes it 3D: every tetrahedron has that apex.
    let apex = match tri.insert([4.5, 4.5, 40.0]).unwrap() {
        Insertion::Inserted(v) => v,
        other => panic!("{other:?}"),
    };
    assert_eq!(tri.dimension(), 3);
    assert!(tri.finite_tets().all(|(_, vs)| vs.contains(&apex)));
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
}

#[test]
fn duplicates_merge_into_the_lowest_index() {
    let base = grid(4, 4, 4);
    let mut points = base.clone();
    points.extend(base.iter().rev());
    points.push([-0.0, 0.0, -0.0]);
    let tri = Triangulation::new(&points).unwrap();
    tri.verify().unwrap();
    for i in 0..base.len() {
        assert_eq!(tri.duplicate_of(i as u32), None);
        assert_eq!(tri.duplicate_of((2 * base.len() - 1 - i) as u32), Some(i as u32));
    }
    assert_eq!(tri.duplicate_of(points.len() as u32 - 1), Some(0));
    assert_eq!(tet_set(&tri), tet_set(&Triangulation::new(&base).unwrap()));
}

#[test]
fn steiner_points_insert_into_a_built_triangulation() {
    let mut tri = Triangulation::new(&grid(4, 4, 4)).unwrap();
    let mut r = Rng(5);
    for _ in 0..200 {
        let p = [4.0 * r.unit() - 0.5, 4.0 * r.unit() - 0.5, 4.0 * r.unit() - 0.5];
        assert!(matches!(tri.insert(p).unwrap(), Insertion::Inserted(_)));
    }
    // On a vertex, on an edge, on a face, on the hull's plane outside it.
    assert_eq!(tri.insert([1.0, 2.0, 3.0]).unwrap(), Insertion::Duplicate(27));
    assert_eq!(tri.duplicate_of(tri.points().len() as u32 - 1), Some(27));
    for p in [[1.5, 2.0, 3.0], [1.5, 1.5, 0.0], [5.0, 5.0, 0.0], [-1.0, 1.0, 3.0]] {
        assert!(matches!(tri.insert(p).unwrap(), Insertion::Inserted(_)), "{p:?}");
    }
    tri.verify().unwrap();
    tri.verify_empty_spheres().unwrap();
}

#[test]
fn locate_finds_the_containing_tetrahedron() {
    let tri = Triangulation::new(&grid(3, 3, 3)).unwrap();
    match tri.locate([0.5, 0.5, 0.5]).unwrap() {
        Some(Location::Inside(t)) => assert!(!tri.is_ghost(t)),
        other => panic!("{other:?}"),
    }
    assert!(matches!(tri.locate([9.0, 0.5, 0.5]).unwrap(), Some(Location::Outside(t)) if tri.is_ghost(t)));
    assert_eq!(tri.locate([2.0, 2.0, 2.0]).unwrap(), Some(Location::Vertex(26)));
    assert!(tri.locate([f64::NAN, 0.0, 0.0]).is_err());
}

#[test]
fn hostile_coordinates_are_rejected() {
    for bad in [[f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0], [1e300, 0.0, 0.0], [1e-300, 0.0, 0.0]] {
        assert!(matches!(Triangulation::new(&[[0.0; 3], bad]), Err(ref e) if e.kind() == Some(crate::Kind::Invalid)));
    }
    let mut tri = Triangulation::new(&grid(2, 2, 2)).unwrap();
    assert!(tri.insert([0.0, f64::NEG_INFINITY, 0.0]).is_err());
    tri.verify().unwrap();
    // Extreme but admissible magnitudes still triangulate exactly.
    let big: Vec<[f64; 3]> = random(4, 200).into_iter().map(|[x, y, z]| [x * 1e40, y * 1e-30 + 1.0, z]).collect();
    Triangulation::new(&big).unwrap().verify().unwrap();
}

#[test]
fn incident_tets_and_neighbors_are_consistent() {
    let tri = Triangulation::new(&random(6, 300)).unwrap();
    for v in 0..300 {
        let t = tri.incident_tet(v).unwrap();
        assert!(tri.tet(t).unwrap().contains(&v));
    }
    for (t, _) in tri.finite_tets() {
        for n in tri.neighbors(t).unwrap() {
            assert!(tri.neighbors(n).unwrap().contains(&t));
        }
    }
}
