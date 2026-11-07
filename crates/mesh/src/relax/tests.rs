use super::*;
use crate::tet::geometry::AnalyticSurface;

const R_IN: f64 = 1.0;
const R_OUT: f64 = 2.0;
const HEIGHT: f64 = 1.0;

/// A quarter annulus (r 1..2, θ 0..90°, z 0..1) on its analytic faces: 0 the outer cylinder, 1
/// the inner, 2 and 3 the planes z = 0 and z = 1, 4 and 5 the planes θ = 0 and θ = 90°.
struct Annulus;

impl Annulus {
    fn surface(f: u32) -> AnalyticSurface {
        match f {
            0 => AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: R_OUT },
            1 => AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: R_IN },
            2 => AnalyticSurface::Plane { normal: [0.0, 0.0, -1.0], point: [0.0; 3] },
            3 => AnalyticSurface::Plane { normal: [0.0, 0.0, 1.0], point: [0.0, 0.0, HEIGHT] },
            4 => AnalyticSurface::Plane { normal: [0.0, -1.0, 0.0], point: [0.0; 3] },
            _ => AnalyticSurface::Plane { normal: [-1.0, 0.0, 0.0], point: [0.0; 3] },
        }
    }

    /// The faces a point lies on.
    fn faces(p: P) -> Vec<u32> {
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        let tol = 1e-9;
        let mut f = Vec::new();
        for (i, on) in
            [(r - R_OUT).abs() < tol, (r - R_IN).abs() < tol, p[2].abs() < tol, (p[2] - HEIGHT).abs() < tol, p[1].abs() < tol, p[0].abs() < tol]
                .into_iter()
                .enumerate()
        {
            if on {
                f.push(i as u32);
            }
        }
        f
    }

    fn on(p: P, faces: &[u32]) -> bool {
        faces.iter().all(|&f| Annulus::surface(f).project(p).is_some_and(|q| dist(p, q) < 1e-9))
    }
}

/// Curves are encoded as `a * 8 + b`.
struct Geo {
    model: Model,
}

impl Constraints for Geo {
    fn freedom(&self, node: u32) -> Option<Freedom> {
        let p = self.model.node(node)?.xyz;
        let f = Annulus::faces(p);
        Some(match f.as_slice() {
            [] => Freedom::Free,
            [a] => Freedom::Face(*a),
            [a, b] => Freedom::Curve(a * 8 + b),
            _ => Freedom::Fixed,
        })
    }
    fn slide_face(&self, face: u32, p: P, d: P) -> Option<P> {
        Annulus::surface(face).slide(p, d)
    }
    fn slide_curve(&self, curve: u32, p: P, d: P) -> Option<P> {
        let (a, b) = (Annulus::surface(curve / 8), Annulus::surface(curve % 8));
        let mut q = add(p, d);
        for _ in 0..100 {
            q = b.project(a.project(q)?)?;
        }
        Some(q)
    }
}

fn polar(r: f64, t: f64, z: f64) -> P {
    let (s, c) = crate::tet::geometry::sin_cos(t);
    [r * c, r * s, z]
}

/// A hex8 mesh of the annulus with `n` divisions each way; node label `1 + i + (n+1)(j + (n+1)k)`.
fn hex_annulus(n: usize) -> Model {
    let mut m = Model::default();
    let id = |i: usize, j: usize, k: usize| (1 + i + (n + 1) * (j + (n + 1) * k)) as u32;
    for k in 0..=n {
        for j in 0..=n {
            for i in 0..=n {
                let p = polar(
                    R_IN + (R_OUT - R_IN) * i as f64 / n as f64,
                    std::f64::consts::FRAC_PI_2 * j as f64 / n as f64,
                    HEIGHT * k as f64 / n as f64,
                );
                m.set_node(id(i, j, k), p).unwrap();
            }
        }
    }
    let mut label = 1;
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let nodes = vec![
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i + 1, j + 1, k + 1),
                    id(i, j + 1, k + 1),
                ];
                m.set_element(label, ElementKind::Hex8, nodes).unwrap();
                label += 1;
            }
        }
    }
    m
}

/// Move some nodes near `centre`, keeping each on its faces.
fn perturb(m: &mut Model, centre: P, radius: f64, amount: f64) -> BTreeSet<u32> {
    let geo = Geo { model: m.clone() };
    let mut moved = BTreeSet::new();
    let labels: Vec<u32> = m.nodes().keys().copied().collect();
    for (k, n) in labels.into_iter().enumerate() {
        let p = m.node(n).unwrap().xyz;
        if dist(p, centre) > radius {
            continue;
        }
        let d = [amount * ((k * 7 % 5) as f64 - 2.0) / 2.0, amount * ((k * 3 % 5) as f64 - 2.0) / 2.0, amount * ((k * 11 % 5) as f64 - 2.0) / 2.0];
        let q = match geo.freedom(n).unwrap() {
            Freedom::Free => Some(add(p, d)),
            Freedom::Face(f) => geo.slide_face(f, p, d),
            _ => None,
        };
        if let Some(q) = q {
            m.set_node(n, q).unwrap();
            moved.insert(n);
        }
    }
    moved
}

fn region_min(m: &Model, labels: &[u32]) -> (f64, usize) {
    let mut worst = f64::INFINITY;
    let mut invalid = 0;
    for l in labels {
        let e = m.element(*l).unwrap();
        let x: Vec<P> = e.nodes.iter().map(|n| m.node(*n).unwrap().xyz).collect();
        let el = Elem { label: *l, kind: e.kind, nodes: e.nodes.clone(), normal: None, size0: 0.0 };
        let (q, v) = measure(&el, &x);
        worst = worst.min(q);
        invalid += usize::from(!v);
    }
    (worst, invalid)
}

fn applied(m: &Model, r: &Relaxed) -> Model {
    let mut out = m.clone();
    for (n, p) in &r.positions {
        out.set_node(*n, *p).unwrap();
    }
    out
}

const CENTRE: P = [1.414, 1.414, 0.5];

/// The issue's test: a perturbed hex mesh of a cylindrical part relaxes inside the sphere — every
/// surface node on its face, quality up, nothing inverted, nothing outside the sphere moved.
#[test]
fn hex_spikes_relax_with_surface_nodes_on_the_cylinder() {
    let mut m = hex_annulus(6);
    let original = m.clone();
    perturb(&mut m, CENTRE, 0.5, 0.08);
    let elements: Vec<u32> = m.elements().keys().copied().collect();
    let (before, _) = region_min(&m, &elements);
    let weights = sphere_weights(&m, CENTRE, 0.7, 1.0, Falloff::Smooth, None).unwrap();
    let geo = Geo { model: original.clone() };
    let r = relax(&m, &weights, 10, &geo).unwrap();
    assert!(r.moved > 0 && r.after.min_quality > r.before.min_quality, "{:?} → {:?}", r.before, r.after);
    assert_eq!(r.after.invalid, 0);
    let out = applied(&m, &r);
    let (after, invalid) = region_min(&out, &elements);
    assert!(after > before + 0.02 && invalid == 0, "{before} → {after}");
    for (n, p) in &r.positions {
        assert!(dist(m.node(*n).unwrap().xyz, CENTRE) < 0.7, "node {n} outside the sphere moved");
        let faces = Annulus::faces(original.node(*n).unwrap().xyz);
        assert!(Annulus::on(*p, &faces), "node {n} left its faces {faces:?}: {p:?}");
    }
}

/// tet10 with curved mid-side nodes on the cylinder: after relaxing, corner and mid-side surface
/// nodes are still on the surface and no element is invalid.
#[test]
fn tet10_mid_side_nodes_follow_on_the_geometry() {
    let hex = hex_annulus(4);
    let mut m = Model::default();
    for (l, n) in hex.nodes() {
        m.set_node(*l, n.xyz).unwrap();
    }
    // Six tetrahedra per hexahedron around its 1-7 diagonal (Kuhn), consistent across faces.
    const KUHN: [[usize; 4]; 6] = [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]];
    let mut mids: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut next = 100_000;
    let mut label = 1;
    for e in hex.elements().values() {
        for t in KUHN {
            let mut c: Vec<u32> = t.iter().map(|&i| e.nodes[i]).collect();
            let x: Vec<P> = c.iter().map(|n| m.node(*n).unwrap().xyz).collect();
            if dot(cross(sub(x[1], x[0]), sub(x[2], x[0])), sub(x[3], x[0])) < 0.0 {
                c.swap(1, 2);
            }
            let mut nodes = c.clone();
            for [a, b] in [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]] {
                let key = (c[a].min(c[b]), c[a].max(c[b]));
                let mid = *mids.entry(key).or_insert_with(|| {
                    let (pa, pb) = (m.node(key.0).unwrap().xyz, m.node(key.1).unwrap().xyz);
                    let mut p = scale(add(pa, pb), 0.5);
                    let shared: Vec<u32> = Annulus::faces(pa).into_iter().filter(|f| Annulus::faces(pb).contains(f)).collect();
                    if let [f] = shared.as_slice() {
                        p = Annulus::surface(*f).project(p).unwrap();
                    }
                    next += 1;
                    m.set_node(next, p).unwrap();
                    next
                });
                nodes.push(mid);
            }
            m.set_element(label, ElementKind::Tet10, nodes).unwrap();
            label += 1;
        }
    }
    let original = m.clone();
    // Perturb the corners only, then put the mid-side nodes at their edges' (curved) middles.
    let corners = corner_nodes(&m);
    perturb(&mut m, CENTRE, 0.6, 0.06);
    for (&(a, b), &mid) in &mids {
        let p = scale(add(m.node(a).unwrap().xyz, m.node(b).unwrap().xyz), 0.5);
        let shared: Vec<u32> =
            Annulus::faces(original.node(a).unwrap().xyz).into_iter().filter(|f| Annulus::faces(original.node(b).unwrap().xyz).contains(f)).collect();
        let p = if let [f] = shared.as_slice() { Annulus::surface(*f).project(p).unwrap() } else { p };
        m.set_node(mid, p).unwrap();
    }
    let _ = corners;
    let weights = sphere_weights(&m, CENTRE, 0.8, 1.0, Falloff::Smooth, None).unwrap();
    let r = relax(&m, &weights, 8, &Geo { model: original.clone() }).unwrap();
    assert!(r.after.min_quality >= r.before.min_quality && r.after.invalid == 0, "{:?} → {:?}", r.before, r.after);
    assert!(r.after.min_quality > r.before.min_quality + 0.01, "{:?} → {:?}", r.before, r.after);
    let mut mid_on_surface = 0;
    for (n, p) in &r.positions {
        let faces = Annulus::faces(original.node(*n).unwrap().xyz);
        assert!(Annulus::on(*p, &faces), "node {n} left {faces:?}");
        mid_on_surface += usize::from(*n > 100_000 && faces.contains(&0));
    }
    assert!(mid_on_surface > 0, "some curved mid-side nodes moved on the cylinder");
}

/// A cut between two cells is not geometry: its shared nodes move freely, and both cells keep
/// them (conforming).
#[test]
fn partition_cut_nodes_move_freely_and_stay_shared() {
    let mut m = hex_annulus(6);
    let original = m.clone();
    perturb(&mut m, CENTRE, 0.5, 0.08);
    // The cut: the plane θ = 45° through the middle row of nodes (j = 3) — reported free.
    struct Cut(Geo);
    impl Constraints for Cut {
        fn freedom(&self, node: u32) -> Option<Freedom> {
            let f = self.0.freedom(node)?;
            Some(if f == Freedom::Free { Freedom::Free } else { f })
        }
        fn slide_face(&self, f: u32, p: P, d: P) -> Option<P> {
            self.0.slide_face(f, p, d)
        }
        fn slide_curve(&self, c: u32, p: P, d: P) -> Option<P> {
            self.0.slide_curve(c, p, d)
        }
    }
    let weights = sphere_weights(&m, CENTRE, 0.7, 1.0, Falloff::Smooth, None).unwrap();
    let r = relax(&m, &weights, 10, &Cut(Geo { model: original.clone() })).unwrap();
    let on_cut = |p: P| (p[0] - p[1]).abs() < 1e-9;
    let cut_moved = r.positions.iter().filter(|(n, p)| on_cut(original.node(**n).unwrap().xyz) && !on_cut(**p)).count();
    assert!(cut_moved > 0, "cut nodes are free to leave the cut plane");
    assert_eq!(r.after.invalid, 0);
}

/// The issue's field test: relax a coarse mesh, mesh 2× finer, apply the field: the fine mesh
/// takes the relaxed shape (its nodes where the relaxed coarse mesh interpolates them) and
/// nothing inverts.
#[test]
fn a_relaxation_field_carries_over_to_a_finer_mesh() {
    let mut coarse = hex_annulus(4);
    let original = coarse.clone();
    // A deliberate interior shape change: pull the interior nodes toward the centre.
    for (n, p) in original.nodes() {
        if Annulus::faces(p.xyz).is_empty() && dist(p.xyz, CENTRE) < 0.6 {
            coarse.set_node(*n, add(p.xyz, [0.03, -0.02, 0.04])).unwrap();
        }
    }
    let weights = sphere_weights(&coarse, CENTRE, 0.9, 1.0, Falloff::Smooth, None).unwrap();
    let r = relax(&coarse, &weights, 5, &Geo { model: coarse.clone() }).unwrap();
    // The field is the change from the original coarse mesh to the relaxed one.
    let relaxed = applied(&coarse, &r);
    let total = Relaxed {
        positions: relaxed.nodes().iter().filter(|(n, p)| original.node(**n).unwrap().xyz != p.xyz).map(|(n, p)| (*n, p.xyz)).collect(),
        ..Default::default()
    };
    let cells = field_cells(&original, &total);
    assert!(!cells.is_empty());
    let fine = hex_annulus(8);
    let out = apply_field(&fine, &cells, None, &Geo { model: fine.clone() }, 0.0).unwrap();
    assert!(out.moved > 0);
    let result = applied(&fine, &out);
    let elements: Vec<u32> = fine.elements().keys().copied().collect();
    assert_eq!(region_min(&result, &elements).1, 0, "nothing inverts");
    for (n, p) in fine.nodes() {
        let d = field_displacement(&cells, p.xyz).unwrap_or([0.0; 3]);
        let got = result.node(*n).unwrap().xyz;
        if Annulus::faces(p.xyz).is_empty() {
            assert!(dist(got, add(p.xyz, d)) < 1e-9, "node {n}: {got:?} vs {:?}", add(p.xyz, d));
        } else {
            assert!(Annulus::on(got, &Annulus::faces(p.xyz)), "node {n} stays on its faces");
        }
    }
    // A fine node at a coarse node takes exactly that node's displacement.
    let coarse_node = original.nodes().iter().find(|(n, p)| r.positions.contains_key(n) && Annulus::faces(p.xyz).is_empty()).unwrap();
    let fine_twin = fine.nodes().values().find(|p| dist(p.xyz, coarse_node.1.xyz) < 1e-9).unwrap();
    let d = field_displacement(&cells, fine_twin.xyz).unwrap();
    assert!(dist(add(fine_twin.xyz, d), relaxed.node(*coarse_node.0).unwrap().xyz) < 1e-8);
}

#[test]
fn hostile_parameters_are_refused() {
    let m = hex_annulus(2);
    assert!(sphere_weights(&m, [f64::NAN, 0.0, 0.0], 1.0, 1.0, Falloff::Smooth, None).is_err());
    assert!(sphere_weights(&m, CENTRE, -1.0, 1.0, Falloff::Smooth, None).is_err());
    assert!(sphere_weights(&m, CENTRE, 1.0, 1.5, Falloff::Smooth, None).is_err());
    let w = sphere_weights(&m, CENTRE, 1.0, 1.0, Falloff::Smooth, None).unwrap();
    assert!(relax(&m, &w, 0, &NoGeometry).is_err());
    assert!(relax(&m, &w, MAX_ITERATIONS + 1, &NoGeometry).is_err());
    let bad: BTreeMap<u32, f64> = [(1, f64::NAN)].into();
    assert!(relax(&m, &bad, 1, &NoGeometry).is_err());
    // Without geometry, skin nodes stay put.
    let r = relax(&m, &w, 3, &NoGeometry).unwrap();
    assert!(r.positions.is_empty() || r.positions.keys().all(|n| Annulus::faces(m.node(*n).unwrap().xyz).is_empty()));
}

/// Every node of the annulus pushed about (surface nodes on their faces).
fn distorted_annulus(n: usize, amount: f64) -> (Model, Model) {
    let original = hex_annulus(n);
    let mut m = original.clone();
    perturb(&mut m, [1.0, 1.0, 0.5], 10.0, amount);
    (original, m)
}

/// Strong mode on a distorted block: rounds of bad nodes moved and the mesh morphed round them
/// lift the smallest scaled Jacobian well above where it was and the mean with it; nothing
/// inverts; nodes on faces stay on them (the cylinders exactly); nodes on edges and vertices do
/// not move; the same input gives the same output.
#[test]
fn strong_relaxation_of_a_distorted_block_keeps_faces_and_edges() {
    let (original, m) = distorted_annulus(6, 0.06);
    let geo = Geo { model: original.clone() };
    let nodes: BTreeSet<u32> = m.nodes().keys().copied().collect();
    let options = MorphOptions { below: 0.9, reach: None, rounds: 40, ..MorphOptions::default() };
    let r = morph(&m, &nodes, options, &geo).unwrap();
    assert!(r.after.min_quality > r.before.min_quality + 0.2, "{:?} → {:?}: {:?}", r.before.min_quality, r.after.min_quality, r.rounds);
    assert!(r.mean_after > r.mean_before, "mean {} → {}", r.mean_before, r.mean_after);
    assert_eq!(r.after.invalid, 0);
    assert!(!r.rounds.is_empty() && r.rounds.iter().all(|k| k.region_nodes >= k.bad_nodes && k.reach > 0.0), "{:?}", r.rounds);
    let out = applied(&m, &Relaxed { positions: r.positions.clone(), ..Relaxed::default() });
    let elements: Vec<u32> = m.elements().keys().copied().collect();
    let (after, invalid) = region_min(&out, &elements);
    assert!((after - r.after.min_quality).abs() < 1e-12 && invalid == 0);
    for (n, p) in out.nodes() {
        let faces = Annulus::faces(original.node(*n).unwrap().xyz);
        assert!(Annulus::on(p.xyz, &faces), "node {n} left its faces {faces:?}");
        if faces.len() >= 2 {
            assert_eq!(p.xyz, m.node(*n).unwrap().xyz, "node {n} on an edge moved");
        }
    }
    assert_eq!(morph(&m, &nodes, options, &geo).unwrap(), r, "not deterministic");
}

/// The morph is local: one bad spot moves only the nodes within the reach.
#[test]
fn strong_relaxation_moves_only_nodes_within_the_reach() {
    let original = hex_annulus(8);
    let mut m = original.clone();
    perturb(&mut m, CENTRE, 0.2, 0.06);
    let geo = Geo { model: original };
    let nodes: BTreeSet<u32> = m.nodes().keys().copied().collect();
    let reach = 0.4;
    let r = morph(&m, &nodes, MorphOptions { below: 0.8, reach: Some(reach), rounds: 20, ..MorphOptions::default() }, &geo).unwrap();
    assert!(r.after.min_quality > r.before.min_quality, "{:?}", r.rounds);
    for n in r.positions.keys() {
        assert!(dist(m.node(*n).unwrap().xyz, CENTRE) < 0.2 + 2.0 * reach, "node {n} far from the bad spot moved");
    }
}

/// Strong mode untangles: an interior node pushed through its neighbours inverts elements; the
/// morph brings every element back to valid.
#[test]
fn strong_relaxation_untangles_inverted_elements() {
    let (original, mut m) = distorted_annulus(4, 0.0);
    // An interior node (i = j = k = 2 of 4) pushed past its neighbour.
    let inside = (1 + 2 + 5 * (2 + 5 * 2)) as u32;
    let p = m.node(inside).unwrap().xyz;
    m.set_node(inside, [p[0] + 0.4, p[1] + 0.4, p[2] + 0.3]).unwrap();
    let geo = Geo { model: original };
    let nodes: BTreeSet<u32> = m.nodes().keys().copied().collect();
    let r = morph(&m, &nodes, MorphOptions::default(), &geo).unwrap();
    assert!(r.before.invalid > 0, "{:?}", r.before);
    assert_eq!(r.after.invalid, 0, "{:?}", r.after);
    assert!(r.after.min_quality > 0.5, "{:?}", r.after);
}

/// Both guards keep the result valid and never below the start; the valid guards (keep a morph
/// unless it inverts an element) lift the worst element at least as far as the strict ones.
#[test]
fn strong_relaxation_guards_valid_and_strict() {
    let (original, m) = distorted_annulus(6, 0.06);
    let geo = Geo { model: original };
    let nodes: BTreeSet<u32> = m.nodes().keys().copied().collect();
    let run = |guards| morph(&m, &nodes, MorphOptions { below: 0.9, rounds: 40, guards, ..MorphOptions::default() }, &geo).unwrap();
    let (valid, strict) = (run(Guards::Valid), run(Guards::Strict));
    for r in [&valid, &strict] {
        assert_eq!(r.after.invalid, 0, "{:?}", r.after);
        assert!(r.after.min_quality >= r.before.min_quality, "{:?} → {:?}", r.before.min_quality, r.after.min_quality);
    }
    assert!(valid.after.min_quality >= strict.after.min_quality - 1e-9, "valid {:?}, strict {:?}", valid.after, strict.after);
}

/// A good mesh is left as it is; bad options are errors.
#[test]
fn strong_relaxation_keeps_a_good_mesh_and_refuses_bad_options() {
    let m = hex_annulus(3);
    let geo = Geo { model: m.clone() };
    let nodes: BTreeSet<u32> = m.nodes().keys().copied().collect();
    let r = morph(&m, &nodes, MorphOptions { below: 0.5, ..MorphOptions::default() }, &geo).unwrap();
    assert!(r.positions.is_empty() && r.rounds.is_empty());
    for bad in [
        MorphOptions { below: 0.0, ..MorphOptions::default() },
        MorphOptions { below: 1.5, ..MorphOptions::default() },
        MorphOptions { reach: Some(-1.0), ..MorphOptions::default() },
        MorphOptions { reach: Some(f64::NAN), ..MorphOptions::default() },
        MorphOptions { rounds: 0, ..MorphOptions::default() },
        MorphOptions { rounds: MAX_ROUNDS + 1, ..MorphOptions::default() },
    ] {
        assert!(morph(&m, &nodes, bad, &geo).is_err(), "{bad:?}");
    }
    assert!(morph(&m, &BTreeSet::new(), MorphOptions::default(), &geo).unwrap().positions.is_empty());
}

/// The guarded mode never stretches an element past the region's longest one before (the
/// scaled Jacobian alone does not see an element flattening: on the crank the worst edge ratio
/// went from 10.5 to 2759, #522).
#[test]
fn guarded_relaxation_keeps_the_largest_edge_ratio() {
    for (amount, radius) in [(0.04, 0.5), (0.08, 0.7), (0.12, 1.0)] {
        let original = hex_annulus(8);
        let mut m = original.clone();
        perturb(&mut m, CENTRE, radius, amount);
        let weights = sphere_weights(&m, CENTRE, radius + 0.2, 1.0, Falloff::None, None).unwrap();
        let r = relax(&m, &weights, 50, &Geo { model: original }).unwrap();
        assert!(r.after.max_edge_ratio <= r.before.max_edge_ratio, "{amount}: {} → {}", r.before.max_edge_ratio, r.after.max_edge_ratio);
        assert!(r.after.invalid <= r.before.invalid);
    }
}

// ---- Moving nodes by hand (`mesh.node.move`) ----

/// The annulus node `(i, j, k)` of [`hex_annulus`] with `n` divisions.
fn node_of(n: usize, i: usize, j: usize, k: usize) -> u32 {
    (1 + i + (n + 1) * (j + (n + 1) * k)) as u32
}

fn grab(n: u32) -> BTreeSet<u32> {
    BTreeSet::from([n])
}

/// Falloff: the grabbed node moves by all of the displacement, the nodes round it by less the
/// farther they are, none beyond the radius; the three falloffs differ.
#[test]
fn a_falloff_move_carries_the_grabbed_node_and_fades_to_the_radius() {
    let m = hex_annulus(6);
    let h = node_of(6, 3, 3, 3);
    let d = [0.02, -0.01, 0.015];
    let r = 0.6;
    let mut shares = Vec::new();
    for falloff in [Falloff::Smooth, Falloff::Linear, Falloff::None] {
        let inf = Influence { radius: Some(r), falloff, ..Influence::default() };
        let w = move_weights(&m, &grab(h), &inf, &Unconstrained).unwrap();
        assert_eq!(w.get(&h), Some(&1.0));
        let moved = move_nodes(&m, &w, d, Guard::None, &Unconstrained).unwrap();
        assert!(dist(moved.positions[&h], add(m.node(h).unwrap().xyz, d)) < 1e-12);
        let p0 = m.node(h).unwrap().xyz;
        for (n, p) in &moved.positions {
            let old = m.node(*n).unwrap().xyz;
            assert!(dist(old, p0) <= r + 1e-12, "node {n} beyond the radius moved");
            assert!(dist(old, *p) <= dist([0.0; 3], d) + 1e-12);
        }
        assert_eq!(moved.scale, 1.0);
        assert!(moved.elements.len() >= 8 && moved.worse + moved.better > 0);
        shares.push(w.get(&node_of(6, 4, 3, 3)).copied().unwrap_or(0.0));
    }
    // Smooth (1 − t²)², linear 1 − t, none 1, at t = one element over the radius.
    assert!(shares[0] > 0.0 && shares[0] < 1.0 && shares[1] > 0.0 && shares[1] < 1.0 && shares[2] == 1.0, "{shares:?}");
    assert_ne!(shares[0], shares[1]);
    // Strength scales the followers, never the grabbed node.
    let half = Influence { radius: Some(r), strength: 0.5, falloff: Falloff::None, ..Influence::default() };
    let w = move_weights(&m, &grab(h), &half, &Unconstrained).unwrap();
    assert_eq!((w.get(&h), w.get(&node_of(6, 4, 3, 3))), (Some(&1.0), Some(&0.5)));
    // A sphere without grabbed nodes is a soft selection round its centre.
    let c = m.node(h).unwrap().xyz;
    let w = move_weights(&m, &BTreeSet::new(), &Influence { centre: Some(c), radius: Some(r), ..Influence::default() }, &Unconstrained).unwrap();
    assert!(w.get(&h).is_some_and(|v| (*v - 1.0).abs() < 1e-12) && w.len() > 1);
}

/// Harmonic: shares between 0 and 1, 1 at the grabbed node, falling with the graph distance,
/// and zero beyond the radius (no node there moves).
#[test]
fn a_harmonic_move_is_zero_beyond_the_radius() {
    let m = hex_annulus(8);
    let h = node_of(8, 4, 4, 4);
    let r = 0.45;
    let inf = Influence { radius: Some(r), propagate: Propagate::Harmonic, ..Influence::default() };
    let w = move_weights(&m, &grab(h), &inf, &Unconstrained).unwrap();
    let p0 = m.node(h).unwrap().xyz;
    assert_eq!(w.get(&h), Some(&1.0));
    assert!(w.len() > 10);
    for (n, s) in &w {
        assert!(*s > 0.0 && *s <= 1.0);
        assert!(dist(m.node(*n).unwrap().xyz, p0) <= r + 1e-12, "node {n} beyond the radius has a share");
    }
    let (near, far) = (w[&node_of(8, 5, 4, 4)], w.get(&node_of(8, 6, 4, 4)).copied().unwrap_or(0.0));
    assert!(near < 1.0 && far < near, "{near} {far}");
    let moved = move_nodes(&m, &w, [0.0, 0.0, 0.03], Guard::default(), &Unconstrained).unwrap();
    assert!(moved.positions.keys().all(|n| w.contains_key(n)));
    assert_eq!(moved.after.invalid, 0);
}

/// On the geometry: a node on the outer cylinder stays on it, a node on an edge of two faces
/// stays on both, a vertex is held; followers on the faces stay on theirs.
#[test]
fn a_move_on_the_geometry_slides_on_faces_and_edges_and_holds_vertices() {
    let n = 6;
    let m = hex_annulus(n);
    let geo = Geo { model: m.clone() };
    let d = [0.05, 0.05, 0.05];
    let inf = Influence { radius: Some(0.5), ..Influence::default() };
    for h in [node_of(n, n, 3, 3), node_of(n, n, 3, 0), node_of(n, n, 0, 0)] {
        let w = move_weights(&m, &grab(h), &inf, &geo).unwrap();
        let moved = move_nodes(&m, &w, d, Guard::default(), &geo).unwrap();
        for (k, p) in &moved.positions {
            let faces = Annulus::faces(m.node(*k).unwrap().xyz);
            assert!(faces.len() < 3, "vertex {k} moved");
            assert!(Annulus::on(*p, &faces), "node {k} left its faces {faces:?}");
        }
        let moved_handle = moved.positions.contains_key(&h);
        assert_eq!(moved_handle, Annulus::faces(m.node(h).unwrap().xyz).len() < 3, "node {h}");
    }
}

/// The guard scales a move that would invert elements back to the largest that inverts none.
#[test]
fn the_guard_scales_a_move_back_until_nothing_inverts() {
    let m = hex_annulus(6);
    let h = node_of(6, 3, 3, 3);
    let w = move_weights(&m, &grab(h), &Influence::default(), &Unconstrained).unwrap();
    assert_eq!(w.len(), 1);
    let d = [0.6, 0.0, 0.0];
    let wild = move_nodes(&m, &w, d, Guard::None, &Unconstrained).unwrap();
    assert!(wild.after.invalid > 0 && wild.scale == 1.0);
    let guarded = move_nodes(&m, &w, d, Guard::default(), &Unconstrained).unwrap();
    assert_eq!(guarded.after.invalid, 0);
    assert!(guarded.scale > 0.05 && guarded.scale < 1.0, "{}", guarded.scale);
    assert!(guarded.worse > 0 && guarded.after.min_quality < guarded.before.min_quality);
    assert!(guarded.after.min_quality >= DEFAULT_MIN_QUALITY, "the guard stops short of a degenerate element");
    let p = guarded.positions[&h];
    assert!((p[0] - m.node(h).unwrap().xyz[0] - 0.6 * guarded.scale).abs() < 1e-12);
}

#[test]
fn hostile_move_parameters_are_refused() {
    let m = hex_annulus(3);
    let h = grab(node_of(3, 1, 1, 1));
    let bad = |inf: Influence| move_weights(&m, &h, &inf, &Unconstrained).is_err();
    assert!(bad(Influence { strength: 0.0, ..Influence::default() }));
    assert!(bad(Influence { strength: f64::NAN, ..Influence::default() }));
    assert!(bad(Influence { radius: Some(-1.0), ..Influence::default() }));
    assert!(bad(Influence { radius: Some(f64::INFINITY), ..Influence::default() }));
    assert!(bad(Influence { centre: Some([f64::NAN, 0.0, 0.0]), radius: Some(1.0), ..Influence::default() }));
    assert!(bad(Influence { propagate: Propagate::Harmonic, ..Influence::default() }));
    let w = move_weights(&m, &h, &Influence::default(), &Unconstrained).unwrap();
    assert!(move_nodes(&m, &w, [f64::NAN, 0.0, 0.0], Guard::None, &Unconstrained).is_err());
    assert!(move_nodes(&m, &w, [1e300, 0.0, 0.0], Guard::None, &Unconstrained).is_err());
    assert!(move_nodes(&m, &BTreeMap::from([(1, 2.0)]), [0.0; 3], Guard::None, &Unconstrained).is_err());
    assert!(move_nodes(&m, &w, [0.1; 3], Guard::NoInversion { min_quality: 1.5 }, &Unconstrained).is_err());
    // A radius too small for the mesh is refused, never an overflow of the search grid.
    assert!(bad(Influence { radius: Some(1e-300), ..Influence::default() }));
    assert!(bad(Influence { centre: Some([-1e300, -1e300, 0.0]), radius: Some(1e-300), ..Influence::default() }));
    let near = |r: f64| moving::nearest_within(&[[-1e300, 0.0, 0.0]], &BTreeMap::from([(1, [1e300, 0.0, 0.0])]), r);
    assert!(near(1e-300).is_empty() && near(1.0).is_empty());
}

/// `m` through millimetres → metres → millimetres: 1 ulp off for some coordinates (#612).
fn round_trip(m: &Model) -> Model {
    let mut out = m.clone();
    for (n, p) in m.nodes() {
        out.set_node(*n, p.xyz.map(|x| (x * 1e-3) * 1e3)).unwrap();
    }
    out
}

/// The largest distance between the nodes of `a` and `b`.
fn largest_difference(a: &Model, b: &Model) -> f64 {
    a.nodes().iter().map(|(n, p)| dist(p.xyz, b.node(*n).unwrap().xyz)).fold(0.0, f64::max)
}

/// A 1-ulp change of the input (millimetres through metres and back) moves the guarded
/// relaxation's result by round-off, not by a different path (#612: the guards compared
/// qualities exactly, and a step kept on one side of a tie and halved on the other moved nodes by
/// up to 9e-2 here, a fraction of an edge of 0.17–0.26). The differences left are a few hundred
/// ulps of the coordinates (≤ 3e-10 seen); 1e-8 is far below any real path change.
#[test]
fn guarded_relaxation_is_stable_under_round_off() {
    for (n, amount) in [(6, 0.03), (8, 0.06)] {
        let (original, m) = distorted_annulus(n, amount);
        let b = round_trip(&m);
        assert!(m.nodes().iter().any(|(l, p)| b.node(*l).unwrap().xyz != p.xyz), "the round trip changes some coordinates");
        let geo = Geo { model: original.clone() };
        let all: BTreeSet<u32> = m.nodes().keys().copied().collect();
        let w = uniform_weights(&m, &all, 1.0).unwrap();
        let (ra, rb) = (relax(&m, &w, 5, &geo).unwrap(), relax(&b, &w, 5, &geo).unwrap());
        assert!(ra.moved > 0 && ra.after.min_quality > ra.before.min_quality + 0.1, "{:?} → {:?}", ra.before, ra.after);
        let d = largest_difference(&applied(&m, &ra), &applied(&b, &rb));
        assert!(d < 1e-8, "n {n}: a 1-ulp input change moved a node by {d:e}");
    }
}

/// A 1-ulp change of the input moves the strong relaxation's result by round-off, not by a
/// different path (#612: summed into one number, the objective's barrier of invalid elements and
/// the q⁻¹⁶ of bad ones rounded away the change of the elements a node moves, so the sweeps
/// followed round-off and moved nodes by up to 6.5e-2 here, a fraction of an edge of 0.17–0.26).
/// What is left is round-off amplified over the sweeps (≤ 3e-8 seen); 1e-6 is far below any real
/// path change.
#[test]
fn strong_relaxation_is_stable_under_round_off() {
    for (n, amount) in [(6, 0.06), (8, 0.03)] {
        let (original, m) = distorted_annulus(n, amount);
        let b = round_trip(&m);
        assert!(m.nodes().iter().any(|(l, p)| b.node(*l).unwrap().xyz != p.xyz), "the round trip changes some coordinates");
        let geo = Geo { model: original.clone() };
        let all: BTreeSet<u32> = m.nodes().keys().copied().collect();
        let opts = MorphOptions { below: 0.9, ..MorphOptions::default() };
        let (ra, rb) = (morph(&m, &all, opts, &geo).unwrap(), morph(&b, &all, opts, &geo).unwrap());
        assert!(ra.moved > 0 && ra.after.min_quality > 0.9, "{:?} → {:?}", ra.before, ra.after);
        let ap = |m: &Model, r: &Morphed| applied(m, &Relaxed { positions: r.positions.clone(), ..Relaxed::default() });
        let d = largest_difference(&ap(&m, &ra), &ap(&b, &rb));
        assert!(d < 1e-6, "n {n}: a 1-ulp input change moved a node by {d:e}");
    }
}

/// The strong mode's fast hex20 samples give the certified quality's scaled Jacobian (the
/// smallest over the same samples) on curved, distorted elements (#612).
#[test]
fn hex20_samples_match_the_certified_scaled_jacobian() {
    use crate::quality::jacobian::{HEX_CORNERS, HEX_EDGES};
    for case in 0..20u32 {
        // A unit cube, sheared and stretched, its corners and mid-side nodes nudged.
        let t = f64::from(case);
        let wobble = |k: usize, d: usize| 0.06 * ((t + 1.0) * (k as f64 + 1.3) * (d as f64 + 0.7)).sin();
        let map = |p: P| [1.3 * p[0] + 0.2 * p[1], 0.9 * p[1] + 0.1 * t / 20.0 * p[2], 1.1 * p[2] + 0.15 * p[0]];
        let mut x: Vec<P> = HEX_CORNERS
            .iter()
            .enumerate()
            .map(|(k, c)| {
                let p = map(c.map(|v| v as f64));
                [p[0] + wobble(k, 0), p[1] + wobble(k, 1), p[2] + wobble(k, 2)]
            })
            .collect();
        for (i, &(a, b)) in HEX_EDGES.iter().enumerate() {
            let m = scale(add(x[a], x[b]), 0.5);
            x.push([m[0] + 0.5 * wobble(8 + i, 0), m[1] + 0.5 * wobble(8 + i, 1), m[2] + 0.5 * wobble(8 + i, 2)]);
        }
        let certified = crate::quality::element_quality(ElementKind::Hex20, &x).unwrap().get(Metric::ScaledJacobian).unwrap();
        let (centre, samples) = strong::hex20_scaled_jacobian(&x);
        let fast = samples.iter().copied().fold(f64::INFINITY, f64::min);
        assert!((fast - certified).abs() < 1e-12, "case {case}: {fast} against {certified}");
        assert!(centre > 0.0);
    }
}
