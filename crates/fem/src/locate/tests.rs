use super::*;

const SOLIDS: [ElementKind; 6] =
    [ElementKind::Tet4, ElementKind::Tet10, ElementKind::Wedge6, ElementKind::Pyramid5, ElementKind::Hex8, ElementKind::Hex20];

fn mapped(kind: ElementKind, xyz: &[[f64; 3]], natural: [f64; 3]) -> [f64; 3] {
    let shape = Basis::of(kind).eval(natural);
    [0, 1, 2].map(|k| xyz.iter().zip(shape.values()).map(|(x, n)| x[k] * n).sum())
}

fn linear(x: [f64; 3]) -> f64 {
    2.0 * x[0] + 3.0 * x[1] - x[2] + 7.0
}

#[test]
fn linear_patch_all_six_solids_and_every_node() {
    for kind in SOLIDS {
        // Nonsymmetric, rotated/sheared affine map detects transposed inverse mistakes.
        let xyz: Vec<_> = Basis::of(kind)
            .nodes()
            .iter()
            .map(|[r, s, t]| [3.0 + 2.0 * r + 0.4 * s, -7.0 + 0.2 * r + 3.0 * s + 0.5 * t, 2.0 + 0.3 * r + 0.7 * s + 4.0 * t])
            .collect();
        let values: Vec<_> = xyz.iter().copied().map(linear).collect();
        for natural in std::iter::once(center(kind)).chain(Basis::of(kind).nodes().iter().copied()) {
            let point = mapped(kind, &xyz, natural);
            let found = locate(kind, &xyz, point, &mut Work::new(10_000)).unwrap().unwrap();
            let actual = found.interpolate(&values).unwrap();
            assert!((actual - linear(point)).abs() < 1e-12 * linear(point).abs().max(1.0), "{kind:?} {natural:?}: {actual}");
        }
        for (i, point) in xyz.iter().enumerate() {
            let found = locate(kind, &xyz, *point, &mut Work::new(10_000)).unwrap().unwrap();
            assert_eq!(found.interpolate(&values), Some(values[i]), "exact node {kind:?} {i}");
        }
    }
}

#[test]
fn curved_quadratic_maps_extend_beyond_all_nodes_but_are_found() {
    for kind in [ElementKind::Tet10, ElementKind::Hex20] {
        let xyz: Vec<_> = Basis::of(kind)
            .nodes()
            .iter()
            .map(|[r, s, t]| match kind {
                ElementKind::Tet10 => [*r, r + s - 2.4 * r * (1.0 - r - s - t), *t],
                _ => [*r, s + r - (1.0 - r * r), *t],
            })
            .collect();
        let xi = if kind == ElementKind::Tet10 { [0.25, 0.01, 0.01] } else { [-0.5, -0.99, 0.0] };
        let point = mapped(kind, &xyz, xi);
        assert!(point[1] < xyz.iter().map(|x| x[1]).fold(f64::INFINITY, f64::min));
        let enclosure = bounds(kind, &xyz).unwrap();
        assert!(enclosure.contains(point), "{kind:?}: {enclosure:?} misses {point:?}");
        let found = locate(kind, &xyz, point, &mut Work::new(10_000)).unwrap().unwrap();
        assert!((found.interpolate(&xyz.iter().copied().map(linear).collect::<Vec<_>>()).unwrap() - linear(point)).abs() < 1e-12);
        // A deterministic interior grid checks enclosure and inversion beyond the one regression.
        for i in 0..7 {
            for j in 0..7 {
                for k in 0..7 {
                    let xi = if kind == ElementKind::Tet10 {
                        [i as f64 / 20.0, j as f64 / 20.0, k as f64 / 20.0]
                    } else {
                        [i as f64 / 3.0 - 1.0, j as f64 / 3.0 - 1.0, k as f64 / 3.0 - 1.0]
                    };
                    let point = mapped(kind, &xyz, xi);
                    assert!(enclosure.contains(point), "{kind:?} {xi:?}");
                    let found = locate(kind, &xyz, point, &mut Work::new(10_000)).unwrap().unwrap();
                    assert!((0..3).all(|a| (found.natural[a] - xi[a]).abs() < 1e-11), "{kind:?} {xi:?} -> {:?}", found.natural);
                }
            }
        }
    }
}

#[test]
fn outside_requires_reference_exclusion_and_exhaustion_is_an_error() {
    for kind in SOLIDS {
        let xyz = Basis::of(kind).nodes();
        assert!(locate(kind, xyz, [9.0; 3], &mut Work::new(10_000)).unwrap().is_none());
        assert!(locate(kind, xyz, center(kind), &mut Work::new(0)).unwrap_err().to_string().contains("work limit"));
    }
    // Inside the tetrahedron's coordinate AABB but outside its volume, and similarly
    // for the wedge/pyramid. Returning null here must follow actual exclusion.
    for (kind, point) in [
        (ElementKind::Tet4, [0.6, 0.6, 0.6]),
        (ElementKind::Tet10, [0.6, 0.6, 0.6]),
        (ElementKind::Wedge6, [0.8, 0.8, 0.0]),
        (ElementKind::Pyramid5, [0.9, 0.9, 0.7]),
    ] {
        assert!(locate(kind, Basis::of(kind).nodes(), point, &mut Work::new(10_000)).unwrap().is_none(), "{kind:?}");
    }
    // Subdivision can conservatively refuse a point arbitrarily near a curved boundary.
    // Its finite cap is explicit and cannot turn an unresolved point into an outside gap.
    let kind = ElementKind::Hex20;
    let xyz: Vec<_> = Basis::of(kind).nodes().iter().map(|[r, s, t]| [*r, s + r - (1.0 - r * r), *t]).collect();
    let point = [-0.5, -2.250000000001, 0.0];
    let result = locate(kind, &xyz, point, &mut Work::new(10_000));
    assert!(result.is_err(), "near-boundary search must disclose uncertainty: {result:?}");
    assert!(result.unwrap_err().to_string().contains("unresolved"));
}

#[test]
fn scaled_translated_and_distorted_elements_preserve_linear_interpolation() {
    for kind in SOLIDS {
        for scale in [1e-12, 1.0, 1e12] {
            let xyz: Vec<_> = Basis::of(kind)
                .nodes()
                .iter()
                .map(|[r, s, t]| [scale * (3.0 + r + 0.12 * s * t), scale * (-2.0 + s + 0.06 * r * t), scale * (4.0 + t)])
                .collect();
            let xi = center(kind);
            let point = mapped(kind, &xyz, xi);
            let found = locate(kind, &xyz, point, &mut Work::new(10_000)).unwrap().unwrap();
            assert!((0..3).all(|a| (found.natural[a] - xi[a]).abs() < 1e-12), "{kind:?} {scale}");
        }
    }
}

#[test]
fn invalid_geometry_data_and_unsupported_families_refuse_without_panics() {
    for kind in SOLIDS {
        let mut xyz = Basis::of(kind).nodes().to_vec();
        assert!(bounds(kind, &xyz[..1]).is_err());
        xyz[0][0] = f64::NAN;
        assert!(bounds(kind, &xyz).is_err());
        assert!(locate(kind, &xyz, [0.0; 3], &mut Work::new(10_000)).is_err());
        let collapsed = vec![[0.0; 3]; Basis::of(kind).count()];
        assert!(locate(kind, &collapsed, [0.0; 3], &mut Work::new(10_000)).is_err());
        let inverted: Vec<_> = Basis::of(kind).nodes().iter().map(|[r, s, t]| [-r, *s, *t]).collect();
        assert!(locate(kind, &inverted, mapped(kind, &inverted, center(kind)), &mut Work::new(10_000)).is_err());
    }
    assert!(bounds(ElementKind::Tri3, Basis::Tri3.nodes()).is_err());
    assert!(bounds(ElementKind::Tet4, &[[f64::MAX, 0.0, 0.0], [-f64::MAX, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).is_err());
    let location = Location { natural: [0.25; 3], shape: Basis::Tet4.eval([0.25; 3]) };
    assert!(location.interpolate(&[1.0, 2.0]).is_none());
    assert!(location.interpolate(&[1.0, f64::NAN, 3.0, 4.0]).is_none());
    assert!(location.interpolate(&[f64::INFINITY; 4]).is_none());
}

#[test]
fn exact_node_does_not_bypass_the_local_jacobian_check() {
    // Positive center determinant 1-3/4, but the node r=1 has determinant -2.
    let kind = ElementKind::Tet10;
    let xyz: Vec<_> = Basis::of(kind).nodes().iter().map(|[r, s, t]| [*r, s * (1.0 - 3.0 * r), *t]).collect();
    let error = locate(kind, &xyz, xyz[1], &mut Work::new(10_000)).unwrap_err();
    assert!(error.to_string().contains("queried node"), "{error}");
}

#[test]
fn bounds_cost_is_debited_before_any_basis_evaluation() {
    let kind = ElementKind::Hex20;
    let xyz = Basis::of(kind).nodes();
    let mut insufficient = Work::new(27);
    assert!(locate(kind, xyz, [9.0; 3], &mut insufficient).unwrap_err().to_string().contains("work limit"));
    assert_eq!(insufficient.remaining(), 27);
    let mut exact = Work::new(28);
    assert!(locate(kind, xyz, [9.0; 3], &mut exact).unwrap().is_none());
    assert_eq!(exact.remaining(), 0);
    let mut no_center = Work::new(28);
    assert!(locate(kind, xyz, [0.0; 3], &mut no_center).unwrap_err().to_string().contains("work limit"));
    assert_eq!(no_center.remaining(), 0);
}

#[test]
fn candidate_roundoff_band_does_not_drop_exact_boundary_stations() {
    let kind = ElementKind::Hex8;
    let xyz: Vec<_> = Basis::of(kind).nodes().iter().map(|p| p.map(|v| (v + 1.0) * 0.0005)).collect();
    let point = [-1e-19, 0.0005, 0.0005];
    let enclosure = bounds(kind, &xyz).unwrap();
    assert!(!enclosure.contains(point));
    assert!(enclosure.contains_near(point));
    assert!(locate(kind, &xyz, point, &mut Work::new(100)).unwrap().is_some());
    assert!(locate(kind, &xyz, [-1e-8, 0.0005, 0.0005], &mut Work::new(100)).unwrap().is_none());
}

#[test]
fn undefined_values_with_exact_zero_weights_do_not_erase_a_defined_node() {
    let location = Location { natural: [0.0; 3], shape: Basis::Tet4.eval([0.0; 3]) };
    assert_eq!(location.interpolate(&[4.0, f64::NAN, f64::NAN, f64::NAN]), Some(4.0));
    assert!(location.interpolate(&[f64::NAN, 4.0, 4.0, 4.0]).is_none());
    let face = Location { natural: [0.25, 0.25, 0.0], shape: Basis::Tet4.eval([0.25, 0.25, 0.0]) };
    assert_eq!(face.interpolate(&[2.0, 4.0, 8.0, f64::NAN]), Some(4.0));
    assert!(face.interpolate(&[2.0, f64::NAN, 8.0, f64::NAN]).is_none());
}

#[test]
fn pyramid_apex_edges_faces_and_near_apex_are_not_a_cube_domain() {
    let kind = ElementKind::Pyramid5;
    let xyz: Vec<_> = Basis::of(kind).nodes().iter().map(|[r, s, t]| [2.0 + r + 0.1 * s, 3.0 + s, 4.0 + 2.0 * t]).collect();
    for xi in [[0.0, 0.0, 1.0], [0.2, -0.2, 0.8], [0.0, 0.1, 0.9], [1e-5, -1e-5, 1.0 - 2e-5], [0.4, 0.2, 0.0]] {
        let point = mapped(kind, &xyz, xi);
        let found = locate(kind, &xyz, point, &mut Work::new(10_000)).unwrap().unwrap();
        assert!((0..3).all(|a| (found.natural[a] - xi[a]).abs() < 1e-12), "{xi:?} -> {:?}", found.natural);
    }
    assert!(!inside(kind, [0.5, 0.5, 0.9], INSIDE));
}

#[test]
fn large_world_origin_does_not_enter_the_newton_jacobian() {
    for kind in SOLIDS {
        let xyz: Vec<_> =
            Basis::of(kind).nodes().iter().map(|[r, s, t]| [1e9 + 2.0 * r + 0.1 * s, -2e9 + 3.0 * s + 0.2 * t, 3e9 + 4.0 * t]).collect();
        let point = mapped(kind, &xyz, center(kind));
        let found = locate(kind, &xyz, point, &mut Work::new(10_000)).unwrap().unwrap();
        let values: Vec<_> = xyz.iter().copied().map(linear).collect();
        assert!((found.interpolate(&values).unwrap() - linear(point)).abs() < 1e-12 * linear(point).abs());
    }
}

#[test]
fn outward_bounds_union_and_interval_extensions_enclose_samples() {
    for kind in SOLIDS {
        let xyz = Basis::of(kind).nodes();
        let a = bounds(kind, xyz).unwrap();
        let shifted: Vec<_> = xyz.iter().map(|x| x.map(|v| v + 3.0)).collect();
        let b = bounds(kind, &shifted).unwrap();
        let both = a.union(b);
        assert!(xyz.iter().chain(&shifted).all(|x| both.contains(*x)));
        let xi = center(kind);
        let interval = interval_weights(kind, xi.map(Interval::point));
        let actual = if kind == ElementKind::Pyramid5 {
            Basis::of(kind).eval([xi[0] * (1.0 - xi[2]), xi[1] * (1.0 - xi[2]), xi[2]])
        } else {
            Basis::of(kind).eval(xi)
        };
        for (n, extension) in actual.values().iter().zip(interval) {
            assert!(*n >= extension.lo && *n <= extension.hi, "{kind:?} {n} {extension:?}");
        }
    }
}
