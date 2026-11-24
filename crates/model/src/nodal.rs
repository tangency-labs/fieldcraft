//! Consistent nodal loads: turning a force, moment or traction on a set into the nodal values a
//! deck or the solver applies, the same way for every consumer.
//!
//! A load over a surface is spread by the element shape functions: node `i` of a face takes
//! `∫ N_i dA` of a unit traction, integrated by Gauss quadrature on the face's own (possibly
//! curved, quadratic) geometry. The weights of a face sum to its area, so a total force split by
//! them keeps its total on any mesh of the same surface. Quadratic faces give the textbook
//! values (corner 0 and mid-side A/3 on a flat six-node triangle; corner −A/12 and mid-side A/3 on
//! a flat eight-node quadrilateral).

use std::collections::BTreeMap;

use crate::{ElementEdge, ElementFace, Error, Extent, Model, Result};

/// Shape functions and their derivatives at one point of a reference face.
struct Shape {
    n: Vec<f64>,
    dxi: Vec<f64>,
    deta: Vec<f64>,
}

/// Triangle: barycentric L0 = 1 − ξ − η, L1 = ξ, L2 = η.
fn triangle(nodes: usize, xi: f64, eta: f64) -> Shape {
    let l = [1.0 - xi - eta, xi, eta];
    let dl = [(-1.0, -1.0), (1.0, 0.0), (0.0, 1.0)];
    let mut s = Shape { n: Vec::with_capacity(nodes), dxi: Vec::with_capacity(nodes), deta: Vec::with_capacity(nodes) };
    if nodes == 3 {
        for (li, (dx, de)) in l.iter().zip(dl) {
            s.n.push(*li);
            s.dxi.push(dx);
            s.deta.push(de);
        }
        return s;
    }
    for (li, (dx, de)) in l.iter().zip(dl) {
        s.n.push(li * (2.0 * li - 1.0));
        s.dxi.push((4.0 * li - 1.0) * dx);
        s.deta.push((4.0 * li - 1.0) * de);
    }
    for (la, lb, (dxa, dea), (dxb, deb)) in [(l[0], l[1], dl[0], dl[1]), (l[1], l[2], dl[1], dl[2]), (l[2], l[0], dl[2], dl[0])] {
        s.n.push(4.0 * la * lb);
        s.dxi.push(4.0 * (dxa * lb + la * dxb));
        s.deta.push(4.0 * (dea * lb + la * deb));
    }
    s
}

/// Quadrilateral on [−1, 1]², corners (−1,−1), (1,−1), (1,1), (−1,1).
fn quadrilateral(nodes: usize, xi: f64, eta: f64) -> Shape {
    const CORNERS: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let mut s = Shape { n: Vec::with_capacity(nodes), dxi: Vec::with_capacity(nodes), deta: Vec::with_capacity(nodes) };
    for (xi_i, eta_i) in CORNERS {
        let (a, b) = (1.0 + xi * xi_i, 1.0 + eta * eta_i);
        if nodes == 4 {
            s.n.push(a * b / 4.0);
            s.dxi.push(xi_i * b / 4.0);
            s.deta.push(eta_i * a / 4.0);
        } else {
            let c = xi * xi_i + eta * eta_i - 1.0;
            s.n.push(a * b * c / 4.0);
            s.dxi.push(xi_i * b * (2.0 * xi * xi_i + eta * eta_i) / 4.0);
            s.deta.push(eta_i * a * (xi * xi_i + 2.0 * eta * eta_i) / 4.0);
        }
    }
    if nodes == 8 {
        // Mid-sides of edges 1-2 (η = −1), 2-3 (ξ = 1), 3-4 (η = 1), 4-1 (ξ = −1).
        for (xi_i, eta_i) in [(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)] {
            if xi_i == 0.0 {
                s.n.push((1.0 - xi * xi) * (1.0 + eta * eta_i) / 2.0);
                s.dxi.push(-xi * (1.0 + eta * eta_i));
                s.deta.push(eta_i * (1.0 - xi * xi) / 2.0);
            } else {
                s.n.push((1.0 + xi * xi_i) * (1.0 - eta * eta) / 2.0);
                s.dxi.push(xi_i * (1.0 - eta * eta) / 2.0);
                s.deta.push(-eta * (1.0 + xi * xi_i));
            }
        }
    }
    s
}

/// Gauss points `(ξ, η, weight)`: a degree-4 six-point rule on the triangle (Dunavant 1985) and
/// 3 × 3 Gauss–Legendre on the quadrilateral; both integrate the quadratic faces' weights exactly
/// when the face is flat.
fn rule(corners: usize) -> Vec<(f64, f64, f64)> {
    if corners == 3 {
        let (a1, w1) = (0.445_948_490_915_965, 0.223_381_589_678_011 / 2.0);
        let (a2, w2) = (0.091_576_213_509_771, 0.109_951_743_655_322 / 2.0);
        let (b1, b2) = (1.0 - 2.0 * a1, 1.0 - 2.0 * a2);
        vec![(a1, a1, w1), (a1, b1, w1), (b1, a1, w1), (a2, a2, w2), (a2, b2, w2), (b2, a2, w2)]
    } else {
        let g = (0.6f64).sqrt();
        let pts = [(-g, 5.0 / 9.0), (0.0, 8.0 / 9.0), (g, 5.0 / 9.0)];
        pts.iter().flat_map(|&(x, wx)| pts.iter().map(move |&(y, wy)| (x, y, wx * wy))).collect()
    }
}

/// Visit the quadrature points of a face (node coordinates in the face's node order): the shape
/// function values and the weighted normal (∂x/∂ξ × ∂x/∂η) w, whose length is the point's share of
/// the area and whose direction is the face's right-hand normal. Every surface load (pressure,
/// traction, heat flux, convection) integrates with it, so they share one rule.
pub fn integrate_face(xyz: &[[f64; 3]], mut visit: impl FnMut(&[f64], [f64; 3])) -> Result<()> {
    if xyz.iter().flatten().any(|c| !c.is_finite()) {
        return Err(Error::Invalid("face coordinates must be finite".into()));
    }
    let (corners, shape): (usize, fn(usize, f64, f64) -> Shape) = match xyz.len() {
        3 | 6 => (3, triangle),
        4 | 8 => (4, quadrilateral),
        n => return Err(Error::Invalid(format!("a face has 3, 4, 6 or 8 nodes, not {n}"))),
    };
    for (xi, eta, w) in rule(corners) {
        let s = shape(xyz.len(), xi, eta);
        let mut tx = [0.0; 3];
        let mut te = [0.0; 3];
        for ((p, dx), de) in xyz.iter().zip(&s.dxi).zip(&s.deta) {
            for k in 0..3 {
                tx[k] += dx * p[k];
                te[k] += de * p[k];
            }
        }
        let c = [tx[1] * te[2] - tx[2] * te[1], tx[2] * te[0] - tx[0] * te[2], tx[0] * te[1] - tx[1] * te[0]];
        visit(&s.n, c.map(|v| v * w));
    }
    Ok(())
}

/// A quadrature point of a reference face or edge: the shape function values, their
/// derivatives along the face's natural coordinates ξ and η (an edge has ξ only; its `deta` is
/// empty) and the weight.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraturePoint {
    pub n: Vec<f64>,
    pub dxi: Vec<f64>,
    pub deta: Vec<f64>,
    pub weight: f64,
}

/// The quadrature points of a face of `nodes` nodes (3, 4, 6 or 8), with the rule and shape
/// functions of [`integrate_face`]: what a load on the deformed face (follower pressure)
/// integrates with.
pub fn face_points(nodes: usize) -> Result<Vec<QuadraturePoint>> {
    let (corners, shape): (usize, fn(usize, f64, f64) -> Shape) = match nodes {
        3 | 6 => (3, triangle),
        4 | 8 => (4, quadrilateral),
        n => return Err(Error::Invalid(format!("a face has 3, 4, 6 or 8 nodes, not {n}"))),
    };
    Ok(rule(corners)
        .into_iter()
        .map(|(xi, eta, weight)| {
            let s = shape(nodes, xi, eta);
            QuadraturePoint { n: s.n, dxi: s.dxi, deta: s.deta, weight }
        })
        .collect())
}

/// The quadrature points of an element edge of `nodes` nodes (2 or 3), with the rule and shape
/// functions of [`edge_normal_weights`].
pub fn edge_points(nodes: usize) -> Result<Vec<QuadraturePoint>> {
    if !matches!(nodes, 2 | 3) {
        return Err(Error::Invalid(format!("an edge has 2 or 3 nodes, not {nodes}")));
    }
    let g = (0.6f64).sqrt();
    Ok([(-g, 5.0 / 9.0), (0.0, 8.0 / 9.0), (g, 5.0 / 9.0)]
        .into_iter()
        .map(|(s, weight)| {
            let (n, dn) = if nodes == 2 {
                (vec![(1.0 - s) / 2.0, (1.0 + s) / 2.0], vec![-0.5, 0.5])
            } else {
                (vec![s * (s - 1.0) / 2.0, s * (s + 1.0) / 2.0, 1.0 - s * s], vec![s - 0.5, s + 0.5, -2.0 * s])
            };
            QuadraturePoint { n, dxi: dn, deta: Vec::new(), weight }
        })
        .collect())
}

/// ∫ Nᵢ n dA for every node of a face: the consistent nodal forces of a unit pressure pulling
/// along the face's right-hand (outward, for solids) normal.
pub fn face_normal_weights(xyz: &[[f64; 3]]) -> Result<Vec<[f64; 3]>> {
    let mut out = vec![[0.0; 3]; xyz.len()];
    integrate_face(xyz, |n, nda| {
        for (o, ni) in out.iter_mut().zip(n) {
            for k in 0..3 {
                o[k] += ni * nda[k];
            }
        }
    })?;
    Ok(out)
}

/// Visit the 3-point Gauss–Legendre points of an element edge (its ends, then the mid-side node
/// of a quadratic edge): the shape function values, the position and the weighted tangent
/// (∂x/∂s) w, whose length is the point's share of the edge's length.
fn integrate_edge(xyz: &[[f64; 3]], mut visit: impl FnMut(&[f64], [f64; 3], [f64; 3])) -> Result<()> {
    if xyz.iter().flatten().any(|c| !c.is_finite()) {
        return Err(Error::Invalid("edge coordinates must be finite".into()));
    }
    if !matches!(xyz.len(), 2 | 3) {
        return Err(Error::Invalid(format!("an edge has 2 or 3 nodes, not {}", xyz.len())));
    }
    let g = (0.6f64).sqrt();
    for (s, w) in [(-g, 5.0 / 9.0), (0.0, 8.0 / 9.0), (g, 5.0 / 9.0)] {
        let (n, dn) = if xyz.len() == 2 {
            ([(1.0 - s) / 2.0, (1.0 + s) / 2.0, 0.0], [-0.5, 0.5, 0.0])
        } else {
            ([s * (s - 1.0) / 2.0, s * (s + 1.0) / 2.0, 1.0 - s * s], [s - 0.5, s + 0.5, -2.0 * s])
        };
        let (mut x, mut t) = ([0.0; 3], [0.0; 3]);
        for ((p, ni), di) in xyz.iter().zip(&n).zip(&dn) {
            for k in 0..3 {
                x[k] += ni * p[k];
                t[k] += di * p[k];
            }
        }
        visit(n.get(..xyz.len()).unwrap_or_default(), x, t.map(|v| v * w));
    }
    Ok(())
}

/// ∫ Nᵢ n e dL for every node of an element edge in the XY plane (node coordinates in the edge's
/// order), with n the normal to the right of the edge's direction (outward for an edge of a
/// counter-clockwise 2D element) and e the out-of-plane extent (thickness, or 2πr when
/// revolved): the consistent nodal forces of a unit pressure pulling outward.
pub fn edge_normal_weights(xyz: &[[f64; 3]], extent: Extent) -> Result<Vec<[f64; 3]>> {
    let mut out = vec![[0.0; 3]; xyz.len()];
    integrate_edge(xyz, |n, x, t| {
        let e = extent.at(x);
        for (o, ni) in out.iter_mut().zip(n) {
            o[0] += ni * t[1] * e;
            o[1] -= ni * t[0] * e;
        }
    })?;
    Ok(out)
}

/// Consistent nodal forces of a uniform pressure on an element edge in the XY plane (positive
/// pushes into a counter-clockwise element): `−p` times [`edge_normal_weights`].
pub fn edge_pressure(xyz: &[[f64; 3]], extent: Extent, magnitude: f64) -> Result<Vec<[f64; 3]>> {
    Ok(edge_normal_weights(xyz, extent)?.into_iter().map(|w| w.map(|c| -magnitude * c)).collect())
}

impl Model {
    /// The nodal weights `∫ N_i e dL` of a set of 2D element edges (e the out-of-plane extent of
    /// each element in `planar`), summed per node. Their sum is the area of the surface the
    /// edges stand for: length times thickness, or the surface of revolution.
    pub fn edge_weights(&self, edges: &[ElementEdge], planar: &BTreeMap<u32, Extent>) -> Result<BTreeMap<u32, f64>> {
        let mut out: BTreeMap<u32, f64> = BTreeMap::new();
        for &g in edges {
            let extent = planar.get(&g.element).copied().ok_or_else(|| Error::Invalid(format!("element {} is not a 2D element", g.element)))?;
            let labels = self.edge_nodes(g).ok_or_else(|| Error::NotFound(format!("element {} edge {}", g.element, g.edge)))?;
            let xyz: Vec<[f64; 3]> =
                labels.iter().map(|n| self.node(*n).map(|p| p.xyz).ok_or_else(|| Error::NotFound(format!("node {n}")))).collect::<Result<_>>()?;
            integrate_edge(&xyz, |n, x, t| {
                let dl = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt() * extent.at(x);
                for (label, ni) in labels.iter().zip(n) {
                    *out.entry(*label).or_default() += ni * dl;
                }
            })
            .map_err(|e| match e {
                Error::Invalid(m) => Error::Invalid(format!("element {} edge {}: {m}", g.element, g.edge)),
                e => e,
            })?;
        }
        Ok(out)
    }

    /// The nodal weights `∫ N_i r e dL` of a set of 2D element edges of an axisymmetric model
    /// (r = x, e = 2πr), summed per node: the consistent nodal hoop forces of a hoop shear
    /// traction equal to the radius. Their sum weighted by the nodes' radii is the polar moment
    /// `∫ r² dA` of the surface the edges stand for.
    pub fn edge_torque_weights(&self, edges: &[ElementEdge], planar: &BTreeMap<u32, Extent>) -> Result<BTreeMap<u32, f64>> {
        let mut out: BTreeMap<u32, f64> = BTreeMap::new();
        for &g in edges {
            let extent = planar.get(&g.element).copied().ok_or_else(|| Error::Invalid(format!("element {} is not a 2D element", g.element)))?;
            let labels = self.edge_nodes(g).ok_or_else(|| Error::NotFound(format!("element {} edge {}", g.element, g.edge)))?;
            let xyz: Vec<[f64; 3]> =
                labels.iter().map(|n| self.node(*n).map(|p| p.xyz).ok_or_else(|| Error::NotFound(format!("node {n}")))).collect::<Result<_>>()?;
            integrate_edge(&xyz, |n, x, t| {
                let dl = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt() * extent.at(x);
                for (label, ni) in labels.iter().zip(n) {
                    *out.entry(*label).or_default() += ni * x[0] * dl;
                }
            })
            .map_err(|e| match e {
                Error::Invalid(m) => Error::Invalid(format!("element {} edge {}: {m}", g.element, g.edge)),
                e => e,
            })?;
        }
        Ok(out)
    }

    /// The nodal weights `∫ N_i dA` of a set of element faces, summed per node. Their sum is the
    /// faces' area.
    pub fn surface_weights(&self, faces: &[ElementFace]) -> Result<BTreeMap<u32, f64>> {
        let mut out: BTreeMap<u32, f64> = BTreeMap::new();
        for &f in faces {
            let labels = self.face_nodes(f).ok_or_else(|| Error::NotFound(format!("element {} face {}", f.element, f.face)))?;
            let xyz: Vec<[f64; 3]> =
                labels.iter().map(|n| self.node(*n).map(|p| p.xyz).ok_or_else(|| Error::NotFound(format!("node {n}")))).collect::<Result<_>>()?;
            integrate_face(&xyz, |n, nda| {
                let da = (nda[0] * nda[0] + nda[1] * nda[1] + nda[2] * nda[2]).sqrt();
                for (label, ni) in labels.iter().zip(n) {
                    *out.entry(*label).or_default() += ni * da;
                }
            })
            .map_err(|e| match e {
                Error::Invalid(m) => Error::Invalid(format!("element {} face {}: {m}", f.element, f.face)),
                e => e,
            })?;
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ElementKind;
    use std::collections::BTreeMap;

    /// A flat 2 × 3 face of a hex8 or hex20 (face 1, z = 0), and its weights.
    fn face_weights(kind: ElementKind) -> Vec<f64> {
        let mut m = Model::default();
        let corners = [[0., 0., 0.], [2., 0., 0.], [2., 3., 0.], [0., 3., 0.], [0., 0., 1.], [2., 0., 1.], [2., 3., 1.], [0., 3., 1.]];
        let mut xyz = corners.to_vec();
        for e in kind.edges() {
            if let [a, b, _] = **e {
                xyz.push([0, 1, 2].map(|i| (corners[a][i] + corners[b][i]) / 2.0));
            }
        }
        for (i, p) in xyz.iter().enumerate() {
            m.set_node(i as u32 + 1, *p).unwrap();
        }
        m.set_element(1, kind, (1..=kind.node_count() as u32).collect()).unwrap();
        let w = m.surface_weights(&[ElementFace { element: 1, face: 1 }]).unwrap();
        m.face_nodes(ElementFace { element: 1, face: 1 }).unwrap().iter().map(|n| w[n]).collect()
    }

    #[test]
    fn weights_are_the_consistent_nodal_loads() {
        let close = |a: &[f64], b: &[f64]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-12);
        assert!(close(&face_weights(ElementKind::Hex8), &[1.5; 4]), "{:?}", face_weights(ElementKind::Hex8));
        let quad8 = face_weights(ElementKind::Hex20);
        assert!(close(&quad8, &[-0.5, -0.5, -0.5, -0.5, 2.0, 2.0, 2.0, 2.0]), "{quad8:?}");

        // Triangles: a tet4 and a tet10 face (face 1, the right triangle of legs 1 in z = 0).
        let mut m = Model::default();
        let p = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let mut xyz = p.to_vec();
        for e in ElementKind::Tet10.edges() {
            if let [a, b, _] = **e {
                xyz.push([0, 1, 2].map(|i| (p[a][i] + p[b][i]) / 2.0));
            }
        }
        for (i, q) in xyz.iter().enumerate() {
            m.set_node(i as u32 + 1, *q).unwrap();
        }
        m.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 4]).unwrap();
        m.set_element(2, ElementKind::Tet10, (1..=10).collect()).unwrap();
        let w = m.surface_weights(&[ElementFace { element: 1, face: 1 }]).unwrap();
        assert!(close(&w.values().copied().collect::<Vec<_>>(), &[1.0 / 6.0; 3]), "{w:?}");
        let w = m.surface_weights(&[ElementFace { element: 2, face: 1 }]).unwrap();
        let tri6: Vec<f64> = m.face_nodes(ElementFace { element: 2, face: 1 }).unwrap().iter().map(|n| w[n]).collect();
        assert!(close(&tri6, &[0.0, 0.0, 0.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]), "{tri6:?}");
    }

    #[test]
    fn edge_weights_are_the_consistent_nodal_loads() {
        let close = |a: &[f64], b: &[f64]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-12);
        // A quad8 2 × 1 in the XY plane: edge 1 (y = 0, x from 0 to 2) of thickness 3 shares
        // its area 6 as 1/6, 1/6, 2/3; revolved about Y it shares the annulus π(2² − 0) = 4π
        // as ∫ N 2π x dx.
        let mut m = Model::default();
        let p = [[0., 0., 0.], [2., 0., 0.], [2., 1., 0.], [0., 1., 0.], [1., 0., 0.], [2., 0.5, 0.], [1., 1., 0.], [0., 0.5, 0.]];
        for (i, q) in p.iter().enumerate() {
            m.set_node(i as u32 + 1, *q).unwrap();
        }
        m.set_element(1, ElementKind::Quad8, (1..=8).collect()).unwrap();
        let edge = [ElementEdge { element: 1, edge: 1 }];
        let w = m.edge_weights(&edge, &BTreeMap::from([(1, Extent::Thickness(3.0))])).unwrap();
        assert!(close(&[w[&1], w[&2], w[&5]], &[1.0, 1.0, 4.0]), "{w:?}");
        let w = m.edge_weights(&edge, &BTreeMap::from([(1, Extent::Revolved)])).unwrap();
        // ∫₀² N 2π x dx with N on [0, 2]: corner at 0 → 0, corner at 2 → 2π·2/3·... summed: 4π.
        let pi = std::f64::consts::PI;
        assert!(close(&[w[&1], w[&2], w[&5]], &[0.0, 4.0 * pi / 3.0, 8.0 * pi / 3.0]), "{w:?}");
        assert!(m.edge_weights(&edge, &BTreeMap::new()).is_err(), "not a 2D element");
        // The outward normal of a counter-clockwise element's bottom edge is −Y.
        let n = edge_normal_weights(&[p[0], p[1], p[4]], Extent::Thickness(1.0)).unwrap();
        assert!(close(&n.iter().map(|v| v[1]).collect::<Vec<_>>(), &[-1.0 / 3.0, -1.0 / 3.0, -4.0 / 3.0]) && n.iter().all(|v| v[0] == 0.0), "{n:?}");
        assert!(edge_normal_weights(&[p[0]], Extent::Revolved).is_err());
    }
}
