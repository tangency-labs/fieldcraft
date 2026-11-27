//! Display externally computed point and cell arrays on the grid that supplied them.
use crate::results::{DEFAULT_BANDS, MAX_DEFORMATION_SCALE};
use crate::view::{Palette, Section, Show, Shown, Style};
pub use fieldcraft_io::external::{Cell, Field, Frame, Grid, Series};
use fieldcraft_model::ElementKind;
use fieldcraft_render::contour::{Bands, banded, bands_color, split_triangle};
use fieldcraft_render::{Aabb3, DisplayMesh, DisplayTri, Scene, SceneLine, ScenePoint, SceneTri, Vec3};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Association {
    Point,
    Cell,
}

impl Association {
    pub fn fields(self, frame: &Frame) -> &BTreeMap<String, Field> {
        match self {
            Self::Point => &frame.point_data,
            Self::Cell => &frame.cell_data,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Point => "Point",
            Self::Cell => "Cell",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Component {
    Value,
    Magnitude,
    X,
    Y,
    Z,
}

impl Component {
    pub fn label(self) -> &'static str {
        match self {
            Self::Value => "Value",
            Self::Magnitude => "Magnitude",
            Self::X => "X",
            Self::Y => "Y",
            Self::Z => "Z",
        }
    }

    pub fn for_components(components: u8) -> &'static [Self] {
        match components {
            1 => &[Self::Value],
            2 => &[Self::Magnitude, Self::X, Self::Y],
            3 => &[Self::Magnitude, Self::X, Self::Y, Self::Z],
            _ => &[],
        }
    }

    fn sample(self, field: &Field, index: usize) -> Option<f64> {
        let width = usize::from(field.components);
        let start = index.checked_mul(width)?;
        let tuple = field.values.get(start..start.checked_add(width)?)?;
        let v = match (self, tuple) {
            (Self::Value, [v]) => *v,
            (Self::Magnitude, [x, y]) => libm::hypot(*x, *y),
            (Self::Magnitude, [x, y, z]) => libm::hypot(libm::hypot(*x, *y), *z),
            (Self::X, [x, _, ..]) => *x,
            (Self::Y, [_, y, ..]) => *y,
            (Self::Z, [_, _, z]) => *z,
            _ => return None,
        };
        v.is_finite().then_some(v)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExternalView {
    pub frame: usize,
    pub association: Association,
    pub field: Option<String>,
    pub component: Component,
    pub displacement: Option<String>,
    pub scale: f64,
    pub show: bool,
    #[serde(default)]
    pub animate: bool,
}

impl ExternalView {
    pub fn initial(series: &Series) -> Self {
        let mut view = Self {
            frame: 0,
            association: Association::Point,
            field: None,
            component: Component::Value,
            displacement: None,
            scale: 1.0,
            show: true,
            animate: false,
        };
        if let Some(frame) = series.frames.first() {
            view.select_default_field(frame, None);
        }
        view
    }

    /// Pick the first drawable contour at the chosen location, or in point, then cell order.
    fn select_default_field(&mut self, frame: &Frame, location: Option<Association>) {
        self.field = None;
        for a in [Association::Point, Association::Cell].into_iter().filter(|&a| location.is_none_or(|chosen| a == chosen)) {
            for (name, f) in a.fields(frame) {
                if let Some(&c) = Component::for_components(f.components).first() {
                    self.association = a;
                    self.field = Some(name.clone());
                    self.component = c;
                    return;
                }
            }
        }
    }

    /// Keep the selected array and component where possible. A requested location limits
    /// fallback to arrays at that location when the new frame has different fields.
    pub(crate) fn inherit_contour(&mut self, frame: &Frame, location: Option<Association>) {
        if self.field.is_some() {
            let components = self.field(frame).map_or(&[][..], |f| Component::for_components(f.components));
            if let Some(&default) = components.first() {
                if !components.contains(&self.component) {
                    self.component = default;
                }
            } else {
                self.select_default_field(frame, location);
            }
        }
    }

    /// Drop an unavailable inherited warp. An explicit scale is checked after applying it.
    pub(crate) fn inherit_displacement(&mut self, frame: &Frame, scale_requested: bool) {
        if let Some(n) = &self.displacement
            && (!frame.point_data.get(n).is_some_and(|f| matches!(f.components, 2 | 3))
                || (!scale_requested && self.validate_displacement(frame).is_err()))
        {
            self.displacement = None;
        }
    }

    pub fn validate(&self, series: &Series) -> Result<(), String> {
        let frame = series.frames.get(self.frame).ok_or_else(|| format!("Frame {} is outside the imported series", self.frame))?;
        if !(0.0..=MAX_DEFORMATION_SCALE).contains(&self.scale) {
            return Err(format!("Deformation scale must be finite and between 0 and {MAX_DEFORMATION_SCALE}"));
        }
        if self.animate {
            if series.frames.len() < 2 {
                return Err("Playback needs at least two imported frames".into());
            }
            for f in series.frames.iter().skip(self.frame) {
                self.validate_frame(f)?;
            }
        } else {
            self.validate_frame(frame)?;
        }
        Ok(())
    }

    fn validate_frame(&self, frame: &Frame) -> Result<(), String> {
        if let Some(name) = &self.field {
            let f = self
                .association
                .fields(frame)
                .get(name)
                .ok_or_else(|| format!("No {} field named {name} at time {}", self.association.label(), frame.time))?;
            if !Component::for_components(f.components).contains(&self.component) {
                return Err(format!("{} is unavailable for field {name} with {} components", self.component.label(), f.components));
            }
        }
        self.validate_displacement(frame)
    }

    fn validate_displacement(&self, frame: &Frame) -> Result<(), String> {
        if let Some(n) = &self.displacement {
            let f = frame
                .point_data
                .get(n)
                .filter(|f| matches!(f.components, 2 | 3))
                .ok_or_else(|| format!("Displacement {n} must be a two- or three-component point field"))?;
            if !frame
                .grid
                .points
                .iter()
                .zip(f.values.chunks_exact(usize::from(f.components)))
                .all(|(&point, vector)| self.deformed_point(frame, point, vector).is_some())
            {
                return Err("Deformed coordinates exceed the model coordinate limit".into());
            }
        }
        Ok(())
    }

    fn deformed_point(&self, frame: &Frame, mut point: [f64; 3], vector: &[f64]) -> Option<[f64; 3]> {
        for (coordinate, &value) in point.iter_mut().zip(vector) {
            let meters = match &frame.displacement_unit {
                Some(unit) => fieldcraft_model::units::Quantity::Length.to_si(value, unit).ok()?,
                None => value,
            };
            *coordinate += self.scale * meters;
            if !coordinate.is_finite() || coordinate.abs() > fieldcraft_model::MAX_COORD {
                return None;
            }
        }
        Some(point)
    }

    /// The contour scale stays fixed as the user scrubs or plays a time series. A frame
    /// missing the selected array or component contributes no samples.
    pub fn bands(&self, series: &Series) -> Option<Bands> {
        let mut limits: Option<(f64, f64)> = None;
        for frame in &series.frames {
            if let Some(field) = self.field(frame).filter(|f| Component::for_components(f.components).contains(&self.component)) {
                for i in 0..field.values.len() / usize::from(field.components).max(1) {
                    if let Some(x) = self.component.sample(field, i) {
                        let (lo, hi) = limits.get_or_insert((x, x));
                        *lo = lo.min(x);
                        *hi = hi.max(x);
                    }
                }
            }
        }
        limits.map(|(lo, hi)| Bands::new(lo, hi, DEFAULT_BANDS))
    }

    /// Hold each frame until the following recorded time is reached during playback.
    /// Compare physical time differences from the selected frame so neighboring timestamps
    /// remain distinct even when they occupy a tiny fraction of the whole series.
    pub fn playback_frame(series: &Series, start: usize, fraction: f64) -> Option<usize> {
        let first = series.frames.first()?.time;
        let last = series.frames.last()?.time;
        let origin = series.frames.get(start)?.time;
        if last <= first || !fraction.is_finite() {
            return None;
        }
        let factor = if (last - first).is_finite() { 1.0 } else { 0.5 };
        if fraction <= 0.0 {
            return Some(start);
        }
        let span = last * factor - first * factor;
        let fraction = fraction.min(1.0);
        let elapsed = span * fraction;
        let later = series.frames.get(start + 1..)?;
        Some(
            start
                + later.partition_point(|f| {
                    let difference = f.time - origin;
                    let delta = if difference.is_finite() { difference * factor } else { f.time * factor - origin * factor };
                    // For spans at most one, division cannot underflow a nonzero delta.
                    // Compare fractions there so a subnormal elapsed time cannot round up
                    // to the next frame, even when the span itself is normal.
                    if span <= 1.0 { delta / span <= fraction } else { delta <= elapsed }
                }),
        )
    }
    fn field<'a>(&self, frame: &'a Frame) -> Option<&'a Field> {
        self.association.fields(frame).get(self.field.as_ref()?)
    }
}

fn element(vtk: u8) -> Option<ElementKind> {
    use ElementKind::*;
    Some(match vtk {
        3 => Rod2,
        5 => Tri3,
        9 => Quad4,
        10 => Tet4,
        12 => Hex8,
        13 => Wedge6,
        14 => Pyramid5,
        21 => Rod3,
        22 => Tri6,
        23 => Quad8,
        24 => Tet10,
        25 => Hex20,
        _ => return None,
    })
}

// Reserve at most half a frame budget for extracted topology and half for scene primitives.
// The source grid stays available if its surface representation will not fit.
const DISPLAY_BYTES: usize = fieldcraft_io::external::MAX_FRAME_BYTES / 2;

fn limited(xyz: &[[f64; 3]], palette: &Palette, reason: &str) -> (Scene, DisplayMesh, Shown) {
    let mut drawn = DisplayMesh::default();
    for &p in xyz {
        drawn.bounds.add(Vec3::from(p));
    }
    let note = format!("External results exceed the {} MiB display limit ({reason}).", DISPLAY_BYTES >> 20);
    let mut scene = Scene { background: Some(palette.background), ..Default::default() };
    scene.labels.push(fieldcraft_render::SceneLabel {
        p: drawn.bounds.center(),
        text: "External result exceeds display limit".into(),
        color: palette.marked,
    });
    (scene, drawn, Shown { failed: vec![note], ..Default::default() })
}

/// Reserve up to three fan triangles per intersected slab (at most five vertices).
/// When that bound does not fit, count the actual fan sizes before refusing the scene.
fn triangle_capacity(tris: &[DisplayTri], bands: Option<&Bands>, values: impl Fn(&DisplayTri) -> [f64; 3], limit: usize) -> Option<usize> {
    if tris.len() > limit {
        return None;
    }
    let Some(bands) = bands else { return Some(tris.len()) };
    let estimate = tris.iter().try_fold(0usize, |total, tri| {
        let bound = match bands.interior_bands(&values(tri)) {
            Some((lo, hi)) if lo < hi => ((hi - lo) as usize).saturating_add(1).saturating_mul(3),
            _ => 1,
        };
        let total = total.saturating_add(bound);
        (total <= limit).then_some(total)
    });
    estimate.or_else(|| {
        tris.iter().try_fold(0usize, |mut total, tri| {
            split_triangle(tri.p, values(tri), bands, &mut |piece, _| total = total.saturating_add(piece.len().saturating_sub(2)));
            (total <= limit).then_some(total)
        })
    })
}

/// Visit exactly the colored segments used for allocation and drawing. An interval can
/// span adjacent floats, so identify its band from its bounding levels.
fn split_line(line: &fieldcraft_render::DisplayLine, values: [Option<f64>; 2], bands: &Bands, mut emit: impl FnMut(Vec3, Vec3, Option<u32>)) {
    let [Some(first), Some(last)] = values else {
        emit(line.a, line.b, None);
        return;
    };
    let (a, b, lo, hi) = if first <= last { (line.a, line.b, first, last) } else { (line.b, line.a, last, first) };
    let Some((first_band, last_band)) = bands.interior_bands(&[lo, hi]) else {
        emit(a, b, None);
        return;
    };
    let mut start = a;
    for band in first_band..last_band {
        let level = bands.level(band + 1);
        let end = a + (b - a) * ((level - lo) / (hi - lo));
        emit(start, end, Some(band));
        start = end;
    }
    emit(start, b, Some(last_band));
}

/// Build the external picture in its source coordinates. Vertex cells remain visible in
/// their own right, line cells carry their endpoint or cell colors, and each volume's
/// faces share the model renderer's boundary extraction and higher-order tessellation.
pub fn scene(
    frame: &Frame,
    view: &ExternalView,
    style: Style,
    show: Show,
    palette: &Palette,
    bands: Option<Bands>,
    section: Option<(Section, Vec3)>,
) -> (Scene, DisplayMesh, Shown) {
    let vectors = view.displacement.as_ref().and_then(|n| frame.point_data.get(n));
    let mut mesh_budget = DISPLAY_BYTES;
    let mut moved = Vec::new();
    let xyz = if let Some(field) = vectors {
        let count = frame.grid.points.len();
        let Some(left) = mesh_budget.checked_sub(count.saturating_mul(size_of::<[f64; 3]>())) else {
            return limited(&frame.grid.points, palette, "deformation coordinates");
        };
        if moved.try_reserve_exact(count).is_err() {
            return limited(&frame.grid.points, palette, "deformation coordinates");
        }
        mesh_budget = left;
        for (&point, vector) in frame.grid.points.iter().zip(field.values.chunks_exact(usize::from(field.components).max(1))) {
            let Some(point) = view.deformed_point(frame, point, vector) else {
                return limited(&frame.grid.points, palette, "deformation coordinates");
            };
            moved.push(point);
        }
        moved.as_slice()
    } else {
        &frame.grid.points
    };
    // Cut by reference-coordinate cell centroid, exactly as for native mesh sections.
    // Keep the original numbering for both arrays, including when hidden cells are interleaved.
    let cut = section.map(|(cut, toward)| {
        let (at, sign) = cut.cut(frame.grid.points.iter().filter_map(|p| p.get(cut.axis).copied()), toward);
        (cut.axis, at, sign)
    });
    let visible = |cell: &Cell| {
        cut.is_none_or(|(axis, at, sign)| {
            let sum: f64 = cell.points.iter().filter_map(|&n| frame.grid.points.get(n as usize)?.get(axis)).sum();
            (sum / (cell.points.len().max(1) as f64) - at) * sign <= 0.0
        })
    };
    let cells = frame.grid.cells.iter().enumerate().filter(|(_, c)| visible(c));
    let midpoints = cells.clone().filter(|(_, c)| c.vtk_type == 21).count();
    let scratch = midpoints.saturating_mul(size_of::<fieldcraft_render::DisplayLine>()).saturating_add(if cut.is_some() { xyz.len() } else { 0 });
    let Some(mesh_budget) = mesh_budget.checked_sub(scratch) else {
        return limited(xyz, palette, "section or quadratic line elements");
    };
    let elements = cells.clone().filter_map(|(id, c)| Some((u32::try_from(id).ok()?, element(c.vtk_type)?, c.points.as_slice())));
    let mut drawn = match DisplayMesh::from_indexed(xyz, elements, mesh_budget) {
        Ok(mesh) => mesh,
        Err(reason) => return limited(xyz, palette, reason),
    };
    if cut.is_some() {
        let mut kept = vec![false; xyz.len()];
        for &n in cells.clone().flat_map(|(_, c)| &c.points) {
            if let Some(value) = kept.get_mut(n as usize) {
                *value = true;
            }
        }
        drawn.nodes.retain(|n| kept.get(n.label as usize) == Some(&true));
        drawn.bounds = Aabb3::EMPTY;
        for node in &drawn.nodes {
            drawn.bounds.add(node.p);
        }
    }
    // A quadratic line passes through its middle node in VTK (end, end, middle).
    if drawn.lines.try_reserve_exact(midpoints).is_err() {
        return limited(xyz, palette, "quadratic line elements");
    }
    for i in 0..drawn.lines.len() {
        if let Some(l) = drawn.lines.get_mut(i).filter(|l| l.kind == ElementKind::Rod3)
            && let Some(&mid) = frame.grid.cells.get(l.element as usize).and_then(|c| c.points.get(2))
            && let Some(&p) = xyz.get(mid as usize)
        {
            let mid_p = Vec3::from(p);
            let end = fieldcraft_render::DisplayLine { a: mid_p, nodes: [mid, l.nodes[1]], ..*l };
            l.b = mid_p;
            l.nodes[1] = mid;
            drawn.lines.push(end);
        }
    }

    let mut shown = Shown { scale: if vectors.is_some() { view.scale } else { 0.0 }, ..Default::default() };
    let field = view.field(frame);
    let sample = |id: usize| field.and_then(|f| view.component.sample(f, id));
    let total = match view.association {
        Association::Point => xyz.len(),
        Association::Cell => frame.grid.cells.len(),
    };
    let mut range: Option<((u32, f64), (u32, f64))> = None;
    if field.is_some() {
        for id in 0..total {
            if let (Ok(i), Some(v)) = (u32::try_from(id), sample(id)) {
                let (min, max) = range.get_or_insert(((i, v), (i, v)));
                if v < min.1 {
                    *min = (i, v);
                }
                if v > max.1 {
                    *max = (i, v);
                }
            }
        }
    }
    let bands = bands.or_else(|| range.map(|(min, max)| Bands::new(min.1, max.1, DEFAULT_BANDS)));
    if let (Some((min, max)), Some(name)) = (range, view.field.as_ref()) {
        shown.contour = Some((name.clone(), min, max));
    }
    // Interpolate in unit-range coordinates: two adjacent f64 samples have no
    // representable physical value between them, but their spatial band cuts do.
    let (factor, min, span) = bands.map_or((1.0, 0.0, 0.0), |b| {
        let factor = if (b.max - b.min).is_infinite() { 0.5 } else { 1.0 };
        let min = b.min * factor;
        (factor, min, b.max * factor - min)
    });
    let bands = bands.map(|b| Bands { min: 0.0, max: if span > 0.0 { 1.0 } else { 0.0 }, ..b });
    let value = |point: u32, cell: usize| {
        let v = sample(if view.association == Association::Point { point as usize } else { cell })?;
        Some(if span > 0.0 { (v * factor - min) / span } else { 0.0 })
    };
    let color =
        |v: Option<f64>, plain| bands.and_then(|b| Some(bands_color(&palette.contour, palette.contour_over, &b, b.band(v?)?))).unwrap_or(plain);
    let mut lines = if style == Style::Shaded { 0 } else { drawn.edges.len() };
    for line in &drawn.lines {
        if let Some(b) = &bands {
            split_line(line, line.nodes.map(|n| value(n, line.element as usize)), b, |_, _, _| lines = lines.saturating_add(1));
        } else {
            lines = lines.saturating_add(1);
        }
    }
    let mut points = cells.clone().filter(|(_, c)| c.vtk_type == 1).count();
    if show.nodes {
        points = points.saturating_add(drawn.nodes.len());
    }
    let fixed = lines.saturating_mul(size_of::<SceneLine>()).saturating_add(points.saturating_mul(size_of::<ScenePoint>()));
    let Some(remaining) = DISPLAY_BYTES.checked_sub(fixed) else {
        return limited(xyz, palette, "contour or edge primitives");
    };
    let filled = if style == Style::Wireframe { &[][..] } else { drawn.tris.as_slice() };
    let tri_values = |tri: &DisplayTri| tri.nodes.map(|id| value(id, tri.element as usize).unwrap_or(f64::NAN));
    let Some(triangles) = triangle_capacity(filled, bands.as_ref(), tri_values, remaining / size_of::<SceneTri>()) else {
        return limited(xyz, palette, "contour or edge primitives");
    };
    let mut scene = Scene { background: Some(palette.background), ..Default::default() };
    // Reserve exact capacities before drawing so geometric Vec growth cannot overshoot the bound.
    if scene.lines.try_reserve_exact(lines).is_err()
        || scene.tris.try_reserve_exact(triangles).is_err()
        || scene.points.try_reserve_exact(points).is_err()
    {
        return limited(xyz, palette, "scene allocation");
    }
    if style != Style::Shaded {
        scene.lines.extend(crate::view::mesh_outlines(&drawn, style, palette));
    }
    for tri in filled {
        if let Some(b) = bands {
            let band_color = |k: Option<u32>| k.map_or(palette.no_result, |k| bands_color(&palette.contour, palette.contour_over, &b, k));
            banded(tri.p, tri_values(tri), tri.normal, &b, &band_color, &mut scene.tris);
        } else {
            scene.tris.push(SceneTri { p: tri.p, normal: tri.normal, color: palette.element(tri.kind) });
        }
    }
    for line in &drawn.lines {
        if let Some(bands) = &bands {
            split_line(line, line.nodes.map(|n| value(n, line.element as usize)), bands, |a, b, band| {
                let color = band.map_or(palette.no_result, |k| bands_color(&palette.contour, palette.contour_over, bands, k));
                scene.lines.push(SceneLine { a, b, color, width: crate::view::LINE_ELEMENT_WIDTH });
            });
        } else {
            scene.lines.push(SceneLine { a: line.a, b: line.b, color: palette.element(line.kind), width: crate::view::LINE_ELEMENT_WIDTH });
        }
    }
    for (id, cell) in cells.filter(|(_, c)| c.vtk_type == 1) {
        if let Some(&n) = cell.points.first()
            && let Some(&p) = xyz.get(n as usize)
        {
            scene.points.push(ScenePoint { p: Vec3::from(p), color: color(value(n, id), palette.node), size: 5.0 });
        }
    }
    if show.nodes {
        scene.points.extend(drawn.nodes.iter().map(|n| ScenePoint { p: n.p, color: palette.node, size: 3.0 }));
    }
    (scene, drawn, shown)
}

/// The independently imported frame currently taking the place of the model in the viewport.
pub fn visible(s: &crate::Session) -> Option<(&Frame, &ExternalView)> {
    let view = s.external_view.as_ref().filter(|v| v.show)?;
    Some((s.external_results()?.frames.get(view.frame)?, view))
}

/// Return focus to the model or its solved steps, retaining the imported data and selection.
pub(crate) fn hide(s: &mut crate::Session) {
    if let Some(v) = &mut s.external_view {
        v.show = false;
        v.animate = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> Series {
        let displacements = Field { components: 2, values: vec![1.0, 0.0, 2.0, 0.0, 0.0, 1.0] };
        let grid = Grid { points: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]], cells: vec![Cell { vtk_type: 5, points: vec![0, 1, 2] }] };
        let frame = Frame {
            displacement_unit: None,
            time: 0.,
            grid,
            point_data: BTreeMap::from([("u".into(), displacements)]),
            cell_data: BTreeMap::from([("u".into(), Field { components: 1, values: vec![13.] })]),
        };
        Series { frames: vec![frame] }
    }

    #[test]
    fn point_vectors_move_only_the_view_and_cell_colors_stay_cell_values() {
        let series = fixture();
        let frame = &series.frames[0];
        let mut view = ExternalView::initial(&series);
        assert_eq!(view.component, Component::Magnitude);
        view.association = Association::Cell;
        view.component = Component::Value;
        view.displacement = Some("u".into());
        view.scale = 0.5;
        view.validate(&series).unwrap();
        let (s, b, shown) = scene(frame, &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, None, None);
        assert_eq!(b.bounds.max, Vec3::new(2.0, 1.5, 0.0));
        assert_eq!(s.tris.len(), 1);
        assert_eq!(s.tris[0].p, [Vec3::new(0.5, 0., 0.), Vec3::new(2., 0., 0.), Vec3::new(0., 1.5, 0.)]);
        assert_eq!(shown.contour, Some(("u".into(), (0, 13.), (0, 13.))));
        assert_eq!(frame.grid.points[1], [1., 0., 0.]);
        view.scale = 0.;
        let (_, ref_bounds, _) = scene(frame, &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, None, None);
        assert_eq!(ref_bounds.bounds.max, Vec3::new(1., 1., 0.));
        view.displacement = Some("missing".into());
        assert!(view.validate(&series).is_err());
    }
    #[test]
    fn a_section_reveals_the_shared_face_and_keeps_original_cell_indices() {
        let grid = Grid {
            points: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 0., -1.]],
            cells: vec![Cell { vtk_type: 10, points: vec![0, 1, 2, 3] }, Cell { vtk_type: 10, points: vec![0, 2, 1, 4] }],
        };
        let frame = Frame {
            displacement_unit: None,
            grid,
            time: 0.,
            point_data: Default::default(),
            cell_data: BTreeMap::from([("T".into(), Field { components: 1, values: vec![5., 30.] })]),
        };
        let series = Series { frames: vec![frame] };
        let view = ExternalView::initial(&series);
        let frame = &series.frames[0];
        let (_, all, _) = scene(frame, &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, None, None);
        assert_eq!(all.tris.len(), 6);
        let cut = Some((Section { axis: 2, at: Some(0.0) }, Vec3::Z));
        let (colored, sliced, _) = scene(frame, &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, None, cut);
        assert_eq!(sliced.tris.len(), 4);
        assert!(sliced.tris.iter().all(|t| t.element == 1));
        assert_eq!(sliced.bounds.max.z, 0.0);
        assert_eq!(sliced.bounds.min.z, -1.0);
        assert!(
            colored.tris.iter().all(|t| t.color == fieldcraft_render::contour::band_color(&Palette::LIGHT.contour, DEFAULT_BANDS - 1, DEFAULT_BANDS))
        );
    }

    #[test]
    fn playback_follows_irregular_times_and_contours_use_one_scale() {
        let mut series = fixture();
        for (time, pressure) in [(1.0, 17.0), (4.0, 43.0)] {
            let mut next = series.frames[0].clone();
            next.time = time;
            next.cell_data.get_mut("u").unwrap().values[0] = pressure;
            series.frames.push(next);
        }
        let mut view = ExternalView::initial(&series);
        view.association = Association::Cell;
        view.component = Component::Value;
        let bands = view.bands(&series).unwrap();
        assert_eq!((bands.min, bands.max), (13., 43.));
        for (frame, fraction) in [(0, 0.), (1, 0.25), (1, 0.5), (1, 0.99), (2, 1.)] {
            assert_eq!(ExternalView::playback_frame(&series, 0, fraction), Some(frame));
        }
        assert_eq!(ExternalView::playback_frame(&series, 1, 0.75), Some(2));
        let (start, _, _) = scene(&series.frames[0], &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, Some(bands), None);
        let (end, _, _) = scene(&series.frames[2], &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, Some(bands), None);
        assert_ne!(start.tris[0].color, end.tris[0].color);
        view.animate = true;
        series.frames[2].cell_data.clear();
        assert!(view.validate(&series).is_err());
    }

    #[test]
    fn playback_measures_neighboring_samples_from_the_selected_time() {
        let make = |times: [f64; 3]| {
            let sample = fixture().frames.remove(0);
            Series { frames: times.map(|time| Frame { time, ..sample.clone() }).into() }
        };
        let close = make([-1e16, 0., 1.]);
        assert_eq!(ExternalView::playback_frame(&close, 1, 0.), Some(1));
        assert_eq!(ExternalView::playback_frame(&close, 1, 0.5e-16), Some(1));
        assert_eq!(ExternalView::playback_frame(&close, 1, 1e-16), Some(2));
        let wide = make([-f64::MAX, 0., f64::MAX]);
        assert_eq!(ExternalView::playback_frame(&wide, 0, 0.25), Some(0));
        assert_eq!(ExternalView::playback_frame(&wide, 0, 0.5), Some(1));
        assert_eq!(ExternalView::playback_frame(&wide, 0, 1.), Some(2));
        let m = f64::from_bits(1);
        let tiny = make([0., m, 2. * m]);
        for (fraction, expected) in [(0.49, 0), (0.5, 1), (0.75, 1), (1., 2)] {
            assert_eq!(ExternalView::playback_frame(&tiny, 0, fraction), Some(expected));
        }
        assert_eq!(ExternalView::playback_frame(&tiny, 1, 0.4), Some(1));
        assert_eq!(ExternalView::playback_frame(&tiny, 1, 0.5), Some(2));
        for end in [f64::MIN_POSITIVE, f64::MIN_POSITIVE.next_up()] {
            let normal = make([0., m, end]);
            let first_step = m / end;
            assert_eq!(ExternalView::playback_frame(&normal, 0, first_step.next_down()), Some(0));
            assert_eq!(ExternalView::playback_frame(&normal, 0, first_step), Some(1));
            assert_eq!(ExternalView::playback_frame(&normal, 0, 1.0_f64.next_down()), Some(1));
            assert_eq!(ExternalView::playback_frame(&normal, 0, 1.), Some(2));
        }
    }

    #[test]
    fn playback_validates_only_the_frames_it_will_visit() {
        use serde_json::json;
        let sample = fixture().frames.remove(0);
        let mut frames = [0., 1., 2.].map(|time| Frame { time, ..sample.clone() });
        frames[0].point_data.clear();
        let mut session = crate::Session::default();
        session.set_external_results(Series { frames: frames.into() }).unwrap();
        session.execute("results.external.view", &json!({"frame": 1, "association": "point", "field": "u"})).unwrap();
        let play = session.execute("results.animate", &json!({})).unwrap();
        assert_eq!(play["view"]["frame"], 1);
        assert_eq!(play["view"]["animate"], true);
        session.execute("results.external.view", &json!({"frame": 2})).unwrap();
        // Playing from the last frame restarts at frame zero, which lacks the selected point field.
        let error = session.execute("results.animate", &json!({})).unwrap_err();
        assert!(error.to_string().contains("No Point field named u at time 0"));
        session.execute("results.external.view", &json!({"field": ""})).unwrap();
        let restart = session.execute("results.animate", &json!({})).unwrap();
        assert_eq!(restart["view"]["frame"], 0);
        assert_eq!(restart["view"]["animate"], true);
    }

    #[test]
    fn line_contours_follow_every_band_at_extreme_field_values() {
        let grid = Grid { points: vec![[0., 0., 0.], [1., 0., 0.]], cells: vec![Cell { vtk_type: 3, points: vec![0, 1] }] };
        let frame = Frame {
            displacement_unit: None,
            grid,
            time: 0.,
            point_data: BTreeMap::from([("T".into(), Field { components: 1, values: vec![-f64::MAX, f64::MAX] })]),
            cell_data: BTreeMap::new(),
        };
        let series = Series { frames: vec![frame] };
        let view = ExternalView::initial(&series);
        let (scene, _, shown) = scene(&series.frames[0], &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, view.bands(&series), None);
        assert!(shown.failed.is_empty());
        assert_eq!(scene.lines.len(), DEFAULT_BANDS as usize);
        for (k, l) in scene.lines.iter().enumerate() {
            assert!((l.a.x - k as f64 / f64::from(DEFAULT_BANDS)).abs() < 1e-12);
            assert_eq!(l.color, fieldcraft_render::contour::band_color(&Palette::LIGHT.contour, k as u32, DEFAULT_BANDS));
        }
    }

    #[test]
    fn line_intervals_use_the_cut_levels_at_rounded_boundaries() {
        let bands = Bands::new(0., 1., DEFAULT_BANDS);
        let line = fieldcraft_render::DisplayLine { a: Vec3::ZERO, b: Vec3::X, nodes: [0, 1], kind: ElementKind::Rod2, element: 0 };
        let intervals = |a, b| {
            let mut pieces = Vec::new();
            split_line(&line, [Some(a), Some(b)], &bands, |start, end, band| pieces.push((start.x, end.x, band.unwrap())));
            pieces
        };
        // No float lies inside this interval; its rounded midpoint lands in band 6.
        assert_eq!(intervals(0.5_f64.next_down(), 0.5), [(0., 1., 5)]);
        assert_eq!(intervals(0.5, 0.5_f64.next_down()), [(1., 0., 5)]);
        // The endpoint's band calculation rounds up, but level 5 cuts the line's interior.
        let lo = bands.level(5).next_down();
        assert_eq!(lo * f64::from(DEFAULT_BANDS), 5.);
        assert_eq!(bands.band(lo), Some(4));
        assert_eq!(intervals(lo, bands.level(5)), [(0., 1., 4)]);
        assert_eq!(bands.band(0.44), Some(5));
        let two = intervals(lo, 0.44);
        assert_eq!(two.iter().map(|p| p.2).collect::<Vec<_>>(), [4, 5]);
        assert!(two.iter().all(|&(start, end, _)| start < end));
        let pieces = intervals(lo, 0.55);
        assert_eq!(pieces.iter().map(|p| p.2).collect::<Vec<_>>(), [4, 5, 6]);
        assert_eq!(pieces[0].0, 0.);
        assert_eq!(pieces[2].1, 1.);
        assert!(pieces.iter().all(|&(start, end, _)| start < end));
    }

    #[test]
    fn triangle_budget_accepts_the_actual_fan_count_at_the_limit() {
        let bands = Bands::new(0., 1., DEFAULT_BANDS);
        let tri = DisplayTri { p: [Vec3::ZERO, Vec3::X, Vec3::Y], nodes: [0, 1, 2], normal: Vec3::Z, kind: ElementKind::Tri3, element: 0, face: 1 };
        let values = [0., 0., 1.];
        let mut drawn = Vec::new();
        banded(tri.p, values, tri.normal, &bands, &|_| Palette::LIGHT.no_result, &mut drawn);
        // Eleven quadrilateral strips and one triangular tip per source triangle.
        assert_eq!(drawn.len(), 23);
        let tris = [tri; 2];
        let exact = tris.len() * drawn.len();
        assert_eq!(triangle_capacity(&tris, Some(&bands), |_| values, exact), Some(exact));
        assert_eq!(triangle_capacity(&tris, Some(&bands), |_| values, exact - 1), None);
    }

    #[test]
    fn one_ulp_contours_cut_lines_and_triangles_at_each_band() {
        let mut series = fixture();
        let upper = 1.0_f64.next_up();
        let frame = &mut series.frames[0];
        frame.grid.cells.push(Cell { vtk_type: 3, points: vec![0, 1] });
        frame.point_data.insert("u".into(), Field { components: 1, values: vec![1.0, upper, 1.0] });
        frame.cell_data.clear();
        let view = ExternalView::initial(&series);
        let (scene, _, shown) = scene(&series.frames[0], &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, view.bands(&series), None);
        assert_eq!(shown.contour, Some(("u".into(), (0, 1.0), (1, upper))));
        assert_eq!(scene.lines.len(), DEFAULT_BANDS as usize);
        for (k, line) in scene.lines.iter().enumerate() {
            let left = k as f64 / f64::from(DEFAULT_BANDS);
            let right = (k + 1) as f64 / f64::from(DEFAULT_BANDS);
            let color = fieldcraft_render::contour::band_color(&Palette::LIGHT.contour, k as u32, DEFAULT_BANDS);
            assert_eq!((line.a.x, line.b.x, line.color), (left, right, color));
            // The analytic triangle spans 0 ≤ y ≤ 1 − x, so a strip's area
            // distinguishes a real band from coincident, zero-area cuts.
            let area: f64 = scene
                .tris
                .iter()
                .filter(|t| t.color == color)
                .map(|t| {
                    let [a, b, c] = t.p;
                    (b - a).cross(c - a).len() * 0.5
                })
                .sum();
            assert!((area - (right - left) * (1.0 - (left + right) * 0.5)).abs() < 1e-14);
        }
    }

    #[test]
    fn quadratic_line_draws_both_segments_through_the_source_midpoint() {
        let points = vec![[0., 0., 0.], [2., 0., 0.], [1., 1., 0.]];
        let grid = Grid { points, cells: vec![Cell { vtk_type: 21, points: vec![0, 1, 2] }] };
        let frame = Frame { displacement_unit: None, grid, time: 0., point_data: Default::default(), cell_data: Default::default() };
        let series = Series { frames: vec![frame] };
        let view = ExternalView::initial(&series);
        let (scene, mesh, shown) = scene(&series.frames[0], &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, None, None);
        assert!(shown.failed.is_empty());
        assert_eq!(mesh.lines.len(), 2);
        assert_eq!(scene.lines.len(), 2);
        assert_eq!(mesh.lines[0].nodes, [0, 2]);
        assert_eq!(mesh.lines[1].nodes, [2, 1]);
        assert_eq!(scene.lines[0].b, Vec3::new(1., 1., 0.));
        assert_eq!(scene.lines[0].b, scene.lines[1].a);
    }

    #[test]
    fn invalid_scales_and_excessive_displacements_are_refused() {
        let mut series = fixture();
        let mut view = ExternalView::initial(&series);
        view.displacement = Some("u".into());
        for invalid in [-1., f64::NAN, f64::INFINITY, MAX_DEFORMATION_SCALE * 2.0] {
            view.scale = invalid;
            assert!(view.validate(&series).is_err());
        }
        view.scale = 2.;
        series.frames[0].point_data.get_mut("u").unwrap().values[0] = f64::MAX;
        assert!(view.validate(&series).is_err());
    }

    #[test]
    fn vector_components_color_point_contours_and_distinguish_cell_values() {
        let series = fixture();
        let mut view = ExternalView::initial(&series);
        view.component = Component::X;
        let (scene, _, shown) = scene(&series.frames[0], &view, Style::Shaded, Show::DEFAULT, &Palette::LIGHT, None, None);
        assert_eq!(shown.contour, Some(("u".into(), (2, 0.), (1, 2.))));
        let colors: std::collections::BTreeSet<_> = scene.tris.iter().map(|t| (t.color.0, t.color.1, t.color.2)).collect();
        assert_eq!(colors.len(), DEFAULT_BANDS as usize);
        view.component = Component::Z;
        assert!(view.validate(&series).is_err());
    }
}
