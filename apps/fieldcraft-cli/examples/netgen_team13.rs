//! Regenerate the FieldCraft export recipe from Johannes Maierhofer's SolidCraft TEAM 13
//! design. Figure 1 and Table 1 of the COMPUMAG TEAM 13 problem define the steel and winding
//! dimensions. The air enclosure, block partitions, and face selections are workshop choices.
//! https://www.compumag.org/wp/wp-content/uploads/2018/06/problem13.pdf
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

const MATERIALS: [&str; 10] = [
    "air",
    "iron",
    "brick_back",
    "brick_front",
    "brick_left",
    "brick_right",
    "corner_left_back",
    "corner_left_front",
    "corner_right_back",
    "corner_right_front",
];
const BOUNDARIES: [&str; 4] = ["outer", "gamma_core", "gamma_coil", "winding_seams"];
const AIR_HALF_WIDTH: f64 = 250.;
const COIL_HALF_STRAIGHT: f64 = 50.;
const COIL_INNER_RADIUS: f64 = 25.;
const COIL_OUTER_RADIUS: f64 = 50.;
const COIL_HALF_HEIGHT: f64 = 50.;
const STEEL_INNER_Z: f64 = 60.;
const STEEL_OUTER_Z: f64 = 63.2;

// Symmetric planes through the faces of the center plate, channels, and winding.
// The additional 62.5, 120 and 25 mm planes divide the mesh around and through the winding.
fn intervals(positive: &[f64]) -> Vec<[f64; 2]> {
    let planes: Vec<_> = positive.iter().rev().map(|x| -x).chain([0.]).chain(positive.iter().copied()).collect();
    planes.iter().zip(planes.iter().skip(1)).map(|(&a, &b)| [a, b]).collect()
}

fn body([x, y, z]: [f64; 3]) -> &'static str {
    let (u, v, w) = (x.abs(), y.abs(), z.abs());
    if u < 1.6 && v < 25. && w < STEEL_OUTER_Z {
        return "iron_center";
    }
    if x * y > 0. && u > 2.1 && u < 125.3 && v > 15. && v < 65. && w < STEEL_OUTER_Z && (u > 122.1 || w > STEEL_INNER_Z) {
        return if x < 0. { "iron_left" } else { "iron_right" };
    }
    let radius = u.max(v) - COIL_HALF_STRAIGHT;
    if w >= COIL_HALF_HEIGHT || radius <= COIL_INNER_RADIUS || radius >= COIL_OUTER_RADIUS {
        return "air";
    }
    match (u > COIL_HALF_STRAIGHT, v > COIL_HALF_STRAIGHT, x > 0., y > 0.) {
        (true, true, true, true) => "corner_right_back",
        (true, true, true, false) => "corner_right_front",
        (true, true, false, true) => "corner_left_back",
        (true, true, false, false) => "corner_left_front",
        (true, false, true, _) => "brick_right",
        (true, false, false, _) => "brick_left",
        (_, _, _, true) => "brick_back",
        (_, _, _, false) => "brick_front",
    }
}

// Map square rings onto circular winding corners. Fade the map to the Cartesian grid
// before the back of the channels and before the horizontal steel flanges.
fn mapped([x, y, z]: [f64; 3]) -> [f64; 3] {
    let cx = x.clamp(-COIL_HALF_STRAIGHT, COIL_HALF_STRAIGHT);
    let cy = y.clamp(-COIL_HALF_STRAIGHT, COIL_HALF_STRAIGHT);
    let (dx, dy) = (x - cx, y - cy);
    let radius = dx.abs().max(dy.abs());
    if radius == 0. {
        return [x, y, z];
    }
    let radial_weight = ((70. - radius) / (70. - COIL_OUTER_RADIUS)).clamp(0., 1.);
    let axial_weight = ((STEEL_INNER_Z - z.abs()) / (STEEL_INNER_Z - COIL_HALF_HEIGHT)).clamp(0., 1.);
    let factor = 1. + axial_weight * radial_weight * (radius / (dx * dx + dy * dy).sqrt() - 1.);
    [cx + dx * factor, cy + dy * factor, z]
}

fn material(owner: &'static str) -> &'static str {
    if owner.starts_with("iron_") { "iron" } else { owner }
}

// Face IDs are persistent names from the accompanying native design. A cutout in each
// channel supplies its three inner faces and splits the ends of its two flanges.
fn cad_face(owner: &str, face: &str, [x, y, z]: [f64; 3]) -> String {
    let feature = match owner {
        "corner_right_back" => 1,
        "corner_left_back" => 4,
        "corner_left_front" => 7,
        "corner_right_front" => 10,
        "brick_back" => 13,
        "brick_front" => 14,
        "brick_left" => 15,
        "brick_right" => 16,
        "iron_center" => 17,
        "iron_left" => 18,
        "iron_right" => 20,
        _ => 22,
    };
    if owner.starts_with("corner_") {
        if face == "-z" || face == "+z" {
            return format!("F{feature}:{}", if face == "-z" { "bottom" } else { "top" });
        }
        let inner = x.abs().max(y.abs()) == COIL_HALF_STRAIGHT + COIL_INNER_RADIUS;
        if owner == "corner_right_back" {
            return format!("merge:right_back_winding_{}_wall", if inner { "inner" } else { "outer" });
        }
        return format!("F{}:side", feature + u32::from(inner));
    }
    if owner == "iron_left" || owner == "iron_right" {
        let cut = feature + 1;
        let zface = if z > 0. { "+z" } else { "-z" };
        if x.abs() == 2.1 {
            return format!("F{feature}:{face}#[F{feature}:{zface},F{cut}:{zface}]");
        }
        if x.abs() == 122.1 {
            return format!("F{cut}:{}", if x > 0. { "+x" } else { "-x" });
        }
        if z.abs() == STEEL_INNER_Z {
            return format!("F{cut}:{zface}");
        }
    }
    format!("F{feature}:{face}")
}

fn main() -> io::Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    writeln!(out, "[")?;
    let mut command = |id: &str, params: Value| -> io::Result<()> {
        let record = json!({"command": id, "params": params});
        writeln!(out, "  {record},")
    };
    command("model.units", json!({"system":"mm-t-s-N-MPa"}))?;
    let mut materials = BTreeMap::<&str, BTreeSet<&str>>::new();
    let mut boundaries = BTreeMap::<&str, BTreeSet<(&str, String)>>::new();
    let x_intervals = intervals(&[1.6, 2.1, 50., 62.5, 75., 100., 120., 122.1, 125.3, AIR_HALF_WIDTH]);
    let y_intervals = intervals(&[15., 25., 50., 62.5, 65., 75., 100., 120., AIR_HALF_WIDTH]);
    let z_intervals = intervals(&[25., COIL_HALF_HEIGHT, STEEL_INNER_Z, STEEL_OUTER_Z, AIR_HALF_WIDTH]);
    let sides = [(2, false, "-z"), (2, true, "+z"), (1, false, "-y"), (0, true, "+x"), (1, true, "+y"), (0, false, "-x")];
    for &[x0, x1] in &x_intervals {
        for &[y0, y1] in &y_intervals {
            for &[z0, z1] in &z_intervals {
                let limits = [[x0, x1], [y0, y1], [z0, z1]];
                let center = limits.map(|[a, b]| a.midpoint(b));
                let owner = body(center);
                let region = material(owner);
                materials.entry(region).or_default().insert(owner);
                let mut faces = BTreeMap::new();
                for (index, &(axis, positive, direction)) in sides.iter().enumerate() {
                    let mut point = center;
                    point[axis] = limits[axis][usize::from(positive)];
                    let boundary = if point[axis].abs() == AIR_HALF_WIDTH {
                        "outer"
                    } else {
                        // Sample across the face without a tolerance smaller than the steel gap.
                        let mut across = point;
                        across[axis] = if positive { across[axis].next_up() } else { across[axis].next_down() };
                        let neighbor = material(body(across));
                        if owner == "air" || neighbor == region {
                            continue;
                        }
                        if neighbor == "air" {
                            if region == "iron" { "gamma_core" } else { "gamma_coil" }
                        } else if region < neighbor {
                            "winding_seams"
                        } else {
                            continue;
                        }
                    };
                    let name = cad_face(owner, direction, point);
                    boundaries.entry(boundary).or_default().insert((owner, name.clone()));
                    faces.insert(index + 1, name);
                }
                let corners =
                    [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0], [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]].map(mapped);
                command(
                    "mesh.block",
                    json!({"body":owner,"corners":corners,"divisions":[1,1,1],"edges":[],
                    "faces":faces,"kind":"tet4","merge":1e-7}),
                )?;
            }
        }
    }
    for (name, owners) in materials {
        let geometry: Vec<_> = owners.iter().map(|owner| json!({"body":owner})).collect();
        command("model.set.create", json!({"name":name,"geometry":geometry}))?;
    }
    for (name, faces) in boundaries {
        let geometry: Vec<_> = faces.iter().map(|(owner, face)| json!({"body":owner,"face":face})).collect();
        command("model.set.create", json!({"name":name,"geometry":geometry}))?;
    }
    let export = json!({"command":"mesh.export.vol","params":{"path":"team13.vol","materials":MATERIALS,"boundaries":BOUNDARIES}});
    writeln!(out, "  {export}\n]")?;
    out.flush()
}
