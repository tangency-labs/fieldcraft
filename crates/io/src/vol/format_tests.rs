//! Read the writer's documented records and check them against analytic box geometry.

use super::*;
use fieldcraft_model::Set;

type Point = [f64; 3];

fn sub([a, b, c]: Point, [x, y, z]: Point) -> Point {
    [a - x, b - y, c - z]
}

fn cross([a, b, c]: Point, [x, y, z]: Point) -> Point {
    [b * z - c * y, c * x - a * z, a * y - b * x]
}

fn dot([a, b, c]: Point, [x, y, z]: Point) -> f64 {
    a * x + b * y + c * z
}

fn near(value: f64, expected: f64) {
    assert!((value - expected).abs() < 1e-12, "expected {expected}, found {value}");
}

fn records<'a>(text: &'a str, section: &str) -> impl Iterator<Item = &'a str> {
    let mut lines = text.lines();
    assert_eq!(lines.find(|&line| line == section), Some(section));
    let count = lines.next().unwrap().parse().unwrap();
    lines.take(count)
}

fn numbers(line: &str) -> Vec<usize> {
    line.split_whitespace().map(|number| number.parse().unwrap()).collect()
}

/// Three adjacent unit cubes, each a distinct material region. A single named boundary spans
/// the two internal planes. The tetrahedral variant uses the body diagonal in each cube.
fn three_cubes(kind: ElementKind) -> (Model, [String; 3]) {
    let mut model = Model::default();
    for x in 0..=3 {
        for (i, [y, z]) in [[0., 0.], [1., 0.], [0., 1.], [1., 1.]].into_iter().enumerate() {
            model.set_node(4 * x + i as u32 + 1, [f64::from(x), y, z]).unwrap();
        }
    }
    let materials = ["left", "middle", "right"].map(str::to_owned);
    let mut interfaces = Vec::new();
    for (i, material) in materials.iter().enumerate() {
        let labels = [1, 5, 6, 2, 3, 7, 8, 4].map(|node| node + 4 * i as u32);
        let cells = if kind == ElementKind::Hex8 {
            vec![labels.to_vec()]
        } else {
            [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]].map(|tet| tet.map(|j| labels[j]).to_vec()).to_vec()
        };
        let mut elements = Vec::new();
        for nodes in cells {
            let label = model.elements().len() as u32 + 1;
            model.set_element(label, kind, nodes).unwrap();
            elements.push(label);
            for face in 1..=kind.faces().len() as u8 {
                let facet = ElementFace { element: label, face };
                let nodes = model.face_nodes(facet).unwrap();
                if [1., 2.].into_iter().any(|x| nodes.iter().all(|&n| model.node(n).unwrap().xyz[0] == x)) {
                    interfaces.push(facet);
                }
            }
        }
        model.add_set(material, Set { elements, ..Set::default() }).unwrap();
    }
    model.add_set("joint", Set { faces: interfaces, ..Set::default() }).unwrap();
    (model, materials)
}

#[test]
fn vol_records_close_each_domain_and_match_box_measures() {
    for kind in [ElementKind::Tet4, ElementKind::Hex8] {
        let (model, materials) = three_cubes(kind);
        let mut output = Vec::new();
        serialize(&mut output, &model, &materials, &["joint".into()]).unwrap();
        let text = String::from_utf8(output).unwrap();
        let points: Vec<Point> = records(&text, "points")
            .map(|line| {
                let mut xyz = line.split_whitespace().map(|coordinate| coordinate.parse().unwrap());
                [xyz.next().unwrap(), xyz.next().unwrap(), xyz.next().unwrap()]
            })
            .collect();
        let point = |n: usize| points[n - 1];
        let mut domain_volumes = [0.; 3];
        for line in records(&text, "volumeelements") {
            let values = numbers(line);
            let [domain, count, nodes @ ..] = values.as_slice() else {
                panic!("bad volume record");
            };
            assert_eq!(*count, nodes.len());
            let a = point(nodes[0]);
            let u = sub(point(nodes[1]), a);
            let (v, w, scale) = match *count {
                4 => (sub(point(nodes[2]), a), sub(point(nodes[3]), a), 1. / 6.),
                8 => (sub(point(nodes[3]), a), sub(point(nodes[4]), a), 1.),
                _ => panic!("unexpected volume element"),
            };
            // The first face points out. Changing the order of its tangents gives the
            // positive volume basis with the last direction pointing into the cell.
            let volume = dot(cross(v, u), w) * scale;
            assert!(volume > 0., "{kind:?} has signed volume {volume}");
            domain_volumes[domain - 1] += volume;
        }
        for volume in domain_volumes {
            near(volume, 1.);
        }
        near(domain_volumes.iter().sum(), 3.);

        let names: BTreeMap<usize, &str> = records(&text, "bcnames")
            .map(|line| {
                let (index, name) = line.split_once(' ').unwrap();
                (index.parse().unwrap(), name)
            })
            .collect();
        let mut domain_areas = [0.; 3];
        let mut exterior_area = 0.;
        let mut interface_area = 0.;
        let mut domain_edges = BTreeMap::<_, i32>::new();
        let mut joints = BTreeMap::new();
        let center = |domain| [domain as f64 - 0.5, 0.5, 0.5];
        for line in records(&text, "surfaceelements") {
            let values = numbers(line);
            let [descriptor, boundary, inside, outside, count, nodes @ ..] = values.as_slice() else {
                panic!("bad surface record");
            };
            assert_eq!(*count, nodes.len());
            let origin = point(nodes[0]);
            let normal = cross(sub(point(nodes[1]), origin), sub(point(nodes[2]), origin));
            let area = dot(normal, normal).sqrt() * if *count == 3 { 0.5 } else { 1. };
            assert!(dot(normal, sub(center(*inside), origin)) < 0.);
            if *outside == 0 {
                assert_eq!(names[boundary], "default");
                exterior_area += area;
            } else {
                assert_eq!(names[boundary], "joint");
                assert!(dot(normal, sub(center(*outside), origin)) > 0.);
                let expected = joints.entry((*inside, *outside)).or_insert((*descriptor, *boundary));
                assert_eq!(*expected, (*descriptor, *boundary));
                interface_area += area;
            }
            for (domain, sign) in [(*inside, 1), (*outside, -1)].into_iter().filter(|(domain, _)| *domain != 0) {
                domain_areas[domain - 1] += area;
                for (&a, &b) in nodes.iter().zip(nodes.iter().cycle().skip(1)) {
                    *domain_edges.entry((domain, a.min(b), a.max(b))).or_default() += if a < b { sign } else { -sign };
                }
            }
        }
        assert!(domain_edges.values().all(|&orientation_sum| orientation_sum == 0));
        assert_eq!(joints.len(), 2);
        assert_ne!(joints[&(1, 2)].0, joints[&(2, 3)].0);
        assert_eq!(joints[&(1, 2)].1, joints[&(2, 3)].1);
        for area in domain_areas {
            near(area, 6.);
        }
        near(exterior_area, 14.);
        near(interface_area, 2.);
    }
}

/// Read the two planar cells and their named segments back from the VOL sections. The
/// triangles split a unit square on its diagonal. The inside domain is left of each segment.
pub(super) fn check_split_square(text: &str) {
    let points: Vec<Point> = records(text, "points")
        .map(|line| {
            let xyz: Vec<f64> = line.split_whitespace().map(|v| v.parse().unwrap()).collect();
            xyz.try_into().unwrap()
        })
        .collect();
    let point = |n: usize| points[n - 1];
    let mut centers = BTreeMap::new();
    for line in records(text, "surfaceelements") {
        let row = numbers(line);
        let &[domain, matching, 0, 0, 3, i, j, k] = row.as_slice() else { panic!("invalid planar triangle record") };
        assert_eq!(domain, matching);
        let [a, b, c] = [i, j, k].map(point);
        near(0.5 * cross(sub(b, a), sub(c, a))[2], 0.5);
        centers.insert(domain, std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.));
    }
    assert_eq!(centers.len(), 2);
    assert_eq!(records(text, "volumeelements").count(), 0);
    let labels: BTreeMap<usize, &str> = records(text, "bcnames")
        .map(|line| {
            let (id, name) = line.split_once(' ').unwrap();
            (id.parse().unwrap(), name)
        })
        .collect();
    let mut exterior_lengths = BTreeMap::<usize, f64>::new();
    let mut descriptor_pairs = BTreeMap::new();
    let mut joint_length = 0.;
    for line in records(text, "edgesegmentsgi2") {
        let fields: Vec<i32> = line.split_whitespace().map(|f| f.parse().unwrap()).collect();
        let &[id, 0, i, j, -1, -1, inside, outside, first, 0, second, 1] = fields.as_slice() else { panic!("invalid segment record") };
        assert_eq!((first, second), (id, id));
        let id = usize::try_from(id).unwrap();
        let [inside, outside] = [inside, outside].map(|v| usize::try_from(v).unwrap());
        assert_eq!(*descriptor_pairs.entry(id).or_insert((inside, outside)), (inside, outside));
        let [a, b] = [i, j].map(|v| point(usize::try_from(v).unwrap()));
        let d = sub(b, a);
        assert!(cross(d, sub(centers[&inside], a))[2] > 0.);
        let length = dot(d, d).sqrt();
        if outside == 0 {
            assert_eq!(labels[&id], "default");
            *exterior_lengths.entry(inside).or_default() += length;
        } else {
            assert!(cross(d, sub(centers[&outside], a))[2] < 0.);
            assert_eq!(labels[&id], "electrode");
            joint_length += length;
        }
    }
    assert_eq!(descriptor_pairs.len(), 3);
    assert_eq!(exterior_lengths.len(), 2);
    for length in exterior_lengths.values() {
        near(*length, 2.);
    }
    near(joint_length, 2_f64.sqrt());
}
