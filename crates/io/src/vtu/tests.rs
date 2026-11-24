use std::io::Write;

use base64::Engine;
use base64::engine::general_purpose::STANDARD;
use flate2::{Compression, write::ZlibEncoder};

use super::{collection, read, single_piece};

const POINTS: &str = r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 0 1 0 0 0 0</DataArray>"#;
const POINT_FIELD: &str = r#"<DataArray type="Float32" Name="E" NumberOfComponents="2" format="ascii">0 1 2 3 4 5 6 7</DataArray>"#;

fn grid(points: &str, field: &str, appended: &[u8], options: &str) -> Vec<u8> {
    let mut file = format!(
        r#"<VTKFile type="UnstructuredGrid" {options}><UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="2"><Points>{points}</Points><Cells>
          <DataArray type="Int32" Name="connectivity">0 1 2 3 1 2</DataArray>
          <DataArray type="Int64" Name="offsets">3 6</DataArray><DataArray type="UInt8" Name="types">5 5</DataArray>
          </Cells><PointData>{field}</PointData><CellData><DataArray type="Float64" Name="S">7 -2</DataArray></CellData></Piece></UnstructuredGrid>"#
    )
    .into_bytes();
    file.extend_from_slice(appended);
    file.extend_from_slice(b"</VTKFile>");
    file
}

#[test]
fn point_ids_and_both_field_associations_survive() {
    let frame = read(&grid(POINTS, POINT_FIELD, b"", ""), 0.25).unwrap();
    assert_eq!(frame.time, 0.25);
    assert_eq!(frame.grid.points.len(), 4);
    assert_eq!(frame.grid.points[0], frame.grid.points[3]);
    assert_eq!(frame.grid.cells[1].points, [3, 1, 2]);
    assert_eq!(frame.point_data["E"].components, 2);
    assert_eq!(frame.point_data["E"].values, [0., 1., 2., 3., 4., 5., 6., 7.]);
    assert_eq!(frame.cell_data["S"].values, [7., -2.]);
}

#[test]
fn appended_raw_32_and_64_bit_headers_and_both_byte_orders() {
    for (wide, little) in [(false, true), (true, false)] {
        let options = if little { r#"byte_order="LittleEndian""# } else { r#"byte_order="BigEndian" header_type="UInt64""# };
        let mut block = b"<AppendedData encoding=\"raw\">\n_".to_vec();
        if wide {
            block.extend_from_slice(&32_u64.to_be_bytes());
        } else {
            block.extend_from_slice(&32_u32.to_le_bytes());
        }
        for x in [0_f32, 1., 2., 3., 4., 5., 6., -7.] {
            block.extend_from_slice(&if little { x.to_le_bytes() } else { x.to_be_bytes() });
        }
        block.extend_from_slice(b"\n</AppendedData>");
        let field = r#"<DataArray type="Float32" Name="E" NumberOfComponents="2" format="appended" offset="0"/>"#;
        let frame = read(&grid(POINTS, field, &block, options), 0.).unwrap();
        assert_eq!(frame.point_data["E"].values.last(), Some(&-7.));
        assert!(read(&grid(POINTS, &field.replace("offset=\"0", "offset=\"1"), &block, options), 0.).is_err());
    }
}

#[test]
fn inline_zlib_with_separately_padded_header_and_multiple_blocks() {
    let raw: Vec<u8> = (0..8).flat_map(|x| (x as f32).to_le_bytes()).collect();
    let blocks: Vec<Vec<u8>> = raw
        .chunks(10)
        .map(|b| {
            let mut z = ZlibEncoder::new(Vec::new(), Compression::default());
            z.write_all(b).unwrap();
            z.finish().unwrap()
        })
        .collect();
    let mut header: Vec<u8> = [4_u32, 10, 2].into_iter().flat_map(u32::to_le_bytes).collect();
    for b in &blocks {
        header.extend_from_slice(&(b.len() as u32).to_le_bytes());
    }
    let data = format!("{}\n{}", STANDARD.encode(header), STANDARD.encode(blocks.concat()));
    let field = format!(r#"<DataArray type="Float32" Name="E" NumberOfComponents="2" format="binary">{data}</DataArray>"#);
    let opts = r#"compressor="vtkZLibDataCompressor""#;
    assert_eq!(read(&grid(POINTS, &field, b"", opts), 0.).unwrap().point_data["E"].values[7], 7.);
    // A declared expansion that disagrees with the grid is rejected before inflation.
    let malformed = field.replace("NumberOfComponents=\"2\"", "NumberOfComponents=\"3\"");
    assert!(read(&grid(POINTS, &malformed, b"", opts), 0.).is_err());
}

#[test]
fn appended_base64_offsets_refer_to_encoded_bytes() {
    let mut raw = 32_u32.to_le_bytes().to_vec();
    raw.extend((0..8).flat_map(|x| (x as f32).to_le_bytes()));
    let inline = STANDARD.encode(&raw);
    let prefix = STANDARD.encode(0_u32.to_le_bytes());
    let appended = format!(r#"<AppendedData encoding="base64">_{prefix}{inline}</AppendedData>"#);
    let field = format!(r#"<DataArray type="Float32" Name="E" NumberOfComponents="2" format="appended" offset="{}"/>"#, prefix.len());
    assert_eq!(read(&grid(POINTS, &field, appended.as_bytes(), ""), 0.).unwrap().point_data["E"].values[7], 7.);
}

#[test]
fn invalid_grids_report_the_cause() {
    let good = String::from_utf8(grid(POINTS, POINT_FIELD, b"", "")).unwrap();
    let oversized = format!(r#"NumberOfPoints="{}""#, super::MAX_POINTS + 1);
    for (from, to, message) in [
        ("3 6", "3 5", "offsets"),
        ("3 1 2", "4 1 2", "missing point"),
        ("5 5", "5 200", "cell type 200"),
        ("0 1 2 3 4 5 6 7", "0 1 2", "too few"),
        ("NumberOfPoints=\"4\"", oversized.as_str(), "NumberOfPoints"),
        ("0 1 2 3 4 5 6 7", "0 1 2 3 inf 5 6 7", "finite"),
    ] {
        assert!(read(good.replace(from, to).as_bytes(), 0.).unwrap_err().to_string().contains(message), "{message}");
    }
}

#[test]
fn collections_preserve_times_paths_and_reject_unsupported_members() {
    let doc = r#"<VTKFile type="Collection"><Collection><DataSet timestep="0" file="first.vtu"/><DataSet timestep="0.125" part="0" file="time/last.vtu"/></Collection></VTKFile>"#;
    assert_eq!(collection(doc.as_bytes()).unwrap(), [(0., "first.vtu".into()), (0.125, "time/last.vtu".into())]);
    for (from, to) in [
        ("0.125", "NaN"),
        ("0.125", "0"),
        ("part=\"0\"", "part=\"1\""),
        ("first.vtu", "../outside.vtu"),
        ("first.vtu", "/outside.vtu"),
        ("first.vtu", r"C:\outside.vtu"),
        ("first.vtu", "recursive.pvd"),
    ] {
        assert!(collection(doc.replace(from, to).as_bytes()).is_err());
    }
}

#[test]
fn declared_tuple_counts_and_exact_integer_conversion() {
    let field = r#"<DataArray type="Int64" Name="labels" NumberOfTuples="4">-5 0 9007199254740992 12</DataArray>"#;
    let frame = read(&grid(POINTS, field, b"", ""), 0.).unwrap();
    assert_eq!(frame.point_data["labels"].values[2], 9_007_199_254_740_992.);
    for (from, to, message) in [
        ("Tuples=\"4\"", "Tuples=\"3\"", "NumberOfTuples"),
        ("9007199254740992", "9007199254740993", "represented exactly"),
        ("9007199254740992", "9223372036854775807", "represented exactly"),
    ] {
        assert!(read(&grid(POINTS, &field.replace(from, to), b"", ""), 0.).unwrap_err().to_string().contains(message));
    }
    let unsigned = field.replace("Int64", "UInt64").replace("-5", "18446744073709551615");
    assert!(read(&grid(POINTS, &unsigned, b"", ""), 0.).unwrap_err().to_string().contains("represented exactly"));
}

#[test]
fn a_collection_frame_respects_its_remaining_storage() {
    let bytes = grid(POINTS, POINT_FIELD, b"", "");
    let charge = read(&bytes, 0.0).unwrap().storage_bytes();
    assert_eq!(super::read_with_budget(&bytes, 0.0, charge).unwrap().storage_bytes(), charge);
    assert!(matches!(super::read_with_budget(&bytes, 0.0, charge - 1), Err(crate::Error::Limit(_))));
}

#[test]
fn length_conversion_precedes_the_si_coordinate_limit_and_survives_embedding() {
    use fieldcraft_model::units::{Quantity, Unit};
    let mm = Unit::parse("mm").unwrap();
    let source = Quantity::Length.express(fieldcraft_model::MAX_COORD, &mm).unwrap();
    let positions = POINTS.replace("0 0 0 1 0 0", &format!("0 0 0 {source} 0 0"));
    let bytes = grid(&positions, POINT_FIELD, b"", "");
    assert!(read(&bytes, 0.0).is_err());
    let frame = super::read_with_unit(&bytes, 0.0, super::MAX_FRAME_BYTES, Some(&mm)).unwrap();
    assert_eq!(frame.grid.points[1][0], fieldcraft_model::MAX_COORD);
    assert_eq!(frame.point_data["E"].values[1], 1.0);
    let series = crate::external::Series { frames: vec![frame] };
    let encoded = crate::external::encode(&series).unwrap();
    assert_eq!(crate::external::decode(&encoded).unwrap(), series);
}

#[test]
fn large_connectivity_and_result_fields_fit_their_retained_storage() {
    let nc = 64;
    let points = (0..20).map(|i| format!("{i} 0 0")).collect::<Vec<_>>().join(" ");
    let cell = (0..20).map(|i| i.to_string()).collect::<Vec<_>>().join(" ") + " ";
    let offsets = (1..=nc).map(|i| (20 * i).to_string()).collect::<Vec<_>>().join(" ");
    for components in [1, 9] {
        let xml = format!(
            r#"<VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="20" NumberOfCells="{nc}">
        <Points><DataArray type="Float64" NumberOfComponents="3">{points}</DataArray></Points>
        <Cells><DataArray type="UInt8" Name="types">{}</DataArray><DataArray type="Int32" Name="offsets">{offsets}</DataArray>
        <DataArray type="Int32" Name="connectivity">{}</DataArray></Cells>
        <CellData><DataArray type="Float64" Name="T" NumberOfComponents="{components}">{}</DataArray></CellData>
        </Piece></UnstructuredGrid></VTKFile>"#,
            "25 ".repeat(nc),
            cell.repeat(nc),
            "1 ".repeat(components * nc)
        );
        let expected = read(xml.as_bytes(), 0.0).unwrap();
        let storage = expected.storage_bytes();
        let doc = super::document(&xml).unwrap();
        let vtk = doc.root_element();
        let historical_geometry = (3 * 20 + (2 + 20) * nc) * size_of::<f64>();
        if components == 1 {
            assert!(historical_geometry > storage);
        }
        assert!(historical_geometry + components * nc * size_of::<f64>() > storage);
        assert_eq!(super::decode_piece(vtk, None, 0.0, storage, storage, None).unwrap(), expected);
        assert!(matches!(super::decode_piece(vtk, None, 0.0, storage - 1, storage, None), Err(crate::Error::Limit(_))));
    }
}

#[test]
fn topology_workspace_counts_the_arrays_live_during_an_index_read() {
    let doc = super::document(r#"<DataArray type="UInt64">1 4294967295</DataArray>"#).unwrap();
    let node = doc.root_element();
    // Two held type indices and two decoded offsets each use a u32 slot.
    let cap = 4 * size_of::<u32>();
    let full = super::array::Decoder::new(node, None, cap).unwrap();
    assert_eq!(full.indices(node, 2, 2).unwrap(), [1, u32::MAX]);
    let short = super::array::Decoder::new(node, None, cap - 1).unwrap();
    assert!(matches!(short.indices(node, 2, 2), Err(crate::Error::Limit(_))));
    assert_eq!(short.indices(node, 2, 0).unwrap(), [1, u32::MAX]);

    // An uncompressed raw block borrows its input; a zlib block also holds the
    // inflated bytes and one sentinel byte while converting its integers.
    let doc = super::document(
        r#"<VTKFile compressor="vtkZLibDataCompressor"><AppendedData encoding="raw"/><DataArray type="UInt32" format="appended" offset="0"/></VTKFile>"#,
    )
    .unwrap();
    let root = doc.root_element();
    let node = super::one(root, "DataArray").unwrap();
    let mut z = ZlibEncoder::new(Vec::new(), Compression::default());
    z.write_all(&[1_u32.to_le_bytes(), u32::MAX.to_le_bytes()].concat()).unwrap();
    let block = z.finish().unwrap();
    let mut payload: Vec<u8> = [1_u32, 8, 8, block.len() as u32].into_iter().flat_map(u32::to_le_bytes).collect();
    payload.extend(block);
    let peak = cap + 8 + 1;
    let full = super::array::Decoder::new(root, Some(&payload), peak).unwrap();
    assert_eq!(full.indices(node, 2, 2).unwrap(), [1, u32::MAX]);
    let short = super::array::Decoder::new(root, Some(&payload), peak - 1).unwrap();
    assert!(matches!(short.indices(node, 2, 2), Err(crate::Error::Limit(_))));
}

#[test]
fn topology_indices_preserve_binary_byte_order_and_integer_checks() {
    for little in [false, true] {
        let mut data = if little { 4_u32.to_le_bytes() } else { 4_u32.to_be_bytes() }.to_vec();
        for index in [1_u16, 513] {
            data.extend(if little { index.to_le_bytes() } else { index.to_be_bytes() });
        }
        let order = if little { "LittleEndian" } else { "BigEndian" };
        let xml = format!(r#"<DataArray byte_order="{order}" type="UInt16" format="binary">{}</DataArray>"#, STANDARD.encode(data));
        let doc = super::document(&xml).unwrap();
        let node = doc.root_element();
        assert_eq!(super::array::Decoder::new(node, None, 64).unwrap().indices(node, 2, 0).unwrap(), [1, 513]);
    }
    for (kind, value, message) in [
        ("Int8", "-1", "negative or oversized"),
        ("UInt64", "4294967296", "negative or oversized"),
        ("UInt64", "9007199254740993", "represented exactly"),
        ("Float64", "1", "requires integer arrays"),
    ] {
        let xml = format!(r#"<DataArray type="{kind}">{value}</DataArray>"#);
        let doc = super::document(&xml).unwrap();
        let node = doc.root_element();
        assert!(super::array::Decoder::new(node, None, 64).unwrap().indices(node, 1, 0).unwrap_err().to_string().contains(message));
    }
}

#[test]
fn a_pvtu_indexes_exactly_one_relative_serial_piece() {
    let doc = r#"<VTKFile type="PUnstructuredGrid"><PUnstructuredGrid GhostLevel="0"><PPoints><PDataArray type="Float64" NumberOfComponents="3"/></PPoints><Piece Source="step/grid.VTU"/></PUnstructuredGrid></VTKFile>"#;
    assert_eq!(single_piece(doc.as_bytes()).unwrap(), "step/grid.VTU");
    for name in ["../escape.vtu", "/outside.vtu", r"C:\temp\x.vtu", "nested.pvtu", "index.pvd"] {
        assert!(single_piece(doc.replace("step/grid.VTU", name).as_bytes()).is_err());
    }
    let extra = doc.replace("</PUnstructuredGrid>", r#"<Piece Source="part1.vtu"/></PUnstructuredGrid>"#);
    assert!(single_piece(extra.as_bytes()).unwrap_err().to_string().contains("multi-piece"));
    assert!(single_piece(doc.replace(r#"<Piece Source="step/grid.VTU"/>"#, "").as_bytes()).is_err());
}

#[test]
fn xml_markup_does_not_change_array_values_or_identify_a_payload() {
    let ascii = String::from_utf8(grid(POINTS, POINT_FIELD, b"", "")).unwrap();
    let xml = format!("<!-- <AppendedData> --><?note <AppendedData>?>{ascii}")
        .replace("0 1 2 3 4 5 6 7", "0 1 2 3 4 5 6 1<!--split--><![CDATA[7]]>")
        .replace(">3 6<", ">3<!--split--> 6<");
    assert_eq!(read(xml.as_bytes(), 0.).unwrap().point_data["E"].values[7], 17.);
    let raw: Vec<u8> = 32_u32.to_le_bytes().into_iter().chain((0..8).flat_map(|n| (n as f32).to_le_bytes())).collect();
    let encoded = STANDARD.encode(&raw);
    let (a, b) = encoded.split_at(5);
    let field = format!(r#"<DataArray type="Float32" Name="E" NumberOfComponents="2" format="binary">{a}<!-- split a quartet -->{b}</DataArray>"#);
    assert_eq!(read(&grid(POINTS, &field, b"", ""), 0.).unwrap().point_data["E"].values[7], 7.);

    // An XML trailer can itself mention closing tags. Appended bytes have no XML encoding.
    let mut appended = br#"<!-- <AppendedData> --><AppendedData hint="a > b" encoding="raw">_"#.to_vec();
    appended.extend(raw);
    appended.extend_from_slice(b"</AppendedData   ><!-- </AppendedData> -->");
    let from_raw = r#"<DataArray type="Float32" Name="E" NumberOfComponents="2" format="appended" offset="0"/>"#;
    assert_eq!(read(&grid(POINTS, from_raw, &appended, ""), 0.).unwrap().point_data["E"].values[7], 7.);
}

#[test]
fn lagrange_cells_keep_their_point_order_and_all_samples_in_project_storage() {
    // The wedge uses 0,1,2 on one triangular face and 3,4,5 directly above them.
    // The quadratic tetrahedron has edge samples on 01,12,20,03,13,23, in that order.
    let tetrahedron = [0, 1, 3, 4, 8, 9, 10, 11, 12, 13];
    let cases: &[(u32, u8, &[u32])] = &[
        (68, 3, &[0, 1]),
        (68, 21, &[0, 1, 8]),
        (69, 5, &[0, 1, 3]),
        (69, 22, &[0, 1, 3, 8, 9, 10]),
        (70, 9, &[0, 1, 2, 3]),
        (71, 10, &tetrahedron[..4]),
        (71, 24, &tetrahedron),
        (72, 12, &[0, 1, 2, 3, 4, 5, 6, 7]),
        (73, 13, &[0, 1, 3, 4, 5, 7]),
    ];
    let points = "0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1  0.5 0 0  0.5 0.5 0  0 0.5 0  0 0 0.5  0.5 0 0.5  0 0.5 0.5";
    let names = |values: &[u32]| values.iter().map(u32::to_string).collect::<Vec<_>>().join(" ");
    let ids = cases.iter().flat_map(|(_, _, nodes)| nodes.iter().copied()).collect::<Vec<_>>();
    let offsets = cases
        .iter()
        .scan(0, |n, (_, _, nodes)| {
            *n += nodes.len() as u32;
            Some(*n)
        })
        .collect::<Vec<_>>();
    let types = cases.iter().map(|c| c.0).collect::<Vec<_>>();
    let xml = format!(
        r#"<VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="14" NumberOfCells="9">
      <Points><DataArray type="Float64" NumberOfComponents="3">{points}</DataArray></Points>
      <Cells><DataArray type="Int32" Name="connectivity">{}</DataArray><DataArray type="Int32" Name="offsets">{}</DataArray>
      <DataArray type="UInt8" Name="types">{}</DataArray></Cells>
      <PointData><DataArray type="UInt8" Name="sample">0 1 2 3 4 5 6 7 8 9 10 11 12 13</DataArray></PointData></Piece></UnstructuredGrid></VTKFile>"#,
        names(&ids),
        names(&offsets),
        names(&types)
    );
    let frame = read(xml.as_bytes(), 0.).unwrap();
    for (cell, (_, fixed, source)) in frame.grid.cells.iter().zip(cases) {
        assert_eq!(cell.vtk_type, *fixed);
        assert_eq!(cell.points, *source);
    }
    let series = crate::external::Series { frames: vec![frame] };
    assert_eq!(crate::external::decode(&crate::external::encode(&series).unwrap()).unwrap(), series);
    for (kind, n) in [(68, 4), (69, 10), (70, 8), (70, 9), (71, 20), (72, 20), (72, 27), (73, 18)] {
        assert!(super::canonical_cell(kind, n).unwrap_err().to_string().contains("unsupported Lagrange"));
    }
}
