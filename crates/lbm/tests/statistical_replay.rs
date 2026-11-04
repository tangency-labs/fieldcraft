//! Own FieldCraft cb976117 force histories reported in #796, not validation fixtures.
use fieldcraft_lbm::statistics::Blocks;

#[test]
fn audited_cylinder_windows_do_not_pass_statistical_acceptance() {
    for (name, history, expected) in [
        ("100–200", include_str!("data/statistics/cylinder-early.csv"), 1.3539),
        ("400–600", include_str!("data/statistics/cylinder-late.csv"), 1.6602),
        ("potential start, 400–600", include_str!("data/statistics/cylinder-potential-late.csv"), 1.4818),
    ] {
        let values: Vec<f64> = history.lines().map(|s| s.parse().unwrap()).collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        assert!((mean - expected).abs() < 0.0001, "{name}: {mean}");
        for len in [25, 50] {
            let mut blocks = Blocks::new(len, 0).unwrap();
            for &v in &values {
                blocks.push(v).unwrap();
            }
            let r = blocks.acceptance();
            println!("audit {name}, block {len}: {r:?}");
            assert!(r.analysed_blocks >= 20);
            assert!(!r.accepted, "{name} at {len}: {r:?}");
        }
    }
}
