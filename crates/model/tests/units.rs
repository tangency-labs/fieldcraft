//! The unit table against hand-computed values.

use fieldcraft_model::units::{DisplaySystem, Preset, Quantity, Unit, UnitError, ValueInput};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1e-300)
}

fn si(q: Quantity, text: &str) -> f64 {
    q.parse(text, &DisplaySystem::new(Preset::MKgS)).unwrap_or_else(|e| panic!("{text}: {e}"))
}

fn shown(q: Quantity, si_value: f64, preset: Preset) -> f64 {
    q.express(si_value, &DisplaySystem::new(preset).unit(q).unwrap()).unwrap()
}

#[test]
fn hand_values() {
    use Quantity as Q;
    // E = 210 GPa = 210 000 MPa = 2.1e11 Pa.
    let e = si(Q::Pressure, "210 GPa");
    assert!(close(e, 2.1e11));
    assert!(close(si(Q::Pressure, "210000 MPa"), e));
    assert!(close(shown(Q::Pressure, e, Preset::MmTonneS), 210_000.0));
    // 1 N/mm² = 1 MPa.
    assert!(close(si(Q::Pressure, "1 N/mm²"), si(Q::Pressure, "1 MPa")));
    // Steel density 7850 kg/m³ = 7.85e-9 t/mm³.
    let rho = si(Q::Density, "7850 kg/m³");
    assert!(close(shown(Q::Density, rho, Preset::MmTonneS), 7.85e-9));
    assert!(close(si(Q::Density, "7.85e-9 t/mm^3"), 7850.0));
    // 20 °C = 293.15 K, but a difference of 20 °C is 20 K.
    assert!(close(si(Q::Temperature, "20 °C"), 293.15));
    assert!(close(si(Q::TemperatureDifference, "20 °C"), 20.0));
    assert!(close(si(Q::Temperature, "68 °F"), 293.15));
    assert!(close(si(Q::TemperatureDifference, "9 °F"), 5.0));
    assert!(close(shown(Q::Temperature, 293.15, Preset::MmTonneS), 20.0));
    assert!(close(shown(Q::Temperature, 293.15, Preset::InLbfS), 68.0));
    // Expansion 12e-6 /°C = 12e-6 /K; 6.667e-6 /°F.
    assert!(close(si(Q::ThermalExpansion, "12e-6 °C⁻¹"), 12e-6));
    assert!(close(shown(Q::ThermalExpansion, 12e-6, Preset::InLbfS), 12e-6 * 5.0 / 9.0));
    // 25 mm/s = 0.025 m/s; g = 9.81 m/s² = 9810 mm/s² = 386.22 in/s².
    assert!(close(si(Q::Velocity, "25 mm/s"), 0.025));
    assert!(close(shown(Q::Acceleration, 9.81, Preset::MmTonneS), 9810.0));
    assert!(close(shown(Q::Acceleration, 9.81, Preset::InLbfS), 9.81 / 0.0254));
    // Steel in in-lbf-s: E = 2.1e11 Pa = 30.458 Msi; ρ = 7850 kg/m³ = 7.345e-4 lbf·s²/in⁴.
    let psi = 4.448_221_615_260_5 / (0.0254 * 0.0254);
    assert!(close(shown(Q::Pressure, e, Preset::InLbfS), 2.1e11 / psi));
    assert!(close(shown(Q::Density, rho, Preset::InLbfS), 7850.0 * 0.0254 / 4.448_221_615_260_5 * 0.0254 * 0.0254 * 0.0254));
    // Steel conductivity 50 W/(m·K) is 50 mW/(mm·K); specific heat 460 J/(kg·K) is 4.6e8 mJ/(t·K).
    assert!(close(shown(Q::ThermalConductivity, 50.0, Preset::MmTonneS), 50.0));
    assert!(close(shown(Q::SpecificHeat, 460.0, Preset::MmTonneS), 4.6e8));
    // 10 W/(m²·K) = 0.01 mW/(mm²·K); 1 kW/m² = 1 mW/mm².
    assert!(close(shown(Q::HeatTransferCoefficient, 10.0, Preset::MmTonneS), 0.01));
    assert!(close(shown(Q::HeatFlux, 1000.0, Preset::MmTonneS), 1.0));
    // Water viscosity 1e-3 Pa·s = 1e-9 MPa·s; 1 N·mm = 1e-3 N·m = 1 mJ.
    assert!(close(shown(Q::DynamicViscosity, 1e-3, Preset::MmTonneS), 1e-9));
    assert!(close(si(Q::Moment, "1 N·mm"), si(Q::Energy, "1 mJ")));
    // 90 deg = π/2 rad; 5 % = 0.05.
    assert!(close(si(Q::Angle, "90 deg"), std::f64::consts::FRAC_PI_2));
    assert!(close(si(Q::Dimensionless, "5 %"), 0.05));
}

#[test]
fn every_preset_unit_has_its_quantitys_dimension() {
    for preset in Preset::ALL {
        let display = DisplaySystem::new(preset);
        for &q in Quantity::ALL {
            let unit = display.unit(q).unwrap_or_else(|e| panic!("{preset:?} {q:?}: {e}"));
            assert_eq!(unit.dimension(), q.dimension(), "{preset:?} {q:?} `{}`", unit.symbol());
        }
        assert_eq!(Preset::from_label(preset.label()), Some(preset));
    }
}

#[test]
fn mm_tonne_second_is_consistent() {
    // A consistent system: force = mass · acceleration and stress = force / area hold with no
    // factor, so the deck writers can write its numbers directly.
    let d = DisplaySystem::new(Preset::MmTonneS);
    let f = |q| d.unit(q).unwrap().factor();
    use Quantity as Q;
    assert!(close(f(Q::Force), f(Q::Mass) * f(Q::Acceleration)));
    assert!(close(f(Q::Pressure), f(Q::Force) / f(Q::Area)));
    assert!(close(f(Q::Energy), f(Q::Force) * f(Q::Length)));
    assert!(close(f(Q::Density), f(Q::Mass) / f(Q::Volume)));
    let i = DisplaySystem::new(Preset::InLbfS);
    let g = |q| i.unit(q).unwrap().factor();
    assert!(close(g(Q::Force), g(Q::Mass) * g(Q::Acceleration)));
    assert!(close(g(Q::Pressure), g(Q::Force) / g(Q::Area)));
}

#[test]
fn round_trips_through_every_preset() {
    let values = [(Quantity::Pressure, 2.1e11), (Quantity::Density, 7850.0), (Quantity::Temperature, 293.15), (Quantity::Length, 0.0254)];
    for preset in Preset::ALL {
        let d = DisplaySystem::new(preset);
        for (q, v) in values {
            let text = q.format(v, &d).unwrap();
            let back = q.parse(&text, &d).unwrap();
            assert!((back - v).abs() <= 1e-5 * v.abs(), "{preset:?} {q:?}: {text} -> {back}");
        }
    }
    let mm = DisplaySystem::default();
    assert_eq!(Quantity::Pressure.format(2.1e11, &mm).unwrap(), "210000 MPa");
    assert_eq!(Quantity::Density.format(7850.0, &mm).unwrap(), "7.85e-9 t/mm³");
    assert_eq!(Quantity::Temperature.format(293.15, &mm).unwrap(), "20 °C");
    assert_eq!(Quantity::Dimensionless.format(0.3, &mm).unwrap(), "0.3");
}

#[test]
fn inputs_from_commands() {
    let mm = DisplaySystem::default();
    let parse = |json: &str| serde_json::from_str::<ValueInput>(json).unwrap().to_si(Quantity::Pressure, &mm);
    assert!(close(parse("210000").unwrap(), 2.1e11)); // bare number: display system (MPa)
    assert!(close(parse(r#""210 GPa""#).unwrap(), 2.1e11));
    assert!(close(parse(r#"{"value": 210, "unit": "GPa"}"#).unwrap(), 2.1e11));
    assert!(close(Quantity::Pressure.parse("210000", &mm).unwrap(), 2.1e11));
    assert!(serde_json::from_str::<ValueInput>(r#"{"value": 1, "unit": "Pa", "extra": 0}"#).is_err());
    let over = DisplaySystem::default().with_override(Quantity::Energy, "J").unwrap();
    assert_eq!(Quantity::Energy.format(2.0, &over).unwrap(), "2 J");
    assert!(DisplaySystem::default().with_override(Quantity::Energy, "MPa").is_err());
}

#[test]
fn hostile_values_are_refused() {
    let mm = DisplaySystem::default();
    use Quantity as Q;
    let wrong = Q::Pressure.parse("210 mm", &mm);
    assert!(matches!(wrong, Err(UnitError::WrongDimension { .. })), "{wrong:?}");
    assert_eq!(wrong.unwrap_err().to_string(), "`mm` is L, but a pressure is L^-1 M T^-2");
    for text in ["", "GPa", "nan GPa", "inf Pa", "1e400 Pa", "1e300 GPa", "210 GPa extra", "210 parsecs", "--1 Pa"] {
        assert!(Q::Pressure.parse(text, &mm).is_err(), "{text:?}");
    }
    assert!(Q::Density.parse("1e300 t/mm³", &mm).is_err());
    assert!(matches!(Q::Temperature.parse("-300 °C", &mm), Err(UnitError::BelowAbsoluteZero(_))));
    assert!(Q::TemperatureDifference.parse("-300 °C", &mm).is_ok());
    assert!(Q::Pressure.to_si(f64::NAN, &Unit::parse("Pa").unwrap()).is_err());
    assert!(Q::Pressure.parse(&"9".repeat(500), &mm).is_err());
    let input = ValueInput::WithUnit(fieldcraft_model::units::WithUnit { value: f64::INFINITY, unit: "Pa".into() });
    assert!(input.to_si(Q::Pressure, &mm).is_err());
}

#[test]
fn angles_and_ratios_are_kept_apart() {
    let mm = DisplaySystem::default();
    use Quantity as Q;
    assert!(matches!(Q::Angle.parse("5 %", &mm), Err(UnitError::WrongKind { .. })));
    assert!(matches!(Q::Dimensionless.parse("1 deg", &mm), Err(UnitError::WrongKind { .. })));
    assert!(matches!(Q::Angle.parse("90 mm", &mm), Err(UnitError::WrongDimension { .. })));
    assert!(Q::Pressure.parse("1 MPa·deg", &mm).is_err());
    assert!(close(Q::Angle.parse("90 °", &mm).unwrap(), std::f64::consts::FRAC_PI_2));
    assert!(close(Q::Angle.parse("1.5", &mm).unwrap(), 1.5)); // display: rad, as rotations always were
    assert!(close(Q::AngularVelocity.parse("60 rpm", &mm).unwrap(), std::f64::consts::TAU));
    assert!(close(Q::AngularVelocity.parse("180 deg/s", &mm).unwrap(), std::f64::consts::PI));
    assert!(matches!(Q::AngularVelocity.parse("1 Hz", &mm), Err(UnitError::WrongKind { .. })));
    assert!(matches!(Q::Frequency.parse("1 rad/s", &mm), Err(UnitError::WrongKind { .. })));
    assert!(close(Q::RotationalStiffness.parse("2 N·m/rad", &mm).unwrap(), 2.0));
    assert!(close(Q::RotationalStiffness.parse("2", &mm).unwrap(), 2e-3)); // N·mm/rad
    assert!(close(Q::Dimensionless.parse("0.3", &mm).unwrap(), 0.3));
    assert!(close(Q::ThermalExpansion.parse("12e-6 1/K", &mm).unwrap(), 12e-6));
    assert!(close(Q::Frequency.parse("0.01 1/s", &mm).unwrap(), 0.01));
    assert!(close(Q::ThermalConductivity.parse("50 W/(m K)", &mm).unwrap(), 50.0));
}

/// A value typed with up to 15 significant digits in a metric preset converts to SI and back to
/// exactly what was typed, and its SI value is the double nearest the exact decimal (7.85e-9
/// t/mm³ is 7850 kg/m³, not 7850.000000000001).
#[test]
fn typed_values_convert_exactly_in_metric_presets() {
    use Quantity as Q;
    for preset in [Preset::MKgS, Preset::MmTonneS, Preset::MmKgMs, Preset::MmGMs] {
        let d = DisplaySystem::new(preset);
        for &q in Quantity::ALL {
            for i in 1..2000 {
                for v in [f64::from(i) / 100.0, -f64::from(i) * 7.0, format!("{i}e-9").parse::<f64>().unwrap()] {
                    // Absolute temperatures stop at absolute zero.
                    let v = if q == Q::Temperature { v.abs() } else { v };
                    let si = d.to_si(q, v).unwrap();
                    assert_eq!(d.express(q, si).unwrap(), v, "{preset:?} {q:?} {v}");
                }
            }
        }
    }
    let mm = DisplaySystem::default();
    assert_eq!(mm.to_si(Q::Density, 7.85e-9).unwrap(), 7850.0);
    assert_eq!(mm.to_si(Q::Length, 3.97).unwrap(), 0.00397);
    assert_eq!(mm.to_si(Q::Pressure, 210000.0).unwrap(), 2.1e11);
    assert_eq!(mm.to_si(Q::Temperature, 20.0).unwrap(), 293.15);
    assert_eq!(mm.to_si(Q::SpecificHeat, 4.6e8).unwrap(), 460.0);
}

/// °F goes through °C in decimal: whole degrees convert exactly both ways (80 °C is 176 °F, a
/// legend range of 68 … 176 °F is 293.15 … 353.15 K), and every whole °F comes back as typed.
#[test]
fn fahrenheit_is_exact_through_celsius() {
    use Quantity as Q;
    let f = DisplaySystem::new(Preset::InLbfS);
    let c = DisplaySystem::new(Preset::MmTonneS);
    assert_eq!(f.express(Q::Temperature, c.to_si(Q::Temperature, 80.0).unwrap()).unwrap(), 176.0);
    assert_eq!(f.to_si(Q::Temperature, 68.0).unwrap(), 293.15);
    assert_eq!(f.to_si(Q::Temperature, 176.0).unwrap(), 353.15);
    for whole in -273..=2000 {
        let t = f64::from(whole);
        let k = c.to_si(Q::Temperature, t).unwrap();
        let shown = f.express(Q::Temperature, k).unwrap();
        assert_eq!(c.to_si(Q::Temperature, c.express(Q::Temperature, f.to_si(Q::Temperature, shown).unwrap()).unwrap()).unwrap(), k, "{t} °C");
        assert_eq!(f.express(Q::Temperature, f.to_si(Q::Temperature, t).unwrap()).unwrap(), t, "{t} °F");
    }
}

/// A temperature override carries to the temperature difference and the per-degree units unless
/// they have their own: mm-t-s with °F shows ΔT in °F and an expansion per °F.
#[test]
fn a_temperature_override_carries_to_the_per_degree_units() {
    use Quantity as Q;
    let d = DisplaySystem::new(Preset::MmTonneS).with_override(Q::Temperature, "°F").unwrap();
    let symbols: Vec<&str> = [Q::TemperatureDifference, Q::ThermalExpansion, Q::ThermalConductivity, Q::SpecificHeat, Q::HeatTransferCoefficient]
        .iter()
        .map(|q| d.symbol(*q))
        .collect();
    assert_eq!(symbols, ["°F", "°F⁻¹", "mW/(mm·°F)", "mJ/(t·°F)", "mW/(mm²·°F)"]);
    // 1.8 °F of difference is 1 K; 1e-5 per K is 5.555…e-6 per °F.
    assert!((d.to_si(Q::TemperatureDifference, 1.8).unwrap() - 1.0).abs() < 1e-15);
    assert!((d.express(Q::ThermalExpansion, 1e-5).unwrap() - 1e-5 / 1.8).abs() < 1e-20);
    // Unless overridden itself; and nothing else changes.
    let own = d.clone().with_override(Q::ThermalExpansion, "1/K").unwrap();
    assert_eq!((own.symbol(Q::ThermalExpansion), own.symbol(Q::TemperatureDifference), own.symbol(Q::Pressure)), ("1/K", "°F", "MPa"));
    // SI shown in °C: per-degree units in °C (the same step as K).
    let si = DisplaySystem::new(Preset::MKgS).with_override(Q::Temperature, "°C").unwrap();
    assert_eq!((si.symbol(Q::ThermalConductivity), si.to_si(Q::TemperatureDifference, 3.0).unwrap()), ("W/(m·°C)", 3.0));
}
