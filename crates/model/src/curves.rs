//! Material properties that follow the temperature: tables of (absolute temperature, SI value)
//! interpolated linearly between their points and held at the end values outside them.
//!
//! A material keeps its constants beside the tables (the library takes them at a temperature the
//! user chooses), so a solver that does not follow a table yet still has a value; the heat
//! transfer solver reads the conductivity, specific heat and density tables, the fluid solver the
//! viscosity table (with heat in the fluid domain). The tables are saved with the model, so a
//! result does not change when the library is updated.

use serde::{Deserialize, Serialize};

use crate::{Error, Result, positive};

/// Most points of a temperature table.
pub const MAX_POINTS: usize = 256;

/// The value of the `[x, y]` points (x increasing) at `x`, linear between them and held at the
/// end values outside them (a NaN `x` reads the first point, no points read 0): the interval is
/// found by bisection, so a lookup costs O(log n).
pub fn sorted_value(points: &[[f64; 2]], x: f64) -> f64 {
    let (Some(first), Some(last)) = (points.first(), points.last()) else { return 0.0 };
    if !(x > first[0]) {
        return first[1];
    }
    if x >= last[0] {
        return last[1];
    }
    // The first point at or beyond x, past the first since x is beyond it.
    let i = points.partition_point(|p| p[0] < x);
    match (points.get(i.wrapping_sub(1)), points.get(i)) {
        (Some(a), Some(b)) => a[1] + (x - a[0]) / (b[0] - a[0]) * (b[1] - a[1]),
        _ => last[1],
    }
}

/// A property as a function of the absolute temperature: `[kelvin, value]` points with strictly
/// increasing temperatures, linear between them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct TemperatureTable(pub Vec<[f64; 2]>);

impl TemperatureTable {
    /// The value at `kelvin`, held at the end values outside the table (a NaN temperature reads
    /// the lowest point).
    pub fn value(&self, kelvin: f64) -> f64 {
        sorted_value(&self.0, kelvin)
    }

    /// Whether `kelvin` lies within the table's temperatures.
    pub fn covers(&self, kelvin: f64) -> bool {
        let [a, b] = self.range();
        (a..=b).contains(&kelvin)
    }

    /// The temperature range `[lowest, highest]` (kelvin).
    pub fn range(&self) -> [f64; 2] {
        let t = |p: Option<&[f64; 2]>| p.map_or(f64::NAN, |p| p[0]);
        [t(self.0.first()), t(self.0.last())]
    }

    /// Check `what`'s table: 2 to [`MAX_POINTS`] points, finite temperatures at or above zero
    /// kelvin rising strictly, every value finite and positive.
    pub fn validate(&self, what: &str) -> Result<()> {
        let n = self.0.len();
        if !(2..=MAX_POINTS).contains(&n) {
            return Err(Error::Invalid(format!("the {what} table needs 2 to {MAX_POINTS} points, got {n}")));
        }
        let mut previous = f64::NEG_INFINITY;
        for &[t, v] in &self.0 {
            if !(t.is_finite() && t >= 0.0 && t > previous) {
                return Err(Error::Invalid(format!(
                    "the {what} table's temperatures must be finite, at or above zero kelvin and strictly increasing, got {t} K after {previous} K"
                )));
            }
            positive(&format!("{what} table value at {t} K"), v)?;
            previous = t;
        }
        Ok(())
    }
}

/// The temperature tables of a material; each replaces its constant in the solvers that follow
/// it (the constant must still be given).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TemperatureDependent {
    /// Isotropic conductivity k(T).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<TemperatureTable>,
    /// Specific heat c(T).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<TemperatureTable>,
    /// Density ρ(T).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub density: Option<TemperatureTable>,
    /// Dynamic viscosity μ(T) of a Newtonian fluid.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub viscosity: Option<TemperatureTable>,
}

impl TemperatureDependent {
    /// The tables by name, in a fixed order.
    pub fn tables(&self) -> impl Iterator<Item = (&'static str, &TemperatureTable)> {
        [("conductivity", &self.conductivity), ("specific heat", &self.specific_heat), ("density", &self.density), ("viscosity", &self.viscosity)]
            .into_iter()
            .filter_map(|(what, t)| t.as_ref().map(|t| (what, t)))
    }

    /// The tables mutably, for unit conversion: the value quantity of each.
    pub(crate) fn tables_mut(&mut self) -> impl Iterator<Item = (crate::units::Quantity, &mut TemperatureTable)> {
        use crate::units::Quantity as Q;
        [
            (Q::ThermalConductivity, &mut self.conductivity),
            (Q::SpecificHeat, &mut self.specific_heat),
            (Q::Density, &mut self.density),
            (Q::DynamicViscosity, &mut self.viscosity),
        ]
        .into_iter()
        .filter_map(|(q, t)| t.as_mut().map(|t| (q, t)))
    }

    /// Whether the heat transfer solver follows a table of this material.
    pub fn conducts(&self) -> bool {
        self.conductivity.is_some() || self.specific_heat.is_some() || self.density.is_some()
    }

    pub fn is_empty(&self) -> bool {
        self.tables().next().is_none()
    }

    /// The tables the temperatures `[low, high]` (kelvin) leave, where their end values are
    /// held; `material` names it. Data in SI: the caller words it in its own display system
    /// ([`TableRange::message`]).
    pub fn range(&self, material: &str, low: f64, high: f64) -> Option<TableRange> {
        let tables: Vec<(String, [f64; 2])> = self
            .tables()
            .filter_map(|(what, t)| {
                let [a, b] = t.range();
                (low < a || high > b).then(|| (what.to_string(), [a, b]))
            })
            .collect();
        (!tables.is_empty()).then(|| TableRange { material: material.to_string(), tables, reached: [low, high] })
    }
}

/// Temperatures a material's tables were read at beyond their ranges, where the end values are
/// held: an extrapolation the source does not cover. All temperatures in kelvin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TableRange {
    pub material: String,
    /// Each table left, by property, with its range.
    pub tables: Vec<(String, [f64; 2])>,
    /// The lowest and highest temperatures the material reached.
    pub reached: [f64; 2],
}

impl TableRange {
    /// The warning, each temperature written by `show` (kelvin in, value and unit out).
    pub fn message(&self, show: impl Fn(f64) -> String) -> String {
        let tables: Vec<String> = self.tables.iter().map(|(what, [a, b])| format!("{what} ({} to {})", show(*a), show(*b))).collect();
        format!(
            "material `{}` reaches {} to {}, outside its {} table range: the end values are held there, an extrapolation the source does not cover",
            self.material,
            show(self.reached[0]),
            show(self.reached[1]),
            tables.join(", ")
        )
    }
}

/// The warning in kelvin.
impl std::fmt::Display for TableRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message(|k| format!("{} K", crate::units::format_number(k))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolates_linearly_and_holds_the_ends() {
        let t = TemperatureTable(vec![[300.0, 1.0], [400.0, 3.0], [500.0, 2.0]]);
        t.validate("k").unwrap();
        for (kelvin, value) in [(300.0, 1.0), (350.0, 2.0), (400.0, 3.0), (450.0, 2.5), (500.0, 2.0)] {
            assert_eq!((t.value(kelvin), t.covers(kelvin)), (value, true));
        }
        assert_eq!((t.value(200.0), t.covers(200.0)), (1.0, false));
        assert_eq!((t.value(900.0), t.covers(900.0)), (2.0, false));
        assert_eq!(t.value(f64::NAN), 1.0);
        assert!(!t.covers(f64::NAN));
        assert_eq!(t.value(f64::INFINITY), 2.0);
        assert_eq!(t.range(), [300.0, 500.0]);
    }

    #[test]
    fn the_bisection_reads_what_the_linear_scan_reads() {
        let points: Vec<[f64; 2]> = (0..MAX_POINTS).map(|i| [200.0 + 3.0 * i as f64 + 0.01 * (i * i) as f64, (0.37 * i as f64).sin()]).collect();
        let mut probes = vec![f64::NAN, f64::NEG_INFINITY, f64::INFINITY, -1.0, 1e9];
        probes.extend(points.iter().map(|p| p[0]));
        probes.extend(points.windows(2).map(|w| 0.5 * (w[0][0] + w[1][0])));
        probes.extend(points.windows(2).map(|w| w[1][0] - 1e-9));
        for n in [0, 1, 2, 3, MAX_POINTS] {
            let table = &points[..n];
            for &x in &probes {
                let (got, want) = (sorted_value(table, x), crate::fluid_zones::curve_value(table, x));
                assert_eq!(got.to_bits(), want.to_bits(), "{n} points at {x}: {got} vs {want}");
            }
        }
    }

    #[test]
    fn malformed_tables_are_refused() {
        for bad in [
            vec![[300.0, 1.0]],
            vec![[300.0, 1.0], [300.0, 2.0]],
            vec![[400.0, 1.0], [300.0, 2.0]],
            vec![[-1.0, 1.0], [300.0, 2.0]],
            vec![[300.0, 1.0], [f64::INFINITY, 2.0]],
            vec![[300.0, 0.0], [400.0, 2.0]],
            vec![[300.0, f64::NAN], [400.0, 2.0]],
            (0..=MAX_POINTS).map(|i| [i as f64, 1.0]).collect(),
        ] {
            assert!(TemperatureTable(bad).validate("k").is_err());
        }
    }

    #[test]
    fn a_range_warning_names_the_tables_left() {
        let d = TemperatureDependent {
            conductivity: Some(TemperatureTable(vec![[300.0, 1.0], [400.0, 1.0]])),
            density: Some(TemperatureTable(vec![[200.0, 1.0], [600.0, 1.0]])),
            ..TemperatureDependent::default()
        };
        assert!(d.range("Steel", 310.0, 390.0).is_none());
        let r = d.range("Steel", 290.0, 390.0).unwrap();
        assert_eq!(r, TableRange { material: "Steel".into(), tables: vec![("conductivity".into(), [300.0, 400.0])], reached: [290.0, 390.0] });
        let w = r.to_string();
        assert!(w.contains("290 K to 390 K") && w.contains("conductivity (300 K to 400 K)") && !w.contains("density"), "{w}");
        // The caller's own formatting, as a command's display system does.
        let w = r.message(|k| format!("{:.2} °C", k - 273.15));
        assert!(w.contains("16.85 °C to 116.85 °C") && w.contains("(26.85 °C to 126.85 °C)"), "{w}");
    }
}
