//! The geometry frame: SolidCraft designs are in millimetres, and the geometry commands and the
//! meshers work in them (mm-t-s-N-MPa) around a model converted from SI. An absolute amount of the
//! frame is written with its unit ([`Typed`]) and expressed in it here.

use std::sync::OnceLock;

use super::{DisplaySystem, Preset, Quantity, Typed};

/// The system the geometry and the meshers work in: millimetres, in mm-t-s-N-MPa.
pub fn geometry_frame() -> &'static DisplaySystem {
    static MM: OnceLock<DisplaySystem> = OnceLock::new();
    MM.get_or_init(|| DisplaySystem::new(Preset::MmTonneS))
}

/// An absolute amount of the geometry frame, written with its unit (`Typed::new(1e-9, "mm")`):
/// its value in the frame's millimetres.
pub fn geometry_amount(q: Quantity, t: Typed) -> f64 {
    t.express(q, geometry_frame()).unwrap_or(t.value)
}

/// The length the geometry frame takes for no length at all: a floor that keeps a size off zero.
pub const MIN_GEOMETRY_LENGTH: Typed = Typed::new(1e-9, "mm"); // unit-factors: a floor against a zero size, in the geometry frame's millimetres

/// [`MIN_GEOMETRY_LENGTH`] in the geometry frame.
pub fn min_geometry_length() -> f64 {
    geometry_amount(Quantity::Length, MIN_GEOMETRY_LENGTH)
}
