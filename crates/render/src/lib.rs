//! Display data for FieldCraft, with no GPU dependency: what a finite element model looks like
//! on screen ([`DisplayMesh`]: the visible faces, element edges, line elements and nodes), an
//! orbit camera with perspective or orthographic projection, picking ([`pick`]: what a ray through
//! a pixel hits), flow fields on a grid of cells ([`flow`]: slices, surface values and
//! streamlines), and a CPU rasterizer for headless snapshots. The desktop viewport draws the
//! same [`Scene`] on the GPU.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

mod camera;
pub mod caption;
pub mod contour;
mod display;
pub mod flow;
mod font;
pub mod glyphs;
pub mod pick;
mod raster;
mod vec;

pub use camera::{ASPECT_RANGE, Camera, FIT_MARGIN, FULL_VIEW, MIN_SIZE, Mat4, Orientation, StandardView};
pub use display::{DisplayEdge, DisplayLine, DisplayMesh, DisplayNode, DisplayTri, SurfaceCurve, SurfaceTri};
pub use raster::{
    Canvas, LABEL_OFFSET, MAX_DIM, MAX_PIXELS, MAX_PRIMITIVES, RenderError, Rgb, Scene, SceneFill, SceneLabel, SceneLine, ScenePoint, SceneTri,
    draw_labels, fit_size, png, render, render_layers, render_layers_png, render_png, shade,
};
pub use vec::{Aabb3, Vec3};
