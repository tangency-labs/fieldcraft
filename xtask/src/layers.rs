//! Dependency layering rules (plan/architecture.md §3).
//!
//! Adapted from SolidCraft's xtask (MIT OR Apache-2.0), see ATTRIBUTION.md.
//!
//! The rule engine works on a small, metadata-independent model so it can be
//! unit-tested; `from_metadata` builds that model from `cargo metadata`.

use serde_json::Value;

/// Where a workspace crate sits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class {
    /// Regular layered crate.
    Layer(u8),
    /// Layer 0, and additionally may depend on no workspace crate at all.
    Standalone,
    /// Test tooling: may depend on anything up to L7 (it sits at [`TESTKIT_LAYER`]
    /// for rule purposes); other crates may use it only as a dev-dependency.
    Testkit,
    /// Binaries and build tooling: exempt from the rules.
    Exempt,
}

impl Class {
    fn layer(self) -> Option<u8> {
        match self {
            Class::Layer(l) => Some(l),
            Class::Standalone => Some(0),
            Class::Testkit => Some(TESTKIT_LAYER),
            Class::Exempt => None,
        }
    }
}

/// The layering table. Names are package names without the `fieldcraft-`
/// prefix.
pub const TABLE: &[(&str, Class)] = &[
    // L0: the solver-neutral model.
    ("model", Class::Layer(0)),
    // L0, no workspace dependencies: deterministic thread pool helpers for the solvers.
    ("parallel", Class::Standalone),
    // L0, no workspace dependencies: checking for, verifying and installing releases.
    ("update", Class::Standalone),
    // L1: geometry and math kernels.
    ("geometry", Class::Layer(1)),
    ("mesh", Class::Layer(1)),
    ("fem", Class::Layer(1)),
    ("decks", Class::Layer(1)),
    ("render", Class::Layer(1)),
    // L2: solvers.
    ("solver", Class::Layer(2)),
    ("explicit", Class::Layer(2)),
    ("lbm", Class::Layer(2)),
    // L3: coupling between solvers (`fsi`, later).
    // L4 and up: the document, files, commands, UI.
    ("doc", Class::Layer(4)),
    ("io", Class::Layer(5)),
    ("engine", Class::Layer(6)),
    ("gpu", Class::Layer(7)),
    ("lbm-gpu", Class::Layer(7)),
    ("mesh-gpu", Class::Layer(7)),
    ("ui-egui", Class::Layer(7)),
    ("mcp", Class::Layer(7)),
    ("testkit", Class::Testkit),
    // apps and tooling
    ("fieldcraft", Class::Exempt),
    ("cli", Class::Exempt),
    ("web", Class::Exempt),
    ("xtask", Class::Exempt),
];

/// Explicit orderings *within* a layer (earlier may be used by later).
pub const INTRA_LAYER_ORDER: &[&[&str]] = &[
    // The window draws through the GPU crate's pipelines.
    &["gpu", "ui-egui"],
    // The MCP server's headless screenshots draw with them too.
    &["gpu", "mcp"],
    // The lattice Boltzmann GPU backend opens its device with the GPU crate's adapter policy.
    &["gpu", "lbm-gpu"],
    // So does the projection smoothing GPU device.
    &["gpu", "mesh-gpu"],
];

fn intra_layer_allowed(from: &str, to: &str) -> bool {
    let (from, to) = (short_name(from), short_name(to));
    INTRA_LAYER_ORDER.iter().any(|chain| match (chain.iter().position(|n| *n == from), chain.iter().position(|n| *n == to)) {
        (Some(f), Some(t)) => t < f,
        _ => false,
    })
}

/// External crates that constitute a UI toolkit / windowing dependency.
/// Entries ending in `*` are prefixes.
pub const UI_CRATES: &[&str] = &["egui", "eframe", "egui-wgpu", "wgpu", "winit", "egui_kittest", "rfd", "bevy*"];

/// The only crate that may depend on SolidCraft (`solidcraft-*`, and the solidcore kernel crates
/// behind it, `solidcore*`): the geometry boundary.
pub const SOLIDCRAFT_GATEWAY: &str = "geometry";

/// SolidCraft and its kernel: the crates only [`SOLIDCRAFT_GATEWAY`] may use.
fn is_solidcraft_crate(name: &str) -> bool {
    name.starts_with("solidcraft-") || name == "solidcore" || name.starts_with("solidcore-")
}

/// The one exception to [`SOLIDCRAFT_GATEWAY`] (owner decision, 2026-11-10): the desktop UI draws
/// with SolidCraft's UI kit and the token contract it draws from, and with nothing else of
/// SolidCraft's. Short crate name, then the SolidCraft crates it may use.
pub const SOLIDCRAFT_UI_KIT: (&str, &[&str]) = ("ui-egui", &["solidcraft-ui-kit", "solidcraft-theme"]);

/// Whether `krate` (short name) may depend on the SolidCraft crate `dep`.
fn may_use_solidcraft(krate: &str, dep: &str) -> bool {
    krate == SOLIDCRAFT_GATEWAY || (krate == SOLIDCRAFT_UI_KIT.0 && SOLIDCRAFT_UI_KIT.1.contains(&dep))
}

/// First layer allowed to use UI crates.
pub const UI_MIN_LAYER: u8 = 7;

/// The layer the testkit sits at for rule purposes: above every layered crate.
pub const TESTKIT_LAYER: u8 = UI_MIN_LAYER + 1;

/// The layer reported for apps and tooling when a layered crate depends on one.
const APPS_LAYER: u8 = TESTKIT_LAYER + 1;

pub fn short_name(pkg: &str) -> &str {
    pkg.strip_prefix("fieldcraft-").unwrap_or(pkg)
}

pub fn classify(pkg: &str) -> Option<Class> {
    let s = short_name(pkg);
    TABLE.iter().find(|(n, _)| *n == s).map(|(_, c)| *c)
}

fn is_ui_crate(name: &str) -> bool {
    UI_CRATES.iter().any(|p| match p.strip_suffix('*') {
        Some(prefix) => name.starts_with(prefix),
        None => name == *p,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepKind {
    Normal,
    Dev,
    Build,
}

#[derive(Debug, Clone)]
pub struct Dep {
    pub name: String,
    pub kind: DepKind,
    /// `true` if the dependency is a workspace member.
    pub workspace: bool,
}

#[derive(Debug, Clone)]
pub struct Crate {
    pub name: String,
    pub deps: Vec<Dep>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Violation {
    Unregistered { krate: String },
    Upward { krate: String, dep: String, from: u8, to: u8, kind: DepKind },
    StandaloneHasWorkspaceDep { krate: String, dep: String },
    TestkitAsNormalDep { krate: String },
    UiBelowUiLayer { krate: String, dep: String, layer: u8 },
    SolidCraftOutsideGeometry { krate: String, dep: String },
}

impl std::fmt::Display for Violation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Violation::Unregistered { krate } => {
                write!(f, "{krate}: unknown workspace crate; register it in xtask/src/layers.rs TABLE (see plan/architecture.md §3)")
            }
            Violation::Upward { krate, dep, from, to, kind } => {
                write!(f, "{krate} (L{from}) -> {dep} (L{to}) [{kind:?}]: may only depend on strictly lower layers")
            }
            Violation::StandaloneHasWorkspaceDep { krate, dep } => {
                write!(f, "{krate}: standalone crate must not depend on workspace crate {dep}")
            }
            Violation::TestkitAsNormalDep { krate } => {
                write!(f, "{krate}: fieldcraft-testkit may only be a dev-dependency")
            }
            Violation::UiBelowUiLayer { krate, dep, layer } => {
                write!(f, "{krate} (L{layer}) depends on UI crate `{dep}`; UI toolkits are only allowed in L{UI_MIN_LAYER}+")
            }
            Violation::SolidCraftOutsideGeometry { krate, dep } => {
                write!(
                    f,
                    "{krate} depends on `{dep}`; only fieldcraft-{SOLIDCRAFT_GATEWAY} may use SolidCraft (and fieldcraft-{} its UI kit: {})",
                    SOLIDCRAFT_UI_KIT.0,
                    SOLIDCRAFT_UI_KIT.1.join(", ")
                )
            }
        }
    }
}

/// Check all rules. Returns violations sorted for stable output.
pub fn check(crates: &[Crate]) -> Vec<Violation> {
    let mut out = Vec::new();
    for c in crates {
        let Some(class) = classify(&c.name) else {
            out.push(Violation::Unregistered { krate: c.name.clone() });
            continue;
        };
        if class == Class::Exempt {
            continue;
        }
        let layer = class.layer().unwrap_or(0);
        for d in &c.deps {
            // Self dev-dependencies (e.g. to enable features in tests) are fine.
            if d.name == c.name {
                continue;
            }
            if d.workspace {
                if class == Class::Standalone {
                    out.push(Violation::StandaloneHasWorkspaceDep { krate: c.name.clone(), dep: d.name.clone() });
                    continue;
                }
                match classify(&d.name) {
                    // Unregistered deps are reported on their own entry.
                    None => {}
                    Some(Class::Testkit) if d.kind != DepKind::Dev => {
                        out.push(Violation::TestkitAsNormalDep { krate: c.name.clone() });
                    }
                    Some(Class::Testkit) => {}
                    Some(dc) => {
                        let to = dc.layer().unwrap_or(u8::MAX);
                        if to >= layer && !(to == layer && intra_layer_allowed(&c.name, &d.name)) {
                            out.push(Violation::Upward {
                                krate: c.name.clone(),
                                dep: d.name.clone(),
                                from: layer,
                                to: to.min(APPS_LAYER),
                                kind: d.kind,
                            });
                        }
                    }
                }
            } else if is_solidcraft_crate(&d.name) && !may_use_solidcraft(short_name(&c.name), &d.name) {
                out.push(Violation::SolidCraftOutsideGeometry { krate: c.name.clone(), dep: d.name.clone() });
            } else if layer < UI_MIN_LAYER && is_ui_crate(&d.name) {
                out.push(Violation::UiBelowUiLayer { krate: c.name.clone(), dep: d.name.clone(), layer });
            }
        }
    }
    out.sort_by_key(|v| v.to_string());
    out.dedup();
    out
}

/// Build the model from `cargo metadata --format-version 1 --no-deps`.
pub fn from_metadata(meta: &Value) -> Result<Vec<Crate>, String> {
    let pkgs = meta["packages"].as_array().ok_or("metadata: no packages array")?;
    let members: Vec<&str> = pkgs.iter().filter_map(|p| p["name"].as_str()).collect();
    let mut out = Vec::new();
    for p in pkgs {
        let name = p["name"].as_str().ok_or("package without name")?.to_owned();
        let mut deps = Vec::new();
        for d in p["dependencies"].as_array().into_iter().flatten() {
            let dname = d["name"].as_str().unwrap_or_default().to_owned();
            let kind = match d["kind"].as_str() {
                Some("dev") => DepKind::Dev,
                Some("build") => DepKind::Build,
                _ => DepKind::Normal,
            };
            let workspace = members.contains(&dname.as_str()) || d["path"].is_string();
            deps.push(Dep { name: dname, kind, workspace });
        }
        out.push(Crate { name, deps });
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(out)
}

pub fn describe(class: Option<Class>) -> String {
    match class {
        Some(Class::Layer(l)) => format!("L{l}"),
        Some(Class::Standalone) => "L0 standalone".into(),
        Some(Class::Testkit) => "testkit".into(),
        Some(Class::Exempt) => "exempt".into(),
        None => "UNREGISTERED".into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(name: &str, deps: &[(&str, DepKind, bool)]) -> Crate {
        Crate { name: name.into(), deps: deps.iter().map(|(n, k, w)| Dep { name: (*n).into(), kind: *k, workspace: *w }).collect() }
    }
    use DepKind::*;

    #[test]
    fn clean_downward_graph_passes() {
        let g = [
            c("fieldcraft-model", &[("serde", Normal, false)]),
            c("fieldcraft-doc", &[("fieldcraft-model", Normal, true)]),
            c("fieldcraft-engine", &[("fieldcraft-doc", Normal, true), ("fieldcraft-testkit", Dev, true)]),
            c("fieldcraft-ui-egui", &[("fieldcraft-engine", Normal, true), ("egui", Normal, false)]),
            c("fieldcraft-cli", &[("fieldcraft-ui-egui", Normal, true)]),
        ];
        assert!(check(&g).is_empty(), "{:?}", check(&g));
    }

    #[test]
    fn upward_dependency_flagged() {
        let v = check(&[c("fieldcraft-doc", &[("fieldcraft-engine", Normal, true)])]);
        assert!(matches!(v[..], [Violation::Upward { from: 4, to: 6, .. }]));
        let v = check(&[c("fieldcraft-fem", &[("fieldcraft-solver", Normal, true)])]);
        assert!(matches!(v[..], [Violation::Upward { from: 1, to: 2, .. }]));
    }

    #[test]
    fn sideways_dependency_flagged() {
        let v = check(&[c("fieldcraft-mesh", &[("fieldcraft-fem", Normal, true)])]);
        assert!(matches!(v[..], [Violation::Upward { from: 1, to: 1, .. }]));
    }

    #[test]
    fn self_dev_dependency_ignored() {
        assert!(check(&[c("fieldcraft-doc", &[("fieldcraft-doc", Dev, true)])]).is_empty());
    }

    #[test]
    fn upward_dev_dependency_flagged() {
        let v = check(&[c("fieldcraft-model", &[("fieldcraft-doc", Dev, true)])]);
        assert!(matches!(v[..], [Violation::Upward { kind: Dev, .. }]));
    }

    #[test]
    fn solvers_sit_above_the_kernels() {
        let g = [
            c("fieldcraft-fem", &[("fieldcraft-model", Normal, true)]),
            c("fieldcraft-solver", &[("fieldcraft-fem", Normal, true), ("fieldcraft-mesh", Normal, true)]),
            c("fieldcraft-engine", &[("fieldcraft-solver", Normal, true)]),
        ];
        assert!(check(&g).is_empty(), "{:?}", check(&g));
        let v = check(&[c("fieldcraft-solver", &[("fieldcraft-engine", Normal, true)])]);
        assert!(matches!(v[..], [Violation::Upward { from: 2, to: 6, .. }]));
        let v = check(&[c("fieldcraft-engine", &[("fieldcraft-cli", Normal, true)])]);
        assert!(matches!(v[..], [Violation::Upward { from: 6, to: APPS_LAYER, .. }]));
    }

    #[test]
    fn ui_crates_below_the_ui_layer_flagged() {
        for dep in ["egui", "eframe", "winit", "egui_kittest", "rfd", "bevy_ecs", "bevy"] {
            let v = check(&[c("fieldcraft-engine", &[(dep, Normal, false)])]);
            assert!(matches!(v[..], [Violation::UiBelowUiLayer { layer: 6, .. }]), "{dep}");
            let v = check(&[c("fieldcraft-solver", &[(dep, Normal, false)])]);
            assert!(matches!(v[..], [Violation::UiBelowUiLayer { layer: 2, .. }]), "{dep}");
        }
        assert!(check(&[c("fieldcraft-engine", &[("egui_extras_not", Normal, false)])]).is_empty());
        assert!(check(&[c("fieldcraft-ui-egui", &[("winit", Normal, false)])]).is_empty());
    }

    #[test]
    fn only_geometry_uses_solidcraft() {
        assert!(check(&[c("fieldcraft-geometry", &[("solidcraft-kernel", Normal, false)])]).is_empty());
        let v = check(&[c("fieldcraft-mesh", &[("solidcraft-geom", Normal, false)])]);
        assert!(matches!(v[..], [Violation::SolidCraftOutsideGeometry { .. }]));
        // The kernel behind it too.
        assert!(check(&[c("fieldcraft-geometry", &[("solidcore-record", Normal, false)])]).is_empty());
        for kernel in ["solidcore", "solidcore-record"] {
            let v = check(&[c("fieldcraft-mesh", &[(kernel, Normal, false)])]);
            assert!(matches!(v[..], [Violation::SolidCraftOutsideGeometry { .. }]), "{kernel}");
        }
    }

    #[test]
    fn only_the_ui_uses_the_ui_kit_and_nothing_else_of_solidcraft() {
        for kit in ["solidcraft-ui-kit", "solidcraft-theme"] {
            assert!(check(&[c("fieldcraft-ui-egui", &[(kit, Normal, false)])]).is_empty(), "{kit}");
            let v = check(&[c("fieldcraft-engine", &[(kit, Normal, false)])]);
            assert!(matches!(v[..], [Violation::SolidCraftOutsideGeometry { .. }]), "{kit}");
        }
        for other in ["solidcraft-kernel", "solidcraft-ui-egui", "solidcore-record"] {
            let v = check(&[c("fieldcraft-ui-egui", &[(other, Normal, false)])]);
            assert!(matches!(v[..], [Violation::SolidCraftOutsideGeometry { .. }]), "{other}");
        }
    }

    #[test]
    fn unregistered_crate_is_error() {
        let v = check(&[c("fieldcraft-mystery", &[])]);
        assert!(matches!(&v[..], [Violation::Unregistered { krate }] if krate == "fieldcraft-mystery"));
        assert!(v[0].to_string().contains("register"));
    }

    #[test]
    fn testkit_only_as_dev_dependency() {
        let v = check(&[c("fieldcraft-doc", &[("fieldcraft-testkit", Normal, true)])]);
        assert!(matches!(v[..], [Violation::TestkitAsNormalDep { .. }]));
        assert!(check(&[c("fieldcraft-doc", &[("fieldcraft-testkit", Dev, true)])]).is_empty());
        assert!(check(&[c("fieldcraft-testkit", &[("fieldcraft-engine", Normal, true)])]).is_empty());
    }

    #[test]
    fn apps_and_xtask_exempt() {
        for app in ["fieldcraft", "fieldcraft-cli", "xtask"] {
            assert!(check(&[c(app, &[("egui", Normal, false), ("fieldcraft-ui-egui", Normal, true)])]).is_empty());
        }
    }

    #[test]
    fn metadata_parsing() {
        let meta: Value = serde_json::from_str(
            r#"{"packages":[
                {"name":"fieldcraft-doc","dependencies":[
                    {"name":"fieldcraft-model","kind":null,"path":"/x/crates/model"},
                    {"name":"serde","kind":null},
                    {"name":"proptest","kind":"dev"}]},
                {"name":"fieldcraft-model","dependencies":[]}
            ]}"#,
        )
        .unwrap();
        let g = from_metadata(&meta).unwrap();
        assert_eq!(g.len(), 2);
        let doc = g.iter().find(|c| c.name == "fieldcraft-doc").unwrap();
        assert!(doc.deps[0].workspace && !doc.deps[1].workspace);
        assert_eq!(doc.deps[2].kind, Dev);
        assert!(check(&g).is_empty());
    }
}
