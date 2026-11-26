//! Shared Craft-family panel host (`solidcraft_ui_kit::panel_runtime`).
//!
//! Commit 1 establishes the 4-panel registry, seeded default layout, and layout persistence store;
//! Commit 2 wires `PanelHost` into `ModelApp::shell`.

#![allow(dead_code)]

pub mod store;

use std::collections::BTreeSet;

use solidcraft_panel_layout::{
    Column, ColumnId, EdgeRule, HeaderActionSpec, InstanceId, Layout, MAX_LOGICAL, PanelFeatures, PanelInstance, PanelKind, PanelRegistry,
    PanelResult, PanelSpec, PanelSurface, PlacementMemory, Section, SectionId, SectionLocation, SideHeightPolicy, Size, WorkspaceId, Zone,
    normalize_column_widths, validate,
};

use crate::files::Prefs;

pub const WORKSPACE_ID: &str = "fieldcraft.main";
pub const PANEL_BROWSER: &str = "browser";
pub const PANEL_HISTORY: &str = "history";
pub const PANEL_RESULTS: &str = "results";
pub const PANEL_INSPECTOR: &str = "inspector";
pub const ALL_PANEL_IDS: [&str; 4] = [PANEL_BROWSER, PANEL_HISTORY, PANEL_RESULTS, PANEL_INSPECTOR];

pub const KIND_BROWSER: &str = "fieldcraft.browser";
pub const KIND_HISTORY: &str = "fieldcraft.history";
pub const KIND_RESULTS: &str = "fieldcraft.results";
pub const KIND_INSPECTOR: &str = "fieldcraft.inspector";

pub const MIN_PANEL_SIZE: f32 = 180.0;
pub const DEFAULT_BROWSER_HEIGHT: f32 = 520.0;
pub const DEFAULT_HISTORY_WIDTH: f32 = 640.0;
pub const DEFAULT_RESULTS_HEIGHT: f32 = 420.0;
pub const DEFAULT_INSPECTOR_HEIGHT: f32 = 360.0;
pub const RIGHT_INSPECTOR_WEIGHT: f32 = 0.45;
pub const RIGHT_RESULTS_WEIGHT: f32 = 0.55;

pub(crate) fn panel_kind(id: &str) -> &'static str {
    match id {
        PANEL_BROWSER => KIND_BROWSER,
        PANEL_HISTORY => KIND_HISTORY,
        PANEL_RESULTS => KIND_RESULTS,
        PANEL_INSPECTOR => KIND_INSPECTOR,
        _ => "",
    }
}

/// Register FieldCraft's four persistent panels.
pub(crate) fn build_registry() -> PanelResult<PanelRegistry> {
    let mut registry = PanelRegistry::default();
    let side_zones: BTreeSet<Zone> = [Zone::Left, Zone::Right, Zone::Floating, Zone::Window].into_iter().collect();
    let history_zones: BTreeSet<Zone> = [Zone::Bottom, Zone::Floating, Zone::Window].into_iter().collect();

    let browser_kind = PanelKind::new(KIND_BROWSER)?;
    registry.register(PanelSpec {
        kind: browser_kind.clone(),
        title_key: "panel.browser".into(),
        icon: "list".into(),
        allowed_zones: side_zones.clone(),
        default_zone: Zone::Left,
        min_size: Size { width: crate::browser::MIN_WIDTH, height: MIN_PANEL_SIZE },
        default_size: Size { width: crate::browser::WIDTH, height: DEFAULT_BROWSER_HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Translucent,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: false, popout: true, tabbable: true, ..Default::default() },
        side_height_policy: SideHeightPolicy::Fill,
        fixed_bottom_height: None,
        header_actions: Vec::new(),
        text_heavy: false,
    })?;
    registry.set_side_widths(&browser_kind, vec![(crate::browser::MIN_WIDTH, crate::browser::MAX_WIDTH)])?;

    let history_actions = [
        ("first", "ui.history_first", "timeline_first"),
        ("back", "ui.history_back", "timeline_back"),
        ("forward", "ui.history_forward", "timeline_forward"),
        ("last", "ui.history_last", "timeline_last"),
    ]
    .into_iter()
    .map(|(id, label_key, icon)| {
        Ok(HeaderActionSpec { id: solidcraft_panel_layout::ControlId::new(id)?, label_key: label_key.into(), icon: icon.into() })
    })
    .collect::<PanelResult<Vec<_>>>()?;

    registry.register(PanelSpec {
        kind: PanelKind::new(KIND_HISTORY)?,
        title_key: "panel.history".into(),
        icon: "ui.views".into(),
        allowed_zones: history_zones,
        default_zone: Zone::Bottom,
        min_size: Size { width: 240.0, height: crate::timeline::HEIGHT },
        default_size: Size { width: DEFAULT_HISTORY_WIDTH, height: crate::timeline::HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Opaque,
        features: PanelFeatures {
            collapsible: true,
            closable: true,
            hot_closable: false,
            popout: true,
            tabbable: false,
            ..Default::default()
        },
        side_height_policy: SideHeightPolicy::Content,
        fixed_bottom_height: Some(crate::timeline::HEIGHT),
        header_actions: history_actions,
        text_heavy: false,
    })?;

    registry.register(PanelSpec {
        kind: PanelKind::new(KIND_INSPECTOR)?,
        title_key: "panel.inspector".into(),
        icon: "list".into(),
        allowed_zones: side_zones.clone(),
        default_zone: Zone::Right,
        min_size: Size { width: MIN_PANEL_SIZE, height: MIN_PANEL_SIZE },
        default_size: Size { width: crate::context::WIDTH, height: DEFAULT_INSPECTOR_HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Opaque,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: true, popout: true, tabbable: true, ..Default::default() },
        side_height_policy: SideHeightPolicy::Fill,
        fixed_bottom_height: None,
        header_actions: Vec::new(),
        text_heavy: false,
    })?;

    registry.register(PanelSpec {
        kind: PanelKind::new(KIND_RESULTS)?,
        title_key: "panel.results".into(),
        icon: "results.deformed".into(),
        allowed_zones: side_zones,
        default_zone: Zone::Right,
        min_size: Size { width: MIN_PANEL_SIZE, height: MIN_PANEL_SIZE },
        default_size: Size { width: crate::context::WIDTH, height: DEFAULT_RESULTS_HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Opaque,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: true, popout: true, tabbable: true, ..Default::default() },
        side_height_policy: SideHeightPolicy::Fill,
        fixed_bottom_height: None,
        header_actions: Vec::new(),
        text_heavy: false,
    })?;

    Ok(registry)
}

/// Seed the initial layout from legacy `Prefs` (`browser_width`, `browser_collapsed`).
/// `history_dock == HistoryDock::Top` is tracked by `FieldCraftLayoutStore` alongside `Zone::Bottom`.
pub(crate) fn seeded_default_layout(registry: &PanelRegistry, prefs: &Prefs) -> PanelResult<Layout> {
    let workspace = WorkspaceId::new(WORKSPACE_ID)?;
    let mut layout = Layout::empty(workspace);

    let left_col = ColumnId::new("col-left")?;
    let right_col = ColumnId::new("col-right")?;
    let browser_w = crate::browser::width(prefs.browser_width);
    let right_w = crate::context::WIDTH;

    layout.columns.push(Column { id: left_col.clone(), zone: Zone::Left, width: browser_w });
    layout.columns.push(Column { id: right_col.clone(), zone: Zone::Right, width: right_w });

    struct Entry {
        id: &'static str,
        section: &'static str,
        zone: Zone,
        column: Option<ColumnId>,
        size: Size,
        collapsed: bool,
        weight: f32,
    }

    let entries = [
        Entry {
            id: PANEL_BROWSER,
            section: "sec-browser",
            zone: Zone::Left,
            column: Some(left_col),
            size: Size { width: browser_w, height: DEFAULT_BROWSER_HEIGHT },
            collapsed: prefs.browser_collapsed,
            weight: 1.0,
        },
        Entry {
            id: PANEL_INSPECTOR,
            section: "sec-inspector",
            zone: Zone::Right,
            column: Some(right_col.clone()),
            size: Size { width: right_w, height: DEFAULT_INSPECTOR_HEIGHT },
            collapsed: false,
            weight: RIGHT_INSPECTOR_WEIGHT,
        },
        Entry {
            id: PANEL_RESULTS,
            section: "sec-results",
            zone: Zone::Right,
            column: Some(right_col),
            size: Size { width: right_w, height: DEFAULT_RESULTS_HEIGHT },
            collapsed: false,
            weight: RIGHT_RESULTS_WEIGHT,
        },
        Entry {
            id: PANEL_HISTORY,
            section: "sec-history",
            zone: Zone::Bottom,
            column: None,
            size: Size { width: DEFAULT_HISTORY_WIDTH, height: crate::timeline::HEIGHT },
            collapsed: false,
            weight: 1.0,
        },
    ];

    for (section_index, entry) in entries.into_iter().enumerate() {
        let id = InstanceId::new(entry.id)?;
        let kind = PanelKind::new(panel_kind(entry.id))?;
        let section = SectionId::new(entry.section)?;
        let location = SectionLocation::Dock { zone: entry.zone };
        let column_width = entry.column.as_ref().map(|_| entry.size.width);

        layout.sections.push(Section {
            id: section.clone(),
            location: location.clone(),
            tabs: vec![id.clone()],
            active: id.clone(),
            size: entry.size,
            collapsed: entry.collapsed,
            column: entry.column.clone(),
            weight: entry.weight,
            preferred_extent: None,
            extent_priority: 0,
        });
        layout.instances.push(PanelInstance {
            id,
            kind,
            requested_open: true,
            memory: PlacementMemory {
                section,
                location,
                section_index,
                tab_index: 0,
                size: entry.size,
                collapsed: entry.collapsed,
                column: entry.column,
                column_index: 0,
                column_width,
                weight: entry.weight,
                preferred_extent: None,
                extent_priority: 0,
            },
        });
    }

    normalize_column_widths(&mut layout, registry)?;
    validate(&layout, registry)?;
    Ok(layout)
}
