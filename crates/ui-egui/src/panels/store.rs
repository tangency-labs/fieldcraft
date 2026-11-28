//! Panel layout persistence, one-time legacy migration, and hostile-payload validation.

use std::collections::BTreeSet;

use solidcraft_panel_layout::{
    ActionSink, Layout, LayoutAction, LayoutEnvelope, LayoutReply, MAX_ENVELOPE_BYTES, MAX_LOGICAL, PanelRegistry, PanelResult, Point, Rule,
    SectionLocation, Size, WorkspaceId, decode, encode, error, preparation_context, validate, with_live_layout,
};
use solidcraft_ui_kit::panel_runtime::{PanelLayoutStore, standalone_header_action};

use super::{ALL_PANEL_IDS, panel_kind, seeded_default_layout};
use crate::files::{MAX_PREFS_BYTES, Prefs, RawJson};

/// Clamp origin coordinates when projecting a saved `Zone::Window` panel to an in-window
/// `Floating` card on a `native_windows == false` session.
const PROJECTED_ORIGIN_MAX: f32 = 400.0;

/// Result of encoding the live layout for `Prefs::panel_layout`.
pub(crate) enum EncodedLayout {
    Ready(RawJson),
    ExceedsBudget,
    Error,
}

/// Persistence bridge between `WorkspacePanelHost` and `Prefs::panel_layout`.
pub(crate) struct FieldCraftLayoutStore {
    envelope: LayoutEnvelope,
    encoded_bytes: Vec<u8>,
    dirty: bool,
    has_projected_native: bool,
    preserved_native_envelope: Option<LayoutEnvelope>,
}

impl FieldCraftLayoutStore {
    /// Load and validate `prefs.panel_layout`, or seed once from legacy `Prefs` fields.
    ///
    /// When `native_windows == false` and a valid saved layout contains `SectionLocation::Window`
    /// placements, those placements are projected non-destructively to clamped in-window
    /// `Floating` coordinates for the live runtime while retaining the original `Window` envelope
    /// until the user commits an explicit layout change.
    pub(crate) fn load_from_prefs(prefs: &Prefs, registry: &PanelRegistry, native_windows: bool) -> PanelResult<(Layout, Self)> {
        let fallback = seeded_default_layout(registry, prefs)?;
        let fallback_envelope = LayoutEnvelope::new(fallback.clone());

        let Some(raw_value) = &prefs.panel_layout else {
            let encoded_bytes = encode(&fallback_envelope, registry)?;
            return Ok((
                fallback,
                Self { envelope: fallback_envelope, encoded_bytes, dirty: false, has_projected_native: false, preserved_native_envelope: None },
            ));
        };

        match try_decode_validated(raw_value, registry, &fallback) {
            Ok(envelope) => {
                let mut live = envelope.layout.clone();
                let mut has_projected_native = false;
                let mut preserved_native_envelope = None;
                if !native_windows && has_window_location(&live) {
                    preserved_native_envelope = Some(envelope.clone());
                    project_window_to_floating(&mut live, registry)?;
                    has_projected_native = true;
                }
                let live_envelope = if has_projected_native { LayoutEnvelope::new(live.clone()) } else { envelope };
                let encoded_bytes = encode(&live_envelope, registry)?;
                Ok((live, Self { envelope: live_envelope, encoded_bytes, dirty: false, has_projected_native, preserved_native_envelope }))
            }
            Err(reason) => {
                log::warn!("panel_layout invalid; resetting to defaults: {reason}");
                let encoded_bytes = encode(&fallback_envelope, registry)?;
                Ok((
                    fallback,
                    Self { envelope: fallback_envelope, encoded_bytes, dirty: false, has_projected_native: false, preserved_native_envelope: None },
                ))
            }
        }
    }

    pub(crate) fn is_dirty(&self) -> bool {
        self.dirty
    }

    pub(crate) fn mark_dirty(&mut self) {
        self.dirty = true;
        self.preserved_native_envelope = None;
        self.has_projected_native = false;
    }

    pub(crate) fn clear_dirty(&mut self) {
        self.dirty = false;
    }

    /// Return and clear any preserved `Zone::Window` envelope when upgrading an unedited
    /// session from `native_windows == false` to desktop native capability.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn take_preserved_native(&mut self) -> Option<LayoutEnvelope> {
        if self.dirty || !self.has_projected_native {
            return None;
        }
        self.has_projected_native = false;
        self.preserved_native_envelope.take()
    }

    /// Return section reservations from the current layout envelope.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn preparation_reservations(
        &self,
        layout: &Layout,
        registry: &PanelRegistry,
    ) -> PanelResult<BTreeSet<solidcraft_panel_layout::SectionId>> {
        let envelope = with_live_layout(&self.envelope, layout.clone(), registry)?;
        Ok(preparation_context(&envelope)?.sections)
    }

    /// Snapshot the pre-close `Zone::Window` layout into `preserved_native_envelope`
    /// while leaving `dirty` untouched so OS shutdown saves cannot overwrite saved placements.
    /// Retains the first preserved snapshot while one exists so sequential native closes keep all
    /// windows in `Zone::Window`.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn preserve_before_deferred_close(&mut self, pre_close_layout: &Layout, registry: &PanelRegistry) {
        if self.preserved_native_envelope.is_some() {
            return;
        }
        if let Ok(envelope) = with_live_layout(&self.envelope, pre_close_layout.clone(), registry) {
            self.preserved_native_envelope = Some(envelope);
            self.has_projected_native = true;
        }
    }

    /// Replace the cached envelope with `layout` (for example after `reset_layout`).
    pub(crate) fn replace_live(&mut self, layout: &Layout, registry: &PanelRegistry, dirty: bool) -> PanelResult<()> {
        let next = LayoutEnvelope::new(layout.clone());
        let bytes = encode(&next, registry)?;
        self.envelope = next;
        self.encoded_bytes = bytes;
        if dirty {
            self.mark_dirty();
        }
        Ok(())
    }

    /// Project any live `Zone::Window` placements to in-window `Floating` non-destructively
    /// while preserving the original `Zone::Window` envelope until an explicit user edit.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn project_live_to_floating(&mut self, runtime_layout: &Layout, registry: &PanelRegistry) -> PanelResult<Layout> {
        let mut live = runtime_layout.clone();
        if !self.dirty && !self.has_projected_native && has_window_location(&live) {
            let envelope = with_live_layout(&self.envelope, live.clone(), registry)?;
            self.preserved_native_envelope = Some(envelope);
            self.has_projected_native = true;
        }
        if has_window_location(&live) {
            project_window_to_floating(&mut live, registry)?;
        }
        let next_env = LayoutEnvelope::new(live.clone());
        self.encoded_bytes = encode(&next_env, registry)?;
        self.envelope = next_env;
        Ok(live)
    }

    /// Serialize the current layout for `Prefs::panel_layout`.
    ///
    /// If `native_windows == false` projected saved `Zone::Window` placements into `Floating`
    /// and the user has not explicitly moved any panel yet (`!self.dirty`), this writes back the
    /// preserved original `Zone::Window` envelope so unrelated preference saves never discard
    /// desktop multi-monitor coordinates.
    pub(crate) fn encode_prefs_layout(&self, runtime_layout: &Layout, registry: &PanelRegistry) -> EncodedLayout {
        let envelope = if !self.dirty
            && self.has_projected_native
            && let Some(preserved) = &self.preserved_native_envelope
        {
            preserved.clone()
        } else {
            with_live_layout(&self.envelope, runtime_layout.clone(), registry).unwrap_or_else(|_| LayoutEnvelope::new(runtime_layout.clone()))
        };
        let Ok(bytes) = encode(&envelope, registry) else {
            return EncodedLayout::Error;
        };
        if bytes.len() > MAX_PREFS_BYTES as usize {
            log::warn!("panel_layout exceeds {} bytes; omitting from prefs", MAX_PREFS_BYTES);
            return EncodedLayout::ExceedsBudget;
        }
        let Ok(text) = String::from_utf8(bytes) else {
            return EncodedLayout::Error;
        };
        match RawJson::from_json_string(text) {
            Ok(raw) => EncodedLayout::Ready(raw),
            Err(_) => EncodedLayout::Error,
        }
    }

    #[cfg(test)]
    pub(crate) fn to_prefs_value(&self, runtime_layout: &Layout, registry: &PanelRegistry) -> Option<RawJson> {
        match self.encode_prefs_layout(runtime_layout, registry) {
            EncodedLayout::Ready(raw) => Some(raw),
            EncodedLayout::ExceedsBudget | EncodedLayout::Error => None,
        }
    }
}

impl PanelLayoutStore for FieldCraftLayoutStore {
    fn remembered_key(&self, workspace: &WorkspaceId) -> Option<&[u8]> {
        (self.envelope.layout.workspace == *workspace).then_some(self.encoded_bytes.as_slice())
    }

    fn current_envelope(&self, state: &LayoutReply, registry: &PanelRegistry) -> PanelResult<LayoutEnvelope> {
        with_live_layout(&self.envelope, state.layout.clone(), registry)
    }

    fn execute_action(&mut self, action: LayoutAction, sink: &mut dyn ActionSink, registry: &PanelRegistry) -> PanelResult<LayoutReply> {
        if standalone_header_action(&action) {
            return sink.dispatch(action);
        }
        let before_rev = sink.inspect().revision;
        let state = sink.inspect();
        let envelope = self.current_envelope(&state, registry)?;
        let preparation = preparation_context(&envelope)?;
        let reply = sink.dispatch_in_context(action, &preparation)?;
        if let Ok(updated) = with_live_layout(&self.envelope, reply.layout.clone(), registry)
            && let Ok(bytes) = encode(&updated, registry)
        {
            self.envelope = updated;
            self.encoded_bytes = bytes;
            if reply.revision != before_rev {
                self.mark_dirty();
            }
        }
        Ok(reply)
    }
}

fn try_decode_validated(raw_value: &RawJson, registry: &PanelRegistry, fallback: &Layout) -> PanelResult<LayoutEnvelope> {
    let bytes = raw_value.as_bytes();
    let cap = (MAX_PREFS_BYTES as usize).min(MAX_ENVELOPE_BYTES);
    if bytes.len() > cap {
        return Err(error(Rule::CollectionLimit, "panel_layout exceeds byte cap"));
    }
    let outcome = decode(bytes, registry, fallback)?;
    if let Some(notice) = outcome.notice {
        return Err(notice);
    }
    validate_loaded_envelope(&outcome.envelope, registry)?;
    Ok(outcome.envelope)
}

fn validate_loaded_envelope(envelope: &LayoutEnvelope, registry: &PanelRegistry) -> PanelResult<()> {
    if !envelope.quarantine.is_empty() {
        return Err(error(Rule::MissingTarget, "panel_layout contains unknown panel kinds"));
    }
    validate(&envelope.layout, registry)?;

    let expected: BTreeSet<&str> = ALL_PANEL_IDS.into_iter().collect();
    let found: BTreeSet<&str> = envelope.layout.instances.iter().map(|i| i.id.0.as_str()).collect();
    if found != expected || envelope.layout.instances.len() != ALL_PANEL_IDS.len() {
        return Err(error(Rule::ExactlyOnce, "panel_layout panel id set mismatch"));
    }

    for col in &envelope.layout.columns {
        if !col.width.is_finite() || col.width <= 0.0 || col.width > MAX_LOGICAL {
            return Err(error(Rule::Bounds, "non-finite or out-of-range column width"));
        }
    }
    for sec in &envelope.layout.sections {
        validate_size(sec.size)?;
        validate_weight(sec.weight)?;
        validate_opt_extent(sec.preferred_extent)?;
        validate_location(&sec.location)?;
    }
    for inst in &envelope.layout.instances {
        if inst.kind.0 != panel_kind(&inst.id.0) {
            return Err(error(Rule::Schema, "panel_layout panel kind mismatch"));
        }
        validate_size(inst.memory.size)?;
        validate_weight(inst.memory.weight)?;
        validate_opt_extent(inst.memory.preferred_extent)?;
        if let Some(w) = inst.memory.column_width
            && (!w.is_finite() || w <= 0.0 || w > MAX_LOGICAL)
        {
            return Err(error(Rule::Bounds, "non-finite or out-of-range memory column width"));
        }
        validate_location(&inst.memory.location)?;
    }
    Ok(())
}

fn validate_size(size: Size) -> PanelResult<()> {
    if !size.width.is_finite()
        || !size.height.is_finite()
        || size.width <= 0.0
        || size.height <= 0.0
        || size.width > MAX_LOGICAL
        || size.height > MAX_LOGICAL
    {
        return Err(error(Rule::Bounds, "non-finite or out-of-range panel size"));
    }
    Ok(())
}

fn validate_weight(weight: f32) -> PanelResult<()> {
    if !weight.is_finite() || !(0.0..=MAX_LOGICAL).contains(&weight) {
        return Err(error(Rule::Bounds, "non-finite or out-of-range panel weight"));
    }
    Ok(())
}

fn validate_opt_extent(extent: Option<f32>) -> PanelResult<()> {
    if let Some(v) = extent
        && (!v.is_finite() || v <= 0.0 || v > MAX_LOGICAL)
    {
        return Err(error(Rule::Bounds, "non-finite or out-of-range preferred extent"));
    }
    Ok(())
}

fn validate_location(loc: &SectionLocation) -> PanelResult<()> {
    match loc {
        SectionLocation::Dock { .. } => Ok(()),
        SectionLocation::Floating { origin } => validate_point(*origin),
        SectionLocation::Window { origin, return_to, .. } => {
            validate_point(*origin)?;
            validate_size(return_to.size)?;
            validate_weight(return_to.weight)?;
            validate_opt_extent(return_to.preferred_extent)?;
            if let Some(w) = return_to.column_width
                && (!w.is_finite() || w <= 0.0 || w > MAX_LOGICAL)
            {
                return Err(error(Rule::Bounds, "non-finite or out-of-range return column width"));
            }
            Ok(())
        }
    }
}

fn validate_point(pt: Point) -> PanelResult<()> {
    if !pt.x.is_finite() || !pt.y.is_finite() || pt.x.abs() > MAX_LOGICAL || pt.y.abs() > MAX_LOGICAL {
        return Err(error(Rule::Bounds, "non-finite or out-of-range point"));
    }
    Ok(())
}

fn has_window_location(layout: &Layout) -> bool {
    layout.sections.iter().any(|s| matches!(s.location, SectionLocation::Window { .. }))
        || layout.instances.iter().any(|i| matches!(i.memory.location, SectionLocation::Window { .. }))
}

fn clamp_floating_origin(origin: Point) -> Point {
    Point {
        x: if origin.x.is_finite() { origin.x.clamp(0.0, PROJECTED_ORIGIN_MAX) } else { 0.0 },
        y: if origin.y.is_finite() { origin.y.clamp(0.0, PROJECTED_ORIGIN_MAX) } else { 0.0 },
    }
}

fn project_window_to_floating(layout: &mut Layout, registry: &PanelRegistry) -> PanelResult<()> {
    for section in &mut layout.sections {
        if let SectionLocation::Window { origin, .. } = section.location {
            section.location = SectionLocation::Floating { origin: clamp_floating_origin(origin) };
            section.column = None;
        }
    }
    for instance in &mut layout.instances {
        if let SectionLocation::Window { origin, .. } = instance.memory.location {
            instance.memory.location = SectionLocation::Floating { origin: clamp_floating_origin(origin) };
            instance.memory.column = None;
            instance.memory.column_index = 0;
            instance.memory.column_width = None;
        }
    }
    validate(layout, registry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::panels::{PANEL_BROWSER, PANEL_HISTORY, PANEL_INSPECTOR, PANEL_RESULTS, build_registry};
    use crate::shell::HistoryDock;
    use serde_json::json;
    use solidcraft_panel_layout::{
        ColumnId, Destination, DockReturn, HostCapabilities, HostId, InstanceId, Mutation, SectionId, Zone, next_section_id,
    };

    #[test]
    fn fresh_load_seeds_default_zones_and_split_weights() {
        let reg = build_registry().unwrap();
        let prefs = Prefs::default();
        let (layout, store) = FieldCraftLayoutStore::load_from_prefs(&prefs, &reg, false).unwrap();
        assert!(!store.is_dirty());
        assert_eq!(layout.instances.len(), 4);
        assert_eq!(layout.containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap().location.zone(), Zone::Left);
        assert_eq!(layout.containing(&InstanceId::new(PANEL_HISTORY).unwrap()).unwrap().location.zone(), Zone::Bottom);
        let insp = layout.containing(&InstanceId::new(PANEL_INSPECTOR).unwrap()).unwrap();
        let res = layout.containing(&InstanceId::new(PANEL_RESULTS).unwrap()).unwrap();
        assert_eq!(insp.location.zone(), Zone::Right);
        assert_eq!(res.location.zone(), Zone::Right);
        assert_eq!(insp.column, res.column);
        assert!((insp.weight - 0.45).abs() < 1e-6);
        assert!((res.weight - 0.55).abs() < 1e-6);
    }

    #[test]
    fn legacy_prefs_seed_browser_width_and_collapsed_once() {
        let reg = build_registry().unwrap();
        let prefs = Prefs { browser_width: 280.0, browser_collapsed: true, history_dock: HistoryDock::Top, browser_float: true, ..Prefs::default() };
        let (layout, _store) = FieldCraftLayoutStore::load_from_prefs(&prefs, &reg, false).unwrap();
        let browser = layout.containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap();
        assert!(browser.collapsed);
        assert!((browser.size.width - 280.0).abs() < 1e-6);
    }

    #[test]
    fn hostile_and_overflow_panel_layouts_recover_to_defaults() {
        let reg = build_registry().unwrap();
        let non_default_prefs = Prefs { browser_width: 340.0, browser_collapsed: true, ..Prefs::default() };
        let fallback_prefs = Prefs { browser_width: 290.0, ..Prefs::default() };
        let expected_fallback = seeded_default_layout(&reg, &fallback_prefs).unwrap();

        let (valid_layout, _) = FieldCraftLayoutStore::load_from_prefs(&non_default_prefs, &reg, false).unwrap();
        let valid_value = serde_json::to_value(LayoutEnvelope::new(valid_layout)).unwrap();
        let valid_raw = RawJson::from_json_string(serde_json::to_string(&valid_value).unwrap()).unwrap();

        // Sanity check: unmodified non-default placement loads as-is instead of falling back.
        let (loaded_valid, valid_store) =
            FieldCraftLayoutStore::load_from_prefs(&Prefs { panel_layout: Some(valid_raw), ..fallback_prefs.clone() }, &reg, false).unwrap();
        assert!(!valid_store.is_dirty());
        let valid_browser = loaded_valid.containing(&InstanceId::new(PANEL_BROWSER).unwrap()).unwrap();
        assert!((valid_browser.size.width - 340.0).abs() < 1e-6);
        assert!(valid_browser.collapsed);

        // 1. Float overflow (1e40 -> f32::INFINITY) on a section size.
        let mut overflow = valid_value.clone();
        overflow["layout"]["sections"][0]["size"]["width"] = json!(1e40);

        // 2. Raw JSON overflow (1e400 -> outside f64 range).
        let mut sentinel_overflow = valid_value.clone();
        sentinel_overflow["layout"]["sections"][0]["size"]["width"] = json!(98765.0);
        let overflow_1e400 = serde_json::to_string(&sentinel_overflow).unwrap().replacen("98765.0", "1e400", 1);
        assert!(overflow_1e400.contains("1e400"));

        // 3. Unknown panel id added.
        let mut unknown_id = valid_value.clone();
        unknown_id["layout"]["instances"][0]["id"] = json!("mystery");

        // 4. Swapped kind between `results` and `inspector`.
        let mut swapped_kinds = valid_value.clone();
        let instances = swapped_kinds["layout"]["instances"].as_array_mut().unwrap();
        let res_idx = instances.iter().position(|i| i["id"] == PANEL_RESULTS).unwrap();
        let insp_idx = instances.iter().position(|i| i["id"] == PANEL_INSPECTOR).unwrap();
        let res_kind = instances[res_idx]["kind"].clone();
        instances[res_idx]["kind"] = instances[insp_idx]["kind"].clone();
        instances[insp_idx]["kind"] = res_kind;

        // 5. Negative size.
        let mut negative = valid_value.clone();
        negative["layout"]["sections"][0]["size"]["height"] = json!(-50.0);

        // 6. Oversized payload (> 64 KiB).
        let mut oversized = valid_value;
        oversized["recovery_extra"] = json!("a".repeat(70_000));

        let mut bad_payloads: Vec<RawJson> = [json!("corrupt"), json!([1, 2, 3]), overflow, unknown_id, swapped_kinds, negative, oversized]
            .into_iter()
            .map(|v| RawJson::from_json_string(serde_json::to_string(&v).unwrap()).unwrap())
            .collect();
        bad_payloads.push(RawJson::from_json_string(overflow_1e400).unwrap());

        for bad in bad_payloads {
            let p = Prefs { panel_layout: Some(bad), ..fallback_prefs.clone() };
            let (recovered, store) = FieldCraftLayoutStore::load_from_prefs(&p, &reg, false).unwrap();
            assert!(!store.is_dirty());
            assert_eq!(recovered, expected_fallback);
        }
    }

    #[test]
    fn persisted_native_window_projects_non_destructively_when_native_unavailable() {
        let reg = build_registry().unwrap();
        let (mut layout, _) = FieldCraftLayoutStore::load_from_prefs(&Prefs::default(), &reg, true).unwrap();
        let browser_id = InstanceId::new(PANEL_BROWSER).unwrap();
        let win_loc = SectionLocation::Window {
            host: HostId::new("win-browser").unwrap(),
            origin: Point { x: 1800.0, y: 120.0 },
            return_to: DockReturn {
                zone: Zone::Left,
                index: 0,
                size: Size { width: 250.0, height: 520.0 },
                collapsed: false,
                section: Some(SectionId::new("sec-browser").unwrap()),
                column: Some(ColumnId::new("col-left").unwrap()),
                column_index: 0,
                column_width: Some(250.0),
                weight: 1.0,
                preferred_extent: None,
                extent_priority: 0,
                join: None,
            },
        };
        if let Some(sec) = layout.sections.iter_mut().find(|s| s.tabs.contains(&browser_id)) {
            sec.location = win_loc.clone();
            sec.column = None;
        }
        if let Some(inst) = layout.instances.iter_mut().find(|i| i.id == browser_id) {
            inst.memory.location = win_loc;
            inst.memory.column = None;
            inst.memory.column_width = None;
        }
        let saved_bytes = encode(&LayoutEnvelope::new(layout), &reg).unwrap();
        let saved_val = RawJson::from_json_string(String::from_utf8(saved_bytes).unwrap()).unwrap();

        let prefs = Prefs { panel_layout: Some(saved_val), ..Prefs::default() };
        let (projected, store) = FieldCraftLayoutStore::load_from_prefs(&prefs, &reg, false).unwrap();

        // Presented as clamped in-window Floating.
        let sec = projected.containing(&browser_id).unwrap();
        assert_eq!(sec.location.zone(), Zone::Floating);
        if let SectionLocation::Floating { origin } = sec.location {
            assert!(origin.x <= PROJECTED_ORIGIN_MAX && origin.y <= PROJECTED_ORIGIN_MAX);
        }

        // Unrelated Prefs write-back preserves original Zone::Window coordinates.
        let preserved = store.to_prefs_value(&projected, &reg).unwrap();
        let decoded = decode(preserved.as_bytes(), &reg, &seeded_default_layout(&reg, &Prefs::default()).unwrap()).unwrap();
        assert_eq!(decoded.envelope.layout.containing(&browser_id).unwrap().location.zone(), Zone::Window);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn preserve_before_deferred_close_keeps_first_snapshot_across_sequential_closes() {
        let reg = build_registry().unwrap();
        let (seeded, mut store) = FieldCraftLayoutStore::load_from_prefs(&Prefs::default(), &reg, true).unwrap();
        let caps = HostCapabilities { persistence: true, native_windows: true, max_native_hosts: 8 };
        let mut service = solidcraft_panel_layout::LayoutService::new(reg.clone(), seeded, caps).unwrap();
        let browser_id = InstanceId::new(PANEL_BROWSER).unwrap();
        let history_id = InstanceId::new(PANEL_HISTORY).unwrap();

        let apply_mut = |svc: &mut solidcraft_panel_layout::LayoutService, m: Mutation| {
            let p = svc.prepare(LayoutAction::Apply { expected_revision: None, mutations: vec![m] }).unwrap();
            svc.commit(p).unwrap().layout
        };
        let pop_mut = |svc: &solidcraft_panel_layout::LayoutService, id: &InstanceId, host_name: &str, x: f32| {
            let cur_layout = svc.layout();
            let cur = cur_layout.containing(id).unwrap().clone();
            let zone = cur.location.zone();
            let sec_id = next_section_id(cur_layout).unwrap();
            let return_idx = cur_layout.sections.iter().position(|s| s.id == cur.id).unwrap_or(0);
            let insert_idx = cur_layout.sections.len().saturating_sub(1);
            Mutation::Move {
                instance: id.clone(),
                destination: Destination::NewSection {
                    section: sec_id,
                    index: insert_idx,
                    location: SectionLocation::Window {
                        host: HostId::new(host_name).unwrap(),
                        origin: Point { x, y: 120.0 },
                        return_to: DockReturn {
                            zone,
                            index: return_idx,
                            size: cur.size,
                            collapsed: false,
                            section: Some(cur.id),
                            column: cur.column,
                            column_index: 0,
                            column_width: None,
                            weight: 1.0,
                            preferred_extent: None,
                            extent_priority: 0,
                            join: None,
                        },
                    },
                    size: cur.size,
                },
            }
        };

        let m1 = pop_mut(&service, &browser_id, "win-browser", 2000.0);
        let _ = apply_mut(&mut service, m1);
        let m2 = pop_mut(&service, &history_id, "win-history", 2300.0);
        let two_windows = apply_mut(&mut service, m2);

        store.preserve_before_deferred_close(&two_windows, &reg);

        // Simulate first window redocked before second closes: second call must NOT overwrite first snapshot.
        let after_first_close = apply_mut(&mut service, Mutation::NativeClose { host: HostId::new("win-browser").unwrap() });
        store.preserve_before_deferred_close(&after_first_close, &reg);

        let preserved = store.to_prefs_value(&after_first_close, &reg).unwrap();
        let decoded = decode(preserved.as_bytes(), &reg, &seeded_default_layout(&reg, &Prefs::default()).unwrap()).unwrap();
        assert_eq!(decoded.envelope.layout.containing(&browser_id).unwrap().location.zone(), Zone::Window);
        assert_eq!(decoded.envelope.layout.containing(&history_id).unwrap().location.zone(), Zone::Window);

        store.mark_dirty();
        assert!(store.take_preserved_native().is_none());
    }
}
