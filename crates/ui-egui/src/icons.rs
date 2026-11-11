//! Icons: SolidCraft's icon set (`solidcraft_look::icons`), extended with FieldCraft's own masters,
//! drawn by the UI kit's renderer (`solidcraft_ui_kit::icons`).
//!
//! SolidCraft's artwork draws what both apps have (open, save, undo, the view and navigation
//! controls, help, the kit's own close and chevrons, the CAD operations FieldCraft shares:
//! materials, sections, datums, mirror, bodies, components, and what the two apps both do with
//! files, contacts, fasteners, joints, probes and playback): FieldCraft's ids for them are aliases of SolidCraft's
//! masters ([`ALIASES`]). FieldCraft's own artwork (`assets/icons`, original, drawn with
//! VectorCraft by `cargo xtask icons`) draws the CAE concepts SolidCraft has no icon for: meshing,
//! elements, loads and boundary conditions, steps, solving and results. `build.rs` embeds those
//! masters in the set's colour roles and on its 64-unit square ([`MASTERS`]); a master of ours wins
//! a name both sets have.
//!
//! Icons are drawn in the theme's colours; the accent roles take FieldCraft's teal
//! ([`crate::theme`]). Icons named `ui.*` belong to the window itself (navigation bar, title bar),
//! the others to catalogue entries.

use egui::Rect;
use solidcraft_ui_kit::icons::{self as kit_icons, IconExtension, IconSet};

include!(concat!(env!("OUT_DIR"), "/icons.rs"));

/// FieldCraft's ids drawn with SolidCraft's artwork: (FieldCraft id, SolidCraft master, accent).
pub const ALIASES: &[(&str, &str, &str)] = &[
    ("ui.home", "home", ""),
    ("ui.eye", "eye", ""),
    ("ui.fit", "fit", ""),
    ("ui.zoom", "zoom", ""),
    ("ui.pan", "pan", ""),
    ("ui.orbit", "orbit", ""),
    ("ui.perspective", "perspective", ""),
    ("ui.style", "style", ""),
    ("ui.moon", "moon", ""),
    ("ui.help", "help", ""),
    ("file.open", "open", ""),
    ("file.save", "save", ""),
    ("edit.undo", "undo", ""),
    ("edit.redo", "redo", ""),
    ("diagnostics.bug_report", "bug_report", ""),
    ("model.material.create", "material", ""),
    ("results.cut", "section", ""),
    ("geometry.datum", "plane", ""),
    ("geometry.datum.axis", "axis", ""),
    ("geometry.datum.point", "point", ""),
    ("geometry.partition", "split_body", ""),
    ("geometry.stitch", "stitch", ""),
    ("mesh.mirror", "mirror", ""),
    ("mesh.node.move", "move", ""),
    ("ui.search", "zoom", ""),
    ("view.show_all", "eye", ""),
    ("file.recover", "restore", ""),
    ("geometry.open", "import", ""),
    ("geometry.import.mesh", "mesh_body", ""),
    ("geometry.bodies", "body", ""),
    ("geometry.face.split", "split_face", ""),
    ("assembly.instance", "component", ""),
    ("interaction.contact.general", "contact_all", ""),
    ("interaction.contact.pair", "contact_new", ""),
    ("interaction.bolt", "fastener", ""),
    ("interaction.connector", "joint", ""),
    ("deck.read", "import", ""),
    ("deck.write.abaqus", "export", ""),
    ("deck.write.nastran", "export", ""),
    ("deck.write.lsdyna", "export", ""),
    ("mesh.export.stl", "export", ""),
    ("results.read", "import", ""),
    ("results.export.vtu", "export", ""),
    ("results.animate", "play", ""),
    ("solve.probe", "measure", ""),
];

/// SolidCraft's icons with FieldCraft's masters and aliases over them.
pub static SET: IconSet = IconSet::extended(solidcraft_look::icons::DATA, IconExtension { masters: MASTERS, aliases: ALIASES, icon_accents: &[] });

/// Make [`SET`] the icons the kit draws (once per process; again does nothing).
pub fn install() {
    kit_icons::install(&SET);
}

/// How an icon is drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Look {
    /// In the theme's colours.
    Normal,
    /// Greyed: coming soon or unavailable.
    Greyed,
    /// On the title bar, which is dark in both themes.
    AppBar,
}

impl Look {
    pub fn of(enabled: bool) -> Look {
        if enabled { Look::Normal } else { Look::Greyed }
    }
}

/// How faint a greyed icon is.
const GREYED: f32 = 0.35;

/// Paint icon `id` into `rect` (a square).
pub fn paint(ui: &egui::Ui, id: &str, rect: Rect, look: Look) {
    paint_faded(ui, id, rect, look, 1.0);
}

/// [`paint`] at `opacity` (0 to 1).
pub fn paint_faded(ui: &egui::Ui, id: &str, rect: Rect, look: Look, opacity: f32) {
    let p = ui.painter();
    match look {
        Look::Normal => kit_icons::paint_icon(p, rect, id, opacity),
        Look::Greyed => kit_icons::paint_icon(p, rect, id, opacity * GREYED),
        Look::AppBar => kit_icons::paint_app_bar(p, rect, id, opacity),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_engine::catalogue::CATALOGUE;

    /// Every catalogue entry and window icon has artwork (ours or SolidCraft's), none falls back
    /// to the placeholder, and every master of ours belongs to an entry or the window.
    #[test]
    fn every_entry_draws_and_every_master_has_an_entry() {
        let window = [
            "ui.home",
            "ui.eye",
            "ui.fit",
            "ui.zoom",
            "ui.pan",
            "ui.orbit",
            "ui.perspective",
            "ui.style",
            "ui.views",
            "ui.search",
            "ui.sun",
            "ui.moon",
            "ui.help",
        ];
        for id in CATALOGUE.iter().map(|c| c.id).chain(window) {
            let r = SET.resolve(id);
            assert!(!r.placeholder, "{id} has no icon");
        }
        for (id, _, _) in MASTERS {
            assert!(id.starts_with("ui.") || CATALOGUE.iter().any(|c| c.id == *id), "assets/icons/{id}.svg belongs to no catalogue entry");
            assert!(!ALIASES.iter().any(|(a, ..)| a == id), "{id} is both an alias and a master");
        }
    }

    /// The kit reads every master of ours in SolidCraft's colour roles.
    #[test]
    fn the_kit_reads_every_master_of_ours() {
        for (id, family, svg) in MASTERS {
            assert_eq!(*family, "", "{id}");
            if let Err(e) = SET.parse(svg) {
                panic!("{id}: the kit cannot read it: {e}");
            }
        }
    }

    /// Every colour in our masters is one of SolidCraft's roles (`build.rs` mapped them all): a
    /// colour the mapping missed would draw the same in every theme.
    #[test]
    fn our_masters_paint_only_in_solidcrafts_roles() {
        const ROLES: [&str; 5] = ["currentColor", "#236ead", "#e2e8ef", "#c06f00", "none"];
        for (id, _, svg) in MASTERS {
            for attr in ["fill=\"", "stroke=\""] {
                for value in svg.split(attr).skip(1).filter_map(|rest| rest.split('"').next()) {
                    assert!(ROLES.contains(&value), "{id}: {attr}{value}\" is no SolidCraft colour role");
                }
            }
        }
    }

    /// Every alias names a master of SolidCraft's.
    #[test]
    fn aliases_name_solidcraft_masters() {
        for (alias, master, _) in ALIASES {
            assert!(solidcraft_look::icons::SET.has_master(master), "{alias} -> {master}");
            assert_eq!(SET.resolve(alias).master, *master, "{alias}");
        }
    }

    /// Our masters are scaled from the artwork's 24 units to the kit's 64, so they fill their
    /// square as SolidCraft's do.
    #[test]
    fn our_masters_draw_at_full_size() {
        let look = kit_icons::Look::themed(true, None, solidcraft_theme::ColorVision::Off);
        let doc = SET.doc_for("mesh.tet.generate", 64).expect("the tet icon");
        let pixmap = kit_icons::render(doc, look, 64).expect("rendered");
        let inked: Vec<usize> = pixmap.data().chunks(4).enumerate().filter(|(_, p)| p[3] > 0).map(|(i, _)| i % 64).collect();
        let (lo, hi) = (inked.iter().min().copied().unwrap_or(64), inked.iter().max().copied().unwrap_or(0));
        assert!(lo < 16 && hi > 48, "the artwork spans columns {lo}..{hi} of 64");
    }
}
