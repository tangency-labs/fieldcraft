//! The community links: the Donate and Discord pills in the title bar, before Report a bug, drawn
//! by the kit's application bar (`solidcraft_ui_kit::app_bar::Pill`) in the services' own
//! colours (`app_bar::brand`). A click opens the link's page through the host's URL opener
//! ([`ModelApp::open_link`]); a click the control channel sends opens nothing. The small × inside
//! a pill hides it, and the preferences remember that ([`Prefs::show_donate_button`],
//! [`Prefs::show_discord_button`], both on by default).
//!
//! Adapted from SolidCraft's `crates/ui-egui/src/links.rs` (MIT OR Apache-2.0, ATTRIBUTION.md).

use solidcraft_ui_kit::app_bar;

use crate::ModelApp;
use crate::files::Prefs;

pub const DONATE_URL: &str = "https://ko-fi.com/jyviko";
pub const DISCORD_URL: &str = "https://discord.gg/FtU4NQNVDH";

/// One link and its pill.
pub struct Link {
    /// The pill's id: the kit publishes it under the bar's pill prefix, so `region` is that
    /// prefix and the id, and `hide_region` its ×.
    pub id: &'static str,
    pub region: &'static str,
    pub hide_region: &'static str,
    pub url: &'static str,
    /// The service's mark in SolidCraft's icon set.
    pub icon: &'static str,
    /// Message ids (`locales/en.json`): the label, the tooltip, the ×'s tooltip, and the
    /// accessible names of the pill and of its ×.
    pub label: &'static str,
    pub tip: &'static str,
    pub hide_tip: &'static str,
    pub name: &'static str,
    pub hide_name: &'static str,
    /// Whether the pill carries its label on a wide bar (Discord's mark says its name).
    pub named: bool,
    pub fill: egui::Color32,
    /// Whether the title bar shows the pill, and the switch that stops it.
    pub shown: fn(&Prefs) -> bool,
    pub hide: fn(&mut Prefs),
}

/// Where the kit's bar publishes the pills: each [`Link::region`] starts with it.
pub const PILL_PREFIX: &str = "app_bar:";

/// The links, left to right.
pub const LINKS: &[Link] = &[
    Link {
        id: "donate",
        region: "app_bar:donate",
        hide_region: "app_bar:donate:hide",
        url: DONATE_URL,
        icon: "kofi",
        label: "links.donate_label",
        tip: "links.donate_tip",
        hide_tip: "links.donate_hide_tip",
        name: "links.donate_name",
        hide_name: "links.donate_hide_name",
        named: true,
        fill: app_bar::brand::KOFI,
        shown: |p| p.show_donate_button,
        hide: |p| p.show_donate_button = false,
    },
    Link {
        id: "discord",
        region: "app_bar:discord",
        hide_region: "app_bar:discord:hide",
        url: DISCORD_URL,
        icon: "discord",
        label: "links.discord_label",
        tip: "links.discord_tip",
        hide_tip: "links.discord_hide_tip",
        name: "links.discord_name",
        hide_name: "links.discord_hide_name",
        named: false,
        fill: app_bar::brand::DISCORD,
        shown: |p| p.show_discord_button,
        hide: |p| p.show_discord_button = false,
    },
];

pub fn by_id(id: &str) -> Option<&'static Link> {
    LINKS.iter().find(|l| l.id == id)
}

/// The pills the preferences leave on, left to right; `wide` gives the named ones their label.
pub fn pills(prefs: &Prefs, wide: bool) -> Vec<app_bar::Pill<'static>> {
    LINKS
        .iter()
        .filter(|l| (l.shown)(prefs))
        .map(|l| app_bar::Pill {
            id: l.id,
            icon: l.icon,
            label: if l.named { solidcraft_ui_text::tr(l.label) } else { "" },
            labelled: l.named && wide,
            tip: solidcraft_ui_text::tr(l.tip),
            hide_tip: solidcraft_ui_text::tr(l.hide_tip),
            name: solidcraft_ui_text::tr(l.name),
            hide_name: solidcraft_ui_text::tr(l.hide_name),
            fill: l.fill,
            ink: egui::Color32::WHITE,
        })
        .collect()
}

impl ModelApp {
    /// A pill's × was clicked: hide it, and remember that.
    pub(crate) fn hide_link(&mut self, id: &str) {
        if let Some(link) = by_id(id) {
            (link.hide)(&mut self.prefs);
            self.prefs.save(self.services.prefs.as_ref());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every link has its text in the catalogue and its mark in the icon set, and its regions are
    /// where the kit publishes its pill.
    #[test]
    fn every_link_has_text_a_mark_and_the_kits_regions() {
        crate::theme::install_look();
        for link in LINKS {
            for id in [link.label, link.tip, link.hide_tip, link.name, link.hide_name] {
                assert_ne!(solidcraft_ui_text::tr(id), id, "{id} has no English text");
            }
            assert!(crate::icons::SET.has_master(link.icon), "{}", link.icon);
            assert_eq!(link.region, format!("{PILL_PREFIX}{}", link.id));
            assert_eq!(link.hide_region, format!("{}:hide", link.region));
        }
    }

    /// Both pills are on by default, and each switch hides only its own.
    #[test]
    fn both_pills_are_on_by_default_and_each_switch_hides_only_its_own() {
        let ids = |p: &Prefs| pills(p, true).iter().map(|p| p.id).collect::<Vec<_>>();
        let mut prefs = Prefs::default();
        assert_eq!(ids(&prefs), ["donate", "discord"]);
        (LINKS[0].hide)(&mut prefs);
        assert_eq!(ids(&prefs), ["discord"]);
        (LINKS[1].hide)(&mut prefs);
        assert!(ids(&prefs).is_empty());
    }
}
