//! The accent's colour checks: text on it and in it reads (WCAG 2.x contrast), it is teal, and
//! it stays apart from the status colours (ok, warning, error) for normal vision and for the three
//! dichromacies.
//!
//! Dichromacy is simulated with Machado, Oliveira and Fernandes, "A Physiologically-based Model
//! for Simulation of Color Vision Deficiency", IEEE TVCG 15(6), 2009 (the severity 1.0 matrices,
//! applied to linear sRGB). Colour difference is CIEDE2000 as given by Sharma, Wu and Dalal, "The
//! CIEDE2000 Color-Difference Formula: Implementation Notes, Supplementary Test Data, and
//! Mathematical Observations", Color Research & Application 30(1), 2005, checked against their
//! test data below. Contrast is the WCAG 2.x relative luminance ratio.

use egui::Color32;

use super::Tokens;

type V3 = [f64; 3];

/// Smallest CIEDE2000 difference between the accent and a status colour, for every kind of
/// vision: 10 is a difference seen at a glance, well past the just-noticeable ~1–2.
const STATUS_FLOOR: f64 = 10.0;

fn linear(c: u8) -> f64 {
    let c = f64::from(c) / 255.0;
    if c <= 0.04045 { c / 12.92 } else { ((c + 0.055) / 1.055).powf(2.4) }
}

fn linear_rgb(c: Color32) -> V3 {
    [linear(c.r()), linear(c.g()), linear(c.b())]
}

fn luminance(c: Color32) -> f64 {
    let [r, g, b] = linear_rgb(c);
    0.2126 * r + 0.7152 * g + 0.0722 * b
}

/// WCAG 2's enhanced contrast for text (success criterion 1.4.6, level AAA).
const WCAG_AAA: f64 = 7.0;

/// WCAG 2.x contrast ratio, 1 to 21.
fn contrast(a: Color32, b: Color32) -> f64 {
    let (la, lb) = (luminance(a), luminance(b));
    (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
}

#[derive(Clone, Copy, Debug)]
enum Vision {
    Normal,
    Protan,
    Deutan,
    Tritan,
}

/// Machado et al. 2009, severity 1.0, on linear RGB.
fn simulate(c: Color32, vision: Vision) -> V3 {
    let m: [V3; 3] = match vision {
        Vision::Normal => return linear_rgb(c),
        Vision::Protan => [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
        Vision::Deutan => [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
        Vision::Tritan => [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
    };
    let v = linear_rgb(c);
    m.map(|row| (row[0] * v[0] + row[1] * v[1] + row[2] * v[2]).clamp(0.0, 1.0))
}

/// CIELAB (D65 white) of a linear sRGB colour.
fn lab([r, g, b]: V3) -> V3 {
    let x = (0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / 0.95047;
    let y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
    let z = (0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / 1.08883;
    let f = |t: f64| if t > 216.0 / 24389.0 { t.cbrt() } else { (24389.0 / 27.0 * t + 16.0) / 116.0 };
    let (fx, fy, fz) = (f(x), f(y), f(z));
    [116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

/// CIEDE2000 (Sharma, Wu and Dalal 2005), kL = kC = kH = 1.
fn de2000([l1, a1, b1]: V3, [l2, a2, b2]: V3) -> f64 {
    let c7 = |c: f64| c.powi(7) / (c.powi(7) + 25f64.powi(7));
    let c_bar = (a1.hypot(b1) + a2.hypot(b2)) / 2.0;
    let g = 0.5 * (1.0 - c7(c_bar).sqrt());
    let (a1p, a2p) = ((1.0 + g) * a1, (1.0 + g) * a2);
    let (c1p, c2p) = (a1p.hypot(b1), a2p.hypot(b2));
    let hue = |b: f64, a: f64| if b == 0.0 && a == 0.0 { 0.0 } else { b.atan2(a).to_degrees().rem_euclid(360.0) };
    let (h1, h2) = (hue(b1, a1p), hue(b2, a2p));
    let chromatic = c1p * c2p != 0.0;
    let dh = if !chromatic {
        0.0
    } else if (h2 - h1).abs() <= 180.0 {
        h2 - h1
    } else if h2 - h1 > 180.0 {
        h2 - h1 - 360.0
    } else {
        h2 - h1 + 360.0
    };
    let (dl, dc) = (l2 - l1, c2p - c1p);
    let d_big_h = 2.0 * (c1p * c2p).sqrt() * (dh / 2.0).to_radians().sin();
    let l_bar = (l1 + l2) / 2.0;
    let c_bar_p = (c1p + c2p) / 2.0;
    let h_bar = if !chromatic {
        h1 + h2
    } else if (h1 - h2).abs() <= 180.0 {
        (h1 + h2) / 2.0
    } else if h1 + h2 < 360.0 {
        (h1 + h2 + 360.0) / 2.0
    } else {
        (h1 + h2 - 360.0) / 2.0
    };
    let cos = |deg: f64| deg.to_radians().cos();
    let t = 1.0 - 0.17 * cos(h_bar - 30.0) + 0.24 * cos(2.0 * h_bar) + 0.32 * cos(3.0 * h_bar + 6.0) - 0.20 * cos(4.0 * h_bar - 63.0);
    let d_theta = 30.0 * (-((h_bar - 275.0) / 25.0).powi(2)).exp();
    let r_c = 2.0 * c7(c_bar_p).sqrt();
    let s_l = 1.0 + 0.015 * (l_bar - 50.0).powi(2) / (20.0 + (l_bar - 50.0).powi(2)).sqrt();
    let s_c = 1.0 + 0.045 * c_bar_p;
    let s_h = 1.0 + 0.015 * c_bar_p * t;
    let r_t = -(2.0 * d_theta).to_radians().sin() * r_c;
    let (l, c, h) = (dl / s_l, dc / s_c, d_big_h / s_h);
    (l * l + c * c + h * h + r_t * c * h).sqrt()
}

/// The accent's hue in OKLab (Ottosson 2020), degrees 0–360.
fn oklab_hue(c: Color32) -> f64 {
    let [r, g, b] = linear_rgb(c);
    let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
    let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
    let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
    let a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
    let b = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
    b.atan2(a).to_degrees().rem_euclid(360.0)
}

fn themes() -> [(&'static str, Tokens); 2] {
    [("light", Tokens::light()), ("dark", Tokens::dark())]
}

/// The implementation reproduces Sharma, Wu and Dalal's published test pairs (Table 1).
#[test]
fn ciede2000_matches_the_published_test_data() {
    let pairs: [(V3, V3, f64); 5] = [
        ([50.0, 2.6772, -79.7751], [50.0, 0.0, -82.7485], 2.0425),
        ([50.0, 0.0, 0.0], [50.0, -1.0, 2.0], 2.3669),
        ([50.0, 2.5, 0.0], [50.0, 0.0, -2.5], 4.3065),
        ([50.0, 2.5, 0.0], [73.0, 25.0, -18.0], 27.1492),
        ([60.2574, -34.0099, 36.2677], [60.4626, -34.1751, 39.4387], 1.2644),
    ];
    for (a, b, want) in pairs {
        let got = de2000(a, b);
        assert!((got - want).abs() < 1e-4, "{a:?} {b:?}: {got} != {want}");
    }
}

/// Text and marks on an accent fill, resting or hovered, read at 4.5:1 (WCAG AA body text).
#[test]
fn text_on_the_accent_reads() {
    for (name, t) in themes() {
        for fill in [t.accent, t.accent_hover] {
            let ratio = contrast(t.on_accent, fill);
            assert!(ratio >= 4.5, "{name}: on_accent on {fill:?} is {ratio:.2}:1");
        }
    }
}

/// Text in the accent reads at 4.5:1 on every ground it is written on.
#[test]
fn accent_text_reads_on_panels() {
    for (name, t) in themes() {
        for ground in [t.panel, t.panel_header, t.toolbar, t.dialog_bg, t.heading_bg, t.subheading_bg, t.accent_soft] {
            let ratio = contrast(t.accent_ink, ground);
            assert!(ratio >= 4.5, "{name}: accent_ink on {ground:?} is {ratio:.2}:1");
        }
    }
}

/// Selected labels and buttons, whose text egui draws in the selection stroke's colour on the
/// selection fill, read at 4.5:1.
#[test]
fn selected_text_reads() {
    for dark in [false, true] {
        let v = super::tests::visuals(dark);
        let ratio = contrast(v.selection.stroke.color, v.selection.bg_fill);
        assert!(ratio >= 4.5, "dark={dark}: selected text is {ratio:.2}:1");
    }
}

/// FieldCraft's accent is teal: neither SolidCraft's blue, a violet nor orange (OKLab hue ~187°; blue is
/// ~255°, violet ~297°, orange ~55°).
#[test]
fn the_accent_is_teal() {
    for (name, t) in themes() {
        for c in [t.accent, t.accent_ink, t.accent_hover] {
            let hue = oklab_hue(c);
            assert!((175.0..=200.0).contains(&hue), "{name}: {c:?} has hue {hue:.0}°");
        }
    }
}

/// The accent never reads as ok, warning or error, for normal vision or any dichromacy.
#[test]
fn the_accent_stays_apart_from_status_colours() {
    for (name, t) in themes() {
        for (status, c) in [("ok", t.ok), ("warning", t.warm), ("error", t.error)] {
            for vision in [Vision::Normal, Vision::Protan, Vision::Deutan, Vision::Tritan] {
                let d = de2000(lab(simulate(t.accent, vision)), lab(simulate(c, vision)));
                assert!(d >= STATUS_FLOOR, "{name}: accent vs {status} under {vision:?} is ΔE00 {d:.1}");
            }
        }
    }
}

/// The soft accent and the hover tint sit between the panel and the accent, hover the fainter, so
/// a selected row stays visible under the pointer.
#[test]
fn soft_accents_are_ordered() {
    for (name, t) in themes() {
        let d = |c: Color32| de2000(lab(linear_rgb(t.panel)), lab(linear_rgb(c)));
        assert!(d(t.hover) < d(t.accent_soft) && d(t.accent_soft) < d(t.accent), "{name}");
        assert!(d(t.hover) > 2.0, "{name}: hover is invisible on the panel");
    }
}

/// The kit's tokens carry the accent wherever SolidCraft's carry its blue: every accent slot is
/// teal, so no kit component draws SolidCraft's accent.
#[test]
fn the_kit_has_no_other_accent() {
    for (name, dark) in [("light", false), ("dark", true)] {
        let k = super::kit(dark);
        for c in [
            k.accent,
            k.selected_ink,
            k.form_primary,
            k.form_primary_fill,
            k.figure_accent,
            k.file_primary_ink,
            k.file_primary_border,
            k.file_rail,
            k.file_rail_hover,
            k.chip,
            k.agent_busy,
            k.sketch_accent,
            k.surface_accent,
            k.mesh_accent,
            k.form_accent,
            k.render_accent,
        ] {
            let hue = oklab_hue(c);
            assert!((175.0..=200.0).contains(&hue), "{name}: {c:?} has hue {hue:.0}°");
        }
    }
}

/// Text the kit writes in the accent reads at 4.5:1 on every kit ground it is written on, and
/// text on an accent fill reads on it.
#[test]
fn the_kits_accent_text_reads() {
    for (name, dark) in [("light", false), ("dark", true)] {
        let k = super::kit(dark);
        for ink in [k.selected_ink, k.form_primary, k.file_primary_ink] {
            for ground in [k.panel, k.panel_header, k.dialog_body, k.dialog_footer, k.control_fill, k.selected_bg, k.accent_soft, k.file_primary] {
                let ratio = contrast(ink, ground);
                assert!(ratio >= 4.5, "{name}: {ink:?} on {ground:?} is {ratio:.2}:1");
            }
        }
        for fill in [k.accent, k.form_primary_fill, k.file_rail_hover] {
            let ratio = contrast(k.on_accent, fill);
            assert!(ratio >= 4.5, "{name}: on_accent on {fill:?} is {ratio:.2}:1");
        }
    }
}

/// The selected document tab stands out from the title bar as much in the light theme as in the
/// dark one (the bar is dark in both), and its name reads on it.
#[test]
fn the_document_tab_stands_out_on_the_title_bar() {
    for (name, dark) in [("light", false), ("dark", true)] {
        let k = super::kit(dark);
        let ratio = contrast(k.app_bar_tab_active, k.app_bar);
        assert!(ratio >= 1.35, "{name}: the active tab on the bar is {ratio:.2}:1");
        let ratio = contrast(k.app_bar_text, k.app_bar_tab_active);
        assert!(ratio >= 4.5, "{name}: the tab's name on it is {ratio:.2}:1");
    }
}

/// OKLab chroma (Ottosson 2020): 0 for a grey.
fn oklab_chroma(c: Color32) -> f64 {
    let [r, g, b] = linear_rgb(c);
    let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
    let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
    let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
    let a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s;
    let b = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s;
    a.hypot(b)
}

/// Half a named hue's sector: the hue circle read as six named hues (red, yellow, green, cyan,
/// blue, magenta) gives each 60°, so a colour within 30° of SolidCraft's accent hue reads as its
/// blue.
const BLUE_HALF_SECTOR: f64 = 30.0;

/// Kit colours that stay SolidCraft's however blue they are: its mark on the application bar,
/// which neither FieldCraft nor the kit draws.
const NOT_AN_ACCENT: [&str; 1] = ["logo"];

/// Every colour slot of the kit's tokens, by the kit's own inventory of them (the theme-file
/// keys, which the kit keeps complete at build time), with its colours.
fn kit_slots(mut t: solidcraft_theme::Tokens) -> Vec<(&'static str, Vec<Color32>)> {
    solidcraft_ui_kit::theme_files::keys()
        .into_iter()
        .map(|key| (key, solidcraft_ui_kit::theme_files::slot(&mut t, key).map(|s| s.to_vec()).unwrap_or_default()))
        .collect()
}

/// The kit slots SolidCraft draws in its own 3D view (`Tokens::with_viewport` takes exactly
/// those from the view's tokens): FieldCraft draws its view with its own colours.
fn view_slots(base: solidcraft_theme::Tokens) -> Vec<&'static str> {
    let mut other = base;
    for key in solidcraft_ui_kit::theme_files::keys() {
        for c in solidcraft_ui_kit::theme_files::slot(&mut other, key).into_iter().flatten() {
            // Every channel complemented: never equal to the colour it replaces.
            *c = Color32::from_rgba_premultiplied(255 - c.r(), 255 - c.g(), 255 - c.b(), 255 - c.a());
        }
    }
    let swapped = kit_slots(base.with_viewport(other));
    kit_slots(base).into_iter().zip(swapped).filter(|((_, a), (_, b))| a != b).map(|((key, _), _)| key).collect()
}

/// No kit slot keeps SolidCraft's blue: every slot SolidCraft gives a colour of its accent's hue
/// (more chromatic than its greys) takes another colour from [`super::kit`]. The slots are walked
/// from the kit's inventory, so a slot added to the kit, or one missed, fails here.
#[test]
fn the_kit_keeps_no_solidcraft_blue() {
    for (name, dark) in [("light", false), ("dark", true)] {
        let sc = if dark { solidcraft_theme::Tokens::dark() } else { solidcraft_theme::Tokens::light() };
        // SolidCraft's greys (grounds, rules, dim text) lean blue; its accents are more chromatic.
        let greys = [
            sc.app_bar,
            sc.toolbar,
            sc.panel,
            sc.panel_header,
            sc.border,
            sc.text,
            sc.text_dim,
            sc.field,
            sc.dialog_bg,
            sc.dialog_body,
            sc.dialog_label,
            sc.form_dim,
            sc.file_dim,
            sc.nav_dim,
            sc.nav_press,
        ];
        let grey_chroma = greys.into_iter().map(oklab_chroma).fold(0.0, f64::max);
        let blue_hue = oklab_hue(sc.accent);
        let is_blue = |c: Color32| {
            let off = (oklab_hue(c) - blue_hue).rem_euclid(360.0);
            oklab_chroma(c) > grey_chroma && off.min(360.0 - off) <= BLUE_HALF_SECTOR
        };
        let view = view_slots(sc);
        let kept: Vec<_> = kit_slots(sc)
            .into_iter()
            .zip(kit_slots(super::kit(dark)))
            .filter(|((key, _), _)| !view.contains(key) && !NOT_AN_ACCENT.contains(key))
            .flat_map(|((key, ours), (_, theirs))| ours.into_iter().zip(theirs).map(move |(s, f)| (key, s, f)))
            .filter(|&(_, s, f)| is_blue(s) && s == f)
            .map(|(key, s, _)| format!("{key} {s:?}"))
            .collect();
        assert!(kept.is_empty(), "{name}: SolidCraft's blue kept in {kept:?}");
    }
}

/// The start page's rail: its text reads on it at WCAG AAA (what [`super::RAIL_DEPTH`] is chosen
/// for) and on a hovered or selected row at AA, and such a row stands out from the rail at least
/// as much as it does in SolidCraft's.
#[test]
fn the_start_rail_reads_and_its_rows_stand_out() {
    for (name, dark) in [("light", false), ("dark", true)] {
        let sc = if dark { solidcraft_theme::Tokens::dark() } else { solidcraft_theme::Tokens::light() };
        let k = super::kit(dark);
        let rail = contrast(k.file_rail_ink, k.file_rail);
        assert!(rail >= WCAG_AAA, "{name}: rail text on the rail is {rail:.2}:1");
        let row = contrast(k.file_rail_ink, k.file_rail_hover);
        assert!(row >= 4.5, "{name}: rail text on a hovered row is {row:.2}:1");
        let row = |t: &solidcraft_theme::Tokens| de2000(lab(linear_rgb(t.file_rail)), lab(linear_rgb(t.file_rail_hover)));
        assert!(row(&k) >= row(&sc), "{name}: a row on the rail is ΔE00 {:.1}, SolidCraft's {:.1}", row(&k), row(&sc));
    }
}

/// The agent's states on the title bar (busy, ready, stopped) stay apart for colour-blind users:
/// busy against ready and against stopped (the kit draws stopped in `error`) clears the status
/// floor under each dichromacy (Machado et al. 2009) in both themes. Its arc differs from their
/// dot too, but shape alone is not enough.
#[test]
fn the_agent_states_stay_apart_for_colour_blind_users() {
    for (name, dark) in [("light", false), ("dark", true)] {
        let k = super::kit(dark);
        for (state, c) in [("ready", k.agent_ok), ("stopped", k.error)] {
            for vision in [Vision::Protan, Vision::Deutan, Vision::Tritan] {
                let d = de2000(lab(simulate(k.agent_busy, vision)), lab(simulate(c, vision)));
                assert!(d >= STATUS_FLOOR, "{name}: busy vs {state} under {vision:?} is ΔE00 {d:.1}");
            }
        }
        // The mark is read on the title bar, at rest and under the pointer: the kit draws it on
        // the bar with its translucent hover fill laid over (`app_bar.blend(fill)`).
        for bar in [k.app_bar, k.app_bar.blend(k.app_bar_hover_button), k.app_bar.blend(k.app_bar_hover_strong)] {
            let ratio = contrast(k.agent_busy, bar);
            assert!(ratio >= 4.5, "{name}: busy on {bar:?} is {ratio:.2}:1");
        }
    }
}
