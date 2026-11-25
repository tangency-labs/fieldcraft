//! A character none of the interface's faces has (SolidCraft's IBM Plex, egui's fallbacks, the
//! Hebrew face of [`crate::theme`]) is drawn with an installed font that has it, not as an empty
//! box: a file, folder or body named in Japanese, Korean or Chinese (#991).
//!
//! The installed fonts are listed ([`fontdb`]), and the characters each one maps read with
//! [`skrifa`], the reader egui draws with, so a face is offered only for a character egui will find
//! in it. That is done once per process, on a thread of its own, and only then is
//! [`SystemFallback`] added to a context, so the window never waits for the listing: adding a
//! provider makes egui ask again for the characters it could not draw, and repaint. egui asks
//! [`SystemFallback::font_for`] at most once per family and character, and appends the font it gets
//! to that family's fallback chain, so the rest of the script follows from it. The first answer
//! from a font file reads that file, on the UI thread (the listing has just read it, so it comes
//! from the system's cache), once for all its faces. The window and an offscreen capture share the
//! one listing.
//!
//! Only the desktop app turns it on (`Services::system_fonts`); tests draw without it, so their
//! text does not depend on the machine's fonts. A wasm build has no installed fonts and no module.
//!
//! Adapted from SolidCraft's `crates/ui-egui/src/fonts/system_fallback.rs` at 7371a01a
//! (tangency-labs/solidcraft#20030, MIT OR Apache-2.0).

use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Condvar, Mutex, PoisonError};
use std::time::Duration;

use egui::epaint::text::Blob;
use egui::{FallbackRequest, FontData, FontFamily, FontInsert, FontPriority, InsertFontFamily};
use skrifa::MetadataProvider as _;

/// The process's one listing of the installed fonts.
static LISTING: Mutex<Listing> = Mutex::new(Listing::Idle);
/// Signalled when the listing is done and every waiting context has the provider.
static LISTED: Condvar = Condvar::new();

/// The longest [`wait`] waits: a cold listing of a few thousand faces takes seconds.
const LISTING_WAIT: Duration = Duration::from_secs(30);

enum Listing {
    Idle,
    /// Being listed; these contexts get the provider when it is done.
    Running(Vec<egui::Context>),
    Ready(Arc<SystemFallback>),
}

/// Adds the installed fonts to `ctx` once they are listed, the listing started by the first
/// context. Once per context: the host calls it when it styles one (`Services::system_fonts`).
pub(crate) fn install(ctx: &egui::Context) {
    if ctx.data_mut(|d| std::mem::replace(d.get_temp_mut_or_default::<bool>(given()), true)) {
        return;
    }
    let mut listing = LISTING.lock().unwrap_or_else(PoisonError::into_inner);
    match &mut *listing {
        Listing::Ready(provider) => ctx.add_font_provider(provider.clone()),
        Listing::Running(waiting) => waiting.push(ctx.clone()),
        Listing::Idle => {
            let started = std::thread::Builder::new().name("font-fallback".into()).spawn(|| {
                let mut db = fontdb::Database::new();
                db.load_system_fonts();
                // In a Flatpak sandbox `/usr/share/fonts` holds the runtime's fonts; the host's
                // system, local and user fonts are mounted here. A folder that is missing is
                // skipped.
                #[cfg(target_os = "linux")]
                for dir in ["/run/host/fonts", "/run/host/local-fonts", "/run/host/user-fonts"] {
                    db.load_fonts_dir(dir);
                }
                // fontdb reads `~/.local/share/fonts`, not `$XDG_DATA_HOME/fonts` when that is
                // set elsewhere.
                #[cfg(target_os = "linux")]
                if let Some(data) = std::env::var_os("XDG_DATA_HOME").map(std::path::PathBuf::from).filter(|d| d.is_absolute()) {
                    db.load_fonts_dir(data.join("fonts"));
                }
                let provider = Arc::new(SystemFallback::over(db));
                let mut listing = LISTING.lock().unwrap_or_else(PoisonError::into_inner);
                // Each repaints, and egui asks again for the characters it missed. Under the lock,
                // so a context that [`wait`]s finds the provider added.
                if let Listing::Running(waiting) = std::mem::replace(&mut *listing, Listing::Ready(provider.clone())) {
                    for ctx in waiting {
                        ctx.add_font_provider(provider.clone());
                    }
                }
                drop(listing);
                LISTED.notify_all();
            });
            match started {
                Ok(_) => *listing = Listing::Running(vec![ctx.clone()]),
                Err(error) => log::warn!("no installed fonts for missing characters: {error}"),
            }
        }
    }
}

/// Waits until the installed fonts are listed, if they are being listed, and every context
/// [`install`] was given has them: for an offscreen capture, a still frame that must not show
/// boxes where the window will show text. At most [`LISTING_WAIT`]; `false` when it is still
/// listing then.
pub(crate) fn wait() -> bool {
    let listing = LISTING.lock().unwrap_or_else(PoisonError::into_inner);
    let (listing, _) = LISTED.wait_timeout_while(listing, LISTING_WAIT, |l| matches!(l, Listing::Running(_))).unwrap_or_else(PoisonError::into_inner);
    !matches!(*listing, Listing::Running(_))
}

/// Marks a context that has been given the installed fonts (or is waiting for them).
fn given() -> egui::Id {
    egui::Id::unique("fc_system_font_fallback")
}

/// An installed face egui can draw, and the characters it has, as sorted, disjoint ranges.
struct Face {
    id: fontdb::ID,
    chars: Box<[(u32, u32)]>,
}

impl Face {
    fn has(&self, c: char) -> bool {
        let c = u32::from(c);
        let i = self.chars.partition_point(|&(_, last)| last < c);
        self.chars.get(i).is_some_and(|&(first, _)| first <= c)
    }
}

/// An [`egui::FontProvider`] that finds a character in the installed fonts.
pub(super) struct SystemFallback {
    db: fontdb::Database,
    faces: Vec<Face>,
    /// The face found for a character, for a monospaced family or not (looked up, never
    /// iterated).
    found: Mutex<HashMap<(char, bool), Option<fontdb::ID>>>,
    /// Each font file a face was handed to egui from, read once: egui asks again after every
    /// `set_fonts` (each theme change and capture), and the faces of a collection share it.
    files: Mutex<HashMap<PathBuf, Blob>>,
}

impl SystemFallback {
    /// Over the faces of `db`, reading the characters each one has. Slow: not on the UI thread.
    fn over(db: fontdb::Database) -> Self {
        let mut faces = Vec::new();
        // The faces of a collection are listed together: their file is read once.
        let mut file: Option<(&Path, Blob)> = None;
        for f in db.faces() {
            let data = match (path_of(&f.source), &file) {
                (Some(path), Some((last, data))) if path == *last => data.clone(),
                (path, _) => {
                    let Some(data) = bytes_of(&f.source) else { continue };
                    file = path.map(|p| (p, data.clone()));
                    data
                }
            };
            if let Some(chars) = chars_of((*data).as_ref(), f.index) {
                faces.push(Face { id: f.id, chars });
            }
        }
        Self { db, faces, found: Mutex::new(HashMap::new()), files: Mutex::new(HashMap::new()) }
    }

    /// The installed faces that draw `c`, best first: one of the family's kind (monospaced or
    /// not), upright, normal width, nearest the regular weight, a Noto face before others, then by
    /// name and file so the order does not depend on the order the folders were read in.
    fn ranked(&self, c: char, monospace: bool) -> Vec<&fontdb::FaceInfo> {
        let mut faces: Vec<_> = self.faces.iter().filter(|f| f.has(c)).filter_map(|f| self.db.face(f.id)).collect();
        faces.sort_by_cached_key(|f| {
            let name = f.families.first().map(|(n, _)| n.as_str()).unwrap_or_default();
            (
                f.monospaced != monospace,
                f.style != fontdb::Style::Normal,
                f.stretch != fontdb::Stretch::Normal,
                f.weight.0.abs_diff(fontdb::Weight::NORMAL.0),
                !name.starts_with("Noto"),
                name.to_owned(),
                path_of(&f.source),
                f.index,
            )
        });
        faces
    }

    /// The bytes of the file that holds `face`, read once per file. `None` when it cannot be read
    /// (removed or changed since it was listed).
    fn bytes(&self, face: &fontdb::FaceInfo) -> Option<Blob> {
        let Some(path) = path_of(&face.source) else { return bytes_of(&face.source) };
        match self.files.lock().unwrap_or_else(PoisonError::into_inner).entry(path.to_owned()) {
            Entry::Occupied(file) => Some(file.get().clone()),
            Entry::Vacant(slot) => Some(slot.insert(bytes_of(&face.source)?).clone()),
        }
    }

    /// `face` with the bytes of its file, if they can be read.
    fn read<'a>(&self, face: &'a fontdb::FaceInfo) -> Option<(&'a fontdb::FaceInfo, Blob)> {
        self.bytes(face).map(|bytes| (face, bytes))
    }

    /// The face that draws `c` for a family of the kind, with its bytes: the one found before
    /// while its file can still be read, else the best whose file can (and it is kept).
    fn face_for(&self, c: char, monospace: bool) -> Option<(&fontdb::FaceInfo, Blob)> {
        let key = (c, monospace);
        let found = self.found.lock().unwrap_or_else(PoisonError::into_inner).get(&key).copied();
        match found {
            Some(None) => return None,
            Some(Some(id)) => {
                if let Some(face) = self.db.face(id).and_then(|f| self.read(f)) {
                    return Some(face);
                }
            }
            None => {}
        }
        let face = self.ranked(c, monospace).into_iter().find_map(|f| self.read(f));
        self.found.lock().unwrap_or_else(PoisonError::into_inner).insert(key, face.as_ref().map(|(f, _)| f.id));
        face
    }
}

impl egui::FontProvider for SystemFallback {
    fn font_for(&self, request: &FallbackRequest<'_>) -> Option<FontInsert> {
        let c = request.cluster.chars().next()?;
        let monospace = matches!(request.family, FontFamily::Monospace);
        let (face, bytes) = self.face_for(c, monospace)?;
        let data = FontData::from_blob(bytes, face.index);
        // The database id too: two files can hold faces of one PostScript name and index (two
        // versions, or subsets) with different characters, and egui keeps the first of a name.
        let name = format!("system:{}:{}:{:?}", face.post_script_name, face.index, face.id);
        log::debug!("{name} draws {c:?}");
        Some(FontInsert::new(&name, data, vec![InsertFontFamily { family: request.family.clone(), priority: FontPriority::Lowest }]))
    }
}

/// The file a face was read from; `None` for one held in memory.
fn path_of(source: &fontdb::Source) -> Option<&Path> {
    match source {
        fontdb::Source::Binary(_) => None,
        fontdb::Source::File(path) => Some(path),
    }
}

/// The bytes of the file that holds a face: read, not mapped, so a file changed meanwhile cannot
/// fault the process.
fn bytes_of(source: &fontdb::Source) -> Option<Blob> {
    match source {
        fontdb::Source::Binary(data) => Some(data.clone()),
        fontdb::Source::File(path) => match std::fs::read(path) {
            Ok(data) => Some(Arc::new(data)),
            Err(error) => {
                log::debug!("font {} not read: {error}", path.display());
                None
            }
        },
    }
}

/// The characters face `index` of `data` maps to a glyph, as egui reads them (skrifa's character
/// map: the Unicode subtable of format 4 or 12 it picks, up to `char::MAX` and the face's glyph
/// count, so a malformed range cannot run away), as sorted, disjoint ranges. `None` when it maps
/// none, as a last-resort face such as macOS's `LastResort.otf` (whose format 13 subtable egui does
/// not read), or egui cannot draw it: without epaint's `color_fonts`, which FieldCraft does not
/// turn on, it draws outlines only, so a face of bitmaps alone (Noto Color Emoji) is left out.
fn chars_of(data: &[u8], index: u32) -> Option<Box<[(u32, u32)]>> {
    let font = skrifa::FontRef::from_index(data, index).ok()?;
    let drawable = [b"glyf", b"CFF ", b"CFF2"].into_iter().any(|t| font.table_data(skrifa::Tag::new(t)).is_some());
    if !drawable {
        return None;
    }
    let mut chars: Vec<u32> = font.charmap().mappings().map(|(c, _)| c).collect();
    chars.sort_unstable();
    chars.dedup();
    let mut ranges: Vec<(u32, u32)> = Vec::new();
    for c in chars {
        match ranges.last_mut() {
            Some((_, last)) if *last + 1 == c => *last = c,
            _ => ranges.push((c, c)),
        }
    }
    (!ranges.is_empty()).then(|| ranges.into_boxed_slice())
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui::FontProvider as _;

    /// Latin capital letter schwa: in Ubuntu Light, not in Hack.
    const SCHWA: char = 'Ə';

    /// egui's two default faces, Hack (monospaced, regular) and Ubuntu Light (proportional,
    /// light).
    fn defaults() -> fontdb::Database {
        let mut db = fontdb::Database::new();
        let defs = egui::FontDefinitions::default();
        for name in ["Hack", "Ubuntu-Light"] {
            if let Some(font) = defs.font_data.get(name) {
                db.load_font_data(font.bytes().to_vec());
            }
        }
        db
    }

    fn family_of(p: &SystemFallback, c: char, monospace: bool) -> Option<String> {
        p.ranked(c, monospace).first().and_then(|f| f.families.first()).map(|(n, _)| n.clone())
    }

    #[test]
    fn a_family_takes_a_face_of_its_kind_and_only_one_that_has_the_character() {
        let p = SystemFallback::over(defaults());
        assert_eq!(p.faces.len(), 2);
        assert_eq!(family_of(&p, 'a', true).as_deref(), Some("Hack"));
        assert_eq!(family_of(&p, 'a', false).as_deref(), Some("Ubuntu"));
        // Only Ubuntu has it, so even a monospaced family takes it.
        assert_eq!(family_of(&p, SCHWA, true).as_deref(), Some("Ubuntu"));
        assert_eq!(family_of(&p, '한', false), None);
    }

    #[test]
    fn a_face_has_the_characters_of_its_ranges_only() {
        let face = Face { id: fontdb::ID::dummy(), chars: Box::new([(0x41, 0x5a), (0x61, 0x7a), (0x259, 0x259)]) };
        for c in ['A', 'Z', 'a', 'm', 'z', 'ə'] {
            assert!(face.has(c), "{c:?}");
        }
        for c in ['\0', '@', '[', '`', '{', 'Ə', '한', char::MAX] {
            assert!(!face.has(c), "{c:?}");
        }
    }

    #[test]
    fn a_character_is_looked_up_once_per_kind_of_family() {
        let p = SystemFallback::over(defaults());
        for _ in 0..2 {
            assert!(p.font_for(&FallbackRequest { cluster: "Ə", family: &FontFamily::Monospace }).is_some());
            assert!(p.font_for(&FallbackRequest { cluster: "한", family: &FontFamily::Proportional }).is_none());
        }
        let found = p.found.lock().unwrap_or_else(PoisonError::into_inner);
        assert_eq!(found.len(), 2);
        assert!(found[&(SCHWA, true)].is_some());
        assert_eq!(found[&('한', false)], None);
    }

    #[test]
    fn a_context_without_the_character_draws_it_with_the_found_face() {
        let mut only_hack = egui::FontDefinitions::empty();
        let defs = egui::FontDefinitions::default();
        let Some(hack) = defs.font_data.get("Hack") else { panic!("egui's default fonts") };
        only_hack.font_data.insert("Hack".into(), hack.clone());
        only_hack.families.insert(FontFamily::Proportional, vec!["Hack".into()]);
        only_hack.families.insert(FontFamily::Monospace, vec!["Hack".into()]);
        let ctx = egui::Context::default();
        ctx.set_fonts(only_hack);
        ctx.add_font_provider(Arc::new(SystemFallback::over(defaults())));
        ctx.run_ui(egui::RawInput::default(), |_| {}).textures_delta.clear();
        assert!(ctx.fonts_mut(|f| f.has_glyph(&egui::FontId::proportional(14.0), SCHWA)));
        assert!(!ctx.fonts_mut(|f| f.has_glyph(&egui::FontId::proportional(14.0), '한')));
    }

    /// A font of a `cmap` with one subtable (`platform`, `encoding`) of format 12 or 13 mapping
    /// `first..=last` from glyph 1, and a glyph table `glyphs` (`glyf`, or `CBDT` for bitmaps):
    /// no more than skrifa needs to read it.
    fn font(platform: u16, encoding: u16, format: u16, first: u32, last: u32, glyphs: &[u8; 4]) -> Vec<u8> {
        let mut cmap = Vec::new();
        for v in [0u16, 1, platform, encoding] {
            cmap.extend(v.to_be_bytes());
        }
        cmap.extend(12u32.to_be_bytes());
        cmap.extend(format.to_be_bytes());
        cmap.extend(0u16.to_be_bytes());
        for v in [28u32, 0, 1, first, last, 1] {
            cmap.extend(v.to_be_bytes());
        }
        let tables: [(&[u8; 4], &[u8]); 2] = [(b"cmap", &cmap), (glyphs, &[0; 4])];
        let mut font = 0x0001_0000u32.to_be_bytes().to_vec();
        font.extend(2u16.to_be_bytes());
        font.extend([0; 6]);
        let mut offset = 12 + 16 * tables.len();
        for (tag, data) in tables {
            font.extend(tag);
            font.extend(0u32.to_be_bytes());
            for v in [offset, data.len()] {
                font.extend(u32::try_from(v).unwrap_or_default().to_be_bytes());
            }
            offset += data.len();
        }
        for (_, data) in tables {
            font.extend(data);
        }
        font
    }

    #[test]
    fn a_face_has_the_characters_egui_finds_in_it() {
        // Format 12, Windows UCS-4: egui reads it.
        assert_eq!(chars_of(&font(3, 10, 12, 0x4e00, 0x4e10, b"glyf"), 0).as_deref(), Some(&[(0x4e00, 0x4e10)][..]));
        // Format 13, as in macOS's LastResort.otf: it claims every character, egui reads none,
        // so the face is left out rather than picked before the face that draws them.
        assert_eq!(chars_of(&font(0, 6, 13, 0, 0x10_ffff, b"glyf"), 0), None);
        assert_eq!(chars_of(&font(3, 10, 13, 0, 0x10_ffff, b"glyf"), 0), None);
        // Bitmaps alone, as Noto Color Emoji: egui without `color_fonts` draws nothing from it.
        assert_eq!(chars_of(&font(3, 10, 12, 0x1f600, 0x1f610, b"CBDT"), 0), None);
    }

    #[test]
    fn a_malformed_range_stops_at_the_last_character() {
        // A group to u32::MAX: read up to char::MAX and the glyph count, not 4.3e9 characters.
        let Some(chars) = chars_of(&font(3, 10, 12, 0, u32::MAX, b"glyf"), 0) else { panic!("characters") };
        assert!(chars.iter().all(|&(first, last)| first <= last && last <= u32::from(char::MAX)));
        assert!(chars.iter().map(|&(first, last)| last - first + 1).sum::<u32>() <= u32::from(u16::MAX));
    }

    #[test]
    fn a_font_file_is_read_once_for_all_answers() {
        let dir = std::env::temp_dir().join(format!("fc-fonts-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let mut db = fontdb::Database::new();
        let defs = egui::FontDefinitions::default();
        let Some(ubuntu) = defs.font_data.get("Ubuntu-Light") else { panic!("egui's default fonts") };
        let path = dir.join("ubuntu.ttf");
        let _ = std::fs::write(&path, ubuntu.bytes());
        db.load_font_file(&path).unwrap_or_else(|e| panic!("{e}"));
        let p = SystemFallback::over(db);
        for _ in 0..2 {
            assert!(p.font_for(&FallbackRequest { cluster: "Ə", family: &FontFamily::Proportional }).is_some());
        }
        assert_eq!(p.files.lock().unwrap_or_else(PoisonError::into_inner).len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_face_whose_file_is_gone_gives_way_to_the_next() {
        let dir = std::env::temp_dir().join(format!("fc-fonts-gone-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let defs = egui::FontDefinitions::default();
        let Some(ubuntu) = defs.font_data.get("Ubuntu-Light") else { panic!("egui's default fonts") };
        let mut db = fontdb::Database::new();
        // Two copies of one face: equal but for their file, so `a.ttf` ranks first.
        let (a, b) = (dir.join("a.ttf"), dir.join("b.ttf"));
        for path in [&a, &b] {
            let _ = std::fs::write(path, ubuntu.bytes());
            db.load_font_file(path).unwrap_or_else(|e| panic!("{e}"));
        }
        let p = SystemFallback::over(db);
        let first = p.ranked('Ə', false).first().map(|f| f.id);
        let _ = std::fs::remove_file(&a);
        let request = FallbackRequest { cluster: "Ə", family: &FontFamily::Proportional };
        for _ in 0..2 {
            assert!(p.font_for(&request).is_some());
        }
        let found = p.found.lock().unwrap_or_else(PoisonError::into_inner).get(&('Ə', false)).copied().flatten();
        assert!(found.is_some() && found != first);
        assert_eq!(p.files.lock().unwrap_or_else(PoisonError::into_inner).keys().collect::<Vec<_>>(), [&b]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_request_without_fonts_finds_none() {
        let p = SystemFallback::over(fontdb::Database::new());
        assert!(p.font_for(&FallbackRequest { cluster: "한", family: &FontFamily::Proportional }).is_none());
        assert!(SystemFallback::over(defaults()).font_for(&FallbackRequest { cluster: "", family: &FontFamily::Proportional }).is_none());
    }
}
