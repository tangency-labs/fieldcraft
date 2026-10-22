//! A caption bar under a picture: a line of text on the left, and on the right the contour
//! legend (its lowest value, a swatch per band, its highest value) and a detail such as the
//! frame's time. Drawn into RGBA pixels in the bitmap font ([`crate::font::cased`]), so a picture
//! drawn on the GPU or the CPU gets the same bar (the clips of `results.export.animation`).

use crate::Rgb;
use crate::font;

/// The bar's colours.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Colors {
    pub background: Rgb,
    /// The caption.
    pub text: Rgb,
    /// The legend's values and the detail.
    pub dim: Rgb,
}

/// The legend: its ends as shown (with their unit) and one colour per band, low to high.
#[derive(Clone, Debug, PartialEq)]
pub struct Legend {
    pub low: String,
    pub high: String,
    pub bands: Vec<Rgb>,
}

/// What a bar shows; nothing to show draws no bar.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Bar {
    pub caption: Option<String>,
    pub legend: Option<Legend>,
    pub detail: Option<String>,
}

impl Bar {
    pub fn is_empty(&self) -> bool {
        self.caption.is_none() && self.legend.is_none() && self.detail.is_none()
    }
}

/// Most bands a legend draws.
const MAX_BANDS: usize = 64;
/// Most characters of one text.
const MAX_CHARS: usize = 200;

/// The bar's height under a picture `h` pixels high: a twelfth of it, from 24 to 120 pixels, even
/// (so a clip of an even picture stays even).
pub fn height(h: usize) -> usize {
    ((h / 12).clamp(24, 120) / 2) * 2
}

/// The pixel size of the font in a bar `bar_h` high: its text about a third of the bar.
fn scale(bar_h: usize) -> usize {
    (bar_h / (3 * font::H)).max(1)
}

/// The width the right side of `bar` takes (the legend and the detail) at font pixel size `px`.
fn right_width(bar: &Bar, px: usize) -> usize {
    let gap = px * font::ADVANCE;
    let detail = bar.detail.as_ref().map_or(0, |d| text_width(d, px) + 2 * gap);
    let legend = bar
        .legend
        .as_ref()
        .map_or(0, |l| text_width(&l.low, px) + text_width(&l.high, px) + 2 * gap + swatch(px) * l.bands.len().min(MAX_BANDS) + 2 * gap);
    detail + legend
}

/// The width of one band's swatch.
fn swatch(px: usize) -> usize {
    font::ADVANCE * px
}

/// The width of `text` in pixels at font pixel size `px`.
fn text_width(text: &str, px: usize) -> usize {
    text.chars().take(MAX_CHARS).count() * font::ADVANCE * px
}

/// RGBA pixels drawn into, `w` wide and `h` high; what falls outside is clipped.
struct Pixels<'a> {
    rgba: &'a mut [u8],
    w: usize,
    h: usize,
}

impl Pixels<'_> {
    /// Fill a `bw` × `bh` box at (`x`, `y`).
    fn fill(&mut self, (x, y): (usize, usize), (bw, bh): (usize, usize), color: Rgb) {
        for yy in y..(y + bh).min(self.h) {
            for xx in x..(x + bw).min(self.w) {
                let i = (yy * self.w + xx) * 4;
                if let Some(p) = self.rgba.get_mut(i..i + 4) {
                    p.copy_from_slice(&[color.0, color.1, color.2, 255]);
                }
            }
        }
    }

    /// `text` with its top left at `at` in `color`, font pixels `px` wide.
    fn text(&mut self, at: (usize, usize), px: usize, text: &str, color: Rgb) {
        for (i, ch) in text.chars().take(MAX_CHARS).enumerate() {
            for (fy, row) in font::cased(ch).iter().enumerate() {
                for fx in (0..font::W).filter(|fx| row & (0x10 >> fx) != 0) {
                    self.fill((at.0 + (i * font::ADVANCE + fx) * px, at.1 + fy * px), (px, px), color);
                }
            }
        }
    }
}

/// `rgba` (`w` × `h` pixels) with `bar` drawn under it (in place, grown): the new pixels and
/// height; the picture as it is when the bar is empty or the sizes disagree.
pub fn under(rgba: Vec<u8>, w: usize, h: usize, bar: &Bar, colors: &Colors) -> (Vec<u8>, usize) {
    if bar.is_empty() || w.checked_mul(h).and_then(|n| n.checked_mul(4)) != Some(rgba.len()) || w == 0 {
        return (rgba, h);
    }
    let bh = height(h);
    let total = h + bh;
    let mut out = rgba;
    out.resize(w * total * 4, 0);
    let mut p = Pixels { rgba: &mut out, w, h: total };
    p.fill((0, h), (w, bh), colors.background);
    // The font as large as the bar allows while the legend and the detail keep to three quarters
    // of the width.
    let px = (1..=scale(bh)).rev().find(|px| right_width(bar, *px) <= w * 3 / 4).unwrap_or(1);
    let y = h + (bh.saturating_sub(font::H * px)) / 2;
    let margin = 2 * px * font::ADVANCE;
    // Right to left: the detail, the legend's high end, its bands, its low end.
    let mut right = w.saturating_sub(margin);
    if let Some(d) = &bar.detail {
        right = right.saturating_sub(text_width(d, px));
        p.text((right, y), px, d, colors.dim);
        right = right.saturating_sub(margin);
    }
    if let Some(l) = &bar.legend {
        let bands = l.bands.get(..l.bands.len().min(MAX_BANDS)).unwrap_or_default();
        let (sw, sh) = (swatch(px), font::H * px);
        right = right.saturating_sub(text_width(&l.high, px));
        p.text((right, y), px, &l.high, colors.dim);
        right = right.saturating_sub(px * font::ADVANCE + sw * bands.len());
        for (k, c) in bands.iter().enumerate() {
            p.fill((right + k * sw, y), (sw, sh), *c);
        }
        right = right.saturating_sub(px * font::ADVANCE + text_width(&l.low, px));
        p.text((right, y), px, &l.low, colors.dim);
        right = right.saturating_sub(margin);
    }
    if let Some(c) = &bar.caption {
        // The caption takes what is left on the left, cut short to fit.
        let room = right.saturating_sub(margin) / (font::ADVANCE * px);
        let shown: String = c.chars().take(room).collect();
        p.text((margin, y), px, &shown, colors.text);
    }
    (out, total)
}

#[cfg(test)]
mod tests {
    use super::*;

    const COLORS: Colors = Colors { background: Rgb(20, 22, 26), text: Rgb(240, 240, 240), dim: Rgb(180, 180, 190) };

    fn count(rgba: &[u8], c: Rgb) -> usize {
        rgba.as_chunks::<4>().0.iter().filter(|p| p[..3] == [c.0, c.1, c.2]).count()
    }

    #[test]
    fn a_bar_goes_under_the_picture_with_its_texts_and_bands() {
        let (w, h) = (320, 180);
        let picture = vec![255u8; w * h * 4];
        let bar = Bar {
            caption: Some("AIO radiator · fan air".into()),
            legend: Some(Legend { low: "25 °C".into(), high: "45 °C".into(), bands: vec![Rgb(1, 2, 200), Rgb(200, 2, 1)] }),
            detail: Some("t = 0.25 s".into()),
        };
        let (out, total) = under(picture.clone(), w, h, &bar, &COLORS);
        assert_eq!(total, h + height(h));
        assert_eq!(total % 2, 0);
        assert_eq!(out.len(), w * total * 4);
        assert_eq!(&out[..w * h * 4], &picture[..], "the picture is untouched");
        for c in [COLORS.background, COLORS.text, COLORS.dim, Rgb(1, 2, 200), Rgb(200, 2, 1)] {
            assert!(count(&out[w * h * 4..], c) > 10, "{c:?}");
        }
        // Nothing to show, or pixels of another size: the picture as it is.
        assert_eq!(under(picture.clone(), w, h, &Bar::default(), &COLORS), (picture.clone(), h));
        assert_eq!(under(picture[..8].to_vec(), w, h, &bar, &COLORS).1, h);
        // A long caption and a narrow picture never draw outside.
        let long = Bar { caption: Some("x".repeat(10_000)), ..bar };
        let (out, total) = under(vec![0u8; 16 * 16 * 4], 16, 16, &long, &COLORS);
        assert_eq!(out.len(), 16 * total * 4);
    }
}
