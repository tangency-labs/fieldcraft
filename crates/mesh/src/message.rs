//! Messages that name lengths, areas, volumes and points as data, not as text: the meshers work
//! in whatever system their caller expresses the model in, so they cannot name a unit. The caller
//! renders a message with its units at its boundary ([`Message::render`]); shown bare
//! ([`std::fmt::Display`]) the numbers are as the mesher had them, without a unit.

use std::fmt;

use fieldcraft_model::units::{Quantity, format_number};

/// One piece of a [`Message`].
#[derive(Debug, Clone, PartialEq)]
pub enum Part {
    Text(String),
    /// A number of a quantity, in the mesher's system.
    Value(Quantity, f64),
    /// A point, in the mesher's system: `(x, y, z)` (or `(x, y)` in a plane).
    Point(Vec<f64>),
}

impl From<&str> for Part {
    fn from(s: &str) -> Self {
        Part::Text(s.to_string())
    }
}

impl From<String> for Part {
    fn from(s: String) -> Self {
        Part::Text(s)
    }
}

/// A length.
pub fn length(v: f64) -> Part {
    Part::Value(Quantity::Length, v)
}

/// An area.
pub fn area(v: f64) -> Part {
    Part::Value(Quantity::Area, v)
}

/// A volume.
pub fn volume(v: f64) -> Part {
    Part::Value(Quantity::Volume, v)
}

/// A point (three coordinates, or two in a plane).
pub fn point<const N: usize>(p: [f64; N]) -> Part {
    Part::Point(p.to_vec())
}

/// A message whose numbers keep their quantity.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Message(pub Vec<Part>);

/// A [`Message`] from its parts: text (`&str`, `String`) and [`length`], [`area`], [`volume`],
/// [`point`].
#[macro_export]
macro_rules! message {
    ($($part:expr),* $(,)?) => {
        $crate::message::Message(vec![$($crate::message::Part::from($part)),*])
    };
}

/// How a caller shows the numbers of a message: each value converted from the mesher's system,
/// and the symbol of its unit (`"in"`; empty for none).
pub trait Units {
    fn value(&self, q: Quantity, v: f64) -> f64;
    fn symbol(&self, q: Quantity) -> String;
}

/// Numbers as the mesher had them, without units.
pub struct Bare;

impl Units for Bare {
    fn value(&self, _: Quantity, v: f64) -> f64 {
        v
    }
    fn symbol(&self, _: Quantity) -> String {
        String::new()
    }
}

impl Message {
    /// The message as text, its numbers shown in `units`: `2.5 in`, `(1, 0, 2) in`.
    pub fn render(&self, units: &dyn Units) -> String {
        let with = |text: String, q: Quantity| {
            let unit = units.symbol(q);
            if unit.is_empty() { text } else { format!("{text} {unit}") }
        };
        let mut out = String::new();
        for p in &self.0 {
            match p {
                Part::Text(t) => out.push_str(t),
                Part::Value(q, v) => out.push_str(&with(format_number(units.value(*q, *v)), *q)),
                Part::Point(p) => {
                    let xyz: Vec<String> = p.iter().map(|x| format_number(units.value(Quantity::Length, *x))).collect();
                    out.push_str(&with(format!("({})", xyz.join(", ")), Quantity::Length));
                }
            }
        }
        out
    }

    /// `self` followed by `other`.
    pub fn then(mut self, other: Message) -> Message {
        self.0.extend(other.0);
        self
    }

    /// The messages joined by `sep`.
    pub fn join(items: impl IntoIterator<Item = Message>, sep: &str) -> Message {
        let mut out = Message::default();
        for (i, m) in items.into_iter().enumerate() {
            if i > 0 {
                out.0.push(Part::Text(sep.to_string()));
            }
            out.0.extend(m.0);
        }
        out
    }

    /// The message as `units` renders it, as text only.
    pub fn shown(&self, units: &dyn Units) -> Message {
        Message(vec![Part::Text(self.render(units))])
    }
}

impl From<String> for Message {
    fn from(s: String) -> Self {
        Message(vec![Part::Text(s)])
    }
}

impl fmt::Display for Message {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.render(&Bare))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_message_renders_its_numbers_with_their_units() {
        let m = message!["weld distance ", length(2.5), " at ", point([1.0, 0.0, 2.0]), ", volume ", volume(8.0)];
        assert_eq!(m.to_string(), "weld distance 2.5 at (1, 0, 2), volume 8");
        struct Halves;
        impl Units for Halves {
            fn value(&self, _: Quantity, v: f64) -> f64 {
                v / 2.0
            }
            fn symbol(&self, q: Quantity) -> String {
                if q == Quantity::Volume { "u³".into() } else { "u".into() }
            }
        }
        assert_eq!(m.render(&Halves), "weld distance 1.25 u at (0.5, 0, 1) u, volume 4 u³");
    }
}
