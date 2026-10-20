//! The window thread's end of a channel (#552). An [`Inbox`] is polled each frame: it has no
//! blocking receive and no iterator that waits, so a `for` loop over it (which would wait for
//! the next message forever) does not compile. The crate's clippy configuration disallows
//! `std::sync::mpsc::Receiver` and the channel constructors everywhere else, so the window
//! receives only through here. (The web app has an unrelated `Inbox` of its own, for messages
//! from the page: `apps/fieldcraft-web/src/web.rs`.)
#![allow(clippy::disallowed_types, clippy::disallowed_methods, reason = "the one place the window holds a `Receiver`, behind polling only")]

use std::sync::mpsc::{Receiver, Sender, TryIter, TryRecvError};

/// What the window receives on: [`Inbox::try_recv`] and [`Inbox::try_iter`] only.
///
/// ```
/// let (tx, inbox) = fieldcraft_ui_egui::inbox::channel::<u32>();
/// tx.send(7).unwrap();
/// assert_eq!(inbox.try_iter().collect::<Vec<_>>(), [7]);
/// assert!(inbox.try_recv().is_err());
/// ```
///
/// Looping over it, which would wait for messages, does not compile:
///
/// ```compile_fail
/// let (_tx, inbox) = fieldcraft_ui_egui::inbox::channel::<u32>();
/// for _ in &inbox {}
/// ```
///
/// ```compile_fail
/// let (_tx, inbox) = fieldcraft_ui_egui::inbox::channel::<u32>();
/// let _ = inbox.into_iter();
/// ```
pub struct Inbox<T>(Receiver<T>);

impl<T> Inbox<T> {
    /// The next message if one is waiting.
    pub fn try_recv(&self) -> Result<T, TryRecvError> {
        self.0.try_recv()
    }

    /// The messages waiting now; it ends when there are none, never waiting for more.
    pub fn try_iter(&self) -> TryIter<'_, T> {
        self.0.try_iter()
    }
}

/// A receiver the host made (the control server, the update thread), for the window to poll.
impl<T> From<Receiver<T>> for Inbox<T> {
    fn from(rx: Receiver<T>) -> Self {
        Inbox(rx)
    }
}

/// A channel whose receiving end is an [`Inbox`].
pub fn channel<T>() -> (Sender<T>, Inbox<T>) {
    let (tx, rx) = std::sync::mpsc::channel();
    (tx, Inbox(rx))
}
