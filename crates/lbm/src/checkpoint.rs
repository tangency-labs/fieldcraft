//! CPU restart state, streamed without cloning population or output arrays.
//!
//! The host owns file confinement, atomic replacement and the identity of the original
//! model/settings. The solver owns this versioned binary state and its shape checks. Numbers
//! retain their IEEE bits; signed populations are valid and are never required to be positive.

use std::io::{Read, Write};

use crate::{Error, Result};

/// An explicit checkpoint destination. A host must validate identity before invoking `restore`
/// and publish a save atomically only after `write` succeeds. No filesystem is used by LBM.
pub trait Store {
    fn every(&self) -> u64;
    /// Conservative solver-state bound checked before loading output/history arrays.
    fn checked_memory(&mut self, _bytes: u64) {}
    fn load(&mut self, restore: &mut dyn FnMut(&mut dyn Read) -> Result<u64>) -> Result<bool>;
    fn save(&mut self, step: u64, write: &mut dyn FnMut(&mut dyn Write) -> Result<()>) -> Result<()>;
}

#[doc(hidden)]
pub trait Archive {
    fn loading(&self) -> bool;
    fn bytes(&mut self, value: &mut [u8]) -> Result<()>;
}

pub(crate) struct Reader<'a>(pub &'a mut dyn Read);
pub(crate) struct Writer<'a>(pub &'a mut dyn Write);

impl Archive for Reader<'_> {
    fn loading(&self) -> bool {
        true
    }
    fn bytes(&mut self, value: &mut [u8]) -> Result<()> {
        self.0.read_exact(value).map_err(|e| Error::Invalid(format!("checkpoint read: {e}")))
    }
}
impl Archive for Writer<'_> {
    fn loading(&self) -> bool {
        false
    }
    fn bytes(&mut self, value: &mut [u8]) -> Result<()> {
        self.0.write_all(value).map_err(|e| Error::Invalid(format!("checkpoint write: {e}")))
    }
}

pub(crate) fn u64(a: &mut dyn Archive, value: &mut u64) -> Result<()> {
    let mut bytes = value.to_le_bytes();
    a.bytes(&mut bytes)?;
    *value = u64::from_le_bytes(bytes);
    Ok(())
}
pub(crate) fn boolean(a: &mut dyn Archive, value: &mut bool) -> Result<()> {
    let mut bytes = [u8::from(*value)];
    a.bytes(&mut bytes)?;
    *value = match bytes {
        [0] => false,
        [1] => true,
        _ => return Err(Error::Invalid("checkpoint boolean is invalid".into())),
    };
    Ok(())
}
pub(crate) fn f64(a: &mut dyn Archive, value: &mut f64) -> Result<()> {
    let mut bits = value.to_bits();
    u64(a, &mut bits)?;
    let read = f64::from_bits(bits);
    if !read.is_finite() {
        return Err(Error::Invalid("checkpoint contains a non-finite number".into()));
    }
    *value = read;
    Ok(())
}
pub(crate) fn fixed(a: &mut dyn Archive, value: u64) -> Result<()> {
    let mut read = value;
    u64(a, &mut read)?;
    if read != value {
        return Err(Error::Invalid("checkpoint topology or settings differ from this flow".into()));
    }
    Ok(())
}
pub(crate) fn count(a: &mut dyn Archive, value: usize, max: usize) -> Result<usize> {
    let mut read = value as u64;
    u64(a, &mut read)?;
    usize::try_from(read).ok().filter(|n| *n <= max).ok_or_else(|| Error::Limit("checkpoint array exceeds its expected bound".into()))
}
pub(crate) fn resize<T: Default>(values: &mut Vec<T>, len: usize) -> Result<()> {
    if len > values.len() {
        values.try_reserve_exact(len - values.len()).map_err(|_| Error::Limit("not enough memory for checkpoint state".into()))?;
    }
    values.resize_with(len, T::default);
    Ok(())
}
pub(crate) fn doubles(a: &mut dyn Archive, values: &mut Vec<f64>, expected: usize) -> Result<()> {
    fixed(a, expected as u64)?;
    if a.loading() {
        resize(values, expected)?;
    }
    if values.len() != expected {
        return Err(Error::Invalid("checkpoint state array has the wrong length".into()));
    }
    for value in values {
        f64(a, value)?;
    }
    Ok(())
}
pub(crate) fn vectors<const N: usize>(a: &mut dyn Archive, values: &mut Vec<[f64; N]>, expected: usize) -> Result<()>
where
    [f64; N]: Default,
{
    fixed(a, expected as u64)?;
    if a.loading() {
        resize(values, expected)?;
    }
    if values.len() != expected {
        return Err(Error::Invalid("checkpoint state array has the wrong length".into()));
    }
    for value in values {
        for v in value {
            f64(a, v)?;
        }
    }
    Ok(())
}

pub(crate) fn singles(a: &mut dyn Archive, values: &mut Vec<f32>, expected: usize) -> Result<()> {
    fixed(a, expected as u64)?;
    if a.loading() {
        resize(values, expected)?;
    }
    if values.len() != expected {
        return Err(Error::Invalid("checkpoint field has the wrong length".into()));
    }
    for v in values {
        let mut bytes = v.to_le_bytes();
        a.bytes(&mut bytes)?;
        *v = f32::from_le_bytes(bytes);
        if !v.is_finite() {
            return Err(Error::Invalid("checkpoint field is non-finite".into()));
        }
    }
    Ok(())
}

/// The time-loop state is borrowed: saving it never copies the saved frames or histories.
pub(crate) struct RunState<'a> {
    pub k: &'a mut u64,
    pub next_frame: &'a mut u64,
    pub finished: &'a mut bool,
    pub frames: &'a mut Vec<crate::solve::Frame>,
    pub times: &'a mut Vec<f64>,
    pub forces: &'a mut Vec<Vec<[f64; 3]>>,
    pub heats: &'a mut Vec<Vec<f64>>,
    pub zones: &'a mut Vec<crate::solve::ZoneHistory>,
    pub monitor: &'a mut Option<crate::steady::Monitor>,
}

impl RunState<'_> {
    pub fn archive(&mut self, a: &mut dyn Archive, setup: &crate::Setup) -> Result<()> {
        fixed(a, 0x315453524643)?; // FCRST1; changing any state field changes this version.
        u64(a, self.k)?;
        if *self.k == 0 || *self.k > setup.steps {
            return Err(Error::Invalid("checkpoint step is outside this run".into()));
        }
        u64(a, self.next_frame)?;
        boolean(a, self.finished)?;
        let total = setup.steps;
        let planned = u64::from(setup.frames).min(total).max(1);
        let frame_step = |j: u64| (j * total + planned / 2) / planned;
        let scheduled = (1..=planned).take_while(|j| frame_step(*j) <= *self.k).count();
        if *self.next_frame != scheduled as u64 + 1 {
            return Err(Error::Invalid("checkpoint frame cursor differs from its step".into()));
        }
        let final_extra = *self.finished && *self.k != frame_step(scheduled as u64);
        let expected_frames = 1 + scheduled + usize::from(final_extra);
        fixed(a, expected_frames as u64)?;
        if a.loading() {
            self.frames.clear();
        }
        resize(self.frames, expected_frames)?;
        for (i, f) in self.frames.iter_mut().enumerate() {
            let step = if final_extra && i + 1 == expected_frames { *self.k } else { frame_step(i as u64) };
            fixed(a, step)?;
            f.step = step;
            fixed(a, (step as f64 * setup.time_step).to_bits())?;
            f.time = step as f64 * setup.time_step;
            fixed(a, setup.grid.len() as u64)?;
            if a.loading() {
                resize(&mut f.velocity, setup.grid.len())?;
            }
            if f.velocity.len() != setup.grid.len() {
                return Err(Error::Invalid("checkpoint velocity field has the wrong length".into()));
            }
            for v in &mut f.velocity {
                for c in v {
                    let mut b = c.to_le_bytes();
                    a.bytes(&mut b)?;
                    *c = f32::from_le_bytes(b);
                    if !c.is_finite() {
                        return Err(Error::Invalid("checkpoint velocity is non-finite".into()));
                    }
                }
            }
            singles(a, &mut f.pressure, setup.grid.len())?;
            singles(a, &mut f.temperature, if setup.thermal.is_some() { setup.grid.len() } else { 0 })?;
            f64(a, &mut f.mass)?;
        }
        let every = total.div_ceil(crate::solve::MAX_HISTORY).max(1);
        let samples = (*self.k / every + u64::from(*self.k == total && !self.k.is_multiple_of(every))) as usize;
        doubles(a, self.times, samples)?;
        for (i, time) in self.times.iter().enumerate() {
            let step = (((i as u64) + 1) * every).min(total);
            if time.to_bits() != (step as f64 * setup.time_step).to_bits() {
                return Err(Error::Invalid("checkpoint history time differs from its step".into()));
            }
        }
        fixed(a, self.forces.len() as u64)?;
        for h in self.forces.iter_mut() {
            vectors(a, h, samples)?;
        }
        fixed(a, self.heats.len() as u64)?;
        for h in self.heats.iter_mut() {
            doubles(a, h, samples)?;
        }
        fixed(a, self.zones.len() as u64)?;
        for h in self.zones.iter_mut() {
            doubles(a, &mut h.flow, samples)?;
            doubles(a, &mut h.pressure, samples)?;
            doubles(a, &mut h.fan_rise, if h.fan { samples } else { 0 })?;
        }
        fixed(a, u64::from(self.monitor.is_some()))?;
        if let Some(m) = self.monitor.as_mut() {
            m.checkpoint(a, setup, *self.k, *self.finished)?;
        } else if *self.finished {
            return Err(Error::Invalid("checkpoint has an early stop without a steady rule".into()));
        }
        Ok(())
    }
}
