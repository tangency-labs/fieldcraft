//! Fan/porous data appended to the existing cell-flags buffer: no additional shader binding.
//! The first `n` words remain flags (0 solid, 1 unzoned fluid, 2+zone owned fluid). Then each
//! zone has 16 words, followed by its curve points as pairs of f32 bits. The shader changes
//! only the fan's current rise (word 12). All offsets are absolute word indices.

use fieldcraft_lbm::{Error, Flow, Result, ZoneKind};

pub(super) const WORDS: usize = fieldcraft_lbm::memory::DEVICE_ZONE_WORDS;

fn index(n: usize) -> Result<u32> {
    u32::try_from(n).map_err(|_| Error::Limit("too much fan/porous data for the GPU".into()))
}

fn value(name: &str, label: &str, v: f64) -> Result<u32> {
    let stored = v as f32;
    if !stored.is_finite() || (v != 0.0 && stored == 0.0) {
        return Err(Error::Invalid(format!("zone `{name}`: {label} cannot be represented in GPU single precision")));
    }
    Ok(stored.to_bits())
}

/// Add the zones to the flags prepared by the caller, preserving CPU cell ownership and any
/// fan rise already reached by a warmed-up CPU flow. Refuse bad f32 curves before device work.
pub(super) fn append(flow: &Flow, flags: &mut Vec<u32>) -> Result<()> {
    let zones = flow.zones();
    if zones.is_empty() {
        return Ok(());
    }
    let n = flow.grid().len();
    if flags.len() != n {
        return Err(Error::Invalid("GPU zone upload needs one flag per cell".into()));
    }
    let headers = zones.len().checked_mul(WORDS).ok_or_else(|| Error::Limit("too many GPU zones".into()))?;
    let curves = zones
        .iter()
        .try_fold(0usize, |len, z| {
            let points = match &z.kind {
                ZoneKind::Fan { curve, .. } => curve.len(),
                _ => 0,
            };
            len.checked_add(points.checked_mul(2)?)
        })
        .ok_or_else(|| Error::Limit("too many GPU fan curve points".into()))?;
    let extra = headers.checked_add(curves).ok_or_else(|| Error::Limit("too much GPU zone data".into()))?;
    index(n.checked_add(extra).ok_or_else(|| Error::Limit("too much GPU zone data".into()))?)?;
    flags.try_reserve_exact(extra).map_err(|_| Error::Limit("not enough memory for GPU zones".into()))?;
    flags.resize(n + headers, 0);
    let samples = flow.zone_samples();
    for (k, zone) in zones.iter().enumerate() {
        let mut data = [0; WORDS];
        data[1] = index(zone.axis)?;
        for axis in 0..3 {
            data[2 + axis] = index(zone.lo[axis])?;
            data[5 + axis] = index(zone.hi[axis])?;
        }
        data[10] = 1.0f32.to_bits();
        data[15] = index(zone.layers())?;
        match &zone.kind {
            ZoneKind::Porous { alpha, beta } => {
                data[8] = value(&zone.name, "Darcy resistance", *alpha)?;
                data[9] = value(&zone.name, "inertial resistance", *beta)?;
            }
            ZoneKind::Fan { sign, curve, relaxation } => {
                data[0] = 1;
                data[10] = value(&zone.name, "direction", *sign)?;
                data[11] = value(&zone.name, "relaxation", *relaxation)?;
                data[12] = value(&zone.name, "current pressure rise", samples.get(k).map_or(0.0, |s| s.fan_rise))?;
                data[13] = index(flags.len())?;
                data[14] = index(curve.len())?;
                let mut previous = None;
                for point in curve {
                    let x = value(&zone.name, "curve flow", point[0])?;
                    let y = value(&zone.name, "curve rise", point[1])?;
                    let flow = f32::from_bits(x);
                    if previous.is_some_and(|p| flow <= p) {
                        return Err(Error::Invalid(format!("zone `{}`: fan curve flows are indistinguishable in GPU single precision", zone.name)));
                    }
                    previous = Some(flow);
                    flags.extend([x, y]);
                }
            }
        }
        let start = n + k * WORDS;
        flags.get_mut(start..start + WORDS).ok_or_else(|| Error::Invalid("missing GPU zone header".into()))?.copy_from_slice(&data);
        let cells = flow.zone_cells().get(k).ok_or_else(|| Error::Invalid(format!("zone `{}` has no prepared cells", zone.name)))?;
        for cell in cells {
            let flag = flags
                .get_mut(*cell as usize)
                .filter(|_| (*cell as usize) < n)
                .ok_or_else(|| Error::Invalid(format!("zone `{}` has a cell outside the grid", zone.name)))?;
            *flag = index(k + 2)?;
        }
    }
    Ok(())
}
