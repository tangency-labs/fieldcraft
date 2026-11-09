//! Linear superposition of reference-load static stress tensors. Same model,
//! recovery location and homogeneous boundary operator are required by the caller.
//! Coefficients multiply the complete reference pattern; they are dimensionless.
//! Rainflow and Miner are applied AFTER tensor superposition and fixed projection.
use super::{Basquin, MAX_HISTORY_POINTS, life::StressMean, miner, rainflow};
use crate::{Error, Result, StaticResult};

#[derive(Debug, serde::Serialize)]
pub struct NodeDamage {
    pub node: u32,
    pub damage: f64,
    pub cycle_count: f64,
    pub duration_seconds: f64,
    pub repeated_record_life_seconds: Option<f64>,
}

pub fn evaluate(
    results: &[&StaticResult],
    coefficients: &[Vec<f64>],
    times: &[f64],
    nodes: &[u32],
    weights: [f64; 6],
    curve: &Basquin,
    correction: StressMean,
) -> Result<Vec<NodeDamage>> {
    curve.validate()?;
    if results.is_empty()
        || results.len() > 64
        || coefficients.len() != results.len()
        || times.len() < 2
        || times.len() > MAX_HISTORY_POINTS
        || times.iter().any(|t| !t.is_finite())
        || times.windows(2).any(|w| w[1] <= w[0])
        || coefficients.iter().any(|c| c.len() != times.len() || c.iter().any(|v| !v.is_finite()))
        || nodes.is_empty()
        || nodes.len() > 100000
        || nodes.iter().collect::<std::collections::BTreeSet<_>>().len() != nodes.len()
        || weights.iter().any(|x| !x.is_finite())
        || weights.iter().all(|x| *x == 0.0)
        || nodes.len().checked_mul(times.len()).and_then(|v| v.checked_mul(results.len())).is_none_or(|n| n > 20_000_000)
    {
        return Err(Error::Invalid("fatigue superposition dimensions, times, projection or work cap invalid".into()));
    }
    let first = results.first().ok_or_else(|| Error::Invalid("empty reference results".into()))?;
    if results.iter().any(|r| {
        r.nodes != first.nodes || r.stress.len() != r.nodes.len() || r.nonlinear.is_some() || r.temperature.is_some() || !r.contacts.is_empty()
    }) {
        return Err(Error::Invalid(
            "fatigue references must be same-node linear mechanical static results without thermal/contact base state".into(),
        ));
    }
    let duration = times.last().copied().unwrap_or(0.0) - times.first().copied().unwrap_or(0.0);
    if !duration.is_finite() {
        return Err(Error::Invalid("history duration overflow".into()));
    }
    let mut out = Vec::with_capacity(nodes.len());
    for &node in nodes {
        let index = first.nodes.binary_search(&node).map_err(|_| Error::Invalid(format!("node {node} absent from reference stress")))?;
        let fields =
            results.iter().map(|r| r.stress.get(index).ok_or_else(|| Error::Invalid("missing stress".into()))).collect::<Result<Vec<_>>>()?;
        let mut signal = Vec::with_capacity(times.len());
        for t in 0..times.len() {
            let mut tensor = [0.0; 6];
            for (field, history) in fields.iter().zip(coefficients) {
                let factor = history.get(t).copied().ok_or_else(|| Error::Invalid("missing load factor".into()))?;
                for (dest, &value) in tensor.iter_mut().zip(field.iter()) {
                    *dest += value * factor;
                }
            }
            let value: f64 = tensor.iter().zip(weights).map(|(x, w)| x * w).sum();
            if tensor.iter().any(|x| !x.is_finite()) || !value.is_finite() {
                return Err(Error::Invalid("superposed stress overflow".into()));
            }
            signal.push(value);
        }
        let cycles = rainflow(&signal)?;
        let damage = miner(&cycles, curve, correction)?;
        let life = duration / damage;
        out.push(NodeDamage {
            node,
            damage,
            cycle_count: cycles.iter().map(|c| c.count).sum(),
            duration_seconds: duration,
            repeated_record_life_seconds: if damage > 0.0 && life.is_finite() { Some(life) } else { None },
        });
    }
    Ok(out)
}
