//! Lee & Moser's mean profiles of the turbulent channel (#1105, #1108), shared by the CPU and the GPU
//! matched-wall studies.
//!
//! M. Lee, R. D. Moser, *Direct numerical simulation of turbulent channel flow up to Re_τ ≈ 5200*, J. Fluid
//! Mech. 774:395–415, 2015, doi:10.1017/jfm.2015.268. Sampled from the public mean profiles at
//! turbulence.oden.utexas.edu (`channel2015/data/LM_Channel_<Re_τ>_mean_prof.dat`, the header's U_mean = 1 and
//! u_τ; the file's U column is U⁺). Each table is the file's U⁺ at 16 points y/δ, linearly interpolated in y
//! between the file's rows, and clamped to the last row (the file ends at y/δ ≈ 0.997, where the centreline's
//! value is taken). The bulk is U_b⁺ = 1/u_τ.

/// One Lee & Moser mean profile at a friction Reynolds number.
pub struct LeeMoserChannel {
    /// The friction Reynolds number of the profile.
    pub re_tau: u32,
    /// The bulk velocity U_b⁺ = 1/u_τ (the file's U_mean = 1).
    pub bulk: f64,
    /// U⁺ at y/δ.
    pub profile: [(f64, f64); 16],
}

/// Re_τ 1000: u_τ = 5.00256 × 10⁻², the file's bulk is 19.990.
pub const LM_1000_CHANNEL: LeeMoserChannel = LeeMoserChannel {
    re_tau: 1000,
    bulk: 1.0 / 5.00256e-2,
    profile: [
        (0.01, 8.437),
        (0.02, 11.968),
        (0.03, 13.455),
        (0.05, 14.839),
        (0.07, 15.624),
        (0.1, 16.456),
        (0.15, 17.450),
        (0.2, 18.200),
        (0.3, 19.320),
        (0.4, 20.168),
        (0.5, 20.858),
        (0.6, 21.434),
        (0.7, 21.909),
        (0.8, 22.272),
        (0.9, 22.509),
        (1.0, 22.593),
    ],
};

/// Re_τ 2000: u_τ = 4.58794 × 10⁻², the file's bulk is 21.796.
pub const LM_2000_CHANNEL: LeeMoserChannel = LeeMoserChannel {
    re_tau: 2000,
    bulk: 1.0 / 4.58794e-2,
    profile: [
        (0.01, 11.924),
        (0.02, 14.241),
        (0.03, 15.233),
        (0.05, 16.417),
        (0.07, 17.228),
        (0.1, 18.129),
        (0.15, 19.195),
        (0.2, 19.976),
        (0.3, 21.121),
        (0.4, 21.971),
        (0.5, 22.659),
        (0.6, 23.237),
        (0.7, 23.716),
        (0.8, 24.081),
        (0.9, 24.315),
        (1.0, 24.396),
    ],
};

/// Re_τ 5200: u_τ = 4.14872 × 10⁻², the file's bulk is 24.102.
pub const LM_5200_CHANNEL: LeeMoserChannel = LeeMoserChannel {
    re_tau: 5200,
    bulk: 1.0 / 4.14872e-2,
    profile: [
        (0.01, 14.882),
        (0.02, 16.500),
        (0.03, 17.481),
        (0.05, 18.771),
        (0.07, 19.641),
        (0.1, 20.569),
        (0.15, 21.628),
        (0.2, 22.384),
        (0.3, 23.479),
        (0.4, 24.291),
        (0.5, 24.946),
        (0.6, 25.492),
        (0.7, 25.942),
        (0.8, 26.283),
        (0.9, 26.500),
        (1.0, 26.575),
    ],
};

/// The profile at `re_tau`, when the studies have one.
pub fn at_re_tau(re_tau: u32) -> Option<&'static LeeMoserChannel> {
    [&LM_1000_CHANNEL, &LM_2000_CHANNEL, &LM_5200_CHANNEL].into_iter().find(|lm| lm.re_tau == re_tau)
}

impl LeeMoserChannel {
    /// U⁺ at `y` (over δ, in (0, 1]), interpolated linearly in ln y.
    pub fn at(&self, y: f64) -> f64 {
        let p = &self.profile;
        let i = p.windows(2).position(|w| y <= w[1].0).unwrap_or(p.len() - 2);
        let ((y0, u0), (y1, u1)) = (p[i], p[i + 1]);
        let t = (y.max(y0).ln() - y0.ln()) / (y1.ln() - y0.ln());
        u0 + t.clamp(0.0, 1.0) * (u1 - u0)
    }
}
