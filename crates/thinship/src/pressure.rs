//! A pressure distribution on a hull's bottom — the dynamic lift of a hull
//! running trimmed — as a free-wave source in Michell's integral.
//!
//! Thin-ship theory represents a hull by sources on its centreplane: its
//! thickness. A lifting bottom is something else, a pressure pushing down
//! on the water, and its waves are the momentum the lift's induced drag
//! carries off. In linear theory the two simply superpose, so the lift's
//! waves enter Michell's integral as one more amplitude beside the hull's,
//! and the cross term between them — the hull climbing its own bow wave, or
//! not — comes out of the same integral, with a multihull's interference.
//! Its own drag is that integral's: nothing like `L·tan τ` is to be added.
//!
//! **The amplitude.** Forcing the steady free surface with a pressure `p`
//! instead of a source sheet `q` changes only the numerator of its
//! response at the dispersion pole `g|k| = U²k_x²`: a source of flux `q` at
//! depth `h` gives `iUk_x q̂ e^{−|k|h}`, a surface pressure `−|k| p̂/ρ`. So
//! the pressure radiates as a source `q̂ = iλ p̂/(ρU)` (`|k|/k_x = λ` on the
//! pole), and a pressure acting on the bottom at depth `z̄` — a vertical
//! dipole there — carries the source's own `e^{−νλ²z̄}`. A thin hull's
//! sources are `q = 2U ∂f/∂x`, Michell's amplitude `∬ ∂f/∂x …`, and the
//! hull's `x` runs bow-ward against the stream, so in the same convention
//!
//! ```text
//! A_p±(λ) = −(iλ / 2ρU²) ∬ p(x, y) e^{−νλ² z̄} e^{iνλ(x − x_c)} e^{±iνλ√(λ²−1) y} dx dy
//! ```
//!
//! for the ±θ Kelvin half-systems. In the limit of a patch infinitely wide
//! it gives Havelock's two-dimensional `R/B = g |P̂(ν)|² / (ρU⁴)` (a test
//! below).
//!
//! The pressure is carried as **strips**, one per station: the lift on a
//! length of hull, spread evenly over the wetted beam at the section's mean
//! depth.

use hullgeom::moments::C64;

/// The lift on one length of hull: a uniform pressure across its wetted
/// beam.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Strip {
    /// Station, in the hull's frame (bow at high `x`) [m].
    pub x: f64,
    /// The upward force on the hull over the strip's length [N] (the
    /// pressure's integral; negative for suction).
    pub force: f64,
    /// Half the wetted beam it spreads over [m].
    pub half_width: f64,
    /// The depth it acts at, below the waterline [m].
    pub depth: f64,
}

/// A hull's dynamic lift as strips of bottom pressure.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PressurePatch {
    pub strips: Vec<Strip>,
}

impl PressurePatch {
    /// The total upward force [N].
    pub fn force(&self) -> f64 {
        self.strips.iter().map(|s| s.force).sum()
    }

    /// The bow-up moment about `x_ref`, in the hull's frame [N·m].
    pub fn moment_bow_up(&self, x_ref: f64) -> f64 {
        self.strips.iter().map(|s| s.force * (s.x - x_ref)).sum()
    }

    /// The widest strip's half-width [m]: the transverse extent the outer
    /// quadrature has to resolve.
    pub fn half_width(&self) -> f64 {
        self.strips.iter().map(|s| s.half_width).fold(0.0, f64::max)
    }

    /// `A_p±(λ)` in Michell's convention (see the module), phase centre
    /// `x_c`, at `ν = g/U²`; `scale` is `1/(2ρU²)`. The patch is symmetric
    /// about the hull's centreplane, so both half-systems carry the same
    /// amplitude.
    pub fn amplitude(&self, nu: f64, lambda: f64, x_c: f64, scale: f64) -> C64 {
        let kx = nu * lambda;
        let kappa = nu * lambda * lambda;
        let ky = nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
        let mut sum = C64::ZERO;
        for s in &self.strips {
            let decay = (-kappa * s.depth.max(0.0)).exp();
            if decay == 0.0 || s.force == 0.0 {
                continue;
            }
            let a = ky * s.half_width;
            let sinc = if a.abs() < 1e-8 { 1.0 } else { a.sin() / a };
            sum = sum + C64::cis(kx * (s.x - x_c)).scale(s.force * decay * sinc);
        }
        // −iλ·scale·sum
        C64::new(sum.im, -sum.re).scale(lambda * scale)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::michell::{run_outer, MemberWave, OuterParams, WaveOptions};
    use std::f64::consts::PI;

    /// A patch on its own, as a wave-resistance member.
    #[derive(Clone)]
    struct Alone<'p> {
        patch: &'p PressurePatch,
        scale: f64,
    }

    impl MemberWave for Alone<'_> {
        fn amps(&mut self, nu: f64, lambda: f64) -> (C64, C64) {
            let a = self.patch.amplitude(nu, lambda, 0.0, self.scale);
            (a, a)
        }
        fn offsets(&self) -> (f64, f64) {
            (0.0, 0.0)
        }
    }

    fn resistance(patch: &PressurePatch, rho: f64, g: f64, u: f64) -> f64 {
        let nu = g / (u * u);
        let x_half = patch.strips.iter().map(|s| s.x.abs()).fold(0.0, f64::max);
        let params = OuterParams {
            nu,
            x_half,
            y_half: patch.half_width(),
            t_max: patch.strips.iter().map(|s| s.depth).fold(0.0, f64::max),
            floor: 0.0,
        };
        let coeff = 4.0 * rho * g * g / (PI * u * u);
        let member = Alone {
            patch,
            scale: 1.0 / (2.0 * rho * u * u),
        };
        run_outer(&params, &WaveOptions::default(), coeff, vec![member]).resistance
    }

    /// A uniform pressure `p₀` over `|x| < a`, very wide: per unit width it
    /// is Havelock's two-dimensional patch, `R/B = g|P̂(ν)|²/(ρU⁴)` with
    /// `P̂(ν) = 2p₀ sin(νa)/ν` — the normalisation of the amplitude, checked
    /// against a result that owes nothing to Michell's integral.
    #[test]
    fn a_wide_patch_is_havelocks_two_dimensional_one() {
        let (rho, g, p0, a, w) = (1025.0, 9.81f64, 1000.0, 1.0, 100.0);
        let u = g.sqrt(); // ν = 1
        let n = 400;
        let dx = 2.0 * a / n as f64;
        let patch = PressurePatch {
            strips: (0..n)
                .map(|i| Strip {
                    x: -a + (i as f64 + 0.5) * dx,
                    force: p0 * dx * 2.0 * w,
                    half_width: w,
                    depth: 0.0,
                })
                .collect(),
        };
        let nu = g / (u * u);
        let p_hat = 2.0 * p0 * (nu * a).sin() / nu;
        let want = 2.0 * w * g * p_hat * p_hat / (rho * u.powi(4));
        let got = resistance(&patch, rho, g, u);
        eprintln!("wide patch: {got:.1} N, two-dimensional {want:.1} N");
        assert!((got - want).abs() < 0.02 * want, "{got} vs {want}");
    }

    /// The patch's amplitude is linear in its forces, decays with depth as
    /// a source at that depth would, and is the same for a strip anywhere
    /// up to its phase.
    #[test]
    fn the_amplitude_is_a_dipole_at_depth() {
        let strip = |x, depth| PressurePatch {
            strips: vec![Strip {
                x,
                force: 500.0,
                half_width: 0.3,
                depth,
            }],
        };
        let (nu, lam, scale) = (0.4, 1.7, 1e-5);
        let a0 = strip(0.0, 0.0).amplitude(nu, lam, 0.0, scale);
        let a1 = strip(0.0, 0.2).amplitude(nu, lam, 0.0, scale);
        let decay = (-nu * lam * lam * 0.2).exp();
        assert!((a1.abs() - decay * a0.abs()).abs() < 1e-12 * a0.abs());
        let shifted = strip(1.5, 0.0).amplitude(nu, lam, 0.0, scale);
        assert!((shifted.abs() - a0.abs()).abs() < 1e-12 * a0.abs());
        // −iλ·scale·F·sinc at x = x_c: purely negative imaginary.
        assert!(a0.re.abs() < 1e-15 && a0.im < 0.0);
    }
}
