//! The two-dimensional, deep-water, frequency-domain **pulsating source**.
//!
//! Coordinates are the section plane: `y` transverse, `z` **up**, fluid in
//! `z < 0`; time dependence `e^{−iωt}`, `ν = ω²/g`. The potential of a unit
//! source at `(η, ζ)` satisfying `∇²G = 2πδ`, the free-surface condition
//! `G_z − νG = 0` on `z = 0`, decay at depth and outgoing waves is
//!
//! ```text
//! G = ln r − ln r₁ − 2 PV∫₀^∞ e^{k(z+ζ)} cos k(y−η) / (k − ν) dk − 2πi e^{ν(z+ζ)} cos ν(y−η)
//! ```
//!
//! with `r₁` the distance to the image `(η, −ζ)`. The principal value is
//! closed form in the exponential integral,
//! `PV∫ = Re P(w)`, `P(w) = e^{w}(E₁(w) + iπ)`, `w = ν(z + ζ + i|y − η|)`
//! (principal branch; the `iπ` makes `P` real and continuous as `w` reaches
//! the negative real axis from above). Far away `G → −2πi e^{ν(z+ζ)}
//! e^{iν|y−η|}`: outgoing.
//!
//! For panel integration `G` is split as `ln r + ln r₁ + W`: the two
//! logarithms integrate in closed form over a straight panel, and
//!
//! ```text
//! W = −2 Re[P(w) + ln w] + 2 ln ν − 2πi e^{νZ} cos νY,     Z = z + ζ, Y = y − η
//! ```
//!
//! is bounded — `P + ln w` is continuous at `w = 0` — so Gauss quadrature
//! handles it even where a source panel meets the waterline.

use hullgeom::C64;
use std::f64::consts::{EULER_GAMMA, PI};

/// Beyond this `|w|` the continued fraction is used; below it the power
/// series (whose terms peak near `e^{|w|}`, so at most ~5 digits are lost).
const SERIES_LIMIT: f64 = 10.0;

/// Beyond this `|w|` the asymptotic series, truncated at its smallest term
/// (error ~`e^{−|w|}`), replaces the continued fraction, which converges
/// slowly near the negative real axis — where short waves put the sources
/// of a deep section.
const ASYMPTOTIC_LIMIT: f64 = 40.0;

/// `(P(w), P(w) + ln w)` for `Im w ≥ 0`, `Re w ≤ 0` (the source's own
/// quadrant; `w = 0` gives `(∞, −γ + iπ)` — only the second is finite).
fn p_and_regular(w: C64) -> (C64, C64) {
    let ipi = C64::new(0.0, PI);
    let r = w.abs();
    if r == 0.0 {
        return (C64::new(f64::INFINITY, 0.0), C64::new(-EULER_GAMMA, PI));
    }
    let ew = w.exp();
    // The series' cancellation error, ~1e-16·e^{|w|} on E₁, reaches P
    // multiplied by |e^w| = e^{Re w}: harmless wherever w leans toward the
    // negative real axis — exactly where the continued fraction crawls.
    let series = r < SERIES_LIMIT || (r < ASYMPTOTIC_LIMIT && w.re < -0.7 * r);
    if series {
        // E₁(w) = −γ − ln w − S(w),  S(w) = Σ_{n≥1} (−w)ⁿ/(n·n!).
        let mut term = C64::ONE;
        let mut s = C64::ZERO;
        for n in 1..400 {
            term = term * (-w).scale(1.0 / n as f64);
            let add = term.scale(1.0 / n as f64);
            s = s + add;
            if add.abs() < 1e-17 * s.abs().max(1e-300) {
                break;
            }
        }
        let ln_w = w.ln();
        // P + ln w = e^w(−γ − S + iπ) + ln w·(1 − e^w): no cancellation at w → 0.
        let reg = ew * (C64::new(-EULER_GAMMA, PI) - s) + ln_w * (C64::ONE - ew);
        (reg - ln_w, reg)
    } else {
        let e = if r > ASYMPTOTIC_LIMIT {
            ew_e1_asymptotic(w)
        } else {
            ew_e1_continued_fraction(w)
        };
        // On the negative real axis the fraction returns the principal
        // value e^w(−Ei(−w)), which is already P there.
        let p = if w.im == 0.0 && w.re < 0.0 {
            e
        } else {
            e + ipi * ew
        };
        (p, p + w.ln())
    }
}

/// `e^{w} E₁(w) ~ Σ (−1)ⁿ n!/w^{n+1}`, summed to its smallest term; valid
/// for `|arg w| < 3π/2`, and on the negative real axis the principal value.
fn ew_e1_asymptotic(w: C64) -> C64 {
    let inv = w.recip();
    let mut term = inv;
    let mut sum = inv;
    let mut last = term.abs();
    for n in 1..200 {
        let next = term * inv.scale(-(n as f64));
        let mag = next.abs();
        if mag >= last || mag < 1e-17 * sum.abs() {
            break;
        }
        sum = sum + next;
        term = next;
        last = mag;
    }
    sum
}

/// `e^{w} E₁(w)` by the even continued fraction
/// `1/(w + 1 − 1²/(w + 3 − 2²/(w + 5 − …)))`, modified Lentz.
fn ew_e1_continued_fraction(w: C64) -> C64 {
    let tiny = 1e-300;
    let mut b = w + C64::ONE;
    let mut c = C64::new(1.0 / tiny, 0.0);
    let mut d = b.recip();
    let mut h = d;
    for i in 1..5000 {
        let an = -((i * i) as f64);
        b = b + C64::new(2.0, 0.0);
        d = (d.scale(an) + b).recip();
        c = b + C64::new(an, 0.0) / c;
        if c.abs() < tiny {
            c = C64::new(tiny, 0.0);
        }
        let del = c * d;
        h = h * del;
        if (del - C64::ONE).abs() < 1e-15 {
            break;
        }
    }
    h
}

/// The bounded wave part `W` of the source potential, and its gradient
/// `(∂W/∂y, ∂W/∂z)` with respect to the **field** point, from the offsets
/// `dy = y − η` and `zs = z + ζ ≤ 0`.
pub fn wave_part(nu: f64, dy: f64, zs: f64) -> (C64, C64, C64) {
    let w = C64::new(nu * zs, nu * dy.abs());
    let (p, reg) = p_and_regular(w);
    let decay = (nu * zs).exp();
    let (s, c) = (nu * dy).sin_cos();
    let rad = C64::new(0.0, -2.0 * PI * decay);
    let val = C64::new(-2.0 * reg.re + 2.0 * nu.ln(), 0.0) + rad.scale(c);
    if !p.re.is_finite() {
        // At the source's own image point the gradient is singular (a
        // logarithm, integrable); report the radiating part only.
        return (val, rad.scale(-nu * s), rad.scale(nu * c));
    }
    let sgn = if dy >= 0.0 { 1.0 } else { -1.0 };
    let dy_part = C64::new(2.0 * nu * sgn * p.im, 0.0) + rad.scale(-nu * s);
    let dz_part = C64::new(-2.0 * nu * p.re, 0.0) + rad.scale(nu * c);
    (val, dy_part, dz_part)
}

/// Far-field factor: as `|y| → ∞` the source's potential tends to
/// `−2πi e^{ν(z+ζ)} e^{iν|y−η|}`. Returned is `e^{νζ} e^{−iνη sgn y}` per
/// unit source, for the side `sgn y`.
pub fn far_factor(nu: f64, eta: f64, zeta: f64, side: f64) -> C64 {
    C64::cis(-nu * eta * side).scale((nu * zeta).exp())
}

#[cfg(test)]
mod tests {
    use super::*;
    use hullgeom::quadrature::gauss_legendre;

    /// `PV∫₀^∞ e^{kZ} cos kY/(k − ν) dk` by quadrature, the pole subtracted
    /// on [0, 2ν].
    fn pv_integral(nu: f64, y: f64, z: f64) -> f64 {
        let f = |k: f64| (k * z).exp() * (k * y).cos();
        let (gx, gw) = gauss_legendre(20);
        let mut total = 0.0;
        let seg = |a: f64, b: f64, g: &dyn Fn(f64) -> f64| -> f64 {
            gx.iter()
                .zip(&gw)
                .map(|(&t, &w)| {
                    let k = 0.5 * (a + b) + 0.5 * (b - a) * t;
                    0.5 * (b - a) * w * g(k)
                })
                .sum()
        };
        let fnu = f(nu);
        let n = 400;
        for i in 0..n {
            let (a, b) = (
                2.0 * nu * i as f64 / n as f64,
                2.0 * nu * (i + 1) as f64 / n as f64,
            );
            total += seg(a, b, &|k| {
                if (k - nu).abs() < 1e-14 {
                    0.0
                } else {
                    (f(k) - fnu) / (k - nu)
                }
            });
        }
        let kmax = 2.0 * nu + 60.0 / (-z);
        let n2 = 4000;
        for i in 0..n2 {
            let a = 2.0 * nu + (kmax - 2.0 * nu) * i as f64 / n2 as f64;
            let b = 2.0 * nu + (kmax - 2.0 * nu) * (i + 1) as f64 / n2 as f64;
            total += seg(a, b, &|k| f(k) / (k - nu));
        }
        total
    }

    #[test]
    fn the_principal_value_is_the_real_part_of_p() {
        for &(nu, y, z) in &[
            (1.0, 0.3, -0.5),
            (1.0, 2.0, -0.2),
            (2.5, 0.0, -0.4),
            (0.7, 5.0, -1.0),
            (3.0, 1.0, -0.05),
            (4.0, 0.2, -3.0),
        ] {
            let w = C64::new(nu * z, nu * f64::abs(y));
            let (p, _) = p_and_regular(w);
            let want = pv_integral(nu, y, z);
            assert!(
                (p.re - want).abs() < 1e-7 * want.abs().max(1.0),
                "ν {nu} y {y} z {z}: {} vs {want}",
                p.re
            );
        }
    }

    #[test]
    fn series_and_continued_fraction_agree() {
        for &(re, im) in &[
            (-9.5, 0.5),
            (-10.5, 3.0),
            (-3.0, 9.8),
            (-7.0, 7.2),
            (-10.2, 0.0),
        ] {
            let w = C64::new(re, im);
            let ew = w.exp();
            // Series result as e^w E1(w) = P − iπ e^w (off the axis).
            let r = w.abs();
            let (p, _) = p_and_regular(w);
            let series_or_cf = p;
            let other = {
                let e = ew_e1_continued_fraction(w);
                if im == 0.0 {
                    e
                } else {
                    e + C64::new(0.0, PI) * ew
                }
            };
            if r < SERIES_LIMIT {
                assert!(
                    (series_or_cf - other).abs() < 1e-9 * other.abs().max(1.0),
                    "{w:?}: {series_or_cf:?} vs {other:?}"
                );
            }
        }
    }

    /// Near the negative real axis the power series serves far past
    /// `SERIES_LIMIT`, its error scaled away by `e^{Re w}`.
    #[test]
    fn the_series_holds_near_the_negative_axis() {
        for &(re, im) in &[(-15.0, 2.0), (-25.0, 8.0), (-35.0, 12.0), (-20.0, 0.3)] {
            let w = C64::new(re, im);
            let (p, _) = p_and_regular(w);
            let cf = ew_e1_continued_fraction(w) + C64::new(0.0, PI) * w.exp();
            assert!((p - cf).abs() < 1e-9 * cf.abs(), "{w:?}: {p:?} vs {cf:?}");
        }
    }

    #[test]
    fn the_asymptotic_series_matches_the_continued_fraction() {
        for &(re, im) in &[
            (-45.0, 5.0),
            (-30.0, 30.0),
            (0.0, 60.0),
            (-60.0, 20.0),
            (-41.0, 0.5),
        ] {
            let w = C64::new(re, im);
            let (a, c) = (ew_e1_asymptotic(w), ew_e1_continued_fraction(w));
            assert!((a - c).abs() < 1e-12 * c.abs(), "{w:?}: {a:?} vs {c:?}");
        }
    }

    #[test]
    fn the_source_radiates_outgoing_waves() {
        let nu = 1.3;
        let (zeta, z) = (-0.4, -0.2);
        for &y in &[30.0, 60.0] {
            let (wv, _, _) = wave_part(nu, y, z + zeta);
            let r = (y * y + (z - zeta) * (z - zeta)).sqrt().ln();
            let r1 = (y * y + (z + zeta) * (z + zeta)).sqrt().ln();
            let g = wv + C64::new(r + r1, 0.0);
            let want = C64::new(0.0, -2.0 * PI) * C64::cis(nu * y).scale((nu * (z + zeta)).exp());
            assert!((g - want).abs() < 2e-2, "y {y}: {g:?} vs {want:?}");
        }
    }

    #[test]
    fn the_gradient_matches_finite_differences() {
        let nu = 1.7;
        for &(dy, zs) in &[(0.4, -0.3), (-1.2, -0.8), (3.0, -0.05)] {
            let (_, gy, gz) = wave_part(nu, dy, zs);
            let h = 1e-6;
            let fy = (wave_part(nu, dy + h, zs).0 - wave_part(nu, dy - h, zs).0).scale(0.5 / h);
            let fz = (wave_part(nu, dy, zs + h).0 - wave_part(nu, dy, zs - h).0).scale(0.5 / h);
            assert!(
                (gy - fy).abs() < 1e-5 * (1.0 + fy.abs()),
                "{dy},{zs}: {gy:?} vs {fy:?}"
            );
            assert!(
                (gz - fz).abs() < 1e-5 * (1.0 + fz.abs()),
                "{dy},{zs}: {gz:?} vs {fz:?}"
            );
        }
    }
}
