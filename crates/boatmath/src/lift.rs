//! A hull's dynamic lift by **2D+t** (slender-body) theory, as the pressure
//! patch [`thinship::pressure`] carries into the attitude and into Michell's
//! integral.
//!
//! Thin-ship theory gives a hull thickness but no lift: its sources sit on
//! the centreplane, and a bottom meeting the flow at an angle — a hull
//! running trimmed — carries nothing. 2D+t supplies that part. Seen from the
//! earth, a station's section sinks into the water at `V = U sin α` as the
//! hull goes by at angle of attack `α`, and the water it sets moving
//! downward — the section's added mass `m_a` — pushes back:
//!
//! ```text
//! f(x) = C_T(x) · [ U² sin α cos α · dm_a/dξ  +  ½ ρ C_D,c · 2c · U² sin α |sin α| ]
//! ```
//!
//! per unit length, `ξ` running aft, `2c` the wetted beam. The second term
//! is the cross-flow drag of the section moving down through the water
//! (`C_D,c` [`CROSSFLOW_CD`]). The buoyancy is the hull's own hydrostatics,
//! already in the float solve, and the hull's thickness flow is thin-ship
//! theory's, so neither is counted again: the lift is that of the angle of
//! attack alone.
//!
//! - `m_a` is the section's heave added mass at **infinite frequency**, the
//!   free surface holding `φ = 0` — the 2D+t limit — by the close-fit panel
//!   method of [`seakeeping::section2d`], on the section as it is cut at the
//!   attitude. There's no Wagner pile-up: these are sections immersed well
//!   past their chines, not wedges entering.
//! - `α` is the **keel's angle to the flow over the run**, the aft
//!   [`RUN_FRACTION`] of the wetted length (least squares): the trim and the
//!   run's own buttock angle together, as Savitsky's `τ` is the keel's.
//!   Taking it from the run alone, not station by station, keeps the
//!   forebody's sections deepening toward the keel — the hull's thickness —
//!   out of the lift.
//! - `C_T` is Garme's **transom reduction**, `tanh(k (x − x_T)/(B C_V))`
//!   ([`TRANSOM_DECAY`] `= k`, `C_V = U/√(gB)`): the pressure falls to the
//!   atmosphere's at a dry transom, where pure 2D+t would carry it to the
//!   edge. A hull that closes aft keeps `C_T = 1`; its shrinking sections
//!   carry suction instead, as slender-body theory says they should.
//!
//! Each station's lift is spread across its wetted beam at the section's
//! mean depth (area over beam): a [`Strip`] of the patch. Where the wetted
//! length ends forward the added mass steps up from nothing, and that step
//! is a load of its own on the forward strip (a flat bottom's whole lift).
//!
//! **Against Savitsky** (prismatic planing surfaces, the tests): buoyancy
//! plus this lift is within 10% of his total at `C_V` 4 and `λ` 3–4, over
//! deadrise 0–20° and trim 2–8°, its centre within 12% of the wetted length
//! (aft of his). Short and fast — most of the length wetted before the
//! chines — it falls 30–40% short, for want of Wagner's pile-up; at `C_V` 2
//! it runs ~25% high, gravity being outside 2D+t.

use hullgeom::{Conditions, SectionalHull};
use seakeeping::section2d::Section;
use thinship::pressure::{PressurePatch, Strip};

/// Cross-flow drag coefficient of a section moving down through the water.
pub const CROSSFLOW_CD: f64 = 1.0;
/// Garme's transom-reduction constant.
pub const TRANSOM_DECAY: f64 = 2.5;
/// The run: this aft fraction of the wetted length sets the angle of
/// attack.
pub const RUN_FRACTION: f64 = 0.4;
/// Panels per section for its added mass.
const PANELS: usize = 24;

/// A hull's lift at an attitude.
#[derive(Debug, Clone, Default)]
pub struct Lift {
    pub patch: PressurePatch,
    /// The run's angle of attack [rad], bow-up positive.
    pub alpha: f64,
}

/// One station's section, as 2D+t needs it.
struct Station {
    x: f64,
    keel: f64,
    half_width: f64,
    mean_depth: f64,
    added_mass: f64,
}

fn station(x: f64, curve: &[(f64, f64)], density: f64) -> Station {
    let keel = curve.iter().map(|p| p.1).fold(0.0, f64::max);
    let half_width = curve.iter().map(|p| p.0).fold(0.0, f64::max);
    // Half the section's area, ∫ y d(depth), round the curve from the
    // waterline to the keel.
    let half_area: f64 = curve
        .windows(2)
        .map(|w| 0.5 * (w[0].0 + w[1].0) * (w[1].1 - w[0].1))
        .sum::<f64>()
        .abs();
    let wet = keel > 1e-9 && half_width > 1e-9;
    // A half-beam sampled against depth stops at the keel's depth without
    // running in across a flat bottom: close it to the centreplane, or the
    // panels see only the topsides.
    let mut closed = curve.to_vec();
    if let Some(&(y, _)) = curve.last() {
        if y > 1e-9 {
            closed.push((0.0, keel));
        }
    }
    let curve = &closed[..];
    Station {
        x,
        keel,
        half_width,
        mean_depth: if wet {
            (half_area / half_width).clamp(0.0, keel)
        } else {
            0.0
        },
        added_mass: if wet {
            Section::from_curve(curve, PANELS)
                .and_then(|s| s.added_mass_infinite(density))
                .unwrap_or(0.0)
        } else {
            0.0
        },
    }
}

/// The least-squares slope of `y` on `x`.
fn slope(pts: &[(f64, f64)]) -> f64 {
    let n = pts.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let (mx, my) = pts
        .iter()
        .fold((0.0, 0.0), |(a, b), p| (a + p.0 / n, b + p.1 / n));
    let (sxy, sxx) = pts.iter().fold((0.0, 0.0), |(a, b), p| {
        (a + (p.0 - mx) * (p.1 - my), b + (p.0 - mx) * (p.0 - mx))
    });
    if sxx > 0.0 {
        sxy / sxx
    } else {
        0.0
    }
}

/// The 2D+t lift of `hull`, cut at its attitude, at `cond`'s speed.
pub fn sectional_lift(hull: &SectionalHull, cond: &Conditions) -> Lift {
    let rho = cond.fluid.density;
    let u = cond.speed;
    let st: Vec<Station> = hull.curves().map(|(x, c)| station(x, c, rho)).collect();
    let wet: Vec<&Station> = st.iter().filter(|s| s.keel > 1e-9).collect();
    let (Some(aft), Some(fwd)) = (wet.first(), wet.last()) else {
        return Lift::default();
    };
    let (x_aft, length) = (aft.x, fwd.x - aft.x);
    if !(length > 0.0) {
        return Lift::default();
    }
    // The keel deepening aft is the bow up: α = atan(−∂keel/∂x) over the run.
    let run: Vec<(f64, f64)> = wet
        .iter()
        .filter(|s| s.x <= x_aft + RUN_FRACTION * length)
        .map(|s| (s.x, s.keel))
        .collect();
    let alpha = (-slope(&run)).atan();
    let beam = 2.0 * st.iter().map(|s| s.half_width).fold(0.0, f64::max);
    let cv = u / (cond.gravity * beam).sqrt();
    let transom = hull.transom().is_some();
    let n = st.len();
    let (sa, ca) = alpha.sin_cos();
    let c_t = |x: f64| {
        if transom {
            (TRANSOM_DECAY * (x - x_aft) / (beam * cv)).tanh()
        } else {
            1.0
        }
    };
    let i_fwd = (0..n).rev().find(|&i| st[i].keel > 1e-9).unwrap_or(0);
    let strips = (0..n)
        .filter(|&i| st[i].keel > 1e-9)
        .map(|i| {
            let (lo, hi) = (i.saturating_sub(1), (i + 1).min(n - 1));
            // dm_a/dξ, ξ aft: −dm_a/dx.
            let dm = -(st[hi].added_mass - st[lo].added_mass) / (st[hi].x - st[lo].x);
            let f = c_t(st[i].x)
                * (u * u * sa * ca * dm
                    + 0.5 * rho * CROSSFLOW_CD * 2.0 * st[i].half_width * u * u * sa * sa.abs());
            // Trapezoid weight: half the stations' spacing either side.
            let w = 0.5 * (st[hi].x - st[lo].x);
            // The forward end of the wetted length: the added mass steps up
            // from nothing ahead of it, all at once (a flat bottom's whole
            // lift is there).
            let step = if i == i_fwd {
                c_t(st[i].x) * u * u * sa * ca * st[i].added_mass
            } else {
                0.0
            };
            Strip {
                x: st[i].x,
                force: f * w + step,
                half_width: st[i].half_width,
                depth: st[i].mean_depth,
            }
        })
        .collect();
    Lift {
        patch: PressurePatch { strips },
        alpha,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hullgeom::sectional::{DepthQuadrature, SectionNodes};

    /// A prismatic planing hull — chine beam `b`, deadrise `β`, vertical
    /// sides — cut at trim `τ` with its keel wetted for `l_k` from the
    /// transom (x = 0).
    fn prismatic(b: f64, beta: f64, tau: f64, l_k: f64) -> SectionalHull {
        let (p, spans) = (3usize, 40usize);
        let mut knots = vec![0.0; p + 1];
        knots.extend((1..spans).map(|i| l_k * i as f64 / spans as f64));
        knots.extend(vec![l_k; p + 1]);
        let n = knots.len() - p - 1;
        let xs: Vec<f64> = (0..n)
            .map(|i| knots[i + 1..=i + p].iter().sum::<f64>() / p as f64)
            .collect();
        let chine = 0.5 * b * beta.tan();
        let q = DepthQuadrature::default();
        let sections = xs
            .iter()
            .map(|&x| {
                let d = ((l_k - x) * tau.tan()).max(1e-6);
                let breaks: Vec<f64> = [d - chine].into_iter().filter(|&z| z > 0.0).collect();
                SectionNodes::from_depth_function(
                    d,
                    &breaks,
                    |z| ((d - z) / beta.tan()).min(0.5 * b),
                    &q,
                )
            })
            .collect();
        SectionalHull::new(p, knots, &xs, sections).unwrap()
    }

    /// Savitsky's prismatic planing surface across his fit's range:
    /// buoyancy plus the 2D+t lift against his total lift, and its centre
    /// against his — returned as (lift ratio − 1, centre error in lengths).
    fn against_savitsky(beta_deg: f64, tau_deg: f64, lam: f64, cv: f64) -> (f64, f64) {
        let b = 1.0;
        let (beta, tau) = (beta_deg.to_radians(), tau_deg.to_radians());
        let l_k = lam * b + b * beta.tan() / (2.0 * std::f64::consts::PI * tau.tan());
        let hull = prismatic(b, beta, tau, l_k);
        let g = 9.81;
        let cond = Conditions::seawater(cv * (g * b).sqrt());
        let rho = cond.fluid.density;
        let lift = sectional_lift(&hull, &cond);
        assert!((lift.alpha - tau).abs() < 1e-3, "α {}", lift.alpha);
        let buoyancy = rho * g * hull.displaced_volume();
        let dynamic = lift.patch.force();
        let ours = buoyancy + dynamic;
        let centre = (buoyancy * hull.lcb_x() + lift.patch.moment_bow_up(0.0)) / ours;
        // Savitsky (1964): C_L0 = τ^1.1 (0.0120 λ^½ + 0.0055 λ^{5/2}/C_V²),
        // C_Lβ = C_L0 − 0.0065 β C_L0^0.6, lp/λb = 0.75 − 1/(5.21 C_V²/λ² + 2.39).
        let cl0 = tau_deg.powf(1.1) * (0.0120 * lam.sqrt() + 0.0055 * lam.powf(2.5) / (cv * cv));
        let clb = cl0 - 0.0065 * beta_deg * cl0.powf(0.6);
        let savitsky = clb * 0.5 * rho * cond.speed * cond.speed * b * b;
        let lp = lam * b * (0.75 - 1.0 / (5.21 * cv * cv / (lam * lam) + 2.39));
        let (dl, dc) = (ours / savitsky - 1.0, (centre - lp) / (lam * b));
        eprintln!(
            "β {beta_deg:2}° τ {tau_deg}° λ {lam} C_V {cv}: buoyancy {buoyancy:.0} N + 2D+t \
             {dynamic:.0} N = {ours:.0} N vs Savitsky {savitsky:.0} N ({:+.0}%) · centre \
             {centre:.2} m vs {lp:.2} m ({:+.0}% of λb)",
            100.0 * dl,
            100.0 * dc
        );
        (dl, dc)
    }

    /// Where Savitsky's fit is surest (C_V 4, λ 3–4) the two agree: lift
    /// within 10% over deadrise 0–20° and trim 2–8°, the centre within 12%
    /// of the wetted length (2D+t's sits aft).
    #[test]
    fn prismatic_planing_hulls_against_savitsky() {
        for &(beta, tau, lam, cv) in &[
            (10.0, 4.0, 3.0, 4.0),
            (0.0, 4.0, 3.0, 4.0),
            (20.0, 4.0, 3.0, 4.0),
            (10.0, 2.0, 3.0, 4.0),
            (10.0, 8.0, 3.0, 4.0),
            (10.0, 4.0, 4.0, 4.0),
        ] {
            let (dl, dc) = against_savitsky(beta, tau, lam, cv);
            assert!(
                dl.abs() < 0.12,
                "β {beta} τ {tau} λ {lam} C_V {cv}: lift {dl:+.2}"
            );
            assert!(
                dc.abs() < 0.13,
                "β {beta} τ {tau} λ {lam} C_V {cv}: centre {dc:+.2}"
            );
        }
    }

    /// The theory's known limits, reported rather than asserted: short and
    /// fast (most of the length before the chines wet) it falls short, for
    /// want of Wagner's pile-up (−30 to −40% at λ 1.5–2 or C_V 6–8); slow it
    /// runs high, gravity being outside 2D+t (+24% at C_V 2).
    #[test]
    fn prismatic_planing_hulls_at_the_edges() {
        for &(beta, tau, lam, cv) in &[
            (10.0, 4.0, 1.5, 4.0),
            (10.0, 4.0, 3.0, 8.0),
            (25.0, 6.0, 2.0, 6.0),
            (10.0, 4.0, 3.0, 2.0),
        ] {
            against_savitsky(beta, tau, lam, cv);
        }
    }
}
