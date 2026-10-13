//! Michell's thin-ship wave resistance integral.
//!
//! With ν = g/U², half-beam f(x, z), z downward from the waterline:
//!
//! ```text
//! R_w = (4 ρ g²)/(π U²) ∫_1^∞ (I² + J²) λ²/√(λ² − 1) dλ
//! I(λ) + i J(λ) = ∬ (∂f/∂x) exp(−ν λ² z) exp(i ν λ x) dx dz
//! ```
//!
//! (Tuck 1989; Dambrine, Pierre & Rousseaux 2016.) Substituting λ = sec θ
//! turns the outer integral into `∫_0^{π/2} (I² + J²) sec³θ dθ`, removing the
//! integrable singularity at λ = 1.
//!
//! This module holds what every hull representation shares: the options and
//! results, the transom closure, and the smooth outer θ-integral, integrated
//! with Gauss–Legendre panels sized to the local oscillation rate and then
//! refined until the requested tolerance is met. The amplitudes come from
//! [`crate::sectional`].
//!
//! In tests it also carries the **exact B-spline kernel** ([`InnerIntegral`]):
//! for a half-breadth spline `f(x, z)` I and J are evaluated exactly per knot
//! span via the closed-form moments in [`hullgeom::moments`]. It is the oracle
//! the sectional kernel is validated against.

#[cfg(test)]
use crate::hull::Hull;
#[cfg(test)]
use hullgeom::moments::exp_moments;
use hullgeom::moments::{osc_moments, C64};
use hullgeom::quadrature::gauss_legendre;
#[cfg(test)]
use hullgeom::Conditions;
#[cfg(test)]
use hullgeom::{Error, Result};
use std::f64::consts::{FRAC_PI_2, PI};

/// How the wave integral closes a hull whose half-breadth does not vanish at
/// the transom.
///
/// Michell's integral is over a **closed** body: the free-wave amplitude
/// `∬(∂f/∂x) …` assumes `f` reaches zero at both ends. A transom leaves a step
/// there, and taking the step at face value models a body that shuts
/// instantaneously — which radiates far too much. What a real ventilated
/// transom does is let the flow leave the edge cleanly and close in a hollow
/// some way downstream.
///
/// The closure follows the standard thin-ship treatment (Doctors & Day;
/// Couser & Molland): append a **virtual appendage** running aft from the
/// transom over a hollow length `L_v`, with the half-beam decaying
/// `f_v(x, z) = f_T(z)·φ(s)`, `s = (x_T − x)/L_v`. The shape is the smoothstep
/// `φ(s) = (1 − s)²(1 + 2s)`, flat at both ends, so the appendage introduces no
/// new discontinuity of its own: `f` is continuous at the transom and closes
/// tangentially. The bare step is the `L_v → 0` limit and is recovered exactly.
///
/// For a hull that already closes aft this is inert — every variant gives a
/// bit-for-bit unchanged result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransomClosure {
    /// Leave the transom open: classical Michell, which assumes the body
    /// closes and therefore drops the transom's contribution entirely.
    None,
    /// Hollow length from the **ballistic** estimate `L_v = c·d_T·Fn_T`, i.e.
    /// `c·U·√(d_T/g)`: water leaving the transom horizontally at `U` falls the
    /// transom depth `d_T` under gravity, and `c = √2` is the free-fall value
    /// ([`BALLISTIC_COEFF`]). Raising `c` lengthens the hollow and softens the
    /// transom's wave-making.
    Ballistic { coeff: f64 },
    /// The ballistic hollow while the transom is wet, a fraction of the
    /// transverse wavelength `2πU²/g` once it runs dry, blended between by
    /// transom Froude number `Fn_T = U/√(g·d_T)` ([`AUTO_WET_FN_T`] to
    /// [`AUTO_DRY_FN_T`]).
    ///
    /// Behind a dry transom the hollow is part of the hull's own wave, so
    /// its length scales with the wave rather than with `d_T`. The ballistic
    /// length there is a small fraction of a wavelength (0.95 m against 17 m
    /// on a catamaran at Fn 0.56, Fn_T ≈ 4), and a hollow that short closes
    /// almost as abruptly as the bare step: there it gave `R_w` 1.6× the
    /// infinite hollow's, while any hollow from half a wavelength up came
    /// within ±6% of it. This is an argument, not yet a validation: no CFD or
    /// tank result with the transom immersed has been compared.
    Auto { coeff: f64, fraction: f64 },
    /// Hollow of a prescribed length [m], independent of speed.
    ///
    /// `length = 0` is the bare step — the transom shutting instantaneously.
    /// It is the limit the closure is built to reproduce, and useful as a
    /// bound, but it is not a physical model: a discontinuous body radiates at
    /// every wave angle, so its amplitude falls off only algebraically in `λ`
    /// and the outer quadrature needs a far wider `λ` range to converge
    /// (orders of magnitude slower than any positive hollow).
    Fixed { length: f64 },
}

/// The free-fall hollow-length coefficient: a particle leaving the transom
/// horizontally at `U` falls `d_T` in `√(2d_T/g)`, travelling `√2·U·√(d_T/g)`.
pub const BALLISTIC_COEFF: f64 = std::f64::consts::SQRT_2;

/// [`TransomClosure::Auto`]'s dry hollow, in transverse wavelengths.
pub const DRY_HOLLOW_FRACTION: f64 = 0.5;

/// The transom Froude number below which [`TransomClosure::Auto`] takes the
/// transom as wet (the ballistic hollow) …
pub const AUTO_WET_FN_T: f64 = 2.5;

/// … and above which it takes it as dry (a fraction of a wavelength).
pub const AUTO_DRY_FN_T: f64 = 3.5;

impl Default for TransomClosure {
    fn default() -> Self {
        TransomClosure::Auto {
            coeff: BALLISTIC_COEFF,
            fraction: DRY_HOLLOW_FRACTION,
        }
    }
}

impl TransomClosure {
    /// Hollow length [m] for a transom of immersion `depth` at `ν = g/U²`.
    /// `None` switches the closure off; a non-positive length collapses to the
    /// bare step, which the amplitude handles as the `L_v → 0` limit.
    pub fn hollow_length(self, depth: f64, nu: f64) -> Option<f64> {
        match self {
            TransomClosure::None => None,
            // L_v = c·U·√(d_T/g) = c·√(d_T/ν), since ν = g/U².
            TransomClosure::Ballistic { coeff } => Some(coeff * (depth / nu).sqrt()),
            TransomClosure::Auto { coeff, fraction } => {
                let wet = coeff * (depth / nu).sqrt();
                let dry = fraction * 2.0 * PI / nu;
                // Fn_T² = U²/(g·d_T) = 1/(ν·d_T); smoothstep across the band.
                let fn_t = 1.0 / (nu * depth).sqrt();
                let t = ((fn_t - AUTO_WET_FN_T) / (AUTO_DRY_FN_T - AUTO_WET_FN_T)).clamp(0.0, 1.0);
                let w = t * t * (3.0 - 2.0 * t);
                Some(wet + w * (dry - wet))
            }
            TransomClosure::Fixed { length } => Some(length),
        }
        .filter(|l| l.is_finite() && *l >= 0.0)
    }
}

/// Options for the outer-integral quadrature.
#[derive(Debug, Clone, Copy)]
pub struct WaveOptions {
    /// Target relative tolerance on the resistance.
    pub rel_tol: f64,
    /// Maximum number of panel-halving refinement passes.
    pub max_refinements: usize,
    /// How to close a transom stern. Inert for a hull that closes aft.
    pub transom: TransomClosure,
}

impl Default for WaveOptions {
    fn default() -> Self {
        WaveOptions {
            rel_tol: 1e-5,
            max_refinements: 4,
            transom: TransomClosure::default(),
        }
    }
}

/// `∫_0^1 φ'(s) e^{iKs} ds` for the smoothstep hollow `φ = (1−s)²(1+2s)`,
/// whose derivative is `6s² − 6s`. At `K = 0` this is `φ(1) − φ(0) = −1`,
/// the bare-step limit.
pub(crate) fn hollow_shape_moment(k: f64, scratch: &mut Vec<C64>) -> C64 {
    osc_moments(k, 1.0, 2, scratch);
    scratch[2].scale(6.0) - scratch[1].scale(6.0)
}

/// Wave resistance result with quadrature diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct WaveResistance {
    /// Wave resistance R_w [N].
    pub resistance: f64,
    /// Estimated relative quadrature error (difference between the last two
    /// refinement passes).
    pub est_rel_error: f64,
    /// Total number of inner-integral evaluations performed.
    pub inner_evaluations: usize,
    /// Largest λ = sec θ reached before truncation.
    pub max_lambda: f64,
}

#[cfg(test)]
use hullgeom::Placement;

#[cfg(test)]
/// Compute Michell wave resistance with explicit quadrature options.
pub(crate) fn wave_resistance_with(
    hull: &Hull,
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    multihull_wave_resistance_with(&[(hull, Placement::default())], cond, opts)
}

#[cfg(test)]
/// Combined wave resistance of several thin hulls (multihull), with default
/// options. See [`multihull_wave_resistance_with`].
pub(crate) fn multihull_wave_resistance(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
) -> Result<WaveResistance> {
    multihull_wave_resistance_with(members, cond, &WaveOptions::default())
}

#[cfg(test)]
/// Combined wave resistance of several thin hulls.
///
/// In thin-ship theory the far-field free-wave amplitudes superpose: hull `j`
/// at longitudinal offset `Δx_j` and transverse position `y_j` contributes
/// `F_j(λ) · exp(i ν (λ Δx_j ± λ √(λ²−1) y_j))` to the wave system
/// propagating at ±θ off the track, and the resistance integrand is the mean
/// of the two systems, `½ (|Σ_j …₊|² + |Σ_j …₋|²)` (the θ < 0 half of the
/// free-wave spectrum; the two differ only for fleets that are not
/// mirror-symmetric about their mean centerplane, e.g. staggered pairs). For
/// two identical hulls separated by `s` this reduces to the classical
/// catamaran interference factor `4 cos²(½ ν s λ √(λ²−1))`.
pub(crate) fn multihull_wave_resistance_with(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    validate_fleet(members, cond, opts)?;
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let rho = cond.fluid.density;

    let (cx_ref, y_ref) = fleet_phase_refs(members);
    let params = fleet_outer_params(members, nu, cx_ref, y_ref);

    let mem: Vec<SourceMember> = members
        .iter()
        .map(|(h, p)| SourceMember {
            inner: InnerIntegral::new(h, nu, opts.transom),
            dx: h.x_center() + p.x - cx_ref,
            dy: p.y - y_ref,
        })
        .collect();

    let coeff = 4.0 * rho * g * g / (PI * u * u);
    Ok(run_outer(&params, opts, coeff, mem))
}

#[cfg(test)]
/// The Michell inner integrals `(I(λ), J(λ))` — the free-wave amplitude
/// functions — evaluated exactly (per-span closed forms) for `λ >= 1`.
///
/// Phases are taken relative to the hull's x-midpoint, so I and J individually
/// depend on that (physically irrelevant) choice of origin while `I² + J²`
/// does not.
pub(crate) fn inner_integrals(hull: &Hull, cond: &Conditions, lambda: f64) -> Result<(f64, f64)> {
    inner_integrals_with(hull, cond, lambda, &WaveOptions::default())
}

#[cfg(test)]
/// [`inner_integrals`] with an explicit transom closure (the rest of
/// [`WaveOptions`] governs the outer quadrature and does not apply here).
pub(crate) fn inner_integrals_with(
    hull: &Hull,
    cond: &Conditions,
    lambda: f64,
    opts: &WaveOptions,
) -> Result<(f64, f64)> {
    cond.validate()?;
    if !(lambda.is_finite() && lambda >= 1.0) {
        return Err(Error::InvalidConditions(format!(
            "lambda must be >= 1, got {lambda}"
        )));
    }
    let nu = cond.gravity / (cond.speed * cond.speed);
    let mut inner = InnerIntegral::new(hull, nu, opts.transom);
    let f = inner.eval(lambda);
    Ok((f.re, f.im))
}

/// Geometry-derived phase-rate parameters for the outer quadrature.
pub(crate) struct OuterParams {
    pub(crate) nu: f64,
    /// Half-extent of the whole fleet about its longitudinal phase centre.
    pub(crate) x_half: f64,
    /// Largest transverse offset from the fleet's phase centre.
    pub(crate) y_half: f64,
    /// Deepest draft in the fleet.
    pub(crate) t_max: f64,
    /// An integral the march needn't resolve below: truncation and
    /// refinement are relative to the larger of it and the integral itself.
    /// Zero for a fleet on its own; a member's share of a larger fleet (a
    /// drive's leg, say) is wanted only to the fleet's precision.
    pub(crate) floor: f64,
}

/// Marching-panel Gauss–Legendre integration of
/// `∫_0^{π/2} |A(sec θ)|² sec³θ dθ`, where `|A|²` is the combined amplitude
/// of the fleet's two (±θ) wave systems ([`superpose`] over `members`).
///
/// **Parallel structure.** The panel schedule depends only on the geometry
/// (`params`) and `frac` — never on integrand values — so it is generated
/// ahead in batches, every panel of a batch is evaluated independently
/// (each worker on its own clone of the members, see [`hullgeom::parallel`]),
/// and the truncation test is then applied to the panels *in θ order*. The
/// accumulation order and every per-panel operation are those of the plain
/// serial march, so the result is bit-for-bit independent of the thread
/// count; the only cost of parallelism is the tail of the final batch past
/// the truncation point, bounded by the batch size. `theta_hint` — where the
/// previous (coarser) pass truncated — lets a refinement pass schedule its
/// whole expected range as one batch; the first pass grows its batches.
///
/// Returns the integral and the θ the march stopped at.
fn integrate_outer<M: MemberWave>(
    params: &OuterParams,
    frac: f64,
    members: &[M],
    theta_hint: Option<f64>,
    evals: &mut usize,
) -> (f64, f64) {
    const GL_N: usize = 16;
    /// Truncate once a full quiet window contributes below this fraction.
    const STOP_REL: f64 = 1e-9;
    /// Width of the quiet window in accumulated phase — several full periods
    /// of the oscillating integrand, so a trough of cos²-type oscillation
    /// (width < π) can never trigger truncation on its own. Phase-based, so
    /// the criterion is independent of the panel-refinement level.
    const STOP_WINDOW_PHASE: f64 = 8.0 * PI;
    const LAMBDA_HARD_CAP: f64 = 1e4;
    const MAX_EVALS_PER_PASS: usize = 4_000_000;
    /// Panels per worker in a first batch: small enough that a pass which
    /// truncates early wastes little, large enough to amortise the thread
    /// spawn; batches double while the march continues.
    const FIRST_BATCH_PER_WORKER: usize = 32;
    const MAX_BATCH: usize = 2048;

    let (gx, gw) = gauss_legendre(GL_N);
    let nu = params.nu;
    let (x_half, y_half, t_max) = (params.x_half, params.y_half, params.t_max);

    // Local phase rate of |A|² in θ: the x-oscillation contributes
    // 2 ν x_half d(sec θ)/dθ, the z-decay envelope 2 ν T d(sec²θ)/dθ, and the
    // transverse separation phase ν y λ√(λ²−1) = ν y sec θ tan θ contributes
    // 2 ν y_half d(sec θ tan θ)/dθ = 2 ν y_half sec θ (sec²θ + tan²θ).
    let rate = |theta: f64| -> f64 {
        let sec = 1.0 / theta.cos();
        let tan = theta.tan();
        2.0 * nu * sec * tan * (x_half + t_max * sec)
            + 2.0 * nu * y_half * sec * (sec * sec + tan * tan)
            + 4.0
    };
    // Near θ = 0 the longitudinal phase grows like ν x_half θ², so cap the
    // first panels at one period of that quadratic phase (the transverse
    // phase is linear near 0 and already covered by rate(0)).
    let cap = (2.0 * PI / (2.0 * nu * x_half).sqrt().max(1.0)).min(0.12);

    /// One scheduled panel: `[theta, theta + dt]` with the local phase rate
    /// that sized it.
    struct Panel {
        theta: f64,
        dt: f64,
        rate: f64,
    }

    // One panel's Gauss–Legendre sum (before the half-width factor).
    let panel_sum = |mem: &mut Vec<M>, p: &Panel| -> f64 {
        let half = p.dt / 2.0;
        let mid = p.theta + half;
        let mut panel = 0.0;
        for (k, &xi) in gx.iter().enumerate() {
            let th = mid + half * xi;
            let sec = 1.0 / th.cos();
            panel += gw[k] * superpose(mem, nu, sec) * sec * sec * sec;
        }
        panel
    };

    // With a single worker there is nothing to batch: march panel by panel
    // on one private copy of the members, exactly as the serial loop did,
    // so no panel past the truncation point is ever evaluated.
    let workers = hullgeom::parallel::threads();
    let serial = workers == 1;
    let mut local: Vec<M> = if serial { members.to_vec() } else { Vec::new() };

    let mut theta = 0.0f64;
    let mut total = 0.0f64;
    let mut window_sum = 0.0f64;
    let mut window_phase = 0.0f64;
    let mut pass_evals = 0usize;
    let mut batch = if serial {
        1
    } else {
        FIRST_BATCH_PER_WORKER * workers
    };
    let mut stopped = false;
    while !stopped && theta < FRAC_PI_2 - 1e-12 {
        // Schedule the next batch, mirroring the march's own stopping rules
        // (hard λ cap, evaluation budget) so no panel is scheduled that the
        // serial loop would not have reached.
        // The first batch of a refinement pass runs straight to the hint.
        let to_hint = theta_hint.filter(|_| theta == 0.0 && !serial);
        let mut panels: Vec<Panel> = Vec::with_capacity(batch);
        let mut th = theta;
        let mut ev = pass_evals;
        while (panels.len() < batch || to_hint.is_some_and(|h| th < h)) && th < FRAC_PI_2 - 1e-12 {
            let local_rate = rate(th);
            let dt = (frac * 2.0 * PI / local_rate)
                .min(frac * cap)
                .min(FRAC_PI_2 - th)
                .max(1e-15);
            panels.push(Panel {
                theta: th,
                dt,
                rate: local_rate,
            });
            th += dt;
            ev += GL_N;
            if 1.0 / th.cos() > LAMBDA_HARD_CAP || ev > MAX_EVALS_PER_PASS {
                break;
            }
        }

        // Evaluate every panel's Gauss–Legendre sum, independently.
        let sums: Vec<f64> = if serial {
            panels.iter().map(|p| panel_sum(&mut local, p)).collect()
        } else {
            hullgeom::parallel::map_indexed(
                panels.len(),
                || members.to_vec(),
                |mem, i| panel_sum(mem, &panels[i]),
            )
        };

        // Accumulate in θ order and apply the truncation rules exactly as the
        // serial march does, panel by panel.
        for (p, sum) in panels.iter().zip(sums) {
            let panel = sum * (p.dt / 2.0);
            total += panel;
            pass_evals += GL_N;
            theta = p.theta + p.dt;

            // Truncation: only past λ = 2, and only when an entire window of
            // accumulated oscillation phase contributed negligibly.
            if 1.0 / theta.cos() > 2.0 {
                window_sum += panel;
                window_phase += p.rate * p.dt;
                if window_phase >= STOP_WINDOW_PHASE {
                    if window_sum.abs()
                        <= STOP_REL * total.abs().max(params.floor) + f64::MIN_POSITIVE
                    {
                        stopped = true;
                        break;
                    }
                    window_sum = 0.0;
                    window_phase = 0.0;
                }
            }
            if 1.0 / theta.cos() > LAMBDA_HARD_CAP || pass_evals > MAX_EVALS_PER_PASS {
                stopped = true;
                break;
            }
        }
        if !serial {
            batch = (batch * 2).min(MAX_BATCH);
        }
    }
    *evals += pass_evals;
    (total, theta)
}

// ---------------------------------------------------------------------------
// Shared multihull driver
//
// Every wave-resistance path integrates the same outer integral of the fleet's
// combined free-wave amplitude; they differ only in how each member's ±θ
// amplitude is formed. `MemberWave` captures that per-member amplitude;
// `superpose` shares the placement-phase superposition and `run_outer` the
// panel refinement and result assembly.
// ---------------------------------------------------------------------------

/// One fleet member's contribution to the far-field free-wave amplitude.
/// `Clone + Send` so the outer quadrature can hand each worker thread its own
/// copy (the scratch buffers are per-member; the hull itself is shared).
pub(crate) trait MemberWave: Clone + Send + Sync {
    /// `(A₊, A₋)` — the amplitudes this member carries into the +θ and −θ wave
    /// systems at `λ` (before the placement phase).
    fn amps(&mut self, nu: f64, lambda: f64) -> (C64, C64);
    /// `(dx, dy)`: longitudinal and transverse offsets from the fleet phase
    /// centre. The placement phase is `exp(iν(λ dx ± λ√(λ²−1) dy))`.
    fn offsets(&self) -> (f64, f64);
}

/// Combined `½(|A₊|² + |A₋|²)` of the fleet's two Kelvin half-systems at `λ`,
/// summing each member's amplitude with its placement phase. Shared by every
/// multihull path — the members differ only in [`MemberWave::amps`].
fn superpose<M: MemberWave>(members: &mut [M], nu: f64, lambda: f64) -> f64 {
    let kx = nu * lambda;
    let ky = nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
    let mut plus = C64::ZERO;
    let mut minus = C64::ZERO;
    for m in members.iter_mut() {
        let (a_plus, a_minus) = m.amps(nu, lambda);
        if a_plus == C64::ZERO && a_minus == C64::ZERO {
            continue;
        }
        let (dx, dy) = m.offsets();
        plus = plus + a_plus * C64::cis(kx * dx + ky * dy);
        minus = minus + a_minus * C64::cis(kx * dx - ky * dy);
    }
    0.5 * (plus.abs_sq() + minus.abs_sq())
}

/// March the outer integral to the requested tolerance and assemble the
/// [`WaveResistance`]; `coeff = 4ρg²/(πU²)` is the Michell prefactor and
/// `members` the fleet whose combined amplitude is integrated. Shared by
/// every wave-resistance entry point.
pub(crate) fn run_outer<M: MemberWave>(
    params: &OuterParams,
    opts: &WaveOptions,
    coeff: f64,
    members: Vec<M>,
) -> WaveResistance {
    let mut evals_total = 0usize;
    let mut frac = 1.0;
    let mut evals = 0usize;
    let (mut integral, mut theta_stop) = integrate_outer(params, frac, &members, None, &mut evals);
    evals_total += evals;
    let mut est_rel = f64::INFINITY;
    for _ in 0..opts.max_refinements {
        frac *= 0.5;
        let mut evals = 0usize;
        let (refined, th) = integrate_outer(params, frac, &members, Some(theta_stop), &mut evals);
        evals_total += evals;
        let scale = refined.abs().max(params.floor).max(f64::MIN_POSITIVE);
        est_rel = (refined - integral).abs() / scale;
        integral = refined;
        theta_stop = th;
        if est_rel <= opts.rel_tol {
            break;
        }
    }
    WaveResistance {
        resistance: coeff * integral,
        est_rel_error: est_rel,
        inner_evaluations: evals_total,
        max_lambda: 1.0 / theta_stop.cos().max(1e-300),
    }
}

#[cfg(test)]
/// Shared validation for the multihull entry points.
fn validate_fleet(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<()> {
    cond.validate()?;
    if !(opts.rel_tol.is_finite() && opts.rel_tol > 0.0) {
        return Err(Error::InvalidConditions(
            "rel_tol must be finite and positive".into(),
        ));
    }
    if members.is_empty() {
        return Err(Error::InvalidConditions(
            "at least one hull is required".into(),
        ));
    }
    if members
        .iter()
        .any(|(_, p)| !(p.x.is_finite() && p.y.is_finite()))
    {
        return Err(Error::InvalidConditions(
            "hull placements must be finite".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
/// Fleet phase references: the mean hull x-centre and mean transverse position.
/// Constant overall phase is irrelevant; centring the oscillatory arguments
/// keeps them small.
fn fleet_phase_refs(members: &[(&Hull, Placement)]) -> (f64, f64) {
    let n = members.len() as f64;
    let cx_ref = members.iter().map(|(h, p)| h.x_center() + p.x).sum::<f64>() / n;
    let y_ref = members.iter().map(|(_, p)| p.y).sum::<f64>() / n;
    (cx_ref, y_ref)
}

#[cfg(test)]
/// Outer-integral panel-sizing parameters for a fleet.
fn fleet_outer_params(
    members: &[(&Hull, Placement)],
    nu: f64,
    cx_ref: f64,
    y_ref: f64,
) -> OuterParams {
    OuterParams {
        nu,
        x_half: members
            .iter()
            .map(|(h, p)| (h.x_center() + p.x - cx_ref).abs() + h.x_half_extent())
            .fold(0.0, f64::max),
        y_half: members
            .iter()
            .map(|(_, p)| (p.y - y_ref).abs())
            .fold(0.0, f64::max),
        t_max: members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
        floor: 0.0,
    }
}

#[cfg(test)]
/// Source (thickness) member: its free-wave amplitude, the same for the ±θ
/// systems before the placement phase.
#[derive(Clone)]
struct SourceMember<'h> {
    inner: InnerIntegral<'h>,
    dx: f64,
    dy: f64,
}

#[cfg(test)]
impl MemberWave for SourceMember<'_> {
    fn amps(&mut self, _nu: f64, lambda: f64) -> (C64, C64) {
        let f = self.inner.eval(lambda);
        (f, f)
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

/// Near-field hull transforms at one wavenumber pair, in the `e^{−i k_x x}`
/// convention with `x` measured from the hull's x-centre:
///
/// ```text
/// q     = ∬ ∂f/∂x          e^{−ik_x x} e^{−κz}     (the Michell amplitude, conjugated)
/// p     = ∬ f              e^{−ik_x x} e^{−κz}
/// q1    = ∬ (x−x_c) ∂f/∂x  e^{−ik_x x} e^{−κz}
/// w     = ∫ ∂f/∂x(x,0)         e^{−ik_x x} dx       (waterline)
/// p_wl  = ∫ f(x,0)             e^{−ik_x x} dx
/// q1_wl = ∫ (x−x_c) ∂f/∂x(x,0) e^{−ik_x x} dx
/// ```
///
/// The common phase `e^{ik_x x_c}` cancels in every product the sinkage and
/// trim integrals form; the moment reference is applied by the caller as
/// `(x − x_ref) = (x − x_c) + (x_c − x_ref)`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SquatTransforms {
    /// The source strength's transform: the hull's `∂f/∂x` plus, on a
    /// transom stern, the closing virtual appendage. What makes the pressure.
    pub q_src: C64,
    /// The rest are the **weights**: the real hull the pressure acts on,
    /// which ends at its transom — `f` drops from `f_T` to 0 there, a step
    /// its `∂f/∂x` carries. The appendage is a hollow in the water whose
    /// surface is free (zero pressure), so it adds sources but no weight.
    pub q: C64,
    pub p: C64,
    pub q1: C64,
    pub w: C64,
    pub p_wl: C64,
    pub q1_wl: C64,
}

#[cfg(test)]
/// A hull's coefficient nets contracted against the z-moments at one decay
/// rate `κ` (see [`InnerIntegral::contract_z`]).
#[derive(Debug, Clone, Default)]
pub(crate) struct ZContracted {
    /// `Σ_{t,b} c_fx[s,t,a,b] Zm_t[b]`, indexed `[s * p + a]`.
    pub g_fx: Vec<f64>,
    /// The same for `f`, indexed `[s * (p + 1) + a]`.
    pub g_f: Vec<f64>,
    /// `∫ f_T(z) e^{−κz} dz`: the transom section's z-factor.
    pub z_t: f64,
    /// False when the whole hull lies below the exponential's support.
    pub any: bool,
}

#[cfg(test)]
impl ZContracted {
    /// Reset without releasing the buffers (they are reused per k-node).
    fn clear(&mut self) {
        self.g_fx.clear();
        self.g_f.clear();
        self.z_t = 0.0;
        self.any = false;
    }
}

impl Default for SquatTransforms {
    fn default() -> Self {
        SquatTransforms {
            q_src: C64::ZERO,
            q: C64::ZERO,
            p: C64::ZERO,
            q1: C64::ZERO,
            w: C64::ZERO,
            p_wl: C64::ZERO,
            q1_wl: C64::ZERO,
        }
    }
}

#[cfg(test)]
impl SquatTransforms {
    fn conj_in_place(&mut self) {
        let conj = |v: C64| C64::new(v.re, -v.im);
        self.q_src = conj(self.q_src);
        self.q = conj(self.q);
        self.p = conj(self.p);
        self.q1 = conj(self.q1);
        self.w = conj(self.w);
        self.p_wl = conj(self.p_wl);
        self.q1_wl = conj(self.q1_wl);
    }
}

/// The weights' share of a transom at `x_T` (`phase = e^{ik_x(x_T − x_c)}`,
/// `dx_t = x_T − x_c`): the hull's half-breadth steps from `f_T` to 0 there,
/// so `∂f/∂x` carries `f_T δ(x − x_T)` — through the depth z-factor `z_t` in
/// the volume transforms and the waterline half-beam `f_t0` in the waterline
/// ones. `p` (the transform of `f` itself) has no share.
pub(crate) fn add_transom_step(
    phase: C64,
    dx_t: f64,
    z_t: f64,
    f_t0: f64,
    t: &mut SquatTransforms,
) {
    t.q = t.q + phase.scale(z_t);
    t.w = t.w + phase.scale(f_t0);
    t.q1 = t.q1 + phase.scale(dx_t * z_t);
    t.q1_wl = t.q1_wl + phase.scale(dx_t * f_t0);
}

/// The closing appendage's `∬ ∂f_v/∂x` per unit z-factor, from the
/// oscillatory moments `m` of `e^{−ik_x L_v s}` on `s ∈ [0, 1]` (orders 0–3):
/// `−phase · ∫φ′`, with `φ = 1 − 3s² + 2s³`, `φ′ = 6s² − 6s`.
pub(crate) fn appendage_source(phase: C64, m: &[C64]) -> C64 {
    C64::ZERO - phase * (m[2].scale(6.0) - m[1].scale(6.0))
}

/// What the near-field (sinkage/trim) quadrature needs from a hull: every
/// `κ`-dependent quantity once per decay rate, then the six transforms at
/// any number of `k_x` from it. Implemented by the sectional kernel
/// ([`crate::sectional::SectionalKernel`]) and, in tests, the exact B-spline
/// oracle (`InnerIntegral`).
pub(crate) trait NearFieldKernel: Clone + Send + Sync {
    type Contracted: Clone + Default + Send + Sync;
    fn contract_z(&mut self, kappa: f64, out: &mut Self::Contracted);
    fn transforms_at(&mut self, zc: &Self::Contracted, kx: f64) -> SquatTransforms;
}

#[cfg(test)]
impl NearFieldKernel for InnerIntegral<'_> {
    type Contracted = ZContracted;
    fn contract_z(&mut self, kappa: f64, out: &mut ZContracted) {
        InnerIntegral::contract_z(self, kappa, out)
    }
    fn transforms_at(&mut self, zc: &ZContracted, kx: f64) -> SquatTransforms {
        InnerIntegral::transforms_at(self, zc, kx)
    }
}

#[cfg(test)]
#[derive(Clone)]
pub(crate) struct InnerIntegral<'h> {
    hull: &'h Hull,
    nu: f64,
    transom: TransomClosure,
    /// Scratch: z-moments including the e^{−κ z0} shift, [n_spans_z][q+1].
    zm: Vec<f64>,
    /// Scratch: raw x-moments for one span.
    xm: Vec<C64>,
    /// Scratch: raw z-moments for one span.
    zm_raw: Vec<f64>,
    /// Scratch: hollow-shape moments for the transom closure.
    sm: Vec<C64>,
    /// Number of leading (shallowest) z-spans still above the decay floor at
    /// the current `κ`; the deeper ones are exactly zero in `zm`.
    z_active: usize,
    p: usize,
    q: usize,
}

#[cfg(test)]
/// Relative size below which a z-span's `e^{−κ z₀}` factor cannot matter.
///
/// The amplitude is a sum over z-spans whose shallowest term carries factor
/// `1`, so dropping every span under `1e-20` perturbs it by `≲ n_spans·1e-20`
/// — four orders below double epsilon. What it buys is the whole point of
/// affording a finely-resolved loft: at large `λ`, `κ = νλ²` confines the
/// integrand to a sliver below the waterline, and the outer quadrature spends
/// most of its evaluations out there.
const Z_DECAY_FLOOR: f64 = 1e-20;

#[cfg(test)]
impl<'h> InnerIntegral<'h> {
    pub(crate) fn new(hull: &'h Hull, nu: f64, transom: TransomClosure) -> Self {
        let p = hull.surface().degree_x();
        let q = hull.surface().degree_z();
        InnerIntegral {
            hull,
            nu,
            transom,
            zm: vec![0.0; hull.spans_z().len() * (q + 1)],
            xm: Vec::with_capacity(p),
            zm_raw: Vec::with_capacity(q + 1),
            sm: Vec::with_capacity(3),
            z_active: 0,
            p,
            q,
        }
    }

    /// Free-wave amplitude of the transom's virtual appendage.
    ///
    /// With `f_v = f_T(z)·φ(s)`, `s = (x_T − x)/L_v`, substituting
    /// `x = x_T − s·L_v` separates the double integral completely:
    ///
    /// ```text
    /// F_v = −e^{i k_x x_T} · [∫ f_T(z) e^{−κz} dz] · ∫_0^1 φ'(s) e^{−i k_x L_v s} ds
    /// ```
    ///
    /// The z-factor is the transom section against the **same** per-span
    /// z-moments the hull itself uses (already filled by `fill_zm`), so the
    /// closure costs one dot product and one 3-term moment call per λ, and
    /// carries no quadrature error either.
    fn transom_term(&mut self, kx: f64) -> C64 {
        let Some(tr) = self.hull.transom() else {
            return C64::ZERO;
        };
        let Some(lv) = self.transom.hollow_length(tr.depth, self.nu) else {
            return C64::ZERO;
        };
        let q = self.q;
        let z_factor: f64 = (0..self.z_active)
            .flat_map(|t| {
                let base = t * (q + 1);
                (0..=q).map(move |b| (base + b, base + b))
            })
            .map(|(i, j)| tr.coeff[i] * self.zm[j])
            .sum();
        let shape = hollow_shape_moment(-kx * lv, &mut self.sm);
        let phase = C64::cis(kx * (tr.x - self.hull.x_center()));
        C64::ZERO - (phase * shape).scale(z_factor)
    }

    /// Every near-field hull transform the sinkage/trim integrals need, at an
    /// arbitrary wavenumber pair `(k_x, κ)` — off the dispersion curve, which
    /// is where the local (principal-value) part of the pressure lives.
    /// Composes [`Self::contract_z`] and [`Self::transforms_at`]; a quadrature
    /// that walks many `k_x` at one `κ` should call those two directly (as
    /// [`crate::squat`] does) rather than re-contracting per point.
    #[allow(dead_code, reason = "convenience wrapper exercised directly by tests")]
    pub(crate) fn eval_transforms(&mut self, kx: f64, kappa: f64) -> SquatTransforms {
        let mut zc = ZContracted::default();
        self.contract_z(kappa, &mut zc);
        self.transforms_at(&zc, kx)
    }

    /// Contract every coefficient net against the z-moments at decay `κ`:
    /// `G[s][a] = Σ_t Σ_b c[s,t,a,b] · Zm_t[b]`, for `∂f/∂x` and for `f`, plus
    /// the transom section's own z-factor. Everything that depends on `κ`
    /// happens here, once; [`Self::transforms_at`] then costs one
    /// `osc_moments` call per x-span for each `k_x`, which is what makes a
    /// 2-D wavenumber quadrature affordable on a hull with hundreds of spans.
    pub(crate) fn contract_z(&mut self, kappa: f64, out: &mut ZContracted) {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let nsx = hull.spans_x().len();
        out.clear();
        out.g_fx.resize(nsx * p, 0.0);
        out.g_f.resize(nsx * (p + 1), 0.0);
        out.any = self.fill_zm(kappa);
        if !out.any {
            return;
        }
        let nz_active = self.z_active;
        let cfx = hull.fx_coeff();
        let cf = hull.f_coeff();
        for s in 0..nsx {
            for a in 0..p {
                let mut g = 0.0;
                for tz in 0..nz_active {
                    let base = ((s * nsz + tz) * p + a) * (q + 1);
                    let zrow = &self.zm[tz * (q + 1)..(tz + 1) * (q + 1)];
                    for b in 0..=q {
                        g += cfx[base + b] * zrow[b];
                    }
                }
                out.g_fx[s * p + a] = g;
            }
            for a in 0..=p {
                let mut g = 0.0;
                for tz in 0..nz_active {
                    let base = ((s * nsz + tz) * (p + 1) + a) * (q + 1);
                    let zrow = &self.zm[tz * (q + 1)..(tz + 1) * (q + 1)];
                    for b in 0..=q {
                        g += cf[base + b] * zrow[b];
                    }
                }
                out.g_f[s * (p + 1) + a] = g;
            }
        }
        if let Some(tr) = hull.transom() {
            out.z_t = (0..nz_active)
                .flat_map(|tz| (0..=q).map(move |b| tz * (q + 1) + b))
                .map(|i| tr.coeff[i] * self.zm[i])
                .sum();
        }
    }

    /// The six transforms at `k_x`, from a contraction at some `κ`. The
    /// waterline rows (`Z = 0` on the first z-span) do not depend on `κ` and
    /// are read straight from the coefficient nets. Results are in the
    /// `e^{−ik_x x}` convention with `x` from the hull's x-centre.
    pub(crate) fn transforms_at(&mut self, zc: &ZContracted, kx: f64) -> SquatTransforms {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let x_center = hull.x_center();
        let mut t = SquatTransforms::default();
        if !zc.any {
            return t;
        }
        let cfx = hull.fx_coeff();
        let cf = hull.f_coeff();
        for (s, sx) in hull.spans_x().iter().enumerate() {
            // Moments up to X^p: f has degree p, and so does (X + d)·∂f/∂x.
            osc_moments(kx, sx.len, p, &mut self.xm);
            let d = sx.start - x_center;
            let phase = C64::cis(kx * d);
            let (mut q_s, mut q1_s, mut w_s, mut q1w_s) =
                (C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO);
            for a in 0..p {
                let g = zc.g_fx[s * p + a];
                let gw = cfx[(s * nsz * p + a) * (q + 1)];
                q_s = q_s + self.xm[a].scale(g);
                w_s = w_s + self.xm[a].scale(gw);
                // (x − x_c)·∂f/∂x = (X + d)·∂f/∂x on this span.
                let shifted = self.xm[a + 1] + self.xm[a].scale(d);
                q1_s = q1_s + shifted.scale(g);
                q1w_s = q1w_s + shifted.scale(gw);
            }
            let (mut p_s, mut pw_s) = (C64::ZERO, C64::ZERO);
            for a in 0..=p {
                let g = zc.g_f[s * (p + 1) + a];
                let gw = cf[(s * nsz * (p + 1) + a) * (q + 1)];
                p_s = p_s + self.xm[a].scale(g);
                pw_s = pw_s + self.xm[a].scale(gw);
            }
            t.q = t.q + phase * q_s;
            t.q1 = t.q1 + phase * q1_s;
            t.w = t.w + phase * w_s;
            t.q1_wl = t.q1_wl + phase * q1w_s;
            t.p = t.p + phase * p_s;
            t.p_wl = t.p_wl + phase * pw_s;
        }
        t.q_src = t.q;
        self.add_transom_transforms(kx, zc.z_t, &mut t);
        // Everything above is in the kernel's e^{+i k_x (x − x_c)} convention;
        // the near-field formulas are written for e^{−i k_x x}. All nets are
        // real, so the conversion is a conjugation.
        t.conj_in_place();
        t
    }

    /// The transom's share of the transforms (kernel convention): the step
    /// `f_T → 0` in the weights (see [`SquatTransforms`]), and the closing
    /// appendage `f_v = f_T(z)·φ(s)`, `s = (x_T − x)/L_v`, in the sources:
    ///   ∬ ∂f_v/∂x → −e^{iκ_x(x_T−x_c)} · Z_T ∫φ′ e^{−ik_xL_v s}
    fn add_transom_transforms(&mut self, kx: f64, z_t: f64, t: &mut SquatTransforms) {
        let Some(tr) = self.hull.transom() else {
            return;
        };
        let dx_t = tr.x - self.hull.x_center();
        let phase = C64::cis(kx * dx_t);
        add_transom_step(phase, dx_t, z_t, tr.half_beam, t);
        let Some(lv) = self.transom.hollow_length(tr.depth, self.nu) else {
            return;
        };
        osc_moments(-kx * lv, 1.0, 3, &mut self.sm);
        t.q_src = t.q_src + appendage_source(phase, &self.sm).scale(z_t);
    }

    /// Source free-wave amplitude `I + iJ` at λ = sec θ. Phases use x relative
    /// to the hull midpoint (a pure phase factor on I + iJ that leaves
    /// |I + iJ|² unchanged, but keeps the oscillatory arguments as small as
    /// possible).
    pub(crate) fn eval(&mut self, lambda: f64) -> C64 {
        let hull = self.hull;
        let kx = self.nu * lambda;
        let kappa = self.nu * lambda * lambda;
        if !self.fill_zm(kappa) {
            return C64::ZERO;
        }
        self.accumulate(kx, hull.fx_coeff()) + self.transom_term(kx)
    }

    /// Fill the per-span z-moments (shifted by the decay to each span start) at
    /// decay rate `kappa = ν λ²`. Returns `false` if every span has underflowed
    /// (the whole hull is below the exponential's support), in which case the
    /// amplitudes are zero.
    fn fill_zm(&mut self, kappa: f64) -> bool {
        let hull = self.hull;
        let q = self.q;
        let mut any = false;
        self.z_active = 0;
        for (t, sz) in hull.spans_z().iter().enumerate() {
            let decay = (-kappa * sz.start).exp();
            if decay < Z_DECAY_FLOOR {
                // `spans_z` runs from the waterline down, so `decay` is
                // monotonically decreasing: once a span is below the floor,
                // every deeper one is too. Zero the rest and stop — the
                // remaining spans are where most of the work at large λ used
                // to go, for a contribution below double precision.
                for z in self.zm[t * (q + 1)..].iter_mut() {
                    *z = 0.0;
                }
                break;
            }
            any = true;
            self.z_active = t + 1;
            exp_moments(kappa, sz.len, q, &mut self.zm_raw);
            for b in 0..=q {
                self.zm[t * (q + 1) + b] = decay * self.zm_raw[b];
            }
        }
        any
    }

    /// Accumulate `∬ (Σ coeff·XᵃZᵇ) e^{−κz} e^{iνλx} dx dz` over all span pairs
    /// using the already-filled z-moments and the per-span x-oscillation
    /// moments, for the `∂f/∂x` coefficient net `coeff`.
    fn accumulate(&mut self, kx: f64, coeff: &[f64]) -> C64 {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let x_center = hull.x_center();
        let mut f = C64::ZERO;
        for (s, sx) in hull.spans_x().iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut self.xm);
            let phase = C64::cis(kx * (sx.start - x_center));
            let mut span_sum = C64::ZERO;
            for (a, &xma) in self.xm.iter().enumerate() {
                // g_a = Σ_t Σ_b c[s,t,a,b] · Zm[t][b]
                let mut g_a = 0.0;
                for t in 0..self.z_active {
                    let base = ((s * nsz + t) * p + a) * (q + 1);
                    let zrow = &self.zm[t * (q + 1)..(t + 1) * (q + 1)];
                    let crow = &coeff[base..base + q + 1];
                    for b in 0..=q {
                        g_a += crow[b] * zrow[b];
                    }
                }
                span_sum = span_sum + xma.scale(g_a);
            }
            f = f + phase * span_sum;
        }
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hulls::wigley;

    /// On the dispersion pair (k_x, κ) = (νλ, νλ²) the near-field `q` is the
    /// Michell amplitude, conjugated: one code path, two entry points.
    #[test]
    fn near_field_q_is_the_conjugate_amplitude_on_the_dispersion_curve() {
        let hull = wigley(10.0, 1.0, 0.625).unwrap();
        let nu = 9.80665 / (3.0 * 3.0);
        let mut inner = InnerIntegral::new(&hull, nu, TransomClosure::default());
        for lambda in [1.0, 1.3, 2.2, 5.0] {
            let amp = inner.eval(lambda);
            let t = inner.eval_transforms(nu * lambda, nu * lambda * lambda);
            let scale = amp.abs().max(1e-300);
            assert!(
                (t.q_src.re - amp.re).abs() < 1e-12 * scale
                    && (t.q_src.im + amp.im).abs() < 1e-12 * scale,
                "λ = {lambda}: q_src = {:?}, amplitude = {:?}",
                t.q_src,
                amp
            );
        }
    }

    /// Integration by parts in x on a closed body: ∫ f e^{−ik_x x} = (1/ik_x) ∫ f_x e^{−ik_x x},
    /// i.e. `i k_x p = q` (and `i k_x p_wl = w` along the waterline). Holds
    /// off the dispersion curve too, and — the real point — for a transom
    /// hull only once the appendage is included, so this pins the appendage's
    /// share of every transform.
    #[test]
    fn near_field_transforms_satisfy_the_closure_identity() {
        let check = |hull: &Hull, closure: TransomClosure, tag: &str| {
            let nu = 9.80665 / (2.5 * 2.5);
            let mut inner = InnerIntegral::new(hull, nu, closure);
            for (kx, kappa) in [(0.7, 0.9), (2.0, 4.0), (3.5, 1.2), (nu * 1.4, nu * 1.96)] {
                let t = inner.eval_transforms(kx, kappa);
                let ikx = C64::new(0.0, kx);
                let lhs = ikx * t.p;
                let lhs_wl = ikx * t.p_wl;
                let s1 = t.q.abs().max(1e-12);
                let s2 = t.w.abs().max(1e-12);
                assert!(
                    (lhs - t.q).abs() < 1e-10 * s1,
                    "{tag} (k_x {kx}, κ {kappa}): i k_x p = {:?} vs q = {:?}",
                    lhs,
                    t.q
                );
                assert!(
                    (lhs_wl - t.w).abs() < 1e-10 * s2,
                    "{tag} (k_x {kx}, κ {kappa}): i k_x p_wl = {:?} vs w = {:?}",
                    lhs_wl,
                    t.w
                );
            }
        };
        check(
            &wigley(10.0, 1.0, 0.625).unwrap(),
            TransomClosure::None,
            "wigley",
        );

        // A wedge open at the stern closes only with its appendage.
        let knots_x = vec![0.0, 0.0, 8.0, 8.0];
        let knots_z = vec![0.0, 0.0, 0.25, 0.25];
        let wedge = Hull::new(
            hullgeom::bspline::BSplineSurface::new(
                1,
                1,
                knots_x,
                knots_z,
                vec![0.4, 0.0, 0.0, 0.0],
            )
            .unwrap(),
        )
        .unwrap();
        assert!(wedge.transom().is_some());
        check(
            &wedge,
            TransomClosure::Fixed { length: 1.3 },
            "wedge+appendage",
        );
        check(&wedge, TransomClosure::default(), "wedge+ballistic");
    }

    #[test]
    fn staggered_fleet_resistance_is_mirror_symmetric() {
        // A staggered pair and its mirror image about y = 0 are the same
        // physical system; the resistance must not change. (Regression: the
        // integrand once used only the +θ wave system, which broke this.)
        let hull = wigley(8.0, 0.8, 0.5).unwrap();
        let cond = hullgeom::Conditions::seawater(2.5);
        let a = [
            (&hull, Placement { x: 0.0, y: 1.4 }),
            (&hull, Placement { x: 1.7, y: -1.4 }),
        ];
        let b = [
            (&hull, Placement { x: 0.0, y: -1.4 }),
            (&hull, Placement { x: 1.7, y: 1.4 }),
        ];
        let ra = multihull_wave_resistance(&a, &cond).unwrap().resistance;
        let rb = multihull_wave_resistance(&b, &cond).unwrap().resistance;
        assert!(
            (ra - rb).abs() <= 1e-9 * ra,
            "staggered {ra} vs mirrored {rb}"
        );
    }

    #[test]
    fn inner_integral_is_translation_invariant_in_modulus() {
        let l = 10.0;
        let hull = wigley(l, 1.0, 0.625).unwrap();
        // Same hull, shifted +37 m in x.
        let s = hull.surface();
        let shifted = hullgeom::bspline::BSplineSurface::new(
            s.degree_x(),
            s.degree_z(),
            s.knots_x().iter().map(|k| k + 37.0).collect(),
            s.knots_z().to_vec(),
            s.control().to_vec(),
        )
        .unwrap();
        let hull2 = Hull::new(shifted).unwrap();
        let nu = 9.80665 / (3.0f64 * 3.0);
        let mut i1 = InnerIntegral::new(&hull, nu, TransomClosure::default());
        let mut i2 = InnerIntegral::new(&hull2, nu, TransomClosure::default());
        for lambda in [1.0, 1.4, 2.5, 6.0] {
            let a = i1.eval(lambda).abs_sq();
            let b = i2.eval(lambda).abs_sq();
            assert!(
                (a - b).abs() <= 1e-9 * a.abs().max(1e-300),
                "lambda={lambda}: {a} vs {b}"
            );
        }
    }
}
