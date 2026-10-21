//! The thin-ship wave kernel on **sectional** hulls
//! ([`hullgeom::sectional`]): the Michell free-wave amplitude and the
//! near-field transforms, built on the hull's per-span interpolant of its
//! depth integrals `Z(x; κ)`, and the resistance and sinkage/trim entry
//! points.
//!
//! In tests, `from_hull` (`crate::hull::FromHull`) builds a sectional hull
//! from a B-spline half-breadth surface — stations at the Greville points of
//! its own x knots, where the interpolant reproduces `Z(x; κ)` exactly —
//! which is the reference harness: any difference from the exact B-spline
//! kernel (`crate::michell::InnerIntegral`, a test oracle) is the depth
//! quadrature's.

use crate::friction::{viscous_resistance_for, ViscousOptions, ViscousResistance};
#[cfg(test)]
use crate::hull::{FromHull, Hull};
use crate::michell::{
    add_transom_step, appendage_source, hollow_shape_moment, run_outer, MemberWave,
    NearFieldKernel, OuterParams, SquatTransforms, TransomClosure, WaveOptions, WaveResistance,
};
use crate::pressure::PressurePatch;
use crate::squat::{integrate_force, DynamicForce, Fleet, Member, SquatOptions};
use crate::MultihullResistance;
use hullgeom::float::{DynamicLoad, DynamicModel, FleetState};
use hullgeom::moments::{osc_moments, C64};
#[cfg(test)]
use hullgeom::sectional::DepthQuadrature;
use hullgeom::sectional::{SectionNodes, SectionalContracted, SectionalHull};
use hullgeom::{Conditions, Error, Placement, Result};
use std::f64::consts::PI;

/// The thin-ship free-wave amplitude of a sectional hull.
pub trait SectionalWave {
    /// The source free-wave amplitude `∬ ∂f/∂x e^{−κz} e^{ik_x(x−x_c)}` at
    /// `λ` (no transom closure).
    fn amplitude(&self, nu: f64, lambda: f64, scratch: &mut SectionalContracted) -> C64;
    /// [`SectionalWave::amplitude`] with the transom closed by `closure`:
    /// the virtual appendage over the hollow behind the transom, its
    /// z-factor the aft end station's own depth integral.
    fn amplitude_closed(
        &self,
        nu: f64,
        lambda: f64,
        closure: TransomClosure,
        scratch: &mut SectionalContracted,
    ) -> C64;
}

impl SectionalWave for SectionalHull {
    /// The source free-wave amplitude `∬ ∂f/∂x e^{−κz} e^{ik_x(x−x_c)}` at
    /// `λ` (the exact B-spline kernel's `InnerIntegral::eval` convention, no transom
    /// closure).
    fn amplitude(&self, nu: f64, lambda: f64, scratch: &mut SectionalContracted) -> C64 {
        self.amplitude_closed(nu, lambda, TransomClosure::None, scratch)
    }

    /// [`SectionalHull::amplitude`] with the transom closed by `closure`: the
    /// virtual appendage `f_v = f_T(z)·φ(s)` over the hollow behind the
    /// transom, exactly as the exact B-spline kernel adds it — its z-factor is the
    /// aft end station's own depth integral.
    fn amplitude_closed(
        &self,
        nu: f64,
        lambda: f64,
        closure: TransomClosure,
        scratch: &mut SectionalContracted,
    ) -> C64 {
        let kx = nu * lambda;
        self.contract(nu * lambda * lambda, scratch);
        if !scratch.any() {
            return C64::ZERO;
        }
        let closing = self.transom_term(kx, nu, closure, scratch.aft());
        let p = self.degree();
        let mut xm = Vec::with_capacity(p + 1);
        let mut f = C64::ZERO;
        for (s, sx) in self.spans().iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut xm);
            let phase = C64::cis(kx * (sx.start - self.x_center()));
            let mut sum = C64::ZERO;
            for (a, &m) in xm.iter().enumerate() {
                sum = sum + m.scale(scratch.g_fx()[s * p + a]);
            }
            f = f + phase * sum;
        }
        f + closing + self.bow_source(kx, scratch.fore())
    }
}

/// The near-field transforms of a sectional hull, for [`crate::squat`] and
/// [`crate::nearfield`] (kernel convention; see [`SquatTransforms`]).
pub(crate) trait NearFieldTransforms {
    fn bow_source(&self, kx: f64, z_b: f64) -> C64;
    fn transom_term(&self, kx: f64, nu: f64, closure: TransomClosure, z_t: f64) -> C64;
    fn add_transom_transforms(
        &self,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
        z_t: f64,
        z_b: f64,
        t: &mut SquatTransforms,
    );
    #[allow(
        dead_code,
        reason = "the open-transom form, exercised by the parity tests"
    )]
    fn transforms_at(&self, zc: &SectionalContracted, kx: f64) -> SquatTransforms;
    fn transforms_at_closed(
        &self,
        zc: &SectionalContracted,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
    ) -> SquatTransforms;
    fn q_closed(
        &self,
        zc: &SectionalContracted,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
        xm: &mut Vec<C64>,
    ) -> C64;
}

impl NearFieldTransforms for SectionalHull {
    /// The fore end's step in the sources (kernel convention): a hull that
    /// ends forward on a finite section — a plumb stem face, or the thin
    /// wide sliver it leaves when trimmed bow-up — has `f` drop to 0 there,
    /// so `∂f/∂x` carries `−f_b δ(x − x_b)`. The face is real hull, so
    /// unlike a transom it is never left open. Zero for a pointed bow.
    fn bow_source(&self, kx: f64, z_b: f64) -> C64 {
        let Some(&x_b) = self.station_xs().last() else {
            return C64::ZERO;
        };
        C64::ZERO - C64::cis(kx * (x_b - self.x_center())).scale(z_b)
    }

    /// The transom appendage's free-wave amplitude (kernel convention): the
    /// exact `InnerIntegral::transom_term` with the transom section's
    /// z-factor `z_t` taken from the aft end station.
    fn transom_term(&self, kx: f64, nu: f64, closure: TransomClosure, z_t: f64) -> C64 {
        let Some(tr) = self.transom() else {
            return C64::ZERO;
        };
        let Some(lv) = closure.hollow_length(tr.depth, nu) else {
            return C64::ZERO;
        };
        let mut sm = Vec::with_capacity(3);
        let shape = hollow_shape_moment(-kx * lv, &mut sm);
        let phase = C64::cis(kx * (tr.x - self.x_center()));
        C64::ZERO - (phase * shape).scale(z_t)
    }

    /// The transom's share of the near-field transforms (kernel convention),
    /// as the exact `InnerIntegral::add_transom_transforms`: the step in the
    /// weights, the closing appendage in the sources (see
    /// [`SquatTransforms`]).
    fn add_transom_transforms(
        &self,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
        z_t: f64,
        z_b: f64,
        t: &mut SquatTransforms,
    ) {
        let (Some(&x_a), Some(&x_b)) = (self.station_xs().first(), self.station_xs().last()) else {
            return;
        };
        // The weights: the hull ends where its end stations are, on whatever
        // section they carry — a transom aft, a plumb stem or its trimmed
        // sliver forward, nothing at a pointed end. Taking the steps from the
        // end stations themselves (not only a detected transom) keeps the
        // force continuous as an end section shrinks to nothing.
        let wl = |s: &SectionNodes| s.waterline();
        let (wl_a, wl_b) = (
            self.section_nodes().first().map_or(0.0, wl),
            self.section_nodes().last().map_or(0.0, wl),
        );
        let dx_a = x_a - self.x_center();
        let phase_a = C64::cis(kx * dx_a);
        add_transom_step(phase_a, dx_a, z_t, wl_a, t);
        let dx_b = x_b - self.x_center();
        add_transom_step(C64::cis(kx * dx_b), dx_b, -z_b, -wl_b, t);
        // The sources: the bow's step always (the stem face is hull), the
        // transom's closing appendage when the closure applies.
        t.q_src = t.q_src + self.bow_source(kx, z_b);
        let Some(tr) = self.transom() else {
            return;
        };
        let Some(lv) = closure.hollow_length(tr.depth, nu) else {
            return;
        };
        let mut m = Vec::with_capacity(4);
        osc_moments(-kx * lv, 1.0, 3, &mut m);
        let phase = C64::cis(kx * (tr.x - self.x_center()));
        t.q_src = t.q_src + appendage_source(phase, &m).scale(z_t);
    }

    /// The six near-field transforms at `k_x` from a contraction at some `κ`
    /// — the same quantities, convention and x origin as the exact B-spline kernel's
    /// `InnerIntegral::transforms_at` (no transom closure).
    #[allow(
        dead_code,
        reason = "the open-transom form, exercised by the parity tests"
    )]
    fn transforms_at(&self, zc: &SectionalContracted, kx: f64) -> SquatTransforms {
        self.transforms_at_closed(zc, kx, 0.0, TransomClosure::None)
    }

    /// [`SectionalHull::transforms_at`] with the transom closed by `closure`
    /// at `ν` (the hollow length depends on it).
    fn transforms_at_closed(
        &self,
        zc: &SectionalContracted,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
    ) -> SquatTransforms {
        let p = self.degree();
        let mut t = SquatTransforms::default();
        if !zc.any() {
            return t;
        }
        let mut xm = Vec::with_capacity(p + 1);
        for (s, sx) in self.spans().iter().enumerate() {
            osc_moments(kx, sx.len, p, &mut xm);
            let d = sx.start - self.x_center();
            let phase = C64::cis(kx * d);
            let (mut q_s, mut q1_s, mut w_s, mut q1w_s) =
                (C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO);
            for a in 0..p {
                let (g, gw) = (
                    zc.g_fx()[s * p + a],
                    self.waterline_slope_power()[s * p + a],
                );
                q_s = q_s + xm[a].scale(g);
                w_s = w_s + xm[a].scale(gw);
                let shifted = xm[a + 1] + xm[a].scale(d);
                q1_s = q1_s + shifted.scale(g);
                q1w_s = q1w_s + shifted.scale(gw);
            }
            let (mut p_s, mut pw_s) = (C64::ZERO, C64::ZERO);
            for a in 0..=p {
                p_s = p_s + xm[a].scale(zc.g_f()[s * (p + 1) + a]);
                pw_s = pw_s + xm[a].scale(self.waterline_power()[s * (p + 1) + a]);
            }
            t.q = t.q + phase * q_s;
            t.q1 = t.q1 + phase * q1_s;
            t.w = t.w + phase * w_s;
            t.q1_wl = t.q1_wl + phase * q1w_s;
            t.p = t.p + phase * p_s;
            t.p_wl = t.p_wl + phase * pw_s;
        }
        t.q_src = t.q;
        self.add_transom_transforms(kx, nu, closure, zc.aft(), zc.fore(), &mut t);
        let conj = |v: C64| C64::new(v.re, -v.im);
        t.q_src = conj(t.q_src);
        t.q = conj(t.q);
        t.p = conj(t.p);
        t.q1 = conj(t.q1);
        t.w = conj(t.w);
        t.p_wl = conj(t.p_wl);
        t.q1_wl = conj(t.q1_wl);
        t
    }

    /// Just the source transform `q` of [`SectionalHull::transforms_at_closed`]
    /// (same convention), for the near field, which needs nothing else.
    fn q_closed(
        &self,
        zc: &SectionalContracted,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
        xm: &mut Vec<C64>,
    ) -> C64 {
        if !zc.any() {
            return C64::ZERO;
        }
        let p = self.degree();
        let mut q = C64::ZERO;
        for (s, sx) in self.spans().iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, xm);
            let mut q_s = C64::ZERO;
            for (a, &m) in xm.iter().enumerate().take(p) {
                q_s = q_s + m.scale(zc.g_fx()[s * p + a]);
            }
            q = q + C64::cis(kx * (sx.start - self.x_center())) * q_s;
        }
        if let Some(tr) = self.transom() {
            if let Some(lv) = closure.hollow_length(tr.depth, nu) {
                let mut m = Vec::with_capacity(4);
                osc_moments(-kx * lv, 1.0, 3, &mut m);
                let shape_dx = m[2].scale(6.0) - m[1].scale(6.0);
                let phase = C64::cis(kx * (tr.x - self.x_center()));
                q = q - (phase * shape_dx).scale(zc.aft());
            }
        }
        q = q + self.bow_source(kx, zc.fore());
        C64::new(q.re, -q.im)
    }
}

/// The sectional hull as a near-field kernel, for [`crate::squat`].
#[derive(Clone)]
pub(crate) struct SectionalKernel<'h> {
    hull: &'h SectionalHull,
    nu: f64,
    closure: TransomClosure,
}

impl NearFieldKernel for SectionalKernel<'_> {
    type Contracted = SectionalContracted;
    fn contract_z(&mut self, kappa: f64, out: &mut SectionalContracted) {
        self.hull.contract(kappa, out)
    }
    fn transforms_at(&mut self, zc: &SectionalContracted, kx: f64) -> SquatTransforms {
        self.hull
            .transforms_at_closed(zc, kx, self.nu, self.closure)
    }
}

/// A hull's whole free-wave amplitude: its thickness, closed by `closure`,
/// and the pressure of its lift (see [`crate::pressure`]), if it carries
/// one; `scale` is `1/(2ρU²)`.
fn amplitude_lifted(
    hull: &SectionalHull,
    patch: Option<&PressurePatch>,
    nu: f64,
    lambda: f64,
    closure: TransomClosure,
    scale: f64,
    scratch: &mut SectionalContracted,
) -> C64 {
    let f = hull.amplitude_closed(nu, lambda, closure, scratch);
    match patch {
        Some(p) => f + p.amplitude(nu, lambda, hull.x_center(), scale),
        None => f,
    }
}

/// The sectional hull as a wave-resistance member, for the shared outer
/// quadrature, with its placement phase offsets.
#[derive(Clone)]
struct SectionalMember<'h> {
    hull: &'h SectionalHull,
    patch: Option<&'h PressurePatch>,
    scale: f64,
    scratch: SectionalContracted,
    dx: f64,
    dy: f64,
    closure: TransomClosure,
}

impl MemberWave for SectionalMember<'_> {
    fn amps(&mut self, nu: f64, lambda: f64) -> (C64, C64) {
        let f = amplitude_lifted(
            self.hull,
            self.patch,
            nu,
            lambda,
            self.closure,
            self.scale,
            &mut self.scratch,
        );
        (f, f)
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

/// One hull's free-wave amplitude `A(λ)` on `[1, λ_max]`, tabulated on
/// panels of Chebyshev points and read back by barycentric interpolation.
/// `A` is smooth in `λ` — it oscillates at the hull's own half-length phase,
/// ~`ν L/2` per unit `λ`, under a decaying envelope — so panels a few
/// radians of that phase wide interpolate it spectrally.
struct AmplitudeTable {
    lo: f64,
    h: f64,
    lam_max: f64,
    /// `NODES` values per panel.
    vals: Vec<C64>,
}

const NODES: usize = 17;

/// Chebyshev (second-kind) points on `[-1, 1]`, ascending.
fn cheb(j: usize) -> f64 {
    -(PI * j as f64 / (NODES - 1) as f64).cos()
}

impl AmplitudeTable {
    fn new(
        hull: &SectionalHull,
        patch: Option<&PressurePatch>,
        nu: f64,
        closure: TransomClosure,
        scale: f64,
        lam_max: f64,
    ) -> Self {
        let lo = 1.0;
        // Panels ~3 radians of the amplitude's own phase wide.
        let h = (3.0 / (nu * 0.5 * hull.length()).max(1e-9)).min(0.5);
        let panels = (((lam_max - lo) / h).ceil() as usize).max(1);
        let per_panel =
            hullgeom::parallel::map_indexed(panels, SectionalContracted::default, |scratch, p| {
                (0..NODES)
                    .map(|j| {
                        let lam = lo + h * (p as f64 + 0.5 * (1.0 + cheb(j)));
                        amplitude_lifted(hull, patch, nu, lam, closure, scale, scratch)
                    })
                    .collect::<Vec<_>>()
            });
        AmplitudeTable {
            lo,
            h,
            lam_max: lo + h * panels as f64,
            vals: per_panel.into_iter().flatten().collect(),
        }
    }

    fn eval(&self, lam: f64) -> C64 {
        if !(lam >= self.lo && lam < self.lam_max) {
            return C64::ZERO;
        }
        let u = (lam - self.lo) / self.h;
        let p = (u.floor() as usize).min(self.vals.len() / NODES - 1);
        let t = 2.0 * (u - p as f64) - 1.0;
        let v = &self.vals[p * NODES..(p + 1) * NODES];
        let (mut num, mut den) = (C64::ZERO, 0.0);
        for (j, &vj) in v.iter().enumerate() {
            let d = t - cheb(j);
            if d == 0.0 {
                return vj;
            }
            let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
            let w = sign * if j == 0 || j == NODES - 1 { 0.5 } else { 1.0 } / d;
            num = num + vj.scale(w);
            den += w;
        }
        num.scale(1.0 / den)
    }
}

/// A fleet member reading its amplitude from a shared table.
#[derive(Clone)]
struct Tabulated<'t> {
    table: &'t AmplitudeTable,
    dx: f64,
    dy: f64,
}

impl MemberWave for Tabulated<'_> {
    fn amps(&mut self, _nu: f64, lambda: f64) -> (C64, C64) {
        let f = self.table.eval(lambda);
        (f, f)
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

/// Michell wave resistance of a single sectional hull: an adaptive outer
/// quadrature over `λ`, the transom (if any) closed by `opts.transom`.
pub fn wave_resistance(
    hull: &SectionalHull,
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    multihull_wave_resistance(&[(hull, Placement::default())], cond, opts)
}

/// Wave resistance of a fleet of sectional hulls, placements and interference
/// included: each member's amplitude carries its placement phase
/// `exp(i ν (λ Δx_j ± λ √(λ²−1) y_j))`, and the integrand is the mean of the
/// two (±θ) wave systems. For two identical hulls separated by `s` this
/// reduces to the classical catamaran factor `4 cos²(½ ν s λ √(λ²−1))`.
pub fn multihull_wave_resistance(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    floored_wave_resistance(members, &[], cond, opts, 0.0)
}

/// [`multihull_wave_resistance`] with each member's lift as a pressure
/// patch beside its thickness (`patches[i]` for member `i`; missing or
/// `None`: no lift).
pub fn multihull_wave_resistance_lifted(
    members: &[(&SectionalHull, Placement)],
    patches: &[Option<&PressurePatch>],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    floored_wave_resistance(members, patches, cond, opts, 0.0)
}

/// [`multihull_wave_resistance`], resolved only to `floor` [N]: the
/// truncation and refinement are relative to the larger of it and the
/// result. A member's own resistance, wanted beside a fleet's, needs only
/// the fleet's precision; a small part's tail (a surface-piercing leg's,
/// say) otherwise runs to wavelengths far below the theory's reach.
fn floored_wave_resistance(
    members: &[(&SectionalHull, Placement)],
    patches: &[Option<&PressurePatch>],
    cond: &Conditions,
    opts: &WaveOptions,
    floor: f64,
) -> Result<WaveResistance> {
    cond.validate()?;
    if members.is_empty() {
        return Err(Error::InvalidConditions("empty fleet".into()));
    }
    let (u, g) = (cond.speed, cond.gravity);
    let nu = g / (u * u);
    let n = members.len() as f64;
    let cx_ref = members.iter().map(|(h, p)| h.x_center() + p.x).sum::<f64>() / n;
    let y_ref = members.iter().map(|(_, p)| p.y).sum::<f64>() / n;
    let coeff = 4.0 * cond.fluid.density * g * g / (PI * u * u);
    let scale = 1.0 / (2.0 * cond.fluid.density * u * u);
    let patch = |i: usize| patches.get(i).copied().flatten();
    let params = OuterParams {
        nu,
        x_half: members
            .iter()
            .map(|(h, p)| (h.x_center() + p.x - cx_ref).abs() + 0.5 * h.length())
            .fold(0.0, f64::max),
        // A patch's own width oscillates in θ as a separation would.
        y_half: members
            .iter()
            .enumerate()
            .map(|(i, (_, p))| (p.y - y_ref).abs() + patch(i).map_or(0.0, |q| q.half_width()))
            .fold(0.0, f64::max),
        t_max: members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
        floor: floor / coeff,
    };
    // Copies of members side by side (a catamaran's demihulls, each with
    // its drive): the fleet's integrand oscillates with the spacing and
    // needs many λ-nodes, but each copy carries its original's smooth
    // amplitude. Integrate one copy of each, side by side collapsed (which
    // finds how far in λ the fleet matters, without the oscillation),
    // tabulate each original's amplitude that far, and let the fleet
    // integral interpolate the tables — exact to the interpolation's ~1e-10.
    let tw = twins(members);
    if tw.iter().any(Option::is_some) && members.iter().any(|(_, p)| p.y != members[0].1.y) {
        let originals: Vec<usize> = (0..members.len()).filter(|&i| tw[i].is_none()).collect();
        let collapsed: Vec<(&SectionalHull, Placement)> = originals
            .iter()
            .map(|&i| {
                (
                    members[i].0,
                    Placement {
                        y: 0.0,
                        ..members[i].1
                    },
                )
            })
            .collect();
        let collapsed_patches: Vec<Option<&PressurePatch>> =
            originals.iter().map(|&i| patch(i)).collect();
        let group = floored_wave_resistance(&collapsed, &collapsed_patches, cond, opts, floor)?;
        let tables: Vec<AmplitudeTable> = originals
            .iter()
            .map(|&i| {
                AmplitudeTable::new(
                    members[i].0,
                    patch(i),
                    nu,
                    opts.transom,
                    scale,
                    1.1 * group.max_lambda,
                )
            })
            .collect();
        let table_of = |i: usize| {
            let o = tw[i].unwrap_or(i);
            &tables[originals.iter().position(|&j| j == o).expect("an original")]
        };
        let mem = members
            .iter()
            .enumerate()
            .map(|(i, (h, p))| Tabulated {
                table: table_of(i),
                dx: h.x_center() + p.x - cx_ref,
                dy: p.y - y_ref,
            })
            .collect();
        return Ok(run_outer(&params, opts, coeff, mem));
    }
    let mem = members
        .iter()
        .enumerate()
        .map(|(i, (h, p))| SectionalMember {
            hull: h,
            patch: patch(i),
            scale,
            scratch: SectionalContracted::default(),
            dx: h.x_center() + p.x - cx_ref,
            dy: p.y - y_ref,
            closure: opts.transom,
        })
        .collect();
    Ok(run_outer(&params, opts, coeff, mem))
}

/// Wave + viscous resistance of a fleet of sectional hulls, broken down into
/// [`MultihullResistance`]. The viscous part uses each hull's own shell area
/// ([`SectionalHull::wetted_surface`]).
pub fn multihull_resistance(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    wave_opts: &WaveOptions,
    viscous_opts: &ViscousOptions,
) -> Result<MultihullResistance> {
    multihull_resistance_lifted(members, &[], cond, wave_opts, viscous_opts)
}

/// [`multihull_resistance`] with each member's lift as a pressure patch
/// beside its thickness (see [`multihull_wave_resistance_lifted`]).
pub fn multihull_resistance_lifted(
    members: &[(&SectionalHull, Placement)],
    patches: &[Option<&PressurePatch>],
    cond: &Conditions,
    wave_opts: &WaveOptions,
    viscous_opts: &ViscousOptions,
) -> Result<MultihullResistance> {
    let wave = floored_wave_resistance(members, patches, cond, wave_opts, 0.0)?;
    // A copy of an earlier member has its resistance alone.
    let tw = twins(members);
    let mut solo: Vec<f64> = Vec::with_capacity(members.len());
    // A member's own resistance needn't be resolved finer than a tenth of
    // the fleet's: that leaves a hull's (always a good share of it) exact,
    // and stops a small part's (a drive's leg) long before its tail runs
    // out to wavelengths far below the theory's reach.
    let floor = 0.1 * wave.resistance.abs();
    for (i, m) in members.iter().enumerate() {
        let r = match tw[i] {
            Some(j) => solo[j],
            None => {
                let own = [patches.get(i).copied().flatten()];
                floored_wave_resistance(&[*m], &own, cond, wave_opts, floor)?.resistance
            }
        };
        solo.push(r);
    }
    let solo_wave_total: f64 = solo.iter().sum();
    let viscous: Vec<ViscousResistance> = members
        .iter()
        .map(|(h, _)| viscous_resistance_for(h.length(), h.wetted_surface(), cond, viscous_opts))
        .collect::<Result<_>>()?;
    let viscous_total: f64 = viscous.iter().map(|v| v.resistance).sum();
    let wetted_surface: f64 = members.iter().map(|(h, _)| h.wetted_surface()).sum();
    let total = wave.resistance + viscous_total;
    let q = 0.5 * cond.fluid.density * cond.speed * cond.speed * wetted_surface;
    Ok(MultihullResistance {
        wave,
        viscous,
        viscous_total,
        total,
        effective_power: total * cond.speed,
        wetted_surface,
        cw: wave.resistance / q,
        cv: viscous_total / q,
        ct: total / q,
        solo_wave_total,
        interference: if solo_wave_total > f64::MIN_POSITIVE {
            wave.resistance / solo_wave_total
        } else {
            1.0
        },
    })
}

/// Near-field vertical force and pitch moment on a single sectional hull
/// about `x_ref` (see [`crate::squat`] for the theory), the transom closed by
/// `opts.wave.transom`.
pub fn dynamic_force(
    hull: &SectionalHull,
    cond: &Conditions,
    x_ref: f64,
    opts: &SquatOptions,
) -> Result<DynamicForce> {
    multihull_dynamic_force(&[(hull, Placement::default())], cond, x_ref, opts)
}

/// Near-field force and moment on a fleet of sectional hulls, demihull
/// interaction included: every member's pressure includes what every other
/// member's source sheet induces on it, through the placement phases.
/// `x_ref` is the pitch pivot in fleet coordinates.
pub fn multihull_dynamic_force(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    x_ref: f64,
    opts: &SquatOptions,
) -> Result<DynamicForce> {
    cond.validate()?;
    if members.is_empty() {
        return Err(Error::InvalidGeometry("empty fleet".into()));
    }
    let nu = cond.gravity / (cond.speed * cond.speed);
    let fleet = Fleet {
        members: members
            .iter()
            .map(|(h, p)| Member {
                inner: SectionalKernel {
                    hull: h,
                    nu,
                    closure: opts.wave.transom,
                },
                cx: h.x_center() + p.x,
                y: p.y,
            })
            .collect(),
        nu,
        x_ref,
        scratch: vec![SquatTransforms::default(); members.len()],
        zc_tmp: vec![SectionalContracted::default(); members.len()],
        twin: twins(members),
    };
    let l_max = members.iter().map(|(h, _)| h.length()).fold(0.0, f64::max);
    let t_max = members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max);
    let volume = members.iter().map(|(h, _)| h.displaced_volume()).sum();
    Ok(integrate_force(fleet, cond, l_max, t_max, volume, opts))
}

/// For each member, an earlier member carrying the very same hull (a clone:
/// bit-identical hydrostatics and stations), which the force quadrature can
/// then transform once for both.
fn twins(members: &[(&SectionalHull, Placement)]) -> Vec<Option<usize>> {
    let print = |h: &SectionalHull| {
        [
            h.displaced_volume(),
            h.lcb_x(),
            h.length(),
            h.x_center(),
            h.wetted_surface(),
            h.draft(),
        ]
        .map(f64::to_bits)
        .into_iter()
        .chain([h.stations() as u64])
        .collect::<Vec<_>>()
    };
    let prints: Vec<_> = members.iter().map(|(h, _)| print(h)).collect();
    (0..members.len())
        .map(|i| (0..i).find(|&j| prints[j] == prints[i]))
        .collect()
}

/// The dynamic load a sectional fleet carries at speed, as the equilibrium
/// solver's model (see [`DynamicModel`]): the near-field
/// force and moment about `x_ref`, with a single-pass quadrature for the
/// solver's slope probes — within ~0.1% of the converged value, and 5–20×
/// cheaper.
pub fn dynamic_load_closure<'a>(
    cond: &'a Conditions,
    x_ref: f64,
    opts: &'a SquatOptions,
) -> SectionalDynamic<'a> {
    SectionalDynamic { cond, x_ref, opts }
}

/// See [`dynamic_load_closure`].
#[derive(Clone, Copy)]
pub struct SectionalDynamic<'a> {
    cond: &'a Conditions,
    x_ref: f64,
    opts: &'a SquatOptions,
}

impl SectionalDynamic<'_> {
    fn eval(&self, fleet: &FleetState<SectionalHull>, opts: &SquatOptions) -> Result<DynamicLoad> {
        if fleet.members.is_empty() {
            return Ok(DynamicLoad::default());
        }
        let members: Vec<(&SectionalHull, Placement)> =
            fleet.members.iter().map(|(h, p)| (h, *p)).collect();
        let d = multihull_dynamic_force(&members, self.cond, self.x_ref, opts)?;
        Ok(DynamicLoad {
            force_up: d.force_up,
            moment_bow_up: d.moment_bow_up,
        })
    }
}

impl DynamicModel<SectionalHull> for SectionalDynamic<'_> {
    fn load(&mut self, fleet: &FleetState<SectionalHull>) -> Result<DynamicLoad> {
        self.eval(fleet, self.opts)
    }

    fn probe(&mut self, fleet: &FleetState<SectionalHull>) -> Option<Result<DynamicLoad>> {
        let coarse = SquatOptions {
            max_refinements: 0,
            ..*self.opts
        };
        Some(self.eval(fleet, &coarse))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::michell::{InnerIntegral, TransomClosure};
    use crate::squat::SquatOptions;

    /// Largest relative disagreement of the six transforms over a (k_x, κ)
    /// grid spanning the near-field quadrature's range, each transform
    /// against its own largest magnitude on the grid.
    fn transform_error(hull: &Hull, sec: &SectionalHull, nu: f64) -> f64 {
        let mut inner = InnerIntegral::new(hull, nu, TransomClosure::None);
        let mut zc = SectionalContracted::default();
        let parts = |t: &SquatTransforms| [t.q, t.p, t.q1, t.w, t.p_wl, t.q1_wl];
        let (mut worst_abs, mut scale) = ([0.0f64; 6], [0.0f64; 6]);
        for ik in 0..40 {
            let kappa = nu * 1e-3 * 10f64.powf(ik as f64 * 0.15); // up to ~1e3·ν
            sec.contract(kappa, &mut zc);
            for ix in 0..40 {
                let kx = nu * (ix as f64 * 0.5 - 10.0);
                let a = parts(&inner.eval_transforms(kx, kappa));
                let b = parts(&sec.transforms_at(&zc, kx));
                for k in 0..6 {
                    worst_abs[k] = worst_abs[k].max((a[k] - b[k]).abs());
                    scale[k] = scale[k].max(a[k].abs());
                }
            }
        }
        (0..6).map(|k| worst_abs[k] / scale[k]).fold(0.0, f64::max)
    }

    fn amplitude_error(hull: &Hull, sec: &SectionalHull, nu: f64) -> f64 {
        let mut inner = InnerIntegral::new(hull, nu, TransomClosure::None);
        let mut zc = SectionalContracted::default();
        let (mut worst, mut scale) = (0.0f64, 0.0f64);
        for i in 0..400 {
            let lambda = 1.0 + i as f64 * 0.25; // κ up to ~1e4·ν
            let a = inner.eval(lambda);
            let b = sec.amplitude(nu, lambda, &mut zc);
            worst = worst.max((a - b).abs());
            scale = scale.max(a.abs());
        }
        worst / scale
    }

    #[test]
    fn reproduces_the_lofted_kernel_on_wigley() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        for nu in [0.3, 1.1, 4.0] {
            let (t, a) = (
                transform_error(&hull, &sec, nu),
                amplitude_error(&hull, &sec, nu),
            );
            assert!(
                t < 1e-9 && a < 1e-9,
                "ν {nu}: transforms {t:.1e}, amplitude {a:.1e}"
            );
        }
    }

    #[test]
    #[ignore = "timing report, not a check"]
    fn cost_against_the_lofted_kernel() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let nu = 9.81 / 25.0;
        let mut inner = InnerIntegral::new(&hull, nu, TransomClosure::None);
        let n = 2000;
        let mut zl = crate::michell::ZContracted::default();
        let t = std::time::Instant::now();
        for i in 0..n {
            inner.contract_z(nu * (1.0 + i as f64 * 0.01), &mut zl);
        }
        let lofted_k = t.elapsed().as_secs_f64() / n as f64;
        let mut zs = SectionalContracted::default();
        let t = std::time::Instant::now();
        for i in 0..n {
            sec.contract(nu * (1.0 + i as f64 * 0.01), &mut zs);
        }
        let sec_k = t.elapsed().as_secs_f64() / n as f64;
        inner.contract_z(nu * 2.0, &mut zl);
        sec.contract(nu * 2.0, &mut zs);
        let m = 20_000;
        let mut acc = C64::ZERO;
        let t = std::time::Instant::now();
        for i in 0..m {
            acc = acc + inner.transforms_at(&zl, nu * (1.0 + i as f64 * 1e-3)).q;
        }
        let lofted_x = t.elapsed().as_secs_f64() / m as f64;
        let t = std::time::Instant::now();
        for i in 0..m {
            acc = acc + sec.transforms_at(&zs, nu * (1.0 + i as f64 * 1e-3)).q;
        }
        let sec_x = t.elapsed().as_secs_f64() / m as f64;
        eprintln!(
            "lofted {}x{} spans: {:.1} us/κ, {:.1} us/kx\n\
             sectional {} stations, {} depth nodes, {} x-spans: {:.1} us/κ, {:.1} us/kx ({acc:?})",
            hull.spans_x().len(),
            hull.spans_z().len(),
            lofted_k * 1e6,
            lofted_x * 1e6,
            sec.stations(),
            sec.depth_nodes(),
            sec.x_spans(),
            sec_k * 1e6,
            sec_x * 1e6
        );
    }

    /// Speeds spanning Fn ≈ 0.15–0.5 for a hull of length `l`.
    fn speeds(l: f64) -> Vec<f64> {
        [0.15, 0.25, 0.35, 0.5]
            .iter()
            .map(|fn_| fn_ * (9.81 * l).sqrt())
            .collect()
    }

    fn untransomed() -> WaveOptions {
        WaveOptions {
            transom: TransomClosure::None,
            ..WaveOptions::default()
        }
    }

    fn compare_end_to_end(hull: &Hull, name: &str, tol: f64) {
        let sec = SectionalHull::from_hull(hull, &DepthQuadrature::default()).unwrap();
        let wave = untransomed();
        let squat = SquatOptions {
            wave,
            ..SquatOptions::default()
        };
        let x_ref = hull.x_center();
        for u in speeds(hull.length()) {
            let cond = Conditions::seawater(u);
            let t = std::time::Instant::now();
            let a = crate::michell::wave_resistance_with(hull, &cond, &wave).unwrap();
            let fa = crate::squat::dynamic_force(hull, &cond, x_ref, &squat).unwrap();
            let t_lofted = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let b = wave_resistance(&sec, &cond, &wave).unwrap();
            let fb = dynamic_force(&sec, &cond, x_ref, &squat).unwrap();
            let t_sec = t.elapsed().as_secs_f64();
            let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(1e-300);
            let (er, ef, em) = (
                rel(a.resistance, b.resistance),
                rel(fa.force_up, fb.force_up),
                rel(fa.moment_bow_up, fb.moment_bow_up),
            );
            eprintln!(
                "{name} U {u:.2}: Rw {:.4e}/{:.4e} ({er:.1e}), Fz {:.4e} ({ef:.1e}), M {:.4e} ({em:.1e}); \
                 {t_lofted:.3} s lofted, {t_sec:.3} s sectional",
                a.resistance, b.resistance, fa.force_up, fa.moment_bow_up
            );
            assert!(er < tol && ef < tol && em < tol, "{name} U {u}");
        }
        let vol = rel_vol(hull, &sec);
        assert!(vol < 1e-10, "{name} volume {vol:.1e}");
    }

    fn rel_vol(hull: &Hull, sec: &SectionalHull) -> f64 {
        (hull.displaced_volume() - sec.displaced_volume()).abs() / hull.displaced_volume()
    }

    #[test]
    fn wave_resistance_and_squat_match_the_lofted_hull_on_wigley() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        compare_end_to_end(&hull, "wigley", 1e-7);
    }

    /// The sectional importer end to end on geometry with an exact answer:
    /// the Wigley hull written out as CAD patches and cut back into
    /// sections. Its sections are smooth all the way to the keel and its
    /// depth integral is quadratic in x, so the rays and the x interpolant
    /// should both be essentially exact.
    #[test]
    fn iges_sections_of_a_wigley_reproduce_the_exact_hull() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let surfs = hullgeom::iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let text = hullgeom::iges::write(&surfs, "wigley").unwrap();
        let imp =
            hullgeom::iges::import_sectional(&text, &hullgeom::iges::SectionalOptions::default())
                .unwrap()
                .hulls
                .remove(0);
        eprintln!("{:?}", imp.report);
        let sec = &imp.hull;
        let vol = rel_vol(&hull, sec);
        let wave = untransomed();
        let squat = SquatOptions {
            wave,
            ..SquatOptions::default()
        };
        let rel = |x: f64, y: f64| (x - y).abs() / x.abs();
        for u in speeds(hull.length()) {
            let cond = Conditions::seawater(u);
            let a = crate::michell::wave_resistance_with(&hull, &cond, &wave).unwrap();
            let b = wave_resistance(sec, &cond, &wave).unwrap();
            let fa = crate::squat::dynamic_force(&hull, &cond, 0.0, &squat).unwrap();
            let fb = dynamic_force(sec, &cond, 0.0, &squat).unwrap();
            let (er, ef, em) = (
                rel(a.resistance, b.resistance),
                rel(fa.force_up, fb.force_up),
                rel(fa.moment_bow_up, fb.moment_bow_up),
            );
            eprintln!("U {u:.2}: Rw {er:.1e}, Fz {ef:.1e}, M {em:.1e}");
            assert!(
                er < 1e-8 && ef < 1e-9 && em < 1e-7,
                "U {u}: {er:.1e} {ef:.1e} {em:.1e}"
            );
        }
        eprintln!("volume {vol:.1e}");
        assert!(vol < 1e-8);
    }

    /// A real CAD hull with every awkward end at once — a wet transom whose
    /// side skins run on past it, a separately-patched stem, a fine entry —
    /// imported by sections at two resolutions: the transom is carried as a
    /// section (not closed to nothing), and resistance has converged.
    #[test]
    fn cad_sections_converge_and_keep_the_transom() {
        let Some(text) = crate::cad_fixture("e12.igs") else {
            return;
        };
        let wave = untransomed();
        let import = |stations: usize, rays: usize| {
            let so = hullgeom::iges::SectionalOptions {
                waterline_z: -0.95,
                stations,
                rays,
                ..Default::default()
            };
            let mut fleet = hullgeom::iges::import_sectional(&text, &so).unwrap();
            assert!(fleet.failed.is_empty(), "{:?}", fleet.failed);
            assert_eq!(fleet.hulls.len(), 1, "the stem is part of the hull");
            fleet.hulls.remove(0)
        };
        let (coarse, fine) = (import(61, 17), import(121, 33));
        // The aft end station is the transom, 8–9 mm deep and 0.37 m wide.
        let (x_aft, aft) = &fine.sections[0];
        assert!((x_aft - 0.1799).abs() < 1e-3, "aft end at {x_aft}");
        let wl_half = aft.first().unwrap().0;
        let depth = aft.last().unwrap().1;
        assert!(
            wl_half > 0.18 && depth > 0.008,
            "transom {wl_half} x {depth}"
        );
        let rel = (coarse.hull.displaced_volume() - fine.hull.displaced_volume()).abs()
            / fine.hull.displaced_volume();
        assert!(rel < 1e-4, "volume {rel:.1e}");
        let l = fine.hull.length();
        for fnum in [0.2, 0.3] {
            let cond = Conditions::seawater(fnum * (9.81 * l).sqrt());
            let a = wave_resistance(&coarse.hull, &cond, &wave)
                .unwrap()
                .resistance;
            let b = wave_resistance(&fine.hull, &cond, &wave)
                .unwrap()
                .resistance;
            assert!((a - b).abs() / b < 2e-3, "Fn {fnum}: {a} vs {b}");
        }
    }

    /// Hydrostatics from sections against the exact Wigley's: the
    /// waterplane and buoyancy come from the interpolants exactly; the
    /// wetted surface from strips between station curves converges to the
    /// graph's true area.
    #[test]
    fn sectional_hydrostatics_match_the_exact_hull() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let surfs = hullgeom::iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let text = hullgeom::iges::write(&surfs, "wigley").unwrap();
        let sec =
            hullgeom::iges::import_sectional(&text, &hullgeom::iges::SectionalOptions::default())
                .unwrap()
                .hulls
                .remove(0)
                .hull;
        let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
        let area = rel(sec.waterplane_area(), hull.waterplane_area());
        let m2 = rel(
            sec.waterplane_second_moment(),
            hull.waterplane_second_moment(),
        );
        let m1 = (sec.waterplane_moment() - hull.waterplane_moment()).abs();
        let lcb = (sec.lcb_x() - hull.lcb_x()).abs();
        let s = rel(sec.wetted_surface(), hull.wetted_surface());
        eprintln!("A_w {area:.1e}, I {m2:.1e}, M {m1:.1e}, lcb {lcb:.1e}, S {s:.1e}");
        assert!(
            area < 1e-9 && m2 < 1e-9 && m1 < 1e-9 && lcb < 1e-9,
            "{area} {m2} {m1} {lcb}"
        );
        assert!(s < 1e-3, "wetted surface {s:.1e}");
    }
}

#[cfg(test)]
mod twin_tests {
    use super::*;

    /// The direct fleet integral: every member evaluates its own amplitude.
    fn direct(
        members: &[(&SectionalHull, Placement)],
        cond: &Conditions,
        opts: &WaveOptions,
    ) -> f64 {
        let (u, g) = (cond.speed, cond.gravity);
        let nu = g / (u * u);
        let n = members.len() as f64;
        let cx_ref = members.iter().map(|(h, p)| h.x_center() + p.x).sum::<f64>() / n;
        let y_ref = members.iter().map(|(_, p)| p.y).sum::<f64>() / n;
        let params = OuterParams {
            floor: 0.0,
            nu,
            x_half: members
                .iter()
                .map(|(h, p)| (h.x_center() + p.x - cx_ref).abs() + 0.5 * h.length())
                .fold(0.0, f64::max),
            y_half: members
                .iter()
                .map(|(_, p)| (p.y - y_ref).abs())
                .fold(0.0, f64::max),
            t_max: members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
        };
        let mem = members
            .iter()
            .map(|(h, p)| SectionalMember {
                hull: h,
                patch: None,
                scale: 0.0,
                scratch: SectionalContracted::default(),
                dx: h.x_center() + p.x - cx_ref,
                dy: p.y - y_ref,
                closure: opts.transom,
            })
            .collect();
        run_outer(
            &params,
            opts,
            4.0 * cond.fluid.density * g * g / (PI * u * u),
            mem,
        )
        .resistance
    }

    /// A catamaran's resistance from one tabulated demihull amplitude is the
    /// direct two-hull integral's, for a closed hull and a transom hull.
    #[test]
    fn a_catamaran_from_a_tabulated_amplitude_is_exact() {
        let wigley = {
            let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
            hullgeom::iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0)
                .unwrap()
                .situate_sectional(
                    0,
                    0.0,
                    &Default::default(),
                    &Default::default(),
                    &Default::default(),
                )
                .unwrap()
                .unwrap()
                .hull
        };
        let e12 = crate::cad_fixture("e12.igs").map(|t| {
            let o = hullgeom::iges::SectionalOptions {
                waterline_z: -0.95,
                ..Default::default()
            };
            hullgeom::iges::import_sectional(&t, &o)
                .unwrap()
                .hulls
                .remove(0)
                .hull
        });
        let opts = WaveOptions::default();
        for h in std::iter::once(&wigley).chain(e12.as_ref()) {
            for (fnum, span) in [
                (0.3, 2.0 * h.length() / 10.0),
                (0.45, 3.5 * h.length() / 10.0),
            ] {
                let cond = Conditions::seawater(fnum * (9.81 * h.length()).sqrt());
                let pair = [
                    (
                        h,
                        Placement {
                            x: 0.0,
                            y: 0.5 * span,
                        },
                    ),
                    (
                        h,
                        Placement {
                            x: 0.0,
                            y: -0.5 * span,
                        },
                    ),
                ];
                let fast = multihull_wave_resistance(&pair, &cond, &opts)
                    .unwrap()
                    .resistance;
                let slow = direct(&pair, &cond, &opts);
                assert!(
                    (fast - slow).abs() < 1e-5 * slow,
                    "L {} Fn {fnum} span {span}: tabulated {fast} vs direct {slow}",
                    h.length()
                );
            }
        }
    }

    /// The same for a catamaran whose demihulls each carry a part of their
    /// own (a drive's surface-piercing leg, here a small Wigley strut near
    /// the stern): one table per original, each copy reading its own.
    #[test]
    fn a_catamaran_with_parts_from_tabulated_amplitudes_is_exact() {
        let wigley = |l: f64, b: f64, t: f64| {
            let surfaces = hullgeom::iges::wigley_surfaces(l, b, t).unwrap();
            hullgeom::iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0)
                .unwrap()
                .situate_sectional(
                    0,
                    0.0,
                    &Default::default(),
                    &Default::default(),
                    &Default::default(),
                )
                .unwrap()
                .unwrap()
                .hull
        };
        let (hull, strut) = (wigley(10.0, 1.0, 0.625), wigley(0.4, 0.05, 0.5));
        let cond = Conditions::seawater(0.35 * (9.81f64 * 10.0).sqrt());
        let opts = WaveOptions::default();
        let span = 3.0;
        let fleet: Vec<(&SectionalHull, Placement)> = [0.5 * span, -0.5 * span]
            .into_iter()
            .flat_map(|y| {
                [
                    (&hull, Placement { x: 0.0, y }),
                    (&strut, Placement { x: -4.0, y }),
                ]
            })
            .collect();
        let fast = multihull_wave_resistance(&fleet, &cond, &opts)
            .unwrap()
            .resistance;
        let slow = direct(&fleet, &cond, &opts);
        assert!(
            (fast - slow).abs() < 1e-5 * slow,
            "tabulated {fast} vs direct {slow}"
        );
        // The struts are in it: without them the resistance differs.
        let bare = multihull_wave_resistance(&[fleet[0], fleet[2]], &cond, &opts)
            .unwrap()
            .resistance;
        assert!((fast - bare).abs() > 1e-3 * bare, "{fast} vs {bare}");
    }
}

#[cfg(test)]
mod transom_tests {
    use super::*;
    use crate::michell::{InnerIntegral, TransomClosure};

    const E12_WL: f64 = -0.95;

    fn e12_text() -> Option<String> {
        crate::cad_fixture("e12.igs")
    }

    /// A smooth hull with a wet transom: cubic along x, full aft (half-beam
    /// 0.3 m at x = 0) and closed at the bow (x = 8), quadratic in depth and
    /// closed at the keel (T = 0.5).
    fn transom_hull() -> Hull {
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 4.0, 8.0, 8.0, 8.0, 8.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5];
        let control: Vec<f64> = [0.3, 0.45, 0.5, 0.35, 0.0]
            .iter()
            .flat_map(|&c| [c, 0.9 * c, 0.0])
            .collect();
        let surface =
            hullgeom::bspline::BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        Hull::new(surface).unwrap()
    }

    /// The harness's transom is `f_T = 0.3(1 − t)² + 0.54 t(1 − t)` at
    /// `z = 0.5t`: its first moment, both sides, is `2 · 0.25 · (0.3/12 +
    /// 0.27/6) = 0.035 m³`. And the dryness the transom drag fades in with
    /// runs from wet to dry across the closure's band.
    #[test]
    fn the_transom_carries_its_hydrostatic_moment() {
        let hull = transom_hull();
        let t = hull.transom().expect("a wet transom");
        assert!((t.depth_moment - 0.035).abs() < 1e-12, "{}", t.depth_moment);
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        assert_eq!(sec.transom().map(|t| t.depth_moment), Some(t.depth_moment));
        use crate::michell::{transom_dryness, AUTO_DRY_FN_T, AUTO_WET_FN_T};
        assert_eq!(transom_dryness(AUTO_WET_FN_T), 0.0);
        assert_eq!(transom_dryness(AUTO_DRY_FN_T), 1.0);
        assert_eq!(transom_dryness(0.5 * (AUTO_WET_FN_T + AUTO_DRY_FN_T)), 0.5);
    }

    fn closures() -> [TransomClosure; 3] {
        [
            TransomClosure::default(),
            TransomClosure::Fixed { length: 0.3 },
            TransomClosure::None,
        ]
    }

    /// The closure machinery against the exact B-spline kernel's: on the
    /// harness (the spline's own stations and transom), the amplitude and all
    /// six near-field transforms with each closure, at two speeds.
    #[test]
    fn transom_closure_matches_the_lofted_kernel() {
        let hull = transom_hull();
        assert!(
            hull.transom().is_some(),
            "the hull should carry a wet transom"
        );
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let l = hull.length();
        for fnum in [0.2, 0.4] {
            let nu = 1.0 / (fnum * fnum * l);
            for closure in closures() {
                let mut inner = InnerIntegral::new(&hull, nu, closure);
                let mut zc = SectionalContracted::default();
                let (mut worst, mut scale) = (0.0f64, 0.0f64);
                for i in 0..400 {
                    let lambda = 1.0 + i as f64 * 0.25;
                    let a = inner.eval(lambda);
                    let b = sec.amplitude_closed(nu, lambda, closure, &mut zc);
                    worst = worst.max((a - b).abs());
                    scale = scale.max(a.abs());
                }
                let amp = worst / scale;
                let parts = |t: &SquatTransforms| [t.q, t.p, t.q1, t.w, t.p_wl, t.q1_wl];
                let (mut wa, mut sa) = ([0.0f64; 6], [0.0f64; 6]);
                for ik in 0..30 {
                    let kappa = nu * 1e-3 * 10f64.powf(ik as f64 * 0.2);
                    sec.contract(kappa, &mut zc);
                    for ix in 0..30 {
                        let kx = nu * (ix as f64 * 0.7 - 10.0);
                        let a = parts(&inner.eval_transforms(kx, kappa));
                        let b = parts(&sec.transforms_at_closed(&zc, kx, nu, closure));
                        for k in 0..6 {
                            wa[k] = wa[k].max((a[k] - b[k]).abs());
                            sa[k] = sa[k].max(a[k].abs());
                        }
                    }
                }
                let tr = (0..6).map(|k| wa[k] / sa[k]).fold(0.0, f64::max);
                eprintln!("Fn {fnum} {closure:?}: amplitude {amp:.1e}, transforms {tr:.1e}");
                assert!(
                    amp < 1e-9 && tr < 1e-9,
                    "Fn {fnum} {closure:?}: {amp:.1e} {tr:.1e}"
                );
            }
        }
    }

    /// End to end through the outer quadrature and the near-field
    /// integrals: resistance, force and moment with the closure on the
    /// harness equal the lofted hull's.
    #[test]
    fn closed_transom_resistance_and_squat_match_the_lofted_hull() {
        let hull = transom_hull();
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let wave = WaveOptions::default();
        let squat = SquatOptions {
            wave,
            ..SquatOptions::default()
        };
        let x_ref = hull.x_center();
        let rel = |a: f64, b: f64| (a - b).abs() / a.abs();
        for fnum in [0.2, 0.4] {
            let cond = Conditions::seawater(fnum * (9.81 * hull.length()).sqrt());
            let a = crate::michell::wave_resistance_with(&hull, &cond, &wave).unwrap();
            let b = wave_resistance(&sec, &cond, &wave).unwrap();
            let fa = crate::squat::dynamic_force(&hull, &cond, x_ref, &squat).unwrap();
            let fb = dynamic_force(&sec, &cond, x_ref, &squat).unwrap();
            let (er, ef, em) = (
                rel(a.resistance, b.resistance),
                rel(fa.force_up, fb.force_up),
                rel(fa.moment_bow_up, fb.moment_bow_up),
            );
            eprintln!("Fn {fnum}: Rw {er:.1e}, Fz {ef:.1e}, M {em:.1e}");
            assert!(
                er < 1e-8 && ef < 1e-8 && em < 1e-8,
                "Fn {fnum}: {er:.1e} {ef:.1e} {em:.1e}"
            );
        }
    }

    /// Sections cut from CAD find e12's transom: the aft end station, about
    /// 0.37 m wide at the waterline and a centimetre deep.
    #[test]
    fn a_cad_transom_is_detected_from_its_end_section() {
        let so = hullgeom::iges::SectionalOptions {
            waterline_z: E12_WL,
            ..Default::default()
        };
        let Some(text) = e12_text() else {
            return;
        };
        let imp = hullgeom::iges::import_sectional(&text, &so)
            .unwrap()
            .hulls
            .remove(0);
        let tr = imp.hull.transom().expect("e12 has a wet transom");
        eprintln!("transom {tr:?}");
        assert!((tr.x - 0.18).abs() < 5e-3, "x {}", tr.x);
        assert!(
            tr.half_beam > 0.17 && tr.half_beam < 0.2,
            "half-beam {}",
            tr.half_beam
        );
        assert!(tr.depth > 0.004 && tr.depth < 0.012, "depth {}", tr.depth);
        assert_eq!(imp.report.transom.as_ref().map(|t| t.x), Some(tr.x));
    }

    /// The transom's closure is a hollow in the water, not hull: it adds
    /// sources but no pressure-bearing surface. So on e12, whose transom is
    /// barely immersed (6.5 mm equivalent depth, but 0.37 m wide at the
    /// waterline), the near-field force and moment hardly depend on the
    /// hollow's length — they once moved by 35% and changed sign when the
    /// appendage's waterplane was counted as hull.
    #[test]
    fn a_shallow_transom_closure_barely_moves_the_squat_force() {
        let Some(text) = e12_text() else {
            return;
        };
        let opts = hullgeom::iges::SectionalOptions {
            waterline_z: -0.95,
            ..Default::default()
        };
        let fleet = hullgeom::iges::import_sectional(&text, &opts).unwrap();
        let h = &fleet.hulls[0].hull;
        assert!(h.transom().is_some());
        let cond = Conditions::seawater(0.3 * (9.81f64 * h.length()).sqrt());
        let force = |closure| {
            let so = crate::squat::SquatOptions {
                wave: WaveOptions {
                    transom: closure,
                    ..Default::default()
                },
                ..Default::default()
            };
            dynamic_force(h, &cond, h.lcb_x(), &so).unwrap()
        };
        let open = force(TransomClosure::None);
        for closure in [
            TransomClosure::default(),
            TransomClosure::Fixed { length: 1.0 },
        ] {
            let d = force(closure);
            let (df, dm) = (
                (d.force_up - open.force_up).abs() / open.force_up.abs(),
                (d.moment_bow_up - open.moment_bow_up).abs() / open.moment_bow_up.abs(),
            );
            assert!(
                df < 0.05 && dm < 0.05,
                "{closure:?}: force {df:.3}, moment {dm:.3} off"
            );
        }
    }
}
