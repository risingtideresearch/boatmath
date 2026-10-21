//! Hull–propeller interaction in potential flow: the **nominal wake** a
//! propeller disc sits in, and the **thrust deduction** its suction costs
//! the hull, both from the thin-ship singularities [`crate::nearfield`]
//! already places.
//!
//! ## Wake
//!
//! The hull's near field (its centreplane source sheets, Rankine pair and
//! wave part) gives the axial perturbation velocity `φ_x` anywhere in the
//! water. With the stream toward −x, `φ_x > 0` slows it, so the wake
//! fraction at a point is `φ_x/U`, and the disc's is its area average over
//! the annulus from the hub to the tip. It is split into the local
//! (Rankine) part, the potential wake proper, and the wave part, the
//! orbital velocity of the hull's own waves. The **frictional** wake — the
//! boundary layer, usually the largest part behind a full single-screw
//! stern — is not in a potential model and is not here.
//!
//! ## Thrust deduction
//!
//! Dickmann's model: upstream of it, an actuator disc's flow is that of a
//! uniform sink sheet of density `2u_a` over the disc, `u_a` its induced
//! velocity from momentum theory, `T = 2ρA u_a (V_A + u_a)`. By Lagally's
//! theorem a source sheet `σ` in an external flow `u` feels
//! `F_x = −ρ∬σ u dA`, so the hull's resistance rises by
//!
//! ```text
//! ΔR = ρ ∬ σ u_p dA,     t = ΔR / T,
//! ```
//!
//! `u_p` the sinks' axial velocity over the hull's sheets. The afterbody's
//! sinks (`σ < 0`, the hull narrowing aft) sit in the flow drawn aft toward
//! the propeller (`u_p < 0`): `ΔR > 0`. The sinks' field here is their
//! Rankine pair (sink minus its negative image, φ = 0 on the free surface —
//! the near field's local part); the waves the propeller itself makes are
//! not included.
//!
//! Frame as the near field: x forward (the stream toward −x), y to port, z
//! **down** from the still water, metres, the fleet's placement frame.

// `!(x > 0.0)` rejects NaN too.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use crate::nearfield::{
    build_nodes, build_sheet, check, rankine_phi_x, Mode, NearFieldOptions, Nodes, Sheet,
};
use hullgeom::quadrature::gauss_legendre;
use hullgeom::sectional::SectionalHull;
use hullgeom::{Conditions, Error, Placement, Result};
use std::f64::consts::PI;

/// A propeller disc, its axis along x.
#[derive(Debug, Clone, Copy)]
pub struct Disc {
    pub x: f64,
    pub y: f64,
    /// Of its centre below the still water.
    pub depth: f64,
    pub radius: f64,
    /// Hub radius over tip radius: the annulus the averages run over.
    pub hub: f64,
}

/// A disc's nominal wake.
#[derive(Debug, Clone)]
pub struct Wake {
    /// The area-mean wake fraction, including potential (local + wave) and
    /// frictional boundary-layer parts.
    pub w: f64,
    /// The local (Rankine) part: the potential wake.
    pub w_local: f64,
    /// The wave part: the hull's own waves' orbital velocity.
    pub w_wave: f64,
    /// The frictional boundary-layer wake.
    pub w_friction: f64,
    /// `(r/R, w)` of each ring, hub to tip.
    pub radial: Vec<(f64, f64)>,
}

/// The thrust deduction of a set of discs, each giving `thrust`.
#[derive(Debug, Clone)]
pub struct Deduction {
    pub t: f64,
    /// The hull's added resistance [N], all hulls.
    pub delta_r: f64,
    /// Each hull's share of it [N].
    pub per_hull: Vec<f64>,
    /// Each disc's induced axial velocity at the disc, from momentum theory.
    pub u_a: f64,
}

/// The ring radii and the angles of a disc's quadrature, and each point's
/// share of the annulus' area: `(y, z, dA, ring)`.
fn disc_points(d: &Disc) -> Vec<(f64, f64, f64, usize)> {
    let (gx, gw) = gauss_legendre(6);
    let n_a = 24;
    let (r0, r1) = (d.hub * d.radius, d.radius);
    let mut out = Vec::with_capacity(gx.len() * n_a);
    for (ring, (&t, &w)) in gx.iter().zip(&gw).enumerate() {
        let r = 0.5 * (r1 - r0) * t + 0.5 * (r1 + r0);
        let dr = 0.5 * (r1 - r0) * w;
        for a in 0..n_a {
            let th = 2.0 * PI * (a as f64 + 0.5) / n_a as f64;
            let da = r * dr * 2.0 * PI / n_a as f64;
            out.push((d.y + r * th.cos(), d.depth + r * th.sin(), da, ring));
        }
    }
    out
}

/// A fleet's singularities at one speed, ready to evaluate the flow at a
/// propeller.
pub struct Interaction {
    cond: Conditions,
    sheets: Vec<Sheet>,
    nodes: Nodes,
    cx: Vec<f64>,
    ys: Vec<f64>,
}

impl Interaction {
    pub fn new(
        members: &[(&SectionalHull, Placement)],
        cond: &Conditions,
        opts: &NearFieldOptions,
    ) -> Result<Interaction> {
        check(members, cond)?;
        let u = cond.speed;
        let nu = cond.gravity / (u * u);
        Ok(Interaction {
            sheets: members
                .iter()
                .map(|(h, p)| build_sheet(h, *p, u, nu, opts))
                .collect(),
            nodes: build_nodes(
                members,
                cond,
                opts.closure,
                Mode::Hull {
                    level: opts.level.max(1),
                },
            ),
            cx: members.iter().map(|(h, p)| h.x_center() + p.x).collect(),
            ys: members.iter().map(|(_, p)| p.y).collect(),
            cond: *cond,
        })
    }

    /// The hulls' `φ_x` at `(x, y, z)`, `z` down: `(local, wave)`.
    pub fn phi_x(&self, x: f64, y: f64, z: f64) -> (f64, f64) {
        (
            rankine_phi_x(&self.sheets, x, y, z),
            self.nodes.phi_x(&self.cx, &self.ys, x, y, z),
        )
    }

    fn check_disc(&self, d: &Disc) -> Result<()> {
        if !(d.radius > 0.0) || !(0.0..1.0).contains(&d.hub) {
            return Err(Error::InvalidInput(format!(
                "propeller disc: radius {} and hub ratio {}",
                d.radius, d.hub
            )));
        }
        if d.depth - d.radius <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "the propeller breaks the surface: its centre is {:.3} m down, its radius {:.3} m",
                d.depth, d.radius
            )));
        }
        if disc_points(d)
            .iter()
            .any(|&(y, z, _, _)| self.sheets.iter().any(|s| inside(s, d.x, y, z)))
        {
            return Err(Error::InvalidInput(format!(
                "the propeller disc at x {:.3} m, y {:.3} m, {:.3} m down cuts a hull",
                d.x, d.y, d.depth
            )));
        }
        Ok(())
    }

    /// Boundary-layer frictional wake fraction at `(x, y, z)`: 1/7th-power-law
    /// turbulent boundary layer delta(x) = 0.37 * x * Re_x^(-1/5) along the
    /// wetted surface upstream of the disc, decaying downstream of the stern.
    pub fn frictional_wake_at(&self, x: f64, y: f64, z: f64) -> f64 {
        let u = self.cond.speed;
        if u <= 0.0 {
            return 0.0;
        }
        let nu = self.cond.fluid.kinematic_viscosity.max(1e-9);
        let mut max_wf: f64 = 0.0;

        for s in &self.sheets {
            let xs = &s.xs[s.real_from..];
            if xs.len() < 2 {
                continue;
            }
            let x_stern = xs[0];
            let x_bow = xs[xs.len() - 1];
            // If the point is forward of the bow, there is no upstream boundary layer.
            if x >= x_bow {
                continue;
            }
            // Run length along the hull upstream of this station (x points forward).
            let x_station = x.clamp(x_stern, x_bow);
            let l_run = x_bow - x_station;
            if l_run <= 0.0 {
                continue;
            }

            let re_x = (u * l_run / nu).max(1e3);
            let delta0 = 0.37 * l_run * re_x.powf(-0.2);

            // Distance downstream of the stern (if aft of it).
            let d_x = (x_stern - x).max(0.0);
            // Turbulent wake spreading and centerline velocity defect decay downstream.
            let delta = delta0 * (1.0 + 0.5 * d_x / delta0.max(1e-4)).sqrt();
            let decay = 1.0 / (1.0 + 0.8 * d_x / delta0.max(1e-4)).sqrt();

            // Find half-beam b at depth z and station x_station on this sheet.
            let b = half_beam_at(s, x_station, z);

            let dy = (y - s.y).abs() - b;
            let dist = if z <= s.depth[s.depth.len() - 1] {
                dy.max(0.0)
            } else {
                let dz = z - s.depth[s.depth.len() - 1];
                (dy.max(0.0).powi(2) + dz.powi(2)).sqrt()
            };

            if dist < delta {
                let local_wf = (1.0 - (dist / delta).powf(1.0 / 7.0)) * decay;
                max_wf = max_wf.max(local_wf);
            }
        }
        max_wf
    }

    /// The nominal wake over a disc.
    pub fn wake(&self, d: &Disc) -> Result<Wake> {
        self.check_disc(d)?;
        let u = self.cond.speed;
        let pts = disc_points(d);
        let mut rings = vec![(0.0f64, 0.0f64); 6];
        let (mut loc, mut wav, mut fric, mut area) = (0.0, 0.0, 0.0, 0.0);
        let vals: Vec<(f64, f64, f64)> = hullgeom::parallel::map_indexed(
            pts.len(),
            || (),
            |_, i| {
                let (y, z, _, _) = pts[i];
                let (l, w) = self.phi_x(d.x, y, z);
                let wf = self.frictional_wake_at(d.x, y, z);
                (l, w, wf)
            },
        );
        for (&(_, _, da, ring), &(l, w, wf)) in pts.iter().zip(&vals) {
            loc += l * da;
            wav += w * da;
            fric += wf * da;
            area += da;
            rings[ring].0 += (l + w + wf * u) * da;
            rings[ring].1 += da;
        }
        let (gx, _) = gauss_legendre(6);
        let radial = gx
            .iter()
            .zip(&rings)
            .map(|(&t, &(s, a))| {
                let r = 0.5 * (1.0 - d.hub) * t + 0.5 * (1.0 + d.hub);
                (r, s / a / u)
            })
            .collect();
        let w_local = loc / area / u;
        let w_wave = wav / area / u;
        let w_friction = fric / area;
        Ok(Wake {
            w: w_local + w_wave + w_friction,
            w_local,
            w_wave,
            w_friction,
            radial,
        })
    }

    /// The thrust deduction of `discs`, each giving `thrust` [N] in an
    /// inflow `v_a` [m/s].
    pub fn thrust_deduction(&self, discs: &[Disc], thrust: f64, v_a: f64) -> Result<Deduction> {
        for d in discs {
            self.check_disc(d)?;
        }
        if !(thrust > 0.0) || !(v_a > 0.0) || discs.is_empty() {
            return Err(Error::InvalidInput(
                "thrust deduction: needs discs, a thrust and an inflow".into(),
            ));
        }
        let rho = self.cond.fluid.density;
        let mut sinks: Vec<(f64, f64, f64, f64)> = Vec::new();
        let mut u_a = 0.0;
        for d in discs {
            let pts = disc_points(d);
            let a: f64 = pts.iter().map(|p| p.2).sum();
            u_a = 0.5 * (-v_a + (v_a * v_a + 2.0 * thrust / (rho * a)).sqrt());
            sinks.extend(
                pts.iter()
                    .map(|&(y, z, da, _)| (d.x, y, z, -2.0 * u_a * da)),
            );
        }
        let per_hull: Vec<f64> = self
            .sheets
            .iter()
            .map(|s| rho * sheet_integral(s, self.cond.speed, |x, y, z| sinks_u(&sinks, x, y, z)))
            .collect();
        let delta_r: f64 = per_hull.iter().sum();
        Ok(Deduction {
            t: delta_r / (thrust * discs.len() as f64),
            delta_r,
            per_hull,
            u_a,
        })
    }

    /// The same added resistance the other way round, by Newton's third
    /// law: the force of the hulls' local field on the sinks, `−ρ Σ m u_h`.
    /// The two agree when there's no closure appendage (whose sheet is not
    /// hull); a check on [`Interaction::thrust_deduction`].
    pub fn added_resistance_reciprocal(&self, discs: &[Disc], thrust: f64, v_a: f64) -> f64 {
        let rho = self.cond.fluid.density;
        let mut total = 0.0;
        for d in discs {
            let pts = disc_points(d);
            let a: f64 = pts.iter().map(|p| p.2).sum();
            let u_a = 0.5 * (-v_a + (v_a * v_a + 2.0 * thrust / (rho * a)).sqrt());
            for &(y, z, da, _) in &pts {
                let m = -2.0 * u_a * da;
                total += -rho * m * rankine_phi_x(&self.sheets, d.x, y, z);
            }
        }
        total
    }
}

fn half_beam_at(s: &Sheet, x: f64, z: f64) -> f64 {
    let xs = &s.xs[s.real_from..];
    if xs.len() < 2 {
        return 0.0;
    }
    let t = s.depth[s.depth.len() - 1];
    let nz1 = s.nz1();
    let i = (xs.partition_point(|&v| v <= x).max(1) - 1).min(xs.len() - 2) + s.real_from;
    let dz = t / (nz1 - 1).max(1) as f64;
    let j = ((z.clamp(0.0, t) / dz.max(1e-12)).floor() as usize).min(nz1 - 2);
    let (fx, fz) = (
        ((x - s.xs[i]) / (s.xs[i + 1] - s.xs[i]).max(1e-12)).clamp(0.0, 1.0),
        ((z.clamp(0.0, t) - s.depth[j]) / dz.max(1e-12)).clamp(0.0, 1.0),
    );
    let f = |i: usize, j: usize| s.f[i * nz1 + j];
    (1.0 - fx) * ((1.0 - fz) * f(i, j) + fz * f(i, j + 1))
        + fx * ((1.0 - fz) * f(i + 1, j) + fz * f(i + 1, j + 1))
}

/// Whether `(x, y, z)` lies inside a sheet's real hull: within its
/// half-beam there, interpolated between stations and depth rows.
fn inside(s: &Sheet, x: f64, y: f64, z: f64) -> bool {
    let xs = &s.xs[s.real_from..];
    let (Some(&lo), Some(&hi)) = (xs.first(), xs.last()) else {
        return false;
    };
    let t = s.depth[s.depth.len() - 1];
    if x < lo || x > hi || !(0.0..=t).contains(&z) {
        return false;
    }
    let half = half_beam_at(s, x, z);
    (y - s.y).abs() < half
}

/// `∬ σ u dA` over the real hull of a sheet, `σ = −2U ∂f/∂x` per cell,
/// `u` at each cell's centre.
fn sheet_integral(s: &Sheet, speed: f64, u: impl Fn(f64, f64, f64) -> f64 + Sync) -> f64 {
    let nz1 = s.nz1();
    let cols: Vec<usize> = (s.real_from..s.xs.len().saturating_sub(1)).collect();
    let per: Vec<f64> = hullgeom::parallel::map_indexed(
        cols.len(),
        || (),
        |_, ci| {
            let i = cols[ci];
            let dx = s.xs[i + 1] - s.xs[i];
            if dx <= 0.0 {
                return 0.0;
            }
            let x = 0.5 * (s.xs[i] + s.xs[i + 1]);
            let mut acc = 0.0;
            for j in 0..nz1 - 1 {
                let fx = (s.f[(i + 1) * nz1 + j] + s.f[(i + 1) * nz1 + j + 1]
                    - s.f[i * nz1 + j]
                    - s.f[i * nz1 + j + 1])
                    / (2.0 * dx);
                if fx == 0.0 {
                    continue;
                }
                let sigma = -2.0 * speed * fx;
                let dz = s.depth[j + 1] - s.depth[j];
                let z = 0.5 * (s.depth[j] + s.depth[j + 1]);
                acc += sigma * u(x, s.y, z) * dx * dz;
            }
            acc
        },
    );
    per.iter().sum()
}

/// The axial velocity at `(x, y, z)` of point sinks `(ξ, η, ζ, m)`, `m` the
/// volume outflow (negative for a sink), each with its negative image above
/// the free surface.
fn sinks_u(sinks: &[(f64, f64, f64, f64)], x: f64, y: f64, z: f64) -> f64 {
    let mut u = 0.0;
    for &(xi, eta, zeta, m) in sinks {
        let (dx, dy) = (x - xi, y - eta);
        let r2 = dx * dx + dy * dy + (z - zeta) * (z - zeta);
        let ri2 = dx * dx + dy * dy + (z + zeta) * (z + zeta);
        if r2 < 1e-18 {
            continue;
        }
        u += m * dx * (r2.powf(-1.5) - ri2.powf(-1.5));
    }
    u / (4.0 * PI)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::michell::TransomClosure;

    fn wigley() -> SectionalHull {
        let w = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let text = hullgeom::iges::write(
            &hullgeom::iges::halfbreadth_surfaces(w.surface(), 0.0, 0.0),
            "wigley",
        )
        .unwrap();
        let f = hullgeom::iges::import_sectional(&text, &Default::default()).unwrap();
        f.hulls[0].hull.clone()
    }

    fn opts() -> NearFieldOptions {
        NearFieldOptions {
            closure: TransomClosure::None,
            ..Default::default()
        }
    }

    /// On the hull the field is the near field's: `2φ_x/U` is the
    /// pressure coefficient `hull_pressure` reports, node for node.
    #[test]
    fn the_field_on_the_hull_is_the_pressure() {
        let h = wigley();
        let cond = Conditions::seawater(0.35 * (9.81f64 * 10.0).sqrt());
        let m = [(&h, Placement::default())];
        let it = Interaction::new(&m, &cond, &opts()).unwrap();
        let hp = &crate::nearfield::hull_pressure(&m, &cond, &opts()).unwrap()[0];
        let nz1 = hp.depth.len();
        let mut worst: f64 = 0.0;
        for ix in (0..hp.x.len()).step_by(7) {
            for j in (0..nz1).step_by(5) {
                let (l, w) = it.phi_x(hp.x[ix], hp.y, hp.depth[j]);
                let cp = 2.0 * (l + w) / cond.speed;
                worst = worst.max((cp - hp.cp[ix * nz1 + j]).abs());
            }
        }
        assert!(worst < 1e-9, "{worst}");
        // At the surface the local part vanishes.
        assert!(it.phi_x(0.0, 2.0, 0.0).0.abs() < 1e-12);
    }

    /// Behind the stern the flow is slowed (a positive wake), symmetric
    /// either side, and it fades with distance.
    #[test]
    fn the_wake_behind_a_stern() {
        let h = wigley();
        let cond = Conditions::seawater(0.3 * (9.81f64 * 10.0).sqrt());
        let m = [(&h, Placement::default())];
        let it = Interaction::new(&m, &cond, &opts()).unwrap();
        let (xa, _) = h.x_range();
        let disc = |dx: f64, y: f64| Disc {
            x: xa - dx,
            y,
            depth: 0.35,
            radius: 0.2,
            hub: 0.2,
        };
        let near = it.wake(&disc(0.1, 0.0)).unwrap();
        let far = it.wake(&disc(3.0, 0.0)).unwrap();
        let (port, stbd) = (
            it.wake(&disc(0.1, 0.3)).unwrap(),
            it.wake(&disc(0.1, -0.3)).unwrap(),
        );
        eprintln!(
            "wake 0.1 m astern: w {:.4} (local {:.4}, wave {:.4}, friction {:.4}); 3 m astern {:.4}",
            near.w, near.w_local, near.w_wave, near.w_friction, far.w
        );
        assert!(near.w_local > 0.0, "{near:?}");
        assert!(near.w_friction > 0.0, "{near:?}");
        assert!(
            far.w_friction < near.w_friction,
            "far {far:?} vs near {near:?}"
        );
        assert!(far.w_local.abs() < near.w_local);
        assert!((port.w - stbd.w).abs() < 1e-9 * near.w.abs().max(1e-12));
        assert!((port.w_friction - stbd.w_friction).abs() < 1e-9 * near.w.abs().max(1e-12));
    }

    /// The added resistance by the sinks' field over the hull and by the
    /// hull's field at the sinks agree (Newton's third law, a check on
    /// both); a propeller astern adds resistance; doubling the thrust
    /// roughly doubles it, and t falls as the disc moves aft.
    #[test]
    fn a_disc_in_the_hull_is_refused() {
        let h = wigley();
        let cond = Conditions::seawater(3.0);
        let m = [(&h, Placement::default())];
        let it = Interaction::new(&m, &cond, &opts()).unwrap();
        let (xa, xb) = h.x_range();
        let disc = |x: f64| Disc {
            x,
            y: 0.0,
            depth: 0.3,
            radius: 0.1,
            hub: 0.2,
        };
        let e = it.wake(&disc(0.5 * (xa + xb))).unwrap_err().to_string();
        assert!(e.contains("cuts a hull"), "{e}");
        assert!(it.wake(&disc(xa - 0.2)).is_ok());
    }

    #[test]
    fn thrust_deduction_both_ways() {
        let h = wigley();
        let cond = Conditions::seawater(0.3 * (9.81f64 * 10.0).sqrt());
        let m = [(&h, Placement::default())];
        let it = Interaction::new(&m, &cond, &opts()).unwrap();
        let (xa, _) = h.x_range();
        let disc = |dx: f64| Disc {
            x: xa - dx,
            y: 0.0,
            depth: 0.35,
            radius: 0.2,
            hub: 0.2,
        };
        let (t, v_a) = (200.0, 2.9);
        let d = it.thrust_deduction(&[disc(0.15)], t, v_a).unwrap();
        let r = it.added_resistance_reciprocal(&[disc(0.15)], t, v_a);
        eprintln!(
            "ΔR {:.3} N by the sheet, {:.3} N at the sinks; t {:.4}, u_a {:.3} m/s",
            d.delta_r, r, d.t, d.u_a
        );
        assert!(d.delta_r > 0.0 && d.t > 0.0, "{d:?}");
        assert!(
            (d.delta_r - r).abs() < 0.03 * r.abs(),
            "{} vs {r}",
            d.delta_r
        );
        let aft = it.thrust_deduction(&[disc(1.0)], t, v_a).unwrap();
        assert!(aft.t < d.t, "{} then {}", d.t, aft.t);
        let twice = it.thrust_deduction(&[disc(0.15)], 2.0 * t, v_a).unwrap();
        assert!(twice.delta_r > 1.5 * d.delta_r && twice.delta_r < 2.0 * d.delta_r);
    }
}
