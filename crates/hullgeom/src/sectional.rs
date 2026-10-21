//! **Sectional** hulls: a hull described by its stations, each a section
//! curve with its own depth quadrature, the whole interpolated along the
//! length by a B-spline.
//!
//! Every hull transform the hydrodynamic crates need is of the form
//! `∬ f(x, z) e^{−κz} e^{ik_x x} dx dz` (or its `∂f/∂x` and
//! `(x − x_c) ∂f/∂x` variants). Write it as
//!
//! ```text
//! ∫ Z(x; κ) e^{ik_x x} dx,      Z(x; κ) = ∫ f(x, z) e^{−κz} dz,
//! ```
//!
//! and do the depth integral **per station**, along the station's own section
//! curve. The keel is then just the end of that curve: the crease a lofted
//! `f(x, z)` has along a rockered keel line, and the square-root closure of a
//! round bilge onto it, are both inside a 1-D integral that absorbs them, and
//! `Z(x; κ)` is smooth in `x` wherever the hull's sections vary smoothly.
//!
//! Everything that depends on `κ` is one pass over the stations' depth nodes
//! (`Σ_j w_j e^{−κ z_j}` per station) followed by one banded solve for the
//! B-spline interpolant of `Z(·; κ)` along `x` ([`SectionalHull::contract`]);
//! each `k_x` at that `κ` then costs one closed-form oscillatory-moment call
//! per x-span. The interpolant is exposed span by span in each span's local
//! power basis ([`SectionalHull::spans`], [`SectionalContracted::g_f`]), so
//! a hydrodynamic crate can build its own transforms on it; the plain
//! Fourier transforms are here ([`SectionalHull::x_transform`]).
//!
//! Hulls are cut from CAD through [`crate::source::HullSource`]; the
//! hydrostatics (volume, centre of buoyancy, waterplane, wetted surface) are
//! exact integrals of the same interpolants.

use crate::bspline::{ders_basis, find_span};
use crate::error::{Error, Result};
use crate::moments::{osc_moments, C64};
use crate::quadrature::gauss_legendre;

/// Below this, `e^{−κz}` cannot matter against the shallowest node's weight
/// (the same floor the exact B-spline kernel drops z-spans at).
const DECAY_EXPONENT_FLOOR: f64 = 46.0; // e^{-46} ≈ 1e-20

/// One station's depth integral `Z(κ) = ∫ f(x_i, z) e^{−κz} dz`, as a
/// quadrature that costs one exponential per node at each `κ`.
#[derive(Debug, Clone)]
pub struct SectionNodes {
    nodes: Nodes,
    /// Half-beam at the waterline, `f(x_i, 0)`.
    waterline: f64,
    /// Depth of the section's lowest point.
    depth: f64,
    /// `(half-beam, depth)` at each quadrature node, in node order: where
    /// the depth integral evaluates the section.
    outline: Vec<(f64, f64)>,
    /// The section curve the quadrature integrates, sampled densely and
    /// evenly from the waterline (or the section's top) round to the keel —
    /// for drawing it and for measuring its girth, where joining the
    /// quadrature nodes (graded toward the waterline, sparse at the keel)
    /// would cut corners.
    curve: Vec<(f64, f64)>,
}

/// Points per section in [`SectionNodes::curve`].
const CURVE_POINTS: usize = 129;

#[derive(Debug, Clone)]
enum Nodes {
    /// `Σ_j w_j e^{−κ z_j}`, nodes sorted by depth: a half-beam sampled as
    /// a function of depth.
    Depth { z: Vec<f64>, w: Vec<f64> },
    /// `e^{−κ z₀} Σ_k w_k F(κ d sin θ_k, R_k)` with `F(a, R) = ∫_0^R e^{−ar} r dr`
    /// (`sin` stores `d sin θ_k`, `w` carries the scales `b·d`):
    /// the section as a region swept by rays from its top centreplane point
    /// `(0, z₀)`, the ray at angle `θ` below the horizontal reaching the shell
    /// at distance `R(θ)`. The depth integral is then an area integral in
    /// polar coordinates whose radial part is closed form, and `R(θ)` is
    /// smooth right into the keel — where a half-beam as a function of depth
    /// closes like a square root on any round bilge.
    Polar {
        z0: f64,
        sin: Vec<f64>,
        r: Vec<f64>,
        w: Vec<f64>,
    },
}

/// `∫_0^R e^{−ar} r dr = R²·(1 − e^{−x}(1 + x))/x²`, `x = aR`, by series where
/// the closed form cancels.
fn polar_radial(a: f64, r: f64) -> f64 {
    let x = a * r;
    if x < 0.1 {
        // Σ_n (−1)ⁿ (n+1) xⁿ/(n+2)!
        let (mut term, mut sum, mut fact) = (1.0f64, 0.5f64, 2.0f64);
        for n in 1..12 {
            term *= -x;
            fact *= (n + 2) as f64;
            sum += (n + 1) as f64 * term / fact;
        }
        r * r * sum
    } else {
        r * r * (1.0 - (-x).exp() * (1.0 + x)) / (x * x)
    }
}

impl SectionNodes {
    /// Nodes for a section given as a half-beam function of depth on
    /// `[0, depth]`: Gauss–Legendre panels between the `breaks` (kinks,
    /// chines, knot lines — anywhere the half-beam is not smooth), further
    /// split geometrically toward the waterline so that `e^{−κz}` stays
    /// resolved when a large `κ` confines it to a sliver under the surface.
    pub fn from_depth_function(
        depth: f64,
        breaks: &[f64],
        half_beam: impl Fn(f64) -> f64,
        opts: &DepthQuadrature,
    ) -> SectionNodes {
        let cuts = graded_cuts(depth, breaks, opts);
        let (gx, gw) = gauss_legendre(opts.points);
        let (mut z, mut w) = (Vec::new(), Vec::new());
        let mut outline = vec![(half_beam(0.0), 0.0)];
        for c in cuts.windows(2) {
            let (a, b) = (c[0], c[1]);
            let half = 0.5 * (b - a);
            for (&t, &wt) in gx.iter().zip(&gw) {
                let zj = a + half * (t + 1.0);
                let h = half_beam(zj);
                z.push(zj);
                w.push(half * wt * h);
                outline.push((h, zj));
            }
        }
        outline.push((half_beam(depth), depth));
        // Evenly in depth, with every break (a chine, a knot line) on it.
        let mut zs: Vec<f64> = (0..CURVE_POINTS)
            .map(|k| depth * k as f64 / (CURVE_POINTS - 1) as f64)
            .chain(breaks.iter().copied().filter(|&b| b > 0.0 && b < depth))
            .collect();
        zs.sort_by(f64::total_cmp);
        let curve = zs.iter().map(|&z| (half_beam(z), z)).collect();
        SectionNodes {
            nodes: Nodes::Depth { z, w },
            waterline: half_beam(0.0),
            depth,
            outline,
            curve,
        }
    }

    /// Nodes for a section swept by rays from `(0, z0)` on the centreplane,
    /// in the section's own proportions: with half-beam scaled by `beam` and
    /// depth by `depth`, `radius(θ)` is the scaled distance to the shell
    /// along the scaled ray `θ` below the horizontal, `θ ∈ [0, π/2]` (`π/2`
    /// runs down the centreplane to the keel). Sweeping the scaled plane
    /// keeps rays spread over a thin section (a fine entry is millimetres
    /// wide and the full draft deep), where rays even in physical angle
    /// would crowd its whole outline into a sliver of angle beside the
    /// keel. Requires the section to be star-shaped about `(0, z0)`. Graded
    /// toward `θ = 0`, where a large `κ` concentrates the integrand against
    /// the surface.
    pub fn from_polar(
        z0: f64,
        beam: f64,
        depth: f64,
        radius: impl Fn(f64) -> f64,
        opts: &DepthQuadrature,
    ) -> SectionNodes {
        let half_pi = std::f64::consts::FRAC_PI_2;
        let cuts = graded_cuts(half_pi, &[], opts);
        let (gx, gw) = gauss_legendre(opts.points);
        let (mut sin, mut r, mut w) = (Vec::new(), Vec::new(), Vec::new());
        let at = |th: f64, rr: f64| (beam * rr * th.cos(), z0 + depth * rr * th.sin());
        let mut outline = vec![at(0.0, radius(0.0))];
        for c in cuts.windows(2) {
            let (a, b) = (c[0], c[1]);
            let half = 0.5 * (b - a);
            for (&t, &wt) in gx.iter().zip(&gw) {
                let th = a + half * (t + 1.0);
                let rr = radius(th);
                // Physical depth along the scaled ray is depth·r·sin θ.
                sin.push(depth * th.sin());
                r.push(rr);
                w.push(half * wt * beam * depth);
                outline.push(at(th, rr));
            }
        }
        outline.push(at(half_pi, radius(half_pi)));
        // Evenly in the scaled ray angle: evenly round a section of any
        // proportions.
        let curve = (0..CURVE_POINTS)
            .map(|k| {
                let th = half_pi * k as f64 / (CURVE_POINTS - 1) as f64;
                at(th, radius(th))
            })
            .collect();
        SectionNodes {
            nodes: Nodes::Polar { z0, sin, r, w },
            waterline: if z0 == 0.0 { beam * radius(0.0) } else { 0.0 },
            depth: z0 + depth * radius(half_pi),
            outline,
            curve,
        }
    }

    /// A station with no wetted section (beyond the hull's ends).
    pub fn empty() -> SectionNodes {
        SectionNodes {
            nodes: Nodes::Depth {
                z: Vec::new(),
                w: Vec::new(),
            },
            waterline: 0.0,
            depth: 0.0,
            outline: Vec::new(),
            curve: Vec::new(),
        }
    }

    /// The section at its quadrature nodes, `(half-beam, depth)` from the
    /// waterline (or the section's top) round to the keel: exactly what the
    /// depth integral integrates, for display.
    pub fn outline(&self) -> &[(f64, f64)] {
        &self.outline
    }

    /// The section curve the quadrature integrates, sampled densely and
    /// evenly, `(half-beam, depth)` from the waterline (or top) to the keel.
    pub fn curve(&self) -> &[(f64, f64)] {
        &self.curve
    }

    /// `Z(κ)`.
    fn integrate(&self, kappa: f64) -> f64 {
        match &self.nodes {
            Nodes::Depth { z, w } => {
                let mut s = 0.0;
                for (&z, &w) in z.iter().zip(w) {
                    let e = kappa * z;
                    if e > DECAY_EXPONENT_FLOOR {
                        break;
                    }
                    s += w * (-e).exp();
                }
                s
            }
            Nodes::Polar { z0, sin, r, w } => {
                let e0 = kappa * z0;
                if e0 > DECAY_EXPONENT_FLOOR {
                    return 0.0;
                }
                let s: f64 = sin
                    .iter()
                    .zip(r)
                    .zip(w)
                    .map(|((&sn, &r), &w)| w * polar_radial(kappa * sn, r))
                    .sum();
                s * (-e0).exp()
            }
        }
    }

    /// `∫ z·f(z) dz = −Z′(0)`: the section's first moment of area about the
    /// waterline, one side.
    fn first_moment(&self) -> f64 {
        match &self.nodes {
            Nodes::Depth { z, w } => z.iter().zip(w).map(|(z, w)| w * z).sum(),
            // Depth along a ray is z₀ + r·sin, so the moment's radial part is
            // ∫ (z₀ + r·sin) r dr = z₀R²/2 + sin·R³/3.
            Nodes::Polar { z0, sin, r, w } => sin
                .iter()
                .zip(r)
                .zip(w)
                .map(|((&sn, &r), &w)| w * (z0 * r * r / 2.0 + sn * r * r * r / 3.0))
                .sum(),
        }
    }

    /// Number of quadrature nodes.
    pub fn len(&self) -> usize {
        match &self.nodes {
            Nodes::Depth { z, .. } => z.len(),
            Nodes::Polar { r, .. } => r.len(),
        }
    }

    /// True when the section has no nodes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Panel edges on `[0, extent]`: the `breaks`, plus a geometric sequence
/// toward 0 with ratio `opts.grading` down to `opts.finest · extent`.
fn graded_cuts(extent: f64, breaks: &[f64], opts: &DepthQuadrature) -> Vec<f64> {
    let mut cuts = vec![0.0, extent];
    let mut t = extent;
    while t > opts.finest * extent {
        t *= opts.grading;
        cuts.push(t);
    }
    cuts.extend(breaks.iter().copied().filter(|&b| b > 0.0 && b < extent));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-14 * extent);
    cuts
}

/// How finely each station's depth integral is resolved.
#[derive(Debug, Clone, Copy)]
pub struct DepthQuadrature {
    /// Gauss–Legendre points per panel.
    pub points: usize,
    /// Ratio between successive panel edges toward the waterline.
    pub grading: f64,
    /// Shallowest panel edge, as a fraction of the section depth.
    pub finest: f64,
}

impl Default for DepthQuadrature {
    fn default() -> Self {
        // 8 points on panels graded by 0.35 toward the waterline hold the
        // transforms to ~1e-11 of the exact B-spline kernel's exact values on the
        // Wigley (one z-span: the deepest panel is ~0.4 T long) and ~1e-12
        // on a CAD import, for κ up to ~1e4·ν. Six points leave ~1e-8.
        DepthQuadrature {
            points: 8,
            grading: 0.35,
            finest: 1e-7,
        }
    }
}

/// One non-empty knot span in one direction.
#[derive(Debug, Clone, Copy)]
pub struct Span {
    /// Physical coordinate of the span start.
    pub start: f64,
    /// Span length (> 0).
    pub len: f64,
}

/// A transom whose area is under this fraction of the hull's maximum section
/// area is not reported: it is hydrodynamically negligible, and at a closing
/// stern indistinguishable from the geometry's own noise.
pub const TRANSOM_AREA_REL: f64 = 1e-3;

/// The aft-end section of a hull whose half-breadth does not close there — a
/// **transom**.
///
/// The geometry contract puts the bow at the high-`x` end, so the transom, if
/// there is one, is the section at the hull's lowest `x`. `f_T(z) = f(x_T, z)`
/// is the aft end station's section.
///
/// Presence alone is not a warning: a transom clear of the water is simply a
/// closed hull as far as the wave integral is concerned. What matters is
/// [`Transom::depth`] and the area ratio against
/// [`SectionalHull::max_section_area`].
#[derive(Debug, Clone)]
pub struct Transom {
    /// Station of the transom — the aft end of the hull's x-domain [m].
    pub x: f64,
    /// Immersion depth [m], as the **equivalent rectangle**:
    /// `A_T / (2 · max_z f_T)` — the depth of a rectangle of the transom's
    /// widest beam carrying the same immersed area. Exact for a rectangular
    /// transom; `T/2` for one tapering linearly to the keel.
    ///
    /// Deliberately not a level crossing of `f_T`: an area measure is
    /// insensitive to a thin tail of beam running down the draft, and it is
    /// also the scale a closure wants, since what sets the hollow is how much
    /// water has to fill in behind the transom, not where its edge sits.
    pub depth: f64,
    /// Immersed transom area, `2∫₀^T f_T(z) dz` [m²] — both sides.
    pub area: f64,
    /// Its first moment about the waterline, `2∫₀^T z·f_T(z) dz` [m³] — both
    /// sides: `ρg` times it is the hydrostatic thrust the water would put on
    /// the transom's face.
    pub depth_moment: f64,
    /// Half-beam at the waterline, `f_T(0)` [m].
    pub half_beam: f64,
}

/// A hull as stations along `x`, each with its section's depth quadrature,
/// and the B-spline space along `x` its depth integrals are interpolated in.
#[derive(Debug, Clone)]
pub struct SectionalHull {
    p: usize,
    n: usize,
    sections: Vec<SectionNodes>,
    /// Banded LU of the collocation matrix (interpolation at the stations).
    lu: BandLu,
    spans: Vec<Span>,
    /// First basis index active on each span.
    span_first: Vec<usize>,
    /// Per span, `(p+1) × (p+1)`: `N_{first+r}^{(a)}(x_s) / a!` at `[a][r]` —
    /// maps B-spline coefficients to the span's local power basis.
    taylor: Vec<f64>,
    /// Waterline `f(x, 0)` and `∂f/∂x(x, 0)` in each span's power basis.
    wl_f: Vec<f64>,
    wl_fx: Vec<f64>,
    x_center: f64,
    length: f64,
    draft: f64,
    volume: f64,
    xs: Vec<f64>,
    lcb_x: f64,
    waterplane: [f64; 3],
    wetted_surface: f64,
    /// The aft end station's section, when it is a real (immersed) one: a
    /// transom, closed by [`TransomClosure`] as on any hull.
    transom: Option<Transom>,
}

/// A sectional hull's depth integrals at one `κ`, in each x-span's local
/// power basis (the layout the exact B-spline kernel's `ZContracted` uses).
#[derive(Debug, Clone, Default)]
pub struct SectionalContracted {
    g_f: Vec<f64>,
    g_fx: Vec<f64>,
    any: bool,
    /// The aft end station's depth integral: the transom section's
    /// z-factor for the closure.
    z_t: f64,
    /// The fore end station's depth integral: a blunt bow's step.
    z_b: f64,
}

impl SectionalHull {
    /// A sectional hull from stations `xs` (strictly increasing) and their
    /// sections, interpolated along `x` by a clamped B-spline of degree `p`
    /// on `knots` (which must satisfy Schoenberg–Whitney against `xs`: one
    /// station per basis function, each inside its support).
    pub fn new(p: usize, knots: Vec<f64>, xs: &[f64], sections: Vec<SectionNodes>) -> Result<Self> {
        let n = knots.len() - p - 1;
        if xs.len() != n || sections.len() != n {
            return Err(Error::InvalidInput(format!(
                "sectional hull: {} stations and {} sections for {n} basis functions",
                xs.len(),
                sections.len()
            )));
        }
        // Collocation: row i holds the basis functions alive at station i.
        let mut lu = BandLu::new(n, p);
        for (i, &x) in xs.iter().enumerate() {
            let span = find_span(&knots, p, n, x);
            let d = ders_basis(&knots, p, span, x, 0);
            for r in 0..=p {
                lu.set(i, span - p + r, d[0][r]);
            }
        }
        if !lu.factor() {
            return Err(Error::InvalidInput(
                "sectional hull: stations do not interpolate the x knots \
                 (Schoenberg–Whitney violated)"
                    .into(),
            ));
        }
        let mut spans = Vec::new();
        let mut span_first = Vec::new();
        let mut taylor = Vec::new();
        for s in p..n {
            let (a, b) = (knots[s], knots[s + 1]);
            if b <= a {
                continue;
            }
            spans.push(Span {
                start: a,
                len: b - a,
            });
            span_first.push(s - p);
            let d = ders_basis(&knots, p, s, a, p);
            let mut fact = 1.0;
            for k in 0..=p {
                if k > 0 {
                    fact *= k as f64;
                }
                for r in 0..=p {
                    taylor.push(d[k][r] / fact);
                }
            }
        }
        let mut hull = SectionalHull {
            p,
            n,
            sections,
            lu,
            spans,
            span_first,
            taylor,
            wl_f: Vec::new(),
            wl_fx: Vec::new(),
            x_center: 0.5 * (knots[p] + knots[n]),
            length: knots[n] - knots[p],
            draft: 0.0,
            volume: 0.0,
            xs: xs.to_vec(),
            lcb_x: 0.0,
            waterplane: [0.0; 3],
            wetted_surface: 0.0,
            transom: None,
        };
        hull.draft = hull.sections.iter().map(|s| s.depth).fold(0.0, f64::max);
        // Displaced volume 2∬ f = 2∫ Z(x; 0) dx and its x-moment, exactly
        // from the interpolant of the section areas.
        let mut z0 = SectionalContracted::default();
        hull.contract(0.0, &mut z0);
        hull.volume = 2.0 * hull.x_moment(&z0.g_f, 0);
        hull.lcb_x = if hull.volume > 0.0 {
            2.0 * hull.x_moment(&z0.g_f, 1) / hull.volume
        } else {
            0.0
        };
        let wl: Vec<f64> = hull.sections.iter().map(|s| s.waterline).collect();
        let (f, fx) = hull.to_power_basis(wl);
        hull.wl_f = f;
        hull.wl_fx = fx;
        // Waterplane: A_w = ∫ 2 f(x, 0) dx and its first and second moments
        // about x = 0, from the waterline interpolant.
        hull.waterplane = [0, 1, 2].map(|k| 2.0 * hull.x_moment(&hull.wl_f, k));
        hull.wetted_surface = hull.strip_area();
        hull.transom = hull.detect_transom();
        Ok(hull)
    }

    /// The transom, if the aft end station carries a real section: immersed area
    /// above `TRANSOM_AREA_REL` of the largest section's — with its depth as
    /// the same equivalent rectangle, `A_T / (2 · max f_T)`. The geometry
    /// contract puts the bow at high x, so this is the lowest-x station.
    ///
    /// Only the aft end is closed. A truncated **bow**
    /// (a forward end station with a real section) is left as it is: its
    /// section steps to nothing past the end of the interpolant, which the
    /// wave kernel reads as an open, bluff forward face — no virtual
    /// appendage is added ahead of it.
    fn detect_transom(&self) -> Option<Transom> {
        let aft = self.sections.first()?;
        let area = 2.0 * aft.integrate(0.0);
        let max_section = self
            .sections
            .iter()
            .map(|s| 2.0 * s.integrate(0.0))
            .fold(0.0, f64::max);
        if !(area > TRANSOM_AREA_REL * max_section) {
            return None;
        }
        let peak = aft
            .outline()
            .iter()
            .map(|p| p.0)
            .fold(aft.waterline, f64::max);
        let depth = if peak > 0.0 {
            (area / (2.0 * peak)).min(self.draft)
        } else {
            0.0
        };
        Some(Transom {
            x: self.xs[0],
            depth,
            area,
            depth_moment: 2.0 * aft.first_moment(),
            half_beam: aft.waterline,
        })
    }

    /// The transom the aft end presents, if any.
    pub fn transom(&self) -> Option<&Transom> {
        self.transom.as_ref()
    }

    /// Largest immersed section area `2∫ f dz` over the stations [m²] — the
    /// reference a transom's area is judged against.
    pub fn max_section_area(&self) -> f64 {
        self.sections
            .iter()
            .map(|s| 2.0 * s.integrate(0.0))
            .fold(0.0, f64::max)
    }

    /// `∫ x^k g(x) dx` over the hull for `g` given in each span's local power
    /// basis (`[s·(p+1) + a]`, local coordinate `t = x − start`): exact, by
    /// expanding `(start + t)^k` binomially.
    fn x_moment(&self, g: &[f64], k: usize) -> f64 {
        let p = self.p;
        let binom = |n: usize, r: usize| -> f64 {
            (0..r).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
        };
        self.spans
            .iter()
            .enumerate()
            .map(|(s, sx)| {
                let mut total = 0.0;
                for a in 0..=p {
                    let c = g[s * (p + 1) + a];
                    for j in 0..=k {
                        // ∫_0^len start^(k−j) t^(j+a) dt
                        let e = j + a + 1;
                        total +=
                            c * binom(k, j) * sx.start.powi((k - j) as i32) * sx.len.powi(e as i32)
                                / e as f64;
                    }
                }
                total
            })
            .sum()
    }

    /// Wetted surface of both sides [m²]: the shell between neighbouring
    /// stations as a strip of triangles joining their outlines node for
    /// node, which carries the slope along x a projected area would miss.
    /// An end station without a section closes its strip to a point.
    fn strip_area(&self) -> f64 {
        let tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            let (u, v) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
            );
            let w = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
        };
        let mut area = 0.0;
        for i in 0..self.xs.len().saturating_sub(1) {
            let (xa, xb) = (self.xs[i], self.xs[i + 1]);
            let (oa, ob) = (self.sections[i].curve(), self.sections[i + 1].curve());
            let n = oa.len().max(ob.len());
            if n < 2 {
                continue;
            }
            // Node k of an outline, stretched over n nodes (a point if empty).
            let at = |o: &[(f64, f64)], x: f64, k: usize| -> [f64; 3] {
                if o.is_empty() {
                    return [x, 0.0, 0.0];
                }
                let j = (k * (o.len() - 1)) / (n - 1);
                [x, o[j].0, o[j].1]
            };
            for k in 0..n - 1 {
                let (a0, a1) = (at(oa, xa, k), at(oa, xa, k + 1));
                let (b0, b1) = (at(ob, xb, k), at(ob, xb, k + 1));
                area += tri(a0, b0, b1) + tri(a0, b1, a1);
            }
        }
        2.0 * area
    }

    /// B-spline interpolant of per-station values, as power coefficients of
    /// the value (`[s·(p+1) + a]`) and its x-derivative (`[s·p + a]`) on
    /// each span.
    fn to_power_basis(&self, mut v: Vec<f64>) -> (Vec<f64>, Vec<f64>) {
        let mut f = vec![0.0; self.spans.len() * (self.p + 1)];
        let mut fx = vec![0.0; self.spans.len() * self.p];
        self.lu.solve(&mut v);
        self.power_into(&v, &mut f, &mut fx);
        (f, fx)
    }

    fn power_into(&self, c: &[f64], f: &mut [f64], fx: &mut [f64]) {
        let p = self.p;
        for (s, &first) in self.span_first.iter().enumerate() {
            let t = &self.taylor[s * (p + 1) * (p + 1)..(s + 1) * (p + 1) * (p + 1)];
            for a in 0..=p {
                let row = &t[a * (p + 1)..(a + 1) * (p + 1)];
                f[s * (p + 1) + a] = row
                    .iter()
                    .zip(&c[first..=first + p])
                    .map(|(m, c)| m * c)
                    .sum();
            }
            for a in 0..p {
                fx[s * p + a] = (a + 1) as f64 * f[s * (p + 1) + a + 1];
            }
        }
    }

    /// Everything that depends on `κ`: each station's depth integral, and
    /// their interpolant along `x` in each span's power basis.
    pub fn contract(&self, kappa: f64, out: &mut SectionalContracted) {
        self.contract_with_tangency(kappa, &[], out);
    }

    /// [`SectionalHull::contract`] with an optional boundary-layer tangency
    /// slope correction `df_dx` (`∂δ*/∂x` per station; empty when disabled)
    /// added to the longitudinal slope `g_fx` over each station's wetted
    /// depth, leaving the geometric depth integrals `g_f`, `z_t`, and `z_b`
    /// unchanged.
    pub fn contract_with_tangency(&self, kappa: f64, df_dx: &[f64], out: &mut SectionalContracted) {
        let mut v: Vec<f64> = self.sections.iter().map(|s| s.integrate(kappa)).collect();
        out.any = v.iter().any(|&x| x != 0.0);
        out.z_t = v.first().copied().unwrap_or(0.0);
        out.z_b = v.last().copied().unwrap_or(0.0);
        out.g_f.resize(self.spans.len() * (self.p + 1), 0.0);
        out.g_fx.resize(self.spans.len() * self.p, 0.0);
        if !out.any {
            out.g_f.fill(0.0);
            out.g_fx.fill(0.0);
            return;
        }
        self.lu.solve(&mut v);
        self.power_into(&v, &mut out.g_f, &mut out.g_fx);
        if df_dx.len() == self.sections.len() && df_dx.iter().any(|&x| x != 0.0) {
            let mut w: Vec<f64> = self
                .sections
                .iter()
                .zip(df_dx)
                .map(|(s, &slope)| {
                    if s.depth <= 0.0 || slope == 0.0 {
                        0.0
                    } else {
                        let kt = kappa * s.depth;
                        let iz = if kt.abs() < 1e-4 {
                            s.depth * (1.0 - 0.5 * kt + kt * kt / 6.0)
                        } else {
                            (1.0 - (-kt).exp()) / kappa
                        };
                        slope * iz
                    }
                })
                .collect();
            self.lu.solve(&mut w);
            let p = self.p;
            for (s, &first) in self.span_first.iter().enumerate() {
                let t = &self.taylor[s * (p + 1) * (p + 1)..(s + 1) * (p + 1) * (p + 1)];
                let mut c_span = vec![0.0; p + 1];
                for a in 0..=p {
                    let row = &t[a * (p + 1)..(a + 1) * (p + 1)];
                    c_span[a] = row
                        .iter()
                        .zip(&w[first..=first + p])
                        .map(|(m, c)| m * c)
                        .sum();
                }
                for a in 0..p {
                    out.g_fx[s * p + a] += c_span[a];
                }
                // Project the top degree-p term onto degree 0..p-1 via the
                // shifted Legendre polynomial on [0, L_s], preserving span
                // moments up to degree p-1.
                if p == 3 && c_span[3] != 0.0 {
                    let l = self.spans[s].len;
                    let c3 = c_span[3];
                    out.g_fx[s * p] += 0.05 * l * l * l * c3;
                    out.g_fx[s * p + 1] -= 0.6 * l * l * c3;
                    out.g_fx[s * p + 2] += 1.5 * l * c3;
                }
            }
        }
    }

    /// Boundary-layer displacement-thickness tangency slope `∂δ*/∂x` at each
    /// station for unit-length Reynolds number `re_per_m = U / ν` [1/m].
    ///
    /// Integrates the 1/7th-power-law turbulent boundary layer (`H = 9/7`,
    /// `dδ*/ds = 0.037 Re_s^(-1/5)`) from the bow (`x_max`) to the stern
    /// (`x_min`), amplified over the closing afterbody (`∂f̄/∂x > 0`) by the
    /// von Kármán momentum-integral adverse-pressure-gradient factor `(H + 2)`.
    pub fn boundary_layer_tangency(&self, re_per_m: f64) -> Vec<f64> {
        let n = self.xs.len();
        if n < 2 || !(re_per_m > 0.0) {
            return vec![0.0; n];
        }
        let f_bar: Vec<f64> = self
            .sections
            .iter()
            .map(|s| {
                if s.depth > 1e-6 {
                    s.integrate(0.0) / s.depth
                } else {
                    0.0
                }
            })
            .collect();
        let f_max = f_bar.iter().copied().fold(0.0f64, f64::max).max(1e-6);
        let x_bow = self.xs[n - 1];
        let l = (x_bow - self.xs[0]).max(1e-6);
        let mut df_dx = vec![0.0; n];
        let mut delta_star = 0.0;
        let s0 = 1e-4 * l;
        df_dx[n - 1] = -0.037 * (re_per_m * s0).max(1e3).powf(-0.2);
        // Shape factor H = 9/7 for 1/7th-power profile -> H + 2 = 23/7.
        const H_PLUS_2: f64 = 23.0 / 7.0;
        for i in (0..n - 1).rev() {
            let ds = (self.xs[i + 1] - self.xs[i]).max(1e-9);
            let s = (x_bow - self.xs[i]).max(1e-4 * l);
            let re_s = (re_per_m * s).max(1e3);
            let g_fp = 0.037 * re_s.powf(-0.2);
            let fx_closing = ((f_bar[i + 1] - f_bar[i]) / ds).max(0.0);
            let alpha = H_PLUS_2 * fx_closing / f_max;
            delta_star = (delta_star + ds * g_fp) / (1.0 - (alpha * ds).min(0.5));
            let d_delta_ds = g_fp + alpha * delta_star;
            df_dx[i] = -d_delta_ds;
        }
        df_dx
    }

    /// Length between the end stations' knots [m].
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Deepest depth node's section depth [m] (the deepest section).
    pub fn draft(&self) -> f64 {
        self.draft
    }

    /// Displaced volume `2∬ f dx dz` [m³].
    pub fn displaced_volume(&self) -> f64 {
        self.volume
    }

    /// Longitudinal centre of buoyancy [m], in the hull's x coordinates.
    pub fn lcb_x(&self) -> f64 {
        self.lcb_x
    }

    /// Waterplane area `∫ 2 f(x, 0) dx` [m²].
    pub fn waterplane_area(&self) -> f64 {
        self.waterplane[0]
    }

    /// First moment of the waterplane about x = 0 [m³].
    pub fn waterplane_moment(&self) -> f64 {
        self.waterplane[1]
    }

    /// Second moment of the waterplane about x = 0 [m⁴].
    pub fn waterplane_second_moment(&self) -> f64 {
        self.waterplane[2]
    }

    /// Transverse second moment of the waterplane about the centreplane,
    /// `I_T = ∫ ⅔ f(x, 0)³ dx` [m⁴]: with the displaced volume, the
    /// transverse metacentric radius `BM_T = I_T/∇` — the roll stiffness.
    /// Exact for the waterline interpolant (Gauss per span).
    pub fn transverse_waterplane_inertia(&self) -> f64 {
        let (gx, gw) = gauss_legendre(8);
        self.spans
            .iter()
            .map(|sp| {
                gx.iter()
                    .zip(&gw)
                    .map(|(&t, &w)| {
                        let x = sp.start + 0.5 * sp.len * (1.0 + t);
                        0.5 * sp.len * w * 2.0 / 3.0 * self.waterline_half_beam(x).powi(3)
                    })
                    .sum::<f64>()
            })
            .sum()
    }

    /// Depth of the centre of buoyancy below the waterline [m]: each
    /// station's section area and its centroid from the section curve,
    /// closed along the centreplane and the waterline, integrated along the
    /// stations.
    pub fn buoyancy_depth(&self) -> f64 {
        let per_station: Vec<(f64, f64, f64)> = self
            .xs
            .iter()
            .zip(&self.sections)
            .map(|(&x, sec)| {
                // Shoelace over (y, depth), the curve closed back through
                // the centreplane at the waterline.
                let c = sec.curve();
                let mut pts: Vec<(f64, f64)> = c.to_vec();
                if let (Some(first), Some(last)) = (c.first(), c.last()) {
                    pts.push((0.0, last.1));
                    pts.push((0.0, first.1.min(0.0)));
                }
                let (mut a, mut m) = (0.0, 0.0);
                for i in 0..pts.len() {
                    let (p0, p1) = (pts[i], pts[(i + 1) % pts.len()]);
                    let cross = p0.0 * p1.1 - p1.0 * p0.1;
                    a += 0.5 * cross;
                    m += cross * (p0.1 + p1.1) / 6.0;
                }
                (x, a.abs(), if a != 0.0 { m / a } else { 0.0 })
            })
            .collect();
        let (mut vol, mut mom) = (0.0, 0.0);
        for w in per_station.windows(2) {
            let dx = w[1].0 - w[0].0;
            vol += 0.5 * dx * (w[0].1 + w[1].1);
            mom += 0.5 * dx * (w[0].1 * w[0].2 + w[1].1 * w[1].2);
        }
        if vol > 0.0 {
            mom / vol
        } else {
            0.0
        }
    }

    /// Wetted surface of both sides [m²] (see [`SectionalHull`]'s strip
    /// construction): the shell itself, not its centreplane projection.
    pub fn wetted_surface(&self) -> f64 {
        self.wetted_surface
    }

    /// Centre of the x domain, the phase reference of the transforms [m].
    pub fn x_center(&self) -> f64 {
        self.x_center
    }

    /// The x domain: the end stations [m].
    pub fn x_range(&self) -> (f64, f64) {
        (
            self.spans.first().map_or(0.0, |s| s.start),
            self.spans.last().map_or(0.0, |s| s.start + s.len),
        )
    }

    /// Waterline half-beam `f(x, 0)` [m] as the kernel interpolates it; 0
    /// outside the hull.
    pub fn waterline_half_beam(&self, x: f64) -> f64 {
        let p = self.p;
        let Some(s) = self
            .spans
            .iter()
            .position(|sp| x >= sp.start && x <= sp.start + sp.len)
        else {
            return 0.0;
        };
        let t = x - self.spans[s].start;
        (0..=p)
            .rev()
            .fold(0.0, |acc, a| acc * t + self.wl_f[s * (p + 1) + a])
            .max(0.0)
    }

    /// The stations: each one's x and its section at the quadrature nodes
    /// (see [`SectionNodes::outline`]).
    pub fn sections(&self) -> impl Iterator<Item = (f64, &[(f64, f64)])> {
        self.xs
            .iter()
            .copied()
            .zip(self.sections.iter().map(|s| s.outline()))
    }

    /// The stations: each one's x and its section curve, sampled densely
    /// (see [`SectionNodes::curve`]).
    pub fn curves(&self) -> impl Iterator<Item = (f64, &[(f64, f64)])> {
        self.xs
            .iter()
            .copied()
            .zip(self.sections.iter().map(|s| s.curve()))
    }

    /// Each station's depth integral `Z(xᵢ; κ)`, and the interpolant the
    /// kernel actually integrates along x, sampled `per_span` times per
    /// x-span: `(stations, curve)` as `(x, Z)` pairs. At `κ = 0` this is the
    /// sectional-area curve (half-areas).
    pub fn depth_integral_curve(
        &self,
        kappa: f64,
        per_span: usize,
    ) -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
        let at_st = self
            .xs
            .iter()
            .zip(&self.sections)
            .map(|(&x, s)| (x, s.integrate(kappa)))
            .collect();
        let mut zc = SectionalContracted::default();
        self.contract(kappa, &mut zc);
        let p = self.p;
        let mut curve = Vec::new();
        for (s, sx) in self.spans.iter().enumerate() {
            for k in 0..per_span {
                let t = sx.len * k as f64 / per_span as f64;
                let v = if zc.any {
                    (0..=p)
                        .rev()
                        .fold(0.0, |acc, a| acc * t + zc.g_f[s * (p + 1) + a])
                } else {
                    0.0
                };
                curve.push((sx.start + t, v));
            }
        }
        // Close the curve at the last station.
        if let (Some(&xe), Some(se)) = (self.xs.last(), self.sections.last()) {
            curve.push((xe, se.integrate(kappa)));
        }
        (at_st, curve)
    }

    /// The x-derivative of the interpolated depth integral, `∂Z/∂x(x; κ)`,
    /// sampled `per_span` times per x-span as `(x, ∂Z/∂x)` — the source
    /// strength the kernel integrates (the Michell amplitude is its Fourier
    /// transform along x).
    pub fn depth_integral_slope_curve(&self, kappa: f64, per_span: usize) -> Vec<(f64, f64)> {
        let mut zc = SectionalContracted::default();
        self.contract(kappa, &mut zc);
        let p = self.p;
        let mut curve = Vec::new();
        for (s, sx) in self.spans.iter().enumerate() {
            for k in 0..=per_span {
                if k == per_span && s + 1 < self.spans.len() {
                    continue;
                }
                let t = sx.len * k as f64 / per_span as f64;
                let v = if zc.any {
                    (0..p)
                        .rev()
                        .fold(0.0, |acc, a| acc * t + zc.g_fx[s * p + a])
                } else {
                    0.0
                };
                curve.push((sx.start + t, v));
            }
        }
        curve
    }

    /// Number of stations.
    pub fn stations(&self) -> usize {
        self.n
    }

    /// Total depth nodes over all stations.
    pub fn depth_nodes(&self) -> usize {
        self.sections.iter().map(SectionNodes::len).sum()
    }

    /// Number of x-spans (what each `k_x` evaluation walks).
    pub fn x_spans(&self) -> usize {
        self.spans.len()
    }
    /// Replace the transom description — for a reference construction whose
    /// own geometry describes its transom better than the end station does
    /// (the exact B-spline test oracle).
    pub fn with_transom(mut self, transom: Option<Transom>) -> Self {
        self.transom = transom;
        self
    }

    /// Degree `p` of the B-spline along x.
    pub fn degree(&self) -> usize {
        self.p
    }

    /// The non-empty x knot spans, in order: the local coordinate of every
    /// power-basis layout below is `t = x − start` on its span.
    pub fn spans(&self) -> &[Span] {
        &self.spans
    }

    /// Waterline `f(x, 0)` in each span's power basis (`[s·(p+1) + a]`).
    pub fn waterline_power(&self) -> &[f64] {
        &self.wl_f
    }

    /// Waterline slope `∂f/∂x(x, 0)` in each span's power basis (`[s·p + a]`).
    pub fn waterline_slope_power(&self) -> &[f64] {
        &self.wl_fx
    }

    /// The stations' x, strictly increasing: the aft end first.
    pub fn station_xs(&self) -> &[f64] {
        &self.xs
    }

    /// The stations' sections, in station order.
    pub fn section_nodes(&self) -> &[SectionNodes] {
        &self.sections
    }

    /// `(∫ g e^{ik_x (x − x_c)} dx, ∫ (x − x_c) g e^{ik_x (x − x_c)} dx)` over
    /// the hull for `g` in each span's power basis of degree `deg`
    /// (`[s·(deg+1) + a]`): closed form per span.
    fn power_transforms(&self, g: &[f64], deg: usize, kx: f64) -> (C64, C64) {
        let mut xm = Vec::with_capacity(deg + 2);
        let (mut total, mut first) = (C64::ZERO, C64::ZERO);
        for (s, sx) in self.spans.iter().enumerate() {
            osc_moments(kx, sx.len, deg + 1, &mut xm);
            let d = sx.start - self.x_center;
            let (mut sum, mut sum1) = (C64::ZERO, C64::ZERO);
            for a in 0..=deg {
                let c = g[s * (deg + 1) + a];
                sum = sum + xm[a].scale(c);
                sum1 = sum1 + (xm[a + 1] + xm[a].scale(d)).scale(c);
            }
            let phase = C64::cis(kx * d);
            total = total + phase * sum;
            first = first + phase * sum1;
        }
        (total, first)
    }

    /// `∫ Z(x; κ) e^{ik_x (x − x_c)} dx` for a contraction at some `κ`, and
    /// its first moment `∫ (x − x_c) Z(x; κ) e^{ik_x (x − x_c)} dx`: the
    /// hull's half-volume weighted by `e^{−κz}` and a longitudinal wave. At
    /// `κ = k`, `k_x = k cos β` these are the depth-decayed incident-wave
    /// integrals a Froude–Krylov force and moment need.
    pub fn x_transform(&self, zc: &SectionalContracted, kx: f64) -> (C64, C64) {
        if !zc.any {
            return (C64::ZERO, C64::ZERO);
        }
        self.power_transforms(&zc.g_f, self.p, kx)
    }

    /// `∫ f(x, 0) e^{ik_x (x − x_c)} dx` and its first moment about `x_c`:
    /// the waterline half-beam's Fourier transform.
    pub fn waterline_transform(&self, kx: f64) -> (C64, C64) {
        self.power_transforms(&self.wl_f, self.p, kx)
    }
}

impl SectionNodes {
    /// Half-beam at the waterline, `f(x_i, 0)` [m].
    pub fn waterline(&self) -> f64 {
        self.waterline
    }

    /// Depth of the section's lowest point [m].
    pub fn depth(&self) -> f64 {
        self.depth
    }

    /// The station's depth integral `Z(κ) = ∫ f(x_i, z) e^{−κz} dz`.
    pub fn depth_integral(&self, kappa: f64) -> f64 {
        self.integrate(kappa)
    }
}

impl SectionalContracted {
    /// `Z(x; κ)` in each span's power basis (`[s·(p+1) + a]`).
    pub fn g_f(&self) -> &[f64] {
        &self.g_f
    }

    /// `∂Z/∂x(x; κ)` in each span's power basis (`[s·p + a]`).
    pub fn g_fx(&self) -> &[f64] {
        &self.g_fx
    }

    /// Whether any station's depth integral is non-zero (a dry hull's is not).
    pub fn any(&self) -> bool {
        self.any
    }

    /// The aft end station's depth integral (a transom section's z-factor).
    pub fn aft(&self) -> f64 {
        self.z_t
    }

    /// The fore end station's depth integral (a blunt bow's step).
    pub fn fore(&self) -> f64 {
        self.z_b
    }
}

/// LU of a banded matrix without pivoting — B-spline collocation matrices
/// are totally positive, for which elimination without pivoting is stable
/// (de Boor). Dense storage, band-limited loops: stations number in the
/// hundreds.
#[derive(Debug, Clone)]
struct BandLu {
    n: usize,
    /// Half-bandwidth: entries with `|i − j| > w` are zero.
    w: usize,
    a: Vec<f64>,
}

impl BandLu {
    fn new(n: usize, p: usize) -> Self {
        BandLu {
            n,
            w: p + 1,
            a: vec![0.0; n * n],
        }
    }

    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] = v;
    }

    fn factor(&mut self) -> bool {
        let (n, w) = (self.n, self.w);
        for k in 0..n {
            let piv = self.a[k * n + k];
            if piv.abs() < 1e-14 {
                return false;
            }
            for i in k + 1..(k + w + 1).min(n) {
                let l = self.a[i * n + k] / piv;
                if l == 0.0 {
                    continue;
                }
                self.a[i * n + k] = l;
                for j in k + 1..(k + w + 1).min(n) {
                    self.a[i * n + j] -= l * self.a[k * n + j];
                }
            }
        }
        true
    }

    fn solve(&self, b: &mut [f64]) {
        let (n, w) = (self.n, self.w);
        for i in 0..n {
            let mut s = b[i];
            for k in i.saturating_sub(w)..i {
                s -= self.a[i * n + k] * b[k];
            }
            b[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = b[i];
            for j in i + 1..(i + w + 1).min(n) {
                s -= self.a[i * n + j] * b[j];
            }
            b[i] = s / self.a[i * n + i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{DepthQuadrature, SectionNodes};
    use crate::iges::{self, HullPose, Platform, SectionalOptions};

    /// A section's first moment of area, `∫ z f dz`, from either quadrature:
    /// `bT²/2` for a rectangle by depth, `bT²/3` for a quarter ellipse by
    /// rays.
    #[test]
    fn section_first_moments_are_exact() {
        let (b, t) = (0.4, 0.3);
        let opts = DepthQuadrature::default();
        let rect = SectionNodes::from_depth_function(t, &[], |_| b, &opts);
        let ellipse = SectionNodes::from_polar(0.0, b, t, |_| 1.0, &opts);
        for (what, s, want) in [
            ("rectangle", rect, b * t * t / 2.0),
            ("ellipse", ellipse, b * t * t / 3.0),
        ] {
            let m = s.first_moment();
            assert!((m - want).abs() < 1e-10 * want, "{what}: {m} vs {want}");
        }
    }

    /// The Wigley hull `f = (B/2)(1 − ξ²)(1 − ζ²)`: `I_T = (4/105) B³ L`
    /// and the centre of buoyancy `3T/8` below the waterline.
    #[test]
    fn wigley_roll_hydrostatics_are_exact() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let surfaces = iges::wigley_surfaces(l, b, t).unwrap();
        let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
        let hull = source
            .situate_sectional(
                0,
                0.0,
                &HullPose::default(),
                &Platform::default(),
                &SectionalOptions::default(),
            )
            .unwrap()
            .unwrap()
            .hull;
        let it = hull.transverse_waterplane_inertia();
        let want = 4.0 / 105.0 * b * b * b * l;
        assert!((it - want).abs() < 2e-3 * want, "I_T {it} vs {want}");
        let zb = hull.buoyancy_depth();
        assert!(
            (zb - 0.375 * t).abs() < 2e-3 * t,
            "z_B {zb} vs {}",
            0.375 * t
        );
    }
}
