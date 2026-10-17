//! Test oracle: a validated B-spline half-breadth surface `y = f(x, z)` plus
//! the span polynomials the exact kernel (`crate::michell::InnerIntegral`)
//! integrates, precomputed once. The sectional kernel is checked against it,
//! through [`FromHull::from_hull`].

use hullgeom::bspline::BSplineSurface;
use hullgeom::quadrature::gauss_legendre;
use hullgeom::sectional::{
    DepthQuadrature, SectionNodes, SectionalHull, Span, Transom, TRANSOM_AREA_REL,
};
use hullgeom::{Error, Result};

/// The oracle's transom: the geometry crate's description plus the exact
/// section polynomial the B-spline kernel integrates.
#[derive(Debug, Clone)]
pub struct HullTransom {
    pub section: Transom,
    /// Per-z-span polynomial coefficients `f_T(z) = Σ_b c[b] (z − z0_sz)^b`,
    /// flattened as `[sz * (q + 1) + b]`.
    pub coeff: Vec<f64>,
}

impl std::ops::Deref for HullTransom {
    type Target = Transom;
    fn deref(&self) -> &Transom {
        &self.section
    }
}

/// The reference construction of a sectional hull from the oracle.
pub(crate) trait FromHull: Sized {
    fn from_hull(hull: &Hull, opts: &DepthQuadrature) -> Result<Self>;
}

impl FromHull for SectionalHull {
    /// Sections of a B-spline hull at the Greville points of its own x
    /// knots, depth-integrated panel by panel between its z knots.
    /// `Z(x; κ)` is then exactly in the interpolating space, so the only
    /// approximation is the depth quadrature. The spline's own transom
    /// description (its depth measure samples the spline, not the
    /// quadrature nodes) replaces the end station's, so the harness checks
    /// the closure machinery exactly.
    fn from_hull(hull: &Hull, opts: &DepthQuadrature) -> Result<Self> {
        let s = hull.surface();
        let (p, n) = (s.degree_x(), s.n_ctrl_x());
        let kx = s.knots_x();
        let xs: Vec<f64> = (0..n)
            .map(|i| kx[i + 1..=i + p].iter().sum::<f64>() / p.max(1) as f64)
            .collect();
        let depth = s.z_domain().1;
        let sections = xs
            .iter()
            .map(|&x| SectionNodes::from_depth_function(depth, s.knots_z(), |z| s.eval(x, z), opts))
            .collect();
        Ok(SectionalHull::new(p, kx.to_vec(), &xs, sections)?
            .with_transom(hull.transom().map(|t| t.section.clone())))
    }
}

/// A validated hull.
///
/// Geometry contract (see crate docs): `y = f(x, z) >= 0` is the local
/// half-beam of a port/starboard-symmetric hull; `x` runs along the hull
/// (arbitrary origin), `z` runs vertically **downward** from the undisturbed
/// waterline at `z = 0`; the surface domain is `x ∈ [x0, x1]`, `z ∈ [0, T]`.
#[derive(Debug, Clone)]
pub struct Hull {
    surface: BSplineSurface,
    spans_x: Vec<Span>,
    spans_z: Vec<Span>,
    /// Local polynomial coefficients of ∂f/∂x per span pair:
    /// `fx(x, z) = Σ_{a,b} c[a][b] (x - x0)^a (z - z0)^b`,
    /// flattened as `[(sx * spans_z.len() + sz) * p + a][(q+1) elements over b]`,
    /// i.e. index `((sx * nsz + sz) * p + a) * (q + 1) + b`,
    /// with `a = 0..p` (p = degree_x) and `b = 0..=q` (q = degree_z).
    fx_coeff: Vec<f64>,
    /// Local polynomial coefficients of `f` itself per span pair, same
    /// indexing as `fx_coeff` but with `a = 0..=p` (one more x power), i.e.
    /// index `((sx * nsz + sz) * (p + 1) + a) * (q + 1) + b`. The near-field
    /// (sinkage/trim) transforms need moments of `f` and `x·∂f/∂x`, not only
    /// of `∂f/∂x`; carrying `f` per span keeps them exact too.
    f_coeff: Vec<f64>,
    length: f64,
    draft: f64,
    x_center: f64,
    wetted_surface: f64,
    displaced_volume: f64,
    lcb_x: f64,
    waterplane_area: f64,
    waterplane_moment: f64,
    waterplane_second_moment: f64,
    max_section_area: f64,
    transom: Option<HullTransom>,
}

impl Hull {
    /// Validate the surface as a hull and precompute span polynomials and
    /// geometric integrals.
    ///
    /// Beyond spline validity this requires:
    /// - the z-domain starts at the waterline: `z0 = 0` (small negative
    ///   tolerance rejected) and has positive draft;
    /// - a non-negative control net (a conservative sufficient condition for
    ///   `f >= 0`, by the B-spline convex-hull property).
    pub fn new(surface: BSplineSurface) -> Result<Hull> {
        let (z0, z1) = surface.z_domain();
        let (x0, x1) = surface.x_domain();
        let draft = z1;
        if z0 < 0.0 || z0 > 1e-9 * (z1 - z0) {
            return Err(Error::InvalidGeometry(format!(
                "z-domain must start at the waterline z = 0 (z measured downward); got z0 = {z0}"
            )));
        }
        let scale = surface
            .control()
            .iter()
            .fold(0.0f64, |m, &v| m.max(v.abs()));
        if surface
            .control()
            .iter()
            .any(|&v| v < -1e-12 * scale.max(1.0))
        {
            return Err(Error::InvalidGeometry(
                "control net contains negative half-beam values; the half-breadth \
                 surface must satisfy f(x, z) >= 0"
                    .into(),
            ));
        }

        let p = surface.degree_x();
        let q = surface.degree_z();
        let xs: Vec<usize> = surface.x_span_indices();
        let zs: Vec<usize> = surface.z_span_indices();
        let spans_x: Vec<Span> = xs
            .iter()
            .map(|&s| Span {
                start: surface.knots_x()[s],
                len: surface.knots_x()[s + 1] - surface.knots_x()[s],
            })
            .collect();
        let spans_z: Vec<Span> = zs
            .iter()
            .map(|&s| Span {
                start: surface.knots_z()[s],
                len: surface.knots_z()[s + 1] - surface.knots_z()[s],
            })
            .collect();

        // Local polynomial coefficients of fx = ∂f/∂x on every span pair.
        let (fx_coeff, f_coeff) = compute_span_coeffs(&surface, &xs, &zs);

        // Geometric integrals by per-span Gauss-Legendre.
        // Volume: integrand is polynomial of degree (p, q) => exact.
        // Wetted surface: smooth integrand, use a generous rule.
        // x·f and z·f raise one degree by one; +2 keeps the rule exact.
        let n_vol = (p.max(q) + 2) / 2 + 2;
        // Wetted-surface integrand is smooth but non-polynomial, so it needs a
        // genuine rule rather than an exact one. What sets a *composite*
        // rule's accuracy is the total node count along each axis, not the
        // count per span — so spend a fixed budget across however many spans
        // the loft produced. A single-span hull keeps the original 24-point
        // rule; a finely-lofted one drops to a handful per span, which is what
        // keeps `Hull::new` affordable once a dense net puts a few thousand
        // span pairs here (the integrand is a polynomial *within* each span,
        // so a few nodes on a short span are already far past converged).
        const WET_NODES_PER_AXIS: usize = 200;
        let n_wet_x = WET_NODES_PER_AXIS
            .div_ceil(spans_x.len().max(1))
            .clamp(4, 24);
        let n_wet_z = WET_NODES_PER_AXIS
            .div_ceil(spans_z.len().max(1))
            .clamp(4, 24);
        let (xv, wv) = gauss_legendre(n_vol);
        let (xwx, wwx) = gauss_legendre(n_wet_x);
        let (xwz, wwz) = gauss_legendre(n_wet_z);
        let mut volume = 0.0;
        let mut volume_mx = 0.0;
        let mut wetted = 0.0;
        for sx in &spans_x {
            for sz in &spans_z {
                let jac = sx.len / 2.0 * (sz.len / 2.0);
                for (i, &xi) in xv.iter().enumerate() {
                    let x = sx.start + sx.len * (xi + 1.0) / 2.0;
                    for (j, &zj) in xv.iter().enumerate() {
                        let z = zj_map(sz, zj);
                        let f = surface.eval(x, z);
                        volume += wv[i] * wv[j] * jac * f;
                        volume_mx += wv[i] * wv[j] * jac * x * f;
                    }
                }
                for (i, &xi) in xwx.iter().enumerate() {
                    let x = sx.start + sx.len * (xi + 1.0) / 2.0;
                    for (j, &zj) in xwz.iter().enumerate() {
                        let z = zj_map(sz, zj);
                        let fx = surface.eval_deriv(x, z, 1, 0);
                        let fz = surface.eval_deriv(x, z, 0, 1);
                        wetted += wwx[i] * wwz[j] * jac * (1.0 + fx * fx + fz * fz).sqrt();
                    }
                }
            }
        }
        // Both sides of the hull.
        volume *= 2.0;
        volume_mx *= 2.0;
        wetted *= 2.0;

        // Waterplane properties: 1-D integrals of the beam b(x) = 2 f(x, 0)
        // (the rule is exact to degree 3p + 3, which covers x²·f at p + 2).
        let (xq, wq) = gauss_legendre(3 * p / 2 + 2);
        let mut wp_area = 0.0;
        let mut wp_mx = 0.0;
        let mut wp_ixx = 0.0;
        for sx in &spans_x {
            let jac = sx.len / 2.0;
            for (i, &xi) in xq.iter().enumerate() {
                let x = sx.start + sx.len * (xi + 1.0) / 2.0;
                let f = surface.eval(x, z0.max(0.0));
                let b = 2.0 * f;
                wp_area += wq[i] * jac * b;
                wp_mx += wq[i] * jac * x * b;
                wp_ixx += wq[i] * jac * x * x * b;
            }
        }

        let max_section_area = max_section_area_of(&surface, &spans_x, &spans_z);
        let transom = detect_transom(&surface, &xs, &zs, &spans_z, draft, max_section_area);

        Ok(Hull {
            surface,
            spans_x,
            spans_z,
            fx_coeff,
            f_coeff,
            length: x1 - x0,
            draft,
            x_center: 0.5 * (x0 + x1),
            wetted_surface: wetted,
            displaced_volume: volume,
            lcb_x: if volume > 0.0 {
                volume_mx / volume
            } else {
                0.0
            },
            waterplane_area: wp_area,
            waterplane_moment: wp_mx,
            waterplane_second_moment: wp_ixx,
            max_section_area,
            transom,
        })
    }

    pub fn surface(&self) -> &BSplineSurface {
        &self.surface
    }

    /// Length of the x-domain [m]. Used as the reference length.
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Depth of the z-domain below the waterline [m].
    pub fn draft(&self) -> f64 {
        self.draft
    }

    /// Wetted surface area of both sides [m²], computed on the centerplane
    /// projection: `S = 2 ∬ √(1 + fx² + fz²) dx dz` (consistent with
    /// thin-ship theory; no girth correction).
    pub fn wetted_surface(&self) -> f64 {
        self.wetted_surface
    }

    /// Displaced volume `∇ = 2 ∬ f dx dz` [m³] (thin-ship approximation).
    pub fn displaced_volume(&self) -> f64 {
        self.displaced_volume
    }

    /// Longitudinal centre of buoyancy `x_B = ∬ x f / ∬ f` [m], in the hull's
    /// x coordinates.
    pub fn lcb_x(&self) -> f64 {
        self.lcb_x
    }

    /// Waterplane area `A_w = ∫ 2 f(x, 0) dx` [m²].
    pub fn waterplane_area(&self) -> f64 {
        self.waterplane_area
    }

    /// First moment of the waterplane about x = 0: `∫ 2 x f(x, 0) dx` [m³].
    pub fn waterplane_moment(&self) -> f64 {
        self.waterplane_moment
    }

    /// Second moment of the waterplane about x = 0: `∫ 2 x² f(x, 0) dx` [m⁴].
    pub fn waterplane_second_moment(&self) -> f64 {
        self.waterplane_second_moment
    }

    /// Longitudinal centre of flotation (waterplane centroid) [m].
    pub fn lcf_x(&self) -> f64 {
        if self.waterplane_area > 0.0 {
            self.waterplane_moment / self.waterplane_area
        } else {
            0.0
        }
    }

    /// Maximum immersed section area `A_X = max_x 2∫₀^T f(x, z) dz` [m²].
    ///
    /// The reference area the transom is judged against: `A_T/A_X` is the
    /// standard measure of how much of a transom-stern vessel this is.
    pub fn max_section_area(&self) -> f64 {
        self.max_section_area
    }

    /// The hull's transom, if its half-breadth does not close at the aft
    /// (low-`x`) end. `None` for a hull that closes there — which is what
    /// classical Michell theory assumes, and the only case this crate's wave
    /// integral currently models.
    pub fn transom(&self) -> Option<&HullTransom> {
        self.transom.as_ref()
    }

    pub(crate) fn spans_x(&self) -> &[Span] {
        &self.spans_x
    }

    pub(crate) fn spans_z(&self) -> &[Span] {
        &self.spans_z
    }

    pub(crate) fn fx_coeff(&self) -> &[f64] {
        &self.fx_coeff
    }

    pub(crate) fn f_coeff(&self) -> &[f64] {
        &self.f_coeff
    }

    /// Half-extent of the hull about its x-midpoint (bandwidth of the
    /// oscillatory inner integral).
    pub(crate) fn x_half_extent(&self) -> f64 {
        self.length / 2.0
    }

    pub(crate) fn x_center(&self) -> f64 {
        self.x_center
    }
}

#[inline]
fn zj_map(sz: &Span, node: f64) -> f64 {
    sz.start + sz.len * (node + 1.0) / 2.0
}

/// Local polynomial coefficients of `∂f/∂x` and of `f` per span pair, in the
/// flattened layouts documented on [`Hull::fx_coeff`] and [`Hull::f_coeff`].
/// Both come from one corner-Taylor pass.
fn compute_span_coeffs(
    surface: &BSplineSurface,
    xs: &[usize],
    zs: &[usize],
) -> (Vec<f64>, Vec<f64>) {
    let p = surface.degree_x();
    let q = surface.degree_z();
    // Factorials up to max degree (degrees are small).
    let mut fact = vec![1.0f64; p.max(q) + 2];
    for i in 1..fact.len() {
        fact[i] = fact[i - 1] * i as f64;
    }
    let nsz = zs.len();
    let mut fx_coeff = vec![0.0f64; xs.len() * nsz * p * (q + 1)];
    let mut f_coeff = vec![0.0f64; xs.len() * nsz * (p + 1) * (q + 1)];
    for (isx, &sx) in xs.iter().enumerate() {
        for (isz, &sz) in zs.iter().enumerate() {
            let d = surface.corner_partials(sx, sz);
            for b in 0..=q {
                // f = Σ D[a][b]/(a! b!) X^a Z^b  =>
                // fx coefficient of X^a Z^b is D[a+1][b]/(a! b!).
                for a in 0..p {
                    fx_coeff[((isx * nsz + isz) * p + a) * (q + 1) + b] =
                        d[a + 1][b] / (fact[a] * fact[b]);
                }
                for a in 0..=p {
                    f_coeff[((isx * nsz + isz) * (p + 1) + a) * (q + 1) + b] =
                        d[a][b] / (fact[a] * fact[b]);
                }
            }
        }
    }
    (fx_coeff, f_coeff)
}

/// Maximum immersed section area `A_X = max_x 2∫₀^T f(x, z) dz`.
///
/// Exact in z (the integrand is polynomial of degree `q` on each z-span); the
/// maximum over x is located by sampling each x-span, which is ample for a
/// quantity that only ever appears as the denominator of a ratio.
fn max_section_area_of(surface: &BSplineSurface, spans_x: &[Span], spans_z: &[Span]) -> f64 {
    // As in `Hull::new`: a fixed number of stations over the whole hull, not
    // per span, so a dense loft does not pay for a scan it does not need.
    const SECTION_STATIONS: usize = 256;
    let samples_per_span = SECTION_STATIONS.div_ceil(spans_x.len().max(1)).clamp(2, 16);
    let q = surface.degree_z();
    let (zn, zw) = gauss_legendre(q / 2 + 2);
    let mut best = 0.0f64;
    for sx in spans_x {
        for i in 0..=samples_per_span {
            let x = sx.start + sx.len * i as f64 / samples_per_span as f64;
            let mut area = 0.0;
            for sz in spans_z {
                let jac = sz.len / 2.0;
                for (j, &zj) in zn.iter().enumerate() {
                    area += zw[j] * jac * surface.eval(x, zj_map(sz, zj));
                }
            }
            best = best.max(2.0 * area);
        }
    }
    best
}

/// Detect a transom: a non-closing half-breadth at the aft (low-`x`) end.
///
/// The section polynomial comes from the same corner-Taylor data the span
/// derivatives use — at `X = 0` on the first x-span, `f_T`'s coefficient of
/// `Z^b` is `D[0][b]/b!` — so the transom is carried in exactly the form the
/// free-wave kernel already integrates exactly.
fn detect_transom(
    surface: &BSplineSurface,
    xs: &[usize],
    zs: &[usize],
    spans_z: &[Span],
    draft: f64,
    max_section_area: f64,
) -> Option<HullTransom> {
    let q = surface.degree_z();
    let x_t = surface.x_domain().0;

    let mut fact = vec![1.0f64; q + 2];
    for i in 1..fact.len() {
        fact[i] = fact[i - 1] * i as f64;
    }
    let mut coeff = vec![0.0f64; zs.len() * (q + 1)];
    for (isz, &sz) in zs.iter().enumerate() {
        let d = surface.corner_partials(xs[0], sz);
        for b in 0..=q {
            coeff[isz * (q + 1) + b] = d[0][b] / fact[b];
        }
    }

    // Area and its first moment: exact per z-span (polynomials of degree q
    // and q + 1), both sides.
    let (zn, zw) = gauss_legendre(q / 2 + 2);
    let (mut area, mut depth_moment) = (0.0, 0.0);
    for sz in spans_z {
        let jac = sz.len / 2.0;
        for (j, &zj) in zn.iter().enumerate() {
            let z = zj_map(sz, zj);
            let f = zw[j] * jac * surface.eval(x_t, z);
            area += f;
            depth_moment += f * z;
        }
    }
    area *= 2.0;
    depth_moment *= 2.0;
    if !(area > TRANSOM_AREA_REL * max_section_area) {
        return None;
    }

    // Equivalent-rectangle immersion depth (see `Transom::depth`). Normalised
    // on the transom's largest half-beam rather than f_T(0), so a section
    // carrying no beam right at the waterline still gets a sane depth.
    let n = 16 * spans_z.len();
    let peak = (0..=n).fold(0.0f64, |m, i| {
        m.max(surface.eval(x_t, draft * i as f64 / n as f64))
    });
    let depth = if peak > 0.0 {
        (area / (2.0 * peak)).min(draft)
    } else {
        0.0
    };

    Some(HullTransom {
        section: Transom {
            x: x_t,
            depth,
            area,
            depth_moment,
            half_beam: surface.eval(x_t, 0.0),
        },
        coeff,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The transom section polynomial is what the virtual-appendage closure
    /// will integrate, so it has to reproduce the surface along `x = x_T`
    /// span for span — including across an interior knot.
    #[test]
    fn transom_coefficients_reproduce_the_section() {
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 3.0, 6.0, 9.0, 9.0, 9.0, 9.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        let (nx, nz) = (6usize, 4usize);
        // A net well clear of zero at the aft end (i = 0) tapering to a closed
        // bow: a genuine transom-sterned hull.
        let mut control = vec![0.0; nx * nz];
        for i in 0..nx {
            for j in 0..nz {
                let taper = 1.0 - i as f64 / (nx - 1) as f64;
                control[i * nz + j] = 0.6 * taper * (1.0 - 0.25 * j as f64);
            }
        }
        let surface = BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        let hull = Hull::new(surface).unwrap();
        let tr = hull.transom().expect("tapered-to-bow hull has a transom");

        let q = hull.surface().degree_z();
        for (isz, sz) in hull.spans_z().iter().enumerate() {
            for k in 0..=8 {
                let z = sz.start + sz.len * k as f64 / 8.0;
                let dz = z - sz.start;
                let row = &tr.coeff[isz * (q + 1)..(isz + 1) * (q + 1)];
                let poly = row.iter().rev().fold(0.0, |acc, &c| acc * dz + c);
                let exact = hull.surface().eval(tr.x, z);
                assert!(
                    (poly - exact).abs() < 1e-12 * exact.abs().max(1.0),
                    "z = {z}: polynomial {poly} vs surface {exact}"
                );
            }
        }
    }
}
