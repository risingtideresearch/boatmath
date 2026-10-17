//! End-to-end validation against closed forms: the Wigley hull, whose
//! Michell inner integrals have an analytic expression, and transom
//! detection and closure on hulls whose transoms have closed forms. Built on
//! the exact B-spline kernel and on sectional hulls cut from the same
//! splines at their own stations.

mod wigley {
    use crate::hull::FromHull;
    use crate::hull::Hull;
    use crate::hulls;
    use crate::sectional;
    use crate::{Conditions, Placement, ViscousOptions, WaveOptions};
    use hullgeom::bspline::BSplineSurface;
    use hullgeom::sectional::{DepthQuadrature, SectionalHull};

    const G: f64 = hullgeom::STANDARD_GRAVITY;

    /// The sectional hull cut from a B-spline one at its own Greville stations:
    /// exact up to the depth quadrature.
    fn sec(hull: &Hull) -> SectionalHull {
        SectionalHull::from_hull(hull, &DepthQuadrature::default()).unwrap()
    }

    /// Analytic Michell inner integrals for the Wigley hull
    /// f = (B/2)(1-(2x/L)²)(1-(z/T)²), x ∈ [-L/2, L/2], z ∈ [0, T]:
    /// I = 0 (fore-aft symmetry about the midpoint), and
    /// J(λ) = -(4B/L²) · X(k) · Z(κ), k = νλ, κ = νλ²,
    /// X = ∫ x sin(kx) dx = 2 (sin(ka)/k² - a cos(ka)/k), a = L/2,
    /// Z = ∫ (1-(z/T)²) e^{-κz} dz = (1-E)/κ - (2 - E(κ²T²+2κT+2))/(κ³T²).
    fn wigley_j(l: f64, b: f64, t: f64, nu: f64, lambda: f64) -> f64 {
        let a = l / 2.0;
        let k = nu * lambda;
        let kappa = nu * lambda * lambda;
        let x_int = 2.0 * ((k * a).sin() / (k * k) - a * (k * a).cos() / k);
        let e = (-kappa * t).exp();
        let z_int = (1.0 - e) / kappa
            - (2.0 - e * (kappa * kappa * t * t + 2.0 * kappa * t + 2.0)) / (kappa.powi(3) * t * t);
        -(4.0 * b / (l * l)) * x_int * z_int
    }

    #[test]
    fn wigley_inner_integrals_match_analytic() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = hulls::wigley(l, b, t).unwrap();
        for fn_ in [0.2, 0.3, 0.45] {
            let u = fn_ * (G * l).sqrt();
            let cond = Conditions::freshwater(u);
            let nu = G / (u * u);
            for lambda in [1.0, 1.05, 1.3, 2.0, 3.7, 8.0, 20.0] {
                let (i, j) = crate::michell::inner_integrals(&hull, &cond, lambda).unwrap();
                let want = wigley_j(l, b, t, nu, lambda);
                let scale = want.abs().max(1e-12);
                assert!(
                    i.abs() < 1e-10 * scale.max(1.0),
                    "Fn={fn_} λ={lambda}: I = {i}, expected 0"
                );
                assert!(
                    (j - want).abs() < 1e-9 * scale,
                    "Fn={fn_} λ={lambda}: J = {j}, want {want}"
                );
            }
        }
    }

    /// Reference outer integral by brute-force Simpson in θ using the analytic J,
    /// with a multihull interference factor applied to the integrand.
    #[allow(clippy::too_many_arguments)]
    fn reference_wave_resistance_factored(
        l: f64,
        b: f64,
        t: f64,
        u: f64,
        rho: f64,
        lambda_max: f64,
        n: usize,
        factor: impl Fn(f64) -> f64,
    ) -> f64 {
        let nu = G / (u * u);
        let theta_max = (1.0 / lambda_max).acos();
        let dt = theta_max / n as f64;
        let f = |theta: f64| -> f64 {
            let sec = 1.0 / theta.cos();
            let j = wigley_j(l, b, t, nu, sec);
            j * j * factor(theta) * sec * sec * sec
        };
        let mut s = f(0.0) + f(theta_max);
        for i in 1..n {
            let w = if i % 2 == 1 { 4.0 } else { 2.0 };
            s += w * f(i as f64 * dt);
        }
        let integral = s * dt / 3.0;
        4.0 * rho * G * G / (std::f64::consts::PI * u * u) * integral
    }

    fn reference_wave_resistance(l: f64, b: f64, t: f64, u: f64, rho: f64) -> f64 {
        reference_wave_resistance_factored(l, b, t, u, rho, 200.0, 2_000_000, |_| 1.0)
    }

    #[test]
    fn wigley_wave_resistance_matches_reference() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = hulls::wigley(l, b, t).unwrap();
        let opts = WaveOptions {
            rel_tol: 1e-7,
            max_refinements: 6,
            ..Default::default()
        };
        for fn_ in [0.25, 0.35] {
            let u = fn_ * (G * l).sqrt();
            let cond = Conditions::freshwater(u);
            let got = sectional::wave_resistance(&sec(&hull), &cond, &opts).unwrap();
            let want = reference_wave_resistance(l, b, t, u, cond.fluid.density);
            let rel = (got.resistance - want).abs() / want;
            assert!(
                rel < 2e-5,
                "Fn={fn_}: got {} N, reference {} N (rel {rel:.2e}, est {:.2e})",
                got.resistance,
                want,
                got.est_rel_error
            );
            assert!(got.resistance > 0.0);
        }
    }

    #[test]
    fn coincident_pair_quadruples_wave_resistance() {
        // Two coincident hulls: amplitudes add, |2F|² = 4|F|².
        let hull = sec(&hulls::wigley(10.0, 1.0, 0.625).unwrap());
        let cond = Conditions::seawater(3.0);
        let w = WaveOptions::default();
        let single = sectional::wave_resistance(&hull, &cond, &w)
            .unwrap()
            .resistance;
        let pair = [(&hull, Placement::default()), (&hull, Placement::default())];
        let both = sectional::multihull_wave_resistance(&pair, &cond, &w)
            .unwrap()
            .resistance;
        assert!(
            (both - 4.0 * single).abs() < 1e-9 * both,
            "pair {both} vs 4x single {}",
            4.0 * single
        );
    }

    #[test]
    fn catamaran_matches_analytic_interference() {
        // Two Wigley demihulls at y = ±s/2: interference factor
        // 4 cos²(½ ν s sec θ tan θ).
        let (l, b, t) = (10.0, 1.0, 0.625);
        let s = 2.0f64;
        let hull = sec(&hulls::wigley(l, b, t).unwrap());
        let fn_ = 0.35;
        let u = fn_ * (G * l).sqrt();
        let cond = Conditions::freshwater(u);
        let nu = G / (u * u);
        let members = [
            (&hull, Placement { x: 0.0, y: s / 2.0 }),
            (
                &hull,
                Placement {
                    x: 0.0,
                    y: -s / 2.0,
                },
            ),
        ];
        let opts = WaveOptions {
            rel_tol: 1e-7,
            max_refinements: 6,
            ..Default::default()
        };
        let got = sectional::multihull_wave_resistance(&members, &cond, &opts)
            .unwrap()
            .resistance;
        let want = reference_wave_resistance_factored(
            l,
            b,
            t,
            u,
            cond.fluid.density,
            100.0,
            4_000_000,
            |theta| {
                let sec = 1.0 / theta.cos();
                let alpha = 0.5 * nu * s * sec * theta.tan();
                4.0 * alpha.cos().powi(2)
            },
        );
        let rel = (got - want).abs() / want;
        assert!(
            rel < 2e-4,
            "catamaran: got {got} N, reference {want} N (rel {rel:.2e})"
        );
    }

    #[test]
    fn tandem_matches_analytic_interference() {
        // Two Wigley hulls in line, staggered by d: interference factor
        // 4 cos²(½ ν d sec θ) from the longitudinal phase.
        let (l, b, t) = (10.0, 1.0, 0.625);
        let d = 6.0f64;
        let hull = sec(&hulls::wigley(l, b, t).unwrap());
        let fn_ = 0.35;
        let u = fn_ * (G * l).sqrt();
        let cond = Conditions::freshwater(u);
        let nu = G / (u * u);
        let members = [
            (&hull, Placement { x: d / 2.0, y: 0.0 }),
            (
                &hull,
                Placement {
                    x: -d / 2.0,
                    y: 0.0,
                },
            ),
        ];
        let opts = WaveOptions {
            rel_tol: 1e-7,
            max_refinements: 6,
            ..Default::default()
        };
        let got = sectional::multihull_wave_resistance(&members, &cond, &opts)
            .unwrap()
            .resistance;
        let want = reference_wave_resistance_factored(
            l,
            b,
            t,
            u,
            cond.fluid.density,
            100.0,
            4_000_000,
            |theta| {
                let alpha = 0.5 * nu * d / theta.cos();
                4.0 * alpha.cos().powi(2)
            },
        );
        let rel = (got - want).abs() / want;
        assert!(
            rel < 2e-4,
            "tandem: got {got} N, reference {want} N (rel {rel:.2e})"
        );
    }

    #[test]
    fn multihull_breakdown_is_consistent() {
        let hull = sec(&hulls::wigley(10.0, 1.0, 0.625).unwrap());
        let cond = Conditions::seawater(3.0);
        let members = [
            (&hull, Placement { x: 0.0, y: 1.5 }),
            (&hull, Placement { x: 0.0, y: -1.5 }),
        ];
        let (w, v) = (WaveOptions::default(), ViscousOptions::default());
        let r = sectional::multihull_resistance(&members, &cond, &w, &v).unwrap();
        assert_eq!(r.viscous.len(), 2);
        let solo = sectional::wave_resistance(&hull, &cond, &w)
            .unwrap()
            .resistance;
        assert!((r.solo_wave_total - 2.0 * solo).abs() < 1e-9 * r.solo_wave_total);
        assert!((r.interference - r.wave.resistance / r.solo_wave_total).abs() < 1e-12);
        assert!((r.total - r.wave.resistance - r.viscous_total).abs() < 1e-9 * r.total);
        assert!((r.effective_power - r.total * 3.0).abs() < 1e-9 * r.effective_power);
        assert!((r.wetted_surface - 2.0 * hull.wetted_surface()).abs() < 1e-9 * r.wetted_surface);
    }

    #[test]
    fn wigley_geometry_integrals() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = hulls::wigley(l, b, t).unwrap();
        // Volume: ∇ = 2 (B/2)(2L/3)(2T/3) = 4 B L T / 9.
        let vol_exact = 4.0 * b * l * t / 9.0;
        assert!(
            (hull.displaced_volume() - vol_exact).abs() < 1e-12 * vol_exact,
            "volume {} vs {vol_exact}",
            hull.displaced_volume()
        );
        // Wetted surface: 2D Simpson reference with analytic derivatives.
        let n = 800usize;
        let (dx, dz) = (l / n as f64, t / n as f64);
        let integrand = |x: f64, z: f64| -> f64 {
            let fx = -4.0 * b / (l * l) * x * (1.0 - (z / t).powi(2));
            let fz = -(b / 2.0) * (1.0 - (2.0 * x / l).powi(2)) * 2.0 * z / (t * t);
            (1.0 + fx * fx + fz * fz).sqrt()
        };
        let mut s = 0.0;
        for i in 0..=n {
            let wx = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let x = -l / 2.0 + i as f64 * dx;
            for j in 0..=n {
                let wz = if j == 0 || j == n {
                    1.0
                } else if j % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                s += wx * wz * integrand(x, j as f64 * dz);
            }
        }
        let wetted_ref = 2.0 * s * dx * dz / 9.0;
        let rel = (hull.wetted_surface() - wetted_ref).abs() / wetted_ref;
        assert!(
            rel < 1e-8,
            "wetted {} vs {wetted_ref} (rel {rel:.2e})",
            hull.wetted_surface()
        );
        assert!((hull.length() - l).abs() < 1e-12);
        assert!((hull.draft() - t).abs() < 1e-12);
        // Hydrostatics closed forms for the Wigley hull (symmetric about x = 0):
        // LCB = LCF = 0, A_w = 2BL/3, second moment = BL³/30.
        assert!(hull.lcb_x().abs() < 1e-10, "lcb {}", hull.lcb_x());
        let aw_exact = 2.0 * b * l / 3.0;
        assert!(
            (hull.waterplane_area() - aw_exact).abs() < 1e-12 * aw_exact,
            "Aw {} vs {aw_exact}",
            hull.waterplane_area()
        );
        assert!(hull.lcf_x().abs() < 1e-10);
        let ixx_exact = b * l.powi(3) / 30.0;
        assert!(
            (hull.waterplane_second_moment() - ixx_exact).abs() < 1e-12 * ixx_exact,
            "Ixx {} vs {ixx_exact}",
            hull.waterplane_second_moment()
        );
    }

    #[test]
    fn multispan_inner_integrals_match_brute_force() {
        // Multi-span cubic × quadratic hull; compare the closed-form inner
        // integrals against direct 2D Gauss-Legendre quadrature of
        // ∬ fx e^{-κz} (cos, sin)(k(x - x_mid)) dx dz at low λ (mild oscillation,
        // so a dense per-span rule is a trustworthy reference).
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 2.5, 5.0, 7.5, 10.0, 10.0, 10.0, 10.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.6, 1.2, 1.2, 1.2];
        let (nx, nz) = (7usize, 4usize);
        let mut control = vec![0.0; nx * nz];
        for i in 0..nx {
            for j in 0..nz {
                // Positive, closing at the ends in x; smooth-ish in z.
                let gi = [0.0, 0.35, 0.8, 1.0, 0.8, 0.35, 0.0][i];
                control[i * nz + j] = gi * (1.5 - 0.31 * j as f64);
            }
        }
        let surf = BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        let hull = Hull::new(surf).unwrap();

        let u = 6.0;
        let cond = Conditions::seawater(u);
        let nu = G / (u * u);
        let x_mid = 5.0;

        // Dense Gauss-Legendre per span (nodes/weights on [-1,1], order 30).
        let order = 30usize;
        let (gn, gw) = gauss_legendre_ref(order);
        let xs = [0.0, 2.5, 5.0, 7.5, 10.0];
        let zs = [0.0, 0.6, 1.2];

        for lambda in [1.0, 1.2, 1.6] {
            let k = nu * lambda;
            let kappa = nu * lambda * lambda;
            let (mut i_ref, mut j_ref) = (0.0f64, 0.0f64);
            for w in xs.windows(2) {
                let (x0, x1) = (w[0], w[1]);
                for v in zs.windows(2) {
                    let (z0, z1) = (v[0], v[1]);
                    let jac = (x1 - x0) / 2.0 * ((z1 - z0) / 2.0);
                    for (a, &na) in gn.iter().enumerate() {
                        let x = x0 + (x1 - x0) * (na + 1.0) / 2.0;
                        for (b, &nb) in gn.iter().enumerate() {
                            let z = z0 + (z1 - z0) * (nb + 1.0) / 2.0;
                            let fx = hull.surface().eval_deriv(x, z, 1, 0);
                            let common = gw[a] * gw[b] * jac * fx * (-kappa * z).exp();
                            i_ref += common * (k * (x - x_mid)).cos();
                            j_ref += common * (k * (x - x_mid)).sin();
                        }
                    }
                }
            }
            let (i_got, j_got) = crate::michell::inner_integrals(&hull, &cond, lambda).unwrap();
            let scale = (i_ref * i_ref + j_ref * j_ref).sqrt().max(1e-12);
            assert!(
                (i_got - i_ref).abs() < 1e-10 * scale && (j_got - j_ref).abs() < 1e-10 * scale,
                "λ={lambda}: got ({i_got}, {j_got}), want ({i_ref}, {j_ref})"
            );
        }
    }

    #[test]
    fn deep_span_cutoff_is_invisible_at_large_lambda() {
        // At large λ the vertical decay κ = νλ² confines the integrand to a
        // sliver under the waterline, and the kernel stops walking z-spans whose
        // e^{−κz₀} has fallen below 1e-20. That is the optimisation that pays for
        // a finely-resolved loft, so pin it against a brute-force reference at
        // λ well past the point where the deep spans drop out.
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 2.5, 5.0, 7.5, 10.0, 10.0, 10.0, 10.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.3, 0.6, 0.9, 1.2, 1.2, 1.2];
        let (nx, nz) = (7usize, 6usize);
        let mut control = vec![0.0; nx * nz];
        for i in 0..nx {
            for j in 0..nz {
                let gi = [0.0, 0.35, 0.8, 1.0, 0.8, 0.35, 0.0][i];
                control[i * nz + j] = gi * (1.5 - 0.22 * j as f64);
            }
        }
        let surf = BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        let hull = Hull::new(surf).unwrap();

        let u = 6.0;
        let cond = Conditions::seawater(u);
        let nu = G / (u * u);
        let x_mid = 5.0;
        let (gn, gw) = gauss_legendre_ref(24);
        let xs = [0.0, 2.5, 5.0, 7.5, 10.0];
        let zs = [0.0, 0.3, 0.6, 0.9, 1.2];

        for lambda in [8.0f64, 20.0, 40.0] {
            let k = nu * lambda;
            let kappa = nu * lambda * lambda;
            // Panel the reference finely enough for the oscillation at this λ,
            // and for the decay: both scale with λ.
            let np_x = ((k * 2.5 / 0.5).ceil() as usize).clamp(4, 400);
            let np_z = ((kappa * 0.3 / 0.5).ceil() as usize).clamp(4, 800);
            let (mut i_ref, mut j_ref) = (0.0f64, 0.0f64);
            for w in xs.windows(2) {
                for px in 0..np_x {
                    let x0 = w[0] + (w[1] - w[0]) * px as f64 / np_x as f64;
                    let x1 = w[0] + (w[1] - w[0]) * (px + 1) as f64 / np_x as f64;
                    for v in zs.windows(2) {
                        for pz in 0..np_z {
                            let z0 = v[0] + (v[1] - v[0]) * pz as f64 / np_z as f64;
                            let z1 = v[0] + (v[1] - v[0]) * (pz + 1) as f64 / np_z as f64;
                            // Once the decay is negligible the rest contributes
                            // nothing to the reference either.
                            if (-kappa * z0).exp() < 1e-25 {
                                break;
                            }
                            let jac = (x1 - x0) / 2.0 * ((z1 - z0) / 2.0);
                            for (a, &na) in gn.iter().enumerate() {
                                let x = x0 + (x1 - x0) * (na + 1.0) / 2.0;
                                for (b, &nb) in gn.iter().enumerate() {
                                    let z = z0 + (z1 - z0) * (nb + 1.0) / 2.0;
                                    let fx = hull.surface().eval_deriv(x, z, 1, 0);
                                    let common = gw[a] * gw[b] * jac * fx * (-kappa * z).exp();
                                    i_ref += common * (k * (x - x_mid)).cos();
                                    j_ref += common * (k * (x - x_mid)).sin();
                                }
                            }
                        }
                    }
                }
            }
            let (i_got, j_got) = crate::michell::inner_integrals(&hull, &cond, lambda).unwrap();
            let scale = (i_ref * i_ref + j_ref * j_ref).sqrt().max(1e-300);
            let err = ((i_got - i_ref).hypot(j_got - j_ref)) / scale;
            assert!(
                err < 1e-8,
                "λ={lambda}: got ({i_got:e}, {j_got:e}), want ({i_ref:e}, {j_ref:e}), rel {err:e}"
            );
        }
    }

    /// Local Gauss-Legendre reference (independent of the crate's internal one).
    fn gauss_legendre_ref(n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut nodes = vec![0.0; n];
        let mut weights = vec![0.0; n];
        for i in 0..n {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut dp = 0.0;
            for _ in 0..200 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    let kf = k as f64;
                    let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                    p0 = p1;
                    p1 = p2;
                }
                dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
                let dx = p1 / dp;
                x -= dx;
                if dx.abs() < 1e-15 {
                    break;
                }
            }
            nodes[i] = -x;
            weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
        }
        (nodes, weights)
    }

    #[test]
    fn wedge_prism_geometry() {
        // Piecewise-linear tent in x, constant in z: a wall-sided "wedge" prism.
        // Volume = B L T / 2; wetted = 2 L T sqrt(1 + (B/L)²).
        let (l, b, t) = (8.0, 1.2, 0.9);
        let knots_x = vec![0.0, 0.0, l / 2.0, l, l];
        let knots_z = vec![0.0, 0.0, t, t];
        let control = vec![0.0, 0.0, b / 2.0, b / 2.0, 0.0, 0.0];
        let hull =
            Hull::new(BSplineSurface::new(1, 1, knots_x, knots_z, control).unwrap()).unwrap();
        let vol = b * l * t / 2.0;
        let wet = 2.0 * l * t * (1.0 + (b / l) * (b / l)).sqrt();
        assert!((hull.displaced_volume() - vol).abs() < 1e-12 * vol);
        assert!((hull.wetted_surface() - wet).abs() < 1e-10 * wet);
    }

    #[test]
    fn froude_sweep_is_finite_positive_and_humped() {
        let hull = sec(&hulls::wigley(10.0, 1.0, 0.625).unwrap());
        let members = [(&hull, Placement::default())];
        let (w, v) = (WaveOptions::default(), ViscousOptions::default());
        let mut cws = Vec::new();
        for i in 0..=20 {
            let fn_ = 0.15 + 0.35 * i as f64 / 20.0;
            let u = fn_ * (G * 10.0f64).sqrt();
            let r = sectional::multihull_resistance(&members, &Conditions::seawater(u), &w, &v)
                .unwrap();
            assert!(
                r.wave.resistance.is_finite() && r.wave.resistance >= 0.0,
                "Fn={fn_}: Rw = {}",
                r.wave.resistance
            );
            assert!(r.viscous_total > 0.0 && r.total > r.wave.resistance);
            cws.push(r.cw);
        }
        // The Cw(Fn) curve must not be monotone (humps and hollows).
        let increasing = cws.windows(2).all(|w| w[1] >= w[0]);
        let decreasing = cws.windows(2).all(|w| w[1] <= w[0]);
        assert!(!increasing && !decreasing, "Cw curve unexpectedly monotone");
    }

    #[test]
    fn rejects_bad_hulls_and_conditions() {
        // Negative half-beam.
        let s = BSplineSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.5, 0.5, -0.5, 0.5],
        )
        .unwrap();
        assert!(Hull::new(s).is_err());
        // z-domain not starting at the waterline.
        let s = BSplineSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.5, 0.5, 1.0, 1.0],
            vec![0.5; 4],
        )
        .unwrap();
        assert!(Hull::new(s).is_err());
        // Nonsense speed.
        let hull = sec(&hulls::wigley(10.0, 1.0, 0.625).unwrap());
        let w = WaveOptions::default();
        assert!(sectional::wave_resistance(&hull, &Conditions::seawater(0.0), &w).is_err());
        assert!(sectional::wave_resistance(&hull, &Conditions::seawater(f64::NAN), &w).is_err());
    }
}

/// Transom detection (does the half-breadth close at the aft end, and if
/// not, how much of a transom is it?) and the virtual-appendage closure.
mod transom {
    use crate::hull::FromHull;
    use crate::hull::Hull;
    use crate::hulls;
    use crate::sectional;
    use crate::SectionalWave;
    use crate::{Conditions, Placement, TransomClosure, WaveOptions};
    use hullgeom::bspline::BSplineSurface;
    use hullgeom::sectional::{DepthQuadrature, SectionalHull};

    /// A wedge closing linearly toward the bow: `f = A (1 − x/L)(1 − z/T)`, so the
    /// aft section is the full `f_T(z) = A(1 − z/T)` and every quantity the
    /// detector reports has a closed form.
    fn wedge(a: f64, length: f64, draft: f64) -> Hull {
        let knots_x = vec![0.0, 0.0, length, length];
        let knots_z = vec![0.0, 0.0, draft, draft];
        // Row-major, z fastest: (x0,z0) (x0,T) (x1,z0) (x1,T).
        let control = vec![a, 0.0, 0.0, 0.0];
        Hull::new(BSplineSurface::new(1, 1, knots_x, knots_z, control).unwrap()).unwrap()
    }

    #[test]
    fn wigley_closes_at_both_ends() {
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        assert!(
            hull.transom().is_none(),
            "the Wigley half-breadth vanishes at x = ±L/2; it has no transom"
        );
    }

    #[test]
    fn wigley_max_section_area_matches_the_closed_form() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = hulls::wigley(l, b, t).unwrap();
        // A(x) = 2 ∫₀^T (B/2)(1 − (2x/L)²)(1 − (z/T)²) dz, maximal at x = 0:
        // A_X = B · (2T/3).
        let expect = 2.0 * b * t / 3.0;
        let got = hull.max_section_area();
        assert!(
            (got - expect).abs() < 1e-9 * expect,
            "A_X = {got}, expected {expect}"
        );
    }

    #[test]
    fn wedge_transom_matches_the_closed_form() {
        let (a, l, t) = (0.4, 8.0, 0.25);
        let hull = wedge(a, l, t);
        let tr = hull.transom().expect("the wedge does not close aft");

        assert!((tr.x - 0.0).abs() < 1e-12, "transom station {}", tr.x);
        assert!(
            (tr.half_beam - a).abs() < 1e-12,
            "f_T(0) = {}, expected {a}",
            tr.half_beam
        );
        // A_T = 2 ∫₀^T A(1 − z/T) dz = A·T.
        let area = a * t;
        assert!(
            (tr.area - area).abs() < 1e-12 * area,
            "A_T = {}, expected {area}",
            tr.area
        );
        // A(x) = A·T·(1 − x/L) is maximal at the transom itself, so A_T/A_X = 1.
        assert!(
            (tr.area / hull.max_section_area() - 1.0).abs() < 1e-9,
            "A_T/A_X = {}",
            tr.area / hull.max_section_area()
        );
        // Equivalent rectangle of a linearly tapering transom: A_T/(2·A) = T/2.
        let depth = t / 2.0;
        assert!(
            (tr.depth - depth).abs() < 1e-12 * t,
            "depth = {}, expected {depth}",
            tr.depth
        );
    }

    #[test]
    fn a_numerically_closed_stern_is_not_a_transom() {
        // A hull with a full midsection and a stern half-beam 1e-5 of it: the loft
        // wiggle case, which must not read as a transom.
        let (l, t) = (8.0, 0.25);
        let knots_x = vec![0.0, 0.0, 0.0, l, l, l];
        let knots_z = vec![0.0, 0.0, t, t];
        let control = vec![1e-5, 0.0, 1.0, 0.0, 0.0, 0.0];
        let hull =
            Hull::new(BSplineSurface::new(2, 1, knots_x, knots_z, control).unwrap()).unwrap();
        let tr = hull.transom();
        assert!(
            tr.is_none(),
            "stern area ratio {:?} should read as closed",
            tr.map(|t| t.area / hull.max_section_area())
        );
    }

    // ---------------------------------------------------------------------------
    // The closure
    // ---------------------------------------------------------------------------

    /// The closure on the sectional kernel, cut from the B-spline hull at its
    /// own stations (the harness: exact up to the depth quadrature).
    fn sec(hull: &Hull) -> SectionalHull {
        SectionalHull::from_hull(hull, &DepthQuadrature::default()).unwrap()
    }

    fn wave_resistance_with(
        hull: &SectionalHull,
        cond: &Conditions,
        opts: &WaveOptions,
    ) -> crate::Result<crate::WaveResistance> {
        sectional::wave_resistance(hull, cond, opts)
    }

    fn opts(transom: TransomClosure) -> WaveOptions {
        WaveOptions {
            transom,
            ..Default::default()
        }
    }

    #[test]
    fn closure_is_inert_on_a_hull_that_closes_aft() {
        // Wigley has no transom, so every closure setting must give bit-for-bit
        // the same resistance — the guarantee that this feature cannot perturb
        // classical Michell.
        let hull = sec(&hulls::wigley(10.0, 1.0, 0.625).unwrap());
        let cond = Conditions::seawater(3.0);
        let base = wave_resistance_with(&hull, &cond, &opts(TransomClosure::None))
            .unwrap()
            .resistance;
        for t in [
            TransomClosure::default(),
            TransomClosure::Ballistic { coeff: 4.0 },
            TransomClosure::Fixed { length: 2.0 },
        ] {
            let got = wave_resistance_with(&hull, &cond, &opts(t))
                .unwrap()
                .resistance;
            assert_eq!(
                base.to_bits(),
                got.to_bits(),
                "{t:?} perturbed a closed hull"
            );
        }
    }

    #[test]
    fn zero_hollow_reproduces_the_analytic_step_term() {
        // As L_v -> 0 the virtual appendage collapses to a delta sheet at the
        // transom, whose amplitude is the closed form
        //     F_step = e^{i nu lambda x_T} * integral f_T(z) e^{-nu lambda^2 z} dz.
        // Here f_T(z) = A(1 - z/T), so the z-integral is elementary.
        let (a, l, t) = (0.4, 8.0, 0.25);
        let hull = sec(&wedge(a, l, t));
        let cond = Conditions::seawater(3.0);
        let nu = cond.gravity / (cond.speed * cond.speed);

        for lambda in [1.0, 1.7, 4.0] {
            let kappa = nu * lambda * lambda;
            // integral_0^T A(1 - z/T) e^{-kz} dz = A[(1-e^{-kT})/k - (1 - (1+kT)e^{-kT})/(k^2 T)]
            let e = (-kappa * t).exp();
            let zint =
                a * ((1.0 - e) / kappa - (1.0 - (1.0 + kappa * t) * e) / (kappa * kappa * t));
            // Phase is measured from the hull's x-centre, as the kernel does.
            let phase = nu * lambda * (0.0 - l / 2.0);
            let (want_re, want_im) = (zint * phase.cos(), zint * phase.sin());

            let closed = amplitude(&hull, nu, lambda, TransomClosure::Fixed { length: 0.0 });
            let open = amplitude(&hull, nu, lambda, TransomClosure::None);
            let (dre, dim) = (closed.0 - open.0, closed.1 - open.1);
            let scale = zint.abs().max(1e-12);
            assert!(
                (dre - want_re).abs() < 1e-9 * scale && (dim - want_im).abs() < 1e-9 * scale,
                "lambda {lambda}: closure added ({dre}, {dim}), analytic step ({want_re}, {want_im})"
            );
        }
    }

    /// The free-wave amplitude `I + iJ` with an explicit closure.
    fn amplitude(
        hull: &SectionalHull,
        nu: f64,
        lambda: f64,
        transom: TransomClosure,
    ) -> (f64, f64) {
        let f = hull.amplitude_closed(nu, lambda, transom, &mut Default::default());
        (f.re, f.im)
    }

    #[test]
    fn a_longer_hollow_radiates_less() {
        // The physical content of the closure: stretching the hollow spreads the
        // same net change in half-beam over more length, so it makes less wave.
        // Monotone in the hollow length, and bracketed by the two limits.
        let hull = sec(&wedge(0.4, 8.0, 0.25));
        let cond = Conditions::seawater(3.0);
        let r = |t| {
            wave_resistance_with(&hull, &cond, &opts(t))
                .unwrap()
                .resistance
        };
        let step = r(TransomClosure::Fixed { length: 0.0 });
        let mut prev = step;
        for len in [0.5, 1.0, 2.0, 4.0] {
            let now = r(TransomClosure::Fixed { length: len });
            assert!(
                now < prev,
                "hollow {len} m gave {now} N, not less than {prev} N"
            );
            prev = now;
        }
        let open = r(TransomClosure::None);
        assert!(
            step > open,
            "an abrupt transom ({step} N) should make more wave than ignoring it ({open} N)"
        );
    }

    #[test]
    fn ballistic_hollow_grows_with_speed() {
        // L_v = c*U*sqrt(d_T/g): the faster the boat, the longer the hollow, so
        // the transom's share of the wave-making falls away with speed.
        let hull = sec(&wedge(0.4, 8.0, 0.25));
        let tr = hull.transom().unwrap();
        let g = 9.80665;
        for u in [2.0, 4.0, 8.0] {
            let cond = Conditions::seawater(u);
            let want = std::f64::consts::SQRT_2 * u * (tr.depth / g).sqrt();
            // Match the ballistic closure against an explicit Fixed hollow of the
            // length the formula predicts.
            let ballistic = TransomClosure::Ballistic {
                coeff: crate::BALLISTIC_COEFF,
            };
            let a = wave_resistance_with(&hull, &cond, &opts(ballistic))
                .unwrap()
                .resistance;
            let b =
                wave_resistance_with(&hull, &cond, &opts(TransomClosure::Fixed { length: want }))
                    .unwrap()
                    .resistance;
            assert!(
                (a - b).abs() <= 1e-9 * a.abs().max(1.0),
                "U = {u}: ballistic {a} N vs explicit L_v = {want} m giving {b} N"
            );
        }
    }

    #[test]
    fn auto_hollow_is_ballistic_when_wet_and_wave_scaled_when_dry() {
        let auto = TransomClosure::default();
        let ballistic = TransomClosure::Ballistic {
            coeff: crate::BALLISTIC_COEFF,
        };
        let g: f64 = 9.80665;
        let d: f64 = 0.2;
        let at = |fn_t: f64| {
            let u = fn_t * (g * d).sqrt();
            (g / (u * u), u)
        };
        // Wet: the ballistic length exactly.
        let (nu, _) = at(2.0);
        assert_eq!(auto.hollow_length(d, nu), ballistic.hollow_length(d, nu));
        // Dry: half a transverse wavelength.
        let (nu, u) = at(4.0);
        let want = std::f64::consts::PI * u * u / g;
        let got = auto.hollow_length(d, nu).unwrap();
        assert!((got - want).abs() <= 1e-12 * want, "{got} vs {want}");
        // Between: continuous and growing with speed.
        let mut last = 0.0;
        for i in 0..=40 {
            let (nu, _) = at(2.0 + 0.05 * i as f64);
            let l = auto.hollow_length(d, nu).unwrap();
            assert!(l > last, "Fn_T {}: {l} after {last}", 2.0 + 0.05 * i as f64);
            last = l;
        }
    }

    #[test]
    fn closure_survives_placement_in_a_fleet() {
        // The transom term carries its own phase, so a fleet of two transom hulls
        // must still be translation-covariant: shifting the whole fleet in x
        // changes no resistance.
        let hull = sec(&wedge(0.4, 8.0, 0.25));
        let cond = Conditions::seawater(3.0);
        let o = opts(TransomClosure::default());
        let at = |dx: f64| {
            let m = [
                (&hull, Placement { x: dx, y: 1.6 }),
                (&hull, Placement { x: dx, y: -1.6 }),
            ];
            sectional::multihull_wave_resistance(&m, &cond, &o)
                .unwrap()
                .resistance
        };
        let (a, b) = (at(0.0), at(25.0));
        assert!(
            (a - b).abs() < 1e-9 * a,
            "fleet resistance moved with a rigid x-shift: {a} vs {b}"
        );
    }
}
