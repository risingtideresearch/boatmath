//! The best Wageningen B-series propeller for a thrust at a speed: ported
//! from propopt's `web/propcore.js`, which ports its `src/propopt.py` and
//! `src/drivetrain.py`.
//!
//! At each shaft speed, and for each blade count, a coarse grid over
//! diameter and blade-area ratio finds the basin and a pattern search
//! polishes it. Pitch isn't searched: the thrust requirement pins it. A
//! sweep over shaft speed then gives the cheapest propeller at each rpm,
//! and the cheapest of those is the answer. Cavitation (Keller, Burrill),
//! a Reynolds correction to the series' Re = 2e6 fit, and an optional
//! second operating point the same propeller must reach are all as in the
//! original; see its comments for the reasoning.
//!
//! Units SI throughout: m, m/s, N, W, rpm.

use serde::{Deserialize, Serialize};

// Wageningen B-series polynomial coefficients (the Re = 2e6 fit). Each row
// is [i, j, k, l, c] meaning c · J^i · (P/D)^j · EAR^k · Z^l.
const KT_TERMS: [[f64; 5]; 39] = [
    [0., 0., 0., 0., 0.00880496],
    [0., 0., 0., 1., 0.0144043],
    [0., 0., 0., 2., -0.000606848],
    [0., 0., 1., 1., -0.0125894],
    [0., 0., 1., 2., 0.000690904],
    [0., 0., 2., 0., -0.0507214],
    [0., 1., 0., 0., 0.166351],
    [0., 1., 0., 1., 0.0143481],
    [0., 2., 0., 0., 0.158114],
    [0., 2., 1., 0., 0.415437],
    [0., 2., 2., 1., -0.00410798],
    [0., 3., 0., 0., -0.133698],
    [0., 3., 0., 1., -0.00841728],
    [0., 3., 1., 1., -0.0317791],
    [0., 3., 1., 2., 0.00421749],
    [0., 3., 2., 2., -0.00146564],
    [0., 6., 0., 0., 0.00638407],
    [1., 0., 0., 0., -0.204554],
    [1., 0., 0., 2., -0.0049819],
    [1., 0., 1., 1., 0.0109689],
    [1., 0., 2., 1., 0.018604],
    [1., 1., 0., 1., 0.0606826],
    [1., 1., 1., 0., -0.481497],
    [1., 2., 0., 2., -0.00163652],
    [1., 3., 0., 1., 0.0168424],
    [1., 6., 0., 2., -0.000328787],
    [1., 6., 2., 0., 0.010465],
    [2., 0., 0., 1., -0.0530054],
    [2., 0., 0., 2., 0.0025983],
    [2., 0., 1., 0., -0.147581],
    [2., 0., 2., 0., 0.0854559],
    [2., 6., 0., 0., -0.00132718],
    [2., 6., 0., 2., 0.000116502],
    [2., 6., 2., 0., -0.00648272],
    [3., 0., 0., 2., -0.000560528],
    [3., 0., 1., 0., 0.168496],
    [3., 0., 2., 0., -0.0504475],
    [3., 3., 0., 1., -0.00102296],
    [3., 6., 1., 2., 5.65229e-05],
];

const KQ_TERMS: [[f64; 5]; 47] = [
    [0., 0., 0., 0., 0.00379368],
    [0., 0., 2., 0., 0.015896],
    [0., 0., 2., 2., -0.0001843],
    [0., 1., 0., 1., 0.00513696],
    [0., 1., 1., 0., -0.0408811],
    [0., 1., 2., 0., -0.0502782],
    [0., 2., 0., 0., 0.00344778],
    [0., 2., 1., 0., 0.188561],
    [0., 2., 1., 1., -0.0269403],
    [0., 2., 1., 2., 0.00155334],
    [0., 2., 2., 1., 0.0126803],
    [0., 3., 1., 0., 0.0161886],
    [0., 3., 2., 0., -0.0397722],
    [0., 3., 2., 2., -0.000425399],
    [0., 6., 0., 1., -0.000313912],
    [0., 6., 1., 1., -0.00142121],
    [0., 6., 1., 2., 0.000302683],
    [0., 6., 2., 0., -0.00350024],
    [0., 6., 2., 1., 0.00334268],
    [0., 6., 2., 2., -0.0004659],
    [1., 0., 0., 1., -0.00370871],
    [1., 0., 1., 2., 0.000269551],
    [1., 0., 2., 0., 0.0471729],
    [1., 0., 2., 1., -0.00383637],
    [1., 1., 0., 0., -0.032241],
    [1., 1., 0., 1., 0.0209449],
    [1., 1., 0., 2., -0.00183491],
    [1., 1., 1., 0., -0.108009],
    [1., 1., 1., 1., 0.00438388],
    [1., 3., 1., 0., 0.00318086],
    [1., 6., 2., 2., 5.54194e-05],
    [2., 0., 0., 0., 0.00886523],
    [2., 0., 1., 1., -0.00723408],
    [2., 0., 1., 2., 0.00083265],
    [2., 1., 0., 1., 0.00474319],
    [2., 1., 1., 0., -0.0885381],
    [2., 2., 2., 0., 0.0417122],
    [2., 3., 2., 1., -0.00318278],
    [3., 0., 0., 1., -0.0106854],
    [3., 0., 1., 0., 0.0558082],
    [3., 0., 1., 1., 0.0035985],
    [3., 0., 2., 0., 0.0196283],
    [3., 1., 2., 0., -0.030055],
    [3., 2., 0., 2., 0.000112451],
    [3., 3., 0., 1., 0.00110903],
    [3., 3., 2., 2., 8.69243e-05],
    [3., 6., 0., 2., -2.97228e-05],
];

/// The fit's validity box.
pub const J_BOUNDS: [f64; 2] = [0.001, 2.0];
pub const PD_BOUNDS: [f64; 2] = [0.6, 1.4];
pub const EAR_BOUNDS: [f64; 2] = [0.3, 1.05];

pub const RHO: f64 = 1000.0; // kg/m³, fresh water
pub const NU: f64 = 1.0e-6; // m²/s, at 20 °C
const RE_FIT: f64 = 2e6;
const RE_CD_FLOOR: f64 = 5e4;
const G: f64 = 9.81;
const P_ATM: f64 = 101325.0;
const P_VAPOUR: f64 = 2339.0;
const BURRILL_A: f64 = 0.26;
const BURRILL_B: f64 = 0.57;
const BURRILL_SIGMA_RANGE: [f64; 2] = [0.1, 2.0];
/// The margin another blade count must beat to take over the curve.
pub const HYSTERESIS: f64 = 0.002;

/// The series' actual blade-area-ratio coverage per blade count.
pub fn ear_range(z: f64) -> Option<[f64; 2]> {
    Some(match z as i64 {
        2 => [0.30, 0.38],
        3 => [0.35, 0.80],
        4 => [0.40, 1.00],
        5 => [0.45, 1.05],
        6 => [0.50, 0.80],
        7 => [0.55, 0.85],
        _ => return None,
    })
}

// ------------------------------------------------------------- polynomials

/// K_T and K_Q at fixed (J, EAR, Z), as degree-6 polynomials in P/D.
#[derive(Clone, Copy, Default)]
struct Pitch {
    a: [f64; 7],
    b: [f64; 7],
}

fn pitch_polynomials(j: f64, ear: f64, z: f64) -> Pitch {
    let mut p = Pitch::default();
    let jp = [1.0, j, j * j, j * j * j];
    let bp = [1.0, ear, ear * ear];
    let zp = [1.0, z, z * z];
    for r in &KT_TERMS {
        p.a[r[1] as usize] += r[4] * jp[r[0] as usize] * bp[r[2] as usize] * zp[r[3] as usize];
    }
    for r in &KQ_TERMS {
        p.b[r[1] as usize] += r[4] * jp[r[0] as usize] * bp[r[2] as usize] * zp[r[3] as usize];
    }
    p
}

fn polyval6(c: &[f64; 7], x: f64) -> f64 {
    (((((c[6] * x + c[5]) * x + c[4]) * x + c[3]) * x + c[2]) * x + c[1]) * x + c[0]
}

pub fn k_t(j: f64, pd: f64, ear: f64, z: f64) -> f64 {
    polyval6(&pitch_polynomials(j, ear, z).a, pd)
}

pub fn k_q(j: f64, pd: f64, ear: f64, z: f64) -> f64 {
    polyval6(&pitch_polynomials(j, ear, z).b, pd)
}

// ------------------------------------------------------------- Reynolds

fn chord075(d: f64, ear: f64, z: f64) -> f64 {
    2.073 * (ear / z) * d
}

fn reynolds075(d: f64, ear: f64, z: f64, v_a: f64, n: f64, shaft_angle_deg: f64) -> f64 {
    let eps = shaft_angle_deg.to_radians();
    let (sin_e, cos_e) = eps.sin_cos();
    let v_r2 = (v_a * cos_e).powi(2)
        + (0.75 * std::f64::consts::PI * n * d).powi(2)
        + 0.5 * (v_a * sin_e).powi(2);
    chord075(d, ear, z) * v_r2.sqrt() / NU
}

fn thickness_ratio075(ear: f64, z: f64) -> f64 {
    (0.0185 - 0.00125 * z) / (2.073 * ear / z)
}

fn section_drag_coefficient(re: f64, t_over_c: f64) -> f64 {
    let r = re.max(1.0);
    2.0 * (1.0 + 2.0 * t_over_c) * (0.044 / r.powf(1.0 / 6.0) - 5.0 / r.powf(2.0 / 3.0))
}

struct Reynolds {
    re: f64,
    d_cd: f64,
    dkt_per_pd: f64,
    dkq: f64,
    cd_formula_ok: bool,
}

fn reynolds_correction(
    d: f64,
    ear: f64,
    z: f64,
    v_a: f64,
    n: f64,
    shaft_angle_deg: f64,
) -> Reynolds {
    let re = reynolds075(d, ear, z, v_a, n, shaft_angle_deg);
    let tc = thickness_ratio075(ear, z);
    let d_cd = section_drag_coefficient(re, tc) - section_drag_coefficient(RE_FIT, tc);
    let cz_d = 2.073 * ear;
    Reynolds {
        re,
        d_cd,
        dkt_per_pd: -0.30 * d_cd * cz_d,
        dkq: 0.25 * d_cd * cz_d,
        cd_formula_ok: re >= RE_CD_FLOOR,
    }
}

// ------------------------------------------------------------- cavitation

fn static_pressure(depth: f64) -> f64 {
    P_ATM + RHO * G * depth
}

fn keller_min_ear(t: f64, d: f64, z: f64, depth: f64, k: f64) -> f64 {
    let dp = (static_pressure(depth) - P_VAPOUR).max(1e-6);
    (1.3 + 0.3 * z) * t / (dp * d * d) + k
}

fn projected_area_ratio(ear: f64, pd: f64) -> f64 {
    ear * (1.067 - 0.229 * pd)
}

fn vr2_07(v_a: f64, n: f64, d: f64, shaft_angle_deg: f64) -> f64 {
    let eps = shaft_angle_deg.to_radians();
    let (sin_e, cos_e) = eps.sin_cos();
    (v_a * cos_e).powi(2) + (0.7 * std::f64::consts::PI * n * d + v_a * sin_e.abs()).powi(2)
}

fn cavitation_number07(v_a: f64, n: f64, d: f64, depth: f64, shaft_angle_deg: f64) -> f64 {
    (static_pressure(depth) - P_VAPOUR)
        / (0.5 * RHO * vr2_07(v_a, n, d, shaft_angle_deg).max(1e-12))
}

fn burrill_tau_c(t: f64, ear: f64, pd: f64, d: f64, v_a: f64, n: f64, shaft_angle_deg: f64) -> f64 {
    let a_p = projected_area_ratio(ear, pd) * (std::f64::consts::PI * d * d / 4.0);
    t / (0.5 * RHO * a_p.max(1e-12) * vr2_07(v_a, n, d, shaft_angle_deg).max(1e-12))
}

fn burrill_tau_allow(sigma: f64) -> f64 {
    BURRILL_A * sigma.max(1e-12).powf(BURRILL_B)
}

fn ear_in_series_range(ear: f64, z: f64) -> bool {
    ear_range(z).is_some_and(|r| ear >= r[0] && ear <= r[1])
}

// ------------------------------------------------------------- inputs

/// What the boat asks of its propellers. `thrust` is the total across all
/// shafts; everything downstream is one shaft's.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Inputs {
    /// Ship speed [m/s].
    pub speed: f64,
    /// Total thrust [N].
    pub thrust: f64,
    pub shafts: u32,
    /// Wake fraction w: the propeller sees V_A = V_s (1 − w).
    pub wake: f64,
    /// Thrust deduction t: only T (1 − t) shows up as tow.
    pub thrust_deduction: f64,
    pub d_min: f64,
    pub d_max: f64,
    pub blades: Vec<u32>,
    /// Shaft immersion [m], for cavitation.
    pub depth: f64,
    /// Keller's margin k.
    pub keller_k: f64,
    pub cavitation: bool,
    pub re_correct: bool,
    pub strict_ear: bool,
    /// A second operating point the same propeller must reach: speed [m/s]
    /// and total thrust [N].
    pub top: Option<TopInputs>,
    /// Shaft inclination angle [deg] (shaft tilted down toward stern / bow up).
    #[serde(default)]
    pub shaft_angle_deg: f64,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub struct TopInputs {
    pub speed: f64,
    pub thrust: f64,
}

/// [`Inputs`] worked out to one shaft.
#[derive(Clone, Debug)]
pub struct Config {
    pub v_s: f64,
    pub v_a: f64,
    pub r_t: f64,
    pub shafts: u32,
    pub t_total: f64,
    pub t_target: f64,
    pub p_e: f64,
    pub eta_h: f64,
    pub d_min: f64,
    pub d_max: f64,
    pub z: Vec<f64>,
    pub depth: f64,
    pub keller_k: f64,
    pub cavitation: bool,
    pub strict_ear: bool,
    pub re_correct: bool,
    pub grid_d: usize,
    pub grid_ear: usize,
    pub top: Option<TopConfig>,
    pub shaft_angle_deg: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TopConfig {
    pub v_s: f64,
    pub v_a: f64,
    pub t_total: f64,
    pub t_target: f64,
}

impl Config {
    pub fn new(i: &Inputs) -> Config {
        let (w, t) = (i.wake, i.thrust_deduction);
        let shafts = i.shafts.max(1);
        let r_t = i.thrust * (1.0 - t);
        Config {
            v_s: i.speed,
            v_a: i.speed * (1.0 - w),
            r_t,
            shafts,
            t_total: i.thrust,
            t_target: i.thrust / shafts as f64,
            p_e: r_t * i.speed,
            eta_h: (1.0 - t) / (1.0 - w),
            d_min: i.d_min.min(i.d_max),
            d_max: i.d_max,
            z: i.blades.iter().map(|&z| z as f64).collect(),
            depth: i.depth,
            keller_k: i.keller_k,
            cavitation: i.cavitation,
            strict_ear: i.strict_ear,
            re_correct: i.re_correct,
            grid_d: 16,
            grid_ear: 12,
            top: i
                .top
                .filter(|t| t.speed > 0.0 && t.thrust > 0.0)
                .map(|t| TopConfig {
                    v_s: t.speed,
                    v_a: t.speed * (1.0 - w),
                    t_total: t.thrust,
                    t_target: t.thrust / shafts as f64,
                }),
            shaft_angle_deg: i.shaft_angle_deg,
        }
    }
}

// ------------------------------------------------------------- answers

/// The same propeller at the second operating point.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct TopPoint {
    pub n: f64,
    pub rpm: f64,
    pub J: f64,
    pub K_T: f64,
    pub K_Q: f64,
    pub T: f64,
    pub Q: f64,
    pub P_shaft: f64,
    pub eta0: f64,
    pub sigma: f64,
    pub tau_c: f64,
    pub tau_allow: f64,
    pub burrillOk: bool,
    pub EAR_min: f64,
    pub kellerOk: bool,
}

/// One blade count's best at a shaft speed, in brief.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct PerZ {
    pub Z: f64,
    pub D: f64,
    pub PD: f64,
    pub EAR: f64,
    pub eta0: f64,
    pub P_shaft: f64,
}

/// One geometry at one shaft speed, its pitch sized to hold the thrust.
/// Field names are propcore.js's.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct Design {
    pub n: f64,
    pub rpm: f64,
    pub D: f64,
    pub PD: f64,
    pub EAR: f64,
    pub Z: f64,
    pub J: f64,
    pub K_T: f64,
    pub K_Q: f64,
    pub T: f64,
    pub Q: f64,
    pub P_shaft: f64,
    pub eta0: f64,
    pub Cth: Option<f64>,
    pub etaIdeal: Option<f64>,
    pub etaOfIdeal: Option<f64>,
    pub top: Option<TopPoint>,
    pub topOk: bool,
    pub Re: f64,
    pub dC_D: f64,
    pub cdFormulaOk: bool,
    pub EAR_min: f64,
    pub kellerOk: bool,
    pub sigma: f64,
    pub tau_c: f64,
    pub tau_allow: f64,
    pub burrillOk: bool,
    pub sigmaCharted: bool,
    pub earInSeries: bool,
    /// Every blade count's best at this shaft speed, best first.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub perZ: Vec<PerZ>,
}

impl Design {
    fn brief(&self) -> PerZ {
        PerZ {
            Z: self.Z,
            D: self.D,
            PD: self.PD,
            EAR: self.EAR,
            eta0: self.eta0,
            P_shaft: self.P_shaft,
        }
    }
}

// ------------------------------------------------------------- one geometry

const PD_BISECT_ITERS: usize = 40;

#[allow(clippy::too_many_arguments)]
fn solve_pitch(p: &Pitch, kt_req: f64, dkt_per_pd: f64, pd_lo: f64, pd_hi: f64) -> f64 {
    let mut a = p.a;
    a[1] += dkt_per_pd;
    let f_lo = polyval6(&a, pd_lo) - kt_req;
    let f_hi = polyval6(&a, pd_hi) - kt_req;
    if f_lo > 0.0 || f_hi < 0.0 {
        return f64::NAN;
    }
    let (mut lo, mut hi) = (pd_lo, pd_hi);
    for _ in 0..PD_BISECT_ITERS {
        let mid = 0.5 * (lo + hi);
        if polyval6(&a, mid) - kt_req < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Gutsche's oblique-flow / inclined-shaft correction: when the shaft angle
/// is nonzero, the inflow velocity decomposes into an axial component
/// V_A cos(eps) and a transverse component V_A sin(eps). As the blade sweeps
/// azimuth theta in [0, 2pi), the effective advance ratio J(theta) oscillates
/// and the local dynamic head varies. Integrating over the revolution gives
/// the effective degree-6 pitch polynomials for mean K_T and K_Q.
fn effective_pitch_polynomials(
    v_a: f64,
    shaft_angle_deg: f64,
    n: f64,
    d: f64,
    ear: f64,
    z: f64,
) -> (Pitch, f64) {
    if shaft_angle_deg == 0.0 {
        let j = v_a / (n * d);
        return (pitch_polynomials(j, ear, z), j);
    }
    let eps = shaft_angle_deg.to_radians();
    let (sin_e, cos_e) = eps.sin_cos();
    let v_ax = v_a * cos_e;
    let j_a = v_ax / (n * d);
    let lambda = (j_a * sin_e.abs() / cos_e.abs().max(1e-6)) / (0.7 * std::f64::consts::PI);
    let denom0 = j_a * j_a + (0.7 * std::f64::consts::PI).powi(2);
    const N: usize = 16;
    let mut p_eff = Pitch::default();
    for k in 0..N {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (N as f64);
        let s_th = theta.sin();
        let denom_th = (1.0 - lambda * s_th).max(0.05);
        let j_th = (j_a / denom_th).clamp(J_BOUNDS[0], J_BOUNDS[1]);
        let w_th = (j_a * j_a + (0.7 * std::f64::consts::PI * denom_th).powi(2)) / denom0;
        let p = pitch_polynomials(j_th, ear, z);
        for i in 0..7 {
            p_eff.a[i] += p.a[i] * w_th / (N as f64);
            p_eff.b[i] += p.b[i] * w_th / (N as f64);
        }
    }
    (p_eff, j_a)
}

/// A fixed propeller at the second operating point: the shaft speed that
/// makes its thrust there, or `None` if it can't inside the series.
fn at_condition(
    c: &Config,
    top: &TopConfig,
    d: f64,
    pd: f64,
    ear: f64,
    z: f64,
) -> Option<TopPoint> {
    let cos_e = c.shaft_angle_deg.to_radians().cos();
    let v_ax = top.v_a * cos_e;
    let thrust_at = |n: f64| -> Option<f64> {
        let (p, j) = effective_pitch_polynomials(top.v_a, c.shaft_angle_deg, n, d, ear, z);
        if !(J_BOUNDS[0]..=J_BOUNDS[1]).contains(&j) {
            return None;
        }
        let dkt = if c.re_correct {
            reynolds_correction(d, ear, z, top.v_a, n, c.shaft_angle_deg).dkt_per_pd * pd
        } else {
            0.0
        };
        let kt = polyval6(&p.a, pd) + dkt;
        (kt > 0.0).then(|| RHO * n * n * d.powi(4) * kt)
    };
    let mut lo = v_ax / (J_BOUNDS[1] * d);
    let mut hi = v_ax / (J_BOUNDS[0] * d);
    let want = top.t_target;
    match thrust_at(hi * 0.999) {
        Some(t) if t >= want => {}
        _ => return None,
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        match thrust_at(mid) {
            Some(t) if t >= want => hi = mid,
            _ => lo = mid,
        }
    }
    let n = 0.5 * (lo + hi);
    let (p, j) = effective_pitch_polynomials(top.v_a, c.shaft_angle_deg, n, d, ear, z);
    let (dkt, dkq) = if c.re_correct {
        let rc = reynolds_correction(d, ear, z, top.v_a, n, c.shaft_angle_deg);
        (rc.dkt_per_pd * pd, rc.dkq)
    } else {
        (0.0, 0.0)
    };
    let kt = polyval6(&p.a, pd) + dkt;
    let kq = polyval6(&p.b, pd) + dkq;
    if !(kq > 0.0) || !(kt > 0.0) {
        return None;
    }
    let q = RHO * n * n * d.powi(5) * kq;
    let p_shaft = 2.0 * std::f64::consts::PI * n * q;
    if !(p_shaft > 0.0) {
        return None;
    }
    let sigma = cavitation_number07(top.v_a, n, d, c.depth, c.shaft_angle_deg);
    let tau_c = burrill_tau_c(want, ear, pd, d, top.v_a, n, c.shaft_angle_deg);
    let tau_allow = burrill_tau_allow(sigma);
    let ear_min = keller_min_ear(want, d, z, c.depth, c.keller_k);
    Some(TopPoint {
        n,
        rpm: 60.0 * n,
        J: j,
        K_T: kt,
        K_Q: kq,
        T: want,
        Q: q,
        P_shaft: p_shaft,
        eta0: want * v_ax / p_shaft,
        sigma,
        tau_c,
        tau_allow,
        burrillOk: tau_c <= tau_allow,
        EAR_min: ear_min,
        kellerOk: ear >= ear_min,
    })
}

/// Thrust loading C_Th and the open-propeller ceiling it sets.
fn thrust_loading(t: f64, d: f64, v_a: f64) -> Option<(f64, f64)> {
    let a0 = std::f64::consts::PI * d * d / 4.0;
    if !(a0 > 0.0) || !(v_a > 0.0) {
        return None;
    }
    let cth = t / (0.5 * RHO * v_a * v_a * a0);
    if !(cth > 0.0) || !cth.is_finite() {
        return None;
    }
    Some((cth, 2.0 / (1.0 + (1.0 + cth).sqrt())))
}

/// One (D, EAR, Z) at shaft speed `n` [rev/s], pitch sized to the thrust.
pub fn evaluate(c: &Config, n: f64, d: f64, ear: f64, z: f64) -> Option<Design> {
    let (p, j) = effective_pitch_polynomials(c.v_a, c.shaft_angle_deg, n, d, ear, z);
    if !(J_BOUNDS[0]..=J_BOUNDS[1]).contains(&j) {
        return None;
    }
    let kt_req = c.t_target / (RHO * n * n * d.powi(4));
    if !(kt_req > 0.0) {
        return None;
    }
    let mut rc = None;
    let (mut dkt_per_pd, mut dkq) = (0.0, 0.0);
    if c.re_correct {
        let r = reynolds_correction(d, ear, z, c.v_a, n, c.shaft_angle_deg);
        dkt_per_pd = r.dkt_per_pd;
        dkq = r.dkq;
        rc = Some(r);
    }
    let pd = solve_pitch(&p, kt_req, dkt_per_pd, PD_BOUNDS[0], PD_BOUNDS[1]);
    if !(pd > 0.0) {
        return None;
    }
    let kq = polyval6(&p.b, pd) + dkq;
    if !(kq > 0.0) {
        return None;
    }
    let q = RHO * n * n * d.powi(5) * kq;
    let p_shaft = 2.0 * std::f64::consts::PI * n * q;
    if !(p_shaft > 0.0) {
        return None;
    }
    let cos_e = c.shaft_angle_deg.to_radians().cos();
    let eta0 = c.t_target * (c.v_a * cos_e) / p_shaft;
    if !(eta0 > 0.0) || !eta0.is_finite() {
        return None;
    }
    let load = thrust_loading(c.t_target, d, c.v_a * cos_e);
    let top = c
        .top
        .as_ref()
        .and_then(|t| at_condition(c, t, d, pd, ear, z));
    let rc = rc.unwrap_or_else(|| reynolds_correction(d, ear, z, c.v_a, n, c.shaft_angle_deg));
    let ear_min = keller_min_ear(c.t_target, d, z, c.depth, c.keller_k);
    let keller_ok = ear >= ear_min;
    let sigma = cavitation_number07(c.v_a, n, d, c.depth, c.shaft_angle_deg);
    let tau_c = burrill_tau_c(c.t_target, ear, pd, d, c.v_a, n, c.shaft_angle_deg);
    let tau_allow = burrill_tau_allow(sigma);
    let burrill_ok = tau_c <= tau_allow;
    if c.cavitation && !(keller_ok && burrill_ok) {
        return None;
    }
    if c.strict_ear && !ear_in_series_range(ear, z) {
        return None;
    }
    let top_ok = match &c.top {
        None => true,
        Some(_) => top
            .as_ref()
            .is_some_and(|t| !c.cavitation || (t.kellerOk && t.burrillOk)),
    };
    Some(Design {
        n,
        rpm: 60.0 * n,
        D: d,
        PD: pd,
        EAR: ear,
        Z: z,
        J: j,
        K_T: kt_req,
        K_Q: kq,
        T: c.t_target,
        Q: q,
        P_shaft: p_shaft,
        eta0,
        Cth: load.map(|l| l.0),
        etaIdeal: load.map(|l| l.1),
        etaOfIdeal: load.filter(|l| l.1 > 0.0).map(|l| eta0 / l.1),
        top,
        topOk: top_ok,
        Re: rc.re,
        dC_D: rc.d_cd,
        cdFormulaOk: rc.cd_formula_ok,
        EAR_min: ear_min,
        kellerOk: keller_ok,
        sigma,
        tau_c,
        tau_allow,
        burrillOk: burrill_ok,
        sigmaCharted: (BURRILL_SIGMA_RANGE[0]..=BURRILL_SIGMA_RANGE[1]).contains(&sigma),
        earInSeries: ear_in_series_range(ear, z),
        perZ: Vec::new(),
    })
}

// ------------------------------------------------------------- best at an rpm

#[allow(clippy::too_many_arguments)]
fn pattern_search(
    c: &Config,
    n: f64,
    z: f64,
    d0: f64,
    ear0: f64,
    (d_lo, d_hi): (f64, f64),
    (e_lo, e_hi): (f64, f64),
    best: Option<Design>,
) -> Option<Design> {
    let (mut d, mut ear) = (d0, ear0);
    let Some(seed) = evaluate(c, n, d, ear, z) else {
        return best;
    };
    let mut best = match best {
        Some(b) if seed.eta0 <= b.eta0 => {
            d = b.D;
            ear = b.EAR;
            b
        }
        _ => seed,
    };
    let (mut sd, mut se) = ((d_hi - d_lo) * 0.08, (e_hi - e_lo) * 0.08);
    let (min_d, min_e) = ((d_hi - d_lo) * 1e-5, (e_hi - e_lo) * 1e-5);
    const AXES: [(f64, f64); 4] = [(1., 0.), (-1., 0.), (0., 1.), (0., -1.)];
    const DIAGS: [(f64, f64); 4] = [(1., 1.), (1., -1.), (-1., 1.), (-1., -1.)];
    let mut it = 0;
    while it < 400 && (sd > min_d || se > min_e) {
        let mut moved = false;
        for set in [&AXES, &DIAGS] {
            for &(ud, ue) in set.iter() {
                let d2 = d_hi.min(d_lo.max(d + ud * sd));
                let e2 = e_hi.min(e_lo.max(ear + ue * se));
                if d2 == d && e2 == ear {
                    continue;
                }
                if let Some(r) = evaluate(c, n, d2, e2, z) {
                    if r.eta0 > best.eta0 {
                        best = r;
                        d = d2;
                        ear = e2;
                        moved = true;
                    }
                }
            }
            if moved {
                break;
            }
        }
        if !moved {
            sd *= 0.5;
            se *= 0.5;
        }
        it += 1;
    }
    Some(best)
}

/// The best geometry at an rpm, with every blade count's best (full
/// designs, best first; the best carries them in brief as `perZ`).
#[derive(Clone, Debug)]
pub struct AtRpm {
    pub best: Design,
    pub per_z: Vec<Design>,
}

/// The best geometry at `rpm`, per blade count and overall. `seeds` is the
/// neighbouring rpm's answer per blade count, `(Z, D, EAR)`, to start from.
pub fn at_rpm(c: &Config, rpm: f64, seeds: Option<&[(f64, f64, f64)]>) -> Option<AtRpm> {
    let n = rpm / 60.0;
    let mut best: Option<usize> = None;
    let mut per_z: Vec<Design> = Vec::new();
    for &z in &c.z {
        let (mut e_lo, mut e_hi) = (EAR_BOUNDS[0], EAR_BOUNDS[1]);
        if c.strict_ear {
            let Some(r) = ear_range(z) else { continue };
            e_lo = e_lo.max(r[0]);
            e_hi = e_hi.min(r[1]);
        }
        let (d_lo, d_hi) = (c.d_min, c.d_max);
        let warm = seeds.and_then(|s| s.iter().find(|w| w.0 == z));
        // Warm start first, so the coarse sweep only insures against it.
        let mut bz = None;
        if let Some(&(_, wd, we)) = warm {
            bz = pattern_search(
                c,
                n,
                z,
                d_hi.min(d_lo.max(wd)),
                e_hi.min(e_lo.max(we)),
                (d_lo, d_hi),
                (e_lo, e_hi),
                None,
            );
        }
        let n_d = if bz.is_some() {
            c.grid_d.div_ceil(2)
        } else {
            c.grid_d
        };
        let n_e = if bz.is_some() {
            c.grid_ear.div_ceil(2)
        } else {
            c.grid_ear
        };
        let mut seed: Option<Design> = None;
        for i in 0..n_d {
            let d = if n_d == 1 {
                d_lo
            } else {
                d_lo + (d_hi - d_lo) * i as f64 / (n_d - 1) as f64
            };
            for k in 0..n_e {
                let ear = if n_e == 1 {
                    e_lo
                } else {
                    e_lo + (e_hi - e_lo) * k as f64 / (n_e - 1) as f64
                };
                if let Some(r) = evaluate(c, n, d, ear, z) {
                    if seed.as_ref().is_none_or(|s| r.eta0 > s.eta0) {
                        seed = Some(r);
                    }
                }
            }
        }
        if let Some(s) = seed {
            bz = pattern_search(c, n, z, s.D, s.EAR, (d_lo, d_hi), (e_lo, e_hi), bz);
        }
        let Some(bz) = bz else { continue };
        if best.is_none_or(|b| bz.eta0 > per_z[b].eta0) {
            best = Some(per_z.len());
        }
        per_z.push(bz);
    }
    let mut best = per_z[best?].clone();
    // Stable, as JS's sort is: ties keep their blade-count order.
    per_z.sort_by(|a, b| b.eta0.total_cmp(&a.eta0));
    best.perZ = per_z.iter().map(Design::brief).collect();
    Some(AtRpm { best, per_z })
}

/// The best geometry at `rpm` (see [`at_rpm`]).
pub fn best_at_rpm(c: &Config, rpm: f64, seeds: Option<&[(f64, f64, f64)]>) -> Option<Design> {
    at_rpm(c, rpm, seeds).map(|a| a.best)
}

// ------------------------------------------------------------- the sweep

/// A point on the curve: the best propeller at this rpm, or none. As JSON,
/// the design's fields with `"ok": true`, or `{"rpm": …, "ok": false}`.
#[derive(Clone, Debug)]
pub struct CurvePoint {
    pub rpm: f64,
    pub ok: bool,
    pub design: Option<Design>,
}

impl Serialize for CurvePoint {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::Error;
        let v = match &self.design {
            Some(d) => {
                let mut v = serde_json::to_value(d).map_err(S::Error::custom)?;
                v["ok"] = serde_json::Value::Bool(true);
                v
            }
            None => serde_json::json!({ "rpm": self.rpm, "ok": false }),
        };
        v.serialize(s)
    }
}

impl<'de> Deserialize<'de> for CurvePoint {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        use serde::de::Error;
        let v = serde_json::Value::deserialize(d)?;
        let rpm = v["rpm"]
            .as_f64()
            .ok_or_else(|| D::Error::custom("a curve point without rpm"))?;
        if v["ok"] == true {
            let design: Design = serde_json::from_value(v).map_err(D::Error::custom)?;
            Ok(CurvePoint {
                rpm,
                ok: true,
                design: Some(design),
            })
        } else {
            Ok(CurvePoint {
                rpm,
                ok: false,
                design: None,
            })
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct Sweep {
    pub feasible: bool,
    pub rpmLo: Option<f64>,
    pub rpmHi: Option<f64>,
    pub feasLo: Option<f64>,
    pub feasHi: Option<f64>,
    /// The cheapest design that also reaches the second operating point.
    pub best: Option<Design>,
    /// With a second point, the cheapest ignoring it.
    pub bestUnconstrained: Option<Design>,
    pub topOkCount: usize,
    pub topFeasible: bool,
    pub points: Vec<CurvePoint>,
}

pub struct SweepOptions {
    pub probe_min: f64,
    pub probe_max: f64,
    pub probes: usize,
    pub samples: usize,
    /// The displayed window runs out from the optimum until shaft power
    /// passes this multiple of it.
    pub power_cap: f64,
}

impl Default for SweepOptions {
    fn default() -> Self {
        SweepOptions {
            probe_min: 20.0,
            probe_max: 40000.0,
            probes: 90,
            samples: 180,
            power_cap: 1.6,
        }
    }
}

fn refine_edge(c: &Config, mut bad: f64, mut good: f64) -> f64 {
    for _ in 0..22 {
        let mid = 0.5 * (bad + good);
        if best_at_rpm(c, mid, None).is_some() {
            good = mid;
        } else {
            bad = mid;
        }
    }
    good
}

fn refine_optimum(c: &Config, best: &Design, span: f64) -> Design {
    const PHI: f64 = 0.6180339887;
    let (mut lo, mut hi) = (best.rpm - span, best.rpm + span);
    let p = |rpm: f64| best_at_rpm(c, rpm, None);
    let mut x1 = hi - PHI * (hi - lo);
    let mut x2 = lo + PHI * (hi - lo);
    let mut f1 = p(x1);
    let mut f2 = p(x2);
    let mut i = 0;
    while i < 24 && hi - lo > 1e-4 {
        let p1 = f1.as_ref().map_or(f64::INFINITY, |f| f.P_shaft);
        let p2 = f2.as_ref().map_or(f64::INFINITY, |f| f.P_shaft);
        if p1 < p2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - PHI * (hi - lo);
            f1 = p(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + PHI * (hi - lo);
            f2 = p(x2);
        }
        i += 1;
    }
    let mut win = best.clone();
    for f in [f1, f2].into_iter().flatten() {
        if f.P_shaft < win.P_shaft {
            win = f;
        }
    }
    win
}

/// Sweep shaft speed: the best propeller at each rpm, and the cheapest.
pub fn sweep(c: &Config, o: &SweepOptions) -> Sweep {
    let (probe_lo, probe_hi, n_probe, n_final) = (o.probe_min, o.probe_max, o.probes, o.samples);
    let none = Sweep {
        feasible: false,
        rpmLo: None,
        rpmHi: None,
        feasLo: None,
        feasHi: None,
        best: None,
        bestUnconstrained: None,
        topOkCount: 0,
        topFeasible: false,
        points: Vec::new(),
    };

    // Pass 1: where is there an answer, and where is it cheapest?
    let mut probe: Vec<(f64, Option<Design>)> = Vec::with_capacity(n_probe);
    let (mut feas_lo, mut feas_hi) = (None, None);
    let mut best_probe: Option<Design> = None;
    for i in 0..n_probe {
        let rpm = probe_lo * (probe_hi / probe_lo).powf(i as f64 / (n_probe - 1) as f64);
        let r = best_at_rpm(c, rpm, None);
        if let Some(r) = &r {
            if feas_lo.is_none() {
                feas_lo = Some(rpm);
            }
            feas_hi = Some(rpm);
            if best_probe.as_ref().is_none_or(|b| r.P_shaft < b.P_shaft) {
                best_probe = Some(r.clone());
            }
        }
        probe.push((rpm, r));
    }
    let Some(best_probe) = best_probe else {
        return none;
    };
    let ratio = (probe_hi / probe_lo).powf(1.0 / (n_probe - 1) as f64);
    let feas_lo = refine_edge(c, feas_lo.unwrap() / ratio, feas_lo.unwrap());
    let feas_hi = refine_edge(c, feas_hi.unwrap() * ratio, feas_hi.unwrap());

    // The displayed window: out from the optimum until power passes the cap.
    let (mut rpm_lo, mut rpm_hi) = (feas_lo, feas_hi);
    let cap = best_probe.P_shaft * o.power_cap;
    for (rpm, r) in probe.iter().rev() {
        if *rpm < best_probe.rpm && r.as_ref().is_some_and(|r| r.P_shaft > cap) {
            rpm_lo = *rpm;
            break;
        }
    }
    for (rpm, r) in probe.iter() {
        if *rpm > best_probe.rpm && r.as_ref().is_some_and(|r| r.P_shaft > cap) {
            rpm_hi = *rpm;
            break;
        }
    }
    let rpm_lo = rpm_lo.max(feas_lo);
    let rpm_hi = rpm_hi.min(feas_hi);

    // Pass 2: the curve, walked outward from the cheapest rpm with warm starts.
    let rpms: Vec<f64> = (0..n_final)
        .map(|i| rpm_lo + (rpm_hi - rpm_lo) * i as f64 / (n_final - 1) as f64)
        .collect();
    let mut mid = 0;
    for i in 1..n_final {
        if (rpms[i] - best_probe.rpm).abs() < (rpms[mid] - best_probe.rpm).abs() {
            mid = i;
        }
    }
    let mut raw: Vec<Option<AtRpm>> = vec![None; n_final];
    let seeds_of = |r: &AtRpm| -> Vec<(f64, f64, f64)> {
        r.best.perZ.iter().map(|z| (z.Z, z.D, z.EAR)).collect()
    };
    for dir in [1isize, -1] {
        let mut seeds: Option<Vec<(f64, f64, f64)>> = None;
        let mut i = if dir == 1 {
            mid as isize
        } else {
            mid as isize - 1
        };
        while i >= 0 && (i as usize) < n_final {
            let r = at_rpm(c, rpms[i as usize], seeds.as_deref());
            seeds = r.as_ref().map(seeds_of);
            raw[i as usize] = r;
            i += dir;
        }
    }

    // Keep the blade count the curve already has while it stays within
    // HYSTERESIS of the best, so a change of Z is a real change of answer.
    let mut points: Vec<Option<CurvePoint>> = vec![None; n_final];
    let (mut best, mut best_any): (Option<Design>, Option<Design>) = (None, None);
    for dir in [1isize, -1] {
        let mut held_z: Option<f64> = None;
        let mut i = if dir == 1 {
            mid as isize
        } else {
            mid as isize - 1
        };
        while i >= 0 && (i as usize) < n_final {
            let k = i as usize;
            let Some(at) = raw[k].clone() else {
                held_z = None;
                points[k] = Some(CurvePoint {
                    rpm: rpms[k],
                    ok: false,
                    design: None,
                });
                i += dir;
                continue;
            };
            let mut r = at.best.clone();
            if let Some(hz) = held_z {
                if r.Z != hz {
                    if let Some(held) = at.per_z.iter().find(|x| x.Z == hz) {
                        if held.eta0 >= r.eta0 * (1.0 - HYSTERESIS) {
                            let mut h = held.clone();
                            h.perZ = r.perZ.clone();
                            r = h;
                        }
                    }
                }
            }
            held_z = Some(r.Z);
            if r.topOk && best.as_ref().is_none_or(|b| r.P_shaft < b.P_shaft) {
                best = Some(r.clone());
            }
            if best_any.as_ref().is_none_or(|b| r.P_shaft < b.P_shaft) {
                best_any = Some(r.clone());
            }
            points[k] = Some(CurvePoint {
                rpm: rpms[k],
                ok: true,
                design: Some(r),
            });
            i += dir;
        }
    }
    let span = (rpm_hi - rpm_lo) / (n_final - 1) as f64;
    if let Some(b) = &best {
        let rf = refine_optimum(c, b, span);
        if rf.topOk || c.top.is_none() {
            best = Some(rf);
        }
    }
    let same = match (&best, &best_any) {
        (Some(a), Some(b)) => a.rpm == b.rpm && a.Z == b.Z && a.P_shaft == b.P_shaft,
        _ => false,
    };
    if let Some(b) = &best_any {
        if !same {
            best_any = Some(refine_optimum(c, b, span));
        }
    }
    let points: Vec<CurvePoint> = points
        .into_iter()
        .map(|p| p.expect("every sample"))
        .collect();
    let top_ok_count = points
        .iter()
        .filter(|p| p.ok && p.design.as_ref().is_some_and(|d| d.topOk))
        .count();
    Sweep {
        feasible: true,
        rpmLo: Some(rpm_lo),
        rpmHi: Some(rpm_hi),
        feasLo: Some(feas_lo),
        feasHi: Some(feas_hi),
        best,
        bestUnconstrained: if c.top.is_some() { best_any } else { None },
        topOkCount: top_ok_count,
        topFeasible: c.top.is_none() || top_ok_count > 0,
        points,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_inputs() -> Inputs {
        Inputs {
            speed: 4.0,
            thrust: 230.0,
            shafts: 2,
            wake: 0.1,
            thrust_deduction: 0.05,
            d_min: 0.04,
            d_max: 0.35,
            blades: vec![3, 4],
            depth: 0.25,
            keller_k: 0.1,
            cavitation: true,
            re_correct: true,
            strict_ear: false,
            top: None,
            shaft_angle_deg: 0.0,
        }
    }

    #[test]
    fn test_zero_shaft_angle_matches_axial() {
        let i0 = base_inputs();
        let c0 = Config::new(&i0);
        let d0 = evaluate(&c0, 20.0, 0.25, 0.50, 3.0).unwrap();

        let mut i_zero = base_inputs();
        i_zero.shaft_angle_deg = 0.0;
        let c_zero = Config::new(&i_zero);
        let d_zero = evaluate(&c_zero, 20.0, 0.25, 0.50, 3.0).unwrap();

        assert_eq!(d0.J, d_zero.J);
        assert_eq!(d0.PD, d_zero.PD);
        assert_eq!(d0.K_T, d_zero.K_T);
        assert_eq!(d0.K_Q, d_zero.K_Q);
        assert_eq!(d0.P_shaft, d_zero.P_shaft);
        assert_eq!(d0.eta0, d_zero.eta0);
        assert_eq!(d0.sigma, d_zero.sigma);
    }

    #[test]
    fn test_oblique_flow_correction_effects() {
        let i0 = base_inputs();
        let c0 = Config::new(&i0);
        let d0 = evaluate(&c0, 20.0, 0.25, 0.50, 3.0).unwrap();

        let mut i_oblique = base_inputs();
        i_oblique.shaft_angle_deg = 8.0; // 8 degrees inclined shaft
        let c_oblique = Config::new(&i_oblique);
        let d_ob = evaluate(&c_oblique, 20.0, 0.25, 0.50, 3.0).unwrap();

        // Cross-flow increases peak relative blade velocity, reducing cavitation index sigma:
        assert!(
            d_ob.sigma < d0.sigma,
            "sigma: {} vs {}",
            d_ob.sigma,
            d0.sigma
        );

        // Axial advance ratio is J * cos(8 deg) < J:
        let cos_8 = 8.0f64.to_radians().cos();
        assert!((d_ob.J - d0.J * cos_8).abs() < 1e-12);

        // Due to shaft inclination, useful thrust power along direction of advance is reduced by cos(eps):
        assert!(d_ob.eta0 < d0.eta0);

        // Effective pitch required to hold thrust changes with cross-flow:
        assert!(d_ob.PD != d0.PD);
    }
}
