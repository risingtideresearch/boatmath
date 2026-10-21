//! Hulls from camber's documents: its saved JSON (format version 2), swept
//! as camber sweeps it and fitted with B-spline patches, so `hull` reads a
//! camber file as it reads an IGES one.
//!
//! camber (risingtideresearch/camber) describes a hull as a section swept
//! along its sheer: a sheer plan (a clamped B-spline of half-breadth in the
//! deck plane, whose own parameter u places everything), stations of (n, z)
//! points at given u lofted between (a Catmull–Rom knotted at the stations'
//! u, the knuckles by PCHIP), each section drawn through its lofted points
//! (a centripetal Catmull–Rom creased by the knuckles) in the vertical plane
//! normal to the plan, and the sheet trimmed by the sheer trim (a PCHIP
//! graph z(x)), the centreline and the transom plane. This is a port of that
//! sweep (`src/core/{bspline,pchip,spline,model,mesh}.ts`), held to camber's
//! own output by golden tests (`tests/camber/make.mts`).
//!
//! The swept hull is sampled on a grid of trimmed sections and interpolated,
//! as camber's STEP export does it (Piegl & Tiller's global interpolation,
//! creased at the knuckle rows), into
//!
//! - a patch over the sections that close on the centreline;
//! - behind a raked transom, where sections end on the transom instead,
//!   another;
//! - the transom (or a plan's open end): a ruled face across from its
//!   starboard edge to the centreline;
//!
//! each starboard and mirrored to port. The halves meet at the keel without
//! a common tangent, so a V keel stays a V; and no patch spans the
//! centreplane, which the solver finds by the patches either side of it.
//!
//! The patches are moved into boatmath's frame: metres, floated at the
//! document's deck trim, with its design waterline at z = 0. camber's
//! keel knuckle (`keelK`) is not read, as camber's own sweep doesn't read it
//! yet.

use hullgeom::iges::NurbsSurface3;
use serde_json::Value;

type V2 = [f64; 2];
type V3 = [f64; 3];

/// Whether a JSON value is a camber hull document.
pub fn is_document(v: &Value) -> bool {
    v.get("sheerPlan").is_some()
}

/// Whether these bytes are a camber hull document.
pub fn is_camber(bytes: &[u8]) -> bool {
    if !crate::native::is_native(bytes) {
        return false;
    }
    serde_json::from_slice::<Value>(bytes).is_ok_and(|v| is_document(&v))
}

// ---------------------------------------------------------------------------
// The document
// ---------------------------------------------------------------------------

struct Station {
    u: f64,
    /// (n, z, k) per point, deck point first.
    pts: Vec<[f64; 3]>,
}

/// A camber hull document, in its own unit and frame: x forward from the
/// sheer plan's start, y the half-breadth, z up from the deck datum.
pub struct Document {
    pub name: String,
    /// Metres per document unit.
    pub unit_m: f64,
    /// The design waterline's depth below the deck datum.
    pub waterline: f64,
    /// The deck trim [rad], bow up.
    pub deck_trim: f64,
    plan: Vec<V2>,
    trim: Vec<V2>,
    transom: [V2; 2],
    stations: Vec<Station>,
}

fn num(v: &Value, at: &str) -> Result<f64, String> {
    v.as_f64()
        .filter(|x| x.is_finite())
        .ok_or_else(|| format!("{at} must be a finite number"))
}

fn arr<'a>(v: &'a Value, at: &str, min: usize) -> Result<&'a Vec<Value>, String> {
    let a = v
        .as_array()
        .ok_or_else(|| format!("{at} must be an array"))?;
    if a.len() < min {
        return Err(format!("{at} must have at least {min} entries"));
    }
    Ok(a)
}

/// A knuckle: absent or not a number reads as smooth, and it's held to [0, 1].
fn knuckle(v: &Value) -> f64 {
    v.as_f64()
        .filter(|x| x.is_finite())
        .map_or(0.0, |k| k.clamp(0.0, 1.0))
}

impl Document {
    /// Read and check a document, as camber's `parseDocument` does. Only
    /// version 2, which camber writes, is read.
    pub fn parse(v: &Value) -> Result<Document, String> {
        if !v.is_object() {
            return Err("a camber document must be an object".into());
        }
        let version = match v.get("version") {
            None => 1.0,
            Some(x) => x.as_f64().unwrap_or(f64::NAN),
        };
        if version.is_nan() || version > 2.0 {
            return Err(format!(
                "camber document version {}: this reads version 2",
                v["version"]
            ));
        }
        if version < 2.0 {
            return Err(
                "a version 1 camber document: open it in camber and save it, \
                        which writes version 2"
                    .into(),
            );
        }
        let unit_m = match v["unit"].as_str() {
            Some("mm") => 0.001,
            Some("cm") => 0.01,
            Some("m") => 1.0,
            Some("in") => 0.0254,
            Some("ft") => 0.3048,
            _ => return Err("unit must be one of mm, cm, m, in, ft".into()),
        };
        let plan = arr(&v["sheerPlan"], "sheerPlan", 2)?
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let at = format!("sheerPlan[{i}]");
                Ok([
                    num(&p["x"], &format!("{at}.x"))?,
                    num(&p["y"], &format!("{at}.y"))?,
                ])
            })
            .collect::<Result<Vec<V2>, String>>()?;
        let trim = arr(&v["sheerTrim"], "sheerTrim", 2)?
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let at = format!("sheerTrim[{i}]");
                Ok([
                    num(&p["x"], &format!("{at}.x"))?,
                    num(&p["z"], &format!("{at}.z"))?,
                ])
            })
            .collect::<Result<Vec<V2>, String>>()?;
        for (name, pts) in [("sheerPlan", &plan), ("sheerTrim", &trim)] {
            for i in 1..pts.len() {
                if pts[i][0] <= pts[i - 1][0] {
                    return Err(format!("{name}[{i}].x must be > the previous point's x"));
                }
            }
        }
        let tr = arr(&v["transom"], "transom", 2)?;
        if tr.len() != 2 {
            return Err("transom must have exactly 2 points (top and bottom)".into());
        }
        let tp = |i: usize| -> Result<V2, String> {
            let at = format!("transom[{i}]");
            Ok([
                num(&tr[i]["x"], &format!("{at}.x"))?,
                num(&tr[i]["z"], &format!("{at}.z"))?,
            ])
        };
        let transom = [tp(0)?, tp(1)?];
        if transom[1][1] >= transom[0][1] {
            return Err("transom[1] (the bottom) must be below transom[0]".into());
        }
        let sts = arr(&v["stations"], "stations", 1)?;
        let n_pts = arr(&sts[0]["points"], "stations[0].points", 2)?.len();
        let mut stations = sts
            .iter()
            .enumerate()
            .map(|(j, s)| {
                let pts = arr(&s["points"], &format!("stations[{j}].points"), 2)?;
                if pts.len() != n_pts {
                    return Err(format!(
                        "stations[{j}].points must have {n_pts} points (every station shares one count)"
                    ));
                }
                Ok(Station {
                    u: num(&s["u"], &format!("stations[{j}].u"))?.clamp(0.0, 1.0),
                    pts: pts
                        .iter()
                        .enumerate()
                        .map(|(i, q)| {
                            let at = format!("stations[{j}].points[{i}]");
                            // The ends are corners, whatever was written.
                            let k = if i == 0 || i == n_pts - 1 { 1.0 } else { knuckle(&q["k"]) };
                            Ok([num(&q["n"], &format!("{at}.n"))?, num(&q["z"], &format!("{at}.z"))?, k])
                        })
                        .collect::<Result<Vec<_>, String>>()?,
                })
            })
            .collect::<Result<Vec<Station>, String>>()?;
        stations.sort_by(|a, b| a.u.total_cmp(&b.u));
        for j in 1..stations.len() {
            if stations[j].u <= stations[j - 1].u {
                return Err(format!(
                    "stations: two share the same u ({})",
                    stations[j].u
                ));
            }
        }
        let opt = |k: &str| v.get(k).and_then(Value::as_f64).filter(|x| x.is_finite());
        let deck_trim_deg = opt("deckTrimDeg")
            .or_else(|| opt("deckRakeDeg"))
            .unwrap_or(0.0);
        Ok(Document {
            name: v["name"].as_str().unwrap_or("").to_string(),
            unit_m,
            waterline: opt("waterline").unwrap_or(0.0),
            deck_trim: deck_trim_deg.to_radians(),
            plan,
            trim,
            transom,
            stations,
        })
    }

    /// The hull's length overall, plan start to bow.
    fn loa(&self) -> f64 {
        self.plan[self.plan.len() - 1][0] - self.plan[0][0]
    }
}

// ---------------------------------------------------------------------------
// Curves (camber's bspline.ts, pchip.ts, spline.ts)
// ---------------------------------------------------------------------------

fn find_span(n: usize, p: usize, u: f64, knots: &[f64]) -> usize {
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[p] {
        return p;
    }
    let (mut lo, mut hi) = (p, n + 1);
    let mut mid = (lo + hi) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

fn de_boor(pts: &[V2], knots: &[f64], p: usize, u: f64) -> V2 {
    let span = find_span(pts.len() - 1, p, u, knots);
    let mut d: Vec<V2> = (0..=p).map(|j| pts[span - p + j]).collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = span - p + j;
            let den = knots[i + p - r + 1] - knots[i];
            let a = if den > 0.0 { (u - knots[i]) / den } else { 0.0 };
            d[j] = [
                (1.0 - a) * d[j - 1][0] + a * d[j][0],
                (1.0 - a) * d[j - 1][1] + a * d[j][1],
            ];
        }
    }
    d[p]
}

fn uniform_knots(count: usize, p: usize) -> Vec<f64> {
    let interior = count - p - 1;
    let mut k = vec![0.0; p + 1];
    k.extend((1..=interior).map(|i| i as f64 / (interior + 1) as f64));
    k.extend(std::iter::repeat_n(1.0, p + 1));
    k
}

/// The sheer plan: a clamped B-spline P(u) = (x, y) on uniform knots, its
/// degree dropping for short polygons, with its exact derivative.
struct Plan {
    pts: Vec<V2>,
    p: usize,
    knots: Vec<f64>,
    d_pts: Vec<V2>,
    d_knots: Vec<f64>,
}

impl Plan {
    fn new(pts: Vec<V2>) -> Plan {
        let n = pts.len();
        let p = 3.min(n - 1);
        let knots = uniform_knots(n, p);
        let d_pts = (0..n - 1)
            .map(|i| {
                let den = knots[i + p + 1] - knots[i + 1];
                if den > 0.0 {
                    let s = p as f64 / den;
                    [
                        s * (pts[i + 1][0] - pts[i][0]),
                        s * (pts[i + 1][1] - pts[i][1]),
                    ]
                } else {
                    [0.0, 0.0]
                }
            })
            .collect();
        let d_knots = knots[1..knots.len() - 1].to_vec();
        Plan {
            pts,
            p,
            knots,
            d_pts,
            d_knots,
        }
    }

    fn at(&self, u: f64) -> V2 {
        de_boor(&self.pts, &self.knots, self.p, u.clamp(0.0, 1.0))
    }

    fn d(&self, u: f64) -> V2 {
        if self.p < 2 {
            self.d_pts[0]
        } else {
            de_boor(&self.d_pts, &self.d_knots, self.p - 1, u.clamp(0.0, 1.0))
        }
    }
}

fn c2_ease_in(q: f64) -> f64 {
    q * q * q * (6.0 + q * (3.0 * q - 8.0))
}

fn c2_ease_in_out(q: f64) -> f64 {
    q * q * q * (10.0 + q * (6.0 * q - 15.0))
}

const MID_EASE: f64 = 0.25;
const END_EASE: f64 = 0.25;

/// camber's C²-eased Fritsch–Carlson slopes.
fn pchip_slopes(xs: &[f64], ys: &[f64]) -> Vec<f64> {
    let n = xs.len();
    if n == 1 {
        return vec![0.0];
    }
    let h: Vec<f64> = (0..n - 1).map(|i| xs[i + 1] - xs[i]).collect();
    let d: Vec<f64> = (0..n - 1).map(|i| (ys[i + 1] - ys[i]) / h[i]).collect();
    if n == 2 {
        return vec![d[0], d[0]];
    }
    let mut m = vec![0.0; n];
    for i in 1..n - 1 {
        let (a, b) = (d[i - 1], d[i]);
        if a * b <= 0.0 {
            m[i] = 0.0;
        } else {
            let w1 = 2.0 * h[i] + h[i - 1];
            let w2 = h[i] + 2.0 * h[i - 1];
            let mut mi = (w1 + w2) / (w1 / a + w2 / b);
            let r = a.abs().min(b.abs()) / a.abs().max(b.abs());
            if r < MID_EASE {
                mi *= c2_ease_in_out(r / MID_EASE);
            }
            m[i] = mi;
        }
    }
    m[0] = pchip_end(h[0], h[1], d[0], d[1]);
    m[n - 1] = pchip_end(h[n - 2], h[n - 3], d[n - 2], d[n - 3]);
    m
}

fn pchip_end(h0: f64, h1: f64, d0: f64, d1: f64) -> f64 {
    let m = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
    let s = 3.0 * d0;
    if s == 0.0 {
        return 0.0;
    }
    let t = m / s;
    if t <= 0.0 {
        return 0.0;
    }
    if t >= 1.0 {
        return s;
    }
    if t < END_EASE {
        return s * END_EASE * c2_ease_in(t / END_EASE);
    }
    if t > 1.0 - END_EASE {
        return s * (1.0 - END_EASE * c2_ease_in((1.0 - t) / END_EASE));
    }
    m
}

fn hermite(xs: &[f64], ys: &[f64], m: &[f64], tt: f64) -> f64 {
    let mut i = 0;
    while i + 2 < xs.len() && tt > xs[i + 1] {
        i += 1;
    }
    let h = xs[i + 1] - xs[i];
    let t = (tt - xs[i]) / h;
    let (t2, t3) = (t * t, t * t * t);
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * ys[i] + h10 * h * m[i] + h01 * ys[i + 1] + h11 * h * m[i + 1]
}

/// A cubic Bézier segment of a Catmull–Rom chain.
type Bez = [V2; 4];

fn add(a: V2, b: V2) -> V2 {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn scale(a: V2, s: f64) -> V2 {
    [a[0] * s, a[1] * s]
}
fn lerp2(a: V2, b: V2, t: f64) -> V2 {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]
}

fn centripetal_params(pts: &[V2]) -> Vec<f64> {
    let mut t = vec![0.0];
    for i in 1..pts.len() {
        let d = sub(pts[i], pts[i - 1]);
        t.push(t[i - 1] + d[0].hypot(d[1]).sqrt().max(1e-6));
    }
    t
}

fn cr_tangent(p0: V2, p1: V2, p2: V2, t0: f64, t1: f64, t2: f64) -> V2 {
    let a = scale(sub(p1, p0), 1.0 / (t1 - t0));
    let b = scale(sub(p2, p0), 1.0 / (t2 - t0));
    let c = scale(sub(p2, p1), 1.0 / (t2 - t1));
    add(sub(a, b), c)
}

/// The Catmull–Rom chain through `vals` on knots `t`, creased by `ks`
/// (camber's `crChain`): the ends count as corners.
fn cr_chain(vals: &[V2], t: &[f64], ks: &[f64]) -> Vec<Bez> {
    let n = vals.len();
    (0..n - 1)
        .map(|j| {
            let (p, q) = (vals[j], vals[j + 1]);
            let chord = sub(q, p);
            let dt = t[j + 1] - t[j];
            let i1 =
                (j > 0).then(|| scale(cr_tangent(vals[j - 1], p, q, t[j - 1], t[j], t[j + 1]), dt));
            let i2 = (j + 2 < n)
                .then(|| scale(cr_tangent(p, q, vals[j + 2], t[j], t[j + 1], t[j + 2]), dt));
            let e1 = i2.map_or(chord, |m| sub(scale(chord, 2.0), m));
            let e2 = i1.map_or(chord, |m| sub(scale(chord, 2.0), m));
            let kp = if j == 0 { 1.0 } else { ks[j] };
            let kq = if j + 2 == n { 1.0 } else { ks[j + 1] };
            let (m1, m2) = (i1.unwrap_or(e1), i2.unwrap_or(e2));
            let (cp, cq) = (kp.clamp(0.0, 1.0), kq.clamp(0.0, 1.0));
            let t1 = lerp2(e1, chord, cq);
            let t2 = lerp2(e2, chord, cp);
            let b1 = add(p, scale(lerp2(m1, t1, cp), 1.0 / 3.0));
            let b2 = sub(q, scale(lerp2(m2, t2, cq), 1.0 / 3.0));
            [p, b1, b2, q]
        })
        .collect()
}

/// The chain at v ∈ [0, segments], knot j at v = j.
fn eval_chain(segs: &[Bez], v: f64) -> V2 {
    let n = segs.len();
    let c = v.clamp(0.0, n as f64);
    let j = (c.floor() as usize).min(n - 1);
    let s = c - j as f64;
    let [b0, b1, b2, b3] = segs[j];
    let u = 1.0 - s;
    let (w0, w1, w2, w3) = (u * u * u, 3.0 * u * u * s, 3.0 * u * s * s, s * s * s);
    [
        w0 * b0[0] + w1 * b1[0] + w2 * b2[0] + w3 * b3[0],
        w0 * b0[1] + w1 * b1[1] + w2 * b2[1] + w3 * b3[1],
    ]
}

// ---------------------------------------------------------------------------
// The sweep (camber's model.ts and mesh.ts)
// ---------------------------------------------------------------------------

/// The stations lofted along u: each point's (n, z) by a Catmull–Rom
/// knotted at the stations' u, its knuckle by PCHIP.
enum Loft {
    One {
        pts: Vec<V2>,
        ks: Vec<f64>,
    },
    Many {
        us: Vec<f64>,
        chains: Vec<Vec<Bez>>,
        k_ys: Vec<Vec<f64>>,
        k_ms: Vec<Vec<f64>>,
    },
}

impl Loft {
    fn new(sts: &[Station]) -> Loft {
        let n = sts[0].pts.len();
        if sts.len() == 1 {
            return Loft::One {
                pts: sts[0].pts.iter().map(|p| [p[0], p[1]]).collect(),
                ks: sts[0].pts.iter().map(|p| p[2]).collect(),
            };
        }
        let us: Vec<f64> = sts.iter().map(|s| s.u).collect();
        let zeros = vec![0.0; sts.len()];
        let chains = (0..n)
            .map(|i| {
                let vals: Vec<V2> = sts.iter().map(|s| [s.pts[i][0], s.pts[i][1]]).collect();
                cr_chain(&vals, &us, &zeros)
            })
            .collect();
        let k_ys: Vec<Vec<f64>> = (0..n)
            .map(|i| sts.iter().map(|s| s.pts[i][2]).collect())
            .collect();
        let k_ms = k_ys.iter().map(|ys| pchip_slopes(&us, ys)).collect();
        Loft::Many {
            us,
            chains,
            k_ys,
            k_ms,
        }
    }

    fn at(&self, u: f64) -> (Vec<V2>, Vec<f64>) {
        match self {
            Loft::One { pts, ks } => (pts.clone(), ks.clone()),
            Loft::Many {
                us,
                chains,
                k_ys,
                k_ms,
            } => {
                let k = us.len();
                let uc = u.clamp(us[0], us[k - 1]);
                let mut j = 0;
                while j + 2 < k && uc > us[j + 1] {
                    j += 1;
                }
                let span = us[j + 1] - us[j];
                let t = j as f64 + (uc - us[j]) / if span != 0.0 { span } else { 1.0 };
                let pts = chains.iter().map(|c| eval_chain(c, t)).collect();
                let ks = k_ys
                    .iter()
                    .zip(k_ms)
                    .map(|(ys, m)| hermite(us, ys, m, uc).clamp(0.0, 1.0))
                    .collect();
                (pts, ks)
            }
        }
    }
}

/// What ends a trimmed section, at its top or its bottom.
#[derive(Debug, Clone, Copy, PartialEq)]
enum End {
    /// The deck edge, untrimmed (top), or the section's own last point
    /// (bottom: an open section).
    Sheet,
    Trim,
    Centreline,
    Transom,
}

/// A section of the hull at u: its frame and its curve.
struct Section {
    p: V2,
    n_hat: V2,
    segs: Vec<Bez>,
    ks: Vec<f64>,
}

impl Section {
    fn vmax(&self) -> f64 {
        self.segs.len() as f64
    }

    fn world(&self, v: f64) -> V3 {
        let [n, z] = eval_chain(&self.segs, v);
        [
            self.p[0] + n * self.n_hat[0],
            self.p[1] + n * self.n_hat[1],
            z,
        ]
    }
}

/// A trimmed half-section: `(S − 1)·R + 1` points from its top to its
/// bottom, the station knots on rows `i·R`.
struct Column {
    pts: Vec<V3>,
    top: End,
    bottom: End,
    /// Rows carrying a knuckle (camber's crease rows).
    creases: Vec<usize>,
}

/// camber's swept hull, read at any u.
pub struct Hull {
    plan: Plan,
    trim_x: Vec<f64>,
    trim_z: Vec<f64>,
    trim_m: Vec<f64>,
    transom: [V2; 2],
    loft: Loft,
}

/// Refine a bracketed sign change of g on [a, b] (camber's `bisectRoot`).
fn bisect_root(g: impl Fn(f64) -> f64, mut a: f64, mut b: f64, ga: f64) -> f64 {
    for _ in 0..40 {
        let m = 0.5 * (a + b);
        if (g(m) < 0.0) == (ga < 0.0) {
            a = m;
        } else {
            b = m;
        }
    }
    0.5 * (a + b)
}

impl Hull {
    pub fn new(doc: &Document) -> Hull {
        let trim_x: Vec<f64> = doc.trim.iter().map(|p| p[0]).collect();
        let trim_z: Vec<f64> = doc.trim.iter().map(|p| p[1]).collect();
        let trim_m = pchip_slopes(&trim_x, &trim_z);
        Hull {
            plan: Plan::new(doc.plan.clone()),
            trim_x,
            trim_z,
            trim_m,
            transom: doc.transom,
            loft: Loft::new(&doc.stations),
        }
    }

    fn trim_z(&self, x: f64) -> f64 {
        let (x0, x1) = (self.trim_x[0], self.trim_x[self.trim_x.len() - 1]);
        hermite(&self.trim_x, &self.trim_z, &self.trim_m, x.clamp(x0, x1))
    }

    fn x_transom(&self, z: f64) -> f64 {
        let [a, b] = self.transom;
        let dz = b[1] - a[1];
        a[0] + (b[0] - a[0]) * ((z - a[1]) / if dz != 0.0 { dz } else { 1.0 })
    }

    /// The three trims at a point: each ≥ 0 where the hull keeps it.
    fn constraints(&self, p: V3) -> [f64; 3] {
        [self.trim_z(p[0]) - p[2], p[1], p[0] - self.x_transom(p[2])]
    }

    fn keep(&self, p: V3) -> f64 {
        let c = self.constraints(p);
        c[0].min(c[1]).min(c[2])
    }

    fn section(&self, u: f64) -> Section {
        let p = self.plan.at(u);
        let d = self.plan.d(u);
        let l = d[0].hypot(d[1]);
        let l = if l != 0.0 { l } else { 1.0 };
        let t = [d[0] / l, d[1] / l];
        let (pts, ks) = self.loft.at(u);
        let segs = cr_chain(&pts, &centripetal_params(&pts), &ks);
        Section {
            p,
            n_hat: [t[1], -t[0]],
            segs,
            ks,
        }
    }

    /// The kept span [vTop, vBot] of a section (camber's `keptSpan`).
    fn kept_span(&self, sec: &Section) -> Option<(f64, f64)> {
        const FN: usize = 96;
        let g = |v: f64| self.keep(sec.world(v));
        let at = |i: usize| sec.vmax() * i as f64 / FN as f64;
        let gs: Vec<f64> = (0..=FN).map(|i| g(at(i))).collect();
        let lo = gs.iter().position(|&x| x >= 0.0)?;
        let hi = gs.iter().rposition(|&x| x >= 0.0)?;
        let v_top = if lo == 0 {
            0.0
        } else {
            bisect_root(g, at(lo - 1), at(lo), gs[lo - 1])
        };
        let v_bot = if hi == FN {
            sec.vmax()
        } else {
            bisect_root(g, at(hi + 1), at(hi), gs[hi + 1])
        };
        (v_bot > v_top + 1e-9).then_some((v_top, v_bot))
    }

    fn exists(&self, u: f64) -> bool {
        self.kept_span(&self.section(u)).is_some()
    }

    /// The first u with a section (camber's `aftLimit`).
    pub fn aft_limit(&self) -> f64 {
        if self.exists(0.0) {
            return 0.0;
        }
        let (mut lo, mut hi) = (0.0, 0.5);
        if !self.exists(hi) {
            return 0.0;
        }
        for _ in 0..24 {
            let m = 0.5 * (lo + hi);
            if self.exists(m) {
                hi = m;
            } else {
                lo = m;
            }
        }
        hi
    }

    /// The last u with a section (camber's `forwardLimit`).
    pub fn forward_limit(&self) -> f64 {
        if self.exists(1.0) {
            return 1.0;
        }
        let (mut lo, mut hi) = (0.5, 1.0);
        if !self.exists(lo) {
            return 1.0;
        }
        for _ in 0..24 {
            let m = 0.5 * (lo + hi);
            if self.exists(m) {
                lo = m;
            } else {
                hi = m;
            }
        }
        lo
    }

    /// Which trim holds a section's end, `v` its parameter there.
    fn end_at(&self, p: V3, at_sheet_end: bool) -> End {
        if at_sheet_end {
            return End::Sheet;
        }
        let c = self.constraints(p);
        let k = (0..3)
            .min_by(|&a, &b| c[a].abs().total_cmp(&c[b].abs()))
            .unwrap();
        [End::Trim, End::Centreline, End::Transom][k]
    }

    /// The trimmed starboard half-section at u (camber's `sweptSection`),
    /// with what ends it; None where there's no hull.
    fn column(&self, u: f64, r: usize) -> Option<Column> {
        let sec = self.section(u);
        let (v_top, v_bot) = self.kept_span(&sec)?;
        let s = sec.segs.len() + 1;
        let margin = (v_bot - v_top) * 1e-3;
        let (anchors, pinned): (Vec<f64>, Vec<bool>) = (0..s)
            .map(|i| {
                let lo = v_top + i as f64 * margin;
                let hi = v_bot - (s - 1 - i) as f64 * margin;
                let a = (i as f64).max(lo).min(hi);
                (a, (a - i as f64).abs() < 1e-9)
            })
            .unzip();
        let mut pts = vec![sec.world(anchors[0])];
        let mut creases = Vec::new();
        for i in 0..s - 1 {
            for k in 1..=r {
                pts.push(
                    sec.world(anchors[i] + (anchors[i + 1] - anchors[i]) * k as f64 / r as f64),
                );
            }
            let j = i + 1;
            if j < s - 1 && pinned[j] && sec.ks[j] > 1e-6 {
                creases.push(j * r);
            }
        }
        let last = pts.len() - 1;
        if pts[last][1].abs() < 1e-6 {
            pts[last][1] = 0.0;
        }
        let top = self.end_at(pts[0], v_top == 0.0);
        let bottom = if pts[last][1] == 0.0 {
            End::Centreline
        } else {
            self.end_at(pts[last], v_bot == sec.vmax())
        };
        Some(Column {
            pts,
            top,
            bottom,
            creases,
        })
    }

    /// The trimmed starboard half-section at u as camber's `sweptSection`
    /// gives it: its points, and whether it ends on the centreline.
    pub fn swept_section(&self, u: f64, r: usize) -> Option<(Vec<V3>, bool)> {
        self.column(u, r).map(|c| {
            let keel = c.bottom == End::Centreline;
            (c.pts, keel)
        })
    }
}

// ---------------------------------------------------------------------------
// Interpolation (camber's step.ts)
// ---------------------------------------------------------------------------

fn basis_funs(span: usize, u: f64, p: usize, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

fn dist3(a: V3, b: V3) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn chord_params(pts: &[V3]) -> Vec<f64> {
    let mut d = vec![0.0];
    let mut total = 0.0;
    for i in 1..pts.len() {
        total += dist3(pts[i], pts[i - 1]).max(1e-9);
        d.push(total);
    }
    let total = if total != 0.0 { total } else { 1.0 };
    d.iter().map(|x| x / total).collect()
}

/// Parameters averaged over a family of point rows.
fn averaged_params(rows: &[Vec<V3>]) -> Vec<f64> {
    let mut acc = vec![0.0; rows[0].len()];
    for r in rows {
        for (a, t) in acc.iter_mut().zip(chord_params(r)) {
            *a += t;
        }
    }
    acc.iter().map(|x| x / rows.len() as f64).collect()
}

/// A clamped knot vector over the data parameters' own range, its interior
/// knots averaged.
fn clamped_knots(ub: &[f64], p: usize) -> Vec<f64> {
    let n = ub.len() - 1;
    let mut k = vec![ub[0]; p + 1];
    for j in 1..=n - p {
        k.push(ub[j..j + p].iter().sum::<f64>() / p as f64);
    }
    k.extend(std::iter::repeat_n(ub[n], p + 1));
    k
}

/// LU (Doolittle, partial pivoting) of a square matrix.
struct Lu {
    lu: Vec<Vec<f64>>,
    piv: Vec<usize>,
}

impl Lu {
    fn new(mut a: Vec<Vec<f64>>) -> Lu {
        let n = a.len();
        let mut piv: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let pr = (k..n)
                .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
                .unwrap();
            if pr != k {
                a.swap(k, pr);
                piv.swap(k, pr);
            }
            let d = if a[k][k] != 0.0 { a[k][k] } else { 1e-12 };
            for i in k + 1..n {
                a[i][k] /= d;
                for j in k + 1..n {
                    a[i][j] -= a[i][k] * a[k][j];
                }
            }
        }
        Lu { lu: a, piv }
    }

    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.lu.len();
        let mut y = vec![0.0; n];
        for i in 0..n {
            let row = &self.lu[i];
            y[i] = b[self.piv[i]] - (0..i).map(|j| row[j] * y[j]).sum::<f64>();
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let row = &self.lu[i];
            let s = y[i] - (i + 1..n).map(|j| row[j] * x[j]).sum::<f64>();
            let d = row[i];
            x[i] = s / if d != 0.0 { d } else { 1e-12 };
        }
        x
    }
}

/// The collocation matrix of a degree-p curve with knots `knots` at `ub`.
fn collocation(ub: &[f64], p: usize, knots: &[f64]) -> Lu {
    let n = ub.len() - 1;
    let mut a = vec![vec![0.0; n + 1]; n + 1];
    for (k, &u) in ub.iter().enumerate() {
        let span = find_span(n, p, u, knots);
        for (t, b) in basis_funs(span, u, p, knots).into_iter().enumerate() {
            a[k][span - p + t] = b;
        }
    }
    Lu::new(a)
}

/// The B-spline surface through a grid `q[i][l]` (i along u, l along v):
/// cubic both ways, interpolated by strips in v split at `creases`, which
/// it joins with a multiplicity-q knot so it can break its tangent there
/// (camber's `interpSurfaceCreased`).
fn interp_surface(q: &[Vec<V3>], creases: &[usize]) -> NurbsSurface3 {
    let nu = q.len() - 1;
    let nv = q[0].len() - 1;
    let (p, qd) = (3.min(nu), 3.min(nv));
    let cols: Vec<Vec<V3>> = (0..=nv)
        .map(|l| q.iter().map(|row| row[l]).collect())
        .collect();
    let ub = averaged_params(&cols);
    let vb = averaged_params(q);
    let knots_u = clamped_knots(&ub, p);
    // Down u for each v-column: intermediate controls r[i][l].
    let au = collocation(&ub, p, &knots_u);
    let mut r = vec![vec![[0.0; 3]; nv + 1]; nu + 1];
    for l in 0..=nv {
        for c in 0..3 {
            let sol = au.solve(&q.iter().map(|row| row[l][c]).collect::<Vec<_>>());
            for i in 0..=nu {
                r[i][l][c] = sol[i];
            }
        }
    }
    // Strips in v between creases, each long enough to interpolate.
    let mut cs: Vec<usize> = creases
        .iter()
        .copied()
        .filter(|&c| c > 0 && c < nv)
        .collect();
    cs.sort_unstable();
    cs.dedup();
    let mut bounds = vec![0];
    for c in cs {
        if c - bounds[bounds.len() - 1] > qd && nv - c > qd {
            bounds.push(c);
        }
    }
    bounds.push(nv);
    let segs: Vec<(usize, usize, Vec<f64>, Lu)> = bounds
        .windows(2)
        .map(|w| {
            let params = &vb[w[0]..=w[1]];
            let knots = clamped_knots(params, qd);
            let lu = collocation(params, qd, &knots);
            (w[0], w[1], knots, lu)
        })
        .collect();
    let mut knots_v: Vec<f64> = Vec::new();
    for (si, (a, _, knots, _)) in segs.iter().enumerate() {
        if si == 0 {
            knots_v = knots.clone();
        } else {
            knots_v.truncate(knots_v.len() - (qd + 1));
            knots_v.extend(std::iter::repeat_n(vb[*a], qd));
            knots_v.extend_from_slice(&knots[qd + 1..]);
        }
    }
    let ncv = knots_v.len() - qd - 1;
    let mut ctrl = vec![[0.0; 3]; (nu + 1) * ncv];
    for i in 0..=nu {
        for c in 0..3 {
            let mut col = 0;
            for (si, (a, b, _, lu)) in segs.iter().enumerate() {
                let rhs: Vec<f64> = (*a..=*b).map(|l| r[i][l][c]).collect();
                let sol = lu.solve(&rhs);
                for x in sol.into_iter().skip(usize::from(si > 0)) {
                    ctrl[i * ncv + col][c] = x;
                    col += 1;
                }
            }
        }
    }
    NurbsSurface3 {
        degree_u: p,
        degree_v: qd,
        knots_u,
        knots_v,
        n_ctrl_u: nu + 1,
        n_ctrl_v: ncv,
        weights: vec![1.0; (nu + 1) * ncv],
        ctrl,
        trim_uv: None,
    }
}

/// The curve through `pts`, cubic (or less, for few points), on chord
/// parameters `ub`: its knots and control points.
fn interp_curve(pts: &[V3], ub: &[f64]) -> (usize, Vec<f64>, Vec<V3>) {
    let n = pts.len() - 1;
    let p = 3.min(n);
    let knots = clamped_knots(ub, p);
    let lu = collocation(ub, p, &knots);
    let mut ctrl = vec![[0.0; 3]; n + 1];
    for c in 0..3 {
        let sol = lu.solve(&pts.iter().map(|q| q[c]).collect::<Vec<_>>());
        for (i, x) in sol.into_iter().enumerate() {
            ctrl[i][c] = x;
        }
    }
    (p, knots, ctrl)
}

/// The ruled surface between a curve through the edge points and one
/// through where each one's ruling meets the centreline: half a transom,
/// or of a plan's open end.
fn ruled_to_centreline(edge: &[(V3, V3)]) -> NurbsSurface3 {
    let outer: Vec<V3> = edge.iter().map(|e| e.0).collect();
    let inner: Vec<V3> = edge.iter().map(|e| e.1).collect();
    let ub = chord_params(&outer);
    let (p, knots_v, mut ctrl) = interp_curve(&outer, &ub);
    ctrl.extend(interp_curve(&inner, &ub).2);
    let n = outer.len();
    NurbsSurface3 {
        degree_u: 1,
        degree_v: p,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v,
        n_ctrl_u: 2,
        n_ctrl_v: n,
        weights: vec![1.0; 2 * n],
        ctrl,
        trim_uv: None,
    }
}

// ---------------------------------------------------------------------------
// Patches
// ---------------------------------------------------------------------------

/// Sub-steps per section segment (girth), and sections along the hull.
const GIRTH_R: usize = 8;
const SECTIONS: usize = 120;

/// A starboard patch and its mirror to port. The solver finds a hull's
/// centreplane by the patches either side of it (one crossing per patch),
/// so no patch spans it.
fn both_sides(stbd: NurbsSurface3) -> [NurbsSurface3; 2] {
    let mut port = stbd.clone();
    for q in port.ctrl.iter_mut() {
        q[1] = -q[1];
    }
    [stbd, port]
}

/// The sections over [u0, u1], closer together toward the ends, where a
/// hull's edges turn hardest (a stem, a forefoot, a transom's corners).
fn columns(hull: &Hull, u0: f64, u1: f64, n: usize) -> Result<Vec<Column>, String> {
    (0..=n)
        .map(|i| {
            let f = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos());
            let u = u0 + (u1 - u0) * f;
            hull.column(u, GIRTH_R)
                .ok_or_else(|| format!("no section at u = {u:.4} along the sheer plan"))
        })
        .collect()
}

/// The u in [a, b] where `is_a` stops holding, given it holds at a and not
/// at b.
fn transition(a: f64, b: f64, is_a: impl Fn(f64) -> bool) -> f64 {
    let (mut lo, mut hi) = (a, b);
    for _ in 0..40 {
        let m = 0.5 * (lo + hi);
        if is_a(m) {
            lo = m;
        } else {
            hi = m;
        }
    }
    hi
}

/// The hull's B-spline patches, in the document's own unit and frame.
///
/// A patch is interpolated through sections that all end the same way, top
/// and bottom: where that changes (the sections' tops leaving the sheer for
/// the centreline at a stem, or the transom for the sheer at its head; their
/// bottoms leaving the transom for the keel at its foot), the hull's edge
/// turns a corner, which a cubic through it would overshoot. So each span
/// between changes, found by bisection, is a patch of its own.
pub fn patches(doc: &Document) -> Result<Vec<NurbsSurface3>, String> {
    let hull = Hull::new(doc);
    let (u_aft, u_fwd) = (hull.aft_limit(), hull.forward_limit());
    if !hull.exists(u_aft) || u_fwd <= u_aft {
        return Err("the hull has no sections: its trims cut it all away".into());
    }
    let span = u_fwd - u_aft;
    let ends = |u: f64| hull.column(u, 1).map(|c| (c.top, c.bottom));

    let mut cuts = vec![u_aft];
    const PROBES: usize = 400;
    for i in 0..PROBES {
        let (a, b) = (
            u_aft + span * i as f64 / PROBES as f64,
            u_aft + span * (i + 1) as f64 / PROBES as f64,
        );
        let at_a = ends(a);
        if ends(b) != at_a {
            cuts.push(transition(a, b, |u| ends(u) == at_a));
        }
    }
    cuts.push(u_fwd);

    let mut out = Vec::new();
    let mut spans: Vec<Vec<Column>> = Vec::new();
    for w in cuts.windows(2) {
        if w[1] - w[0] < 1e-9 * span {
            continue;
        }
        let n = ((SECTIONS as f64 * (w[1] - w[0]) / span).ceil() as usize).max(12);
        let cols = columns(&hull, w[0], w[1], n)?;
        if let Some(c) = cols.iter().find(|c| c.bottom == End::Sheet) {
            let x = c.pts[c.pts.len() - 1][0];
            return Err(format!(
                "the section at x = {x:.4} reaches neither the centreline nor the transom: \
                 an open bottom"
            ));
        }
        let creases: Vec<usize> = cols.iter().flat_map(|c| c.creases.clone()).collect();
        let grid: Vec<Vec<V3>> = cols.iter().map(|c| c.pts.clone()).collect();
        out.extend(both_sides(interp_surface(&grid, &creases)));
        spans.push(cols);
    }

    // The aft end: its starboard edge, from the sheer down to the keel,
    // ruled across to the centreline, and its mirror. That's the first
    // section, and the transom's edge either side of it: the tops of the
    // sections that end on it there (up to its head), then the bottoms (down
    // to its foot).
    let tol = 1e-9 * doc.loa();
    let first = hull
        .column(u_aft, GIRTH_R)
        .expect("the aft limit has a section");
    let leading = |on: fn(&Column) -> End| -> Vec<&Column> {
        spans
            .iter()
            .take_while(|cols| on(&cols[0]) == End::Transom)
            .flat_map(|cols| cols.iter().skip(1))
            .collect()
    };
    // Each edge point is ruled to the centreline in its own face: straight
    // across on the transom (a plane level in y), along the section's plane
    // on the first section (normal to the sheer plan, so a plan's open end
    // is a shallow V).
    let n_hat = hull.section(u_aft).n_hat;
    let across = |p: V3| (p, [p[0], 0.0, p[2]]);
    let along = |p: V3| {
        // p + t·n̂ at y = 0; n̂ points inboard.
        let t = if n_hat[1] < 0.0 {
            -p[1] / n_hat[1]
        } else {
            0.0
        };
        (p, [p[0] + t * n_hat[0], 0.0, p[2]])
    };
    let mut edge: Vec<(V3, V3)> = Vec::new();
    let mut push = |e: (V3, V3)| {
        if edge.last().is_none_or(|q| dist3(q.0, e.0) > tol) {
            edge.push(e);
        }
    };
    for c in leading(|c| c.top).iter().rev() {
        push(across(c.pts[0]));
    }
    first.pts.iter().for_each(|&p| push(along(p)));
    for c in leading(|c| c.bottom) {
        push(across(c.pts[c.pts.len() - 1]));
    }
    // A hull that comes to a point aft has no end to close.
    if edge.len() >= 4 && edge.iter().any(|e| e.0[1] > tol) {
        out.extend(both_sides(ruled_to_centreline(&edge)));
    }
    Ok(out)
}

/// A camber document's hull as boatmath geometry: its patches in metres,
/// floated at the document's deck trim with the design waterline at z = 0.
/// `waterline` [m], if given, is the design waterline's height in the
/// document's frame instead (its deck datum at 0), and `unit_m` a scale to
/// metres in place of the document's unit.
pub fn geometry(
    doc: &Document,
    waterline: Option<f64>,
    unit_m: Option<f64>,
) -> Result<Value, String> {
    let s = unit_m.unwrap_or(doc.unit_m);
    let h = waterline.unwrap_or(-doc.waterline * s);
    let (sin, cos) = doc.deck_trim.sin_cos();
    let mut patches = patches(doc)?;
    for p in patches.iter_mut().flat_map(|p| p.ctrl.iter_mut()) {
        let (x, y, z) = (p[0] * s, p[1] * s, p[2] * s);
        *p = [x * cos - z * sin, y, x * sin + z * cos - h];
    }
    for p in &patches {
        p.validate().map_err(|e| format!("camber hull: {e}"))?;
    }
    Ok(serde_json::json!({
        "kind": "nurbs",
        "hulls": [ { "patches": patches.iter().map(crate::native::patch_value).collect::<Vec<_>>() } ],
    }))
}
