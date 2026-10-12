use std::rc::Rc;

use crate::dim::{Dim, Rational};
use crate::interp::{register_source, Interp, SOURCES};
use crate::ops::{self, q_add, q_div, q_mul, q_sqrt, q_sub, R};
use crate::samples::{self, correlation, ranks, sample_count, Samples};
use crate::units::UnitExpr;
use crate::value::{fmt_num, fmt_value, Env, Quantity, Value};

const BUILTINS: &[&str] = &[
    // distributions
    "normal",
    "uniform",
    "triangular",
    "pert",
    "lognormal",
    "bernoulli",
    // math
    "sqrt",
    "abs",
    "exp",
    "ln",
    "log10",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "atan2",
    "floor",
    "ceil",
    "round",
    "min",
    "max",
    // statistics
    "mean",
    "median",
    "sd",
    "quantile",
    "P",
    "sensitivity",
    // lists and vectors
    "len",
    "sum",
    "map",
    "filter",
    "concat",
    "reverse",
    "range",
    "linspace",
    "dot",
    "cross",
    "norm",
    // geometry
    "area",
    "centroid",
    "section",
    "polyline_length",
    // weights
    "combine",
    "inertia",
    "weight_table",
    // misc
    "print",
];

pub fn install(env: &Env) {
    let mut s = env.borrow_mut();
    for b in BUILTINS {
        s.vars.insert(b.to_string(), Value::Builtin(b));
    }
    let q = |x: f64, u: &str| {
        let ue = parse_unit(u);
        let (scale, dim) = ue.resolve().unwrap();
        Value::Num(Quantity::new(Samples::Exact(x * scale), dim, Some(ue)))
    };
    s.vars.insert(
        "pi".into(),
        Value::Num(Quantity::plain(std::f64::consts::PI)),
    );
    s.vars.insert("g0".into(), q(9.80665, "m/s^2"));
    s.vars.insert("rho_sw".into(), q(1.025, "t/m^3"));
    s.vars.insert("rho_fw".into(), q(1.000, "t/m^3"));
    s.vars.insert("rho_steel".into(), q(7.85, "t/m^3"));
    s.vars.insert("rho_al".into(), q(2.70, "t/m^3"));
}

/// Parse a simple unit string like "t/m^3" or "t·m/cm" (used for builtin constants and outputs).
pub fn parse_unit(s: &str) -> UnitExpr {
    let mut u = UnitExpr::default();
    let mut sign = 1;
    for part in s.split_inclusive(['/', '·']) {
        let (body, next) = match part.chars().last() {
            Some('/') => (&part[..part.len() - 1], -1),
            Some('·') => (&part[..part.len() - '·'.len_utf8()], 1),
            _ => (part, sign),
        };
        let (name, exp) = match body.split_once('^') {
            Some((n, e)) => (n, e.parse::<i32>().unwrap()),
            None => (body, 1),
        };
        u = u.mul(&UnitExpr(vec![(
            name.to_string(),
            Rational::int(exp * sign),
        )]));
        sign = next;
    }
    u
}

pub fn q_unit(s: Samples, unit: &str) -> Value {
    let u = parse_unit(unit);
    let (_, dim) = u.resolve().unwrap();
    Value::Num(Quantity::new(s, dim, Some(u)))
}

// ---------------------------------------------------------------- argument helpers

pub fn arity(args: &[Value], lo: usize, hi: usize) -> R<()> {
    if args.len() < lo || args.len() > hi {
        if lo == hi {
            return Err(format!("expected {} argument(s), got {}", lo, args.len()));
        }
        return Err(format!(
            "expected {} to {} arguments, got {}",
            lo,
            hi,
            args.len()
        ));
    }
    Ok(())
}

pub fn num(v: &Value) -> R<Quantity> {
    match v {
        Value::Num(q) => Ok(q.clone()),
        _ => Err(format!("expected a number, got {}", v.type_name())),
    }
}

pub fn list(v: &Value) -> R<Rc<Vec<Value>>> {
    match v {
        Value::List(xs) => Ok(xs.clone()),
        _ => Err(format!("expected a list, got {}", v.type_name())),
    }
}

pub fn exact_f(v: &Value) -> R<f64> {
    num(v)?
        .s
        .exact()
        .ok_or_else(|| "expected an exact number".to_string())
}

pub fn same_dims(qs: &[&Quantity]) -> R<()> {
    for q in qs {
        if q.dim != qs[0].dim {
            return Err(format!(
                "arguments must share a dimension ({} vs {})",
                qs[0].dim.describe(),
                q.dim.describe()
            ));
        }
    }
    Ok(())
}

pub fn dimensionless(q: &Quantity, what: &str) -> R<()> {
    if !q.dim.is_none() {
        return Err(format!(
            "{} must be dimensionless, got {}",
            what,
            q.dim.describe()
        ));
    }
    Ok(())
}

/// Map a scalar function over a number, a list, or each sample.
fn map_num(v: &Value, f: &dyn Fn(&Quantity) -> R<Quantity>) -> R<Value> {
    match v {
        Value::Num(q) => Ok(Value::Num(f(q)?)),
        Value::List(xs) => Ok(Value::list(
            xs.iter().map(|x| map_num(x, f)).collect::<R<_>>()?,
        )),
        _ => Err(format!("expected a number, got {}", v.type_name())),
    }
}

fn plain_fn(q: &Quantity, f: fn(f64) -> f64, name: &str) -> R<Quantity> {
    dimensionless(q, &format!("the argument of {}", name))?;
    Ok(Quantity::new(
        q.s.map(f),
        Dim::NONE,
        Some(UnitExpr::default()),
    ))
}

/// Points as a list of [x, y] or [x, y, z] coordinate lists sharing one dimension.
fn points(v: &Value) -> R<(Vec<Vec<Samples>>, Quantity)> {
    let pts = list(v)?;
    if pts.len() < 3 {
        return Err("a polygon needs at least 3 points".into());
    }
    let mut out = vec![];
    let mut proto: Option<Quantity> = None;
    for p in pts.iter() {
        let cs = list(p)?;
        let mut row = vec![];
        for c in cs.iter() {
            let q = num(c)?;
            if let Some(pr) = &proto {
                same_dims(&[pr, &q])?;
            } else {
                proto = Some(q.clone());
            }
            row.push(q.s.clone());
        }
        if row.len()
            != out
                .first()
                .map(|r: &Vec<Samples>| r.len())
                .unwrap_or(row.len())
        {
            return Err("all points must have the same number of coordinates".into());
        }
        out.push(row);
    }
    let d = out[0].len();
    if d != 2 && d != 3 {
        return Err("points must be 2D [x, y] or 3D [x, y, z]".into());
    }
    Ok((out, proto.unwrap()))
}

/// Evaluate `f` per sample over all `inputs`, producing `nout` sample vectors.
pub fn per_sample(
    inputs: &[&Samples],
    nout: usize,
    mut f: impl FnMut(usize) -> R<Vec<f64>>,
) -> R<Vec<Samples>> {
    if inputs.iter().all(|s| s.is_exact()) {
        return Ok(f(0)?.into_iter().map(Samples::Exact).collect());
    }
    let n = sample_count();
    let mut outs = vec![Vec::with_capacity(n); nout];
    for i in 0..n {
        for (o, x) in outs.iter_mut().zip(f(i)?) {
            o.push(x);
        }
    }
    Ok(outs
        .into_iter()
        .map(|v| Samples::Dist(Rc::new(v)))
        .collect())
}

/// Quantity with the dimension of `proto` raised to `p`, in proto's units.
fn like(proto: &Quantity, p: i32, s: Samples) -> Value {
    let r = Rational::int(p);
    Value::Num(Quantity::new(
        s,
        proto.dim.pow(r),
        proto.unit.as_ref().map(|u| u.pow(r)),
    ))
}

// ---------------------------------------------------------------- distributions

/// 90% interval [a, b]: lognormal when both ends are positive, otherwise normal.
pub fn ci90(a: &Quantity, b: &Quantity) -> R<Quantity> {
    same_dims(&[a, b])?;
    let z95 = samples::norm_ppf(0.95);
    let (a0, b0) = (a.s.mean(), b.s.mean());
    if b0 < a0 {
        return Err("the lower bound of 'to' must not exceed the upper bound".into());
    }
    let positive = a0 > 0.0;
    let z = samples::draw(|r| r.normal());
    let s = Samples::from_fn(&[&z], |i| {
        let (lo, hi) = (a.s.get(i), b.s.get(i));
        if positive {
            let mu = (lo.ln() + hi.ln()) / 2.0;
            let sigma = (hi.ln() - lo.ln()) / (2.0 * z95);
            (mu + sigma * z.get(i)).exp()
        } else {
            (lo + hi) / 2.0 + (hi - lo) / (2.0 * z95) * z.get(i)
        }
    });
    register_source(&s);
    Ok(Quantity::new(s, a.dim, a.unit.clone().or(b.unit.clone())))
}

fn dist(name: &str, args: &[Value]) -> R<Value> {
    let qs = args.iter().map(num).collect::<R<Vec<_>>>()?;
    let s = match name {
        "normal" => {
            arity(args, 2, 2)?;
            same_dims(&[&qs[0], &qs[1]])?;
            let (m, sd) = (qs[0].s.clone(), qs[1].s.clone());
            let z = samples::draw(|r| r.normal());
            Samples::from_fn(&[&z], |i| m.get(i) + sd.get(i) * z.get(i))
        }
        "uniform" => {
            arity(args, 2, 2)?;
            same_dims(&[&qs[0], &qs[1]])?;
            let (a, b) = (qs[0].s.clone(), qs[1].s.clone());
            let u = samples::draw(|r| r.uniform());
            Samples::from_fn(&[&u], |i| a.get(i) + (b.get(i) - a.get(i)) * u.get(i))
        }
        "triangular" => {
            arity(args, 3, 3)?;
            same_dims(&[&qs[0], &qs[1], &qs[2]])?;
            let (a, c, b) = (qs[0].s.clone(), qs[1].s.clone(), qs[2].s.clone());
            let u = samples::draw(|r| r.uniform());
            Samples::from_fn(&[&u], |i| {
                let (a, c, b, u) = (a.get(i), c.get(i), b.get(i), u.get(i));
                let fc = (c - a) / (b - a);
                if u < fc {
                    a + (u * (b - a) * (c - a)).sqrt()
                } else {
                    b - ((1.0 - u) * (b - a) * (b - c)).sqrt()
                }
            })
        }
        "pert" => {
            arity(args, 3, 3)?;
            same_dims(&[&qs[0], &qs[1], &qs[2]])?;
            let (a, m, b) = (qs[0].s.mean(), qs[1].s.mean(), qs[2].s.mean());
            if !(a <= m && m <= b && a < b) {
                return Err("pert(min, mode, max) needs min ≤ mode ≤ max".into());
            }
            let al = 1.0 + 4.0 * (m - a) / (b - a);
            let be = 1.0 + 4.0 * (b - m) / (b - a);
            samples::draw(|r| a + r.beta(al, be) * (b - a))
        }
        "lognormal" => {
            arity(args, 2, 2)?;
            return Ok(Value::Num(ci90(&qs[0], &qs[1])?));
        }
        "bernoulli" => {
            arity(args, 1, 1)?;
            dimensionless(&qs[0], "probability")?;
            let p = qs[0].s.clone();
            let u = samples::draw(|r| r.uniform());
            let s = Samples::from_fn(&[&u], |i| (u.get(i) < p.get(i)) as u8 as f64);
            register_source(&s);
            return Ok(Value::Bool(s));
        }
        _ => unreachable!(),
    };
    register_source(&s);
    Ok(Value::Num(Quantity::new(s, qs[0].dim, qs[0].unit.clone())))
}

// ---------------------------------------------------------------- statistics

fn stat(name: &str, args: &[Value]) -> R<Value> {
    let q = num(&args[0])?;
    let x = match name {
        "mean" => q.s.mean(),
        "median" => q.s.quantile(0.5),
        "sd" => q.s.sd(),
        "quantile" => {
            arity(args, 2, 2)?;
            let p = exact_f(&args[1])?;
            q.s.quantile(p)
        }
        _ => unreachable!(),
    };
    Ok(Value::Num(Quantity::new(
        Samples::Exact(x),
        q.dim,
        q.unit.clone(),
    )))
}

fn sensitivity(q: &Quantity) -> R<Value> {
    let Samples::Dist(out) = &q.s else {
        return Err("value is exact; nothing to attribute".into());
    };
    let ro = ranks(out);
    let mut rows: Vec<(String, f64, f64)> = vec![];
    SOURCES.with(|src| {
        let mut seen = std::collections::HashMap::new();
        for s in src.borrow().iter() {
            let r = if Rc::ptr_eq(&s.samples, out) {
                1.0
            } else {
                correlation(&ro, &ranks(&s.samples))
            };
            let n = seen.entry(s.label.clone()).or_insert(0);
            *n += 1;
            let label = if *n > 1 {
                format!("{}#{}", s.label, n)
            } else {
                s.label.clone()
            };
            rows.push((label, r, r * r));
        }
    });
    let total: f64 = rows.iter().map(|r| r.2).sum();
    rows.retain(|r| r.2 / total > 0.005);
    rows.sort_by(|a, b| b.2.total_cmp(&a.2));
    rows.truncate(15);
    let fields = rows
        .into_iter()
        .map(|(l, r, r2)| {
            let pct = Quantity::new(
                Samples::Exact(r2 / total),
                Dim::NONE,
                Some(UnitExpr::single("%")),
            );
            let dir = if r >= 0.0 { "↑" } else { "↓" };
            (format!("{} {}", dir, l), Value::Num(pct))
        })
        .collect();
    Ok(Value::Record(Rc::new(fields)))
}

// ---------------------------------------------------------------- lists and vectors

fn sum(xs: &[Value]) -> R<Value> {
    let mut it = xs.iter();
    let mut acc = it.next().ok_or("cannot sum an empty list")?.clone();
    for x in it {
        acc = ops::binop("+", &acc, x)?;
    }
    Ok(acc)
}

pub fn vec3(v: &Value, what: &str) -> R<Vec<Quantity>> {
    let xs = list(v)?;
    if xs.len() != 3 {
        return Err(format!(
            "{} must be a 3-vector, got {} components",
            what,
            xs.len()
        ));
    }
    xs.iter().map(num).collect()
}

fn dot(a: &Value, b: &Value) -> R<Value> {
    let (a, b) = (list(a)?, list(b)?);
    if a.len() != b.len() {
        return Err("dot: vectors must have the same length".into());
    }
    let prods = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| ops::binop("*", x, y))
        .collect::<R<Vec<_>>>()?;
    sum(&prods)
}

fn cross(a: &Value, b: &Value) -> R<Value> {
    let (a, b) = (vec3(a, "cross argument")?, vec3(b, "cross argument")?);
    let c = |i: usize, j: usize| q_sub(&q_mul(&a[i], &b[j]), &q_mul(&a[j], &b[i]));
    Ok(Value::list(vec![
        Value::Num(c(1, 2)?),
        Value::Num(c(2, 0)?),
        Value::Num(c(0, 1)?),
    ]))
}

fn norm(a: &Value) -> R<Quantity> {
    let d = dot(a, a)?;
    Ok(q_sqrt(&num(&d)?))
}

// ---------------------------------------------------------------- geometry

/// Area, centroid and second moments of a simple 2D polygon (shoelace formulas).
/// Returns [A, Cx, Cy, Ixx0, Iyy0, Ixy0] with moments about the coordinate axes.
fn poly2(p: &[(f64, f64)]) -> [f64; 6] {
    let n = p.len();
    let mut r = [0.0; 6];
    for i in 0..n {
        let (x0, y0) = p[i];
        let (x1, y1) = p[(i + 1) % n];
        let a = x0 * y1 - x1 * y0;
        r[0] += a / 2.0;
        r[1] += (x0 + x1) * a / 6.0;
        r[2] += (y0 + y1) * a / 6.0;
        r[3] += (y0 * y0 + y0 * y1 + y1 * y1) * a / 12.0;
        r[4] += (x0 * x0 + x0 * x1 + x1 * x1) * a / 12.0;
        r[5] += (x0 * y1 + 2.0 * x0 * y0 + 2.0 * x1 * y1 + x1 * y0) * a / 24.0;
    }
    let sign = r[0].signum();
    r[1] /= r[0];
    r[2] /= r[0];
    r[0] *= sign;
    r[3] *= sign;
    r[4] *= sign;
    r[5] *= sign;
    r
}

/// Area and centroid of a planar polygon in 3D: [A, Cx, Cy, Cz].
fn poly3(p: &[[f64; 3]]) -> [f64; 4] {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cr = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let mut va = [0.0; 3];
    for i in 1..p.len() - 1 {
        let c = cr(sub(p[i], p[0]), sub(p[i + 1], p[0]));
        for k in 0..3 {
            va[k] += c[k] / 2.0;
        }
    }
    let area = (va[0] * va[0] + va[1] * va[1] + va[2] * va[2]).sqrt();
    let nrm = [va[0] / area, va[1] / area, va[2] / area];
    let mut c = [0.0; 3];
    for i in 1..p.len() - 1 {
        let t = cr(sub(p[i], p[0]), sub(p[i + 1], p[0]));
        let a = (t[0] * nrm[0] + t[1] * nrm[1] + t[2] * nrm[2]) / 2.0;
        for k in 0..3 {
            c[k] += a * (p[0][k] + p[i][k] + p[i + 1][k]) / 3.0;
        }
    }
    [area, c[0] / area, c[1] / area, c[2] / area]
}

fn geometry(name: &str, args: &[Value]) -> R<Value> {
    arity(args, 1, 1)?;
    let (pts, proto) = points(&args[0])?;
    let inputs: Vec<&Samples> = pts.iter().flatten().collect();
    if pts[0].len() == 2 {
        let r = per_sample(&inputs, 6, |i| {
            let p: Vec<(f64, f64)> = pts.iter().map(|c| (c[0].get(i), c[1].get(i))).collect();
            Ok(poly2(&p).to_vec())
        })?;
        let [a, cx, cy, ixx0, iyy0, ixy0] = <[Samples; 6]>::try_from(r).unwrap();
        let point = || {
            Value::list(vec![
                like(&proto, 1, cx.clone()),
                like(&proto, 1, cy.clone()),
            ])
        };
        return Ok(match name {
            "area" => like(&proto, 2, a),
            "centroid" => point(),
            "section" => {
                let ixx = Samples::from_fn(&[&a, &cy, &ixx0], |i| {
                    ixx0.get(i) - a.get(i) * cy.get(i).powi(2)
                });
                let iyy = Samples::from_fn(&[&a, &cx, &iyy0], |i| {
                    iyy0.get(i) - a.get(i) * cx.get(i).powi(2)
                });
                let ixy = Samples::from_fn(&[&a, &cx, &cy, &ixy0], |i| {
                    ixy0.get(i) - a.get(i) * cx.get(i) * cy.get(i)
                });
                Value::record(vec![
                    ("area", like(&proto, 2, a.clone())),
                    ("centroid", point()),
                    ("ixx", like(&proto, 4, ixx)),
                    ("iyy", like(&proto, 4, iyy)),
                    ("ixy", like(&proto, 4, ixy)),
                    ("ixx_axis", like(&proto, 4, ixx0)),
                    ("iyy_axis", like(&proto, 4, iyy0)),
                ])
            }
            _ => unreachable!(),
        });
    }
    if name == "section" {
        return Err("section properties need a 2D polygon".into());
    }
    let r = per_sample(&inputs, 4, |i| {
        let p: Vec<[f64; 3]> = pts
            .iter()
            .map(|c| [c[0].get(i), c[1].get(i), c[2].get(i)])
            .collect();
        Ok(poly3(&p).to_vec())
    })?;
    Ok(match name {
        "area" => like(&proto, 2, r[0].clone()),
        _ => Value::list(r[1..].iter().map(|s| like(&proto, 1, s.clone())).collect()),
    })
}

fn polyline_length(v: &Value) -> R<Value> {
    let pts = list(v)?;
    let mut segs = vec![];
    for w in pts.windows(2) {
        let d = ops::binop("-", &w[1], &w[0])?;
        segs.push(Value::Num(norm(&d)?));
    }
    sum(&segs)
}

// ---------------------------------------------------------------- weights

pub fn item_parts(item: &Value) -> R<(Quantity, Vec<Quantity>)> {
    let m = item.field("mass").ok_or("each item needs a 'mass' field")?;
    let m = num(m)?;
    if m.dim != Dim::mass() {
        return Err(format!(
            "item mass must be a mass, got {}",
            m.dim.describe()
        ));
    }
    let cg = vec3(
        item.field("cg")
            .ok_or("each item needs a 'cg' field [x, y, z]")?,
        "cg",
    )?;
    for c in &cg {
        if c.dim != Dim::length() {
            return Err(format!(
                "cg components must be lengths, got {}",
                c.dim.describe()
            ));
        }
    }
    Ok((m, cg))
}

/// Total mass and centre of gravity of a list of items `{mass, cg}`.
pub fn combine(items: &[Value], name: &str) -> R<Value> {
    if items.is_empty() {
        return Err("no items to combine".into());
    }
    let mut total: Option<Quantity> = None;
    let mut mom: Vec<Option<Quantity>> = vec![None, None, None];
    for it in items {
        let (m, cg) = item_parts(it)?;
        total = Some(match total {
            None => m.clone(),
            Some(t) => q_add(&t, &m)?,
        });
        for k in 0..3 {
            let mk = q_mul(&m, &cg[k]);
            mom[k] = Some(match mom[k].take() {
                None => mk,
                Some(acc) => q_add(&acc, &mk)?,
            });
        }
    }
    let total = total.unwrap();
    let cg: Vec<Value> = mom
        .into_iter()
        .map(|mk| {
            let mut c = q_div(&mk.unwrap(), &total);
            if c.unit.is_none()
                || c.unit
                    .as_ref()
                    .is_some_and(|u| u.resolve().map(|r| r.1) != Ok(Dim::length()))
            {
                c.unit = Some(UnitExpr::single("m"));
            }
            Value::Num(c)
        })
        .collect();
    Ok(Value::record(vec![
        ("name", Value::Str(name.into())),
        ("mass", Value::Num(total)),
        ("cg", Value::list(cg)),
    ]))
}

/// Mass moments of inertia of point-mass items about their combined CG, and radii of gyration.
fn inertia(items: &[Value]) -> R<Value> {
    let c = combine(items, "")?;
    let mt = num(c.field("mass").unwrap())?;
    let cg = vec3(c.field("cg").unwrap(), "cg")?;
    let mut acc: [Option<Quantity>; 3] = [None, None, None];
    for it in items {
        let (m, p) = item_parts(it)?;
        let d: Vec<Quantity> = (0..3).map(|k| q_sub(&p[k], &cg[k])).collect::<R<_>>()?;
        let sq: Vec<Quantity> = d.iter().map(|x| q_mul(x, x)).collect();
        let terms = [
            q_add(&sq[1], &sq[2])?,
            q_add(&sq[0], &sq[2])?,
            q_add(&sq[0], &sq[1])?,
        ];
        for k in 0..3 {
            let t = q_mul(&m, &terms[k]);
            acc[k] = Some(match acc[k].take() {
                None => t,
                Some(a) => q_add(&a, &t)?,
            });
        }
    }
    let unit = |s: &Samples| q_unit(s.clone(), "t·m^2");
    let mut fields = vec![];
    let names = [("ixx", "kxx"), ("iyy", "kyy"), ("izz", "kzz")];
    let mut ks = vec![];
    for k in 0..3 {
        let i = acc[k].take().unwrap();
        fields.push((names[k].0, unit(&i.s)));
        let mut g = q_sqrt(&q_div(&i, &mt));
        g.unit = Some(UnitExpr::single("m"));
        ks.push((names[k].1, Value::Num(g)));
    }
    fields.extend(ks);
    Ok(Value::record(fields))
}

fn weight_table(it: &mut Interp, items: &[Value]) -> R<Value> {
    let total = combine(items, "TOTAL")?;
    let mut rows: Vec<Value> = items.to_vec();
    rows.push(total);
    let cell = |q: &Quantity| {
        let (_, scale) = q.display_unit();
        let m = q.s.quantile(0.5) / scale;
        if q.s.is_exact() {
            fmt_num(m)
        } else {
            format!("{} ±{}", fmt_num(m), fmt_num(q.s.sd() / scale))
        }
    };
    let mut lines = vec![];
    for (i, r) in rows.iter().enumerate() {
        let (m, cg) = item_parts(r)?;
        let name = match r.field("name") {
            Some(Value::Str(s)) => s.to_string(),
            _ => format!("item {}", i + 1),
        };
        let mut mm = m.clone();
        if i == rows.len() - 1 {
            mm.unit = num(items[0].field("mass").unwrap())?.unit;
        }
        lines.push([name, cell(&mm), cell(&cg[0]), cell(&cg[1]), cell(&cg[2])]);
    }
    let (mu, _) = num(items[0].field("mass").unwrap())?.display_unit();
    let (lu, _) = num(&list(rows[0].field("cg").unwrap())?[0])?.display_unit();
    let header = [
        "item".to_string(),
        format!("mass ({})", mu),
        format!("LCG ({})", lu),
        format!("TCG ({})", lu),
        format!("VCG ({})", lu),
    ];
    let mut w = [0; 5];
    for row in std::iter::once(&header).chain(lines.iter()) {
        for k in 0..5 {
            w[k] = w[k].max(row[k].chars().count());
        }
    }
    let fmt_row = |r: &[String; 5]| {
        let mut s = format!("{:<w$}", r[0], w = w[0]);
        for k in 1..5 {
            s += &format!("  {:>w$}", r[k], w = w[k]);
        }
        s
    };
    let total_w: usize = w.iter().sum::<usize>() + 8;
    let mut out = vec![fmt_row(&header), "─".repeat(total_w)];
    for (i, l) in lines.iter().enumerate() {
        if i == lines.len() - 1 {
            out.push("─".repeat(total_w));
        }
        out.push(fmt_row(l));
    }
    it.emit(&out.join("\n"));
    Ok(Value::Unit)
}

// ---------------------------------------------------------------- dispatch

pub fn call(it: &mut Interp, name: &str, args: Vec<Value>) -> R<Value> {
    match name {
        "normal" | "uniform" | "triangular" | "pert" | "lognormal" | "bernoulli" => {
            dist(name, &args)
        }

        "sqrt" => {
            arity(&args, 1, 1)?;
            map_num(&args[0], &|q| Ok(q_sqrt(q)))
        }
        "abs" => {
            arity(&args, 1, 1)?;
            map_num(&args[0], &|q| {
                Ok(Quantity::new(q.s.map(f64::abs), q.dim, q.unit.clone()))
            })
        }
        "floor" | "ceil" | "round" => {
            arity(&args, 1, 1)?;
            let f: fn(f64) -> f64 = match name {
                "floor" => f64::floor,
                "ceil" => f64::ceil,
                _ => f64::round,
            };
            map_num(&args[0], &|q| plain_fn(q, f, name))
        }
        "exp" | "ln" | "log10" | "sin" | "cos" | "tan" | "asin" | "acos" | "atan" => {
            arity(&args, 1, 1)?;
            let f: fn(f64) -> f64 = match name {
                "exp" => f64::exp,
                "ln" => f64::ln,
                "log10" => f64::log10,
                "sin" => f64::sin,
                "cos" => f64::cos,
                "tan" => f64::tan,
                "asin" => f64::asin,
                "acos" => f64::acos,
                _ => f64::atan,
            };
            map_num(&args[0], &|q| plain_fn(q, f, name))
        }
        "atan2" => {
            arity(&args, 2, 2)?;
            let (y, x) = (num(&args[0])?, num(&args[1])?);
            same_dims(&[&y, &x])?;
            Ok(Value::Num(Quantity::new(
                y.s.zip(&x.s, f64::atan2),
                Dim::NONE,
                Some(UnitExpr::default()),
            )))
        }
        "min" | "max" => {
            let vals: Vec<Value> = if args.len() == 1 {
                list(&args[0])?.to_vec()
            } else {
                args
            };
            let qs = vals.iter().map(num).collect::<R<Vec<_>>>()?;
            if qs.is_empty() {
                return Err("needs at least one value".into());
            }
            same_dims(&qs.iter().collect::<Vec<_>>())?;
            let is_min = name == "min";
            let mut acc = qs[0].s.clone();
            for q in &qs[1..] {
                acc = acc.zip(&q.s, |a, b| if is_min { a.min(b) } else { a.max(b) });
            }
            Ok(Value::Num(Quantity::new(
                acc,
                qs[0].dim,
                qs[0].unit.clone(),
            )))
        }

        "mean" | "median" | "sd" | "quantile" => {
            arity(&args, 1, 2)?;
            stat(name, &args)
        }
        "P" => {
            arity(&args, 1, 1)?;
            match &args[0] {
                Value::Bool(s) => Ok(Value::Num(Quantity::plain(s.mean()))),
                v => Err(format!("expected a condition, got {}", v.type_name())),
            }
        }
        "sensitivity" => {
            arity(&args, 1, 1)?;
            sensitivity(&num(&args[0])?)
        }

        "len" => {
            arity(&args, 1, 1)?;
            Ok(Value::Num(Quantity::plain(list(&args[0])?.len() as f64)))
        }
        "sum" => {
            arity(&args, 1, 1)?;
            sum(&list(&args[0])?)
        }
        "map" => {
            arity(&args, 2, 2)?;
            let xs = list(&args[0])?;
            let out = xs
                .iter()
                .map(|x| it.call(&args[1], vec![x.clone()]))
                .collect::<R<Vec<_>>>()?;
            Ok(Value::list(out))
        }
        "filter" => {
            arity(&args, 2, 2)?;
            let xs = list(&args[0])?;
            let mut out = vec![];
            for x in xs.iter() {
                match it.call(&args[1], vec![x.clone()])? {
                    Value::Bool(Samples::Exact(b)) => {
                        if b != 0.0 {
                            out.push(x.clone())
                        }
                    }
                    Value::Bool(_) => {
                        return Err("filter condition must be certain, not uncertain".into())
                    }
                    v => {
                        return Err(format!(
                            "filter condition must be a bool, got {}",
                            v.type_name()
                        ))
                    }
                }
            }
            Ok(Value::list(out))
        }
        "concat" => {
            let mut out = vec![];
            for a in &args {
                out.extend(list(a)?.iter().cloned());
            }
            Ok(Value::list(out))
        }
        "reverse" => {
            arity(&args, 1, 1)?;
            Ok(Value::list(list(&args[0])?.iter().rev().cloned().collect()))
        }
        "range" => {
            arity(&args, 1, 2)?;
            let (a, b) = if args.len() == 1 {
                (0.0, exact_f(&args[0])?)
            } else {
                (exact_f(&args[0])?, exact_f(&args[1])?)
            };
            Ok(Value::list(
                (a as i64..b as i64)
                    .map(|i| Value::Num(Quantity::plain(i as f64)))
                    .collect(),
            ))
        }
        "linspace" => {
            arity(&args, 3, 3)?;
            let (a, b) = (num(&args[0])?, num(&args[1])?);
            same_dims(&[&a, &b])?;
            let n = exact_f(&args[2])? as usize;
            if n < 2 {
                return Err("needs at least 2 points".into());
            }
            Ok(Value::list(
                (0..n)
                    .map(|i| {
                        let f = i as f64 / (n - 1) as f64;
                        Value::Num(Quantity::new(
                            a.s.zip(&b.s, |x, y| x + f * (y - x)),
                            a.dim,
                            a.unit.clone(),
                        ))
                    })
                    .collect(),
            ))
        }
        "dot" => {
            arity(&args, 2, 2)?;
            dot(&args[0], &args[1])
        }
        "cross" => {
            arity(&args, 2, 2)?;
            cross(&args[0], &args[1])
        }
        "norm" => {
            arity(&args, 1, 1)?;
            Ok(Value::Num(norm(&args[0])?))
        }

        "area" | "centroid" | "section" => geometry(name, &args),
        "polyline_length" => {
            arity(&args, 1, 1)?;
            polyline_length(&args[0])
        }

        "combine" => {
            arity(&args, 1, 2)?;
            let label = match args.get(1) {
                Some(Value::Str(s)) => s.to_string(),
                Some(v) => return Err(format!("name must be a string, got {}", v.type_name())),
                None => "total".into(),
            };
            combine(&list(&args[0])?, &label)
        }
        "inertia" => {
            arity(&args, 1, 1)?;
            inertia(&list(&args[0])?)
        }
        "weight_table" => {
            arity(&args, 1, 1)?;
            weight_table(it, &list(&args[0])?)
        }

        "print" => {
            let parts: Vec<String> = args
                .iter()
                .map(|a| match a {
                    Value::Str(s) => s.to_string(),
                    v => fmt_value(v, 0),
                })
                .collect();
            it.emit(&parts.join(" "));
            Ok(Value::Unit)
        }
        _ => Err(format!("unknown builtin {}", name)),
    }
}
