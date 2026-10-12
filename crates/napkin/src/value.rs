use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::ast::Expr;
use crate::dim::{Dim, Rational};
use crate::samples::Samples;
use crate::units::{si_unit_for, UnitExpr};

/// A number with a dimension, a preferred display unit, and its samples (in SI).
#[derive(Clone, Debug)]
pub struct Quantity {
    pub dim: Dim,
    /// Display hint: the unit expression this value was written in. `None` means
    /// "display in SI"; `Some(empty)` is a plain number.
    pub unit: Option<UnitExpr>,
    pub s: Samples,
}

impl Quantity {
    pub fn plain(x: f64) -> Quantity {
        Quantity {
            dim: Dim::NONE,
            unit: Some(UnitExpr::default()),
            s: Samples::Exact(x),
        }
    }

    pub fn new(s: Samples, dim: Dim, unit: Option<UnitExpr>) -> Quantity {
        Quantity { dim, unit, s }
    }

    /// The unit used for display, and its scale relative to SI.
    pub fn display_unit(&self) -> (UnitExpr, f64) {
        if let Some(u) = &self.unit {
            if let Ok((scale, dim)) = u.resolve() {
                if dim == self.dim {
                    return (u.clone(), scale);
                }
            }
        }
        let u = si_unit_for(&self.dim);
        let scale = u.resolve().map(|(s, _)| s).unwrap_or(1.0);
        (u, scale)
    }
}

/// A value owned by a host extension (a hull, say): opaque to the language,
/// handed back to the extension's native functions.
pub trait Native: std::any::Any {
    fn type_name(&self) -> &'static str;
    fn describe(&self) -> String;
    fn as_any(&self) -> &dyn std::any::Any;
}

pub type NativeFnImpl = dyn Fn(&mut crate::interp::Interp, Vec<Value>) -> Result<Value, String>;

/// A function provided by a host extension.
pub struct NativeFn {
    pub name: String,
    pub f: Box<NativeFnImpl>,
}

pub struct Closure {
    pub name: String,
    pub params: Vec<String>,
    pub body: Expr,
    pub env: Env,
}

#[derive(Clone)]
pub enum Value {
    Num(Quantity),
    Bool(Samples),
    Str(Rc<str>),
    List(Rc<Vec<Value>>),
    Record(Rc<Vec<(String, Value)>>),
    Func(Rc<Closure>),
    Builtin(&'static str),
    NativeFn(Rc<NativeFn>),
    Native(Rc<dyn Native>),
    Unit,
}

impl Value {
    pub fn type_name(&self) -> &'static str {
        match self {
            Value::Num(q) if q.dim.is_none() => "number",
            Value::Num(_) => "quantity",
            Value::Bool(_) => "bool",
            Value::Str(_) => "string",
            Value::List(_) => "list",
            Value::Record(_) => "record",
            Value::Func(_) | Value::Builtin(_) | Value::NativeFn(_) => "function",
            Value::Native(n) => n.type_name(),
            Value::Unit => "nothing",
        }
    }

    pub fn field(&self, name: &str) -> Option<&Value> {
        match self {
            Value::Record(fs) => fs.iter().find(|(n, _)| n == name).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn record(fields: Vec<(&str, Value)>) -> Value {
        Value::Record(Rc::new(
            fields
                .into_iter()
                .map(|(n, v)| (n.to_string(), v))
                .collect(),
        ))
    }

    pub fn list(items: Vec<Value>) -> Value {
        Value::List(Rc::new(items))
    }
}

#[derive(Default)]
pub struct Scope {
    pub vars: HashMap<String, Value>,
    pub parent: Option<Env>,
}

pub type Env = Rc<RefCell<Scope>>;

pub fn lookup(env: &Env, name: &str) -> Option<Value> {
    let s = env.borrow();
    match s.vars.get(name) {
        Some(v) => Some(v.clone()),
        None => s.parent.as_ref().and_then(|p| lookup(p, name)),
    }
}

// ---------------------------------------------------------------- display

pub fn fmt_num(x: f64) -> String {
    // Round-off residue (e.g. the TCG of a symmetric shape) prints as zero.
    if x.abs() < 1e-12 {
        return "0".into();
    }
    if !x.is_finite() {
        return format!("{}", x);
    }
    let a = x.abs();
    if !(1e-4..1e7).contains(&a) {
        return format!("{:.3e}", x);
    }
    let digits = 4i32;
    let mag = a.log10().floor() as i32;
    let decimals = (digits - 1 - mag).max(0) as usize;
    let s = format!("{:.*}", decimals, x);
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

fn with_unit(x: f64, u: &UnitExpr) -> String {
    if u.is_empty() {
        fmt_num(x)
    } else if u.0.len() == 1 && u.0[0].0 == "%" && u.0[0].1 == Rational::ONE {
        format!("{}%", fmt_num(x))
    } else {
        format!("{} {}", fmt_num(x), u)
    }
}

const SPARK: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];

pub fn sparkline(sorted: &[f64], lo: f64, hi: f64, bins: usize) -> String {
    if hi <= lo {
        return String::new();
    }
    let mut counts = vec![0usize; bins];
    for &x in sorted {
        if x >= lo && x <= hi {
            let b = (((x - lo) / (hi - lo)) * bins as f64) as usize;
            counts[b.min(bins - 1)] += 1;
        }
    }
    let max = *counts.iter().max().unwrap_or(&1) as f64;
    counts
        .iter()
        .map(|&c| {
            if c == 0 {
                ' '
            } else {
                SPARK[((c as f64 / max) * 7.0).round() as usize]
            }
        })
        .collect()
}

/// Compact form: `12.5 t` or `12.5 t ± 0.6`.
/// True when a sample vector is constant but for round-off (a total whose
/// parts' uncertainties cancel, say), so it reads better as one number.
pub fn effectively_exact(s: &Samples) -> bool {
    match s {
        Samples::Exact(_) => true,
        Samples::Dist(v) => {
            let (lo, hi) = v
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
                    (a.min(x), b.max(x))
                });
            // (or round-off about zero: display rounds below 1e-12 to 0 anyway)
            hi - lo <= 1e-9 * hi.abs().max(lo.abs()) || hi - lo < 1e-12
        }
    }
}

pub fn fmt_quantity_short(q: &Quantity) -> String {
    let (u, scale) = q.display_unit();
    match &q.s {
        Samples::Exact(x) => with_unit(x / scale, &u),
        Samples::Dist(v) if effectively_exact(&q.s) => with_unit(v[0] / scale, &u),
        s => {
            let m = s.quantile(0.5) / scale;
            let sd = s.sd() / scale;
            if u.is_empty() {
                format!("{} ± {}", fmt_num(m), fmt_num(sd))
            } else {
                format!("{} ± {} {}", fmt_num(m), fmt_num(sd), u)
            }
        }
    }
}

/// Full form with 90% interval and histogram.
pub fn fmt_quantity(q: &Quantity) -> String {
    let (u, scale) = q.display_unit();
    match &q.s {
        Samples::Exact(x) => with_unit(x / scale, &u),
        Samples::Dist(v) if effectively_exact(&q.s) => with_unit(v[0] / scale, &u),
        s => {
            let v = s.sorted();
            let p = |q: f64| crate::samples::quantile_sorted(&v, q) / scale;
            let (p5, p50, p95) = (p(0.05), p(0.5), p(0.95));
            let sd = s.sd() / scale;
            let lo = crate::samples::quantile_sorted(&v, 0.005);
            let hi = crate::samples::quantile_sorted(&v, 0.995);
            format!(
                "{}  ±{}  [90%: {} .. {}]  {}",
                with_unit(p50, &u),
                fmt_num(sd),
                fmt_num(p5),
                fmt_num(p95),
                sparkline(&v, lo, hi, 16)
            )
        }
    }
}

pub fn fmt_value(v: &Value, indent: usize) -> String {
    let pad = " ".repeat(indent);
    match v {
        Value::Num(q) => fmt_quantity(q),
        Value::Bool(s) => match s {
            Samples::Exact(x) => (if *x != 0.0 { "true" } else { "false" }).into(),
            s => format!("P(true) = {}", fmt_num(s.mean())),
        },
        Value::Str(s) => format!("\"{}\"", s),
        Value::List(items) => {
            let simple = items
                .iter()
                .all(|i| matches!(i, Value::Num(_) | Value::Str(_) | Value::Bool(_)));
            if simple && items.len() <= 6 {
                let parts: Vec<String> = items
                    .iter()
                    .map(|i| match i {
                        Value::Num(q) => fmt_quantity_short(q),
                        other => fmt_value(other, 0),
                    })
                    .collect();
                format!("[{}]", parts.join(", "))
            } else {
                let mut s = "[\n".to_string();
                for i in items.iter() {
                    s += &format!("{}  {},\n", pad, fmt_value(i, indent + 2));
                }
                s + &pad + "]"
            }
        }
        Value::Record(fs) => {
            let w = fs.iter().map(|(n, _)| n.chars().count()).max().unwrap_or(0);
            let mut s = "{\n".to_string();
            for (n, v) in fs.iter() {
                s += &format!("{}  {:w$} : {}\n", pad, n, fmt_value(v, indent + 2), w = w);
            }
            s + &pad + "}"
        }
        Value::Func(c) => format!("<fn {}({})>", c.name, c.params.join(", ")),
        Value::Builtin(n) => format!("<builtin {}>", n),
        Value::NativeFn(f) => format!("<fn {}>", f.name),
        Value::Native(n) => n.describe(),
        Value::Unit => String::new(),
    }
}
