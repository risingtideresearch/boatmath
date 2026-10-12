//! Arithmetic on quantities and (recursively) on lists and records of them.

use std::rc::Rc;

use crate::dim::{Dim, Rational};
use crate::samples::Samples;
use crate::units::UnitExpr;
use crate::value::{Quantity, Value};

pub type R<T> = Result<T, String>;

fn same_dim(op: &str, a: &Quantity, b: &Quantity) -> R<()> {
    if a.dim != b.dim {
        return Err(format!(
            "cannot {} {} and {}",
            op,
            a.dim.describe(),
            b.dim.describe()
        ));
    }
    Ok(())
}

fn add_hint(a: &Quantity, b: &Quantity) -> Option<UnitExpr> {
    match (&a.unit, &b.unit) {
        (Some(u), _) if !u.is_empty() || a.dim.is_none() => Some(u.clone()),
        (_, Some(u)) => Some(u.clone()),
        (u, None) => u.clone(),
    }
}

fn mul_hint(a: &Option<UnitExpr>, b: &Option<UnitExpr>, dim: &Dim) -> Option<UnitExpr> {
    if dim.is_none() {
        return Some(UnitExpr::default());
    }
    match (a, b) {
        (Some(x), Some(y)) => Some(x.mul(y).simplify()),
        _ => None,
    }
}

pub fn q_add(a: &Quantity, b: &Quantity) -> R<Quantity> {
    same_dim("add", a, b)?;
    Ok(Quantity::new(
        a.s.zip(&b.s, |x, y| x + y),
        a.dim,
        add_hint(a, b),
    ))
}

pub fn q_sub(a: &Quantity, b: &Quantity) -> R<Quantity> {
    same_dim("subtract", a, b)?;
    Ok(Quantity::new(
        a.s.zip(&b.s, |x, y| x - y),
        a.dim,
        add_hint(a, b),
    ))
}

pub fn q_mul(a: &Quantity, b: &Quantity) -> Quantity {
    let dim = a.dim.mul(&b.dim);
    Quantity::new(
        a.s.zip(&b.s, |x, y| x * y),
        dim,
        mul_hint(&a.unit, &b.unit, &dim),
    )
}

pub fn q_div(a: &Quantity, b: &Quantity) -> Quantity {
    let dim = a.dim.div(&b.dim);
    let binv = b.unit.as_ref().map(|u| u.inv());
    Quantity::new(
        a.s.zip(&b.s, |x, y| x / y),
        dim,
        mul_hint(&a.unit, &binv, &dim),
    )
}

pub fn q_neg(a: &Quantity) -> Quantity {
    Quantity::new(a.s.map(|x| -x), a.dim, a.unit.clone())
}

pub fn q_pow(a: &Quantity, b: &Quantity) -> R<Quantity> {
    if !b.dim.is_none() {
        return Err(format!(
            "exponent must be dimensionless, got {}",
            b.dim.describe()
        ));
    }
    if a.dim.is_none() {
        return Ok(Quantity::new(
            a.s.zip(&b.s, |x, y| x.powf(y)),
            Dim::NONE,
            Some(UnitExpr::default()),
        ));
    }
    let e =
        b.s.exact()
            .ok_or("a quantity with units can only be raised to an exact power")?;
    let r = Rational::from_f64(e).ok_or("unit exponent must be a simple fraction")?;
    Ok(q_powr(a, r))
}

pub fn q_powr(a: &Quantity, r: Rational) -> Quantity {
    let e = r.to_f64();
    let s = if r == Rational::int(2) {
        a.s.map(|x| x * x)
    } else {
        a.s.map(|x| x.powf(e))
    };
    let dim = a.dim.pow(r);
    let unit = if dim.is_none() {
        Some(UnitExpr::default())
    } else {
        a.unit.as_ref().map(|u| u.pow(r))
    };
    Quantity::new(s, dim, unit)
}

pub fn q_sqrt(a: &Quantity) -> Quantity {
    q_powr(a, Rational::new(1, 2))
}

pub fn compare(op: &str, a: &Quantity, b: &Quantity) -> R<Samples> {
    // Comparing against a literal zero is allowed for any dimension.
    let zero = |q: &Quantity| q.dim.is_none() && q.s.exact() == Some(0.0);
    if a.dim != b.dim && !zero(a) && !zero(b) {
        return Err(format!(
            "cannot compare {} with {}",
            a.dim.describe(),
            b.dim.describe()
        ));
    }
    let f: fn(f64, f64) -> bool = match op {
        "<" => |x, y| x < y,
        ">" => |x, y| x > y,
        "<=" => |x, y| x <= y,
        ">=" => |x, y| x >= y,
        "==" => |x, y| x == y,
        "!=" => |x, y| x != y,
        _ => unreachable!(),
    };
    Ok(a.s.zip(&b.s, |x, y| f(x, y) as u8 as f64))
}

/// Apply a binary op, broadcasting over lists.
pub fn binop(op: &str, a: &Value, b: &Value) -> R<Value> {
    use Value::*;
    match (a, b) {
        (Num(x), Num(y)) => Ok(match op {
            "+" => Num(q_add(x, y)?),
            "-" => Num(q_sub(x, y)?),
            "*" => Num(q_mul(x, y)),
            "/" => Num(q_div(x, y)),
            "^" => Num(q_pow(x, y)?),
            "<" | ">" | "<=" | ">=" | "==" | "!=" => Bool(compare(op, x, y)?),
            _ => return Err(format!("operator '{}' does not apply to numbers", op)),
        }),
        (Bool(x), Bool(y)) => Ok(match op {
            "and" => Bool(x.zip(y, |p, q| (p != 0.0 && q != 0.0) as u8 as f64)),
            "or" => Bool(x.zip(y, |p, q| (p != 0.0 || q != 0.0) as u8 as f64)),
            "==" => Bool(x.zip(y, |p, q| (p == q) as u8 as f64)),
            "!=" => Bool(x.zip(y, |p, q| (p != q) as u8 as f64)),
            _ => return Err(format!("operator '{}' does not apply to booleans", op)),
        }),
        (Str(x), Str(y)) => match op {
            "+" => Ok(Str(format!("{}{}", x, y).into())),
            "==" => Ok(Bool(Samples::Exact((x == y) as u8 as f64))),
            "!=" => Ok(Bool(Samples::Exact((x != y) as u8 as f64))),
            _ => Err(format!("operator '{}' does not apply to strings", op)),
        },
        (List(xs), List(ys)) => {
            if xs.len() != ys.len() {
                return Err(format!(
                    "list lengths differ ({} vs {})",
                    xs.len(),
                    ys.len()
                ));
            }
            let out: R<Vec<_>> = xs
                .iter()
                .zip(ys.iter())
                .map(|(x, y)| binop(op, x, y))
                .collect();
            Ok(List(Rc::new(out?)))
        }
        (List(xs), y @ Num(_)) => Ok(List(Rc::new(
            xs.iter().map(|x| binop(op, x, y)).collect::<R<_>>()?,
        ))),
        (x @ Num(_), List(ys)) => Ok(List(Rc::new(
            ys.iter().map(|y| binop(op, x, y)).collect::<R<_>>()?,
        ))),
        _ => Err(format!(
            "operator '{}' does not apply to {} and {}",
            op,
            a.type_name(),
            b.type_name()
        )),
    }
}

pub fn neg(v: &Value) -> R<Value> {
    match v {
        Value::Num(q) => Ok(Value::Num(q_neg(q))),
        Value::List(xs) => Ok(Value::List(Rc::new(xs.iter().map(neg).collect::<R<_>>()?))),
        _ => Err(format!("cannot negate a {}", v.type_name())),
    }
}

/// Per-sample select: where `c` is true take `a`, else `b`.
pub fn select(c: &Samples, a: &Value, b: &Value) -> R<Value> {
    use Value::*;
    match (a, b) {
        (Num(x), Num(y)) => {
            if x.dim != y.dim {
                return Err(format!(
                    "the branches of an uncertain 'if' must have the same dimension ({} vs {})",
                    x.dim.describe(),
                    y.dim.describe()
                ));
            }
            let s = Samples::from_fn(&[c, &x.s, &y.s], |i| {
                if c.get(i) != 0.0 {
                    x.s.get(i)
                } else {
                    y.s.get(i)
                }
            });
            Ok(Num(Quantity::new(s, x.dim, add_hint(x, y))))
        }
        (Bool(x), Bool(y)) => Ok(Bool(Samples::from_fn(&[c, x, y], |i| {
            if c.get(i) != 0.0 {
                x.get(i)
            } else {
                y.get(i)
            }
        }))),
        (List(xs), List(ys)) if xs.len() == ys.len() => Ok(List(Rc::new(
            xs.iter()
                .zip(ys.iter())
                .map(|(x, y)| select(c, x, y))
                .collect::<R<_>>()?,
        ))),
        (Record(xs), Record(ys)) if xs.len() == ys.len() => Ok(Record(Rc::new(
            xs.iter()
                .zip(ys.iter())
                .map(|((n, x), (m, y))| {
                    if n != m {
                        return Err(
                            "the branches of an uncertain 'if' must have the same fields"
                                .to_string(),
                        );
                    }
                    Ok((n.clone(), select(c, x, y)?))
                })
                .collect::<R<_>>()?,
        ))),
        (Str(x), Str(y)) if x == y => Ok(a.clone()),
        _ => Err(format!(
            "the branches of an uncertain 'if' must have the same shape ({} vs {})",
            a.type_name(),
            b.type_name()
        )),
    }
}
