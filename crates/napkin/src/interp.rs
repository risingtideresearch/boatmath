use std::cell::RefCell;
use std::rc::Rc;

use crate::ast::{Expr, Stmt, StmtKind};
use crate::builtins;

use crate::ops::{self, R};
use crate::samples::{self, Samples};
use crate::units::UnitExpr;
use crate::value::{lookup, Closure, Env, NativeFn, Quantity, Scope, Value};

/// Every primitive distribution is registered as an input source so `sensitivity`
/// can attribute output variance to it. A source is labelled by where it was written:
/// the enclosing `let`, record fields, list positions, or an item's `name`.
pub struct Source {
    pub label: String,
    pub samples: Rc<Vec<f64>>,
}

thread_local! {
    pub static SOURCES: RefCell<Vec<Source>> = const { RefCell::new(vec![]) };
    static PATH: RefCell<Vec<String>> = const { RefCell::new(vec![]) };
}

pub fn register_source(s: &Samples) {
    if let Samples::Dist(v) = s {
        let label = PATH.with(|p| {
            let mut out = String::new();
            for seg in p.borrow().iter() {
                if let Some(name) = seg.strip_prefix('@') {
                    // a named item: restart the label at its name
                    out = name.to_string();
                } else if seg.starts_with('[') || out.is_empty() {
                    out.push_str(seg);
                } else {
                    out.push('.');
                    out.push_str(seg);
                }
            }
            if out.is_empty() {
                "(unnamed)".to_string()
            } else {
                out
            }
        });
        SOURCES.with(|src| {
            src.borrow_mut().push(Source {
                label,
                samples: v.clone(),
            })
        });
    }
}

/// Evaluate `f` with `seg` appended to the source-label path.
fn labelled<T>(seg: String, f: impl FnOnce() -> T) -> T {
    PATH.with(|p| p.borrow_mut().push(seg));
    let r = f();
    PATH.with(|p| p.borrow_mut().pop());
    r
}

/// A record literal's `name: "..."` field, if it has one.
fn record_name(e: &Expr) -> Option<&str> {
    match e {
        Expr::Record(fs) => fs.iter().find_map(|(n, v)| match v {
            Expr::Str(s) if n == "name" => Some(s.as_str()),
            _ => None,
        }),
        _ => None,
    }
}

pub struct Interp {
    pub global: Env,
    /// Source text, used to echo top-level expressions.
    pub lines: Vec<String>,
    /// Where printed output goes (stdout by default).
    pub output: Box<dyn FnMut(&str)>,
}

pub fn err_at<T>(line: usize, e: String) -> R<T> {
    if e.starts_with("line ") {
        Err(e)
    } else {
        Err(format!("line {}: {}", line, e))
    }
}

impl Interp {
    /// A fresh interpreter drawing `n` samples per distribution from `seed`.
    pub fn new(n: usize, seed: u64) -> Interp {
        samples::configure(n, seed);
        SOURCES.with(|s| s.borrow_mut().clear());
        let global = Rc::new(RefCell::new(Scope::default()));
        builtins::install(&global);
        Interp {
            global,
            lines: vec![],
            output: Box::new(|s| println!("{}", s)),
        }
    }

    pub fn emit(&mut self, s: &str) {
        (self.output)(s)
    }

    /// Make a host function callable from napkin code as `name(...)`.
    pub fn register(
        &mut self,
        name: &str,
        f: impl Fn(&mut Interp, Vec<Value>) -> R<Value> + 'static,
    ) {
        let v = Value::NativeFn(Rc::new(NativeFn {
            name: name.to_string(),
            f: Box::new(f),
        }));
        self.global.borrow_mut().vars.insert(name.to_string(), v);
    }

    pub fn run(&mut self, src: &str) -> R<()> {
        self.lines = src.lines().map(|s| s.to_string()).collect();
        let toks = crate::lexer::lex(src)?;
        let prog = crate::parser::Parser::new(toks).program()?;
        let env = self.global.clone();
        for st in &prog {
            let v = match self.exec(st, &env) {
                Ok(v) => v,
                Err(e) => return err_at(st.line, e),
            };
            if let (StmtKind::Expr(e), Some(v)) = (&st.kind, v) {
                if matches!(v, Value::Unit) {
                    continue;
                }
                let label = self.source_label(st.line, e);
                let body = crate::value::fmt_value(&v, 0);
                let text = match label {
                    Some(l) if !body.contains('\n') => format!("{} = {}", l, body),
                    Some(l) => format!("{} =\n{}", l, body),
                    None => body,
                };
                self.emit(&text);
            }
        }
        Ok(())
    }

    fn source_label(&self, line: usize, e: &Expr) -> Option<String> {
        if matches!(e, Expr::Call(f, _) if matches!(f.as_ref(), Expr::Var(n) if n == "print")) {
            return None;
        }
        let text = self.lines.get(line - 1)?;
        let text = text.split('#').next().unwrap_or("").trim();
        if text.is_empty() {
            None
        } else {
            Some(text.to_string())
        }
    }

    pub fn exec(&mut self, st: &Stmt, env: &Env) -> R<Option<Value>> {
        let res = match &st.kind {
            StmtKind::Let(name, e) => {
                let v = labelled(name.clone(), || self.eval(e, env))?;
                env.borrow_mut().vars.insert(name.clone(), v);
                None
            }
            StmtKind::Fn(name, params, body) => {
                let c = Closure {
                    name: name.clone(),
                    params: params.clone(),
                    body: body.clone(),
                    env: env.clone(),
                };
                env.borrow_mut()
                    .vars
                    .insert(name.clone(), Value::Func(Rc::new(c)));
                None
            }
            StmtKind::Expr(e) => Some(self.eval(e, env)?),
        };
        Ok(res)
    }

    pub fn eval(&mut self, e: &Expr, env: &Env) -> R<Value> {
        match e {
            Expr::Num(x, None) => Ok(Value::Num(Quantity::plain(*x))),
            Expr::Num(x, Some(u)) => {
                let (scale, dim) = u.resolve()?;
                Ok(Value::Num(Quantity::new(Samples::Exact(x * scale), dim, Some(u.clone()))))
            }
            Expr::Str(s) => Ok(Value::Str(s.as_str().into())),
            Expr::Var(n) => match n.as_str() {
                "true" => Ok(Value::Bool(Samples::Exact(1.0))),
                "false" => Ok(Value::Bool(Samples::Exact(0.0))),
                _ => lookup(env, n).ok_or_else(|| {
                    if crate::units::lookup_unit(n).is_some() {
                        format!("unknown variable '{}' ('{}' is a unit — write it after a number, e.g. 1 {})", n, n, n)
                    } else {
                        format!("unknown variable '{}'", n)
                    }
                }),
            },
            Expr::List(items) => {
                let mut vs = vec![];
                for (i, e) in items.iter().enumerate() {
                    let seg = match record_name(e) {
                        Some(n) => format!("@{}", n),
                        None => format!("[{}]", i),
                    };
                    vs.push(labelled(seg, || self.eval(e, env))?);
                }
                Ok(Value::list(vs))
            }
            Expr::Record(fields) => {
                let mut out = vec![];
                for (n, e) in fields {
                    out.push((n.clone(), labelled(n.clone(), || self.eval(e, env))?));
                }
                Ok(Value::Record(Rc::new(out)))
            }
            Expr::Block(stmts) => {
                let scope = Rc::new(RefCell::new(Scope { vars: Default::default(), parent: Some(env.clone()) }));
                let mut last = Value::Unit;
                for st in stmts {
                    last = match self.exec(st, &scope) {
                        Ok(v) => v.unwrap_or(Value::Unit),
                        Err(e) => return err_at(st.line, e),
                    };
                }
                Ok(last)
            }
            Expr::Unary(op, x) => {
                let v = self.eval(x, env)?;
                match *op {
                    "-" => ops::neg(&v),
                    "not" => match v {
                        Value::Bool(s) => Ok(Value::Bool(s.map(|x| (x == 0.0) as u8 as f64))),
                        _ => Err(format!("'not' needs a bool, got {}", v.type_name())),
                    },
                    _ => unreachable!(),
                }
            }
            Expr::Binary(op, a, b) => {
                let a = self.eval(a, env)?;
                // Short-circuit exact booleans.
                if let Value::Bool(Samples::Exact(x)) = &a {
                    if (*op == "and" && *x == 0.0) || (*op == "or" && *x != 0.0) {
                        return Ok(a);
                    }
                }
                let b = self.eval(b, env)?;
                ops::binop(op, &a, &b)
            }
            Expr::PlusMinus(a, b) => {
                let (a, b) = (self.num(a, env)?, self.num(b, env)?);
                if a.dim != b.dim && !a.dim.is_none() {
                    return Err(format!("'±' needs both sides in the same dimension ({} vs {})", a.dim.describe(), b.dim.describe()));
                }
                let z = samples::draw(|r| r.normal());
                let s = Samples::from_fn(&[&z], |i| a.s.get(i) + b.s.get(i) * z.get(i));
                register_source(&s);
                let unit = if a.dim.is_none() && !b.dim.is_none() { b.unit.clone() } else { a.unit.clone() };
                Ok(Value::Num(Quantity::new(s, b.dim, unit)))
            }
            Expr::PlusMinusRel(a, p) => {
                let a = self.num(a, env)?;
                let z = samples::draw(|r| r.normal());
                let s = Samples::from_fn(&[&z], |i| a.s.get(i) * (1.0 + p * z.get(i)));
                register_source(&s);
                Ok(Value::Num(Quantity::new(s, a.dim, a.unit.clone())))
            }
            Expr::To(a, b) => {
                let (a, b) = (self.num(a, env)?, self.num(b, env)?);
                Ok(Value::Num(builtins::ci90(&a, &b)?))
            }
            Expr::ApplyUnit(x, u) => {
                let v = self.eval(x, env)?;
                apply_unit(&v, u)
            }
            Expr::Convert(x, u) => {
                let v = self.eval(x, env)?;
                convert(&v, u)
            }
            Expr::Call(f, args) => {
                let fv = self.eval(f, env)?;
                let argv = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
                self.call(&fv, argv)
            }
            Expr::Field(x, name) => {
                let v = self.eval(x, env)?;
                if let Value::List(items) = &v {
                    // items.mass -> list of each item's mass
                    if items.iter().all(|i| i.field(name).is_some()) && !items.is_empty() {
                        return Ok(Value::list(items.iter().map(|i| i.field(name).unwrap().clone()).collect()));
                    }
                }
                v.field(name).cloned().ok_or_else(|| format!("{} has no field '{}'", v.type_name(), name))
            }
            Expr::Index(x, i) => {
                let v = self.eval(x, env)?;
                let i = self.eval(i, env)?;
                index(&v, &i)
            }
            Expr::If(c, a, b) => {
                let c = self.eval(c, env)?;
                match c {
                    Value::Bool(Samples::Exact(x)) => {
                        if x != 0.0 {
                            self.eval(a, env)
                        } else {
                            self.eval(b, env)
                        }
                    }
                    Value::Bool(s) => {
                        let va = self.eval(a, env)?;
                        let vb = self.eval(b, env)?;
                        ops::select(&s, &va, &vb)
                    }
                    other => Err(format!("'if' needs a bool condition, got {}", other.type_name())),
                }
            }
            Expr::Lambda(params, body) => Ok(Value::Func(Rc::new(Closure {
                name: "lambda".into(),
                params: params.clone(),
                body: (**body).clone(),
                env: env.clone(),
            }))),
        }
    }

    fn num(&mut self, e: &Expr, env: &Env) -> R<Quantity> {
        match self.eval(e, env)? {
            Value::Num(q) => Ok(q),
            v => Err(format!("expected a number, got {}", v.type_name())),
        }
    }

    pub fn call(&mut self, f: &Value, args: Vec<Value>) -> R<Value> {
        match f {
            Value::Func(c) => {
                if c.params.len() != args.len() {
                    return Err(format!(
                        "{} takes {} argument(s), got {}",
                        c.name,
                        c.params.len(),
                        args.len()
                    ));
                }
                let scope = Rc::new(RefCell::new(Scope {
                    vars: Default::default(),
                    parent: Some(c.env.clone()),
                }));
                for (p, a) in c.params.iter().zip(args) {
                    scope.borrow_mut().vars.insert(p.clone(), a);
                }
                self.eval(&c.body, &scope)
                    .map_err(|e| format!("in {}: {}", c.name, e))
            }
            Value::Builtin(name) => {
                builtins::call(self, name, args).map_err(|e| format!("{}: {}", name, e))
            }
            Value::NativeFn(nf) => (nf.f)(self, args).map_err(|e| format!("{}: {}", nf.name, e)),
            other => Err(format!("cannot call a {}", other.type_name())),
        }
    }
}

fn index(v: &Value, i: &Value) -> R<Value> {
    match (v, i) {
        (Value::List(xs), Value::Num(q)) if q.dim.is_none() => {
            let k = q.s.exact().ok_or("list index must be exact")?;
            let n = xs.len() as i64;
            let mut k = k as i64;
            if k < 0 {
                k += n;
            }
            xs.get(k as usize)
                .cloned()
                .ok_or_else(|| format!("index {} out of range (length {})", k, n))
        }
        (Value::Record(_), Value::Str(s)) => v
            .field(s)
            .cloned()
            .ok_or_else(|| format!("no field '{}'", s)),
        _ => Err(format!(
            "cannot index a {} with a {}",
            v.type_name(),
            i.type_name()
        )),
    }
}

pub fn apply_unit(v: &Value, u: &UnitExpr) -> R<Value> {
    let (scale, dim) = u.resolve()?;
    match v {
        Value::Num(q) if q.dim.is_none() => Ok(Value::Num(Quantity::new(
            q.s.map(|x| x * scale),
            dim,
            Some(u.clone()),
        ))),
        Value::Num(q) => Err(format!(
            "cannot attach unit {} to a value that already has dimension {}",
            u,
            q.dim.describe()
        )),
        Value::List(xs) => Ok(Value::list(
            xs.iter().map(|x| apply_unit(x, u)).collect::<R<_>>()?,
        )),
        _ => Err(format!("cannot attach a unit to a {}", v.type_name())),
    }
}

pub fn convert(v: &Value, u: &UnitExpr) -> R<Value> {
    let (_, dim) = u.resolve()?;
    match v {
        Value::Num(q) => {
            if q.dim != dim {
                return Err(format!(
                    "cannot convert {} to {} ({})",
                    q.dim.describe(),
                    u,
                    dim.describe()
                ));
            }
            Ok(Value::Num(Quantity::new(
                q.s.clone(),
                q.dim,
                Some(u.clone()),
            )))
        }
        Value::List(xs) => Ok(Value::list(
            xs.iter().map(|x| convert(x, u)).collect::<R<_>>()?,
        )),
        Value::Record(fs) => Ok(Value::Record(Rc::new(
            fs.iter()
                .map(|(n, x)| match x {
                    // Convert only the fields with a matching dimension.
                    Value::Num(q) if q.dim != dim => Ok((n.clone(), x.clone())),
                    _ => Ok((
                        n.clone(),
                        convert(x, u).or_else(|_| Ok::<_, String>(x.clone()))?,
                    )),
                })
                .collect::<R<_>>()?,
        ))),
        _ => Err(format!("cannot convert a {}", v.type_name())),
    }
}
