//! Unit table and unit expressions.
//!
//! A unit is a scale factor relative to SI plus a dimension. Unit expressions such as
//! `kg/m^3` or `t·m` are products of named units raised to rational powers; values keep
//! the expression they were written in as a display hint.

use crate::dim::{Dim, Rational};

#[derive(Clone, Debug)]
pub struct UnitDef {
    pub scale: f64,
    pub dim: Dim,
}

const PREFIXES: &[(&str, f64)] = &[
    ("Y", 1e24),
    ("Z", 1e21),
    ("E", 1e18),
    ("P", 1e15),
    ("T", 1e12),
    ("G", 1e9),
    ("M", 1e6),
    ("k", 1e3),
    ("h", 1e2),
    ("da", 1e1),
    ("d", 1e-1),
    ("c", 1e-2),
    ("m", 1e-3),
    ("u", 1e-6),
    ("µ", 1e-6),
    ("n", 1e-9),
    ("p", 1e-12),
];

fn d(exps: [i32; 7]) -> Dim {
    let mut out = Dim::NONE;
    for i in 0..7 {
        out.0[i] = Rational::int(exps[i]);
    }
    out
}

/// (name, scale to SI, dimension exponents [m kg s A K mol cd], accepts SI prefixes)
fn table() -> Vec<(&'static str, f64, [i32; 7], bool)> {
    const L: [i32; 7] = [1, 0, 0, 0, 0, 0, 0];
    const M: [i32; 7] = [0, 1, 0, 0, 0, 0, 0];
    const T: [i32; 7] = [0, 0, 1, 0, 0, 0, 0];
    const NONE: [i32; 7] = [0; 7];
    const FORCE: [i32; 7] = [1, 1, -2, 0, 0, 0, 0];
    const ENERGY: [i32; 7] = [2, 1, -2, 0, 0, 0, 0];
    const POWER: [i32; 7] = [2, 1, -3, 0, 0, 0, 0];
    const PRESSURE: [i32; 7] = [-1, 1, -2, 0, 0, 0, 0];
    const VOLUME: [i32; 7] = [3, 0, 0, 0, 0, 0, 0];
    const AREA: [i32; 7] = [2, 0, 0, 0, 0, 0, 0];
    const SPEED: [i32; 7] = [1, 0, -1, 0, 0, 0, 0];
    const FREQ: [i32; 7] = [0, 0, -1, 0, 0, 0, 0];
    vec![
        // SI base and derived
        ("m", 1.0, L, true),
        ("g", 1e-3, M, true),
        ("s", 1.0, T, true),
        ("A", 1.0, [0, 0, 0, 1, 0, 0, 0], true),
        ("K", 1.0, [0, 0, 0, 0, 1, 0, 0], true),
        ("mol", 1.0, [0, 0, 0, 0, 0, 1, 0], true),
        ("cd", 1.0, [0, 0, 0, 0, 0, 0, 1], false),
        ("N", 1.0, FORCE, true),
        ("J", 1.0, ENERGY, true),
        ("W", 1.0, POWER, true),
        ("Pa", 1.0, PRESSURE, true),
        ("Hz", 1.0, FREQ, true),
        ("V", 1.0, [2, 1, -3, -1, 0, 0, 0], true),
        ("L", 1e-3, VOLUME, true),
        ("bar", 1e5, PRESSURE, true),
        // time
        ("min", 60.0, T, false),
        ("h", 3600.0, T, false),
        ("hr", 3600.0, T, false),
        ("day", 86400.0, T, false),
        // mass (naval)
        ("t", 1000.0, M, true),
        ("tonne", 1000.0, M, false),
        ("LT", 1016.0469088, M, false),
        ("ST", 907.18474, M, false),
        ("lb", 0.45359237, M, false),
        // length
        ("ft", 0.3048, L, false),
        ("in", 0.0254, L, false),
        ("yd", 0.9144, L, false),
        ("nmi", 1852.0, L, false),
        ("mi", 1609.344, L, false),
        ("ha", 1e4, AREA, false),
        ("gal", 3.785411784e-3, VOLUME, false),
        // speed
        ("kn", 1852.0 / 3600.0, SPEED, false),
        ("knot", 1852.0 / 3600.0, SPEED, false),
        // force / pressure / power
        ("tf", 9806.65, FORCE, false),
        ("lbf", 4.4482216152605, FORCE, false),
        ("kip", 4448.2216152605, FORCE, false),
        ("psi", 6894.757293168, PRESSURE, false),
        ("hp", 745.69987158227, POWER, false),
        ("Wh", 3600.0, ENERGY, true),
        // dimensionless
        ("rad", 1.0, NONE, false),
        ("deg", std::f64::consts::PI / 180.0, NONE, false),
        ("%", 0.01, NONE, false),
    ]
}

pub fn lookup_unit(name: &str) -> Option<UnitDef> {
    let tbl = table();
    if let Some((_, s, e, _)) = tbl.iter().find(|u| u.0 == name) {
        return Some(UnitDef {
            scale: *s,
            dim: d(*e),
        });
    }
    for (p, ps) in PREFIXES {
        if let Some(rest) = name.strip_prefix(p) {
            if let Some((_, s, e, _)) = tbl.iter().find(|u| u.0 == rest && u.3) {
                return Some(UnitDef {
                    scale: s * ps,
                    dim: d(*e),
                });
            }
        }
    }
    None
}

/// A product of named units with rational exponents, e.g. kg·m^2/s^2.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct UnitExpr(pub Vec<(String, Rational)>);

impl UnitExpr {
    pub fn single(name: &str) -> UnitExpr {
        UnitExpr(vec![(name.to_string(), Rational::ONE)])
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Resolve to (scale, dim). Fails on unknown unit names.
    pub fn resolve(&self) -> Result<(f64, Dim), String> {
        let mut scale = 1.0;
        let mut dim = Dim::NONE;
        for (name, e) in &self.0 {
            let u = lookup_unit(name).ok_or_else(|| format!("unknown unit '{}'", name))?;
            scale *= u.scale.powf(e.to_f64());
            dim = dim.mul(&u.dim.pow(*e));
        }
        Ok((scale, dim))
    }

    pub fn mul(&self, o: &UnitExpr) -> UnitExpr {
        let mut out = self.0.clone();
        for (n, e) in &o.0 {
            if let Some(p) = out.iter_mut().find(|(m, _)| m == n) {
                p.1 = p.1 + *e;
            } else {
                out.push((n.clone(), *e));
            }
        }
        out.retain(|(_, e)| !e.is_zero());
        UnitExpr(out)
    }

    /// Merge units that measure the same thing into the first one written, so
    /// `m^2·mm·t/m^3` becomes `t` rather than `mm·t/m`.
    pub fn simplify(&self) -> UnitExpr {
        let mut out: Vec<(String, Rational, Option<Dim>)> = vec![];
        for (n, e) in &self.0 {
            let d = lookup_unit(n).map(|u| u.dim).filter(|d| !d.is_none());
            match out
                .iter_mut()
                .find(|(m, _, dm)| m == n || (d.is_some() && *dm == d))
            {
                Some(p) => p.1 = p.1 + *e,
                None => out.push((n.clone(), *e, d)),
            }
        }
        UnitExpr(
            out.into_iter()
                .filter(|(_, e, _)| !e.is_zero())
                .map(|(n, e, _)| (n, e))
                .collect(),
        )
    }

    pub fn pow(&self, p: Rational) -> UnitExpr {
        let mut out: Vec<_> = self.0.iter().map(|(n, e)| (n.clone(), *e * p)).collect();
        out.retain(|(_, e)| !e.is_zero());
        UnitExpr(out)
    }

    pub fn inv(&self) -> UnitExpr {
        self.pow(Rational::int(-1))
    }
}

impl std::fmt::Display for UnitExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let num: Vec<_> = self.0.iter().filter(|(_, e)| e.num > 0).collect();
        let den: Vec<_> = self.0.iter().filter(|(_, e)| e.num < 0).collect();
        let fmt_part = |(n, e): &(String, Rational)| {
            if *e == Rational::ONE {
                n.clone()
            } else {
                format!("{}^{}", n, e)
            }
        };
        let mut s = num
            .iter()
            .map(|p| fmt_part(p))
            .collect::<Vec<_>>()
            .join("·");
        if !den.is_empty() {
            if s.is_empty() {
                s.push('1');
            }
            for (n, e) in den {
                s.push('/');
                s.push_str(&fmt_part(&(n.clone(), -*e)));
            }
        }
        write!(f, "{}", s)
    }
}

/// Preferred names for common SI dimensions when a value has no display hint.
pub fn si_unit_for(dim: &Dim) -> UnitExpr {
    for name in ["N", "J", "W", "Pa"] {
        if lookup_unit(name).unwrap().dim == *dim {
            return UnitExpr::single(name);
        }
    }
    let mut parts = vec![];
    for i in 0..7 {
        let e = dim.0[i];
        if !e.is_zero() {
            parts.push((crate::dim::BASE_NAMES[i].to_string(), e));
        }
    }
    UnitExpr(parts)
}
