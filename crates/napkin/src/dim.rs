//! Physical dimensions as vectors of rational exponents over the SI base dimensions.

use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Rational {
    pub num: i32,
    pub den: i32,
}

fn gcd(a: i32, b: i32) -> i32 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.max(1)
}

impl Rational {
    pub const ZERO: Rational = Rational { num: 0, den: 1 };
    pub const ONE: Rational = Rational { num: 1, den: 1 };

    pub fn new(num: i32, den: i32) -> Rational {
        assert!(den != 0, "zero denominator");
        let s = if den < 0 { -1 } else { 1 };
        let g = gcd(num, den);
        Rational {
            num: s * num / g,
            den: s * den / g,
        }
    }

    pub fn int(n: i32) -> Rational {
        Rational { num: n, den: 1 }
    }

    pub fn is_zero(self) -> bool {
        self.num == 0
    }

    pub fn to_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }

    /// Recover a small rational from a float exponent (e.g. 0.5 -> 1/2).
    pub fn from_f64(x: f64) -> Option<Rational> {
        for den in 1..=12 {
            let n = x * den as f64;
            if (n - n.round()).abs() < 1e-9 && n.abs() < 1e6 {
                return Some(Rational::new(n.round() as i32, den));
            }
        }
        None
    }
}

impl Add for Rational {
    type Output = Rational;
    fn add(self, o: Rational) -> Rational {
        Rational::new(self.num * o.den + o.num * self.den, self.den * o.den)
    }
}
impl Sub for Rational {
    type Output = Rational;
    fn sub(self, o: Rational) -> Rational {
        self + (-o)
    }
}
impl Neg for Rational {
    type Output = Rational;
    fn neg(self) -> Rational {
        Rational {
            num: -self.num,
            den: self.den,
        }
    }
}
impl Mul for Rational {
    type Output = Rational;
    fn mul(self, o: Rational) -> Rational {
        Rational::new(self.num * o.num, self.den * o.den)
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

pub const BASE_NAMES: [&str; 7] = ["m", "kg", "s", "A", "K", "mol", "cd"];
pub const BASE_DIM_NAMES: [&str; 7] = [
    "length",
    "mass",
    "time",
    "current",
    "temperature",
    "amount",
    "luminosity",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Dim(pub [Rational; 7]);

impl Dim {
    pub const NONE: Dim = Dim([Rational::ZERO; 7]);

    pub fn base(i: usize) -> Dim {
        let mut d = Dim::NONE;
        d.0[i] = Rational::ONE;
        d
    }

    pub fn length() -> Dim {
        Dim::base(0)
    }
    pub fn mass() -> Dim {
        Dim::base(1)
    }

    pub fn is_none(&self) -> bool {
        self.0.iter().all(|r| r.is_zero())
    }

    pub fn mul(&self, o: &Dim) -> Dim {
        let mut d = *self;
        for i in 0..7 {
            d.0[i] = d.0[i] + o.0[i];
        }
        d
    }

    pub fn div(&self, o: &Dim) -> Dim {
        let mut d = *self;
        for i in 0..7 {
            d.0[i] = d.0[i] - o.0[i];
        }
        d
    }

    pub fn pow(&self, p: Rational) -> Dim {
        let mut d = *self;
        for i in 0..7 {
            d.0[i] = d.0[i] * p;
        }
        d
    }

    /// Human-readable description, e.g. "length^2·mass".
    pub fn describe(&self) -> String {
        if self.is_none() {
            return "dimensionless".into();
        }
        let mut parts = vec![];
        for i in 0..7 {
            let e = self.0[i];
            if !e.is_zero() {
                if e == Rational::ONE {
                    parts.push(BASE_DIM_NAMES[i].to_string());
                } else {
                    parts.push(format!("{}^{}", BASE_DIM_NAMES[i], e));
                }
            }
        }
        parts.join("·")
    }
}
