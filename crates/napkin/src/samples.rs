//! Sample vectors. Every uncertain number is N draws indexed by a shared "world" index,
//! so sample i of every variable comes from the same simulated world and correlations
//! between derived quantities are preserved automatically.

use std::cell::{Cell, RefCell};
use std::rc::Rc;

thread_local! {
    static N: Cell<usize> = const { Cell::new(10_000) };
    static RNG: RefCell<Rng> = RefCell::new(Rng::new(0x9E37_79B9_7F4A_7C15));
}

pub fn sample_count() -> usize {
    N.with(|n| n.get())
}

pub fn configure(n: usize, seed: u64) {
    N.with(|c| c.set(n));
    RNG.with(|r| *r.borrow_mut() = Rng::new(seed));
}

/// xoshiro256** seeded through splitmix64.
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Rng {
            s: [next(), next(), next(), next()],
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let r = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }

    /// Uniform in the open interval (0, 1).
    pub fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    pub fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    /// Marsaglia–Tsang gamma(shape, 1).
    pub fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let u = self.uniform();
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.uniform();
            if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
                return d * v;
            }
        }
    }

    pub fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        x / (x + y)
    }
}

/// Draw a fresh independent sample vector.
pub fn draw(mut f: impl FnMut(&mut Rng) -> f64) -> Samples {
    let n = sample_count();
    RNG.with(|r| {
        let mut r = r.borrow_mut();
        Samples::Dist(Rc::new((0..n).map(|_| f(&mut r)).collect()))
    })
}

#[derive(Clone, Debug)]
pub enum Samples {
    Exact(f64),
    Dist(Rc<Vec<f64>>),
}

impl Samples {
    pub fn is_exact(&self) -> bool {
        matches!(self, Samples::Exact(_))
    }

    pub fn get(&self, i: usize) -> f64 {
        match self {
            Samples::Exact(x) => *x,
            Samples::Dist(v) => v[i],
        }
    }

    pub fn exact(&self) -> Option<f64> {
        match self {
            Samples::Exact(x) => Some(*x),
            Samples::Dist(_) => None,
        }
    }

    pub fn map(&self, f: impl Fn(f64) -> f64) -> Samples {
        match self {
            Samples::Exact(x) => Samples::Exact(f(*x)),
            Samples::Dist(v) => Samples::Dist(Rc::new(v.iter().map(|&x| f(x)).collect())),
        }
    }

    pub fn zip(&self, o: &Samples, f: impl Fn(f64, f64) -> f64) -> Samples {
        match (self, o) {
            (Samples::Exact(a), Samples::Exact(b)) => Samples::Exact(f(*a, *b)),
            (Samples::Exact(a), Samples::Dist(v)) => {
                Samples::Dist(Rc::new(v.iter().map(|&b| f(*a, b)).collect()))
            }
            (Samples::Dist(v), Samples::Exact(b)) => {
                Samples::Dist(Rc::new(v.iter().map(|&a| f(a, *b)).collect()))
            }
            (Samples::Dist(v), Samples::Dist(w)) => Samples::Dist(Rc::new(
                v.iter().zip(w.iter()).map(|(&a, &b)| f(a, b)).collect(),
            )),
        }
    }

    /// Build from a per-sample function over several inputs; exact if all inputs are exact.
    pub fn from_fn(inputs: &[&Samples], f: impl Fn(usize) -> f64) -> Samples {
        if inputs.iter().all(|s| s.is_exact()) {
            return Samples::Exact(f(0));
        }
        Samples::Dist(Rc::new((0..sample_count()).map(f).collect()))
    }

    pub fn to_vec(&self) -> Vec<f64> {
        match self {
            Samples::Exact(x) => vec![*x; sample_count()],
            Samples::Dist(v) => v.as_ref().clone(),
        }
    }

    pub fn mean(&self) -> f64 {
        match self {
            Samples::Exact(x) => *x,
            Samples::Dist(v) => v.iter().sum::<f64>() / v.len() as f64,
        }
    }

    pub fn sd(&self) -> f64 {
        match self {
            Samples::Exact(_) => 0.0,
            Samples::Dist(v) => {
                let m = self.mean();
                (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt()
            }
        }
    }

    pub fn sorted(&self) -> Vec<f64> {
        let mut v = self.to_vec();
        v.sort_by(|a, b| a.total_cmp(b));
        v
    }

    pub fn quantile(&self, q: f64) -> f64 {
        match self {
            Samples::Exact(x) => *x,
            Samples::Dist(_) => quantile_sorted(&self.sorted(), q),
        }
    }
}

pub fn quantile_sorted(v: &[f64], q: f64) -> f64 {
    let pos = q.clamp(0.0, 1.0) * (v.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    v[lo] + (v[hi] - v[lo]) * (pos - lo as f64)
}

/// Ranks (average rank for ties) used by Spearman correlation.
pub fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0;
        for k in i..=j {
            r[idx[k]] = avg;
        }
        i = j + 1;
    }
    r
}

pub fn correlation(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    if saa == 0.0 || sbb == 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

/// Inverse standard normal CDF (Acklam's approximation).
pub fn norm_ppf(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    let pl = 0.02425;
    if p < pl {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - pl {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        -norm_ppf(1.0 - p)
    }
}
