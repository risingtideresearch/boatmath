//! Expensive physics on uncertain inputs: solve on a small grid, interpolate per sample.
//!
//! A float or a GZ curve costs a tenth of a second to seconds, so it can't be
//! solved once per sample. The quantities they produce are smooth in their load
//! inputs (mass, LCG, VCG), so each is solved at Chebyshev–Lobatto nodes spanning
//! the inputs' sampled range, one axis per uncertain input, and interpolated
//! (barycentric, tensor product) for every sample. An exact input is a single
//! node, so a fully exact call is one solve. A couple of samples are re-solved
//! exactly afterwards to measure the interpolation error.

use hullgeom::parallel;
use napkin::samples::sample_count;
use napkin::Samples;
use std::rc::Rc;

/// One input axis: its samples and its nodes.
struct Axis {
    nodes: Vec<f64>,
    weights: Vec<f64>,
}

impl Axis {
    fn new(s: &Samples, n: usize) -> Axis {
        let (lo, hi) = match s {
            Samples::Exact(x) => (*x, *x),
            Samples::Dist(v) => v
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
                    (a.min(x), b.max(x))
                }),
        };
        if hi - lo <= 1e-12 * hi.abs().max(1.0) || n < 2 {
            return Axis {
                nodes: vec![0.5 * (lo + hi)],
                weights: vec![1.0],
            };
        }
        let (c, h) = (0.5 * (lo + hi), 0.5 * (hi - lo));
        let nodes = (0..n)
            .map(|j| c + h * (std::f64::consts::PI * j as f64 / (n - 1) as f64).cos())
            .collect();
        let weights = (0..n)
            .map(|j| {
                let w = if j % 2 == 0 { 1.0 } else { -1.0 };
                if j == 0 || j == n - 1 {
                    0.5 * w
                } else {
                    w
                }
            })
            .collect();
        Axis { nodes, weights }
    }

    /// Piecewise-linear basis values at x: the two bracketing nodes.
    fn linear_basis(&self, x: f64, out: &mut Vec<f64>) {
        out.clear();
        out.resize(self.nodes.len(), 0.0);
        if self.nodes.len() == 1 {
            out[0] = 1.0;
            return;
        }
        // Lobatto nodes run from the top of the range down.
        let n = self.nodes.len();
        let k = (0..n - 1)
            .find(|&k| x >= self.nodes[k + 1])
            .unwrap_or(n - 2);
        let (a, b) = (self.nodes[k], self.nodes[k + 1]);
        let f = ((x - b) / (a - b)).clamp(0.0, 1.0);
        out[k] = f;
        out[k + 1] = 1.0 - f;
    }

    /// Barycentric basis values at x (sum to one).
    fn basis(&self, x: f64, out: &mut Vec<f64>) {
        out.clear();
        if self.nodes.len() == 1 {
            out.push(1.0);
            return;
        }
        if let Some(j) = self.nodes.iter().position(|&n| n == x) {
            out.extend((0..self.nodes.len()).map(|k| (k == j) as u8 as f64));
            return;
        }
        let mut sum = 0.0;
        for (n, w) in self.nodes.iter().zip(&self.weights) {
            let t = w / (x - n);
            out.push(t);
            sum += t;
        }
        for t in out.iter_mut() {
            *t /= sum;
        }
    }
}

pub struct Result {
    /// One sample vector per output.
    pub outputs: Vec<Samples>,
    /// Worst held-out error over the outputs, as a fraction of each output's
    /// spread over the grid; `None` when nothing needed interpolating.
    pub check: Option<f64>,
    pub solves: usize,
}

/// Evaluate `f(point) -> outputs` over the uncertain `inputs`, `nodes` per axis.
/// Outputs listed in `kinked` (the transom immersing, a deck edge going under)
/// are interpolated piecewise-linearly, which can't ring across a kink; the
/// rest by polynomial. The held-out check covers the smooth outputs.
pub fn evaluate(
    inputs: &[&Samples],
    nodes: &[usize],
    kinked: &[usize],
    f: impl Fn(&[f64]) -> std::result::Result<Vec<f64>, String> + Sync,
) -> std::result::Result<Result, String> {
    let axes: Vec<Axis> = inputs
        .iter()
        .zip(nodes)
        .map(|(s, &n)| Axis::new(s, n))
        .collect();
    let dims: Vec<usize> = axes.iter().map(|a| a.nodes.len()).collect();
    let total: usize = dims.iter().product();
    let point = |mut k: usize| -> Vec<f64> {
        // Last axis fastest.
        let mut p = vec![0.0; axes.len()];
        for a in (0..axes.len()).rev() {
            p[a] = axes[a].nodes[k % dims[a]];
            k /= dims[a];
        }
        p
    };
    let grid: Vec<Vec<f64>> = parallel::map_indexed(total, || (), |_, k| f(&point(k)))
        .into_iter()
        .collect::<std::result::Result<_, _>>()?;
    let nout = grid[0].len();
    if total == 1 {
        return Ok(Result {
            outputs: grid[0].iter().map(|&x| Samples::Exact(x)).collect(),
            check: None,
            solves: 1,
        });
    }

    // An output some node couldn't produce (NaN) is interpolated locally too,
    // so only the samples near that node lose it.
    let kinked: Vec<usize> = (0..nout)
        .filter(|o| kinked.contains(o) || grid.iter().any(|v| !v[*o].is_finite()))
        .collect();
    let kinked = &kinked[..];
    let smooth: Vec<usize> = (0..nout).filter(|o| !kinked.contains(o)).collect();
    let weights = |bases: &[Vec<f64>], k: usize| {
        let mut w = 1.0;
        let mut kk = k;
        for a in (0..axes.len()).rev() {
            w *= bases[a][kk % dims[a]];
            kk /= dims[a];
        }
        w
    };
    let interp = |i: usize, bases: &mut [Vec<f64>], lin: &mut [Vec<f64>], out: &mut [f64]| {
        for (a, ax) in axes.iter().enumerate() {
            ax.basis(inputs[a].get(i), &mut bases[a]);
            ax.linear_basis(inputs[a].get(i), &mut lin[a]);
        }
        out.iter_mut().for_each(|o| *o = 0.0);
        for (k, vals) in grid.iter().enumerate() {
            let w = weights(bases, k);
            if w != 0.0 {
                for &o in smooth.iter() {
                    out[o] += w * vals[o];
                }
            }
            if !kinked.is_empty() {
                let w = weights(lin, k);
                if w != 0.0 {
                    for &o in kinked {
                        out[o] += w * vals[o];
                    }
                }
            }
        }
    };

    let n = sample_count();
    let mut cols = vec![Vec::with_capacity(n); nout];
    let mut bases = vec![vec![]; axes.len()];
    let mut lin = vec![vec![]; axes.len()];
    let mut row = vec![0.0; nout];
    for i in 0..n {
        interp(i, &mut bases, &mut lin, &mut row);
        for (c, &x) in cols.iter_mut().zip(&row) {
            c.push(x);
        }
    }

    // Held-out check at two samples' inputs.
    // Each output's spread over the grid, floored at a millionth of the
    // largest output so a quantity that is ~0 everywhere doesn't read as noisy.
    let scale = grid
        .iter()
        .flatten()
        .filter(|x| x.is_finite())
        .fold(0.0f64, |a, x| a.max(x.abs()));
    let spread: Vec<f64> = (0..nout)
        .map(|o| {
            let (lo, hi) = grid
                .iter()
                .filter(|v| v[o].is_finite())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
                    (a.min(v[o]), b.max(v[o]))
                });
            (hi - lo).max(1e-6 * scale).max(1e-300)
        })
        .collect();
    let mut worst: f64 = 0.0;
    for &i in &[n / 3, (2 * n) / 3] {
        let p: Vec<f64> = inputs.iter().map(|s| s.get(i)).collect();
        let exact = f(&p)?;
        for &o in &smooth {
            if exact[o].is_finite() && cols[o][i].is_finite() {
                worst = worst.max((exact[o] - cols[o][i]).abs() / spread[o]);
            }
        }
    }
    Ok(Result {
        outputs: cols
            .into_iter()
            .map(|c| Samples::Dist(Rc::new(c)))
            .collect(),
        check: Some(worst),
        solves: total + 2,
    })
}
