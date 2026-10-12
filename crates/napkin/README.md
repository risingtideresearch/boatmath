# napkin

A small language for engineering calculations where every number carries
**units** and **uncertainty**. This crate is the language. `napkin-boat` adds
boatmath's hulls and builds the `napkin` command.

```
let h    = 12 ± 0.5 m
let mass = 2 to 5 kg
let E    = mass * g0 * h -> kJ

E                  # E = 0.37 kJ  ±0.1114  [90%: 0.2341 .. 0.5914]  ▂▄▆▇█▇▆▅▄▃▂▂▁▁▁▁
P(E > 0.4 kJ)      # 0.391
sensitivity(E)     # ↑ mass 97.5%, ↑ h 2.5%
h + mass           # line 7: cannot add length and mass
```

## How uncertainty works

Each uncertain value is a vector of N samples (10,000 by default). Sample *i*
of every variable belongs to the same simulated "world", so arithmetic works
sample by sample and correlations are kept: `h - h` is exactly zero, and a CG
computed from item masses is correlated with the total mass. Exact values stay
scalars and cost nothing.

Every primitive distribution is recorded as an input with a label taken from
where it was written: `shell.mass`, `Crew.cg[0]`, or an item's `name`.
`sensitivity(x)` ranks those inputs by squared rank correlation with `x`. That
share is an estimate of how much of the variance each input explains, and it's
most reliable when the model is roughly monotone.

## Language

### Numbers and units

```
12 m        9.81 m/s^2      7.85 t/m^3     5 kN·m     m^-1    1.5e-3 m
[48, 0, 4.2] m              # a unit after ] or ) applies to every element
(12 ± 0.5) mm
x -> ft                     # convert for display (the dimension must match)
```

- A unit expression ends at whitespace, so `5 m/s` is a speed and `5 m / t`
  divides by a variable `t`. Use `*` to multiply by a variable: `2 * x`, not
  `2 x`.
- Supported units:
  - SI with prefixes: m, g, s, N, J, W, Pa, Hz, L, …
  - Mass: t, LT, ST, lb.
  - Length, area, volume: ft, in, yd, nmi, mi, ha, gal.
  - Other: kn, tf, lbf, kip, psi, hp, bar, min, h, day, deg, rad, %.
- Results keep the units you wrote. Units of the same kind merge, so
  `m^2 · mm · t/m^3` displays as `t`.
- Dimensions are checked when the program runs. Adding or comparing values with
  different dimensions is an error with a line number.

### Uncertainty

| syntax | meaning |
|---|---|
| `12 ± 0.5 m` (or `+-`) | normal, standard deviation 0.5 m |
| `200 t ± 5%` | normal, standard deviation 5% of the value |
| `2 to 5 t` | 90% interval: lognormal if both ends are positive, normal otherwise |
| `normal(m, sd)`, `uniform(a, b)` | |
| `triangular(a, mode, b)`, `pert(a, mode, b)` | three-point estimates |
| `bernoulli(p)` | uncertain bool, e.g. an optional item |

- **Statistics:** `mean`, `median`, `sd`, `quantile(x, q)`, `P(condition)`,
  `sensitivity(x)`.
- **Comparisons:** a comparison gives an uncertain bool. An `if` on an uncertain
  bool runs both branches and picks per sample, so
  `if bernoulli(0.3) then 12 t else 0 t` is a mixture.

### Values and functions

```
let x = 3 m
fn area_of(b, h) = b * h / 2
fn f(x) {                       # block body: the last expression is the result
  let k = 2.6
  k * sqrt(x)
}
|x| x * 2                       # lambda
[1, 2, 3]   xs[0]   xs[-1]      # lists; arithmetic broadcasts element-wise
{name: "Fuel", mass: 30 t}      # records; items.mass maps the field over a list
if c then a else b    and  or  not    == != < > <= >=
```

- **Lists:** `len sum map filter concat reverse range linspace`, and
  `trapz(xs, ys)` for ∫ y dx by the trapezoid rule.
- **Vectors:** `dot cross norm`.
- **Math:** `sqrt abs exp ln log10 sin cos tan asin acos atan atan2 floor ceil
  round min max`.
- **Constants:** `pi g0 rho_sw rho_fw rho_steel rho_al`.

### Geometry

Points are `[x, y]` or `[x, y, z]` lists, and any coordinate can be uncertain.

- `area(pts)`, `centroid(pts)`: a 2D polygon, or a planar polygon in 3D.
- `section(pts)`: for a 2D polygon, returns `area`, `centroid`, centroidal
  `ixx iyy ixy`, and `ixx_axis iyy_axis` about the coordinate axes.
- `polyline_length(pts)`.

### Weights

An item is any record with `mass` and `cg: [x, y, z]`, plus an optional `name`.

- `combine(items, name?)`: total mass and CG. The result is itself an item, so
  groups nest (lightship → full load).
- `inertia(items)`: mass moments `ixx iyy izz` about the combined CG and radii
  of gyration `kxx kyy kzz`. Items are treated as point masses.
- `weight_table(items)`: prints mass, LCG, TCG and VCG as median ± sd, with a
  total row.

## Embedding

```rust
let mut it = napkin::Interp::new(10_000, seed);
it.register("my_fn", |it, args| { /* … */ Ok(napkin::Value::Unit) });
it.run(source)?;
```

- **Native functions:** `Interp::register` adds a function written in Rust.
  `napkin::builtins` has the argument helpers (`num`, `list`, `arity`,
  `item_parts`) and `per_sample`.
- **Opaque values:** a value such as a hull implements `value::Native`.
- **Output:** printing goes through `Interp::output`, and `run_captured`
  collects it into a string.

## Layout

| file | contents |
|---|---|
| `lexer.rs`, `parser.rs`, `ast.rs` | syntax |
| `dim.rs`, `units.rs` | rational dimension vectors and the unit table |
| `samples.rs` | sample vectors, RNG (xoshiro256**), quantiles, rank correlation |
| `interp.rs`, `ops.rs`, `value.rs` | evaluation, broadcasting arithmetic, display, the extension API |
| `builtins.rs` | distributions, statistics, geometry, weights |
