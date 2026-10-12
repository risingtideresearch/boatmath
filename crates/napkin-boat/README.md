# napkin-boat

The `napkin` command: the [napkin](../napkin) language with hulls from
`hullgeom`. A weight study's uncertain mass and CG run through the free-trim
float and the GZ curve, and come out as distributions with probabilities and
the inputs that drive them.

```
cargo run --release -p napkin-boat -- crates/napkin-boat/examples/graf_weights.nap
```

```
let cat  = hull("graf-demihull.igs", {span: 1.6 m})
let ship = combine(items)              # uncertain mass and CG
let s    = stability(cat, ship)

s.max_gz           # 0.622 m  ±0.03113  [90%: 0.5596 .. 0.6607]
P(s.vanishing < 65 deg)
sensitivity(s.max_gz)                  # ↓ Crew.mass 39%, ↓ Crew.cg[1] 34%, ↑ Crew.cg[0] 17%, …
```

## Functions

| function | returns |
|---|---|
| `hull(path, {waterline, span, units, stations, rays, centerplane})` | a platform from IGES or STL, cut as boatmath cuts it. `waterline` is the design waterline's CAD height. `span` doubles a single hull into a catamaran. `units` is required for STL (e.g. `1 mm`). |
| `wigley(L, B, T, {span})` | the exact Wigley hull (no freeboard, so it can't be heeled) |
| `hydrostatics(hull)` | at the design waterline: `length draft volume displacement lcb awp wetted` |
| `float(hull, item, {rho, grid})` | the free-trim hydrostatic float: `sinkage trim draft volume displacement lcb awp wetted` |
| `stability(hull, item, {rho, grid, max_heel, step})` | `gm list max_gz heel_at_max vanishing area_30 area_40 area_total heels gz` |

- **The load:** `item` is any record with `mass` and `cg: [x, y, z]`, such as
  the result of `combine(items)`.
- **Density:** `rho` defaults to sea water at 15 °C (1025.9 kg/m³), as boatmath
  uses.
- **Heels:** default 0° to 180° in 2° steps.

### Frames

A CG is given in the hull file's own frame, in metres:
- **x:** as in the file, bow toward +x.
- **z:** up. hullgeom measures heights from the design waterline, and the
  conversion is `vcg = z − waterline`.
- **y:** from the platform's centreline, which is a monohull's centreplane or
  the midpoint between a catamaran's demihulls.

### Off-centre CGs

hullgeom's loads have no TCG. An off-centre CG shortens the righting arm by
|TCG|·cos φ when heeling toward it, so `stability` reports the curve heeled to
that, the worse, side. This relies on boatmath's hulls being port/starboard
symmetric. The curve's figures are then read from the angle of list, where the
curve crosses up through zero, instead of from upright. With TCG = 0 they are
exactly hullgeom's `summarise`.

## How uncertain loads stay cheap

A float takes about 0.1 s and a GZ curve up to a second, so they can't be
solved once per sample. Instead they are solved on a grid and interpolated
(`surrogate.rs`).

- **The grid:** Chebyshev–Lobatto nodes span each uncertain input's sampled
  range.
  - `float`: mass × LCG, 5 × 5 nodes.
  - `stability`: mass × LCG × VCG, 5 × 5 × 3 nodes. VCG reaches the free-trim
    balance only through the trim, so three nodes capture it.
- **Interpolation:** tensor-product barycentric for each sample. An exact input
  is a single node, so an exact load costs one solve.
- **Kinked outputs:** waterplane area, wetted surface and draft kink as the
  transom immerses, so they are interpolated piecewise-linearly instead, which
  can't ring across a kink. Heels that some node couldn't be balanced at are
  handled locally the same way.
- **Error check:** two samples are re-solved exactly and compared with the
  interpolation. A held-out error above 1% of an output's spread prints a
  warning; pass `{grid: 7}` for more nodes.
- **Graf catamaran:** on the weight study above, the float interpolates to
  within 0.1% of spread. GZ is within 0.6%, about 2 mm, worst past the peak
  where the deck edge goes under.
- **Extreme loads:** at the most extreme corner of a grid (heaviest and highest
  at once) hullgeom may fail to balance large heels. The curve is cut short
  there, and the samples it affects are counted in a warning, with their
  vanishing angle treated as a lower bound.

## Tests

`tests/hulls.rs` checks against:
- **The Wigley float:** V = 4/9·LBT.
- **A box barge's wall-sided GZ curve and list angle:** the box is written as
  STL by the test. hullgeom's ray sectioning is good to about 1% of volume and
  5 mm of GM on its sharp corners.
- **Interpolated against direct solves:** the surrogate's results are compared
  with exact solves.
- **The Graf catamaran against `boatmath statics`:** GM, max GZ, vanishing angle
  and trim agree to the digits shown. This test skips when `graf-demihull.igs`
  isn't at the repository root (or in `$BOATMATH_CAD_DIR`).

## Not yet

- **Wind heeling** and IMO weather criteria.
- **Free-surface corrections** for tanks.
- **Seakeeping and resistance:** strip theory and Michell resistance, through
  the same surrogate.
- **Reading boatmath's JSONL records,** so napkin can sit at the end of a
  pipeline.
- **Moving formulas into napkin:** ITTC-57 friction, the resistance build-up,
  the roll period, and `CaseParams` defaults from a weight study.
