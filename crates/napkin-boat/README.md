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
| `shell(hull, {x: [a, b], y: [..], z: [..], hull})` | `area centroid` of the shell plating inside a box (each axis optional) |
| `slice(hull, x, {below, hull})` | the transverse section at `x`: `area centroid girth outlines`, cut off at height `below` if given |

- **The load:** `item` is any record with `mass` and `cg: [x, y, z]`, such as
  the result of `combine(items)`.
- **Density:** `rho` defaults to sea water at 15 °C (1025.9 kg/m³), as boatmath
  uses.
- **Heels:** default 0° to 180° in 2° steps.
- **`hydrostatics`** also gives `waterline`, the design waterline's height in
  the file frame, for writing regions relative to it.
- **`hull: k`** picks one hull of the platform, numbered from 0. For a
  catamaran, hull 0 is the one at +y.

## Weights from geometry

`shell` and `slice` measure the hull's surface, so a weight can come from the
shape rather than a guess. Each returns a centroid, so an item is
`{mass: area × areal weight, cg: centroid}`.
`examples/graf_structure.nap` builds the Graf catamaran's structure this way,
then floats and heels it:

```
let bottom = shell(cat, {z: [-1 m, wl + (100 ± 20 mm)]})       # plating below an uncertain line
fn bulkhead(name, x, areal) {
  let s = slice(cat, x)                                        # both demihulls, keel to sheer
  {name: name, mass: s.area * areal, cg: s.centroid}
}
let frames = map(linspace(0.4 m, 6.0 m, 15), |x| {
  let s = slice(cat, x)
  {mass: s.girth * (0.9 kg/m ± 10%), cg: s.centroid}           # girth: shell length around
})
let tank = trapz(xs, map(xs, |x| slice(cat, x, {below: wl}).area))   # a volume, by slices
```

- **What's measured:** both work on the triangle mesh hullgeom keeps of each
  hull. For STL that's the file's own facets; for IGES the patches are
  tessellated at about 1/400 of the hull's size. A one-sided (half-hull) file
  is mirrored.
- **Accuracy:** the facets are flat, so on curved plating the mesh falls short
  of the surface by roughly (h/R)²/8. On the Graf demihull:
  - the shell below the waterline is 5.777 m², against hullgeom's wetted surface
    of 5.794 m² (−0.3%);
  - 131 slices integrated along x give 0.7476 m³, against its displaced volume
    of 0.7505 m³ (−0.4%).

  On flat-sided shapes (a box) and the Wigley form, areas, centroids and girths
  are exact.
- **Open sections:** a section without a deck is closed straight across its
  top, as hullgeom's GZ does. `girth` counts only the shell, not that closing
  line or the cut at `below`.
- **Uncertain bounds:** region bounds, a station or a `below` height may be
  uncertain. They go through the same grid as `float`, and the held-out check
  ignores errors within the mesh's own 0.3%.
- **`outlines`:** each closed loop as `[x, y, z]` points, ready for napkin's
  `area`, `centroid` and `section`. Given only when the station is certain.

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
- **Shell and slices:** a box's shell and slices exactly, the Wigley form's
  section areas and centroids, and the Graf demihull's slices and shell against
  hullgeom's volume and wetted surface.
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
