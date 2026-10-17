# boatmath

Physics on boats, in Rust. Hulls from CAD (IGES, STL) or from
[camber](https://github.com/risingtideresearch/camber) are cut into
sections, then floated, heeled, run through calm water and through waves,
and fitted with propellers and motors. The work is done by command-line
tools that pass JSON records down a pipe (`boatmath`), and by a web app
with a queue (`boatmath-web`). Both compute through the same library
(`boatmath`).

```sh
boatmath hull e12.igs --waterline -0.95 \
  | boatmath case --span 3.0 \
  | boatmath study --froude 0.3:0.7:0.05 \
  | boatmath run \
  | boatmath plot -x speed -y forces.rt -o rt.svg
```

## Workspace

| crate | what it holds | depends on |
|---|---|---|
| `hullgeom` | hulls cut from IGES/STL patches or triangles into sections (`SectionalHull`), hydrostatics, hydrostatic and dynamic equilibrium (`float`), GZ curves (`stability`), B-splines, closed-form moments, `Conditions` | — |
| `thinship` | thin-ship theory: wave resistance by Michell's integral (`thinship::michell`), the near field and dynamic sinkage and trim (`squat`), the far-field wave spectrum, the propeller's wake and thrust deduction (`propulsion`), ITTC-57 friction | `hullgeom` |
| `seakeeping` | linear strip-theory seakeeping: the five rigid-body modes, added resistance, irregular seas | `hullgeom` |
| `propeller` | B-series propeller search and motor matching (a port of propopt's web cores, with its motor database) | — |
| `boatmath` | what the CLI and the web app share: hulls, cases and studies; the calm-water and wave computations on a platform; native JSON geometry; camber documents (`boatmath::camber`); drives on mounts; CAD export | all of the above |
| `boatmath-cli` | the `boatmath` command: Unix-style tools on JSON-lines streams | `boatmath`, `propeller` |
| `boatmath-web` | the browser front end, a store and a queue of studies | `boatmath` |

The geometry crate knows no flow theory. `thinship` gives the equilibrium
solver its speed-dependent load through `float::DynamicModel`
(`thinship::sectional::dynamic_load_closure`). The seakeeping crate builds
on the sectional hull's transforms without going through `thinship`.

## Where to read more

- **The CLI:** [`docs/boatmath-cli.md`](docs/boatmath-cli.md) covers every
  record and command, mounts and drives, propellers, CAD export and camber
  import.
- **The theory:** each crate's module docs (`cargo doc --open`). `thinship`'s
  docs have the Michell integral and how it's evaluated on sections;
  [`docs/michell-calculation.tex`](docs/michell-calculation.tex) has
  derivation notes, written for the earlier tensor-product B-spline hulls. [`docs/dynamic-squat-derivation.md`](docs/dynamic-squat-derivation.md)
  derives the dynamic sinkage and trim.
- **Seakeeping:** [`docs/seakeeping-findings.md`](docs/seakeeping-findings.md)
  covers what was built and how it validates against Journée's Wigley hulls
  and Vugts' cylinders.
- **Deployment:** [`deploy/README.md`](deploy/README.md) runs the web app as
  a service.

## What to trust

- **Calm water** is thin-ship (Michell) theory: a slender hull
  (`|∂f/∂x| ≪ 1`), no wave-breaking, deep water, an infinite fluid.
  Multihulls interfere exactly within that theory. A transom stern is
  closed by a virtual appendage (`TransomClosure`) whose hollow length is
  a modelling choice (by default the ballistic hollow while the transom is
  wet, half a transverse wavelength once it runs dry), so transom-sterned
  results carry that uncertainty.
  Thin-ship theory overstates `C_W` on full hulls.
- **Viscous resistance** is the ITTC-57 line on each hull's wetted area,
  times a form factor, plus a roughness allowance, both 0 unless a study
  sets them. Neither is derived from the hull.
- **A dry transom's drag** (the hydrostatic thrust it loses) is outside
  thin-ship theory; `study --transom-drag` adds it. Don't back the form factor out of a measured `C_T` using this
  `C_W`: the fit would absorb thin-ship theory's error.
- **Seakeeping** is strip theory: trust heave; treat pitch at speed and
  added resistance near resonance as indicative (see the findings doc).
  There's no hull-to-hull wave interaction in waves.
- **Wake and thrust deduction** come from the thin-ship singularities, so
  they miss the frictional wake and the suction on a flat run or transom.
  Behind a transom, read `t` as a lower bound.

## Web front end

`boatmath-web` keeps a permanent record at three levels and works through a
queue:

- **Hulls** (`/hulls`): an uploaded IGES or STL file with how it is cut (its
  design waterline, stations and rays, units), shown in 3D with its
  sections and hydrostatics.
- **Cases** (`/cases`): a platform on a hull and its load: the hull on its
  own or doubled into a catamaran at a span, its mass (carried by sinking,
  or by scaling the hull), LCG, VCG, radii of gyration and roll damping. Its
  statics are computed when it's made: the float at rest, GM_T and the roll
  period, and the GZ curve with its peak, angle of vanishing stability and
  areas.
- **Studies** (`/studies`): a speed on a case, in calm water (the near-field
  pressure, the free surface, resistance, sinkage and trim at speed) or in
  waves from one heading (responses over a wavelength sweep, added
  resistance, an optional irregular sea, an animated seaway).

**Queue** (`/queue`) shows the running study and what waits; **Plot**
(`/plot`) plots any result against any parameter and exports CSV. Every
result records the solver version it was computed with (the last commit to
touch the solver's code). A result from another version is marked stale
and can be run again.

```text
cargo run --release -p boatmath-web -- --data boatmath-data   # http://127.0.0.1:8080/
BOATMATH_WEB_DIR=crates/boatmath-web/src/web boatmath-web     # pages read from disk, to edit them live
```

## References

- J. H. Michell, *The wave resistance of a ship*, Phil. Mag. 45 (1898).
- E. O. Tuck, *The wave resistance formula of J.H. Michell (1898) and its
  significance to recent research in ship hydrodynamics*, J. Austral. Math.
  Soc. B 30 (1989).
- E. O. Tuck, D. C. Scullen & L. Lazauskas, *Ship-wave patterns in the
  spirit of Michell*, IUTAM Symposium (2001).
- J. Dambrine, M. Pierre, G. Rousseaux, *A theoretical and numerical
  determination of optimal ship forms based on Michell's wave resistance*,
  ESAIM: COCV (2016), arXiv:1410.2800.
- N. Salvesen, E. O. Tuck & O. Faltinsen, *Ship motions and sea loads*,
  Trans. SNAME 78 (1970).
- J. M. J. Journée, *Experiments and calculations on four Wigley hullforms*,
  Delft report 0909 (1992).
- ITTC Recommended Procedures: *1957 ITTC Performance Prediction Method*.
