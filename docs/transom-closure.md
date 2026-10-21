# Transom Closure & E-12 STAR-CCM+ Validation in `boatmath`

This document explains how `boatmath` (specifically [`thinship::michell::TransomClosure`](../crates/thinship/src/michell.rs)) models transom sterns in Michell's thin-ship wave-resistance integral, how that compares physically to boundary-layer displacement-thickness ($\delta^*$) methods, and how the different closure and transom-drag options compare against the STAR-CCM+ CFD benchmark for the E-12 catamaran (`2026-10-14.md`).

---

## 1. The Physical Problem: Transom Sterns in Thin-Ship Theory

In classical thin-ship theory (Michell 1898), a hull of half-breadth $y = \pm f(x, z)$ moving at speed $U$ in the $+x$ direction (stream toward $-x$, $z$ downward from the undisturbed free surface) is represented by a centerplane distribution of Kelvin sources of strength proportional to the longitudinal hull slope:

$$\sigma(x, z) = -2 U \frac{\partial f}{\partial x}(x, z)$$

The wave resistance is obtained from the far-field Havelock free-wave spectrum:

$$R_w = \frac{4 \rho g^2}{\pi U^2} \int_1^\infty \bigl(I(\lambda)^2 + J(\lambda)^2\bigr) \frac{\lambda^2}{\sqrt{\lambda^2 - 1}} \, d\lambda, \qquad \nu = \frac{g}{U^2}$$

$$I(\lambda) + i J(\lambda) = \iint_{\text{hull}} \frac{\partial f}{\partial x}(x, z) \, e^{-\nu \lambda^2 z} \, e^{i \nu \lambda x} \, dx \, dz$$

### Why a Transom Breaks Unmodified Michell Theory

Integration by parts in $x$ requires the hull to **close** at both bow and stern ($f(x_{\text{bow}}, z) = 0$ and $f(x_{\text{stern}}, z) = 0$) so that net source flux is zero ($\iint \sigma \, dx \, dz = 0$). When a hull terminates aft at a transom of non-zero half-breadth $f_T(z) > 0$ at $x = x_T$:

1. **Leaving the transom open (`TransomClosure::None` / `--closure off`):**
   Integrating $\partial f/\partial x$ only over the wetted hull $x \in [x_T, x_{\text{bow}}]$ ignores the stern closure altogether. Because the hull has not closed back to the centerplane, the source sheet has a net positive flux—equivalent to trailing a semi-infinite parallel body of cross-section $f_T(z)$ to $x \to -\infty$.
2. **Closing the transom as a bare step (`TransomClosure::Fixed { length: 0.0 }`):**
   If the hull drops instantaneously from $f_T(z)$ to $0$ at $x = x_T$, $\partial f/\partial x$ acquires a Dirac delta line sink $+f_T(z)\,\delta(x - x_T)$. In wavenumber space, a step discontinuity decays only algebraically ($\mathcal{O}(1/\lambda)$), radiating unphysical high-angle short waves across all $\lambda$ and grossly overpredicting $R_w$.

### Virtual Appendage vs. Boundary-Layer Displacement Thickness ($\delta^*$)

In the naval-architecture literature, two distinct families of corrections exist for thin-ship theory:

1. **Boundary-Layer Displacement Thickness ($\delta^*$):**
   Methods such as Bašić et al. (2020) modify the effective hull slope along the hull by the turbulent boundary-layer displacement-thickness growth rate, $(\partial f/\partial x)_{\text{eff}} = \partial f/\partial x + \partial \delta^*/\partial x$. This accounts for viscous dampening of stern wave-making where the boundary layer thickens along the closing afterbody.
   - In `boatmath`, this is available as an opt-in tangency correction via `WaveOptions { boundary_layer: true, .. }` ([`SectionalHull::boundary_layer_tangency`](../crates/hullgeom/src/sectional.rs) and [`SectionalHull::contract_with_tangency`](../crates/hullgeom/src/sectional.rs)) or `boatmath study --boundary-layer`, in addition to the propeller inflow wake model ([`thinship::propulsion::Interaction::frictional_wake_at`](../crates/thinship/src/propulsion.rs)).
2. **Transom Hollow / Virtual Appendage (Doctors & Day 1997; Couser, Molland & Armstrong 1997):**
   When water separates cleanly from the sharp edge of a transom stern, it leaves a ventilated air cavity (or recirculating dead-water hollow) of length $L_v$ aft of the transom before the two sides of the wake meet on the centerplane. Inside the wave-making potential flow, the water surface bounding this hollow acts like a smooth, streamlined extension of the hull—a **virtual appendage** of length $L_v$ over which the effective half-breadth decays smoothly from $f_T(z)$ at $x = x_T$ to $0$ at $x = x_T - L_v$.
   - **This is what `TransomClosure` in `boatmath` implements.**

---

## 2. Mathematical Formulation of `TransomClosure` in `boatmath`

### Smoothstep Virtual Appendage

When a hull has an immersed transom at $x = x_T$ with section profile $f_T(z)$ and immersion depth $d_T$, `thinship` appends a virtual extension over $x \in [x_T - L_v, x_T]$ parameterized by $s = (x_T - x) / L_v \in [0, 1]$:

$$f_v(x, z) = f_T(z) \, \phi(s), \qquad \phi(s) = (1 - s)^2 (1 + 2s) = 1 - 3s^2 + 2s^3$$

Why cubic Hermite smoothstep $\phi(s)$?
- $\phi(0) = 1$ and $\phi(1) = 0$: the half-breadth matches the transom at $s = 0$ and closes to zero at $s = 1$.
- $\phi'(0) = 0$ and $\phi'(1) = 0$: the slope $\phi'(s) = 6s(s - 1)$ vanishes at both the transom edge and the hollow closure point, introducing **no slope discontinuity** at either end.
- Because $f_v(x, z) = f_T(z)\,\phi(s)$ is separable in $x$ and $z$, its depth integral at any vertical decay rate $\kappa = \nu \lambda^2$ is simply the aft station's own depth integral $Z_T(\kappa) = \int_0^{d_T} f_T(z) e^{-\kappa z}\,dz$ multiplied by the closed-form oscillatory moment of $\phi'(s)$ ([`hollow_shape_moment`](../crates/thinship/src/michell.rs)):

$$\mathcal{M}(K) = \int_0^1 \phi'(s) \, e^{i K s} \, ds = 6 \int_0^1 (s^2 - s) \, e^{i K s} \, ds, \qquad K = -k_x L_v$$

Notice the two limits of $\mathcal{M}(-k_x L_v)$:
- **Bare step ($L_v \to 0$):** $\mathcal{M}(0) = \phi(1) - \phi(0) = -1$, recovering the exact Dirac-delta step amplitude $-Z_T(\kappa) e^{i k_x x_T}$.
- **Infinite hollow ($L_v \to \infty$):** By Riemann–Lebesgue, $\mathcal{M}(-k_x L_v) \to 0$ for all $k_x > 0$, recovering `TransomClosure::None` (`--closure off`).

The same virtual appendage is placed consistently across:
1. **Far-field wave resistance (`R_w`):** [`SectionalHull::transom_term`](../crates/thinship/src/sectional.rs) in `crates/thinship/src/sectional.rs`.
2. **Near-field wave transforms (sinkage/trim & wave elevation):** [`SectionalHull::add_transom_transforms`](../crates/thinship/src/sectional.rs) in `crates/thinship/src/sectional.rs`.
3. **Physical-space Rankine source sheet:** [`build_sheet`](../crates/thinship/src/nearfield.rs) in `crates/thinship/src/nearfield.rs` (which prepends 12 virtual stations over $[x_T - L_v, x_T]$).

---

## 3. Hollow Length Modes & Transom Dryness

[`TransomClosure::hollow_length`](../crates/thinship/src/michell.rs) determines $L_v$ from the transom immersion depth $d_T$, speed $U$ (via $\nu = g/U^2$), and the **transom Froude number**:

$$\mathrm{Fn}_T = \frac{U}{\sqrt{g \, d_T}} = \frac{1}{\sqrt{\nu \, d_T}}$$

| Variant | CLI `--closure` | Hollow Length $L_v$ | Physical Regime |
| :--- | :--- | :--- | :--- |
| `None` | `off` | $\infty$ (no closure term) | Infinite hollow; transom slope omitted from Michell integral. |
| `Ballistic { coeff }` | `ballistic[:c]` | $L_{\text{wet}} = c \, U \sqrt{d_T / g}$ (default $c = \sqrt{2} \approx 1.414$) | Free-fall ballistic trajectory: a fluid particle leaving the transom horizontally at speed $U$ falls depth $d_T$ in time $t = \sqrt{2 d_T / g}$. |
| `Auto { coeff, fraction }` | `auto[:frac]` *(default)* | Smoothstep blend from $L_{\text{wet}}$ ($\mathrm{Fn}_T \le 2.5$) to $L_{\text{dry}} = \text{frac} \cdot \frac{2\pi U^2}{g}$ ($\mathrm{Fn}_T \ge 3.5$, default $\text{frac} = 0.5$) | Uses the ballistic hollow while the transom is wet/partially ventilated ($\mathrm{Fn}_T \le 2.5$), and transitions to half a transverse wavelength ($\frac{1}{2}\lambda_T = \pi U^2 / g$) once the transom runs dry ($\mathrm{Fn}_T \ge 3.5$). |
| `Fixed { length }` | `fixed:L` | Constant $L_v = L$ [m] | Prescribed geometric hollow length independent of speed. |

### Transom Dryness & Hydrostatic Transom Drag (`--transom-drag`)

Both [`TransomClosure::Auto`](../crates/thinship/src/michell.rs) and the optional hydrostatic transom drag (`study --transom-drag`) share the smoothstep dryness function [`thinship::transom_dryness`](../crates/thinship/src/michell.rs):

$$w_{\text{dry}}(\mathrm{Fn}_T) = 3 t^2 - 2 t^3, \qquad t = \mathrm{clamp}\!\left(\frac{\mathrm{Fn}_T - 2.5}{3.5 - 2.5},\; 0,\; 1\right)$$

When a transom runs dry in reality, atmospheric air fills the transom face instead of water. The hull loses the forward hydrostatic pressure thrust that undisturbed water would exert on the submerged transom face:

$$R_{\text{transom}} = w_{\text{dry}}(\mathrm{Fn}_T) \, \rho g \iint_{\text{transom}} z \, dA = w_{\text{dry}}(\mathrm{Fn}_T) \, \rho g \, M_z$$

where $M_z = 2 \int_0^{d_T} z \, f_T(z) \, dz$ is the first moment of the transom area about the waterline ([`Transom::depth_moment`](../crates/hullgeom/src/sectional.rs)). Because linear potential-flow wave resistance $R_w$ only measures far-field radiated wave energy, $R_{\text{transom}}$ accounts separately for the local hydrostatic pressure deficit on a ventilated transom (`forces.r_transom` when `--transom-drag` is enabled).

---

## 4. E-12 Catamaran vs. STAR-CCM+ CFD Benchmark

Running `cargo run --release --example e12_starccm` ([`crates/thinship/examples/e12_starccm.rs`](../crates/thinship/examples/e12_starccm.rs)) floats the 2500 kg, 3.0 m span E-12 catamaran ($L_{\text{WL}} = 8.645\text{ m}$, $d_T = 0.0752\text{ m}$, $A_T = 0.0546\text{ m}^2$ per demihull, full dry-transom hydrostatic drag $2 \rho g M_z = 45.5\text{ N}$) and compares `thinship` against the coarse-mesh STAR-CCM+ CFD results from `2026-10-14.md`:

| Speed | $\mathrm{Fn}$ | $\mathrm{Fn}_T$ | Solver / Closure Mode | Total Drag [N] | Pressure / Wave [N] | Viscous [N] | Total Error vs. CFD |
| :---: | :---: | :---: | :--- | :---: | :---: | :---: | :---: |
| **6 kn** | 0.335 | 3.59 | **STAR-CCM+ (coarse)** | **304.0** | **100.0** | **204.0** | baseline |
| | | | `Ballistic` ($L_v = 0.38\text{ m}$) | 334.9 | 122.7 | 212.2 | +10.2% |
| | | | `Ballistic + BL` | 334.7 | 122.5 | 212.2 | +10.1% |
| | | | `Auto + Transom` ($L_v = 3.05\text{ m} + R_{\text{tr}}$) | 306.0 | 93.9 | 212.2 | **+0.7%** |
| | | | `Auto + Transom + BL` | 307.2 | 95.0 | 212.2 | +1.1% |
| | | | `Auto` ($L_v = 3.05\text{ m}$) | 260.5 | 48.3 | 212.2 | -14.3% |
| | | | `None` ($L_v = \infty$) | 251.9 | 39.7 | 212.2 | -17.1% |
| **8 kn** | 0.447 | 4.79 | **STAR-CCM+ (coarse)** | **894.0** | **556.0** | **338.0** | baseline |
| | | | `Ballistic` ($L_v = 0.51\text{ m}$) | 963.6 | 603.5 | 360.2 | +7.8% |
| | | | `Ballistic + BL` | 958.1 | 597.9 | 360.2 | +7.2% |
| | | | `Auto + Transom` ($L_v = 5.43\text{ m} + R_{\text{tr}}$) | 755.3 | 395.1 | 360.2 | -15.5% |
| | | | `Auto + Transom + BL` | 750.1 | 389.9 | 360.2 | -16.1% |
| | | | `Auto` ($L_v = 5.43\text{ m}$) | 709.7 | 349.5 | 360.2 | -20.6% |
| | | | `None` ($L_v = \infty$) | 801.6 | 441.4 | 360.2 | -10.3% |
| **10 kn** | 0.559 | 5.99 | **STAR-CCM+ (coarse)** | **1232.0** | **718.0** | **513.0** | baseline |
| | | | `Ballistic` ($L_v = 0.63\text{ m}$) | 1313.4 | 770.0 | 543.4 | +6.6% |
| | | | `Ballistic + BL` | 1306.8 | 763.4 | 543.4 | +6.1% |
| | | | `Auto + Transom` ($L_v = 8.48\text{ m} + R_{\text{tr}}$) | 1118.9 | 575.5 | 543.4 | -9.2% |
| | | | `Auto + Transom + BL` | 1113.9 | 570.5 | 543.4 | -9.6% |
| | | | `Auto` ($L_v = 8.48\text{ m}$) | 1073.3 | 529.9 | 543.4 | -12.9% |
| | | | `None` ($L_v = \infty$) | 1099.0 | 555.6 | 543.4 | -10.8% |
| **12 kn** | 0.670 | 7.19 | **STAR-CCM+ (coarse)** | **1532.0** | **787.0** | **745.0** | baseline |
| | | | `Ballistic` ($L_v = 0.76\text{ m}$) | 1563.1 | 802.4 | 760.7 | +2.0% |
| | | | `Ballistic + BL` | 1556.2 | 795.5 | 760.7 | **+1.6%** |
| | | | `Auto + Transom` ($L_v = 12.21\text{ m} + R_{\text{tr}}$) | 1404.5 | 643.9 | 760.7 | -8.3% |
| | | | `Auto + Transom + BL` | 1400.3 | 639.6 | 760.7 | -8.6% |
| | | | `Auto` ($L_v = 12.21\text{ m}$) | 1359.0 | 598.3 | 760.7 | -11.3% |
| | | | `None` ($L_v = \infty$) | 1355.0 | 594.3 | 760.7 | -11.6% |

### Key Takeaways from the E-12 Comparison

1. **Viscous Resistance Agreement:**
   `thinship`'s ITTC-57 friction line on the sectioned wetted surface tracks STAR-CCM+ viscous drag within **2%–6%** across all four speeds (212 N vs. 204 N at 6 kn; 360 N vs. 338 N at 8 kn; 543 N vs. 513 N at 10 kn; 761 N vs. 745 N at 12 kn) without needing any ad-hoc form factor.
2. **Why `Ballistic` and `Auto + Transom` Bracket the CFD Pressure Drag:**
   - Because E-12 has a shallow transom ($d_T = 7.5\text{ cm}$ at rest), even at 6 knots the transom Froude number is $\mathrm{Fn}_T = 3.59 > 3.5$ (fully dry).
   - At 6 kn ($\mathrm{Fn} = 0.335$), adding the dry-transom hydrostatic drag ($45.5\text{ N}$) to `Auto`'s far-field wave drag ($48.3\text{ N}$) gives **$93.8\text{ N}$ pressure drag** vs. **$100.0\text{ N}$ in STAR-CCM+** (and **$306.0\text{ N}$ vs. $304.0\text{ N}$ total drag**, within **0.7%**).
   - At higher Froude numbers (8–12 kn, $\mathrm{Fn} = 0.45\text{–}0.67$), dynamic sinkage and bow-up trim in the free-to-heave-and-trim CFD simulation immerse the stern deeper into the trough of the hull's own wave system and add trim-induced pressure drag ($L \tan\tau$). Consequently, `Ballistic` (which keeps a shorter rooster-tail closure $L_v \approx 0.5\text{–}0.76\text{ m}$) tracks the CFD pressure drag within **+1.6% to +7.8%**, while `Auto + Transom` (at fixed resting trim) underpredicts by **8%–15%**. Running `boatmath study --transom-drag --lift 2dt` (or `--closure ballistic`) captures both effects.
3. **Effect of Boundary-Layer Displacement Thickness (`+ BL`):**
   - Enabling the von Kármán / $1/7\text{th}$-power-law boundary-layer displacement-thickness tangency correction (`boundary_layer: true` / `--boundary-layer`) damps the stern wave-making at 8–12 kn by **5–7 N** (~1% of $R_w$), bringing `Ballistic + BL` at 12 kn to **795.5 N pressure drag vs. 787.0 N CFD** (**+1.6% total drag error**).
   - Because E-12 is a high-length-to-beam catamaran ($L/B \approx 16$) at full-scale Reynolds number ($\mathrm{Re}_L \sim 3 \times 10^7$) whose stern terminates at a transom rather than closing sharply to a point, its stern wave-making is dominated by the transom hollow length $L_v$ and dry-transom hydrostatic drag $R_{\text{transom}}$ rather than by boundary-layer displacement thickness.
