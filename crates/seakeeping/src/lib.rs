//! # seakeeping
//!
//! Linear **strip-theory seakeeping** for sectional hulls: heave and pitch in
//! regular waves, in the frequency domain, following Salvesen, Tuck &
//! Faltinsen (1970). It is built on [`hullgeom`] alone — the
//! thin-ship wave kernel (`thinship`) plays no part in it; the two share the
//! hull geometry and the equilibrium attitude it is linearised about.
//!
//! ## Conventions
//!
//! - Hull frame as in [`hullgeom`]: `x` forward (the bow is the
//!   high-`x` end), sections measure `z` **down** from the waterline.
//! - Heave `η₃` is positive **up**; pitch `η₅` is positive **bow up**, about a
//!   pivot `x_ref` (usually the longitudinal centre of gravity).
//! - Waves travel at heading `β` to the ship's course: `β = 0` following
//!   seas, `β = π` head seas. The incident elevation, in the frame moving
//!   with the ship at speed `U`, is `ζ = Re[ζ_a e^{i k (x − x_ref) cos β}
//!   e^{−i ω_e t}]` with the encounter frequency `ω_e = ω − k U cos β`
//!   ([`waves`]): complex amplitudes here are per unit wave amplitude, with
//!   their phase relative to a crest at `x_ref`.
//!
//! ## What is here
//!
//! - [`restoring`]: hydrostatic restoring coefficients `C₃₃, C₃₅, C₅₅`
//!   from the waterplane;
//! - [`froude_krylov`]: the incident-wave heave force and pitch moment, in
//!   closed form per x-span from the sectional hull's depth integrals
//!   ([`hullgeom::SectionalHull::x_transform`]);
//! - [`green`], [`section2d`]: the 2-D pulsating source and Frank's
//!   close-fit method on each station's section curve — added mass,
//!   damping, radiated waves, and the diffraction force by Haskind;
//! - [`green`], [`section2d`] also solve each section's **sway and roll**
//!   (the antisymmetric problems) with their diffraction and
//!   Froude–Krylov forces;
//! - [`strip`]: the platform's five-mode response — sway, heave, roll,
//!   pitch and yaw about G — in regular waves at forward speed: each
//!   station's motion follows from the platform's through a kinematic map
//!   (hull offsets and G's height included, so multihulls, and asymmetric
//!   ones such as proas, need nothing special), and its force is the
//!   forward-speed operator round its complex added mass, which gives
//!   Salvesen–Tuck–Faltinsen's coefficients with their transom terms;
//! - [`strip::added_resistance`]: mean added resistance by Gerritsma &
//!   Beukelman's radiated energy (no short-wave correction), and
//!   [`strip::added_resistance_maruo`] by Maruo's far-field momentum,
//!   from the Kochin function of the stations' sources;
//! - [`sea`]: Bretschneider and JONSWAP spectra, significant motions and
//!   accelerations, mean added resistance in a sea.
//!
//! Roll damping is potential flow's alone, which leaves a monohull's roll
//! resonance several times too high (the real damping is viscous);
//! [`strip::StripOptions::roll_damping`] adds a fraction of critical in its
//! place. A multihull's roll is mostly its hulls' heave, and does not have
//! this problem. The lateral modes are checked by long beam waves (the hull
//! follows the water: sway 1, roll the slope), zero-speed reciprocity of the
//! five-mode coefficients, and mirror symmetry; and, section by section,
//! against Vugts' (1970) cylinder experiments in beam waves as plotted in
//! Journée's SEAWAY validation report: sway, heave and roll added mass,
//! damping and wave loads match SEAWAY's curves to about 1–6% and the tank
//! to about 1–12% (roll's small added mass and damping, where viscosity
//! matters, 15–40%). The three-dimensional lateral response has no
//! experimental check yet.
//!
//! Still to come: a short-wave added-resistance correction; empirical
//! (Ikeda-type) roll damping.
//!
//! Checks: the source's principal value against quadrature, the section
//! damping against the energy its waves carry, the exact infinite-frequency
//! added mass of a semicircle (and the approach to it), the Haskind diffraction
//! force against the solved diffraction problem, long-wave limits of the
//! full response.
//!
//! Against experiment — Journée's four Wigley hulls in head waves (DUT-SHL
//! Report 0909, 1992; `validation` tests, `journee_wigley` example): the
//! heave added mass and damping agree to ~10–15% over the mid frequencies,
//! the wave force and moment on the restrained hull to ~5–15%, and heave in
//! head waves closely, resonance peaks included; zero-speed pitch comes
//! out 10–20% low. At speed the pitch added inertia is ~30% low and the
//! pitch damping grows with U² where the tank shows none, so the pitch
//! resonance falls at shorter waves than measured, and the radiated-energy
//! added resistance (which goes with the square of the motions) peaks early
//! and two to six times too high at Fn 0.3–0.4. The far-field (Maruo)
//! estimate, in which the sections' waves interfere and the forward-scattered
//! wave carries no momentum, is within a factor of two of the measured peaks
//! there but about half the measurements at Fn 0.2 and on the beamy L/B 5
//! hulls, and neither method ranks the four hulls reliably. Journée found the same discrepancies with his own
//! Frank- and Ursell-based strip codes (MARIND 2001); they are strip
//! theory's, on these hulls.

#![allow(
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::neg_cmp_op_on_partial_ord
)]

pub mod froude_krylov;
pub mod green;
pub mod linalg;
pub mod platform;
pub mod restoring;
pub mod sea;
pub mod section2d;
pub mod strip;
#[cfg(test)]
mod validation;
pub mod waves;

pub use froude_krylov::{froude_krylov, WaveLoad};
pub use restoring::{restoring, Restoring};
