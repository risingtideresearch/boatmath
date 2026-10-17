//! # thinship
//!
//! Thin-ship theory: **wave resistance** (Michell's integral, [`michell`]), **dynamic sinkage and
//! trim**, the far-field **wave pattern**, and **viscous resistance**
//! (ITTC-57) for **sectional hulls** — hulls described by their stations,
//! cut from IGES or STL geometry by [`hullgeom`], which also owns
//! the hydrostatics and the equilibrium solver. This crate is the flow
//! theory on top: [`sectional::dynamic_load_closure`] hands the equilibrium
//! solver its speed-dependent load.
//!
//! ## Theory
//!
//! With `ν = g/U²` (Tuck 1989; Dambrine, Pierre & Rousseaux 2016):
//!
//! ```text
//! R_w = (4 ρ g²)/(π U²) ∫_1^∞ (I² + J²) λ²/√(λ² − 1) dλ
//! I(λ) + i J(λ) = ∬ (∂f/∂x) exp(−ν λ² z) exp(i ν λ x) dx dz
//! ```
//!
//! The depth integral is done per station along its own section curve and
//! interpolated along `x` by a B-spline, whose oscillatory integral is closed
//! form per span ([`sectional`]); only the smooth outer integral is
//! quadratured, with panels sized to the local oscillation rate and refined
//! to a requested tolerance.
//!
//! Viscous resistance uses the ITTC-57 correlation line with an optional form
//! factor and roughness allowance, referencing each hull's wetted shell area.
//!
//! ## Assumptions and limits
//!
//! - Thin-ship (Michell) linearisation: slender hull, `|∂f/∂x| ≪ 1`; no
//!   wave-breaking; deep water; infinite fluid extent. Hulls are
//!   port/starboard symmetric about their centreplanes.
//! - A **transom stern** is closed by a virtual appendage — see
//!   [`TransomClosure`] — whose hollow length defaults to the ballistic
//!   estimate while the transom is wet and half a transverse wavelength once
//!   it runs dry; the bow has no such treatment.
//!
//! ## Example
//!
//! ```
//! use thinship::{sectional, Conditions};
//! use hullgeom::iges::{self, HullPose, Platform, SectionalOptions};
//!
//! let surfaces = iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
//! let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
//! let cut = source
//!     .situate_sectional(
//!         0,
//!         0.0,
//!         &HullPose::default(),
//!         &Platform::default(),
//!         &SectionalOptions::default(),
//!     )
//!     .unwrap()
//!     .expect("the hull is wet");
//! let cond = Conditions::freshwater(3.0); // 3 m/s
//! let r = sectional::multihull_resistance(
//!     &[(&cut.hull, cut.placement)],
//!     &cond,
//!     &Default::default(),
//!     &Default::default(),
//! )
//! .unwrap();
//! assert!(r.wave.resistance > 0.0 && r.viscous_total > 0.0);
//! assert!(r.total > r.wave.resistance);
//! ```

#[cfg(test)]
mod cut_tests;
mod friction;
#[cfg(test)]
mod hull;
#[cfg(test)]
mod hulls;
mod michell;
pub mod nearfield;
pub mod propulsion;
pub mod sectional;
pub mod spectrum;
pub mod squat;
#[cfg(test)]
mod validation;

// The geometry vocabulary this crate's signatures speak.
pub use friction::{
    ittc57_cf, roughness_delta_cf, roughness_reynolds, schlichting_rough_cf,
    viscous_resistance_for, Roughness, ViscousOptions, ViscousResistance,
};
pub use hullgeom::{
    Conditions, Error, Fluid, Placement, Result, SectionalHull, Transom, C64, STANDARD_GRAVITY,
};
pub use michell::{
    transom_dryness, TransomClosure, WaveOptions, WaveResistance, AUTO_DRY_FN_T, AUTO_WET_FN_T,
    BALLISTIC_COEFF, DRY_HOLLOW_FRACTION,
};
pub use sectional::SectionalWave;
pub use spectrum::{FreeWaveSpectrum, WaveGrid};

/// Combined resistance breakdown for a multihull.
#[derive(Debug, Clone)]
pub struct MultihullResistance {
    /// Combined wave resistance (with interference).
    pub wave: WaveResistance,
    /// Per-member viscous resistance, in input order.
    pub viscous: Vec<ViscousResistance>,
    /// Σ of the members' viscous resistances [N].
    pub viscous_total: f64,
    /// Combined wave + viscous resistance [N].
    pub total: f64,
    /// Effective (towing) power P_E = R_t · U [W].
    pub effective_power: f64,
    /// Σ of the members' wetted surfaces [m²] (reference area for Cw/Cv/Ct).
    pub wetted_surface: f64,
    pub cw: f64,
    pub cv: f64,
    pub ct: f64,
    /// Σ of each member's standalone wave resistance [N].
    pub solo_wave_total: f64,
    /// Wave interference factor: combined R_w / Σ standalone R_w
    /// (< 1 favourable, > 1 unfavourable, 1 when far apart).
    pub interference: f64,
}

/// Test support: the CAD files the real-hull tests run on live at the
/// repository root and are not committed (`.gitignore`: `/*.igs`), so a test
/// that needs one skips, saying so, when it is absent.
#[cfg(test)]
pub(crate) fn cad_fixture(name: &str) -> Option<String> {
    let path = format!("{}/../../{name}", env!("CARGO_MANIFEST_DIR"));
    match std::fs::read_to_string(&path) {
        Ok(text) => Some(text),
        Err(_) => {
            eprintln!("skipping: CAD fixture {name} is not present (not in the repository)");
            None
        }
    }
}
