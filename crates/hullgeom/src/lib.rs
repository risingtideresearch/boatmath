//! # hullgeom
//!
//! Hull **geometry** for the hydrodynamic crates of this workspace: hulls cut
//! from CAD into **sections**, their hydrostatics, and the platform's
//! hydrostatic equilibrium. It knows nothing of any flow theory — the
//! thin-ship wave kernel (`thinship`) and strip-theory seakeeping
//! (`seakeeping`) are both built on it.
//!
//! ## Where hulls come from
//!
//! A [`SectionalHull`] is cut from source geometry kept in its own frame, so
//! it can be re-posed and re-cut as often as a sweep or an equilibrium solve
//! needs ([`source::HullSource`]):
//!
//! - **IGES** B-spline surfaces ([`iges::source_fleet`], or surfaces already
//!   in hand through [`iges::source_fleet_from_surfaces`] — e.g. the exact
//!   Wigley test hull, [`iges::wigley_surfaces`]),
//! - **STL** triangle meshes ([`stl::mesh_fleet`]).
//!
//! Frames: CAD sources are `x` along the hull, `y` transverse, `z` **up**;
//! a cut hull's sections measure `z` downward from the effective waterline.
//! All quantities are SI. Hulls are port/starboard symmetric about their
//! centreplanes; the bow is the high-`x` end.
//!
//! ## What a sectional hull offers
//!
//! - hydrostatics: displaced volume, centre of buoyancy, waterplane area and
//!   its first and second moments, wetted surface, transom;
//! - the per-station section curves ([`SectionalHull::curves`]);
//! - the depth-decayed section integrals `Z(x; κ) = ∫ f(x, z) e^{−κz} dz`
//!   interpolated along `x` in closed form per span
//!   ([`SectionalHull::contract`]), and their Fourier transforms along `x`
//!   ([`SectionalHull::x_transform`]) — the building block of every
//!   wavenumber-domain hull integral.
//!
//! ## Equilibrium
//!
//! [`float`] floats a fleet of hulls to carry a load, re-cutting each hull at
//! every trial pose; a hydrodynamic crate can fold a speed-dependent load
//! into the same solve through [`float::DynamicModel`].

#![allow(
    clippy::type_complexity,
    clippy::needless_range_loop,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::suspicious_arithmetic_impl
)]

pub mod bspline;
mod conditions;
mod error;
pub mod float;
pub mod hulls;
pub mod iges;
pub mod moments;
pub mod parallel;
pub mod quadrature;
pub mod sectional;
pub mod source;
pub mod stability;
pub mod stl;

pub use conditions::{Conditions, Fluid, STANDARD_GRAVITY};
pub use error::{Error, Result};
pub use moments::C64;
pub use sectional::{SectionalHull, Transom};

/// Position of one hull of a multihull, in the global (fleet) frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Placement {
    /// Longitudinal shift **added** to the hull's own x coordinates [m]
    /// (0 keeps the coordinates from the hull's file).
    pub x: f64,
    /// Transverse position of the hull's centerplane [m].
    pub y: f64,
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
