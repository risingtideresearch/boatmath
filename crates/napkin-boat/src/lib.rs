//! # napkin-boat
//!
//! napkin with boatmath's hulls: `hull("file.igs", {waterline: …, span: …})`
//! loads a platform from CAD through `hullgeom`, and `float` and `stability`
//! carry a napkin weight study's uncertain mass and CG through the free-trim
//! float and the GZ curve (see [`surrogate`] for how that stays cheap).

pub mod hull;
pub mod surrogate;

/// An interpreter with the hull functions installed.
pub fn interp(samples: usize, seed: u64) -> napkin::Interp {
    let mut it = napkin::Interp::new(samples, seed);
    hull::install(&mut it);
    it
}
