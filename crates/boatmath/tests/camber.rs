//! The port against camber: each `tests/camber/<name>.golden.json` (made by
//! `tests/camber/make.mts` from camber's own sweep) holds the trimmed
//! half-sections of `<name>.json` at a spread of the plan's parameter, the
//! hull's ends, and camber's hydrostatics; boatmath must give the same.

use boatmath::camber::{Document, Hull};
use boatmath::{hull_summary, loft, LoftRequest};
use serde_json::Value;

const NAMES: [&str; 5] = ["default", "cruiser", "flat-bottom", "inverted-bow", "graf"];

fn load(file: &str) -> Value {
    let path = format!("{}/tests/camber/{file}", env!("CARGO_MANIFEST_DIR"));
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn doc(name: &str) -> Document {
    Document::parse(&load(&format!("{name}.json"))).unwrap()
}

/// The sweep, section by section, to bisection precision.
#[test]
fn sections_are_cambers() {
    for name in NAMES {
        let g = load(&format!("{name}.golden.json"));
        let hull = Hull::new(&doc(name));
        let r = g["R"].as_u64().unwrap() as usize;
        let mut worst = 0.0f64;
        for s in g["sections"].as_array().unwrap() {
            let u = s["u"].as_f64().unwrap();
            let mine = hull.swept_section(u, r);
            assert_eq!(
                mine.is_none(),
                s["empty"].as_bool().unwrap(),
                "{name} u {u}: empty"
            );
            let Some((pts, keel)) = mine else { continue };
            assert_eq!(keel, s["keel"].as_bool().unwrap(), "{name} u {u}: keel");
            let theirs = s["pts"].as_array().unwrap();
            assert_eq!(pts.len(), theirs.len(), "{name} u {u}: rows");
            for (p, q) in pts.iter().zip(theirs) {
                for c in 0..3 {
                    worst = worst.max((p[c] - q[c].as_f64().unwrap()).abs());
                }
            }
        }
        // Lengths in the document's unit (mm for most, inches for graf),
        // against a hull of a few hundred to several thousand.
        assert!(worst < 1e-6, "{name}: off by {worst:e}");
        for (k, f) in [
            ("aft_limit", hull.aft_limit()),
            ("forward_limit", hull.forward_limit()),
        ] {
            let want = g[k].as_f64().unwrap();
            assert!((f - want).abs() < 1e-12, "{name}: {k} {f}, camber's {want}");
        }
    }
}

/// The patches cut by boatmath's solver float as camber's hull does: its
/// displaced volume, and its waterline's length and beam, at the design
/// waterline. camber's figures are converged (`make.mts` samples finely);
/// its default sampling reads a little low. Only a hull camber calls closed
/// counts: one whose sections end on a transom (camber leaves the transom
/// out of its volume), or whose sheer is under water, isn't a fair test.
/// flat-bottom's open stern is a shallow V across its fanned end section,
/// 1.2% of its length, which the solver's sections (square to x) end with
/// a whole section near the V's point: so its volume is held to 0.5% and
/// its length to 1.5%. A hull of e12's size whose end is a transom (rvh,
/// alumc) comes within 0.07% of camber's volume.
#[test]
fn the_patches_displace_what_camber_does() {
    for name in ["flat-bottom"] {
        let d = doc(name);
        let h = &load(&format!("{name}.golden.json"))["hydro"];
        assert!(h["closed"].as_bool().unwrap());
        let geometry = boatmath::camber::geometry(&d, None, None).unwrap();
        let bytes = serde_json::to_vec(&geometry).unwrap();
        let sections = loft("geometry.json", bytes, &LoftRequest::default()).unwrap();
        let s = &hull_summary(&sections)["hulls"][0];
        let (m, m3) = (d.unit_m, d.unit_m.powi(3));
        let near = |k: &str, mine: f64, theirs: f64, tol: f64| {
            let e = (mine / theirs - 1.0).abs();
            assert!(e < tol, "{name}: {k} {mine}, camber's {theirs} ({e:.4})");
        };
        near(
            "volume",
            s["displaced_volume"].as_f64().unwrap(),
            h["vol"].as_f64().unwrap() * m3,
            0.005,
        );
        {
            near(
                "length",
                s["length"].as_f64().unwrap(),
                h["lwl"].as_f64().unwrap() * m,
                0.015,
            );
            near(
                "beam",
                s["beam"].as_f64().unwrap(),
                h["bwl"].as_f64().unwrap() * m,
                0.01,
            );
        }
    }
}

#[test]
fn a_version_1_document_is_refused() {
    let v: Value = serde_json::json!({ "length": 1000, "sheerPlan": [] });
    let e = Document::parse(&v).err().unwrap();
    assert!(e.contains("version 1"), "{e}");
}

/// Every patch lies on one side of the centreplane, and within the hull's
/// length: the solver finds a
/// hull's centreplane from the patches either side of it, one crossing per
/// patch, and a patch spanning it reads as a one-sided hull, which can't be
/// cut once it's off the centreline (a catamaran's demihull).
#[test]
fn no_patch_spans_the_centreplane() {
    for name in NAMES {
        let d = doc(name);
        let v = load(&format!("{name}.json"));
        let plan = v["sheerPlan"].as_array().unwrap();
        let bow = plan[plan.len() - 1]["x"].as_f64().unwrap() * d.unit_m;
        // Interpolated surfaces stray a little between their sections:
        // 0.05% of the length, against a camber hull's own.
        let tol = 5e-4 * bow;
        let g = boatmath::camber::geometry(&d, None, None).unwrap();
        for p in g["hulls"][0]["patches"].as_array().unwrap() {
            let s = boatmath::native::patch_from_value(p).unwrap();
            let ((u0, u1), (v0, v1)) = (s.u_domain(), s.v_domain());
            let n = 60;
            let pts: Vec<[f64; 3]> = (0..=n)
                .flat_map(|i| (0..=n).map(move |j| (i, j)))
                .map(|(i, j)| {
                    let u = u0 + (u1 - u0) * i as f64 / n as f64;
                    let v = v0 + (v1 - v0) * j as f64 / n as f64;
                    s.point(u, v)
                })
                .collect();
            // Nor does it overshoot the bow (the deck is level here).
            let x_max = pts.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max);
            assert!(
                x_max <= bow + tol,
                "{name}: a patch reaches x {x_max}, past the bow {bow}"
            );
            let ys: Vec<f64> = pts.iter().map(|q| q[1]).collect();
            let lo = ys.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            assert!(
                lo >= -tol || hi <= tol,
                "{name}: a patch spans y {lo}..{hi}"
            );
        }
    }
}
