//! Hull geometry as JSON: the patches or triangles the solver re-poses and
//! cuts, already split into hulls, in metres, with the design waterline at
//! z = 0.
//!
//! ```json
//! { "kind": "nurbs",
//!   "hulls": [ { "patches": [ { "degree": [3, 3], "n_ctrl": [8, 6],
//!                               "knots_u": [...], "knots_v": [...],
//!                               "ctrl": [[x, y, z], ...], "trim_uv": null } ] } ] }
//! { "kind": "mesh",
//!   "hulls": [ { "vertices": [[x, y, z], ...], "triangles": [[0, 1, 2], ...] } ] }
//! ```
//!
//! x and y are the file's; z is up. Patches are polynomial B-splines, so
//! they carry no weights; `ctrl` is row-major, v fastest. An IGES file
//! becomes `nurbs`, an STL `mesh`. Its import settings (waterline, units,
//! scale) are applied once, here, and are not needed again. A trim, a
//! sinkage or another scale maps the control points (or vertices) exactly,
//! so the solver re-poses this as it would the file.

use crate::files::{open_source_bytes, Kind, LoadSettings, SourceFile};
use crate::LoftRequest;
use hullgeom::iges::{HullPose, NurbsSurface3, Platform, SourceFleet};
use hullgeom::stl::MeshFleet;
use serde_json::{json, Value};

/// The file's geometry, its import settings applied, as JSON.
pub fn from_file(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Value, String> {
    let settings = LoadSettings {
        waterline_z: req.waterline.unwrap_or(0.0),
        waterline: req.waterline,
        units: req.units,
        ..LoadSettings::default()
    };
    if let Some(doc) = camber_document(&bytes) {
        let doc = crate::camber::Document::parse(&doc?).map_err(|e| format!("{name}: {e}"))?;
        let g = crate::camber::geometry(&doc, req.waterline, req.units)?;
        let pose = req.pose();
        return if pose.scale == 1.0 && pose.scale_yz == 1.0 {
            Ok(g)
        } else {
            posed(&g, &pose)
        };
    }
    if is_native(&bytes) {
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| format!("{name}: {e}"))?;
        if v.get("kind").is_some() && v.get("hulls").is_some() {
            let pose = req.pose();
            return if pose.scale == 1.0 && pose.scale_yz == 1.0 {
                Ok(v)
            } else {
                posed(&v, &pose)
            };
        }
    }
    let wl = settings.waterline_z;
    let pose = req.pose();
    let at_rest = Platform::default();
    match open_source_bytes(name, bytes.clone(), &settings)?.kind {
        Kind::Iges => {
            let text = String::from_utf8(bytes).map_err(|e| e.to_string())?;
            let fleet =
                hullgeom::iges::source_fleet(&text, wl).map_err(|e| format!("{name}: {e}"))?;
            let hulls = (0..fleet.len())
                .map(|i| {
                    let mut surfs = fleet
                        .posed_surfaces(i, wl, &pose, &at_rest)
                        .map_err(|e| e.to_string())?;
                    for p in surfs.iter_mut().flat_map(|s| s.ctrl.iter_mut()) {
                        p[2] -= wl;
                    }
                    Ok(json!({ "patches": surfs.iter().map(patch_json).collect::<Vec<_>>() }))
                })
                .collect::<Result<Vec<_>, String>>()?;
            Ok(json!({ "kind": "nurbs", "hulls": hulls }))
        }
        Kind::Stl => {
            let units = req.units.ok_or("an STL needs its units")?;
            let fleet =
                hullgeom::stl::mesh_fleet(&bytes, units, wl).map_err(|e| format!("{name}: {e}"))?;
            let hulls = (0..fleet.len())
                .map(|i| {
                    // Posed about the waterline, z up from it.
                    let (v, t) = fleet
                        .posed_tessellation(i, wl, &pose, &at_rest)
                        .map_err(|e| e.to_string())?;
                    Ok(json!({ "vertices": v, "triangles": t }))
                })
                .collect::<Result<Vec<_>, String>>()?;
            Ok(json!({ "kind": "mesh", "hulls": hulls }))
        }
    }
}

/// A camber hull document, if these bytes are one (a JSON object with a
/// sheer plan).
fn camber_document(bytes: &[u8]) -> Option<Result<Value, String>> {
    if !is_native(bytes) {
        return None;
    }
    match serde_json::from_slice::<Value>(bytes) {
        Ok(v) if crate::camber::is_document(&v) => Some(Ok(v)),
        Ok(_) => None,
        Err(e) => Some(Err(format!("not JSON: {e}"))),
    }
}

/// A patch as the geometry's JSON writes it.
pub fn patch_value(s: &NurbsSurface3) -> Value {
    patch_json(s)
}

fn patch_json(s: &NurbsSurface3) -> Value {
    json!({
        "degree": [s.degree_u, s.degree_v],
        "n_ctrl": [s.n_ctrl_u, s.n_ctrl_v],
        "knots_u": s.knots_u,
        "knots_v": s.knots_v,
        "ctrl": s.ctrl,
        "trim_uv": s.trim_uv,
    })
}

fn field<T: serde::de::DeserializeOwned>(v: &Value, k: &str, at: &str) -> Result<T, String> {
    serde_json::from_value(v[k].clone()).map_err(|e| format!("{at}: {k}: {e}"))
}

/// A patch from the geometry's JSON.
pub fn patch_from_value(v: &Value) -> Result<NurbsSurface3, String> {
    patch(v, "patch")
}

fn patch(v: &Value, at: &str) -> Result<NurbsSurface3, String> {
    let [degree_u, degree_v]: [usize; 2] = field(v, "degree", at)?;
    let [n_ctrl_u, n_ctrl_v]: [usize; 2] = field(v, "n_ctrl", at)?;
    let ctrl: Vec<[f64; 3]> = field(v, "ctrl", at)?;
    let s = NurbsSurface3 {
        degree_u,
        degree_v,
        knots_u: field(v, "knots_u", at)?,
        knots_v: field(v, "knots_v", at)?,
        n_ctrl_u,
        n_ctrl_v,
        weights: vec![1.0; ctrl.len()],
        ctrl,
        trim_uv: field(v, "trim_uv", at)?,
    };
    s.validate().map_err(|e| format!("{at}: {e}"))?;
    Ok(s)
}

/// Whether these bytes are a JSON geometry rather than a hull file.
pub fn is_native(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .find(|b| !b.is_ascii_whitespace())
        .is_some_and(|&b| b == b'{')
}

/// A JSON geometry, ready to cut: its design waterline is at z = 0.
pub fn open(bytes: &[u8]) -> Result<SourceFile, String> {
    let v: Value = serde_json::from_slice(bytes).map_err(|e| format!("geometry: {e}"))?;
    open_value(&v)
}

enum Fleet {
    Nurbs(SourceFleet),
    Mesh(MeshFleet),
}

fn fleet(v: &Value) -> Result<Fleet, String> {
    let hulls = v["hulls"].as_array().ok_or("geometry: no hulls")?;
    match v["kind"].as_str() {
        Some("nurbs") => {
            let hulls = hulls
                .iter()
                .enumerate()
                .map(|(i, h)| {
                    h["patches"]
                        .as_array()
                        .ok_or_else(|| format!("geometry hull {}: no patches", i + 1))?
                        .iter()
                        .enumerate()
                        .map(|(j, p)| patch(p, &format!("hull {} patch {}", i + 1, j + 1)))
                        .collect::<Result<Vec<_>, String>>()
                })
                .collect::<Result<Vec<_>, String>>()?;
            Ok(Fleet::Nurbs(
                SourceFleet::from_hulls(hulls).map_err(|e| e.to_string())?,
            ))
        }
        Some("mesh") => {
            let meshes = hulls
                .iter()
                .enumerate()
                .map(|(i, h)| {
                    let at = format!("hull {}", i + 1);
                    Ok((field(h, "vertices", &at)?, field(h, "triangles", &at)?))
                })
                .collect::<Result<Vec<_>, String>>()?;
            Ok(Fleet::Mesh(
                MeshFleet::from_meshes(meshes).map_err(|e| e.to_string())?,
            ))
        }
        other => Err(format!(
            "geometry kind {other:?}: expected \"nurbs\" or \"mesh\""
        )),
    }
}

/// The geometry with a pose applied to every hull — a scale, say, about the
/// design waterline — as a new geometry.
pub fn posed(v: &Value, pose: &HullPose) -> Result<Value, String> {
    let at_rest = Platform::default();
    Ok(match fleet(v)? {
        Fleet::Nurbs(f) => {
            let hulls = (0..f.len())
                .map(|i| {
                    let surfs = f
                        .posed_surfaces(i, 0.0, pose, &at_rest)
                        .map_err(|e| e.to_string())?;
                    Ok(json!({ "patches": surfs.iter().map(patch_json).collect::<Vec<_>>() }))
                })
                .collect::<Result<Vec<_>, String>>()?;
            json!({ "kind": "nurbs", "hulls": hulls })
        }
        Fleet::Mesh(f) => {
            let hulls = (0..f.len())
                .map(|i| {
                    let (v, t) = f
                        .posed_tessellation(i, 0.0, pose, &at_rest)
                        .map_err(|e| e.to_string())?;
                    Ok(json!({ "vertices": v, "triangles": t }))
                })
                .collect::<Result<Vec<_>, String>>()?;
            json!({ "kind": "mesh", "hulls": hulls })
        }
    })
}

pub fn open_value(v: &Value) -> Result<SourceFile, String> {
    if crate::camber::is_document(v) {
        let doc = crate::camber::Document::parse(v).map_err(|e| format!("geometry: {e}"))?;
        let g = crate::camber::geometry(&doc, None, None)?;
        return open_value(&g);
    }
    let (kind, source): (Kind, Box<dyn hullgeom::source::HullSource>) = match fleet(v)? {
        Fleet::Nurbs(f) => (Kind::Iges, Box::new(f)),
        Fleet::Mesh(f) => (Kind::Stl, Box::new(f)),
    };
    Ok(SourceFile {
        path: "geometry".into(),
        kind,
        source,
        waterline_z: 0.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Wigley hull as IGES, its design waterline at CAD height `wl`.
    fn wigley(wl: f64) -> Vec<u8> {
        let mut surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        for p in surfaces.iter_mut().flat_map(|s| s.ctrl.iter_mut()) {
            p[2] += wl;
        }
        hullgeom::iges::write(&surfaces, "wigley")
            .unwrap()
            .into_bytes()
    }

    fn hull0(v: &Value) -> (f64, f64, f64) {
        let h = &v["hulls"][0];
        let f = |k: &str| h[k].as_f64().unwrap();
        (f("displaced_volume"), f("length"), f("draft"))
    }

    /// The JSON geometry cuts into the same hull as the file it came from,
    /// with the file's waterline and a scale applied once.
    #[test]
    fn a_file_and_its_geometry_cut_alike() {
        let req = LoftRequest {
            waterline: Some(1.5),
            scale: Some(1.2),
            ..LoftRequest::default()
        };
        let file = crate::loft("w.igs", wigley(1.5), &req).unwrap();
        let g = from_file("w.igs", wigley(1.5), &req).unwrap();
        assert_eq!(g["kind"], "nurbs");
        let bytes = serde_json::to_vec(&g).unwrap();
        assert!(is_native(&bytes));
        let native = crate::loft("hull.json", bytes, &LoftRequest::default()).unwrap();
        let (a, b) = (hull0(&file), hull0(&native));
        let exact = 4.0 / 9.0 * 10.0 * 0.625 * 1.2f64.powi(3);
        assert!((a.0 - exact).abs() < 1e-6 * exact, "{a:?}");
        for (x, y) in [(a.0, b.0), (a.1, b.1), (a.2, b.2)] {
            assert!((x - y).abs() < 1e-9 * x.abs().max(1.0), "{a:?} vs {b:?}");
        }
    }

    /// The Wigley as an ASCII STL in millimetres, its waterline at `wl` [m].
    fn wigley_stl(wl: f64) -> Vec<u8> {
        let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let fleet = SourceFleet::from_hulls(vec![surfaces.to_vec()]).unwrap();
        let (v, t) = fleet
            .posed_tessellation(0, 0.0, &HullPose::default(), &Platform::default())
            .unwrap();
        let mut s = String::from("solid w\n");
        for tri in t {
            s.push_str("facet normal 0 0 0\nouter loop\n");
            for k in tri {
                let p = v[k as usize];
                s.push_str(&format!(
                    "vertex {} {} {}\n",
                    p[0] * 1e3,
                    p[1] * 1e3,
                    (p[2] + wl) * 1e3
                ));
            }
            s.push_str("endloop\nendfacet\n");
        }
        s.push_str("endsolid w\n");
        s.into_bytes()
    }

    #[test]
    fn an_stl_and_its_geometry_cut_alike() {
        let req = LoftRequest {
            waterline: Some(0.4),
            units: Some(1e-3),
            ..LoftRequest::default()
        };
        let file = crate::loft("w.stl", wigley_stl(0.4), &req).unwrap();
        let g = from_file("w.stl", wigley_stl(0.4), &req).unwrap();
        assert_eq!(g["kind"], "mesh");
        let native = crate::loft(
            "hull.json",
            serde_json::to_vec(&g).unwrap(),
            &LoftRequest::default(),
        )
        .unwrap();
        let (a, b) = (hull0(&file), hull0(&native));
        assert!(a.0 > 0.0, "{a:?}");
        for (x, y) in [(a.0, b.0), (a.1, b.1), (a.2, b.2)] {
            assert!((x - y).abs() < 1e-9 * x.abs().max(1.0), "{a:?} vs {b:?}");
        }
    }

    #[test]
    fn bad_geometry_is_refused() {
        let mut g = from_file("w.igs", wigley(0.0), &LoftRequest::default()).unwrap();
        g["hulls"][0]["patches"][0]["knots_u"] = json!([0.0, 1.0]);
        let e = open_value(&g).err().unwrap();
        assert!(e.contains("knots"), "{e}");
        assert!(open_value(&json!({ "kind": "voxels", "hulls": [] })).is_err());
    }
}
