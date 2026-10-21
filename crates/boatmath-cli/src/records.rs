//! Making records. The definitions — hulls from files, scaled hulls, cases
//! on hulls, studies on cases — and the statics computed from them. Each
//! record's id is that of its inputs, so making the same thing twice gives
//! the same record.

use crate::cache::Cache;
use crate::stream::{id_of, sha256, short, Stream};
use boatmath::params::{CaseParams, StudyParams};
use boatmath::{hull_summary, loft, statics, wave_gz, LoftRequest, SOLVER_VERSION};
use serde_json::{json, Value};
use std::path::Path;

/// The type a record must be, or an error naming what it is instead.
pub fn expect<'a>(r: &'a Value, kind: &str) -> Result<&'a str, String> {
    let id = r["id"].as_str().unwrap_or("?");
    match r["type"].as_str() {
        Some(t) if t == kind => Ok(id),
        Some(t) => Err(format!("expected a {kind}, got a {t} ({})", short(id))),
        None => Err(format!("expected a {kind} record, got {}", truncate(r))),
    }
}

fn truncate(v: &Value) -> String {
    let s = v.to_string();
    if s.len() > 60 {
        format!("{}…", &s[..60])
    } else {
        s
    }
}

/// How a hull is cut into sections: settings of the cut, not of the
/// geometry, so they stay with the hull record.
#[derive(Clone, Copy, Debug, Default, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CutSettings {
    /// Stations along the hull.
    pub stations: Option<usize>,
    /// Rays across each section.
    pub rays: Option<usize>,
    /// Centreplane override [m, y].
    pub centerplane: Option<f64>,
}

/// What the solver needs of a hull: its geometry, as the bytes `boatmath`
/// reads, and how it is cut.
pub struct HullSource {
    pub file_name: String,
    pub bytes: Vec<u8>,
    pub import: LoftRequest,
}

pub fn hull_source(hull: &Value) -> Result<HullSource, String> {
    let cut: CutSettings =
        serde_json::from_value(hull["cut"].clone()).map_err(|e| format!("hull cut: {e}"))?;
    if !hull["geometry"].is_object() {
        return Err("a hull without geometry".into());
    }
    Ok(HullSource {
        file_name: "geometry.json".into(),
        bytes: serde_json::to_vec(&hull["geometry"]).map_err(|e| e.to_string())?,
        import: LoftRequest {
            stations: cut.stations,
            rays: cut.rays,
            centerplane: cut.centerplane,
            ..LoftRequest::default()
        },
    })
}

/// A hull record, checked to cut. Its hydrostatics are `statics`'.
fn hull_record(
    name: &str,
    source: Value,
    geometry: Value,
    cut: CutSettings,
    parent: Option<&str>,
) -> Result<Value, String> {
    let r = json!({
        "type": "hull",
        "id": id_of(&json!({ "geometry": geometry, "cut": cut })),
        "name": name,
        "source": source,
        "cut": cut,
        "parent": parent,
        "geometry": geometry,
    });
    let src = hull_source(&r)?;
    loft(&src.file_name, src.bytes, &src.import)?;
    Ok(r)
}

/// A hull from a file: its geometry, with the file's import settings
/// (waterline, units) applied once. `source` records where it came from.
pub fn hull_from_file(
    path: &Path,
    name: Option<&str>,
    import: LoftRequest,
) -> Result<Value, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "hull".into());
    let stem = path
        .file_stem()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| file_name.clone());
    let abs = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
    let source = json!({
        "path": abs.to_string_lossy(),
        "sha256": sha256(&bytes),
        "waterline": import.waterline,
        "units": import.units,
    });
    let geometry = boatmath::native::from_file(&file_name, bytes, &import)
        .map_err(|e| format!("{file_name}: {e}"))?;
    let cut = CutSettings {
        stations: import.stations,
        rays: import.rays,
        centerplane: import.centerplane,
    };
    hull_record(name.unwrap_or(&stem), source, geometry, cut, None)
}

/// How to scale a hull, about its design waterline.
#[derive(Clone, Copy, Debug)]
pub enum Scaling {
    /// The whole hull by this factor.
    By(f64),
    /// Beam and draft by this factor, the length kept.
    Beam(f64),
    /// To carry this mass [kg] at the design waterline: uniformly, or with
    /// `keep_length` beam and draft only.
    Mass { kg: f64, keep_length: bool },
}

/// A hull's design displacement [kg], in seawater.
fn design_mass(hull: &Value) -> Result<f64, String> {
    let src = hull_source(hull)?;
    let cut = loft(&src.file_name, src.bytes, &src.import)?;
    let v: f64 = cut["hulls"]
        .as_array()
        .ok_or("no hulls")?
        .iter()
        .filter_map(|h| h["displaced_volume"].as_f64())
        .sum();
    Ok(hullgeom::Fluid::SEAWATER_15C.density * v)
}

/// A hull made by scaling another, about its design waterline. The new
/// geometry is the old with its control points (or vertices) scaled.
pub fn scaled_hull(hull: &Value, how: Scaling, name: Option<&str>) -> Result<Value, String> {
    let id = expect(hull, "hull")?;
    let base = hull["name"].as_str().unwrap_or("hull");
    let positive = |k: f64, what: &str| {
        if k > 0.0 && k.is_finite() {
            Ok(k)
        } else {
            Err(format!("{what} {k}: expected a positive number"))
        }
    };
    let (by, yz, label) = match how {
        Scaling::By(k) => (Some(positive(k, "--by")?), None, format!("{base} ×{k}")),
        Scaling::Beam(k) => (
            None,
            Some(positive(k, "--beam")?),
            format!("{base} ×{k} beam"),
        ),
        Scaling::Mass { kg, keep_length } => {
            let m0 = design_mass(hull)?;
            let r = positive(kg, "--mass")? / m0;
            if keep_length {
                (
                    None,
                    Some(r.sqrt()),
                    format!("{base} at {kg} kg, length kept"),
                )
            } else {
                (Some(r.cbrt()), None, format!("{base} at {kg} kg"))
            }
        }
    };
    let pose = LoftRequest {
        scale: by,
        scale_yz: yz,
        ..LoftRequest::default()
    }
    .pose();
    let geometry = boatmath::native::posed(&hull["geometry"], &pose)?;
    let cut: CutSettings =
        serde_json::from_value(hull["cut"].clone()).map_err(|e| format!("hull cut: {e}"))?;
    let source = json!({ "scaled": { "by": by, "beam": yz } });
    hull_record(name.unwrap_or(&label), source, geometry, cut, Some(id))
}

/// A case's params as the CLI writes them: the web app's `CaseParams`
/// without `mass_by` (a CLI case always carries its mass by sinking; to
/// carry it at the design waterline, scale the hull first).
fn case_params_json(p: &CaseParams) -> Value {
    let mut v = serde_json::to_value(p).expect("params serialize");
    if let Some(o) = v.as_object_mut() {
        o.remove("mass_by");
    }
    v
}

/// A case's params, for the solver.
pub fn case_params(case: &Value) -> Result<CaseParams, String> {
    let v = case["params"].clone();
    serde_json::from_value(v).map_err(|e| format!("case params: {e}"))
}

/// A case on a hull: a definition only.
pub fn case(hull: &Value, params: CaseParams, name: Option<&str>) -> Result<Value, String> {
    let hull_id = expect(hull, "hull")?;
    let params = params.canonical()?;
    let p = case_params_json(&params);
    Ok(json!({
        "type": "case",
        "id": id_of(&json!({ "hull": hull_id, "params": p })),
        "name": name.unwrap_or(""),
        "hull": hull_id,
        "params": p,
        "parent": null,
    }))
}

/// A study on a case: only the request; `run` computes it.
pub fn study(case: &Value, params: StudyParams) -> Result<Value, String> {
    let case_id = expect(case, "case")?;
    Ok(study_record(case_id, params.canonical()?))
}

pub fn study_record(case_id: &str, params: StudyParams) -> Value {
    json!({
        "type": "study",
        "id": id_of(&json!({ "case": case_id, "params": params })),
        "case": case_id,
        "kind": params.kind(),
        "params": params,
    })
}

/// What the solver needs of a case: its hull, its parameters, and its
/// geometry — the hull's, with the case's mount (if it has one) appended
/// as members of their own — and where each drive's propeller sits.
pub struct CaseSource<'a> {
    pub hull: &'a Value,
    pub params: CaseParams,
    pub src: HullSource,
    pub mount: Option<boatmath::mount::Mount>,
    /// Each drive's propeller, in the hull's design frame (one per hull).
    pub thrusts: Vec<boatmath::mount::Thrust>,
}

/// A case's mount, if it has one.
pub fn case_mount(case: &Value) -> Result<Option<boatmath::mount::Mount>, String> {
    match case["params"].get("mount") {
        None | Some(Value::Null) => Ok(None),
        Some(m) => serde_json::from_value(m.clone())
            .map(Some)
            .map_err(|e| format!("case mount: {e}")),
    }
}

pub fn case_source<'a>(s: &'a Stream, case: &Value) -> Result<CaseSource<'a>, String> {
    let hull = s.follow(case, "hull")?;
    let params = case_params(case)?;
    let mut src = hull_source(hull)?;
    let mount = params
        .mount
        .clone()
        .or_else(|| case_mount(case).ok().flatten());
    let mut thrusts = Vec::new();
    if let Some(m) = &mount {
        let bare = loft(&src.file_name, src.bytes.clone(), &src.import)?;
        let (g, t) = boatmath::mount::mounted(&hull["geometry"], &bare, m)?;
        src.bytes = serde_json::to_vec(&g).map_err(|e| e.to_string())?;
        // A catamaran doubles its hull, and the drives go with each copy,
        // each as far off its centreplane as on the single hull.
        let yc = bare["hulls"][0]["placement"]["y"].as_f64().unwrap_or(0.0);
        thrusts = match params.span {
            Some(span) => t
                .iter()
                .flat_map(|d| {
                    let off = d.y - yc;
                    [0.5 * span, -0.5 * span].map(|y| boatmath::mount::Thrust { y: y + off, ..*d })
                })
                .collect(),
            None => t,
        };
    }
    Ok(CaseSource {
        hull,
        params,
        src,
        mount,
        thrusts,
    })
}

/// A case with a drive: the case's own parameters and `mount`, its parent
/// the bare case (a case already mounted is remounted from its parent).
pub fn mounted_case(
    case: &Value,
    mount: &boatmath::mount::Mount,
    name: Option<&str>,
) -> Result<Value, String> {
    let case_id = expect(case, "case")?;
    let hull_id = case["hull"].as_str().ok_or("case: no hull")?;
    let mut p = case["params"].clone();
    let o = p.as_object_mut().ok_or("case: params not an object")?;
    let remount = o.contains_key("mount");
    o.insert(
        "mount".into(),
        serde_json::to_value(mount).map_err(|e| e.to_string())?,
    );
    let parent = match (remount, case["parent"].as_str()) {
        (true, Some(pid)) => pid,
        _ => case_id,
    };
    let name = name
        .map(str::to_string)
        .unwrap_or_else(|| case["name"].as_str().unwrap_or("").to_string());
    Ok(json!({
        "type": "case",
        "id": id_of(&json!({ "hull": hull_id, "params": p })),
        "name": name,
        "hull": hull_id,
        "params": p,
        "parent": parent,
    }))
}

/// The sections record of a case cut at `(sinkage [m], trim [rad])`. Its id
/// is that of the case and the attitude, so a case at rest and the studies
/// held there share one.
pub fn sections(
    s: &Stream,
    case: &Value,
    attitude: (f64, f64),
    cache: &Cache,
) -> Result<Value, String> {
    let case_id = expect(case, "case")?;
    let id = id_of(&json!({ "case": case_id, "attitude": [attitude.0, attitude.1] }));
    if let Some(r) = s
        .get("sections", &id)
        .cloned()
        .or_else(|| cache.get("sections", &id))
    {
        return Ok(r);
    }
    let CaseSource { params, src, .. } = case_source(s, case)?;
    let mut r =
        boatmath::sections::sections(&src.file_name, src.bytes, &src.import, &params, attitude)?;
    let o = r.as_object_mut().expect("an object");
    o.insert("type".into(), json!("sections"));
    o.insert("id".into(), json!(id));
    o.insert("case".into(), json!(case_id));
    o.insert("solver_version".into(), json!(SOLVER_VERSION));
    cache.put(&r);
    Ok(r)
}

/// A hull's statics: its hydrostatics at the design waterline.
pub fn hull_statics(hull: &Value, cache: &Cache) -> Result<Value, String> {
    let id = expect(hull, "hull")?;
    if let Some(r) = cache.get("statics", id) {
        return Ok(r);
    }
    let src = hull_source(hull)?;
    let summary = hull_summary(&loft(&src.file_name, src.bytes, &src.import)?);
    let r = json!({
        "type": "statics",
        "id": id,
        "hull": id,
        "hulls": summary["hulls"],
        "notes": summary["notes"],
        "solver_version": SOLVER_VERSION,
    });
    cache.put(&r);
    Ok(r)
}

/// A case's statics, and the sections record of its float at rest (first):
/// the equilibrium, its hydrostatics, roll stability and the GZ curve.
pub fn case_statics(s: &Stream, case: &Value, cache: &Cache) -> Result<Vec<Value>, String> {
    let id = expect(case, "case")?;
    if let Some(r) = cache.get("statics", id) {
        let mut out = Vec::new();
        if let Some(sid) = r["sections"].as_str() {
            match cache.get("sections", sid) {
                Some(sec) => out.push(sec),
                None => return case_statics(s, case, &Cache::default()),
            }
        }
        out.push(r);
        return Ok(out);
    }
    let CaseSource { params, src, .. } = case_source(s, case)?;
    let mut r = json!({ "type": "statics", "id": id, "case": id });
    let mut out = Vec::new();
    match statics(&src.file_name, src.bytes, &src.import, &params) {
        Ok(v) => {
            let Value::Object(o) = v else {
                return Err("statics: not an object".into());
            };
            for (k, v) in o {
                // The meshes are for drawing; a picture re-poses the geometry.
                if k != "meshes" && k != "body_meshes" {
                    r[k] = v;
                }
            }
            let rest = (
                r["at_rest"]["sinkage"].as_f64(),
                r["at_rest"]["trim_rad"].as_f64(),
            );
            if let (Some(z), Some(t)) = rest {
                match sections(s, case, (z, t), cache) {
                    Ok(sec) => {
                        r["sections"] = sec["id"].clone();
                        out.push(sec);
                    }
                    Err(e) => eprintln!("boatmath: case {}: sections: {e}", short(id)),
                }
            }
        }
        Err(e) => r["error"] = json!(e),
    }
    r["solver_version"] = json!(SOLVER_VERSION);
    cache.put(&r);
    out.push(r);
    Ok(out)
}

/// A case's quasi-static GZ in a regular wave of given length and height [m].
pub fn gz_waves(
    s: &Stream,
    case: &Value,
    length: f64,
    height: f64,
    cache: &Cache,
) -> Result<Value, String> {
    let case_id = expect(case, "case")?;
    let id = id_of(&json!({ "case": case_id, "length": length, "height": height }));
    if let Some(r) = cache.get("gz_waves", &id) {
        return Ok(r);
    }
    let CaseSource { params, src, .. } = case_source(s, case)?;
    let mut r = json!({
        "type": "gz_waves",
        "id": id,
        "case": case_id,
        "length": length,
        "height": height,
    });
    match wave_gz(
        &src.file_name,
        src.bytes,
        &src.import,
        &params,
        length,
        height,
    ) {
        Ok(v) => {
            let Value::Object(o) = v else {
                return Err("wave_gz: not an object".into());
            };
            for (k, v) in o {
                r[k] = v;
            }
        }
        Err(e) => r["error"] = json!(e),
    }
    r["solver_version"] = json!(SOLVER_VERSION);
    cache.put(&r);
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn case_params_with_and_without_mount() {
        let p_bare = CaseParams::default().canonical().unwrap();
        let json_bare = case_params_json(&p_bare);
        assert!(json_bare.get("mount").is_none());

        let case_val = json!({
            "type": "case",
            "id": "abc",
            "hull": "h123",
            "params": json_bare,
        });
        let roundtrip = case_params(&case_val).unwrap();
        assert_eq!(roundtrip.mount, None);

        let mount = boatmath::mount::build(
            Some("oceanvolt-sd8"),
            &boatmath::mount::Options {
                x: Some(2.0),
                ..Default::default()
            },
        )
        .unwrap();

        let mut p_mounted = p_bare.clone();
        p_mounted.mount = Some(mount.clone());
        let json_mounted = case_params_json(&p_mounted);
        assert!(json_mounted.get("mount").is_some());

        let case_mounted_val = json!({
            "type": "case",
            "id": "def",
            "hull": "h123",
            "params": json_mounted,
        });
        let roundtrip_mounted = case_params(&case_mounted_val).unwrap();
        assert_eq!(roundtrip_mounted.mount, Some(mount.clone()));

        let mounted_record = mounted_case(&case_val, &mount, None).unwrap();
        assert_eq!(mounted_record["params"], json_mounted);
    }
}
