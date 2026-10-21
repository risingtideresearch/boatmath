//! The JSON API: hulls, their cases, studies on those,
//! results, and the queue. Bodies and answers are JSON unless noted.
//!
//! ```text
//! GET  /api/version                     {"solver_version": ..}
//! GET  /api/whoami                      {"login", "name"} from the tailnet, or null
//!
//! Behind `tailscale serve` (the server bound to loopback), uploads and
//! studies are labelled with the tailnet's name for the asker, and a `by`
//! given in the request is ignored.
//!
//! GET  /api/hulls                       every hull, newest first
//! POST /api/hulls?file=F&name=&notes=&by=&parent=&waterline=&centerplane=&stations=&rays=&units=
//!                                       body = the file's bytes; cut, saved →
//!                                       {"id", "created", "hull"}
//! GET  /api/hulls/:id                   one hull
//! POST /api/hulls/:id                   {"name"?, "notes"?}
//! GET  /api/hulls/:id/sections          the cut (as /api/loft answers)
//! GET  /api/hulls/:id/geometry          the display geometry (as /api/geometry)
//! GET  /api/hulls/:id/file              the uploaded file
//! POST /api/hulls/:id/scale             {"mode": "xyz" | "yz", "length"? | "displacement"?,
//!                                        "name"?, "notes"?, "by"?, "dry_run"?}
//!                                       a new hull: this one scaled about its design
//!                                       waterline to a waterline length [m] (xyz) or
//!                                       a displacement there [kg] → {"id", "created",
//!                                       "hull", "factor"}; a dry run answers the factor
//!                                       and the dimensions it gives
//!
//! GET  /api/cases?hull=               cases, with their statics' summary
//! POST /api/cases                     {"hull", "params": CaseParams, "name"?, "by"?,
//!                                        "dry_run"?}
//!                                       → {"id", "created", "case"}; the statics
//!                                       are computed before it is saved. A dry study
//!                                       checks it and answers {"params", "existing"}
//! GET  /api/cases/:id                 one (its statics computed again if stale)
//! POST /api/cases/:id                 {"name"?, "notes"?}
//! GET  /api/cases/:id/gz_waves?length=&height=
//!                                       its GZ curves in a regular wave held still
//!                                       around it: beam and quartering, the worst
//!                                       crest of each, and the energy to peak GZ
//! GET  /api/cases/:id/statics         the whole statics (GZ curve, meshes at rest)
//!
//! POST /api/studies                        {"case_ids": [..], "studies": [StudyParams..],
//!                                        "priority"?, "by"?, "dry_run"?}
//!                                       every study on every case (a study in
//!                                       waves brings its calm-water study) →
//!                                       {"studies": [..]}; a dry run
//!                                       queues nothing and says which exist
//! GET  /api/studies?hull=&case=&ids=a,b&status=a,b&kind=&limit=&sort=finished
//!                                       studies with their case and result
//! GET  /api/studies/:id                    one study
//! GET  /api/studies/:id/result             its whole answer (fields and all)
//! POST /api/studies/:id/cancel             stop it (queued, or running)
//! POST /api/studies/:id/requeue            run it again
//! POST /api/studies/:id/priority?p=N       higher runs sooner
//! POST /api/requeue_stale?hull=         requeue every study with a stale result
//! GET  /api/queue                       {"running", "studies", "counts"}
//! ```

use boatmath_web::params::{CaseParams, StudyParams};
use boatmath_web::store::{NewHull, Store, StudyFilter, SOLVER_VERSION};
use boatmath_web::worker::Worker;
use boatmath_web::{geometry, loft, statics, wave_gz, LoftRequest};
use serde_json::{json, Value};
use std::sync::Arc;
use tiny_http::Method;

pub struct App {
    pub store: Arc<Store>,
    pub worker: Arc<Worker>,
}

pub enum Reply {
    Json(Value),
    /// A stored blob, gzipped, and its content type.
    Gz(Vec<u8>, &'static str),
}

/// An answer that is not a success: a status and a message.
pub type Fail = (u16, String);

fn bad(e: impl std::fmt::Display) -> Fail {
    (400, e.to_string())
}

fn missing(what: &str) -> Fail {
    (404, format!("no such {what}"))
}

fn internal(e: String) -> Fail {
    (500, e)
}

/// Who is asking, as `tailscale serve` says (its `Tailscale-User-*`
/// headers): what uploads and studies are labelled with, in place of the name
/// the page sends.
pub struct Who {
    pub login: String,
    pub name: String,
}

/// The label for work asked for: the tailnet's name for the asker when
/// there is one, else what the request gave.
fn by<'a>(who: Option<&'a Who>, given: &'a str) -> &'a str {
    match who {
        Some(w) if !w.name.is_empty() => &w.name,
        Some(w) => &w.login,
        None => given.trim(),
    }
}

/// Route an `/api/…` request; `None` if no route matches.
pub fn route(
    app: &App,
    method: &Method,
    path: &str,
    q: &[(String, String)],
    who: Option<&Who>,
    body: &mut dyn FnMut() -> Result<Vec<u8>, String>,
) -> Option<Result<Reply, Fail>> {
    let seg: Vec<&str> = path.trim_matches('/').split('/').collect();
    if seg.first() != Some(&"api") {
        return None;
    }
    let get = |k: &str| {
        q.iter()
            .find(|(n, v)| n == k && !v.trim().is_empty())
            .map(|(_, v)| v.trim())
    };
    let id = |s: &str| s.parse::<i64>().map_err(|_| bad(format!("bad id {s:?}")));
    let opt_id = |k: &str| get(k).map(id).transpose();
    let json_body = |body: &mut dyn FnMut() -> Result<Vec<u8>, String>| -> Result<Value, Fail> {
        serde_json::from_slice(&body().map_err(bad)?).map_err(|e| bad(format!("body: {e}")))
    };
    let s = &app.store;
    let list = |v: Result<Vec<Value>, String>| v.map(|v| Reply::Json(json!(v))).map_err(internal);
    let one = |v: Result<Option<Value>, String>, what: &str| {
        v.map_err(internal)?.map(Reply::Json).ok_or(missing(what))
    };
    use Method::{Get, Post};
    let r = match (method, &seg[1..]) {
        (Get, ["version"]) => Ok(Reply::Json(json!({ "solver_version": SOLVER_VERSION }))),
        (Get, ["whoami"]) => Ok(Reply::Json(match who {
            Some(w) => json!({ "login": w.login, "name": w.name }),
            None => Value::Null,
        })),

        (Get, ["hulls"]) => list(s.hulls()),
        (Post, ["hulls"]) => upload(app, q, who, body),
        (Get, ["hulls", h]) => id(h).and_then(|h| one(s.hull(h), "hull")),
        (Post, ["hulls", h, "scale"]) => {
            id(h).and_then(|h| json_body(body).and_then(|v| scale_hull(app, h, &v, who)))
        }
        (Post, ["hulls", h]) => id(h).and_then(|h| {
            let v = json_body(body)?;
            let text = |k| v[k].as_str().map(str::trim);
            if s.update_hull(h, text("name"), text("notes"))
                .map_err(internal)?
            {
                one(s.hull(h), "hull")
            } else {
                Err(missing("hull"))
            }
        }),
        (Get, ["hulls", h, which @ ("file" | "sections" | "geometry")]) => id(h).and_then(|h| {
            let blob = s
                .hull_blob(h, which)
                .map_err(internal)?
                .ok_or(missing("hull"))?;
            let kind = if *which == "file" {
                "application/octet-stream"
            } else {
                "application/json"
            };
            Ok(Reply::Gz(s.blob_gz(&blob).map_err(internal)?, kind))
        }),

        (Get, ["cases"]) => opt_id("hull").and_then(|h| list(s.cases(h))),
        (Post, ["cases"]) => json_body(body).and_then(|v| {
            let hull = v["hull"].as_i64().ok_or(bad("hull (an id) is required"))?;
            let p = case_params(&v["params"])?;
            if v["dry_run"].as_bool() == Some(true) {
                // Checked, and found if it exists; nothing made.
                check_case(app, hull, &p)?;
                return Ok(Reply::Json(json!({
                    "params": p,
                    "existing": s.find_case(hull, &p).map_err(internal)?,
                })));
            }
            let by = by(who, v["by"].as_str().unwrap_or(""));
            let (id, created, _) = make_case(app, hull, &p, v["name"].as_str(), by)?;
            Ok(Reply::Json(json!({
                "id": id,
                "created": created,
                "case": s.case(id).map_err(internal)?,
            })))
        }),
        (Get, ["cases", c]) => id(c).and_then(|c| {
            let cfg = s.case(c).map_err(internal)?.ok_or(missing("case"))?;
            if cfg["stale"] == true {
                // Statics are cheap: computed again when looked at.
                recompute(app, c)?;
                return one(s.case(c), "case");
            }
            Ok(Reply::Json(cfg))
        }),
        (Post, ["cases", c]) => id(c).and_then(|c| {
            let v = json_body(body)?;
            let text = |k| v[k].as_str().map(str::trim);
            if s.update_case(c, text("name"), text("notes"))
                .map_err(internal)?
            {
                one(s.case(c), "case")
            } else {
                Err(missing("case"))
            }
        }),
        (Get, ["cases", c, "gz_waves"]) => id(c).and_then(|c| {
            let num = |k: &str| {
                get(k)
                    .ok_or(bad(format!("{k} is required")))?
                    .parse::<f64>()
                    .map_err(|_| bad(format!("{k}: expected a number")))
            };
            let (length, height) = (num("length")?, num("height")?);
            let (hull, p, _) = s.case_source(c).map_err(internal)?.ok_or(missing("case"))?;
            let src = s
                .hull_source(hull)
                .map_err(internal)?
                .ok_or(missing("hull"))?;
            let bytes = s.blob(&src.file_blob).map_err(internal)?;
            wave_gz(&src.file_name, bytes, &src.import, &p, length, height)
                .map(Reply::Json)
                .map_err(bad)
        }),
        (Get, ["cases", c, "statics"]) => id(c).and_then(|c| {
            let (_, _, blob) = s.case_source(c).map_err(internal)?.ok_or(missing("case"))?;
            let blob = blob.ok_or((409, "its statics could not be computed".to_string()))?;
            Ok(Reply::Gz(
                s.blob_gz(&blob).map_err(internal)?,
                "application/json",
            ))
        }),

        (Post, ["studies"]) => json_body(body).and_then(|v| add_studies(app, &v, who)),
        (Get, ["studies"]) => (|| {
            let f = StudyFilter {
                hull: opt_id("hull")?,
                case: opt_id("case")?,
                ids: get("ids").map(String::from),
                status: get("status").map(String::from),
                kind: get("kind").map(String::from),
                limit: opt_id("limit")?,
                by_finish: get("sort") == Some("finished"),
            };
            list(s.studies(&f))
        })(),
        (Get, ["studies", r]) => id(r).and_then(|r| one(s.study(r), "study")),
        (Get, ["studies", r, "result"]) => id(r).and_then(|r| {
            let blob = s
                .result_blob(r)
                .map_err(internal)?
                .ok_or(missing("result"))?;
            Ok(Reply::Gz(
                s.blob_gz(&blob).map_err(internal)?,
                "application/json",
            ))
        }),
        (Post, ["studies", r, "cancel"]) => id(r).and_then(|r| {
            let before = s
                .cancel_queued(r)
                .map_err(internal)?
                .ok_or(missing("study"))?;
            let stopping = before == "running" && app.worker.cancel(r);
            Ok(Reply::Json(json!({ "was": before, "stopping": stopping })))
        }),
        (Post, ["studies", r, "requeue"]) => id(r).and_then(|r| {
            let n = s.requeue(Some(&[r]), None).map_err(internal)?;
            app.worker.poke();
            Ok(Reply::Json(json!({ "requeued": n })))
        }),
        (Post, ["studies", r, "priority"]) => id(r).and_then(|r| {
            let p = get("p").ok_or(bad("p is required"))?;
            let p = p
                .parse::<i64>()
                .map_err(|_| bad(format!("p: bad number {p:?}")))?;
            if s.set_priority(r, p).map_err(internal)? {
                Ok(Reply::Json(json!({ "id": r, "priority": p })))
            } else {
                Err(missing("study"))
            }
        }),
        (Post, ["requeue_stale"]) => opt_id("hull").and_then(|h| {
            let n = s.requeue(None, h).map_err(internal)?;
            app.worker.poke();
            Ok(Reply::Json(json!({ "requeued": n })))
        }),
        (Get, ["queue"]) => (|| {
            Ok(Reply::Json(json!({
                "running": app.worker.status(),
                "studies": s.queue().map_err(internal)?,
                "counts": s.counts().map_err(internal)?,
                "solver_version": SOLVER_VERSION,
            })))
        })(),
        _ => return None,
    };
    Some(r)
}

/// `POST /api/hulls/:id/scale`: a new hull, the same file cut with a scale
/// added to its settings — about its design waterline, so the waterline
/// stays: the whole hull by `k` to a waterline length (`L/L₀`) or a
/// displacement (`(Δ/Δ₀)^⅓`), or its beam and draft only to a displacement
/// (`(Δ/Δ₀)^½`; the length is kept, so a length cannot be the target). A
/// hull already scaled is scaled again: the factors multiply. A file of
/// several hulls is refused (each would scale about its own centre, not
/// the gaps between them).
fn scale_hull(app: &App, id: i64, v: &Value, who: Option<&Who>) -> Result<Reply, Fail> {
    let s = &app.store;
    let hull = s.hull(id).map_err(internal)?.ok_or(missing("hull"))?;
    let src = s
        .hull_source(id)
        .map_err(internal)?
        .ok_or(missing("hull"))?;
    let hs = hull["summary"]["hulls"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    if hs.len() != 1 {
        return Err(bad(format!(
            "scaling makes a new hull of a file with one hull; this one holds {}",
            hs.len()
        )));
    }
    let h0 = &hs[0];
    let dim = |k: &str| h0[k].as_f64().filter(|x| *x > 0.0);
    let (l0, b0, t0, v0) = match (
        dim("length"),
        dim("beam"),
        dim("draft"),
        dim("displaced_volume"),
    ) {
        (Some(l), Some(b), Some(t), Some(vol)) => (l, b, t, vol),
        _ => return Err(bad("this hull's dimensions are unknown")),
    };
    let rho = hullgeom::Fluid::SEAWATER_15C.density;
    let d0 = rho * v0;
    let target = |k: &str| v[k].as_f64().filter(|x| x.is_finite());
    let yz = match v["mode"].as_str().unwrap_or("xyz") {
        "xyz" => false,
        "yz" => true,
        m => return Err(bad(format!("mode {m:?}: expected xyz or yz"))),
    };
    let (factor, what) = match (yz, target("length"), target("displacement")) {
        (_, Some(_), Some(_)) => return Err(bad("give a length or a displacement, not both")),
        (false, Some(l), None) => (l / l0, format!("L {l} m")),
        (true, Some(_), None) => {
            return Err(bad(
                "scaling beam and draft keeps the length; give a displacement",
            ))
        }
        (false, None, Some(d)) => ((d / d0).cbrt(), format!("{d} kg")),
        (true, None, Some(d)) => ((d / d0).sqrt(), format!("{d} kg, beam and draft")),
        (_, None, None) => return Err(bad("give a length [m] or a displacement [kg]")),
    };
    if !(factor > 0.05 && factor < 20.0) {
        return Err(bad(format!(
            "that is a factor of {factor:.4}; expected between 0.05 and 20"
        )));
    }
    let (kl, kb, kv) = if yz {
        (1.0, factor, factor * factor)
    } else {
        (factor, factor, factor.powi(3))
    };
    let dims = json!({
        "factor": factor,
        "mode": if yz { "yz" } else { "xyz" },
        "length": l0 * kl, "beam": b0 * kb, "draft": t0 * kb,
        "displacement": d0 * kv,
        "from": { "length": l0, "beam": b0, "draft": t0, "displacement": d0 },
    });
    if v["dry_run"].as_bool() == Some(true) {
        return Ok(Reply::Json(dims));
    }
    let mut import = src.import;
    if yz {
        import.scale_yz = Some(import.scale_yz.unwrap_or(1.0) * factor);
    } else {
        import.scale = Some(import.scale.unwrap_or(1.0) * factor);
    }
    let bytes = s.blob(&src.file_blob).map_err(internal)?;
    let sections = loft(&src.file_name, bytes.clone(), &import).map_err(bad)?;
    let geom = geometry(&src.file_name, bytes.clone(), &import).map_err(bad)?;
    let parent_name = hull["name"].as_str().unwrap_or("hull");
    let default_name = format!("{parent_name} scaled to {what}");
    let name = v["name"]
        .as_str()
        .map(str::trim)
        .filter(|n| !n.is_empty())
        .unwrap_or(&default_name);
    let (nid, created) = s
        .add_hull(&NewHull {
            name,
            notes: v["notes"].as_str().unwrap_or("").trim(),
            uploaded_by: by(who, v["by"].as_str().unwrap_or("")),
            file_name: &src.file_name,
            bytes: &bytes,
            import,
            parent_id: Some(id),
            sections: &sections,
            geometry: &geom,
        })
        .map_err(internal)?;
    Ok(Reply::Json(json!({
        "id": nid,
        "created": created,
        "factor": factor,
        "expected": dims,
        "hull": s.hull(nid).map_err(internal)?,
    })))
}

/// `POST /api/hulls`: cut the upload (refusing what cannot be cut) and save
/// it with its sections and display geometry.
fn upload(
    app: &App,
    q: &[(String, String)],
    who: Option<&Who>,
    body: &mut dyn FnMut() -> Result<Vec<u8>, String>,
) -> Result<Reply, Fail> {
    let get = |k: &str| q.iter().find(|(n, _)| n == k).map_or("", |(_, v)| v.trim());
    let file = get("file");
    if file.is_empty() {
        return Err(bad("file (the file's name, for its format) is required"));
    }
    let import = LoftRequest::from_query(q).map_err(bad)?;
    let parent_id = match get("parent") {
        "" => None,
        p => Some(
            p.parse::<i64>()
                .map_err(|_| bad(format!("parent: bad id {p:?}")))?,
        ),
    };
    let bytes = body().map_err(bad)?;
    let t0 = std::time::Instant::now();
    let sections = loft(file, bytes.clone(), &import).map_err(bad)?;
    let geom = geometry(file, bytes.clone(), &import).map_err(bad)?;
    let stem = file.rsplit_once('.').map_or(file, |(s, _)| s);
    let name = match get("name") {
        "" => stem,
        n => n,
    };
    let (id, created) = app
        .store
        .add_hull(&NewHull {
            name,
            notes: get("notes"),
            uploaded_by: by(who, get("by")),
            file_name: file,
            bytes: &bytes,
            import,
            parent_id,
            sections: &sections,
            geometry: &geom,
        })
        .map_err(internal)?;
    eprintln!(
        "hull {id} ({file}): {} ({:.2} s)",
        if created { "saved" } else { "already saved" },
        t0.elapsed().as_secs_f64()
    );
    Ok(Reply::Json(json!({
        "id": id,
        "created": created,
        "hull": app.store.hull(id).map_err(internal)?,
    })))
}

fn case_params(v: &Value) -> Result<CaseParams, Fail> {
    let v = if v.is_null() { json!({}) } else { v.clone() };
    serde_json::from_value::<CaseParams>(v)
        .map_err(|e| e.to_string())
        .and_then(CaseParams::canonical)
        .map_err(bad)
}

/// The checks a case must pass before its statics are tried: a
/// catamaran doubles one hull, its demihulls clear of each other.
fn check_case(app: &App, hull: i64, p: &CaseParams) -> Result<(), Fail> {
    let h = app
        .store
        .hull(hull)
        .map_err(internal)?
        .ok_or(missing("hull"))?;
    let hulls = h["summary"]["hulls"].as_array().map_or(0, Vec::len);
    if let Some(span) = p.span {
        if hulls != 1 {
            return Err(bad(format!(
                "a catamaran doubles a single hull; this file holds {hulls}"
            )));
        }
        let beam = h["summary"]["hulls"][0]["beam"].as_f64().unwrap_or(0.0);
        if span <= beam {
            return Err(bad(format!(
                "span {span} m: the demihulls overlap (beam {beam:.3} m)"
            )));
        }
    }
    Ok(())
}

/// A case's statics, from its hull's file.
fn compute_statics(app: &App, hull: i64, p: &CaseParams) -> Result<Result<Value, String>, Fail> {
    let src = app
        .store
        .hull_source(hull)
        .map_err(internal)?
        .ok_or(missing("hull"))?;
    let bytes = app.store.blob(&src.file_blob).map_err(internal)?;
    let t0 = std::time::Instant::now();
    let r = statics(&src.file_name, bytes, &src.import, p);
    eprintln!(
        "statics on hull {hull}: {} ({:.2} s)",
        r.as_ref().map_or_else(|e| e.as_str(), |_| "ok"),
        t0.elapsed().as_secs_f64()
    );
    Ok(r)
}

/// The case on `hull` with these parameters: found, or made (its
/// statics computed first). `(id, created, statics error)`.
fn make_case(
    app: &App,
    hull: i64,
    p: &CaseParams,
    name: Option<&str>,
    by: &str,
) -> Result<(i64, bool, Option<String>), Fail> {
    if let Some(id) = app.store.find_case(hull, p).map_err(internal)? {
        let error = app
            .store
            .case(id)
            .map_err(internal)?
            .and_then(|c| c["error"].as_str().map(String::from));
        return Ok((id, false, error));
    }
    check_case(app, hull, p)?;
    let st = compute_statics(app, hull, p)?;
    let name = name.map(str::trim).unwrap_or("");
    let (id, created) = app
        .store
        .add_case(hull, p, name, by, st.as_ref().map_err(String::as_str))
        .map_err(internal)?;
    Ok((id, created, st.err()))
}

fn recompute(app: &App, id: i64) -> Result<(), Fail> {
    let (hull, p, _) = app
        .store
        .case_source(id)
        .map_err(internal)?
        .ok_or(missing("case"))?;
    let st = compute_statics(app, hull, &p)?;
    app.store
        .set_statics(id, st.as_ref().map_err(String::as_str))
        .map_err(internal)
}

/// `POST /api/studies`: every study on every case named (each made
/// already, its statics computed), all checked before anything is queued; a
/// study in waves brings the calm-water study it is taken about (with the grid
/// of a calm study asked for alongside, if any).
fn add_studies(app: &App, v: &Value, who: Option<&Who>) -> Result<Reply, Fail> {
    let s = &app.store;
    let cases: Vec<i64> = v["case_ids"]
        .as_array()
        .filter(|a| !a.is_empty())
        .ok_or(bad("case_ids: expected a list of case ids"))?
        .iter()
        .map(|c| {
            let id = c.as_i64().ok_or(bad(format!("case_ids: bad id {c}")))?;
            let cfg = s.case(id).map_err(internal)?.ok_or(missing("case"))?;
            if let Some(e) = cfg["error"].as_str() {
                return Err(bad(format!(
                    "case #{id} has no statics, so nothing can float on it: {e}"
                )));
            }
            Ok(id)
        })
        .collect::<Result<_, Fail>>()?;
    let asked: Vec<StudyParams> = v["studies"]
        .as_array()
        .filter(|a| !a.is_empty())
        .ok_or(bad("studies: expected a list of studies"))?
        .iter()
        .enumerate()
        .map(|(i, r)| {
            serde_json::from_value::<StudyParams>(r.clone())
                .map_err(|e| e.to_string())
                .and_then(StudyParams::canonical)
                .map_err(|e| bad(format!("study {}: {e}", i + 1)))
        })
        .collect::<Result<_, _>>()?;
    // The studies to make on each case: those asked for, each study in
    // waves after the calm-water study it needs (`implied` when not asked).
    let mut studies: Vec<(StudyParams, bool)> = Vec::new();
    for r in &asked {
        if r.waves.is_some() {
            let grid = asked
                .iter()
                .find(|o| o.waves.is_none() && o.attitude_key() == r.attitude_key())
                .map_or(640, |o| o.grid);
            let calm = r.calm(grid);
            if !studies
                .iter()
                .any(|(o, _)| o.waves.is_none() && o.attitude_key() == calm.attitude_key())
            {
                studies.push((calm, !asked.contains(&r.calm(grid))));
            }
        }
        if !studies.iter().any(|(o, _)| o == r) {
            studies.push((r.clone(), false));
        }
    }
    let total = cases.len() * studies.len();
    if total > 400 {
        return Err(bad(format!(
            "that is {total} studies; ask for at most 400 at a time"
        )));
    }

    if v["dry_run"].as_bool() == Some(true) {
        let mut rs = Vec::new();
        for &c in &cases {
            for (r, implied) in &studies {
                let found = s.find_study(c, r).map_err(internal)?;
                rs.push(json!({
                    "case_id": c,
                    "params": r,
                    "kind": r.kind(),
                    "implied": implied,
                    "existing": found.map(|(id, status, stale)| json!({
                        "id": id, "status": status, "stale": stale,
                    })),
                }));
            }
        }
        return Ok(Reply::Json(json!({ "studies": rs })));
    }

    let by = by(who, v["by"].as_str().unwrap_or(""));
    let priority = v["priority"].as_i64().unwrap_or(0);
    let mut rs = Vec::new();
    for &cid in &cases {
        for (r, implied) in &studies {
            let (rid, created) = s.add_study(cid, r, priority, by).map_err(internal)?;
            rs.push(json!({ "id": rid, "created": created, "case_id": cid, "implied": implied }));
        }
    }
    app.worker.poke();
    Ok(Reply::Json(json!({ "studies": rs })))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_store(tag: &str) -> Arc<Store> {
        let dir = std::env::temp_dir().join(format!(
            "boatmath-api-{tag}-{}-{}",
            std::process::id(),
            boatmath_web::store::now()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        Arc::new(Store::open(&dir).unwrap())
    }

    #[test]
    fn camber_document_uploads_and_views() {
        let store = temp_store("camber");
        let worker = Worker::new(Arc::clone(&store));
        let app = App { store, worker };

        let camber_path = format!(
            "{}/../boatmath/tests/camber/default.json",
            env!("CARGO_MANIFEST_DIR")
        );
        let bytes = std::fs::read(&camber_path).unwrap();

        // 1. Loft through api / shared loft
        let sections = loft("default.json", bytes.clone(), &LoftRequest::default()).unwrap();
        assert_eq!(sections["name"], "default.json");
        let geom = geometry("default.json", bytes.clone(), &LoftRequest::default()).unwrap();
        assert!(!geom["meshes"].as_array().unwrap().is_empty());

        // 2. Upload via api::route POST /api/hulls
        let mut body = {
            let b = bytes.clone();
            move || Ok(b.clone())
        };
        let q = vec![
            ("file".to_string(), "default.json".to_string()),
            ("name".to_string(), "camber_default".to_string()),
        ];
        let reply = route(&app, &Method::Post, "/api/hulls", &q, None, &mut body)
            .unwrap()
            .unwrap();
        let json_reply = match reply {
            Reply::Json(v) => v,
            _ => panic!("expected json"),
        };
        assert_eq!(json_reply["created"], true);
        let hull_id = json_reply["id"].as_i64().unwrap();

        // 3. GET /api/hulls/:id
        let reply = route(
            &app,
            &Method::Get,
            &format!("/api/hulls/{hull_id}"),
            &[],
            None,
            &mut || Ok(Vec::new()),
        )
        .unwrap()
        .unwrap();
        let h = match reply {
            Reply::Json(v) => v,
            _ => panic!("expected json"),
        };
        assert_eq!(h["name"], "camber_default");

        // 4. Create a case on this camber hull
        let case = CaseParams::default().canonical().unwrap();
        let statics = statics("default.json", bytes, &LoftRequest::default(), &case).unwrap();
        let (case_id, _) = app
            .store
            .add_case(hull_id, &case, "design", "test", Ok(&statics))
            .unwrap();
        assert!(case_id > 0);
    }
}
