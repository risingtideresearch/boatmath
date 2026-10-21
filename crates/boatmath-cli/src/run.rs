//! Running studies. A study whose result is already in the stream (or the
//! cache) isn't computed again. A study in waves is held at the attitude of
//! the calm-water study at its speed (as the web app's queue does), so
//! those calm-water studies run first, and join the stream if they weren't
//! in it.
//!
//! Each calm-water result is written with its sections and its field just
//! before it, so the stream stays whole.

use crate::cache::Cache;
use crate::records::{case_source, expect, sections, study_record, CaseSource};
use crate::stream::{short, Out, Stream};
use boatmath::params::{default_grid, StudyParams};
use boatmath::SOLVER_VERSION;
use boatmath::{
    flow_with_progress, wave_scalars, waves_with_progress, FlowRequest, Progress, ThrustLine,
};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::Mutex;
use std::time::Instant;

pub struct Options {
    pub jobs: usize,
    pub quiet: bool,
}

/// Run `f` on each item (with its index) on up to `jobs` threads.
fn pool<T: Sync>(items: &[T], jobs: usize, f: impl Fn(usize, &T) + Sync) {
    if items.is_empty() {
        return;
    }
    let next = Mutex::new(0usize);
    let workers = jobs.clamp(1, items.len().max(1));
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| {
                hullgeom::parallel::set_threads((hullgeom::parallel::available() / workers).max(1));
                loop {
                    let i = {
                        let mut n = next.lock().unwrap();
                        let idx = *n;
                        *n += 1;
                        idx
                    };
                    match items.get(i) {
                        Some(item) => f(i, item),
                        None => break,
                    }
                }
            });
        }
    });
}

fn params(study: &Value) -> Result<StudyParams, String> {
    serde_json::from_value(study["params"].clone()).map_err(|e| format!("study params: {e}"))
}

/// The calm-water study a study in waves is held at the attitude of.
fn calm_of(study: &Value) -> Result<Option<Value>, String> {
    let p = params(study)?;
    Ok(p.waves
        .is_some()
        .then(|| study_record(study["case"].as_str().unwrap_or(""), p.calm(default_grid()))))
}

#[derive(Clone, Debug)]
struct SolvedAttitude {
    case: String,
    params: StudyParams,
    attitude: (f64, f64),
}

/// What's known as the run goes: results computed (or found) so far, each
/// case's attitude at rest, and solved dynamic attitudes for warm-starting.
#[derive(Default)]
struct Known {
    results: HashMap<String, Value>,
    at_rest: HashMap<String, (f64, f64)>,
    solved: Vec<SolvedAttitude>,
}

impl Known {
    fn record(&mut self, study: &Value, result: &Value) {
        let id = study["id"].as_str().unwrap_or("");
        self.results.insert(id.to_string(), result.clone());

        let case_id = match study["case"].as_str() {
            Some(c) => c,
            None => return,
        };
        let p = match params(study) {
            Ok(p) => p,
            Err(_) => return,
        };
        if !p.dynamic || p.waves.is_some() {
            return;
        }
        if result["forces"]["solved"].as_bool() != Some(true) {
            return;
        }
        if let Some(att) = pair(&result["forces"], "sinkage", "trim_rad") {
            self.record_attitude(case_id, &p, att);
        }
    }

    fn record_attitude(&mut self, case_id: &str, p: &StudyParams, att: (f64, f64)) {
        if let Some(existing) = self
            .solved
            .iter_mut()
            .find(|s| s.case == case_id && s.params == *p)
        {
            existing.attitude = att;
        } else {
            self.solved.push(SolvedAttitude {
                case: case_id.to_string(),
                params: p.clone(),
                attitude: att,
            });
        }
    }

    /// The equilibrium of the nearest speed already solved on this
    /// case otherwise alike (`same_but_speed`), to start the solve from.
    fn warm_start(&self, case_id: &str, p: &StudyParams) -> Option<(f64, f64)> {
        if !p.dynamic || p.waves.is_some() {
            return None;
        }
        self.solved
            .iter()
            .filter(|s| {
                s.case == case_id && s.params.waves.is_none() && p.same_but_speed(&s.params)
            })
            .min_by(|a, b| {
                (a.params.froude - p.froude)
                    .abs()
                    .total_cmp(&(b.params.froude - p.froude).abs())
            })
            .map(|s| s.attitude)
    }
}

/// A result already in hand: in the stream, or in the cache with its field
/// and sections. Returns the records to write, the result last.
fn found(s: &Stream, cache: &Cache, id: &str) -> Option<Vec<Value>> {
    if let Some(r) = s
        .get("result", id)
        .filter(|r| r["solver_version"] == SOLVER_VERSION)
    {
        return Some(vec![r.clone()]);
    }
    let r = cache.get("result", id)?;
    let mut out = Vec::new();
    for key in ["sections", "field"] {
        if let Some(pid) = r[key].as_str() {
            out.push(cache.get(key, pid)?);
        }
    }
    out.push(r);
    Some(out)
}

/// Run the studies in the stream, writing it through and then each result
/// (after the records it refers to) as it's ready. Returns how many failed.
pub fn run<W: Write + Send>(
    s: &Stream,
    out: &mut Out<W>,
    cache: &Cache,
    opts: &Options,
) -> Result<usize, String> {
    out.pass(s)?;
    // Results are written in study order, each as soon as every one before
    // it is done, so a stream comes out the same however the jobs fall.
    struct Queue<'o, W: Write + Send> {
        out: &'o mut Out<W>,
        next: usize,
        ready: HashMap<usize, Vec<Value>>,
    }
    let queue = Mutex::new(Queue {
        out,
        next: 0,
        ready: HashMap::new(),
    });
    let emit = |rs: &[Value]| -> Result<(), String> {
        let mut q = queue.lock().unwrap();
        for r in rs {
            q.out.emit(r)?;
        }
        Ok(())
    };
    let finish = |i: usize, rs: Vec<Value>| -> Result<(), String> {
        let mut q = queue.lock().unwrap();
        q.ready.insert(i, rs);
        loop {
            let n = q.next;
            let Some(rs) = q.ready.remove(&n) else { break };
            for r in &rs {
                q.out.emit(r)?;
            }
            q.next += 1;
        }
        Ok(())
    };
    let known = Mutex::new(Known::default());
    let failed = Mutex::new(0usize);

    // First the calm-water studies, those asked for and those the studies
    // in waves are held about; then the studies in waves.
    let mut seen: HashSet<String> = HashSet::new();
    let mut first: Vec<Value> = Vec::new();
    let mut later: Vec<Value> = Vec::new();
    for study in s.of_type("study") {
        let id = study["id"].as_str().unwrap_or("").to_string();
        if !seen.insert(id.clone()) {
            continue;
        }
        match calm_of(study)? {
            None => first.push(study.clone()),
            Some(calm) => {
                let cid = calm["id"].as_str().unwrap_or("").to_string();
                if seen.insert(cid) {
                    // Into the stream, before its result.
                    emit(std::slice::from_ref(&calm))?;
                    first.push(calm);
                }
                later.push(study.clone());
            }
        }
    }

    // Pre-populate known with any studies already in stream or cache.
    for study in &first {
        let id = study["id"].as_str().unwrap_or("");
        if let Some(rs) = found(s, cache, id) {
            if let Some(r) = rs.last() {
                known.lock().unwrap().record(study, r);
            }
        }
    }

    let attempt = |seq: usize, study: &Value| {
        let id = study["id"].as_str().unwrap_or("");
        let records = match found(s, cache, id) {
            Some(rs) => Ok(rs),
            None => compute(s, cache, &known, study, opts.quiet),
        };
        let rs = match records {
            Ok(rs) => {
                if let Some(r) = rs.last() {
                    known.lock().unwrap().record(study, r);
                }
                rs
            }
            Err(e) => {
                eprintln!("boatmath: study {}: {e}", short(id));
                *failed.lock().unwrap() += 1;
                Vec::new()
            }
        };
        if let Err(e) = finish(seq, rs) {
            eprintln!("boatmath: {e}");
        }
    };
    pool(&first, opts.jobs, |i, st| attempt(i, st));
    let n = first.len();
    pool(&later, opts.jobs, |i, st| attempt(n + i, st));
    Ok(failed.into_inner().unwrap())
}

/// Each hull's entry less its display mesh.
fn without_meshes(hulls: &Value) -> Value {
    let mut h = hulls.clone();
    for e in h.as_array_mut().into_iter().flatten() {
        if let Some(o) = e.as_object_mut() {
            o.remove("mesh");
        }
    }
    h
}

fn pair(v: &Value, a: &str, b: &str) -> Option<(f64, f64)> {
    Some((v[a].as_f64()?, v[b].as_f64()?))
}

/// A case's attitude at rest: from its statics in the stream, or solved.
fn at_rest(s: &Stream, known: &Mutex<Known>, case: &Value) -> Result<(f64, f64), String> {
    let id = case["id"].as_str().unwrap_or("");
    if let Some(&a) = known.lock().unwrap().at_rest.get(id) {
        return Ok(a);
    }
    let a = match s
        .get("statics", id)
        .and_then(|st| pair(&st["at_rest"], "sinkage", "trim_rad"))
    {
        Some(a) => a,
        None => {
            let CaseSource { params, src, .. } = case_source(s, case)?;
            boatmath::at_rest(&src.file_name, src.bytes, &src.import, &params)?
        }
    };
    known.lock().unwrap().at_rest.insert(id.to_string(), a);
    Ok(a)
}

/// Compute one study: the records it makes, its result last.
fn compute(
    s: &Stream,
    cache: &Cache,
    known: &Mutex<Known>,
    study: &Value,
    quiet: bool,
) -> Result<Vec<Value>, String> {
    let id = expect(study, "study")?;
    let p = params(study)?;
    let case = s.follow(study, "case")?;
    let CaseSource {
        hull,
        params: cp,
        src,
        thrusts,
        ..
    } = case_source(s, case)?;

    let t0 = Instant::now();
    let label = format!(
        "{} {} Fn {}{}",
        short(id),
        hull["name"].as_str().unwrap_or(""),
        p.froude,
        p.waves
            .as_ref()
            .map_or(String::new(), |w| format!(" waves {}°", w.heading))
    );
    let mut last = "";
    let mut report = |pr: &Progress| {
        if !quiet && pr.stage != last {
            eprintln!(
                "boatmath: {label}: {} ({:.0}%)",
                pr.stage.to_lowercase(),
                100.0 * pr.fraction
            );
            last = pr.stage;
        }
        true
    };

    let mut records = Vec::new();
    let mut r = match &p.waves {
        None => {
            let hold = if p.dynamic {
                None
            } else {
                Some(at_rest(s, known, case)?)
            };
            let warm = if p.dynamic {
                known
                    .lock()
                    .unwrap()
                    .warm_start(case["id"].as_str().unwrap_or(""), &p)
            } else {
                None
            };
            // With a drive and floating at speed, the attitude is the
            // self-propelled one: towed first, for the resistance the
            // thrust must meet, then with that thrust along the drives'
            // line. (Held at rest, the attitude is given.)
            let angle = thrusts.first().map_or(0.0, |t| t.angle);
            let mut propulsion = None;
            let (thrust, warm) = if p.dynamic && !thrusts.is_empty() {
                let towed = FlowRequest {
                    cut: src.import,
                    case: cp.clone(),
                    study: p.clone(),
                    warm,
                    hold: None,
                    thrust: None,
                    field: false,
                };
                let v0 =
                    flow_with_progress(&src.file_name, src.bytes.clone(), &towed, &mut report)?;
                let rt0 = v0["forces"]["rt"]
                    .as_f64()
                    .ok_or("the towed pass has no R_t")?;
                let t = rt0 / angle.cos();
                propulsion = Some(json!({
                    "thrust": t,
                    "shaft_angle_deg": angle.to_degrees(),
                    "drives": thrusts,
                    "towed": {
                        "rt": rt0,
                        "sinkage": v0["forces"]["sinkage"],
                        "trim_deg": v0["forces"]["trim_deg"],
                    },
                }));
                let towed_att = pair(&v0["forces"], "sinkage", "trim_rad");
                (
                    Some(ThrustLine {
                        thrust: t,
                        at: thrusts.iter().map(|d| (d.x, d.z)).collect(),
                        angle,
                    }),
                    towed_att.or(warm),
                )
            } else {
                (None, warm)
            };
            let req = FlowRequest {
                cut: src.import,
                case: cp,
                study: p.clone(),
                warm,
                hold,
                thrust,
                field: true,
            };
            let v = flow_with_progress(&src.file_name, src.bytes, &req, &mut report)?;
            if p.dynamic && v["forces"]["solved"].as_bool() == Some(true) {
                if let Some(att) = pair(&v["forces"], "sinkage", "trim_rad") {
                    known.lock().unwrap().record_attitude(
                        case["id"].as_str().unwrap_or(""),
                        &p,
                        att,
                    );
                }
            }
            // The meshes are for drawing; the pictures re-pose the hull's
            // geometry at the attitude instead of keeping them.
            let field = json!({
                "type": "field",
                "id": id,
                "study": id,
                "hulls": without_meshes(&v["hulls"]),
                "surface": v["surface"],
                "solver_version": SOLVER_VERSION,
            });
            let sec = match pair(&v["forces"], "sinkage", "trim_rad") {
                Some(a) => Some(sections(s, case, a, cache)?),
                None => None,
            };
            let r = json!({
                "kind": "calm",
                "froude": v["froude"],
                "speed": v["speed"],
                "transverse_wavelength": v["transverse_wavelength"],
                "seconds": v["seconds"],
                "forces": v["forces"],
                "self_propelled": propulsion,
                "field": id,
                "sections": sec.as_ref().map(|x| x["id"].clone()),
            });
            if let Some(sec) = sec {
                records.push(sec);
            }
            cache.put(&field);
            records.push(field);
            r
        }
        Some(_) => {
            let calm = calm_of(study)?.expect("a study in waves");
            let calm_id = calm["id"].as_str().unwrap_or("");
            let cr = known
                .lock()
                .unwrap()
                .results
                .get(calm_id)
                .cloned()
                .ok_or_else(|| format!("its calm-water study {} has no result", short(calm_id)))?;
            let hold = pair(&cr["forces"], "sinkage", "trim_rad")
                .ok_or("its calm-water study has no attitude")?;
            let v = waves_with_progress(
                &src.file_name,
                src.bytes,
                &src.import,
                &cp,
                &p,
                hold,
                &mut report,
            )?;
            json!({
                "kind": "waves",
                "froude": v["froude"],
                "speed": v["speed"],
                "seconds": v["seconds"],
                "attitude": v["attitude"],
                "calm": calm_id,
                "sections": cr["sections"],
                "peaks": wave_scalars(&v),
                "seakeeping": v["seakeeping"],
            })
        }
    };
    let o = r.as_object_mut().expect("an object");
    o.insert("type".into(), json!("result"));
    o.insert("id".into(), json!(id));
    o.insert("study".into(), json!(id));
    o.insert("solver_version".into(), json!(SOLVER_VERSION));
    cache.put(&r);
    records.push(r);
    if !quiet {
        eprintln!(
            "boatmath: {label}: done in {:.1} s",
            t0.elapsed().as_secs_f64()
        );
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use boatmath::params::Closure;

    fn test_study_params(froude: f64, dynamic: bool, grid: usize) -> StudyParams {
        StudyParams {
            froude,
            dynamic,
            closure: Closure::default(),
            grid,
            astern: None,
            waves: None,
        }
    }

    #[test]
    fn warm_start_picks_nearest_froude_on_same_case() {
        let mut known = Known::default();
        let p_fn2 = test_study_params(0.2, true, 40);
        let p_fn4 = test_study_params(0.4, true, 40);

        let st_c1_fn2 = json!({
            "type": "study",
            "id": "s1",
            "case": "c1",
            "params": p_fn2,
        });
        let res_c1_fn2 = json!({
            "type": "result",
            "id": "s1",
            "study": "s1",
            "forces": {
                "solved": true,
                "sinkage": 0.05,
                "trim_rad": 0.01,
            }
        });

        // Record c1 at Fn 0.2
        known.record(&st_c1_fn2, &res_c1_fn2);

        // Querying c1 at Fn 0.3 should find Fn 0.2
        let p_fn3 = test_study_params(0.3, true, 40);
        assert_eq!(known.warm_start("c1", &p_fn3), Some((0.05, 0.01)));

        // Querying c2 should return None
        assert_eq!(known.warm_start("c2", &p_fn3), None);

        // Querying non-dynamic study should return None
        let p_fn3_static = test_study_params(0.3, false, 40);
        assert_eq!(known.warm_start("c1", &p_fn3_static), None);

        // Record c1 at Fn 0.4
        let st_c1_fn4 = json!({
            "type": "study",
            "id": "s2",
            "case": "c1",
            "params": p_fn4,
        });
        let res_c1_fn4 = json!({
            "type": "result",
            "id": "s2",
            "study": "s2",
            "forces": {
                "solved": true,
                "sinkage": 0.08,
                "trim_rad": 0.02,
            }
        });
        known.record(&st_c1_fn4, &res_c1_fn4);

        // Fn 0.35 is closer to Fn 0.4 (diff 0.05 vs 0.15)
        let p_fn35 = test_study_params(0.35, true, 40);
        assert_eq!(known.warm_start("c1", &p_fn35), Some((0.08, 0.02)));

        // Fn 0.25 is closer to Fn 0.2 (diff 0.05 vs 0.15)
        let p_fn25 = test_study_params(0.25, true, 40);
        assert_eq!(known.warm_start("c1", &p_fn25), Some((0.05, 0.01)));
    }

    #[test]
    fn run_prepopulates_and_streams_in_order() {
        let p1 = test_study_params(0.2, true, 40);
        let p2 = test_study_params(0.3, true, 40);

        let case = json!({
            "type": "case",
            "id": "c1",
            "hull": "h1",
            "mass": 1000.0,
        });
        let st1 = json!({
            "type": "study",
            "id": "s1",
            "case": "c1",
            "params": p1,
        });
        let st2 = json!({
            "type": "study",
            "id": "s2",
            "case": "c1",
            "params": p2,
        });
        let res1 = json!({
            "type": "result",
            "id": "s1",
            "study": "s1",
            "solver_version": SOLVER_VERSION,
            "forces": {
                "solved": true,
                "sinkage": 0.05,
                "trim_rad": 0.01,
            }
        });
        let res2 = json!({
            "type": "result",
            "id": "s2",
            "study": "s2",
            "solver_version": SOLVER_VERSION,
            "forces": {
                "solved": true,
                "sinkage": 0.07,
                "trim_rad": 0.015,
            }
        });

        let mut s = Stream::default();
        s.push(case);
        s.push(st1);
        s.push(st2);
        s.push(res1);
        s.push(res2);

        let mut out_bytes = Vec::new();
        let mut out = Out::new(&mut out_bytes);
        let cache = Cache::default();
        let failed = run(
            &s,
            &mut out,
            &cache,
            &Options {
                jobs: 2,
                quiet: true,
            },
        )
        .unwrap();

        assert_eq!(failed, 0);

        let mut out_stream = Stream::default();
        for line in String::from_utf8(out_bytes).unwrap().lines() {
            if !line.trim().is_empty() {
                let rec: Value = serde_json::from_str(line).unwrap();
                out_stream.push(rec);
            }
        }

        let results: Vec<&Value> = out_stream.of_type("result").collect();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0]["id"], "s1");
        assert_eq!(results[1]["id"], "s2");
    }
}
