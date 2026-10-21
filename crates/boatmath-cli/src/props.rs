//! Propellers and motors as records: a `prop` is the best B-series
//! propeller for an operating point, with its rpm curve; a `drive` is one
//! motor driving it, at its best reduction and point on that curve; a
//! `motor` is a database entry. See the `propeller` crate for the models.

use crate::cache::Cache;
use crate::records::{expect, study_record};
use crate::stream::{id_of, short, Stream};
use boatmath::params::StudyParams;
use boatmath::SOLVER_VERSION;
use propeller::bseries::{sweep, Config, CurvePoint, Inputs, SweepOptions};
use propeller::motor::{scatter, Database, Dot, Motor, Options, RankBy, Rating};
use serde_json::{json, Value};

/// A calm-water result's operating point: its speed [m/s] and the thrust
/// [N] its resistance asks for, `R_t / (1 − t)`, and its hulls.
pub fn operating_point(r: &Value, thrust_deduction: f64) -> Result<(f64, f64), String> {
    let id = expect(r, "result")?;
    if r["kind"] != "calm" {
        return Err(format!("result {}: not a calm-water result", short(id)));
    }
    let speed = r["speed"].as_f64().ok_or("a result without a speed")?;
    let rt = r["forces"]["rt"].as_f64().ok_or("a result without R_t")?;
    Ok((speed, rt / (1.0 - thrust_deduction)))
}

/// The number of hulls a result's platform floats on: a shaft each.
pub fn hull_count(s: &Stream, r: &Value) -> usize {
    let case = s
        .follow(r, "study")
        .and_then(|st| s.follow(st, "case"))
        .ok();
    match case {
        Some(c) if !c["params"]["span"].is_null() => 2,
        Some(c) => s
            .follow(c, "hull")
            .ok()
            .and_then(|h| h["geometry"]["hulls"].as_array().map(Vec::len))
            .unwrap_or(1),
        None => 1,
    }
}

/// A result's drives, from its case's mount: where each propeller sits at
/// the result's attitude, the shafts' angle to the flow (their own plus the
/// running trim), and the drive's own propeller, if it comes with one.
pub struct Drives {
    /// Each propeller's centre `(x, y, depth)` at the attitude [m].
    pub at: Vec<(f64, f64, f64)>,
    /// The thrust line's angle to the flow [rad].
    pub angle: f64,
    pub prop_diameter: Option<f64>,
    /// The thrust line, as the prop record keeps it.
    pub record: Value,
}

pub fn drives_of(s: &Stream, r: &Value) -> Result<Option<Drives>, String> {
    let case = s.follow(s.follow(r, "study")?, "case")?;
    let cs = crate::records::case_source(s, case)?;
    let Some(m) = &cs.mount else { return Ok(None) };
    let f = &r["forces"];
    let attitude = (
        f["sinkage"].as_f64().unwrap_or(0.0),
        f["trim_rad"].as_f64().unwrap_or(0.0),
    );
    let mounts = boatmath::propulsion::mounts(
        &cs.src.file_name,
        cs.src.bytes.clone(),
        &cs.src.import,
        &cs.params,
        attitude,
    )?;
    let at: Vec<(f64, f64, f64)> = cs
        .thrusts
        .iter()
        .map(|t| {
            let d = mounts.disc_at(t.x, t.y, t.z, 0.0);
            (d.x, d.y, d.depth)
        })
        .collect();
    let angle = m.shaft_angle_deg.to_radians() + attitude.1;
    let record = json!({
        "mount": m,
        "drives": at.iter().map(|(x, y, d)| json!({ "x": x, "y": y, "depth": d })).collect::<Vec<_>>(),
        "angle_deg": angle.to_degrees(),
        // The cross flow loads the blades cyclically and brings cavitation on
        // early; the B-series sees only the axial part.
        "oblique": angle.to_degrees().abs() > 3.0,
    });
    Ok(Some(Drives {
        at,
        angle,
        prop_diameter: m.prop_diameter,
        record,
    }))
}

/// The result of the same study at Froude number `froude`: the second
/// operating point, which must be in the stream.
pub fn result_at_froude(s: &Stream, r: &Value, froude: f64) -> Result<Value, String> {
    let study = s.follow(r, "study")?;
    let mut p: StudyParams = serde_json::from_value(study["params"].clone())
        .map_err(|e| format!("study params: {e}"))?;
    p.froude = froude;
    let other = study_record(study["case"].as_str().unwrap_or(""), p.canonical()?);
    let id = other["id"].as_str().unwrap_or("");
    s.get("result", id)
        .filter(|r| r["solver_version"] == SOLVER_VERSION)
        .cloned()
        .ok_or_else(|| {
            format!(
                "no result at Fn {froude} for case {} in the stream: run that study too \
                 (boatmath study --froude {froude},… | boatmath run)",
                short(study["case"].as_str().unwrap_or(""))
            )
        })
}

/// The best propeller for `inputs`, and its curve, as a record (from the
/// cache when it's there).
pub fn prop(
    cache: &Cache,
    inputs: &Inputs,
    result: Option<&str>,
    power_cap: f64,
    interaction: Option<&Value>,
    thrust_line: Option<&Value>,
) -> Result<Value, String> {
    let id = id_of(&json!({
        "inputs": inputs, "result": result, "power_cap": power_cap,
        "interaction": interaction.map(|i| &i["settings"]),
        "thrust_line": thrust_line,
    }));
    if let Some(old) = cache.get("prop", &id) {
        return Ok(old);
    }
    let c = Config::new(inputs);
    let s = sweep(
        &c,
        &SweepOptions {
            power_cap,
            ..SweepOptions::default()
        },
    );
    let qpc = s.best.as_ref().map(|b| b.eta0 * c.eta_h);
    let r = json!({
        "type": "prop",
        "id": id,
        "result": result,
        "inputs": inputs,
        // Per shaft, as the propeller sees it; and the whole boat's.
        "shaft": { "V_A": c.v_a, "T": c.t_target, "top_T": c.top.map(|t| t.t_target) },
        "boat": { "R_T": c.r_t, "P_E": c.p_e, "eta_H": c.eta_h, "QPC": qpc },
        "feasible": s.feasible,
        "best": s.best,
        "best_unconstrained": s.bestUnconstrained,
        "top_feasible": s.topFeasible,
        "window": { "rpm_lo": s.rpmLo, "rpm_hi": s.rpmHi },
        "feasible_rpm": { "lo": s.feasLo, "hi": s.feasHi },
        "curve": s.points,
        "interaction": interaction,
        "thrust_line": thrust_line,
        "solver_version": SOLVER_VERSION,
    });
    cache.put(&r);
    Ok(r)
}

/// Which hull–propeller coefficients to estimate, and where the
/// propellers sit.
pub struct Auto {
    pub wake: bool,
    pub deduction: bool,
}

/// A result's propeller with its wake fraction and/or thrust deduction
/// from potential flow (`thinship::propulsion`): the propeller found with a
/// guess, its disc placed on each hull at the result's attitude, w and t
/// computed there, and the propeller found again, until they settle.
pub fn auto_prop(
    s: &Stream,
    cache: &Cache,
    r: &Value,
    mut inputs: Inputs,
    auto: &Auto,
    power_cap: f64,
) -> Result<Value, String> {
    use boatmath::params::StudyParams;
    let result_id = expect(r, "result")?;
    let study = s.follow(r, "study")?;
    let case = s.follow(study, "case")?;
    let cs = crate::records::case_source(s, case)?;
    let (src, cp) = (cs.src, cs.params);
    let sp: StudyParams = serde_json::from_value(study["params"].clone())
        .map_err(|e| format!("study params: {e}"))?;
    let f = &r["forces"];
    let attitude = (
        f["sinkage"].as_f64().unwrap_or(0.0),
        f["trim_rad"].as_f64().unwrap_or(0.0),
    );
    let speed = r["speed"].as_f64().ok_or("a result without a speed")?;
    let rt = f["rt"].as_f64().ok_or("a result without R_t")?;
    let placed = boatmath::propulsion::placed(
        &src.file_name,
        src.bytes,
        &src.import,
        &cp,
        attitude,
        speed,
        sp.closure.transom(),
    )?;
    let drives = drives_of(s, r)?
        .ok_or("--wake auto and --thrust-deduction auto need a drive: give the case a mount (`boatmath mount`)")?;
    let thrusts = crate::records::case_source(s, case)?.thrusts;
    let cos_e = drives.angle.cos();
    let (mut w, mut t) = (
        if auto.wake { 0.0 } else { inputs.wake },
        if auto.deduction {
            0.0
        } else {
            inputs.thrust_deduction
        },
    );
    let mut steps = Vec::new();
    let mut last = None;
    let mut converged = false;
    for _ in 0..8 {
        inputs.wake = w;
        inputs.thrust_deduction = t;
        inputs.thrust = rt / ((1.0 - t) * cos_e);
        inputs.shaft_angle_deg = drives.angle.to_degrees();
        let c = Config::new(&inputs);
        let s = sweep(
            &c,
            &SweepOptions {
                power_cap,
                ..SweepOptions::default()
            },
        );
        let best = s
            .best
            .ok_or("no B-series propeller can do it (try a larger --d-max, or --no-cavitation)")?;
        let discs = placed.mounts.mount_discs(&thrusts, 0.5 * best.D);
        let wakes = discs
            .iter()
            .map(|d| placed.interaction.wake(d))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| e.to_string())?;
        let w_new = wakes.iter().map(|k| k.w).sum::<f64>() / wakes.len() as f64;
        let v_a = speed * (1.0 - if auto.wake { w_new } else { w });
        let ded = placed
            .interaction
            .thrust_deduction(&discs, inputs.thrust / discs.len() as f64, v_a)
            .map_err(|e| e.to_string())?;
        steps.push(json!({ "D": best.D, "w": w_new, "t": ded.t }));
        let (w_next, t_next) = (
            if auto.wake { w_new } else { w },
            if auto.deduction { ded.t } else { t },
        );
        if !(-0.5..0.9).contains(&w_next) || !(-0.5..0.9).contains(&t_next) {
            return Err(format!(
                "the propeller's wake {w_next:.3} and thrust deduction {t_next:.3} are out of range: is the mount's propeller inside the hull?"
            ));
        }
        let done = (w_next - w).abs() < 5e-4 && (t_next - t).abs() < 5e-4;
        w = w_next;
        t = t_next;
        last = Some((wakes, ded, discs));
        if done {
            converged = true;
            break;
        }
    }
    let (wakes, ded, discs) = last.expect("one step at least");
    inputs.wake = w;
    inputs.thrust_deduction = t;
    inputs.thrust = rt / ((1.0 - t) * cos_e);
    inputs.shaft_angle_deg = drives.angle.to_degrees();
    let n = wakes.len() as f64;
    let interaction = json!({
        "model": "potential",
        "settings": { "wake": auto.wake, "thrust_deduction": auto.deduction },
        "w": wakes.iter().map(|k| k.w).sum::<f64>() / n,
        "w_local": wakes.iter().map(|k| k.w_local).sum::<f64>() / n,
        "w_wave": wakes.iter().map(|k| k.w_wave).sum::<f64>() / n,
        "w_friction": wakes.iter().map(|k| k.w_friction).sum::<f64>() / n,
        "radial": wakes[0].radial,
        "t": ded.t,
        "delta_r": ded.delta_r,
        "u_a": ded.u_a,
        "discs": discs.iter().map(|d| json!({ "x": d.x, "y": d.y, "depth": d.depth, "radius": d.radius, "hub": d.hub })).collect::<Vec<_>>(),
        "steps": steps,
        "converged": converged,
    });
    prop(
        cache,
        &inputs,
        Some(result_id),
        power_cap,
        Some(&interaction),
        Some(&drives.record),
    )
}

// ------------------------------------------------------------- motors

/// Which motors to consider, as the web page's filters.
#[derive(Default, Clone)]
pub struct Filter {
    pub vendors: Vec<String>,
    pub mapped_only: bool,
    pub single_only: bool,
    pub max_mass: Option<f64>,
    pub max_od: Option<f64>,
    pub max_price: Option<f64>,
}

impl Filter {
    /// A size or price limit is hard: a motor that doesn't publish the
    /// number can't meet it.
    pub fn admits(&self, m: &Motor) -> bool {
        (self.vendors.is_empty()
            || self
                .vendors
                .iter()
                .any(|v| v.eq_ignore_ascii_case(&m.vendor)))
            && (!self.mapped_only || matches!(m.tier.as_deref(), Some("A") | Some("A-")))
            && (!self.single_only || m.stack.is_none_or(|s| s == 0.0))
            && self
                .max_mass
                .is_none_or(|l| m.mass_kg.is_some_and(|v| v <= l))
            && self.max_od.is_none_or(|l| m.od_mm.is_some_and(|v| v <= l))
            && self
                .max_price
                .is_none_or(|l| m.price_usd.is_some_and(|v| v <= l))
    }
}

/// A motor's database entry as a record.
pub fn motor_record(m: &Motor) -> Value {
    let mut r = m.record.clone();
    r["type"] = json!("motor");
    r
}

pub struct MatchOptions {
    pub filter: Filter,
    pub rank_by: RankBy,
    pub all_windings: bool,
    /// A fixed reduction (1: direct drive), else the best for each motor.
    pub ratio: Option<f64>,
    pub ratio_max: f64,
    pub gear_eta: f64,
    pub controller_loss: bool,
    pub peak: bool,
    /// Require the prop's second operating point, through the same gearbox.
    pub check_top: bool,
}

/// Every admitted motor at its best point on `prop`'s curve, ranked, as
/// `drive` records; and the ids of those that can't drive it.
pub fn drives(
    db: &Database,
    prop: &Value,
    o: &MatchOptions,
) -> Result<(Vec<Value>, Vec<String>), String> {
    let prop_id = expect(prop, "prop")?;
    let curve: Vec<CurvePoint> = serde_json::from_value(prop["curve"].clone())
        .map_err(|e| format!("prop {}: curve: {e}", short(prop_id)))?;
    let has_top = !prop["inputs"]["top"].is_null();
    let mo = Options {
        rating: if o.peak { Rating::Peak } else { Rating::Cont },
        gear_eta: o.gear_eta,
        controller_eta: if o.controller_loss {
            db.controller_eta
        } else {
            1.0
        },
        ratio: o.ratio,
        ratio_max: o.ratio_max,
        check_top: o.check_top && has_top,
        ..Options::default()
    };
    let sc = scatter(
        db,
        &curve,
        &mo,
        |m| o.filter.admits(m),
        o.rank_by,
        o.all_windings,
    );
    let settings = json!({
        "filter": {
            "vendors": o.filter.vendors, "mapped_only": o.filter.mapped_only,
            "single_only": o.filter.single_only, "max_mass": o.filter.max_mass,
            "max_od": o.filter.max_od, "max_price": o.filter.max_price,
        },
        "ratio": o.ratio, "ratio_max": o.ratio_max, "gear_eta": o.gear_eta,
        "controller_eta": mo.controller_eta, "rating": if o.peak { "peak" } else { "cont" },
        "check_top": mo.check_top,
    });
    let mut out = Vec::new();
    for (rank, d) in sc.dots.iter().enumerate() {
        let r = drive_record(prop_id, d, rank + 1, &settings);
        out.push(r);
    }
    Ok((out, sc.unreachable.iter().map(|m| m.id.clone()).collect()))
}

fn drive_record(prop_id: &str, d: &Dot, rank: usize, settings: &Value) -> Value {
    let p = d.point.design.as_ref();
    let run = &d.run;
    json!({
        "type": "drive",
        "id": id_of(&json!({ "prop": prop_id, "motor": d.motor.id, "settings": settings })),
        "prop": prop_id,
        "motor": d.motor.id,
        "rank": rank,
        "settings": settings,
        // The propeller it runs: the curve's point at this shaft speed.
        "rpm": d.point.rpm,
        "P_shaft": p.map(|p| p.P_shaft),
        "propeller": p.map(|p| json!({ "D": p.D, "PD": p.PD, "EAR": p.EAR, "Z": p.Z, "eta0": p.eta0 })),
        // The motor, through its reduction.
        "P_elec": run.P_elec,
        "eta": run.eta,
        "etaWithController": run.etaWithController,
        "etaDrive": run.etaDrive,
        "ratio": run.ratio,
        "gearEta": run.gearEta,
        "motorRpm": run.rpm,
        "motorQ": run.Q,
        "V_bus": run.V_bus,
        "I_dc": run.I_dc,
        "I_arms": run.I_arms,
        "load": d.load,
        "extrapolated": run.extrapolated,
        "top": run.top,
        "tier": d.motor.tier,
    })
}
