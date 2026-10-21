//! The port against propopt's web cores: each golden file (made by
//! `tests/golden/make.js`) holds the original's sweep and motor ranking for
//! one configuration, and this crate must give the same.

use propeller::bseries::{sweep, Config, Inputs, SweepOptions, TopInputs};
use propeller::motor::{scatter, Database, Options, RankBy};
use serde_json::Value;

fn load(name: &str) -> Value {
    let path = format!("{}/tests/golden/{name}.json", env!("CARGO_MANIFEST_DIR"));
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn inputs(o: &Value) -> Inputs {
    let f = |k: &str| o[k].as_f64().unwrap();
    let top = match (o["V_top"].as_f64(), o["T_top"].as_f64()) {
        (Some(s), Some(t)) => Some(TopInputs {
            speed: s,
            thrust: t,
        }),
        _ => None,
    };
    Inputs {
        speed: f("V_s"),
        thrust: f("T"),
        shafts: f("shafts") as u32,
        wake: f("wake"),
        thrust_deduction: f("thrustDeduction"),
        d_min: f("D_min"),
        d_max: f("D_max"),
        blades: o["Z"]
            .as_array()
            .unwrap()
            .iter()
            .map(|z| z.as_u64().unwrap() as u32)
            .collect(),
        depth: f("depth"),
        keller_k: f("kellerK"),
        cavitation: o["cavitation"].as_bool().unwrap(),
        re_correct: o["reCorrect"].as_bool().unwrap(),
        strict_ear: o["strictEAR"].as_bool().unwrap(),
        top,
        shaft_angle_deg: 0.0,
    }
}

fn close(what: &str, got: f64, want: f64, rel: f64) {
    let tol = rel * want.abs().max(1e-9);
    assert!(
        (got - want).abs() <= tol,
        "{what}: got {got}, want {want} (off by {:.2e} relative)",
        (got - want).abs() / want.abs().max(1e-300)
    );
}

fn design(what: &str, got: &Value, want: &Value) {
    for k in [
        "rpm", "D", "PD", "EAR", "Z", "eta0", "P_shaft", "Q", "J", "K_T", "K_Q",
    ] {
        close(
            &format!("{what} {k}"),
            got[k].as_f64().unwrap(),
            want[k].as_f64().unwrap(),
            1e-6,
        );
    }
    for k in ["kellerOk", "burrillOk", "topOk", "earInSeries"] {
        assert_eq!(got[k], want[k], "{what} {k}");
    }
}

fn check(name: &str) {
    let g = load(name);
    let c = Config::new(&inputs(&g["opts"]));
    let s = sweep(
        &c,
        &SweepOptions {
            power_cap: g["powerCap"].as_f64().unwrap(),
            ..SweepOptions::default()
        },
    );
    let s = serde_json::to_value(&s).unwrap();
    let w = &g["sweep"];
    assert_eq!(s["feasible"], w["feasible"]);
    for k in ["rpmLo", "rpmHi", "feasLo", "feasHi"] {
        close(
            &format!("{name} {k}"),
            s[k].as_f64().unwrap(),
            w[k].as_f64().unwrap(),
            1e-9,
        );
    }
    design(&format!("{name} best"), &s["best"], &w["best"]);
    if !w["bestUnconstrained"].is_null() {
        design(
            &format!("{name} bestUnconstrained"),
            &s["bestUnconstrained"],
            &w["bestUnconstrained"],
        );
    }
    assert_eq!(s["topOkCount"], w["topOkCount"], "{name} topOkCount");
    let (sp, wp) = (
        s["points"].as_array().unwrap(),
        w["points"].as_array().unwrap(),
    );
    assert_eq!(sp.len(), wp.len());
    for (i, (a, b)) in sp.iter().zip(wp).enumerate() {
        assert_eq!(a["ok"], b["ok"], "{name} point {i} ok");
        if b["ok"] == true {
            design(&format!("{name} point {i}"), a, b);
        }
    }

    // Motors, ranked by electrical power, one per family and every winding.
    let db = Database::vendored();
    let pts: Vec<propeller::bseries::CurvePoint> =
        serde_json::from_value(s["points"].clone()).unwrap();
    let o = Options {
        controller_eta: db.controller_eta,
        check_top: c.top.is_some(),
        ..Options::default()
    };
    for (key, all) in [("dots", false), ("all", true)] {
        let sc = scatter(db, &pts, &o, |_| true, RankBy::Power, all);
        let want = g["motors"][key].as_array().unwrap();
        let got = if all { &sc.all } else { &sc.dots };
        let ids: Vec<&str> = got.iter().map(|d| d.motor.id.as_str()).collect();
        let wids: Vec<&str> = want.iter().map(|d| d["motor"].as_str().unwrap()).collect();
        assert_eq!(ids, wids, "{name} {key}: motor order");
        for (d, w) in got.iter().zip(want) {
            let id = &d.motor.id;
            close(
                &format!("{name} {id} P_elec"),
                d.run.P_elec,
                w["P_elec"].as_f64().unwrap(),
                1e-6,
            );
            close(
                &format!("{name} {id} ratio"),
                d.run.ratio,
                w["ratio"].as_f64().unwrap(),
                1e-6,
            );
            close(
                &format!("{name} {id} rpm"),
                d.point.rpm,
                w["rpm"].as_f64().unwrap(),
                1e-9,
            );
            close(
                &format!("{name} {id} eta"),
                d.run.eta,
                w["eta"].as_f64().unwrap(),
                1e-6,
            );
        }
    }
    let unreachable: Vec<&str> = scatter(db, &pts, &o, |_| true, RankBy::Power, false)
        .unreachable
        .iter()
        .map(|m| m.id.as_str())
        .collect();
    let wu: Vec<&str> = g["motors"]["unreachable"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(unreachable, wu, "{name} unreachable");
}

#[test]
fn web_default() {
    check("web_default");
}

#[test]
fn top_speed() {
    check("top_speed");
}

#[test]
fn small_twin() {
    check("small_twin");
}

#[test]
fn the_vendored_database_reads() {
    let db = Database::vendored();
    assert_eq!(db.motors.len(), 159);
    assert_eq!(db.maps.len(), 2);
    assert!((db.controller_eta - 0.97).abs() < 1e-12);
}
