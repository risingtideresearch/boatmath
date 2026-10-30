use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

fn boatmath_bin() -> &'static str {
    env!("CARGO_BIN_EXE_boatmath")
}

fn run_cmd(args: &[&str], stdin_bytes: Option<&[u8]>) -> std::process::Output {
    let mut cmd = Command::new(boatmath_bin());
    cmd.args(args).stdout(Stdio::piped()).stderr(Stdio::piped());
    if stdin_bytes.is_some() {
        cmd.stdin(Stdio::piped());
    } else {
        cmd.stdin(Stdio::null());
    }
    let mut child = cmd.spawn().expect("spawn boatmath");
    if let Some(input) = stdin_bytes {
        child.stdin.as_mut().unwrap().write_all(input).unwrap();
    }
    child.wait_with_output().unwrap()
}

fn parse_jsonl(bytes: &[u8]) -> Vec<Value> {
    serde_json::Deserializer::from_slice(bytes)
        .into_iter::<Value>()
        .map(|r| r.expect("valid JSON record"))
        .collect()
}

#[test]
fn outright_tonnage_modes_units_and_24m_limit() {
    // 1. Outright monohull with units and canonical ID equivalence vs --hull-dims.
    let out1 = run_cmd(
        &[
            "tonnage", "--tml", "12m", "--tmb", "360cm", "--tmd", "1800mm",
        ],
        None,
    );
    assert!(
        out1.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out1.stderr)
    );
    let recs1 = parse_jsonl(&out1.stdout);
    assert_eq!(recs1.len(), 1);
    assert_eq!(recs1[0]["type"], "tonnage");
    assert_eq!(recs1[0]["gt"], 12.44);
    assert_eq!(recs1[0]["nt"], 9.33);
    assert_eq!(recs1[0]["form"], "4A");

    let out2 = run_cmd(&["tonnage", "--hull-dims", "12x3.6x1.8"], None);
    assert!(out2.status.success());
    let recs2 = parse_jsonl(&out2.stdout);
    assert_eq!(recs1[0]["id"], recs2[0]["id"]);

    // 2. Outright asymmetric trimaran with bridge and superstructure.
    let tri = run_cmd(
        &[
            "tonnage",
            "--kind",
            "sail",
            "--hull-dims",
            "12x1.8x1.6",
            "--hull-dims",
            "9x0.6x0.7",
            "--hull-dims",
            "9x0.6x0.7",
            "--bridge",
            "5x3x0.9",
            "--space",
            "4x2.5x1.9",
        ],
        None,
    );
    assert!(
        tri.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&tri.stderr)
    );
    let tri_rec = &parse_jsonl(&tri.stdout)[0];
    assert_eq!(tri_rec["hull_count"], 3);
    assert_eq!(tri_rec["above_deck_exempt"], false);
    assert!(tri_rec["bridge_tonnage"].as_f64().unwrap() > 0.0);
    assert!(tri_rec["above_deck_tonnage"].as_f64().unwrap() > 0.0);

    // 3. >= 24 m statutory length check.
    let too_long = run_cmd(
        &["tonnage", "--tml", "24.5", "--tmb", "5.5", "--tmd", "2.8"],
        None,
    );
    assert!(!too_long.status.success());
    assert!(String::from_utf8_lossy(&too_long.stderr).contains("24 m"));

    let overhang_ok = run_cmd(
        &[
            "tonnage",
            "--tml",
            "24.5",
            "--tmb",
            "5.5",
            "--tmd",
            "2.8",
            "--reg-length",
            "23.4",
        ],
        None,
    );
    assert!(
        overhang_ok.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&overhang_ok.stderr)
    );
}

#[test]
fn stream_hull_scale_catamaran_mount_table_plot_pick_and_cache() {
    let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.2, 1.0).unwrap();
    let iges_text = hullgeom::iges::write(&surfaces, "wigley").unwrap();
    let pid = std::process::id();
    let dir = std::env::temp_dir().join(format!("boatmath_tonnage_test_{pid}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let igs_path = dir.join("wigley.igs");
    let cache_dir = dir.join("cache");
    std::fs::write(&igs_path, iges_text).unwrap();

    // 1. hull -> scale (0.8, 1.0, 1.2) -> tonnage (with exempt small cabin) -> table & plot
    let hull_out = run_cmd(
        &["hull", igs_path.to_str().unwrap(), "--waterline", "-0.4"],
        None,
    );
    assert!(hull_out.status.success());

    let scale_out = run_cmd(&["scale", "--by", "0.8,1.0,1.2"], Some(&hull_out.stdout));
    assert!(scale_out.status.success());

    // Note: original hull (10 m) + 3 scaled hulls (8 m, 10 m, 12 m) = 4 hulls in stream.
    let ton_out = run_cmd(
        &[
            "--cache",
            cache_dir.to_str().unwrap(),
            "tonnage",
            "--kind",
            "sail",
            "--space",
            "4x0.9x1.4",
        ],
        Some(&scale_out.stdout),
    );
    assert!(
        ton_out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ton_out.stderr)
    );

    // Second run hits the cache and produces identical records.
    let ton_cached = run_cmd(
        &[
            "--cache",
            cache_dir.to_str().unwrap(),
            "tonnage",
            "--kind",
            "sail",
            "--space",
            "4x0.9x1.4",
        ],
        Some(&scale_out.stdout),
    );
    assert!(ton_cached.status.success());
    assert_eq!(
        parse_jsonl(&ton_out.stdout),
        parse_jsonl(&ton_cached.stdout)
    );

    let table_out = run_cmd(
        &[
            "table",
            "hull.name",
            "tml",
            "tmb",
            "tmd",
            "gt",
            "nt",
            "form",
        ],
        Some(&ton_out.stdout),
    );
    assert!(table_out.status.success());
    let table_text = String::from_utf8(table_out.stdout).unwrap();
    assert!(table_text.contains("hull.name\ttml\ttmb\ttmd\tgt\tnt\tform"));
    assert!(table_text.contains("4A"));

    let plot_out = run_cmd(
        &["plot", "-x", "tml", "-y", "gt", "-y", "nt"],
        Some(&ton_out.stdout),
    );
    assert!(plot_out.status.success());
    assert!(String::from_utf8_lossy(&plot_out.stdout).contains("<svg"));

    // Pick keeps tonnage and its hull + parent hull ancestry.
    let pick_out = run_cmd(
        &["pick", "tonnage", "--where", "tml=12"],
        Some(&ton_out.stdout),
    );
    assert!(pick_out.status.success());
    let picked = parse_jsonl(&pick_out.stdout);
    assert!(picked
        .iter()
        .any(|r| r["type"] == "tonnage" && r["tml"] == 12.0));
    assert!(picked.iter().any(|r| r["type"] == "hull"));

    // 2. hull -> case --span 3.0 -> mount -> tonnage (catamaran ignores mount, counts 2 demihulls)
    let case_out = run_cmd(&["case", "--span", "3.0"], Some(&hull_out.stdout));
    assert!(case_out.status.success());
    let mount_out = run_cmd(
        &["mount", "--stock", "oceanvolt-sd8", "--x", "2.0"],
        Some(&case_out.stdout),
    );
    assert!(mount_out.status.success());
    // Note: mount passes bare case + mounted case (2 cases, both span 3.0).
    let cat_ton = run_cmd(
        &["tonnage", "--bridge", "5x2.2x0.8"],
        Some(&mount_out.stdout),
    );
    assert!(
        cat_ton.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&cat_ton.stderr)
    );
    let cat_recs: Vec<Value> = parse_jsonl(&cat_ton.stdout)
        .into_iter()
        .filter(|r| r["type"] == "tonnage")
        .collect();
    assert_eq!(cat_recs.len(), 2);
    // Bare catamaran case and mounted catamaran case MUST have identical GT & NT!
    assert_eq!(cat_recs[0]["hull_count"], 2);
    assert_eq!(cat_recs[0]["gt"], cat_recs[1]["gt"]);
    assert_eq!(cat_recs[0]["nt"], cat_recs[1]["nt"]);

    // 3. Partial-stream failure when one scaled hull exceeds 24 m statutory length:
    // 10 m * 2.6 = 26 m -> reg_length = 0.96 * 26 = 24.96 m >= 24 m.
    let mixed_scale = run_cmd(&["scale", "--by", "2.6"], Some(&hull_out.stdout));
    assert!(mixed_scale.status.success());
    let mixed_ton = run_cmd(&["tonnage"], Some(&mixed_scale.stdout));
    assert!(!mixed_ton.status.success());
    let mixed_recs: Vec<Value> = parse_jsonl(&mixed_ton.stdout)
        .into_iter()
        .filter(|r| r["type"] == "tonnage")
        .collect();
    assert_eq!(mixed_recs.len(), 2);
    assert!(mixed_recs[0]["gt"].as_f64().is_some());
    assert!(mixed_recs[1]["error"].as_str().unwrap().contains("24 m"));

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn camber_and_stl_hull_tonnage() {
    // Camber v2 JSON document (cruiser.json is a 0.97 m model in mm; scale ×10 to 9.70 m)
    let cruiser = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../boatmath/tests/camber/cruiser.json");
    let hull_out = run_cmd(&["hull", cruiser.to_str().unwrap()], None);
    assert!(
        hull_out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&hull_out.stderr)
    );
    let scaled_out = run_cmd(&["scale", "--by", "10"], Some(&hull_out.stdout));
    assert!(scaled_out.status.success());
    let ton_out = run_cmd(&["tonnage", "--kind", "sail"], Some(&scaled_out.stdout));
    assert!(
        ton_out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ton_out.stderr)
    );
    let ton: Vec<Value> = parse_jsonl(&ton_out.stdout)
        .into_iter()
        .filter(|r| r["type"] == "tonnage")
        .collect();
    assert_eq!(ton.len(), 2);
    let l0 = ton[0]["tml"].as_f64().unwrap();
    let l1 = ton[1]["tml"].as_f64().unwrap();
    assert!((0.95..1.30).contains(&l0), "l0 = {l0}");
    assert!((l1 - 10.0 * l0).abs() <= 0.05, "l0 = {l0}, l1 = {l1}");
    assert!(ton[1]["gt"].as_f64().unwrap() > 1.0);
    assert_eq!(ton[1]["form"], "4A");

    // ASCII STL box barge in millimetres (10000 x 3000 x 1500 mm, waterline at z = 0)
    let (x0, x1, y0, y1, z0, z1) = (0.0, 10000.0, -1500.0, 1500.0, -750.0, 750.0);
    let quads = [
        ([x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0]), // bottom
        ([x0, y0, z1], [x0, y1, z1], [x1, y1, z1], [x1, y0, z1]), // top
        ([x0, y0, z0], [x0, y0, z1], [x1, y0, z1], [x1, y0, z0]), // stbd
        ([x0, y1, z0], [x1, y1, z0], [x1, y1, z1], [x0, y1, z1]), // port
        ([x0, y0, z0], [x0, y1, z0], [x0, y1, z1], [x0, y0, z1]), // aft
        ([x1, y0, z0], [x1, y0, z1], [x1, y1, z1], [x1, y1, z0]), // fore
    ];
    let mut stl = String::from("solid barge\n");
    for (a, b, c, d) in quads {
        for tri in [[a, b, c], [a, c, d]] {
            stl.push_str("  facet normal 0 0 0\n    outer loop\n");
            for v in tri {
                stl.push_str(&format!("      vertex {} {} {}\n", v[0], v[1], v[2]));
            }
            stl.push_str("    endloop\n  endfacet\n");
        }
    }
    stl.push_str("endsolid barge\n");
    let stl_path = std::env::temp_dir().join(format!("barge_tonnage_{}.stl", std::process::id()));
    std::fs::write(&stl_path, stl).unwrap();

    let stl_hull = run_cmd(&["hull", stl_path.to_str().unwrap(), "--units", "mm"], None);
    assert!(
        stl_hull.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&stl_hull.stderr)
    );
    let stl_ton = run_cmd(
        &["tonnage", "--kind", "barge", "--propulsion", "power"],
        Some(&stl_hull.stdout),
    );
    assert!(
        stl_ton.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&stl_ton.stderr)
    );
    let rec = parse_jsonl(&stl_ton.stdout)
        .into_iter()
        .find(|r| r["type"] == "tonnage")
        .unwrap();
    // 10.00 x 3.00 x 1.50 * 0.20 = 9.00 GT, 9.00 * 0.75 = 6.75 NT
    assert_eq!(rec["tml"], 10.0);
    assert_eq!(rec["tmb"], 3.0);
    assert_eq!(rec["tmd"], 1.5);
    assert_eq!(rec["gt"], 9.0);
    assert_eq!(rec["nt"], 6.75);

    let _ = std::fs::remove_file(&stl_path);
}
