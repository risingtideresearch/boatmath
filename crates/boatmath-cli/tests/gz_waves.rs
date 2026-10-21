use std::io::Write;
use std::process::{Command, Stdio};

fn boatmath_bin() -> &'static str {
    env!("CARGO_BIN_EXE_boatmath")
}

#[test]
fn test_cli_gz_waves_table_and_plot() {
    let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 1.0).unwrap();
    let iges_text = hullgeom::iges::write(&surfaces, "wigley").unwrap();
    let dir = std::env::temp_dir();
    let path = dir.join("test_cli_wigley.igs");
    std::fs::write(&path, iges_text).unwrap();

    // 1. boatmath hull <path> --waterline -0.375
    let hull_out = Command::new(boatmath_bin())
        .arg("hull")
        .arg(&path)
        .arg("--waterline")
        .arg("-0.375")
        .output()
        .expect("boatmath hull");
    assert!(
        hull_out.status.success(),
        "hull stderr: {}",
        String::from_utf8_lossy(&hull_out.stderr)
    );

    // 2. boatmath case --vcg -0.3
    let mut case_child = Command::new(boatmath_bin())
        .arg("case")
        .arg("--vcg")
        .arg("-0.3")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("boatmath case spawn");
    case_child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(&hull_out.stdout)
        .unwrap();
    let case_out = case_child.wait_with_output().unwrap();
    assert!(
        case_out.status.success(),
        "case stderr: {}",
        String::from_utf8_lossy(&case_out.stderr)
    );

    // 3. boatmath gz-waves --length 10m --height 0.2m
    let mut gz_child = Command::new(boatmath_bin())
        .arg("gz-waves")
        .arg("--length")
        .arg("10m")
        .arg("--height")
        .arg("0.2m")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("boatmath gz-waves spawn");
    gz_child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(&case_out.stdout)
        .unwrap();
    let gz_out = gz_child.wait_with_output().unwrap();
    assert!(
        gz_out.status.success(),
        "gz-waves stderr: {}",
        String::from_utf8_lossy(&gz_out.stderr)
    );
    let gz_stream = String::from_utf8(gz_out.stdout).unwrap();
    assert!(
        gz_stream.contains("\"type\":\"gz_waves\"") || gz_stream.contains("\"type\": \"gz_waves\"")
    );

    // 4. boatmath table case.params.mass length height calm.max_gz
    let mut table_child = Command::new(boatmath_bin())
        .arg("table")
        .args(["case.params.mass", "length", "height", "calm.max_gz"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("boatmath table spawn");
    table_child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(gz_stream.as_bytes())
        .unwrap();
    let table_out = table_child.wait_with_output().unwrap();
    assert!(
        table_out.status.success(),
        "table stderr: {}",
        String::from_utf8_lossy(&table_out.stderr)
    );
    let table_str = String::from_utf8(table_out.stdout).unwrap();
    assert!(table_str.contains("mass"));
    assert!(table_str.contains("length"));
    assert!(table_str.contains("height"));
    assert!(table_str.contains("max_gz"));

    // 5. boatmath plot -x heading_deg -y worst.max_gz --explode headings
    let mut plot_child = Command::new(boatmath_bin())
        .arg("plot")
        .args([
            "-x",
            "heading_deg",
            "-y",
            "worst.max_gz",
            "--explode",
            "headings",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("boatmath plot spawn");
    plot_child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(gz_stream.as_bytes())
        .unwrap();
    let plot_out = plot_child.wait_with_output().unwrap();
    assert!(
        plot_out.status.success(),
        "plot stderr: {}",
        String::from_utf8_lossy(&plot_out.stderr)
    );
    let plot_str = String::from_utf8(plot_out.stdout).unwrap();
    assert!(plot_str.contains("<svg"));

    let _ = std::fs::remove_file(path);
}
