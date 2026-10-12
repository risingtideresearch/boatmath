//! napkin-boat against hulls with known answers: a Wigley hull (exact patches)
//! for the float, and a box barge (an STL written here) for GZ, whose
//! wall-sided curve is closed form. The CAD hulls the README uses are not in
//! the repository; the test on one skips when it is absent.

fn run(src: &str) -> String {
    let mut it = napkin_boat::interp(4000, 3);
    let (out, r) = napkin::run_captured(&mut it, src);
    if let Err(e) = r {
        panic!("error: {}\noutput:\n{}", e, out);
    }
    out
}

fn num_after(out: &str, prefix: &str) -> f64 {
    let at = out
        .find(prefix)
        .unwrap_or_else(|| panic!("{:?} not in output:\n{}", prefix, out));
    let tok: String = out[at + prefix.len()..]
        .chars()
        .take_while(|c| c.is_ascii_digit() || "-.e".contains(*c))
        .collect();
    tok.parse()
        .unwrap_or_else(|_| panic!("no number after {:?} in:\n{}", prefix, out))
}

fn close(out: &str, prefix: &str, want: f64, tol: f64) {
    let got = num_after(out, prefix);
    assert!(
        (got - want).abs() <= tol,
        "{}: got {}, want {} ± {}\n{}",
        prefix,
        got,
        want,
        tol,
        out
    );
}

/// A closed box, `l × b × d`, keel at z = 0, as ASCII STL.
fn box_stl(l: f64, b: f64, d: f64) -> String {
    let h = b / 2.0;
    let v = |i: usize| {
        [
            if i & 1 == 0 { 0.0 } else { l },
            if i & 2 == 0 { -h } else { h },
            if i & 4 == 0 { 0.0 } else { d },
        ]
    };
    // Each face as two triangles over corner indices, wound outward.
    let quads = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut s = String::from("solid box\n");
    for q in quads {
        for t in [[q[0], q[1], q[2]], [q[0], q[2], q[3]]] {
            s += "facet normal 0 0 0\nouter loop\n";
            for i in t {
                let p = v(i);
                s += &format!("vertex {} {} {}\n", p[0], p[1], p[2]);
            }
            s += "endloop\nendfacet\n";
        }
    }
    s + "endsolid box\n"
}

fn box_file() -> String {
    // One file per call: tests run in parallel.
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let k = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("napkin-boat-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("box{k}.stl"));
    std::fs::write(&path, box_stl(10.0, 2.0, 1.5)).unwrap();
    path.to_string_lossy().into_owned()
}

#[test]
fn wigley_floats_at_its_design_waterline() {
    // V = 4/9 · L B T for the Wigley form.
    let out = run("let h = wigley(10 m, 1 m, 0.625 m)\n\
         let d = hydrostatics(h)\n\
         d.volume\n\
         let f = float(h, {mass: d.displacement, cg: [0, 0, 0.2] m})\n\
         f.sinkage -> mm\n\
         f.trim\n");
    close(&out, "d.volume = ", 4.0 / 9.0 * 10.0 * 0.625, 0.003);
    close(&out, "f.sinkage -> mm = ", 0.0, 0.01);
    close(&out, "f.trim = ", 0.0, 1e-4);
}

#[test]
fn uncertain_load_is_interpolated_accurately() {
    // A lighter load (the Wigley surfaces end at the waterline, so it can't
    // sink): the surrogate's median must agree with a direct solve at the
    // median load to well under a mm.
    let out = run("let h = wigley(10 m, 1 m, 0.625 m)\n\
         let m = 0.8 * hydrostatics(h).displacement\n\
         let direct = float(h, {mass: m, cg: [0.05, 0, 0] m})\n\
         let f = float(h, {mass: m ± 3%, cg: [0.05 ± 0.02, 0, 0] m})\n\
         direct.sinkage -> mm\n\
         median(f.sinkage) -> mm\n\
         direct.trim\n\
         median(f.trim)\n");
    assert!(!out.contains("warning"), "unexpected warning:\n{}", out);
    let direct = num_after(&out, "direct.sinkage -> mm = ");
    close(&out, "median(f.sinkage) -> mm = ", direct, 0.5);
    let trim = num_after(&out, "direct.trim = ");
    close(&out, "median(f.trim) = ", trim, 0.01);
}

#[test]
fn box_barge_gz_is_wall_sided() {
    // L 10, B 2, D 1.5, T 0.5 (design waterline z = 0.5), KG 0.8:
    // KB = 0.25, BM = B²/12T = 2/3, GM = 0.1167; wall-sided until the deck
    // edge goes under at 45°: GZ = sin φ (GM + ½ BM tan² φ). hullgeom's ray
    // sectioning is good to about 1% in volume and 5 mm in GM on this box's
    // sharp corners; a frame or VCG mistake here would be tenths of a metre.
    let path = box_file();
    let src = format!(
        "let b = hull(\"{path}\", {{units: 1 m, waterline: 0.5 m}})\n\
         let d = hydrostatics(b)\n\
         let s = stability(b, {{mass: d.displacement, cg: [5, 0, 0.8] m}}, {{max_heel: 60 deg}})\n\
         d.volume\n\
         s.gm\n\
         s.gz[10]\n\
         s.gz[15]\n"
    );
    let out = run(&src);
    close(&out, "d.volume = ", 10.0, 0.1);
    let (gm, bm) = (0.25 + 2.0 / 3.0 - 0.8, 2.0 / 3.0);
    let wall = |deg: f64| {
        let p = deg.to_radians();
        p.sin() * (gm + 0.5 * bm * p.tan().powi(2))
    };
    close(&out, "s.gm = ", gm, 0.005);
    close(&out, "s.gz[10] = ", wall(20.0), 0.005);
    close(&out, "s.gz[15] = ", wall(30.0), 0.008);
}

#[test]
fn an_off_centre_cg_lists_the_box() {
    // Small angles: tan(list) ≈ TCG / GM.
    let path = box_file();
    let src = format!(
        "let b = hull(\"{path}\", {{units: 1 m, waterline: 0.5 m}})\n\
         let s = stability(b, {{mass: hydrostatics(b).displacement, cg: [5, 0.01, 0.8] m}}, {{max_heel: 60 deg}})\n\
         s.list\n"
    );
    let out = run(&src);
    let gm: f64 = 0.25 + 2.0 / 3.0 - 0.8;
    close(&out, "s.list = ", (0.01 / gm).atan().to_degrees(), 0.3);
}

#[test]
fn gz_vcg_is_interpolated_accurately() {
    // An uncertain KG spreads GM by exactly its own spread.
    let path = box_file();
    let src = format!(
        "let b = hull(\"{path}\", {{units: 1 m, waterline: 0.5 m}})\n\
         let s = stability(b, {{mass: hydrostatics(b).displacement, cg: [5, 0, 0.8 ± 0.01] m}}, {{max_heel: 40 deg}})\n\
         sd(s.gm) -> mm\n\
         median(s.gm) -> mm\n"
    );
    let out = run(&src);
    assert!(!out.contains("warning"), "unexpected warning:\n{}", out);
    close(&out, "sd(s.gm) -> mm = ", 10.0, 0.5);
    close(
        &out,
        "median(s.gm) -> mm = ",
        1000.0 * (0.25 + 2.0 / 3.0 - 0.8),
        5.0,
    );
}

/// The Graf catamaran agrees with `boatmath statics` (case: span 1.6 m,
/// 1540 kg at LCG 2.66 m, VCG 0.30 m): trim −0.0211°, GM 2.4356 m,
/// max GZ 0.6893 m, vanishing 75.39°.
#[test]
fn graf_catamaran_matches_boatmath_statics() {
    let root = std::env::var("BOATMATH_CAD_DIR")
        .unwrap_or_else(|_| format!("{}/../..", env!("CARGO_MANIFEST_DIR")));
    let path = format!("{}/graf-demihull.igs", root);
    if !std::path::Path::new(&path).exists() {
        eprintln!("skipping: {} is not present (not in the repository)", path);
        return;
    }
    let out = run(&format!(
        "let cat = hull(\"{path}\", {{span: 1.6 m}})\n\
         let load = {{mass: 1540 kg, cg: [2.66, 0, 0.30] m}}\n\
         float(cat, load).trim\n\
         let s = stability(cat, load)\n\
         s.gm -> mm\n\
         s.max_gz -> mm\n\
         s.vanishing\n"
    ));
    close(&out, "float(cat, load).trim = ", -0.02112, 2e-4);
    close(&out, "s.gm -> mm = ", 2435.6, 1.0);
    close(&out, "s.max_gz -> mm = ", 689.3, 0.5);
    close(&out, "s.vanishing = ", 75.39, 0.01);
}

#[test]
fn box_shell_and_slices_are_exact() {
    // A closed 10 × 2 × 1.5 m box: shell 2(LB + LD + BD) = 76 m²; below
    // z = 0.5 the bottom, sides and ends make 20 + 10 + 2 = 32 m².
    let path = box_file();
    let out = run(&format!(
        "let b = hull(\"{path}\", {{units: 1 m, waterline: 0.5 m}})\n\
         shell(b).area\n\
         shell(b, {{z: [0 m, 0.5 m]}}).area\n\
         let s = slice(b, 5 m)\n\
         s.area\n\
         s.girth\n\
         let w = slice(b, 5 m, {{below: 0.5 m}})\n\
         w.area\n\
         w.girth\n\
         w.centroid[2]\n\
         let u = shell(b, {{z: [0 m, 0.5 ± 0.05 m]}})\n\
         mean(u.area)\n\
         sd(u.area)\n"
    ));
    close(&out, "shell(b).area = ", 76.0, 1e-9);
    close(&out, "shell(b, {z: [0 m, 0.5 m]}).area = ", 32.0, 1e-9);
    close(&out, "s.area = ", 3.0, 1e-9);
    close(&out, "s.girth = ", 7.0, 1e-9);
    close(&out, "w.area = ", 1.0, 1e-9);
    close(&out, "w.girth = ", 3.0, 1e-9);
    close(&out, "w.centroid[2] = ", 0.25, 1e-9);
    // The wetted band grows by the girth of sides and ends, 24 m, per metre.
    close(&out, "mean(u.area) = ", 32.0, 0.05);
    close(&out, "sd(u.area) = ", 24.0 * 0.05, 0.03);
}

#[test]
fn wigley_slices_match_the_form() {
    // y = B/2 (1 − (2x/L)²)(1 − (z/T)²): a section's area is ⅔ B T (1 − (2x/L)²),
    // its centroid ⅜ T below the waterline.
    let out = run("let w = wigley(10 m, 1 m, 0.625 m)\n\
         slice(w, 0 m).area -> cm^2\n\
         slice(w, 2.5 m).area -> cm^2\n\
         slice(w, 0 m).centroid[2] -> mm\n\
         shell(w).area\n\
         hydrostatics(w).wetted\n");
    close(
        &out,
        "slice(w, 0 m).area -> cm^2 = ",
        1e4 * 2.0 / 3.0 * 0.625,
        2.0,
    );
    close(
        &out,
        "slice(w, 2.5 m).area -> cm^2 = ",
        1e4 * 2.0 / 3.0 * 0.625 * 0.75,
        2.0,
    );
    close(
        &out,
        "slice(w, 0 m).centroid[2] -> mm = ",
        -1000.0 * 0.375 * 0.625,
        0.5,
    );
    let wetted = num_after(&out, "hydrostatics(w).wetted = ");
    close(&out, "shell(w).area = ", wetted, 0.01 * wetted);
}

/// Slices integrated along the Graf demihull give its displaced volume, and
/// its shell below the waterline its wetted surface, to the mesh's ~0.4%.
#[test]
fn graf_slices_and_shell_agree_with_hullgeom() {
    let root = std::env::var("BOATMATH_CAD_DIR")
        .unwrap_or_else(|_| format!("{}/../..", env!("CARGO_MANIFEST_DIR")));
    let path = format!("{}/graf-demihull.igs", root);
    if !std::path::Path::new(&path).exists() {
        eprintln!("skipping: {} is not present (not in the repository)", path);
        return;
    }
    let out = run(&format!(
        "let h = hull(\"{path}\")\n\
         let d = hydrostatics(h)\n\
         let xs = linspace(0 m, 6.49 m, 131)\n\
         trapz(xs, map(xs, |x| slice(h, x, {{below: d.waterline}}).area)) -> m^3\n\
         d.volume\n\
         shell(h, {{z: [-1 m, d.waterline]}}).area\n\
         d.wetted\n"
    ));
    let v = num_after(&out, "d.volume = ");
    close(&out, "}).area)) -> m^3 = ", v, 0.006 * v);
    let s = num_after(&out, "d.wetted = ");
    close(
        &out,
        "shell(h, {z: [-1 m, d.waterline]}).area = ",
        s,
        0.006 * s,
    );
}
