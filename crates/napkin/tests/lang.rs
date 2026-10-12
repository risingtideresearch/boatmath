//! End-to-end tests: run napkin source through the binary and check its output.

fn run(src: &str) -> (String, String) {
    let mut it = napkin::Interp::new(20_000, 7);
    std::env::set_current_dir(env!("CARGO_MANIFEST_DIR")).unwrap();
    let (out, r) = napkin::run_captured(&mut it, src);
    (out, r.err().unwrap_or_default())
}

fn ok(src: &str) -> String {
    let (out, err) = run(src);
    assert!(err.is_empty(), "unexpected error: {}", err);
    out
}

/// The leading number after `prefix` in the output, e.g. "P(x > 1 m) = " -> 0.75.
fn num_after(out: &str, prefix: &str) -> f64 {
    let rest = &out[out
        .find(prefix)
        .unwrap_or_else(|| panic!("{:?} not in output:\n{}", prefix, out))
        + prefix.len()..];
    let tok: String = rest
        .chars()
        .take_while(|c| c.is_ascii_digit() || "-.e".contains(*c))
        .collect();
    tok.parse().unwrap()
}

fn close(out: &str, prefix: &str, want: f64, tol: f64) {
    let got = num_after(out, prefix);
    assert!(
        (got - want).abs() <= tol,
        "{}: got {}, want {} ± {}",
        prefix,
        got,
        want,
        tol
    );
}

fn has(out: &str, want: &str) {
    assert!(
        out.contains(want),
        "expected {:?} in output:\n{}",
        want,
        out
    );
}

#[test]
fn units_and_conversions() {
    let out = ok("1 t + 500 kg\n12 kn -> m/s\n2 m * 3 m\n10 mm * 2 m^2 * 7.85 t/m^3\n30 deg -> rad\n1 LT -> kg");
    has(&out, "1 t + 500 kg = 1.5 t");
    has(&out, "= 6.173 m/s");
    has(&out, "2 m * 3 m = 6 m^2");
    has(&out, "= 0.157 t");
    has(&out, "= 0.5236 rad");
    has(&out, "= 1016 kg");
}

#[test]
fn dimension_errors_are_reported() {
    let (_, err) = run("let a = 2 m\nlet b = 3 kg\na + b");
    has(&err, "line 3");
    has(&err, "cannot add length and mass");
    let (_, err) = run("5 m -> kg");
    has(&err, "cannot convert length to kg");
}

#[test]
fn shared_samples_preserve_correlation() {
    let out = ok("let h = 12 ± 0.5 m\nh - h\nh / h");
    has(&out, "h - h = 0 m");
    has(&out, "h / h = 1");
}

#[test]
fn distributions_have_expected_moments() {
    let out = ok("mean(10 ± 2 m)\nsd(10 ± 2 m)\nquantile(2 to 8 t, 0.05)\nquantile(2 to 8 t, 0.95)\nmean(200 t ± 5%)");
    close(&out, "mean(10 ± 2 m) = ", 10.0, 0.05);
    close(&out, "sd(10 ± 2 m) = ", 2.0, 0.05);
    close(&out, "quantile(2 to 8 t, 0.05) = ", 2.0, 0.1);
    close(&out, "quantile(2 to 8 t, 0.95) = ", 8.0, 0.3);
    close(&out, "mean(200 t ± 5%) = ", 200.0, 0.5);
}

#[test]
fn uncertain_if_selects_per_sample() {
    let out = ok("let x = uniform(0 m, 1 m)\nP(x > 0.25 m)\nmean(if x > 0.5 m then 1 t else 0 t)");
    close(&out, "P(x > 0.25 m) = ", 0.75, 0.01);
    close(&out, "then 1 t else 0 t) = ", 0.5, 0.01);
}

#[test]
fn section_properties_of_a_triangle() {
    let out = ok("section([[0, 0], [4, 0], [0, 3]] m)");
    has(&out, "area     : 6 m^2");
    has(&out, "centroid : [1.333 m, 1 m]");
    has(&out, "ixx      : 3 m^4");
    has(&out, "iyy      : 5.333 m^4");
    has(&out, "ixy      : -2 m^4");
}

#[test]
fn planar_polygon_in_3d() {
    let out = ok(
        "area([[0,0,0],[2,0,0],[2,3,0],[0,3,0]] m)\ncentroid([[0,0,0],[2,0,1],[2,3,1],[0,3,0]] m)",
    );
    has(&out, "= 6 m^2");
    has(&out, "= [1 m, 1.5 m, 0.5 m]");
}

#[test]
fn combined_cg_of_point_masses() {
    let out = ok("let items = [{mass: 2 t, cg: [0, 0, 0] m}, {mass: 1 t, cg: [3, 0, 6] m}]\ncombine(items).cg\ninertia(items).ixx");
    has(&out, "combine(items).cg = [1 m, 0 m, 2 m]");
    // Ixx about CG: 2·(0+4) + 1·(0+16) = 24 t·m²
    has(&out, "inertia(items).ixx = 24 t·m^2");
}

#[test]
fn sensitivity_names_inputs_by_where_they_are_written() {
    let out = ok("let items = [{name: \"A\", mass: 10 t ± 1 t, cg: [0, 0, 1] m}, {name: \"B\", mass: 1 t, cg: [0, 0, 5 ± 0.01] m}]\nsensitivity(combine(items).cg[2])");
    has(&out, "↓ A.mass");
}
