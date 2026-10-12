//! Hulls from CAD, floated and heeled under napkin's uncertain loads.
//!
//! Frames: an item's `cg` is in the hull file's frame, in metres: `x` as in
//! the file (bow toward +x), `z` up, and `y` from the platform's centreline
//! (a monohull's centreplane, or midway between a catamaran's demihulls).
//! hullgeom measures heights from the design waterline; the conversion is
//! `vcg = z − waterline`.

use std::rc::Rc;

use hullgeom::float::{solve_equilibrium_sectional, LoadCase};
use hullgeom::iges::{HullPose, Platform, SectionalOptions};
use hullgeom::source::{HullSource, SourceHull};
use hullgeom::stability::{gz_curve, FullSection, GzOptions, StabilityHull, StabilityLoad};
use hullgeom::{Fluid, SectionalHull};
use napkin::builtins::{arity, item_parts, num, q_unit};
use napkin::ops::R;
use napkin::value::Native;
use napkin::{Dim, Interp, Quantity, Samples, UnitExpr, Value};

use crate::geometry::{self, Bounds, Mesh};
use crate::surrogate;

/// A platform of hulls cut from one file: the file's hulls where they lie, or
/// its one hull doubled into a catamaran.
pub struct BoatHull {
    name: String,
    source: Box<dyn HullSource>,
    /// CAD height of the design waterline [m].
    waterline_z: f64,
    opts: SectionalOptions,
    layout: Vec<(usize, HullPose)>,
    /// Each layout entry cut at the design waterline.
    design: Vec<(SectionalHull, hullgeom::Placement)>,
    /// Whole sections, keel to sheer, for GZ.
    whole: Vec<StabilityHull>,
    /// Each layout entry as triangles in the file frame, for measuring.
    meshes: Vec<Mesh>,
}

impl Native for BoatHull {
    fn type_name(&self) -> &'static str {
        "hull"
    }
    fn describe(&self) -> String {
        format!(
            "<hull {}: {} hull{}, {:.3} m long, Δ {:.4} t at its design waterline>",
            self.name,
            self.design.len(),
            if self.design.len() == 1 { "" } else { "s" },
            self.length(),
            self.volume() * Fluid::SEAWATER_15C.density / 1000.0
        )
    }
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

impl BoatHull {
    fn volume(&self) -> f64 {
        self.design.iter().map(|(h, _)| h.displaced_volume()).sum()
    }

    fn length(&self) -> f64 {
        self.design
            .iter()
            .map(|(h, _)| h.length())
            .fold(0.0, f64::max)
    }

    fn sources(&self) -> Vec<SourceHull<'_>> {
        self.layout
            .iter()
            .map(|&(index, pose)| SourceHull {
                source: self.source.as_ref(),
                index,
                waterline_z: self.waterline_z,
                pose,
            })
            .collect()
    }
}

// ---------------------------------------------------------------- options

fn opt<'a>(opts: Option<&'a Value>, name: &str) -> Option<&'a Value> {
    opts.and_then(|o| o.field(name))
}

fn opt_len(opts: Option<&Value>, name: &str) -> R<Option<f64>> {
    match opt(opts, name) {
        None => Ok(None),
        Some(v) => {
            let q = num(v)?;
            if q.dim != Dim::length() {
                return Err(format!(
                    "{} must be a length, got {}",
                    name,
                    q.dim.describe()
                ));
            }
            q.s.exact()
                .map(Some)
                .ok_or_else(|| format!("{} must be exact", name))
        }
    }
}

fn opt_count(opts: Option<&Value>, name: &str) -> R<Option<usize>> {
    match opt(opts, name) {
        None => Ok(None),
        Some(v) => {
            let q = num(v)?;
            match q.s.exact() {
                Some(x) if q.dim.is_none() && x >= 1.0 => Ok(Some(x as usize)),
                _ => Err(format!("{} must be an exact count", name)),
            }
        }
    }
}

fn check_opts(opts: Option<&Value>, allowed: &[&str]) -> R<()> {
    if let Some(v) = opts {
        let Value::Record(fs) = v else {
            return Err(format!("options must be a record, got {}", v.type_name()));
        };
        for (n, _) in fs.iter() {
            if !allowed.contains(&n.as_str()) {
                return Err(format!(
                    "unknown option '{}' (expected one of: {})",
                    n,
                    allowed.join(", ")
                ));
            }
        }
    }
    Ok(())
}

fn hull_arg(v: &Value) -> R<Rc<dyn Native>> {
    match v {
        Value::Native(n) if n.as_any().is::<BoatHull>() => Ok(n.clone()),
        _ => Err(format!("expected a hull, got {}", v.type_name())),
    }
}

fn density(opts: Option<&Value>) -> R<f64> {
    match opt(opts, "rho") {
        None => Ok(Fluid::SEAWATER_15C.density),
        Some(v) => {
            let q = num(v)?;
            if q.dim != Dim::mass().div(&Dim::length().pow(napkin::Rational::int(3))) {
                return Err("rho must be a density".into());
            }
            q.s.exact().ok_or_else(|| "rho must be exact".to_string())
        }
    }
}

// ---------------------------------------------------------------- loading

const HULL_OPTS: &[&str] = &[
    "waterline",
    "span",
    "units",
    "stations",
    "rays",
    "centerplane",
];

fn build(name: &str, source: Box<dyn HullSource>, opts: Option<&Value>) -> R<BoatHull> {
    let d = SectionalOptions::default();
    let waterline_z = opt_len(opts, "waterline")?.unwrap_or(0.0);
    let sopts = SectionalOptions {
        waterline_z,
        centerplane: opt_len(opts, "centerplane")?,
        stations: opt_count(opts, "stations")?.unwrap_or(d.stations),
        rays: opt_count(opts, "rays")?.unwrap_or(d.rays),
        ..d
    };
    // The file's hulls that are wet at the design waterline, as boatmath cuts them.
    let mut cut = vec![];
    for i in 0..source.len() {
        if let Some(h) = source
            .situate_sectional(
                i,
                waterline_z,
                &HullPose::default(),
                &Platform::default(),
                &sopts,
            )
            .map_err(|e| format!("hull {}: {}", i + 1, e))?
        {
            cut.push((i, h));
        }
    }
    if cut.is_empty() {
        return Err(format!(
            "{}: no hull is wet at waterline z = {} m",
            name, waterline_z
        ));
    }
    let layout: Vec<(usize, HullPose)> = match opt_len(opts, "span")? {
        None => cut.iter().map(|(i, _)| (*i, HullPose::default())).collect(),
        Some(span) => {
            if cut.len() != 1 {
                return Err(format!(
                    "a catamaran doubles a single hull; this file holds {}",
                    cut.len()
                ));
            }
            let (i, h) = &cut[0];
            let beam = 2.0
                * h.hull
                    .curves()
                    .flat_map(|(_, c)| c.iter().map(|p| p.0))
                    .fold(0.0, f64::max);
            if span <= beam {
                return Err(format!(
                    "span {} m: the demihulls overlap (beam {:.3} m)",
                    span, beam
                ));
            }
            let yc = h.placement.y;
            [0.5 * span, -0.5 * span]
                .map(|y| {
                    (
                        *i,
                        HullPose {
                            dy: y - yc,
                            ..HullPose::default()
                        },
                    )
                })
                .to_vec()
        }
    };
    let design = layout
        .iter()
        .map(|(i, pose)| {
            let h = &cut.iter().find(|(j, _)| j == i).unwrap().1;
            (
                h.hull.clone(),
                hullgeom::Placement {
                    x: h.placement.x,
                    y: h.placement.y + pose.dy,
                },
            )
        })
        .collect::<Vec<_>>();
    let l_ref = design.iter().map(|(h, _)| h.length()).fold(0.0, f64::max);
    let whole = whole_sections(source.as_ref(), waterline_z, &layout, &sopts, l_ref)?;
    let meshes = layout
        .iter()
        .zip(&design)
        .map(|((i, pose), (_, place))| {
            let (v, t) = source
                .posed_tessellation(*i, waterline_z, pose, &Platform::default())
                .map_err(|e| e.to_string())?;
            // Posed heights are from the waterline; the file frame's are not.
            let v = v
                .into_iter()
                .map(|p| [p[0], p[1], p[2] + waterline_z])
                .collect();
            Ok(Mesh::new(v, t, place.y))
        })
        .collect::<R<Vec<_>>>()?;
    Ok(BoatHull {
        name: name.to_string(),
        source,
        waterline_z,
        opts: sopts,
        layout,
        design,
        whole,
        meshes,
    })
}

/// Each hull's whole sections, keel to sheer, at the design pose: the hull cut
/// with the water raised past its top (as boatmath's statics does).
fn whole_sections(
    src: &dyn HullSource,
    wl: f64,
    layout: &[(usize, HullPose)],
    opts: &SectionalOptions,
    l_ref: f64,
) -> R<Vec<StabilityHull>> {
    let mut out = vec![];
    for (i, pose) in layout {
        let (v, _) = src
            .posed_tessellation(*i, wl, pose, &Platform::default())
            .map_err(|e| e.to_string())?;
        let top = v.iter().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max);
        if !top.is_finite() {
            return Err("the hull has no geometry above its keel".into());
        }
        let over = top + 1e-4 * l_ref.max(1e-3);
        let raised = Platform {
            sinkage: over,
            trim: 0.0,
            pivot_x: 0.0,
        };
        let imp = src
            .situate_sectional(*i, wl, pose, &raised, opts)
            .map_err(|e| e.to_string())?
            .ok_or("the hull could not be cut whole")?;
        out.push(StabilityHull {
            sections: imp
                .hull
                .curves()
                .map(|(x, curve)| FullSection {
                    x: x + imp.placement.x,
                    half: curve.iter().map(|&(h, d)| (h, over - d)).collect(),
                })
                .collect(),
            y: imp.placement.y,
        });
    }
    Ok(out)
}

fn load_file(args: &[Value]) -> R<Value> {
    arity(args, 1, 2)?;
    let Value::Str(path) = &args[0] else {
        return Err("expected a file path".into());
    };
    let opts = args.get(1);
    check_opts(opts, HULL_OPTS)?;
    let wl = opt_len(opts, "waterline")?.unwrap_or(0.0);
    let bytes = std::fs::read(path.as_ref()).map_err(|e| format!("cannot read {}: {}", path, e))?;
    let lower = path.to_ascii_lowercase();
    let source: Box<dyn HullSource> = if lower.ends_with(".stl")
        || std::str::from_utf8(&bytes).is_err()
    {
        let units = match opt(opts, "units") {
            Some(v) => {
                let q = num(v)?;
                if q.dim != Dim::length() {
                    return Err("units must be a length, e.g. {units: 1 mm}".into());
                }
                q.s.exact().ok_or("units must be exact")?
            }
            None => {
                return Err(format!(
                    "{}: STL files carry no units; pass {{units: 1 mm}} (or m, in, ft)",
                    path
                ))
            }
        };
        Box::new(
            hullgeom::stl::mesh_fleet(&bytes, units, wl).map_err(|e| format!("{}: {}", path, e))?,
        )
    } else {
        let text = String::from_utf8(bytes).unwrap();
        Box::new(hullgeom::iges::source_fleet(&text, wl).map_err(|e| format!("{}: {}", path, e))?)
    };
    Ok(Value::Native(Rc::new(build(path, source, opts)?)))
}

fn wigley(args: &[Value]) -> R<Value> {
    arity(args, 3, 4)?;
    let mut dims = [0.0; 3];
    for k in 0..3 {
        let q = num(&args[k])?;
        if q.dim != Dim::length() {
            return Err("wigley(length, beam, draft) takes lengths".into());
        }
        dims[k] = q.s.exact().ok_or("wigley dimensions must be exact")?;
    }
    let opts = args.get(3);
    check_opts(opts, &["span", "stations", "rays"])?;
    let surfaces =
        hullgeom::iges::wigley_surfaces(dims[0], dims[1], dims[2]).map_err(|e| e.to_string())?;
    let source = hullgeom::iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0)
        .map_err(|e| e.to_string())?;
    Ok(Value::Native(Rc::new(build(
        "wigley",
        Box::new(source),
        opts,
    )?)))
}

// ---------------------------------------------------------------- results

fn q(s: Samples, unit: &str) -> Value {
    q_unit(s, unit)
}

fn deg(s: Samples) -> Value {
    Value::Num(Quantity::new(s, Dim::NONE, Some(UnitExpr::single("deg"))))
}

fn hydrostatics(args: &[Value]) -> R<Value> {
    arity(args, 1, 1)?;
    let hn = hull_arg(&args[0])?;
    let h = hn.as_any().downcast_ref::<BoatHull>().unwrap();
    let rho = Fluid::SEAWATER_15C.density;
    let vol = h.volume();
    let lcb = h
        .design
        .iter()
        .map(|(s, p)| s.displaced_volume() * (s.lcb_x() + p.x))
        .sum::<f64>()
        / vol;
    let awp: f64 = h.design.iter().map(|(s, _)| s.waterplane_area()).sum();
    let wetted: f64 = h.design.iter().map(|(s, _)| s.wetted_surface()).sum();
    let draft = h.design.iter().map(|(s, _)| s.draft()).fold(0.0, f64::max);
    let ex = Samples::Exact;
    Ok(Value::record(vec![
        ("waterline", q(ex(h.waterline_z), "m")),
        ("length", q(ex(h.length()), "m")),
        ("draft", q(ex(draft), "m")),
        ("volume", q(ex(vol), "m^3")),
        ("displacement", q(ex(rho * vol), "t")),
        ("lcb", q(ex(lcb), "m")),
        ("awp", q(ex(awp), "m^2")),
        ("wetted", q(ex(wetted), "m^2")),
    ]))
}

/// An item's mass, LCG, TCG and VCG (above the design waterline), checked.
fn load_of(h: &BoatHull, item: &Value) -> R<(Samples, Samples, Samples, Samples)> {
    let (m, cg) = item_parts(item)?;
    let vcg = cg[2].s.map(|z| z - h.waterline_z);
    Ok((m.s, cg[0].s.clone(), cg[1].s.clone(), vcg))
}

fn warn_check(it: &mut Interp, what: &str, r: &surrogate::Result) {
    if let Some(e) = r.check {
        if e > 1e-2 {
            it.emit(&format!(
                "warning: {} interpolated from {} solves; held-out error {:.2}% of its spread — pass {{grid: 7}} for more nodes",
                what,
                r.solves,
                100.0 * e
            ));
        }
    }
}

fn grid_opt(opts: Option<&Value>) -> R<usize> {
    Ok(opt_count(opts, "grid")?.unwrap_or(5).clamp(2, 15))
}

/// Free-trim hydrostatic float of the item's load.
fn float(it: &mut Interp, args: &[Value]) -> R<Value> {
    arity(args, 2, 3)?;
    let hn = hull_arg(&args[0])?;
    let h = hn.as_any().downcast_ref::<BoatHull>().unwrap();
    let opts = args.get(2);
    check_opts(opts, &["rho", "grid"])?;
    let rho = density(opts)?;
    let n = grid_opt(opts)?;
    let (m, lcg, _, _) = load_of(h, &args[1])?;
    let sources = h.sources();
    // Waterplane area, wetted surface and draft kink as the transom immerses.
    let r = surrogate::evaluate(&[&m, &lcg], &[n, n], &[4, 5, 6], 1e-3, |p| {
        let eq = solve_equilibrium_sectional(
            &sources,
            &LoadCase {
                mass: p[0],
                lcg: Some(p[1]),
            },
            rho,
            &h.opts,
        )
        .map_err(|e| {
            format!(
                "floating {:.4} t at LCG {:.4} m: {}",
                p[0] / 1000.0,
                p[1],
                e
            )
        })?;
        let members = &eq.fleet.members;
        Ok(vec![
            eq.sinkage,
            eq.trim,
            eq.volume,
            eq.lcb,
            eq.waterplane_area,
            members.iter().map(|(s, _)| s.wetted_surface()).sum(),
            members.iter().map(|(s, _)| s.draft()).fold(0.0, f64::max),
        ])
    })?;
    warn_check(it, "float", &r);
    let o = r.outputs;
    let disp = o[2].map(|v| v * rho);
    Ok(Value::record(vec![
        ("sinkage", q(o[0].clone(), "m")),
        ("trim", deg(o[1].clone())),
        ("draft", q(o[6].clone(), "m")),
        ("volume", q(o[2].clone(), "m^3")),
        ("displacement", q(disp, "t")),
        ("lcb", q(o[3].clone(), "m")),
        ("awp", q(o[4].clone(), "m^2")),
        ("wetted", q(o[5].clone(), "m^2")),
    ]))
}

/// The GZ curve over as many of the heels as the load can be balanced at. At
/// an extreme load (a grid corner: heaviest and highest at once) the balance
/// can fail at a large heel; the curve is then cut short just before it.
fn gz_reaching(
    hulls: &[StabilityHull],
    load: &StabilityLoad,
    opts: &GzOptions,
) -> R<hullgeom::stability::GzCurve> {
    let mut opts = opts.clone();
    loop {
        match gz_curve(hulls, load, &opts) {
            Ok(c) => return Ok(c),
            Err(e) => {
                let msg = e.to_string();
                // "GZ at 108.0°: …"
                let failed = msg
                    .split("GZ at ")
                    .nth(1)
                    .and_then(|r| r.split('°').next())
                    .and_then(|d| d.parse::<f64>().ok())
                    .map(f64::to_radians);
                let keep = match failed {
                    Some(h) if h > 0.0 => opts.heels.iter().filter(|&&x| x < h - 1e-9).count(),
                    _ => return Err(msg),
                };
                if keep < 2 || keep >= opts.heels.len() {
                    return Err(msg);
                }
                opts.heels.truncate(keep);
            }
        }
    }
}

/// The GZ curve's figures for one sample, read as hullgeom's `summarise` reads
/// them, except that a curve starting negative (an off-centre CG) is read from
/// its angle of list, where it crosses up through zero, rather than from upright.
struct Figures {
    list: f64,
    max_gz: f64,
    heel_at_max: f64,
    vanishing: f64,
    area_30: f64,
    area_40: f64,
    area_total: f64,
    /// The curve stopped (the load couldn't be balanced further) before it
    /// came back through zero: the vanishing angle is a lower bound.
    cut_short: bool,
}

fn figures(heels: &[f64], gz: &[f64], gm: f64) -> Figures {
    // Only as far as the curve was solved.
    let solved = gz.iter().take_while(|g| g.is_finite()).count();
    let cut_short = solved < heels.len();
    let (heels, gz) = (&heels[..solved], &gz[..solved]);
    let crossing = |k: usize| heels[k] + gz[k] / (gz[k] - gz[k + 1]) * (heels[k + 1] - heels[k]);
    // Where the curve starts: upright, or the angle of list.
    let mut list = if gz.first().is_some_and(|g| *g >= 0.0) {
        Some(0.0)
    } else {
        None
    };
    if list.is_none() {
        list = (0..gz.len().saturating_sub(1))
            .find(|&k| gz[k] < 0.0 && gz[k + 1] >= 0.0)
            .map(crossing);
    }
    let unstable = gm.is_nan() || gm <= 0.0 || list.is_none();
    let list = list.unwrap_or(0.0);
    let (mut max_gz, mut heel_at_max) = (f64::NEG_INFINITY, list);
    for (&h, &g) in heels.iter().zip(gz) {
        if h >= list && g > max_gz {
            max_gz = g;
            heel_at_max = h;
        }
    }
    let mut vanishing = *heels.last().unwrap_or(&0.0);
    let mut found = false;
    if unstable {
        vanishing = list;
    } else {
        for k in 0..heels.len().saturating_sub(1) {
            if heels[k] > list && gz[k] > 0.0 && gz[k + 1] <= 0.0 {
                vanishing = crossing(k);
                found = true;
                break;
            }
        }
    }
    // ∫ GZ dφ from the start of the curve, trapezoidal over the points.
    let area = |upto: f64| {
        let mut a = 0.0;
        for k in 0..heels.len().saturating_sub(1) {
            let (h0, h1) = (heels[k].max(list), heels[k + 1].min(upto));
            if h1 <= h0 || heels[k + 1] <= heels[k] {
                continue;
            }
            let at =
                |x: f64| gz[k] + (gz[k + 1] - gz[k]) * (x - heels[k]) / (heels[k + 1] - heels[k]);
            a += 0.5 * (h1 - h0) * (at(h0) + at(h1));
        }
        a
    };
    Figures {
        list,
        max_gz: if max_gz.is_finite() { max_gz } else { 0.0 },
        heel_at_max,
        vanishing,
        area_30: area(30f64.to_radians()),
        area_40: area(40f64.to_radians()),
        area_total: if unstable { 0.0 } else { area(vanishing) },
        cut_short: cut_short && !found && !unstable,
    }
}

/// Large-angle stability with free trim: GM and the GZ curve with its figures.
fn stability(it: &mut Interp, args: &[Value]) -> R<Value> {
    arity(args, 2, 3)?;
    let hn = hull_arg(&args[0])?;
    let h = hn.as_any().downcast_ref::<BoatHull>().unwrap();
    let opts = args.get(2);
    check_opts(opts, &["rho", "grid", "max_heel", "step"])?;
    let rho = density(opts)?;
    let n = grid_opt(opts)?;
    let angle = |name: &str, default: f64| -> R<f64> {
        match opt(opts, name) {
            None => Ok(default.to_radians()),
            Some(v) => {
                let q = num(v)?;
                if !q.dim.is_none() {
                    return Err(format!("{} must be an angle, e.g. 60 deg", name));
                }
                q.s.exact().ok_or_else(|| format!("{} must be exact", name))
            }
        }
    };
    let (max_heel, step) = (angle("max_heel", 180.0)?, angle("step", 2.0)?);
    if !(step > 0.0) {
        return Err("step must be positive".into());
    }
    let nh = (max_heel / step + 1e-9).floor() as usize + 1;
    let heels: Vec<f64> = (0..nh).map(|k| k as f64 * step).collect();
    let (m, lcg, tcg, vcg) = load_of(h, &args[1])?;
    let gopts = GzOptions {
        heels: heels.clone(),
        density: rho,
        ..GzOptions::default()
    };
    // VCG enters the free-trim balance only through the trim, so 3 nodes suffice.
    let r = surrogate::evaluate(&[&m, &lcg, &vcg], &[n, n, 3], &[], 1e-3, |p| {
        let load = StabilityLoad {
            mass: p[0],
            lcg: p[1],
            vcg: p[2],
        };
        let c = gz_reaching(&h.whole, &load, &gopts).map_err(|e| {
            format!(
                "GZ for {:.4} t at LCG {:.4} m, VCG {:.4} m above the DWL: {}",
                p[0] / 1000.0,
                p[1],
                p[2],
                e
            )
        })?;
        let mut out = vec![c.gm];
        out.extend(c.points.iter().map(|q| q.gz));
        // Heels this load couldn't be balanced at stay NaN.
        out.resize(1 + heels.len(), f64::NAN);
        Ok(out)
    })?;
    warn_check(it, "GZ", &r);
    let o = r.outputs;
    // hullgeom's load has no TCG. An off-centre CG shortens the arm by
    // |TCG|·cos φ heeling toward it, so the curve is the platform's (port and
    // starboard symmetric) heeled to that, the worse, side.
    let gz_cols: Vec<Samples> = (0..nh)
        .map(|k| {
            let c = heels[k].cos();
            Samples::from_fn(&[&o[k + 1], &tcg], |i| {
                o[k + 1].get(i) - tcg.get(i).abs() * c
            })
        })
        .collect();
    let mut inputs: Vec<&Samples> = gz_cols.iter().collect();
    inputs.push(&o[0]);
    let n_out = 8;
    let figs = napkin::builtins::per_sample(&inputs, n_out, |i| {
        let gz: Vec<f64> = gz_cols.iter().map(|s| s.get(i)).collect();
        let f = figures(&heels, &gz, o[0].get(i));
        Ok(vec![
            f.max_gz,
            f.heel_at_max,
            f.vanishing,
            f.area_30,
            f.area_40,
            f.area_total,
            f.cut_short as u8 as f64,
            f.list,
        ])
    })?;
    match &figs[6] {
        Samples::Exact(x) if *x != 0.0 => it.emit(
            "warning: the GZ curve couldn't be balanced as far as it returns to zero; \
             the vanishing angle is a lower bound",
        ),
        Samples::Dist(v) if v.iter().any(|x| *x != 0.0) => it.emit(&format!(
            "warning: for {} of {} samples (the most extreme loads) the GZ curve couldn't be balanced \
             as far as it returns to zero; their vanishing angle is a lower bound",
            v.iter().filter(|x| **x != 0.0).count(),
            v.len()
        )),
        _ => {}
    }
    Ok(Value::record(vec![
        ("gm", q(o[0].clone(), "m")),
        ("list", deg(figs[7].clone())),
        ("max_gz", q(figs[0].clone(), "m")),
        ("heel_at_max", deg(figs[1].clone())),
        ("vanishing", deg(figs[2].clone())),
        ("area_30", q(figs[3].clone(), "m·rad")),
        ("area_40", q(figs[4].clone(), "m·rad")),
        ("area_total", q(figs[5].clone(), "m·rad")),
        (
            "heels",
            Value::list(heels.iter().map(|&h| deg(Samples::Exact(h))).collect()),
        ),
        (
            "gz",
            Value::list(gz_cols.into_iter().map(|s| q(s, "m")).collect()),
        ),
    ]))
}

// ---------------------------------------------------------------- measuring

/// How closely the flat-faceted mesh follows a curved hull (see `geometry`).
const MESH_ACCURACY: f64 = 3e-3;

/// The meshes an option `hull: k` picks (all of them by default). For a
/// catamaran, hull 0 is the one at +y.
fn pick<'a>(h: &'a BoatHull, opts: Option<&Value>) -> R<Vec<&'a Mesh>> {
    match opt_count0(opts, "hull")? {
        None => Ok(h.meshes.iter().collect()),
        Some(k) if k < h.meshes.len() => Ok(vec![&h.meshes[k]]),
        Some(k) => Err(format!(
            "hull {} out of range: this platform has {} (numbered from 0)",
            k,
            h.meshes.len()
        )),
    }
}

fn opt_count0(opts: Option<&Value>, name: &str) -> R<Option<usize>> {
    match opt(opts, name) {
        None => Ok(None),
        Some(v) => match num(v)?.s.exact() {
            Some(x) if x >= 0.0 && x.fract() == 0.0 => Ok(Some(x as usize)),
            _ => Err(format!("{} must be an exact whole number", name)),
        },
    }
}

/// A length's samples, checked.
fn length_samples(v: &Value, what: &str) -> R<Samples> {
    let q = num(v)?;
    if q.dim != Dim::length() {
        return Err(format!(
            "{} must be a length, got {}",
            what,
            q.dim.describe()
        ));
    }
    Ok(q.s)
}

fn point(c: [Samples; 3]) -> Value {
    let [x, y, z] = c;
    Value::list(vec![q(x, "m"), q(y, "m"), q(z, "m")])
}

/// Area and centroid of the shell within a box: `{x: [a, b], y: [..], z: [..]}`,
/// each axis optional and each bound possibly uncertain.
fn shell(it: &mut Interp, args: &[Value]) -> R<Value> {
    arity(args, 1, 2)?;
    let hn = hull_arg(&args[0])?;
    let h = hn.as_any().downcast_ref::<BoatHull>().unwrap();
    let opts = args.get(1);
    check_opts(opts, &["x", "y", "z", "hull", "grid"])?;
    let meshes = pick(h, opts)?;
    // Each finite bound is an input; (axis, end) says where it goes.
    let mut inputs: Vec<Samples> = vec![];
    let mut slots: Vec<(usize, usize)> = vec![];
    for (k, name) in ["x", "y", "z"].iter().enumerate() {
        if let Some(v) = opt(opts, name) {
            let Value::List(ends) = v else {
                return Err(format!("{} must be a range [from, to]", name));
            };
            if ends.len() != 2 {
                return Err(format!("{} must be a range [from, to]", name));
            }
            for (e, end) in ends.iter().enumerate() {
                inputs.push(length_samples(end, name)?);
                slots.push((k, e));
            }
        }
    }
    let n = grid_opt(opts)?;
    let refs: Vec<&Samples> = inputs.iter().collect();
    let r = surrogate::evaluate(&refs, &vec![n; refs.len()], &[], MESH_ACCURACY, |p| {
        let mut b: Bounds = [[f64::NEG_INFINITY, f64::INFINITY]; 3];
        for (&(k, e), &v) in slots.iter().zip(p) {
            b[k][e] = v;
        }
        let (a, c) = geometry::shell(&meshes, &b);
        Ok(vec![a, c[0], c[1], c[2]])
    })?;
    warn_check(it, "shell", &r);
    let [a, x, y, z] = <[Samples; 4]>::try_from(r.outputs).unwrap();
    Ok(Value::record(vec![
        ("area", q(a, "m^2")),
        ("centroid", point([x, y, z])),
    ]))
}

/// The transverse section at `x`: its area, centroid and girth (the length of
/// shell around it), cut off at `below` if given. An open section (a hull with
/// no deck) is closed straight across its top.
fn slice(it: &mut Interp, args: &[Value]) -> R<Value> {
    arity(args, 2, 3)?;
    let hn = hull_arg(&args[0])?;
    let h = hn.as_any().downcast_ref::<BoatHull>().unwrap();
    let opts = args.get(2);
    check_opts(opts, &["below", "hull", "grid"])?;
    let meshes = pick(h, opts)?;
    let xs = length_samples(&args[1], "the station")?;
    let below = match opt(opts, "below") {
        Some(v) => Some(length_samples(v, "below")?),
        None => None,
    };
    let outlines = |x: f64, zc: Option<f64>| -> Vec<geometry::Outline> {
        let all = geometry::slice(&meshes, x);
        match zc {
            Some(z) => all.iter().filter_map(|o| o.below(z)).collect(),
            None => all,
        }
    };
    // Without a cut, its axis is a placeholder single node.
    let zc = below.clone().unwrap_or(Samples::Exact(f64::NAN));
    let n = grid_opt(opts)?;
    let has_below = below.is_some();
    let r = surrogate::evaluate(&[&xs, &zc], &[n, n], &[], MESH_ACCURACY, |p| {
        let zc = if has_below { Some(p[1]) } else { None };
        let (mut a, mut my, mut mz, mut girth) = (0.0, 0.0, 0.0, 0.0);
        for o in outlines(p[0], zc) {
            let (ai, c, g) = o.measure();
            if ai > 0.0 {
                a += ai;
                my += ai * c[0];
                mz += ai * c[1];
            }
            girth += g;
        }
        let (cy, cz) = if a > 0.0 {
            (my / a, mz / a)
        } else {
            (f64::NAN, f64::NAN)
        };
        Ok(vec![a, cy, cz, girth])
    })?;
    warn_check(it, "slice", &r);
    let [a, y, z, g] = <[Samples; 4]>::try_from(r.outputs).unwrap();
    let mut fields = vec![
        ("x", q(xs.clone(), "m")),
        ("area", q(a, "m^2")),
        ("centroid", point([xs.clone(), y, z])),
        ("girth", q(g, "m")),
    ];
    // The outline itself, when the station is certain: [x, y, z] points per
    // closed loop, ready for napkin's own geometry functions.
    if let (Some(x), true) = (xs.exact(), below.as_ref().is_none_or(|b| b.is_exact())) {
        let zc = below.as_ref().and_then(|b| b.exact());
        let loops = outlines(x, zc)
            .into_iter()
            .map(|o| {
                Value::list(
                    o.pts
                        .iter()
                        .map(|p| {
                            point([
                                Samples::Exact(x),
                                Samples::Exact(p[0]),
                                Samples::Exact(p[1]),
                            ])
                        })
                        .collect(),
                )
            })
            .collect();
        fields.push(("outlines", Value::list(loops)));
    }
    Ok(Value::record(fields))
}

/// Register the hull functions on an interpreter.
pub fn install(it: &mut Interp) {
    it.register("hull", |_, a| load_file(&a));
    it.register("wigley", |_, a| wigley(&a));
    it.register("hydrostatics", |_, a| hydrostatics(&a));
    it.register("float", |it, a| float(it, &a));
    it.register("stability", |it, a| stability(it, &a));
    it.register("shell", |it, a| shell(it, &a));
    it.register("slice", |it, a| slice(it, &a));
}
