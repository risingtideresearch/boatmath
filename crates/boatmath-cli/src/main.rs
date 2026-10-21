//! `boatmath` — Unix-style tools for doing physics on boats. Records are
//! JSON, one per line; each command reads a stream, writes it through, and
//! adds what it made, so the stream at the end of a pipeline holds the
//! whole computation:
//!
//! ```sh
//! boatmath hull guillemot.igs --waterline 0.12 \
//!   | boatmath case --mass 150,200 \
//!   | boatmath study --froude 0.2:0.6:0.05 \
//!   | boatmath run -j 8 \
//!   | boatmath table study.case.params.mass study.params.froude forces.rt
//! ```
//!
//! See docs/boatmath-cli.md.

mod cache;
mod cad;
mod list;
mod path;
mod pick;
mod plot;
mod props;
mod records;
mod run;
mod stream;
mod units;
mod views;

use boatmath::params::{CaseParams, Closure, Sea, StudyParams, Waves};
use boatmath::LoftRequest;
use cache::Cache;
use clap::{Parser, Subcommand};
use serde_json::Value;
use std::io::Write;
use std::path::PathBuf;
use stream::{short, Out, Stream};

#[derive(Parser)]
#[command(
    name = "boatmath",
    version,
    about = "Unix-style tools for doing physics on boats: JSON records in a pipe"
)]
struct Cli {
    /// Memoize expensive steps here (default $BOATMATH_CACHE; none without).
    #[arg(long, global = true)]
    cache: Option<PathBuf>,
    /// Read records from these files instead of stdin; repeat for several.
    #[arg(long = "input", short = 'i', global = true)]
    inputs: Vec<PathBuf>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Hulls from files (IGES, STL, camber's JSON): a hull record per file and
    /// setting, its geometry inline as JSON (B-spline patches, or triangles
    /// from an STL), in metres with the design waterline at z = 0.
    Hull {
        files: Vec<PathBuf>,
        #[arg(long)]
        name: Option<String>,
        /// Design waterline height in the file's frame [m] (IGES; a camber
        /// document has its own); a LIST.
        #[arg(long, allow_hyphen_values = true)]
        waterline: Option<String>,
        /// Centreplane override [m] (IGES).
        #[arg(long, allow_hyphen_values = true)]
        centerplane: Option<f64>,
        /// Stations along the hull (IGES).
        #[arg(long)]
        stations: Option<usize>,
        /// Rays across each section (IGES).
        #[arg(long)]
        rays: Option<usize>,
        /// Units of an STL (or in place of a camber document's): mm, cm, m,
        /// in, ft, or metres per unit.
        #[arg(long)]
        units: Option<String>,
    },
    /// Hulls scaled from the stream's, about the design waterline.
    Scale {
        /// Scale the whole hull by this factor; a LIST.
        #[arg(long)]
        by: Option<String>,
        /// Scale beam and draft only by this factor; a LIST.
        #[arg(long)]
        beam: Option<String>,
        /// Scale to carry this mass [kg] at the design waterline: uniformly,
        /// or with --keep-length beam and draft only; a LIST.
        #[arg(long)]
        mass: Option<String>,
        #[arg(long)]
        keep_length: bool,
        #[arg(long)]
        name: Option<String>,
    },
    /// Cases on the stream's hulls: a platform and load, one per hull and
    /// combination of the LISTs. A definition only; see `statics`.
    Case {
        /// Make a catamaran with this centre span [m]; a LIST.
        #[arg(long)]
        span: Option<String>,
        /// Load [kg], carried by sinking; default the design displacement; a LIST.
        #[arg(long)]
        mass: Option<String>,
        /// Longitudinal centre of gravity [m]; default the LCB; a LIST.
        #[arg(long, allow_hyphen_values = true)]
        lcg: Option<String>,
        /// Centre of gravity above the design waterline [m]; a LIST.
        #[arg(long, allow_hyphen_values = true)]
        vcg: Option<String>,
        /// Roll radius of gyration [m]; a LIST.
        #[arg(long)]
        kxx: Option<String>,
        /// Pitch radius of gyration, a fraction of the length; a LIST.
        #[arg(long)]
        kyy: Option<String>,
        /// Yaw radius of gyration [m]; a LIST.
        #[arg(long)]
        kzz: Option<String>,
        /// Roll damping, a fraction of critical; a LIST.
        #[arg(long)]
        roll_damping: Option<String>,
        #[arg(long)]
        name: Option<String>,
    },
    /// The stream's cases with a drive: a leg (and pod) hung from the hull,
    /// or a shaft and its strut, cut as members of their own, carrying the
    /// propeller. A stock drive
    /// (`boatmath mounts`) gives its dimensions; the flags give or override
    /// them. Lengths take units (m, mm, in, …).
    Mount {
        /// A stock drive by name.
        #[arg(long)]
        stock: Option<String>,
        /// saildrive (a leg through the bottom to a pod), outboard (a leg
        /// from the transom), pod (a pod on a short strut), or shaft (an
        /// inclined shaft from the hull bottom, on a P-strut).
        #[arg(long)]
        kind: Option<String>,
        /// The leg's mid-chord forward of the hull's aft end (an outboard's
        /// astern of it, negative); a shaft drive's propeller.
        #[arg(long, allow_hyphen_values = true)]
        x: Option<String>,
        /// Out from the hull's centreplane.
        #[arg(long, default_value = "0", allow_hyphen_values = true)]
        y: String,
        /// Two drives on each hull, at ±--y (twin screws).
        #[arg(long)]
        pair: bool,
        /// The shaft below the keel at the leg (saildrive, pod) or at the
        /// propeller (shaft), or below the transom's bottom (outboard).
        #[arg(long)]
        shaft_depth: Option<String>,
        /// The leg's (a shaft drive's strut's) chord.
        #[arg(long)]
        chord: Option<String>,
        /// The leg's (strut's) thickness.
        #[arg(long)]
        thickness: Option<String>,
        /// A shaft drive's shaft diameter.
        #[arg(long)]
        shaft_diameter: Option<String>,
        /// A shaft drive's strut, its mid-chord ahead of the propeller;
        /// default half its chord, 0.1 m and 0.15 of --prop-diameter.
        #[arg(long)]
        strut_ahead: Option<String>,
        #[arg(long)]
        pod_length: Option<String>,
        #[arg(long)]
        pod_diameter: Option<String>,
        /// The pod's nose ahead of the leg's mid-chord.
        #[arg(long, allow_hyphen_values = true)]
        nose_ahead: Option<String>,
        /// The propeller aft of the pod's nose.
        #[arg(long, allow_hyphen_values = true)]
        prop_from_nose: Option<String>,
        /// No pod: the leg alone.
        #[arg(long)]
        no_pod: bool,
        /// The propeller ahead of the pod, pulling.
        #[arg(long)]
        tractor: bool,
        /// The drive's own propeller; `prop` defaults --d-max to it.
        #[arg(long)]
        prop_diameter: Option<String>,
        /// The shaft's angle, bow up of the baseline [deg].
        #[arg(long, default_value = "0", allow_hyphen_values = true)]
        shaft_angle: f64,
        #[arg(long)]
        name: Option<String>,
    },
    /// The stock drives `mount --stock` takes, a record per drive.
    Mounts,
    /// Statics of the stream's hulls (hydrostatics at the design waterline)
    /// and cases (the float at rest, roll stability, the GZ curve).
    Statics,
    /// Quasi-static GZ in a regular wave (platform::wave_gz).
    #[command(name = "gz-waves")]
    GzWaves {
        /// Wavelength(s) [m]; a LIST.
        #[arg(long)]
        length: String,
        /// Wave height(s) [m]; a LIST.
        #[arg(long)]
        height: String,
    },
    /// Studies on the stream's cases: one per case and combination of the LISTs.
    Study {
        /// Length Froude numbers; a LIST.
        #[arg(long)]
        froude: String,
        /// Hold the attitude at rest instead of floating at speed.
        #[arg(long)]
        hold: bool,
        /// Transom closure: auto[:FRACTION] (the default), ballistic[:COEFF],
        /// fixed:LENGTH, or off.
        #[arg(long)]
        closure: Option<String>,
        /// Free-surface grid columns (calm water).
        #[arg(long)]
        grid: Option<usize>,
        /// How far astern the free surface reaches, in lengths (calm water;
        /// default 1.5).
        #[arg(long)]
        astern: Option<f64>,
        /// In waves from these headings [deg, 180 head seas]; a LIST.
        #[arg(long, allow_hyphen_values = true)]
        waves: Option<String>,
        /// Wavelengths over the length (with --waves); a LIST.
        #[arg(long)]
        lambdas: Option<String>,
        /// An irregular sea (with --waves): bretschneider:hs=H,tp=T or
        /// jonswap:hs=H,tp=T[,gamma=G].
        #[arg(long)]
        sea: Option<String>,
    },
    /// Compute the stream's studies; each result is written as it is ready.
    Run {
        /// Studies at once; default the number of cores.
        #[arg(short, long)]
        jobs: Option<usize>,
        /// No progress on stderr.
        #[arg(short, long)]
        quiet: bool,
    },
    /// The best B-series propeller for each calm-water result in the stream
    /// (its speed, and the thrust R_t / ((1 − t) cos ε) along its case's
    /// drive), or for --speed and --thrust.
    Prop {
        /// Total thrust (N, kN, kgf, lbf), instead of reading results.
        #[arg(long)]
        thrust: Option<String>,
        /// Ship speed (m/s, kn, mph, km/h), with --thrust.
        #[arg(long)]
        speed: Option<String>,
        /// Shafts the thrust is split across; default the case's drives, or
        /// one per hull.
        #[arg(long)]
        shafts: Option<u32>,
        /// Wake fraction w: the propeller sees V (1 − w). `auto`: from
        /// potential flow at the result's attitude, on its case's drive.
        #[arg(long, default_value = "0")]
        wake: String,
        /// Thrust deduction t; `auto` as for --wake.
        #[arg(long, default_value = "0")]
        thrust_deduction: String,
        /// Largest diameter (m, mm, in, …); default the drive's own
        /// propeller's, if it comes with one.
        #[arg(long)]
        d_max: Option<String>,
        #[arg(long, default_value = "40mm")]
        d_min: String,
        /// Blade counts to consider; a LIST.
        #[arg(long, default_value = "2,3,4,5,6,7")]
        blades: String,
        /// Shaft immersion, for cavitation; default the drive's. Needed
        /// without one.
        #[arg(long)]
        depth: Option<String>,
        /// Keller's margin k: 0.2 single screw, 0.1 twin, 0 fast.
        #[arg(long, default_value_t = 0.2)]
        keller_k: f64,
        /// Allow cavitating designs.
        #[arg(long)]
        no_cavitation: bool,
        /// Leave the series at its Re = 2e6 fit.
        #[arg(long)]
        no_re_correct: bool,
        /// Keep each blade count to the series' own blade-area range.
        #[arg(long)]
        strict_ear: bool,
        /// A second point the same propeller must reach: the same study's
        /// result at this Froude number, in the stream.
        #[arg(long)]
        top_froude: Option<f64>,
        /// Or the second point's speed and total thrust.
        #[arg(long)]
        top_speed: Option<String>,
        #[arg(long)]
        top_thrust: Option<String>,
        /// Keep the whole feasible curve, not just to 1.6× the least power.
        #[arg(long)]
        full_range: bool,
    },
    /// The motors that can drive each prop in the stream, each at its best
    /// reduction and point on the prop's curve, ranked: a drive record each.
    Match {
        /// Rank by electrical power, mass or price.
        #[arg(long, default_value = "power")]
        rank_by: String,
        /// Every winding, not just the best of each family.
        #[arg(long)]
        all_windings: bool,
        /// Direct drive only (a ratio of 1).
        #[arg(long)]
        direct_only: bool,
        /// A fixed reduction ratio, instead of each motor's best.
        #[arg(long)]
        ratio: Option<f64>,
        /// The largest reduction considered.
        #[arg(long, default_value_t = 25.0)]
        ratio_max: f64,
        /// Efficiency of one reduction stage.
        #[arg(long, default_value_t = propeller::motor::GEAR_ETA_DEFAULT)]
        gear_eta: f64,
        /// Leave out the controller's loss (motor efficiency only).
        #[arg(long)]
        no_controller_loss: bool,
        /// Hold the design point on peak ratings, not continuous.
        #[arg(long)]
        peak: bool,
        /// Don't require the prop's second point.
        #[arg(long)]
        no_top: bool,
        #[command(flatten)]
        filter: MotorFilter,
        /// Motors from this file instead of the vendored database: a
        /// motors.js, its JSON, or motor records (as `boatmath motors`
        /// writes, filtered with jq, say).
        #[arg(long)]
        motors: Option<PathBuf>,
    },
    /// The motor database, a record per motor.
    Motors {
        #[command(flatten)]
        filter: MotorFilter,
    },
    /// Filter the stream by record, keeping the ancestors of what it keeps.
    Pick {
        /// Keep records of this type (and what they refer to).
        kind: Option<String>,
        /// Only the first N of them.
        #[arg(long)]
        first: Option<usize>,
        /// Only those where PATH reads VALUE; repeat for several.
        #[arg(long = "where")]
        wheres: Vec<String>,
        /// Leave out records of this type; repeat for several.
        #[arg(long)]
        drop: Vec<String>,
    },
    /// The stream as a table, a column per FIELD path (`forces.rt`,
    /// `study.case.params.span`, `|heave|`, ...); tab-separated.
    Table {
        #[arg(required = true)]
        fields: Vec<String>,
        /// A row per record of this type (default: the type the first field
        /// reads on).
        #[arg(long = "type")]
        kind: Option<String>,
        /// A row per element of the array at this path, its fields read
        /// from the element first and then the record.
        #[arg(long)]
        explode: Option<String>,
        /// Comma-separated, quoted where needed.
        #[arg(long)]
        csv: bool,
        /// No header row.
        #[arg(long)]
        no_header: bool,
    },
    /// The stream as an SVG line plot.
    Plot {
        /// The x field.
        #[arg(short)]
        x: String,
        /// A y field; repeat for several.
        #[arg(short, required = true)]
        y: Vec<String>,
        /// A series per value of this field (or fields, repeated).
        #[arg(long)]
        by: Vec<String>,
        /// A point per record of this type (see `table`).
        #[arg(long = "type")]
        kind: Option<String>,
        /// A point per element of the array at this path (see `table`).
        #[arg(long)]
        explode: Option<String>,
        #[arg(long)]
        title: Option<String>,
        #[arg(long)]
        xlabel: Option<String>,
        #[arg(long)]
        ylabel: Option<String>,
        /// Write here instead of stdout.
        #[arg(short)]
        o: Option<PathBuf>,
    },
    /// The wake of the stream's calm-water results, from above, as SVG.
    Wake {
        /// Write here; with several results, a pattern naming each by
        /// `{id8}` or `{froude}` (wake-{froude}.svg).
        #[arg(short)]
        o: Option<String>,
        /// Colour scale ±RANGE [m]; default the 99th percentile of |ζ|.
        #[arg(long)]
        range: Option<f64>,
        #[arg(long)]
        title: Option<String>,
    },
    /// The pressure on the hulls of the stream's calm-water results, from
    /// below, as SVG.
    Pressure {
        /// As for `wake`.
        #[arg(short)]
        o: Option<String>,
        /// Colour scale ±RANGE in C_p; default the 99th percentile of |C_p|.
        #[arg(long)]
        range: Option<f64>,
        #[arg(long)]
        title: Option<String>,
    },
    /// The hull of the stream's calm-water results in profile at its
    /// attitude, with the wave along its side, as SVG.
    Profile {
        /// As for `wake`.
        #[arg(short)]
        o: Option<String>,
        /// Draw the wave this many times its height.
        #[arg(long, default_value_t = 1.0)]
        wave_scale: f64,
        #[arg(long)]
        title: Option<String>,
    },
    /// The stream's props (or, with none, its calm-water results) as IGES
    /// for CAD: the hulls at their attitude, the free surface, and a prop's
    /// discs, each on its own level (1 hulls, 2 water, 3 props). Its
    /// statics too: a case at rest and heeled along its GZ curve.
    Cad {
        /// With statics: also heel to each of these angles [deg]; a LIST.
        #[arg(long)]
        heel: Option<String>,
        /// Write here; with several records, a pattern naming each by
        /// `{id8}` or `{froude}`.
        #[arg(short)]
        o: Option<String>,
        /// Leave out the free surface.
        #[arg(long)]
        no_water: bool,
        /// Draw the waves this many times their height.
        #[arg(long, default_value_t = 1.0)]
        wave_scale: f64,
        /// Draw a prop's propellers as plain discs, not their blades.
        #[arg(long)]
        prop_discs: bool,
        /// Its propellers turn anticlockwise seen from astern.
        #[arg(long)]
        left_handed: bool,
    },
}

#[derive(clap::Args)]
struct MotorFilter {
    /// Only this maker's motors; repeat for several.
    #[arg(long)]
    vendor: Vec<String>,
    /// Only motors with a published efficiency map (tiers A, A-).
    #[arg(long)]
    mapped_only: bool,
    /// No stacked (multi-machine) motors.
    #[arg(long)]
    single_only: bool,
    /// Heaviest motor [kg].
    #[arg(long)]
    max_mass: Option<f64>,
    /// Largest outside diameter (m, mm, in, …).
    #[arg(long)]
    max_od: Option<String>,
    /// Dearest motor [USD].
    #[arg(long)]
    max_price: Option<f64>,
}

impl MotorFilter {
    fn filter(&self) -> Result<props::Filter, String> {
        Ok(props::Filter {
            vendors: self.vendor.clone(),
            mapped_only: self.mapped_only,
            single_only: self.single_only,
            max_mass: self.max_mass,
            max_od: self
                .max_od
                .as_deref()
                .map(units::length)
                .transpose()?
                .map(|m| m * 1e3),
            max_price: self.max_price,
        })
    }
}

/// A motor database from a file: propopt's motors.js or its JSON, or motor
/// records (their maps and controller efficiency the vendored database's).
fn motor_database(path: &PathBuf) -> Result<propeller::motor::Database, String> {
    use propeller::motor::Database;
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let start = text.find('{').unwrap_or(0);
    if let Ok(v) = serde_json::from_str::<Value>(text[start..].trim().trim_end_matches(';')) {
        if v.get("motors").is_some() {
            return Database::from_value(v);
        }
    }
    if text.trim_start().starts_with("const") || text.contains("MOTOR_DB") {
        return Database::from_js(&text);
    }
    let motors: Vec<Value> = serde_json::Deserializer::from_str(&text)
        .into_iter::<Value>()
        .collect::<Result<_, _>>()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    let base = &Database::vendored().raw;
    Database::from_value(serde_json::json!({
        "motors": motors,
        "maps": base["maps"],
        "controller_eta": base["controller_eta"],
    }))
}

/// The views a picture command draws.
enum View {
    Wake(Option<f64>),
    Pressure(Option<f64>),
    Profile(f64),
}

/// A file name from a pattern: `{id8}` and `{froude}` filled in.
fn named(pattern: &str, id: &str, froude: f64) -> String {
    pattern
        .replace("{id8}", &id[..8.min(id.len())])
        .replace("{froude}", &format!("{froude}"))
}

/// Several outputs need a pattern to name each.
fn check_pattern(n: usize, o: &Option<String>, what: &str) -> Result<(), String> {
    if n > 1 && !o.as_deref().is_some_and(|p| p.contains('{')) {
        return Err(format!(
            "{n} {what}: give -o a pattern with {{id8}} or {{froude}} to name each"
        ));
    }
    Ok(())
}

/// Draw `view` of each calm-water result in the stream, to `o` (or stdout).
fn draw(
    s: &Stream,
    out: &mut impl Write,
    view: View,
    o: Option<String>,
    title: Option<String>,
) -> Result<usize, String> {
    let results: Vec<&Value> = s
        .of_type("result")
        .filter(|r| r["kind"] == "calm")
        .collect();
    if results.is_empty() {
        return Err("no calm-water results in the stream".into());
    }
    check_pattern(results.len(), &o, "results")?;
    let res = path::Resolver::new(s);
    let mut failed = 0;
    for r in results {
        let id = r["id"].as_str().unwrap_or("");
        let picture = (|| -> Result<String, String> {
            let field = views::Field::parse(s.follow(r, "field")?)?;
            let name = path::cell(&res.get(r, "study.case.hull.name"));
            let froude = r["froude"].as_f64().unwrap_or(0.0);
            let f = &r["forces"];
            let num = |k: &str| f[k].as_f64().map_or("?".into(), |v| format!("{v:.1}"));
            let title = title
                .clone()
                .unwrap_or_else(|| format!("{name} · Fn {froude}"));
            let sub = format!("R_t {} N · R_w {} N · {}", num("rt"), num("rw"), short(id));
            Ok(match view {
                View::Wake(range) => views::wake(&field, &title, &sub, range),
                View::Pressure(range) => views::pressure(&field, &title, &sub, range),
                View::Profile(scale) => {
                    let sec = s.follow(r, "sections")?;
                    let attitude = (
                        sec["attitude"]["sinkage"].as_f64().unwrap_or(0.0),
                        sec["attitude"]["trim_rad"].as_f64().unwrap_or(0.0),
                    );
                    let case = s.follow(s.follow(r, "study")?, "case")?;
                    let records::CaseSource { params, src, .. } = records::case_source(s, case)?;
                    let meshes = boatmath::sections::meshes(
                        &src.file_name,
                        src.bytes,
                        &src.import,
                        &params,
                        attitude,
                    )?;
                    // The first hull, and the drive's parts on it.
                    let first = meshes
                        .iter()
                        .find(|m| m.role == "hull")
                        .ok_or("the case has no hull")?;
                    let parts: Vec<&[[f64; 3]]> = meshes
                        .iter()
                        .filter(|m| m.role != "hull" && m.copy == first.copy)
                        .map(|m| m.mesh.0.as_slice())
                        .collect();
                    views::profile(&field, &first.mesh.0, &parts, &title, &sub, attitude, scale)?
                }
            })
        })();
        match picture {
            Ok(svg) => match &o {
                Some(p) => {
                    let path = named(p, id, r["froude"].as_f64().unwrap_or(0.0));
                    std::fs::write(&path, svg).map_err(|e| format!("{path}: {e}"))?;
                }
                None => out.write_all(svg.as_bytes()).map_err(stdout_err)?,
            },
            Err(e) => {
                eprintln!("boatmath: result {}: {e}", short(id));
                failed += 1;
            }
        }
    }
    Ok(failed)
}

/// The records a table or plot has a row for: of `kind`, or else of the
/// most derived type the first field (or the exploded array) reads on.
fn rows_of<'a>(
    s: &'a Stream,
    kind: &Option<String>,
    first: &str,
    explode: Option<&str>,
) -> Vec<&'a Value> {
    let first = explode.unwrap_or(first);
    if let Some(k) = kind {
        return s.of_type(k).collect();
    }
    let res = path::Resolver::new(s);
    const ORDER: [&str; 11] = [
        "drive", "prop", "result", "gz_waves", "statics", "study", "case", "hull", "motor",
        "sections", "field",
    ];
    for k in ORDER {
        let rs: Vec<&Value> = s.of_type(k).collect();
        if rs
            .iter()
            .any(|r| res.get(r, first).is_some_and(|v| !v.is_null()))
        {
            return rs;
        }
    }
    Vec::new()
}

/// A field path's last segment, to name a series by: `span` for
/// `study.case.params.span`, `|heave|` for `|seakeeping.heave|`.
fn short_field(f: &str) -> String {
    match f.strip_prefix('|').and_then(|p| p.strip_suffix('|')) {
        Some(inner) => format!("|{}|", short_field(inner)),
        None => f.rsplit('.').next().unwrap_or(f).to_string(),
    }
}

/// A failed write to stdout. A reader that stops early (`| head`) closes
/// the pipe; like any Unix filter, stop quietly then.
pub(crate) fn stdout_err(e: std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::BrokenPipe {
        std::process::exit(0);
    }
    format!("stdout: {e}")
}

fn csv_cell(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

fn list(s: &Option<String>, what: &str) -> Result<Option<Vec<f64>>, String> {
    s.as_deref()
        .map(|s| list::parse(s).map_err(|e| format!("--{what}: {e}")))
        .transpose()
}

fn parse_lengths(s: &str, what: &str) -> Result<Vec<f64>, String> {
    let num = |t: &str| units::length(t).map_err(|e| format!("--{what}: {e}"));
    let mut out = Vec::new();
    for part in s.split(',').filter(|p| !p.trim().is_empty()) {
        let bits: Vec<&str> = part.split(':').collect();
        match bits[..] {
            [x] => out.push(num(x)?),
            [a, b, step] => {
                let (a, b, step) = (num(a)?, num(b)?, num(step)?);
                if step <= 0.0 || b < a {
                    return Err(format!(
                        "--{what} {part}: expected start:stop:step, start ≤ stop, step > 0"
                    ));
                }
                let n = ((b - a) / step + 1e-9).floor() as usize;
                if n > 10_000 {
                    return Err(format!("--{what} {part}: more than 10000 values"));
                }
                let digits = (-step.log10().floor()).max(0.0) as i32 + 6;
                let scale = 10f64.powi(digits);
                out.extend((0..=n).map(|i| ((a + i as f64 * step) * scale).round() / scale));
            }
            _ => {
                return Err(format!(
                    "--{what} {part}: expected a number or start:stop:step"
                ))
            }
        }
    }
    if out.is_empty() {
        return Err(format!("--{what} {s:?}: an empty list"));
    }
    Ok(out)
}

fn closure(s: &str) -> Result<Closure, String> {
    let (kind, arg) = s.split_once(':').unwrap_or((s, ""));
    let num = |a: &str| {
        a.parse::<f64>()
            .map_err(|_| format!("--closure {s:?}: expected a number after the colon"))
    };
    Ok(match kind {
        "ballistic" if arg.is_empty() => Closure::Ballistic { coeff: None },
        "ballistic" => Closure::Ballistic {
            coeff: Some(num(arg)?),
        },
        "auto" if arg.is_empty() => Closure::Auto {
            coeff: None,
            fraction: None,
        },
        "auto" => Closure::Auto {
            coeff: None,
            fraction: Some(num(arg)?),
        },
        "fixed" => Closure::Fixed { length: num(arg)? },
        "off" => Closure::Off,
        _ => {
            return Err(format!(
                "--closure {s:?}: expected auto[:FRACTION], ballistic[:COEFF], fixed:LENGTH or off"
            ))
        }
    })
}

fn sea(s: &str) -> Result<Sea, String> {
    let (kind, args) = s.split_once(':').unwrap_or((s, ""));
    let (mut hs, mut tp, mut gamma) = (None, None, None);
    for kv in args.split(',').filter(|kv| !kv.is_empty()) {
        let (k, v) = kv
            .split_once('=')
            .ok_or_else(|| format!("--sea: {kv:?}: expected key=value"))?;
        let v: f64 = v
            .parse()
            .map_err(|_| format!("--sea: {kv:?}: expected a number"))?;
        match k {
            "hs" => hs = Some(v),
            "tp" => tp = Some(v),
            "gamma" => gamma = Some(v),
            _ => return Err(format!("--sea: {k:?}: expected hs, tp or gamma")),
        }
    }
    let (hs, tp) = match (hs, tp) {
        (Some(h), Some(t)) => (h, t),
        _ => return Err(format!("--sea {s:?}: needs hs= and tp=")),
    };
    match kind {
        "bretschneider" | "ittc" if gamma.is_none() => Ok(Sea::Bretschneider { hs, tp }),
        "jonswap" => Ok(Sea::Jonswap { hs, tp, gamma }),
        _ => Err(format!(
            "--sea {s:?}: expected bretschneider:… or jonswap:…"
        )),
    }
}

fn main() {
    let cli = Cli::parse();
    match go(cli) {
        Ok(0) => {}
        Ok(_) => std::process::exit(1),
        Err(e) => {
            eprintln!("boatmath: {e}");
            std::process::exit(1);
        }
    }
}

/// Run a command; returns how many of its records failed.
fn go(cli: Cli) -> Result<usize, String> {
    let cache = Cache::open(cli.cache.clone());
    let read = || Stream::read(&cli.inputs);
    let mut out = Out::new(std::io::stdout());
    let mut failed = 0;
    let fail = |what: &str, id: &str, e: String| {
        eprintln!("boatmath: {what} {}: {e}", short(id));
        1
    };
    match cli.command {
        Command::Hull {
            files,
            name,
            waterline,
            centerplane,
            stations,
            rays,
            units,
        } => {
            if files.is_empty() {
                return Err("hull: name a hull file".into());
            }
            let units = units.as_deref().map(boatmath::parse_units).transpose()?;
            for f in &files {
                for w in list::product(&[list(&waterline, "waterline")?]) {
                    let import = LoftRequest {
                        waterline: w[0],
                        centerplane,
                        stations,
                        rays,
                        units,
                        ..LoftRequest::default()
                    };
                    match records::hull_from_file(f, name.as_deref(), import) {
                        Ok(r) => out.emit(&r)?,
                        Err(e) => {
                            eprintln!("boatmath: {e}");
                            failed += 1;
                        }
                    }
                }
            }
        }
        Command::Scale {
            by,
            beam,
            mass,
            keep_length,
            name,
        } => {
            use records::Scaling;
            let s = read()?;
            out.pass(&s)?;
            let mut how: Vec<Scaling> = Vec::new();
            for k in list(&by, "by")?.unwrap_or_default() {
                how.push(Scaling::By(k));
            }
            for k in list(&beam, "beam")?.unwrap_or_default() {
                how.push(Scaling::Beam(k));
            }
            for kg in list(&mass, "mass")?.unwrap_or_default() {
                how.push(Scaling::Mass { kg, keep_length });
            }
            if how.is_empty() {
                return Err("scale: give --by, --beam or --mass".into());
            }
            if keep_length && mass.is_none() {
                return Err("--keep-length goes with --mass".into());
            }
            for h in s.of_type("hull") {
                for &w in &how {
                    match records::scaled_hull(h, w, name.as_deref()) {
                        Ok(r) => out.emit(&r)?,
                        Err(e) => failed += fail("hull", h["id"].as_str().unwrap_or(""), e),
                    }
                }
            }
        }
        Command::Case {
            span,
            mass,
            lcg,
            vcg,
            kxx,
            kyy,
            kzz,
            roll_damping,
            name,
        } => {
            let s = read()?;
            out.pass(&s)?;
            let axes = [
                list(&span, "span")?,
                list(&mass, "mass")?,
                list(&lcg, "lcg")?,
                list(&vcg, "vcg")?,
                list(&kxx, "kxx")?,
                list(&kyy, "kyy")?,
                list(&kzz, "kzz")?,
                list(&roll_damping, "roll-damping")?,
            ];
            for h in s.of_type("hull") {
                for v in list::product(&axes) {
                    let params = CaseParams {
                        span: v[0],
                        mass: v[1],
                        lcg: v[2],
                        vcg: v[3],
                        kxx: v[4],
                        kyy: v[5],
                        kzz: v[6],
                        roll_damping: v[7].unwrap_or(0.0),
                        ..CaseParams::default()
                    };
                    out.emit(&records::case(h, params, name.as_deref())?)?;
                }
            }
        }
        Command::Mount {
            stock,
            kind,
            x,
            y,
            pair,
            shaft_depth,
            chord,
            thickness,
            shaft_diameter,
            strut_ahead,
            pod_length,
            pod_diameter,
            nose_ahead,
            prop_from_nose,
            no_pod,
            tractor,
            prop_diameter,
            shaft_angle,
            name,
        } => {
            let len = |s: &Option<String>| s.as_deref().map(units::length).transpose();
            let kind = kind
                .map(|k| {
                    serde_json::from_value::<boatmath::mount::Kind>(serde_json::json!(k))
                        .map_err(|_| format!("--kind {k}: saildrive, outboard, pod or shaft"))
                })
                .transpose()?;
            let opts = boatmath::mount::Options {
                kind,
                x: len(&x)?,
                y: units::length(&y)?,
                shaft_depth: len(&shaft_depth)?,
                chord: len(&chord)?,
                thickness: len(&thickness)?,
                pod_length: len(&pod_length)?,
                pod_diameter: len(&pod_diameter)?,
                nose_ahead: len(&nose_ahead)?,
                prop_from_nose: len(&prop_from_nose)?,
                no_pod,
                tractor: tractor.then_some(true),
                prop_diameter: len(&prop_diameter)?,
                shaft_angle_deg: shaft_angle,
                shaft_diameter: len(&shaft_diameter)?,
                strut_ahead: len(&strut_ahead)?,
                pair,
            };
            let mount = boatmath::mount::build(stock.as_deref(), &opts)?;
            let s = read()?;
            out.pass(&s)?;
            for c in s.of_type("case") {
                match records::mounted_case(c, &mount, name.as_deref()) {
                    Ok(r) => out.emit(&r)?,
                    Err(e) => failed += fail("case", c["id"].as_str().unwrap_or(""), e),
                }
            }
        }
        Command::Mounts => {
            for m in boatmath::mount::stock_table()["mounts"]
                .as_array()
                .into_iter()
                .flatten()
            {
                let mut r = m.clone();
                if let Some(o) = r.as_object_mut() {
                    o.insert("type".into(), serde_json::json!("stock_mount"));
                    o.insert("id".into(), m["name"].clone());
                }
                out.emit(&r)?;
            }
        }
        Command::Statics => {
            let s = read()?;
            out.pass(&s)?;
            for h in s.of_type("hull") {
                match records::hull_statics(h, &cache) {
                    Ok(r) => out.emit(&r)?,
                    Err(e) => failed += fail("hull", h["id"].as_str().unwrap_or(""), e),
                }
            }
            for c in s.of_type("case") {
                let id = c["id"].as_str().unwrap_or("");
                match records::case_statics(&s, c, &cache) {
                    Ok(rs) => {
                        for r in &rs {
                            if let Some(e) = r["error"].as_str() {
                                failed += fail("case", id, e.to_string());
                            }
                            out.emit(r)?;
                        }
                    }
                    Err(e) => failed += fail("case", id, e),
                }
            }
        }
        Command::GzWaves { length, height } => {
            let lengths = parse_lengths(&length, "length")?;
            let heights = parse_lengths(&height, "height")?;
            let s = read()?;
            out.pass(&s)?;
            for c in s.of_type("case") {
                let id = c["id"].as_str().unwrap_or("");
                for &l in &lengths {
                    for &h in &heights {
                        match records::gz_waves(&s, c, l, h, &cache) {
                            Ok(r) => {
                                if let Some(e) = r["error"].as_str() {
                                    failed += fail("case", id, e.to_string());
                                }
                                out.emit(&r)?;
                            }
                            Err(e) => failed += fail("case", id, e),
                        }
                    }
                }
            }
        }
        Command::Study {
            froude,
            hold,
            closure: cl,
            grid,
            astern,
            waves,
            lambdas,
            sea: sea_arg,
        } => {
            let froudes = list::parse(&froude).map_err(|e| format!("--froude: {e}"))?;
            let closure = cl.as_deref().map(closure).transpose()?.unwrap_or_default();
            let lambdas = list(&lambdas, "lambdas")?;
            let sea = sea_arg.as_deref().map(sea).transpose()?;
            if waves.is_none() && (lambdas.is_some() || sea.is_some()) {
                return Err("--lambdas and --sea need --waves".into());
            }
            let headings = list(&waves, "waves")?;
            let s = read()?;
            out.pass(&s)?;
            for c in s.of_type("case") {
                for v in list::product(&[headings.clone(), Some(froudes.clone())]) {
                    let params = StudyParams {
                        froude: v[1].expect("a froude number"),
                        dynamic: !hold,
                        closure,
                        grid: grid.unwrap_or_else(boatmath::params::default_grid),
                        astern,
                        waves: v[0].map(|heading| Waves {
                            heading,
                            lambdas: lambdas.clone(),
                            sea,
                        }),
                    };
                    out.emit(&records::study(c, params)?)?;
                }
            }
        }
        Command::Run { jobs, quiet } => {
            let jobs =
                jobs.unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
            let s = read()?;
            failed += run::run(&s, &mut out, &cache, &run::Options { jobs, quiet })?;
        }
        Command::Prop {
            thrust,
            speed,
            shafts,
            wake,
            thrust_deduction,
            d_max,
            d_min,
            blades,
            depth,
            keller_k,
            no_cavitation,
            no_re_correct,
            strict_ear,
            top_froude,
            top_speed,
            top_thrust,
            full_range,
        } => {
            let blades = list::parse(&blades)
                .map_err(|e| format!("--blades: {e}"))?
                .into_iter()
                .map(|z| {
                    let z = z.round();
                    if (2.0..=7.0).contains(&z) {
                        Ok(z as u32)
                    } else {
                        Err(format!("--blades: {z}: the series has 2 to 7 blades"))
                    }
                })
                .collect::<Result<Vec<u32>, String>>()?;
            let coefficient = |s: &str, what: &str| -> Result<Option<f64>, String> {
                if s == "auto" {
                    return Ok(None);
                }
                match s.parse::<f64>() {
                    Ok(v) if (0.0..1.0).contains(&v) => Ok(Some(v)),
                    _ => Err(format!("--{what} {s:?}: a fraction, 0 to 1, or auto")),
                }
            };
            let (wake_given, ded_given) = (
                coefficient(&wake, "wake")?,
                coefficient(&thrust_deduction, "thrust-deduction")?,
            );
            let d_min = units::length(&d_min)?;
            let d_max = d_max.as_deref().map(units::length).transpose()?;
            let depth = depth.as_deref().map(units::length).transpose()?;
            let auto = (wake_given.is_none() || ded_given.is_none()).then(|| props::Auto {
                wake: wake_given.is_none(),
                deduction: ded_given.is_none(),
            });
            let (wake, thrust_deduction) = (wake_given.unwrap_or(0.0), ded_given.unwrap_or(0.0));
            let explicit_top = match (&top_speed, &top_thrust) {
                (Some(s), Some(t)) => Some(propeller::bseries::TopInputs {
                    speed: units::speed(s)?,
                    thrust: units::force(t)?,
                }),
                (None, None) => None,
                _ => return Err("--top-speed and --top-thrust go together".into()),
            };
            let base = |speed: f64, thrust: f64, shafts: u32, top, d_max: f64, depth: f64| {
                propeller::bseries::Inputs {
                    speed,
                    thrust,
                    shafts,
                    wake,
                    thrust_deduction,
                    d_min,
                    d_max,
                    blades: blades.clone(),
                    depth,
                    keller_k,
                    cavitation: !no_cavitation,
                    re_correct: !no_re_correct,
                    strict_ear,
                    top,
                }
            };
            let cap = if full_range { f64::INFINITY } else { 1.6 };
            let report = |r: &Value| {
                let id = r["id"].as_str().unwrap_or("");
                if let Some(i) = r.get("interaction").filter(|i| !i.is_null()) {
                    // What was used: a coefficient given outright is not the
                    // one the interaction computed alongside.
                    let used = |k: &str| r["inputs"][k].as_f64().unwrap_or(0.0);
                    let w = if i["settings"]["wake"] == true {
                        format!(
                            "w {:.3} (local {:.3}, wave {:.3})",
                            used("wake"),
                            i["w_local"].as_f64().unwrap_or(0.0),
                            i["w_wave"].as_f64().unwrap_or(0.0)
                        )
                    } else {
                        format!("w {:.3} given", used("wake"))
                    };
                    let t = if i["settings"]["thrust_deduction"] == true {
                        format!("t {:.3}", used("thrust_deduction"))
                    } else {
                        format!("t {:.3} given", used("thrust_deduction"))
                    };
                    eprintln!(
                        "boatmath: prop {}: {w}, {t}, after {} steps{}",
                        short(id),
                        i["steps"].as_array().map_or(0, Vec::len),
                        if i["converged"] == true {
                            ""
                        } else {
                            " (not converged)"
                        },
                    );
                }
                let b = &r["best"];
                if b.is_null() {
                    eprintln!(
                        "boatmath: prop {}: no B-series propeller can do it (try a larger --d-max, or --no-cavitation)",
                        short(id)
                    );
                    return 1;
                }
                eprintln!(
                    "boatmath: prop {}: {:.0} W at {:.0} rpm, {}-blade, D {:.0} mm, P/D {:.2}, EAR {:.2}, η₀ {:.3}",
                    short(id),
                    b["P_shaft"].as_f64().unwrap_or(0.0),
                    b["rpm"].as_f64().unwrap_or(0.0),
                    b["Z"].as_f64().unwrap_or(0.0),
                    1e3 * b["D"].as_f64().unwrap_or(0.0),
                    b["PD"].as_f64().unwrap_or(0.0),
                    b["EAR"].as_f64().unwrap_or(0.0),
                    b["eta0"].as_f64().unwrap_or(0.0),
                );
                0
            };
            match (&thrust, &speed) {
                (Some(t), Some(v)) => {
                    if top_froude.is_some() {
                        return Err(
                            "--top-froude needs results in the stream; give --top-speed and --top-thrust"
                                .into(),
                        );
                    }
                    if auto.is_some() {
                        return Err(
                            "--wake auto and --thrust-deduction auto need results in the stream"
                                .into(),
                        );
                    }
                    let inputs = base(
                        units::speed(v)?,
                        units::force(t)?,
                        shafts.unwrap_or(1),
                        explicit_top,
                        d_max.ok_or("give --d-max")?,
                        depth.ok_or("give --depth (the shaft's immersion)")?,
                    );
                    let r = props::prop(&cache, &inputs, None, cap, None, None)?;
                    failed += report(&r);
                    out.emit(&r)?;
                }
                (None, None) => {
                    let s = read()?;
                    out.pass(&s)?;
                    for r in s.of_type("result").filter(|r| r["kind"] == "calm") {
                        let id = r["id"].as_str().unwrap_or("");
                        let one = (|| -> Result<Value, String> {
                            // The case's drive, if it has one: where the
                            // propellers sit, and the thrust line's angle.
                            let drives = props::drives_of(&s, r)?;
                            let cos_e = drives.as_ref().map_or(1.0, |d| d.angle.cos());
                            let (v, t) = props::operating_point(r, thrust_deduction)?;
                            let top = match top_froude {
                                Some(f) => {
                                    let o = props::result_at_froude(&s, r, f)?;
                                    let (vt, tt) = props::operating_point(&o, thrust_deduction)?;
                                    Some(propeller::bseries::TopInputs {
                                        speed: vt,
                                        thrust: tt / cos_e,
                                    })
                                }
                                None => explicit_top,
                            };
                            let n = shafts.unwrap_or_else(|| match &drives {
                                Some(d) => d.at.len() as u32,
                                None => props::hull_count(&s, r) as u32,
                            });
                            let d_max = d_max
                                .or(drives.as_ref().and_then(|d| d.prop_diameter))
                                .ok_or("give --d-max (the drive has no propeller of its own)")?;
                            let depth = match (depth, &drives) {
                                (Some(d), _) => d,
                                (None, Some(d)) => {
                                    d.at.iter().map(|a| a.2).sum::<f64>() / d.at.len().max(1) as f64
                                }
                                (None, None) => {
                                    return Err("give --depth, or give the case a mount".into())
                                }
                            };
                            let inputs = base(v, t / cos_e, n, top, d_max, depth);
                            match &auto {
                                Some(a) => props::auto_prop(&s, &cache, r, inputs, a, cap),
                                None => props::prop(
                                    &cache,
                                    &inputs,
                                    Some(id),
                                    cap,
                                    None,
                                    drives.as_ref().map(|d| &d.record),
                                ),
                            }
                        })();
                        match one {
                            Ok(p) => {
                                failed += report(&p);
                                out.emit(&p)?;
                            }
                            Err(e) => failed += fail("result", id, e),
                        }
                    }
                }
                _ => return Err("--thrust and --speed go together".into()),
            }
        }
        Command::Match {
            rank_by,
            all_windings,
            direct_only,
            ratio,
            ratio_max,
            gear_eta,
            no_controller_loss,
            peak,
            no_top,
            filter,
            motors,
        } => {
            use propeller::motor::RankBy;
            let rank_by = match rank_by.as_str() {
                "power" => RankBy::Power,
                "mass" => RankBy::Mass,
                "price" => RankBy::Price,
                o => return Err(format!("--rank-by {o:?}: expected power, mass or price")),
            };
            let owned;
            let db = match &motors {
                Some(p) => {
                    owned = motor_database(p)?;
                    &owned
                }
                None => propeller::motor::Database::vendored(),
            };
            let o = props::MatchOptions {
                filter: filter.filter()?,
                rank_by,
                all_windings,
                ratio: if direct_only { Some(1.0) } else { ratio },
                ratio_max,
                gear_eta,
                controller_loss: !no_controller_loss,
                peak,
                check_top: !no_top,
            };
            let s = read()?;
            out.pass(&s)?;
            for p in s.of_type("prop") {
                let (drives, unreachable) = props::drives(db, p, &o)?;
                eprintln!(
                    "boatmath: prop {}: {} motors can drive it, {} can't",
                    short(p["id"].as_str().unwrap_or("")),
                    drives.len(),
                    unreachable.len()
                );
                for d in &drives {
                    out.emit(d)?;
                }
            }
        }
        Command::Motors { filter } => {
            let f = filter.filter()?;
            for m in &propeller::motor::Database::vendored().motors {
                if f.admits(m) {
                    out.emit(&props::motor_record(m))?;
                }
            }
        }
        Command::Pick {
            kind,
            first,
            wheres,
            drop,
        } => {
            let s = read()?;
            let p = pick::Pick {
                kind,
                first,
                wheres: wheres
                    .iter()
                    .map(|w| pick::parse_where(w))
                    .collect::<Result<_, _>>()?,
                drop,
            };
            for i in pick::pick(&s, &p) {
                out.emit(&s.records[i])?;
            }
        }
        Command::Table {
            fields,
            kind,
            explode,
            csv,
            no_header,
        } => {
            let s = read()?;
            let res = path::Resolver::new(&s);
            let line = |cells: Vec<String>| {
                if csv {
                    cells
                        .iter()
                        .map(|c| csv_cell(c))
                        .collect::<Vec<_>>()
                        .join(",")
                } else {
                    cells
                        .iter()
                        .map(|c| c.replace(['\t', '\n'], " "))
                        .collect::<Vec<_>>()
                        .join("\t")
                }
            };
            let w = out.raw();
            if !no_header {
                writeln!(w, "{}", line(fields.clone())).map_err(stdout_err)?;
            }
            for r in rows_of(&s, &kind, &fields[0], explode.as_deref()) {
                for row in res.rows(r, explode.as_deref()) {
                    let cells = fields
                        .iter()
                        .map(|f| path::cell(&res.get_in(&row, r, f)))
                        .collect();
                    writeln!(w, "{}", line(cells)).map_err(stdout_err)?;
                }
            }
        }
        Command::Plot {
            x,
            y,
            by,
            kind,
            explode,
            title,
            xlabel,
            ylabel,
            o,
        } => {
            let s = read()?;
            let res = path::Resolver::new(&s);
            let mut series: Vec<plot::Series> = Vec::new();
            for r in rows_of(&s, &kind, &x, explode.as_deref()) {
                for row in res.rows(r, explode.as_deref()) {
                    let Some(xv) = res.get_in(&row, r, &x).and_then(|v| v.as_f64()) else {
                        continue;
                    };
                    let group: Vec<String> = by
                        .iter()
                        .map(|f| {
                            format!("{} {}", short_field(f), path::cell(&res.get_in(&row, r, f)))
                        })
                        .collect();
                    for yf in &y {
                        let Some(yv) = res.get_in(&row, r, yf).and_then(|v| v.as_f64()) else {
                            continue;
                        };
                        let mut name = group.clone();
                        if y.len() > 1 {
                            name.insert(0, short_field(yf));
                        }
                        let name = name.join(", ");
                        match series.iter_mut().find(|s| s.name == name) {
                            Some(s) => s.points.push((xv, yv)),
                            None => series.push(plot::Series {
                                name,
                                points: vec![(xv, yv)],
                            }),
                        }
                    }
                }
            }
            let svg = plot::svg(&plot::Plot {
                title,
                x_label: xlabel.unwrap_or_else(|| x.clone()),
                y_label: ylabel.unwrap_or_else(|| y.join(", ")),
                series,
            })?;
            match o {
                Some(p) => std::fs::write(&p, svg).map_err(|e| format!("{}: {e}", p.display()))?,
                None => out.raw().write_all(svg.as_bytes()).map_err(stdout_err)?,
            }
        }
        Command::Wake { o, range, title } => {
            failed += draw(&read()?, out.raw(), View::Wake(range), o, title)?;
        }
        Command::Pressure { o, range, title } => {
            failed += draw(&read()?, out.raw(), View::Pressure(range), o, title)?;
        }
        Command::Profile {
            o,
            wave_scale,
            title,
        } => {
            failed += draw(&read()?, out.raw(), View::Profile(wave_scale), o, title)?;
        }
        Command::Cad {
            heel,
            o,
            no_water,
            wave_scale,
            prop_discs,
            left_handed,
        } => {
            let opts = cad::Options {
                water: !no_water,
                wave_scale,
                blades: !prop_discs,
                right_handed: !left_handed,
            };
            let s = read()?;
            // A prop draws its result too; without props, each result.
            let mut targets: Vec<&Value> = s
                .of_type("prop")
                .filter(|p| !p["result"].is_null())
                .collect();
            if targets.is_empty() {
                targets = s
                    .of_type("result")
                    .filter(|r| r["kind"] == "calm")
                    .collect();
            }
            let heels = list(&heel, "heel")?.unwrap_or_default();
            targets.extend(s.of_type("statics"));
            if targets.is_empty() {
                return Err("no calm-water results, props on them or statics in the stream".into());
            }
            check_pattern(targets.len(), &o, "records")?;
            for r in targets {
                let id = r["id"].as_str().unwrap_or("");
                let model = if r["type"] == "statics" {
                    cad::statics_model(&s, r, &heels)
                } else {
                    cad::model(&s, r, &opts)
                };
                match model {
                    Ok(text) => match &o {
                        Some(p) => {
                            let froude = r["froude"]
                                .as_f64()
                                .or_else(|| {
                                    r["result"]
                                        .as_str()
                                        .and_then(|rid| s.get("result", rid))
                                        .and_then(|res| res["froude"].as_f64())
                                })
                                .unwrap_or(0.0);
                            let path = named(p, id, froude);
                            std::fs::write(&path, text).map_err(|e| format!("{path}: {e}"))?;
                            eprintln!("boatmath: wrote {path}");
                        }
                        None => out.raw().write_all(text.as_bytes()).map_err(stdout_err)?,
                    },
                    Err(e) => failed += fail(r["type"].as_str().unwrap_or("record"), id, e),
                }
            }
        }
    }
    Ok(failed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_lengths_parses_units_and_ranges() {
        assert_eq!(parse_lengths("10", "length").unwrap(), vec![10.0]);
        assert_eq!(parse_lengths("10m", "length").unwrap(), vec![10.0]);
        assert_eq!(parse_lengths("500mm", "length").unwrap(), vec![0.5]);
        assert_eq!(
            parse_lengths("5, 10, 15", "length").unwrap(),
            vec![5.0, 10.0, 15.0]
        );
        assert_eq!(
            parse_lengths("5m, 10m, 15m", "length").unwrap(),
            vec![5.0, 10.0, 15.0]
        );
        assert_eq!(
            parse_lengths("5:15:5", "length").unwrap(),
            vec![5.0, 10.0, 15.0]
        );
        assert_eq!(
            parse_lengths("5m:15m:5m", "length").unwrap(),
            vec![5.0, 10.0, 15.0]
        );
    }

    #[test]
    fn gz_waves_on_wigley_case() {
        let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 1.0).unwrap();
        let iges_text = hullgeom::iges::write(&surfaces, "wigley").unwrap();
        let dir = std::env::temp_dir();
        let path = dir.join("test_wigley_gz_waves.igs");
        std::fs::write(&path, iges_text).unwrap();

        let req = LoftRequest {
            waterline: Some(-0.375),
            ..LoftRequest::default()
        };
        let hull = records::hull_from_file(&path, Some("wigley"), req).unwrap();
        let cp = CaseParams {
            vcg: Some(-0.3),
            ..CaseParams::default()
        };
        let case = records::case(&hull, cp, Some("c1")).unwrap();
        let mut s = Stream::default();
        s.push(hull.clone());
        s.push(case.clone());

        let cache = Cache::default();
        let gz = records::gz_waves(&s, &case, 10.0, 0.2, &cache).unwrap();
        assert!(gz["error"].is_null(), "{}", gz["error"]);
        assert_eq!(gz["type"], "gz_waves");
        assert_eq!(gz["case"], case["id"]);
        assert_eq!(gz["length"], 10.0);
        assert_eq!(gz["height"], 0.2);
        assert_eq!(gz["solver_version"], boatmath::SOLVER_VERSION);

        let headings = gz["headings"].as_array().unwrap();
        assert_eq!(headings.len(), 3);
        assert_eq!(headings[0]["name"], "beam");
        assert_eq!(headings[1]["name"], "stern quartering");
        assert_eq!(headings[2]["name"], "bow quartering");
        assert!(gz["calm"]["max_gz"].as_f64().unwrap() > 0.0);

        s.push(gz.clone());

        // Ancestor link test
        let parents = stream::parents(&gz);
        assert!(parents
            .iter()
            .any(|(t, id)| t == "case" && id == case["id"].as_str().unwrap()));

        // Table / resolver test: explode headings
        let res = path::Resolver::new(&s);
        let rows = rows_of(&s, &None, "name", Some("headings"));
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0]["id"], gz["id"]);
        let exploded = res.rows(rows[0], Some("headings"));
        assert_eq!(exploded.len(), 3);
        for row in &exploded {
            assert!(res.get_in(row, rows[0], "name").is_some());
            assert_eq!(res.get_in(row, rows[0], "length").unwrap(), 10.0);
            assert_eq!(res.get_in(row, rows[0], "height").unwrap(), 0.2);
            assert!(res.get_in(row, rows[0], "case.params.mass").is_some());
        }

        // Plot / resolver test: explode headings for heading_deg vs worst.max_gz
        let p_rows = rows_of(&s, &None, "heading_deg", Some("headings"));
        assert_eq!(p_rows.len(), 1);
        let p_exploded = res.rows(p_rows[0], Some("headings"));
        let mut series = Vec::new();
        for row in &p_exploded {
            let x = res
                .get_in(row, p_rows[0], "heading_deg")
                .and_then(|v| v.as_f64())
                .unwrap();
            let y = res
                .get_in(row, p_rows[0], "worst.max_gz")
                .and_then(|v| v.as_f64())
                .unwrap();
            series.push((x, y));
        }
        assert_eq!(series.len(), 3);
        let svg = plot::svg(&plot::Plot {
            title: Some("Wigley GZ in Waves".into()),
            x_label: "Heading [deg]".into(),
            y_label: "Max GZ [m]".into(),
            series: vec![plot::Series {
                name: "worst".into(),
                points: series,
            }],
        })
        .unwrap();
        assert!(svg.contains("<svg"));

        // Clean up
        let _ = std::fs::remove_file(path);
    }
}
