//! `boatmath` — hulls, cases and studies on the `thinship`, `seakeeping` and `hullgeom` solvers: the
//! computations the web app (`boatmath-web`) and the CLI share.
//!
//! A hull is seen the way the physics sees it: cut into sections.
//! IGES and STL hulls are cut straight from their patches or triangles —
//! read by [`files`]. The page is sent each
//! station's section curve (what the depth integral integrates), the CAD
//! ray hits it was interpolated from, the
//! depth-integral curve the kernel interpolates along x, the hydrostatics,
//! and the transom with what the closure needs to draw its virtual appendage
//! at any speed.

use crate::files::{open_source_bytes, Kind, LoadSettings};
use hullgeom::iges::{HullPose, Platform, SectionalImport};
use hullgeom::SectionalHull;
use serde_json::{json, Value};
use thinship::{Conditions, Placement, WaveOptions, STANDARD_GRAVITY};

pub mod cad;
pub mod camber;
pub mod files;
pub mod mount;
pub mod native;
pub mod params;
pub mod platform;
pub mod propulsion;
pub mod sections;

/// The diverging colour map the wake and pressure views share.
pub use files::diverging;
pub use files::parse_units;
use params::{CaseParams, StudyParams};
pub use platform::{at_rest, heeled, statics, wave_gz, waves_with_progress};

/// The solver version results are stamped with: the last commit to touch
/// the solver's code (the geometry, thin-ship, seakeeping and CLI crates,
/// and this one), `-dirty` if they have uncommitted changes. Set by
/// `build.rs`; `BOATMATH_SOLVER_VERSION` at build time overrides it.
pub const SOLVER_VERSION: &str = env!("BOATMATH_SOLVER_VERSION");

/// Largest upload accepted [bytes]. Big enough for a finely tessellated STL.
pub const MAX_UPLOAD: usize = 128 << 20;

/// Import options the page can set.
#[derive(Default, Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LoftRequest {
    /// Design waterline height in the file's frame [m] (IGES).
    pub waterline: Option<f64>,
    /// Centreplane override [m] (IGES).
    pub centerplane: Option<f64>,
    /// Stations along the hull (IGES).
    pub stations: Option<usize>,
    /// Rays across each section (IGES).
    pub rays: Option<usize>,
    /// Scale to metres (STL, which carries no units): mm, m, in, ... or a
    /// number.
    pub units: Option<f64>,
    /// A hull made by scaling another: the whole hull by this factor, about
    /// its design waterline (length, beam and draft alike) …
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale: Option<f64>,
    /// … and/or its beam and draft only, the length kept.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale_yz: Option<f64>,
}

impl LoftRequest {
    /// The design pose the settings put each hull in: its scale, if any.
    pub fn pose(&self) -> HullPose {
        HullPose {
            scale: self.scale.unwrap_or(1.0),
            scale_yz: self.scale_yz.unwrap_or(1.0),
            ..HullPose::default()
        }
    }
}

impl LoftRequest {
    /// Parse from `key=value` query pairs; unknown keys are ignored.
    pub fn from_query(pairs: &[(String, String)]) -> Result<LoftRequest, String> {
        let mut r = LoftRequest::default();
        let num = |k: &str, v: &str| {
            v.trim()
                .parse::<f64>()
                .map_err(|_| format!("{k}: expected a number, got {v:?}"))
        };
        let count = |k: &str, v: &str, min: usize| match v.trim().parse::<usize>() {
            Ok(n) if n >= min => Ok(n),
            _ => Err(format!(
                "{k}: expected a count of at least {min}, got {v:?}"
            )),
        };
        for (k, v) in pairs {
            if v.trim().is_empty() {
                continue;
            }
            match k.as_str() {
                "waterline" => r.waterline = Some(num(k, v)?),
                "centerplane" => r.centerplane = Some(num(k, v)?),
                "stations" => r.stations = Some(count(k, v, 8)?),
                "rays" => r.rays = Some(count(k, v, 5)?),
                "units" => r.units = Some(files::parse_units(v.trim())?),
                "scale" => r.scale = Some(num(k, v)?).filter(|&k| k != 1.0),
                "scale_yz" => r.scale_yz = Some(num(k, v)?).filter(|&k| k != 1.0),
                _ => {}
            }
        }
        Ok(r)
    }
}

/// Cut an uploaded hull into sections and describe it as the JSON the page
/// draws.
pub fn loft(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Value, String> {
    let t0 = std::time::Instant::now();
    let cut = cut(name, bytes, req)?;
    let hulls: Vec<Value> = cut
        .hulls
        .iter()
        .map(|h| sectioned_json(h, cut.kind))
        .collect();
    Ok(json!({
        "name": name,
        "seconds": t0.elapsed().as_secs_f64(),
        "notes": cut.notes,
        "hulls": hulls,
    }))
}

/// The whole of an upload's geometry for display, topsides and all, in the
/// file's own frame (z up, metres): each hull's tessellation as base64
/// little-endian arrays, and the geometry's z range — what the waterline
/// can be set within. Its units and scale are the settings'; a scale is
/// about the design waterline.
pub fn geometry(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Value, String> {
    let settings = LoadSettings {
        units: req.units,
        ..LoadSettings::default()
    };
    let file = if native::is_native(&bytes) {
        native::open(&bytes)?
    } else {
        files::open_geometry(name, bytes, &settings)?
    };
    let wl = req.waterline.unwrap_or(0.0);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut meshes = Vec::new();
    for i in 0..file.source.len() {
        let (mut v, t) = file
            .source
            .posed_tessellation(i, wl, &req.pose(), &Platform::default())
            .map_err(|e| e.to_string())?;
        // Posed about the waterline, z up from it: back into the file's frame.
        for p in v.iter_mut() {
            p[2] += wl;
            lo = lo.min(p[2]);
            hi = hi.max(p[2]);
        }
        let flat: Vec<f32> = v.iter().flat_map(|p| p.map(|c| c as f32)).collect();
        let idx: Vec<u32> = t.iter().flatten().copied().collect();
        meshes.push(json!({
            "vertices": b64(bytemuck_f32(&flat)),
            "triangles": b64(bytemuck_u32(&idx)),
        }));
    }
    Ok(json!({ "z_range": [lo, hi], "meshes": meshes }))
}

/// An upload cut into sections at its design pose, with the source kept so
/// its hulls can be re-cut at another attitude.
struct Cut {
    kind: Kind,
    file: files::SourceFile,
    /// Each cut hull's index in the file.
    index: Vec<usize>,
    hulls: Vec<SectionalImport>,
    opts: hullgeom::iges::SectionalOptions,
    notes: Vec<String>,
}

fn cut(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Cut, String> {
    let d = LoadSettings::default();
    let settings = LoadSettings {
        waterline_z: req.waterline.unwrap_or(0.0),
        centerplane: req.centerplane,
        stations: req.stations.unwrap_or(d.stations),
        rays: req.rays.unwrap_or(d.rays),
        units: req.units,
    };
    let file = if native::is_native(&bytes) {
        native::open(&bytes)?
    } else {
        open_source_bytes(name, bytes, &settings)?
    };
    let opts = settings.sectional(file.waterline_z);
    let mut notes = Vec::new();
    let mut hulls = Vec::new();
    let mut index = Vec::new();
    for i in 0..file.source.len() {
        match file.source.situate_sectional(
            i,
            file.waterline_z,
            &req.pose(),
            &Platform::default(),
            &opts,
        ) {
            Ok(Some(h)) => {
                hulls.push(h);
                index.push(i);
            }
            Ok(None) => notes.push(format!("hull {} is dry at this waterline", i + 1)),
            Err(e) => notes.push(format!("hull {} not sectioned: {e}", i + 1)),
        }
    }
    if hulls.is_empty() {
        return Err(notes
            .first()
            .cloned()
            .unwrap_or_else(|| "no hull found".into()));
    }
    Ok(Cut {
        kind: file.kind,
        file,
        index,
        hulls,
        opts,
        notes,
    })
}

/// The solver's dynamic-load model with a progress report on every full
/// evaluation (the solver's cheap slope probes pass straight through), so a
/// cancelled request also stops the Newton loop.
struct Counted<'a, R> {
    inner: thinship::sectional::SectionalDynamic<'a>,
    report: R,
    /// A load added to every evaluation: the drives' thrust.
    extra: hullgeom::float::DynamicLoad,
}

impl<R> hullgeom::float::DynamicModel<SectionalHull> for Counted<'_, R>
where
    R: FnMut(&hullgeom::float::DynamicLoad, f64) -> thinship::Result<()>,
{
    fn load(
        &mut self,
        fleet: &hullgeom::float::FleetState<SectionalHull>,
    ) -> thinship::Result<hullgeom::float::DynamicLoad> {
        let mut d = self.inner.load(fleet)?;
        d.force_up += self.extra.force_up;
        d.moment_bow_up += self.extra.moment_bow_up;
        let vol: f64 = fleet
            .members
            .iter()
            .map(|(h, _)| h.displaced_volume())
            .sum();
        (self.report)(&d, vol)?;
        Ok(d)
    }

    fn probe(
        &mut self,
        fleet: &hullgeom::float::FleetState<SectionalHull>,
    ) -> Option<thinship::Result<hullgeom::float::DynamicLoad>> {
        let extra = self.extra;
        self.inner.probe(fleet).map(|r| {
            r.map(|mut d| {
                d.force_up += extra.force_up;
                d.moment_bow_up += extra.moment_bow_up;
                d
            })
        })
    }
}

/// How a case's mass is carried (see [`FlowRequest::mass_by`]).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default, serde::Serialize, serde::Deserialize)]
pub enum MassBy {
    /// The hull floats deeper or shallower: the waterline moves.
    #[default]
    #[serde(rename = "sinking")]
    Sinking,
    /// The hull is scaled uniformly, `k = (m/m₀)^(1/3)`, about its design
    /// waterline: the waterline stays.
    #[serde(rename = "scale")]
    ScaleXyz,
    /// Beam and draft are scaled, `k = (m/m₀)^(1/2)`, the length kept: the
    /// waterline stays.
    #[serde(rename = "scale_yz")]
    ScaleYz,
}

/// A calm-water study on a case: the hull's cut, the case
/// (layout and load), and the study's speed and model settings.
pub struct FlowRequest {
    pub cut: LoftRequest,
    pub case: CaseParams,
    pub study: StudyParams,
    /// A previous solution `(sinkage [m], trim [rad])` to start the
    /// equilibrium from — a neighbouring speed's.
    pub warm: Option<(f64, f64)>,
    /// Hold the platform at this `(sinkage [m], trim [rad])` instead of
    /// solving for it: a study at rest's attitude (see [`statics`]).
    pub hold: Option<(f64, f64)>,
    /// The drives' thrust, for a self-propelled attitude: it enters the
    /// equilibrium with the hull's own forces.
    pub thrust: Option<ThrustLine>,
    /// Compute the hull pressure and the free surface; without, only the
    /// attitude and the forces.
    pub field: bool,
}

/// The drives' thrust on a platform: what a self-propelled attitude adds to
/// a towed one.
#[derive(Clone, Debug)]
pub struct ThrustLine {
    /// The total thrust along the shafts [N], shared equally.
    pub thrust: f64,
    /// Each drive's propeller centre, `(x, z)` in the hull's design frame
    /// (z up from the design waterline) [m].
    pub at: Vec<(f64, f64)>,
    /// The shafts' angle to the baseline, bow up [rad].
    pub angle: f64,
}

impl ThrustLine {
    /// Its vertical force [N] and bow-up moment [N m] about `pivot_x`, with
    /// the resistance it balances acting along `z_resistance`: the couple of
    /// the thrust below the resistance trims the bow up; an inclined shaft's
    /// vertical share lifts where the propeller is.
    pub fn load(&self, pivot_x: f64, z_resistance: f64) -> hullgeom::float::DynamicLoad {
        let n = self.at.len().max(1) as f64;
        let (s, c) = self.angle.sin_cos();
        let mut d = hullgeom::float::DynamicLoad {
            force_up: 0.0,
            moment_bow_up: 0.0,
        };
        for &(x, z) in &self.at {
            let t = self.thrust / n;
            d.force_up += t * s;
            d.moment_bow_up += t * c * (z_resistance - z) + t * s * (x - pivot_x);
        }
        d
    }
}

/// The height (z up from the design waterline) the resistance is taken to
/// act along: the centroid of the hulls' immersed centreplane area.
fn resistance_line(design: &[(SectionalHull, Placement)], roles: &[(String, Option<f64>)]) -> f64 {
    let (mut a, mut m) = (0.0, 0.0);
    for ((h, _), (role, _)) in design.iter().zip(roles) {
        if role != "hull" {
            continue;
        }
        let st: Vec<(f64, f64)> = h
            .curves()
            .map(|(x, c)| (x, c.iter().map(|p| p.1).fold(0.0, f64::max)))
            .collect();
        for w in st.windows(2) {
            let (dx, d) = (w[1].0 - w[0].0, 0.5 * (w[0].1 + w[1].1));
            a += d * dx;
            m += -0.5 * d * d * dx;
        }
    }
    if a > 0.0 {
        m / a
    } else {
        0.0
    }
}

/// Where a computation has got to.
#[derive(Debug, Clone)]
pub struct Progress {
    /// The stage under way, e.g. "Floating at speed".
    pub stage: &'static str,
    /// Stage number (1-based) and how many there are.
    pub step: usize,
    pub steps: usize,
    /// What the stage is doing, e.g. the force evaluation count.
    pub detail: String,
    /// The whole computation's estimated fraction done, `0..=1`. Stages are
    /// weighted by their typical cost; the equilibrium's iteration count is
    /// not known in advance, so its share approaches its end asymptotically.
    pub fraction: f64,
}

/// The stages of a flow and their typical shares of its time: floating the
/// hull at speed dominates when it is solved; otherwise the free surface
/// and the dynamic force (then computed at the forces stage) do.
fn flow_stages(dynamic: bool) -> Vec<(&'static str, f64)> {
    // Measured on e12 (Fn 0.3, dynamic: 0.04 / 8.4 / 0.6 / 3.6 / 7.2 s):
    // the forces stage also encodes the hull meshes into the answer.
    if dynamic {
        vec![
            ("Cutting sections", 0.01),
            ("Floating at speed", 0.5),
            ("Hull pressure", 0.04),
            ("Free surface", 0.2),
            ("Forces", 0.25),
        ]
    } else {
        vec![
            ("Cutting sections", 0.02),
            ("At rest's attitude", 0.0),
            ("Hull pressure", 0.08),
            ("Free surface", 0.35),
            ("Forces", 0.55),
        ]
    }
}

/// Reports progress through the caller's callback, and turns its refusal
/// (the study was cancelled) into an error that stops the computation.
pub(crate) struct Tracker<'a> {
    pub(crate) stages: Vec<(&'static str, f64)>,
    pub(crate) report: &'a mut dyn FnMut(&Progress) -> bool,
}

/// The error a computation stops with when its progress callback asks it to.
pub const CANCELLED: &str = "cancelled";

impl Tracker<'_> {
    /// Stage `i` (0-based) is `within` (0..1) of the way through.
    pub(crate) fn at(&mut self, i: usize, within: f64, detail: String) -> Result<(), String> {
        let total: f64 = self.stages.iter().map(|s| s.1).sum::<f64>().max(1e-12);
        let before: f64 = self.stages[..i].iter().map(|s| s.1).sum();
        let p = Progress {
            stage: self.stages[i].0,
            step: i + 1,
            steps: self.stages.len(),
            detail,
            fraction: ((before + self.stages[i].1 * within.clamp(0.0, 1.0)) / total).min(1.0),
        };
        if (self.report)(&p) {
            Ok(())
        } else {
            Err(CANCELLED.into())
        }
    }
}

/// The steady flow at one speed: each hull's near-field pressure, the free
/// surface around the fleet (local field and waves), and the forces.
pub fn flow(name: &str, bytes: Vec<u8>, req: &FlowRequest) -> Result<Value, String> {
    flow_with_progress(name, bytes, req, &mut |_| true)
}

/// A case set up on its hull's cut: the platform's hulls as
/// `(hull index, pose)` — the file's hulls where they are, or its one hull
/// doubled into a catamaran — cut at the design waterline (scaled, if the
/// load is carried by scaling), and the load.
pub(crate) struct Setup {
    pub(crate) cut: Cut,
    pub(crate) layout: Vec<(usize, HullPose)>,
    pub(crate) design: Vec<(SectionalHull, Placement)>,
    /// The load [kg] and its LCG [m, fleet x]: by default the design
    /// displacement at its LCB, so at rest the fleet floats exactly at its
    /// design waterline.
    pub(crate) mass: f64,
    pub(crate) lcg: f64,
    /// The longest hull's length [m], the Froude number's.
    pub(crate) l_ref: f64,
    /// Each layout entry's role (`hull`, or a drive's `leg` or `pod`) and,
    /// for a part, its viscous form factor.
    pub(crate) roles: Vec<(String, Option<f64>)>,
    /// Each layout entry's cross-flow drag, for an inclined body (a shaft):
    /// see [`mount::CrossFlow`].
    pub(crate) cross: Vec<Option<mount::CrossFlow>>,
}

pub(crate) fn setup(
    name: &str,
    bytes: Vec<u8>,
    cut_req: &LoftRequest,
    c: &CaseParams,
) -> Result<Setup, String> {
    // The members' roles, by source index: a JSON geometry may carry a
    // drive's parts as members of their own; a file's are all hulls.
    let (source_roles, source_cross) = if native::is_native(&bytes) {
        serde_json::from_slice::<Value>(&bytes)
            .map(|g| (mount::roles(&g), mount::cross_flows(&g)))
            .unwrap_or_default()
    } else {
        Default::default()
    };
    let role_of = |i: usize| -> (String, Option<f64>) {
        source_roles
            .get(i)
            .cloned()
            .unwrap_or_else(|| ("hull".to_string(), None))
    };
    let cut = cut(name, bytes, cut_req)?;
    let rho = hullgeom::Fluid::SEAWATER_15C.density;
    // Each hull in the pose its cut settings give it (a scaled hull's scale).
    let base = cut_req.pose();
    let mut layout: Vec<(usize, HullPose)> = match c.span {
        None => cut.index.iter().map(|&i| (i, base)).collect(),
        Some(span) => {
            let hulls: Vec<usize> = (0..cut.hulls.len())
                .filter(|&k| role_of(cut.index[k]).0 == "hull")
                .collect();
            if hulls.len() != 1 {
                return Err(format!(
                    "a catamaran doubles a single hull; this file holds {}",
                    hulls.len()
                ));
            }
            let h = &cut.hulls[hulls[0]];
            let (yc, beam) = (h.placement.y, files::max_beam(&h.hull));
            if span <= beam {
                return Err(format!(
                    "span {span} m: the demihulls overlap (beam {beam:.3} m)"
                ));
            }
            // The pose's `dy` puts the copies' centreplanes at ±span/2, and
            // a drive's parts go with their hull.
            cut.index
                .iter()
                .flat_map(|&i| {
                    [0.5 * span, -0.5 * span].map(|y| (i, HullPose { dy: y - yc, ..base }))
                })
                .collect()
        }
    };
    let placed = |cut: &Cut, i: usize, pose: &HullPose| {
        let h = &cut.hulls[cut.index.iter().position(|j| *j == i).expect("a cut hull")];
        (
            h.hull.clone(),
            Placement {
                x: h.placement.x,
                y: h.placement.y + pose.dy,
            },
        )
    };
    let mut design: Vec<(SectionalHull, Placement)> = layout
        .iter()
        .map(|(i, pose)| placed(&cut, *i, pose))
        .collect();
    let vol: f64 = design.iter().map(|(h, _)| h.displaced_volume()).sum();
    // A mass carried by scaling: every hull scaled alike until the fleet
    // displaces it at the design waterline, and re-cut there — the design
    // data (and the LCB the LCG defaults to) are the scaled hull's.
    if let (Some(m), by) = (c.mass, c.mass_by) {
        if by != MassBy::Sinking {
            let ratio = m / (rho * vol);
            for (_, pose) in layout.iter_mut() {
                // On top of the hull's own scale, if it has one.
                match by {
                    MassBy::ScaleXyz => pose.scale = base.scale * ratio.cbrt(),
                    _ => pose.scale_yz = base.scale_yz * ratio.sqrt(),
                }
            }
            design = layout
                .iter()
                .map(|(i, pose)| {
                    cut.file
                        .source
                        .situate_sectional(
                            *i,
                            cut.file.waterline_z,
                            pose,
                            &Platform::default(),
                            &cut.opts,
                        )
                        .map_err(|e| e.to_string())?
                        .map(|h| (h.hull, h.placement))
                        .ok_or_else(|| "the scaled hull is dry".to_string())
                })
                .collect::<Result<Vec<_>, String>>()?;
        }
    }
    let vol: f64 = design.iter().map(|(h, _)| h.displaced_volume()).sum();
    let lcb = design
        .iter()
        .map(|(h, pl)| h.displaced_volume() * (h.lcb_x() + pl.x))
        .sum::<f64>()
        / vol.max(f64::MIN_POSITIVE);
    let roles: Vec<(String, Option<f64>)> = layout.iter().map(|(i, _)| role_of(*i)).collect();
    let cross: Vec<Option<mount::CrossFlow>> = layout
        .iter()
        .map(|(i, _)| source_cross.get(*i).copied().flatten())
        .collect();
    // The Froude number's length is the hulls', not a drive's parts'.
    let l_ref = design
        .iter()
        .zip(&roles)
        .filter(|(_, r)| r.0 == "hull")
        .map(|((h, _), _)| h.length())
        .fold(0.0f64, f64::max);
    Ok(Setup {
        mass: c.mass.unwrap_or(rho * vol),
        lcg: c.lcg.unwrap_or(lcb),
        l_ref,
        cut,
        layout,
        design,
        roles,
        cross,
    })
}

impl Setup {
    /// The platform at `(sinkage, trim)` about the LCG.
    pub(crate) fn platform(&self, sinkage: f64, trim: f64) -> Platform {
        Platform {
            sinkage,
            trim,
            pivot_x: self.lcg,
        }
    }

    /// The fleet cut at a platform state, no solve. A copy of an earlier
    /// layout entry, shifted sideways, is that cut moved.
    pub(crate) fn situate(
        &self,
        platform: &Platform,
    ) -> Result<Vec<(SectionalHull, Placement)>, String> {
        let mut out: Vec<(SectionalHull, Placement)> = Vec::new();
        for (k, (i, pose)) in self.layout.iter().enumerate() {
            if let Some(j) = self.layout[..k].iter().position(|(o, q)| {
                o == i && HullPose { dy: 0.0, ..*q } == HullPose { dy: 0.0, ..*pose }
            }) {
                let (h, p) = out[j].clone();
                out.push((
                    h,
                    Placement {
                        y: p.y + pose.dy - self.layout[j].1.dy,
                        ..p
                    },
                ));
                continue;
            }
            let h = self
                .cut
                .file
                .source
                .situate_sectional(
                    *i,
                    self.cut.file.waterline_z,
                    pose,
                    platform,
                    &self.cut.opts,
                )
                .map_err(|e| e.to_string())?
                .ok_or("the hull is dry at this attitude")?;
            out.push((h.hull, h.placement));
        }
        Ok(out)
    }

    /// The whole of each hull at a platform state, topsides included, for
    /// display: the source's tessellation, x forward, y across, z up from
    /// the water.
    pub(crate) fn meshes(&self, platform: &Platform) -> Result<Vec<Value>, String> {
        let mut meshes = Vec::new();
        for (i, pose) in &self.layout {
            let (v, t) = self
                .cut
                .file
                .source
                .posed_tessellation(*i, self.cut.file.waterline_z, pose, platform)
                .map_err(|e| e.to_string())?;
            let flat: Vec<f32> = v.iter().flat_map(|p| p.map(|c| c as f32)).collect();
            let idx: Vec<u32> = t.iter().flatten().copied().collect();
            meshes.push(json!({
                "vertices": b64(bytemuck_f32(&flat)),
                "triangles": b64(bytemuck_u32(&idx)),
            }));
        }
        Ok(meshes)
    }
}

/// [`flow`], reporting each stage (and each force evaluation of the
/// equilibrium solve) through `report`; returning `false` from it stops the
/// computation with the error [`CANCELLED`].
pub fn flow_with_progress(
    name: &str,
    bytes: Vec<u8>,
    req: &FlowRequest,
    report: &mut dyn FnMut(&Progress) -> bool,
) -> Result<Value, String> {
    use hullgeom::float::{solve_equilibrium_sectional_dynamic, LoadCase};
    use hullgeom::source::SourceHull;
    use thinship::nearfield::{free_surface, hull_pressure, NearFieldOptions};
    let study = &req.study;
    let t0 = std::time::Instant::now();
    let mut track = Tracker {
        stages: flow_stages(study.dynamic),
        report,
    };
    track.at(0, 0.0, name.to_string())?;
    let s = setup(name, bytes, &req.cut, &req.case)?;
    track.at(
        0,
        1.0,
        format!(
            "{} hull{}",
            s.layout.len(),
            if s.layout.len() == 1 { "" } else { "s" }
        ),
    )?;
    let (mass, lcg, l_ref) = (s.mass, s.lcg, s.l_ref);
    let closure = study.closure.transom();
    let cond = Conditions::seawater(study.froude * (STANDARD_GRAVITY * l_ref).sqrt());
    let rho = cond.fluid.density;
    let wave = WaveOptions {
        transom: closure,
        ..WaveOptions::default()
    };
    let squat = thinship::squat::SquatOptions {
        wave,
        ..Default::default()
    };

    // The attitude: the dynamic equilibrium at this speed, or held (at
    // rest's, when not floated at speed).
    let mut platform = Platform::default();
    let mut solved = None;
    let at_attitude: Vec<(SectionalHull, Placement)> = if let Some((sinkage, trim)) = req.hold {
        platform = s.platform(sinkage, trim);
        track.at(1, 0.5, "at the held attitude".into())?;
        s.situate(&platform)?
    } else if study.dynamic {
        let sources: Vec<SourceHull> = s
            .layout
            .iter()
            .map(|&(index, pose)| SourceHull {
                source: s.cut.file.source.as_ref(),
                index,
                waterline_z: s.cut.file.waterline_z,
                pose,
            })
            .collect();
        track.at(1, 0.0, "starting the Newton solve".into())?;
        // Count the solve's force evaluations (each is most of an
        // iteration's cost) and report the latest lift; the count is open
        // ended, so the stage's share fills asymptotically.
        let inner = thinship::sectional::dynamic_load_closure(&cond, lcg, &squat);
        let mut evaluations = 0usize;
        let weight = mass * cond.gravity;
        let track_ref = &mut track;
        let report = move |d: &hullgeom::float::DynamicLoad, vol: f64| {
            evaluations += 1;
            track_ref
                .at(
                    1,
                    1.0 - (-(evaluations as f64) / 10.0).exp(),
                    format!(
                        "force evaluation {evaluations} · lift {:.1}% of weight · \
                         displacing {:.1}% of the load",
                        100.0 * d.force_up / weight,
                        100.0 * rho * vol / mass,
                    ),
                )
                .map_err(thinship::Error::InvalidInput)
        };
        let extra = req
            .thrust
            .as_ref()
            .map(|t| t.load(lcg, resistance_line(&s.design, &s.roles)))
            .unwrap_or(hullgeom::float::DynamicLoad {
                force_up: 0.0,
                moment_bow_up: 0.0,
            });
        let counted = Counted {
            inner,
            report,
            extra,
        };
        let eq = solve_equilibrium_sectional_dynamic(
            &sources,
            &LoadCase {
                mass,
                lcg: Some(lcg),
            },
            rho,
            cond.gravity,
            &s.cut.opts,
            counted,
            req.warm,
        )
        .map_err(|e| {
            if e.to_string().contains(CANCELLED) {
                return CANCELLED.to_string();
            }
            let hint = if e.to_string().contains("waterplane") {
                " — if the geometry ends at the waterline (no topsides), it cannot sink or \
                 trim; hold it at rest's attitude instead"
            } else {
                ""
            };
            format!("dynamic equilibrium at Fn {}: {e}{hint}", study.froude)
        })?;
        platform = s.platform(eq.sinkage, eq.trim);
        solved = Some((
            eq.sinkage,
            eq.trim,
            eq.iterations,
            eq.dynamic,
            eq.lift_fraction,
        ));
        eq.fleet.members
    } else {
        s.design.clone()
    };
    let members: Vec<(&SectionalHull, Placement)> =
        at_attitude.iter().map(|(h, p)| (h, *p)).collect();
    let t_attitude = t0.elapsed().as_secs_f64();

    let nf = NearFieldOptions {
        closure,
        ..NearFieldOptions::default()
    };
    track.at(2, 0.0, format!("{} rows deep per hull", nf.depth_rows))?;
    let pressures = if req.field {
        hull_pressure(&members, &cond, &nf).map_err(|e| e.to_string())?
    } else {
        Vec::new()
    };

    // The free surface: from ahead of the bows to ~1.5 lengths astern (or
    // `astern`).
    let (mut xa, mut xb, mut yh) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
    for (h, pl) in &members {
        let (a, b) = h.x_range();
        xa = xa.min(a + pl.x);
        xb = xb.max(b + pl.x);
        yh = yh.max(pl.y.abs() + 0.5 * files::max_beam(h));
    }
    let (x0, x1) = (xa - study.astern.unwrap_or(1.5) * l_ref, xb + 0.4 * l_ref);
    let yh = (yh + 0.45 * l_ref).max(0.3 * (x1 - x0));
    let nx = study.grid;
    // An odd row count, so y = 0 is a row and a symmetric fleet mirrors.
    let nh = ((nx as f64) * yh / (x1 - x0)).round().clamp(8.0, 600.0) as usize + 1;
    let ny = 2 * nh - 1;
    track.at(
        3,
        0.0,
        format!("{nx} × {ny} grid over {:.0} × {:.0} m", x1 - x0, 2.0 * yh),
    )?;
    // A fleet symmetric about y = 0 — one hull on the centreline, or
    // mirrored pairs — has a symmetric field: compute the half with y ≥ 0.
    let symmetric = members.iter().all(|(h, pl)| {
        let tol = 1e-6 * l_ref;
        pl.y.abs() < tol
            || members.iter().any(|(o, q)| {
                (q.y + pl.y).abs() < tol
                    && (q.x - pl.x).abs() < tol
                    && (o.displaced_volume() - h.displaced_volume()).abs()
                        <= 1e-6 * h.displaced_volume()
                    && (o.length() - h.length()).abs() < tol
            })
    });
    let g = if !req.field {
        None
    } else {
        Some(if let Some(span) = req.case.span {
            // A catamaran's field is its demihull's, twice, shifted by ±span/2
            // (exact in thin-ship theory): one hull's field on a y-grid whose
            // spacing divides span/2, summed at offset rows. Both demihulls are
            // the same cut at the same attitude, and each is symmetric about its
            // own centreplane, so the sum is symmetric too — half is computed.
            let dxg = (x1 - x0) / (nx - 1) as f64;
            let m = ((0.5 * span) / dxg).round().max(1.0) as usize;
            let dy = 0.5 * span / m as f64;
            let nh = (yh / dy).ceil() as usize + 1;
            let (yh, ny) = ((nh - 1) as f64 * dy, 2 * nh - 1);
            let rows = nh + m;
            let (h0, p0) = members[0];
            let one = [(h0, Placement { x: p0.x, y: 0.0 })];
            let half = free_surface(
                &one,
                &cond,
                &nf,
                x0,
                x1,
                0.0,
                (rows - 1) as f64 * dy,
                nx,
                rows,
            )
            .map_err(|e| e.to_string())?;
            let row = |j: usize| &half.zeta[j * nx..(j + 1) * nx];
            let mut zeta = Vec::with_capacity(nx * ny);
            for iy in 0..ny {
                let k = iy.abs_diff(nh - 1);
                let (a, b) = (row(k.abs_diff(m)), row(k + m));
                zeta.extend(a.iter().zip(b).map(|(a, b)| a + b));
            }
            thinship::WaveGrid {
                y0: -yh,
                y1: yh,
                ny,
                zeta,
                ..half
            }
        } else if symmetric {
            let half = free_surface(&members, &cond, &nf, x0, x1, 0.0, yh, nx, nh)
                .map_err(|e| e.to_string())?;
            let mut zeta = Vec::with_capacity(nx * ny);
            for iy in 0..ny {
                let k = iy.abs_diff(nh - 1);
                zeta.extend_from_slice(&half.zeta[k * nx..(k + 1) * nx]);
            }
            thinship::WaveGrid {
                y0: -yh,
                ny,
                zeta,
                ..half
            }
        } else {
            free_surface(&members, &cond, &nf, x0, x1, -yh, yh, nx, ny)
                .map_err(|e| e.to_string())?
        })
    };
    track.at(
        4,
        0.0,
        if solved.is_some() {
            "resistance, then the hull meshes".into()
        } else {
            "resistance, dynamic force, then the hull meshes".into()
        },
    )?;
    let t_field = t0.elapsed().as_secs_f64() - t_attitude;

    let viscous = thinship::ViscousOptions {
        form_factor: study.form_factor.unwrap_or(0.0),
        roughness: study
            .roughness
            .map_or(thinship::Roughness::None, thinship::Roughness::DeltaCf),
    };
    let res = thinship::sectional::multihull_resistance(&members, &cond, &wave, &viscous)
        .map_err(|e| e.to_string())?;

    // Forces: the solved balance, or at a held attitude the dynamic load
    // there, unbalanced (how far the attitude is from its own equilibrium at
    // speed), or at the design attitude the force and the first-order
    // sinkage and trim it implies.
    let forces = match (solved, req.hold) {
        (Some((sinkage, trim, iterations, d, lift)), _) => json!({
            "fz": d.force_up,
            "moment": d.moment_bow_up,
            "lift_fraction": lift,
            "sinkage": sinkage,
            "trim_deg": trim.to_degrees(),
            "trim_rad": trim,
            "iterations": iterations,
            "solved": true,
        }),
        (None, Some((sinkage, trim))) => {
            let d = thinship::sectional::multihull_dynamic_force(&members, &cond, lcg, &squat)
                .map_err(|e| e.to_string())?;
            json!({
                "fz": d.force_up,
                "moment": d.moment_bow_up,
                // Of the weight, as the solved studies report it.
                "lift_fraction": d.force_up / (mass * cond.gravity),
                "sinkage": sinkage,
                "trim_deg": trim.to_degrees(),
                "trim_rad": trim,
                "solved": false,
                "held": true,
            })
        }
        (None, None) => {
            let (mut aw, mut mw, mut iw) = (0.0, 0.0, 0.0);
            for (h, pl) in &members {
                let a = h.waterplane_area();
                aw += a;
                mw += h.waterplane_moment() + pl.x * a;
                iw += h.waterplane_second_moment()
                    + 2.0 * pl.x * h.waterplane_moment()
                    + pl.x * pl.x * a;
            }
            let lcf = mw / aw.max(f64::MIN_POSITIVE);
            let i_l = iw - mw * mw / aw.max(f64::MIN_POSITIVE);
            let d = thinship::sectional::multihull_dynamic_force(&members, &cond, lcf, &squat)
                .map_err(|e| e.to_string())?;
            json!({
                "fz": d.force_up,
                "moment": d.moment_bow_up,
                "lift_fraction": d.lift_fraction,
                "sinkage": -d.force_up / (rho * cond.gravity * aw),
                "trim_deg": (d.moment_bow_up / (rho * cond.gravity * i_l)).to_degrees(),
                "solved": false,
            })
        }
    };

    let meshes = if req.field {
        s.meshes(&platform)?
    } else {
        Vec::new()
    };
    let hulls: Vec<Value> = pressures
        .iter()
        .zip(meshes)
        .map(|(p, mesh)| {
            json!({
                "x": round(&p.x),
                "depth": round(&p.depth),
                "half_beam": round(&p.half_beam),
                "cp": round(&p.cp),
                "y": p.y,
                "force_up": p.force_up,
                "mesh": mesh,
            })
        })
        .collect();
    // A drive's parts are friction on their own lengths, as the hull's
    // members are, but at their own form factors (Hoerner's, for a foil and
    // for a body of revolution) rather than the hull's; an inclined shaft
    // adds the drag of the flow across it, at its angle to the flow (its
    // own plus the trim). The roughness allowance is the skin's, so every
    // part carries it.
    let (mut rv, mut r_app) = (0.0, 0.0);
    let parts_known = res.viscous.len() == s.roles.len();
    let q = 0.5 * rho * cond.speed * cond.speed;
    for (k, v) in res.viscous.iter().enumerate() {
        match s.roles.get(k).filter(|_| parts_known) {
            Some((role, Some(kp))) if role != "hull" => {
                let cross = s.cross.get(k).copied().flatten().map_or(0.0, |c| {
                    q * c.area * c.cd * (c.angle + platform.trim).sin().abs().powi(3)
                });
                let r = ((1.0 + kp) * v.cf + v.roughness_cf) * q * v.wetted_surface + cross;
                rv += r;
                r_app += r;
            }
            _ => rv += v.resistance,
        }
    }
    // A transom running dry: the water no longer closes in behind it, so the
    // hull loses the hydrostatic thrust ρg∫∫z dA the transom's face would
    // carry — faded in across the closure's own wet-to-dry band of Fn_T.
    // Thin-ship theory has no term for it, closure or not.
    let r_transom: f64 = if study.transom_drag {
        members
            .iter()
            .filter_map(|(h, _)| h.transom())
            .filter(|t| t.depth > 0.0)
            .map(|t| {
                let fn_t = cond.speed / (cond.gravity * t.depth).sqrt();
                rho * cond.gravity * t.depth_moment * thinship::transom_dryness(fn_t)
            })
            .sum()
    } else {
        0.0
    };
    let rt = res.wave.resistance + rv + r_transom;
    let q_s = 0.5
        * rho
        * cond.speed
        * cond.speed
        * members.iter().map(|(h, _)| h.wetted_surface()).sum::<f64>();
    let mut forces = forces;
    if s.roles.iter().any(|r| r.0 != "hull") {
        forces["r_appendages_viscous"] = json!(r_app);
    }
    if study.transom_drag {
        forces["r_transom"] = json!(r_transom);
    }
    for (k, v) in [
        ("rw", res.wave.resistance),
        ("rv", rv),
        ("rt", rt),
        ("pe", rt * cond.speed),
        ("cw", res.cw),
        ("ct", rt / q_s.max(f64::MIN_POSITIVE)),
        ("interference", res.interference),
        ("mass", mass),
        ("lcg", lcg),
        (
            "displaced_volume",
            members.iter().map(|(h, _)| h.displaced_volume()).sum(),
        ),
        (
            "length",
            members.iter().map(|(h, _)| h.length()).fold(0.0, f64::max),
        ),
    ] {
        forces[k] = json!(v);
    }
    Ok(json!({
        "froude": study.froude,
        "speed": cond.speed,
        "transverse_wavelength": 2.0 * std::f64::consts::PI * cond.speed * cond.speed / cond.gravity,
        "seconds": t0.elapsed().as_secs_f64(),
        "timing": { "attitude": t_attitude, "field": t_field },
        "hulls": hulls,
        "surface": g.as_ref().map(|g| json!({
            "x0": g.x0, "x1": g.x1, "y0": g.y0, "y1": g.y1, "nx": g.nx, "ny": g.ny,
            "zeta": b64(bytemuck_f32(&g.zeta.iter().map(|&z| z as f32).collect::<Vec<_>>())),
        })),
        "forces": forces,
    }))
}

pub(crate) fn bytemuck_f32(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

pub(crate) fn bytemuck_u32(v: &[u32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Standard base64, for little-endian arrays the page decodes into typed
/// arrays (a third the size of JSON numbers, and no parsing).
pub(crate) fn b64(bytes: Vec<u8>) -> String {
    const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for c in bytes.chunks(3) {
        let n = (c[0] as u32) << 16
            | (*c.get(1).unwrap_or(&0) as u32) << 8
            | *c.get(2).unwrap_or(&0) as u32;
        for k in 0..4 {
            if k <= c.len() {
                out.push(A[(n >> (18 - 6 * k) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// Values rounded to 6 significant figures: a third of the JSON.
fn round(v: &[f64]) -> Vec<Value> {
    v.iter()
        .map(|&x| {
            if x == 0.0 || !x.is_finite() {
                json!(0.0)
            } else {
                let e = 5 - x.abs().log10().floor() as i32;
                let k = 10f64.powi(e);
                json!((x * k).round() / k)
            }
        })
        .collect()
}

fn sectioned_json(imp: &SectionalImport, kind: Kind) -> Value {
    let r = &imp.report;
    let sides = if r.two_sided {
        format!("both sides averaged about y = {:.4} m", r.centerplane)
    } else {
        format!("one side about y = {:.4} m", r.centerplane)
    };
    let what = match kind {
        Kind::Iges => format!("IGES, {} patches", r.patches),
        Kind::Stl => format!("STL, {} triangles", r.patches),
    };
    let mut lines = vec![format!(
        "source: {what}, cut at {} stations over x {:.4}..{:.4} m ({sides})",
        r.stations, r.x_range.0, r.x_range.1
    )];
    if r.dropped_stations > 0 {
        lines.push(format!(
            "{} interior stations could not be sectioned (bridged by the interpolant)",
            r.dropped_stations
        ));
    }
    if r.max_asymmetry > 1e-3 * r.draft.max(1e-9) {
        lines.push(format!(
            "port and starboard differ by up to {:.3e} m (averaged)",
            r.max_asymmetry
        ));
    }
    hull_json(&imp.hull, imp.placement, Some(&imp.sections), lines)
}

/// One hull, drawn as the physics uses it (see the module docs).
fn hull_json(
    hull: &SectionalHull,
    placement: Placement,
    rays: Option<&Vec<(f64, Vec<(f64, f64)>)>>,
    lines: Vec<String>,
) -> Value {
    // Each station as the curve the quadrature integrates (not a polyline
    // through its quadrature nodes, which are graded toward the waterline).
    let stations: Vec<(f64, Vec<(f64, f64)>)> =
        hull.curves().map(|(x, o)| (x, o.to_vec())).collect();
    // A see-through surface between stations, for orientation only (the
    // kernel interpolates each station's depth integral along x, not a
    // surface): rows joined at equal fractions of girth, so neighbouring
    // sections with different node spacing still meet cleanly.
    let rows = 48usize;
    let (mut mx, mut mz, mut my) = (Vec::new(), Vec::new(), Vec::new());
    for (x, o) in &stations {
        for (y, z) in by_girth(o, rows) {
            mx.push(*x);
            mz.push(z);
            my.push(y);
        }
    }
    let keel: Vec<(f64, f64)> = stations
        .iter()
        .map(|(x, o)| (*x, o.last().map_or(0.0, |p| p.1)))
        .collect();
    let beam = 2.0
        * stations
            .iter()
            .flat_map(|(_, o)| o.first().map(|p| p.0))
            .fold(0.0, f64::max);
    // The depth-integral curves at κ = 0 (sectional area) and at a short
    // wave's decay rate: λ = 2 at Fn 0.15, κ = νλ², ν = g/U² = 1/(0.0225 L).
    let kappas = [0.0, 4.0 / (0.0225 * hull.length())];
    let curves: Vec<Value> = kappas
        .iter()
        .map(|&kappa| {
            let (st, c) = hull.depth_integral_curve(kappa, 8);
            json!({ "kappa": kappa, "stations": st, "curve": c })
        })
        .collect();
    // The transom and its section — the aft end station's outline, whose
    // depth integral is the closure's depth factor — so the page can draw
    // the virtual appendage at whatever speed and closure it is asked about.
    let transom = hull.transom().map(|t| {
        let aft = stations
            .iter()
            .min_by(|a, b| (a.0 - t.x).abs().total_cmp(&(b.0 - t.x).abs()));
        json!({
            "x": t.x,
            "depth": t.depth,
            "half_beam": t.half_beam,
            "area": t.area,
            "area_ratio": t.area / hull.max_section_area().max(1e-300),
            "outline": aft.map(|s| by_girth(&s.1, rows)).unwrap_or_default(),
            // Z_T at each curve's κ: the kernel extends Z over the hollow as
            // Z_T · φ(s).
            "z_t": kappas
                .iter()
                .map(|&k| hull.depth_integral_curve(k, 1).0.first().map_or(0.0, |p| p.1))
                .collect::<Vec<_>>(),
        })
    });
    json!({
        "placement": { "x": placement.x, "y": placement.y },
        "diagnostics": lines,
        "length": hull.length(),
        "beam": beam,
        "draft": hull.draft(),
        "displaced_volume": hull.displaced_volume(),
        "wetted_surface": hull.wetted_surface(),
        "lcb_x": hull.lcb_x(),
        "waterplane_area": hull.waterplane_area(),
        "stations": stations,
        "rays": rays,
        "mesh": { "nx": stations.len(), "nz": rows + 1, "x": mx, "z": mz, "y": my },
        "keel": keel,
        "area": curves,
        "transom": transom,
    })
}

/// A section outline resampled at `rows + 1` points evenly spaced in girth
/// (a point for an empty section).
fn by_girth(o: &[(f64, f64)], rows: usize) -> Vec<(f64, f64)> {
    if o.len() < 2 {
        return vec![(0.0, 0.0); rows + 1];
    }
    let mut cum = vec![0.0];
    for w in o.windows(2) {
        let d = ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt();
        cum.push(cum.last().unwrap() + d);
    }
    let total = *cum.last().unwrap();
    (0..=rows)
        .map(|k| {
            let g = total * k as f64 / rows as f64;
            let j = cum.partition_point(|&c| c < g).clamp(1, o.len() - 1);
            let t = if cum[j] > cum[j - 1] {
                (g - cum[j - 1]) / (cum[j] - cum[j - 1])
            } else {
                0.0
            };
            (
                o[j - 1].0 + t * (o[j].0 - o[j - 1].0),
                o[j - 1].1 + t * (o[j].1 - o[j - 1].1),
            )
        })
        .collect()
}

/// What the hull list shows of a cut: per hull, the principal dimensions and
/// hydrostatics.
pub fn hull_summary(sections: &Value) -> Value {
    let hulls: Vec<Value> = sections["hulls"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|h| {
                    let mut s = json!({});
                    for k in [
                        "length",
                        "beam",
                        "draft",
                        "displaced_volume",
                        "wetted_surface",
                        "lcb_x",
                        "waterplane_area",
                    ] {
                        s[k] = h[k].clone();
                    }
                    s["transom"] = json!(!h["transom"].is_null());
                    s
                })
                .collect()
        })
        .unwrap_or_default();
    json!({ "hulls": hulls, "notes": sections["notes"] })
}

/// A response measure: its name, and how to read it off a wavelength's point.
type Measure<'a> = (&'a str, &'a dyn Fn(&Value) -> Option<f64>);

/// What a study in waves is plotted by: each response's peak over the
/// wavelengths and where it is, the added resistance's, the roll stability
/// at the attitude, and the irregular sea's statistics.
pub fn wave_scalars(v: &Value) -> Value {
    let sk = &v["seakeeping"];
    let h = &sk["headings"][0];
    let pts = h["points"].as_array().cloned().unwrap_or_default();
    let abs = |z: &Value| {
        z.as_array().map(|a| {
            a[0].as_f64()
                .unwrap_or(0.0)
                .hypot(a[1].as_f64().unwrap_or(0.0))
        })
    };
    let peak = |f: &dyn Fn(&Value) -> Option<f64>| {
        pts.iter()
            .filter_map(|p| Some((f(p)?, p["lambda"].as_f64()?)))
            .fold(None, |best: Option<(f64, f64)>, (y, l)| match best {
                Some((b, _)) if b >= y => best,
                _ => Some((y, l)),
            })
    };
    let mut out = json!({
        "froude": v["froude"],
        "speed": v["speed"],
        "heading": h["heading"],
        "sinkage": v["attitude"]["sinkage"],
        "trim_rad": v["attitude"]["trim_rad"],
        "trim_deg": v["attitude"]["trim_deg"],
        "gm_t": sk["gm_t"],
        "roll_period": sk["roll_period"],
        "failed_points": pts.iter().filter(|p| !p["error"].is_null()).count(),
    });
    let heave = |p: &Value| abs(&p["heave"]);
    let pitch = |p: &Value| abs(&p["pitch"]);
    let roll = |p: &Value| abs(&p["roll"]);
    let added = |p: &Value| p["raw_gb"].as_f64();
    let measures: [Measure; 4] = [
        ("heave", &heave),
        ("pitch", &pitch),
        ("roll", &roll),
        ("added_resistance", &added),
    ];
    for (k, f) in measures {
        if let Some((y, l)) = peak(f) {
            out[format!("{k}_peak")] = json!(y);
            out[format!("{k}_peak_lambda")] = json!(l);
        }
    }
    let sea = &h["sea"];
    if sea.is_object() && sea["error"].is_null() {
        for k in [
            "heave",
            "pitch_deg",
            "accel_bow",
            "accel_lcg",
            "raw_gb",
            "raw_far",
        ] {
            out[format!("sea_{k}")] = sea[k].clone();
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exact Wigley as IGES, shown by sections.
    #[test]
    fn a_wigley_is_shown_in_sections() {
        let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = hullgeom::iges::write(&surfaces, "wigley").unwrap();
        let v = loft("w.igs", text.into_bytes(), &LoftRequest::default()).unwrap();
        let h = &v["hulls"][0];
        let vol = h["displaced_volume"].as_f64().unwrap();
        // The Wigley's exact volume, 4/9 L B T.
        let exact = 4.0 / 9.0 * 10.0 * 0.625;
        assert!((vol - exact).abs() < 1e-6 * exact, "{vol}");
        assert!(h["transom"].is_null());
        assert!(!h["stations"].as_array().unwrap().is_empty());
    }

    fn wigley() -> Vec<u8> {
        let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        hullgeom::iges::write(&surfaces, "wigley")
            .unwrap()
            .into_bytes()
    }

    /// A calm-water request from a study's and a case's JSON.
    fn request(study: &str, case: &str) -> FlowRequest {
        FlowRequest {
            cut: LoftRequest::default(),
            case: serde_json::from_str::<CaseParams>(case)
                .unwrap()
                .canonical()
                .unwrap(),
            study: serde_json::from_str::<StudyParams>(study)
                .unwrap()
                .canonical()
                .unwrap(),
            warm: None,
            hold: None,
            thrust: None,
            field: true,
        }
    }

    /// The flow at one speed: pressure per hull, the surface grid asked for,
    /// and a dynamic lift that sinks the hull. This Wigley ends at its
    /// waterline (no topsides to sink into), so at the design attitude.
    #[test]
    fn a_wigley_flow_has_pressure_waves_and_lift() {
        let req = request(r#"{"froude": 0.35, "grid": 60, "dynamic": false}"#, "{}");
        let v = flow("w.igs", wigley(), &req).unwrap();
        let h = &v["hulls"][0];
        let (nx, nz) = (
            h["x"].as_array().unwrap().len(),
            h["depth"].as_array().unwrap().len(),
        );
        assert_eq!(h["cp"].as_array().unwrap().len(), nx * nz);
        assert_eq!(v["surface"]["nx"], 60);
        let f = &v["forces"];
        assert!(
            f["fz"].as_f64().unwrap() < 0.0 && f["sinkage"].as_f64().unwrap() > 0.0,
            "{f}"
        );
        assert!(f["rw"].as_f64().unwrap() > 0.0);
    }

    /// Progress walks every stage in order, its fraction never goes back,
    /// and a refusal from the callback stops the flow.
    #[test]
    fn flow_progress_advances_and_can_be_cancelled() {
        let (bytes, req) = (
            wigley(),
            request(r#"{"froude": 0.35, "grid": 40, "dynamic": false}"#, "{}"),
        );
        let mut seen: Vec<(usize, f64)> = Vec::new();
        flow_with_progress("w.igs", bytes.clone(), &req, &mut |p| {
            seen.push((p.step, p.fraction));
            true
        })
        .unwrap();
        assert!(seen.windows(2).all(|w| w[1].1 >= w[0].1), "{seen:?}");
        assert_eq!(seen.first().map(|s| s.0), Some(1));
        assert_eq!(seen.last().map(|s| s.0), Some(5));
        let mut calls = 0;
        let e = flow_with_progress("w.igs", bytes, &req, &mut |_| {
            calls += 1;
            calls < 3
        })
        .unwrap_err();
        assert_eq!(e, CANCELLED);
        assert_eq!(calls, 3);
    }

    /// A catamaran's free surface, summed from one demihull's field shifted
    /// by ±span/2, is the two-hull field computed directly.
    #[test]
    fn a_catamaran_field_is_its_demihulls_summed() {
        let text = String::from_utf8(wigley()).unwrap();
        let req = request(
            r#"{"froude": 0.35, "grid": 60, "dynamic": false}"#,
            r#"{"span": 3.1}"#,
        );
        let v = flow("w.igs", text.clone().into_bytes(), &req).unwrap();
        let sfc = &v["surface"];
        let get = |k: &str| sfc[k].as_f64().unwrap();
        let (nx, ny) = (get("nx") as usize, get("ny") as usize);
        let bytes = unb64(sfc["zeta"].as_str().unwrap());
        let summed: Vec<f64> = bytes
            .chunks(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64)
            .collect();
        // The same two hulls, directly, on the same grid.
        let hull = hullgeom::iges::source_fleet(&text, 0.0)
            .unwrap()
            .situate_sectional(
                0,
                0.0,
                &HullPose::default(),
                &Platform::default(),
                &Default::default(),
            )
            .unwrap()
            .unwrap()
            .hull;
        let fleet = [
            (&hull, Placement { x: 0.0, y: 1.55 }),
            (&hull, Placement { x: 0.0, y: -1.55 }),
        ];
        let cond = Conditions::seawater(0.35 * (STANDARD_GRAVITY * hull.length()).sqrt());
        let direct = thinship::nearfield::free_surface(
            &fleet,
            &cond,
            &Default::default(),
            get("x0"),
            get("x1"),
            get("y0"),
            get("y1"),
            nx,
            ny,
        )
        .unwrap();
        let peak = direct.zeta.iter().fold(0.0f64, |m, z| m.max(z.abs()));
        let worst = summed
            .iter()
            .zip(&direct.zeta)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(
            worst < 1e-4 * peak,
            "summed vs direct: {worst} of peak {peak}"
        );
        assert_eq!(v["hulls"].as_array().unwrap().len(), 2);
    }

    fn unb64(s: &str) -> Vec<u8> {
        const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let val = |c: u8| A.iter().position(|&a| a == c).unwrap() as u32;
        let mut out = Vec::new();
        for q in s.as_bytes().chunks(4) {
            let n = q.iter().take_while(|&&c| c != b'=').count();
            let v = q[..n].iter().fold(0u32, |acc, &c| acc << 6 | val(c)) << (6 * (4 - n));
            out.extend(&[(v >> 16) as u8, (v >> 8) as u8, v as u8][..n - 1]);
        }
        out
    }

    /// A mass carried by scaling keeps the design waterline: the scaled
    /// Wigley displaces exactly that mass there, its length scaled by the
    /// cube root of the mass ratio (x, y, z) or kept (beam and draft only).
    #[test]
    fn a_mass_carried_by_scaling_keeps_the_waterline() {
        let rho = Conditions::seawater(1.0).fluid.density;
        let v0 = 4.0 / 9.0 * 10.0 * 0.625;
        let mass = 1.2 * rho * v0;
        for (by, length) in [("scale", 10.0 * 1.2f64.cbrt()), ("scale_yz", 10.0)] {
            let req = request(
                r#"{"froude": 0.3, "grid": 40, "dynamic": false}"#,
                &format!(r#"{{"mass": {mass}, "mass_by": "{by}"}}"#),
            );
            let v = flow("w.igs", wigley(), &req).unwrap();
            let f = &v["forces"];
            let (vol, len) = (
                f["displaced_volume"].as_f64().unwrap(),
                f["length"].as_f64().unwrap(),
            );
            assert!(
                (vol - 1.2 * v0).abs() < 1e-5 * v0,
                "{by}: volume {vol} vs {}",
                1.2 * v0
            );
            assert!(
                (len - length).abs() < 1e-6 * length,
                "{by}: length {len} vs {length}"
            );
        }
    }

    /// A case's statics: at the design load it floats at its
    /// design waterline, and the GZ curve's slope at upright (whole
    /// sections, clipped) is the GM_T the waterplane gives. A Wigley form a
    /// metre deep, cut 0.375 m below its top (so it has freeboard to heel
    /// into); deep and narrow, it rolls over with G at the waterline, so G
    /// is put below it.
    #[test]
    fn a_case_has_its_statics() {
        let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 1.0).unwrap();
        let deep = hullgeom::iges::write(&surfaces, "wigley")
            .unwrap()
            .into_bytes();
        let cut = LoftRequest {
            waterline: Some(-0.375),
            ..LoftRequest::default()
        };
        let c = serde_json::from_str::<CaseParams>(r#"{"vcg": -0.3}"#)
            .unwrap()
            .canonical()
            .unwrap();
        let v = statics("w.igs", deep, &cut, &c).unwrap();
        let rest = &v["at_rest"];
        assert!(rest["sinkage"].as_f64().unwrap().abs() < 1e-4, "{rest}");
        assert!(rest["trim_rad"].as_f64().unwrap().abs() < 1e-5, "{rest}");
        let gm_t = v["roll"]["gm_t"].as_f64().unwrap();
        let gz = &v["gz"];
        let gm = gz["gm"].as_f64().unwrap_or_else(|| panic!("{gz}"));
        assert!(gm_t > 0.0, "{}", v["roll"]);
        assert!(
            (gm - gm_t).abs() < 0.02 * gm_t,
            "GZ slope {gm} vs GM_T {gm_t}"
        );
        assert_eq!(gz["heel_deg"].as_array().unwrap().len(), 91);
    }

    /// A study in waves carries its responses: in long head waves a Wigley
    /// heaves and pitches with the water, and the far field is reported
    /// only where it applies.
    #[test]
    fn a_study_in_waves_carries_its_responses() {
        let case = serde_json::from_str::<CaseParams>(r#"{"vcg": -0.3}"#)
            .unwrap()
            .canonical()
            .unwrap();
        let study = |heading: f64| {
            serde_json::from_str::<StudyParams>(&format!(
                r#"{{"froude": 0.2, "dynamic": false, "waves": {{"heading": {heading}, "lambdas": [6]}}}}"#
            ))
            .unwrap()
            .canonical()
            .unwrap()
        };
        let go = |r: &StudyParams| {
            waves_with_progress(
                "w.igs",
                wigley(),
                &LoftRequest::default(),
                &case,
                r,
                (0.0, 0.0),
                &mut |_| true,
            )
            .unwrap()
        };
        let head = go(&study(180.0));
        let sk = &head["seakeeping"];
        assert!(sk["gm_t"].as_f64().unwrap() > 0.0, "{}", sk);
        let p = &sk["headings"][0]["points"][0];
        let abs = |z: &Value| z[0].as_f64().unwrap().hypot(z[1].as_f64().unwrap());
        let k = p["k"].as_f64().unwrap();
        assert!(
            (abs(&p["heave"]) - 1.0).abs() < 0.1,
            "heave {}",
            abs(&p["heave"])
        );
        assert!(
            (abs(&p["pitch"]) / k - 1.0).abs() < 0.15,
            "pitch {}",
            abs(&p["pitch"]) / k
        );
        assert!(p["raw_far"].is_f64());
        let bow = go(&study(120.0));
        assert!(bow["seakeeping"]["headings"][0]["points"][0]["raw_far"].is_null());
    }

    /// A hull scaled in its settings is cut scaled about its design
    /// waterline: the whole Wigley's length by k and volume by k³, its beam
    /// and draft only its volume by k² and its length kept.
    #[test]
    fn a_scaled_hull_is_cut_scaled() {
        let exact = 4.0 / 9.0 * 10.0 * 0.625;
        for (req, length, volume) in [
            (
                LoftRequest {
                    scale: Some(1.2),
                    ..LoftRequest::default()
                },
                12.0,
                exact * 1.2f64.powi(3),
            ),
            (
                LoftRequest {
                    scale_yz: Some(0.9),
                    ..LoftRequest::default()
                },
                10.0,
                exact * 0.81,
            ),
        ] {
            let v = loft("w.igs", wigley(), &req).unwrap();
            let h = &v["hulls"][0];
            let (l, vol) = (
                h["length"].as_f64().unwrap(),
                h["displaced_volume"].as_f64().unwrap(),
            );
            assert!((l - length).abs() < 1e-6 * length, "{req:?}: length {l}");
            assert!(
                (vol - volume).abs() < 1e-6 * volume,
                "{req:?}: volume {vol}"
            );
        }
        // Unset, the scale is not in the stored settings: an unscaled hull's
        // settings are what they were before hulls could be scaled.
        let plain = serde_json::to_string(&LoftRequest::default()).unwrap();
        assert!(!plain.contains("scale"), "{plain}");
    }

    #[test]
    fn an_stl_needs_its_units() {
        let e = loft("x.stl", vec![0; 200], &LoftRequest::default()).unwrap_err();
        assert!(e.contains("units"), "{e}");
    }
}
