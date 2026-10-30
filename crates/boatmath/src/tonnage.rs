//! Gross and net tonnage for vessels less than 24 m in statutory length,
//! under Canada's *Vessel Registration and Tonnage Regulations* (SOR/2007-126,
//! Division 3, ss. 16–17) and Part 3 of Transport Canada's *Standard for the
//! Tonnage Measurement of Vessels* (TP 13430).
//!
//! All linear dimensions enter in metres rounded to the nearest centimetre
//! (`0.01 m`), and intermediate and final tonnages are rounded to two decimal
//! places with a third decimal digit of `5` or more rounded up (TP 13430
//! s. 3.3.1). On 3D geometry, `TML` is measured parallel to the design
//! waterline across hull members only (propulsion mounts excluded, s. 3.2.8),
//! and `TMB` and `TMD` are taken on the transverse section at the middle of
//! `TML` (ss. 3.2.6–3.2.7), with `TMD` measured to the side edge of a cambered
//! upper deck rather than its crown (s. 3.2.7.1, Fig. 3.3).

use crate::params::CaseParams;
use crate::{setup, LoftRequest};
use hullgeom::iges::Platform;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

/// Round to two decimal places (`0.01`), increasing the second decimal digit
/// by one when the third is `5` or more (TP 13430 s. 3.3.1).
///
/// A tiny relative guard (`1e-11`) prevents binary floating-point downward
/// representation bias on exact decimal half-way values such as `1.005` or
/// `3.30 * 0.95 = 3.135`.
pub fn round2(x: f64) -> f64 {
    if !x.is_finite() || x == 0.0 {
        return 0.0;
    }
    let eps = 1e-11 * x.abs().max(1.0);
    if x >= 0.0 {
        ((x * 100.0) + 0.5 + eps).floor() / 100.0
    } else {
        -(((-x * 100.0) + 0.5 + eps).floor() / 100.0)
    }
}

/// Hull form class for the Gross Tonnage Coefficient (`GTC`, TP 13430 s. 3.2.3).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum HullKind {
    /// Vessel not designed for sailing (normal hull form): `GTC = 0.16`.
    #[default]
    #[serde(rename = "power")]
    Power,
    /// Vessel designed for sailing (with or without auxiliary engine): `GTC = 0.08`.
    #[serde(rename = "sail")]
    Sail,
    /// Barge hull form (rectangular box shape): `GTC = 0.20`.
    #[serde(rename = "barge")]
    Barge,
}

impl HullKind {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "power" => Ok(Self::Power),
            "sail" => Ok(Self::Sail),
            "barge" => Ok(Self::Barge),
            other => Err(format!("--kind {other:?}: expected power, sail or barge")),
        }
    }

    pub fn gtc(self) -> f64 {
        match self {
            Self::Power => 0.16,
            Self::Sail => 0.08,
            Self::Barge => 0.20,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Power => "power",
            Self::Sail => "sail",
            Self::Barge => "barge",
        }
    }
}

/// Propulsion class for the Net Tonnage Coefficient (`NTC`, TP 13430 s. 3.2.4).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Propulsion {
    /// Power-driven vessel: `NTC = 0.75`.
    #[serde(rename = "power")]
    Power,
    /// Vessel propelled by sail (alone or with auxiliary engine): `NTC = 0.95`.
    #[serde(rename = "sail")]
    Sail,
    /// Non-propelled vessel: `NTC = 1.00`.
    #[serde(rename = "none")]
    None,
}

impl Propulsion {
    pub fn default_for(kind: HullKind) -> Self {
        match kind {
            HullKind::Power => Self::Power,
            HullKind::Sail => Self::Sail,
            HullKind::Barge => Self::None,
        }
    }

    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "power" => Ok(Self::Power),
            "sail" => Ok(Self::Sail),
            "none" | "unpropelled" => Ok(Self::None),
            other => Err(format!(
                "--propulsion {other:?}: expected power, sail or none"
            )),
        }
    }

    pub fn ntc(self) -> f64 {
        match self {
            Self::Power => 0.75,
            Self::Sail => 0.95,
            Self::None => 1.00,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Power => "power",
            Self::Sail => "sail",
            Self::None => "none",
        }
    }
}

/// Mean length, mean breadth and mean height [m] of an enclosed superstructure,
/// deckhouse or break above the upper deck (s. 3.5.3.1), or of an enclosed
/// multi-hull bridge structure below the upper deck (s. 3.5.2.2).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Space {
    pub ml: f64,
    pub mb: f64,
    pub mh: f64,
}

impl Space {
    fn canonical(self, what: &str) -> Result<Self, String> {
        for (k, v) in [("ml", self.ml), ("mb", self.mb), ("mh", self.mh)] {
            if !(v > 0.0 && v.is_finite()) {
                return Err(format!("{what} {k} must be positive"));
            }
        }
        let s = Space {
            ml: round2(self.ml),
            mb: round2(self.mb),
            mh: round2(self.mh),
        };
        if s.ml <= 0.0 || s.mb <= 0.0 || s.mh <= 0.0 {
            return Err(format!("{what} dimensions must be at least 0.01 m"));
        }
        Ok(s)
    }

    /// `(ML * MB * MH) / 2.83` rounded per s. 3.3.1.
    pub fn tonnage(self) -> f64 {
        round2((self.ml * self.mb * self.mh) / 2.83)
    }
}

/// Principal dimensions of a single hull [m], optionally accompanied by the
/// 3D coordinates they were measured from (z up from the design waterline).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct HullDimensions {
    pub tml: f64,
    pub tmb: f64,
    pub tmd: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub x_aft: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub x_fore: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub x_mid: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub keel_z: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deck_z: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub l_wl85: Option<f64>,
}

impl HullDimensions {
    pub fn from_dims(tml: f64, tmb: f64, tmd: f64) -> Self {
        Self {
            tml,
            tmb,
            tmd,
            x_aft: None,
            x_fore: None,
            x_mid: None,
            keel_z: None,
            deck_z: None,
            l_wl85: None,
        }
    }
}

/// Options and overrides for a Part 3 tonnage calculation.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TonnageParams {
    pub kind: HullKind,
    /// Defaults from `kind` when `None` (`power` -> `power`, `sail` -> `sail`,
    /// `barge` -> `none`).
    pub propulsion: Option<Propulsion>,
    /// Override (or outright) Tonnage Measurement Length [m].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tml: Option<f64>,
    /// Override (or outright) Tonnage Measurement Breadth [m].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tmb: Option<f64>,
    /// Override (or outright) Tonnage Measurement Depth [m].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tmd: Option<f64>,
    /// Override upper-deck height at side at mid-TML [m, z up from design waterline].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deck_z: Option<f64>,
    /// Override bottom-of-keel height at mid-TML [m, z up from design waterline].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub keel_z: Option<f64>,
    /// Override statutory length L [m] (96% of 85%-depth waterline, or stem to
    /// rudder stock; VRTR s. 6 / TP 13430 Part 1 s. 1.1).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reg_length: Option<f64>,
    /// Enclosed spaces (and side-to-side breaks) above the upper deck (s. 3.5.3.1).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub spaces: Vec<Space>,
    /// Number of tiers of superstructures/deckhouses above the upper deck.
    pub tiers: usize,
    /// Enclosed bridge structures connecting hulls below the upper deck on a
    /// multi-hull vessel (s. 3.5.2.2).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub bridges: Vec<Space>,
}

impl TonnageParams {
    pub fn canonical(mut self) -> Result<Self, String> {
        self.propulsion = Some(
            self.propulsion
                .unwrap_or_else(|| Propulsion::default_for(self.kind)),
        );
        for (name, val) in [
            ("tml", self.tml),
            ("tmb", self.tmb),
            ("tmd", self.tmd),
            ("reg-length", self.reg_length),
        ] {
            if let Some(x) = val {
                if !(x > 0.0 && x.is_finite()) {
                    return Err(format!("{name} must be positive"));
                }
                if round2(x) <= 0.0 {
                    return Err(format!("{name} must be at least 0.01 m"));
                }
            }
        }
        for (name, val) in [("deck-z", self.deck_z), ("keel-z", self.keel_z)] {
            if let Some(x) = val {
                if !x.is_finite() {
                    return Err(format!("{name}: not a number"));
                }
            }
        }
        if self.tmd.is_some() && (self.deck_z.is_some() || self.keel_z.is_some()) {
            return Err("give --tmd, or --deck-z/--keel-z, not both".into());
        }
        if let (Some(dz), Some(kz)) = (self.deck_z, self.keel_z) {
            if dz <= kz {
                return Err(format!("deck-z ({dz} m) must be above keel-z ({kz} m)"));
            }
        }
        self.tml = self.tml.map(round2);
        self.tmb = self.tmb.map(round2);
        self.tmd = self.tmd.map(round2);
        self.reg_length = self.reg_length.map(round2);
        self.spaces = self
            .spaces
            .into_iter()
            .map(|s| s.canonical("space"))
            .collect::<Result<Vec<_>, _>>()?;
        self.bridges = self
            .bridges
            .into_iter()
            .map(|b| b.canonical("bridge"))
            .collect::<Result<Vec<_>, _>>()?;
        if !self.spaces.is_empty() && self.tiers == 0 {
            self.tiers = 1;
        }
        Ok(self)
    }
}

/// Compute Transport Canada Part 3 gross and net tonnage from measured or
/// outright hull dimensions.
pub fn compute(
    hulls: &[HullDimensions],
    params: &TonnageParams,
    mut notes: Vec<String>,
) -> Result<Value, String> {
    if hulls.is_empty() {
        return Err("no hulls to measure".into());
    }
    let p = params.clone().canonical()?;
    let propulsion = p
        .propulsion
        .unwrap_or_else(|| Propulsion::default_for(p.kind));
    if hulls.len() == 1 && !p.bridges.is_empty() {
        return Err("--bridge is for multi-hull vessels (TP 13430 s. 3.5.2.2)".into());
    }

    let gtc = p.kind.gtc();
    let ntc = propulsion.ntc();

    let mut hull_entries = Vec::with_capacity(hulls.len());
    let mut hulls_tonnage_sum = 0.0;
    let mut max_tml = 0.0f64;
    let mut max_tmb = 0.0f64;
    let mut max_tmd = 0.0f64;

    for (idx, h) in hulls.iter().enumerate() {
        let tml = round2(p.tml.unwrap_or(h.tml));
        let tmb = round2(p.tmb.unwrap_or(h.tmb));
        let tmd = round2(p.tmd.unwrap_or(h.tmd));
        for (what, v) in [("TML", tml), ("TMB", tmb), ("TMD", tmd)] {
            if !(v > 0.0 && v.is_finite()) {
                return Err(format!(
                    "hull {}: {what} ({v} m) must be at least 0.01 m",
                    idx + 1
                ));
            }
        }
        max_tml = max_tml.max(tml);
        max_tmb = max_tmb.max(tmb);
        max_tmd = max_tmd.max(tmd);
        let ht = round2(tml * tmb * tmd * gtc);
        hulls_tonnage_sum = round2(hulls_tonnage_sum + ht);

        let mut obj = json!({
            "tml": tml,
            "tmb": tmb,
            "tmd": tmd,
            "tonnage": ht,
        });
        if let Some(o) = obj.as_object_mut() {
            for (k, v) in [
                ("x_aft", h.x_aft),
                ("x_fore", h.x_fore),
                ("x_mid", h.x_mid),
                ("keel_z", h.keel_z),
                ("deck_z", h.deck_z),
                ("l_wl85", h.l_wl85),
            ] {
                if let Some(val) = v {
                    o.insert(k.into(), json!(round_geom(val)));
                }
            }
        }
        hull_entries.push(obj);
    }

    // Statutory length L check: Part 3 applies strictly to vessels < 24.0 m.
    let reg_length: Option<f64> = if let Some(rl) = p.reg_length {
        Some(rl)
    } else if p.tml.is_none() && hulls.iter().all(|h| h.l_wl85.is_some()) {
        let l85 = hulls.iter().filter_map(|h| h.l_wl85).fold(0.0f64, f64::max);
        Some(round2(0.96 * l85))
    } else {
        None
    };

    if let Some(rl) = reg_length {
        if rl >= 24.0 {
            return Err(format!(
                "statutory length {rl:.2} m is 24 m or more: TP 13430 Part 3 applies only to vessels under 24 m (use Part 2 / ITC 69)"
            ));
        }
    } else if max_tml >= 24.0 {
        return Err(format!(
            "TML {max_tml:.2} m is 24 m or more: TP 13430 Part 3 applies only to vessels under 24 m in statutory length (pass --reg-length if 96% of the 85%-depth waterline is under 24 m)"
        ));
    }

    // Multi-hull enclosed bridge structures below upper deck (s. 3.5.2.2).
    let mut bridge_entries = Vec::with_capacity(p.bridges.len());
    let mut bridge_tonnage = 0.0;
    for b in &p.bridges {
        let bt = b.tonnage();
        bridge_tonnage = round2(bridge_tonnage + bt);
        bridge_entries.push(json!({
            "ml": b.ml,
            "mb": b.mb,
            "mh": b.mh,
            "tonnage": bt,
        }));
    }

    // Per s. 3.5.2.2, enclosed bridge structures connecting hulls below the
    // upper deck form part of the tonnage of the hull.
    let hull_tonnage = round2(hulls_tonnage_sum + bridge_tonnage);

    // Small monohull exemption for spaces above the upper deck (s. 3.5.3.2):
    // monohull, TML <= 15.00 m, at most 1 tier of superstructures/deckhouses
    // (including breaks), and combined ML <= 70% of TML. Evaluated in integer
    // centimetres so exact boundary cases such as 70% of 10.10 m = 7.07 m never
    // suffer binary float rounding drift.
    let tml_cm = (max_tml * 100.0).round() as i64;
    let sum_ml_cm: i64 = p.spaces.iter().map(|s| (s.ml * 100.0).round() as i64).sum();
    let exempt_eligible =
        hulls.len() == 1 && tml_cm <= 1500 && p.tiers <= 1 && sum_ml_cm * 10 <= tml_cm * 7;
    let above_deck_exempt = exempt_eligible && !p.spaces.is_empty();

    let mut space_entries = Vec::with_capacity(p.spaces.len());
    let mut raw_spaces_tonnage = 0.0;
    for s in &p.spaces {
        let st = s.tonnage();
        raw_spaces_tonnage = round2(raw_spaces_tonnage + st);
        space_entries.push(json!({
            "ml": s.ml,
            "mb": s.mb,
            "mh": s.mh,
            "tonnage": st,
            "exempt": exempt_eligible,
        }));
    }

    let above_deck_tonnage = if exempt_eligible {
        if !p.spaces.is_empty() {
            notes.push(
                "spaces above upper deck omitted under TP 13430 s. 3.5.3.2 (monohull, TML <= 15 m, <= 1 tier, total ML <= 70% TML)"
                    .into(),
            );
        }
        0.0
    } else {
        raw_spaces_tonnage
    };

    let gt = round2(hull_tonnage + above_deck_tonnage);
    let nt = round2(gt * ntc);

    // Owner-measurement eligibility (TP 13430 s. 3.4.2, Form 4A / Form 4B).
    let form: Option<&'static str> = if exempt_eligible {
        Some("4A")
    } else if hulls.len() > 1
        && tml_cm <= 1500
        && p.bridges.is_empty()
        && p.spaces.is_empty()
        && p.tiers == 0
    {
        Some("4B")
    } else {
        None
    };

    Ok(json!({
        "regulation": "TP 13430 Part 3",
        "kind": p.kind.as_str(),
        "propulsion": propulsion.as_str(),
        "gtc": gtc,
        "ntc": ntc,
        "hull_count": hulls.len(),
        "tml": max_tml,
        "tmb": max_tmb,
        "tmd": max_tmd,
        "reg_length": reg_length,
        "gt": gt,
        "nt": nt,
        "hull_tonnage": hull_tonnage,
        "bridge_tonnage": bridge_tonnage,
        "above_deck_tonnage": above_deck_tonnage,
        "above_deck_exempt": above_deck_exempt,
        "owner_measurable": form.is_some(),
        "form": form,
        "hulls": hull_entries,
        "bridges": bridge_entries,
        "spaces": space_entries,
        "tiers": p.tiers,
        "notes": notes,
    }))
}

fn round_geom(x: f64) -> f64 {
    if !x.is_finite() || x == 0.0 {
        0.0
    } else {
        (x * 1e6).round() / 1e6
    }
}

/// Measure a hull file (or native JSON geometry) at its design pose and compute
/// Transport Canada Part 3 tonnage.
pub fn hull_tonnage(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    params: &TonnageParams,
) -> Result<Value, String> {
    case_tonnage(name, bytes, cut, &CaseParams::default(), params)
}

/// Measure a case's hull platform (doubling a catamaran's demihull when `span`
/// is set and ignoring propulsion mount parts per TP 13430 s. 3.2.8) and
/// compute Transport Canada Part 3 tonnage.
pub fn case_tonnage(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    case: &CaseParams,
    params: &TonnageParams,
) -> Result<Value, String> {
    let p = params.clone().canonical()?;
    let s = setup(name, bytes, cut, case)?;
    let mut notes = s.cut.notes.clone();
    let mut measured = Vec::new();
    let wl = s.cut.file.waterline_z;
    let at_rest = Platform::default();

    for (m, &(src_idx, ref pose)) in s.layout.iter().enumerate() {
        if s.roles.get(m).map(|r| r.0.as_str()) != Some("hull") {
            continue;
        }
        let cut_pos = s
            .cut
            .index
            .iter()
            .position(|&j| j == src_idx)
            .ok_or("missing hull in cut")?;
        let two_sided = s.cut.hulls[cut_pos].report.two_sided;
        let y_c = s.design[m].1.y;
        let (verts, tris) = s
            .cut
            .file
            .source
            .posed_tessellation(src_idx, wl, pose, &at_rest)
            .map_err(|e| e.to_string())?;
        let (dims, hull_notes) =
            measure_mesh(&verts, &tris, two_sided, y_c, &p, measured.len() + 1)?;
        for n in hull_notes {
            if !notes.contains(&n) {
                notes.push(n);
            }
        }
        measured.push(dims);
    }

    if measured.is_empty() {
        return Err("the platform has no hull geometry".into());
    }
    compute(&measured, &p, notes)
}

/// Measure one hull's posed tessellation (`x` forward parallel to the design
/// waterline, `y` transverse, `z` up from the design waterline).
fn measure_mesh(
    verts: &[[f64; 3]],
    tris: &[[u32; 3]],
    two_sided: bool,
    y_c: f64,
    params: &TonnageParams,
    hull_num: usize,
) -> Result<(HullDimensions, Vec<String>), String> {
    let mut x_aft = f64::INFINITY;
    let mut x_fore = f64::NEG_INFINITY;
    let mut z_mesh_max = f64::NEG_INFINITY;

    for v in verts {
        if v[0].is_finite() && v[1].is_finite() && v[2].is_finite() {
            x_aft = x_aft.min(v[0]);
            x_fore = x_fore.max(v[0]);
            z_mesh_max = z_mesh_max.max(v[2]);
        }
    }
    if !(x_fore > x_aft && x_aft.is_finite() && x_fore.is_finite()) {
        return Err(format!("hull {hull_num}: degenerate longitudinal extent"));
    }

    let raw_tml = x_fore - x_aft;
    let x_mid = 0.5 * (x_aft + x_fore);
    let eps = 1e-12 * raw_tml.max(1.0);

    let segs = slice_mesh_x(verts, tris, x_mid, eps);
    if segs.is_empty() {
        return Err(format!(
            "hull {hull_num}: no cross-section geometry at mid-TML (x = {x_mid:.3} m)"
        ));
    }

    // TMB (s. 3.2.6, Fig. 3.2): maximum horizontal breadth at mid-TML.
    let raw_tmb = if two_sided {
        let mut y_lo = f64::INFINITY;
        let mut y_hi = f64::NEG_INFINITY;
        for &(a, b) in &segs {
            y_lo = y_lo.min(a[0]).min(b[0]);
            y_hi = y_hi.max(a[0]).max(b[0]);
        }
        y_hi - y_lo
    } else {
        let mut h_max = 0.0f64;
        for &(a, b) in &segs {
            h_max = h_max.max((a[0] - y_c).abs()).max((b[0] - y_c).abs());
        }
        2.0 * h_max
    };

    if !(raw_tmb > 0.0 && raw_tmb.is_finite()) {
        return Err(format!(
            "hull {hull_num}: zero breadth at mid-TML (x = {x_mid:.3} m)"
        ));
    }

    // Bottom of keel and upper deck at side at mid-TML (s. 3.2.7, Fig. 3.3).
    let mut z_bot = f64::INFINITY;
    let mut z_top = f64::NEG_INFINITY;
    for &(a, b) in &segs {
        z_bot = z_bot.min(a[1]).min(b[1]);
        z_top = z_top.max(a[1]).max(b[1]);
    }

    let mut notes = Vec::new();
    let keel_z = match params.keel_z {
        Some(kz) => kz,
        None => {
            if params.tmd.is_none() {
                if let Some((z_aft, z_fore)) = adjacent_keel_z(verts, tris, x_mid, raw_tml, eps) {
                    let z_proj = 0.5 * (z_aft + z_fore);
                    // Only flag when BOTH neighbouring stations are lower than
                    // mid-TML and the projected keel line sits > 2 cm deeper.
                    if z_aft < z_bot - 0.01 && z_fore < z_bot - 0.01 && z_bot - z_proj > 0.02 {
                        notes.push(format!(
                            "hull {hull_num}: keel at mid-TML (z = {z_bot:.3} m) sits above adjacent stations ({z_proj:.3} m); pass --keel-z if bridging a localized keel cutout (SSB 17/2025)"
                        ));
                    }
                }
            }
            z_bot
        }
    };

    let deck_z = match params.deck_z {
        Some(dz) => dz,
        None => detect_deck_side_z(&segs, two_sided, y_c, z_bot, z_top),
    };

    let raw_tmd = deck_z - keel_z;
    if params.tmd.is_none() && !(raw_tmd > 0.0 && raw_tmd.is_finite()) {
        return Err(format!(
            "hull {hull_num}: deck height ({deck_z:.3} m) is not above keel ({keel_z:.3} m) at mid-TML"
        ));
    }

    // Waterline at 85% of depth from top/bottom of keel for statutory length L
    // (TP 13430 Part 1 s. 1.1 / VRTR s. 6).
    let depth_for_wl = params.tmd.unwrap_or(raw_tmd);
    let z_wl85 = keel_z + 0.85 * depth_for_wl;
    let l_wl85 = if z_wl85 <= z_mesh_max - 1e-6 * depth_for_wl.max(1.0) {
        waterline_length_z(verts, tris, z_wl85, eps).unwrap_or(raw_tml)
    } else {
        raw_tml
    };

    Ok((
        HullDimensions {
            tml: raw_tml,
            tmb: raw_tmb,
            tmd: raw_tmd.max(0.0),
            x_aft: Some(x_aft),
            x_fore: Some(x_fore),
            x_mid: Some(x_mid),
            keel_z: Some(keel_z),
            deck_z: Some(deck_z),
            l_wl85: Some(l_wl85),
        },
        notes,
    ))
}

fn adjacent_keel_z(
    verts: &[[f64; 3]],
    tris: &[[u32; 3]],
    x_mid: f64,
    tml: f64,
    eps: f64,
) -> Option<(f64, f64)> {
    let dx = 0.05 * tml;
    let min_z = |x: f64| {
        let s = slice_mesh_x(verts, tris, x, eps);
        if s.is_empty() {
            None
        } else {
            Some(
                s.iter()
                    .fold(f64::INFINITY, |m, &(a, b)| m.min(a[1]).min(b[1])),
            )
        }
    };
    Some((min_z(x_mid - dx)?, min_z(x_mid + dx)?))
}

/// Intersect triangles with the transverse plane `x = x_s`, returning `(y, z)`
/// segments on the shell while skipping transverse diaphragm triangles that lie
/// wholly in `x = x_s`.
fn slice_mesh_x(
    verts: &[[f64; 3]],
    tris: &[[u32; 3]],
    x_s: f64,
    eps: f64,
) -> Vec<([f64; 2], [f64; 2])> {
    let mut out = Vec::new();
    for tri in tris {
        let p = [
            verts[tri[0] as usize],
            verts[tri[1] as usize],
            verts[tri[2] as usize],
        ];
        let d = [p[0][0] - x_s, p[1][0] - x_s, p[2][0] - x_s];
        let on: Vec<usize> = (0..3).filter(|&i| d[i].abs() <= eps).collect();
        match on.len() {
            3 => continue, // transverse bulkhead / end-cap face in the plane
            2 => {
                let (a, b) = (p[on[0]], p[on[1]]);
                out.push(([a[1], a[2]], [b[1], b[2]]));
            }
            1 => {
                let k = on[0];
                let (a, b) = ((k + 1) % 3, (k + 2) % 3);
                if (d[a] > 0.0) != (d[b] > 0.0) {
                    let t = -d[a] / (d[b] - d[a]);
                    let q = [
                        p[a][1] + t * (p[b][1] - p[a][1]),
                        p[a][2] + t * (p[b][2] - p[a][2]),
                    ];
                    out.push(([p[k][1], p[k][2]], q));
                }
            }
            0 => {
                if (d[0] > 0.0) == (d[1] > 0.0) && (d[1] > 0.0) == (d[2] > 0.0) {
                    continue;
                }
                let mut pts = Vec::with_capacity(2);
                for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                    if (d[a] > 0.0) != (d[b] > 0.0) {
                        let t = -d[a] / (d[b] - d[a]);
                        pts.push([
                            p[a][1] + t * (p[b][1] - p[a][1]),
                            p[a][2] + t * (p[b][2] - p[a][2]),
                        ]);
                    }
                }
                if pts.len() == 2 {
                    out.push((pts[0], pts[1]));
                }
            }
            _ => unreachable!(),
        }
    }
    out
}

/// Longitudinal span `x_max - x_min` on the horizontal plane `z = z_s`.
fn waterline_length_z(verts: &[[f64; 3]], tris: &[[u32; 3]], z_s: f64, eps: f64) -> Option<f64> {
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut any = false;
    let mut mark = |x: f64| {
        x_lo = x_lo.min(x);
        x_hi = x_hi.max(x);
        any = true;
    };
    for tri in tris {
        let p = [
            verts[tri[0] as usize],
            verts[tri[1] as usize],
            verts[tri[2] as usize],
        ];
        let d = [p[0][2] - z_s, p[1][2] - z_s, p[2][2] - z_s];
        for i in 0..3 {
            if d[i].abs() <= eps {
                mark(p[i][0]);
            }
        }
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            if (d[a] > eps && d[b] < -eps) || (d[a] < -eps && d[b] > eps) {
                let t = -d[a] / (d[b] - d[a]);
                mark(p[a][0] + t * (p[b][0] - p[a][0]));
            }
        }
    }
    (any && x_hi > x_lo).then_some(x_hi - x_lo)
}

/// Upper-deck height at the side of the vessel (TP 13430 s. 3.2.7.1, Fig. 3.3).
///
/// Open-top hulls have no deck closing onto `y_c` near the top, so `z_top` is
/// the sheer/gunwale height (`s. 3.2.7.3`). Closed-top hulls with a cambered
/// deck start near `(h = 0, z_crown)` and slope gently outward (`-dz/dh <= 0.40`)
/// to a sheer chine at `h >= 0.65 * h_max`, where the shell turns sharply down
/// into the hull side (`-dz/dh > 0.8`) or reaches `h_max`; circular pontoons
/// curve past `-dz/dh = 0.40` well before `0.65 R`, keeping their full diameter.
fn detect_deck_side_z(
    segs: &[([f64; 2], [f64; 2])],
    two_sided: bool,
    y_c: f64,
    z_bot: f64,
    z_top: f64,
) -> f64 {
    let total_h = z_top - z_bot;
    if total_h <= 1e-9 {
        return z_top;
    }
    if two_sided {
        let pos = side_deck_z(segs, y_c, 1.0, z_bot, z_top);
        let neg = side_deck_z(segs, y_c, -1.0, z_bot, z_top);
        match (pos, neg) {
            (Some(a), Some(b)) => 0.5 * (a + b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => z_top,
        }
    } else {
        // One-sided mesh lies on whichever side of y_c carries half-breadth.
        let mut h_pos = 0.0f64;
        let mut h_neg = 0.0f64;
        for &(a, b) in segs {
            h_pos = h_pos.max(a[0] - y_c).max(b[0] - y_c);
            h_neg = h_neg.max(y_c - a[0]).max(y_c - b[0]);
        }
        let dir = if h_pos >= h_neg { 1.0 } else { -1.0 };
        side_deck_z(segs, y_c, dir, z_bot, z_top).unwrap_or(z_top)
    }
}

fn side_deck_z(
    segs: &[([f64; 2], [f64; 2])],
    y_c: f64,
    dir: f64,
    z_bot: f64,
    z_top: f64,
) -> Option<f64> {
    let total_h = z_top - z_bot;
    // Collect segments clipped to the half-plane `s = (y - y_c) * dir >= 0`,
    // oriented with `0 <= h0 <= h1`.
    let mut half_segs: Vec<([f64; 2], [f64; 2])> = Vec::new();
    let mut h_max = 0.0f64;
    for &(a, b) in segs {
        let mut sa = (a[0] - y_c) * dir;
        let mut za = a[1];
        let mut sb = (b[0] - y_c) * dir;
        let mut zb = b[1];
        if sa < -1e-12 && sb < -1e-12 {
            continue;
        }
        if sa < 0.0 && sb > 0.0 {
            let t = -sa / (sb - sa);
            za = za + t * (zb - za);
            sa = 0.0;
        } else if sb < 0.0 && sa > 0.0 {
            let t = -sb / (sa - sb);
            zb = zb + t * (za - zb);
            sb = 0.0;
        }
        sa = sa.max(0.0);
        sb = sb.max(0.0);
        h_max = h_max.max(sa).max(sb);
        if sa <= sb {
            half_segs.push(([sa, za], [sb, zb]));
        } else {
            half_segs.push(([sb, zb], [sa, za]));
        }
    }
    if h_max <= 1e-9 {
        return None;
    }

    let h_tol = 1e-4 * h_max;
    let non_vert: Vec<([f64; 2], [f64; 2])> = half_segs
        .iter()
        .copied()
        .filter(|&(a, b)| b[0] - a[0] > 1e-7 * h_max && a[1].max(b[1]) >= z_bot + 0.5 * total_h)
        .collect();

    // Must have an upper deck meeting the centreplane in the upper 40% of the section.
    let start = non_vert
        .iter()
        .filter(|&&(a, _)| a[0] <= 0.05 * h_max && a[1] >= z_bot + 0.6 * total_h)
        .max_by(|a, b| a.0[1].total_cmp(&b.0[1]))?;

    let mut h_cur = start.0[0];
    let mut z_cur = start.0[1];
    let mut steep_break = false;
    let mut at_edge = false;

    for _ in 0..=non_vert.len() {
        // Among segments covering h_cur and extending outboard past it, follow
        // the uppermost one.
        let next = non_vert
            .iter()
            .filter(|&&(a, b)| a[0] <= h_cur + h_tol && b[0] > h_cur + h_tol)
            .max_by(|&&(a1, b1), &&(a2, b2)| {
                let t1 = ((h_cur - a1[0]) / (b1[0] - a1[0])).clamp(0.0, 1.0);
                let t2 = ((h_cur - a2[0]) / (b2[0] - a2[0])).clamp(0.0, 1.0);
                let z1 = a1[1] + t1 * (b1[1] - a1[1]);
                let z2 = a2[1] + t2 * (b2[1] - a2[1]);
                z1.total_cmp(&z2)
            });

        let Some(&(a, b)) = next else {
            at_edge = true;
            break;
        };
        let slope = (a[1] - b[1]) / (b[0] - a[0]);
        if (-0.05..=0.40).contains(&slope) {
            h_cur = b[0];
            z_cur = b[1];
        } else {
            steep_break = slope > 0.8;
            break;
        }
    }

    if h_cur >= 0.65 * h_max && (at_edge || steep_break) {
        Some(z_cur)
    } else {
        Some(z_top)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LoftRequest;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "got {a}, expected {b}");
    }

    #[test]
    fn round2_follows_tp13430_s_3_3_1() {
        close(round2(0.0), 0.0);
        close(round2(1.004), 1.00);
        close(round2(1.005), 1.01);
        close(round2(2.675), 2.68);
        // 3.30 * 0.95 = 3.135 exact in decimal, 3.1349999999999998 in IEEE-754 f64.
        close(round2(3.30 * 0.95), 3.14);
        close(round2(12.3449), 12.34);
        close(round2(12.345), 12.35);
    }

    #[test]
    fn monohull_coefficients_and_rounding() {
        // 12.00 x 3.60 x 1.80 = 77.76 m^3:
        // power: GT = 77.76 * 0.16 = 12.4416 -> 12.44, NT = 12.44 * 0.75 = 9.33
        let dims = [HullDimensions::from_dims(12.0, 3.6, 1.8)];
        let p = compute(&dims, &TonnageParams::default(), vec![]).unwrap();
        close(p["gt"].as_f64().unwrap(), 12.44);
        close(p["nt"].as_f64().unwrap(), 9.33);
        assert_eq!(p["form"], "4A");
        assert_eq!(p["owner_measurable"], true);

        // sail: GT = 77.76 * 0.08 = 6.2208 -> 6.22, NT = 6.22 * 0.95 = 5.909 -> 5.91
        let s = compute(
            &dims,
            &TonnageParams {
                kind: HullKind::Sail,
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        close(s["gt"].as_f64().unwrap(), 6.22);
        close(s["nt"].as_f64().unwrap(), 5.91);

        // barge (default unpropelled): GT = 77.76 * 0.20 = 15.552 -> 15.55, NT = 15.55
        let b = compute(
            &dims,
            &TonnageParams {
                kind: HullKind::Barge,
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        close(b["gt"].as_f64().unwrap(), 15.55);
        close(b["nt"].as_f64().unwrap(), 15.55);

        // self-propelled barge: GTC = 0.20, NTC = 0.75 -> NT = 15.55 * 0.75 = 11.6625 -> 11.66
        let spb = compute(
            &dims,
            &TonnageParams {
                kind: HullKind::Barge,
                propulsion: Some(Propulsion::Power),
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        close(spb["gt"].as_f64().unwrap(), 15.55);
        close(spb["nt"].as_f64().unwrap(), 11.66);
    }

    #[test]
    fn small_monohull_exemption_s_3_5_3_2_boundaries() {
        // Exact 70% boundary at TML = 10.10 m (where 0.7 * 10.10 = 7.069999999999999 in f64):
        // ML = 7.07 m MUST be exempt!
        let dims = [HullDimensions::from_dims(10.10, 3.20, 1.50)];
        let exempt = compute(
            &dims,
            &TonnageParams {
                spaces: vec![Space {
                    ml: 7.07,
                    mb: 2.40,
                    mh: 1.90,
                }],
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        assert_eq!(exempt["above_deck_exempt"], true);
        close(exempt["above_deck_tonnage"].as_f64().unwrap(), 0.0);
        assert_eq!(exempt["form"], "4A");

        // Exceeding 70% by 1 cm (ML = 7.08 m on 10.10 m): NOT exempt!
        let over_len = compute(
            &dims,
            &TonnageParams {
                spaces: vec![Space {
                    ml: 7.08,
                    mb: 2.40,
                    mh: 1.90,
                }],
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        assert_eq!(over_len["above_deck_exempt"], false);
        let expected_space = round2((7.08 * 2.40 * 1.90) / 2.83);
        close(
            over_len["above_deck_tonnage"].as_f64().unwrap(),
            expected_space,
        );
        assert!(over_len["form"].is_null());

        // Two tiers on <= 15 m monohull: NOT exempt!
        let two_tiers = compute(
            &dims,
            &TonnageParams {
                spaces: vec![Space {
                    ml: 4.0,
                    mb: 2.0,
                    mh: 1.9,
                }],
                tiers: 2,
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        assert_eq!(two_tiers["above_deck_exempt"], false);
        assert!(two_tiers["above_deck_tonnage"].as_f64().unwrap() > 0.0);

        // TML = 15.01 m (> 15.00 m): NOT exempt even with small single-tier cabin!
        let over_15 = compute(
            &[HullDimensions::from_dims(15.01, 4.20, 2.10)],
            &TonnageParams {
                spaces: vec![Space {
                    ml: 3.0,
                    mb: 2.5,
                    mh: 2.0,
                }],
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        assert_eq!(over_15["above_deck_exempt"], false);
        assert!(over_15["above_deck_tonnage"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn multihull_bridge_and_no_3_5_3_2_exemption() {
        let hulls = [
            HullDimensions::from_dims(11.50, 1.40, 1.80),
            HullDimensions::from_dims(11.50, 1.40, 1.80),
        ];
        // Bare catamaran <= 15 m: Form 4B
        let bare = compute(&hulls, &TonnageParams::default(), vec![]).unwrap();
        let single_ht = round2(11.50 * 1.40 * 1.80 * 0.16);
        close(bare["gt"].as_f64().unwrap(), round2(2.0 * single_ht));
        assert_eq!(bare["form"], "4B");

        // Catamaran with bridge structure and small deckhouse: s. 3.5.3.2 applies
        // to monohulls only, so both bridge and deckhouse are included!
        let bridged = compute(
            &hulls,
            &TonnageParams {
                kind: HullKind::Sail,
                bridges: vec![Space {
                    ml: 5.50,
                    mb: 3.20,
                    mh: 0.85,
                }],
                spaces: vec![Space {
                    ml: 4.20,
                    mb: 3.40,
                    mh: 1.95,
                }],
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        let sail_ht = round2(11.50 * 1.40 * 1.80 * 0.08);
        let bt = round2((5.50 * 3.20 * 0.85) / 2.83);
        let st = round2((4.20 * 3.40 * 1.95) / 2.83);
        close(
            bridged["hull_tonnage"].as_f64().unwrap(),
            round2(2.0 * sail_ht + bt),
        );
        close(bridged["above_deck_tonnage"].as_f64().unwrap(), st);
        close(
            bridged["gt"].as_f64().unwrap(),
            round2(2.0 * sail_ht + bt + st),
        );
        assert_eq!(bridged["above_deck_exempt"], false);
        assert!(bridged["form"].is_null());
    }

    #[test]
    fn statutory_length_24m_limit() {
        // Outright TML = 24.0 m without --reg-length is refused.
        let err = compute(
            &[HullDimensions::from_dims(24.0, 5.0, 2.5)],
            &TonnageParams::default(),
            vec![],
        )
        .unwrap_err();
        assert!(err.contains("24 m"), "{err}");

        // TML = 24.5 m with --reg-length 23.4 m (< 24 m) is accepted under Part 3!
        let ok = compute(
            &[HullDimensions::from_dims(24.5, 5.5, 2.8)],
            &TonnageParams {
                reg_length: Some(23.4),
                ..Default::default()
            },
            vec![],
        )
        .unwrap();
        close(ok["tml"].as_f64().unwrap(), 24.5);
        close(ok["reg_length"].as_f64().unwrap(), 23.4);

        // --reg-length 24.0 m is refused.
        assert!(compute(
            &[HullDimensions::from_dims(23.8, 5.5, 2.8)],
            &TonnageParams {
                reg_length: Some(24.0),
                ..Default::default()
            },
            vec![],
        )
        .is_err());
    }

    /// Build a closed triangular mesh of a prism along `x in [0, l]` with given
    /// closed 2D cross-section polygon `[(y, z), ...]` (counter-clockwise).
    fn prism_native_mesh(l: f64, section: &[[f64; 2]]) -> Vec<u8> {
        let n = section.len();
        let mut verts = Vec::with_capacity(2 * n + 2);
        for &[y, z] in section {
            verts.push([0.0, y, z]);
        }
        for &[y, z] in section {
            verts.push([l, y, z]);
        }
        // End-cap centres
        let yc = section.iter().map(|p| p[0]).sum::<f64>() / n as f64;
        let zc = section.iter().map(|p| p[1]).sum::<f64>() / n as f64;
        verts.push([0.0, yc, zc]);
        verts.push([l, yc, zc]);
        let c0 = (2 * n) as u32;
        let c1 = (2 * n + 1) as u32;

        let mut tris = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            let (a0, b0) = (i as u32, j as u32);
            let (a1, b1) = ((n + i) as u32, (n + j) as u32);
            tris.push([a0, a1, b1]);
            tris.push([a0, b1, b0]);
            tris.push([c0, b0, a0]);
            tris.push([c1, a1, b1]);
        }
        serde_json::to_vec(&json!({
            "kind": "mesh",
            "hulls": [{ "vertices": verts, "triangles": tris }]
        }))
        .unwrap()
    }

    #[test]
    fn box_barge_cambered_deck_and_cylinder_geometry() {
        // 1. Rectangular box barge: 12 m x 4 m x 2 m (draft 1 m, freeboard 1 m).
        let box_bytes =
            prism_native_mesh(12.0, &[[-2.0, -1.0], [2.0, -1.0], [2.0, 1.0], [-2.0, 1.0]]);
        let b = hull_tonnage(
            "barge.json",
            box_bytes,
            &LoftRequest::default(),
            &TonnageParams {
                kind: HullKind::Barge,
                ..Default::default()
            },
        )
        .unwrap();
        close(b["tml"].as_f64().unwrap(), 12.0);
        close(b["tmb"].as_f64().unwrap(), 4.0);
        close(b["tmd"].as_f64().unwrap(), 2.0);
        close(b["reg_length"].as_f64().unwrap(), round2(0.96 * 12.0));
        close(b["gt"].as_f64().unwrap(), round2(12.0 * 4.0 * 2.0 * 0.20));

        // 2. Hull with tumblehome (max beam 3.0 m at z = 0.4 m, narrowing to 2.8 m
        // at sheer z = 1.0 m) AND a cambered deck rising 0.18 m to z = 1.18 m at
        // centreline, keel at z = -0.80 m.
        // Per TP 13430 ss. 3.2.6 & 3.2.7.1:
        //   TMB = 3.00 m (maximum breadth at mid-TML, including tumblehome)
        //   TMD = 1.00 - (-0.80) = 1.80 m (to deck at side, NOT 1.98 m crown!)
        let mut sec = vec![
            [0.0, -0.80],
            [1.20, -0.40],
            [1.50, 0.40], // tumblehome max beam: half-breadth 1.50 -> TMB = 3.00
            [1.40, 1.00], // starboard sheer edge: z = 1.00
        ];
        // Parabolic deck camber from h = 1.40 down to -1.40 peaking at z = 1.18
        for k in (1..20).rev() {
            let y = -1.40 + 2.80 * (k as f64 / 20.0);
            let u = y / 1.40;
            let z = 1.18 - 0.18 * u * u;
            sec.push([y, z]);
        }
        sec.extend_from_slice(&[[-1.40, 1.00], [-1.50, 0.40], [-1.20, -0.40]]);
        let cambered_bytes = prism_native_mesh(10.0, &sec);
        let cam = hull_tonnage(
            "cambered.json",
            cambered_bytes,
            &LoftRequest::default(),
            &TonnageParams::default(),
        )
        .unwrap();
        close(cam["tmb"].as_f64().unwrap(), 3.00);
        close(cam["tmd"].as_f64().unwrap(), 1.80);

        // 3. Cylindrical pontoon (radius 0.50 m, centered at z = 0.05 m so bottom
        // is at -0.45 m and top at +0.55 m): TMD must be full diameter 1.00 m!
        let n_cyl = 72;
        let cyl_sec: Vec<[f64; 2]> = (0..n_cyl)
            .map(|i| {
                let th = 2.0 * std::f64::consts::PI * (i as f64 / n_cyl as f64);
                [0.50 * th.cos(), 0.05 + 0.50 * th.sin()]
            })
            .collect();
        let cyl_bytes = prism_native_mesh(8.0, &cyl_sec);
        let cyl = hull_tonnage(
            "pontoon.json",
            cyl_bytes,
            &LoftRequest::default(),
            &TonnageParams::default(),
        )
        .unwrap();
        close(cyl["tmb"].as_f64().unwrap(), 1.00);
        close(cyl["tmd"].as_f64().unwrap(), 1.00);
    }

    #[test]
    fn full_vs_half_hull_wigley_catamaran_and_mount_exclusion() {
        let full_surfs = hullgeom::iges::wigley_surfaces(10.0, 1.2, 1.0).unwrap();
        let full_igs = hullgeom::iges::write(&full_surfs, "wigley_full")
            .unwrap()
            .into_bytes();
        // Half-hull IGES: keep only the starboard patch (ctrl y >= 0).
        let half_surfs: Vec<_> = full_surfs
            .iter()
            .filter(|s| s.ctrl.iter().all(|p| p[1] >= -1e-9))
            .cloned()
            .collect();
        assert_eq!(half_surfs.len(), 1);
        let half_igs = hullgeom::iges::write(&half_surfs, "wigley_half")
            .unwrap()
            .into_bytes();

        let cut = LoftRequest {
            waterline: Some(-0.4),
            ..Default::default()
        };
        let t_full = hull_tonnage("full.igs", full_igs.clone(), &cut, &Default::default()).unwrap();
        let t_half = hull_tonnage("half.igs", half_igs.clone(), &cut, &Default::default()).unwrap();

        close(t_full["tml"].as_f64().unwrap(), 10.0);
        close(t_full["tmb"].as_f64().unwrap(), 1.2);
        close(t_full["tmd"].as_f64().unwrap(), 1.0);
        assert_eq!(t_full["gt"], t_half["gt"]);
        assert_eq!(t_full["tmb"], t_half["tmb"]);
        assert_eq!(t_full["tmd"], t_half["tmd"]);

        // Catamaran via CaseParams::span on both full and half-hull IGES:
        let cat_case = CaseParams {
            span: Some(3.0),
            ..Default::default()
        };
        let cat_full = case_tonnage(
            "full.igs",
            full_igs.clone(),
            &cut,
            &cat_case,
            &Default::default(),
        )
        .unwrap();
        let cat_half =
            case_tonnage("half.igs", half_igs, &cut, &cat_case, &Default::default()).unwrap();
        assert_eq!(cat_full["hull_count"], 2);
        assert_eq!(cat_full["gt"], cat_half["gt"]);
        close(
            cat_full["gt"].as_f64().unwrap(),
            round2(2.0 * t_full["gt"].as_f64().unwrap()),
        );

        // Outboard mount extending astern (x = -0.35 m) and saildrive hanging
        // below the keel MUST NOT increase TML or TMD (TP 13430 s. 3.2.8)!
        let outboard = crate::mount::build(
            None,
            &crate::mount::Options {
                kind: Some(crate::mount::Kind::Outboard),
                x: Some(-0.35),
                shaft_depth: Some(0.45),
                chord: Some(0.18),
                thickness: Some(0.03),
                pod_length: Some(0.40),
                pod_diameter: Some(0.10),
                nose_ahead: Some(0.15),
                prop_from_nose: Some(0.40),
                ..Default::default()
            },
        )
        .unwrap();
        let saildrive = crate::mount::build(
            Some("oceanvolt-sd8"),
            &crate::mount::Options {
                x: Some(2.0),
                ..Default::default()
            },
        )
        .unwrap();
        for m in [outboard, saildrive] {
            let mounted_case = CaseParams {
                mount: Some(m),
                ..Default::default()
            };
            let t_mounted = case_tonnage(
                "full.igs",
                full_igs.clone(),
                &cut,
                &mounted_case,
                &Default::default(),
            )
            .unwrap();
            assert_eq!(t_mounted["tml"], t_full["tml"]);
            assert_eq!(t_mounted["tmd"], t_full["tmd"]);
            assert_eq!(t_mounted["gt"], t_full["gt"]);
        }
    }
}
