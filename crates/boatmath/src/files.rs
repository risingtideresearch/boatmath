//! Reading hull files (IGES, STL) into source geometry the solvers re-pose
//! and re-cut ([`hullgeom::source::HullSource`]), and the small parsers that
//! go with it.

use hullgeom::iges::{source_fleet, SectionalOptions};
use hullgeom::source::HullSource;
use hullgeom::SectionalHull;

/// What a file held.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Iges,
    Stl,
}

/// Import settings shared by every command that reads a hull.
#[derive(Clone, Copy)]
pub struct LoadSettings {
    /// CAD height of the design waterline.
    pub waterline_z: f64,
    /// Explicit design waterline, if given (Camber documents use their own
    /// design waterline when this is None).
    pub waterline: Option<f64>,
    pub centerplane: Option<f64>,
    pub stations: usize,
    pub rays: usize,
    /// Scale to metres for STL, which carries no units.
    pub units: Option<f64>,
}

impl Default for LoadSettings {
    fn default() -> Self {
        let d = SectionalOptions::default();
        LoadSettings {
            waterline_z: 0.0,
            waterline: None,
            centerplane: None,
            stations: d.stations,
            rays: d.rays,
            units: None,
        }
    }
}

impl LoadSettings {
    /// Sectioning options for a file whose design waterline is at
    /// `waterline_z`.
    pub fn sectional(&self, waterline_z: f64) -> SectionalOptions {
        SectionalOptions {
            waterline_z,
            centerplane: self.centerplane,
            stations: self.stations,
            rays: self.rays,
            ..Default::default()
        }
    }
}

/// A parsed file: its hulls as re-cuttable geometry.
pub struct SourceFile {
    pub path: String,
    pub kind: Kind,
    pub source: Box<dyn HullSource>,
    /// CAD height of the design waterline in this file's frame.
    pub waterline_z: f64,
}

/// A hull file's contents as source geometry; `path` names the file, and
/// its format is sniffed by extension and content.
pub fn open_source_bytes(
    path: &str,
    bytes: Vec<u8>,
    settings: &LoadSettings,
) -> Result<SourceFile, String> {
    let lower = path.to_ascii_lowercase();
    // STL: by extension or binary layout (binary STL is not UTF-8).
    if lower.ends_with(".stl") || looks_binary_stl(&bytes) || std::str::from_utf8(&bytes).is_err() {
        let scale = settings.units.ok_or_else(|| {
            format!(
                "{path}: STL files carry no units; pass --units mm|cm|m|in|ft \
                 (or a scale to metres)"
            )
        })?;
        let src = hullgeom::stl::mesh_fleet(&bytes, scale, settings.waterline_z)
            .map_err(|e| format!("{path}: STL import failed: {e}"))?;
        return Ok(SourceFile {
            path: path.into(),
            kind: Kind::Stl,
            source: Box::new(src),
            waterline_z: settings.waterline_z,
        });
    }
    let text = String::from_utf8(bytes).expect("checked utf8");
    let first = text
        .lines()
        .find(|l| !l.trim().is_empty())
        .unwrap_or("")
        .trim_end();
    if first.starts_with("michell-hull") || first.starts_with("michell-offsets") {
        return Err(format!(
            "{path}: .hull control nets, offsets tables and sample grids are no longer \
             read; import the CAD geometry (IGES, STL or Camber JSON) instead"
        ));
    }
    if text.trim_start().starts_with('{') || lower.ends_with(".json") {
        let v: serde_json::Value =
            serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))?;
        if crate::camber::is_document(&v) {
            let doc = crate::camber::Document::parse(&v).map_err(|e| format!("{path}: {e}"))?;
            let g = crate::camber::geometry(&doc, settings.waterline, settings.units)?;
            let mut sf = crate::native::open_value(&g)?;
            sf.path = path.into();
            sf.waterline_z = settings.waterline_z;
            return Ok(sf);
        }
        if v.get("kind").is_some() && v.get("hulls").is_some() {
            let mut sf = crate::native::open_value(&v)?;
            sf.path = path.into();
            sf.waterline_z = settings.waterline_z;
            return Ok(sf);
        }
        return Err(format!(
            "{path}: expected a Camber document or native geometry JSON"
        ));
    }
    let looks_iges = lower.ends_with(".igs")
        || lower.ends_with(".iges")
        || first.len() >= 73 && matches!(first.as_bytes()[72], b'S' | b'G');
    if looks_iges {
        let src = source_fleet(&text, settings.waterline_z)
            .map_err(|e| format!("{path}: IGES import failed: {e}"))?;
        return Ok(SourceFile {
            path: path.into(),
            kind: Kind::Iges,
            source: Box::new(src),
            waterline_z: settings.waterline_z,
        });
    }
    Err(format!(
        "cannot determine the format of {path}: expected an IGES, STL or Camber/native JSON file"
    ))
}

/// A file opened for display rather than for cutting: its hulls clustered
/// with the water notionally above everything, so the whole geometry shows
/// whatever the waterline (even one that misses the hull entirely). The
/// returned source's `waterline_z` is 0: tessellations come back in the
/// file's own frame.
pub fn open_geometry(
    path: &str,
    bytes: Vec<u8>,
    settings: &LoadSettings,
) -> Result<SourceFile, String> {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".stl") || looks_binary_stl(&bytes) || std::str::from_utf8(&bytes).is_err() {
        let scale = settings.units.ok_or_else(|| {
            format!("{path}: STL files carry no units; give their scale (mm, m, in, …)")
        })?;
        let tris = hullgeom::stl::parse_stl(&bytes, scale).map_err(|e| format!("{path}: {e}"))?;
        let top = tris
            .iter()
            .flatten()
            .fold(f64::NEG_INFINITY, |m, v| m.max(v[2]));
        let src = hullgeom::stl::mesh_fleet(&bytes, scale, top + 1.0)
            .map_err(|e| format!("{path}: STL import failed: {e}"))?;
        return Ok(SourceFile {
            path: path.into(),
            kind: Kind::Stl,
            source: Box::new(src),
            waterline_z: 0.0,
        });
    }
    let text = String::from_utf8(bytes).expect("checked utf8");
    if text.trim_start().starts_with('{') || lower.ends_with(".json") {
        let v: serde_json::Value =
            serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))?;
        if crate::camber::is_document(&v) {
            let doc = crate::camber::Document::parse(&v).map_err(|e| format!("{path}: {e}"))?;
            let g = crate::camber::geometry(&doc, settings.waterline, settings.units)?;
            let mut sf = crate::native::open_value(&g)?;
            sf.path = path.into();
            return Ok(sf);
        }
        if v.get("kind").is_some() && v.get("hulls").is_some() {
            let mut sf = crate::native::open_value(&v)?;
            sf.path = path.into();
            return Ok(sf);
        }
        return Err(format!(
            "{path}: expected a Camber document or native geometry JSON"
        ));
    }
    let file =
        hullgeom::iges::parse(&text).map_err(|e| format!("{path}: IGES import failed: {e}"))?;
    let top = file
        .surfaces
        .iter()
        .flat_map(|s| s.ctrl.iter())
        .fold(f64::NEG_INFINITY, |m, p| m.max(p[2]));
    if !top.is_finite() {
        return Err(format!("{path}: no surfaces to show"));
    }
    let src =
        source_fleet(&text, top + 1.0).map_err(|e| format!("{path}: IGES import failed: {e}"))?;
    Ok(SourceFile {
        path: path.into(),
        kind: Kind::Iges,
        source: Box::new(src),
        waterline_z: 0.0,
    })
}

/// Largest beam [m] of a sectional hull: twice its widest section.
pub fn max_beam(hull: &SectionalHull) -> f64 {
    2.0 * hull
        .curves()
        .flat_map(|(_, c)| c.iter().map(|p| p.0))
        .fold(0.0f64, f64::max)
}

/// Parse a units name (or raw scale) to metres-per-unit.
pub fn parse_units(s: &str) -> Result<f64, String> {
    match s.trim() {
        "mm" => Ok(0.001),
        "cm" => Ok(0.01),
        "m" => Ok(1.0),
        "in" => Ok(0.0254),
        "ft" => Ok(0.3048),
        other => other
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v > 0.0)
            .ok_or_else(|| format!("--units {other:?}: expected mm|cm|m|in|ft or a scale")),
    }
}

/// Binary-STL detection: exact size match on the triangle count.
pub fn looks_binary_stl(bytes: &[u8]) -> bool {
    bytes.len() >= 84 && {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        bytes.len() == 84 + 50 * n
    }
}

/// Diverging blue–gray–red colormap on t ∈ [−1, 1] (trough → crest, gray at
/// undisturbed water), interpolated in linear-light RGB between fixed
/// anchors (cool and warm poles around a neutral midpoint).
pub fn diverging(t: f64) -> [u8; 3] {
    const ANCHORS: &[(f64, [u8; 3])] = &[
        (-1.0, [0x0d, 0x36, 0x6b]),
        (-0.75, [0x1c, 0x5c, 0xab]),
        (-0.5, [0x2a, 0x78, 0xd6]),
        (-0.3, [0x55, 0x98, 0xe7]),
        (-0.12, [0x9e, 0xc5, 0xf4]),
        (0.0, [0xf0, 0xef, 0xec]),
        (0.12, [0xf5, 0xb8, 0xab]),
        (0.3, [0xee, 0x8f, 0x77]),
        (0.5, [0xe3, 0x49, 0x48]),
        (0.75, [0xb0, 0x2a, 0x2a]),
        (1.0, [0x6b, 0x14, 0x14]),
    ];
    let t = t.clamp(-1.0, 1.0);
    let t = if t.is_nan() { 0.0 } else { t };
    let mut i = 0;
    while i + 2 < ANCHORS.len() && t > ANCHORS[i + 1].0 {
        i += 1;
    }
    let (t0, c0) = ANCHORS[i];
    let (t1, c1) = ANCHORS[i + 1];
    let s = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    let mut out = [0u8; 3];
    for k in 0..3 {
        let a = srgb_to_linear(c0[k]);
        let b = srgb_to_linear(c1[k]);
        out[k] = linear_to_srgb(a + (b - a) * s);
    }
    out
}

fn srgb_to_linear(v: u8) -> f64 {
    let x = v as f64 / 255.0;
    if x <= 0.04045 {
        x / 12.92
    } else {
        ((x + 0.055) / 1.055).powf(2.4)
    }
}

fn linear_to_srgb(x: f64) -> u8 {
    let x = x.clamp(0.0, 1.0);
    let v = if x <= 0.003_130_8 {
        12.92 * x
    } else {
        1.055 * x.powf(1.0 / 2.4) - 0.055
    };
    (v * 255.0).round() as u8
}
