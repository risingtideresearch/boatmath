//! Pictures of a calm-water result, as SVG: the wake seen from above, the
//! pressure on the hull seen from below, and the hull in profile at its
//! attitude with the wave along its side. All three read the result's
//! `field` blob: the free surface `ζ(x, y)` and each hull's pressure
//! coefficient on its centreplane sheet (stations × depths, with the
//! half-breadth at each node). The profile's hull is its geometry re-posed
//! at the attitude of the result's sections.
//!
//! Frame: x forward (the hull advances toward +x), y to port, z up from the
//! still water, metres. Lengths are drawn to scale; where a picture would
//! be too thin to read, its short axis is stretched and the axis says by
//! how much.

use crate::plot::{esc, label, sig, ticks};
use boatmath::diverging;
use serde_json::Value;
use std::fmt::Write;

// ------------------------------------------------------------------ data

pub struct HullField {
    x: Vec<f64>,
    depth: Vec<f64>,
    /// `[ix * depth.len() + j]`.
    half_beam: Vec<f64>,
    cp: Vec<f64>,
    /// Centreplane [m].
    y: f64,
}

impl HullField {
    /// Half-breadth at the waterline at `x`, 0 off the hull.
    fn waterline(&self, x: f64) -> f64 {
        let nd = self.depth.len();
        let n = self.x.len();
        if n < 2 || x < self.x[0] || x > self.x[n - 1] {
            return 0.0;
        }
        let i = self.x.partition_point(|&v| v <= x).clamp(1, n - 1);
        let (a, b) = (self.x[i - 1], self.x[i]);
        let t = if b > a { (x - a) / (b - a) } else { 0.0 };
        (1.0 - t) * self.half_beam[(i - 1) * nd] + t * self.half_beam[i * nd]
    }

    fn x_range(&self) -> (f64, f64) {
        (self.x[0], self.x[self.x.len() - 1])
    }
}

pub struct Surface {
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    nx: usize,
    ny: usize,
    /// `[iy * nx + ix]`, row 0 at `y0`.
    zeta: Vec<f32>,
}

impl Surface {
    /// `ζ` at `(x, y)`, bilinear, or `None` off the grid.
    fn at(&self, x: f64, y: f64) -> Option<f64> {
        let fx = (x - self.x0) / (self.x1 - self.x0) * (self.nx - 1) as f64;
        let fy = (y - self.y0) / (self.y1 - self.y0) * (self.ny - 1) as f64;
        if !(0.0..=(self.nx - 1) as f64).contains(&fx)
            || !(0.0..=(self.ny - 1) as f64).contains(&fy)
        {
            return None;
        }
        let (i, j) = (
            (fx.floor() as usize).min(self.nx - 2),
            (fy.floor() as usize).min(self.ny - 2),
        );
        let (tx, ty) = (fx - i as f64, fy - j as f64);
        let z = |i: usize, j: usize| self.zeta[j * self.nx + i] as f64;
        Some(
            (1.0 - ty) * ((1.0 - tx) * z(i, j) + tx * z(i + 1, j))
                + ty * ((1.0 - tx) * z(i, j + 1) + tx * z(i + 1, j + 1)),
        )
    }
}

pub struct Field {
    hulls: Vec<HullField>,
    surface: Surface,
}

fn b64_decode(s: &str) -> Result<Vec<u8>, String> {
    let val = |c: u8| -> Result<u32, String> {
        Ok(match c {
            b'A'..=b'Z' => c - b'A',
            b'a'..=b'z' => c - b'a' + 26,
            b'0'..=b'9' => c - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            _ => return Err(format!("base64: unexpected {:?}", c as char)),
        } as u32)
    };
    let bytes: Vec<u8> = s
        .bytes()
        .filter(|&c| c != b'=' && !c.is_ascii_whitespace())
        .collect();
    let mut out = Vec::with_capacity(bytes.len() * 3 / 4);
    for chunk in bytes.chunks(4) {
        let mut n = 0u32;
        for (k, &c) in chunk.iter().enumerate() {
            n |= val(c)? << (18 - 6 * k);
        }
        let take = chunk.len().saturating_sub(1);
        for k in 0..take {
            out.push((n >> (16 - 8 * k)) as u8);
        }
    }
    Ok(out)
}

fn b64_encode(bytes: &[u8]) -> String {
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

pub(crate) fn f32s(v: &Value, what: &str) -> Result<Vec<f32>, String> {
    let bytes = b64_decode(v.as_str().ok_or_else(|| format!("field: no {what}"))?)?;
    Ok(bytes
        .as_chunks::<4>()
        .0
        .iter()
        .copied()
        .map(f32::from_le_bytes)
        .collect())
}

fn floats(v: &Value, what: &str) -> Result<Vec<f64>, String> {
    v.as_array()
        .ok_or_else(|| format!("field: no {what}"))?
        .iter()
        .map(|x| {
            x.as_f64()
                .ok_or_else(|| format!("field: {what}: not a number"))
        })
        .collect()
}

impl Field {
    /// A calm-water result's field, from its blob's JSON.
    pub fn parse(v: &Value) -> Result<Field, String> {
        let s = &v["surface"];
        let n = |k: &str| s[k].as_f64().ok_or_else(|| format!("field: surface {k}"));
        let surface = Surface {
            x0: n("x0")?,
            x1: n("x1")?,
            y0: n("y0")?,
            y1: n("y1")?,
            nx: n("nx")? as usize,
            ny: n("ny")? as usize,
            zeta: f32s(&s["zeta"], "zeta")?,
        };
        if surface.nx < 2 || surface.ny < 2 || surface.zeta.len() != surface.nx * surface.ny {
            return Err("field: a malformed surface grid".into());
        }
        let hulls = v["hulls"]
            .as_array()
            .ok_or("field: no hulls")?
            .iter()
            .map(|h| {
                let hf = HullField {
                    x: floats(&h["x"], "x")?,
                    depth: floats(&h["depth"], "depth")?,
                    half_beam: floats(&h["half_beam"], "half_beam")?,
                    cp: floats(&h["cp"], "cp")?,
                    y: h["y"].as_f64().unwrap_or(0.0),
                };
                let cells = hf.x.len() * hf.depth.len();
                if hf.x.len() < 2 || hf.half_beam.len() != cells || hf.cp.len() != cells {
                    return Err("field: a malformed pressure grid".to_string());
                }
                Ok(hf)
            })
            .collect::<Result<Vec<_>, String>>()?;
        if hulls.is_empty() {
            return Err("field: no hulls".into());
        }
        Ok(Field { hulls, surface })
    }
}

// ------------------------------------------------------------------ png

fn png_rgb(w: usize, h: usize, rgb: &[u8]) -> Vec<u8> {
    use flate2::write::ZlibEncoder;
    use std::io::Write as _;
    let mut raw = Vec::with_capacity(h * (3 * w + 1));
    for row in rgb.chunks(3 * w) {
        raw.push(0);
        raw.extend_from_slice(row);
    }
    let mut z = ZlibEncoder::new(Vec::new(), flate2::Compression::default());
    z.write_all(&raw).expect("in memory");
    let idat = z.finish().expect("in memory");
    let mut out = b"\x89PNG\r\n\x1a\n".to_vec();
    let mut chunk = |kind: &[u8], data: &[u8]| {
        out.extend_from_slice(&(data.len() as u32).to_be_bytes());
        let mut crc = flate2::Crc::new();
        crc.update(kind);
        crc.update(data);
        out.extend_from_slice(kind);
        out.extend_from_slice(data);
        out.extend_from_slice(&crc.sum().to_be_bytes());
    };
    let mut ihdr = Vec::new();
    ihdr.extend_from_slice(&(w as u32).to_be_bytes());
    ihdr.extend_from_slice(&(h as u32).to_be_bytes());
    ihdr.extend_from_slice(&[8, 2, 0, 0, 0]);
    chunk(b"IHDR", &ihdr);
    chunk(b"IDAT", &idat);
    chunk(b"IEND", &[]);
    out
}

// ------------------------------------------------------------------ frame

/// The 99th percentile of `|v|`: a colour scale that a single spike (at a
/// bow, say) doesn't wash out.
fn robust_max(v: impl Iterator<Item = f64>) -> f64 {
    let mut a: Vec<f64> = v.filter(|x| x.is_finite()).map(f64::abs).collect();
    if a.is_empty() {
        return 1.0;
    }
    a.sort_by(f64::total_cmp);
    let m = a[((a.len() - 1) as f64 * 0.99) as usize];
    if m > 0.0 {
        m
    } else {
        a[a.len() - 1].max(1e-9)
    }
}

/// A picture's data box and where it lands on the page.
struct Frame {
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    left: f64,
    top: f64,
    pw: f64,
    ph: f64,
    /// Draw +y downward (a view from below, with port at the bottom).
    flip: bool,
    /// How much the short axis is stretched beyond true scale.
    stretch: f64,
}

impl Frame {
    /// Lengths to scale at `pw` wide, unless that leaves the picture less
    /// than `min_h` tall: then y is stretched to it.
    fn new(x: (f64, f64), y: (f64, f64), top: f64, min_h: f64, flip: bool) -> Frame {
        let (left, pw) = (64.0, 840.0);
        let true_h = pw * (y.1 - y.0) / (x.1 - x.0);
        let ph = true_h.clamp(min_h, 900.0);
        Frame {
            x0: x.0,
            x1: x.1,
            y0: y.0,
            y1: y.1,
            left,
            top,
            pw,
            ph,
            flip,
            stretch: ph / true_h,
        }
    }

    fn sx(&self, x: f64) -> f64 {
        self.left + (x - self.x0) / (self.x1 - self.x0) * self.pw
    }

    fn sy(&self, y: f64) -> f64 {
        let t = (y - self.y0) / (self.y1 - self.y0);
        self.top + if self.flip { t } else { 1.0 - t } * self.ph
    }

    fn width(&self) -> f64 {
        self.left + self.pw + 24.0
    }

    fn bottom(&self) -> f64 {
        self.top + self.ph
    }
}

fn open(w: f64, h: f64, title: &str, subtitle: &str) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w:.0} {h:.0}" width="{w:.0}" height="{h:.0}" font-family="system-ui, -apple-system, sans-serif" font-size="12">
<style>:root{{--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--grid:#e4e3df;--hull:#c9c7c0;--part:#8f8d86;--water:#dbe8f7;--wave:#2a78d6}}
@media (prefers-color-scheme: dark){{:root{{--surface:#1a1a19;--ink:#ffffff;--ink2:#c3c2b7;--grid:#3a3a37;--hull:#6b6a64;--part:#9a9890;--water:#1f2f45;--wave:#3987e5}}}}
text{{fill:var(--ink2)}} .title{{fill:var(--ink);font-size:15px;font-weight:600}}</style>
<rect width="{w:.0}" height="{h:.0}" fill="var(--surface)"/>
<text class="title" x="64" y="24">{}</text><text x="64" y="42">{}</text>
"#,
        esc(title),
        esc(subtitle)
    );
    s
}

/// Ticks and labels on the frame's left and bottom edges.
fn axes(s: &mut String, f: &Frame, x_label: &str, y_label: &str) {
    let xt = ticks(f.x0, f.x1, 8);
    let xs = xt.get(1).zip(xt.first()).map_or(1.0, |(b, a)| b - a);
    for &x in &xt {
        let _ = writeln!(
            s,
            r#"<line x1="{x:.1}" x2="{x:.1}" y1="{b}" y2="{}" stroke="var(--ink2)"/><text x="{x:.1}" y="{}" text-anchor="middle">{}</text>"#,
            f.bottom() + 5.0,
            f.bottom() + 19.0,
            label(x, xs),
            x = f.sx(x),
            b = f.bottom()
        );
    }
    let yt = ticks(f.y0, f.y1, ((f.ph / 60.0).round() as usize).max(2));
    let ys = yt.get(1).zip(yt.first()).map_or(1.0, |(b, a)| b - a);
    for &y in &yt {
        let _ = writeln!(
            s,
            r#"<line x1="{}" x2="{}" y1="{y:.1}" y2="{y:.1}" stroke="var(--ink2)"/><text x="{}" y="{:.1}" text-anchor="end">{}</text>"#,
            f.left - 5.0,
            f.left,
            f.left - 8.0,
            f.sy(y) + 4.0,
            label(y, ys),
            y = f.sy(y)
        );
    }
    let _ = writeln!(
        s,
        r#"<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="var(--ink2)"/>"#,
        f.left, f.top, f.pw, f.ph
    );
    let _ = writeln!(
        s,
        r#"<text x="{}" y="{}" text-anchor="middle">{}</text>"#,
        f.left + f.pw / 2.0,
        f.bottom() + 36.0,
        esc(x_label)
    );
    let stretched = if f.stretch > 1.05 {
        format!(" (×{} stretched)", sig(f.stretch))
    } else {
        String::new()
    };
    let _ = writeln!(
        s,
        r#"<text transform="translate(18 {}) rotate(-90)" text-anchor="middle">{}{}</text>"#,
        f.top + f.ph / 2.0,
        esc(y_label),
        stretched
    );
}

/// A diverging colour bar for `±range`, under the frame at the right.
fn colorbar(s: &mut String, f: &Frame, range: f64, unit: &str, scale: f64, lo: &str, hi: &str) {
    let (w, x) = (240.0, f.left + f.pw - 240.0);
    let y = f.bottom() + 52.0;
    let _ = write!(s, r#"<defs><linearGradient id="cb">"#);
    for k in 0..=20 {
        let t = -1.0 + 0.1 * k as f64;
        let c = diverging(t);
        let _ = write!(
            s,
            r#"<stop offset="{:.3}" stop-color="rgb({},{},{})"/>"#,
            k as f64 / 20.0,
            c[0],
            c[1],
            c[2]
        );
    }
    let _ = writeln!(s, "</linearGradient></defs>");
    let _ = writeln!(
        s,
        r#"<rect x="{x}" y="{y}" width="{w}" height="10" rx="2" fill="url(#cb)"/>"#
    );
    let r = range * scale;
    for (t, txt) in [
        (0.0, format!("−{} {lo}", sig(r))),
        (0.5, "0".to_string()),
        (1.0, format!("+{} {hi}", sig(r))),
    ] {
        let anchor = ["start", "middle", "end"][(t * 2.0) as usize];
        let _ = writeln!(
            s,
            r#"<text x="{:.1}" y="{}" text-anchor="{anchor}">{}</text>"#,
            x + t * w,
            y + 24.0,
            esc(&txt)
        );
    }
    let _ = writeln!(
        s,
        r#"<text x="{:.1}" y="{}" text-anchor="end">{}</text>"#,
        x - 10.0,
        y + 9.0,
        esc(unit)
    );
}

/// The hull's waterline outline from above (or below), as SVG points.
fn waterline_points(h: &HullField, f: &Frame) -> String {
    let nd = h.depth.len();
    let side = |s: f64| {
        h.x.iter()
            .enumerate()
            .map(move |(i, &x)| (x, h.y + s * h.half_beam[i * nd]))
    };
    side(1.0)
        .chain(side(-1.0).collect::<Vec<_>>().into_iter().rev())
        .map(|(x, y)| format!("{:.1},{:.1}", f.sx(x), f.sy(y)))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Metres shown in millimetres when they're small.
fn length_unit(range: f64) -> (&'static str, f64) {
    if range < 0.5 {
        ("mm", 1e3)
    } else {
        ("m", 1.0)
    }
}

// ------------------------------------------------------------------ views

/// The wake from above: `ζ(x, y)` coloured troughs to crests, the hulls'
/// waterlines over it.
pub fn wake(field: &Field, title: &str, subtitle: &str, range: Option<f64>) -> String {
    let g = &field.surface;
    let a = range.unwrap_or_else(|| robust_max(g.zeta.iter().map(|&z| z as f64)));
    let mut rgb = Vec::with_capacity(3 * g.nx * g.ny);
    for row in 0..g.ny {
        let iy = g.ny - 1 - row;
        for ix in 0..g.nx {
            rgb.extend_from_slice(&diverging(g.zeta[iy * g.nx + ix] as f64 / a));
        }
    }
    let png = png_rgb(g.nx, g.ny, &rgb);
    let f = Frame::new((g.x0, g.x1), (g.y0, g.y1), 60.0, 120.0, false);
    let mut s = open(f.width(), f.bottom() + 96.0, title, subtitle);
    let _ = writeln!(
        s,
        r#"<image x="{}" y="{}" width="{}" height="{}" preserveAspectRatio="none" href="data:image/png;base64,{}"/>"#,
        f.left,
        f.top,
        f.pw,
        f.ph,
        b64_encode(&png)
    );
    for h in &field.hulls {
        let _ = writeln!(
            s,
            r#"<polygon points="{}" fill="var(--hull)" stroke="var(--ink)" stroke-width="1"/>"#,
            waterline_points(h, &f)
        );
    }
    axes(&mut s, &f, "x [m], forward →", "y [m], port ↑");
    let (unit, k) = length_unit(a);
    colorbar(&mut s, &f, a, &format!("ζ [{unit}]"), k, "trough", "crest");
    s.push_str("</svg>\n");
    s
}

/// The pressure on the hulls from below: `C_p` on each hull's surface,
/// projected onto the plan, port at the bottom (as seen from beneath).
pub fn pressure(field: &Field, title: &str, subtitle: &str, range: Option<f64>) -> String {
    let a =
        range.unwrap_or_else(|| robust_max(field.hulls.iter().flat_map(|h| h.cp.iter().copied())));
    let (mut xa, mut xb, mut ya, mut yb) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for h in &field.hulls {
        let (a, b) = h.x_range();
        let m = h.half_beam.iter().fold(0.0f64, |m, v| m.max(*v));
        xa = xa.min(a);
        xb = xb.max(b);
        ya = ya.min(h.y - m);
        yb = yb.max(h.y + m);
    }
    let (px, py) = (0.03 * (xb - xa), 0.15 * (yb - ya).max(1e-6));
    let f = Frame::new((xa - px, xb + px), (ya - py, yb + py), 60.0, 160.0, true);
    let mut s = open(f.width(), f.bottom() + 96.0, title, subtitle);
    for h in &field.hulls {
        let nd = h.depth.len();
        let at = |i: usize, j: usize| i * nd + j;
        // Shallow rows first: from below, the deeper hull is nearer.
        for j in 0..nd - 1 {
            for i in 0..h.x.len() - 1 {
                let corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)];
                if corners.iter().all(|&(i, j)| h.half_beam[at(i, j)] <= 0.0) {
                    continue;
                }
                let cp = corners.iter().map(|&(i, j)| h.cp[at(i, j)]).sum::<f64>() / 4.0;
                let c = diverging(cp / a);
                for side in [1.0, -1.0] {
                    let pts: Vec<String> = corners
                        .iter()
                        .map(|&(i, j)| {
                            format!(
                                "{:.1},{:.1}",
                                f.sx(h.x[i]),
                                f.sy(h.y + side * h.half_beam[at(i, j)])
                            )
                        })
                        .collect();
                    let _ = writeln!(
                        s,
                        r#"<polygon points="{}" fill="rgb({r},{g},{b})" stroke="rgb({r},{g},{b})" stroke-width="0.5"/>"#,
                        pts.join(" "),
                        r = c[0],
                        g = c[1],
                        b = c[2]
                    );
                }
            }
        }
        let _ = writeln!(
            s,
            r#"<polygon points="{}" fill="none" stroke="var(--ink)" stroke-width="1"/>"#,
            waterline_points(h, &f)
        );
    }
    axes(&mut s, &f, "x [m], forward →", "y [m], port ↓ (from below)");
    colorbar(&mut s, &f, a, "C_p", 1.0, "suction", "pressure");
    s.push_str("</svg>\n");
    s
}

/// The first hull in profile at its attitude: its silhouette (from `verts`,
/// its mesh there), the still water, and the wave along its side (`ζ` at
/// its waterline half-breadth, on its centreplane beyond its ends), scaled
/// by `wave_scale`.
pub fn profile(
    field: &Field,
    verts: &[[f64; 3]],
    parts: &[&[[f64; 3]]],
    title: &str,
    subtitle: &str,
    attitude: (f64, f64),
    wave_scale: f64,
) -> Result<String, String> {
    let h = &field.hulls[0];
    if verts.is_empty() {
        return Err("the hull has no mesh".into());
    }
    let bins = silhouette(verts, 240);
    let (xa, xb) = (bins[0].0, bins[bins.len() - 1].0);
    let l = xb - xa;
    // A drive's parts, each its own silhouette: finer, being short.
    let part_bins: Vec<Vec<(f64, f64, f64)>> = parts
        .iter()
        .filter(|p| !p.is_empty())
        .map(|p| silhouette(p, 60))
        .collect();
    let (x0, x1) = part_bins
        .iter()
        .flatten()
        .fold((xa - 0.25 * l, xb + 0.25 * l), |(a, b), p| {
            (a.min(p.0 - 0.05 * l), b.max(p.0 + 0.05 * l))
        });
    let wave: Vec<(f64, f64)> = (0..=400)
        .filter_map(|k| {
            let x = x0 + (x1 - x0) * k as f64 / 400.0;
            let y = h.y + h.waterline(x);
            field.surface.at(x, y).map(|z| (x, wave_scale * z))
        })
        .collect();
    let zmin = bins
        .iter()
        .chain(part_bins.iter().flatten())
        .map(|b| b.1)
        .chain(wave.iter().map(|w| w.1))
        .fold(0.0, f64::min);
    let zmax = bins
        .iter()
        .chain(part_bins.iter().flatten())
        .map(|b| b.2)
        .chain(wave.iter().map(|w| w.1))
        .fold(0.0, f64::max);
    let pad = 0.15 * (zmax - zmin).max(1e-3);
    let f = Frame::new((x0, x1), (zmin - pad, zmax + pad), 60.0, 160.0, false);
    let mut s = open(f.width(), f.bottom() + 60.0, title, subtitle);
    // The water below the still waterline.
    let _ = writeln!(
        s,
        r#"<rect x="{}" y="{:.1}" width="{}" height="{:.1}" fill="var(--water)"/>"#,
        f.left,
        f.sy(0.0),
        f.pw,
        f.bottom() - f.sy(0.0)
    );
    let outline = |bins: &[(f64, f64, f64)]| -> String {
        bins.iter()
            .map(|b| (b.0, b.2))
            .chain(bins.iter().rev().map(|b| (b.0, b.1)))
            .map(|(x, z)| format!("{:.1},{:.1}", f.sx(x), f.sy(z)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    // The parts first, so the hull covers the leg's root inside it.
    for p in &part_bins {
        let _ = writeln!(
            s,
            r#"<polygon points="{}" fill="var(--part)" stroke="var(--ink)" stroke-width="1" stroke-linejoin="round"/>"#,
            outline(p)
        );
    }
    let _ = writeln!(
        s,
        r#"<polygon points="{}" fill="var(--hull)" stroke="var(--ink)" stroke-width="1" stroke-linejoin="round"/>"#,
        outline(&bins)
    );
    let _ = writeln!(
        s,
        r#"<line x1="{}" x2="{}" y1="{y:.1}" y2="{y:.1}" stroke="var(--ink2)" stroke-dasharray="4 3"/>"#,
        f.left,
        f.left + f.pw,
        y = f.sy(0.0)
    );
    let line: Vec<String> = wave
        .iter()
        .map(|&(x, z)| format!("{:.1},{:.1}", f.sx(x), f.sy(z)))
        .collect();
    let _ = writeln!(
        s,
        r#"<polyline points="{}" fill="none" stroke="var(--wave)" stroke-width="2" stroke-linejoin="round"/>"#,
        line.join(" ")
    );
    let (sinkage, trim) = attitude;
    let mut note = format!(
        "sinkage {} mm · trim {}° {}",
        sig(sinkage * 1e3),
        sig(trim.to_degrees().abs()),
        if trim >= 0.0 { "bow up" } else { "bow down" }
    );
    if wave_scale != 1.0 {
        let _ = write!(note, " · wave ×{}", sig(wave_scale));
    }
    if field.hulls.len() > 1 + parts.len() {
        note.push_str(" · first hull of the platform");
    }
    let _ = writeln!(
        s,
        r#"<text x="{}" y="{}" text-anchor="end">{}</text>"#,
        f.left + f.pw - 6.0,
        f.top + 16.0,
        esc(&note)
    );
    axes(&mut s, &f, "x [m], forward →", "z [m], up");
    s.push_str("</svg>\n");
    Ok(s)
}

/// A mesh's side view: the highest and lowest point in each of `nb` slices
/// of x, `(x, z_low, z_high)` for each slice it reaches.
fn silhouette(verts: &[[f64; 3]], nb: usize) -> Vec<(f64, f64, f64)> {
    let (xa, xb) = verts
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
            (a.min(v[0]), b.max(v[0]))
        });
    let l = (xb - xa).max(1e-9);
    let mut lo = vec![f64::INFINITY; nb];
    let mut hi = vec![f64::NEG_INFINITY; nb];
    for v in verts {
        let k = (((v[0] - xa) / l * nb as f64) as usize).min(nb - 1);
        lo[k] = lo[k].min(v[2]);
        hi[k] = hi[k].max(v[2]);
    }
    (0..nb)
        .filter(|&k| lo[k].is_finite())
        .map(|k| (xa + (k as f64 + 0.5) / nb as f64 * l, lo[k], hi[k]))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_round_trips() {
        for n in 0..8 {
            let b: Vec<u8> = (0..n as u8).map(|i| i.wrapping_mul(37)).collect();
            assert_eq!(b64_decode(&b64_encode(&b)).unwrap(), b);
        }
    }
}
