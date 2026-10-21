//! A calm-water result as CAD surfaces: the hulls at their attitude (their
//! own patches, re-posed exactly), the free surface through every point of
//! the wave grid, and propellers (B-series blades, or a disc where none has
//! been chosen). Frame: x forward, y to port, z up,
//! the still water at z = 0, metres.

use crate::params::CaseParams;
use crate::{setup, LoftRequest};
use hullgeom::iges::{interpolate_grid, NurbsSurface3, SectionalImport};
use thinship::propulsion::Disc;

/// Interpolate a sectional loft into B-spline surface patches for a mesh hull (from an STL).
fn sectional_surfaces(imp: &SectionalImport) -> Result<Vec<NurbsSurface3>, String> {
    let ns = imp.polar.len();
    if ns < 4 {
        return Err("a sectional loft needs at least 4 stations".into());
    }
    let n_h = 17;
    let n_theta = 2 * n_h - 1;
    let thetas: Vec<f64> = (0..n_h)
        .map(|k| {
            0.5 * std::f64::consts::FRAC_PI_2
                * (1.0 - (std::f64::consts::PI * k as f64 / (n_h - 1) as f64).cos())
        })
        .collect();
    let y_c = imp.placement.y;
    let mut pts = Vec::with_capacity(ns * n_theta);
    for sec in &imp.polar {
        let x = sec.x;
        if sec.radii.is_empty() {
            for _ in 0..n_theta {
                pts.push([x, y_c, 0.0]);
            }
        } else {
            // Port side: waterline down to keel.
            for &th in &thetas {
                let (y, depth) = sec.point(th);
                pts.push([x, y_c + y, -depth]);
            }
            // Starboard side: keel up to waterline.
            for &th in thetas.iter().rev().skip(1) {
                let (y, depth) = sec.point(th);
                pts.push([x, y_c - y, -depth]);
            }
        }
    }
    let s = interpolate_grid(ns, n_theta, &pts).map_err(|e| e.to_string())?;
    Ok(vec![s])
}

/// Each of the case's members — its hulls and a drive's parts — as its
/// role (`hull`, `leg`, `pod`) and surface patches at `(sinkage [m], trim
/// [rad])`.
pub fn hull_surfaces(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    attitude: (f64, f64),
) -> Result<Vec<(String, Vec<NurbsSurface3>)>, String> {
    let s = setup(name, bytes, cut, c)?;
    let platform = s.platform(attitude.0, attitude.1);
    s.layout
        .iter()
        .zip(&s.roles)
        .map(|((i, pose), (role, _))| {
            let patches =
                match s
                    .cut
                    .file
                    .source
                    .posed_surfaces(*i, s.cut.file.waterline_z, pose, &platform)
                {
                    Some(p) => p.map_err(|e| e.to_string())?,
                    None => {
                        let imp = s
                            .cut
                            .file
                            .source
                            .situate_sectional(
                                *i,
                                s.cut.file.waterline_z,
                                pose,
                                &platform,
                                &s.cut.opts,
                            )
                            .map_err(|e| e.to_string())?
                            .ok_or("a mesh hull (from an STL) that does not enter the water")?;
                        sectional_surfaces(&imp)?
                    }
                };
            Ok((role.clone(), patches))
        })
        .collect()
}

/// The free surface `z = scale · ζ(x, y)` through every point of its grid.
#[allow(clippy::too_many_arguments)]
pub fn water_surface(
    x0: f64,
    x1: f64,
    nx: usize,
    y0: f64,
    y1: f64,
    ny: usize,
    zeta: &[f64],
    scale: f64,
) -> Result<NurbsSurface3, String> {
    let z: Vec<f64> = zeta.iter().map(|v| scale * v).collect();
    hullgeom::iges::graph_surface(x0, x1, nx, y0, y1, ny, &z).map_err(|e| e.to_string())
}

/// A propeller disc as a flat annulus, hub to tip, square to the shaft:
/// a rational quadratic circle swept radially, exact.
pub fn disc_surface(d: &Disc) -> NurbsSurface3 {
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let square = [
        (1.0, 0.0, 1.0),
        (1.0, 1.0, h),
        (0.0, 1.0, 1.0),
        (-1.0, 1.0, h),
        (-1.0, 0.0, 1.0),
        (-1.0, -1.0, h),
        (0.0, -1.0, 1.0),
        (1.0, -1.0, h),
        (1.0, 0.0, 1.0),
    ];
    let (yc, zc) = (d.y, -d.depth);
    let radii = [d.hub * d.radius, d.radius];
    let mut ctrl = Vec::new();
    let mut weights = Vec::new();
    for &(cy, cz, w) in &square {
        for &r in &radii {
            ctrl.push([d.x, yc + r * cy, zc + r * cz]);
            weights.push(w);
        }
    }
    NurbsSurface3 {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        n_ctrl_u: 9,
        n_ctrl_v: 2,
        ctrl,
        weights,
        trim_uv: None,
    }
}

// --------------------------------------------------------------- propellers

/// A propeller as `prop` chose it: a Wageningen B-series screw.
#[derive(Clone, Copy, Debug)]
pub struct Propeller {
    pub blades: usize,
    pub diameter: f64,
    pub pitch_ratio: f64,
    /// Expanded blade-area ratio `A_E/A_0`.
    pub area_ratio: f64,
    /// Turning clockwise seen from astern.
    pub right_handed: bool,
}

/// The B-series blade, by `r/R` (Kuiper 1992, *The Wageningen Propeller
/// Series*, table of dimensions).
/// Chord `c Z / (D A_E/A_0)`, for 4 to 7 blades …
const B_CHORD: [f64; 9] = [1.662, 1.882, 2.050, 2.152, 2.187, 2.144, 1.970, 1.582, 0.0];
/// … and for 3.
const B_CHORD_3: [f64; 9] = [1.633, 1.832, 2.000, 2.120, 2.186, 2.168, 2.127, 1.657, 0.0];
/// The leading edge ahead of the generator line, a fraction of the chord.
const B_LE: [f64; 9] = [
    0.616, 0.611, 0.599, 0.583, 0.558, 0.526, 0.481, 0.400, 0.400,
];
/// The greatest thickness aft of the leading edge, a fraction of the chord.
const B_TMAX_AT: [f64; 9] = [
    0.350, 0.350, 0.351, 0.355, 0.389, 0.443, 0.479, 0.500, 0.500,
];
/// Greatest thickness `t/D = A − B Z`.
const B_T_A: [f64; 9] = [
    0.0526, 0.0464, 0.0402, 0.0340, 0.0278, 0.0216, 0.0154, 0.0092, 0.0030,
];
const B_T_B: [f64; 9] = [
    0.0040, 0.0035, 0.0030, 0.0025, 0.0020, 0.0015, 0.0010, 0.0005, 0.0,
];
/// The four-bladed series' pitch, reduced toward the root.
const B4_PITCH: [f64; 9] = [0.822, 0.887, 0.950, 0.992, 1.0, 1.0, 1.0, 1.0, 1.0];
/// The generator line's rake aft.
const B_RAKE_DEG: f64 = 15.0;

/// A table (at `r/R` = 0.2, 0.3, …, 1.0) at `r` in `[0.2, 0.9]`:
/// Catmull–Rom through its points.
fn b_table(t: &[f64; 9], r: f64) -> f64 {
    let u = ((r - 0.2) / 0.1).clamp(0.0, 7.0);
    let i = (u.floor() as usize).min(6);
    let f = u - i as f64;
    let p = |k: isize| t[(i as isize + k).clamp(0, 7) as usize];
    let (p0, p1, p2, p3) = (p(-1), p(0), p(1), p(2));
    0.5 * (2.0 * p1
        + (p2 - p0) * f
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * f * f
        + (3.0 * p1 - p0 - 3.0 * p2 + p3) * f * f * f)
}

/// A propeller at a disc, its shaft at `angle` (bow up, to the still
/// water): each blade a surface through its expanded sections wrapped on
/// their pitch helices and raked aft, and the hub a body of revolution.
/// The outline, thickness, rake and pitch are the series'; the sections
/// follow the series' tabulated ordinates (Kuiper 1992 / Oosterveld & van Oossanen
/// 1975) with aerofoil sections at the root and segmental sections outboard.
pub fn propeller_surfaces(
    d: &Disc,
    angle: f64,
    p: &Propeller,
) -> Result<Vec<NurbsSurface3>, String> {
    use std::f64::consts::PI;
    let (dia, z) = (p.diameter, p.blades);
    if z < 2 || !(dia > 0.0 && p.pitch_ratio > 0.0 && p.area_ratio > 0.0) {
        return Err(format!(
            "a propeller of {z} blades, D {dia}: not a propeller"
        ));
    }
    let big_r = 0.5 * dia;
    let hand = if p.right_handed { 1.0 } else { -1.0 };
    // Square to the shaft: along its axis (forward), across (port) and up.
    let (ca, sa) = (angle.cos(), angle.sin());
    let centre = [d.x, d.y, -d.depth];
    let place = |ax: f64, y: f64, up: f64| -> [f64; 3] {
        [
            centre[0] + ax * ca - up * sa,
            centre[1] + y,
            centre[2] + ax * sa + up * ca,
        ]
    };
    let chord_t = if z == 3 { &B_CHORD_3 } else { &B_CHORD };
    // The hub covers the blades' roots, and whatever it sits on (a pod).
    let r_hub = (0.105 * dia).max(1.02 * d.hub * d.radius);
    let r0 = (0.97 * r_hub / big_r).max(0.2);
    if r0 >= 0.8 {
        return Err(format!(
            "the propeller ({dia:.3} m) is too small for its hub ({:.3} m)",
            2.0 * r_hub
        ));
    }

    // Radii, closer toward the tip, where the outline turns.
    let (n_r, n_c) = (21, 17);
    let radii: Vec<f64> = (0..n_r)
        .map(|i| r0 + (1.0 - r0) * (0.5 * PI * i as f64 / (n_r - 1) as f64).sin())
        .collect();
    let fracs: Vec<f64> = (0..n_c)
        .map(|k| 0.5 * (1.0 - (PI * k as f64 / (n_c - 1) as f64).cos()))
        .collect();
    let mut out = Vec::new();
    let mut root_x = (f64::INFINITY, f64::NEG_INFINITY);
    // The first blade's tip points straight up, toward the hull, so its
    // clearance reads off the drawing.
    for blade in 0..z {
        let theta0 = 0.5 * PI + 2.0 * PI * blade as f64 / z as f64;
        let mut pts = Vec::with_capacity(n_r * (2 * n_c - 1));
        for (ir, &r) in radii.iter().enumerate() {
            // The outline closes elliptically over the last tenth.
            let tip = if r > 0.9 {
                ((1.0 - r) / 0.1).max(0.0).sqrt()
            } else {
                1.0
            };
            let rr = r.min(0.9);
            let c = dia * p.area_ratio / z as f64 * b_table(chord_t, rr) * tip;
            let le = b_table(&B_LE, rr);
            let at = b_table(&B_TMAX_AT, rr);
            let t = dia * (b_table(&B_T_A, rr) - b_table(&B_T_B, rr) * z as f64).max(0.0) * tip;
            let pitch = p.pitch_ratio * dia * if z == 4 { b_table(&B4_PITCH, rr) } else { 1.0 };
            let ra = r * big_r;
            let phi = (pitch / (2.0 * PI * ra)).atan();
            let rake = -(ra - 0.2 * big_r) * B_RAKE_DEG.to_radians().tan();
            // Around the section: the back from the trailing edge to the
            // leading edge, then the face back to the trailing edge.
            let back = fracs.iter().rev().map(|&f| (f, true));
            let face = fracs.iter().skip(1).map(|&f| (f, false));
            for (f, on_back) in back.chain(face) {
                // Wageningen B-series section ordinates (Kuiper 1992 / Oosterveld & van Oossanen 1975):
                // Outer sections (r/R ≥ 0.5) are segmental: flat face, parabolic back (V2 = 1 - P²).
                // Inner sections (r/R < 0.5) have aerofoil wash-back on the face near the leading edge (V1).
                let p_coord = if f < at {
                    (at - f) / at
                } else {
                    (f - at) / (1.0 - at)
                };
                let v2 = (1.0 - p_coord * p_coord).max(0.0);
                let v1 = if f < at && rr < 0.5 {
                    0.15 * ((0.5 - rr) / 0.3) * p_coord * p_coord
                } else {
                    0.0
                };
                let eta = if on_back { (v1 + v2) * t } else { v1 * t };
                // Along the chord toward the leading edge from the
                // generator line, on the helix; the back faces forward.
                let xi = (le - f) * c;
                let arc = xi * phi.cos() - eta * phi.sin();
                let ax = xi * phi.sin() + eta * phi.cos() + rake;
                if ir == 0 {
                    root_x = (root_x.0.min(ax), root_x.1.max(ax));
                }
                let th = theta0 + hand * arc / ra;
                pts.push(place(ax, ra * th.cos(), ra * th.sin()));
            }
        }
        out.push(
            hullgeom::iges::interpolate_grid(n_r, 2 * n_c - 1, &pts).map_err(|e| e.to_string())?,
        );
    }

    // The hub, over the roots and a little more: rounded forward into the
    // shaft, a cone aft.
    let (x_fwd, x_aft) = (root_x.1 + 0.04 * dia, root_x.0 - 0.06 * dia);
    let (n_u, n_v) = (15, 17);
    let mut pts = Vec::with_capacity(n_u * n_v);
    for iu in 0..n_u {
        let s = 0.5 * (1.0 - (PI * iu as f64 / (n_u - 1) as f64).cos());
        let ax = x_fwd + (x_aft - x_fwd) * s;
        let r = r_hub
            * if s < 0.15 {
                (1.0 - ((0.15 - s) / 0.15).powi(2)).max(0.0).sqrt()
            } else if s > 0.6 {
                let u = (s - 0.6) / 0.4;
                (1.0 - u * u).max(0.0).sqrt() * (1.0 - 0.5 * u)
            } else {
                1.0
            };
        for iv in 0..n_v {
            let th = 2.0 * PI * iv as f64 / (n_v - 1) as f64;
            pts.push(place(ax, r * th.cos(), r * th.sin()));
        }
    }
    out.push(hullgeom::iges::interpolate_grid(n_u, n_v, &pts).map_err(|e| e.to_string())?);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The blades' expanded area is the series' `A_E/A_0` of the disc, by
    /// the outline the surfaces are built from; each blade spans hub to tip
    /// and the hub sits on the shaft line.
    #[test]
    fn a_b_series_propeller_has_its_area_and_reach() {
        let d = Disc {
            x: -4.0,
            y: 0.5,
            depth: 0.6,
            radius: 0.2,
            hub: 0.2,
        };
        for z in [2, 3, 4, 5] {
            let p = Propeller {
                blades: z,
                diameter: 0.4,
                pitch_ratio: 1.0,
                area_ratio: 0.55,
                right_handed: true,
            };
            let s = propeller_surfaces(&d, 0.1, &p).unwrap();
            assert_eq!(s.len(), z + 1);
            // The outline's area, hub to tip, against A_E/A_0 · πR²
            // (the series' A_E counts from 0.2R; the hub's part is small).
            let n = 4000;
            let chord_t = if z == 3 { &B_CHORD_3 } else { &B_CHORD };
            let area: f64 = (0..n)
                .map(|i| {
                    let r = 0.2 + 0.8 * (i as f64 + 0.5) / n as f64;
                    let tip = if r > 0.9 {
                        ((1.0 - r) / 0.1).sqrt()
                    } else {
                        1.0
                    };
                    0.4 * 0.55 / z as f64 * b_table(chord_t, r.min(0.9)) * tip * 0.8 * 0.2
                        / n as f64
                })
                .sum::<f64>()
                * z as f64;
            let disc = std::f64::consts::PI * 0.2 * 0.2;
            assert!(
                (area / disc - 0.55).abs() < 0.04,
                "Z {z}: A_E/A_0 {}",
                area / disc
            );
            // Every control point within the disc's radius of the shaft line
            // (and a little: the interpolant overshoots).
            let axis = [0.1f64.cos(), 0.0, 0.1f64.sin()];
            for q in s.iter().flat_map(|s| s.ctrl.iter()) {
                let v = [q[0] - d.x, q[1] - d.y, q[2] + d.depth];
                let along = v[0] * axis[0] + v[2] * axis[2];
                let off = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - along * along)
                    .max(0.0)
                    .sqrt();
                assert!(off < 1.05 * 0.2 && along.abs() < 0.2, "{q:?}");
            }
            // The first blade's tip is straight up, square to the shaft.
            let top = s[0]
                .ctrl
                .iter()
                .max_by(|a, b| a[2].total_cmp(&b[2]))
                .unwrap();
            let v = [top[0] - d.x, top[1] - d.y, top[2] + d.depth];
            let up = [-0.1f64.sin(), 0.0, 0.1f64.cos()];
            let reach = v[0] * up[0] + v[2] * up[2];
            assert!(
                v[1].abs() < 0.01 && (reach - 0.2).abs() < 0.01,
                "Z {z}: tip {v:?}"
            );
        }
    }

    /// An STL hull has its sectional loft interpolated into B-spline surface
    /// patches for CAD export.
    #[test]
    fn an_stl_hull_exports_cad_surfaces() {
        use hullgeom::iges::{HullPose, Platform, SourceFleet};
        let surfaces = hullgeom::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let fleet = SourceFleet::from_hulls(vec![surfaces.to_vec()]).unwrap();
        let (v, t) = fleet
            .posed_tessellation(0, 0.0, &HullPose::default(), &Platform::default())
            .unwrap();
        let mut text = String::from("solid w\n");
        let wl = 0.4;
        for tri in t {
            text.push_str("facet normal 0 0 0\nouter loop\n");
            for k in tri {
                let p = v[k as usize];
                text.push_str(&format!(
                    "vertex {} {} {}\n",
                    p[0] * 1e3,
                    p[1] * 1e3,
                    (p[2] + wl) * 1e3
                ));
            }
            text.push_str("endloop\nendfacet\n");
        }
        text.push_str("endsolid w\n");
        let bytes = text.into_bytes();
        let req = LoftRequest {
            waterline: Some(wl),
            units: Some(1e-3),
            ..LoftRequest::default()
        };
        let cp = CaseParams::default();
        let members = hull_surfaces("w.stl", bytes, &req, &cp, (0.0, 0.0)).unwrap();
        assert_eq!(members.len(), 1);
        assert_eq!(members[0].0, "hull");
        assert!(!members[0].1.is_empty());
        let surf = &members[0].1[0];
        assert!(surf.n_ctrl_u >= 4);
        assert!(surf.n_ctrl_v >= 4);
        for p in &surf.ctrl {
            assert!(p[0] >= -5.5 && p[0] <= 5.5, "x out of range: {}", p[0]);
            assert!(p[1] >= -0.7 && p[1] <= 0.7, "y out of range: {}", p[1]);
            assert!(p[2] >= -0.85 && p[2] <= 0.15, "z out of range: {}", p[2]);
        }
    }
}
