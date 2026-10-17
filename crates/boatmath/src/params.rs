//! What the store keeps below a hull, in canonical form, so that asking for
//! the same thing twice finds the first:
//!
//! - a **case** ([`CaseParams`]): the platform the hull makes —
//!   on its own or doubled into a catamaran — and its load (mass, centre of
//!   gravity, radii of gyration, roll damping). Its statics (the float at
//!   rest, GM, the GZ curve) are computed when it is made.
//! - a **study** ([`StudyParams`]) on a case: a speed and the model's
//!   settings, in calm water or in waves from one heading. A study in waves
//!   is taken about the attitude of the calm-water study at its speed
//!   ([`StudyParams::calm`]), which it waits for.
//!
//! The hull's import settings (waterline, stations, units, ...) belong to
//! the hull. `canonical` fills in every default and checks every range, so
//! two requests that mean the same computation serialize — and hash — alike.

use crate::MassBy;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use thinship::{TransomClosure, BALLISTIC_COEFF, DRY_HOLLOW_FRACTION};

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn sha(s: &str) -> String {
    hex(&Sha256::digest(s.as_bytes()))
}

fn finite(name: &str, v: Option<f64>) -> Result<(), String> {
    match v {
        Some(x) if !x.is_finite() => Err(format!("{name}: not a number")),
        _ => Ok(()),
    }
}

fn positive(name: &str, v: Option<f64>) -> Result<(), String> {
    match v {
        Some(x) if !(x > 0.0 && x.is_finite()) => Err(format!("{name} must be positive")),
        _ => Ok(()),
    }
}

// ---------------------------------------------------------------- case

/// A platform on a hull and its load.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CaseParams {
    /// Double the (single) hull into a catamaran with this centre span, the
    /// distance between the demihulls' centreplanes [m].
    #[serde(default)]
    pub span: Option<f64>,
    /// Load [kg]; default the design displacement.
    #[serde(default)]
    pub mass: Option<f64>,
    /// How a given mass is carried: by sinking, or by scaling the hull.
    #[serde(default)]
    pub mass_by: MassBy,
    /// Longitudinal centre of gravity [m, fleet x]; default the LCB.
    #[serde(default)]
    pub lcg: Option<f64>,
    /// Height of the centre of gravity above the design waterline [m];
    /// default 0 (at the waterline).
    #[serde(default)]
    pub vcg: Option<f64>,
    /// Radii of gyration: roll `k_xx` [m], pitch `k_yy` as a fraction of the
    /// length, yaw `k_zz` [m]; defaults as `seakeeping::platform`.
    #[serde(default)]
    pub kxx: Option<f64>,
    #[serde(default)]
    pub kyy: Option<f64>,
    #[serde(default)]
    pub kzz: Option<f64>,
    /// Roll damping, a fraction of critical (bilge keels, appendages).
    #[serde(default)]
    pub roll_damping: f64,
}

impl Default for CaseParams {
    fn default() -> Self {
        CaseParams {
            span: None,
            mass: None,
            mass_by: MassBy::Sinking,
            lcg: None,
            vcg: None,
            kxx: None,
            kyy: None,
            kzz: None,
            roll_damping: 0.0,
        }
    }
}

impl CaseParams {
    pub fn canonical(mut self) -> Result<CaseParams, String> {
        positive("span", self.span)?;
        positive("mass", self.mass)?;
        finite("lcg", self.lcg)?;
        finite("vcg", self.vcg)?;
        positive("kxx", self.kxx)?;
        positive("kyy", self.kyy)?;
        positive("kzz", self.kzz)?;
        if !(self.roll_damping >= 0.0 && self.roll_damping < 1.0) {
            return Err("roll damping: expected a fraction of critical, 0 to 1".into());
        }
        if self.mass.is_none() && self.mass_by != MassBy::Sinking {
            // Scaling to the design displacement is no scaling at all.
            self.mass_by = MassBy::Sinking;
        }
        Ok(self)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("params serialize")
    }

    /// Hex SHA-256 of the canonical JSON: the case's identity on
    /// its hull.
    pub fn hash(&self) -> String {
        sha(&self.to_json())
    }
}

// ---------------------------------------------------------------- study

/// The transom closure, as stored.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Closure {
    /// The free-fall hollow; `coeff` defaults to the solver's.
    Ballistic {
        #[serde(default)]
        coeff: Option<f64>,
    },
    /// The free-fall hollow while the transom is wet, `fraction` of a
    /// transverse wavelength once it runs dry; both default to the solver's.
    Auto {
        #[serde(default)]
        coeff: Option<f64>,
        #[serde(default)]
        fraction: Option<f64>,
    },
    /// A hollow of fixed length [m].
    Fixed {
        length: f64,
    },
    Off,
}

impl Default for Closure {
    fn default() -> Self {
        Closure::Auto {
            coeff: None,
            fraction: None,
        }
    }
}

impl Closure {
    pub fn transom(self) -> TransomClosure {
        match self {
            Closure::Ballistic { coeff } => TransomClosure::Ballistic {
                coeff: coeff.unwrap_or(BALLISTIC_COEFF),
            },
            Closure::Auto { coeff, fraction } => TransomClosure::Auto {
                coeff: coeff.unwrap_or(BALLISTIC_COEFF),
                fraction: fraction.unwrap_or(DRY_HOLLOW_FRACTION),
            },
            Closure::Fixed { length } => TransomClosure::Fixed { length },
            Closure::Off => TransomClosure::None,
        }
    }

    fn canonical(self) -> Result<Closure, String> {
        let or = |v: Option<f64>, default: f64, what: &str| match v {
            Some(c) if c.is_finite() => Ok(c.max(0.0)),
            Some(c) => Err(format!("closure {what} {c}: not a number")),
            None => Ok(default),
        };
        Ok(match self {
            Closure::Ballistic { coeff } => Closure::Ballistic {
                coeff: Some(or(coeff, BALLISTIC_COEFF, "coeff")?),
            },
            Closure::Auto { coeff, fraction } => Closure::Auto {
                coeff: Some(or(coeff, BALLISTIC_COEFF, "coeff")?),
                fraction: Some(or(fraction, DRY_HOLLOW_FRACTION, "fraction")?),
            },
            Closure::Fixed { length } if length.is_finite() => Closure::Fixed {
                length: length.max(0.0),
            },
            Closure::Fixed { length } => {
                return Err(format!("closure length {length}: not a number"))
            }
            Closure::Off => Closure::Off,
        })
    }
}

/// An irregular sea, as stored.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Sea {
    /// Bretschneider / ITTC, significant height `hs` [m], peak period `tp` [s].
    Bretschneider { hs: f64, tp: f64 },
    /// JONSWAP with peak enhancement `gamma` (default 3.3).
    Jonswap {
        hs: f64,
        tp: f64,
        #[serde(default)]
        gamma: Option<f64>,
    },
}

impl Sea {
    pub fn spectrum(self) -> seakeeping::sea::Spectrum {
        use seakeeping::sea::Spectrum;
        match self {
            Sea::Bretschneider { hs, tp } => Spectrum::Bretschneider { hs, tp },
            Sea::Jonswap { hs, tp, gamma } => Spectrum::Jonswap {
                hs,
                tp,
                gamma: gamma.unwrap_or(3.3),
            },
        }
    }

    fn canonical(self) -> Result<Sea, String> {
        let (hs, tp) = match self {
            Sea::Bretschneider { hs, tp } | Sea::Jonswap { hs, tp, .. } => (hs, tp),
        };
        positive("sea hs", Some(hs))?;
        positive("sea tp", Some(tp))?;
        Ok(match self {
            Sea::Jonswap { hs, tp, gamma } => {
                let g = gamma.unwrap_or(3.3);
                if !(1.0..=10.0).contains(&g) {
                    return Err(format!("sea gamma {g}: expected 1 to 10"));
                }
                Sea::Jonswap {
                    hs,
                    tp,
                    gamma: Some(g),
                }
            }
            s => s,
        })
    }
}

/// Regular waves from one heading, over a range of wavelengths, and
/// optionally one irregular sea from the same heading.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Waves {
    /// [deg]: 180 head seas, 90 beam, 0 following.
    pub heading: f64,
    /// Wavelengths over the longest hull's length; default 0.5 to 3 by 0.125.
    #[serde(default)]
    pub lambdas: Option<Vec<f64>>,
    #[serde(default)]
    pub sea: Option<Sea>,
}

/// The wavelengths a sweep takes by default: λ/L from 0.5 to 3 by 0.125.
pub fn default_lambdas() -> Vec<f64> {
    (0..=20).map(|i| 0.5 + 0.125 * i as f64).collect()
}

/// A speed on a case, in calm water or in waves.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StudyParams {
    /// Length Froude number on the longest hull.
    pub froude: f64,
    /// Float at the dynamic equilibrium at speed; else hold the attitude at
    /// rest.
    #[serde(default = "default_dynamic")]
    pub dynamic: bool,
    #[serde(default)]
    pub closure: Closure,
    /// Free-surface grid columns (a calm-water study's field; display only).
    #[serde(default = "default_grid")]
    pub grid: usize,
    /// How far astern the field reaches, in lengths behind the sterns
    /// (`None`: 1.5).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub astern: Option<f64>,
    /// In waves; `None` in calm water.
    #[serde(default)]
    pub waves: Option<Waves>,
    /// The hulls' form factor `k`, as `(1 + k)·C_F` (`None`: 0, the flat
    /// plate). A drive's parts keep their own.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub form_factor: Option<f64>,
    /// Roughness allowance `ΔC_F`, on every wetted part (`None`: smooth).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub roughness: Option<f64>,
    /// Charge a transom running dry the hydrostatic thrust it loses.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub transom_drag: bool,
    /// The hulls' dynamic lift beside thin-ship theory's thickness, in the
    /// attitude and in the waves (`None`: none).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lift: Option<LiftModel>,
}

/// How a study's hulls lift (see [`crate::lift`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LiftModel {
    /// Slender-body 2D+t on the sections at the attitude.
    #[serde(rename = "2dt")]
    TwoDT,
}

pub fn default_grid() -> usize {
    640
}

fn default_dynamic() -> bool {
    true
}

impl StudyParams {
    /// `"calm"` or `"waves"`.
    pub fn kind(&self) -> &'static str {
        if self.waves.is_some() {
            "waves"
        } else {
            "calm"
        }
    }

    /// Checked, with every default filled in (so equal computations compare
    /// equal). A study in waves has no field of its own, so no grid.
    pub fn canonical(mut self) -> Result<StudyParams, String> {
        if !(self.froude > 0.0 && self.froude < 5.0) {
            return Err(format!("froude {}: expected 0 < Fn < 5", self.froude));
        }
        self.closure = self.closure.canonical()?;
        self.grid = self.grid.clamp(40, 1200);
        // A zero is the default: dropped, so the study is the one without it.
        for (name, v, hi) in [
            ("form factor", &mut self.form_factor, 2.0),
            ("roughness", &mut self.roughness, 0.01),
        ] {
            if let Some(x) = *v {
                if !(0.0..=hi).contains(&x) {
                    return Err(format!("{name} {x}: expected 0 to {hi}"));
                }
                *v = Some(x).filter(|&x| x > 0.0);
            }
        }
        if let Some(a) = self.astern {
            if !(a > 0.0 && a <= 20.0) {
                return Err(format!("astern {a}: expected 0 < lengths ≤ 20"));
            }
        }
        if let Some(w) = &mut self.waves {
            self.grid = default_grid();
            self.astern = None;
            if !w.heading.is_finite() {
                return Err("heading: not a number".into());
            }
            w.heading = w.heading.rem_euclid(360.0);
            let mut l = w.lambdas.take().unwrap_or_else(default_lambdas);
            if l.is_empty() || l.len() > 200 || l.iter().any(|&x| !(x > 0.0 && x.is_finite())) {
                return Err("wavelengths: expected 1 to 200 positive λ/L".into());
            }
            l.sort_by(f64::total_cmp);
            l.dedup();
            w.lambdas = Some(l);
            w.sea = w.sea.map(Sea::canonical).transpose()?;
        }
        Ok(self)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(self).expect("params serialize")
    }

    /// Hex SHA-256 of the canonical JSON: the study's identity on its
    /// case.
    pub fn hash(&self) -> String {
        sha(&self.to_json())
    }

    /// The calm-water study a study in waves is taken about: the same speed,
    /// attitude and closure (and the grid asked for with it).
    pub fn calm(&self, grid: usize) -> StudyParams {
        StudyParams {
            waves: None,
            grid,
            astern: None,
            ..self.clone()
        }
    }

    /// What a study's attitude depends on — the calm-water studies that share it
    /// (any grid) and the studies in waves taken about it: its speed, attitude
    /// and closure.
    pub fn attitude_key(&self) -> String {
        let key = (self.froude, self.dynamic, self.closure);
        sha(&match self.lift {
            // Without lift, the key it always had.
            None => serde_json::to_string(&key),
            Some(lift) => serde_json::to_string(&(key, lift)),
        }
        .expect("serialize"))
    }

    /// The same study at another speed: what a warm start may borrow from.
    pub fn same_but_speed(&self, other: &StudyParams) -> bool {
        self.dynamic == other.dynamic && self.closure == other.closure && self.lift == other.lift
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn study(s: &str) -> StudyParams {
        serde_json::from_str::<StudyParams>(s)
            .unwrap()
            .canonical()
            .unwrap()
    }

    fn case(s: &str) -> CaseParams {
        serde_json::from_str::<CaseParams>(s)
            .unwrap()
            .canonical()
            .unwrap()
    }

    /// Defaults spelt out and left out are the same thing.
    #[test]
    fn equal_requests_hash_alike() {
        let a = study(r#"{"froude": 0.3}"#);
        let b =
            study(r#"{"froude": 0.3, "grid": 640, "dynamic": true, "closure": {"type": "auto"}}"#);
        assert_eq!(a.hash(), b.hash());
        assert_ne!(a.hash(), study(r#"{"froude": 0.31}"#).hash());
        let ballistic = study(r#"{"froude": 0.3, "closure": {"type": "ballistic"}}"#);
        assert_ne!(
            a.hash(),
            ballistic.hash(),
            "a stored ballistic study stays ballistic"
        );
        assert_eq!(a.hash(), study(r#"{"froude": 0.3, "astern": null}"#).hash());
        assert_ne!(a.hash(), study(r#"{"froude": 0.3, "astern": 4}"#).hash());
        let w = study(r#"{"froude": 0.3, "grid": 200, "waves": {"heading": -180}}"#);
        assert_eq!(w.kind(), "waves");
        assert_eq!(w.grid, 640, "a study in waves has no grid of its own");
        assert_eq!(w.waves.as_ref().unwrap().heading, 180.0);
        assert_eq!(
            w.waves.as_ref().unwrap().lambdas.as_ref().unwrap().len(),
            21
        );
        assert_eq!(w.attitude_key(), a.attitude_key());
        assert_eq!(w.calm(640).hash(), a.hash());

        let c = case(r#"{"mass_by": "scale"}"#);
        assert_eq!(c.hash(), case("{}").hash(), "no mass, nothing to scale");
        assert_ne!(case(r#"{"span": 3}"#).hash(), case("{}").hash());
    }

    #[test]
    fn bad_requests_are_refused() {
        for s in [
            r#"{"froude": 0}"#,
            r#"{"froude": 0.3, "waves": {"heading": 180, "lambdas": []}}"#,
            r#"{"froude": 0.3, "waves": {"heading": 180, "sea": {"type": "jonswap", "hs": 0, "tp": 5}}}"#,
        ] {
            let p: StudyParams = serde_json::from_str(s).unwrap();
            assert!(p.canonical().is_err(), "{s}");
        }
        for s in [
            r#"{"mass": -1}"#,
            r#"{"span": 0}"#,
            r#"{"roll_damping": 2}"#,
        ] {
            let p: CaseParams = serde_json::from_str(s).unwrap();
            assert!(p.canonical().is_err(), "{s}");
        }
        assert!(serde_json::from_str::<StudyParams>(r#"{"froude": 0.3, "speed": 2}"#).is_err());
    }
}
