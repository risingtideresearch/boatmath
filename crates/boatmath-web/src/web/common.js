// What every page shares: small helpers, the API client, and the top bar.

export const $ = (id) => document.getElementById(id);
export const G = 9.81;
export const KN = 0.514444;
// The solver's fluid: salt water at 15 °C (hullgeom::Fluid::SEAWATER_15C).
export const RHO = 1025.9;

export const fmt = (v, d = 3) => v == null || !Number.isFinite(Number(v)) ? "—" : Number(v).toFixed(d);
export const fmtBytes = (n) => n < 1024 ? `${n} B` : n < 1 << 20 ? `${(n / 1024).toFixed(0)} kB` : `${(n / (1 << 20)).toFixed(1)} MB`;
export const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

// A unix time as a short local date and time; today's as the time alone.
export function fmtTime(t) {
  if (t == null) return "—";
  const d = new Date(t * 1000), now = new Date();
  const time = d.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
  return d.toDateString() === now.toDateString() ? time
    : `${d.toLocaleDateString([], { day: "numeric", month: "short", ...(d.getFullYear() === now.getFullYear() ? {} : { year: "numeric" }) })} ${time}`;
}
export function fmtDuration(s) {
  if (s == null) return "—";
  return s < 60 ? `${s.toFixed(s < 10 ? 1 : 0)} s` : s < 3600 ? `${Math.floor(s / 60)} min ${Math.round(s % 60)} s` : `${(s / 3600).toFixed(1)} h`;
}

// Base64 to an ArrayBuffer (the server's little-endian typed arrays).
export function b64(s) {
  const bin = atob(s), out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out.buffer;
}

// The name this browser's uploads and requests are labelled with.
export function whoami() {
  try { return localStorage.getItem("boatmath.by") || ""; } catch { return ""; }
}
function setWhoami(v) {
  try { localStorage.setItem("boatmath.by", v); } catch { /* private window */ }
}

// A JSON API call: the parsed answer, or an Error with the server's message.
// With a `json` or raw `body` it is a POST unless `method` says otherwise.
export async function api(path, { method, body, json } = {}) {
  const opts = {};
  if (json !== undefined) { opts.body = JSON.stringify(json); opts.headers = { "Content-Type": "application/json" }; }
  else if (body !== undefined) opts.body = body;
  opts.method = method || (opts.body !== undefined ? "POST" : "GET");
  const r = await fetch(path, opts);
  const text = await r.text();
  let data;
  try { data = text ? JSON.parse(text) : null; } catch { data = { error: text || r.statusText }; }
  if (!r.ok) throw new Error((data && data.error) || `${r.status} ${r.statusText}`);
  return data;
}

// The id in a path like /hulls/12 (the last number in it).
export function pathId() {
  const m = location.pathname.match(/\/(\d+)\/?$/);
  return m ? Number(m[1]) : null;
}

// A hull's summary in words: principal dimensions and displacement.
export function hullDims(summary) {
  const hs = summary?.hulls || [];
  if (!hs.length) return "";
  const h = hs[0];
  const vol = hs.reduce((s, h) => s + (h.displaced_volume || 0), 0);
  const n = hs.length > 1 ? `${hs.length} hulls · ` : "";
  return `${n}L ${fmt(h.length, 2)} · B ${fmt(h.beam, 2)} · T ${fmt(h.draft, 3)} m · Δ ${fmt(RHO * vol, 0)} kg`;
}

// The reference length speeds are made Froude numbers on: the longest hull.
export const lRef = (summary) => Math.max(0, ...(summary?.hulls || []).map((h) => h.length || 0));
// A Froude number on length `l` as knots.
export const knots = (fn, l) => fn * Math.sqrt(G * l) / KN;

// A case's parameters in words.
export function caseLabel(p) {
  const how = { sinking: "", scale: " (hull scaled)", scale_yz: " (beam and draft scaled)" }[p.mass_by] ?? "";
  const parts = [
    p.span != null ? `catamaran, span ${fmt(p.span, 2)} m` : "monohull",
    p.mass == null ? "design mass" : `${fmt(p.mass, 1)} kg${how}`,
  ];
  if (p.lcg != null) parts.push(`LCG ${fmt(p.lcg, 3)}`);
  if (p.vcg != null) parts.push(`VCG ${fmt(p.vcg, 3)}`);
  if (p.kxx != null) parts.push(`kxx ${fmt(p.kxx, 2)} m`);
  if (p.kyy != null) parts.push(`kyy ${fmt(p.kyy, 2)} L`);
  if (p.kzz != null) parts.push(`kzz ${fmt(p.kzz, 2)} m`);
  if (p.roll_damping) parts.push(`roll damping ${fmt(100 * p.roll_damping, 0)}%`);
  return parts.join(" · ");
}

// A case's name, or its parameters when it has none.
export const caseName = (c) => c.name || caseLabel(c.params);

// Headings in words: head, bow, beam, quartering, following seas.
export function headingText(deg) {
  const d = ((deg % 360) + 360) % 360, a = d > 180 ? 360 - d : d;
  const what = a >= 165 ? "head" : a > 105 ? "bow" : a >= 75 ? "beam" : a > 15 ? "quartering" : "following";
  return `${fmt(d, 0)}° ${what} seas`;
}

// A study's parameters in words; `l` (the reference length) adds the speed
// in knots.
export function studyLabel(p, l) {
  const speed = `Fn ${fmt(p.froude, 3)}${l ? ` (${fmt(knots(p.froude, l), 2)} kn)` : ""}`;
  const c = p.closure || {};
  const closure = c.type === "off" ? "no closure"
    : c.type === "fixed" ? `hollow ${fmt(c.length, 2)} m`
    : c.type === "ballistic" ? "ballistic hollow"
    : c.type === "auto" && c.fraction != null && c.fraction !== 0.5 ? `dry hollow ${fmt(c.fraction, 2)} λ`
    : "";
  const w = p.waves;
  const sea = !w ? "calm water" : [headingText(w.heading), w.sea ? `${w.sea.type === "jonswap" ? "JONSWAP" : "Bretschneider"} Hs ${fmt(w.sea.hs, 2)} m Tp ${fmt(w.sea.tp, 1)} s` : ""].filter(Boolean).join(", ");
  return [speed, sea, p.dynamic === false ? "held at rest's attitude" : "", closure].filter(Boolean).join(" · ");
}

const TABS = [
  ["hulls", "/hulls", "Hulls"],
  ["cases", "/cases", "Cases"],
  ["studies", "/studies", "Studies"],
  ["queue", "/queue", "Queue"],
  ["plot", "/plot", "Plot"],
];

// A field's values: blank (the fallback), `a`, `a, b, …`, or `a:b:step`.
export function values(input, fallback) {
  const t = input.value.trim();
  if (!t) return [fallback];
  const out = [];
  for (const part of t.split(/[,\s]+/).filter(Boolean)) {
    const r = part.split(":").map(Number);
    if (r.length === 3 && r.every(Number.isFinite) && r[2] > 0 && r[1] >= r[0]) {
      for (let k = 0; r[0] + k * r[2] <= r[1] + 1e-9 * Math.abs(r[2]); k++) out.push(+(r[0] + k * r[2]).toPrecision(12));
    } else if (r.length === 1 && Number.isFinite(r[0])) out.push(r[0]);
    else throw new Error(`${input.closest("label").firstChild.textContent.trim()}: cannot read “${part}”`);
  }
  return out;
}

// A study's outcome in brief: R_t in calm water, the heave peak in waves.
export function outcome(r) {
  const s = r.result?.scalars;
  if (!s) return "";
  return r.kind === "calm" ? `R<sub>t</sub> ${fmt(s.rt, 1)} N` : `heave ${fmt(s.heave_peak, 2)} m/m`;
}

// A status pill; a done study with a stale result says so.
export const pill = (r) => `<span class="pill ${r.status}${r.stale ? " stale" : ""}">${r.stale && r.status === "done" ? "stale" : r.status}</span>`;

// The top bar: the pages, the queue's length, and who is asking.
export function nav(on) {
  const bar = document.createElement("nav");
  bar.className = "top";
  bar.innerHTML = `<a class="brand" href="/hulls">boatmath</a>`
    + TABS.map(([k, href, label]) => `<a class="tab${k === on ? " on" : ""}" href="${href}" data-tab="${k}">${label}${k === "queue" ? `<span class="badge" id="qbadge" hidden></span>` : ""}</a>`).join("")
    + `<label class="who">you <input id="whoami" placeholder="your name" autocomplete="name"></label>`;
  document.body.prepend(bar);
  const who = bar.querySelector("#whoami");
  who.value = whoami();
  who.addEventListener("change", () => setWhoami(who.value.trim()));
  // Behind tailscale serve the server knows who is asking, and labels work
  // with that: show it, in place of a name to type.
  api("/api/whoami").then((w) => {
    if (!w) return;
    const name = w.name || w.login;
    setWhoami(name);
    const label = bar.querySelector(".who");
    label.textContent = "";
    label.title = w.login;
    label.append("you ", Object.assign(document.createElement("b"), { textContent: name }));
  }).catch(() => {});
  // The queue's length, kept up to date while the page is open.
  const badge = bar.querySelector("#qbadge");
  async function poll() {
    try {
      const q = await api("/api/queue");
      const n = q.studies.length;
      badge.hidden = !n;
      badge.textContent = n;
    } catch { /* offline: leave it */ }
  }
  poll();
  setInterval(() => { if (!document.hidden) poll(); }, 5000);
}

// A status line: a message, or an error in red.
export function status(msg, err = false) {
  const s = $("status");
  if (!s) return;
  s.textContent = msg;
  s.className = err ? "err" : "";
}

// Each hull's hydrostatics, as the cut reports them.
export function hullCards(hulls, notes = []) {
  return notes.map((n) => `<div class="note">${esc(n)}</div>`).join("") + hulls.map((h, i) => {
    const title = hulls.length > 1 ? `Hull ${i + 1}` : "Hull";
    const where = h.placement.y || h.placement.x ? `y ${fmt(h.placement.y)} m` : "";
    const t = h.transom;
    return `<div class="hull">
      <h3>${title}<span class="tag">${esc(where)}</span></h3>
      <dl>
        <dt>Length</dt><dd>${fmt(h.length)} m</dd>
        <dt>Beam (WL)</dt><dd>${fmt(h.beam)} m</dd>
        <dt>Draft</dt><dd>${fmt(h.draft)} m</dd>
        <dt>Displaced volume</dt><dd>${fmt(h.displaced_volume, 4)} m³</dd>
        <dt>Displacement</dt><dd>${fmt(RHO * h.displaced_volume, 1)} kg</dd>
        <dt>Wetted surface</dt><dd>${fmt(h.wetted_surface)} m²</dd>
        <dt>Waterplane area</dt><dd>${fmt(h.waterplane_area)} m²</dd>
        <dt>LCB</dt><dd>x ${fmt(h.lcb_x)} m</dd>
        <dt>Stations</dt><dd>${h.stations.length}</dd>
        ${t ? `<dt>Transom</dt><dd>x ${fmt(t.x)} m, ${(100 * t.area_ratio).toFixed(1)}% of max section</dd>` : ""}
      </dl>
      <details><summary>Sectioning diagnostics</summary><pre>${esc(h.diagnostics.join("\n"))}</pre></details>
    </div>`;
  }).join("");
}
