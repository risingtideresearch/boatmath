// The golden data `camber.rs` is held to: camber's own sweep, run on each
// fixture document. Run from a camber checkout (its tsx resolves the imports):
//
//   cd ../camber && npx tsx <this dir>/make.mts <this dir>
//
// It writes, beside each `<name>.json` hull document, `<name>.golden.json`:
// the trimmed starboard half-section (`sweptSection`) at a spread of the
// plan's parameter u, the hull's aft and forward limits, and camber's
// hydrostatics at its own design waterline. `default.json`, `cruiser.json`
// and `flat-bottom.json` are written first: camber's default hull and two of
// its v1 examples, read through its own v1 → v2 conversion. `graf.json` is
// the Graf catamaran demihull (python/tools/graf), in inches: few stations,
// a curved sheer, a knuckled chine ledge and a bow made by the centreline
// trim.

import { readFileSync, writeFileSync } from "node:fs";
import { join, resolve } from "node:path";

const camber = process.cwd();
const dir = resolve(process.argv[2] ?? ".");
const src = (p: string) => join(camber, "src", p);

const { buildJson, parseDocument } = await import(src("core/json.ts"));
const { defaultHull } = await import(src("core/hull.ts"));
const { convertV1ToV2 } = await import(src("legacy/v1/convert.ts"));
const { assemble } = await import(src("core/runtime.ts"));
const { sweptSection, forwardLimit, aftLimit, computeHullSampling } =
  await import(src("core/mesh.ts"));
const { hydrostatics } = await import(src("core/hydro.ts"));

writeFileSync(join(dir, "default.json"), buildJson(defaultHull()) + "\n");
for (const [name, example] of [
  ["cruiser", "round-bilge-cruiser.json"],
  ["flat-bottom", "flat-bottom.json"],
]) {
  const v1 = JSON.parse(readFileSync(join(camber, "examples", example), "utf8"));
  writeFileSync(
    join(dir, `${name}.json`),
    JSON.stringify(convertV1ToV2(v1), null, 2) + "\n",
  );
}

const R = 4;
for (const name of ["default", "cruiser", "flat-bottom", "inverted-bow", "graf"]) {
  const text = readFileSync(join(dir, `${name}.json`), "utf8");
  const model = assemble(parseDocument(text).state);
  const aft = aftLimit(model),
    fwd = forwardLimit(model);
  const us = [...Array.from({ length: 41 }, (_, i) => i / 40), aft, fwd];
  const sections = us.map((u) => {
    const s = sweptSection(model, u, R, true);
    return { u, empty: s.empty, keel: s.keel, pts: s.pts };
  });
  const h = hydrostatics(model, computeHullSampling(model, 1200, 40));
  const golden = {
    R,
    aft_limit: aft,
    forward_limit: fwd,
    sections,
    hydro: h && {
      vol: h.vol,
      lwl: h.lwl,
      bwl: h.bwl,
      draft: h.draft,
      lcb: h.lcb,
      waterplane_area: h.waterplaneArea,
      closed: h.closed,
    },
  };
  writeFileSync(join(dir, `${name}.golden.json`), JSON.stringify(golden) + "\n");
  console.log(name, "aft", aft, "fwd", fwd, "vol", h?.vol, "closed", h?.closed);
}
