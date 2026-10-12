# Graf's catamaran, from US patent 5,188,049

Lofts a demihull of Lawrence Graf's displacement catamaran (US 5,188,049,
1993) from the patent's drawings, as IGES or as a camber document. The
patent's own example: 20 ft 8 in waterline, 24 in waterline beam, a 16 in
by 36 in tunnel, 2,700 lb.

`graf-sections.json` is the drawing, put in numbers: the ten half-sections
of Fig. 5, traced from a 300 dpi rendering of the patent PDF, at the
stations Fig. 2 marks, scaled to the 24 in beam (and checked against the
patent's 1.6° keel rise). `graf-loft.py` makes the hull from them; its
docstring has what is the patent's and what is a choice.

```sh
pip install numpy scipy
python graf-loft.py graf-sections.json > graf-demihull.igs
python graf-loft.py graf-sections.json --camber > graf-camber.json      # 7 stations x 5 points
python graf-loft.py graf-sections.json --camber-u 0.0,0.4492,0.8142,0.8708,1.0 \
    --camber-points 2 > graf-camber-5.json                              # 5 x 6
```

The camber documents' stations are fitted to the IGES hull, not cut from
it: monotone from deck to keel, the chine ledge level and dry. Against the
IGES, both displace the same to 0.04% and float on the same LCB to 0.1 in;
their surfaces are within 0.54 in (7 stations) and 0.61 in (5), worst at
the stem. The five-station one is a fixture of `crates/boatmath/tests/camber.rs`
(`graf.json`).

The catamaran at the patent's example, floating free at 6-18 kn:

```sh
boatmath hull graf-demihull.igs --waterline -0.001016 \
  | boatmath case --span 1.705 --mass 1225 \
  | boatmath study --froude 0.39268,0.52358,0.65447,0.78537,0.91626,1.04716,1.17805 \
  | boatmath run > graf-cat-calm.jsonl
```

The waterline is 1 mm under the stern chine, so the chine ledge is dry. The
sections as drawn displace more than the patent's 2,700 lb at the chine;
at 2,700 lb the boat floats 1.9 in higher, where its wetted surface
(106 ft²) is close to the patent's 110.
