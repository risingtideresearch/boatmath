"""Graf catamaran demihull (US 5,188,049, 1993) -> IGES, from the patent's drawings.

    python3 graf-loft.py graf-sections.json > graf-demihull.igs
    python3 graf-loft.py graf-sections.json --camber > graf-camber.json

The second form writes the hull as a camber document instead (below, at
camber_document); --camber-u and --camber-points pick its stations and the
points on each, --camber-dense writes the hull's own cuts at 38 stations,
and --camber-search prints the trade-off of station count against error.

graf-sections.json holds the ten half-sections of Fig. 5, traced from the
300 dpi drawing and put in inches. The drawing has no scale, so it is scaled
so that its widest waterline half-breadth (station 7) is 12 in: the patent's
example has a 24 in waterline beam. Its height is scaled the same, which the
keel confirms. The patent puts the keel's rise at 1.6 deg from 37-39% of the
waterline length. Read off Fig. 5 at the stations of Fig. 2, the keel rises
1.3 / 1.8 / 2.9 / 4.0 in at stations 7-10, against 1.3 / 1.9 / 3.1 / 4.3 in
from those numbers. The stations sit where Fig. 2 puts them, along the
example's 20 ft 8 in waterline.

Each station is three pieces: the round bottom from the keel to the chine,
the chine ledge (flat, as in the patent's preferred embodiment) and a
straight flared topside, carried up to the bridge deck at the example's
16 in tunnel height. Forward, all three close onto a stem raked 25 deg from
the vertical (the patent's beta); the chine meets it 8.9 in above the
waterline, where Fig. 5 shows it reaching the centreline. Aft, a flat transom
closes the hull at the last station.

The keel is the patent's: the deepest point of Fig. 5 (13.4 in) back to 38%
of the waterline length, then rising at 1.6 deg. Stations 1 and 2 are on the
forefoot and keep their own keel heights.

Surfaces are bicubic interpolants through the sections (u along the hull,
v around the section), one patch per piece and side, mirrored to port, in
metres, x forward from the transom, z up from the design waterline. The
waterline is the stern chine, which the patent says runs along it.
"""
import json
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator, make_interp_spline

IN = 0.0254
LWL = 248.0           # 20 ft 8 in
DECK = 16.0           # tunnel height: waterline to the underside of the bridge deck
STEM_RAKE = np.radians(25.0)
STEM_CHINE_Z = 8.9    # where the chine meets the stem, Fig. 5
NV = 21               # points around each bottom section

secs = json.load(open(sys.argv[1]))
order = [str(k) for k in range(1, 11)]


def stem_s(z):
    """The stem: inches aft of the waterline's forward end at height z."""
    return -z * np.tan(STEM_RAKE)


def resample(y, z, n):
    """A section (chine -> keel) as n points keel -> chine, even in arc length."""
    y, z = np.asarray(y)[::-1], np.asarray(z)[::-1]
    t = np.r_[0, np.cumsum(np.hypot(np.diff(y), np.diff(z)))]
    ts = np.linspace(0, t[-1], n)
    return np.interp(ts, t, y), np.interp(ts, t, z)


# Rows of section curves, stem first, each as points (s aft, y, z) in inches.
bottom, ledge, topside = [], [], []
k1 = secs['1']
foot = k1['keel']
zs = np.linspace(foot, STEM_CHINE_Z, NV)
bottom.append([(stem_s(z), 0.0, z) for z in zs])
ledge.append([(stem_s(STEM_CHINE_Z), 0.0, STEM_CHINE_Z)] * 2)
topside.append([(stem_s(z), 0.0, z) for z in (STEM_CHINE_Z, DECK)])
for k in order:
    v = secs[k]
    y, z = resample(v['y'], v['z'], NV)
    s = v['s']
    bottom.append([(s, a, b) for a, b in zip(y, z)])
    zc, yi, yo = v['chine_z'], v['chine_in'], v['chine_out']
    ledge.append([(s, yi, zc), (s, yo, zc)])
    flare = np.tan(np.radians(v['flare_deg']))
    topside.append([(s, yo, zc), (s, yo + (DECK - zc) * flare, DECK)])


def densify(rows, step=6.0):
    """More sections between the stations, monotone (PCHIP) in s, so the cubic
    loft through them can't overshoot: through the stations alone it would
    dip 0.3 in below the keel behind the forefoot."""
    P = np.array(rows[1:], dtype=float)              # stations 1..10
    s = P[:, 0, 0]
    fine = np.unique(np.r_[s, np.arange(s[0], s[-1], step)])
    f = PchipInterpolator(s, P[:, :, 1:], axis=0)(fine)
    dense = np.concatenate([np.broadcast_to(fine[:, None, None], f.shape[:2] + (1,)), f], axis=2)
    return [rows[0]] + [list(map(tuple, r)) for r in dense]


bottom, ledge, topside = densify(bottom), densify(ledge), densify(topside)

# One u parameter for all three pieces, so their shared edges agree: chord
# length along the chine, which every piece touches.
chine = np.array([row[-1] for row in bottom])
u = np.r_[0, np.cumsum(np.linalg.norm(np.diff(chine, axis=0), axis=1))]
u /= u[-1]


def surface(rows, degree_v):
    """Interpolate rows (sections) of points: cubic along u, degree_v around."""
    P = np.array(rows, dtype=float)                  # (nu, nv, 3)
    nu, nv, _ = P.shape
    if degree_v == 1:
        kv = np.r_[0.0, np.linspace(0, 1, nv), 1.0]
        Cv = P
    else:
        v = np.linspace(0, 1, nv)
        sp = make_interp_spline(v, P, k=3, axis=1)
        kv, Cv = sp.t, np.moveaxis(sp.c, 0, 1)      # scipy puts v first
    sp = make_interp_spline(u, Cv, k=3, axis=0)
    return dict(ku=sp.t, kv=kv, du=3, dv=degree_v, ctrl=sp.c)


def to_metres(srf, mirror=False):
    c = srf['ctrl'].copy()
    x = (LWL - c[..., 0]) * IN                       # forward from the transom
    y = c[..., 1] * IN * (-1 if mirror else 1)
    z = c[..., 2] * IN
    return dict(srf, ctrl=np.stack([x, y, z], axis=-1))


# Transom: flat, from the centreline to the last section's outline.
last = secs['10']
yb, zb = resample(last['y'], last['z'], NV)
outline = list(zip(yb, zb))
outline.append((last['chine_out'], last['chine_z']))
fl = np.tan(np.radians(last['flare_deg']))
outline.append((last['chine_out'] + (DECK - last['chine_z']) * fl, DECK))
outline = np.array(outline)
centre = np.c_[np.zeros(len(outline)), outline[:, 1]]
transom = dict(du=1, dv=1, ku=np.array([0, 0, 1, 1.0]),
               kv=np.r_[0.0, np.linspace(0, 1, len(outline)), 1.0],
               ctrl=np.array([[(LWL, a, b) for a, b in centre],
                              [(LWL, a, b) for a, b in outline]]))

pieces = dict(bottom=surface(bottom, 3), ledge=surface(ledge, 1), topside=surface(topside, 1))
patches = []
for srf in (*pieces.values(), transom):
    patches += [to_metres(srf), to_metres(srf, mirror=True)]


# ---------------------------------------------------------------- camber out
def vline(srf, v, n=4000):
    """The curve S(u, v) along the hull, as (s, y, z) in inches, stem first."""
    from scipy.interpolate import BSpline
    ku, kv, du, dv = srf['ku'], srf['kv'], srf['du'], srf['dv']
    Bu = BSpline.design_matrix(np.linspace(0, 1, n), ku, du).toarray()
    Bv = BSpline.design_matrix(np.array([v]), kv, dv).toarray()[0]
    return Bu @ np.einsum('j,ijk->ik', Bv, srf['ctrl'])


def entry_angle(line, over=6.0):
    """A v-line's mean half-angle over its first `over` inches from the stem."""
    ls = line[np.argsort(line[:, 0])]
    near = ls[ls[:, 0] < ls[0, 0] + over]
    return np.arctan2(near[-1, 1] - near[0, 1], near[-1, 0] - near[0, 0])


def uniform_knots(n, p):
    """camber's sheer-plan knots: clamped, interior knots evenly spaced."""
    inner = n - p - 1
    return np.r_[np.zeros(p + 1), np.arange(1, inner + 1) / (inner + 1), np.ones(p + 1)]


def fit_plan(xy, n_ctrl=8, p=3, rounds=30):
    """A clamped uniform-knot B-spline through the ends of the polyline `xy`,
    least squares in between, its parameters refined by projecting each
    point onto the curve. Returns (control points, max deviation)."""
    from scipy.interpolate import BSpline
    knots = uniform_knots(n_ctrl, p)
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    t = d / d[-1]
    dense = np.linspace(0, 1, 20001)
    for _ in range(rounds):
        B = BSpline.design_matrix(t, knots, p).toarray()
        # Ends pinned, and the second control point level with the first, so
        # the plan leaves the transom square to it: solve for the rest.
        rhs = xy - np.outer(B[:, 0], xy[0]) - np.outer(B[:, -1], xy[-1])
        cx, *_ = np.linalg.lstsq(B[:, 1:-1], rhs[:, 0], rcond=None)
        ry = rhs[:, 1] - B[:, 1] * xy[0, 1]
        cy, *_ = np.linalg.lstsq(B[:, 2:-1], ry, rcond=None)
        inner = np.c_[cx, np.r_[xy[0, 1], cy]]
        ctrl = np.vstack([xy[0], inner, xy[-1]])
        curve = BSpline(knots, ctrl, p)(dense)
        nearest = np.argmin(np.linalg.norm(xy[:, None, :] - curve[None, ::20, :], axis=2), axis=1) * 20
        for i, j in enumerate(nearest):           # refine each to the dense grid
            lo, hi = max(0, j - 20), min(len(dense), j + 21)
            k = lo + np.argmin(np.linalg.norm(curve[lo:hi] - xy[i], axis=1))
            t[i] = dense[k]
        t[0], t[-1] = 0.0, 1.0
    dev = np.min(np.linalg.norm(xy[:, None, :] - curve[None, :, :], axis=2), axis=1)
    return ctrl, knots, dev.max()


def camber_setup(n_bottom=9):
    """What a camber (v2) document of the hull is made from. Its sheer plan is the deck edge,
    from the transom to the stem head, fitted with camber's clamped
    uniform-knot B-spline, leaving the transom square to it. Each station is the lofted hull cut in the
    vertical plane through the plan at its u, square to the plan (camber's
    section plane), as (n, z): n inboard from the plan along its normal, z
    from the deck datum (the bridge deck, 16 in up). The chine's two edges
    are knuckles. Forward of the stem, where a section plane passes ahead of
    a v-line's end on the stem, the v-line is carried on straight, to port
    at its own entry angle and rising as it rose into the stem, so camber's loft crosses the centreline at the
    stem and its centreline trim cuts the rest away."""
    from scipy.interpolate import BSpline
    to_x = lambda s: LWL - s                                        # noqa: E731
    # The bottom's points between the inner chine (v = 1) and the keel (v = 0).
    if n_bottom == 9:
        v_bottom = [1.0, 0.9, 0.78, 0.66, 0.54, 0.42, 0.31, 0.21, 0.12, 0.05, 0.0]
    else:
        v_bottom = list(np.linspace(1.0, 0.0, n_bottom + 2))
    lines = ([vline(pieces['topside'], 1.0), vline(pieces['topside'], 0.0)]
             + [vline(pieces['bottom'], v) for v in v_bottom])
    knuckled = {1, 2}                       # outer and inner chine
    entries = [entry_angle(line) for line in lines]
    entries[-1] = entries[-2]               # the keel is on the centreline: its neighbour's
    # Each v-line in camber's frame (x forward, y, z), stem first, carried 3 in
    # on behind the transom along its own direction there.
    xl = []
    for line in lines:
        c = np.c_[to_x(line[:, 0]), line[:, 1], line[:, 2]]   # stem first
        tail = c[-1] - c[-30]
        tail = tail / max(-tail[0], 1e-9)
        xl.append(np.vstack([c, c[-1] + np.outer(np.linspace(0.1, 3.0, 30), tail)]))

    # The plan: the deck edge (topside v = 1), transom to stem head.
    deck = xl[0][:, :2]
    deck = deck[deck[:, 0] >= 0.0]
    deck = deck[::-1]                       # aft to forward, as camber's u runs
    deck = np.vstack([deck[:-1:8], deck[-1]])   # thinned, for the fit
    deck[-1, 1] = 0.0                       # the stem head, on the centreline
    ctrl, knots, dev = fit_plan(deck)
    plan = BSpline(knots, ctrl, 3)
    dplan = plan.derivative()
    print(f"# sheer plan: {len(ctrl)} control points, max deviation from the deck "
          f"edge {dev:.4f} in", file=sys.stderr)

    us = np.linspace(0, 1, 20001)
    px = plan(us)[:, 0]

    def u_at(s):
        return float(np.interp(to_x(s), px, us))

    def rise(c, over=4.0):
        """dz/dx of a v-line over its last `over` inches up to the stem."""
        k = np.argmax(c[:, 0] < c[0, 0] - over)
        return (c[0, 2] - c[k, 2]) / (c[0, 0] - c[k, 0])

    rises = [rise(c) for c in xl]

    def cut(c, P, t, entry, gz):
        """Where v-line c crosses the plane through P normal to t: (n, z, on
        the hull)."""
        nh = np.array([t[1], -t[0]])
        f = (c[:, 0] - P[0]) * t[0] + (c[:, 1] - P[1]) * t[1]
        k = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
        if len(k):
            i = k[0]
            a = f[i] / (f[i] - f[i + 1])
            q, on = c[i] + a * (c[i + 1] - c[i]), True
        else:
            # Ahead of the line's stem end: carry it on, straight, to port
            # at its entry angle and rising as it rose into the stem.
            xs, zs = c[0, 0], c[0, 2]
            d = ((P[0] - xs) * t[0] + P[1] * t[1]) / (t[0] - t[1] * np.tan(entry))
            q, on = np.array([xs + d, -d * np.tan(entry), zs + d * gz]), False
        return float((q[0] - P[0]) * nh[0] + (q[1] - P[1]) * nh[1]), float(q[2]), on

    def section(u):
        """camber's station at plan parameter u: (n, z) per point, z from the
        waterline (inches)."""
        P = plan(u)
        d = dplan(u)
        t = d / np.hypot(*d)
        cuts = [cut(c, P, t, entries[i], rises[i]) for i, c in enumerate(xl)]
        nz = [np.array(c[:2]) for c in cuts]
        # Where the keel is on the hull it is on the centreline: carry it on
        # 1 in along its last segment, so the section plainly crosses the
        # centreline, between stations too, and camber trims it there.
        if cuts[-1][2]:
            d = nz[-1] - nz[-2]
            nz[-1] = nz[-1] + 1.0 * d / np.hypot(*d)
        return [tuple(map(float, q)) for q in nz]

    return dict(ctrl=ctrl, plan=plan, u_at=u_at, section=section, knuckled=knuckled,
                n_points=len(xl))


def camber_chain(vals, t):
    """camber's Catmull-Rom chain through `vals` (S, ...) on knots t, with no
    knuckles but at its ends (Loft::new's chains): Bezier segments (S-1, 4, ...)."""
    def tangent(p0, p1, p2, t0, t1, t2):
        return (p1 - p0) / (t1 - t0) - (p2 - p0) / (t2 - t0) + (p2 - p1) / (t2 - t1)
    n, segs = len(vals), []
    for j in range(n - 1):
        p, q = vals[j], vals[j + 1]
        chord, dt = q - p, t[j + 1] - t[j]
        i1 = tangent(vals[j - 1], p, q, t[j - 1], t[j], t[j + 1]) * dt if j > 0 else None
        i2 = tangent(p, q, vals[j + 2], t[j], t[j + 1], t[j + 2]) * dt if j + 2 < n else None
        m1 = i1 if i1 is not None else 2 * chord - i2 if i2 is not None else chord
        m2 = i2 if i2 is not None else 2 * chord - i1 if i1 is not None else chord
        segs.append(np.stack([p, p + m1 / 3, q - m2 / 3, q]))
    return np.array(segs)


def camber_loft(st_u, vals, u):
    """camber's stations (Loft::at): the points at each u, lofted from `vals`
    (S, ...) at st_u. Linear in vals."""
    st_u = np.asarray(st_u)
    segs = camber_chain(vals, st_u)
    uc = np.clip(u, st_u[0], st_u[-1])
    j = np.clip(np.searchsorted(st_u, uc) - 1, 0, len(st_u) - 2)
    s = (uc - st_u[j]) / (st_u[j + 1] - st_u[j])
    w = np.stack([(1 - s) ** 3, 3 * (1 - s) ** 2 * s, 3 * (1 - s) * s * s, s ** 3], 1)
    return np.einsum('ab,ab...->a...', w, segs[j])


def camber_truth(setup, extra_u=()):
    """The hull cut at a fine spread of u (closer at the bow), as camber's
    stations would be: u, (u, P, 2) points and the weight each point gets
    in a fit: 1 on the hull (and the keel just past it), 0.05 past the
    centreline, where only where it crosses matters."""
    u = np.unique(np.r_[np.linspace(0, 1, 1201), np.linspace(0.95, 1, 201), list(extra_u)])
    pts = np.array([setup['section'](v) for v in u])
    d = setup['plan'].derivative()(u)
    n_centre = setup['plan'](u)[:, 1] / (d[:, 0] / np.hypot(d[:, 0], d[:, 1]))
    weight = np.where(pts[:, :, 0] <= n_centre[:, None] + 1.1, 1.0, 0.05)
    return u, pts, weight


def fit_stations(setup, station_u, truth=None, dry=True):
    """Station points whose camber loft is closest (least squares) to the
    hull, rather than the hull's own cuts at the stations, which a few
    stations loft far from. With `dry`, the chine's two edges are held at or
    above the waterline everywhere along the loft: camber's Catmull-Rom
    dips below its stations where the chine comes down onto the waterline."""
    from scipy.optimize import minimize
    u, pts, weight = truth or camber_truth(setup)
    S = len(station_u)
    A = camber_loft(station_u, np.eye(S), u)                      # (U, S)
    V = np.zeros((S,) + pts.shape[1:])
    for p in range(pts.shape[1]):
        w = np.sqrt(weight[:, p])
        for c in range(2):
            V[:, p, c] = np.linalg.lstsq(A * w[:, None], pts[:, p, c] * w, rcond=None)[0]
            if dry and c == 1 and p in setup['knuckled'] and (A @ V[:, p, c]).min() < 0:
                b = pts[:, p, c]
                res = minimize(lambda v: np.sum((A @ v - b) ** 2), V[:, p, c],
                               jac=lambda v: 2 * A.T @ (A @ v - b), method='SLSQP',
                               constraints=[{'type': 'ineq', 'fun': lambda v: A @ v,
                                             'jac': lambda v: A}])
                V[:, p, c] = res.x
    return V


def camber_error(setup, station_u, truth, V=None):
    """How far camber's loft of these stations strays from the hull: the
    distance of each lofted point from the cut it stands for, on the hull."""
    u, pts, weight = truth
    if V is None:
        V = fit_stations(setup, station_u, truth, dry=False)
    d = np.linalg.norm(camber_loft(station_u, V, u) - pts, axis=2)
    return d[weight == 1.0]


def camber_search(setup, min_gap=4.0):
    """Fewer stations: from the dense set, drop the station whose loss
    lofts closest, then slide each of the rest along the plan while that
    helps, keeping them `min_gap` in apart; down to four. Prints each count's
    error and stations."""
    dense = sorted({0.0, 1.0} | {round(setup['u_at'](s), 7) for s in DENSE_S})
    truth = camber_truth(setup, dense)
    grid = truth[0]
    x = setup['plan'](grid)[:, 0]
    err = lambda sel: camber_error(setup, grid[sel], truth).max()   # noqa: E731
    sel = [int(np.argmin(abs(grid - v))) for v in dense]
    steps = np.r_[-40, -20, -10, -5, -2, -1, 1, 2, 5, 10, 20, 40]
    while len(sel) > 4:
        gaps = np.diff(x[sel])
        tight = [k for k in range(1, len(sel) - 1) if min(gaps[k - 1], gaps[k]) < min_gap]
        _, k = min((err(sel[:k] + sel[k + 1:]), k) for k in (tight or range(1, len(sel) - 1)))
        sel = sel[:k] + sel[k + 1:]
        best, moved = err(sel), True
        while moved:
            moved = False
            for k in range(1, len(sel) - 1):
                for c in sorted(set(np.clip(sel[k] + steps, sel[k - 1] + 1, sel[k + 1] - 1))):
                    if min(x[c] - x[sel[k - 1]], x[sel[k + 1]] - x[c]) < min_gap:
                        continue
                    e = err(sel[:k] + [int(c)] + sel[k + 1:])
                    if e < best - 1e-6:
                        best, sel, moved = e, sel[:k] + [int(c)] + sel[k + 1:], True
        d = camber_error(setup, grid[sel], truth)
        print(f"{len(sel):3d} stations: max {d.max():.3f} in, mean {d.mean():.3f} in; "
              f"u = {[round(float(v), 4) for v in grid[sel]]}", file=sys.stderr, flush=True)


def camber_sections(setup, st_u, V, ks, u, per_seg=16, nv=80, curves=False):
    """camber's swept hull, starboard, at plan parameters u: its stations
    lofted (camber_loft), each section drawn through its points (camber's
    centripetal Catmull-Rom, creased by the knuckles ks) in the plane square
    to the plan, and trimmed at the centreline. (len(u), nv, 3) points, even
    in arc length down each section, in the IGES frame (inches, x forward
    from the transom, z up from the waterline); NaN where there's no hull.
    With `curves`, also each section's whole (n, z) curve, sampled, and how
    many of its samples are kept."""
    pts = camber_loft(st_u, V, u)                                   # (U, P, 2)
    U, n = pts.shape[:2]
    t = np.c_[np.zeros(U), np.cumsum(np.sqrt(np.linalg.norm(np.diff(pts, axis=1), axis=2)).clip(1e-6), axis=1)]

    def tangent(a, b, c, ta, tb, tc):
        return ((b - a) / (tb - ta)[:, None] - (c - a) / (tc - ta)[:, None]
                + (c - b) / (tc - tb)[:, None])
    segs = np.zeros((U, n - 1, 4, 2))
    for j in range(n - 1):                                          # camber's cr_chain
        p, q = pts[:, j], pts[:, j + 1]
        chord, dt = q - p, (t[:, j + 1] - t[:, j])[:, None]
        i1 = tangent(pts[:, j - 1], p, q, t[:, j - 1], t[:, j], t[:, j + 1]) * dt if j > 0 else None
        i2 = tangent(p, q, pts[:, j + 2], t[:, j], t[:, j + 1], t[:, j + 2]) * dt if j + 2 < n else None
        e1 = 2 * chord - i2 if i2 is not None else chord
        e2 = 2 * chord - i1 if i1 is not None else chord
        cp = 1.0 if j == 0 else float(np.clip(ks[j], 0, 1))
        cq = 1.0 if j + 2 == n else float(np.clip(ks[j + 1], 0, 1))
        m1 = i1 if i1 is not None else e1
        m2 = i2 if i2 is not None else e2
        t1, t2 = e1 + (chord - e1) * cq, e2 + (chord - e2) * cp
        segs[:, j] = np.stack([p, p + (m1 + (t1 - m1) * cp) / 3, q - (m2 + (t2 - m2) * cq) / 3, q], 1)
    s = np.linspace(0, 1, per_seg + 1)[:-1]
    w = np.stack([(1 - s) ** 3, 3 * (1 - s) ** 2 * s, 3 * (1 - s) * s * s, s ** 3], 1)
    cur = np.concatenate([np.einsum('sk,uqkc->uqsc', w, segs).reshape(U, -1, 2), pts[:, -1:]], 1)
    P, d = setup['plan'](u), setup['plan'].derivative()(u)
    d = d / np.hypot(d[:, :1], d[:, 1:])
    x = P[:, :1] + cur[..., 0] * d[:, 1:]
    y = P[:, 1:] - cur[..., 0] * d[:, :1]
    out = np.full((U, nv, 3), np.nan)
    kept_to = np.zeros(U, int)
    for i in range(U):
        # camber keeps a section from its first point on the hull to its last
        # (kept_span), so a section that strays past the centreline and back
        # keeps the stray; it ends where it last crosses.
        if y[i, 0] < 0:
            continue
        kept = np.flatnonzero(y[i] >= 0)
        m = kept_to[i] = kept[-1] + 1
        q = np.c_[x[i, :m], y[i, :m], cur[i, :m, 1]]
        if m < y.shape[1]:                                          # to the centreline
            a = y[i, m - 1] / (y[i, m - 1] - y[i, m])
            q = np.vstack([q, q[-1] + a * (np.r_[x[i, m], y[i, m], cur[i, m, 1]] - q[-1])])
        L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]
        r = np.linspace(0, L[-1], nv)
        out[i] = np.stack([np.interp(r, L, q[:, k]) for k in range(3)], 1)
    return (out, cur, kept_to) if curves else out


def grid_normals(G):
    n = np.cross(np.gradient(G, axis=0), np.gradient(G, axis=1))
    return n / np.linalg.norm(n, axis=-1, keepdims=True).clip(1e-12)


def hull_cloud(nu=4000):
    """The lofted hull's starboard side (bottom, ledge, topside) as points and
    normals, in the IGES frame (inches)."""
    from scipy.interpolate import BSpline
    P, N = [], []
    for name, nv in (('bottom', 120), ('ledge', 8), ('topside', 60)):
        srf = pieces[name]
        Bu = BSpline.design_matrix(np.linspace(0, 1, nu), srf['ku'], srf['du']).toarray()
        Bv = BSpline.design_matrix(np.linspace(0, 1, nv), srf['kv'], srf['dv']).toarray()
        G = np.stack([Bu @ srf['ctrl'][..., k] @ Bv.T for k in range(3)], -1)
        G[..., 0] = LWL - G[..., 0]
        P.append(G.reshape(-1, 3))
        N.append(grid_normals(G).reshape(-1, 3))
    return np.vstack(P), np.vstack(N)


def refine_stations(setup, st_u, V0, max_nfev=100):
    """Station points (and the bottom points' knuckles) moved so camber's
    swept surface itself lies on the hull, by nonlinear least squares from
    V0: each swept point's distance along the hull's normal, and each hull
    point's along the swept surface's.

    Each station is built so it reads as a section should, deck to keel:
    the deck point at the deck, then each point inboard of the last and no
    higher (steps through softplus), with the chine's two edges at one
    height, so the ledge is level, as in the patent. Penalties hold the
    curves between the points to the same (camber's Catmull-Rom can
    overshoot): the whole of each station's curve, and the kept part of
    the sections lofted between them. More penalties hold the chine at or
    above the waterline, the keel past the centreline and the surface off
    the port side. Returns (V, ks, deviation), deviation being the two-way
    distances."""
    from scipy.optimize import least_squares
    from scipy.spatial import cKDTree
    R, RN = hull_cloud()
    tree = cKDTree(R)
    sub = np.random.default_rng(0).choice(len(R), 20000, replace=False)
    Rs = R[sub]
    S, P = V0.shape[:2]
    u = np.unique(np.r_[np.linspace(0, 1, 400), np.linspace(0.93, 1, 133), st_u])
    at_station = np.isin(u, st_u)
    d = setup['plan'].derivative()(u)
    n_centre = setup['plan'](u)[:, 1] / (d[:, 0] / np.hypot(d[:, 0], d[:, 1]))
    k_bottom = list(range(3, P - 1))
    ks0 = np.array([1.0 if i in setup['knuckled'] else 0.0 for i in range(P)])
    soft = lambda x: np.logaddexp(0.0, x)                                  # noqa: E731
    unsoft = lambda y: y + np.log(-np.expm1(-y))                           # noqa: E731

    def build(p):
        """Stations (S, P, 2) from parameters (S, 2P - 2): deck n, the
        ledge's drop below the deck, the chine edges' steps in, then each
        later point's step in and drop."""
        V = np.zeros((S, P, 2))
        V[:, 0] = np.c_[p[:, 0], np.full(S, DECK)]
        zc = DECK - soft(p[:, 1])
        V[:, 1] = np.c_[V[:, 0, 0] + soft(p[:, 2]), zc]
        V[:, 2] = np.c_[V[:, 1, 0] + soft(p[:, 3]), zc]
        for i in range(3, P):
            a, b = p[:, 2 * i - 2], p[:, 2 * i - 1]
            V[:, i] = np.c_[V[:, i - 1, 0] + soft(a), V[:, i - 1, 1] - soft(b)]
        return V

    def seed(V):
        """The parameters nearest V that build a monotone station."""
        p = np.zeros((S, 2 * P - 2))
        p[:, 0] = V[:, 0, 0]
        zc = np.clip(V[:, 1:3, 1].mean(1), 0.05, DECK - 0.05)
        p[:, 1] = unsoft(DECK - zc)
        n_prev, z_prev = V[:, 0, 0], zc
        for i in range(1, P):
            n_i = np.maximum(V[:, i, 0], n_prev + 0.02)
            if i <= 2:
                p[:, i + 1] = unsoft(n_i - n_prev)
            else:
                z_i = np.minimum(V[:, i, 1], z_prev - 0.02)
                p[:, 2 * i - 2], p[:, 2 * i - 1] = unsoft(n_i - n_prev), unsoft(z_prev - z_i)
                z_prev = z_i
            n_prev = n_i
        return p

    npar = S * (2 * P - 2)

    def unpack(x):
        ks = ks0.copy()
        ks[k_bottom] = x[npar:]
        return build(x[:npar].reshape(S, 2 * P - 2)), ks

    def distances(G):
        """Signed distances: each swept point's from the hull (0 where there's
        no section), and each hull point's from the swept surface."""
        g = G.reshape(-1, 3)
        ok = ~np.isnan(g[:, 0])
        _, i = tree.query(g[ok])
        there = np.zeros(len(g))
        there[ok] = np.einsum('ij,ij->i', g[ok] - R[i], RN[i])
        _, j = cKDTree(g[ok]).query(Rs)
        back = np.einsum('ij,ij->i', Rs - g[ok][j], grid_normals(np.nan_to_num(G)).reshape(-1, 3)[ok][j])
        return there, back

    def unmonotone(cur, kept_to):
        """How far each step down each curve rises or turns outboard: over
        the whole curve at the stations, over the kept part between."""
        dn, dz = np.diff(cur[..., 0], axis=1), np.diff(cur[..., 1], axis=1)
        idx = np.arange(dn.shape[1])
        use = at_station[:, None] | (idx[None, :] < (kept_to - 1)[:, None])
        return (np.maximum(0, dz) + np.maximum(0, -dn)) * use

    def residuals(x):
        V, ks = unpack(x)
        G, cur, kept_to = camber_sections(setup, st_u, V, ks, u, curves=True)
        there, back = distances(G)
        port = np.nan_to_num(np.maximum(0, -G[..., 1])).ravel()
        keel = camber_loft(st_u, V, u)[:, -1, 0]
        chine = camber_loft(st_u, V[:, 1, 1], u)
        return np.r_[there, back, 10 * port, 300 * unmonotone(cur, kept_to).ravel(),
                     10 * np.maximum(0, n_centre + 0.3 - keel),
                     30 * np.maximum(0, 0.01 - chine)]

    x0 = np.r_[seed(V0).ravel(), np.zeros(len(k_bottom))]
    lo = np.r_[np.full(npar, -np.inf), np.zeros(len(k_bottom))]
    hi = np.r_[np.full(npar, np.inf), np.ones(len(k_bottom))]
    res = least_squares(residuals, x0, bounds=(lo, hi), diff_step=1e-4, max_nfev=max_nfev)
    V, ks = unpack(res.x)
    G = camber_sections(setup, st_u, V, ks, u)
    there, back = distances(G)
    return V, ks, np.abs(np.r_[there[~np.isnan(G[..., 0]).ravel()], back])


def station_monotone(setup, st_u, V, ks):
    """The worst rise and the worst outboard turn (inches) down any station's
    whole curve, deck to keel, sampled between its points."""
    _, cur, _ = camber_sections(setup, st_u, V, ks, np.asarray(st_u, float), curves=True)
    dn, dz = np.diff(cur[..., 0], axis=1), np.diff(cur[..., 1], axis=1)
    return float(np.maximum(0, dz).max()), float(np.maximum(0, -dn).max())


# Stations for the editable document, from camber_search (--camber-search),
# their points fitted (fit_stations). Seven is the fewest that keep the form:
# against the IGES, volume, wetted area and LCB agree to 0.1% and 0.02 in at
# five stations or more, and the beam to 0.2 in, but below seven the stem
# goes (0.7 in out near the waterline at seven, 1.2 at six) and the waterline
# lengthens (0.6 in at seven, 0.9 at six). The transom, the run, and the rest
# at the bow, where the chine climbs and the stem rakes.
CAMBER_U = [0.0, 0.3467, 0.8108, 0.8958, 0.9595, 0.9805, 1.0]
# Bottom points per station, between the inner chine and the keel (each
# station also has the deck, the chine's two edges and the keel). Nine are
# the v-lines' own points, fitted; fewer are refined so camber's swept surface
# lies on the hull, every section monotone (refine_stations). At seven
# stations one is enough: the surface within 0.54 in (99% within 0.34),
# volume, wetted area and LCB to 0.14% and 0.02 in, beam +0.05 in, waterline
# +0.45 in, draft +0.26 in. Two come to 0.45 in (99% within 0.31); more gain
# little. At five stations two are best (0.61 in, 99% within 0.35).
CAMBER_POINTS = 1


# The 39 stations the first export used: dense at the bow and where the chine
# comes down onto the waterline (behind station 7), so camber's Catmull-Rom
# loft can't swing between them.
DENSE_S = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2.5, 4, 6, 9, 12.5, 17.9, 24, 30, 34.2, 45,
           55.3, 66, 76.9, 87, 94.2, 98.2, 120, 140.9, 150, 156, 160, 162.7, 166, 170,
           185, 205.8, 227, LWL]


def camber_document(setup, station_u, fit=True, V=None, ks=None):
    """The hull as a camber (v2) document with stations at these plan u:
    their points V (and knuckles ks) as given, fitted to the hull
    (fit_stations), or the hull's own cuts."""
    ctrl, knuckled = setup['ctrl'], setup['knuckled']
    station_u = sorted(station_u)
    if V is None and fit:
        V = fit_stations(setup, station_u)
    if ks is None:
        ks = [1 if i in knuckled else 0 for i in range(setup['n_points'])]
    stations = []
    for j, u in enumerate(station_u):
        nz = V[j] if V is not None else setup['section'](u)
        pts = [{"k": round(float(ks[i]), 4), "n": round(float(n), 5), "z": round(float(z) - DECK, 5)}
               for i, (n, z) in enumerate(nz)]
        stations.append({"u": round(float(u), 7), "keelK": 0, "points": pts})
    plan_x = float(ctrl[-1, 0])
    return {
        "name": "Graf catamaran demihull (US 5,188,049)",
        "version": 2,
        "unit": "in",
        "sheerPlan": [{"x": round(float(x), 5), "y": round(float(y), 5)} for x, y in ctrl],
        "sheerTrim": [{"k": 0, "x": 0, "z": 0.5}, {"k": 0, "x": plan_x, "z": 0.5}],
        # camber's transom plane is just behind the plan's start, so it trims
        # nothing: the hull ends square at the start, and the end is closed
        # flat across. (A transom cutting square-on sections leaves a sliver
        # of hull behind it that fits badly.)
        "transom": [{"x": -0.5, "z": 1.0}, {"x": -0.5, "z": -30.0}],
        # 1 mm below the stern chine, so the chine ledge stays dry
        "waterline": DECK + 0.04,
        "deckTrimDeg": 0,
        "stations": stations,
    }


if any(a.startswith("--camber") for a in sys.argv):
    n_bottom = (int(sys.argv[sys.argv.index("--camber-points") + 1])
                if "--camber-points" in sys.argv else CAMBER_POINTS)
    setup = camber_setup(n_bottom)
    if "--camber-search" in sys.argv:
        camber_search(setup)
        sys.exit()
    if "--camber-dense" in sys.argv:
        us = {0.0, 1.0} | {round(setup['u_at'](s), 7) for s in DENSE_S}
        doc = camber_document(setup, us, fit=False)
    else:
        us = ([float(v) for v in sys.argv[sys.argv.index("--camber-u") + 1].split(",")]
              if "--camber-u" in sys.argv else CAMBER_U)    # e.g. from --camber-search
        if n_bottom == 9:                   # the v-lines' own points, fitted
            doc = camber_document(setup, us)
        else:                               # fewer: refined to the surface
            V, ks, dev = refine_stations(setup, us, fit_stations(setup, us))
            print(f"# {len(us)} stations, {n_bottom} bottom points: surface within "
                  f"{dev.max():.3f} in (99% within {np.percentile(dev, 99):.3f})", file=sys.stderr)
            doc = camber_document(setup, us, V=V, ks=ks)
    json.dump(doc, sys.stdout, indent=1)
    sys.exit()


# ---------------------------------------------------------------- IGES out
def fmt(v):
    return f"{v:.9g}" if not float(v).is_integer() else f"{v:.1f}"


def hollerith(s):
    return f"{len(s)}H{s}"


def lines(text, width, tag):
    out, cur = [], ""
    for tok in text.split(","):
        piece = tok + ("," if not tok.endswith(";") else "")
        if len(cur) + len(piece) > width:
            out.append(cur)
            cur = ""
        cur += piece
    if cur:
        out.append(cur)
    return out


start = ["Graf catamaran demihull, US 5,188,049 (1993), lofted from Figs. 2 and 5"]
glob = ",".join([
    "1H,", "1H;", hollerith("graf-demihull"), hollerith("graf-demihull.igs"),
    hollerith("boatmath graf-loft.py"), hollerith("1"), "32", "38", "6", "308", "15",
    hollerith("graf-demihull"), "1.0", "6", hollerith("M"), "1", "1.0",
    hollerith("20261012.000000"), "1.0E-6", "10.0", hollerith(""), hollerith(""),
    "11", "0", hollerith("20261012.000000")]) + ";"

S = [f"{l:<72}S{i + 1:7d}" for i, l in enumerate(start)]
G = [f"{l:<72}G{i + 1:7d}" for i, l in enumerate(lines(glob, 72, "G"))]
D, P = [], []
for n, srf in enumerate(patches):
    nu, nv = srf['ctrl'].shape[:2]
    ctrl = srf['ctrl']
    params = [128, nu - 1, nv - 1, srf['du'], srf['dv'], 0, 0, 1, 0, 0]
    params += list(srf['ku']) + list(srf['kv'])
    params += [1.0] * (nu * nv)
    for j in range(nv):                           # u fastest
        for i in range(nu):
            params += list(ctrl[i, j])
    params += [0.0, 1.0, 0.0, 1.0]
    text = ",".join(fmt(p) for p in params) + ";"
    de = 2 * n + 1
    pl = lines(text, 64, "P")
    pstart = len(P) + 1
    for l in pl:
        P.append(f"{l:<64}{de:8d}P{len(P) + 1:7d}")
    D.append(f"{128:8d}{pstart:8d}{0:8d}{0:8d}{0:8d}{0:8d}{0:8d}{0:8d}{'00000000':>8}D{de:7d}")
    D.append(f"{128:8d}{0:8d}{0:8d}{len(pl):8d}{0:8d}{'':8}{'':8}{'':8}{0:8d}D{de + 1:7d}")
T = [f"S{len(S):7d}G{len(G):7d}D{len(D):7d}P{len(P):7d}{'':40}T{1:7d}"]
sys.stdout.write("\n".join(S + G + D + P + T) + "\n")

deck_hb = max(r[-1][1] for r in topside)
print(f"# {len(patches)} patches; widest at the deck {2 * deck_hb:.2f} in, so a 36 in "
      f"tunnel puts the hull centres {36 + 2 * deck_hb:.1f} in apart", file=sys.stderr)
