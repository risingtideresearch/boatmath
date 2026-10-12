//! Measuring a hull's surface: areas and centroids of regions of its shell,
//! and transverse slices through it (bulkheads, frames, tank sections).
//!
//! Both work on the triangle mesh hullgeom keeps of each hull: an STL's own
//! facets, or an IGES hull's patches tessellated at about 1/400 of its size.
//! Facets are flat, so on a curved patch the area is that of the inscribed
//! mesh, short of the surface's by roughly (h/R)²/8 for facets of size h on a
//! radius R: 0.1–0.3% on a small hull's bilge, much less elsewhere.
//!
//! Coordinates are the hull file's (napkin's CG frame): x along the hull,
//! z up in the file's frame, y from the platform's centreline.

/// One hull of the platform as triangles, in the file frame.
pub struct Mesh {
    pub verts: Vec<[f64; 3]>,
    pub tris: Vec<[u32; 3]>,
    pub lo: [f64; 3],
    pub hi: [f64; 3],
}

impl Mesh {
    /// A mesh from posed vertices; a one-sided (half-hull) shell is mirrored
    /// about its centreplane `yc` so the mesh is the whole hull.
    pub fn new(mut verts: Vec<[f64; 3]>, mut tris: Vec<[u32; 3]>, yc: f64) -> Mesh {
        let (ylo, yhi) = verts
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), p| {
                (a.min(p[1]), b.max(p[1]))
            });
        let tol = 1e-6 * (yhi - ylo).abs().max(1e-9);
        if ylo >= yc - tol || yhi <= yc + tol {
            let n = verts.len() as u32;
            let mirrored: Vec<[f64; 3]> = verts
                .iter()
                .map(|p| [p[0], 2.0 * yc - p[1], p[2]])
                .collect();
            verts.extend(mirrored);
            let more: Vec<[u32; 3]> = tris
                .iter()
                .map(|t| [t[0] + n, t[2] + n, t[1] + n])
                .collect();
            tris.extend(more);
        }
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &verts {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        Mesh {
            verts,
            tris,
            lo,
            hi,
        }
    }

    fn extent(&self) -> f64 {
        (0..3).map(|k| self.hi[k] - self.lo[k]).fold(0.0, f64::max)
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Clip a convex polygon to `p[axis] ≥ v` (`keep_above`) or `≤ v`.
fn clip(poly: &[[f64; 3]], axis: usize, v: f64, keep_above: bool) -> Vec<[f64; 3]> {
    let inside = |p: &[f64; 3]| {
        if keep_above {
            p[axis] >= v
        } else {
            p[axis] <= v
        }
    };
    let mut out = Vec::with_capacity(poly.len() + 2);
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        let (ia, ib) = (inside(&a), inside(&b));
        if ia {
            out.push(a);
        }
        if ia != ib {
            let t = (v - a[axis]) / (b[axis] - a[axis]);
            out.push([
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ]);
        }
    }
    out
}

/// A region of space: per axis, `[lo, hi]` (infinite for unbounded).
pub type Bounds = [[f64; 2]; 3];

/// Area and area centroid of the shell inside `b`, over `meshes`.
pub fn shell(meshes: &[&Mesh], b: &Bounds) -> (f64, [f64; 3]) {
    let (mut area, mut m) = (0.0, [0.0; 3]);
    for mesh in meshes {
        if (0..3).any(|k| mesh.hi[k] < b[k][0] || mesh.lo[k] > b[k][1]) {
            continue;
        }
        for t in &mesh.tris {
            let tri = [
                mesh.verts[t[0] as usize],
                mesh.verts[t[1] as usize],
                mesh.verts[t[2] as usize],
            ];
            // Whole triangles inside need no clipping; whole ones outside are skipped.
            let mut all_in = true;
            let mut out = false;
            for k in 0..3 {
                let (lo, hi) = (
                    tri.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min),
                    tri.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max),
                );
                if hi < b[k][0] || lo > b[k][1] {
                    out = true;
                }
                if lo < b[k][0] || hi > b[k][1] {
                    all_in = false;
                }
            }
            if out {
                continue;
            }
            let mut poly = tri.to_vec();
            if !all_in {
                for k in 0..3 {
                    if b[k][0].is_finite() {
                        poly = clip(&poly, k, b[k][0], true);
                    }
                    if b[k][1].is_finite() {
                        poly = clip(&poly, k, b[k][1], false);
                    }
                    if poly.len() < 3 {
                        break;
                    }
                }
            }
            for i in 1..poly.len().saturating_sub(1) {
                let a = 0.5 * norm(cross(sub(poly[i], poly[0]), sub(poly[i + 1], poly[0])));
                area += a;
                for k in 0..3 {
                    m[k] += a * (poly[0][k] + poly[i][k] + poly[i + 1][k]) / 3.0;
                }
            }
        }
    }
    if area > 0.0 {
        (area, [m[0] / area, m[1] / area, m[2] / area])
    } else {
        (0.0, [f64::NAN; 3])
    }
}

/// One closed outline of a slice in its (y, z) plane. `shell[i]` says whether
/// the edge from point i to i+1 is the hull's surface (else it closes an open
/// section across its top, as a flush deck at the sheer would).
#[derive(Clone, Debug)]
pub struct Outline {
    pub pts: Vec<[f64; 2]>,
    pub shell: Vec<bool>,
}

/// The hulls' sections at station `x`: closed outlines.
pub fn slice(meshes: &[&Mesh], x: f64) -> Vec<Outline> {
    let mut out = vec![];
    for mesh in meshes {
        if x < mesh.lo[0] || x > mesh.hi[0] {
            continue;
        }
        out.extend(slice_mesh(mesh, x));
    }
    out
}

fn slice_mesh(mesh: &Mesh, x: f64) -> Vec<Outline> {
    // Segments where the plane crosses each triangle. A vertex exactly on the
    // plane counts as ahead of it, so every crossing is of an edge's interior.
    let side = |i: u32| mesh.verts[i as usize][0] >= x;
    let at = |i: u32, j: u32| {
        let (a, b) = (mesh.verts[i as usize], mesh.verts[j as usize]);
        let t = (x - a[0]) / (b[0] - a[0]);
        [a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])]
    };
    let mut segs: Vec<[[f64; 2]; 2]> = vec![];
    for t in &mesh.tris {
        let s = [side(t[0]), side(t[1]), side(t[2])];
        if s[0] == s[1] && s[1] == s[2] {
            continue;
        }
        let mut pts = vec![];
        for k in 0..3 {
            let (i, j) = (t[k], t[(k + 1) % 3]);
            if s[k] != s[(k + 1) % 3] {
                pts.push(at(i, j));
            }
        }
        if pts.len() == 2 && (pts[0][0] - pts[1][0]).abs() + (pts[0][1] - pts[1][1]).abs() > 0.0 {
            segs.push([pts[0], pts[1]]);
        }
    }
    if segs.is_empty() {
        return vec![];
    }
    // Weld endpoints within a tolerance: patches meet along seams whose
    // tessellations differ slightly.
    let tol = 1e-4 * mesh.extent();
    let cell = |p: [f64; 2]| ((p[0] / tol).floor() as i64, (p[1] / tol).floor() as i64);
    let mut grid: std::collections::HashMap<(i64, i64), Vec<usize>> = Default::default();
    let mut nodes: Vec<[f64; 2]> = vec![];
    let mut node_of = |p: [f64; 2], nodes: &mut Vec<[f64; 2]>| -> usize {
        let (cx, cy) = cell(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(ids) = grid.get(&(cx + dx, cy + dy)) {
                    for &id in ids {
                        let q = nodes[id];
                        if (q[0] - p[0]).hypot(q[1] - p[1]) <= tol {
                            return id;
                        }
                    }
                }
            }
        }
        nodes.push(p);
        grid.entry((cx, cy)).or_default().push(nodes.len() - 1);
        nodes.len() - 1
    };
    let mut adj: Vec<Vec<usize>> = vec![];
    let mut edges: Vec<(usize, usize)> = vec![];
    for s in &segs {
        let a = node_of(s[0], &mut nodes);
        let b = node_of(s[1], &mut nodes);
        if a == b {
            continue;
        }
        adj.resize(nodes.len(), vec![]);
        if adj[a].contains(&b) {
            continue;
        }
        adj[a].push(b);
        adj[b].push(a);
        edges.push((a, b));
    }
    adj.resize(nodes.len(), vec![]);
    // Walk chains: from open ends first, then whatever loops remain.
    let mut used = std::collections::HashSet::new();
    let key = |a: usize, b: usize| (a.min(b), a.max(b));
    let mut chains: Vec<(Vec<usize>, bool)> = vec![];
    let starts: Vec<usize> = (0..nodes.len())
        .filter(|&n| adj[n].len() == 1)
        .chain(0..nodes.len())
        .collect();
    for s in starts {
        while let Some(&next) = adj[s].iter().find(|&&m| !used.contains(&key(s, m))) {
            let mut chain = vec![s];
            let (mut prev, mut cur) = (s, next);
            used.insert(key(prev, cur));
            loop {
                chain.push(cur);
                match adj[cur]
                    .iter()
                    .find(|&&m| m != prev && !used.contains(&key(cur, m)))
                {
                    Some(&n) => {
                        used.insert(key(cur, n));
                        prev = cur;
                        cur = n;
                    }
                    None => break,
                }
            }
            let closed = chain.len() > 2 && chain.first() == chain.last();
            if closed {
                chain.pop();
            }
            chains.push((chain, closed));
        }
    }
    // Join open chains end to end across small gaps, nearest first.
    let gap = 1e-2 * mesh.extent();
    let mut open: Vec<Vec<usize>> = chains
        .iter()
        .filter(|c| !c.1)
        .map(|c| c.0.clone())
        .collect();
    let mut loops: Vec<Outline> = chains
        .iter()
        .filter(|c| c.1)
        .map(|c| Outline {
            pts: c.0.iter().map(|&n| nodes[n]).collect(),
            shell: vec![true; c.0.len()],
        })
        .collect();
    loop {
        let mut best: Option<(f64, usize, usize, bool, bool)> = None;
        for i in 0..open.len() {
            for j in 0..open.len() {
                if i == j {
                    continue;
                }
                for (ei, a) in [(false, open[i][0]), (true, *open[i].last().unwrap())] {
                    for (ej, b) in [(false, open[j][0]), (true, *open[j].last().unwrap())] {
                        let d = (nodes[a][0] - nodes[b][0]).hypot(nodes[a][1] - nodes[b][1]);
                        if d < gap && best.is_none_or(|bb| d < bb.0) {
                            best = Some((d, i, j, ei, ej));
                        }
                    }
                }
            }
        }
        let Some((_, i, j, ei, ej)) = best else { break };
        let (mut a, mut b) = (open[i].clone(), open[j].clone());
        if !ei {
            a.reverse();
        }
        if ej {
            b.reverse();
        }
        a.extend(b);
        let (hi, lo) = (i.max(j), i.min(j));
        open.remove(hi);
        open.remove(lo);
        open.push(a);
    }
    // Close each open section straight across its ends.
    for c in open {
        if c.len() < 3 {
            continue;
        }
        let mut shell = vec![true; c.len()];
        *shell.last_mut().unwrap() = false;
        loops.push(Outline {
            pts: c.iter().map(|&n| nodes[n]).collect(),
            shell,
        });
    }
    loops
}

impl Outline {
    /// The part of the outline at or below height `z` (a section to a deck,
    /// or under the waterline); the cut is a closing edge, not shell.
    pub fn below(&self, z: f64) -> Option<Outline> {
        let n = self.pts.len();
        let (mut pts, mut shell) = (vec![], vec![]);
        for i in 0..n {
            let (a, b) = (self.pts[i], self.pts[(i + 1) % n]);
            let (ia, ib) = (a[1] <= z, b[1] <= z);
            if ia {
                pts.push(a);
                shell.push(self.shell[i]);
            }
            if ia != ib {
                let t = (z - a[1]) / (b[1] - a[1]);
                pts.push([a[0] + t * (b[0] - a[0]), z]);
                // Leaving: the cut runs along z until it re-enters. Entering:
                // the rest of this edge is the original edge.
                shell.push(!ia && self.shell[i]);
            }
        }
        if pts.len() < 3 {
            None
        } else {
            Some(Outline { pts, shell })
        }
    }

    /// Area, centroid (y, z), and girth: the length of the shell edges.
    pub fn measure(&self) -> (f64, [f64; 2], f64) {
        let n = self.pts.len();
        let (mut a, mut cy, mut cz, mut girth) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            let (p, q) = (self.pts[i], self.pts[(i + 1) % n]);
            let c = p[0] * q[1] - q[0] * p[1];
            a += c;
            cy += (p[0] + q[0]) * c;
            cz += (p[1] + q[1]) * c;
            if self.shell[i] {
                girth += (q[0] - p[0]).hypot(q[1] - p[1]);
            }
        }
        a *= 0.5;
        if a.abs() < 1e-300 {
            return (0.0, [f64::NAN; 2], girth);
        }
        (a.abs(), [cy / (6.0 * a), cz / (6.0 * a)], girth)
    }
}
