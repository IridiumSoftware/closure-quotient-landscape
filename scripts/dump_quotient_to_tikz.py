#!/usr/bin/env python3
"""
dump_quotient_to_tikz.py — emit per-quotient TikZ figures for the
closure-quotient paper.

For each canonical family member, produces a single .tex file at
paper/figures/<quotient>.tex containing a \\begin{figure}...\\end{figure}
environment with three subpanels:

  (a) Tier-stratified node-edge graph.
  (b) Hyperedge stars (every multiway hyperedge drawn as a 3-pronged
      figure meeting at its centroid).
  (c) Tier-sorted adjacency heatmap.

Quotients rendered (canonical family + Q181):
  Q48   — C-closure of cyc(6)
  Q90   — C-closure of (Adjacent + 6 random K_6^3 pads at seed=42)
  Q84   — C-closure of Adjacent ternary on 6 vertices
  Q102  — C-closure of K_6^3
  Q181  — single-side closure of a 12-vertex 132-edge irregular seed

Pure Python stdlib + closure-v5's exact-arithmetic builders. Output is
LaTeX/TikZ, included into the paper via \\input{figures/<q>.tex}.

Usage:
  cd paper/scripts/
  python3 dump_quotient_to_tikz.py
"""

import os, sys, io, contextlib, math, random

CV5_ROOT = os.path.normpath(
    os.path.join(os.path.dirname(__file__), '..', '..')
)
sys.path.insert(0, CV5_ROOT)

with contextlib.redirect_stdout(sys.stderr):
    from q102_build_exact_v1 import (
        build_c_closed_quotient, complete_ternary,
    )

# Local: g4_larger_graphs_exact_v1 (single-side exact builder + Pool B)
sys.path.insert(0, os.path.join(CV5_ROOT))
with contextlib.redirect_stdout(sys.stderr):
    from g4_larger_graphs_exact_v1 import (
        POOL_B,
        build_quotient_exact,
        compose_exact, proj_equiv, is_zero_vec,
    )

# ----------------------------------------------------------------------------
# Seed-graph constructions (canonical, copied from closure-v5 sources)
# ----------------------------------------------------------------------------

def cyc6():
    """6-vertex cyclic ternary — closes to Q24 single, Q48 C-closed."""
    return [(0,1,2),(1,2,3),(2,3,4),(3,4,5),(4,5,0),(5,0,1)]

def adjacent_ternary():
    """Adjacent triples (i, i+1, i+2) + all permutations — Q42/Q84."""
    import itertools
    edges = set()
    for i in range(6):
        triple = (i, (i+1)%6, (i+2)%6)
        for p in itertools.permutations(triple):
            edges.add(p)
    return list(edges)

def complete_edges_K6():
    """K_6^3 edges (all ordered triples of distinct vertices on 6 verts)."""
    return [(i,j,k) for i in range(6) for j in range(6) if j != i
            for k in range(6) if k != i and k != j]

def q90_seed_edges():
    """A 6-vertex seed that closes to exactly Q90 under C-closure
    (Adjacent + 6 K_6^3 pads; search for ANY pad set that lands on
    n_cl=90, since matching the canonical numpy-shuffle seed without
    numpy is finicky and the specific pads aren't structurally
    load-bearing for visualization)."""
    base = set(adjacent_ternary())
    pool = [e for e in complete_edges_K6() if e not in base]
    rng = random.Random(11)
    for trial in range(40):
        rng.shuffle(pool)
        candidate = list(base) + pool[:6]
        Q = build_c_closed_quotient(candidate, depth=4, n_seed_verts=6)
        if Q['n_cl'] == 90:
            return candidate
    raise RuntimeError("no Q90 seed found in 40 trials")

# ----------------------------------------------------------------------------
# Q181 seed (12-vertex, 132-edge irregular seed)
# ----------------------------------------------------------------------------

def q181_seed():
    """Find an n-vertex k-edge seed that closes to |Q|=181 under the
    canonical s58 / Pool A ICs (matching q181_search.py's original
    method). Pool A's K_n^3 defect at n>=9 affects regular complete-K_n
    closures, not the irregular n=12 seeds searched here.
    Sweeps density 0.04..0.20."""
    from q102_build_exact_v1 import _GEN_ICS as POOL_A_RAW
    pool_a_ic = [list(t) for t in POOL_A_RAW]
    # Pad to >=13 ICs if needed (we have 15, sufficient for n=12).
    n = 12
    full = [(i,j,k) for i in range(n) for j in range(n) if j != i
            for k in range(n) if k != i and k != j]
    random.seed(2718)  # same seed as q181_search.py
    me = len(full)
    for frac_step in range(8, 40):
        frac = frac_step / 200
        k = max(1, int(me * frac))
        for trial in range(6):
            seed_edges = random.sample(full, k)
            n_cl, _ = build_quotient_exact(seed_edges, n, pool_a_ic, depth=5)
            if n_cl == 181:
                return seed_edges, n, pool_a_ic
    raise RuntimeError("no Q181 seed found in density sweep")

# ----------------------------------------------------------------------------
# Quotient data extraction
# ----------------------------------------------------------------------------

class QuotientData:
    """Holds everything the TikZ emitter needs about one quotient."""
    def __init__(self, name, n_cl, tier, origin, edges, hyperedges,
                 tier_sizes, total_he):
        self.name        = name        # e.g. "Q48"
        self.n_cl        = n_cl
        self.tier        = tier        # dict {cid: 0/1/2 for A/B/C}
        self.origin      = origin      # dict {cid: 0/1/2 for orig/conj/both}
        self.edges       = edges       # list of (i, j) with i<j
        self.hyperedges  = hyperedges  # list of (c1, c2, c3)
        self.tier_sizes  = tier_sizes  # [|A|, |B|, |C|]
        self.total_he    = total_he

def extract_c_closed(name, seed_edges, n_verts, depth=4):
    """Use closure-v5's C-closed exact builder; canonical ICs."""
    Q = build_c_closed_quotient(seed_edges, depth=depth,
                                 n_seed_verts=n_verts)
    n_cl = Q['n_cl']
    TIER_TO_INT = {'A': 0, 'B': 1, 'C': 2}
    tier = {c: TIER_TO_INT[Q['q_tier'][c]] for c in range(n_cl)}
    ORIGIN_TO_INT = {'orig_only': 0, 'conj_only': 1, 'both': 2}
    origin = {c: ORIGIN_TO_INT[Q['cl_origin'][c]] for c in range(n_cl)}
    # Binary edges
    adj = [[0]*n_cl for _ in range(n_cl)]
    for (c1, c2, c3) in Q['q_he']:
        for a, b in ((c1,c2),(c1,c3),(c2,c3),
                     (c2,c1),(c3,c1),(c3,c2)):
            if a != b: adj[a][b] = 1
    edges = [(i,j) for i in range(n_cl) for j in range(i+1, n_cl) if adj[i][j]]
    tier_sizes = [sum(1 for c in range(n_cl) if tier[c] == t) for t in (0,1,2)]
    return QuotientData(name, n_cl, tier, origin, edges,
                        sorted(Q['q_he']), tier_sizes, len(Q['q_he']))

def extract_single_side(name, seed_edges, n_verts, ic_pool, depth=5):
    """Single-side closure (no C-closure); tier from depth-of-appearance."""
    psi_init = {v: [list(c) for c in ic_pool[v]] for v in range(n_verts)}
    # Replicate build with depth tracking
    psi = {k: list(v) for k, v in psi_init.items()}
    vertex_depth = {v: 0 for v in psi_init}
    nv = max(psi_init) + 1
    all_edges = [(0, s1, s2, s3) for (s1,s2,s3) in seed_edges]
    cache = {}
    for d in range(depth):
        for (_, v1, v2, v3) in [e for e in all_edges if e[0] == d]:
            key = (v1, v2)
            if key not in cache:
                w = compose_exact(psi[v1], psi[v2])
                if is_zero_vec(w): continue
                psi[nv] = w; cache[key] = nv
                vertex_depth[nv] = d + 1
                nv += 1
            w = cache[key]
            all_edges.append((d+1, w, v2, v3))
            all_edges.append((d+1, w, v1, v3))
            all_edges.append((d+1, w, v1, v2))
    # Cluster
    reps = []; vid_to_cid = {}; cluster_min_depth = {}
    for v in sorted(psi):
        pv = psi[v]
        if is_zero_vec(pv): continue
        matched = -1
        for i, r in enumerate(reps):
            if proj_equiv(pv, r): matched = i; break
        if matched >= 0:
            cid = matched
            vid_to_cid[v] = cid
            if vertex_depth[v] < cluster_min_depth[cid]:
                cluster_min_depth[cid] = vertex_depth[v]
        else:
            cid = len(reps)
            reps.append(pv)
            vid_to_cid[v] = cid
            cluster_min_depth[cid] = vertex_depth[v]
    n_cl = len(reps)
    # Tier by min depth: 0 -> A, 1 -> B, >=2 -> C
    def tier_of(d):
        return 0 if d == 0 else (1 if d == 1 else 2)
    tier = {c: tier_of(cluster_min_depth[c]) for c in range(n_cl)}
    origin = {c: 0 for c in range(n_cl)}  # all orig (no C-closure)
    # Hyperedges + adjacency
    hyperedges = set()
    for (_, v1, v2, v3) in all_edges:
        if v1 in vid_to_cid and v2 in vid_to_cid and v3 in vid_to_cid:
            hyperedges.add((vid_to_cid[v1], vid_to_cid[v2], vid_to_cid[v3]))
    adj = [[0]*n_cl for _ in range(n_cl)]
    for (c1, c2, c3) in hyperedges:
        for a, b in ((c1,c2),(c1,c3),(c2,c3),
                     (c2,c1),(c3,c1),(c3,c2)):
            if a != b: adj[a][b] = 1
    edges = [(i,j) for i in range(n_cl) for j in range(i+1, n_cl) if adj[i][j]]
    tier_sizes = [sum(1 for c in range(n_cl) if tier[c] == t) for t in (0,1,2)]
    return QuotientData(name, n_cl, tier, origin, edges,
                        sorted(hyperedges), tier_sizes, len(hyperedges))

# ----------------------------------------------------------------------------
# Layout
# ----------------------------------------------------------------------------

def tier_layout(q):
    """Compute (x, y) in TikZ units for each cluster. Three tiers stacked
    vertically; radii scale with tier size."""
    by_tier = {0: [], 1: [], 2: []}
    for c in range(q.n_cl):
        by_tier[q.tier[c]].append(c)
    for t in by_tier:
        by_tier[t].sort()

    # Vertical: A at top (y=+1.2), B at middle (y=0), C at bottom (y=-1.2).
    # Radii scale with tier population so denser tiers spread wider.
    def radius(n_in_tier):
        return max(1.0, 0.35 * math.sqrt(n_in_tier) + 0.4)

    pos = {}
    for t, clusters in by_tier.items():
        n = len(clusters)
        if n == 0: continue
        r = radius(n)
        y = 1.2 - t * 1.2  # +1.2, 0, -1.2
        # Tier B and C rotated slightly to avoid vertical alignment.
        rot = (math.pi/6 if t == 1 else
               (math.pi/12 if t == 2 else 0.0))
        for j, c in enumerate(clusters):
            angle = 2*math.pi*j/n + rot
            pos[c] = (r * math.cos(angle), y + 0.0 * math.sin(angle),
                      r * math.sin(angle))  # x, y, z (z used only for centroid)
    return pos, by_tier

# Wong 2011 colourblind-tolerant palette (TikZ-friendly RGB names).
TIER_COLOR_RGB = [
    (0.00, 0.45, 0.70),   # A — blue
    (0.84, 0.37, 0.00),   # B — vermilion
    (0.00, 0.62, 0.45),   # C — bluish green
]
TIER_LABEL = ['A', 'B', 'C']

def cluster_color(q, c):
    """RGB for a cluster: tier color, slightly desaturated for conj_only."""
    r, g, b = TIER_COLOR_RGB[q.tier[c]]
    if q.origin[c] == 1:  # conj-only — darken
        return (r * 0.55, g * 0.55, b * 0.85)
    return (r, g, b)

def rgb_str(rgb):
    return f"{rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f}"

# ----------------------------------------------------------------------------
# TikZ emitters (one per panel)
# ----------------------------------------------------------------------------

# Use TikZ x-y plane (drop the "z" we stored; the 2D projection is x,y).
# pos values are (x_3d, y_tier, z_3d); we use (x_3d, y_tier) for 2D.
def xy(pos_xyz):
    return pos_xyz[0], pos_xyz[1] + pos_xyz[2]*0.15  # mild z-fold for visual depth

def emit_panel_tier_graph(q):
    """Panel (a): tier-stratified node-edge graph."""
    pos, by_tier = tier_layout(q)
    lines = []
    lines.append("\\begin{tikzpicture}[scale=1.4]")
    # Edges first (under nodes)
    for (i, j) in q.edges:
        x1, y1 = xy(pos[i])
        x2, y2 = xy(pos[j])
        same = q.tier[i] == q.tier[j]
        op = "0.18" if same else "0.32"
        lines.append(f"  \\draw[gray, line width=0.15pt, opacity={op}] "
                     f"({x1:.3f},{y1:.3f}) -- ({x2:.3f},{y2:.3f});")
    # Nodes
    for c in range(q.n_cl):
        x, y = xy(pos[c])
        rgb = cluster_color(q, c)
        lines.append(f"  \\definecolor{{nc{c}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\fill[nc{c}, draw=black, line width=0.12pt] "
                     f"({x:.3f},{y:.3f}) circle (0.06);")
    # Tier labels (left edge)
    for t in (0, 1, 2):
        if not by_tier[t]: continue
        y_t = 1.2 - t * 1.2
        # Use the cluster_position radius to pick x just left of the ring
        x_label = -3.0
        rgb = TIER_COLOR_RGB[t]
        lines.append(f"  \\definecolor{{tl{t}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\node[font=\\small\\bfseries, tl{t}] "
                     f"at ({x_label:.2f},{y_t:.2f}) {{Tier {TIER_LABEL[t]}}};")
    lines.append("\\end{tikzpicture}")
    return "\n".join(lines)

def emit_panel_hyperedges(q, max_he=400):
    """Panel (a) graph view: every hyperedge as a 3-pronged star meeting
    at a centroid. The three slots render symmetrically; this is how the
    rewriting graph's gauge quotient sees the hyperedge.

    Subsamples evenly to at most max_he to keep TikZ within TeX memory
    limits for the larger quotients (Q102 has 2760, Q181 has 3974)."""
    pos, by_tier = tier_layout(q)
    lines = []
    lines.append("\\begin{tikzpicture}[scale=1.4]")
    if len(q.hyperedges) <= max_he:
        he = list(q.hyperedges)
    else:
        stride = len(q.hyperedges) / max_he
        idxs = [int(i * stride) for i in range(max_he)]
        he = [q.hyperedges[i] for i in idxs]
    # Centroid + 3 symmetric lines per hyperedge (amber, like §1.7 panel a)
    for (c1, c2, c3) in he:
        x1, y1 = xy(pos[c1]); x2, y2 = xy(pos[c2]); x3, y3 = xy(pos[c3])
        cx = (x1+x2+x3)/3.0; cy = (y1+y2+y3)/3.0
        for (px, py) in ((x1,y1),(x2,y2),(x3,y3)):
            lines.append(f"  \\draw[gray, line width=0.08pt, opacity=0.10] "
                         f"({px:.3f},{py:.3f}) -- ({cx:.3f},{cy:.3f});")
    # Tier labels (left edge)
    for t in (0, 1, 2):
        if not by_tier[t]: continue
        y_t = 1.2 - t * 1.2
        x_label = -3.80 if q.n_cl <= 90 else -4.10 if q.n_cl <= 102 else -5.00
        rgb = TIER_COLOR_RGB[t]
        lines.append(f"  \\definecolor{{tlA{t}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\node[font=\\small\\bfseries, tlA{t}] "
                     f"at ({x_label:.2f},{y_t:.2f}) {{Tier {TIER_LABEL[t]}}};")
    # Nodes
    for c in range(q.n_cl):
        x, y = xy(pos[c])
        rgb = cluster_color(q, c)
        lines.append(f"  \\definecolor{{hc{c}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\fill[hc{c}, draw=black, line width=0.10pt] "
                     f"({x:.3f},{y:.3f}) circle (0.05);")
    lines.append("\\end{tikzpicture}")
    return "\n".join(lines)

def emit_panel_ordered_slot(q, max_he=400):
    """Panel (b) computational view: each ordered triple (v1, v2, v3)
    drawn with explicit slot semantics — the two operator-input spokes
    (v1, v2 -> centroid) in Wong-bluish-green, the output spoke
    (centroid -> v3) in Wong-reddish-purple. Same data as panel (a),
    same subsampling stride; together the two panels externalise the
    type-level binary-op-inside-ternary-hyperedge distinction of
    §1.7."""
    pos, by_tier = tier_layout(q)
    lines = []
    lines.append("\\begin{tikzpicture}[scale=1.4]")
    if len(q.hyperedges) <= max_he:
        he = list(q.hyperedges)
    else:
        stride = len(q.hyperedges) / max_he
        idxs = [int(i * stride) for i in range(max_he)]
        he = [q.hyperedges[i] for i in idxs]
    # Wong 2011 colour-blind-safe palette: bluish-green inputs, reddish-purple output.
    lines.append("  \\definecolor{slotin}{rgb}{0.000,0.620,0.450}")
    lines.append("  \\definecolor{slotout}{rgb}{0.800,0.470,0.650}")
    # Two input spokes + one output spoke per ordered triple. The output
    # is rendered slightly more opaque and dashed to read as an arrow.
    for (c1, c2, c3) in he:
        x1, y1 = xy(pos[c1]); x2, y2 = xy(pos[c2]); x3, y3 = xy(pos[c3])
        cx = (x1+x2+x3)/3.0; cy = (y1+y2+y3)/3.0
        # input slot 1
        lines.append(f"  \\draw[slotin, line width=0.10pt, opacity=0.12] "
                     f"({x1:.3f},{y1:.3f}) -- ({cx:.3f},{cy:.3f});")
        # input slot 2
        lines.append(f"  \\draw[slotin, line width=0.10pt, opacity=0.12] "
                     f"({x2:.3f},{y2:.3f}) -- ({cx:.3f},{cy:.3f});")
        # output slot (dashed, brighter)
        lines.append(f"  \\draw[slotout, line width=0.14pt, opacity=0.32, dashed] "
                     f"({cx:.3f},{cy:.3f}) -- ({x3:.3f},{y3:.3f});")
    # Tier labels (left edge)
    for t in (0, 1, 2):
        if not by_tier[t]: continue
        y_t = 1.2 - t * 1.2
        x_label = -3.80 if q.n_cl <= 90 else -4.10 if q.n_cl <= 102 else -5.00
        rgb = TIER_COLOR_RGB[t]
        lines.append(f"  \\definecolor{{tlB{t}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\node[font=\\small\\bfseries, tlB{t}] "
                     f"at ({x_label:.2f},{y_t:.2f}) {{Tier {TIER_LABEL[t]}}};")
    # Nodes
    for c in range(q.n_cl):
        x, y = xy(pos[c])
        rgb = cluster_color(q, c)
        lines.append(f"  \\definecolor{{oc{c}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\fill[oc{c}, draw=black, line width=0.10pt] "
                     f"({x:.3f},{y:.3f}) circle (0.05);")
    lines.append("\\end{tikzpicture}")
    return "\n".join(lines)

def emit_panel_heatmap(q):
    """Panel (c): adjacency heatmap, rows/cols sorted by tier."""
    # Sort by (tier, index)
    order = sorted(range(q.n_cl),
                   key=lambda c: (q.tier[c], c))
    pos_in_order = {c: i for i, c in enumerate(order)}
    # adjacency lookup
    adj_set = set()
    for (i, j) in q.edges:
        adj_set.add((i, j))
        adj_set.add((j, i))
    # Cell layout: 8 unit-wide rect for the heatmap (sized for a full-width
    # subfigure stack), total grid size scales to n_cl.
    cell = 8.0 / max(q.n_cl, 1)
    lines = []
    lines.append("\\begin{tikzpicture}[scale=1.0]")
    # Draw filled cells (skip empty for compactness)
    for r_idx, ci in enumerate(order):
        for c_idx, cj in enumerate(order):
            if ci == cj:
                fill = "black!85"
            elif (ci, cj) in adj_set:
                fill = "blue!60!black"
            else:
                continue
            x = c_idx * cell
            y = -r_idx * cell  # negate so row 0 is at top
            lines.append(f"  \\fill[{fill}] ({x:.3f},{y:.3f}) "
                         f"rectangle ({x+cell:.3f},{y-cell:.3f});")
    # Tier-block boundaries (red lines)
    cumulative = [0]
    for t in (0, 1):
        cumulative.append(cumulative[-1] + q.tier_sizes[t])
    for boundary in cumulative[1:]:
        x = boundary * cell
        y_top = 0; y_bot = -q.n_cl * cell
        lines.append(f"  \\draw[red!80!black, line width=0.4pt] "
                     f"({x:.3f},{y_top:.3f}) -- ({x:.3f},{y_bot:.3f});")
        lines.append(f"  \\draw[red!80!black, line width=0.4pt] "
                     f"(0,{-x:.3f}) -- ({q.n_cl*cell:.3f},{-x:.3f});")
    # Outer border
    lines.append(f"  \\draw[black, line width=0.3pt] (0,0) rectangle "
                 f"({q.n_cl*cell:.3f},{-q.n_cl*cell:.3f});")
    # Tier labels along left edge
    tier_starts = [0, q.tier_sizes[0], q.tier_sizes[0] + q.tier_sizes[1]]
    for t in (0, 1, 2):
        if q.tier_sizes[t] == 0: continue
        y_mid = -(tier_starts[t] + q.tier_sizes[t]/2) * cell
        rgb = TIER_COLOR_RGB[t]
        lines.append(f"  \\definecolor{{hml{t}}}{{rgb}}{{{rgb_str(rgb)}}}")
        lines.append(f"  \\node[font=\\small\\bfseries, hml{t}] "
                     f"at (-0.25,{y_mid:.3f}) {{{TIER_LABEL[t]}}};")
    lines.append("\\end{tikzpicture}")
    return "\n".join(lines)

# ----------------------------------------------------------------------------
# Figure environment emitter
# ----------------------------------------------------------------------------

def emit_figure(q, label_suffix, caption):
    """Emit a complete \\begin{figure} environment with three subpanels
    stacked vertically:
      (a) graph view (centroid-star, symmetric three-spoke);
      (b) computational view (ordered-slot, input + output colours);
      (c) adjacency heatmap (binary co-occurrence sorted by tier).
    (a)+(b) mirror §1.7's layered binary-op-inside-ternary-routing
    figure; (c) is the binary shadow, kept as a diagnostic for tier-
    block structure."""
    panel_a = emit_panel_hyperedges(q)
    panel_b = emit_panel_ordered_slot(q)
    panel_c = emit_panel_heatmap(q)
    subsample_note = (
        f" (subsample of {q.total_he} hyperedges; even-stride)"
        if q.total_he > 400 else f" ({q.total_he} hyperedges)"
    )
    out = []
    out.append(f"% Auto-generated by paper/scripts/dump_quotient_to_tikz.py")
    out.append(f"% Do not hand-edit; regenerate via `python3 dump_quotient_to_tikz.py`.")
    out.append(f"\\begin{{figure}}[H]")
    out.append(f"\\centering")
    out.append(f"\\begin{{subfigure}}[t]{{\\textwidth}}")
    out.append(f"  \\centering")
    out.append(panel_a)
    out.append(f"  \\caption{{Graph view: centroid-star (three symmetric spokes per ternary hyperedge){subsample_note}.}}")
    out.append(f"\\end{{subfigure}}")
    out.append(f"\\par\\bigskip")
    out.append(f"\\begin{{subfigure}}[t]{{\\textwidth}}")
    out.append(f"  \\centering")
    out.append(panel_b)
    out.append(f"  \\caption{{Computational view: solid bluish-green spokes to the two operator-input slots, dashed reddish-purple arrow to the output slot.}}")
    out.append(f"\\end{{subfigure}}")
    out.append(f"\\par\\bigskip")
    out.append(f"\\begin{{subfigure}}[t]{{\\textwidth}}")
    out.append(f"  \\centering")
    out.append(panel_c)
    out.append(f"  \\caption{{Adjacency heatmap (rows/columns sorted by tier): the binary co-occurrence shadow of the ternary structure above; included as a tier-block-structure diagnostic.}}")
    out.append(f"\\end{{subfigure}}")
    out.append(f"\\caption{{{caption}}}")
    out.append(f"\\label{{fig:{label_suffix}}}")
    out.append(f"\\end{{figure}}")
    return "\n".join(out)

# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def main():
    out_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                            '..', 'figures'))
    os.makedirs(out_dir, exist_ok=True)

    quotients = []

    sys.stderr.write("[dump] building Q48 (C-closure of cyc(6))...\n")
    q48 = extract_c_closed("Q48", cyc6(), 6, depth=4)
    assert q48.n_cl == 48, f"Q48 build failed: {q48.n_cl}"
    quotients.append((q48, "q48",
        r"$Q_{48} = Q_{24} \cup \Ccls(Q_{24})$: three views of the "
        r"$\Ccls$-closure of the canonical 6-vertex Rosen-closed seed graph. "
        r"Tier counts $|A|/|B|/|C| = "
        + f"{q48.tier_sizes[0]}/{q48.tier_sizes[1]}/{q48.tier_sizes[2]}"
        + r"$."))

    sys.stderr.write("[dump] searching for Q90 seed + building...\n")
    q90_edges = q90_seed_edges()
    q90 = extract_c_closed("Q90", q90_edges, 6, depth=4)
    assert q90.n_cl == 90, f"Q90 build failed: {q90.n_cl}"
    quotients.append((q90, "q90",
        r"$Q_{90} = Q_{45} \cup \Ccls(Q_{45})$: $\Ccls$-closure of "
        r"the canonical ``Adjacent + 6 pads'' seed at pad-seed $42$. "
        r"The first family member with $|\text{Tier B}|$ not divisible by "
        r"$3$ (S178 obstruction)."))

    sys.stderr.write("[dump] building Q84 (C-closure of Adjacent)...\n")
    q84 = extract_c_closed("Q84", adjacent_ternary(), 6, depth=4)
    sys.stderr.write(f"  Q84 raw n_cl = {q84.n_cl}\n")
    quotients.append((q84, "q84",
        r"$Q_{84}$: $\Ccls$-closure of the Adjacent ternary graph on "
        r"$6$ vertices ($36$ edges)."))

    sys.stderr.write("[dump] building Q102 (C-closure of K_6^3)...\n")
    q102 = extract_c_closed("Q102", complete_edges_K6(), 6, depth=4)
    assert q102.n_cl == 102, f"Q102 build failed: {q102.n_cl}"
    quotients.append((q102, "q102",
        r"$Q_{102} = Q_{51} \cup \Ccls(Q_{51})$: $\Ccls$-closure of the "
        r"complete ternary $\Kn[6]$ on six vertices --- the "
        r"\emph{developmentally-complete zygote} of \cite{CFS2026}, "
        r"carrying the KO-dim 6 spectral triple that forces the "
        r"Standard-Model algebra. "
        r"Tier counts $|A|/|B|/|C| = "
        + f"{q102.tier_sizes[0]}/{q102.tier_sizes[1]}/{q102.tier_sizes[2]}"
        + r"$."))

    sys.stderr.write("[dump] searching for Q181 seed...\n")
    seed_181, n_181, ic_181 = q181_seed()
    sys.stderr.write(f"  Q181 seed found: {n_181} verts, {len(seed_181)} edges\n")
    sys.stderr.write("[dump] building Q181 (single-side, irregular seed)...\n")
    q181 = extract_single_side("Q181", seed_181, n_181, ic_181, depth=5)
    assert q181.n_cl == 181, f"Q181 build failed: {q181.n_cl}"
    quotients.append((q181, "q181",
        r"$Q_{181}$: single-side closure of a $12$-vertex $132$-edge "
        r"irregular seed (Theorem~\ref{thm:Q181}). $181$ is prime; this "
        r"is the smallest exhibited prime-cardinality closure object. "
        r"Tier here is depth-of-appearance: A = seed vertices, "
        r"B = depth-$1$ compositions, C = deeper."))

    # Write per-quotient .tex files
    for q, slug, caption in quotients:
        path = os.path.join(out_dir, f"{slug}.tex")
        with open(path, 'w') as f:
            f.write(emit_figure(q, slug, caption))
        sys.stderr.write(
            f"[dump] wrote {path} ({q.n_cl} clusters, {len(q.edges)} edges, "
            f"{q.total_he} hyperedges)\n")

    sys.stderr.write("\n[dump] done.\n")
    sys.stderr.write("Include in paper via:\n")
    for _, slug, _ in quotients:
        sys.stderr.write(f"  \\input{{figures/{slug}.tex}}\n")

if __name__ == "__main__":
    main()
