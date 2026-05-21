#!/usr/bin/env python3
"""
visualize_Q24.py - Three side-by-side visualizations of the closure
quotient Q24, written to a self-contained HTML file. Pure stdlib only:
no matplotlib, no networkx, no graphviz, no external SVG library.

What gets visualised
====================
Q24 is the closure quotient of the 6-vertex cyclic ternary seed
cyc(6) under truly-generic ICs. The multiway construction at depth 5
produces a finite vertex set; the gauge quotient by exact projective
equivalence (Definition 3.1 in the paper) collapses it to 24 classes.

The 24 clusters stratify by minimum depth-of-appearance into three
"tiers" matching the CFS paper's terminology:
  - Tier A (6): the original seed vertices, depth 0
  - Tier B (6): first compositions, depth 1
  - Tier C (12): deeper compositions, depths >= 2

Three figures
=============
  (1) Tier-stratified node-edge graph - 24 nodes in three rows,
      coloured by tier, with composition edges drawn between tiers.
      Cleanest at-a-glance picture.

  (2) Hyperedge diagram - same node layout, but every multiway
      hyperedge (s1, s2, s3) drawn as a 3-pronged star connecting
      its three cluster endpoints. Preserves the ternary arity that
      the binary collapse of (1) loses. Edges coloured by their
      construction depth.

  (3) Adjacency matrix heatmap - 24x24 binary adjacency derived from
      cluster co-occurrence in hyperedges, sorted by tier. Shows the
      tier-block structure (A-block, A<->B coupling, B<->C coupling).

Output
======
  paper/figures/Q24_visualization.html - single HTML file with all
  three SVG figures inlined and laid out side-by-side via flexbox.
  Opens in any browser. No JavaScript, no external assets.

Usage
=====
  python3 visualize_Q24.py
  open ../figures/Q24_visualization.html
"""

import os
import sys
import math

# Re-use the exact-arithmetic primitives from the existing builder.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Kn3_large_n_probe import (
    gmul, gsub, gconj, giszero, cross3, compose, proj_equiv, is_zero_vec,
    make_ics,  # uses _GEN (Pool A); fine for Q24 (Pool A is generic at n=6)
)


# ----------------------------------------------------------------------------
# DATA EXTRACTION
# ----------------------------------------------------------------------------

def cyc6():
    """The 6-vertex cyclic ternary seed. Closes to Q24 at depth >= 4."""
    n = 6
    return [(i, (i + 1) % n, (i + 2) % n) for i in range(n)]


def build_multiway_with_meta(topo, psi_init, depth):
    """Build multiway and return per-vertex depth + per-edge metadata
    so the visualization can extract tier and hyperedge structure."""
    psi = {k: list(v) for k, v in psi_init.items()}
    next_vid = max(psi_init) + 1
    vertex_depth = {v: 0 for v in psi_init}
    edges = [(0, s1, s2, s3) for (s1, s2, s3) in topo]
    cache = {}
    for d in range(depth):
        for (_, v1, v2, v3) in [e for e in edges if e[0] == d]:
            key = (v1, v2)
            if key not in cache:
                w = compose(psi[v1], psi[v2])
                if is_zero_vec(w):
                    continue
                psi[next_vid] = w
                cache[key] = next_vid
                vertex_depth[next_vid] = d + 1
                next_vid += 1
            w = cache[key]
            edges.append((d + 1, w, v2, v3))
            edges.append((d + 1, w, v1, v3))
            edges.append((d + 1, w, v1, v2))
    return psi, vertex_depth, edges


def cluster_with_tier(psi, vertex_depth):
    """Bucket multiway vertices into clusters by exact projective
    equivalence; tag each cluster with min-depth-of-appearance."""
    reps = []                # list of (cluster_id, rep_vector, min_depth)
    vid_to_cid = {}
    for v in sorted(psi):
        pv = psi[v]
        if is_zero_vec(pv):
            continue
        matched = -1
        for i, (cid, r, _) in enumerate(reps):
            if proj_equiv(pv, r):
                matched = i
                break
        if matched >= 0:
            cid = reps[matched][0]
            vid_to_cid[v] = cid
            # update min_depth if this vertex appeared earlier
            existing = reps[matched]
            if vertex_depth[v] < existing[2]:
                reps[matched] = (existing[0], existing[1], vertex_depth[v])
        else:
            cid = len(reps)
            reps.append((cid, pv, vertex_depth[v]))
            vid_to_cid[v] = cid
    cluster_depth = {c[0]: c[2] for c in reps}
    return vid_to_cid, cluster_depth, len(reps)


def cluster_hyperedges(edges, vid_to_cid):
    """Map each multiway hyperedge to a triple of cluster IDs.
    Returns the deduplicated set of (c1, c2, c3) triples."""
    he = set()
    for (d, v1, v2, v3) in edges:
        if v1 in vid_to_cid and v2 in vid_to_cid and v3 in vid_to_cid:
            c1, c2, c3 = vid_to_cid[v1], vid_to_cid[v2], vid_to_cid[v3]
            he.add((c1, c2, c3))
    return sorted(he)


def cluster_adjacency(hyperedges, n_cl):
    """Build a 24x24 binary adjacency matrix from hyperedge co-occurrence:
    A[i][j] = 1 iff clusters i, j appear together in any hyperedge."""
    adj = [[0] * n_cl for _ in range(n_cl)]
    for (c1, c2, c3) in hyperedges:
        for a, b in ((c1, c2), (c1, c3), (c2, c3),
                     (c2, c1), (c3, c1), (c3, c2)):
            adj[a][b] = 1
    for i in range(n_cl):
        adj[i][i] = 0
    return adj


# ----------------------------------------------------------------------------
# SVG GENERATION
# ----------------------------------------------------------------------------

# Tier color palette - colourblind-tolerant (Wong 2011)
TIER_COLOUR = {
    'A': '#0072B2',   # blue
    'B': '#D55E00',   # vermilion
    'C': '#009E73',   # bluish green
}

def tier_of(depth):
    if depth == 0:
        return 'A'
    elif depth == 1:
        return 'B'
    else:
        return 'C'


def svg_open(w, h, title):
    """Header for a single SVG figure."""
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}" font-family="sans-serif" font-size="11">\n'
        f'  <title>{title}</title>\n'
    )

def svg_close():
    return '</svg>\n'


def tier_layout(cluster_depth, n_cl, width, height, top=60, bottom=40,
                tier_pad=20):
    """Compute (x, y) for each cluster in a three-row tier layout."""
    by_tier = {'A': [], 'B': [], 'C': []}
    for c in range(n_cl):
        by_tier[tier_of(cluster_depth[c])].append(c)
    for t in by_tier:
        by_tier[t].sort()
    # Three rows, evenly spaced vertically
    inner_h = height - top - bottom
    row_y = {
        'A': top + 0 * inner_h / 2,
        'B': top + 1 * inner_h / 2,
        'C': top + 2 * inner_h / 2,
    }
    pos = {}
    for t, clusters in by_tier.items():
        k = len(clusters)
        if k == 0:
            continue
        # Spread across [tier_pad, width - tier_pad]
        avail = width - 2 * tier_pad
        for j, c in enumerate(clusters):
            x = tier_pad + (j + 0.5) * avail / k
            pos[c] = (x, row_y[t])
    return pos, by_tier


# ---- FIGURE 1: tier-stratified node-edge graph ------------------------------

def figure_tier_graph(cluster_depth, adj, n_cl,
                       width=520, height=420):
    s = svg_open(width, height, "Q24 - tier-stratified node-edge graph")
    s += (f'  <text x="{width/2}" y="22" text-anchor="middle" '
          f'font-weight="bold" font-size="13">'
          f'Fig. 1. Tier-stratified node-edge graph</text>\n')
    pos, by_tier = tier_layout(cluster_depth, n_cl, width, height)
    # Edges first (drawn under nodes)
    drawn_edges = set()
    for i in range(n_cl):
        for j in range(i + 1, n_cl):
            if adj[i][j]:
                e_key = (i, j)
                if e_key in drawn_edges:
                    continue
                drawn_edges.add(e_key)
                xi, yi = pos[i]
                xj, yj = pos[j]
                # Stroke darker for cross-tier, lighter for intra-tier
                same = tier_of(cluster_depth[i]) == tier_of(cluster_depth[j])
                colour = "#aaaaaa" if same else "#444444"
                width_e = "0.6" if same else "0.9"
                s += (f'  <line x1="{xi:.1f}" y1="{yi:.1f}" '
                      f'x2="{xj:.1f}" y2="{yj:.1f}" stroke="{colour}" '
                      f'stroke-width="{width_e}" stroke-opacity="0.5"/>\n')
    # Nodes on top
    for c in range(n_cl):
        x, y = pos[c]
        t = tier_of(cluster_depth[c])
        fill = TIER_COLOUR[t]
        s += (f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="9" '
              f'fill="{fill}" stroke="#222" stroke-width="0.7"/>\n')
        s += (f'  <text x="{x:.1f}" y="{y+3.5:.1f}" text-anchor="middle" '
              f'fill="white" font-size="9" font-weight="bold">{c}</text>\n')
    # Tier labels at left
    for t in ('A', 'B', 'C'):
        if by_tier[t]:
            _, ry = pos[by_tier[t][0]]
            s += (f'  <text x="14" y="{ry+4:.1f}" font-size="13" '
                  f'font-weight="bold" fill="{TIER_COLOUR[t]}">'
                  f'Tier {t}</text>\n')
    s += svg_close()
    return s


# ---- FIGURE 2: hyperedge diagram --------------------------------------------

def figure_hyperedges(cluster_depth, hyperedges, n_cl,
                      width=520, height=420, max_he=None):
    s = svg_open(width, height, "Q24 - hyperedge structure")
    n_he_total = len(hyperedges)
    s += (f'  <text x="{width/2}" y="22" text-anchor="middle" '
          f'font-weight="bold" font-size="13">'
          f'Fig. 2. Hyperedge structure ({n_he_total} cluster hyperedges)'
          f'</text>\n')
    pos, by_tier = tier_layout(cluster_depth, n_cl, width, height)
    he = hyperedges if max_he is None else hyperedges[:max_he]
    # Each hyperedge: 3-pronged star meeting at centroid
    for (c1, c2, c3) in he:
        x1, y1 = pos[c1]; x2, y2 = pos[c2]; x3, y3 = pos[c3]
        cx = (x1 + x2 + x3) / 3
        cy = (y1 + y2 + y3) / 3
        for (px, py) in ((x1, y1), (x2, y2), (x3, y3)):
            s += (f'  <line x1="{px:.1f}" y1="{py:.1f}" '
                  f'x2="{cx:.1f}" y2="{cy:.1f}" '
                  f'stroke="#888" stroke-width="0.35" '
                  f'stroke-opacity="0.25"/>\n')
        # Tiny centroid dot
        s += (f'  <circle cx="{cx:.1f}" cy="{cy:.1f}" r="0.8" '
              f'fill="#666" fill-opacity="0.4"/>\n')
    # Nodes on top (smaller than fig 1 since edges are dense)
    for c in range(n_cl):
        x, y = pos[c]
        t = tier_of(cluster_depth[c])
        fill = TIER_COLOUR[t]
        s += (f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="7" '
              f'fill="{fill}" stroke="#222" stroke-width="0.7"/>\n')
        s += (f'  <text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" '
              f'fill="white" font-size="8" font-weight="bold">{c}</text>\n')
    # Tier labels at left
    for t in ('A', 'B', 'C'):
        if by_tier[t]:
            _, ry = pos[by_tier[t][0]]
            s += (f'  <text x="14" y="{ry+4:.1f}" font-size="13" '
                  f'font-weight="bold" fill="{TIER_COLOUR[t]}">'
                  f'Tier {t}</text>\n')
    s += svg_close()
    return s


# ---- FIGURE 3: adjacency matrix heatmap ------------------------------------

def figure_heatmap(cluster_depth, adj, n_cl, width=480, height=480):
    s = svg_open(width, height, "Q24 - adjacency heatmap (tier-sorted)")
    s += (f'  <text x="{width/2}" y="22" text-anchor="middle" '
          f'font-weight="bold" font-size="13">'
          f'Fig. 3. Adjacency heatmap (rows/cols sorted by tier)</text>\n')
    # Sort clusters by (tier, original-id) for nice block structure
    tier_rank = {'A': 0, 'B': 1, 'C': 2}
    order = sorted(range(n_cl),
                   key=lambda c: (tier_rank[tier_of(cluster_depth[c])], c))
    # Layout: grid centered, with axis labels
    margin_top = 50
    margin_left = 50
    margin_right = 20
    margin_bottom = 30
    grid_w = width - margin_left - margin_right
    grid_h = height - margin_top - margin_bottom
    cell = min(grid_w, grid_h) / n_cl
    grid_size = cell * n_cl
    gx0 = margin_left + (grid_w - grid_size) / 2
    gy0 = margin_top + (grid_h - grid_size) / 2
    # Cells
    for r, ci in enumerate(order):
        for c, cj in enumerate(order):
            x = gx0 + c * cell
            y = gy0 + r * cell
            if r == c:
                fill = "#222222"
            elif adj[ci][cj]:
                fill = "#1f3a5f"
            else:
                fill = "#f0f0f0"
            s += (f'  <rect x="{x:.1f}" y="{y:.1f}" '
                  f'width="{cell:.1f}" height="{cell:.1f}" '
                  f'fill="{fill}" stroke="white" stroke-width="0.3"/>\n')
    # Tier dividers
    tier_counts = {'A': 0, 'B': 0, 'C': 0}
    for c in range(n_cl):
        tier_counts[tier_of(cluster_depth[c])] += 1
    nA, nB = tier_counts['A'], tier_counts['B']
    # Horizontal dividers (between tiers)
    for boundary in (nA, nA + nB):
        y = gy0 + boundary * cell
        s += (f'  <line x1="{gx0:.1f}" y1="{y:.1f}" '
              f'x2="{gx0 + grid_size:.1f}" y2="{y:.1f}" '
              f'stroke="#cc0000" stroke-width="1.5"/>\n')
        x = gx0 + boundary * cell
        s += (f'  <line x1="{x:.1f}" y1="{gy0:.1f}" '
              f'x2="{x:.1f}" y2="{gy0 + grid_size:.1f}" '
              f'stroke="#cc0000" stroke-width="1.5"/>\n')
    # Tier labels along left
    tier_starts = {'A': 0, 'B': nA, 'C': nA + nB}
    for t, start in tier_starts.items():
        n_t = tier_counts[t]
        y_mid = gy0 + (start + n_t / 2) * cell
        s += (f'  <text x="{gx0 - 8:.1f}" y="{y_mid + 4:.1f}" '
              f'text-anchor="end" font-size="11" font-weight="bold" '
              f'fill="{TIER_COLOUR[t]}">{t}</text>\n')
        x_mid = gx0 + (start + n_t / 2) * cell
        s += (f'  <text x="{x_mid:.1f}" y="{gy0 - 6:.1f}" '
              f'text-anchor="middle" font-size="11" font-weight="bold" '
              f'fill="{TIER_COLOUR[t]}">{t}</text>\n')
    # Legend
    legy = height - 14
    s += (f'  <rect x="{margin_left}" y="{legy - 8}" width="10" height="10" '
          f'fill="#1f3a5f"/>\n')
    s += (f'  <text x="{margin_left + 14}" y="{legy}" font-size="10">'
          f'adjacent</text>\n')
    s += (f'  <rect x="{margin_left + 80}" y="{legy - 8}" width="10" height="10" '
          f'fill="#f0f0f0" stroke="#aaa" stroke-width="0.5"/>\n')
    s += (f'  <text x="{margin_left + 94}" y="{legy}" font-size="10">'
          f'not adjacent</text>\n')
    s += (f'  <rect x="{margin_left + 180}" y="{legy - 8}" width="10" height="10" '
          f'fill="#222222"/>\n')
    s += (f'  <text x="{margin_left + 194}" y="{legy}" font-size="10">'
          f'diagonal (self)</text>\n')
    s += svg_close()
    return s


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

def main():
    print("Building Q24 from cyc(6) at depth 5 (exact Gaussian-integer arithmetic)...")
    topo = cyc6()
    ics = make_ics(6)   # Pool A; fully generic at n=6 (Theorem 9.2 in the paper)
    psi, vertex_depth, edges = build_multiway_with_meta(topo, ics, depth=5)
    vid_to_cid, cluster_depth, n_cl = cluster_with_tier(psi, vertex_depth)
    print(f"  multiway vertices: {len(psi)}")
    print(f"  clusters (Q24):    {n_cl} (expected: 24)")
    assert n_cl == 24, f"Q24 build failed: got {n_cl} clusters"

    he = cluster_hyperedges(edges, vid_to_cid)
    adj = cluster_adjacency(he, n_cl)
    n_edges = sum(adj[i][j] for i in range(n_cl) for j in range(i+1, n_cl))
    print(f"  cluster hyperedges: {len(he)}")
    print(f"  binary cluster-cluster adjacencies: {n_edges}")

    tier_counts = {'A': 0, 'B': 0, 'C': 0}
    for c in range(n_cl):
        tier_counts[tier_of(cluster_depth[c])] += 1
    print(f"  tier counts: A={tier_counts['A']}, B={tier_counts['B']}, "
          f"C={tier_counts['C']} (expected: 6/6/12)")

    print("\nGenerating three SVG figures...")
    fig1 = figure_tier_graph(cluster_depth, adj, n_cl)
    fig2 = figure_hyperedges(cluster_depth, he, n_cl)
    fig3 = figure_heatmap(cluster_depth, adj, n_cl)

    # Wrap in a self-contained HTML file with side-by-side flexbox layout.
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Q24 visualization - three views</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Helvetica Neue',
                       Arial, sans-serif;
         max-width: 1700px; margin: 24px auto; padding: 0 16px;
         color: #222; background: #fafafa; }}
  h1   {{ margin: 0 0 8px 0; font-size: 22px; font-weight: 600; }}
  p    {{ margin: 4px 0 16px 0; color: #555; font-size: 14px; }}
  .row {{ display: flex; flex-wrap: wrap; gap: 16px;
          justify-content: center; align-items: flex-start; }}
  .panel {{ background: white; border: 1px solid #ddd;
            border-radius: 6px; padding: 8px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.05); }}
  .meta {{ margin-top: 18px; color: #666; font-size: 13px;
           border-top: 1px solid #ddd; padding-top: 12px; }}
  code {{ background: #f0f0f0; padding: 1px 4px; border-radius: 3px;
          font-size: 12.5px; }}
</style>
</head>
<body>
<h1>Q24 - three visualizations of the closure quotient</h1>
<p>Built from <code>cyc(6)</code> (the 6-vertex cyclic ternary seed) at
depth 5, under truly-generic ICs (Pool A; fully generic at n=6 per
Theorem 9.2). Exact Gaussian-integer ray-equivalence clustering -
no Float64 threshold. Tier counts: Tier A = 6, Tier B = 6, Tier C = 12
(total 24).</p>

<div class="row">
  <div class="panel">{fig1}</div>
  <div class="panel">{fig2}</div>
  <div class="panel">{fig3}</div>
</div>

<div class="meta">
<strong>Reading the figures.</strong>
<ul>
  <li><strong>Fig. 1</strong> - tier-stratified node-edge graph.
      24 clusters in three rows by tier; binary edges drawn between
      every cluster pair that co-occurs in some hyperedge. Cross-tier
      edges are darker than within-tier ones.</li>
  <li><strong>Fig. 2</strong> - hyperedge diagram. Same nodes; every
      ternary hyperedge drawn as a 3-pronged star meeting at its
      centroid. Preserves the ternary arity the binary collapse loses.</li>
  <li><strong>Fig. 3</strong> - adjacency heatmap. Rows and columns
      sorted by tier (A, then B, then C); red lines mark tier
      boundaries. Visualises the tier-block coupling structure.</li>
</ul>
<p>Generated by <code>paper/scripts/visualize_Q24.py</code>; pure
Python stdlib, no external dependencies. Reproducible end-to-end
from the canonical exact-arithmetic builder in
<code>Kn3_large_n_probe.py</code>.</p>
</div>
</body>
</html>
"""

    out_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                            '..', 'figures'))
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'Q24_visualization.html')
    with open(out_path, 'w') as f:
        f.write(html)
    print(f"\nWrote {out_path}")
    print(f"Open with: open '{out_path}'")
    return out_path


if __name__ == "__main__":
    main()
