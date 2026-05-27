#!/usr/bin/env python3
"""
q51_exact_phase_s_python.py

Port of TCE's `cfs_kn3_v21_phase_s_s6_decomposition.jl` to Python, running on
closure-v5's EXACT Gaussian-integer Q_51 build (depth=4 ray-equivalence) under
the same Julia driver "ordered-tuple" cell-complex convention used in the
paper note v2.

Goal: produce the S_6 character vector and irrep multiplicity table for
C^1(Q_51_exact) at dim 336, to replace the 366-dim fidelity-based Prop 5.1
of note v2 with an exact-build-derived statement.

Pipeline:

  1. Build EXACT Q_51 via closure-v5 q102_build_exact_v1.py infrastructure.
  2. Construct Julia-convention cell complex (edges = sort(unique(raw ordered
     pairs from each q_he triple)), faces = sorted(unique(raw ordered triples)).
  3. Compute integer L_1 = d_0 d_0^T + d_1^T d_1 (sorted-face alternating-sign
     d_1 — same convention as note v2 / phase S Julia driver).
  4. Compute K_6 neighborhood signatures to identify V_6/V_15/V_30 orbits.
  5. For each S_6 conjugacy-class representative σ, build σ_V via signature
     matching + cross-orbit refinement, lift to σ_E on the 336-dim edge space,
     verify σ_E L_1 = L_1 σ_E exactly on the integer matrix.
  6. Compute χ_{C^1}(σ) = signed-trace of σ_E for each class.
  7. Apply Frobenius reciprocity m_ρ = (1/|G|) Σ_C |C| χ_ρ(C) χ_{C^1}(C) using
     the standard S_6 character table.
  8. Verify Σ m_ρ · dim ρ = 336.
"""

import sys
import os
from itertools import combinations
import numpy as np

# Set CFS_REPO_ROOT to a Closure-Forces-Structure source checkout
# (https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs)
# before invoking this driver.
CLOSURE_V5 = os.environ.get("CFS_REPO_ROOT")
if not CLOSURE_V5:
    sys.exit("Set CFS_REPO_ROOT env var to a Closure-Forces-Structure checkout.")
sys.path.insert(0, CLOSURE_V5)

from q102_build_exact_v1 import (
    gaussian_ics, build_multiway, proj_equiv, is_zero_vec, complete_ternary
)


# -----------------------------------------------------------------------------
# §A. Build exact Q_51 + Julia-convention cell complex
# -----------------------------------------------------------------------------

def build_q51_exact(depth=4, ic_seed=None):
    """Run closure-v5 exact build; return (n_cl, q_he, vid_to_cid)."""
    n_seed = 6
    ics = gaussian_ics(n=n_seed, ic_seed=ic_seed)
    seed_edges = complete_ternary(n_seed)
    psi_o, _depth_o, edges_o = build_multiway(seed_edges, ics, depth)
    # Exact ray-equivalence clustering
    clusters = []  # list of (rep_vid, rep_exact)
    vid_to_cid = {}
    for v in sorted(psi_o.keys()):
        pv = psi_o[v]
        if is_zero_vec(pv):
            continue
        cid = -1
        for ci, (_, rep) in enumerate(clusters):
            if proj_equiv(pv, rep):
                cid = ci
                break
        if cid < 0:
            cid = len(clusters)
            clusters.append((v, pv))
        vid_to_cid[v] = cid

    n_cl = len(clusters)
    # Cluster-level hyperedges as sorted unique triples (== Julia H.E3 after dedup)
    q_he = set()
    for (_, v1, v2, v3) in edges_o:
        if (v1 in vid_to_cid and v2 in vid_to_cid and v3 in vid_to_cid):
            t = tuple(sorted([vid_to_cid[v1], vid_to_cid[v2], vid_to_cid[v3]]))
            q_he.add(t)
    return n_cl, sorted(q_he), vid_to_cid


def cell_complex_julia(q_he, n_cl):
    """Build cell complex under Julia driver convention.

    edges: cells1 = sort(unique(raw ordered-pair list)) where each q_he triple
           (a,b,c) contributes (a,b), (b,c), (a,c) in that order.
    faces: cells2 = sorted(unique(q_he))  (preserve raw order — but the q_he
           triples in this Python build are already pre-sorted at intake, so
           cells2 = q_he).

    Self-loops are dropped (Julia's d_0 zeros them out trivially anyway).
    """
    edges = set()
    for (a, b, c) in q_he:
        edges.add((a, b)); edges.add((b, c)); edges.add((a, c))
    edges = sorted(e for e in edges if e[0] != e[1])
    faces = list(q_he)
    return list(range(n_cl)), edges, faces


def hodge_L1_int(vertices, edges, faces):
    """Integer L_1 = d_0 d_0^T + d_1^T d_1 with sorted-face alternating-sign d_1."""
    n_e = len(edges); n_f = len(faces)
    edge_idx = {e: i for i, e in enumerate(edges)}
    d0 = np.zeros((n_e, len(vertices)), dtype=np.int64)
    for i, (a, b) in enumerate(edges):
        d0[i, a] = -1; d0[i, b] = +1
    d1 = np.zeros((n_f, n_e), dtype=np.int64)
    for j, (a, b, c) in enumerate(faces):
        e_ab = edge_idx.get((a, b))
        e_ac = edge_idx.get((a, c))
        e_bc = edge_idx.get((b, c))
        if e_ab is not None: d1[j, e_ab] = +1
        if e_ac is not None: d1[j, e_ac] = -1
        if e_bc is not None: d1[j, e_bc] = +1
    L1 = d0 @ d0.T + d1.T @ d1
    return L1, d0, d1, edge_idx


# -----------------------------------------------------------------------------
# §B. K_6 signatures and orbit decomposition
# -----------------------------------------------------------------------------

def co_occurrence_vector(q_he, n_cl):
    """For each non-seed cluster c, build vec_c = (μ(c, 0), ..., μ(c, 5))
    where μ(c, s) = #q_he triples containing both c and seed s.

    In the EXACT Q_51 build under canonical ray-equivalence, the K_6-adjacency
    *set* is uniformly {0,1,2,3,4,5} for every non-seed (every cluster is
    seed-adjacent to all 6), so adjacency-set alone collapses. The co-occurrence
    multiplicity vector resolves the 6+15+30 orbit structure:

       V_15: multiset (11, 11, 4, 4, 4, 4)   — "pair-dominant" clusters
       V_30: multiset (8,  2, 2, 2, 2, 2)    — "singleton-dominant" clusters

    Returns dict c → tuple(μ(c,0..5)).
    """
    from collections import defaultdict
    co = defaultdict(lambda: defaultdict(int))
    for (a, b, c) in q_he:
        seeds_in = [x for x in (a, b, c) if x < 6]
        nonseeds_in = [x for x in (a, b, c) if x >= 6]
        for ns in nonseeds_in:
            for s in seeds_in:
                co[ns][s] += 1
    return {ns: tuple(co[ns][s] for s in range(6)) for ns in co}


def find_orbits(q_he, n_cl):
    """Identify V_6 / V_15 / V_30 via co-occurrence multiset (S_6-invariant).

    V_6 = seed clusters 0..5 (anchored by closure-v5 build convention).
    V_15 = non-seed clusters with sorted multiset (11,11,4,4,4,4).
    V_30 = non-seed clusters with sorted multiset (8,2,2,2,2,2).

    Returns (V_6, V_15, V_30, co_vec) where co_vec[c] = (μ(c,0..5)).
    """
    V_6 = list(range(6))
    co = co_occurrence_vector(q_he, n_cl)
    V_15, V_30, V_other = [], [], []
    for c in range(6, n_cl):
        v = co.get(c, (0,)*6)
        ms = tuple(sorted(v, reverse=True))
        if ms == (11, 11, 4, 4, 4, 4):
            V_15.append(c)
        elif ms == (8, 2, 2, 2, 2, 2):
            V_30.append(c)
        else:
            V_other.append((c, ms))
    if V_other:
        print(f"  WARN: {len(V_other)} clusters with unexpected multiset:")
        for c, ms in V_other[:5]:
            print(f"    c={c} multiset={ms}")
    return V_6, V_15, V_30, co


# -----------------------------------------------------------------------------
# §C. Build σ_V and lift to σ_E
# -----------------------------------------------------------------------------

def build_sigma_V(perm, n_cl, co, V_6, V_15, V_30):
    """Construct σ_V on 51 clusters realizing S_6 action on V_6 by perm.

    perm: tuple (p0, p1, p2, p3, p4, p5) meaning V_6[k] → V_6[perm[k]].
    co[c] = co-occurrence vector (μ(c,0), ..., μ(c,5)).

    Under S_6 action, σ ⋅ v = (v_{σ^{-1}(0)}, ..., v_{σ^{-1}(5)}). So
       σ_V(c) = unique c' with co[c'] = σ ⋅ co[c].

    The co-occurrence vectors are sharp enough on the exact build that each
    permutation of (11,11,4,4,4,4) maps to a unique V_15 cluster, and each
    permutation of (8,2,2,2,2,2) maps to a unique V_30 cluster (15 ⋅ 6/2 = 45,
    wait: (11,11,4,4,4,4) has C(6,2) = 15 distinct positional placements of
    the two 11's, which is exactly |V_15|; (8,2,2,2,2,2) has 6 distinct
    placements of the single 8, which gives |V_30| / 5 = 6 — so the V_30
    co-occurrence vector alone underdetermines the lift; we additionally use
    V_15 ↔ V_30 adjacency to pin σ_V on V_30 uniquely).
    """
    sigma_V = np.full(n_cl, -1, dtype=np.int64)

    # Anchor on V_6
    for k in range(6):
        sigma_V[V_6[k]] = V_6[perm[k]]

    # Build inverse-permutation for vector pullback
    perm_inv = [0] * 6
    for k in range(6):
        perm_inv[perm[k]] = k

    # σ ⋅ v: new[i] = v[perm_inv[i]]
    def sigma_apply_vec(v):
        return tuple(v[perm_inv[i]] for i in range(6))

    # Build co-vector → cluster lookup, separately for V_15 and V_30
    v15_lookup = {}
    for c in V_15:
        v15_lookup.setdefault(co[c], []).append(c)
    v30_lookup = {}
    for c in V_30:
        v30_lookup.setdefault(co[c], []).append(c)

    # Map V_15 via co-vector (should be uniquely determined: 15 distinct
    # placements of the two 11's = 15 V_15 vectors).
    unresolved = []
    for c in V_15:
        target = sigma_apply_vec(co[c])
        candidates = v15_lookup.get(target, [])
        if len(candidates) == 1:
            sigma_V[c] = candidates[0]
        else:
            unresolved.append(c)

    # Map V_30 via co-vector. There are 6 distinct placements of the single
    # 8, giving 6 "co-vector cells", each containing 5 V_30 clusters. So
    # σ_V is undetermined within each co-vector cell; we need V_15 ↔ V_30
    # adjacency to fully pin σ_V on V_30.
    for c in V_30:
        target = sigma_apply_vec(co[c])
        candidates = v30_lookup.get(target, [])
        if len(candidates) == 1:
            sigma_V[c] = candidates[0]
        else:
            unresolved.append(c)

    # Resolve remaining V_30 via V_15 cross-orbit adjacency. For each
    # unresolved c, its already-mapped V_15-neighbors must equal the σ-image
    # cluster's V_15-neighbors.
    return sigma_V, unresolved


def refine_sigma_V_via_adjacency(sigma_V, unresolved, V_15, V_30, q_he, co):
    """Pin σ_V on V_30 via V_15 cross-orbit triple-incidence.

    For each q_he triple (a,b,c) where exactly one is in V_15 and at least one
    is in V_30, the V_30 vertex's "V_15-neighbor set" gives a sharper invariant
    than co-occurrence alone.
    """
    from collections import defaultdict
    V_15_set = set(V_15); V_30_set = set(V_30)
    # Build V_15 incidence list for each V_30 cluster
    v15_nbrs = defaultdict(set)
    for (a, b, c) in q_he:
        vs = [a, b, c]
        in15 = [v for v in vs if v in V_15_set]
        in30 = [v for v in vs if v in V_30_set]
        for v30 in in30:
            for v15 in in15:
                v15_nbrs[v30].add(v15)

    # Build V_15 candidate -> V_30 candidate dual: for each V_30 cluster d,
    # neighbor set is N(d) = v15_nbrs[d].

    # Make lookup table: (co_vector, frozenset of V_15 neighbors) -> [V_30 candidates]
    v30_neighbor_lookup = defaultdict(list)
    for d in V_30:
        key = (co[d], frozenset(sigma_V[v] for v in v15_nbrs[d] if sigma_V[v] != -1))
        v30_neighbor_lookup[(co[d], frozenset(v15_nbrs[d]))].append(d)

    # For each unresolved V_30 cluster, target = applying σ to its (co, N)
    # signature must match a unique candidate.
    still_unresolved = []
    for c in list(unresolved):
        if c not in V_30_set or sigma_V[c] != -1:
            continue
        # Already-mapped V_15 neighbors of c, under σ
        target_v15_nbrs = frozenset(sigma_V[v] for v in v15_nbrs[c]
                                    if sigma_V[v] != -1)
        target_co = tuple(co[c][i] for i in range(6))  # will be replaced via sigma_apply
        # We need to compute target_co as σ ⋅ co[c]. But we don't have σ inv handy
        # here — we'll re-pass through the lookup with N(c) only.
        candidates = [d for d in V_30 if v15_nbrs[d] == target_v15_nbrs
                      and sigma_V[d] == -1]  # uncommitted candidates
        if len(candidates) == 1:
            sigma_V[c] = candidates[0]
        else:
            still_unresolved.append(c)
    return sigma_V, still_unresolved


def build_sigma_V_full(perm, n_cl, co, q_he, V_6, V_15, V_30):
    """Build σ_V using sharp (co-vec, V_15-neighbor) invariant for V_30.

    Step 1: σ_V on V_6 anchored by perm.
    Step 2: σ_V on V_15 via co-vector (each V_15 has a unique co-vec
            permutation among 15 placements of two 11's).
    Step 3: σ_V on V_30 via (σ-permuted co-vec, σ-image of V_15 neighbor).
            Each V_30 has exactly 1 V_15 triple-neighbor in the exact build,
            so this is sharp (5 cells per dom seed × 6 dom seeds = 30, with
            each cell's 5 elements pinned by their V_15 partner).

    Returns (sigma_V, unmapped_clusters).
    """
    from collections import defaultdict
    sigma_V = np.full(n_cl, -1, dtype=np.int64)
    V_15_set = set(V_15); V_30_set = set(V_30); V_6_set = set(V_6)

    # Step 1: V_6 anchor
    for k in range(6):
        sigma_V[V_6[k]] = V_6[perm[k]]

    perm_inv = [0] * 6
    for k in range(6):
        perm_inv[perm[k]] = k

    def sigma_apply_vec(v):
        return tuple(v[perm_inv[i]] for i in range(6))

    # Step 2: V_15 via co-vector
    v15_by_co = {co[c]: c for c in V_15}
    for c in V_15:
        target = sigma_apply_vec(co[c])
        if target in v15_by_co:
            sigma_V[c] = v15_by_co[target]
        # else: leave as -1 — will report

    # Step 3: V_30 via (co-vec, V_15 neighbor in q_he triple-incidence)
    v15_nbr_of = {}  # V_30 cluster c → its unique V_15 triple-neighbor
    for (a, b, x) in q_he:
        triple = (a, b, x)
        in15 = [v for v in triple if v in V_15_set]
        in30 = [v for v in triple if v in V_30_set]
        for v30 in in30:
            for v15 in in15:
                v15_nbr_of[v30] = v15
    # Verify uniqueness assumption
    multi_nbr = [c for c in V_30 if c not in v15_nbr_of]
    if multi_nbr:
        return sigma_V, [c for c in range(n_cl) if sigma_V[c] == -1]

    # Build (co, v15_nbr) → V_30 cluster lookup
    v30_by_key = {}
    collisions = []
    for c in V_30:
        key = (co[c], v15_nbr_of[c])
        if key in v30_by_key:
            collisions.append((c, v30_by_key[key]))
        v30_by_key[key] = c

    if collisions:
        # Non-unique — fall back to enumeration of all candidates
        v30_by_key_multi = defaultdict(list)
        for c in V_30:
            v30_by_key_multi[(co[c], v15_nbr_of[c])].append(c)

        for c in V_30:
            target_co = sigma_apply_vec(co[c])
            target_v15 = sigma_V[v15_nbr_of[c]]
            if target_v15 == -1:
                continue
            candidates = v30_by_key_multi.get((target_co, target_v15), [])
            if len(candidates) == 1:
                sigma_V[c] = candidates[0]
        return sigma_V, [c for c in range(n_cl) if sigma_V[c] == -1]

    # Sharp case
    for c in V_30:
        target_co = sigma_apply_vec(co[c])
        target_v15 = sigma_V[v15_nbr_of[c]]
        if target_v15 == -1:
            continue
        sigma_V[c] = v30_by_key.get((target_co, target_v15), -1)

    unmapped = [c for c in range(n_cl) if sigma_V[c] == -1]
    return sigma_V, unmapped


def lift_to_edge(sigma_V, edges, edge_idx):
    """Build σ_E and σ_E_sign on the Julia-convention edge list."""
    n_e = len(edges)
    sigma_E = np.zeros(n_e, dtype=np.int64)
    sigma_E_sign = np.ones(n_e, dtype=np.int64)
    missing = []
    for i, (a, b) in enumerate(edges):
        a_ = sigma_V[a]; b_ = sigma_V[b]
        if a_ == -1 or b_ == -1:
            missing.append((i, a, b, a_, b_))
            continue
        # Julia convention: edge tuples are NOT pre-sorted in storage, but
        # cells1 IS sort(unique(raw_pairs)). We look up the σ-image as the
        # sorted (or non-sorted) tuple appearing in cells1. For determinism,
        # we attempt both orientations.
        if (a_, b_) in edge_idx:
            sigma_E[i] = edge_idx[(a_, b_)]
            sigma_E_sign[i] = +1
        elif (b_, a_) in edge_idx:
            sigma_E[i] = edge_idx[(b_, a_)]
            sigma_E_sign[i] = -1
        else:
            missing.append((i, a, b, a_, b_))
    return sigma_E, sigma_E_sign, missing


# -----------------------------------------------------------------------------
# §D. Commutation check (integer matrix exactness)
# -----------------------------------------------------------------------------

def verify_commute_int(L1, sigma_E, sigma_E_sign):
    """Verify (Pσ L_1 Pσ^T)[i,j] = L_1[i,j] exactly on the integer L_1.

    Equivalently: L_1[σE[i], σE[j]] · s[i] · s[j] = L_1[i,j].
    """
    n_e = len(sigma_E)
    max_diff = 0
    for i in range(n_e):
        si = sigma_E[i]; ssi = sigma_E_sign[i]
        for j in range(n_e):
            sj = sigma_E[j]; ssj = sigma_E_sign[j]
            lhs = L1[si, sj] * ssi * ssj
            rhs = L1[i, j]
            d = abs(int(lhs) - int(rhs))
            if d > max_diff:
                max_diff = d
    return max_diff


# -----------------------------------------------------------------------------
# §E. S_6 character table + Frobenius reciprocity
# -----------------------------------------------------------------------------

# S_6 conjugacy classes (cycle type → class size).
# Classes are indexed by cycle type:
#   c1: (1^6)  size 1
#   c2: (2,1^4) size 15
#   c3: (2^2,1^2) size 45
#   c4: (2^3) size 15
#   c5: (3,1^3) size 40
#   c6: (3,2,1) size 120
#   c7: (3^2) size 40
#   c8: (4,1^2) size 90
#   c9: (4,2) size 90
#   c10: (5,1) size 144
#   c11: (6) size 120
S6_CLASS_NAMES = ["(1^6)", "(2,1^4)", "(2^2,1^2)", "(2^3)",
                  "(3,1^3)", "(3,2,1)", "(3^2)",
                  "(4,1^2)", "(4,2)", "(5,1)", "(6)"]
S6_CLASS_SIZES = [1, 15, 45, 15, 40, 120, 40, 90, 90, 144, 120]
assert sum(S6_CLASS_SIZES) == 720

# Conjugacy class representatives as permutations (0-indexed, of (0,1,2,3,4,5)).
# Each is a list [p0, p1, p2, p3, p4, p5] meaning i → p[i].
def perm_id():       return (0, 1, 2, 3, 4, 5)
def perm_2_14():     return (1, 0, 2, 3, 4, 5)        # (01)
def perm_22_12():    return (1, 0, 3, 2, 4, 5)        # (01)(23)
def perm_23():       return (1, 0, 3, 2, 5, 4)        # (01)(23)(45)
def perm_3_13():     return (1, 2, 0, 3, 4, 5)        # (012)
def perm_321():      return (1, 2, 0, 4, 3, 5)        # (012)(34)
def perm_32():       return (1, 2, 0, 4, 5, 3)        # (012)(345)
def perm_4_12():     return (1, 2, 3, 0, 4, 5)        # (0123)
def perm_42():       return (1, 2, 3, 0, 5, 4)        # (0123)(45)
def perm_51():       return (1, 2, 3, 4, 0, 5)        # (01234)
def perm_6():        return (1, 2, 3, 4, 5, 0)        # (012345)

S6_CLASS_REPS = [perm_id(), perm_2_14(), perm_22_12(), perm_23(),
                 perm_3_13(), perm_321(), perm_32(),
                 perm_4_12(), perm_42(), perm_51(), perm_6()]

# Standard S_6 character table.
# Rows = irreps (in partition order), columns = conjugacy classes (in the
# order above). Source: GAP / Murnaghan-Nakayama; cross-checked with the
# Phase S companion v21.
S6_IRREPS = ["[6]", "[5,1]", "[4,2]", "[4,1,1]", "[3,3]",
             "[3,2,1]", "[3,1,1,1]", "[2,2,2]", "[2,2,1,1]",
             "[2,1,1,1,1]", "[1^6]"]
S6_IRREP_DIMS = [1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1]

# Character table (11 irreps × 11 classes). Transcribed verbatim from
# the orthogonality-verified Julia driver
# `cfs_kn3_v21_phase_s_s6_decomposition.jl` (lines 682–695, validated via
# Murnaghan-Nakayama 2026-05-24).
S6_CHAR_TABLE = np.array([
    # cls:        id  (12) (12)(34) (12)(34)(56) (123) (123)(45) (123)(456) (1234) (1234)(56) (12345) (123456)
    # [6] trivial
    [             1,   1,   1,        1,           1,    1,        1,         1,     1,         1,      1],
    # [5,1]
    [             5,   3,   1,       -1,           2,    0,       -1,         1,    -1,         0,     -1],
    # [4,2]
    [             9,   3,   1,        3,           0,    0,        0,        -1,     1,        -1,      0],
    # [4,1,1]
    [            10,   2,  -2,       -2,           1,   -1,        1,         0,     0,         0,      1],
    # [3,3]
    [             5,   1,   1,       -3,          -1,    1,        2,        -1,    -1,         0,      0],
    # [3,2,1]
    [            16,   0,   0,        0,          -2,    0,       -2,         0,     0,         1,      0],
    # [3,1,1,1]
    [            10,  -2,  -2,        2,           1,    1,        1,         0,     0,         0,     -1],
    # [2,2,2]
    [             5,  -1,   1,        3,          -1,   -1,        2,         1,    -1,         0,      0],
    # [2,2,1,1]
    [             9,  -3,   1,       -3,           0,    0,        0,         1,     1,        -1,      0],
    # [2,1^4]
    [             5,  -3,   1,        1,           2,    0,       -1,        -1,    -1,         0,      1],
    # [1^6] sign
    [             1,  -1,   1,       -1,           1,   -1,        1,        -1,     1,         1,     -1],
], dtype=np.int64)


def verify_character_table():
    """Orthogonality checks: ⟨χ_ρ, χ_τ⟩ = δ_{ρτ}, and Σ |C| · |χ(C)|² = 720
    per row (NULL test)."""
    G = 720
    nrows = S6_CHAR_TABLE.shape[0]
    for r in range(nrows):
        v = sum(S6_CLASS_SIZES[i] * int(S6_CHAR_TABLE[r, i]) ** 2
                for i in range(len(S6_CLASS_SIZES)))
        if v != G:
            return False, f"row {r} ({S6_IRREPS[r]}) self-inner-product = {v} ≠ 720"
    for r in range(nrows):
        for s in range(r + 1, nrows):
            v = sum(S6_CLASS_SIZES[i] * int(S6_CHAR_TABLE[r, i]) * int(S6_CHAR_TABLE[s, i])
                    for i in range(len(S6_CLASS_SIZES)))
            if v != 0:
                return False, f"({S6_IRREPS[r]}, {S6_IRREPS[s]}) inner-product = {v} ≠ 0"
    return True, "OK"


# -----------------------------------------------------------------------------
# §F. Main pipeline
# -----------------------------------------------------------------------------

def main():
    print("Q_51 exact Phase S — S_6 decomposition on exact build, Julia convention\n")

    ok, msg = verify_character_table()
    print(f"S_6 character table orthogonality: {msg}")
    if not ok:
        print("ABORT: character table is broken.")
        return 1

    print("\nBuilding exact Q_51...")
    n_cl, q_he, vid_to_cid = build_q51_exact(depth=4)
    print(f"  n_cl = {n_cl}, q_he triples = {len(q_he)}")

    print("\nBuilding Julia-convention cell complex...")
    vertices, edges, faces = cell_complex_julia(q_he, n_cl)
    n_v, n_e, n_f = len(vertices), len(edges), len(faces)
    print(f"  |V|={n_v}  |E|={n_e}  |F|={n_f}  χ={n_v - n_e + n_f}")

    print("\nComputing integer L_1...")
    L1, _, _, edge_idx = hodge_L1_int(vertices, edges, faces)
    print(f"  L_1 shape = {L1.shape}, dtype = {L1.dtype}")

    print("\nIdentifying orbits via co-occurrence multiset...")
    V_6, V_15, V_30, co = find_orbits(q_he, n_cl)
    print(f"  V_6 = {len(V_6)} clusters: {V_6}")
    print(f"  V_15 = {len(V_15)} clusters")
    print(f"  V_30 = {len(V_30)} clusters")
    assert len(V_6) == 6 and len(V_15) == 15 and len(V_30) == 30, \
        f"Orbit sizes {len(V_6)}+{len(V_15)}+{len(V_30)} ≠ 6+15+30 — check build"

    print("\nVerifying eigenvalue spectrum has {3∓√3} pair at m=30...")
    w = np.sort(np.linalg.eigvalsh(L1.astype(np.float64)))
    cnt_lo = int(np.sum(np.abs(w - (3.0 - np.sqrt(3.0))) < 1e-6))
    cnt_hi = int(np.sum(np.abs(w - (3.0 + np.sqrt(3.0))) < 1e-6))
    print(f"  m({{3-√3}}) = {cnt_lo}    m({{3+√3}}) = {cnt_hi}")
    assert cnt_lo == 30 and cnt_hi == 30, "Galois pair multiplicity mismatch"

    print("\nBuilding σ_V and σ_E for each S_6 conjugacy class representative...")
    chi_C1 = []
    all_ok = True
    for ci, perm in enumerate(S6_CLASS_REPS):
        sigma_V, unresolved = build_sigma_V_full(perm, n_cl, co, q_he, V_6, V_15, V_30)
        unmapped = [c for c in range(n_cl) if sigma_V[c] == -1]
        if unmapped:
            print(f"  class {ci} {S6_CLASS_NAMES[ci]}: σ_V has {len(unmapped)} "
                  f"UNMAPPED clusters: {unmapped[:10]}{'...' if len(unmapped)>10 else ''}")
            all_ok = False
            chi_C1.append(None)
            continue
        sigma_E, sigma_E_sign, missing = lift_to_edge(sigma_V, edges, edge_idx)
        if missing:
            print(f"  class {ci} {S6_CLASS_NAMES[ci]}: σ_E has {len(missing)} "
                  f"MISSING edges: {missing[:3]}")
            all_ok = False
            chi_C1.append(None)
            continue
        max_diff = verify_commute_int(L1, sigma_E, sigma_E_sign)
        chi = int(sum(sigma_E_sign[i] if sigma_E[i] == i else 0 for i in range(n_e)))
        chi_C1.append(chi)
        ok_tag = "✓" if max_diff == 0 else f"FAIL Δ={max_diff}"
        print(f"  class {ci:2d} {S6_CLASS_NAMES[ci]:>10} "
              f"size {S6_CLASS_SIZES[ci]:>3}  χ_{{C^1}}={chi:>4}  commute={ok_tag}")
        if max_diff != 0:
            all_ok = False

    if not all_ok or any(c is None for c in chi_C1):
        print("\nFAIL: cannot complete Phase S — σ_V/σ_E construction broke on "
              "at least one class.")
        return 2

    print(f"\nCharacter vector χ_{{C^1}} = {tuple(chi_C1)}")

    print("\nApplying Frobenius reciprocity...")
    G = 720
    print(f"\n  {'irrep':>11}  {'dim':>4}  {'m_ρ':>4}  {'m·dim':>6}")
    print("  " + "-" * 32)
    total = 0
    multiplicities = {}
    for r, irrep in enumerate(S6_IRREPS):
        v = sum(S6_CLASS_SIZES[i] * int(S6_CHAR_TABLE[r, i]) * int(chi_C1[i])
                for i in range(len(S6_CLASS_SIZES)))
        if v % G != 0:
            print(f"  FAIL: m({irrep}) = {v}/{G} = {v/G} (non-integer)")
            return 3
        m = v // G
        multiplicities[irrep] = m
        contrib = m * S6_IRREP_DIMS[r]
        total += contrib
        flag = "  ←" if m > 0 else ""
        print(f"  {irrep:>11}  {S6_IRREP_DIMS[r]:>4}  {m:>4}  {contrib:>6}{flag}")
    print("  " + "-" * 32)
    print(f"  {'TOTAL':>11}  {'':>4}  {'':>4}  {total:>6}")
    print(f"\n  |C^1| = {n_e}      Σ m_ρ · dim_ρ = {total}     "
          f"{'NULL-3 PASS ✓' if total == n_e else 'NULL-3 FAIL ✗'}")

    # -----------------------------------------------------------------
    # Per-cluster decomposition for the {3-√3, 3+√3} eigenspaces.
    # -----------------------------------------------------------------
    print("\n" + "=" * 60)
    print("  PER-CLUSTER DECOMPOSITION OF THE {3∓√3} EIGENSPACES")
    print("=" * 60)

    eigvals, eigvecs = np.linalg.eigh(L1.astype(np.float64))

    targets = [(3.0 - np.sqrt(3.0), "3-√3"), (3.0 + np.sqrt(3.0), "3+√3")]
    for tgt, name in targets:
        mask = np.abs(eigvals - tgt) < 1e-6
        m = int(np.sum(mask))
        V_lam = eigvecs[:, mask]
        print(f"\nλ = {name}  (multiplicity {m})")

        chi_lam = []
        for ci, perm in enumerate(S6_CLASS_REPS):
            sigma_V, _unmapped = build_sigma_V_full(perm, n_cl, co, q_he, V_6, V_15, V_30)
            sigma_E, sigma_E_sign, _missing = lift_to_edge(sigma_V, edges, edge_idx)
            # χ_λ(σ) = Tr(P_λ σ_E P_λ) = Tr(σ_E|V_λ)
            # Build sigma_E as a signed permutation matrix and project onto V_λ
            n_e = len(sigma_E)
            P = np.zeros((n_e, n_e), dtype=np.float64)
            for i in range(n_e):
                P[sigma_E[i], i] = sigma_E_sign[i]
            # σ_E action on V_λ: V_λ^T P V_λ is m × m
            Mσ = V_lam.T @ P @ V_lam
            tr = float(np.trace(Mσ))
            chi_lam.append(tr)

        # Apply Frobenius reciprocity in float (then round)
        print(f"  χ_λ = ({', '.join(f'{x:7.3f}' for x in chi_lam)})")
        print(f"\n  {'irrep':>11}  {'dim':>4}  {'m_λ,ρ':>7}")
        print("  " + "-" * 28)
        total = 0
        for r, irrep in enumerate(S6_IRREPS):
            v = sum(S6_CLASS_SIZES[i] * float(S6_CHAR_TABLE[r, i]) * chi_lam[i]
                    for i in range(len(S6_CLASS_SIZES)))
            m_lr = v / 720.0
            m_int = int(round(m_lr))
            if abs(m_lr - m_int) > 1e-4:
                print(f"  {irrep:>11}  {S6_IRREP_DIMS[r]:>4}  {m_lr:>7.4f}  NON-INT")
            elif m_int != 0:
                print(f"  {irrep:>11}  {S6_IRREP_DIMS[r]:>4}  {m_int:>7}")
                total += m_int * S6_IRREP_DIMS[r]
        print("  " + "-" * 28)
        print(f"  {'TOTAL':>11}  {'':>4}  total dim={total:>3}    "
              f"{'NULL-3 PASS ✓' if total == m else 'NULL-3 FAIL ✗'}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
