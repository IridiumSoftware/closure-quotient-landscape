#!/usr/bin/env python3
"""
q51_exact_vs_fidelity_verification.py

Compare two builds of "Q_51" and verify which (if either) the paper note
v2 has been using:

  (A) EXACT: K_6^3 multiway rewriting + exact Gaussian-integer ray-equivalence
      clustering (Q_51 per closure-v5 spec S157 / S178, single-chirality
      orig-side of the exact build_c_closed_quotient).

  (B) FIDELITY: K_6^3 multiway rewriting with Haar Float64 colours +
      fidelity > 0.999 threshold clustering, depth=5, no C-closure (the
      TCE build_q51 currently used in the cross-substrate note).

For each, compute the Hodge 1-Laplacian L_1 on the q_he-induced cell
complex and report:
  * |V|, |E|, |F|, Euler χ
  * L_1 spectrum (eigenvalue clusters with multiplicities)
  * Whether {3-√3, 3+√3} appears

Then report the structural diff between the two builds (cluster-merge
table, hyperedge symmetric difference).

USAGE: set CFS_REPO_ROOT to a Closure-Forces-Structure source checkout
  (https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs)
  then run this driver from anywhere:
    export CFS_REPO_ROOT=/path/to/closure-forces-structure
    python3 q51_exact_vs_fidelity_verification.py
"""

import sys, os
# Add the CFS source root to sys.path so we can import its exact build.
CLOSURE_V5 = os.environ.get("CFS_REPO_ROOT")
if not CLOSURE_V5:
    sys.exit("Set CFS_REPO_ROOT env var to a Closure-Forces-Structure checkout.")
sys.path.insert(0, CLOSURE_V5)

import numpy as np
from collections import defaultdict, Counter

# Import the exact build infrastructure
from q102_build_exact_v1 import (
    gaussian_ics, build_multiway, proj_equiv, is_zero_vec,
    _to_c128, complete_ternary, _gconj, haar_C3, fidelity
)


def build_q51_exact(depth=4, ic_seed=None):
    """Q_51 via exact Gaussian-integer arithmetic (orig-side / no C-closure).

    Per closure-v5 q102_build_exact_v1.py docstring: 'K6^3 single-side -> 51'.
    Returns dict with n_cl, q_psi_exact, q_he, vid_to_cid.
    """
    n_seed = 6
    ics = gaussian_ics(n=n_seed, ic_seed=ic_seed)
    seed_edges = complete_ternary(n_seed)   # K_6^3 = 120 ordered triples
    psi_o, depth_o, edges_o = build_multiway(seed_edges, ics, depth)

    # Exact ray-equivalence clustering (NO threshold)
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
    q_psi_exact = {i: clusters[i][1] for i in range(n_cl)}

    # Cluster-level hyperedges (sorted triples, unique)
    q_he = set()
    for (_, v1, v2, v3) in edges_o:
        if (v1 in vid_to_cid and v2 in vid_to_cid and v3 in vid_to_cid):
            t = tuple(sorted([vid_to_cid[v1], vid_to_cid[v2], vid_to_cid[v3]]))
            q_he.add(t)
    q_he = sorted(q_he)

    return dict(n_cl=n_cl, q_psi_exact=q_psi_exact, q_he=q_he,
                vid_to_cid=vid_to_cid, label="EXACT")


def build_q51_fidelity(depth=5, threshold=0.999, ic_seed=0):
    """TCE-style: Haar Float64 + fidelity-clustering + depth=5 + no C-closure."""
    rng = np.random.default_rng(ic_seed)
    # Mirror the Julia haar_C3 (C^3 unit vector)
    psi_init = {v: haar_C3(rng) for v in range(6)}

    seed_edges = complete_ternary(6)
    # Build multiway in Float64 (mirrors TCE build_quotient_color)
    psi = {k: np.array(v, dtype=np.complex128) for k, v in psi_init.items()}
    next_vid = max(psi.keys()) + 1
    edges = [(0, t[0], t[1], t[2]) for t in seed_edges]
    cache = {}
    for d in range(depth):
        d_edges = [e for e in edges if e[0] == d]
        for (_, v1, v2, v3) in d_edges:
            key = (v1, v2)
            if key not in cache:
                w_psi = np.cross(np.conj(psi[v1]), np.conj(psi[v2]))
                # The closure-v5 compose is conj(a × b)
                # (TCE uses compose_colour which is compose conj(a×b) per the seed)
                if np.linalg.norm(w_psi) < 1e-10:
                    continue
                w_psi = w_psi / np.linalg.norm(w_psi)
                psi[next_vid] = w_psi
                cache[key] = next_vid
                next_vid += 1
            w = cache[key]
            edges.append((d+1, w, v2, v3))
            edges.append((d+1, w, v1, v3))
            edges.append((d+1, w, v1, v2))

    # Fidelity clustering (threshold=0.999)
    all_vids = sorted(psi.keys())
    clusters = []  # (rep_vid, rep_vec)
    vid_to_cid = {}
    for v in all_vids:
        pv = psi[v]
        if np.linalg.norm(pv) < 1e-10:
            continue
        cid = -1
        for ci, (_, rep) in enumerate(clusters):
            if fidelity(pv, rep) > threshold:
                cid = ci
                break
        if cid < 0:
            cid = len(clusters)
            clusters.append((v, pv))
        vid_to_cid[v] = cid

    n_cl = len(clusters)
    q_psi = {i: clusters[i][1] for i in range(n_cl)}

    # Cluster-level hyperedges
    q_he = set()
    for (_, v1, v2, v3) in edges:
        if (v1 in vid_to_cid and v2 in vid_to_cid and v3 in vid_to_cid):
            t = tuple(sorted([vid_to_cid[v1], vid_to_cid[v2], vid_to_cid[v3]]))
            q_he.add(t)
    q_he = sorted(q_he)

    return dict(n_cl=n_cl, q_psi=q_psi, q_he=q_he,
                vid_to_cid=vid_to_cid, label="FIDELITY")


def cell_complex(q_he, n_v):
    """Extract (vertices, edges, faces) from q_he triples — PYTHON convention.

    Canonicalizes faces (sorted unordered triples) and edges (sorted unordered
    pairs). Drops self-loops and degenerate triples. This is the convention used
    by the Python verification path.
    """
    edges = set()
    faces = []
    for t in q_he:
        a, b, c = sorted(t)
        # Self-edges (a==a etc.) and triples with repeated vertices give
        # degenerate 2-cells; we still record them per the q_he convention.
        edges.add((a, b)); edges.add((a, c)); edges.add((b, c))
        faces.append((a, b, c))
    edges = sorted(set(e for e in edges if e[0] != e[1]))   # exclude self-loops
    faces = sorted(set(f for f in faces if len(set(f)) == 3))  # exclude degenerates
    return list(range(n_v)), edges, faces


def cell_complex_julia(q_he, n_v):
    """Julia convention: edges built from RAW ORDERED triples without
    pre-canonicalization, faces preserve raw triples (no dedup of unordered
    permutations). Used by gpg_b1_embedding_test_lib.jl::CellComplex.

    Mimics the Julia code:
        edges = NTuple{2,Int}[]
        for (a,b,c) in H.E3
            push!(edges, (a,b)); push!(edges, (b,c)); push!(edges, (a,c))
        end
        cells1 = sort(unique(edges))  # tuple equality is ORDER-SENSITIVE
        cells2 = H.E3                 # faces kept as-is (no dedup)
    """
    edges = []
    faces = []
    for t in q_he:
        a, b, c = t  # raw order — NOT sorted
        edges.append((a, b))
        edges.append((b, c))
        edges.append((a, c))
        faces.append((a, b, c))
    # Order-sensitive uniqueness (matches Julia tuple equality)
    edges = sorted(set(edges))
    # Drop self-loops (a==b etc.) — Julia code retains these but they contribute
    # zero rows to d_0 and d_1 in any cell complex. For honest comparison we
    # check both with and without; keep all for now to match Julia output.
    faces_unique = sorted(set(faces))  # only deduplicate exact-tuple repeats
    return list(range(n_v)), edges, faces_unique


def hodge_L1_julia(vertices, edges, faces):
    """L_1 with the Julia driver's d_1 convention:
       for each face (a, b, c) [raw order], d_1[j, e_ab] = +1,
       d_1[j, e_ac] = -1, d_1[j, e_bc] = +1.

    Each (a,b) edge index is looked up under ORDERED tuple key — if (b,a) is
    in the edge list and (a,b) is not, the lookup returns None and that face
    row gets zero for that edge.
    """
    n_e = len(edges)
    n_f = len(faces)
    edge_idx = {e: i for i, e in enumerate(edges)}

    d0 = np.zeros((n_e, len(vertices)), dtype=np.int64)
    for i, (a, b) in enumerate(edges):
        if a == b:
            continue  # self-loop, drop
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
    return L1.astype(np.float64), d0, d1


def hodge_L1(vertices, edges, faces):
    """L_1 = d_0 d_0^T + d_1^T d_1 with alternating-sign d_1 on sorted faces."""
    n_e = len(edges)
    n_f = len(faces)
    edge_idx = {e: i for i, e in enumerate(edges)}

    # d_0: n_e × n_v
    d0 = np.zeros((n_e, len(vertices)), dtype=np.int64)
    for i, (a, b) in enumerate(edges):
        d0[i, a] = -1; d0[i, b] = +1

    # d_1: n_f × n_e  with +1 (a,b), -1 (a,c), +1 (b,c) on sorted face (a<b<c)
    d1 = np.zeros((n_f, n_e), dtype=np.int64)
    for j, (a, b, c) in enumerate(faces):
        e_ab = edge_idx.get((a, b))
        e_ac = edge_idx.get((a, c))
        e_bc = edge_idx.get((b, c))
        if e_ab is not None: d1[j, e_ab] = +1
        if e_ac is not None: d1[j, e_ac] = -1
        if e_bc is not None: d1[j, e_bc] = +1

    L1 = d0 @ d0.T + d1.T @ d1
    return L1.astype(np.float64), d0, d1


def cluster_eigvals(L1, tol=1e-6):
    """Return list of (eigenvalue, multiplicity)."""
    w = np.sort(np.linalg.eigvalsh(L1))
    clusters = []
    cur = w[0]; cnt = 1
    for x in w[1:]:
        if abs(x - cur) < tol:
            cnt += 1
        else:
            clusters.append((cur, cnt))
            cur = x; cnt = 1
    clusters.append((cur, cnt))
    return clusters


def check_galois_pair(clusters, tol=1e-4):
    """Search for the {3-√3, 3+√3} pair."""
    sqrt3 = np.sqrt(3.0)
    target_lo = 3.0 - sqrt3   # ≈ 1.2679
    target_hi = 3.0 + sqrt3   # ≈ 4.7321
    found_lo = None; found_hi = None
    for (lam, m) in clusters:
        if abs(lam - target_lo) < tol:
            found_lo = (lam, m)
        if abs(lam - target_hi) < tol:
            found_hi = (lam, m)
    return found_lo, found_hi


def summarise(build):
    """Build → cell complex → spectrum → report."""
    print("=" * 60)
    print(f"  {build['label']} Q_51 build")
    print("=" * 60)
    print(f"  n_cl = {build['n_cl']}")
    print(f"  q_he triples = {len(build['q_he'])}")

    vs, es, fs = cell_complex(build['q_he'], build['n_cl'])
    chi = len(vs) - len(es) + len(fs)
    print(f"  |V|={len(vs)}  |E|={len(es)}  |F|={len(fs)}  χ={chi}")

    L1, _, _ = hodge_L1(vs, es, fs)
    clusters = cluster_eigvals(L1)
    print(f"  L_1 dim = {L1.shape[0]}")
    print(f"  Distinct eigenvalue clusters = {len(clusters)}")
    for k, (lam, m) in enumerate(clusters):
        sym = ""
        if abs(lam - (3 - np.sqrt(3))) < 1e-4: sym = "  ← 3 − √3"
        elif abs(lam - (3 + np.sqrt(3))) < 1e-4: sym = "  ← 3 + √3"
        print(f"    k={k:>2} λ={lam:14.8f}  m={m}{sym}")

    lo, hi = check_galois_pair(clusters)
    if lo and hi:
        print(f"  ✓ {{3∓√3}} PAIR FOUND: m_lo={lo[1]}  m_hi={hi[1]}")
    else:
        print(f"  ✗ {{3∓√3}} PAIR NOT FOUND   (lo={lo}, hi={hi})")

    return dict(n_v=len(vs), n_e=len(es), n_f=len(fs), chi=chi,
                clusters=clusters, q_he=build['q_he'], n_cl=build['n_cl'])


def structural_diff(b_exact, b_fid):
    """Compare the two builds at the q_he level."""
    print("\n" + "=" * 60)
    print("  STRUCTURAL DIFFERENCE")
    print("=" * 60)

    # Can't directly compare q_he (different cluster IDs in each build).
    # Compare on shape only.
    ne_e = b_exact['n_e']; ne_f = b_exact['n_f']
    nf_e = b_fid['n_e']; nf_f = b_fid['n_f']
    print(f"  n_cl:        exact={b_exact['n_cl']:>4}   fidelity={b_fid['n_cl']:>4}"
          f"   Δ={b_exact['n_cl']-b_fid['n_cl']:+d}")
    print(f"  |E|:         exact={ne_e:>4}   fidelity={nf_e:>4}"
          f"   Δ={ne_e-nf_e:+d}")
    print(f"  |F|:         exact={ne_f:>4}   fidelity={nf_f:>4}"
          f"   Δ={ne_f-nf_f:+d}")

    # Spectrum-level comparison
    ec = b_exact['clusters']; fc = b_fid['clusters']
    print(f"\n  eigenvalue clusters: exact={len(ec)}  fidelity={len(fc)}")
    if len(ec) == len(fc):
        max_diff = max(abs(ec[i][0] - fc[i][0]) for i in range(len(ec)))
        max_mult_diff = max(abs(ec[i][1] - fc[i][1]) for i in range(len(ec)))
        print(f"  matching by index: max |Δλ|={max_diff:.2e}  "
              f"max |Δm|={max_mult_diff}")
    return None


def summarise_julia(build):
    """Build → Julia-convention cell complex → spectrum → report."""
    print("=" * 60)
    print(f"  {build['label']} Q_51 build  (JULIA convention)")
    print("=" * 60)
    print(f"  n_cl = {build['n_cl']}")
    print(f"  q_he triples = {len(build['q_he'])}")

    vs, es, fs = cell_complex_julia(build['q_he'], build['n_cl'])
    chi = len(vs) - len(es) + len(fs)
    print(f"  |V|={len(vs)}  |E|={len(es)}  |F|={len(fs)}  χ={chi}")

    L1, _, _ = hodge_L1_julia(vs, es, fs)
    clusters = cluster_eigvals(L1)
    print(f"  L_1 dim = {L1.shape[0]}")
    print(f"  Distinct eigenvalue clusters = {len(clusters)}")
    for k, (lam, m) in enumerate(clusters):
        sym = ""
        if abs(lam - (3 - np.sqrt(3))) < 1e-4: sym = "  ← 3 − √3"
        elif abs(lam - (3 + np.sqrt(3))) < 1e-4: sym = "  ← 3 + √3"
        print(f"    k={k:>2} λ={lam:14.8f}  m={m}{sym}")

    lo, hi = check_galois_pair(clusters)
    if lo and hi:
        print(f"  ✓ {{3∓√3}} PAIR FOUND: m_lo={lo[1]}  m_hi={hi[1]}")
    else:
        print(f"  ✗ {{3∓√3}} PAIR NOT FOUND   (lo={lo}, hi={hi})")

    return dict(n_v=len(vs), n_e=len(es), n_f=len(fs), chi=chi,
                clusters=clusters, q_he=build['q_he'], n_cl=build['n_cl'])


if __name__ == "__main__":
    print("Q_51: EXACT vs FIDELITY comparison\n")
    print("Building EXACT Q_51 (Gaussian-integer ray-equivalence, depth=4)...")
    e = build_q51_exact(depth=4)
    re = summarise(e)

    print("\nBuilding FIDELITY Q_51 (Haar Float64 + threshold=0.999, depth=5)...")
    f = build_q51_fidelity(depth=5, threshold=0.999, ic_seed=0)
    rf = summarise(f)

    print("\n" + "=" * 60)
    print("  JULIA-CONVENTION RE-RUN (raw-triple edges, |E|=366 expected for fidelity)")
    print("=" * 60)
    print("\nEXACT under Julia convention:")
    re_j = summarise_julia(e)
    print("\nFIDELITY under Julia convention:")
    rf_j = summarise_julia(f)

    structural_diff(re, rf)

    # Print key verdict
    print("\n" + "=" * 60)
    print("  VERDICT")
    print("=" * 60)
    lo_e, hi_e = check_galois_pair(re['clusters'])
    lo_f, hi_f = check_galois_pair(rf['clusters'])
    print(f"  EXACT    {{3∓√3}} present: {bool(lo_e and hi_e)}"
          + (f"  mults ({lo_e[1]}, {hi_e[1]})" if lo_e and hi_e else ""))
    print(f"  FIDELITY {{3∓√3}} present: {bool(lo_f and hi_f)}"
          + (f"  mults ({lo_f[1]}, {hi_f[1]})" if lo_f and hi_f else ""))

    if (lo_e and hi_e) and (lo_f and hi_f) \
       and lo_e[1] == lo_f[1] and hi_e[1] == hi_f[1]:
        print("\n  → BUILDS AGREE on the Galois pair structure for Q_51.")
        print("    Paper note's claim survives, but should cite EXACT build for honesty.")
    elif (lo_e and hi_e) and not (lo_f and hi_f):
        print("\n  → FIDELITY MISSES the pair; EXACT carries it.")
        print("    Paper note must switch to EXACT build to make its claim.")
    elif not (lo_e and hi_e) and (lo_f and hi_f):
        print("\n  → FIDELITY artefact: the pair appears in fidelity build but NOT in exact.")
        print("    Paper note's main observation collapses; need to re-write.")
    else:
        print("\n  → BOTH builds give Galois pair, but multiplicities differ.")
        print("    Paper note needs updated multiplicities + an EXACT build switch.")
