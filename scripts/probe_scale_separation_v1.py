#!/usr/bin/env python3
# Probe — Q_51 = K_6^3 carries Rosen "scale separation / coarse-graining".
#
# Rosen (M,R)-system role: scale separation = identity survives
# coarse-graining to lower scales. In CFS: Q_102 -> Q_51 by sector-strip
# (forget the conjugate sector) gives a structurally exact projection;
# the L_1 spectrum (Hodge Laplacian on the 1-cochain space) of Q_102
# decomposes as 2 × spectrum(Q_51) by J-symmetry under fpf ICs; K_n^3
# closure is depth-stable across the family.
#
# Hypothesis:
#   (i)  Sector-strip identity: Q_102's 51 orig_only clusters match
#        Q_51's clusters cluster-for-cluster (projective gauge eq).
#   (ii) L_1 spectrum: Q_102 = 2 × Q_51 — each Q_51 eigenvalue doubles
#        multiplicity in Q_102 (block structure under J-symmetry).
#   (iii) Depth-stability: K_n^3 -> n(3n-1)/2 stable across d=3,4,5
#         for n=5,6,7 (the canonical-density regime).
#
# Reading α V_k decomposition machinery for Q_51 already lives at
# scripts/q51_exact_phase_s_python.py (S_6 char decomp, dim-336 sum
# verified). For Q_102, the V_k = 2 × Q_51 V_k follows from (i)+(ii)
# via the J-symmetry sector-strip argument.

import os
import sys
from pathlib import Path

# Set CFS_REPO_ROOT to a Closure-Forces-Structure source checkout
# (https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs)
# before invoking this probe.
CLOSURE_V5 = os.environ.get("CFS_REPO_ROOT")
if not CLOSURE_V5:
    sys.exit("Set CFS_REPO_ROOT env var to a Closure-Forces-Structure checkout.")
SCRIPTS = Path(CLOSURE_V5) / "paper" / "scripts"
sys.path.insert(0, CLOSURE_V5)
sys.path.insert(0, str(SCRIPTS))

import numpy as np

from q102_build_exact_v1 import (
    build_c_closed_quotient, complete_ternary, build_J,
    gaussian_ics, build_multiway, proj_equiv, is_zero_vec, _gconj,
)
from q51_exact_phase_s_python import (
    build_q51_exact, cell_complex_julia, hodge_L1_int,
)


# --- Helper: compute |single-side Q_n^3| at depth d --------------------------
def single_side_count(seed_edges, n_seed_verts, depth, ic_seed=None):
    ics = gaussian_ics(n=n_seed_verts, ic_seed=ic_seed)
    psi, _, _ = build_multiway(seed_edges, ics, depth)
    # exact ray-equiv clustering
    clusters = []
    for v in sorted(psi.keys()):
        pv = psi[v]
        if is_zero_vec(pv): continue
        found = False
        for rep in clusters:
            if proj_equiv(pv, rep):
                found = True
                break
        if not found:
            clusters.append(pv)
    return len(clusters)


# -----------------------------------------------------------------------------
# Step 1 — Build Q_51 exact + L_1
# -----------------------------------------------------------------------------
print("=" * 78)
print("Probe — Q_51 = K_6^3 carries scale separation (sector-strip identity)")
print("=" * 78)
print()

n51, q51_he, _q51_v_to_c = build_q51_exact(depth=4)
print(f"Q_51: {n51} clusters, {len(q51_he)} q_he triples")
assert n51 == 51, f"Q_51 build failed: {n51} != 51"

# Cell complex + L_1
verts51, edges51, faces51 = cell_complex_julia(q51_he, 51)
L51, _d0_51, _d1_51, _ei51 = hodge_L1_int(verts51, edges51, faces51)
print(f"   cell complex: |V|=51, |E|={len(edges51)}, |F|={len(faces51)}")
print(f"   L_1 shape: {L51.shape}, integer entries")
spec51 = np.linalg.eigvalsh(L51.astype(np.float64))
spec51_sorted = np.sort(spec51)
print(f"   L_1 spectrum range: [{spec51_sorted[0]:.6f}, {spec51_sorted[-1]:.6f}]")
n_zero_51 = int(np.sum(np.abs(spec51_sorted) < 1e-9))
print(f"   |ker L_1(Q_51)|: {n_zero_51}")
print()

# -----------------------------------------------------------------------------
# Step 2 — Build Q_102 exact, label by sector
# -----------------------------------------------------------------------------
print("Building Q_102 (C-closed K_6^3)...")
Q102 = build_c_closed_quotient(complete_ternary(6), depth=4, n_seed_verts=6)
assert Q102['n_cl'] == 102, f"Q_102 build failed: {Q102['n_cl']} != 102"

orig_only = sorted(c for c in range(102) if Q102['cl_origin'][c] == 'orig_only')
conj_only = sorted(c for c in range(102) if Q102['cl_origin'][c] == 'conj_only')
both       = sorted(c for c in range(102) if Q102['cl_origin'][c] == 'both')

print(f"Q_102: orig_only={len(orig_only)}, conj_only={len(conj_only)}, both={len(both)}")
print(f"   Expected: 51, 51, 0 (J fixed-point-free under canonical Gaussian ICs)")

# -----------------------------------------------------------------------------
# Step 3 — J fixed-point-free check
# -----------------------------------------------------------------------------
J_mat, jmap = build_J(Q102)
fpf = all(jmap.get(c, c) != c for c in range(102))
print(f"J fpf on Q_102: {'YES' if fpf else 'NO'}")
# Confirm J maps orig <-> conj
jmap_orig_to_conj = all(jmap[c] in conj_only for c in orig_only)
jmap_conj_to_orig = all(jmap[c] in orig_only for c in conj_only)
print(f"J maps orig <-> conj: {'YES' if (jmap_orig_to_conj and jmap_conj_to_orig) else 'NO'}")
print()

# -----------------------------------------------------------------------------
# Step 4 — Sector-strip identity check: orig_only cluster reps match Q_51's
# -----------------------------------------------------------------------------
print("Sector-strip identity: orig_only clusters of Q_102 match Q_51?")
# Get Q_51's exact cluster reps from build_q51_exact (we have q_he but not
# directly the reps; rebuild from build_q51_exact's psi)
# Easier path: re-build Q_51 with same canonical ICs, get reps, and match
# against Q_102's orig_only cluster reps.
ics = gaussian_ics(n=6, ic_seed=None)
psi51, _, _ = build_multiway(complete_ternary(6), ics, 4)
q51_reps = []
for v in sorted(psi51.keys()):
    pv = psi51[v]
    if is_zero_vec(pv): continue
    found = False
    for rep in q51_reps:
        if proj_equiv(pv, rep):
            found = True; break
    if not found:
        q51_reps.append(pv)
assert len(q51_reps) == 51

q102_orig_reps = [Q102['q_psi_exact'][c] for c in orig_only]
assert len(q102_orig_reps) == 51

# For each Q_51 rep, find a matching Q_102 orig_only rep (proj_equiv)
matched = 0
unmatched_q51 = []
for r51 in q51_reps:
    if any(proj_equiv(r51, r102) for r102 in q102_orig_reps):
        matched += 1
    else:
        unmatched_q51.append(r51)
print(f"   {matched}/51 Q_51 cluster reps proj-equivalent to a Q_102 orig_only rep")
strip_verdict = "PASS" if matched == 51 else "FAIL"
print(f"   Sector-strip identity: {strip_verdict}")
print()

# -----------------------------------------------------------------------------
# Step 5 — L_1 spectrum: Q_102 vs 2 × Q_51
# -----------------------------------------------------------------------------
print("L_1 spectrum comparison (Q_102 vs 2 × Q_51):")
q102_he_sorted = sorted(Q102['q_he'])
verts102, edges102, faces102 = cell_complex_julia(q102_he_sorted, 102)
L102, _, _, _ = hodge_L1_int(verts102, edges102, faces102)
print(f"   Q_102 cell complex: |V|=102, |E|={len(edges102)}, |F|={len(faces102)}")
print(f"   L_1(Q_102) shape: {L102.shape}")
spec102 = np.sort(np.linalg.eigvalsh(L102.astype(np.float64)))

# Compare doubled-Q_51 vs Q_102: each Q_51 eigenvalue should appear with
# doubled multiplicity in Q_102 IF block-diagonal under J.
spec51_doubled = np.sort(np.concatenate([spec51_sorted, spec51_sorted]))
if len(spec51_doubled) == len(spec102):
    diff = np.abs(spec51_doubled - spec102).max()
    print(f"   |2×spec(Q_51) − spec(Q_102)|_inf = {diff:.6e}")
    spec_match = diff < 1e-6
    print(f"   Block-diagonal match: {'YES' if spec_match else 'NO'}")
else:
    spec_match = False
    print(f"   dim mismatch: 2×|spec(Q_51)|={len(spec51_doubled)} vs "
          f"|spec(Q_102)|={len(spec102)} — block-diagonal claim FAILS "
          f"(edges/faces dimensions differ: cross-sector hyperedges present)")
print()

# -----------------------------------------------------------------------------
# Step 6 — Depth-stability of K_n^3 closure for n=5,6,7
# -----------------------------------------------------------------------------
print("K_n^3 depth-stability (single-side):")
print(f"{'n':>3} {'pred n(3n-1)/2':>14}", end="")
for d in (3, 4, 5):
    print(f" {'d='+str(d):>7}", end="")
print(" stable")
print("-" * 60)

depth_stable_all = True
for n in (5, 6, 7):
    pred = n * (3*n - 1) // 2
    print(f"{n:>3} {pred:>14}", end="")
    counts = []
    for d in (3, 4, 5):
        q = single_side_count(complete_ternary(n), n, d, ic_seed=None)
        counts.append(q)
        print(f" {q:>7}", end="")
    is_stable = (counts[-1] == pred and counts[-2] == pred)
    if not is_stable: depth_stable_all = False
    print(f"   {'Y' if is_stable else 'N'}")
print()

# -----------------------------------------------------------------------------
# Verdict
# -----------------------------------------------------------------------------
print("=" * 78)
print("Verdict:")
print(f"  (i)  Sector-strip identity (Q_102 orig_only = Q_51): "
      f"{strip_verdict}")
print(f"  (ii) Spectrum 2×Q_51 = Q_102: "
      f"{'PASS' if spec_match else 'FAIL'}")
print(f"  (iii) K_n^3 depth-stable n=5,6,7 at d=4,5: "
      f"{'PASS' if depth_stable_all else 'FAIL'}")
print()
overall = "PASS" if (strip_verdict == "PASS"
                     and spec_match
                     and depth_stable_all) else "PARTIAL/FAIL"
print(f"Hypothesis 'Q_51 carries scale separation via sector-strip + "
      f"depth-stability': {overall}")
