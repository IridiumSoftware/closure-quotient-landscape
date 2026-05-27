#!/usr/bin/env python3
# Probe 5b — J-decomposition of Q_102's L_1.
#
# Resolves the Probe 5 finding: Q_102's L_1 (840-dim) is NOT simply
# 2 × Q_51's L_1 (would be 630-dim) — cross-sector hyperedges contribute
# 210 additional edges to the cell complex. The richer structure is
# resolved by decomposing L_1 under the charge-conjugation involution J:
#
#   L_1(Q_102) = L_+ ⊕ L_-     on   V_+ ⊕ V_-
#
# where V_+ is the J-symmetric subspace and V_- is the J-antisymmetric
# subspace. By construction J_edge L_1 = L_1 J_edge (verified).
#
# Predicted: V_+ spectrum ≡ Q_51's L_1 spectrum (J-symmetric carries the
# cluster-level scale-separation identity). V_- spectrum is the
# cross-sector "mixing" structure with no Q_51 counterpart.
#
# Q_51 sector-strip carries the cluster-level identity (Probe 5 (i));
# V_+ spectrum match confirms the spectral-level identity. The "scale
# separation" claim is sharpened to: identity persists in the
# J-symmetric block of Q_102's L_1.

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
    build_c_closed_quotient, build_J, complete_ternary,
    gaussian_ics, build_multiway, proj_equiv, is_zero_vec,
)
from q51_exact_phase_s_python import (
    build_q51_exact, cell_complex_julia, hodge_L1_int,
)


# --- Build Q_51 and its L_1 -------------------------------------------------
n51, q51_he, _ = build_q51_exact(depth=4)
assert n51 == 51
verts51, edges51, faces51 = cell_complex_julia(q51_he, 51)
L51, _, _, _ = hodge_L1_int(verts51, edges51, faces51)
spec51 = np.sort(np.linalg.eigvalsh(L51.astype(np.float64)))
print(f"Q_51:  |V|=51, |E|={len(edges51)}, |F|={len(faces51)}, "
      f"L_1 dim={L51.shape[0]}")
print(f"   spectrum range: [{spec51[0]:.6f}, {spec51[-1]:.6f}]")
print()

# --- Build Q_102 + J map + L_1 ----------------------------------------------
Q102 = build_c_closed_quotient(complete_ternary(6), depth=4, n_seed_verts=6)
assert Q102['n_cl'] == 102
_J_mat, jmap = build_J(Q102)
# Verify J is fixed-point-free and an involution
fpf = all(jmap[c] != c for c in range(102))
involution = all(jmap[jmap[c]] == c for c in range(102))
print(f"Q_102: |Q|=102, J fpf={fpf}, J involution={involution}")

q102_he_sorted = sorted(Q102['q_he'])
verts102, edges102, faces102 = cell_complex_julia(q102_he_sorted, 102)
L102, _, _, edge_idx = hodge_L1_int(verts102, edges102, faces102)
n_e = L102.shape[0]
print(f"        |V|=102, |E|={n_e}, |F|={len(faces102)}")
print()

# --- Lift J to edge space (oriented permutation) ----------------------------
#
# Edges in cell_complex_julia preserve original orientation from the q_he
# triple (a, b, c) -> (a, b), (b, c), (a, c). J acts on clusters via
# jmap. The forward image of edge (a, b) is (jmap[a], jmap[b]), which
# IS present in the cell complex by J-symmetry of q_he (orig and conj
# sectors generate mirrored triples). J preserves orientation, so the
# action on the oriented-edge basis is an unsigned permutation.
J_edge = np.zeros((n_e, n_e), dtype=np.int64)
missing = 0
for i, (a, b) in enumerate(edges102):
    ja, jb = jmap[a], jmap[b]
    j = edge_idx.get((ja, jb))
    if j is None:
        missing += 1
        continue
    J_edge[j, i] = +1

if missing:
    print(f"WARNING: {missing} edges had no J-image in the cell complex")

# Verify J_edge is an involution: J_edge^2 = I
J_sq = J_edge @ J_edge
ok_inv = np.allclose(J_sq, np.eye(n_e), atol=0)
print(f"J_edge involution (J^2 = I): {ok_inv}")
print(f"J_edge orthogonal           : {np.allclose(J_edge @ J_edge.T, np.eye(n_e), atol=0)}")

# Verify [J_edge, L_1] = 0
commute = np.array_equal(J_edge @ L102, L102 @ J_edge)
print(f"[J_edge, L_1] = 0           : {commute}")
print()

if not commute:
    # If J_edge does not commute exactly with L_1, the J-decomposition is
    # not block-diagonal under J. Report and exit.
    delta = J_edge @ L102 - L102 @ J_edge
    max_diff = np.abs(delta).max()
    print(f"Commutator max|·| = {max_diff} — J-symmetry NOT exact.")
    print("Proceeding with eigenspace decomp anyway for diagnostic purposes.")
    print()

# --- Compute V_+ and V_- bases via J_edge eigenvectors -----------------------
J_edge_f = J_edge.astype(np.float64)
eig_vals_J, eig_vecs_J = np.linalg.eigh(J_edge_f)
# Round eigenvalues to ±1
pos = np.where(np.abs(eig_vals_J - 1.0) < 1e-8)[0]
neg = np.where(np.abs(eig_vals_J + 1.0) < 1e-8)[0]
print(f"V_+ dim (J=+1): {len(pos)}  V_- dim (J=-1): {len(neg)}  "
      f"sum = {len(pos)+len(neg)} (should equal {n_e})")

V_plus  = eig_vecs_J[:, pos]   # n_e × dim(V_+)
V_minus = eig_vecs_J[:, neg]   # n_e × dim(V_-)

# --- Restricted L_1 spectra --------------------------------------------------
L_plus  = V_plus.T  @ L102.astype(np.float64) @ V_plus
L_minus = V_minus.T @ L102.astype(np.float64) @ V_minus
spec_plus  = np.sort(np.linalg.eigvalsh((L_plus  + L_plus.T)  / 2))
spec_minus = np.sort(np.linalg.eigvalsh((L_minus + L_minus.T) / 2))
print(f"L_+ shape: {L_plus.shape}, L_- shape: {L_minus.shape}")
print(f"L_+ spectrum range: [{spec_plus[0]:.6f}, {spec_plus[-1]:.6f}]")
print(f"L_- spectrum range: [{spec_minus[0]:.6f}, {spec_minus[-1]:.6f}]")
print()

# --- Compare V_+ spectrum to Q_51's L_1 spectrum ----------------------------
print("Comparison: V_+ spectrum vs Q_51's L_1 spectrum")
if len(spec_plus) == len(spec51):
    diff = np.abs(spec_plus - spec51).max()
    print(f"  dimensions match: {len(spec_plus)}")
    print(f"  |spec(L_+) − spec(Q_51 L_1)|_inf = {diff:.6e}")
    plus_matches = diff < 1e-6
    print(f"  V_+ spectrum ≡ Q_51 L_1 spectrum: "
          f"{'YES' if plus_matches else 'NO'}")
else:
    print(f"  dim(V_+) = {len(spec_plus)} ≠ dim(L_1(Q_51)) = {len(spec51)}")
    plus_matches = False

print()

# --- Verdict ----------------------------------------------------------------
print("=" * 76)
print("Verdict: J-decomposition resolves Probe 5(ii)")
print("=" * 76)
print(f"  L_1(Q_102) splits into V_+ ⊕ V_- under J-action.")
print(f"  dim(V_+) = {len(pos)}, dim(V_-) = {len(neg)}, total = {n_e}.")
print()
print(f"  V_+ ≡ Q_51's L_1 spectrum: {'YES (cluster identity → spectral identity)' if plus_matches else 'NO'}")
print()
if plus_matches:
    print("  Refined scale-separation claim:")
    print("    - Sector-strip identity holds at cluster level (Probe 5(i))")
    print("    - J-symmetric block of L_1(Q_102) ≡ L_1(Q_51) spectrally (here)")
    print("    - J-antisymmetric block carries cross-sector mixing structure")
    print("  Together: Q_51 carries scale separation via the J-symmetric block of Q_102.")
