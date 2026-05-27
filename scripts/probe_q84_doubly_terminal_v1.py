#!/usr/bin/env python3
# Probe — Q_84 is doubly-terminal (F-fixed AND C-fixed), §6 omission?
#
# Paper §6 currently enumerates Q_48, Q_90, Q_102 as the doubly-terminal
# members of the canonical family. Probe 1 (probe_boundary_regulation_v1)
# reported Q_84 also doubly-terminal, which would mean §6's enumeration
# is incomplete. This stand-alone probe verifies Q_84's status using:
#
#   (a) The canonical closure-v5 build_c_closed_quotient (NOT probe-local
#       helpers — the authoritative builder used in reground_q_family.py
#       and the paper proper).
#   (b) F-fix at depths 4, 5, 6, 7 (extended range vs Probe 1).
#   (c) C-fix via iterated C-closure of Q_84's cluster reps.
#   (d) Cross-check against Q_90's analogous test (a known doubly-terminal).
#
# If Q_84 is doubly-terminal at all depths and matches Q_90's pattern,
# §6 enumeration should be updated to include Q_84 in v2.
#
# Note: Paper §5 Depth-stability already includes Q_84 in its list (line
# 1110 of the .tex), so the depth-stability claim is established — the
# question here is specifically whether Q_84 is ALSO F-fixed in the
# stronger sense of Theorem double-terminal.

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

import itertools

from q102_build_exact_v1 import (
    build_c_closed_quotient, build_multiway, adjacent_ternary, complete_ternary,
    gaussian_ics, proj_equiv, is_zero_vec, compose_exact, _gconj,
)


def q45_seed_edges(seed=42):
    """Q_90's seed: adjacent + 6 random K_6^3 pads (numpy RNG 42)."""
    import numpy as np
    base = set(adjacent_ternary())
    pool = [e for e in complete_ternary(6) if e not in base]
    rng = np.random.default_rng(seed); rng.shuffle(pool)
    return list(base) + pool[:6]


# --- F-closure (multiway from cluster reps + q_he topology) ----------------
def F_closure_size(reps, q_he, depth):
    seed = {i: [list(c) for c in reps[i]] for i in range(len(reps))}
    psi, _, _ = build_multiway(list(q_he), seed, depth)
    new_reps = []
    for v in sorted(psi.keys()):
        pv = psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in new_reps):
            new_reps.append(pv)
    return len(new_reps)


# --- C-iterate (one-layer conjugate union from cluster reps + q_he) --------
def C_iterate_size(reps, q_he):
    qp = {i: [list(c) for c in reps[i]] for i in range(len(reps))}
    def expand(seed):
        psi = dict(seed); nv = max(psi) + 1
        for (s1, s2, _t) in q_he:
            if s1 not in psi or s2 not in psi: continue
            w = compose_exact(psi[s1], psi[s2])
            if is_zero_vec(w): continue
            psi[nv] = w; nv += 1
        return psi
    po = expand(dict(qp))
    pc = expand({v: [_gconj(x) for x in qp[v]] for v in qp})
    off = max(po) + 1
    allp = dict(po)
    for v, p in pc.items():
        allp[v + off] = p
    new_reps = []
    for v in sorted(allp.keys()):
        pv = allp[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in new_reps):
            new_reps.append(pv)
    return len(new_reps)


# --- Test a single C-closed quotient at multiple depths ---------------------
def test_doubly_terminal(name, seed_edges, n_verts, expected_n):
    print(f"\n=== {name} (seed: {len(seed_edges)} edges on {n_verts} verts) ===")
    Q = build_c_closed_quotient(seed_edges, depth=4, n_seed_verts=n_verts)
    n = Q['n_cl']
    print(f"  |{name}| = {n}  (expected {expected_n}) -> "
          f"{'OK' if n == expected_n else 'MISMATCH'}")
    assert n == expected_n
    reps_dict = Q['q_psi_exact']
    reps = [reps_dict[i] for i in range(n)]
    q_he = Q['q_he']
    print(f"  cell complex: {len(q_he)} hyperedge triples in q_he")

    # F-fix across depths
    print(f"  F-closure across depths:")
    f_all_match = True
    for d in (4, 5, 6, 7):
        sz = F_closure_size(reps, q_he, d)
        ok = (sz == n)
        if not ok: f_all_match = False
        print(f"    depth={d}: |F(Q)|={sz}  {'(F-fixed)' if ok else '(NOT F-fixed)'}")

    # C-fix via iterated C-closure
    cc_size = C_iterate_size(reps, q_he)
    c_fixed = (cc_size == n)
    print(f"  C-iterate: |Q ∪ C(Q)|={cc_size}  "
          f"{'(C-idempotent)' if c_fixed else '(NOT C-idempotent)'}")

    return f_all_match and c_fixed


print("=" * 76)
print("Probe — Q_84 doubly-terminal stand-alone verification")
print("=" * 76)
print("""
§6 of the paper (line 1127) lists Q_48, Q_90, Q_102 as doubly-terminal.
Probe 1 reported Q_84 also doubly-terminal. This probe re-verifies with
the canonical closure-v5 builder + extended depth range (d=4..7).
""")

q84_ok = test_doubly_terminal("Q_84", adjacent_ternary(), 6, 84)
q90_ok = test_doubly_terminal("Q_90", q45_seed_edges(42), 6, 90)
q48_ok = test_doubly_terminal("Q_48", [(i,(i+1)%6,(i+2)%6) for i in range(6)], 6, 48)
q102_ok = test_doubly_terminal("Q_102", complete_ternary(6), 6, 102)

print("\n" + "=" * 76)
print("Stand-alone verification verdict:")
print(f"  Q_48  doubly-terminal: {'PASS' if q48_ok else 'FAIL'}")
print(f"  Q_84  doubly-terminal: {'PASS' if q84_ok else 'FAIL'}   <-- omitted from §6")
print(f"  Q_90  doubly-terminal: {'PASS' if q90_ok else 'FAIL'}")
print(f"  Q_102 doubly-terminal: {'PASS' if q102_ok else 'FAIL'}")
print()
if q84_ok and q48_ok and q90_ok and q102_ok:
    print("Conclusion: §6 should include Q_84 in Theorem double-terminal.")
    print("The four doubly-terminal members on 6-vertex seeds are:")
    print("  Q_48  = C-closed(cyc(6))         smallest")
    print("  Q_84  = C-closed(adjacent-36)    new in v2")
    print("  Q_90  = C-closed(q45-seed)       (adjacent + 6 pads)")
    print("  Q_102 = C-closed(K_6^3)          maximal")
elif not q84_ok:
    print("Probe 1's Q_84 finding does NOT survive stand-alone verification.")
    print("Investigate: probe-local helpers may have given a false positive.")
