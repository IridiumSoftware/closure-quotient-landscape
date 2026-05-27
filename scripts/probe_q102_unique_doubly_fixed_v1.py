#!/usr/bin/env python3
# Probe — Q_102 is uniquely the identity carrier in the canonical family.
#
# Rosen (M,R)-system role: identity = the recoverable structure that
# persists under admissible transport. In CFS: Q_102 is the unique
# canonical-family member satisfying ALL FIVE identity criteria
# simultaneously.
#
# Criteria:
#   (1) F-fixed       — closure-growth endofunctor idempotent
#   (2) C-fixed       — charge-conjugation idempotent
#   (3) =C(K_6^3)     — reached by one C-application from validity anchor Q_51
#   (4) Max-6v        — maximal cardinality on canonical 6-vertex seeds
#   (5) Strip→Q_51    — sector-strip (orig_only) recovers Q_51 cluster-by-cluster
#
# Predicted discrimination:
#    Q | F | C | C(K_6^3) | Max-6v | Strip→Q_51
#   Q_12   Y   N      N         N         N
#   Q_24   Y   N      N         N         N
#   Q_45   Y   N      N         N         N
#   Q_48   Y   Y      N         N         N   (strips to Q_24, not Q_51)
#   Q_51   Y   N      N         N         N   (no sector struct; self)
#   Q_84   Y   Y      N         N         N   (strips to ≠ Q_51)
#   Q_90   Y   Y      N         N         N   (strips to Q_45)
#   Q_102  Y   Y      Y         Y         Y   <-- UNIQUE
#
# Q_181 omitted: 12-vertex seed, fails Max-6v.

import os
import sys
from pathlib import Path

# Set CFS_REPO_ROOT to a Closure-Forces-Structure source checkout
# (https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs)
# before invoking this probe.
CLOSURE_V5 = os.environ.get("CFS_REPO_ROOT")
if not CLOSURE_V5:
    sys.exit("Set CFS_REPO_ROOT env var to a Closure-Forces-Structure checkout.")
sys.path.insert(0, CLOSURE_V5)

import itertools
import numpy as np

from q102_build_exact_v1 import (
    build_c_closed_quotient, complete_ternary, adjacent_ternary,
    gaussian_ics, build_multiway, proj_equiv, is_zero_vec, _gconj,
)


# --- Seed constructors ------------------------------------------------------
def cyc_seed(n):
    return [(i, (i+1)%n, (i+2)%n) for i in range(n)]

def kn3_seed(n):
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

def adjacent_edges_perm():
    s=set()
    for i in range(6):
        t=(i,(i+1)%6,(i+2)%6)
        for p in itertools.permutations(t): s.add(p)
    return sorted(s)

def q45_seed_edges(seed=42):
    base=set(adjacent_edges_perm())
    pool=[e for e in kn3_seed(6) if e not in base]
    rng=np.random.default_rng(seed); rng.shuffle(pool)
    return list(base)+pool[:6]


# --- Single-side closure -----------------------------------------------------
def single_side_clusters(seed_edges, n_verts, depth=4, ic_seed=None):
    """Return list of exact cluster representatives, single-side."""
    ics = gaussian_ics(n=n_verts, ic_seed=ic_seed)
    psi, _, _ = build_multiway(seed_edges, ics, depth)
    reps = []
    for v in sorted(psi.keys()):
        pv = psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in reps):
            reps.append(pv)
    return reps


# --- C-closed quotient (returns dict with cl_origin etc.) -------------------
def c_closed_full(seed_edges, n_verts, depth=4):
    return build_c_closed_quotient(seed_edges, depth=depth, n_seed_verts=n_verts)


# --- F-fix check via running multiway on the cluster reps as a fresh seed ---
def F_closure_size_from_reps(reps, q_he, depth):
    """|F(Q)|: re-run multiway from Q's own reps + q_he topology."""
    seed = {i: list(reps[i]) for i in range(len(reps))}
    psi, _, _ = build_multiway(list(q_he), seed, depth)
    new_reps = []
    for v in sorted(psi.keys()):
        pv = psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in new_reps):
            new_reps.append(pv)
    return len(new_reps)


# --- Per-Q evaluation --------------------------------------------------------
def evaluate_Q(name, mode, seed, n_verts, q51_reps):
    """Compute the 5 criteria for a single Q.

    mode: 'single' = single-side closure of `seed`
          'cclosed' = C-closed closure of `seed`
    seed: the seed edge list
    n_verts: number of vertices in the seed
    q51_reps: 51 cluster reps of Q_51 (for criterion 5 comparison)
    """
    if mode == 'single':
        reps = single_side_clusters(seed, n_verts, depth=4)
        n = len(reps)
        # Build q_he from these clusters: rebuild and quotient
        ics = gaussian_ics(n=n_verts, ic_seed=None)
        psi, _, edges = build_multiway(seed, ics, 4)
        # Cluster psi
        v_to_c = {}
        for v in sorted(psi.keys()):
            pv = psi[v]
            if is_zero_vec(pv): continue
            for ci, r in enumerate(reps):
                if proj_equiv(pv, r):
                    v_to_c[v] = ci; break
        q_he = set()
        for (_, v1, v2, v3) in edges:
            if v1 in v_to_c and v2 in v_to_c and v3 in v_to_c:
                # Preserve raw ordering (unsorted) — matches q102_build_exact's q_he convention
                q_he.add((v_to_c[v1], v_to_c[v2], v_to_c[v3]))
        # C-fix: |Q ∪ C(Q)|
        cQ = c_closed_full(seed, n_verts)
        c_size = cQ['n_cl']
        c_fix = (c_size == n)
        # Strip→Q_51: orig_only sector ≡ Q_51?
        orig_only_count = sum(1 for c in cQ['cl_origin'].values() if c == 'orig_only')
        strip_to_q51 = False  # single-side has no sector structure
        is_C_K6 = False
        # F-fix
        f_size = F_closure_size_from_reps(reps, q_he, 5)
        f_fix = (f_size == n)

    else:  # 'cclosed'
        cQ = c_closed_full(seed, n_verts)
        n = cQ['n_cl']
        reps_dict = cQ['q_psi_exact']
        reps = [reps_dict[i] for i in range(n)]
        q_he = cQ['q_he']
        # F-fix: re-multiway on the reps + topology
        f_size = F_closure_size_from_reps(reps, q_he, 5)
        f_fix = (f_size == n)
        # C-fix: c_closed of an already-c-closed quotient — by construction
        # we can re-c-close the cluster reps as seed, but easier: C-closure
        # of a c-closed quotient is c-idempotent for the closed family
        # (paper Prop C-idempotent). Verify empirically by re-c-closing the
        # cluster representatives' implied topology.
        # We'll trust the structural property: c_closed(seed) is c-idempotent
        # iff |c_closed(q_he, ic_from_reps)| == n. For brevity:
        # Use a quick C-iteration check: |Q ∪ C(Q)| where C is applied to
        # the cluster reps directly.
        cc_size = c_iterate_from_reps(reps, q_he)
        c_fix = (cc_size == n)
        # Strip→Q_51: orig_only reps proj-equiv to Q_51's reps?
        orig_only_idxs = sorted(c for c in range(n) if cQ['cl_origin'][c] == 'orig_only')
        orig_reps = [reps[c] for c in orig_only_idxs]
        if len(orig_reps) == 51 == len(q51_reps):
            matched = sum(1 for r51 in q51_reps
                          if any(proj_equiv(r51, ro) for ro in orig_reps))
            strip_to_q51 = (matched == 51)
        else:
            strip_to_q51 = False
        # is C(K_6^3)?
        is_C_K6 = (sorted(seed) == sorted(kn3_seed(6)))

    return {
        'name': name, 'n': n,
        'F': f_fix, 'C': c_fix,
        'is_C_K6': is_C_K6, 'max_6v': False,  # filled after all evaluated
        'strip': strip_to_q51,
    }


def c_iterate_from_reps(reps, he):
    """|Q ∪ C(Q)| via one-layer expansion from cluster reps + q_he."""
    qp = {i: list(reps[i]) for i in range(len(reps))}
    def expand(seed):
        psi = dict(seed); nv = max(psi)+1
        for (s1, s2, _t) in he:
            if s1 not in psi or s2 not in psi: continue
            # use compose from q102_build_exact (via build_multiway helper)
            # — but compose isn't directly importable, so do it inline:
            from q102_build_exact_v1 import compose_exact
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


# --- Main pipeline ----------------------------------------------------------
print("=" * 78)
print("Probe — Q_102 is unique identity carrier in the canonical family")
print("=" * 78)
print()

# Build Q_51 reps for sector-strip comparison
q51_reps = single_side_clusters(complete_ternary(6), 6, depth=4)
assert len(q51_reps) == 51, f"Q_51 build failed: {len(q51_reps)}"
print(f"Q_51 reference: 51 cluster reps built.")
print()

# Family entries
family = [
    ("Q_12",  'single',  cyc_seed(3),            3),
    ("Q_24",  'single',  cyc_seed(6),            6),
    ("Q_45",  'single',  q45_seed_edges(42),     6),
    ("Q_48",  'cclosed', cyc_seed(6),            6),
    ("Q_51",  'single',  complete_ternary(6),    6),
    ("Q_84",  'cclosed', adjacent_edges_perm(),  6),
    ("Q_90",  'cclosed', q45_seed_edges(42),     6),
    ("Q_102", 'cclosed', complete_ternary(6),    6),
]

results = []
for name, mode, seed, n_verts in family:
    r = evaluate_Q(name, mode, seed, n_verts, q51_reps)
    results.append(r)
    print(f"  evaluated {name} (mode={mode}): |Q|={r['n']}, "
          f"F={r['F']}, C={r['C']}, is_C(K6³)={r['is_C_K6']}, strip→Q_51={r['strip']}")

# Fill in Max-6v (largest cardinality in canonical 6-vertex-seed family)
six_v = [r for r in results if r['n']]
max_n = max(r['n'] for r in six_v)
for r in results:
    r['max_6v'] = (r['n'] == max_n)

# --- Discrimination table ---------------------------------------------------
print()
print("Identity-carrier discrimination table:")
print()
print(f"{'Q':>7} {'|Q|':>5} {'(1) F-fix':>10} {'(2) C-fix':>10} "
      f"{'(3) C(K6³)':>11} {'(4) Max-6v':>11} {'(5) Strip→Q_51':>15} {'ALL 5':>6}")
print("-" * 80)
for r in results:
    all_five = r['F'] and r['C'] and r['is_C_K6'] and r['max_6v'] and r['strip']
    print(f"{r['name']:>7} {r['n']:>5} "
          f"{'Y' if r['F'] else 'N':>10} "
          f"{'Y' if r['C'] else 'N':>10} "
          f"{'Y' if r['is_C_K6'] else 'N':>11} "
          f"{'Y' if r['max_6v'] else 'N':>11} "
          f"{'Y' if r['strip'] else 'N':>15} "
          f"{'YES' if all_five else '':>6}")
print()

unique_carriers = [r for r in results if (r['F'] and r['C'] and r['is_C_K6']
                                            and r['max_6v'] and r['strip'])]
print(f"Members satisfying all 5 criteria: "
      f"{[r['name'] for r in unique_carriers]}")
if len(unique_carriers) == 1 and unique_carriers[0]['name'] == 'Q_102':
    print(f"\nHypothesis 'Q_102 uniquely carries identity (5-criterion)': PASS")
elif len(unique_carriers) == 1:
    print(f"\nHypothesis 'Q_102 uniquely carries identity': "
          f"UNEXPECTED — unique row is {unique_carriers[0]['name']}, not Q_102")
elif len(unique_carriers) == 0:
    print(f"\nHypothesis 'Q_102 uniquely carries identity': "
          f"FAIL — no member satisfies all 5")
else:
    print(f"\nHypothesis 'Q_102 UNIQUELY carries identity': "
          f"FAIL — {len(unique_carriers)} members satisfy all 5")
