#!/usr/bin/env python3
"""
probe_pool_a_defect.py — characterize the Pool A defect at K_n^3 for n>=9.

From audit_threebug_companion_v1.md §1.3 (closure-v5 v335): the
canonical s58 IC set (Pool A) under-counts K_n^3 vs pentagonal
by shortfall 5/12/12/17/19/19/19 at n=9..15 and acquires exactly 1
J-fixed cluster at n>=10. Plateau at +19/f=1 by n=13. The suspected
cause is Pool A's axis-aligned components (Gaussian-integer triples
with one or more components having zero imaginary part — e.g. (1,0),
(2,0), (3,0)). This probe characterizes the cause.

Diagnostics:
  T1 — Pool A's axis-aligned vs strict-generic split (inventory).
  T2 — Strict-generic subset of Pool A (ICs 8-14), padded with Pool B
       entries to reach n=9: does it restore pentagonal?
  T3 — Axis-aligned subset of Pool A (ICs 0-7), padded with Pool B to
       reach n=9: does it reproduce Pool A's defect?
  T4 — Pool A with axis-aligned components replaced by strict-generic
       in-place: does it restore pentagonal?
  T5 — Identify the J-fixed cluster representative at K_10^3 under
       Pool A (the one that appears at n=10, persists thereafter).
  T6 — Shuffled Pool A assignment: does the defect survive IC-index
       permutation? (Tests whether the defect is in the SET vs the
       per-vertex assignment.)

Output: prints findings to stdout; the matching findings doc is
Kn3_pool_a_defect_findings_v1.md (written separately).
"""

import sys, random

# Re-use exact arithmetic from the existing probe corpus.
def gmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c)
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
def gconj(x): return (x[0], -x[1])
def giszero(x): return x == (0, 0)
def cross3(a, b):
    return [gsub(gmul(a[1], b[2]), gmul(a[2], b[1])),
            gsub(gmul(a[2], b[0]), gmul(a[0], b[2])),
            gsub(gmul(a[0], b[1]), gmul(a[1], b[0]))]
def compose(a, b): return [gconj(x) for x in cross3(a, b)]
def proj_equiv(a, b): return all(giszero(x) for x in cross3(a, b))
def is_zero_vec(v): return all(giszero(x) for x in v)

POOL_A = [
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
 [(3,1),(2,0),(1,-1)],[(1,1),(3,-1),(2,0)],[(2,-1),(1,2),(3,1)],
 [(1,2),(3,1),(2,-1)],[(2,-1),(1,3),(1,1)],[(3,1),(2,-1),(1,2)],
 [(1,1),(2,2),(3,-1)],[(2,1),(1,-2),(2,1)],[(1,-1),(3,2),(1,1)],
]
POOL_B = [
 [(4,1),(1,3),(2,-1)],[(1,2),(3,1),(4,-3)],[(2,3),(4,-1),(1,2)],
 [(3,-2),(2,1),(4,3)],[(4,3),(1,-2),(3,1)],[(2,-3),(4,1),(1,2)],
 [(1,4),(3,-1),(2,3)],[(4,-1),(2,3),(1,2)],[(3,2),(1,-3),(4,1)],
 [(2,1),(4,-3),(3,2)],[(1,3),(2,-4),(4,1)],[(3,-1),(4,2),(2,3)],
 [(4,2),(3,-1),(1,3)],[(1,-3),(2,4),(3,-1)],[(2,3),(1,-4),(3,2)],
 [(3,1),(2,-3),(4,1)],
]

def Kn(n):
    return [(i, j, k) for i in range(n) for j in range(n) if j != i
            for k in range(n) if k != i and k != j]

def make_ics(pool, n):
    return {v: [list(c) for c in pool[v]] for v in range(n)}

def make_ics_from_list(triples):
    return {v: [list(c) for c in triples[v]] for v in range(len(triples))}

def build_multiway(topo, psi_init, depth):
    psi = {k: list(v) for k, v in psi_init.items()}
    nv = max(psi_init) + 1
    edges = [(0, s1, s2, s3) for (s1, s2, s3) in topo]
    cache = {}
    for d in range(depth):
        for (_, v1, v2, v3) in [e for e in edges if e[0] == d]:
            key = (v1, v2)
            if key not in cache:
                w = compose(psi[v1], psi[v2])
                if is_zero_vec(w): continue
                psi[nv] = w; cache[key] = nv; nv += 1
            w = cache[key]
            edges.append((d+1, w, v2, v3))
            edges.append((d+1, w, v1, v3))
            edges.append((d+1, w, v1, v2))
    return psi

def count_single(topo, ics, depth):
    psi = build_multiway(topo, ics, depth)
    reps = []
    for v in sorted(psi):
        pv = psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in reps): reps.append(pv)
    return len(reps), reps

def count_c_closed_with_fixed(topo, ics, depth):
    """Return (single, c_closed, f, j_fixed_reps)."""
    ics_o = ics
    ics_c = {v: [list(gconj(tuple(c))) for c in vec] for v, vec in ics.items()}
    psi_o = build_multiway(topo, ics_o, depth)
    psi_c = build_multiway(topo, ics_c, depth)
    # single side
    reps_o = []
    for v in sorted(psi_o):
        pv = psi_o[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in reps_o): reps_o.append(pv)
    single = len(reps_o)
    # combined
    off = max(psi_o) + 1
    allp = {}
    for v, p in psi_o.items(): allp[v] = p
    for v, p in psi_c.items(): allp[v + off] = p
    reps_c = []
    for v in sorted(allp):
        pv = allp[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in reps_c): reps_c.append(pv)
    c_closed = len(reps_c)
    f = 2 * single - c_closed
    # find J-fixed reps (those projectively equal to their own conjugate)
    j_fixed = [r for r in reps_c
               if proj_equiv(r, [gconj(tuple(c)) for c in r])]
    return single, c_closed, f, j_fixed

def axis_aligned_components(triple):
    """Return indices (0,1,2) where a component is axis-aligned
    (re == 0 or im == 0)."""
    out = []
    for i, c in enumerate(triple):
        if c[0] == 0 or c[1] == 0:
            out.append(i)
    return out

print("=" * 70)
print("PROBE: Pool A defect characterization at K_n^3 for n>=9")
print("=" * 70)
print()

# Validity gate (sanity).
gate, _ = count_single(Kn(6), make_ics(POOL_A, 6), depth=5)
assert gate == 51, f"validity gate failed: {gate}"
gate, _ = count_single(Kn(6), make_ics(POOL_B, 6), depth=5)
assert gate == 51, f"Pool B gate failed: {gate}"

# ---- T1: axis-aligned inventory ----------------------------------------
print("T1 — Pool A axis-aligned inventory")
print("-" * 70)
print(f"{'ic':>3}  {'triple':<32}  axis-aligned components")
strict_idx = []
axis_idx = []
for i, ic in enumerate(POOL_A):
    aa = axis_aligned_components(ic)
    if aa:
        axis_idx.append(i)
        print(f"{i:>3}  {str(ic):<32}  {aa}")
    else:
        strict_idx.append(i)
print()
print(f"  axis-aligned ICs  : {axis_idx}  (count={len(axis_idx)})")
print(f"  strict-generic ICs: {strict_idx}  (count={len(strict_idx)})")
print()

# ---- T2: strict-generic subset of Pool A, padded to n=9 ---------------
print("T2 — Strict-generic subset of Pool A (ICs 8..14), padded with Pool B")
print("-" * 70)
# ICs 8..14 are strict-generic per the T1 inventory below; verify.
strict_subset = [POOL_A[i] for i in strict_idx]
print(f"  using Pool A strict-generic ICs at indices {strict_idx}")
print(f"  padding with Pool B[0..{9 - len(strict_subset) - 1}] to reach n=9")
test_ics = strict_subset + POOL_B[: 9 - len(strict_subset)]
single, cc, f, jf = count_c_closed_with_fixed(
    Kn(9), make_ics_from_list(test_ics), depth=5)
print(f"  K_9^3: single={single} (pentagonal 117), C-closed={cc}, "
      f"J-fixed f={f}")
print(f"  -> matches pentagonal? {'YES' if single == 117 else 'NO'}")
print()

# ---- T3: axis-aligned subset of Pool A, padded to n=9 -----------------
print("T3 — Axis-aligned subset of Pool A (ICs 0..7), padded with Pool B")
print("-" * 70)
axis_subset = [POOL_A[i] for i in axis_idx]
print(f"  using Pool A axis-aligned ICs at indices {axis_idx}")
print(f"  padding with Pool B[0..{9 - len(axis_subset) - 1}] to reach n=9")
if 9 - len(axis_subset) > 0:
    test_ics = axis_subset + POOL_B[: 9 - len(axis_subset)]
else:
    test_ics = axis_subset[:9]
single, cc, f, jf = count_c_closed_with_fixed(
    Kn(9), make_ics_from_list(test_ics), depth=5)
print(f"  K_9^3: single={single} (Pool A canonical 112), C-closed={cc}, "
      f"J-fixed f={f}")
print(f"  -> reproduces Pool A defect? {'YES' if single < 117 else 'NO'}")
print()

# ---- T4: Pool A with axis-aligned components replaced in-place --------
print("T4 — Pool A with axis-aligned components replaced by strict-generic")
print("-" * 70)
# Replace each axis-aligned component with a strict-generic alternative
# from a small deterministic pool.
REPLACE_POOL = [(2,1),(1,2),(3,1),(1,3),(2,3),(3,2),(4,1),(1,4),(2,-1),(-1,2)]
def heal_triple(triple, used_replacements):
    healed = []
    for c in triple:
        if c[0] == 0 or c[1] == 0:
            # pick first replacement not already used in this triple
            for r in REPLACE_POOL:
                if r not in used_replacements and r not in healed:
                    healed.append(r)
                    used_replacements.append(r)
                    break
            else:
                healed.append(REPLACE_POOL[0])
        else:
            healed.append(c)
    return healed

pool_a_healed = []
for ic in POOL_A:
    pool_a_healed.append(heal_triple(ic, []))
# Verify all healed components are strict-generic
for i, ic in enumerate(pool_a_healed):
    assert axis_aligned_components(ic) == [], \
        f"healed IC {i} still has axis-aligned: {ic}"
# Validity gate the healed pool at n=6
single_6, _ = count_single(Kn(6), make_ics_from_list(pool_a_healed[:6]), depth=5)
print(f"  Pool A healed: n=6 validity gate = {single_6} (must be 51) -> "
      f"{'PASS' if single_6 == 51 else 'FAIL'}")
if single_6 == 51:
    single, cc, f, jf = count_c_closed_with_fixed(
        Kn(9), make_ics_from_list(pool_a_healed[:9]), depth=5)
    print(f"  K_9^3 under healed Pool A: single={single} (pentagonal 117), "
          f"C-closed={cc}, J-fixed f={f}")
    print(f"  -> matches pentagonal? {'YES' if single == 117 else 'NO'}")
print()

# ---- T5: J-fixed cluster identification at K_10^3 under Pool A --------
print("T5 — J-fixed cluster representative at K_10^3 under Pool A")
print("-" * 70)
single, cc, f, jf = count_c_closed_with_fixed(
    Kn(10), make_ics(POOL_A, 10), depth=5)
print(f"  K_10^3 under Pool A: single={single} (vs pentagonal 145), "
      f"C-closed={cc} (vs 2|Q|=2*single={2*single}), J-fixed f={f}")
print(f"  -> J-fixed clusters: {len(jf)}")
for i, r in enumerate(jf):
    print(f"    J-fixed[{i}]: {r}")
    # Characterize: what makes this rep projectively real?
    # A vec is projectively real iff there exists lambda in C* such that
    # lambda * r has only real components. Equivalently: the cross of r
    # with its conjugate is zero (which is what proj_equiv checks).
    # We can compute the "lambda" by finding the phase that aligns the
    # first nonzero component with the real axis.
    first_nonzero = next((c for c in r if c != (0, 0)), None)
    if first_nonzero:
        a, b = first_nonzero
        norm_sq = a*a + b*b
        # lambda = conj(c) / |c|^2 normalized; in exact arithmetic the
        # "projective representative scaled to real" requires rational
        # scaling; we report the un-scaled rep and note the phase.
        print(f"      first nonzero component: ({a},{b})  "
              f"(|c|^2 = {norm_sq})")
print()

# ---- T6: Shuffled Pool A assignment ----------------------------------
print("T6 — Shuffled Pool A assignment: does the defect persist under")
print("     IC-index permutation? (Tests whether defect is in the SET vs")
print("     the per-vertex assignment.)")
print("-" * 70)
random.seed(17)
defect_counts = []
n_trials = 5
for trial in range(n_trials):
    perm = list(range(15))
    random.shuffle(perm)
    shuffled = [POOL_A[i] for i in perm[:9]]
    single, cc, f, jf = count_c_closed_with_fixed(
        Kn(9), make_ics_from_list(shuffled), depth=5)
    defect_counts.append((single, f))
    print(f"  trial {trial+1}: perm[:9]={perm[:9]}  K_9^3 single={single}, "
          f"J-fixed f={f}")
all_same_single = len({c[0] for c in defect_counts}) == 1
all_below_pent = all(c[0] < 117 for c in defect_counts)
print(f"  all trials below pentagonal 117? {'YES' if all_below_pent else 'NO'}")
print(f"  all trials same single-count?    {'YES' if all_same_single else 'NO'}")
print()

print("=" * 70)
print("End of probe — see findings doc Kn3_pool_a_defect_findings_v1.md")
print("=" * 70)
