#!/usr/bin/env python3
# Probe — Q_45 -> Q_90 C-doubling carries Rosen "repair / efficient-cause".
#
# Rosen (M,R)-system role: the system produces its own efficient cause to
# restore closure after perturbation. In CFS: the charge-conjugation map
# (Q -> Q ∪ C(Q)) is the internal corrective that extends a "half-closure"
# (one IC sector) to a doubly-terminal closure, without external input.
# Operates only under "metabolically-sufficient" ICs (J fixed-point-free).
#
# Hypothesis. Under truly-generic ICs (J fpf, f=0), the C-doubling map
# fires cleanly: |Q ∪ C(Q)| = 2|Q|. Q_45 -> Q_90 is the canonical
# instance. Under degenerate ICs (real-valued: ψ̄ = ψ exactly), the
# repair operator collapses: |Q ∪ C(Q)| = |Q|, f = |Q|. The IC
# regime determines whether repair fires.
#
# Three IC regimes:
#   1. POOL_A (canonical s58 entries 0..5):   generic at this seed
#   2. POOL_B (strict-generic):                 sanity cross-check
#   3. POOL_REAL (purely real Gaussian ints):   degenerate control
#
# Predicted discrimination:
#   POOL_A      : |Q|=45, |Q∪C(Q)|=90, f=0   (clean doubling — repair fires)
#   POOL_B      : (some |Q|), |Q∪C(Q)|=2|Q|, f=0 (clean doubling at this seed)
#   POOL_REAL   : |Q|=X,    |Q∪C(Q)|=X,   f=X  (total collapse — repair fails)
#
# Validity gate: K_6^3 C-closed = 102 under POOL_A. Real ICs are
# expected to fail or shift the K_6^3 gate; this is the structural
# signature of the degenerate regime.

import itertools
import numpy as np

# --- Gaussian-integer arithmetic --------------------------------------------
def gmul(x,y):
    a,b=x; c,d=y; return (a*c-b*d, a*d+b*c)
def gsub(x,y): return (x[0]-y[0], x[1]-y[1])
def gconj(x): return (x[0], -x[1])
def giszero(x): return x==(0,0)
def cross3(a,b):
    return [gsub(gmul(a[1],b[2]),gmul(a[2],b[1])),
            gsub(gmul(a[2],b[0]),gmul(a[0],b[2])),
            gsub(gmul(a[0],b[1]),gmul(a[1],b[0]))]
def compose(a,b): return [gconj(x) for x in cross3(a,b)]
def proj_equiv(a,b): return all(giszero(x) for x in cross3(a,b))
def is_zero_vec(v): return all(giszero(x) for x in v)

# --- IC pools ---------------------------------------------------------------
POOL_A = [
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
]
POOL_B = [
 [(4,1),(1,3),(2,-1)],[(1,2),(3,1),(4,-3)],[(2,3),(4,-1),(1,2)],
 [(3,-2),(2,1),(4,3)],[(4,3),(1,-2),(3,1)],[(2,-3),(4,1),(1,2)],
]
# POOL_REAL: every component (a, 0) — purely real Gaussian integer. Under
# the involution J: ψ̄(v) = (a, -0) = (a, 0) = ψ(v) identically. Predicted
# behaviour: total collapse, f = |Q|.
POOL_REAL = [
 [(2,0),(1,0),(3,0)],[(1,0),(2,0),(4,0)],[(3,0),(1,0),(2,0)],
 [(1,0),(3,0),(2,0)],[(2,0),(3,0),(1,0)],[(3,0),(2,0),(1,0)],
]

def make_ics(pool, n):
    return {v: [list(c) for c in pool[v]] for v in range(n)}

# --- Seed constructors ------------------------------------------------------
def adjacent_edges_perm():
    s=set()
    for i in range(6):
        t=(i,(i+1)%6,(i+2)%6)
        for p in itertools.permutations(t): s.add(p)
    return sorted(s)

def kn3_seed(n):
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

def q45_seed_edges(seed=42):
    base=set(adjacent_edges_perm())
    pool=[e for e in kn3_seed(6) if e not in base]
    rng=np.random.default_rng(seed); rng.shuffle(pool)
    return list(base)+pool[:6]

# --- Multiway + cluster -----------------------------------------------------
def build_multiway(topo, psi_init, depth):
    psi={k:list(v) for k,v in psi_init.items()}
    nv=max(psi_init)+1
    edges=[(0,s1,s2,s3) for (s1,s2,s3) in topo]; cache={}
    for d in range(depth):
        for (_,v1,v2,v3) in [e for e in edges if e[0]==d]:
            key=(v1,v2)
            if key not in cache:
                w=compose(psi[v1],psi[v2])
                if is_zero_vec(w): continue
                psi[nv]=w; cache[key]=nv; nv+=1
            if key in cache:
                w=cache[key]
                edges.append((d+1,w,v2,v3))
                edges.append((d+1,w,v1,v3))
                edges.append((d+1,w,v1,v2))
    return psi, edges

def cluster(psi):
    reps=[]; vtc={}
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        m=-1
        for ci,r in enumerate(reps):
            if proj_equiv(pv,r): m=ci; break
        if m>=0: vtc[v]=m
        else: reps.append(pv); vtc[v]=len(reps)-1
    return reps,vtc

def single_side_count(topo, ics, depth):
    psi, _ = build_multiway(topo, ics, depth)
    reps, _ = cluster(psi)
    return len(reps)

def c_closed_count(topo, ics, depth):
    psi_o, _ = build_multiway(topo, ics, depth)
    psi_c, _ = build_multiway(
        topo,
        {v: [gconj(x) for x in vec] for v, vec in ics.items()},
        depth,
    )
    off = max(psi_o) + 1
    allp = dict(psi_o)
    for v, p in psi_c.items():
        allp[v + off] = p
    reps, _ = cluster(allp)
    return len(reps)

# --- Validity gate check (POOL_A) -------------------------------------------
g = c_closed_count(kn3_seed(6), make_ics(POOL_A, 6), 5)
print(f"Validity gate: K_6^3 C-closed (POOL_A, depth=5) = {g}  "
      f"(must be 102) -> {'PASS' if g == 102 else 'FAIL'}")
assert g == 102
print()

# --- Probe -------------------------------------------------------------------
print("=" * 80)
print("Probe — Q_45 -> Q_90 carries repair / efficient-cause production")
print("=" * 80)
print()
print("The C-doubling map Q -> Q ∪ C(Q) is the Rosen-style internal")
print("corrective. Under generic ICs (J fpf, f=0) it fires cleanly:")
print("|Q ∪ C(Q)| = 2|Q|. Under degenerate ICs (real-valued: ψ̄ ≡ ψ)")
print("it collapses: |Q ∪ C(Q)| = |Q|.")
print()

print(f"{'IC pool':<32} {'|Q_45|':>7} {'|Q∪C(Q)|':>10} {'2|Q|':>7} {'f':>4} {'verdict':>22}")
print("-" * 86)

# Reuse q45_seed_edges with the canonical numpy RNG seed 42
seed_q45 = q45_seed_edges(42)

rows = []
for (label, pool) in [("POOL_A (canonical s58 6-IC)", POOL_A),
                       ("POOL_B (strict-generic)", POOL_B),
                       ("POOL_REAL (purely real)", POOL_REAL)]:
    ics = make_ics(pool, 6)
    n_s = single_side_count(seed_q45, ics, 4)
    n_c = c_closed_count(seed_q45, ics, 4)
    f = 2 * n_s - n_c
    if f == 0:
        verdict = "CLEAN DOUBLING"
    elif f == n_s:
        verdict = "TOTAL COLLAPSE"
    elif 0 < f < n_s:
        verdict = f"PARTIAL COLLAPSE (f={f})"
    else:
        verdict = f"UNEXPECTED (f={f})"
    rows.append((label, n_s, n_c, 2*n_s, f, verdict))
    print(f"{label:<32} {n_s:>7} {n_c:>10} {2*n_s:>7} {f:>4} {verdict:>22}")

print()

# Verdict
generic_row = rows[0]  # POOL_A
real_row = rows[2]     # POOL_REAL

repair_fires_generic = (generic_row[4] == 0 and generic_row[1] == 45)
repair_fails_real = (real_row[4] == real_row[1])

print(f"Repair fires under POOL_A (Q_45 -> Q_90, f=0): "
      f"{'YES' if repair_fires_generic else 'NO'}")
print(f"Repair fails under POOL_REAL (total collapse, f=|Q|): "
      f"{'YES' if repair_fails_real else 'NO'}")
print()

verdict = "PASS" if (repair_fires_generic and repair_fails_real) else "FAIL"
print(f"Hypothesis 'C-doubling = Rosen-style internal repair "
      f"conditioned on generic ICs': {verdict}")
