#!/usr/bin/env python3
# Probe — Q_12 = generic floor carries Rosen "constraint maintenance".
#
# Rosen (M,R)-system role: constraint maintenance = the system's
# structural constraints that prevent collapse to triviality. In CFS:
# under truly-generic ICs, every Rosen-closed ternary seed with >=3
# vertices closes to at least |Q| = 12 (paper Proposition 5.X Generic
# floor). Sub-12 closures (Q_3, Q_5, Q_9) arise only under non-generic
# ICs — i.e., the genericity-of-ICs IS the constraint that maintains
# the floor.
#
# Hypothesis. Generic ICs on small Rosen-closed seeds: |Q| >= 12.
# Real-valued (degenerate) ICs on the same seeds: |Q| can drop
# below 12, exhibiting Q_3, Q_5, or Q_9.
#
# Predicted discrimination:
#   Generic ICs (POOL_A, POOL_B): every small seed -> |Q| >= 12
#   Degenerate ICs (POOL_REAL):    same seeds -> some |Q| < 12
#
# Validity gate: K_6^3 C-closed = 102 under POOL_A.

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
POOL_REAL = [
 [(2,0),(1,0),(3,0)],[(1,0),(2,0),(4,0)],[(3,0),(1,0),(2,0)],
 [(1,0),(3,0),(2,0)],[(2,0),(3,0),(1,0)],[(3,0),(2,0),(1,0)],
]
# POOL_DEGEN_PARALLEL: two ICs are projectively parallel (the third
# is generic). Predicted to give a deeper sub-floor collapse than
# real-valued (J-collapse alone keeps |Q| at least at single-side).
POOL_DEGEN_PARALLEL = [
 [(1,1),(1,1),(1,1)],[(1,1),(1,1),(1,1)],[(1,1),(1,1),(1,1)],
 [(2,2),(2,2),(2,2)],[(2,2),(2,2),(2,2)],[(2,2),(2,2),(2,2)],
]

def make_ics(pool, n):
    return {v: [list(c) for c in pool[v]] for v in range(n)}

# --- Seed constructors ------------------------------------------------------
def cyc_seed(n):
    """cyc(n) per paper Def 2.X."""
    return [(i, (i+1)%n, (i+2)%n) for i in range(n)]

def kn3_seed(n):
    """K_n^3 — all distinct ordered triples."""
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

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

# --- Validity gate ----------------------------------------------------------
g = single_side_count(kn3_seed(6), make_ics(POOL_A, 6), 5)
# Use single-side gate; full C-closed gate is verified in Probe 1
print(f"Validity gate: K_6^3 single-side (POOL_A, depth=5) = {g}  "
      f"(must be 51) -> {'PASS' if g == 51 else 'FAIL'}")
assert g == 51
print()

# --- Probe ------------------------------------------------------------------
print("=" * 78)
print("Probe — Q_12 = generic floor (constraint maintenance)")
print("=" * 78)
print()
print("Under truly-generic ICs, every Rosen-closed ternary seed with")
print(">=3 vertices closes to at least |Q| = 12 (paper Prop 5.X).")
print("Sub-12 closures arise only under constraint-violating ICs.")
print()

seeds = [
    ("cyc(3)",         cyc_seed(3),    3),
    ("cyc(4)",         cyc_seed(4),    4),
    ("cyc(5)",         cyc_seed(5),    5),
    ("K_3^3 (perms)",  kn3_seed(3),    3),
    ("K_4^3",          kn3_seed(4),    4),
]

pools = [
    ("POOL_A   (canonical generic)",   POOL_A),
    ("POOL_B   (strict-generic)",       POOL_B),
    ("POOL_REAL (real-valued, degen)", POOL_REAL),
    ("POOL_DEGEN_PARALLEL (parallel)",  POOL_DEGEN_PARALLEL),
]

print(f"{'seed':<18}", end="")
for label, _ in pools:
    print(f"{label:>34}", end="")
print()
print("-" * (18 + 34 * len(pools)))

for sname, seed, n_v in seeds:
    print(f"{sname:<18}", end="")
    for label, pool in pools:
        ics = make_ics(pool, n_v)
        n_q = single_side_count(seed, ics, 4)
        marker = "" if n_q >= 12 else "  <<<<"
        cell = f"|Q|={n_q}{marker}"
        print(f"{cell:>34}", end="")
    print()

print()

# --- Verdict ----------------------------------------------------------------
generic_above_floor = True
degen_below_floor = False

for sname, seed, n_v in seeds:
    for label, pool in pools[:2]:  # POOL_A, POOL_B
        ics = make_ics(pool, n_v)
        n_q = single_side_count(seed, ics, 4)
        if n_q < 12:
            generic_above_floor = False
    for label, pool in pools[2:]:  # POOL_REAL, POOL_DEGEN_PARALLEL
        ics = make_ics(pool, n_v)
        n_q = single_side_count(seed, ics, 4)
        if 0 < n_q < 12:
            degen_below_floor = True

print(f"All generic-IC closures >= 12 (floor maintained): "
      f"{'YES' if generic_above_floor else 'NO'}")
print(f"Some degenerate-IC closure < 12 (constraint violated): "
      f"{'YES' if degen_below_floor else 'NO'}")
print()
verdict = "PASS" if (generic_above_floor and degen_below_floor) else "FAIL"
print(f"Hypothesis 'Q_12 = constraint-maintaining floor; "
      f"sub-12 = constraint-violation': {verdict}")
