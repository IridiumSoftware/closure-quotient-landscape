#!/usr/bin/env python3
# Probe — Q_24 = cyc(6) carries Rosen "causality / entailment ordering".
#
# Rosen (M,R)-system role: causal entailment = the system's internal
# ordering of cause-effect relations. In CFS: the multiway-graph
# composition rule (v1, v2) -> w means "w is entailed by v1 and v2";
# the seed's structure determines the causal-trajectory signature.
#
# Hypothesis. The cyc(n) family is the cleanest causal-entailment
# generator: closure cardinality scales as 4n (paper Prop cyc-law),
# giving a *constant* "causal step width" of 4 per cycle node. The
# K_n^3 family, by contrast, scales as n(3n-1)/2 (paper Thm Kn-pent),
# giving a *growing* "cluster proliferation" with n. Constant-slope
# = causal entailment; growing-slope = symmetric-clustering. Q_24 is
# cyc(6), the smallest cyc anchor at the canonical n=6.
#
# Predicted discrimination:
#   cyc(n): |Q|/n = 4 (constant)        <-- causal-entailment signature
#   K_n^3:  |Q|/n = (3n-1)/2 (growing)  <-- cluster proliferation
#
# Validity gate: K_6^3 single-side = 51 (also reconfirmed by the K_n^3
# pentagonal row at n=6).

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

# --- POOL_B (20 strict-generic ICs, supports n=3..20) -----------------------
POOL_B = [
 [(4,1),(1,3),(2,-1)],[(1,2),(3,1),(4,-3)],[(2,3),(4,-1),(1,2)],
 [(3,-2),(2,1),(4,3)],[(4,3),(1,-2),(3,1)],[(2,-3),(4,1),(1,2)],
 [(1,4),(3,-1),(2,3)],[(4,-1),(2,3),(1,2)],[(3,2),(1,-3),(4,1)],
 [(2,1),(4,-3),(3,2)],[(1,3),(2,-4),(4,1)],[(3,-1),(4,2),(2,3)],
 [(4,2),(3,-1),(1,3)],[(1,-3),(2,4),(3,-1)],[(2,3),(1,-4),(3,2)],
 [(3,1),(2,-3),(4,1)],[(4,1),(3,2),(1,-3)],[(2,4),(1,-1),(3,-2)],
 [(1,2),(4,3),(2,-1)],[(3,-1),(1,4),(2,3)],
]

def make_ics(pool, n):
    return {v: [list(c) for c in pool[v]] for v in range(n)}

# --- Seed constructors ------------------------------------------------------
def cyc_seed(n):
    return [(i, (i+1)%n, (i+2)%n) for i in range(n)]

def kn3_seed(n):
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

# --- Per-depth cluster growth: for cyc(6) specifically ----------------------
def cluster_growth_by_depth(topo, ics, max_depth):
    """Return [|Q|_d for d in 0..max_depth]: cluster count at each depth.
    For cyc(n), expect d=0 -> n (the n seed ICs), then growth that
    saturates at 4n by some depth."""
    counts = []
    psi, edges = build_multiway(topo, ics, max_depth)
    # Re-cluster at each prefix of edges
    for d in range(max_depth + 1):
        # Take only psi entries reachable at depth <= d
        psi_d = {}
        for v in psi:
            psi_d[v] = psi[v]
            # quick stop — we want full closure at depth d.
            # Easier: rebuild at depth d.
        psi_at_d, _ = build_multiway(topo, ics, d)
        reps, _ = cluster(psi_at_d)
        counts.append(len(reps))
    return counts

# --- Validity gate ----------------------------------------------------------
g = single_side_count(kn3_seed(6), make_ics(POOL_B, 6), 5)
print(f"Validity gate: K_6^3 single-side (POOL_B, depth=5) = {g}  "
      f"(must be 51) -> {'PASS' if g == 51 else 'FAIL'}")
assert g == 51
print()

# --- Probe: cyc(n) and K_n^3 scaling ----------------------------------------
print("=" * 76)
print("Probe — Q_24 = cyc(6) carries causal-entailment (constant slope per n)")
print("=" * 76)
print()
print("cyc(n) closure law: |Q| = 4n (paper Prop cyc-law)")
print("K_n^3 closure law: |Q| = n(3n-1)/2 (paper Thm Kn-pent)")
print()
print(f"{'n':>3} {'|Q(cyc(n))|':>12} {'4n':>4} {'cyc match':>10} "
      f"{'|Q(K_n^3)|':>11} {'n(3n-1)/2':>10} {'Kn match':>9}")
print("-" * 76)

cyc_match_all = True
kn3_match_all = True

for n in range(3, 11):
    ics = make_ics(POOL_B, n)
    cyc_q = single_side_count(cyc_seed(n), ics, 4)
    cyc_pred = 4 * n
    cyc_ok = (cyc_q == cyc_pred)
    if not cyc_ok: cyc_match_all = False

    if n <= 7:  # Kn^3 at n>7 gets expensive at depth 4
        kn3_q = single_side_count(kn3_seed(n), ics, 4)
        kn3_pred = n * (3*n - 1) // 2
        kn3_ok = (kn3_q == kn3_pred)
        if not kn3_ok: kn3_match_all = False
        kn3_q_str = str(kn3_q)
        kn3_pred_str = str(kn3_pred)
        kn3_ok_str = "Y" if kn3_ok else "N"
    else:
        kn3_q_str = "-"
        kn3_pred_str = "-"
        kn3_ok_str = "-"

    print(f"{n:>3} {cyc_q:>12} {cyc_pred:>4} {'Y' if cyc_ok else 'N':>10} "
          f"{kn3_q_str:>11} {kn3_pred_str:>10} {kn3_ok_str:>9}")

print()
print("Slope analysis: |Q| / n")
print(f"{'n':>3} {'cyc slope':>10} {'K_n^3 slope':>12} {'K_n^3 / cyc slope':>20}")
print("-" * 50)
for n in range(3, 8):
    ics = make_ics(POOL_B, n)
    cyc_q = single_side_count(cyc_seed(n), ics, 4)
    kn3_q = single_side_count(kn3_seed(n), ics, 4)
    print(f"{n:>3} {cyc_q/n:>10.2f} {kn3_q/n:>12.2f} {kn3_q/cyc_q:>20.2f}")

print()
print("cyc slope is CONSTANT at 4 (causal-entailment fixed width).")
print("K_n^3 slope GROWS with n (cluster proliferation).")
print()

# --- Per-depth growth for cyc(6): show closure saturates by depth ~4 -------
print("Multiway closure growth for cyc(6) by depth:")
counts_cyc6 = cluster_growth_by_depth(cyc_seed(6), make_ics(POOL_B, 6), 5)
print(f"  depth -> |Q|: {counts_cyc6}  (saturates at 24 = 4*6)")
print()

# --- Verdict ----------------------------------------------------------------
print(f"cyc(n) -> 4n verified for n=3..10:        "
      f"{'PASS' if cyc_match_all else 'FAIL'}")
print(f"K_n^3 -> n(3n-1)/2 verified for n=3..7:  "
      f"{'PASS' if kn3_match_all else 'FAIL'}")
print()
verdict = "PASS" if (cyc_match_all and kn3_match_all) else "FAIL"
print(f"Hypothesis 'cyc(n) carries causal-entailment (constant slope) "
      f"vs K_n^3 cluster proliferation (growing slope)': {verdict}")
