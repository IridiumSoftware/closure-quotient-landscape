#!/usr/bin/env python3
# K_n^3 extended probe — pin the n>=9 behavior.
#
# Findings from Kn3_large_n_probe.py: K_n^3 single-side matches pentagonal
# for n=6..8, then drops below at n=9..12 (112/133/164/193 vs 117/145/176/210).
# Depth is not the issue — d=3..d=8 all give the same count.
#
# Two open questions:
#   (1) IC-independence — does a SECOND generic IC pool give the same n>=9
#       numbers, or does the count vary by IC pool? (If IC-dependent at
#       large n, generic-regime invariance BREAKS at n=9.)
#   (2) Extend to n=13..15 — what does the sequence do beyond n=12?
#
# Both probes use depth=5 (well into stable regime per main probe).
# Also: report C-closed counts at n=9..15 to check whether the |Q∪C(Q)|=2|Q|
# rule still holds (or whether large-n breaks that too).

import time, sys, random

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

# Pool A — the canonical 15 generic ICs from prime_count_probe.py.
POOL_A = [
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
 [(3,1),(2,0),(1,-1)],[(1,1),(3,-1),(2,0)],[(2,-1),(1,2),(3,1)],
 [(1,2),(3,1),(2,-1)],[(2,-1),(1,3),(1,1)],[(3,1),(2,-1),(1,2)],
 [(1,1),(2,2),(3,-1)],[(2,1),(1,-2),(2,1)],[(1,-1),(3,2),(1,1)],
]

# Pool B — independent set of 16 generic mixed real/imag triples,
# deliberately different small integer values and signs from Pool A.
# Each triple: 3 Gaussian-integer components, each component has nonzero
# real and nonzero imag part (strict genericity, avoids axis-aligned
# degeneracies).
POOL_B = [
 [(4,1),(1,3),(2,-1)],[(1,2),(3,1),(4,-3)],[(2,3),(4,-1),(1,2)],
 [(3,-2),(2,1),(4,3)],[(4,3),(1,-2),(3,1)],[(2,-3),(4,1),(1,2)],
 [(1,4),(3,-1),(2,3)],[(4,-1),(2,3),(1,2)],[(3,2),(1,-3),(4,1)],
 [(2,1),(4,-3),(3,2)],[(1,3),(2,-4),(4,1)],[(3,-1),(4,2),(2,3)],
 [(4,2),(3,-1),(1,3)],[(1,-3),(2,4),(3,-1)],[(2,3),(1,-4),(3,2)],
 [(3,1),(2,-3),(4,1)],
]

def make_ics(pool, n):
    assert n <= len(pool), f"need more ICs for n={n} (pool has {len(pool)})"
    return {v: [list(c) for c in pool[v]] for v in range(n)}

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
            w=cache[key]
            edges.append((d+1,w,v2,v3))
            edges.append((d+1,w,v1,v3))
            edges.append((d+1,w,v1,v2))
    return psi

def count_single(topo, pool, n_verts, depth):
    psi=build_multiway(topo, make_ics(pool, n_verts), depth)
    reps=[]
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv,r) for r in reps): reps.append(pv)
    return len(reps)

def count_c_closed(topo, pool, n_verts, depth):
    init_o = make_ics(pool, n_verts)
    init_c = {v:[gconj(x) for x in vec] for v,vec in init_o.items()}
    psi_o = build_multiway(topo, init_o, depth)
    psi_c = build_multiway(topo, init_c, depth)
    off = max(psi_o)+1
    allp = {}
    for v,p in psi_o.items(): allp[v]=p
    for v,p in psi_c.items(): allp[v+off]=p
    reps=[]
    for v in sorted(allp):
        pv=allp[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv,r) for r in reps): reps.append(pv)
    return len(reps)

def Kn(n):
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

def pentagonal(n): return n*(3*n-1)//2

# Validity gate (Pool A, n=6, must give 51)
gate = count_single(Kn(6), POOL_A, 6, depth=5)
print(f"validity gate (Pool A, K_6^3): {gate}  (must be 51) -> "
      f"{'PASS' if gate==51 else 'FAIL'}")
if gate != 51: sys.exit("gate failed")

# Validity gate (Pool B, n=6, must also give 51 — IC-independence at small n)
gate_b = count_single(Kn(6), POOL_B, 6, depth=5)
print(f"validity gate (Pool B, K_6^3): {gate_b}  (must be 51) -> "
      f"{'PASS' if gate_b==51 else 'FAIL — Pool B not generic at n=6'}")
if gate_b != 51: sys.exit("Pool B gate failed")

print()
print("Question 1 + 2: K_n^3 single-side, both pools, n=6..15, depth=5")
print(f"{'n':>3}  {'|K_n^3|':>8}  {'pent':>5}  {'Pool A':>7}  {'Pool B':>7}  "
      f"{'agree?':>7}  {'pent?':>7}  {'shortfall':>10}  {'t_A':>7}  {'t_B':>7}")
print("-"*100)

results = []
for n in range(6, 16):
    topo = Kn(n)
    pent = pentagonal(n)
    t0 = time.time()
    cA = count_single(topo, POOL_A, n, depth=5)
    tA = time.time() - t0
    t0 = time.time()
    cB = count_single(topo, POOL_B, n, depth=5)
    tB = time.time() - t0
    agree = "YES" if cA == cB else "NO"
    pent_ok = "YES" if cA == pent else "no"
    shortfall = pent - cA
    print(f"{n:>3}  {len(topo):>8}  {pent:>5}  {cA:>7}  {cB:>7}  "
          f"{agree:>7}  {pent_ok:>7}  {shortfall:>10}  {tA:>7.2f}  {tB:>7.2f}")
    results.append((n, cA, cB, pent))

# C-closed check — does |Q∪C(Q)| = 2|Q| still hold at large n?
print()
print("Question 3: C-closed at large n — does the 2|Q| rule still hold?")
print(f"{'n':>3}  {'single Pool A':>13}  {'C-closed Pool A':>15}  "
      f"{'2*single':>9}  {'matches?':>9}")
print("-"*70)
for n in range(6, 16):
    topo = Kn(n)
    cs = count_single(topo, POOL_A, n, depth=5)
    cc = count_c_closed(topo, POOL_A, n, depth=5)
    matches = "YES" if cc == 2*cs else "NO"
    print(f"{n:>3}  {cs:>13}  {cc:>15}  {2*cs:>9}  {matches:>9}")

# Summary line for the findings doc
print()
print("=" * 70)
print("Sequence summary (Pool A, single-side, K_n^3 for n=6..15):")
seq = [r[1] for r in results]
print(f"  {seq}")
print(f"Pentagonal predictions:")
print(f"  {[r[3] for r in results]}")
print(f"Pool A == Pool B for all n? {'YES' if all(r[1]==r[2] for r in results) else 'NO'}")
