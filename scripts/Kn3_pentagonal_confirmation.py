#!/usr/bin/env python3
# K_n^3 pentagonal-law confirmation — third independent IC pool + range push.
#
# Prior finding (Kn3_extended_probe.py): under Pool A (the canonical s58
# 9-vector set + q102 IC-set-2), K_n^3 single-side drops below pentagonal
# at n>=9; under Pool B (strict-genericity, every component nonzero real
# AND nonzero imag), it matches pentagonal at n=6..15.
#
# To make "K_n^3 -> n(3n-1)/2 for all n" a paper-level claim, we need:
#   (a) an independent third pool (Pool C) agreeing with Pool B, ruling
#       out a Pool-B fluke
#   (b) extension beyond n=15 to verify the law continues
#   (c) a J-fixed-point report per pool — Pool A acquires exactly one J
#       fixed point at n>=10 (a separate sharp finding about Pool A,
#       not about K_n^3)
#
# Pool C: a fresh strict-generic IC set with NO numerical overlap with
# Pool B (different bounded-integer ranges and irregular sign patterns).

import time, sys

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

# Pool A — canonical s58 set extended (15 entries, used in every prior
# CFS probe). NOT FULLY GENERIC at n>=9 (this paper's finding).
POOL_A = [
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
 [(3,1),(2,0),(1,-1)],[(1,1),(3,-1),(2,0)],[(2,-1),(1,2),(3,1)],
 [(1,2),(3,1),(2,-1)],[(2,-1),(1,3),(1,1)],[(3,1),(2,-1),(1,2)],
 [(1,1),(2,2),(3,-1)],[(2,1),(1,-2),(2,1)],[(1,-1),(3,2),(1,1)],
]

# Pool B — strict-generic, mixed real/imag everywhere. Coefficients in
# {-4,-3,-2,-1,1,2,3,4}, every component has both parts nonzero.
POOL_B = [
 [(4,1),(1,3),(2,-1)],[(1,2),(3,1),(4,-3)],[(2,3),(4,-1),(1,2)],
 [(3,-2),(2,1),(4,3)],[(4,3),(1,-2),(3,1)],[(2,-3),(4,1),(1,2)],
 [(1,4),(3,-1),(2,3)],[(4,-1),(2,3),(1,2)],[(3,2),(1,-3),(4,1)],
 [(2,1),(4,-3),(3,2)],[(1,3),(2,-4),(4,1)],[(3,-1),(4,2),(2,3)],
 [(4,2),(3,-1),(1,3)],[(1,-3),(2,4),(3,-1)],[(2,3),(1,-4),(3,2)],
 [(3,1),(2,-3),(4,1)],[(4,1),(3,2),(1,-3)],[(2,4),(1,-1),(3,-2)],
 [(1,2),(4,3),(2,-1)],[(3,-1),(1,4),(2,3)],
]

# Pool C — independent strict-generic set with DIFFERENT coefficient
# range {-7,..,-3,3,..,7} (no overlap with Pool B's {-4,..,4} except by
# coincidence), and deliberately different sign patterns.
POOL_C = [
 [(5,3),(7,-4),(3,5)],[(6,-5),(4,3),(7,-3)],[(3,-4),(7,5),(5,3)],
 [(7,3),(3,-5),(4,7)],[(4,-3),(6,5),(7,-3)],[(5,-7),(3,4),(6,3)],
 [(7,-3),(5,4),(3,-7)],[(3,5),(7,-3),(4,6)],[(5,4),(7,-3),(6,5)],
 [(7,-5),(3,4),(6,-3)],[(4,7),(5,-3),(3,5)],[(6,3),(7,-5),(4,3)],
 [(3,7),(5,-4),(7,3)],[(7,3),(4,-7),(5,4)],[(3,-5),(7,4),(6,-3)],
 [(5,7),(4,-3),(6,5)],[(7,-4),(3,6),(5,-3)],[(4,5),(7,-6),(3,7)],
 [(6,-3),(5,7),(7,-4)],[(3,4),(5,-7),(6,3)],
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

def count_c_closed_and_f(topo, pool, n_verts, depth):
    """Return (single_side, c_closed, j_fixed_count f).
    Uses the identity |Q∪C(Q)| = 2*single - f.
    """
    init_o = make_ics(pool, n_verts)
    init_c = {v:[gconj(x) for x in vec] for v,vec in init_o.items()}
    psi_o = build_multiway(topo, init_o, depth)
    psi_c = build_multiway(topo, init_c, depth)
    # single-side count from psi_o alone
    reps_o = []
    for v in sorted(psi_o):
        pv = psi_o[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in reps_o): reps_o.append(pv)
    single = len(reps_o)
    # combined
    off = max(psi_o)+1
    allp = {}
    for v,p in psi_o.items(): allp[v]=p
    for v,p in psi_c.items(): allp[v+off]=p
    reps_c = []
    for v in sorted(allp):
        pv = allp[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv, r) for r in reps_c): reps_c.append(pv)
    c_closed = len(reps_c)
    f = 2*single - c_closed
    return single, c_closed, f

def Kn(n):
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

def pentagonal(n): return n*(3*n-1)//2

# Validity gates at n=6 (must give 51 for any generic pool)
for name, pool in (("A", POOL_A), ("B", POOL_B), ("C", POOL_C)):
    g = count_single(Kn(6), pool, 6, depth=5)
    status = "PASS" if g == 51 else "FAIL"
    print(f"validity gate (Pool {name}, K_6^3): {g}  -> {status}")
    if g != 51: sys.exit(f"Pool {name} gate failed")
print()

print("K_n^3: single-side counts and J fixed-point counts, all three pools")
print(f"{'n':>3}  {'pent':>5}  {'A_ss':>5}  {'A_f':>4}  {'B_ss':>5}  {'B_f':>4}  "
      f"{'C_ss':>5}  {'C_f':>4}  {'B==C':>5}  {'B==pent':>8}  {'C==pent':>8}")
print("-"*100)

max_n = min(len(POOL_A), len(POOL_B), len(POOL_C), 20)
results = []
for n in range(6, max_n+1):
    topo = Kn(n)
    pent = pentagonal(n)
    # Pool A
    if n <= len(POOL_A):
        a_ss, a_cc, a_f = count_c_closed_and_f(topo, POOL_A, n, 5)
        a_ss_s, a_f_s = str(a_ss), str(a_f)
    else:
        a_ss, a_f = None, None
        a_ss_s, a_f_s = "n/a", "n/a"
    # Pool B
    b_ss, b_cc, b_f = count_c_closed_and_f(topo, POOL_B, n, 5)
    # Pool C
    c_ss, c_cc, c_f = count_c_closed_and_f(topo, POOL_C, n, 5)
    bc_eq = "YES" if b_ss == c_ss else "NO"
    b_pent = "YES" if b_ss == pent else "no"
    c_pent = "YES" if c_ss == pent else "no"
    print(f"{n:>3}  {pent:>5}  {a_ss_s:>5}  {a_f_s:>4}  {b_ss:>5}  {b_f:>4}  "
          f"{c_ss:>5}  {c_f:>4}  {bc_eq:>5}  {b_pent:>8}  {c_pent:>8}")
    results.append((n, pent, a_ss, a_f, b_ss, b_f, c_ss, c_f))

print()
print("=" * 100)
all_b_pent = all(r[4] == r[1] for r in results)
all_c_pent = all(r[6] == r[1] for r in results)
all_bc_agree = all(r[4] == r[6] for r in results)
all_bf_zero = all(r[5] == 0 for r in results)
all_cf_zero = all(r[7] == 0 for r in results)
print(f"Pool B single-side = pentagonal for all n=6..{max_n}? {all_b_pent}")
print(f"Pool C single-side = pentagonal for all n=6..{max_n}? {all_c_pent}")
print(f"Pool B == Pool C for all n? {all_bc_agree}")
print(f"Pool B J-fixed = 0 for all n? {all_bf_zero}")
print(f"Pool C J-fixed = 0 for all n? {all_cf_zero}")
print()
print("Pool A degeneracy summary (shortfall from pentagonal, then J-fixed):")
for r in results:
    if r[2] is None: continue
    print(f"  n={r[0]:>2}: shortfall={r[1]-r[2]:>3}, J-fixed f={r[3]}")
