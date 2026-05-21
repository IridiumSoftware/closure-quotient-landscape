#!/usr/bin/env python3
# K_n^3 large-n probe — pin the closure law beyond n=6.
#
# Question: under generic ICs, what is single_side(K_n^3) for n >= 7?
# closure_quotient_landscape_findings_v1.md observed counts 112/133/164/193
# at n=9..12 at depth=4, sitting below the pentagonal formula
# n(3n-1)/2 = {117, 145, 176, 210}. Two hypotheses:
#   (a) depth=4 under-closes at large n; higher depth restores pentagonal.
#   (b) the actual asymptotic genuinely differs from pentagonal at large n.
#
# Strategy: for each n in 6..12, run single_side(K_n^3) across depths
#   3, 4, 5, 6, 7, 8 and check for depth-stability. Report measured value
#   vs pentagonal prediction. Stop a row early once two consecutive depths
#   give identical counts (depth-stable).
#
# Exact Gaussian-integer arithmetic. Validity-gated on K_6^3 -> 51 (any
# depth >= ~3 should reach this). EXPLORATION + characterization.

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

# 15 hand-verified GENERIC ICs (same set as prime_count_probe.py — the s58
# 9-vector set + q102 IC-set-2 extension). All mixed real/imag — the basin
# that gives the canonical IC-independent counts.
_GEN=[
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
 [(3,1),(2,0),(1,-1)],[(1,1),(3,-1),(2,0)],[(2,-1),(1,2),(3,1)],
 [(1,2),(3,1),(2,-1)],[(2,-1),(1,3),(1,1)],[(3,1),(2,-1),(1,2)],
 [(1,1),(2,2),(3,-1)],[(2,1),(1,-2),(2,1)],[(1,-1),(3,2),(1,1)],
]
def make_ics(n):
    assert n <= len(_GEN), f"need more generic ICs for n={n} (have {len(_GEN)})"
    return {v: [list(c) for c in _GEN[v]] for v in range(n)}

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

def count_single(topo, n_verts, depth):
    psi=build_multiway(topo, make_ics(n_verts), depth)
    reps=[]
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv,r) for r in reps): reps.append(pv)
    return len(reps)

def Kn(n):
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

def pentagonal(n): return n*(3*n-1)//2

# validity gate
gate = count_single(Kn(6), 6, depth=5)
print(f"validity gate: K_6^3 single-side at d=5 = {gate}  "
      f"(must be 51) -> {'PASS' if gate==51 else 'FAIL'}")
if gate != 51:
    sys.exit("validity gate failed — scan void")

print()
print("K_n^3 single-side depth scan (exact, generic ICs)")
print(f"{'n':>3}  {'|K_n^3|':>8}  {'pent':>5}  {'d=3':>6} {'d=4':>6} {'d=5':>6} "
      f"{'d=6':>6} {'d=7':>6} {'d=8':>6}  {'verdict':<28}  {'time_s':>8}")
print("-"*120)

DEPTHS = [3,4,5,6,7,8]
for n in range(6, 13):
    topo=Kn(n)
    pent=pentagonal(n)
    counts={}; row_t0=time.time(); verdict=""
    last=None; stable_depth=None
    for d in DEPTHS:
        t0=time.time()
        try:
            c = count_single(topo, n, depth=d)
        except Exception as ex:
            counts[d] = f"ERR:{ex}"; break
        counts[d] = c
        if last is not None and c == last and stable_depth is None:
            stable_depth = d
        last = c
        elapsed = time.time()-t0
        # honest cutoff: if a single depth-step takes more than 5 minutes,
        # don't push further at this n.
        if elapsed > 300:
            verdict = f"timeout at d={d}"
            break
    if stable_depth is not None:
        final = counts[stable_depth]
        if final == pent:
            verdict = f"converges to pent at d={stable_depth}"
        else:
            verdict = f"converges to {final} (NOT pent {pent})"
    elif all(isinstance(counts.get(d),int) for d in DEPTHS):
        verdict = verdict or f"still growing through d={DEPTHS[-1]}"
    total_t = time.time()-row_t0
    cells = [f"{counts.get(d,'-'):>6}" for d in DEPTHS]
    print(f"{n:>3}  {len(topo):>8}  {pent:>5}  " + " ".join(cells) +
          f"  {verdict:<28}  {total_t:>8.2f}")

print()
print("KEY: pent = n(3n-1)/2 (the small-n closed form). 'converges to pent'")
print("means depth was the only issue. 'converges to X (NOT pent N)' means")
print("the asymptotic differs and needs its own characterization.")
