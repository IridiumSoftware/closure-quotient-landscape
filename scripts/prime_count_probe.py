#!/usr/bin/env python3
# Thread 3 — prime-count probe. Can a CFS closure produce a prime
# cardinality (specifically 181), or only composites?
# Regular families are composite by derivation (see report); this scans
# IRREGULAR single-side seeds. Exact Gaussian-integer arithmetic.

import random, math

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
def isprime(n):
    if n<2: return False
    for p in range(2,int(math.isqrt(n))+1):
        if n%p==0: return False
    return True

# 15 hand-verified GENERIC ICs (the s58 9-vector set + q102 IC-set-2),
# all mixed real/imag — the basin that gives the canonical IC-independent
# counts (K6^3 -> 51 single-side). Validity-gated below.
_GEN=[
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
 [(3,1),(2,0),(1,-1)],[(1,1),(3,-1),(2,0)],[(2,-1),(1,2),(3,1)],
 [(1,2),(3,1),(2,-1)],[(2,-1),(1,3),(1,1)],[(3,1),(2,-1),(1,2)],
 [(1,1),(2,2),(3,-1)],[(2,1),(1,-2),(2,1)],[(1,-1),(3,2),(1,1)],
]
def make_ics(n):
    assert n <= len(_GEN), f"need more generic ICs for n={n}"
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
            edges.append((d+1,w,v2,v3)); edges.append((d+1,w,v1,v3)); edges.append((d+1,w,v1,v2))
    return psi

def count_single(topo, n_verts, depth=4):
    psi=build_multiway(topo, make_ics(n_verts), depth)
    reps=[]
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv,r) for r in reps): reps.append(pv)
    return len(reps)

def Kn(n): return [(i,j,k) for i in range(n) for j in range(n) if j!=i
                   for k in range(n) if k!=i and k!=j]

random.seed(17)
print("THREAD 3 — prime-count probe (irregular single-side closure counts)")
print("regular families: cyc(n)=4n, K_n^3=n(3n-1)/2 — composite for n>=3 (derived)")
# VALIDITY GATE — generic ICs must reproduce the known anchor K6^3 -> 51.
gate = count_single(Kn(6), 6, depth=5)
print(f"validity gate: K6^3 single-side = {gate}  (must be 51) -> "
      f"{'PASS' if gate==51 else 'FAIL — ICs not generic, scan void'}")
if gate != 51:
    raise SystemExit("validity gate failed")
print()
all_counts=set(); primes_found=[]
for n in range(3,13):
    full=Kn(n); maxe=len(full)
    counts=set()
    # sweep seed densities: fractions of K_n^3
    fracs=[0.15,0.3,0.45,0.6,0.75,0.9,1.0]
    for f in fracs:
        k=max(1,int(maxe*f))
        for trial in range(2):
            topo = full if f>=1.0 else random.sample(full,k)
            try:
                c=count_single(topo,n,depth=4)
                counts.add(c)
            except Exception as ex:
                pass
            if f>=1.0: break
    cl=sorted(counts)
    pr=[c for c in cl if isprime(c)]
    all_counts |= counts
    primes_found += pr
    flag = "  <-- 181!" if 181 in counts else ""
    print(f"n={n:>2}  K_n^3 edges={maxe:>4}  counts seen: {cl}  primes:{pr}{flag}")

print()
print(f"all distinct counts seen: {sorted(all_counts)}")
print(f"PRIME counts found: {sorted(set(primes_found))}")
print(f"181 reachable in this scan: {181 in all_counts}")
