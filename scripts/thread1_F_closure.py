#!/usr/bin/env python3
# Thread 1 — F-closure of Q90: the genuine closure-growth endofunctor
# F = M/~ (full multiway BFS + quotient) applied to Q90 as a fresh seed,
# vs the idempotent one-layer C-closure (the "Q180" retry).
# Exact Gaussian-integer arithmetic; validity-gated.
import numpy as np, itertools

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

ICS={0:[(2,1),(1,0),(3,-1)],1:[(1,0),(2,1),(1,-2)],2:[(1,-1),(3,0),(2,1)],
     3:[(3,0),(1,-1),(1,2)],4:[(1,2),(2,-1),(1,0)],5:[(2,0),(1,1),(3,0)]}

def adjacent_edges():
    s=set()
    for i in range(6):
        t=(i,(i+1)%6,(i+2)%6)
        for p in itertools.permutations(t): s.add(p)
    return sorted(s)
def complete_edges():
    return [(i,j,k) for i in range(6) for j in range(6) if j!=i
            for k in range(6) if k!=i and k!=j]
def g0_edges(): return [(i,(i+1)%6,(i+2)%6) for i in range(6)]
def q45_seed_edges(seed=42):
    base=set(adjacent_edges())
    pool=[e for e in complete_edges() if e not in base]
    rng=np.random.default_rng(seed); rng.shuffle(pool)
    return list(base)+pool[:6]

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

def c_closed(topo, depth):
    psi_o,edg_o=build_multiway(topo, ICS, depth)
    psi_c,edg_c=build_multiway(topo, {v:[gconj(x) for x in vec] for v,vec in ICS.items()}, depth)
    off=max(psi_o)+1
    allp=dict(psi_o)
    for v,p in psi_c.items(): allp[v+off]=p
    reps,vtc=cluster(allp)
    he=set()
    for (_,v1,v2,v3) in edg_o+[(d,a+off,b+off,c+off) for (d,a,b,c) in edg_c]:
        if v1 in vtc and v2 in vtc and v3 in vtc:
            he.add((vtc[v1],vtc[v2],vtc[v3]))
    return len(reps),reps,he

def F_closure(reps, he, depth):
    """genuine closure-growth F = M/~ : full multiway BFS over the object's
       own hyperedge topology, then quotient."""
    seed={i:reps[i] for i in range(len(reps))}
    psi,_=build_multiway(list(he), seed, depth)
    r,_=cluster(psi)
    return len(r)

def c_close_iterate(reps, he):
    """one-layer conjugate-union C-closure (the 'Q180' route)."""
    qp={i:reps[i] for i in range(len(reps))}
    def expand(seed):
        psi=dict(seed); nv=max(psi)+1
        for (s1,s2,_t) in he:
            if s1 not in psi or s2 not in psi: continue
            w=compose(psi[s1],psi[s2])
            if is_zero_vec(w): continue
            psi[nv]=w; nv+=1
        return psi
    po=expand(dict(qp)); pc=expand({v:[gconj(x) for x in qp[v]] for v in qp})
    off=max(po)+1; allp=dict(po)
    for v,p in pc.items(): allp[v+off]=p
    r,_=cluster(allp)
    return len(r)

# --- validity gate ---
g,_,_=c_closed(complete_edges(),5)
print(f"validity gate: K6^3 C-closed = {g}  (must be 102) -> {'PASS' if g==102 else 'FAIL'}")
assert g==102
print()

for name,topo,dep in (("Q48 (C-closed G0)",g0_edges(),5),
                       ("Q90 (C-closed q45)",q45_seed_edges(42),5),
                       ("Q102 (C-closed K6^3)",complete_edges(),5)):
    ncl,reps,he=c_closed(topo,dep)
    print(f"{name}: n_cl={ncl}, hyperedges={len(he)}")
    f4=F_closure(reps,he,4); f5=F_closure(reps,he,5); f6=F_closure(reps,he,6)
    cc=c_close_iterate(reps,he)
    print(f"   F-closure (genuine M/~)  depth4/5/6 : {f4} / {f5} / {f6}")
    print(f"   iterated C-closure (the 'Q180' route): {cc}")
    print(f"   -> F-fixed point? {'YES' if f5==ncl else 'NO ('+str(f5)+')'}    "
          f"C-closure idempotent? {'YES' if cc==ncl else 'NO ('+str(cc)+')'}")
    print()
