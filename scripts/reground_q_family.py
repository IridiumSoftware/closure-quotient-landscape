#!/usr/bin/env python3
# Re-grounding: recompute the Q45/Q84/Q90/Q98/Q180 family under EXACT
# Gaussian-integer arithmetic + exact ray-equivalence — replacing the
# fidelity-0.999 Float64 clustering the closure-v5 family was built with.
# Topologies are reproduced faithfully (incl. numpy RNG seed 42).

import numpy as np

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

# exact ICs, keyed 0..5 (landscape code is 0-indexed)
ICS = {0:[(2,1),(1,0),(3,-1)], 1:[(1,0),(2,1),(1,-2)], 2:[(1,-1),(3,0),(2,1)],
       3:[(3,0),(1,-1),(1,2)], 4:[(1,2),(2,-1),(1,0)], 5:[(2,0),(1,1),(3,0)]}

# --- topologies, faithful to quotient_landscape_autopoiesis_v1.py (0-indexed) ---
def adjacent_edges():
    import itertools
    s=set()
    for i in range(6):
        t=(i,(i+1)%6,(i+2)%6)
        for p in itertools.permutations(t): s.add(p)
    return sorted(s)
def complete_edges():
    return [(i,j,k) for i in range(6) for j in range(6) if j!=i
            for k in range(6) if k!=i and k!=j]
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

def cluster(psi):  # 0-indexed
    reps=[]; vtc={}
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        m=-1
        for ci,r in enumerate(reps):
            if proj_equiv(pv,r): m=ci; break
        if m>=0: vtc[v]=m
        else: reps.append(pv); vtc[v]=len(reps)-1
    return reps, vtc

def single_side(topo, depth):
    psi,_=build_multiway(topo, ICS, depth)
    reps,_=cluster(psi)
    return len(reps)

def c_closed(topo, depth):
    psi_o,edg_o=build_multiway(topo, ICS, depth)
    psi_c,edg_c=build_multiway(topo, {v:[gconj(x) for x in vec] for v,vec in ICS.items()}, depth)
    off=max(psi_o)+1
    allp=dict(psi_o)
    for v,p in psi_c.items(): allp[v+off]=p
    reps,vtc=cluster(allp)
    edges_off=[(d,v1+off,v2+off,v3+off) for (d,v1,v2,v3) in edg_c]
    he=set()
    for (_,v1,v2,v3) in edg_o+edges_off:
        if v1 in vtc and v2 in vtc and v3 in vtc:
            he.add((vtc[v1],vtc[v2],vtc[v3]))
    return len(reps), reps, he

def c_close_iterate(reps, he):
    # exact port of c_close_quotient (q90_q180_explore_v3.py)
    q_psi={i:reps[i] for i in range(len(reps))}
    def expand(seed):
        psi=dict(seed); nv=max(psi)+1
        for (s1,s2,_t) in he:
            if s1 not in psi or s2 not in psi: continue
            w=compose(psi[s1],psi[s2])
            if is_zero_vec(w): continue
            psi[nv]=w; nv+=1
        return psi
    psi_o=expand(dict(q_psi))
    psi_c=expand({v:[gconj(x) for x in q_psi[v]] for v in q_psi})
    off=max(psi_o)+1
    allp=dict(psi_o)
    for v,p in psi_c.items(): allp[v+off]=p
    reps2,_=cluster(allp)
    return len(reps2)

print("EXACT re-grounding of the closure-quotient family")
print("(published counts used fidelity-0.999 Float64 clustering)")
print()
print(f"{'construction':<34}{'published':>11}{'exact':>9}{'drift':>8}")
print("-"*62)
adj=adjacent_edges(); k6=complete_edges(); q45=q45_seed_edges(42)

rows=[]
ss_adj=single_side(adj,5);                 rows.append(("adjacent-36  single-side", None, ss_adj))
cc_adj,_,_=c_closed(adj,5);                rows.append(("adjacent-36  C-closed  (Q84)", 84, cc_adj))
ss_q45=single_side(q45,5);                 rows.append(("q45-seed     single-side (Q45)", 45, ss_q45))
cc_q45,reps90,he90=c_closed(q45,5);        rows.append(("q45-seed     C-closed  (Q90)", 90, cc_q45))
ss_k6=single_side(k6,5);                   rows.append(("K6^3         single-side (Q51)", 51, ss_k6))
cc_k6,_,_=c_closed(k6,5);                  rows.append(("K6^3         C-closed  (Q98)", 98, cc_k6))
q180=c_close_iterate(reps90,he90);         rows.append(("iterated C-closure of Q90 (Q180)", 180, q180))

for name,pub,ex in rows:
    if pub is None:
        print(f"{name:<34}{'—':>11}{ex:>9}{'':>8}")
    else:
        d=ex-pub
        print(f"{name:<34}{pub:>11}{ex:>9}{('+' if d>=0 else '')+str(d):>8}")
print()
print(f"q45 seed: {len(q45)} edges (36 adjacent + 6 RNG-42 extras)")
