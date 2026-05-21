#!/usr/bin/env python3
# Thread 2 — does C-closure carry a self-conjugate (J-fixed) core?
# Direct check: compute the J involution (charge conjugation) on each
# C-closed object and count its fixed points. A J-fixed cluster c is one
# with conj(rep_c) ~ rep_c. #J-fixed = f  =>  C-closed count = 2n - f.
# Also: does a degenerate (real) IC set create self-conjugate clusters?
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
# degenerate: purely REAL ICs (im=0) -> conj(psi)=psi
ICS_REAL={v:[(c[0],0) for c in vec] for v,vec in ICS.items()}

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
    return psi

def cluster(psi):
    reps=[]
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv,r) for r in reps): reps.append(pv)
    return reps

def single_side(topo, ics, depth):
    return cluster(build_multiway(topo, ics, depth))

def c_closed_reps(topo, ics, depth):
    psi_o=build_multiway(topo, ics, depth)
    psi_c=build_multiway(topo, {v:[gconj(x) for x in vec] for v,vec in ics.items()}, depth)
    allp=dict(psi_o)
    off=max(psi_o)+1
    for v,p in psi_c.items(): allp[v+off]=p
    return cluster(allp)

def j_fixed_count(reps):
    """J = charge conjugation. count clusters with conj(rep_c) ~ rep_c."""
    n=len(reps); jmap=[]
    fixed=0; complete=True
    for c in range(n):
        cr=[gconj(x) for x in reps[c]]
        partner=-1
        for c2 in range(n):
            if proj_equiv(cr, reps[c2]): partner=c2; break
        jmap.append(partner)
        if partner==-1: complete=False
        elif partner==c: fixed+=1
    return fixed, complete

print("THREAD 2 — J-fixed-point (self-conjugate cluster) check")
g=len(c_closed_reps(complete_edges(),ICS,5))
print(f"validity gate: K6^3 C-closed = {g} (must be 102) -> {'PASS' if g==102 else 'FAIL'}")
assert g==102
print()
print(f"{'object':<24}{'single':>8}{'C-closed':>10}{'J-fixed':>9}{'2n-f check':>12}")
print("-"*64)
for name,topo in (("G0  -> Q24/Q48",g0_edges()),
                   ("q45 -> Q45/Q90",q45_seed_edges(42)),
                   ("K6^3-> Q51/Q102",complete_edges())):
    ss=len(single_side(topo,ICS,5))
    reps=c_closed_reps(topo,ICS,5)
    cc=len(reps)
    f,complete=j_fixed_count(reps)
    chk = "OK" if cc==2*ss-f else "MISMATCH"
    print(f"{name:<24}{ss:>8}{cc:>10}{f:>9}{f'2*{ss}-{f}={2*ss-f} {chk}':>12}")
print()
print("--- degenerate control: purely REAL ICs (conj(psi)=psi) ---")
for name,topo in (("G0 real-IC",g0_edges()),("K6^3 real-IC",complete_edges())):
    ss=len(single_side(topo,ICS_REAL,5))
    reps=c_closed_reps(topo,ICS_REAL,5)
    cc=len(reps); f,_=j_fixed_count(reps)
    print(f"  {name:<16} single={ss}  C-closed={cc}  J-fixed={f}  "
          f"(real ICs -> conj=identity, sectors collapse)")
