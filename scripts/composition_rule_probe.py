#!/usr/bin/env python3
# Composition rule — combine two seeds A,B and ask whether the closure
# count Q(A o B) is a derivable function of Q(A), Q(B).
# Modes: disjoint sum (no shared verts), vertex-glued (g shared), bridged.
# Exact Gaussian-integer arithmetic; validity-gated.

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

# many distinct generic ICs (validated: K6^3 -> 51)
_GEN=[
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
 [(3,1),(2,0),(1,-1)],[(1,1),(3,-1),(2,0)],[(2,-1),(1,2),(3,1)],
 [(1,2),(3,1),(2,-1)],[(2,-1),(1,3),(1,1)],[(3,1),(2,-1),(1,2)],
 [(1,1),(2,2),(3,-1)],[(2,1),(1,-2),(2,1)],[(1,-1),(3,2),(1,1)],
]
def ics(n): return {v:[list(c) for c in _GEN[v]] for v in range(n)}

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

def closure_count(topo, nverts, depth=5):
    psi=build_multiway(topo, ics(nverts), depth)
    reps=[]
    for v in sorted(psi):
        pv=psi[v]
        if is_zero_vec(pv): continue
        if not any(proj_equiv(pv,r) for r in reps): reps.append(pv)
    return len(reps)

def cyc(n, base=0):  # cyclic triples on vertices base..base+n-1
    return [(base+i, base+(i+1)%n, base+(i+2)%n) for i in range(n)]

# validity gate
g=closure_count([(i,(i+1)%6,(i+2)%6) for i in range(6)],6,5)
print(f"gate: cyc(6) -> {g} (must be 24) -> {'PASS' if g==24 else 'FAIL'}")
assert g==24
print()

def glued(nA,nB,g):
    # B shares g vertices: B verts 0..g-1 = A verts nA-g..nA-1
    def mv(x): return (nA-g+x) if x<g else (nA-g+x)
    Bg=[tuple((nA-g+x) for x in ((i)%nB,(i+1)%nB,(i+2)%nB)) for i in range(nB)]
    return closure_count(cyc(nA,0)+Bg, nA+nB-g, 5)

print("COMPOSITION PROBE — cyc(nA) o cyc(nB)")
print("delta = Q(AoB) - Q(A) - Q(B)  (interface term)")
print(f"{'A':>4}{'B':>4}{'QA':>5}{'QB':>5} |"
      f"{'disj':>6}{'glu1':>6}{'glu2':>6}{'glu3':>6}{'br1':>6}{'br2':>6}  | deltas")
print("-"*82)
for nA,nB in [(3,3),(3,4),(4,4),(3,5),(4,5),(5,5),(3,6),(4,6)]:
    qA=closure_count(cyc(nA),nA,5); qB=closure_count(cyc(nB),nB,5)
    s=qA+qB
    qD=closure_count(cyc(nA,0)+cyc(nB,nA),nA+nB,5)
    qG1=glued(nA,nB,1); qG2=glued(nA,nB,2); qG3=glued(nA,nB,3)
    qB1=closure_count(cyc(nA,0)+cyc(nB,nA)+[(nA-1,nA,nA+1)],nA+nB,5)
    qB2=closure_count(cyc(nA,0)+cyc(nB,nA)+[(nA-1,nA,nA+1),(nA-2,nA+2,nA+3)],nA+nB,5)
    ds=[v-s for v in (qD,qG1,qG2,qG3,qB1,qB2)]
    print(f"{nA:>4}{nB:>4}{qA:>5}{qB:>5} |"
          f"{qD:>6}{qG1:>6}{qG2:>6}{qG3:>6}{qB1:>6}{qB2:>6}  | "
          f"disj{ds[0]:+d} g1{ds[1]:+d} g2{ds[2]:+d} g3{ds[3]:+d} br1{ds[4]:+d} br2{ds[5]:+d}")
print()
print("rule: Q(AoB) = Q(A) + Q(B) + delta(interface), delta universal (A,B-independent)")
