#!/usr/bin/env python3
# Exploratory scan — exact closure cluster counts for small seeds.
# Maps the small-object landscape of the CFS closure model. Exact
# Gaussian-integer arithmetic (no fidelity threshold). EXPLORATION, not a
# pre-registered test — measuring, not claiming.

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

ICS = {
 1:[(2,1),(1,0),(3,-1)], 2:[(1,0),(2,1),(1,-2)], 3:[(1,-1),(3,0),(2,1)],
 4:[(3,0),(1,-1),(1,2)], 5:[(1,2),(2,-1),(1,0)], 6:[(2,0),(1,1),(3,0)],
}

def build_multiway(topo, psi_init, depth):
    psi={k:list(v) for k,v in psi_init.items()}
    nv=max(psi_init)+1
    edges=[(0,s1,s2,s3) for (s1,s2,s3) in topo]
    cache={}
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
        m=0
        for ci,r in enumerate(reps):
            if proj_equiv(pv,r): m=ci+1; break
        if m>0: vtc[v]=m
        else: reps.append(pv); vtc[v]=len(reps)
    return reps, vtc

def single_side(topo, psi_init, depth):
    psi,_=build_multiway(topo,psi_init,depth)
    reps,_=cluster(psi)
    return len(reps)

def c_closed(topo, psi_init, depth):
    psi_o,_=build_multiway(topo,psi_init,depth)
    psi_c,_=build_multiway(topo,{v:[gconj(x) for x in vec] for v,vec in psi_init.items()},depth)
    off=max(psi_o)+1
    allp={}
    for v,p in psi_o.items(): allp[v]=p
    for v,p in psi_c.items(): allp[v+off]=p
    reps,_=cluster(allp)
    return len(reps)

def Kn(n):  # complete ternary on n verts
    return [(i,j,k) for i in range(1,n+1) for j in range(1,n+1) if j!=i
            for k in range(1,n+1) if k!=i and k!=j]
def cyc(n):  # adjacent cyclic triples (G0-style)
    return [((i)%n+1,(i+1)%n+1,(i+2)%n+1) for i in range(n)]

print("EXACT small-closure landscape — single-side / C-closed cluster counts")
print("(known anchors: G0=cyc(6) d5 single->24 ; K6^3 d5 C-closed->102)")
print()
print(f"{'seed':<14}{'verts':>6}{'edges':>7}  | single-side d4/d5/d6  | C-closed d4/d5/d6")
print("-"*78)
for n in (3,4,5,6):
    ic={i:ICS[i] for i in range(1,n+1)}
    for name,topo in (("cyc(%d)"%n,cyc(n)),("K%d^3"%n,Kn(n))):
        ss=[single_side(topo,ic,d) for d in (4,5,6)]
        cc=[c_closed(topo,ic,d) for d in (4,5,6)]
        print(f"{name:<14}{n:>6}{len(topo):>7}  | "
              f"{ss[0]:>5}{ss[1]:>5}{ss[2]:>5}        | {cc[0]:>5}{cc[1]:>5}{cc[2]:>5}")
print()
print("looking for: small generators (3,5,9) and the factor family (45,90,180)")
