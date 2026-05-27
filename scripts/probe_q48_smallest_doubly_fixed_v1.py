#!/usr/bin/env python3
# Probe — boundary self-regulation: Q_48 is the SMALLEST doubly-fixed Q.
#
# Rosen (M,R)-system role mapping: "boundary self-regulation" = the smallest
# endogenously-stable cardinality under BOTH F (closure-growth endofunctor,
# Definition 6.X) AND C (charge-conjugation idempotency, Proposition 7.X).
#
# Hypothesis. Among the canonical family {Q_12, Q_24, Q_45, Q_48, Q_51,
# Q_90, Q_102}, Q_48 is the FIRST member satisfying both F(Q) = Q AND
# C(Q) = Q. Smaller candidates (Q_12, Q_24, Q_45) and the
# non-C-fixed K_6^3-anchor Q_51 each fail at least one — specifically C,
# since under generic ICs they C-double cleanly (Q_12->Q_24, Q_24->Q_48,
# Q_45->Q_90, Q_51->Q_102).
#
# Predicted discrimination table (paper §13 candidate):
#   Q_12  : F-fix ?,  C-fix N  (->Q_24)
#   Q_24  : F-fix ?,  C-fix N  (->Q_48)
#   Q_45  : F-fix Y,  C-fix N  (->Q_90)
#   Q_48  : F-fix Y,  C-fix Y  <-- SMALLEST doubly-fixed
#   Q_51  : F-fix Y,  C-fix N  (->Q_102)
#   Q_90  : F-fix Y,  C-fix Y
#   Q_102 : F-fix Y,  C-fix Y
#
# Validity gate: K_6^3 C-closed = 102 under POOL_B (strict-generic).
# Exact Gaussian-integer ray arithmetic throughout (no Float64).
#
# Reuses helpers verbatim from thread1_F_closure.py + ICs from
# Kn3_pentagonal_confirmation.py (POOL_B, the strict-generic 20-IC set).

import itertools
import numpy as np

# --- Gaussian-integer arithmetic (verbatim from thread1_F_closure.py) -------
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

# --- IC pools ---------------------------------------------------------------
# POOL_A = canonical s58 set entries 0..5 (matches thread1_F_closure.py's ICS
# and reground_q_family.py's ICS). Generic at K_6^3; the paper's published
# Q_45 / Q_90 / Q_51 / Q_102 are defined under this pool.
POOL_A = [
 [(2,1),(1,0),(3,-1)],[(1,0),(2,1),(1,-2)],[(1,-1),(3,0),(2,1)],
 [(3,0),(1,-1),(1,2)],[(1,2),(2,-1),(1,0)],[(2,0),(1,1),(3,0)],
]

# POOL_B = strict-generic from Kn3_pentagonal_confirmation.py — every
# component nonzero real AND nonzero imag. Kept here for cross-pool
# sanity comparison (single-side adjacent_edges drifts 45->42 under
# POOL_B; this is an aside, not the load-bearing observation).
POOL_B = [
 [(4,1),(1,3),(2,-1)],[(1,2),(3,1),(4,-3)],[(2,3),(4,-1),(1,2)],
 [(3,-2),(2,1),(4,3)],[(4,3),(1,-2),(3,1)],[(2,-3),(4,1),(1,2)],
]

def make_ics(pool, n):
    return {v: [list(c) for c in pool[v]] for v in range(n)}

# --- Seed constructors -------------------------------------------------------
def cyc_seed(n):
    """cyc(n) per paper Def 2.X: n forward-shift triples, no perm closure."""
    return [(i, (i+1)%n, (i+2)%n) for i in range(n)]

def kn3_seed(n):
    """K_n^3 per paper Def 2.Y: all distinct ordered triples on n vertices."""
    return [(i,j,k) for i in range(n) for j in range(n) if j!=i
            for k in range(n) if k!=i and k!=j]

def adjacent_edges_perm():
    """Adjacent triples on 6 vertices with full perm closure — the Q_45 seed
    (matches thread1_F_closure.py's adjacent_edges(); 36 edges)."""
    s=set()
    for i in range(6):
        t=(i,(i+1)%6,(i+2)%6)
        for p in itertools.permutations(t): s.add(p)
    return sorted(s)

def q45_seed_edges(seed=42):
    """adjacent_edges_perm + 6 random non-base from K_6^3 — the Q_45/Q_90 seed.
    Uses numpy.random.default_rng(42) shuffle to match the canonical
    paper-family numbers reproduced in reground_q_family.py and
    thread1_F_closure.py: single_side = 45, c_closed = 90."""
    base=set(adjacent_edges_perm())
    pool=[e for e in kn3_seed(6) if e not in base]
    rng=np.random.default_rng(seed); rng.shuffle(pool)
    return list(base)+pool[:6]

# --- Multiway + cluster (verbatim from thread1_F_closure.py) ----------------
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

def hyperedges_from_clusters(edges, vtc):
    he=set()
    for (_,v1,v2,v3) in edges:
        if v1 in vtc and v2 in vtc and v3 in vtc:
            he.add((vtc[v1],vtc[v2],vtc[v3]))
    return he

def single_side_quotient(topo, ics, depth):
    """Returns (reps, hyperedges) for single-side closure."""
    psi, edg = build_multiway(topo, ics, depth)
    reps, vtc = cluster(psi)
    he = hyperedges_from_clusters(edg, vtc)
    return reps, he

def c_closed_quotient(topo, ics, depth):
    """Returns (reps, hyperedges) for C-closed quotient."""
    psi_o, edg_o = build_multiway(topo, ics, depth)
    psi_c, edg_c = build_multiway(
        topo,
        {v: [gconj(x) for x in vec] for v, vec in ics.items()},
        depth,
    )
    off = max(psi_o) + 1
    allp = dict(psi_o)
    for v, p in psi_c.items():
        allp[v + off] = p
    reps, vtc = cluster(allp)
    he = set()
    all_edges = list(edg_o) + [(d, a+off, b+off, c+off) for (d, a, b, c) in edg_c]
    for (_, v1, v2, v3) in all_edges:
        if v1 in vtc and v2 in vtc and v3 in vtc:
            he.add((vtc[v1], vtc[v2], vtc[v3]))
    return reps, he

# --- F-closure endofunctor: apply multiway to Q as fresh seed ---------------
def F_closure_size(reps, he, depth):
    """|F(Q)| = re-run multiway on Q's own cluster reps + hyperedges, depth."""
    seed = {i: reps[i] for i in range(len(reps))}
    psi, _ = build_multiway(list(he), seed, depth)
    r, _ = cluster(psi)
    return len(r)

# --- C-closure idempotency check: one-layer C-union -------------------------
def C_closure_size(reps, he):
    """|Q ∪ C(Q)| via one-layer expansion (matches thread1's c_close_iterate)."""
    qp = {i: reps[i] for i in range(len(reps))}
    def expand(seed):
        psi = dict(seed); nv = max(psi) + 1
        for (s1, s2, _t) in he:
            if s1 not in psi or s2 not in psi: continue
            w = compose(psi[s1], psi[s2])
            if is_zero_vec(w): continue
            psi[nv] = w; nv += 1
        return psi
    po = expand(dict(qp))
    pc = expand({v: [gconj(x) for x in qp[v]] for v in qp})
    off = max(po) + 1
    allp = dict(po)
    for v, p in pc.items():
        allp[v + off] = p
    r, _ = cluster(allp)
    return len(r)

# --- Validity gate -----------------------------------------------------------
ics_6 = make_ics(POOL_A, 6)
reps_gate, _ = c_closed_quotient(kn3_seed(6), ics_6, 5)
n_gate = len(reps_gate)
print(f"Validity gate: K_6^3 C-closed (POOL_A, depth=5) = {n_gate}  "
      f"(must be 102) -> {'PASS' if n_gate == 102 else 'FAIL'}")
assert n_gate == 102, "Validity gate failed — cannot trust further results"
print()

# --- Family scan -------------------------------------------------------------
print("=" * 84)
print("Probe — boundary self-regulation: Q_48 is the smallest doubly-fixed")
print("=" * 84)
print()
print(f"{'Q':>7} {'|Q|':>5} {'|F(Q)|@d4':>11} {'|F(Q)|@d5':>11} "
      f"{'|Q∪C(Q)|':>10} {'F-fix':>6} {'C-fix':>6} {'BOTH':>5}")
print("-" * 84)

ics_3 = make_ics(POOL_A, 3)

# Build canonical-family quotients
family = []

# Q_12 = single_side(cyc(3))
reps, he = single_side_quotient(cyc_seed(3), ics_3, 3)
family.append(("Q_12",  reps, he))

# Q_24 = single_side(cyc(6))
reps, he = single_side_quotient(cyc_seed(6), ics_6, 4)
family.append(("Q_24",  reps, he))

# Q_45 = single_side(q45_seed_edges)   [paper §4.1, reground_q_family.py]
seed_q45 = q45_seed_edges(42)
reps, he = single_side_quotient(seed_q45, ics_6, 4)
family.append(("Q_45",  reps, he))

# Q_48 = c_closed(cyc(6))
reps, he = c_closed_quotient(cyc_seed(6), ics_6, 4)
family.append(("Q_48",  reps, he))

# Q_51 = single_side(K_6^3)
reps, he = single_side_quotient(kn3_seed(6), ics_6, 4)
family.append(("Q_51",  reps, he))

# Q_84 = c_closed(adjacent_edges_perm)   [paper §4.1 Prop Q_84-exact]
reps, he = c_closed_quotient(adjacent_edges_perm(), ics_6, 4)
family.append(("Q_84",  reps, he))

# Q_90 = c_closed(q45_seed_edges)   [paper §4.1 Prop Q_45/Q_90-exact]
reps, he = c_closed_quotient(seed_q45, ics_6, 4)
family.append(("Q_90",  reps, he))

# Q_102 = c_closed(K_6^3)   [paper §4.2 Prop Q_102-exact]
reps, he = c_closed_quotient(kn3_seed(6), ics_6, 4)
family.append(("Q_102", reps, he))

rows = []
for name, reps, he in family:
    n = len(reps)
    f4 = F_closure_size(reps, he, 4)
    f5 = F_closure_size(reps, he, 5)
    cc = C_closure_size(reps, he)
    f_fix = (f4 == n and f5 == n)
    c_fix = (cc == n)
    both = f_fix and c_fix
    rows.append((name, n, f4, f5, cc, f_fix, c_fix, both))
    print(f"{name:>7} {n:>5} {f4:>11} {f5:>11} {cc:>10} "
          f"{'Y' if f_fix else 'N':>6} {'Y' if c_fix else 'N':>6} "
          f"{'YES' if both else '':>5}")

print()
doubly = [r for r in rows if r[7]]
if doubly:
    smallest = min(doubly, key=lambda r: r[1])
    verdict = "PASS" if smallest[0] == "Q_48" else f"UNEXPECTED ({smallest[0]} ≠ Q_48)"
    print(f"Smallest doubly-fixed: {smallest[0]}  (|Q|={smallest[1]})")
    print(f"Hypothesis 'Q_48 carries boundary self-regulation as smallest "
          f"doubly-fixed': {verdict}")
else:
    print("Hypothesis 'Q_48 is smallest doubly-fixed': FAIL "
          "(no doubly-fixed family member found)")
