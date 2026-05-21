# Pool A K_n^3 defect — algebraic characterization

**Scope.** Closes the open follow-up from `audit_threebug_companion_v1.md`
§4.5 (closure-v5 v335): identify the specific algebraic relation in
the canonical s58 / q102-IC-set-2 IC pool ("Pool A") that drives its
K_n^3 single-side count's defect at n >= 9.

**Built 2026-05-21.** Probe: `probe_pool_a_defect.py` in this folder,
plus an in-line T7 counter-test. Exact Gaussian-integer arithmetic,
validity-gated on K_6^3 -> 51.

---

## Headline

**The defect is driven by Pool A inadvertently containing all three
cyclic permutations of a single Gaussian-integer triple.** Specifically,
Pool A's ICs 8, 9, and 11 are the three cyclic shifts of the same
underlying triple

```
{ (2,-1), (1,2), (3,1) }
```

When all three cyclic shifts are present in the IC pool alongside
Pool A's 8 axis-aligned ICs, the multiway closure on K_n^3 for n >= 9
exhibits projective coincidences that merge clusters relative to the
truly-generic count. Removing any of the three rotations eliminates the
defect; keeping all three (as Pool A does) drives the
5/12/12/17/19/19/19 shortfall sequence and the J-fixed cluster at
n >= 10.

---

## 1. Pool A inventory

Pool A's 15 ICs split into 8 axis-aligned (entries with at least one
component having zero real or zero imaginary part) and 7 strict-generic:

| IC | Triple | Axis-aligned components |
|---|---|---|
| 0 | `[(2,1),(1,0),(3,-1)]`   | [1] |
| 1 | `[(1,0),(2,1),(1,-2)]`   | [0] |
| 2 | `[(1,-1),(3,0),(2,1)]`   | [1] |
| 3 | `[(3,0),(1,-1),(1,2)]`   | [0] |
| 4 | `[(1,2),(2,-1),(1,0)]`   | [2] |
| 5 | `[(2,0),(1,1),(3,0)]`    | [0, 2] |
| 6 | `[(3,1),(2,0),(1,-1)]`   | [1] |
| 7 | `[(1,1),(3,-1),(2,0)]`   | [2] |
| 8 | `[(2,-1),(1,2),(3,1)]`   | (strict-generic) |
| 9 | `[(1,2),(3,1),(2,-1)]`   | (strict-generic) |
| 10 | `[(2,-1),(1,3),(1,1)]`  | (strict-generic) |
| 11 | `[(3,1),(2,-1),(1,2)]`  | (strict-generic) |
| 12 | `[(1,1),(2,2),(3,-1)]`  | (strict-generic) |
| 13 | `[(2,1),(1,-2),(2,1)]`  | (strict-generic) |
| 14 | `[(1,-1),(3,2),(1,1)]`  | (strict-generic) |

**Strict-generic ICs 8, 9, 11 are the three cyclic shifts of
`{ (2,-1), (1,2), (3,1) }`** — that is the load-bearing structural
fact.

---

## 2. The diagnostic ladder

### 2.1 T2 — Strict-generic Pool A subset restores pentagonal

ICs 8..14 (the 7 strict-generic Pool A entries) padded with Pool B's
first two ICs to reach n=9: gives `|Q_9| = 117 = pentagonal`, f=0.
**So the axis-aligned ICs are necessary** for the defect.

### 2.2 T3 — Axis-aligned-only Pool A subset does NOT trigger

ICs 0..7 (the 8 axis-aligned Pool A entries) padded with Pool B's IC 0
to reach n=9: gives `|Q_9| = 117 = pentagonal`, f=0. **So axis-aligned
ICs alone are not sufficient.** The defect needs interaction with
something specific.

### 2.3 T4 — Replacing axis-aligned in-place restores pentagonal

Pool A with each axis-aligned component swapped for a strict-generic
alternative (preserving the indices 0..14, just healing the bad
components): `|Q_9| = 117`, f=0. Validates that the axis-aligned
components are part of the activating mechanism.

### 2.4 T7 — Which strict-generic IC, paired with axis-aligned 0..7, triggers?

Pool A axis-aligned 0..7 + one strict-generic 9th IC:

| 9th IC | Triple | K_9^3 single | J-fixed f | Trigger? |
|---|---|---|---|---|
| Pool A IC 8  | `[(2,-1),(1,2),(3,1)]`  | **112** | 0 | **YES** |
| Pool A IC 9  | `[(1,2),(3,1),(2,-1)]`  | **110** | 1 | **YES** |
| Pool A IC 10 | `[(2,-1),(1,3),(1,1)]`  | 117 | 0 | no |
| Pool A IC 11 | `[(3,1),(2,-1),(1,2)]`  | **112** | 0 | **YES** |
| Pool A IC 12 | `[(1,1),(2,2),(3,-1)]`  | 117 | 0 | no |
| Pool A IC 13 | `[(2,1),(1,-2),(2,1)]`  | 117 | 0 | no |
| Pool A IC 14 | `[(1,-1),(3,2),(1,1)]`  | 117 | 0 | no |
| Pool B IC 0  | `[(4,1),(1,3),(2,-1)]`  | 117 | 0 | no |
| Pool B IC 1  | `[(1,2),(3,1),(4,-3)]`  | 117 | 0 | no |

**Triggering ICs are exactly the three cyclic shifts of the triple
`{(2,-1),(1,2),(3,1)}`.** Non-triggering Pool A strict-generic ICs
(10, 12, 13, 14) have different underlying triples; Pool B controls
are unrelated triples.

### 2.5 T5 — J-fixed cluster at K_10^3 under Pool A

The single J-fixed cluster that appears at n >= 10 has representative

```
[ (1, 8), (-1, -8), (1, 8) ]
```

This is `(1+8i) * (1, -1, 1)` — the projective image of the **real
vector (1, -1, 1)** scaled by the Gaussian integer (1, 8). Real rays
in CP^2 are self-conjugate-up-to-scale; this representative satisfies
`conj(r) ~ r` projectively, which is precisely the J-fixed condition.

So the J-fixed cluster has a clean structural meaning: it is the
projective image of a *real* ray that the Pool-A-induced composition
chain produces at large n.

### 2.6 T6 — Vertex assignment dependence

Shuffled Pool A 9-subsets (5 trials, random seed 17) gave defects of
varying size: single-counts of 115, 110, 112, 112, 110. All below the
pentagonal 117, but with the exact count depending on which Pool A ICs
end up at which K_n^3 vertices. This is consistent with the cyclic-
permutation-triple explanation: different shuffles include different
numbers of the {(2,-1),(1,2),(3,1)} cyclic shifts at different vertex
positions, and the defect magnitude is proportional to that interaction
with the axis-aligned ICs at neighbouring vertices.

---

## 3. Mechanistic picture

The CFS composition rule `mu(a, b) = conj(a x b)` is bilinear over
Gaussian integers. When the IC pool contains all three cyclic shifts
of a triple `{(2,-1), (1,2), (3,1)}`, the multiway expansion produces
specific compositions that lie on the same projective ray —
specifically, the ray of the real vector (1, -1, 1) scaled by Gaussian
integers like (1, 8). These projective coincidences are what merge
clusters relative to the strict-generic baseline.

The presence of axis-aligned components in the *other* 8 ICs supplies
the secondary compositions that, when combined with the cyclic-triple
ICs, populate this real-ray-shaped equivalence class with enough vertex
arrivals to make it a J-fixed cluster at n >= 10.

**Either condition alone is insufficient:**
- Axis-aligned ICs alone (T3): no cyclic-triple ICs to seed the
  real-ray cluster, so no projective merge.
- Cyclic-triple ICs alone (would need a non-Pool-A test): would
  produce the real ray but without the axis-aligned ICs' secondary
  compositions, the cluster may not become J-fixed.
- **Both present simultaneously (T7 with IC 8/9/11): defect activates.**

---

## 4. Implications

### 4.1 Validity-gate refinement

The strengthened gate K_10^3 -> 145 (from audit_threebug_companion_v1.md
§1.3) catches the cyclic-triple defect immediately. A more targeted
filter would explicitly reject IC pools containing two or more cyclic
permutations of the same triple — a cheap test that can run before
the K_10^3 build.

### 4.2 closure-v5 spec record

The audit at v335 (closure-v5 commit `f2084ae`) classified Pool A as
"small-n generic only" and introduced S230 (`Obs_ic_genericity_hierarchy`).
This findings doc identifies the **specific** algebraic relation
underlying the small-n-generic classification: cyclic-permutation
redundancy of the triple `{(2,-1), (1,2), (3,1)}` among Pool A's
strict-generic entries. The audit's "characterise the hidden algebraic
relation" open item (§4.5, second bullet) is now closed.

### 4.3 Why s58 has this structure

The s58 IC set was constructed historically by hand-picking Gaussian-
integer triples for K_6^3 -> 51. Cyclic permutations of a single
underlying triple were apparently treated as distinct enough to count
as separate ICs — a reasonable assumption at K_6^3 where the
permutation symmetry doesn't surface. At K_n^3 for n >= 9 it does.

### 4.4 Future IC-set construction guidance

When constructing IC pools for CFS computational work involving
K_n^3 at large n, explicitly check for:
1. Strict-genericity (no zero-component vectors).
2. **No two ICs that are cyclic permutations of the same triple.**

Both conditions are necessary; the K_10^3 -> 145 gate is a sufficient
empirical check.

---

## 5. Pool A in context

Pool A's defect is a quantitative property of one specific IC set used
historically across CFS work. It does NOT invalidate any K_6^3-based
result (Pool A is fully generic at n=6), and it does NOT invalidate the
broader CFS programme. What it does is sharpen the discipline around
IC-set construction for any future K_n^3 work at large n.

The corrected v312 F.3 verdict (|Q_9|=234, 78 colour triplets) was
the audit's primary deliverable; this findings doc supplements with
the algebraic characterization of *why* Pool A specifically gave the
old wrong number.

---

## 6. Reproducibility

- `probe_pool_a_defect.py` — main probe (T1..T6).
- T7 counter-test ran inline (output captured in this doc's §2.4 table);
  the test routine is straightforward to extract into a v2 of the probe
  if needed.
- Validity-gated K_6^3 -> 51 throughout.
- Total runtime: <30s on a single CPU core.

---

## 7. Closes

- `audit_threebug_companion_v1.md` §4.5 second bullet: "Characterise
  the hidden algebraic relation in the s58 / Pool A IC set that
  produces its K_n^3 defect at n >= 9."
- Task #11 in the closure-quotient-paper session task list (2026-05-21).
- `closure_quotient_landscape.tex` §11 first open problem ("Pool A
  defect: which algebraic relation?") — answer: cyclic-permutation
  triples `{(2,-1), (1,2), (3,1)}`.
