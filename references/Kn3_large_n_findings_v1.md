# K_n^3 large-n closure law — findings

**Scope.** Pin the closure law for `K_n^3` under generic ICs beyond n = 6,
addressing the open item in `closure_quotient_landscape_findings_v1.md`
§8: *"Depth 4 may under-close `K_n^3` for n ≥ 9 — the large-n single-side
counts (112/133/164/193 at n=9..12) sit below the pentagonal formula;
the exact `K_n^3` law beyond n = 6 is not pinned."*

**Built 2026-05-20.** Aaron driving; exact Gaussian-integer arithmetic
(no fidelity threshold), three independent IC pools, depth-stability
verified up to depth 8. Scripts in this folder:
`Kn3_large_n_probe.py`, `Kn3_extended_probe.py`,
`Kn3_pentagonal_confirmation.py`, `Kn3_push_to_n20.py`.

Use the closure-v5 venv: `Closure v5/.venv/bin/python`.

---

## Summary

1. **K_n^3 single-side = n(3n−1)/2 (pentagonal) for all measured n=6..20**
   under truly generic ICs. Confirmed by two independent strict-generic
   IC pools (Pool B and Pool C, defined below) — exact agreement at
   every n.
2. **Depth was not the issue.** Counts are stable from depth 3 onward at
   every measured n; extending to depth 8 changes nothing.
3. **The original measurement (112/133/164/193 at n=9..12) was a
   Pool-A artifact**, not a property of `K_n^3`. The canonical 15-vector
   s58/q102 IC set used in every prior CFS probe — Pool A in this
   report — is not fully generic at n ≥ 9. Under Pool A the single-side
   count drops below pentagonal and the charge-conjugation involution
   J acquires exactly one fixed cluster at n ≥ 10.
4. **The IC-genericity story sharpens** from a binary (generic vs
   degenerate) to a hierarchy: there are IC sets that pass the
   K_6^3 → 51 validity gate but fail to be fully generic at larger n.
5. **The closure_quotient_landscape_findings_v1.md §8 open item is
   closed** with a stronger result than its original framing implied.

---

## 1. The three IC pools

All three pools are strict-generic *at n = 6* — every pool gives
K_6^3 → 51 (the canonical validity gate). They diverge in their
behavior at larger n.

| pool | size | coefficient range | construction | provenance |
|---|---|---|---|---|
| **A** | 15 | coords in {−2,−1,0,1,2,3} | the canonical s58 9-vector set + q102 IC-set-2 extension | every prior CFS probe |
| **B** | 20 | coords in {−4,..,−1,1,..,4} | strict: every component has nonzero real AND nonzero imag | constructed for this report |
| **C** | 20 | coords in {−7,..,−3,3,..,7} | strict-generic, deliberately disjoint coefficient range from B | constructed for this report |

Pool A includes coordinates with zero real or zero imag part
(e.g. `(1,0)`, `(3,0)`, `(2,0)`). Pools B and C exclude these by
construction — every component of every triple is a Gaussian integer
with both parts nonzero.

---

## 2. The pentagonal law confirmed

`K_n^3` single-side closure counts at depth 5 (depth-stable), three pools:

| n | \|K_n^3\| | pentagonal n(3n−1)/2 | Pool A | Pool B | Pool C | B==C==pent? |
|---|---|---|---|---|---|---|
| 6 | 120 | 51 | 51 | 51 | 51 | YES |
| 7 | 210 | 70 | 70 | 70 | 70 | YES |
| 8 | 336 | 92 | 92 | 92 | 92 | YES |
| 9 | 504 | 117 | **112** | 117 | 117 | YES |
| 10 | 720 | 145 | **133** | 145 | 145 | YES |
| 11 | 990 | 176 | **164** | 176 | 176 | YES |
| 12 | 1320 | 210 | **193** | 210 | 210 | YES |
| 13 | 1716 | 247 | **228** | 247 | 247 | YES |
| 14 | 2184 | 287 | **268** | 287 | 287 | YES |
| 15 | 2730 | 330 | **311** | 330 | 330 | YES |
| 16 | 3360 | 376 | — | 376 | 376 | YES |
| 17 | 4080 | 425 | — | 425 | 425 | YES |
| 18 | 4896 | 477 | — | 477 | 477 | YES |
| 19 | 5814 | 532 | — | 532 | 532 | YES |
| 20 | 6840 | 590 | — | 590 | 590 | YES |

(Pool A only goes to n=15 — pool size limit.)

**Pool B and Pool C agree exactly at every n=6..20** despite being
constructed from disjoint coefficient ranges and different sign
patterns. This is the strongest available empirical evidence for
IC-independence under genuine genericity.

**The K_n^3 → n(3n−1)/2 law holds across all 15 values of n tested.**
A closed-form algebraic proof is not in this document; the experimental
case is at the limit of what depth-5 exhaustive cluster enumeration can
deliver in seconds of compute (n=20 ran in ~2 seconds per pool).

---

## 3. Depth-stability

The first probe (`Kn3_large_n_probe.py`) ran each `K_n^3` at depths
3, 4, 5, 6, 7, 8 for n = 6..12 and the count was identical at every
depth in every row — Pool A's "depression" was depth-stable, not a
depth artifact.

For Pool B and Pool C: depth=5 reaches the same count as depth=3 at
every n tested. The closure-quotient family is depth-stable at depth ≤ 5
for K_n^3 in this n range, in agreement with the
closure_quotient_landscape_findings §2 claim of depth-stability by
depth 4 for every small object.

---

## 4. The J fixed-point check

Charge-conjugation involution J: a cluster `c` is J-fixed iff
`conj(rep(c)) ~ rep(c)` projectively. The identity
`|Q ∪ C(Q)| = 2·|Q| − f` recovers the J-fixed count `f` from the
single-side count and the C-closed count.

| n | Pool A f | Pool B f | Pool C f |
|---|---|---|---|
| 6 | 0 | 0 | 0 |
| 7 | 0 | 0 | 0 |
| 8 | 0 | 0 | 0 |
| 9 | 0 | 0 | 0 |
| 10 | **1** | 0 | 0 |
| 11 | **1** | 0 | 0 |
| 12 | **1** | 0 | 0 |
| 13 | **1** | 0 | 0 |
| 14 | **1** | 0 | 0 |
| 15 | **1** | 0 | 0 |

**Under Pools B and C, J is fixed-point-free at every measured n.**
This is the sharper version of the
closure_quotient_landscape_findings §5 claim. The earlier statement
("J fixed-point-free under generic ICs") was correct in spirit but
under-specified: J is fixed-point-free under *truly* generic ICs;
under Pool A specifically, J acquires one fixed point at n ≥ 10.

Consequently the C-closure clean-doubling law `|Q ∪ C(Q)| = 2·|Q|`
holds for Pools B and C at every n tested, but **fails for Pool A at
n ≥ 10** (where the formula gives 2|Q| + 1 instead).

---

## 5. The Pool A degeneracy — a hidden algebraic relation

Pool A's shortfall (`pent − single_side`) and J-fixed count:

| n | shortfall | J-fixed |
|---|---|---|
| 6 | 0 | 0 |
| 7 | 0 | 0 |
| 8 | 0 | 0 |
| 9 | 5 | 0 |
| 10 | 12 | 1 |
| 11 | 12 | 1 |
| 12 | 17 | 1 |
| 13 | 19 | 1 |
| 14 | 19 | 1 |
| 15 | 19 | 1 |

The shortfall climbs irregularly through n = 9..12 (5, 12, 12, 17) and
**plateaus at 19 for n = 13..15**. The J-fixed count saturates at
exactly 1 from n = 10 onward.

**Interpretation.** Pool A contains a hidden small-coefficient
algebraic relation among its 15 triples (likely related to its
inclusion of axis-aligned components like `(1,0)`, `(3,0)`, `(2,0)`
and the appearance of repeated coordinate values across triples).
The relation is invisible at K_n^3 for n ≤ 8 (no triple subset large
enough to surface it), produces an irregular cluster collapse at
n = 9..12, and a stable +19 defect with one J-fixed cluster at n ≥ 13.
Characterizing exactly which algebraic relation drives this is **out
of scope** for the closure-quotient paper but worth recording as an
open follow-up.

The headline correction: **Pool A is not the canonical generic basin
the closure_quotient_landscape_findings §6 narrative assumed.** It is
generic enough to satisfy K_6^3 → 51 but exhibits a residual
non-genericity that surfaces only at larger n.

---

## 6. Implications for the closure-quotient paper

For §5.2 of the paper, the K_n^3 law is now:

> **K_n^3 single-side closure count = n(3n−1)/2 under truly generic
> ICs**, confirmed exactly for n = 6..20 across two independent
> strict-generic IC pools (Pool B and Pool C, total 30 distinct
> generic triples drawn from disjoint coefficient ranges). Pentagonal
> numbers: 12, 22, 35, 51, 70, 92, 117, 145, 176, 210, 247, 287, 330,
> 376, 425, 477, 532, 590 at n = 3..20.

For §3 (methods): the K_6^3 → 51 validity gate is necessary but not
sufficient at large n; an IC set passing the gate may still exhibit
hidden non-genericity at n ≥ 9. A stronger gate would be
K_10^3 → 145 or K_12^3 → 210 — these catch Pool A's hidden
relation that K_6^3 does not.

For §9 (IC genericity as master variable): the binary
generic/degenerate split needs refinement. There is a **genericity
hierarchy** — "small-n generic" (passes K_6^3=51) vs "fully generic"
(passes K_n^3 = pentagonal for all n). The canonical s58 set sits in
the gap, and this is itself a structural finding about IC sets for
the CFS model.

For the closure-v5 spec: any closure-v5 claim that depended on
Pool A at K_n^3 for n ≥ 9 should be flagged for re-validation.
Prior closure-v5 work that uses only K_6^3 (most of it) is unaffected.

---

## 7. Honest scope notes

- **Truly-generic-IC characterization is empirical, not algebraic.**
  Pools B and C are constructed by hand (strict-genericity: every
  component nonzero real and nonzero imag, no obvious small-integer
  relations). No proof that "every IC satisfying this property gives
  pentagonal" — only that two independent such pools agree at all
  measured n.
- **Range is bounded.** n = 6..20 is the measured range; the law
  is conjecturally pentagonal at all n ≥ 6, but the paper claims only
  what is measured.
- **The Pool A defect characterization (plateau at 19/1 by n = 13) is
  also empirical** and could shift if n were pushed further. Pool A
  has 15 entries so n = 15 is its ceiling without extension.

---

## 8. Open follow-up (out of scope for the closure-quotient paper)

Identify the specific algebraic relation in Pool A that produces the
n ≥ 9 defect. Candidate diagnostics:
- Probe Pool A subsets of size 9..12 to find a minimal degenerate
  subtuple.
- Check whether the defect is generated by Pool A's axis-aligned
  triples (the ones with zero imag in a component).
- Check whether the J-fixed cluster that appears at n = 10 has a
  specific form (e.g. a particular Gaussian-integer rep) that
  identifies the relation.

This is its own short investigation; not a paper blocker.
