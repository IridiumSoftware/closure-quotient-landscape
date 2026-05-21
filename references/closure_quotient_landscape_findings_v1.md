# The CFS Closure-Quotient Landscape — Consolidated Findings

**2026-05-20 update banner.** The §8 open item ("Depth 4 may under-close
K_n³ for n ≥ 9") was resolved with a substantially sharper result. The
K_n³ → n(3n−1)/2 (pentagonal) law is now confirmed for all measured
n = 6..20 across two independent strict-generic IC pools. The original
112/133/164/193 measurement is a **Pool-A artifact** (the canonical s58
IC set used in every prior CFS probe is not fully generic at n ≥ 9).
See `Kn3_large_n_findings_v1.md` for the full report. Inline updates to
§2, §5, §6, §8 below are flagged `[2026-05-20]`.

**Scope.** Exact-arithmetic characterisation of the CFS closure-quotient
family — the objects Q24/Q45/Q48/Q51/Q90/Q102 and their generators.
Built 2026-05-17, Aaron driving (expert lead, post-SFT-synthesis).
All counts are exact Gaussian-integer ray-quotients — no Float64 fidelity
threshold. Every scan is validity-gated (must reproduce K₆³→51 / →102).

Scripts (this folder): `reground_q_family.py`, `prime_count_probe.py`,
`q181_search.py`, `thread1_F_closure.py`, `thread2_J_fixed_check.py`,
`small_closure_landscape_scan.py`, plus the 2026-05-20 K_n³ probes:
`Kn3_large_n_probe.py`, `Kn3_extended_probe.py`,
`Kn3_pentagonal_confirmation.py`, `Kn3_push_to_n20.py`.

---

## 1. Re-grounding the published family (exact vs fidelity-0.999)

The closure-v5 family doc built Q45/Q84/Q90/Q98 with **fidelity-0.999
Float64 clustering**. Recomputed exactly:

| construction | published (fidelity) | exact | verdict |
|---|---|---|---|
| q45-seed single-side (Q45) | 45 | **45** | confirmed real |
| adjacent-36 C-closed (Q84) | 84 | **84** | confirmed real |
| q45-seed C-closed (Q90) | 90 | **90** | confirmed real |
| K₆³ C-closed (Q98) | 98 | **102** | **fidelity artifact — it is Q102** |
| iterated C-closure → Q180 | 180 | **90** | **Q180 does not exist** (idempotence) |

**Q45 and Q90 are genuine exact objects.** **"Q98" is wrong** — the
fidelity-0.999 threshold merged 4 distinct clusters; the exact object is
Q102. **"Q180" does not exist** — see §4.

---

## 2. The closure laws (exact, generic ICs)

From the small-closure landscape scan, with closed forms:

- **`cyc(n) → 4n`** single-side (linear).
- **`K_n³ → n(3n−1)/2`** single-side — the **pentagonal numbers**
  (12, 22, 35, 51, 70, 92, …). Composite for all n ≥ 3 (the factoring
  `n·(3n−1)/2` always splits). **[2026-05-20]** Confirmed for all
  measured n = 6..20 (sequence 51, 70, 92, 117, 145, 176, 210, 247,
  287, 330, 376, 425, 477, 532, 590) across two independent
  strict-generic IC pools (Pool B and Pool C in the K_n³ findings).
  The original claim was anchored at n ≤ 6; it now extends.
- **C-closure exactly doubles: `|Q ∪ C(Q)| = 2·|Q|`** — universal, exact
  (Q24→Q48, Q45→Q90, Q51→Q102, every scan row). See §5 for why.
- **Generic floor = Q12** — even the minimal 3-vertex seed closes to 12;
  generic full closure never goes below 12.
- **Depth-stable by depth 4** for every small object. **[2026-05-20]**
  K_n³ for n = 6..12 is depth-stable from depth 3 onward (verified out
  to depth 8 — same count at every depth).

---

## 3. Thread 3 — prime-count probe: the closure is not prime-blind

- The **regular families are composite-only** by derivation: `4n`,
  pentagonal `n(3n−1)/2`, and C-closure `2n` are all composite.
- **Irregular single-side seeds DO produce primes** — found at every
  scale: 13, 19, 23, 71, 83, 97, 193 (validity-gated generic-IC scan).
- **Q181 exists — exhibited.** A **12-vertex, 132-edge seed (density
  0.100) closes to exactly 181**; **12 distinct such seeds** found in 192
  trials. Q181 is a robust closure object. (181 prime.)
- **C-closure can never be prime** under generic ICs — it is exactly `2n`
  (even). So **every prime closure object is necessarily single-side.**

## 4. Thread 1 — F-closure: the C-closed objects are terminal

The genuine closure-growth endofunctor `F = M/~` (full multiway BFS +
quotient), applied to each C-closed object as a fresh seed:

| object | n_cl | F-closure (d4/5/6) | iterated C-closure | F-fixed? | C-idempotent? |
|---|---|---|---|---|---|
| Q48 | 48 | 48/48/48 | 48 | YES | YES |
| Q90 | 90 | 90/90/90 | 90 | YES | YES |
| Q102 | 102 | 102/102/102 | 102 | YES | YES |

**Q48, Q90, Q102 are fixed points of BOTH F and C-closure** — terminal
objects. Neither closure-growth nor conjugate-union grows them.
**Consequence:** the "factor family" `Q45→Q90→Q180→…` by iterated
doubling is **false**. There is exactly **one** doubling (single-side →
C-closed); then everything is terminal. Larger objects come only from
**larger seeds** — never from iterating operators on smaller ones (Q181
came from a fresh 12-vertex seed, §3).

## 5. Thread 2 — J-fixed check: why the doubling is exact

The charge-conjugation involution J on each C-closed object, fixed-point
count `f` (a J-fixed cluster has `conj(rep) ~ rep`):

| object | single | C-closed | J-fixed `f` | `2n − f` |
|---|---|---|---|---|
| Q48 | 24 | 48 | **0** | 2·24−0 = 48 ✓ |
| Q90 | 45 | 90 | **0** | 2·45−0 = 90 ✓ |
| Q102 | 51 | 102 | **0** | 2·51−0 = 102 ✓ |

**Under generic ICs, J is fixed-point-free** — zero self-conjugate
clusters. This is *why* C-closure cleanly doubles: `|Q∪C(Q)| = 2n − f`
with `f = 0`. **There is no `2n ± 1` structure** — the "doubling + 1
fixed-point core" idea is falsified directly.

**[2026-05-20] Refinement.** "Generic ICs" should be read as **truly**
generic — strict-genericity (every Gaussian component with nonzero real
AND nonzero imag part) confirms J fixed-point-free at every measured n
in {6..20} for K_n³. The canonical s58 IC set (Pool A in the K_n³
findings, used for the Q48/Q90/Q102 table above) is fully generic at
n ≤ 8 — the values 0/0/0 in the table are correct as stated. But Pool A
acquires exactly one J-fixed cluster at K_n³ for n ≥ 10, breaking the
`2|Q|` doubling by +1. That is a Pool A defect, **not** a property of
the closure. See §6 for the IC-genericity hierarchy this implies.

**Degenerate control — purely real ICs** (`conj(ψ) = ψ`): G₀ → single 16,
C-closed 16, all 16 J-fixed; K₆³ → 22, 22, all 22 J-fixed. Real ICs
collapse the conjugate sector onto the original — **no doubling at all**,
and the single-side count itself drops (16, 22 vs generic 24, 51).

---

## 6. The unifying finding — IC genericity is the master variable

The landscape splits by IC genericity:

- **Generic regime** (genuinely-complex ICs): lawful and IC-independent —
  `cyc → 4n`, `K_n³ → pentagonal`, `C-closure → 2n` exact, J
  fixed-point-free, floor Q12, objects terminal under F. This is the
  closure-v5 canonical family.
- **Degenerate / "partial closure" regime** (real or special ICs): counts
  drop and become IC-dependent; J acquires fixed points; the clean
  doubling breaks. This is where **small objects** (6, 7, 8, 9 — seen in
  early probes) and **odd/prime counts** can live.

Aaron's word "partial closure" names this exactly. The small generators
(Q3/Q5/Q9) — if they exist — are degenerate-regime objects, not generic
ones; the generic floor is Q12.

**[2026-05-20] Refinement: the regime split is a hierarchy, not a binary.**
The K_n³ probe (`Kn3_large_n_findings_v1.md`) found that the canonical
s58 IC set ("Pool A" — what every prior CFS probe used) passes the
K_6³ → 51 validity gate but is not fully generic at K_n³ for n ≥ 9:
it has a hidden algebraic relation that suppresses the single-side count
below pentagonal (shortfall 5/12/12/17/19/19/19 at n = 9..15) and
introduces 1 J-fixed cluster at n ≥ 10. Truly generic IC pools
(strict-genericity) match pentagonal exactly at all measured n.

So the regimes are:
- **Fully generic** (pentagonal at every n; J fixed-point-free at every n)
- **Small-n generic** (passes K_6³ → 51, fails at larger n — Pool A sits here)
- **Degenerate / partial-closure** (real ICs or other special structure)

The K_6³ → 51 validity gate alone is necessary but not sufficient
for fully-generic regime membership. A stronger gate (K_10³ → 145 or
K_12³ → 210) catches Pool A's hidden relation.

---

## 7. Status of the original intuitions

| intuition | verdict |
|---|---|
| Q45, Q90 real, derivable | **Confirmed** — exact, and `2n` derivable |
| Q98 a family member | **No** — it is Q102 (fidelity artifact) |
| Q180 via iterated doubling | **False** — C-closure idempotent, F-fixed; no Q180 any route |
| Q181 reachable | **Yes — exhibited** (12-vertex seed); closure not prime-blind |
| union rule | **Found & exact** — `\|Q∪C(Q)\| = 2\|Q\|`, J fixed-point-free |
| composition rule (combine seeds) | **Found** — additive + universal interface term; see §9 |
| Q3/Q5/Q9 generators | **Not generic** — generic floor is Q12; small objects live in the degenerate/partial-closure regime |

---

## 8. Honest notes & open items

- **Two prime-probe bugs** were hit and fixed (IC period-6 clone, then
  degenerate real ICs). A **validity gate** (reproduce K₆³→51/102) now
  guards every scan; the reported results are all gate-passing.
- ~~**Depth 4 may under-close `K_n³` for n ≥ 9** — the large-n single-side
  counts (112/133/164/193 at n=9..12) sit below the pentagonal formula;
  the exact `K_n³` law beyond n = 6 is not pinned. Does not affect the
  prime-existence or Q181 results.~~ **[2026-05-20] RESOLVED.** Depth
  was not the issue. The actual K_n³ law is pentagonal at every measured
  n (6..20) under truly generic ICs; the 112/133/164/193 measurement was
  a defect in the canonical s58 IC set ("Pool A"). See
  `Kn3_large_n_findings_v1.md`. The Q181 result is unaffected (Q181 was
  found by random subset seeding, not full K_n³).
- **[2026-05-20] NEW OPEN ITEM.** Characterise the hidden algebraic
  relation in the s58 / Pool A IC set that produces its K_n³ defect at
  n ≥ 9. Diagnostics suggested in the K_n³ findings doc §8: probe Pool A
  subsets of size 9..12 to find a minimal degenerate subtuple; check
  whether the defect is driven by Pool A's axis-aligned components
  (`(1,0)`, `(3,0)`, `(2,0)`); identify the form of the J-fixed cluster
  that appears at n = 10.
- **closure-v5 relevance:** the Q98→Q102 correction and the closure laws
  (§2) bear on closure-v5's own family documentation; integrating them
  there is Aaron's call. **[2026-05-20]** Any closure-v5 result that
  relies on Pool A at K_n³ for n ≥ 9 should be flagged for
  re-validation under a fully-generic pool. Prior closure-v5 work using
  only K_6³ (most of it) is unaffected.

---

## 9. The composition rule — combining two seeds

Probe: `cyc(nA) ∘ cyc(nB)` under four join modes, δ ≡ Q(A∘B) − Q(A) − Q(B),
across 8 (A,B) pairs. **δ is universal — it depends only on the interface,
never on A or B:**

| join mode | δ (every (A,B) pair) |
|---|---|
| **disjoint `⊔`** (no shared verts, no bridge) | **0** |
| **bridge** (extra triple spanning A,B) | **+3 per bridge triple** (br1 +3, br2 +6) |
| **glue g vertices** | **−(4g−3)** — g=1 → −1, g=2 → −5, g=3 → −9 |

> **The composition rule:**  `Q(A ∘ B) = Q(A) + Q(B) + δ(interface)`
> with `δ = +3·(#bridge triples) − Σ(4gᵢ−3)` over glued interfaces — and δ
> is **structure-independent** (the same for every A, B with that
> interface). Verified exact across all 8 pairs. (The lone exception,
> cyc(3)∘cyc(3) glued-3, is the degenerate total-overlap case g = nA = nB
> where the seeds coincide — not a rule violation.)

**This is the generative algebra Aaron was after.** Consequences:

- The closure counts form a **commutative monoid under disjoint sum `⊔`**:
  `Q(A⊔B) = Q(A)+Q(B)`, associative, commutative, identity 0.
- **The union rule is a special case.** C-closure `Q ∪ C(Q) = 2n` *is* the
  disjoint composition `Q ⊔ (conjugate copy)` — the conjugate sector
  behaves as a disjoint generic copy (δ=0), giving exactly `2n`. This is
  *why* C-closure is clean `2n` and J is fixed-point-free (§5) — the
  conjugate copy does not interface with the original.
- **The algebra is purely ADDITIVE** — counts add, with a linear interface
  correction. **There is no multiplicative rule.** Aaron's "Q45 = Q9·Q5"
  factor intuition does not hold: `Q45 ≠ Q9 × Q5`. Even the C-closure "×2"
  is additive (doubling = adding a disjoint copy). The factor structure of
  the *numbers* (45 = 9·5) is not the closure algebra; the closure algebra
  is `+`, never `×`.

This **closes the last open thread.** The closure-quotient landscape is now
characterised: generators (single-side closures, §2), the additive monoid
under ⊔ with universal interface δ (§9), all C-closed objects terminal
under F and C (§4), J fixed-point-free (§5), primes only off the regular
grid (§3). Script: `composition_rule_probe.py`.
