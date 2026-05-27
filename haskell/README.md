# cfs-kernel — Haskell reference port of the CFS rewriting kernel

A hermetic, type-checked re-implementation of the closure-quotient
rewriting kernel that the closure-quotient-landscape paper rests on.

## Why this exists

The Python reference (`q102_build_exact_v1.py` at the closure-v5 root) is
authoritative for the empirical results but loose at the type level. An
LLM-swarm review of the manuscript flagged the language as conflating
"ternary hypergraph" with "ternary operator." This module settles that
distinction at the type level:

- `CFS.Cochain.compose :: Cochain -> Cochain -> Cochain` — the rewriting
  /operator/ is binary by signature (conjugated cross-product on `Z[i]^3`).
- `CFS.Multiway.Hyperedge = (Vertex, Vertex, Vertex)` — the rewriting
  /graph/ is genuinely ternary; each step takes a triple `(v1, v2, v3)`,
  applies the binary operator to `(psi v1, psi v2)`, and schedules three
  new hyperedges that propagate `v3`. The third slot is not an operator
  input but is load-bearing for the multiway expansion (and hence for
  `|Q|`).

Re-implementing in Haskell with no shared code with the Python and
matching its cardinalities is the strongest cross-validation: agreement
on `|Q_24|`, `|Q_48|`, `|Q_51|`, `|Q_84|`, `|Q_90|`, `|Q_102|` plus
J-doubling and F-closure terminality means the rewriting rule is
captured correctly in both implementations.

## What's in scope

- The rewriting kernel: Gaussian-integer arithmetic, conjugated cross-
  product, exact projective ray equivalence, ternary-hyperedge multiway
  expansion, exact-ray cluster reduction.
- The canonical Gaussian-integer ICs (the 15-element `_GEN_ICS` set
  verified to put `K_6^3 -> Q_102`).
- The canonical seed topologies (`G_0`, adjacent, complete-ternary,
  adjacent + 6 numpy-PCG64-seed-42 pads → Q_45 / Q_90).
- J: the charge-conjugation involution on Q.
- F-closure terminality: re-running the multiway on `Q`'s clusters and
  hyperedges at depths 4, 5, 6 to confirm `|F^d Q| = |Q|`.

## What's out of scope

- Hodge L_1 spectra and the cyclotomic-Galois decomposition — see the
  paper's Python/Julia drivers.
- Visualization figure generators — see `paper/scripts/dump_quotient_to_tikz.py`.
- Q_181 — its seed is found by random density sweep, not a pinned topology;
  not part of this kernel scope.

## Build + test

```sh
cabal build
cabal test --test-show-details=direct
```

Pinned to GHC 9.6.7; `cabal.project.freeze` pins the dependency graph
(`base`, `containers` only).

## Cross-implementation verification (current state)

Every check below is run against the canonical `_GEN_ICS` initial
conditions, depth 5, with both the Python reference and this Haskell
port. Two independent implementations producing identical numbers is
the load-bearing fact.

| Check | Expected | Status |
|---|---:|---|
| `|Q_24|` (G_0 single-side) | 24 | PASS |
| `|Q_42|` (adjacent single-side) | 42 | PASS |
| `|Q_45|` (adj+6pads single-side) | 45 | PASS |
| `|Q_51|` (K_6^3 single-side) | 51 | PASS |
| `|Q_48|` (G_0 C-closed) | 48 | PASS |
| `|Q_84|` (adjacent C-closed) | 84 | PASS |
| `|Q_90|` (adj+6pads C-closed) | 90 | PASS |
| `|Q_102|` (K_6^3 C-closed) | 102 | PASS |
| J total + involution + fixed-point-free on Q_{48,84,90,102} | all | PASS |
| F-closure terminality `|F^d Q_{48,84,90,102}|` at d∈{4,5,6} | unchanged | PASS |

24 / 24 checks pass.

## Evidence type

- `type-checked`: the rewriting kernel's operator arity (binary) and
  hyperedge arity (ternary) are encoded in the types and verified by GHC.
- `example-tested` (cross-implementation): the cardinality agreement is
  established by running two independent implementations and comparing.

In the closure-v5 spec hierarchy this combination is sufficient for
`:proved` status on the structural-arity claim, and `:verified` on the
cardinality claim.
