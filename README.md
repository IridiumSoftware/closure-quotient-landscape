# The CFS Closure-Quotient Landscape

> *Exact Arithmetic, Validity Gating, and a Generative Algebra*
>
> Aaron Green — 2026

The rendered paper is **[`closure_quotient_landscape.pdf`](closure_quotient_landscape.pdf)**
(co-authored Aaron Green + Brian Crabtree). LaTeX source is
`closure_quotient_landscape.tex`. Repository is the canonical
self-contained build artefact: TeX source, bibliography (inline
`\thebibliography`), figures, cross-substrate Galois-pair diagnostic
artifacts (`cross_substrate_galois_pair/`), and cited reference
documents all live here.

This is **v2**: the original closure-quotient landscape (Part 1) merged
with a 9-substrate discrete Hodge corpus (Part 2) demonstrating the
cyclotomic-Galois pattern across the Platonic family, the torus
$T^2_{3,3}$, $\mathbb{RP}^2$, the Klein bottle, the heptagonal
bipyramid (cubic extension), and the soccer ball. The v1 (closure-
quotient only) paper is preserved in git history.

## Headline

A *generative algebra* on the CFS closure-quotient family: for closure
quotients $A$ and $B$,
$$
\lvert Q(A \circ B) \rvert \;=\; \lvert Q(A) \rvert + \lvert Q(B) \rvert + \delta(\text{interface})
$$
with $\delta$ a structure-independent interface correction
($\delta_\sqcup = 0$ disjoint, $\delta_\text{bridge} = +3$ per bridge
triple, $\delta_\text{glue-}g = -(4g-3)$). The counts form a
commutative monoid under disjoint sum; the C-closure doubling rule
falls out as the special case where the conjugate sector is
structurally disjoint. See §10 of the paper for the headline; §4
re-grounds the published family (Q98 → Q102, Q180 nonexistent); §5
the closure laws; §8 exhibits $Q_{181}$; §9 the IC-genericity
hierarchy.

## Repository layout

```
closure_quotient_landscape.tex   — LaTeX source (single file, inline bibliography)
closure_quotient_landscape.pdf   — rendered output (regenerable via `make`)
Makefile                         — `latexmk -pdf` build + `verify-scripts` target
figures/                         — auto-generated TikZ figures (one .tex per quotient)
references/                      — canonical static copies of cited findings docs
scripts/                         — probe scripts (Python stdlib only for the core,
                                   plus an env-var-gated probe set that imports the
                                   CFS exact-arithmetic builders — see below)
cross_substrate_galois_pair/     — Galois-pair diagnostic run outputs cited in §28
haskell/                         — independent Haskell port of the rewriting kernel
                                   (cabal, type-checked binary-op/ternary-hyperedge
                                   distinction, 24/24 cross-impl cardinality + J +
                                   F-closure checks pass; see haskell/README.md)
LICENSE                          — CC-BY-4.0
TCL.txt + TCL.txt.asc            — Triadic Closure License v1.3 (commons substrate)
```

## Build

```sh
make                              # latexmk -pdf (preferred)
# or, if no latexmk:
tectonic closure_quotient_landscape.tex
```

Both engines produce identical PDF output.

## Reproducibility

### Cross-implementation validation (Haskell)

An independent Haskell port of the rewriting kernel lives in `haskell/`.
It encodes the operator/hyperedge arity distinction at the type level
(`compose :: Cochain -> Cochain -> Cochain` is binary by signature;
`type Hyperedge = (Vertex, Vertex, Vertex)` is ternary by construction)
and matches the Python on every canonical cardinality:

```sh
cd haskell
cabal test --test-show-details=direct
```

24 of 24 checks pass: |Q_24|, |Q_42|, |Q_45|, |Q_51| (single-side);
|Q_48|, |Q_84|, |Q_90|, |Q_102| (C-closed); J fixed-point-free
involution; F-closure terminality at depths 4, 5, 6 for the four
doubly-terminal members.

### Python probe corpus

The probe corpus is in `scripts/` (CPython 3.10+). It splits into two
groups.

**Self-contained core** — pure Python standard library, no
`CFS_REPO_ROOT`, no external builders. These run straight from a fresh
clone (the visualization generators and `q181_search.py` import only
sibling modules in the same directory); every measurement probe runs a
`K_6^3 → 51` validity gate before reporting:

```sh
cd scripts/
for s in *.py; do python3 "$s" || echo "FAIL: $s"; done
```

On a fresh clone the self-contained scripts all pass. The env-gated
group below instead prints a one-line `Set CFS_REPO_ROOT` notice and
returns non-zero — a deliberate clean exit, which the loop above reports
as `FAIL` only because of the exit code. Core-group runtime is under
five minutes on a single CPU core. The figure `.tex` files in `figures/`
are pre-generated; the visualization generators only *regenerate* them.

**Env-gated** — the `probe_*_v1.py` set, `q51_exact_phase_s_python.py`,
`q51_exact_vs_fidelity_verification.py`, and the figure regenerator
`dump_quotient_to_tikz.py` import the canonical exact-arithmetic builders
that live in the broader CFS programme working tree (these also use
`numpy`). Those builders are published in the public CFS mirror: point
`CFS_REPO_ROOT` at the `scripts/quotient_construction/` directory of a
[Closure-Forces-Structure---SM-Rosen-Hypergraphs](https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs)
checkout pinned at commit
[`a2c1c0c`](https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs/commit/a2c1c0c2103c40f83599bb91f193bd78be12b2cd)
— the commit that adds `q102_build_exact_v1.py` and
`g4_larger_graphs_exact_v1.py`:

```sh
git clone https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs cfs
git -C cfs checkout a2c1c0c
export CFS_REPO_ROOT="$PWD/cfs/scripts/quotient_construction"
```

Without `CFS_REPO_ROOT` set they exit cleanly with a clear message. The
captured stdout from each run is committed alongside the script as
`*.stdout.log` (and the figures as `figures/*.tex`) so the cited evidence
survives without re-running. Verified at the pinned commit:
`q51_exact_phase_s_python.py` reproduces the committed Galois-pair result
(`|C¹|=315`, `{3∓√3}` at multiplicity 30, `[4,2]⊕[3,3]⊕[3,2,1]`).

## Citation

```bibtex
@unpublished{Green2026ClosureQuotient,
  author = {Aaron Green},
  title  = {The {CFS} Closure-Quotient Landscape:
            Exact Arithmetic, Validity Gating, and a Generative Algebra},
  year   = {2026},
  note   = {researchers.one preprint;
            \url{github.com/IridiumSoftware/closure-quotient-landscape}},
}
```

## Related projects

- **`closure-forces-structure`** — the parent CFS programme (Rosen
  closure → Standard-Model algebra; this paper's substrate). Public
  mirror at
  [`Closure-Forces-Structure---SM-Rosen-Hypergraphs`](https://github.com/IridiumSoftware/Closure-Forces-Structure---SM-Rosen-Hypergraphs).
- **`triadic-coordination-engine`** — the engine track that the
  paper's OpenGL visualization work (Phase 0/1 Q24/Q102 in 3D) lives
  under.

## License

The paper itself is **CC-BY-4.0** — see [`LICENSE`](LICENSE). Free to
share and adapt with attribution.

The broader research programme this repository sits inside is released
to the commons under the **Triadic Closure License (TCL) v1.3** — see
[`TCL.txt`](TCL.txt) (GPG-signed detached signature at
[`TCL.txt.asc`](TCL.txt.asc)). The TCL is the commons-substrate
license; it governs derived works that extend the programme rather
than just cite the paper. Canonical TCL URL:
<https://github.com/IridiumSoftware/possibilistic-security>.
