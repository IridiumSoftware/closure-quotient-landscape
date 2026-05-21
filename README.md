# The CFS Closure-Quotient Landscape

> *Exact Arithmetic, Validity Gating, and a Generative Algebra*
>
> Aaron Green — 2026

The rendered paper is **[`closure_quotient_landscape.pdf`](closure_quotient_landscape.pdf)**
(~50 pp, 12 pt single-spaced). LaTeX source is
`closure_quotient_landscape.tex`. Repository is the canonical
self-contained build artefact: TeX source, bibliography (inline
`\thebibliography`), figures, and cited reference documents all live
here.

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
RONE_SUBMISSION.md               — submission metadata for researchers.one
figures/                         — auto-generated TikZ figures (one .tex per quotient)
references/                      — canonical static copies of cited findings docs
scripts/                         — probe scripts (Python stdlib only, no installs)
LICENSE                          — CC-BY-4.0
```

## Build

```sh
make                              # latexmk -pdf (preferred)
# or, if no latexmk:
tectonic closure_quotient_landscape.tex
```

Both engines produce identical PDF output.

## Reproducibility

The probe corpus is in `scripts/`. All scripts use only the Python
standard library (no `numpy`, `matplotlib`, etc.) and run on
CPython 3.10+. Every probe runs a `K_6^3 → 51` validity gate before
reporting measurements.

```sh
cd scripts/
for s in *.py; do python3 "$s" || echo "FAIL: $s"; done
```

Total runtime under five minutes on a single CPU core.

Some scripts (`dump_q24_to_c.py`, `dump_q102_to_c.py`,
`dump_quotient_to_tikz.py`) import canonical exact-arithmetic builders
that live in the broader CFS programme working tree
(`github.com/IridiumSoftware/closure-forces-structure` and its public
mirror); the figure `.tex` files in `figures/` are pre-generated, so
running them is only needed to *regenerate* the figures from updated
data.

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

CC-BY-4.0 — see [`LICENSE`](LICENSE). Free to share and adapt with
attribution.
