# Cross-Substrate Galois Pair — Diagnostic Artifacts

Companion artifacts cited by the main paper
(`closure_quotient_landscape.tex`,
bibitem `tce-substrate-provenance`) documenting that the discrete
Hodge 1-Laplacian on two independently-constructed 2-complexes —
the CFS $K_n^3$ quotient $Q_{51}$ and the triangulated $3\times 3$
flat torus $T^2_{3,3}$ — both contain the eigenvalue pair
$\{3 - \sqrt 3,\, 3 + \sqrt 3\}$ (roots of $x^2 - 6x + 6$).
The substrate-provenance discipline forbids upgrading the
appearance to a derivation claim across substrates; it is an
algebraic observation indexed by the $\bZ_3$-cyclotomic
substructure in both automorphism groups.

## Artifacts

- `q51_exact_vs_fidelity_run.txt` — four-build comparison
  (exact Gaussian-integer ray-equivalence vs.\ Float64
  fidelity-$0.999$ clustering; canonical-simplicial vs.\
  ordered-tuple chain conventions). All four agree on the
  $\{3 \mp \sqrt 3\}$ pair at multiplicity 30.
- `q51_exact_phase_s_run.txt` — $S_6$ character-theoretic
  decomposition on the exact build:
  $\sigma_E L_1 = L_1 \sigma_E$ at machine zero on all 11
  conjugacy classes; per-eigenspace decomposition
  $[4,2] \oplus [3,3] \oplus [3,2,1]$ on the $\{3 \mp \sqrt 3\}$
  eigenspaces.

## Drivers

- `scripts/q51_exact_vs_fidelity_verification.py` — four-build
  verifier.
- `scripts/q51_exact_phase_s_python.py` — Phase S decomposition
  on the exact build.

Both depend only on the Python standard library plus `numpy` and
`scipy`. Tested under Python 3.13.
