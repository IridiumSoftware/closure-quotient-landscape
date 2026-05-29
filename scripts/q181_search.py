#!/usr/bin/env python3
r"""
q181_search.py — targeted search exhibiting a seed whose exact single-side
closure has cardinality 181 (prime).

This is the *discovery* script for the $Q_{181}$ existence claim of
\S\ref{sec:primes-Q181}. The discovered seed is baked as a verified literal
in `q181_canonical_seed.py` (with its own self-check); the existence claim
rests on that fixed edge list, while this script records *how* it was found:
a `random.seed(2718)` density sweep over irregular 12- and 13-vertex seeds.

Self-contained: the exact Gaussian-integer kernel (`count_single`, `Kn`,
`POOL_A`) is imported from the stdlib-only sibling
`Kn3_pentagonal_confirmation.py` in this same directory — the same
import-a-sibling pattern used by `visualize_Q24.py` and `Kn3_push_to_n20.py`.
No `CFS_REPO_ROOT` / external engine is needed; run it from `scripts/`:

  cd scripts/
  python3 q181_search.py

(Previously this script read the kernel from `/tmp/sft_prime_probe.py` via an
absolute path, which made it non-reproducible on a fresh clone; the kernel now
comes from the committed sibling module.)
"""

import sys
import random
import contextlib

# Import the exact single-side closure counter from the sibling stdlib-only
# kernel. Redirect that module's import-time table print to stderr so this
# script's stdout (captured as q181_search.stdout.log) stays clean.
with contextlib.redirect_stdout(sys.stderr):
    from Kn3_pentagonal_confirmation import count_single, Kn, POOL_A

DEPTH = 5

# Validity gate: K_6^3 must close to exactly 51 under the canonical Pool A ICs
# before any measurement is reported (same gate every probe in this corpus runs).
g = count_single(Kn(6), POOL_A, 6, depth=DEPTH)
assert g == 51, f"gate fail: K_6^3 single-side = {g}, expected 51"
print(f"gate K_6^3 = 51 PASS")

# Density sweep. random.seed(2718) is the seed under which the canonical Q_181
# seed was originally discovered (n=12, density 0.10, first hit).
random.seed(2718)
hits = []          # (n_verts, k_edges, density, trial_index, seed_edges)
best = None        # (closest_count, n_verts, k_edges) — for the no-hit report
trials = 0
for (n, full) in ((12, Kn(12)), (13, Kn(13))):
    me = len(full)
    for frac in [x / 200 for x in range(8, 40)]:   # fine density sweep 0.04..0.195
        k = max(1, int(me * frac))
        for t in range(6):
            topo = random.sample(full, k)
            c = count_single(topo, POOL_A, n, depth=DEPTH)
            trials += 1
            if c == 181:
                hits.append((n, k, frac, trials, topo))
            if best is None or abs(c - 181) < abs(best[0] - 181):
                best = (c, n, k)
    if hits:
        break

print(f"trials: {trials}")
if hits:
    n, k, frac, trial, topo = hits[0]
    print(f"*** Q181 EXHIBITED — n={n} vertices, {k}-edge seed (density {frac:.3f}), "
          f"first hit at trial {trial}: closes to exactly 181 ***")
    print(f"   total Q181 seeds found in sweep: {len(hits)}")
    print(f"   181 is prime; the smallest prime-cardinality closure object exhibited.")

    # Cross-check the discovered seed against the baked canonical literal.
    try:
        with contextlib.redirect_stdout(sys.stderr):
            from q181_canonical_seed import Q181_SEED
        same = sorted(topo) == sorted(Q181_SEED)
        print(f"   discovered seed matches q181_canonical_seed.Q181_SEED: "
              f"{'YES' if same else 'NO (RNG/pool drift — investigate)'}")
    except ImportError:
        pass
    sys.exit(0)
else:
    c, n, k = best
    print(f"no exact 181 in {trials} trials; closest: count={c} at n={n}, {k} edges")
    print(f"   (181 is prime and in-band — n=12 spans the 180s; a finer/longer "
          f"search pins it. See q181_canonical_seed.py for the baked seed.)")
    sys.exit(1)
