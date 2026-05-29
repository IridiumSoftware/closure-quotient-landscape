Q181_SEED = [
    (0, 2, 1),
    (0, 4, 9),
    (0, 4, 10),
    (0, 8, 3),
    (0, 9, 1),
    (0, 9, 7),
    (0, 9, 8),
    (0, 10, 4),
    (0, 11, 5),
    (1, 0, 9),
    (1, 2, 10),
    (1, 4, 7),
    (1, 6, 4),
    (1, 6, 8),
    (1, 7, 4),
    (1, 7, 5),
    (1, 8, 0),
    (1, 10, 4),
    (2, 0, 6),
    (2, 1, 9),
    (2, 3, 1),
    (2, 3, 5),
    (2, 3, 6),
    (2, 6, 1),
    (2, 6, 10),
    (2, 7, 1),
    (2, 9, 1),
    (2, 10, 5),
    (3, 0, 10),
    (3, 0, 11),
    (3, 2, 10),
    (3, 4, 5),
    (3, 4, 7),
    (3, 5, 6),
    (3, 7, 5),
    (3, 10, 1),
    (4, 0, 5),
    (4, 0, 8),
    (4, 1, 5),
    (4, 1, 10),
    (4, 1, 11),
    (4, 2, 9),
    (4, 3, 1),
    (4, 5, 2),
    (4, 5, 9),
    (4, 7, 8),
    (4, 10, 11),
    (4, 11, 3),
    (4, 11, 7),
    (5, 1, 4),
    (5, 2, 6),
    (5, 2, 7),
    (5, 3, 4),
    (5, 4, 10),
    (5, 6, 4),
    (5, 8, 10),
    (5, 8, 11),
    (5, 9, 6),
    (5, 9, 8),
    (5, 9, 11),
    (5, 11, 0),
    (6, 0, 5),
    (6, 0, 8),
    (6, 0, 9),
    (6, 2, 0),
    (6, 2, 7),
    (6, 4, 9),
    (6, 5, 1),
    (6, 7, 0),
    (6, 7, 3),
    (6, 7, 4),
    (6, 9, 7),
    (6, 10, 7),
    (7, 0, 3),
    (7, 2, 1),
    (7, 2, 9),
    (7, 2, 11),
    (7, 3, 4),
    (7, 5, 10),
    (7, 6, 9),
    (7, 6, 10),
    (7, 6, 11),
    (7, 10, 0),
    (7, 11, 1),
    (8, 1, 3),
    (8, 1, 5),
    (8, 1, 9),
    (8, 2, 3),
    (8, 4, 7),
    (8, 4, 9),
    (8, 5, 4),
    (8, 6, 1),
    (8, 6, 7),
    (8, 7, 9),
    (8, 9, 10),
    (8, 11, 9),
    (9, 0, 5),
    (9, 0, 8),
    (9, 1, 3),
    (9, 2, 5),
    (9, 3, 0),
    (9, 4, 10),
    (9, 5, 7),
    (9, 5, 10),
    (9, 5, 11),
    (9, 6, 7),
    (9, 7, 10),
    (9, 8, 0),
    (9, 10, 1),
    (9, 11, 7),
    (10, 1, 4),
    (10, 1, 8),
    (10, 1, 11),
    (10, 2, 7),
    (10, 2, 9),
    (10, 3, 1),
    (10, 3, 6),
    (10, 7, 2),
    (10, 8, 0),
    (10, 8, 2),
    (10, 9, 11),
    (10, 11, 1),
    (10, 11, 6),
    (11, 1, 8),
    (11, 2, 6),
    (11, 3, 9),
    (11, 4, 2),
    (11, 4, 3),
    (11, 5, 2),
    (11, 6, 0),
    (11, 9, 10),
    (11, 10, 6),
]

Q181_N_VERTS = 12
# Discovered by the random density sweep in q181_search.py (random.seed(2718),
# density 0.10, trial 0). Baked here as a literal so any consumer can verify
# |Q_181| = 181 without re-running the search. The discovery method is the
# search; the existence claim rests on this fixed edge list.

if __name__ == "__main__":
    import sys, os
    # Honor CFS_REPO_ROOT like the other env-gated scripts (probe_*_v1,
    # q51_exact_*): point it at the scripts/quotient_construction/ directory of
    # a Closure-Forces-Structure checkout (see README). The Q_181 existence
    # claim rests on the Q181_SEED literal above; this __main__ re-verifies it
    # against the canonical exact builder.
    CFS_ROOT = os.environ.get("CFS_REPO_ROOT")
    if not CFS_ROOT:
        sys.exit("Set CFS_REPO_ROOT to a Closure-Forces-Structure checkout's "
                 "scripts/quotient_construction/ to verify (see README).")
    sys.path.insert(0, CFS_ROOT)
    try:
        from q102_build_exact_v1 import (
            build_multiway, proj_equiv, is_zero_vec, _GEN_ICS
        )
    except ImportError as e:
        sys.exit(f"Could not import q102_build_exact_v1 from "
                 f"CFS_REPO_ROOT={CFS_ROOT}: {e}")
    ics = {i: [list(t) for t in _GEN_ICS[i]] for i in range(Q181_N_VERTS)}
    psi, _, _ = build_multiway(Q181_SEED, ics, 5)
    reps = []
    for pv in psi.values():
        if is_zero_vec(pv):
            continue
        if not any(proj_equiv(pv, r) for r in reps):
            reps.append(pv)
    n = len(reps)
    ok = "PASS" if n == 181 else "FAIL"
    print(f"Q_181 canonical seed: |seed|={len(Q181_SEED)}, n_verts={Q181_N_VERTS}")
    print(f"|Q| at depth 5 = {n}  expect 181  [{ok}]")
    sys.exit(0 if ok == "PASS" else 1)
