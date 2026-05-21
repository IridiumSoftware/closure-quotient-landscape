#!/usr/bin/env python3
# Extend Pool B and Pool C to n=20, depth=5, K_n^3 single-side only.
# Goal: confirm pentagonal n(3n-1)/2 holds beyond n=15.
import time, sys
sys.path.insert(0, ".")
from Kn3_pentagonal_confirmation import (
    POOL_B, POOL_C, count_single, Kn, pentagonal,
)

print(f"|POOL_B|={len(POOL_B)}, |POOL_C|={len(POOL_C)}")
max_n = min(len(POOL_B), len(POOL_C), 20)

print()
print(f"{'n':>3}  {'|K_n^3|':>8}  {'pent':>6}  {'B_ss':>6}  {'C_ss':>6}  "
      f"{'pent==B==C':>11}  {'t_B':>7}  {'t_C':>7}")
print("-"*80)
for n in range(15, max_n+1):
    topo = Kn(n)
    pent = pentagonal(n)
    t0 = time.time(); cB = count_single(topo, POOL_B, n, 5); tB = time.time()-t0
    t0 = time.time(); cC = count_single(topo, POOL_C, n, 5); tC = time.time()-t0
    ok = "YES" if (cB == cC == pent) else "NO"
    print(f"{n:>3}  {len(topo):>8}  {pent:>6}  {cB:>6}  {cC:>6}  "
          f"{ok:>11}  {tB:>7.2f}  {tC:>7.2f}")
