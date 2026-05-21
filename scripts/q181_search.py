#!/usr/bin/env python3
# Targeted search — exhibit a seed whose exact single-side closure = 181.
# Reuses the validity-gated generic-IC machinery from sft_prime_probe.py.
import random, importlib.util, sys
spec=importlib.util.spec_from_file_location("pp","/tmp/sft_prime_probe.py")
# import functions without running the scan: load module guarded
src=open("/tmp/sft_prime_probe.py").read().split("random.seed(17)")[0]
ns={}
exec(src, ns)
count_single=ns['count_single']; Kn=ns['Kn']; isprime=ns['isprime']

# gate
assert count_single(Kn(6),6,depth=5)==51, "gate fail"
print("gate K6^3=51 PASS")

random.seed(2718)
hits=[]; best=None
full12=Kn(12); full13=Kn(13)
trials=0
for (n,full) in ((12,full12),(13,full13)):
    me=len(full)
    for frac in [x/200 for x in range(8,40)]:   # fine density sweep
        k=max(1,int(me*frac))
        for t in range(6):
            topo=random.sample(full,k)
            c=count_single(topo,n,depth=5)
            trials+=1
            if c==181:
                hits.append((n,k,frac));
            if best is None or abs(c-181)<abs(best[0]-181):
                best=(c,n,k)
    if hits: break

print(f"trials: {trials}")
if hits:
    n,k,frac=hits[0]
    print(f"*** Q181 EXHIBITED — n={n} vertices, {k}-edge seed (density {frac:.3f}) "
          f"closes to exactly 181 ***")
    print(f"   total Q181 seeds found: {len(hits)}")
else:
    c,n,k=best
    print(f"no exact 181 in {trials} trials; closest: count={c} at n={n}, {k} edges")
    print(f"   (181 is prime and in-band — n=12 spans 180s; a finer/longer search pins it)")
