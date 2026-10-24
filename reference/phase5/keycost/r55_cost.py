# ABOUTME: R5.5 cost check: per stratum of the sampled candidates, V×τ, slice law size, codegen
# ABOUTME: operations, ns per call, compile and key time. Args: data/kt_r55.tsv
import sys
from collections import defaultdict

def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]

by = defaultdict(list)
for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if f[7] != "ok":
        by[f[0] + " ERR"].append(f)
        continue
    by[f[0]].append(dict(cx=int(f[3]) + int(f[4]), cxA=int(f[3]), D=int(f[8]), KT=int(f[9]),
                         np=int(f[10]), tlaw=1e3 * float(f[11]), tkey=1e3 * float(f[12]),
                         tcomp=float(f[19]), ns=float(f[20]), al=int(f[22]), ops=int(f[23]),
                         n=int(f[17])))
print("stratum: n | V×τ sum med [min, max] | D med max | ops med max | ns med p90 max | "
      "compile s med max | key ms med max | np med | allocs>0")
for k in ("plain_near", "allo_under", "allo_new", "plain_over"):
    rs = by[k]
    g = lambda key: [r[key] for r in rs]
    print(f"{k}: {len(rs)} | {q(g('cx'), .5)} [{min(g('cx'))}, {max(g('cx'))}] | {q(g('D'), .5)} "
          f"{max(g('D'))} | {q(g('ops'), .5)} {max(g('ops'))} | {q(g('ns'), .5):.1f} "
          f"{q(g('ns'), .9):.1f} {max(g('ns')):.1f} | {q(g('tcomp'), .5):.3f} {max(g('tcomp')):.3f} | "
          f"{q(g('tkey'), .5):.2f} {max(g('tkey')):.2f} | {q(g('np'), .5)} | "
          f"{sum(1 for r in rs if r['al'] > 0)}")
for k in by:
    if k.endswith("ERR"):
        print(k, len(by[k]), by[k][0][7][:120])
