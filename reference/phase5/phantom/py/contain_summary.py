# ABOUTME: Tabulates the numeric containment sample by stratum: verdicts at 1e-12 (variant inside the
# ABOUTME: comparator; comparator inside the variant at fixed L or with L free) and the worst residuals.
import collections, sys
sys.path.insert(0, ".")
from load import *
FS = {rxn: facts(rxn) for rxn in ("bibi", "pingpong", "uniuni_R", "bibi_R")}
tags = {}
for l in open(D + "/sel_tags.tsv"):
    rxn, g, cmp, tag, n = l.rstrip("\n").split("\t"); tags[(rxn, g, cmp)] = (tag, n)
TH = 1e-12
def verdict(ab, bafix, bafree):
    a, b, c = ab < TH, bafix < TH, bafree < TH
    return ("same at fixed L" if a and b else "same, comparator reached as L->0" if a and c else
            "variant strictly inside" if a else "variant contains comparator" if (b or c) else "differs")
res = collections.defaultdict(list)
for l in open(D + "/contain.tsv"):
    x = l.rstrip("\n").split("\t")
    res[tags[tuple(x[:3])]].append((float(x[3]), float(x[4]), float(x[5]), x[0], x[1]))
print("stratum (population size) | pairs | verdicts | worst variant-in-comparator | worst comparator-in-variant (fixed L)")
for (t, n) in sorted(res):
    rs = res[(t, n)]
    c = collections.Counter(verdict(*r[:3]) for r in rs)
    same = [r for r in rs if r[0] < TH]
    print(f"{t} ({n}) | {len(rs)} | {dict(c)} | max over inside: {max((r[0] for r in same), default=float('nan')):.1e} | "
          f"max over fixed-L same: {max((r[1] for r in rs if r[0] < TH and r[1] < TH), default=float('nan')):.1e}")
    for r in rs:
        v = verdict(*r[:3])
        if v not in ("same at fixed L",):
            fl = FS[r[3]][int(r[4])]
            print(f"    {r[3]} #{r[4]} {fl['move']} flag={fl['flag']} cert={fl['cert']} r_var={fl['r_var']} r_0={fl['r_0']}: {v} "
                  f"(variant-in {r[0]:.1e}, in-variant fixed L {r[1]:.1e}, L free {r[2]:.1e})")
