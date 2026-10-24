# ABOUTME: For each spec-rule phantom whose L->0 key a relative holds, the round of the earliest such
# ABOUTME: relative relative to the phantom's round, and how the relative relates (ancestor or child).
import collections, sys
sys.path.insert(0, ".")
from load import *
c = collections.Counter()
for rxn in ("bibi", "pingpong"):
    F, R, E = facts(rxn), recs(rxn), edges(rxn)
    for r in F.values():
        if not (r["kind"] == "allo1" and r["flag"]): continue
        anc, g = [], r["parent"]
        while g != 0:
            anc.append(g); g = int(R[g]["parent"])
        hits = []
        for a in anc:
            if F[a]["kind"] == "plain" and F[a]["key"] == r["key0"]:
                hits.append((F[a]["round"], "ancestor"))
            for ch in E.get(a, []):
                x = F.get(ch)
                if x and x["kind"] == "plain" and x["status"] == "ok" and x["key"] == r["key0"]:
                    hits.append((x["round"], "ancestor's child"))
        if not hits:
            c[(rxn, r["move"], "none")] += 1; continue
        rd, how = min(hits)
        c[(rxn, r["move"], how, "earlier round" if rd < r["round"] else "same round")] += 1
for k in sorted(c): print(k, c[k])
