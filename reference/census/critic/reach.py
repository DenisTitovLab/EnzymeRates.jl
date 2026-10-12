# ABOUTME: Upper bound on the round-2 laws lost if round-1 tag-sibling duplicates (D) or L-phantom
# ABOUTME: candidates were never emitted, from the census's first-parent attribution of children.
import csv
from collections import defaultdict
C = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census"
rows = list(csv.DictReader(open(f"{C}/critic/small.tsv"), delimiter="\t"))
for run in ("bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200", "uni_bi_rxn"):
    R = [r for r in rows if r["run"] == run]
    by = {int(r["gid"]): r for r in R}
    r1 = [r for r in R if r["round"] == "1"]
    SD = {int(r["gid"]) for r in r1 if r["fitted"] == "true" and r["D"] == "true" and r["Z"] != "true"}
    SL = {int(r["gid"]) for r in r1 if r["fitted"] == "true" and r["Psub"] == "L_phantom" and r["D"] != "true"}
    SW = {int(r["gid"]) for r in r1 if r["fitted"] == "true" and r["W"] == "true"}
    r2 = [r for r in R if r["round"] == "2"]
    parents2 = {int(r["parent"]) for r in r2}
    early = {r["key"] for r in R if r["round"] in ("0", "1") and r["key"]}
    ok = lambda r: (r["fitted"] == "true" and r["W"] != "true" and r["D"] != "true" and r["P"] != "true"
                    and r["Z"] != "true" and r["N"] != "true")
    def lost(drop):
        avail = set(early) | {r["key"] for r in r2 if int(r["parent"]) not in drop}
        orph = [r for r in r2 if int(r["parent"]) in drop]
        lk = {r["key"] for r in orph} - avail
        lrec = [r for r in orph if r["key"] in lk]
        lok = {r["key"] for r in lrec if ok(r)}
        lfit = {r["key"] for r in lrec if r["fitted"] == "true" and r["W"] != "true"}
        return len(orph), len(lk), len(lfit), len(lok)
    keys2 = {r["key"] for r in r2} - early
    okk2 = {r["key"] for r in r2 if ok(r)}
    print(f"{run}: round-1 parents expanded {len(parents2 & {int(r['gid']) for r in r1})}; new round-2 laws (keys) {len(keys2)}, of them OK-class {len(okk2)}")
    for lbl, drop in (("D tag siblings", SD), ("L-phantom not D", SL), ("D + L-phantom", SD | SL), ("W", SW)):
        dp = drop & parents2
        o, lk, lf, lo = lost(drop)
        print(f"  drop {lbl}: {len(drop)} round-1 candidates ({len(dp)} expanded); their round-2 children {o}; "
              f"laws found only through them {lk} (non-W {lf}, OK-class {lo}) = {100*lk/max(len(keys2),1):.1f}% of new laws")
