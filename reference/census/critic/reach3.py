# ABOUTME: Classes and moves of the bi-bi round-2 laws reachable only through L-phantom parents
# ABOUTME: (complete expansion of all round-1 parents, so no sampling confound).
import csv
from collections import Counter
C = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census"
rows = [r for r in csv.DictReader(open(f"{C}/critic/small.tsv"), delimiter="\t") if r["run"] == "bi_bi_rxn.r2_100000"]
r1 = [r for r in rows if r["round"] == "1"]; r2 = [r for r in rows if r["round"] == "2"]
early = {r["key"] for r in rows if r["round"] in ("0", "1")}
SL = {int(r["gid"]) for r in r1 if r["fitted"] == "true" and r["Psub"] == "L_phantom" and r["D"] != "true"}
avail = early | {r["key"] for r in r2 if int(r["parent"]) not in SL}
orph = [r for r in r2 if int(r["parent"]) in SL]
lk = {r["key"] for r in orph} - avail
lr = [r for r in orph if r["key"] in lk and r["fitted"] == "true"]
print("lost fitted:", len(lr), "by move:", Counter(r["move"] for r in lr))
print("by class:", Counter((r["Psub"] if r["P"] == "true" else "OK" if all(r[c] != "true" for c in ("W","D","Z","N")) else "other") for r in lr))
okr = [r for r in lr if r["P"] != "true" and all(r[c] != "true" for c in ("W","D","Z","N"))]
for r in okr: print("  OK lost:", r["gid"], r["move"], "parent", r["parent"], "np", r["np"])
