# ABOUTME: Ping-pong check of the tag-sibling reachability bound, restricted to dropped siblings
# ABOUTME: whose kept representative was also expanded (the 1,200-parent sample covers few pairs).
import csv
C = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census"
rows = [r for r in csv.DictReader(open(f"{C}/critic/small.tsv"), delimiter="\t") if r["run"] == "bi_bi_pp_rxn.r2_1200"]
r1 = [r for r in rows if r["round"] == "1"]; r2 = [r for r in rows if r["round"] == "2"]
parents2 = {int(r["parent"]) for r in r2}
early = {r["key"] for r in rows if r["round"] in ("0", "1")}
SD = {int(r["gid"]): int(r["Dof"]) for r in r1 if r["fitted"] == "true" and r["D"] == "true" and r["Z"] != "true"}
both = {y for y, x in SD.items() if y in parents2 and x in parents2}
only_y = {y for y, x in SD.items() if y in parents2 and x not in parents2}
avail = early | {r["key"] for r in r2 if int(r["parent"]) not in SD}
def lost(sel):
    orph = [r for r in r2 if int(r["parent"]) in sel]
    lk = {r["key"] for r in orph} - avail
    okk = {r["key"] for r in orph if r["key"] in lk and r["fitted"] == "true" and
           all(r[c] != "true" for c in ("W", "D", "P", "Z", "N"))}
    return len(orph), len(lk), len(okk)
print("pp 1200: dropped siblings expanded with their representative also expanded:", len(both), "children/lost keys/OK lost", lost(both))
print("pp 1200: dropped siblings expanded without their representative:", len(only_y), "children/lost keys/OK lost", lost(only_y))
from collections import Counter
orph = [r for r in r2 if int(r["parent"]) in both]
lk = {r["key"] for r in orph} - avail
lr = [r for r in orph if r["key"] in lk and r["fitted"] == "true"]
print("lost fitted records by move:", Counter(r["move"] for r in lr))
print("lost fitted records by class:", Counter(",".join(c for c in ("W", "P", "Z", "N") if r[c] == "true") or "OK" for r in lr))
