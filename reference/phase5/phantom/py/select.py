# ABOUTME: Picks the numeric-containment sample: flagged variants with and without a certificate,
# ABOUTME: rank-equal unflagged ones, rank-gaining numeric candidates, enumeration-parent comparisons.
import random, sys
sys.path.insert(0, ".")
from load import *
rng = random.Random(20261014)
rows = []
allv = []
for rxn in ("bibi", "pingpong", "bibi_R", "uniuni_R"):
    F = facts(rxn)
    allv += [(rxn, r) for r in F.values() if r["kind"] == "allo1"]
def pick(cond, k, cmp, tag):
    xs = [(rxn, r) for rxn, r in allv if cond(rxn, r)]
    xs = xs if len(xs) <= k else rng.sample(xs, k)
    for rxn, r in xs:
        rows.append((rxn, r["gid"], cmp, tag, len([1 for q in allv if cond(*q)])))
pick(lambda x, r: r["flag"] and r["cert"] == 0, 100, "LP", "S1_flag_nocert")
pick(lambda x, r: r["flag"] and r["cert"] == 1, 40, "LP", "S2_flag_cert")
pick(lambda x, r: (not r["flag"]) and r["r_var"] == r["r_0"] and r["cert"] == 0, 60, "LP", "S3_rankeq_noflag_nocert")
pick(lambda x, r: (not r["flag"]) and r["r_var"] == r["r_0"] and r["cert"] == 1, 20, "LP", "S3b_rankeq_noflag_cert")
pick(lambda x, r: r["num"] and r["r_var"] > r["r_0"], 30, "LP", "S4_num_rankgain")
pick(lambda x, r: r["num"] and r["move"] != "to_allosteric" and r["round"] > 0, 40, "parent", "S5_num_alloparent")
pick(lambda x, r: r["num"] and r["move"] == "to_allosteric", 20, "parent", "S6_num_plainparent")
pick(lambda x, r: x == "uniuni_R" and r["gid"] == 7 or x == "bibi_R" and r["gid"] == 483, 5, "LP", "S7_negative_control")
with open(D + "/sel.tsv", "w") as f:
    for rxn, g, cmp, tag, n in rows:
        f.write(f"{rxn}\t{g}\t{cmp}\n")
with open(D + "/sel_tags.tsv", "w") as f:
    for rxn, g, cmp, tag, n in rows:
        f.write(f"{rxn}\t{g}\t{cmp}\t{tag}\t{n}\n")
import collections
print(collections.Counter((t, n) for _, _, _, t, n in rows))
