# ABOUTME: Samples R5.5 candidates for the slice cost check (keytime.jl "code" rows): allosteric laws
# ABOUTME: newly over 337, allosteric ones still under, plain ones near and above 337.
import random

D = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/keycost/data/"
files = ["r55_bibi_largest.tsv", "r55_bibiR1_largest.tsv", "r55_pingpong1_largest.tsv",
         "r55_bibi_uniform.tsv"]
seen, rows = set(), []
for fn in files:
    for line in open(D + fn):
        f = line.rstrip("\n").split("\t")
        if f[1] in seen:
            continue
        seen.add(f[1])
        rows.append((fn.split("_")[1], f))
strata = {
    "allo_new": lambda f: f[2] == "1" and int(f[4]) <= 337 < int(f[4]) + int(f[5]),
    "allo_under": lambda f: f[2] == "1" and int(f[4]) + int(f[5]) <= 337 and int(f[4]) > 120,
    "plain_near": lambda f: f[2] == "0" and 200 < int(f[4]) <= 337,
    "plain_over": lambda f: f[2] == "0" and 337 < int(f[4]) <= 720,
}
quota = {"allo_new": 160, "allo_under": 60, "plain_near": 60, "plain_over": 60}
rng = random.Random(20261014)
with open(D + "r55_cost.tsv", "w") as out:
    for name, pred in strata.items():
        pool = [(r, f) for r, f in rows if pred(f)]
        pick = pool if len(pool) <= quota[name] else rng.sample(pool, quota[name])
        print(name, len(pool), len(pick))
        for r, f in pick:
            out.write("\t".join([name, r + ":r" + f[0], f[1], f[7], f[8], f[4], f[5], "", "", "",
                                 f[9]]) + "\n")
