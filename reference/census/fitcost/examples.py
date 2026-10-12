# ABOUTME: Fit-cost calibration, step 4: per-class summary tables, printed examples (step lists with
# ABOUTME: their fit cost) and the loaded-heap vs fresh-process comparison of the same candidates.
import csv
import math
import sys
from collections import defaultdict

import numpy as np

D = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost/data"


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


text = {}
for name in ("selection.tsv", "selection_big.tsv"):
    for r in read(f"{D}/{name}"):
        text[(r["rxn"], r["gid"])] = r["text"]

TAG = sys.argv[1] if len(sys.argv) > 1 else "fresh2"
fresh = [r for r in read(f"{D}/fits_{TAG}.tsv") if r["error"] == ""]
for r in fresh:
    r["P"] = int(r["rank"]) < int(r["np"])
    if "pkg_compile" in r:              # identify's own fit, cold
        r["t_fit_rep"] = r["t_fit"]
        r["t_fit"] = r["t_fit_pkg"]
        r["fit_compile"] = r["pkg_compile"]
    r["tc"] = sum(float(r[k]) for k in ("t_pass1", "t_problem", "t_first_loss", "t_fit"))
    r["G"] = int(r["gens_total"]) / int(r["n_restarts"])

print("per reaction and class (fresh process): n | np median | ops median | G/20 median | "
      "t_fit median, mean, max s | t_candidate mean s | fixed share of t_candidate (mean)")
for rxn in ("bi_bi_rxn", "ter_bi_rxn"):
    for cls, sel in (("full rank (OK+REF+N)", lambda r: not r["P"]), ("P (rank-deficient)", lambda r: r["P"]),
                     ("P plain", lambda r: r["P"] and r["allo"] == "false"),
                     ("P allosteric", lambda r: r["P"] and r["allo"] == "true"),
                     ("full allosteric", lambda r: not r["P"] and r["allo"] == "true"),
                     ("full plain", lambda r: not r["P"] and r["allo"] == "false")):
        sub = [r for r in fresh if r["rxn"] == rxn and sel(r)]
        if not sub:
            continue
        tf = [float(r["t_fit"]) for r in sub]
        fixed = [(r["tc"] - float(r["t_fit"]) + max(float(r["fit_compile"]), 0)) / r["tc"] for r in sub]
        print(f"  {rxn:10s} {cls:22s} n={len(sub):3d} | np {np.median([int(r['np']) for r in sub]):4.1f} | "
              f"ops {np.median([int(r['pkg_ops']) for r in sub]):6.0f} | G/20 {np.median([r['G'] for r in sub]):6.0f} | "
              f"t_fit {np.median(tf):5.2f} {np.mean(tf):5.2f} {max(tf):5.2f} | t_cand {np.mean([r['tc'] for r in sub]):5.2f} | "
              f"fixed {np.mean(fixed):.2f}")

print("\nexamples (fresh process):")
def show(r):
    print(f"  {r['rxn']} #{r['gid']} {r['cls']} np {r['np']} rank {r['rank']} ops {r['pkg_ops']} Dterms {r['dterms']}: "
          f"t_fit {float(r['t_fit']):.2f} s, G/20 {r['G']:.0f}, stop S{r['n_success']}/St{r['n_stalled']}/T{r['n_maxtime']}, "
          f"first call {float(r['t_first_loss']):.3f} s, best loss {float(r['best_loss']):.4g}")
    print(f"      {text[(r['rxn'], r['gid'])]}")

for rxn in ("bi_bi_rxn", "ter_bi_rxn"):
    for label, sel in (("full rank, small law", lambda r: not r["P"] and int(r["pkg_ops"]) < 60),
                       ("full rank, large law", lambda r: not r["P"] and int(r["pkg_ops"]) > 1500),
                       ("rank-deficient, plain", lambda r: r["P"] and r["allo"] == "false"),
                       ("rank-deficient, allosteric", lambda r: r["P"] and r["allo"] == "true"),
                       ("flagged never-saturating", lambda r: r["cls"] == "N")):
        sub = sorted([r for r in fresh if r["rxn"] == rxn and sel(r)], key=lambda r: float(r["t_fit"]))
        if not sub:
            continue
        print(f" {rxn}, {label} ({len(sub)} fitted; cheapest, median, most expensive):")
        picks = [sub[0], sub[len(sub) // 2], sub[-1]] if len(sub) >= 3 else sub
        for r in picks:
            show(r)

# same candidates, loaded heap (run 1, after the selection phase in one process) vs fresh process
run1 = {(r["rxn"], r["gid"]): r for r in read(f"{D.replace('/data', '/data_run1')}/fits.tsv")}
pairs = [(run1[(r["rxn"], r["gid"])], r) for r in fresh if (r["rxn"], r["gid"]) in run1]
same_g = sum(a["gens_total"] == b["gens_total"] for a, b in pairs)
same_loss = sum(abs(float(a["best_loss"]) - float(b["best_loss"])) <= 1e-12 * max(1, float(a["best_loss"])) for a, b in pairs)
ratio = np.array([float(b["t_fit"]) / float(a["t_fit"]) for a, b in pairs])
gc_a = sum(float(a["fit_gctime"]) for a, b in pairs); t_a = sum(float(a["t_fit"]) for a, b in pairs)
gc_b = sum(float(b["fit_gctime"]) for a, b in pairs); t_b = sum(float(b["t_fit"]) for a, b in pairs)
p1a = np.array([float(a["t_pass1"]) for a, b in pairs]); p1b = np.array([float(b["t_pass1"]) for a, b in pairs])
print(f"\nsame {len(pairs)} candidates, run 1 (one process after the selection phase, maxrss 2.0 GB; replicate loop "
      f"specialized per mechanism) vs this fresh run (t_fit = fit_rate_equation, cold):")
print(f"  identical total generations {same_g}/{len(pairs)}, identical best loss {same_loss}/{len(pairs)}")
print(f"  t_fit fresh/run1: median {np.median(ratio):.2f}, range {ratio.min():.2f}..{ratio.max():.2f}; "
      f"summed t_fit {t_a:.1f} → {t_b:.1f} s; GC {gc_a:.1f} s ({100*gc_a/t_a:.0f}%) → {gc_b:.1f} s ({100*gc_b/t_b:.0f}%)")
print(f"  pass 1 s: run 1 median {np.median(p1a):.4f} (generated functions already compiled) → fresh median "
      f"{np.median(p1b):.4f}, max {p1b.max():.3f}")
