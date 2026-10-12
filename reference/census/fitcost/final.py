# ABOUTME: Fit-cost calibration, step 5: the additive structural cost model (fixed costs plus
# ABOUTME: generations × per-generation cost) on the corrected fresh run, its binned scatter and class sums.
import csv
import math
from collections import defaultdict

import numpy as np

D = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost/data"
WARM = {("bi_bi_rxn", "1"), ("ter_bi_rxn", "5")}


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def lstsq(X, y):
    X, y = np.asarray(X, float), np.asarray(y, float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    se = np.sqrt(np.diag((r @ r / max(len(y) - X.shape[1], 1)) * np.linalg.pinv(X.T @ X)))
    return b, se, 1 - (r @ r) / ((y - y.mean()) ** 2).sum(), r


lam = lambda n: 4 + math.floor(3 * math.log(n))
F = [r for r in read(f"{D}/fits_fresh2.tsv") if r["error"] == "" and (r["rxn"], r["gid"]) not in WARM]
n = np.array([int(r["np"]) for r in F], float)
ops = np.array([int(r["pkg_ops"]) for r in F], float)
P = np.array([int(r["rank"]) < int(r["np"]) for r in F], float)
ter = np.array([r["rxn"] == "ter_bi_rxn" for r in F], float)
G = np.array([int(r["gens_total"]) for r in F], float)                  # all 20 restarts
tl = np.array([float(r["t_loss_med"]) for r in F])
t_rep = np.array([float(r["t_fit"]) for r in F])                       # warm replicate
t_pkg = np.array([float(r["t_fit_pkg"]) for r in F])                   # identify's fit, cold
fixed = np.array([float(r["t_pass1"]) + float(r["t_problem"]) + float(r["t_first_loss"]) for r in F])
tc = fixed + t_pkg
L = np.array([lam(int(k)) for k in n], float)

# 1. per-generation cost from the warm replicate: pg = c0 + c2 n^2 + k·λ·t_loss
pg = t_rep / G
(c0, c2, k), _, r2pg, _ = lstsq(np.c_[np.ones(len(F)), n ** 2, L * tl], pg)
# 2. loss time from law size
(a0, a1), _, r2l, _ = lstsq(np.c_[np.ones(len(F)), ops], tl)
# 3. generations from np, P, reaction
(g0, g1, g2, g3), gse, r2g, rg = lstsq(np.c_[np.ones(len(F)), np.log(n), P, ter], np.log(G / 20))
sg = rg.std()
# 4. cold extra inside fit_rate_equation and the fixed costs, from law size
(e0, e1), _, r2e, _ = lstsq(np.c_[np.ones(len(F)), ops], t_pkg - t_rep)
(f0, f1), _, r2f, _ = lstsq(np.c_[np.ones(len(F)), ops], fixed)

print(f"structural model (fresh process, {len(F)} candidates, 48 data points, 20 restarts):")
print(f"  per generation  = {c0*1e6:.2f} us + {c2*1e6:.4f} us·np² + {k:.2f}·λ(np)·t_loss   (R^2 {r2pg:.3f})")
print(f"  t_loss          = {a0*1e6:.3f} us + {a1*1e9:.3f} ns·ops                    (R^2 {r2l:.3f})")
print(f"  generations/20  = exp({g0:.2f}) · np^{g1:.2f}(±{gse[1]:.2f}) · exp({g2:.2f}(±{gse[2]:.2f})·P) · exp({g3:.2f}(±{gse[3]:.2f})·ter_bi)"
      f"   (R^2 {r2g:.3f}, residual sd {sg:.2f})")
print(f"  cold extra in fit_rate_equation = {e0*1e3:.1f} ms + {e1*1e6:.2f} us·ops        (R^2 {r2e:.3f})")
print(f"  fixed (pass 1 + FittingProblem + first loss! call) = {f0*1e3:.1f} ms + {f1*1e6:.2f} us·ops (R^2 {r2f:.3f})")

def predict(n_, ops_, P_, ter_, mean=False):
    g = 20 * math.exp(g0 + g1 * math.log(n_) + g2 * P_ + g3 * ter_ + (sg ** 2 / 2 if mean else 0))
    per = c0 + c2 * n_ ** 2 + k * lam(n_) * (a0 + a1 * ops_)
    return f0 + f1 * ops_ + e0 + e1 * ops_ + g * per

pred = np.array([predict(n[i], ops[i], P[i], ter[i]) for i in range(len(F))])
lr = np.log(tc / pred)
print(f"\n  t_candidate actual/predicted (median-G prediction): median {math.exp(np.median(lr)):.2f}, "
      f"IQR {math.exp(np.percentile(lr, 25)):.2f}..{math.exp(np.percentile(lr, 75)):.2f}, "
      f"5-95% {math.exp(np.percentile(lr, 5)):.2f}..{math.exp(np.percentile(lr, 95)):.2f}, "
      f"range {math.exp(lr.min()):.2f}..{math.exp(lr.max()):.2f}; log-space R^2 "
      f"{1 - (lr @ lr) / ((np.log(tc) - np.log(tc).mean()) ** 2).sum():.3f}")
predm = np.array([predict(n[i], ops[i], P[i], ter[i], mean=True) for i in range(len(F))])
print(f"  summed: actual {tc.sum():.1f} s, predicted with mean G {predm.sum():.1f} s, with median G {pred.sum():.1f} s")
print("\n  binned scatter (by predicted t_candidate, deciles): predicted range → actual median [IQR], n")
edges = np.percentile(pred, np.linspace(0, 100, 11))
for i in range(10):
    m = (pred >= edges[i]) & (pred <= edges[i + 1] if i == 9 else pred < edges[i + 1])
    print(f"    {edges[i]:.3f}–{edges[i+1]:.3f} s → {np.median(tc[m]):.3f} [{np.percentile(tc[m],25):.3f}–"
          f"{np.percentile(tc[m],75):.3f}] s, n={m.sum()}")

print("\n  where the time goes (sample sums, s): fixed", f"{fixed.sum():.1f}", "| cold extra in fit",
      f"{(t_pkg - t_rep).sum():.1f}", "| optimization (warm replicate)", f"{t_rep.sum():.1f}",
      f"| of which λ·t_loss·G {np.sum(L * tl * G):.1f}, GC {sum(float(r['fit_gctime']) for r in F):.1f}")

print("\n  sample sums by reaction and class (t_candidate):")
for rxn in ("bi_bi_rxn", "ter_bi_rxn"):
    m = np.array([r["rxn"] == rxn for r in F])
    for lab, mm in (("full rank", m & (P == 0)), ("P", m & (P == 1))):
        print(f"    {rxn} {lab:9s}: n={mm.sum():3d}, sum {tc[mm].sum():6.1f} s, mean {tc[mm].mean():.3f} s, "
              f"median {np.median(tc[mm]):.3f} s; mean np {n[mm].mean():.2f}, mean ops {ops[mm].mean():.0f}")

print("\n  predicted mean t_candidate per np (P=0 / P=1), median law size of the sample at that np:")
for rxn, tv in (("bi_bi_rxn", 0), ("ter_bi_rxn", 1)):
    row = []
    for k_ in range(5, 12):
        m = (n == k_) & (ter == tv)
        o = np.median(ops[m]) if m.any() else 150
        row.append(f"np {k_}: {predict(k_, o, 0, tv, True):.2f}/{predict(k_, o, 1, tv, True):.2f}")
    print(f"    {rxn}: " + "; ".join(row))
