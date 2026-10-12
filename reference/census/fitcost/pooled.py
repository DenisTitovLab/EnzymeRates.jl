# ABOUTME: Fit-cost calibration, step 7: the structural cost model fitted on both datasets together
# ABOUTME: (fresh2 and alt), its scatter per dataset, and the model's P time shares on the pools.
import csv
import math

import numpy as np

D = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost/data"
WARM = {("bi_bi_rxn", "1"), ("ter_bi_rxn", "5")}
lam = lambda n: 4 + math.floor(3 * math.log(n))


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def lstsq(X, y):
    X, y = np.asarray(X, float), np.asarray(y, float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    se = np.sqrt(np.diag((r @ r / max(len(y) - X.shape[1], 1)) * np.linalg.pinv(X.T @ X)))
    return b, se, 1 - (r @ r) / ((y - y.mean()) ** 2).sum(), r


rows = []
for tag in ("fresh2", "alt"):
    for r in read(f"{D}/fits_{tag}.tsv"):
        if r["error"] == "" and (r["rxn"], r["gid"]) not in WARM:
            r["tag"] = tag
            rows.append(r)
n = np.array([int(r["np"]) for r in rows], float)
ops = np.array([int(r["pkg_ops"]) for r in rows], float)
P = np.array([int(r["rank"]) < int(r["np"]) for r in rows], float)
ter = np.array([r["rxn"] == "ter_bi_rxn" for r in rows], float)
alt = np.array([r["tag"] == "alt" for r in rows], float)
G = np.array([int(r["gens_total"]) for r in rows], float)
tl = np.array([float(r["t_loss_med"]) for r in rows])
trep = np.array([float(r["t_fit"]) for r in rows])
tpkg = np.array([float(r["t_fit_pkg"]) for r in rows])
fixed = np.array([float(r["t_pass1"]) + float(r["t_problem"]) + float(r["t_first_loss"]) for r in rows])
tc = fixed + tpkg
L = np.array([lam(int(k)) for k in n], float)

(c0, c2, k), _, r2pg, _ = lstsq(np.c_[np.ones(len(rows)), n ** 2, L * tl], trep / G)
(a0, a1), _, r2l, _ = lstsq(np.c_[np.ones(len(rows)), ops], tl)
Xg = np.c_[np.ones(len(rows)), np.log(n), P, ter, alt, ter * alt]
bg, sg_, r2g, rg = lstsq(Xg, np.log(G / 20))
sg = rg.std()
(x0, x1), _, _, _ = lstsq(np.c_[np.ones(len(rows)), ops], fixed + (tpkg - trep))
print(f"pooled structural model ({len(rows)} fits = 331 candidates × 2 datasets; 48 points; 20 restarts):")
print(f"  t_candidate = F + 20·G·g")
print(f"  F (pass 1 + FittingProblem + first loss! call + compile inside fit_rate_equation) = {x0*1e3:.0f} ms + {x1*1e6:.0f} us·ops")
print(f"  g (one generation) = {c0*1e6:.1f} us + {c2*1e6:.3f} us·np² + {k:.2f}·λ(np)·t_loss,  λ = 4 + ⌊3 ln np⌋   (R^2 {r2pg:.3f})")
print(f"  t_loss (48 points) = {a0*1e6:.2f} us + {a1*1e9:.2f} ns·ops   (R^2 {r2l:.3f}); per point {a0/48*1e9:.1f} ns + {a1/48*1e9:.3f} ns·ops")
names = ["const", "log np", "P", "ter_bi", "alt data", "ter_bi×alt"]
print("  log G = " + " + ".join(f"{b:.2f}(±{s:.2f})·{nm}" for b, s, nm in zip(bg, sg_, names))
      + f"   (G = generations per restart; R^2 {r2g:.3f}, residual sd {sg:.2f})")

def predict(i, mean=True):
    lg = bg @ Xg[i] + (sg ** 2 / 2 if mean else 0)
    return x0 + x1 * ops[i] + 20 * math.exp(lg) * (c0 + c2 * n[i] ** 2 + k * L[i] * (a0 + a1 * ops[i]))

pred = np.array([predict(i, mean=False) for i in range(len(rows))])
for tag, m in (("main data", alt == 0), ("alt data", alt == 1)):
    lr = np.log(tc[m] / pred[m])
    pm = np.array([predict(i) for i in np.where(m)[0]])
    print(f"  {tag}: actual/predicted median {math.exp(np.median(lr)):.2f}, IQR {math.exp(np.percentile(lr,25)):.2f}.."
          f"{math.exp(np.percentile(lr,75)):.2f}, 5–95% {math.exp(np.percentile(lr,5)):.2f}..{math.exp(np.percentile(lr,95)):.2f}; "
          f"summed actual {tc[m].sum():.1f} s vs predicted (mean G) {pm.sum():.1f} s")

# P shares in the pools (main-data coefficients: alt = 0)
pop = read(f"{D}/population_selection_big.tsv")
print("\n  P share of the pools' predicted time (pooled model, main-data coefficients; mean G):")
print("  pool round | records | P count share | summed-np share | D-terms share | ops share | model share | model total (serial s)")
for rxn in ("bi_bi_rxn", "ter_bi_rxn"):
    for rnd in ("0", "1", "2"):
        sub = [p for p in pop if p["rxn"] == rxn and p["round"] == rnd]
        isP = np.array([p["P"] == "true" for p in sub])
        nn = np.array([int(p["np"]) for p in sub], float)
        oo = np.array([int(p["pkg_ops"]) for p in sub], float)
        dd = np.array([int(p["dterms"]) for p in sub], float)
        tt = np.array([x0 + x1 * oo[i] + 20 * math.exp(bg[0] + bg[1] * math.log(nn[i]) + bg[2] * isP[i]
                       + bg[3] * (rxn == "ter_bi_rxn") + sg ** 2 / 2)
                       * (c0 + c2 * nn[i] ** 2 + k * lam(int(nn[i])) * (a0 + a1 * oo[i])) for i in range(len(sub))])
        sh = lambda w: 100 * w[isP].sum() / w.sum()
        print(f"  {rxn} r{rnd} | {len(sub)} | {100*isP.mean():.1f}% | {sh(nn):.1f}% | {sh(dd):.1f}% | {sh(oo):.1f}% | "
              f"{sh(tt):.1f}% | {tt.sum():.0f}")

# binned scatter of the pooled model over both datasets
print("\n  binned scatter, both datasets (deciles of predicted t_candidate, median-G prediction):")
edges = np.percentile(pred, np.linspace(0, 100, 11))
for i in range(10):
    m = (pred >= edges[i]) & ((pred <= edges[i + 1]) if i == 9 else (pred < edges[i + 1]))
    print(f"    predicted {edges[i]:.3f}–{edges[i+1]:.3f} s → actual median {np.median(tc[m]):.3f} s "
          f"[IQR {np.percentile(tc[m],25):.3f}–{np.percentile(tc[m],75):.3f}], n={m.sum()}")
lr = np.log(tc / pred)
print(f"  all {len(rows)}: actual/predicted median {math.exp(np.median(lr)):.2f}, IQR {math.exp(np.percentile(lr,25)):.2f}.."
      f"{math.exp(np.percentile(lr,75)):.2f}, 5–95% {math.exp(np.percentile(lr,5)):.2f}..{math.exp(np.percentile(lr,95)):.2f}, "
      f"min {math.exp(lr.min()):.2f}, max {math.exp(lr.max()):.2f}; log-space R^2 "
      f"{1 - ((lr - 0) @ (lr - 0)) / ((np.log(tc) - np.log(tc).mean()) ** 2).sum():.3f}")
