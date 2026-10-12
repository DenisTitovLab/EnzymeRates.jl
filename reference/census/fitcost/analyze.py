# ABOUTME: Fit-cost calibration, step 3: reads the fit, restart and CMA-ES overhead tables and fits a
# ABOUTME: cost model (time per fit vs parameter count, law size, rank deficiency) with its scatter.
# usage: analyze.py <data dir> <tag of fits_/restarts_/cma_overhead_ files, or ""> <population file>
import csv
import math
import sys
from collections import Counter, defaultdict

import numpy as np

D = sys.argv[1]
TAG = sys.argv[2] if len(sys.argv) > 2 else ""
POP = sys.argv[3] if len(sys.argv) > 3 else "population.tsv"
tagged = lambda base: f"{base}_{TAG}.tsv" if TAG else f"{base}.tsv"
NPTS = 48


def read(name):
    with open(f"{D}/{name}") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def lstsq(X, y):
    """Least squares with standard errors and R^2."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    dof = max(len(y) - X.shape[1], 1)
    s2 = res @ res / dof
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X)))
    r2 = 1 - (res @ res) / max(((y - y.mean()) ** 2).sum(), 1e-300)
    return beta, se, r2, res


def q(v, ps=(0, 25, 50, 75, 90, 100)):
    v = np.asarray(v, float)
    return "/".join(f"{np.percentile(v, p):.3g}" for p in ps) if len(v) else "-"


lam_of = lambda n: 4 + math.floor(3 * math.log(n))
allrows = read(tagged("fits"))
fits = [r for r in allrows if r["error"] == ""]
rs = read(tagged("restarts"))
ov = read(tagged("cma_overhead"))
pop = read(POP)
print(f"fits: {len(fits)} ok, {len(allrows) - len(fits)} errors")
for r in allrows:
    if r["error"] != "":
        print(f"  ERROR {r['rxn']} #{r['gid']} {r['cls']}: {r['error']}")

for f in fits:
    f["np_"] = int(f["np"])
    f["rk"] = int(f["rank"])
    f["P"] = f["rk"] < f["np_"]
    f["ops"] = int(f["pkg_ops"])
    f["tl"] = num(f["t_loss_med"])
    f["tf"] = num(f["t_fit"])
    f["comp"] = max(num(f["fit_compile"]), 0.0)
    f["gc"] = num(f["fit_gctime"])
    f["topt"] = f["tf"] - f["comp"]                  # the optimization itself
    f["g"] = int(f["gens_total"])
    f["nr"] = int(f["n_restarts"])
    f["lam"] = int(f["lambda"])
    f["fe"] = int(f["fevals_total"])
    f["pg"] = f["topt"] / f["g"]                      # seconds per generation
    f["p1"] = num(f["t_pass1"])
    f["prob"] = num(f["t_problem"]) if f.get("t_problem", "") not in ("", None) else 0.0
    f["jit"] = num(f["t_first_loss"])
    # identify's own fit of the candidate (cold) when measured, else the replicate's
    f["tfp"] = num(f["t_fit_pkg"]) if "pkg_compile" in f and f["t_fit_pkg"] != "" else f["tf"]
    f["pcomp"] = max(num(f["pkg_compile"]), 0.0) if "pkg_compile" in f else f["comp"]
    f["tc"] = f["p1"] + f["prob"] + f["jit"] + f["tfp"]  # per-candidate total
    f["ter"] = f["rxn"] == "ter_bi_rxn"
    f["allo_"] = f["allo"] == "true"

# ── CMA-ES overhead per generation (sphere) ────────────────────────────────────────────
byn = defaultdict(list)
for o in ov:
    byn[int(o["n"])].append(num(o["t"]) / num(o["gens"]))
ns = sorted(byn)
cn = np.array([np.median(byn[n]) for n in ns])
print("\nCMA-ES + Optimization.jl cost per generation with a trivial objective (sphere, 400 generations, median of 3):")
print("  n:      " + " ".join(f"{n:6d}" for n in ns))
print("  lambda: " + " ".join(f"{lam_of(n):6d}" for n in ns))
print("  us/gen: " + " ".join(f"{c * 1e6:6.1f}" for c in cn))
Bs, _, r2s, _ = lstsq(np.c_[np.ones(len(ns)), np.square(ns)], cn)
print(f"  c(n) ≈ {Bs[0]*1e6:.2f} + {Bs[1]*1e6:.4f}·n² us (R^2 {r2s:.3f})")

# ── loss evaluation time vs law size ────────────────────────────────────────────────────
ops = np.array([f["ops"] for f in fits], float)
tl = np.array([f["tl"] for f in fits])
B1, se1, r21, _ = lstsq(np.c_[np.ones(len(ops)), ops], tl)
print(f"\nloss! time per call ({NPTS} points) vs package-law size (ops = call nodes of num, den, dependents):")
print(f"  t_loss = {B1[0]*1e6:.3f} us + {B1[1]*1e9:.3f} ns·ops (R^2 {r21:.3f}); per data point "
      f"{B1[0]/NPTS*1e9:.1f} ns + {B1[1]/NPTS*1e9:.4f} ns·ops")
print(f"  ops {q(ops)}; t_loss us {q(tl*1e6)}")
dt = np.array([num(f["dterms"]) for f in fits])
okd = dt > 0
B1d, _, r21d, _ = lstsq(np.c_[np.ones(okd.sum()), dt[okd]], tl[okd])
Bo, _, r2o, _ = lstsq(np.c_[np.ones(okd.sum()), dt[okd]], ops[okd])
print(f"  vs slice D terms: t_loss = {B1d[0]*1e6:.3f} us + {B1d[1]*1e9:.2f} ns·Dterms (R^2 {r21d:.3f}); "
      f"ops ≈ {Bo[0]:.0f} + {Bo[1]:.2f}·Dterms (R^2 {r2o:.3f})")

# ── per-generation cost ─────────────────────────────────────────────────────────────────
pg = np.array([f["pg"] for f in fits])
gcg = np.array([f["gc"] / f["g"] for f in fits])
obj = np.array([f["lam"] * f["tl"] for f in fits])
n2 = np.array([f["np_"] ** 2 for f in fits], float)
print(f"\nper generation (optimization time without the in-fit compile / generations), us: {q(pg*1e6)}")
print(f"  objective part λ·t_loss us: {q(obj*1e6)}; its share of the generation: {q(obj/pg)}")
print(f"  GC part us: {q(gcg*1e6)}; GC share of the optimization time: {q(np.array([f['gc']/f['topt'] for f in fits]))}; "
      f"fits with no GC {sum(f['gc'] == 0 for f in fits)} of {len(fits)}")
rest = pg - obj - gcg
Bp, sep, r2p, _ = lstsq(np.c_[np.ones(len(fits)), n2], rest)
print(f"  remainder (CMA-ES and wrapper) us: {q(rest*1e6)}; ≈ {Bp[0]*1e6:.2f} + {Bp[1]*1e6:.4f}·n² (R^2 {r2p:.3f})")
Bpg, sepg, r2pg, respg = lstsq(np.c_[np.ones(len(fits)), n2, obj], pg)
print(f"  per-generation model pg ≈ {Bpg[0]*1e6:.1f} us + {Bpg[1]*1e6:.3f} us·n² + {Bpg[2]:.2f}·λ·t_loss "
      f"(R^2 {r2pg:.3f}; GC included)")
print(f"  allocation per generation KB: {q(np.array([num(f['fit_bytes'])/f['g']/1e3 for f in fits]))}")

# ── generations ─────────────────────────────────────────────────────────────────────────
print("\ngenerations per restart (G/20) by reaction, np and class (median [min..max], fits):")
for rxn in sorted(set(f["rxn"] for f in fits)):
    for cls in (False, True):
        by = defaultdict(list)
        for f in fits:
            if f["rxn"] == rxn and f["P"] == cls:
                by[f["np_"]].append(f["g"] / f["nr"])
        print(f"  {rxn} {'P   ' if cls else 'full'}: " + "; ".join(
            f"np {k}: {np.median(v):.0f} [{min(v):.0f}..{max(v):.0f}] n={len(v)}" for k, v in sorted(by.items())))
yg = np.array([math.log(f["g"] / f["nr"]) for f in fits])
X = np.c_[np.ones(len(fits)), [math.log(f["np_"]) for f in fits], [f["P"] for f in fits],
          [f["ter"] for f in fits], [f["allo_"] for f in fits], np.log(ops)]
Bg, seg, r2g, resg = lstsq(X, yg)
names = ["const", "log np", "P", "ter_bi", "allosteric", "log ops"]
print("  log(G/20) = " + " + ".join(f"{b:.2f}(±{s:.2f})·{n}" for b, s, n in zip(Bg, seg, names))
      + f";  R^2 {r2g:.3f}, residual sd {resg.std():.2f}")
Bg2, seg2, r2g2, resg2 = lstsq(X[:, :3], yg)
print(f"  np and P only: G/20 ≈ {math.exp(Bg2[0]):.2f}·np^{Bg2[1]:.2f}(±{seg2[1]:.2f})·{math.exp(Bg2[2]):.2f}^P"
      f" (P coef {Bg2[2]:.2f}±{seg2[2]:.2f}); R^2 {r2g2:.3f}, residual sd {resg2.std():.2f}")
Xd = np.c_[np.ones(len(fits)), [math.log(f["np_"]) for f in fits], [f["np_"] - f["rk"] for f in fits],
           [f["ter"] for f in fits]]
Bd, sed, r2d, _ = lstsq(Xd, yg)
print(f"  deficit (np − rank) per unit, with log np and ter_bi: {Bd[2]:.2f}±{sed[2]:.2f} (R^2 {r2d:.3f})")

# ── stop reasons, maxtime ───────────────────────────────────────────────────────────────
print("\nrestart stop reasons by class:")
key = {(f["rxn"], f["gid"]): ("P" if f["P"] else "full") for f in fits}
cnt, tmax, gmax = defaultdict(Counter), defaultdict(float), defaultdict(int)
for r in rs:
    k = key.get((r["rxn"], r["gid"]))
    if k is None:
        continue
    cnt[k][r["reason"]] += 1
    tmax[k] = max(tmax[k], num(r["t"]))
    gmax[k] = max(gmax[k], int(r["gens"]))
for k in sorted(cnt):
    print(f"  {k}: {dict(cnt[k])}; longest restart {tmax[k]:.2f} s ({gmax[k]} generations); "
          f"maxtime is 600 s per restart")
print(f"  fits truncated by the per-fit cap: {sum(f['truncated'] == 'true' for f in fits)}; "
      f"restarts ending at maxtime: {sum(int(f['n_maxtime']) for f in fits)}; at maxiters: "
      f"{sum(int(f['n_maxiters']) for f in fits)}")

# ── the per-fit model ─────────────────────────────────────────────────────────────────
tf = np.array([f["tfp"] for f in fits])
tc = np.array([f["tc"] for f in fits])
print("\nMODEL B (direct): log t = a + b·log np + c·log ops + d·P   [t_fit = fit_rate_equation, 20 restarts, cold]")
Xb = np.c_[np.ones(len(fits)), [math.log(f["np_"]) for f in fits], np.log(ops), [f["P"] for f in fits]]
BB = resB = None
for label, y in (("t_fit", tf), ("t_candidate (pass 1 + problem + first call + fit)", tc)):
    Bb, seb, r2b, resb = lstsq(Xb, np.log(y))
    print(f"  {label}: a {Bb[0]:.2f}±{seb[0]:.2f}, b {Bb[1]:.2f}±{seb[1]:.2f}, c {Bb[2]:.2f}±{seb[2]:.2f}, "
          f"d {Bb[3]:.2f}±{seb[3]:.2f}; R^2 {r2b:.3f}; residual sd {resb.std():.2f} (×{math.exp(resb.std()):.2f})")
    print(f"     t ≈ {math.exp(Bb[0]):.3g} s · np^{Bb[1]:.2f} · ops^{Bb[2]:.2f} · {math.exp(Bb[3]):.2f}^P")
    if BB is None:
        BB, resB = Bb, resb
dts = np.array([max(num(f["dterms"]), 1.0) for f in fits])
Xc = np.c_[np.ones(len(fits)), [math.log(f["np_"]) for f in fits], np.log(dts), [f["P"] for f in fits]]
Bc, sec, r2c, resc = lstsq(Xc, np.log(tc))
print(f"  t_candidate with slice D terms for size: t ≈ {math.exp(Bc[0]):.3g} s · np^{Bc[1]:.2f}(±{sec[1]:.2f}) · "
      f"Dterms^{Bc[2]:.2f}(±{sec[2]:.2f}) · {math.exp(Bc[3]):.2f}^P (P coef {Bc[3]:.2f}±{sec[3]:.2f}); "
      f"R^2 {r2c:.3f}, residual sd {resc.std():.2f}")
Xe = np.c_[Xb, [f["allo_"] for f in fits], [f["ter"] for f in fits]]
Be, see, r2e, _ = lstsq(Xe, np.log(tf))
print(f"  t_fit with allosteric and ter_bi indicators added: allo {Be[4]:.2f}±{see[4]:.2f}, ter_bi "
      f"{Be[5]:.2f}±{see[5]:.2f} (log units); P {Be[3]:.2f}±{see[3]:.2f}; R^2 {r2e:.3f}")
print(f"  t_fit s {q(tf)}; t_candidate s {q(tc)}")
for c in sorted(set(f["cls"] for f in fits)):
    sub = [f for f in fits if f["cls"] == c]
    print(f"  class {c}: n={len(sub)} t_fit median {np.median([f['tfp'] for f in sub]):.2f} mean "
          f"{np.mean([f['tfp'] for f in sub]):.2f} s; t_candidate mean {np.mean([f['tc'] for f in sub]):.2f} s; "
          f"np median {np.median([f['np_'] for f in sub]):.0f}")
print("  largest t_fit residuals of model B (actual/fitted):")
for i in np.argsort(-np.abs(resB))[:6]:
    f = fits[i]
    print(f"    {f['rxn']} #{f['gid']} {f['cls']} np {f['np_']} rk {f['rk']} ops {f['ops']}: "
          f"t_fit {f['tfp']:.2f} s, G/20 {f['g']/f['nr']:.0f}, ×{math.exp(resB[i]):.2f}")

# ── matched comparison P vs full rank ───────────────────────────────────────────────────
print("\nrank-deficient vs full-rank at equal reaction and np (medians):")
for rxn in sorted(set(f["rxn"] for f in fits)):
    for n in sorted(set(f["np_"] for f in fits if f["rxn"] == rxn)):
        a = [f for f in fits if f["rxn"] == rxn and f["np_"] == n and f["P"]]
        b = [f for f in fits if f["rxn"] == rxn and f["np_"] == n and not f["P"]]
        if a and b:
            print(f"  {rxn} np {n}: P n={len(a)} G/20 {np.median([f['g']/f['nr'] for f in a]):.0f} "
                  f"t_fit {np.median([f['tfp'] for f in a]):.2f} s | full n={len(b)} G/20 "
                  f"{np.median([f['g']/f['nr'] for f in b]):.0f} t_fit {np.median([f['tfp'] for f in b]):.2f} s")

# ── fixed per-candidate costs ───────────────────────────────────────────────────────────
WARM = {("bi_bi_rxn", "1"), ("ter_bi_rxn", "5")}   # fitted during the warm-up: no cold compile
cold = [f for f in fits if (f["rxn"], f["gid"]) not in WARM]
p1 = np.array([f["p1"] for f in cold])
jit = np.array([f["jit"] for f in cold])
comp = np.array([f["pcomp"] for f in cold])
prob = np.array([f["prob"] for f in cold])
ops = np.array([f["ops"] for f in cold], float)
tc = np.array([f["tc"] for f in cold])
print(f"\nfixed costs per candidate, s ({len(cold)} candidates; the warm-up mechanisms excluded):")
print(f"  pass 1 (_eq_complexity, compile_mechanism, fitted_params, rate_equation_string, type string, key): {q(p1)}")
Bq, _, r2q, _ = lstsq(np.c_[np.ones(len(ops)), ops], p1)
print(f"     ≈ {Bq[0]*1e3:.2f} ms + {Bq[1]*1e6:.2f} us·ops (R^2 {r2q:.3f})")
print(f"  FittingProblem: {q(prob)}")
print(f"  first loss! call (rate_equation derivation + JIT of rate_equation and loss!): {q(jit)}")
Bj, _, r2j, _ = lstsq(np.c_[np.ones(len(ops)), ops], jit)
print(f"     ≈ {Bj[0]*1e3:.1f} ms + {Bj[1]*1e6:.1f} us·ops (R^2 {r2j:.3f})")
print(f"  compile inside fit_rate_equation (@timed compile_time; FunctionWrapper thunk, _kcat_forward): {q(comp)}")
fixed = p1 + prob + jit + comp
print(f"  fixed share of the per-candidate total: {q(fixed / tc)}; summed over the sample "
      f"{fixed.sum():.1f} of {tc.sum():.1f} s ({100*fixed.sum()/tc.sum():.0f}%)")

# ── validation against fit_rate_equation ───────────────────────────────────────────────
print("\nreplicated restart loop vs EnzymeRates.fit_rate_equation on the same problem and seed:")
both = [f for f in fits if f["t_fit_pkg"] != ""]
same = sum(num(f["loss_pkg"]) == num(f["best_loss"]) for f in both)
print(f"  {len(both)} candidates with both; identical best loss in {same}")
for f in both[:4]:
    print(f"  {f['rxn']} #{f['gid']} {f['cls']}: replicate {f['tf']:.2f} s, loss {num(f['best_loss']):.6g} | "
          f"package {num(f['t_fit_pkg']):.2f} s, loss {num(f['loss_pkg']):.6g}, retcode {f['ret_pkg']}")

# ── the model applied to the pools' eligible candidates ─────────────────────────────────
Bt, _, _, _ = lstsq(Xb, np.log(np.array([f["tc"] for f in fits])))
print("\nmodel B (t_candidate) applied to every eligible pool record (predicted seconds; P shares):")
for rxn in sorted(set(p["rxn"] for p in pop)):
    for rnd in sorted(set(p["round"] for p in pop if p["rxn"] == rxn)):
        sub = [p for p in pop if p["rxn"] == rxn and p["round"] == rnd]
        t = np.array([math.exp(Bt[0] + Bt[1] * math.log(int(p["np"])) + Bt[2] * math.log(int(p["pkg_ops"]))
                               + Bt[3] * (p["P"] == "true")) for p in sub])
        isP = np.array([p["P"] == "true" for p in sub])
        npsum = np.array([int(p["np"]) for p in sub])
        print(f"  {rxn} round {rnd}: {len(sub)} records, predicted {t.sum():.0f} s (mean {t.mean():.2f} s); "
              f"P {isP.sum()} = {100*isP.mean():.1f}% of count, {100*t[isP].sum()/t.sum():.1f}% of predicted time, "
              f"{100*npsum[isP].sum()/npsum.sum():.1f}% of summed params")
