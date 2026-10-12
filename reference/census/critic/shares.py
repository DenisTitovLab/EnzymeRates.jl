# ABOUTME: Weights the census classes by the fitcost track's pooled cost model (mean G, main-data
# ABOUTME: coefficients) and prints each class's share of predicted fit time per run and round.
import csv, math, sys
from collections import defaultdict
import numpy as np

C = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census"
FD = f"{C}/fitcost/data"
lam = lambda n: 4 + math.floor(3 * math.log(n))
read = lambda p: list(csv.DictReader(open(p), delimiter="\t"))

def lstsq(X, y):
    b, *_ = np.linalg.lstsq(np.asarray(X, float), np.asarray(y, float), rcond=None)
    return b, np.asarray(y, float) - np.asarray(X, float) @ b

# --- refit the pooled structural model exactly as fitcost/pooled.py does --------------------
WARM = {("bi_bi_rxn", "1"), ("ter_bi_rxn", "5")}
rows = []
for tag in ("fresh2", "alt"):
    for r in read(f"{FD}/fits_{tag}.tsv"):
        if r["error"] == "" and (r["rxn"], r["gid"]) not in WARM:
            r["tag"] = tag; rows.append(r)
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
L = np.array([lam(int(k)) for k in n], float)
(c0, c2, kk), _ = lstsq(np.c_[np.ones(len(rows)), n ** 2, L * tl], trep / G)
(a0, a1), _ = lstsq(np.c_[np.ones(len(rows)), ops], tl)
Xg = np.c_[np.ones(len(rows)), np.log(n), P, ter, alt, ter * alt]
bg, rg = lstsq(Xg, np.log(G / 20)); sg = rg.std()
(x0, x1), _ = lstsq(np.c_[np.ones(len(rows)), ops], fixed + (tpkg - trep))
print(f"model: F = {x0*1e3:.0f} ms + {x1*1e6:.0f} us*ops; g = {c0*1e6:.1f} us + {c2*1e6:.3f} us*np^2 + "
      f"{kk:.2f}*lam*t_loss; t_loss = {a0*1e6:.2f} us + {a1*1e9:.2f} ns*ops; logG = {bg[0]:.2f} + "
      f"{bg[1]:.2f} log np + {bg[2]:.2f} P + {bg[3]:.2f} ter; sd {sg:.2f}")

def t_pred(np_, dterms, isP, isTer):
    o = min(2.203 * dterms ** 1.337, 12.7 * dterms)   # ops from slice D terms (fit on the pools)
    lg = bg[0] + bg[1] * math.log(np_) + bg[2] * isP + bg[3] * isTer + sg ** 2 / 2
    return x0 + x1 * o + 20 * math.exp(lg) * (c0 + c2 * np_ ** 2 + kk * lam(np_) * (a0 + a1 * o))

# --- load census records --------------------------------------------------------------------
recs = []
for r in read(f"{C}/critic/small.tsv"):
    if r["fitted"] != "true": continue
    wk = r["Wkind"]
    recs.append(dict(run=r["run"], rnd=int(r["round"]), np=int(r["np"]), dt=int(r["dterms"]),
        P=r["P"] == "true", PL=r["Psub"] == "L_phantom", W=r["W"] == "true",
        Wisl="onlyA:island" in wk, Wg1=("onlyA:G1" in wk),
        Wface=any(x in wk for x in ("split", "trap", "bottomless")),
        D=r["D"] == "true", Z=r["Z"] == "true", N=r["N"] == "true", ter=0))
for r in read(f"{C}/critic/large.tsv"):
    if r["fitted"] != "true": continue
    lab = set(x for x in r["labels"].split(",") if x)
    recs.append(dict(run=r["run"], rnd=int(r["round"]), np=int(r["np"]), dt=int(r["dterms"]),
        P="Ptoday" in lab, PL=("Ptoday" in lab and "L_phantom" in r["explain"]),
        W=any(x.startswith("W") for x in lab), Wisl="W_admit" in lab, Wg1=False,
        Wface=bool(lab & {"W_split", "W_trap", "W_bottomless"}),
        D="D" in lab, Z="Z" in lab, N="N" in lab, ter=1))
for r in recs:
    r["t"] = t_pred(r["np"], max(r["dt"], 1), r["P"], r["ter"])

def pc(a, b): return f"{100*a/b:.1f}" if b else "-"
groups = defaultdict(list)
for r in recs: groups[(r["run"], r["rnd"])].append(r)
hdr = ("run r | n | T(s) | W | isl | D | D\\W | P | PL | PL\\(W,D) | Poth\\(W,D) | Z | N | any | "
       "mean s: OK, P, D | (count %: W, D\\W, PL\\(W,D), any)")
print(hdr)
for key in sorted(groups, key=lambda k: (k[0], k[1])):
    g = groups[key]; T = sum(r["t"] for r in g); nn = len(g)
    s = lambda f: sum(r["t"] for r in g if f(r))
    c = lambda f: sum(1 for r in g if f(r))
    anyc = lambda r: r["W"] or r["D"] or r["P"] or r["Z"] or r["N"]
    ok = [r["t"] for r in g if not anyc(r)]
    pl = [r["t"] for r in g if r["P"]]
    dd = [r["t"] for r in g if r["D"]]
    m = lambda v: f"{np.mean(v):.2f}" if v else "-"
    print(f"{key[0]} r{key[1]} | {nn} | {T:.0f} | {pc(s(lambda r: r['W']),T)} | {pc(s(lambda r: r['Wisl']),T)} | "
          f"{pc(s(lambda r: r['D']),T)} | {pc(s(lambda r: r['D'] and not r['W']),T)} | {pc(s(lambda r: r['P']),T)} | "
          f"{pc(s(lambda r: r['PL']),T)} | {pc(s(lambda r: r['PL'] and not (r['W'] or r['D'])),T)} | "
          f"{pc(s(lambda r: r['P'] and not r['PL'] and not (r['W'] or r['D'])),T)} | {pc(s(lambda r: r['Z']),T)} | "
          f"{pc(s(lambda r: r['N']),T)} | {pc(s(anyc),T)} | {m(ok)}, {m(pl)}, {m(dd)} | "
          f"({pc(c(lambda r: r['W']),nn)}, {pc(c(lambda r: r['D'] and not r['W']),nn)}, "
          f"{pc(c(lambda r: r['PL'] and not (r['W'] or r['D'])),nn)}, {pc(c(anyc),nn)})")
