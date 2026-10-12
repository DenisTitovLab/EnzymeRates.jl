# ABOUTME: Second pass over the cost-weighted census shares: island-only (not D, not P) time, and
# ABOUTME: the same incremental shares when each loss call costs 10x (480 data points instead of 48).
import runpy, math, io, contextlib
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ns = runpy.run_path("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/critic/shares.py")
recs = ns["recs"]; bg, sg, x0, x1, c0, c2, kk, a0, a1, lam = (ns[k] for k in ("bg","sg","x0","x1","c0","c2","kk","a0","a1","lam"))
def t480(np_, dt, isP, isTer):
    o = min(2.203 * dt ** 1.337, 12.7 * dt)
    lg = bg[0] + bg[1] * math.log(np_) + bg[2] * isP + bg[3] * isTer + sg ** 2 / 2
    return x0 + x1 * o + 20 * math.exp(lg) * (c0 + c2 * np_ ** 2 + kk * lam(np_) * 10 * (a0 + a1 * o))
from collections import defaultdict
g = defaultdict(list)
for r in recs:
    r["t480"] = t480(r["np"], max(r["dt"], 1), r["P"], r["ter"]); g[(r["run"], r["rnd"])].append(r)
print("run r | W-only (not D,P) time% 48pt | 480pt: W | D\\W | PL\\(W,D) | Poth\\(W,D) | any | T48 s | T480 s")
for k in sorted(g):
    if k[1] == 0 or k[0].startswith("uni_uni") or (k[0] in ("bi_bi_rxn","bi_bi_pp_rxn") and k[1]==2) : pass
    G = g[k]; T = sum(r["t"] for r in G); T4 = sum(r["t480"] for r in G)
    s = lambda f, key="t480": 100 * sum(r[key] for r in G if f(r)) / (T4 if key == "t480" else T)
    print(f"{k[0]} r{k[1]} | {s(lambda r: r['W'] and not r['D'] and not r['P'], 't'):.1f} | {s(lambda r: r['W']):.1f} | "
          f"{s(lambda r: r['D'] and not r['W']):.1f} | {s(lambda r: r['PL'] and not (r['W'] or r['D'])):.1f} | "
          f"{s(lambda r: r['P'] and not r['PL'] and not (r['W'] or r['D'])):.1f} | "
          f"{s(lambda r: r['W'] or r['D'] or r['P'] or r['Z'] or r['N']):.1f} | {T:.0f} | {T4:.0f}")
