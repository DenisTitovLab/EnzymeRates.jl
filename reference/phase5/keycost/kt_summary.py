# ABOUTME: Summaries of keytime.jl output: key and law time by law size (R5.6), and codegen ns per
# ABOUTME: term under three term counts (R5.7). Args: <kt.tsv> [<kt.tsv> ...].
import math
import sys
from collections import defaultdict


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else float("nan")


def load(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        r = dict(set=f[0], label=f[1], id=f[2], cxA=int(f[5]), cxI=int(f[6]), status=f[7])
        if f[7] == "ok":
            r.update(D=int(f[8]), KT=int(f[9]), np=int(f[10]), tlaw=float(f[11]),
                     tkey=float(f[12]), leaves=int(f[13]), nodes=int(f[14]), key=f[15],
                     qmax=int(f[16]), n=int(f[17]))
            if len(f) > 18:
                r.update(tprog=float(f[18]), tcomp=float(f[19]), ns=float(f[20]),
                         nsi=float(f[21]), allocs=int(f[22]), ops=int(f[23]), body=f[24])
        out.append(r)
    return out


def slope(xs, ys):
    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)


for path in sys.argv[1:]:
    rows = load(path)
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"{path}: {len(rows)} rows, {len(ok)} laws, {len(rows) - len(ok)} errors")
    errs = defaultdict(int)
    for r in rows:
        if r["status"] != "ok":
            errs[r["status"][:90]] += 1
    for e, c in errs.items():
        print(f"  error x{c}: {e}")
    bins = [(0, 50), (51, 337), (338, 1000), (1001, 2000), (2001, 5000), (5001, 10**9)]
    print("  by key-part terms (KT): n | key ms median p90 max | law ms median max | "
          "key/law median max | leaves max | nodes max")
    for lo, hi in bins:
        b = [r for r in ok if lo <= r["KT"] <= hi]
        if not b:
            continue
        kt = [1e3 * r["tkey"] for r in b]
        lt = [1e3 * r["tlaw"] for r in b]
        ratio = [r["tkey"] / r["tlaw"] for r in b]
        print(f"  KT [{lo}, {hi}]: {len(b)} | {q(kt, .5):.2f} {q(kt, .9):.2f} {max(kt):.2f} | "
              f"{q(lt, .5):.1f} {max(lt):.1f} | {q(ratio, .5):.2f} {max(ratio):.2f} | "
              f"{max(r['leaves'] for r in b)} | {max(r['nodes'] for r in b)}")
    big = [r for r in ok if r["KT"] > 50]
    if len(big) > 5:
        print(f"  log-log slope of key time on KT (KT > 50): {slope([r['KT'] for r in big], [r['tkey'] for r in big]):.2f}; "
              f"of law time: {slope([r['KT'] for r in big], [r['tlaw'] for r in big]):.2f}")
    print(f"  leaves > 1 in {sum(r['leaves'] > 1 for r in ok)} laws; max leaves "
          f"{max((r['leaves'] for r in ok), default=0)}; exponent denominators > 2 in "
          f"{sum(r['qmax'] > 2 for r in ok)}, = 2 in {sum(r['qmax'] == 2 for r in ok)}")
    allo = [r for r in ok if r["n"] > 0]
    if allo:
        kt = [1e3 * r["tkey"] for r in allo]
        print(f"  allosteric: {len(allo)} laws, key ms median {q(kt, .5):.2f} max {max(kt):.2f}; "
              f"max KT {max(r['KT'] for r in allo)}")
    top = sorted(ok, key=lambda r: -r["tkey"])[:5]
    for r in top:
        print(f"  slowest key: {r['set']} {r['label']} {r['id']} KT {r['KT']} D {r['D']} cx {r['cxA']}+{r['cxI']} "
              f"key {1e3 * r['tkey']:.1f} ms law {1e3 * r['tlaw']:.1f} ms leaves {r['leaves']}")
    code = [r for r in ok if "ns" in r]
    if code:
        print(f"  codegen: {len(code)} laws; allocations > 0 in {sum(r['allocs'] > 0 for r in code)}; "
              f"ns barrier max {max(r['ns'] for r in code):.2f}; over 120 ns: "
              f"{sum(r['ns'] > 120 for r in code)}")
        for name, key in (("D", "D"), ("V×τ", "cx"), ("ops", "ops")):
            val = (lambda r: r["cxA"] + r["cxI"]) if key == "cx" else (lambda r, k=key: r[k])
            sel = [r for r in code if val(r) >= 300]
            if not sel:
                continue
            per = [r["ns"] / val(r) for r in sel]
            print(f"    ns per {name:4s} over {len(sel)} laws with {name} >= 300: min {min(per):.4f} "
                  f"median {q(per, .5):.4f} max {max(per):.4f} (max/min {max(per) / min(per):.2f})")
        sel = [r for r in code if r["D"] >= 50]
        if len(sel) > 5:
            print(f"    log-log slope of ns on D (D >= 50): {slope([r['D'] for r in sel], [r['ns'] for r in sel]):.2f}; "
                  f"on ops: {slope([r['ops'] for r in sel], [r['ns'] for r in sel]):.2f}; on V×τ: "
                  f"{slope([r['cxA'] + r['cxI'] for r in sel], [r['ns'] for r in sel]):.2f}")
