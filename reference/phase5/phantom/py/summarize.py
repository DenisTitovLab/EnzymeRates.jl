# ABOUTME: Prints the L-phantom study's tables: pre-filter and rank counts per population, rank method
# ABOUTME: agreement, the flagged variants' law identity and keys, numeric classes, R5.3 fit sources.
import collections, statistics, sys
sys.path.insert(0, ".")
from load import *

POPS = [("bibi", 1), ("bibi", 2), ("pingpong", 1), ("uniuni_R", 0), ("uniuni_R", 1),
        ("uniuni_R", 2), ("bibi_R", 0), ("bibi_R", 1)]
FS = {rxn: facts(rxn) for rxn in ("bibi", "pingpong", "uniuni_R", "bibi_R")}
RS = {rxn: recs(rxn) for rxn in FS}
ES = {rxn: edges(rxn) for rxn in FS}

def V(rxn, rd):
    return [r for r in FS[rxn].values() if r["kind"] == "allo1" and r["round"] == rd]

print("== 1. Populations (one-subunit allosteric candidates that pass the caps)")
print("pop round | n admit pf rankeq flag flag&cert flag&!cert | rankeq&!pf (cert/no) | num num&flag | flag_rel!=flag")
tot = collections.Counter()
for rxn, rd in POPS:
    c = collections.Counter()
    for r in V(rxn, rd):
        c["n"] += 1; c["adm"] += r["admit"] == "ok"; c["pf"] += r["pf_bind"] and r["pf_sites"]
        c["rankeq"] += r["r_var"] == r["r_0"]; c["flag"] += r["flag"]
        c["fc"] += r["flag"] and r["cert"] == 1; c["fnc"] += r["flag"] and r["cert"] == 0
        m = r["r_var"] == r["r_0"] and not (r["pf_bind"] and r["pf_sites"])
        c["mis"] += m; c["misc"] += m and r["cert"] == 1; c["misnc"] += m and r["cert"] == 0
        c["num"] += r["num"]; c["numf"] += r["num"] and r["flag"]
        c["rel"] += r["flag"] != r["flag_rel"]
        c["flag_noadm"] += r["flag_noadm"]
    tot.update(c)
    print(f"{rxn} r{rd} | {c['n']} {c['adm']} {c['pf']} {c['rankeq']} {c['flag']} {c['fc']} {c['fnc']} | "
          f"{c['mis']} ({c['misc']}/{c['misnc']}) | {c['num']} {c['numf']} | {c['rel']} | flag ignoring admission {c['flag_noadm']}")
print("total", dict(tot))

print("\n== 2. Rank methods: modular (two 61/62-bit primes) vs BigFloat (192-bit, cut 1e-30)")
allv = [r for rxn in FS for r in FS[rxn].values() if r["kind"] == "allo1"]
print("variants", len(allv), "| r_var mod!=big", sum(r["r_var"] != r["rb_var"] for r in allv),
      "| r_0 mod!=big", sum(r["r_0"] != r["rb_0"] for r in allv))
kept = [float(r["rb_kept"]) for r in allv] + [float(r["rb0_kept"]) for r in allv]
drop = [float(r["rb_drop"]) for r in allv] + [float(r["rb0_drop"]) for r in allv]
print("BigFloat gap: smallest kept pivot ratio %.1e, largest dropped %.1e" % (min(kept), max(drop)))
print("np_slice != package _independent_param_count:", sum(r["np_slice"] != r["np_pkg"] for r in allv))
tm = sorted(float(r["t_mod"]) for r in allv)
tb = sorted(float(r["t_big"]) for r in allv)
q = lambda xs, p: xs[min(len(xs) - 1, int(p * len(xs)))]
print("ms per variant, modular (4 ranks: variant, variant-L, both modes of L->0; 2 primes each): median %.1f p95 %.1f max %.1f" % (q(tm, .5), q(tm, .95), tm[-1]))
print("ms per variant, BigFloat (variant + L->0): median %.1f p95 %.1f max %.1f" % (q(tb, .5), q(tb, .95), tb[-1]))
print("relative-mode flag differs from absolute:", sum(r["flag"] != r["flag_rel"] for r in allv))
print("pf_bind != pf_empty (iso free forms):", sum(r["pf_bind"] != r["pf_empty"] for r in allv))

print("\n== 3. Flagged variants: L->0 law vs the plain catalytic law and the enumeration parent")
fl = [(rxn, r) for rxn in FS for r in FS[rxn].values() if r["kind"] == "allo1" and r["flag"]]
c = collections.Counter()
for rxn, r in fl:
    c["n"] += 1
    c["edges_same"] += r["edges_same"]; c["slot_same"] += r["slot_same"]
    c["ident_slot"] += r["ident_slot"]; c["n_inact_indep=0"] += r["n_inact_indep"] == 0
    c["key0==key_act"] += r["key0"] == r["key_act"]
    c["np=np0+1"] += r["np_slice"] == r["np0"] + 1
    c["cert"] += r["cert"]
    if r["cert"] == 1:
        c["cert_err<=1e-60"] += float(r["cert_err"]) <= 1e-60
    p = FS[rxn].get(r["parent"])
    if p is not None and p["kind"] == "plain":
        c["plain parent"] += 1
        c["key0==parent key"] += p["key"] == r["key0"]
    else:
        c["allosteric parent"] += 1
print(dict(c))
ce = [float(r["cert_err"]) for rxn in FS for r in FS[rxn].values() if r["kind"] == "allo1" and r.get("cert") == 1]
print("certificate check, 256-bit, 5 points: %d certificates, worst relative error %.1e" % (len(ce), max(ce)))
print("flagged by move:", dict(collections.Counter((rxn, r["round"], r["move"]) for rxn, r in fl)))

print("\n== 4. Numeric L-phantom candidates (L in the span of the variant's other columns), census's label")
for rxn, rd in POPS:
    c = collections.Counter()
    for r in V(rxn, rd):
        if not r["num"]: continue
        pf = r["pf_bind"] and r["pf_sites"]
        cls = "B rank>L0" if r["r_var"] > r["r_0"] else ("C rank<L0" if r["r_var"] < r["r_0"] else ("A1 eq+cert" if r["cert"] else "A2 eq,nocert"))
        c[(cls, "pf" if pf else "no-pf", "flag" if r["flag"] else "-")] += 1
    if c: print(rxn, rd, dict(sorted(c.items())))

print("\n== 5. R5.3: where each flagged phantom's L->0 key is already fitted")
def plain_keys(rxn, gids):
    F = FS[rxn]
    return {F[g]["key"] for g in gids if g in F and F[g]["kind"] == "plain" and F[g]["status"] == "ok"}
bykey = {}
for rxn in FS:
    d = collections.defaultdict(list)
    for r in FS[rxn].values():
        if r["kind"] == "plain" and r["status"] == "ok":
            d[r["key"]].append(r["round"])
    bykey[rxn] = d
cat = collections.Counter()
unf = collections.defaultdict(list)
for rxn, r in fl:
    F, R, E = FS[rxn], RS[rxn], ES[rxn]
    k, rd = r["key0"], r["round"]
    anc, g = [], r["parent"]
    while g != 0:
        anc.append(g)
        g = int(R[g]["parent"])
    rel = set(anc) | {c for a in anc for c in E.get(a, [])}
    p = F.get(r["parent"])
    if p is not None and p["kind"] == "plain" and p["key"] == k:
        where = "1 enumeration parent"
    elif k in plain_keys(rxn, rel):
        where = "2 ancestor or ancestor's child"
    elif any(x < rd for x in bykey[rxn].get(k, [])):
        where = "3 other candidate, earlier round"
    elif any(x == rd for x in bykey[rxn].get(k, [])):
        where = "4 other candidate, same round"
    else:
        where = "5 never fitted"
        unf[rxn].append(r)
    cat[(rxn, rd, r["move"], where)] += 1
for k in sorted(cat): print(k, cat[k])
for rxn, rs in unf.items():
    print("never fitted:", rxn, len(rs), "phantoms,", len({r["key0"] for r in rs}), "distinct L->0 keys; np of phantoms",
          dict(collections.Counter(r["np_slice"] for r in rs)))
    for r in rs[:3]:
        print("   e.g.", r["gid"], r["move"], "parent", r["parent"], RS[rxn][r["gid"]]["text"][:400])
