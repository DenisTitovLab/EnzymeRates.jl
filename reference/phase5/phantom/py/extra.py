# ABOUTME: Extra L-phantom tables: the certificate rule (no binds-nothing pre-filter) against the spec
# ABOUTME: rule, its fit sources, and round-2 candidates reachable only through phantom parents.
import collections, sys
sys.path.insert(0, ".")
from load import *
FS = {rxn: facts(rxn) for rxn in ("bibi", "pingpong", "uniuni_R", "bibi_R")}
RS = {rxn: recs(rxn) for rxn in FS}
ES = {rxn: edges(rxn) for rxn in FS}
def cflag(r):  # certificate rule: admitted and the rescaling certificate exists
    return r["kind"] == "allo1" and r["admit"] == "ok" and r["cert"] == 1
print("== A. Certificate rule (admitted + certificate) vs spec rule (admitted + pre-filter + rank)")
for rxn in FS:
    for rd in sorted({r["round"] for r in FS[rxn].values()}):
        vs = [r for r in FS[rxn].values() if r["kind"] == "allo1" and r["round"] == rd]
        if not vs: continue
        c = collections.Counter()
        for r in vs:
            c["n"] += 1; c["spec"] += r["flag"]; c["cert"] += cflag(r)
            c["both"] += r["flag"] and cflag(r); c["spec only"] += r["flag"] and not cflag(r)
            c["cert only"] += cflag(r) and not r["flag"]
            if cflag(r):
                c["cert: rank==L0"] += r["r_var"] == r["r_0"]
                c["cert: ident_slot"] += r["ident_slot"]; c["cert: key0==key_act"] += r["key0"] == r["key_act"]
                c["cert: np==np0+1"] += r["np_slice"] == r["np0"] + 1
                c["cert: n_inact_indep=0"] += r["n_inact_indep"] == 0
        print(rxn, rd, dict(c))
print("\n== B. Where the certificate rule's phantoms find their L->0 key fitted (beam-safe relatives first)")
bykey = {}
for rxn in FS:
    d = collections.defaultdict(list)
    for r in FS[rxn].values():
        if r["kind"] == "plain" and r["status"] == "ok":
            d[r["key"]].append(r["round"])
    bykey[rxn] = d
cat = collections.Counter(); keys5 = collections.defaultdict(set); n5 = collections.Counter()
for rxn in FS:
    F, R, E = FS[rxn], RS[rxn], ES[rxn]
    for r in F.values():
        if not cflag(r): continue
        k, rd = r["key0"], r["round"]
        anc, g = [], r["parent"]
        while g != 0:
            anc.append(g); g = int(R[g]["parent"])
        rel = set(anc) | {c for a in anc for c in E.get(a, [])}
        relkeys = {F[x]["key"] for x in rel if x in F and F[x]["kind"] == "plain" and F[x]["status"] == "ok"}
        p = F.get(r["parent"])
        if p is not None and p["kind"] == "plain" and p["key"] == k: w = "1 parent"
        elif k in relkeys: w = "2 relative"
        elif any(x < rd for x in bykey[rxn].get(k, [])): w = "3 earlier round"
        elif any(x == rd for x in bykey[rxn].get(k, [])): w = "4 same round"
        else:
            w = "5 never"; keys5[(rxn, rd)].add(k); n5[(rxn, rd)] += 1
        cat[(rxn, rd, r["move"], w)] += 1
for k in sorted(cat): print(k, cat[k])
for k in sorted(n5): print("never fitted", k, n5[k], "phantoms,", len(keys5[k]), "distinct L->0 keys")
print("\n== C. Round-2 bi-bi candidates whose sampled parents are all spec-rule phantoms")
F, E = FS["bibi"], ES["bibi"]
parents = collections.defaultdict(set)
for p, cs in E.items():
    for c in cs: parents[c].add(p)
c = collections.Counter()
for r in F.values():
    if r["round"] != 2: continue
    ps = parents[r["gid"]]
    allph = all(F[p]["kind"] == "allo1" and F[p]["flag"] for p in ps)
    anyph = any(F[p]["kind"] == "allo1" and F[p]["flag"] for p in ps)
    kind = "phantom" if r.get("flag") else r["kind"]
    c[(kind, "only phantom parents" if allph else "some phantom parent" if anyph else "no phantom parent")] += 1
for k in sorted(c): print(k, c[k])
nonph_only = [r for r in F.values() if r["round"] == 2 and r["kind"] == "allo1" and not r["flag"]
              and all(F[p]["kind"] == "allo1" and F[p]["flag"] for p in parents[r["gid"]])]
print("non-phantom allosteric children reachable only through phantom parents:", len(nonph_only),
      "| their moves", dict(collections.Counter(r["move"] for r in nonph_only)))
print("sampled round-1 parents that are phantoms:", sum(1 for p in {p for r in F.values() if r["round"] == 2 for p in parents[r["gid"]]} if F[p]["kind"] == "allo1" and F[p]["flag"]))
