# ABOUTME: Summarizes the verifier's phantom recheck (ph_*.tsv): flag counts, rank agreement, the
# ABOUTME: certificate on flagged variants, fits removed per rule and each phantom's fit source.
import sys, glob, collections
rows = {}
for p in sorted(glob.glob(sys.argv[1] + "/ph_*.tsv")):
    for l in open(p):
        d = dict(x.split("=", 1) for x in l.rstrip("\n").split("\t") if "=" in x)
        rows[int(d["gid"])] = d
B = lambda d, k: d.get(k) == "true"
print("rows", len(rows), "errors", sum(d["status"] != "ok" for d in rows.values()))
for d in rows.values():
    if d["status"] != "ok": print("  error", d["gid"], d.get("err", "")[:150])
v = [d for d in rows.values() if d["round"] == "2" and d["allo"] == "true" and d["n"] == "1" and d["status"] == "ok"]
for d in v:
    d["pf"] = B(d, "pf_bind") and B(d, "pf_sites")
    d["rankeq"] = d["r_var"] == d["r_0"]
    d["flag"] = B(d, "admit") and d["pf"] and d["rankeq"]
    d["certb"] = B(d, "cert")
c = collections.Counter()
for d in v:
    c["n"] += 1; c["admit"] += B(d, "admit"); c["pf"] += d["pf"]; c["rankeq"] += d["rankeq"]
    c["flag"] += d["flag"]; c["flag&cert"] += d["flag"] and d["certb"]; c["flag&!cert"] += d["flag"] and not d["certb"]
    c["rankeq&!pf"] += d["rankeq"] and not d["pf"]; c["rankeq&!pf&cert"] += d["rankeq"] and not d["pf"] and d["certb"]
    c["admit&cert"] += B(d, "admit") and d["certb"]; c["admit&cert&!pf"] += B(d, "admit") and d["certb"] and not d["pf"]
    c["cert&rank>0"] += d["certb"] and int(d["r_var"]) > int(d["r_0"])
    c["pf_bind!=pf_empty"] += d["pf_bind"] != d["pf_empty"]
    c["draws disagree"] += not B(d, "draws_agree")
    c["rank>L0"] += int(d["r_var"]) > int(d["r_0"]); c["rank<L0"] += int(d["r_var"]) < int(d["r_0"])
print("round-2 one-subunit variants:", dict(c))
kept = min(float(d["kept"]) for d in v); dropped = max(float(d["dropped"]) for d in v)
print(f"BigFloat pivots: smallest kept {kept:.1e}, largest dropped {dropped:.1e}")
fl = [d for d in v if d["flag"]]
print("flagged:", len(fl), "| key_act==key0", sum(d["key_act"] == d["key0"] for d in fl),
      "| np==np0+1", sum(int(d["np"]) == int(d["np0"]) + 1 for d in fl),
      "| cert", sum(d["certb"] for d in fl),
      "| worst cert map error", max((float(d["cert_err"]) for d in fl if "cert_err" in d), default=None))
print("flagged without certificate, why:", collections.Counter(d["cert_why"] for d in fl if not d["certb"]),
      [(d["gid"], d["move"]) for d in fl if not d["certb"]])
ce = [float(d["cert_err"]) for d in v if "cert_err" in d]
print("certificates checked numerically:", len(ce), "worst", max(ce, default=None))
print("flagged by move:", collections.Counter(d["move"] for d in fl))
# fits removed in round 2 (key dedup first), as phantom/py/fits.py defines it
seen = set()
for rd in ("0", "1", "2"):
    rs = [d for d in rows.values() if d["round"] == rd and d["status"] == "ok"]
    new = {d["key"] for d in rs} - seen
    out = {}
    for name, rule in (("spec", lambda d: d.get("flag", False)),
                       ("cert", lambda d: d.get("certb", False) and B(d, "admit"))):
        ph = {d["key"] for d in rs if rule(d)} & new
        non = {d["key"] for d in rs if not rule(d)}
        out[name] = (len(ph - non), len(ph & non))
    print(f"round {rd}: {len(rs)} candidates, {len(new)} new keys; fits removed: spec rule {out['spec'][0]}, certificate rule {out['cert'][0]} (shared with non-phantoms {out['spec'][1]}, {out['cert'][1]})")
    seen |= {d["key"] for d in rs}
# R5.3: where each flagged phantom's L->0 key is fitted
par = {g: int(d["parent"]) for g, d in rows.items()}
kids = collections.defaultdict(list)
for g, p in par.items(): kids[p].append(g)
plain = {g: d for g, d in rows.items() if d["allo"] == "false" and d["status"] == "ok"}
src = collections.Counter(); never = []
for d in fl:
    g, k0 = int(d["gid"]), d["key0"]
    anc, p = [], par[g]
    while p: anc.append(p); p = par.get(p, 0)
    rel = set(anc) | {x for a in anc for x in kids[a]}
    p0 = par[g]
    if p0 in plain and plain[p0]["key"] == k0: cls = "1 enumeration parent"
    elif any(x in plain and plain[x]["key"] == k0 and plain[x]["round"] in ("0", "1") for x in rel): cls = "2 earlier-round relative"
    elif any(x["key"] == k0 and x["round"] in ("0", "1") for x in plain.values()): cls = "3 earlier round, not a relative"
    elif any(x["key"] == k0 for x in plain.values()): cls = "4 same round only"
    else: cls = "5 never"; never.append((g, d["move"], d["np"], k0[:12]))
    src[(d["move"], cls)] += 1
for k, n in sorted(src.items()): print("  source", k, n)
print("never fitted:", len(never), "distinct L->0 keys", len({x[3] for x in never}), never[:10])
