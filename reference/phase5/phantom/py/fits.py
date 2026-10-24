# ABOUTME: Fits per round with identify's key cache, and the fits each phantom rule removes: distinct
# ABOUTME: keys of phantoms not fitted earlier, against all distinct new keys of the round.
import collections, sys
sys.path.insert(0, ".")
from load import *
for rxn in ("bibi", "pingpong"):
    F = facts(rxn)
    ak = {}
    for l in open(f"{D}/allokeys_{rxn}.tsv"):
        g, k = l.split(); ak[int(g)] = k
    seen = set()
    for rd in sorted({r["round"] for r in F.values()}):
        rs = [r for r in F.values() if r["round"] == rd and (r["kind"] == "allo1" or r.get("status") == "ok")]
        key = lambda r: r["key"] if r["kind"] == "plain" else ak[r["gid"]]
        new = {key(r) for r in rs} - seen
        out = {}
        for name, rule in (("spec", lambda r: r["kind"] == "allo1" and r["flag"]),
                           ("cert", lambda r: r["kind"] == "allo1" and r["admit"] == "ok" and r["cert"] == 1)):
            ph = {key(r) for r in rs if rule(r)} & new
            nonph = {key(r) for r in rs if not rule(r)}
            out[name] = (len(ph - nonph), len(ph & nonph))
        adm = [r for r in rs if r["kind"] == "plain" or r["admit"] == "ok"]
        newadm = {key(r) for r in adm} - seen
        print(f"{rxn} round {rd}: {len(rs)} candidates, {len(new)} new keys ({len(newadm)} admitted); "
              f"fits removed: spec rule {out['spec'][0]}, certificate rule {out['cert'][0]} "
              f"(keys shared with a non-phantom: {out['spec'][1]}, {out['cert'][1]})")
        seen |= {key(r) for r in rs}
