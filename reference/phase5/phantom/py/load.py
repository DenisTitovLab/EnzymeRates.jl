# ABOUTME: Loads the L-phantom study's records, per-candidate facts and edges into dicts, and derives
# ABOUTME: the spec rule's flag, the numeric flag and the plain key sets the summaries read.
import os, collections
D = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/phantom/data"

def unesc(s):
    return s.replace("\\t", "\t").replace("\\n", "\n").replace("\\\\", "\\")

def recs(rxn):
    out = {}
    with open(os.path.join(D, f"recs_{rxn}.tsv")) as f:
        head = f.readline().rstrip("\n").split("\t")
        for l in f:
            r = dict(zip(head, [unesc(x) for x in l.rstrip("\n").split("\t")]))
            out[int(r["gid"])] = r
    return out

def facts(rxn):
    out = {}
    with open(os.path.join(D, f"facts_{rxn}.tsv")) as f:
        for l in f:
            x = l.rstrip("\n").split("\t")
            g = int(x[0])
            r = {"gid": g, "kind": x[1], "round": int(x[2]), "move": x[3], "parent": int(x[4]),
                 "np_pkg": int(x[5])}
            if x[1] == "plain":
                r.update(status=x[6], np_slice=x[7], key=x[8])
            else:
                for kv in x[6:]:
                    k, _, v = kv.partition("=")
                    r[k] = v
                for k in ("r_var", "r_noL", "r_0", "r_var_rel", "r_0_rel", "rb_var", "rb_0",
                          "np_slice", "np0", "pf_bind", "pf_empty", "pf_sites", "cert",
                          "ident_slot", "slot_same", "edges_same", "n_inact_indep", "nsites"):
                    if k in r and r[k] != "":
                        r[k] = int(r[k])
                r["flag"] = (r.get("admit") == "ok" and r["pf_bind"] == 1 and r["pf_sites"] == 1
                             and r["r_var"] == r["r_0"])
                r["flag_noadm"] = (r["pf_bind"] == 1 and r["pf_sites"] == 1
                                   and r["r_var"] == r["r_0"])
                r["flag_rel"] = (r.get("admit") == "ok" and r["pf_bind"] == 1 and r["pf_sites"] == 1
                                 and r["r_var_rel"] == r["r_0_rel"])
                r["num"] = r["r_var"] < r["np_slice"] and r["r_noL"] == r["r_var"]
            out[g] = r
    return out

def edges(rxn):
    ch = collections.defaultdict(list)
    with open(os.path.join(D, f"edges_{rxn}.tsv")) as f:
        for l in f:
            p, c, mv = l.rstrip("\n").split("\t")
            ch[int(p)].append(int(c))
    return ch
