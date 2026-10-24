# ABOUTME: Picks the law study's census sample (specs, then per reaction init and smallest-id children)
# ABOUTME: from a census corpus directory (argv 1) and writes set, id, roles, regulators and oriented text as one TSV.
import sys

B = sys.argv[1]
OUT = sys.argv[2]

def rows(name):
    with open(f"{B}/{name}", encoding="utf-8") as f:
        return [l.rstrip("\n").split("\t") for l in f]

def unescape(s):
    out, i = [], 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s):
            out.append({"n": "\n", "t": "\t", "\\": "\\"}.get(s[i + 1], s[i + 1]))
            i += 2
        else:
            out.append(s[i])
            i += 1
    return "".join(out)

def names(text, key):
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(key + ":"):
            body = line[len(key) + 1:]
            return [x.split("[")[0].split("::")[0].strip() for x in body.split(",") if x.strip()]
    return []

reactions = {r[0]: unescape(r[4]) for r in rows("reactions.tsv")}
roles = {}
regs = {}
for rx, text in reactions.items():
    subs, prods = names(text, "substrates"), names(text, "products")
    mods = names(text, "competitive_inhibitors")
    roles[rx] = "S=%s;P=%s;M=%s" % (",".join(subs), ",".join(prods), ",".join(mods))
    regs[rx] = ",".join(names(text, "allosteric_regulators"))

init = {}
for r in rows("init.tsv"):
    init.setdefault(r[0], []).append(r[2])
kids = {}
for r in rows("children.tsv"):
    kids.setdefault((r[0], int(r[1])), set()).add(r[4])

# (reaction, init quota (None = all), [(parent level, quota)])
PLAN = [("bibi", None, [(0, 361)]),
        ("pingpong", None, [(0, 409), (1, 200)]),
        ("twosites", None, [(0, 60), (1, 100)]),
        ("uniuni", None, [(0, 60), (1, 100)]),
        ("terbi", 200, [(0, 100)]),
        ("prefix", 100, [(0, 100)])]

want = []          # (set, id, roles, regs)
for r in rows("laws.tsv"):
    if r[1] == "spec":
        want.append(("spec:" + r[0][len("spec:"):], r[2], r[3], ""))
for rx, iq, levels in PLAN:
    seen = set()
    ids = sorted(set(init[rx]))
    for i in ids[:iq] if iq is not None else ids:
        seen.add(i)
        want.append((rx + ":init", i, roles[rx], regs[rx]))
    for lev, q in levels:
        pool = sorted(kids.get((rx, lev), set()) - seen)[:q]
        for i in pool:
            seen.add(i)
            want.append((f"{rx}:L{lev + 1}", i, roles[rx], regs[rx]))

need = {w[1] for w in want}
text = {}
with open(f"{B}/mechs.tsv", encoding="utf-8") as f:
    for line in f:
        i, _, rest = line.partition("\t")
        if i in need:
            text[i] = rest.rstrip("\n").split("\t")[1]
with open(OUT, "w", encoding="utf-8") as f:
    f.write("set\tid\troles\tregs\toriented\n")
    for s, i, ro, rg in want:
        f.write(f"{s}\t{i}\t{ro}\t{rg}\t{text[i]}\n")
counts = {}
for s, *_ in want:
    counts[s.split(":")[0] if s.startswith("spec:") else s] = counts.get(s.split(":")[0] if s.startswith("spec:") else s, 0) + 1
print(len(want), "rows", counts)
