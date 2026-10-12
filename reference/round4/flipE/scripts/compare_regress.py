# ABOUTME: Task 4: compares the regression TSVs of the prototype and the unmodified package:
# ABOUTME: mechanism identity, fitted parameters, equation-string hashes and rate values.
import math
import sys
from collections import Counter


def load(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f]


proto, base = load(sys.argv[1]), load(sys.argv[2])
print(f"rows: prototype {len(proto)}, baseline {len(base)}")
assert len(proto) == len(base)
same = Counter()
max_abs, max_rel, worst = 0.0, 0.0, None
for p, b in zip(proto, base):
    for col in ("key_hash", "n_free", "tags", "n_pruned_groups", "fitted_params", "eqstring_hash"):
        same[col] += p[col] == b[col]
    for col in ("v_generic", "v_products0", "v_low", "kcat"):
        x, y = float(p[col]), float(b[col])
        if math.isnan(x) and math.isnan(y):
            same[col + " (both NaN)"] += 1
            continue
        d = abs(x - y)
        r = d / max(abs(y), 1e-300)
        if d > max_abs:
            max_abs = d
        if r > max_rel:
            max_rel, worst = r, (p["i"], col, x, y)
        same[col + " bitwise"] += x == y
for col, n in same.items():
    print(f"  {col}: {n}/{len(proto)} equal")
print(f"max |prototype - baseline| over the values: {max_abs:.3e}; max relative: {max_rel:.3e}"
      + (f" (row {worst[0]}, {worst[1]}: {worst[2]!r} vs {worst[3]!r})" if worst else ""))
print("free forms per mechanism:", dict(Counter(r["n_free"] for r in proto)))
print("tag sets:", dict(Counter(r["tags"] for r in proto)))
print("mechanisms whose inactive state strands a non-OnlyA group:",
      sum(int(r["n_pruned_groups"]) > 0 for r in proto))
