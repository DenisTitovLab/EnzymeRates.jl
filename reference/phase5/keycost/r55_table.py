# ABOUTME: R5.5 table: one line per reaction, sampling mode and round with the candidates each
# ABOUTME: complexity-filter rule newly drops. Args: data/r55_<rxn>_<mode>.tsv ...
import sys
from collections import defaultdict

F = 337
print("file round | cand allo | admitted-today | new:sum new:distinct new:max | new:sum>674 "
      "plain(337,674] | maxA-admitted maxSum-admitted | allo-with-A>168")
for path in sys.argv[1:]:
    by = defaultdict(list)
    for line in open(path):
        f = line.split("\t")
        by[int(f[0])].append((f[2] == "1", int(f[4]), int(f[5]), f[6] == "1"))
    tag = path.split("r55_")[1].replace(".tsv", "")
    for rd in sorted(by):
        rs = by[rd]
        adm = [r for r in rs if r[1] <= F]
        new = lambda g: sum(1 for r in adm if g(r) > F)
        s = new(lambda r: r[1] + r[2])
        d = new(lambda r: r[1] + (0 if r[3] else r[2]))
        m = new(lambda r: max(r[1], r[2]))
        s2 = sum(1 for r in adm if r[1] + r[2] > 2 * F)
        pm = sum(1 for r in rs if not r[0] and F < r[1] <= 2 * F)
        ma = max((r[1] for r in adm), default=0)
        ms = max((r[1] + r[2] for r in adm), default=0)
        big = sum(1 for r in adm if r[0] and r[1] > 168)
        print(f"{tag} r{rd} | {len(rs)} {sum(r[0] for r in rs)} | {len(adm)} | {s} {d} {m} | "
              f"{s2} {pm} | {ma} {ms} | {big}")
