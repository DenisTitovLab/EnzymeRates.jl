# ABOUTME: R5.5 counts: per reaction and round, candidates the complexity filter newly drops when
# ABOUTME: both conformations count, under four candidate rules. Args: data/r55_<rxn>.tsv ...
import sys
from collections import defaultdict

F = 337


def rows(path):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        yield dict(round=int(f[0]), allo=f[2] == "1", n=int(f[3]), a=int(f[4]), i=int(f[5]),
                   same=f[6] == "1")


RULES = {
    "sum": lambda r: r["a"] + r["i"],
    "distinct": lambda r: r["a"] + (0 if r["same"] else r["i"]),
    "max": lambda r: max(r["a"], r["i"]),
}

for path in sys.argv[1:]:
    by = defaultdict(list)
    for r in rows(path):
        by[r["round"]].append(r)
    print(path)
    for rd in sorted(by):
        rs = by[rd]
        al = [r for r in rs if r["allo"]]
        today = [r for r in rs if r["a"] <= F]
        print(f"  round {rd}: {len(rs)} candidates, {len(al)} allosteric "
              f"(n=1 {sum(r['n'] == 1 for r in al)}, n=2 {sum(r['n'] == 2 for r in al)}); "
              f"within {F} today {len(today)} (allosteric {sum(r['allo'] for r in today)})")
        for name, f in RULES.items():
            new = [r for r in today if f(r) > F]
            print(f"    {name:8s} > {F}: newly over {len(new)} = {100 * len(new) / max(len(today), 1):.1f}% "
                  f"of today's admitted; n=1 {sum(r['n'] == 1 for r in new)}, n=2 "
                  f"{sum(r['n'] == 2 for r in new)}; same-law {sum(r['same'] for r in new)}")
        s2 = [r for r in today if r["a"] + r["i"] > 2 * F]
        plain_mid = [r for r in rs if not r["allo"] and F < r["a"] <= 2 * F]
        allo_mid = [r for r in rs if r["allo"] and F < r["a"] + r["i"] <= 2 * F and r["a"] > F]
        print(f"    sum > {2 * F}: newly over {len(s2)}; newly admitted plain (a in ({F}, {2 * F}]) "
              f"{len(plain_mid)}, allosteric with a > {F} and sum <= {2 * F}: {len(allo_mid)}")
        ia = [r["i"] / r["a"] for r in al if r["a"] > 0]
        if ia:
            ia.sort()
            q = lambda p: ia[min(len(ia) - 1, int(p * len(ia)))]
            print(f"    inactive/active V×τ over allosteric: median {q(0.5):.2f}, p10 {q(0.1):.2f}, "
                  f"p90 {q(0.9):.2f}, max {ia[-1]:.2f}; inactive > active in "
                  f"{sum(r['i'] > r['a'] for r in al)}; inactive 0 in {sum(r['i'] == 0 for r in al)}; "
                  f"same-law {sum(r['same'] for r in al)}")
        big = sorted(r["a"] for r in today)
        print(f"    today's admitted V×τ: max {big[-1] if big else 0}; in (168, 337]: "
              f"{sum(168 < x <= 337 for x in big)}")
