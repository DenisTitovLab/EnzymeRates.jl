# ABOUTME: Independent verification of topology counts
# ABOUTME: for simple ter-ter, pyruvate carboxylase, PDH.

from itertools import combinations, permutations
from math import factorial

def fubini(n):
    if n <= 1: return 1
    total = 0
    for k in range(1, n+1):
        s = 0
        for j in range(k+1):
            s += ((-1)**j * factorial(k)
                  // (factorial(j) * factorial(k-j))
                  * (k-j)**n)
        total += s
    return total

def can_pp(acc, need):
    for a, c in need.items():
        if acc.get(a, 0) < c:
            return False
    return True

def subtract(acc, need):
    res = dict(acc)
    for a, c in need.items():
        res[a] = res.get(a, 0) - c
        if res[a] <= 0:
            del res[a]
    return res

def add_atoms(acc, more):
    res = dict(acc)
    for a, c in more.items():
        res[a] = res.get(a, 0) + c
    return res

def prod_atoms_for_set(ps, patoms):
    result = {}
    for p in ps:
        for a, c in patoms[p].items():
            result[a] = result.get(a, 0) + c
    return result

def count_topologies(name, subs, prods, satoms, patoms):
    print(f"\n{'='*60}")
    print(f"{name}")
    print(f"{'='*60}")

    # 1. SEQUENTIAL (3->3)
    seq_topos = fubini(len(subs)) * fubini(len(prods))
    print(f"\nSequential: 1 pattern, {seq_topos} topologies")

    # 2. BI-UNI (2 subs first, 1 sub second)
    bi_uni_count = 0
    bi_uni_topos = 0
    for pair in combinations(subs, 2):
        rsub = [s for s in subs if s not in pair][0]
        acc = {}
        for s in pair:
            acc = add_atoms(acc, satoms[s])

        for k in [1, 2]:
            for ps in combinations(prods, k):
                need = prod_atoms_for_set(ps, patoms)
                if not can_pp(acc, need):
                    continue
                residual = subtract(acc, need)
                rprods = [p for p in prods if p not in ps]
                n_pe = len(ps) + (1 if residual else 0)
                if n_pe > 3 or len(pair) > 3:
                    continue
                # Second iso
                acc2 = add_atoms(dict(residual),
                                 satoms[rsub])
                need2 = prod_atoms_for_set(rprods, patoms)
                if not can_pp(acc2, need2):
                    continue
                res2 = subtract(acc2, need2)
                if res2:  # must be exact
                    continue
                bi_uni_count += 1
                ords = (fubini(len(pair)) * fubini(len(ps))
                        * fubini(1)
                        * fubini(len(rprods)))
                bi_uni_topos += ords
                print(f"  bi-uni: {'+'.join(pair)}->"
                      f"{'+'.join(ps)}, "
                      f"{rsub}->{'+'.join(rprods)}"
                      f"  ({ords} orderings)")

    print(f"  Bi-uni total: {bi_uni_count} patterns, "
          f"{bi_uni_topos} topologies")

    # 3. UNI-BI (1 sub first, 2 subs second)
    uni_bi_count = 0
    uni_bi_topos = 0
    for s1 in subs:
        rsubs = [s for s in subs if s != s1]
        acc1 = dict(satoms[s1])

        for k in [1, 2]:
            for ps in combinations(prods, k):
                need = prod_atoms_for_set(ps, patoms)
                if not can_pp(acc1, need):
                    continue
                residual = subtract(acc1, need)
                rprods = [p for p in prods if p not in ps]
                n_pe = len(ps) + (1 if residual else 0)
                if n_pe > 3:
                    continue
                acc2 = dict(residual)
                for s in rsubs:
                    acc2 = add_atoms(acc2, satoms[s])
                need2 = prod_atoms_for_set(rprods, patoms)
                if not can_pp(acc2, need2):
                    continue
                res2 = subtract(acc2, need2)
                if res2:
                    continue
                n_pe2 = len(rprods)
                if n_pe2 > 3 or len(rsubs) > 3:
                    continue
                uni_bi_count += 1
                ords = (fubini(1) * fubini(len(ps))
                        * fubini(len(rsubs))
                        * fubini(len(rprods)))
                uni_bi_topos += ords
                print(f"  uni-bi: {s1}->{'+'.join(ps)}, "
                      f"{'+'.join(rsubs)}->"
                      f"{'+'.join(rprods)}"
                      f"  ({ords} orderings)")

    print(f"  Uni-bi total: {uni_bi_count} patterns, "
          f"{uni_bi_topos} topologies")

    # 4. HEXA-UNI (1 sub each, 3 iso steps)
    # Each ordering of subs creates different iso FORMS
    # (E_A vs Estar_A), so each valid ordering is a
    # distinct iso pattern.
    hexa_count = 0
    for perm in permutations(subs):
        for pperm in permutations(prods):
            acc = {}
            valid = True
            for i in range(3):
                acc = add_atoms(acc, satoms[perm[i]])
                need = patoms[pperm[i]]
                if not can_pp(acc, need):
                    valid = False
                    break
                acc = subtract(acc, need)
            if valid and not acc:
                hexa_count += 1
    print(f"\n  Hexa-uni: {hexa_count} patterns, "
          f"{hexa_count} topologies (1 ordering each)")

    total = seq_topos + bi_uni_topos + uni_bi_topos + hexa_count
    pp = bi_uni_topos + uni_bi_topos + hexa_count
    print(f"\n  TOTAL: {total} topologies")
    print(f"    Sequential: {seq_topos}")
    print(f"    Ping-pong: {pp}")
    return total, seq_topos, pp

# Simple ter-ter
count_topologies("Simple ter-ter: A[C]+B[N]+D[X]",
    ['A', 'B', 'D'], ['P', 'Q', 'R'],
    {'A': {'C':1}, 'B': {'N':1}, 'D': {'X':1}},
    {'P': {'C':1}, 'Q': {'N':1}, 'R': {'X':1}})

# Pyruvate carboxylase
count_topologies("Pyruvate carboxylase",
    ['Pyr', 'HCO3', 'ATP'], ['OAA', 'ADP', 'Pi'],
    {'Pyr': {'C':3,'H':3,'O':3},
     'HCO3': {'H':1,'C':1,'O':3},
     'ATP': {'C':10,'H':16,'N':5,'O':13,'P':3}},
    {'OAA': {'C':4,'H':3,'O':5},
     'ADP': {'C':10,'H':15,'N':5,'O':10,'P':2},
     'Pi': {'H':2,'P':1,'O':4}})

# Pyruvate dehydrogenase
count_topologies("Pyruvate dehydrogenase",
    ['Pyr', 'NAD', 'CoA'], ['AcCoA', 'NADH', 'CO2'],
    {'Pyr': {'C':3,'H':3,'O':3},
     'NAD': {'C':21,'H':28,'N':7,'O':14,'P':2},
     'CoA': {'C':21,'H':36,'N':7,'O':16,'P':3,'S':1}},
    {'AcCoA': {'C':23,'H':38,'N':7,'O':17,'P':3,'S':1},
     'NADH': {'C':21,'H':29,'N':7,'O':14,'P':2},
     'CO2': {'C':1,'O':2}})
