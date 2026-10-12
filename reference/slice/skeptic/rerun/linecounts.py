# ABOUTME: Code lines (no blanks, comments or docstrings) of the slice per file and per
# ABOUTME: component, before and after integration; usage: python3 linecounts.py
import os, re, sys

SLICE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC, BEFORE = os.path.join(SLICE, "src"), os.path.join(SLICE, "integrate", "src_before")


def flags(lines):
    """Per line: True for a code line (proto/p5b/codelines.py's rule, from line 1)."""
    out, indoc = [], False
    for raw in lines:
        s = raw.strip()
        q = s.count('"""')
        if indoc:
            if q % 2 == 1:
                indoc = False
            out.append(False)
            continue
        if s.startswith('"""'):
            if q % 2 == 1:
                indoc = True
            out.append(False)
            continue
        if q % 2 == 1:
            indoc = True
        out.append(bool(s) and not s.startswith("#"))
    return out


def spans(lines, name):
    """Lines of every top-level definition of `name` (function ... end, or name(...) =
    up to the next blank line)."""
    pat_f = re.compile(r"^(function|struct) " + re.escape(name) + r"\b")
    pat_s = re.compile(r"^" + re.escape(name) + r"\(")
    out = set()
    for i, l in enumerate(lines):
        if pat_f.match(l):
            j = i
            while lines[j].rstrip() != "end":
                j += 1
            out.update(range(i, j + 1))
        elif pat_s.match(l) and "=" in l:
            j = i
            while j + 1 < len(lines) and lines[j + 1].strip() != "":
                j += 1
            out.update(range(i, j + 1))
    return out


def count(path, names=None):
    lines = open(path, encoding="utf-8").read().split("\n")
    f = flags(lines)
    if names is None:
        return sum(f)
    idx = set().union(*(spans(lines, n) for n in names))
    return sum(f[i] for i in idx)


files = sorted(x for x in os.listdir(SRC) if x.endswith(".jl"))
print("per file (code lines): before integration -> after")
for x in files:
    b = count(os.path.join(BEFORE, x)) if os.path.exists(os.path.join(BEFORE, x)) else 0
    print(f"  {x:16s} {b:4d} -> {count(os.path.join(SRC, x)):4d}")
tb = sum(count(os.path.join(BEFORE, x)) for x in os.listdir(BEFORE) if x.endswith(".jl"))
print(f"  {'total':16s} {tb:4d} -> {sum(count(os.path.join(SRC, x)) for x in files):4d}")

KCAT_ALLO = ["kcat_mask", "kcat_terms", "kcat_candidates", "kcat"]
allo_all = count(os.path.join(SRC, "allosteric.jl"))
allo_kcat = count(os.path.join(SRC, "allosteric.jl"), KCAT_ALLO)
c = lambda x: count(os.path.join(SRC, x))
rows = [
    ("Engine, polynomials, size guard", 193, "215-235",
     c("engine.jl") + c("poly.jl"), "engine.jl + poly.jl"),
    ("Solve, folds, labels, all_parameters", 202, "150-200",
     c("thermo.jl") + c("api.jl"), "thermo.jl + api.jl"),
    ("Phantom handling", 0, "0-70", 0, "none (decision 12)"),
    ("Allosteric law (mask in thermo.jl)", 96, "65-90", allo_all - allo_kcat,
     "allosteric.jl without kcat"),
    (":OnlyA admission and kcat", 155, "10-270",
     c("admission.jl") + c("kcat.jl") + allo_kcat, "admission.jl + kcat.jl + allosteric kcat"),
    ("Codegen and rate_equation", 172, "125-160", c("codegen.jl"), "codegen.jl"),
    ("Printing and normalization", 72, "70-90",
     c("print.jl") + c("print_allo.jl") + c("normalize.jl"),
     "print.jl + print_allo.jl + normalize.jl"),
    ("Textbook bridge in src", 0, "0-65", 0, "none (in the tests)"),
    ("Dedup key", 6, "75-87", c("key.jl"), "key.jl"),
    ("Public API (pipeline; no rescale)", 30, "30-40", c("pipeline.jl"), "pipeline.jl"),
]
print("\nper component, against today's 926 and the critic's projection")
for name, today, proj, slice_, files_ in rows:
    print(f"  {name:38s} today {today:4d}  projected {proj:8s}  slice {slice_:4d}  ({files_})")
print(f"  {'total':38s} today {sum(r[1] for r in rows):4d}  projected {'~970':8s}  "
      f"slice {sum(r[3] for r in rows):4d}")
print("\noutside the derivation")
print(f"  types.jl (type, constructor, orientation, RE forest and blocks, checks, Sig): "
      f"{c('types.jl')} (today's types.jl 746 at the probes' baseline, 770 on main)")
print(f"  names.jl (naming chokepoint): {c('names.jl')}")
print(f"  enum.jl (three moves and the expand driver): {c('enum.jl')}")
print(f"  convert.jl (Pop bridge, test only): {c('convert.jl')}; Slice.jl: {c('Slice.jl')}")
