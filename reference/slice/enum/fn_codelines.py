# ABOUTME: Counts code lines (no blanks, comments or docstrings) of named top-level Julia
# ABOUTME: functions, every method of each name. Usage: python3 fn_codelines.py FILE NAME...
import re, sys

def code_flags(lines):
    """Per line: True for a code line (the codelines.py rule, scanning from line 1)."""
    out, indoc = [], False
    for raw in lines:
        s = raw.strip()
        q = s.count('"""')
        if indoc:
            if q % 2 == 1: indoc = False
            out.append(False); continue
        if s.startswith('"""'):
            if q % 2 == 1: indoc = True
            out.append(False); continue
        if q % 2 == 1: indoc = True
        out.append(bool(s) and not s.startswith("#"))
    return out

def spans(lines, name):
    """Line spans of every top-level definition of `name`: `function name(` up to the
    next column-0 `end`, or `name(...) =` up to the next blank line (overlapping spans
    of consecutive one-line methods are merged by the caller)."""
    out = []
    pat_f = re.compile(r'^function ' + re.escape(name) + r'\(')
    pat_s = re.compile(r'^' + re.escape(name) + r'\(.*\)\s*(where.*)?=')
    pat_s2 = re.compile(r'^' + re.escape(name) + r'\(')
    for i, l in enumerate(lines):
        if pat_f.match(l):
            j = i
            while lines[j].rstrip() != "end": j += 1
            out.append((i, j))
        elif pat_s.match(l) or (pat_s2.match(l) and "=" in l):
            j = i
            while j + 1 < len(lines) and lines[j + 1].strip() != "": j += 1
            out.append((i, j))
    return out

path, names = sys.argv[1], sys.argv[2:]
lines = open(path).read().split("\n")
flags = code_flags(lines)
total = 0
for n in names:
    ss = spans(lines, n)
    # a one-line definition's span runs to the next blank line, so count each line once
    c = sum(flags[i] for i in sorted({i for a, b in ss for i in range(a, b + 1)}))
    total += c
    print(f"{n}: {c} ({len(ss)} definitions at lines {[a + 1 for a, _ in ss]})")
print("TOTAL", total)
