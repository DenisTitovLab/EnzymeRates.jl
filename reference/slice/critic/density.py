# ABOUTME: Code-line density (non-blank characters per code line, statements per line) of Julia
# ABOUTME: files, using the codelines.py rule for what counts as a code line; usage: density.py FILES
import sys
def code_lines(path):
    out, indoc = [], False
    for raw in open(path, encoding="utf-8").read().split("\n"):
        s = raw.strip()
        q = s.count('"""')
        if indoc:
            if q % 2 == 1: indoc = False
            continue
        if s.startswith('"""'):
            if q % 2 == 1: indoc = True
            continue
        if q % 2 == 1: indoc = True
        if s and not s.startswith("#"):
            out.append(s.split(" #")[0].rstrip())
    return out
tot_l = tot_c = 0
for f in sys.argv[1:]:
    ls = code_lines(f)
    c = sum(len(l.replace(" ", "")) for l in ls)
    tot_l += len(ls); tot_c += c
    print(f"{f.split('/')[-1]:55s} lines {len(ls):5d} chars/line {c/max(len(ls),1):5.1f}")
print(f"{'TOTAL':55s} lines {tot_l:5d} chars/line {tot_c/max(tot_l,1):5.1f}")
