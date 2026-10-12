# ABOUTME: Counts code lines of Julia files: non-blank, not a comment, not inside a docstring.
# ABOUTME: Usage: python3 codelines.py file.jl ...
import sys
def count(path):
    n, indoc = 0, False
    for line in open(path):
        s = line.strip()
        if indoc:
            if s.endswith('"""'):
                indoc = False
            continue
        if s.startswith('"""'):
            if not (len(s) > 3 and s.endswith('"""')):
                indoc = True
            continue
        if not s or s.startswith('#'):
            continue
        n += 1
    return n
tot = 0
for p in sys.argv[1:]:
    c = count(p); tot += c
    print(f"{c:5d}  {p}")
print(f"{tot:5d}  total")
