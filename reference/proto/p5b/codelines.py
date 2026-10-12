# ABOUTME: Counts Julia code lines (no blanks, comments or docstrings) in a file or a line range.
# ABOUTME: Usage: python3 codelines.py FILE [FIRST LAST]
import sys
path = sys.argv[1]
lines = open(path).read().split("\n")
lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (1, len(lines))
indoc, n = False, 0
for i, raw in enumerate(lines, 1):
    s = raw.strip()
    q = s.count('"""')
    if indoc:
        if q % 2 == 1: indoc = False
        continue
    if s.startswith('"""'):
        if q % 2 == 1: indoc = True
        continue
    if q % 2 == 1:          # a code line that opens a docstring-like block
        indoc = True
    if lo <= i <= hi and s and not s.startswith("#"):
        n += 1
print(n)
