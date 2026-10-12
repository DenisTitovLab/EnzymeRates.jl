# ABOUTME: Counts Julia code lines (non-blank, non-comment, outside docstrings) of files or line ranges.
# ABOUTME: Usage: python3 loc.py file[:start-end] ...
import sys
def code_lines(lines):
    n = 0; in_doc = False
    for ln in lines:
        s = ln.strip()
        if in_doc:
            if s.endswith('"""') : in_doc = False
            continue
        if s.startswith('"""'):
            if not (len(s) > 3 and s.endswith('"""')):
                in_doc = True
            continue
        if not s or s.startswith('#'):
            continue
        n += 1
    return n
tot = 0
for arg in sys.argv[1:]:
    if ':' in arg.split('/')[-1]:
        f, r = arg.rsplit(':', 1); a, b = map(int, r.split('-'))
        lines = open(f).read().split('\n')[a-1:b]
    else:
        f = arg; lines = open(f).read().split('\n')
    c = code_lines(lines); tot += c
    print(f"{c:5d}  {arg}")
print(f"{tot:5d}  total")
