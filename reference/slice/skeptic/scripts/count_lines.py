# ABOUTME: Independent code-line counter for Julia files: skips blank lines, full-line comments,
# ABOUTME: #= =# block comments and docstrings (a triple-quoted string at statement start).
import sys, re

def count(path, verbose=False):
    lines = open(path, encoding="utf-8").read().split("\n")
    n, i, kept = 0, 0, []
    in_block = 0
    while i < len(lines):
        s = lines[i].strip()
        if in_block:
            in_block += s.count("#=") - s.count("=#")
            i += 1
            continue
        if s.startswith("#="):
            in_block = s.count("#=") - s.count("=#")
            i += 1
            continue
        if not s or s.startswith("#"):
            i += 1
            continue
        if s.startswith('"""'):
            # docstring: skip to the line holding its closing quotes
            rest = s[3:]
            if '"""' in rest:
                i += 1
                continue
            i += 1
            while i < len(lines) and '"""' not in lines[i]:
                i += 1
            i += 1
            continue
        n += 1
        kept.append(i + 1)
        i += 1
    return n, kept

if __name__ == "__main__":
    tot = 0
    for p in sys.argv[1:]:
        c, kept = count(p)
        tot += c
        print(f"{p.split('/')[-1]:16s} {c:5d}")
    print(f"{'total':16s} {tot:5d}")
