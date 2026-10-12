# ABOUTME: Skeptic P5a: counts Julia code lines (no blanks, comments or docstrings) for whole files
# ABOUTME: or for named top-level functions, to check the probe's code-size claims.
import re, sys
def code_lines(text):
    out, in_doc, in_block = [], False, 0
    for line in text.split("\n"):
        s = line.strip()
        if in_block:
            if "=#" in s: in_block -= 1
            continue
        if s.startswith("#="):
            if "=#" not in s: in_block += 1
            continue
        if in_doc:
            if s.endswith('"""') or s == '"""': in_doc = False
            continue
        if s.startswith('"""'):
            if not (len(s) > 3 and s.endswith('"""')): in_doc = True
            continue
        if not s or s.startswith("#"): continue
        out.append(line)
    return out
def functions(text):
    """top-level definitions: name -> source block (from its line to the matching 'end' or the one-liner)"""
    lines = text.split("\n"); res = {}; i = 0
    while i < len(lines):
        l = lines[i]
        m = re.match(r'^(?:@generated\s+)?function\s+([A-Za-z_!0-9]+)', l) or re.match(r'^([A-Za-z_!0-9]+)\(.*\)\s*=', l)
        if m:
            name = m.group(1); j = i
            if l.startswith("function") or l.startswith("@generated"):
                while j < len(lines) and lines[j] != "end": j += 1
            else:
                while j + 1 < len(lines) and lines[j + 1].startswith(" "): j += 1
            res.setdefault(name, []).append("\n".join(lines[i:j + 1])); i = j + 1
        else: i += 1
    return res
if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "files":
        tot = 0
        for f in sys.argv[2:]:
            n = len(code_lines(open(f).read())); tot += n; print(f"{n:6d}  {f.split('/')[-1]}")
        print(f"{tot:6d}  total")
    else:
        f = sys.argv[2]; fs = functions(open(f).read()); tot = 0
        for name in sys.argv[3:]:
            blocks = fs.get(name, [])
            n = sum(len(code_lines(b)) for b in blocks); tot += n
            print(f"{n:6d}  {name} ({len(blocks)} method(s))")
        print(f"{tot:6d}  total")
