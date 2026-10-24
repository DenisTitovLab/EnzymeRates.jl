# ABOUTME: Splits a scan_all.jl log at its variant markers and runs hoistscan.py's scan on each
# ABOUTME: segment, printing one line per variant (escape analyses, hoists, dangerous hoists).
import sys, os, tempfile, subprocess
HS = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/hoistscan.py"
for path in sys.argv[1:]:
    lines = open(path, errors="replace").read().split("\n")
    segs, cur, name = [], [], "preamble"
    for l in lines:
        if l.startswith("=== VARIANT "):
            segs.append((name, cur)); cur, name = [], l.split()[-1]
            continue
        if l.startswith("=== END VARIANT"):
            continue
        cur.append(l)
    segs.append((name, cur))
    with tempfile.TemporaryDirectory() as d:
        for i, (n, seg) in enumerate(segs):
            if n == "preamble" and not any("remark:" in x for x in seg):
                continue
            p = os.path.join(d, f"{os.path.basename(path)}.{i}.{n}")
            open(p, "w").write("\n".join(seg))
            out = subprocess.run([sys.executable, "-I", HS, p], capture_output=True, text=True).stdout
            rounds = [x for x in seg if x.startswith(("derivation round", "identify round"))]
            print(f"{os.path.basename(path)} [{n}]: " + " ".join(out.split("\n")[1].split()) +
                  f" | rounds {len(rounds)} | GC error lines {sum('GC error' in x for x in seg)}")
            for x in out.split("\n")[2:]:
                if x.strip(): print("    " + x.strip())
