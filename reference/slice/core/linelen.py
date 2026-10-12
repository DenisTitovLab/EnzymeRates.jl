# ABOUTME: Lists lines longer than 92 characters (not bytes) in the given files.
# ABOUTME: Used to keep the slice sources within the package's line-length limit.
import sys
for p in sys.argv[1:]:
    for i, l in enumerate(open(p, encoding="utf-8").read().split("\n"), 1):
        if len(l) > 92:
            print(f"{p}:{i}: {len(l)}")
