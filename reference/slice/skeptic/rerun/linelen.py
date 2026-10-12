# ABOUTME: Lists lines longer than 92 characters (Unicode characters, not bytes) in the
# ABOUTME: given Julia files; usage: python3 linelen.py FILE...
import sys
for path in sys.argv[1:]:
    for i, line in enumerate(open(path, encoding="utf-8").read().split("\n"), 1):
        if len(line) > 92:
            print(f"{path}:{i}: {len(line)}")
