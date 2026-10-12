#!/bin/bash
# ABOUTME: Counts code lines (non-blank, non-comment, outside docstrings) of a file or a line range.
# ABOUTME: Usage: codelines.sh file [first last]
f=$1; a=${2:-1}; b=${3:-999999}
awk -v a="$a" -v b="$b" '
NR < a || NR > b { next }
{
  line = $0; t = line; gsub(/^[ \t]+|[ \t]+$/, "", t)
  n = gsub(/"""/, "\"\"\"", line)
  if (indoc) { if (n % 2 == 1) indoc = 0; next }
  if (n % 2 == 1) { indoc = 1; next }
  if (n >= 2 && t ~ /^"""/) next
  if (t == "" || t ~ /^#/) next
  c++
}
END { print c + 0 }' "$f"
