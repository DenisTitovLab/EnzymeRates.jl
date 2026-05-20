#!/usr/bin/env bash
# ABOUTME: Per-test-file runtime report for stage closeouts.
# ABOUTME: Records each test file's @elapsed include time; if a
# ABOUTME: baseline exists, flags files whose runtime grew >= 2x.
# Informational ONLY — does not fail CI.
set -uo pipefail

BASELINE_FILE="${1:-docs/superpowers/refactor-test-timings-main-baseline.txt}"
OUT_FILE="${2:-/tmp/stage-test-timings.txt}"

echo "# Per-test-file runtime report (informational)"  > "$OUT_FILE"
echo "# Format: <test_file> <elapsed_seconds>"        >> "$OUT_FILE"

# Time each test file in its own subprocess for clean isolation.
# Skip runtests.jl (wrapper) and any helper files that don't
# contain @testset blocks themselves.
for f in test/test_*.jl; do
    base=$(basename "$f")
    # Wrap in a try so a failing test doesn't abort the loop —
    # the test-integrity gate handles failure separately.
    # All test stdout/stderr goes to /dev/null; the timing is written
    # to a side-channel file we read back. Otherwise test output (println,
    # CMA-ES iteration prints, Test Summary tables) pollutes the report.
    timing_file=$(mktemp)
    julia --project -e "
        using Pkg
        Pkg.activate(\"test\")
        try
            t = @elapsed include(\"$f\")
            open(\"$timing_file\", \"w\") do io
                print(io, round(t; digits=2))
            end
        catch e
            open(\"$timing_file\", \"w\") do io
                print(io, \"FAIL\")
            end
        end
    " >/dev/null 2>&1
    t=$(cat "$timing_file")
    rm -f "$timing_file"
    echo "$base $t" >> "$OUT_FILE"
done

echo ""
echo "=== Test-runtime report ==="
cat "$OUT_FILE"

if [ -f "$BASELINE_FILE" ]; then
    echo ""
    echo "=== Regression check vs $BASELINE_FILE ==="
    echo "    (2x+ flagged for INVESTIGATE; informational only, not failing)"
    awk '
        NR==FNR {
            b[$1] = $2
            next
        }
        $1 in b && $2 != "FAIL" && b[$1] != "FAIL" {
            r = $2 / b[$1]
            tag = (r >= 2.0 ? " <-- INVESTIGATE (>=2x baseline)" : "")
            printf "  %-50s %.2fx  (%.2fs vs %.2fs baseline)%s\n", $1, r, $2, b[$1], tag
        }
        $1 in b && ($2 == "FAIL" || b[$1] == "FAIL") {
            printf "  %-50s   FAIL state changed (was %s, now %s)\n", $1, b[$1], $2
        }
    ' "$BASELINE_FILE" "$OUT_FILE"
else
    echo ""
    echo "(no baseline file at $BASELINE_FILE — first run is establishing the baseline)"
fi
