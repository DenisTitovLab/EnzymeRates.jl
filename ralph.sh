#!/bin/bash
set -euo pipefail

MAX_ITER=${1:-10}
FILE="src/mechanism_enumeration.jl"
LOG_DIR=".ralph-logs"
PROMPT_FILE="PROMPT.md"

mkdir -p "$LOG_DIR"

initial_lines=$(wc -l < "$FILE" | tr -d ' ')
echo "=== Ralph Loop ==="
echo "File: $FILE ($initial_lines lines)"
echo "Max iterations: $MAX_ITER"
echo "Logs: $LOG_DIR/"
echo ""

trap 'echo ""; echo "Interrupted at iteration $i. Check $LOG_DIR/summary.log for progress."; exit 1' INT TERM

for i in $(seq 1 "$MAX_ITER"); do
    before=$(wc -l < "$FILE" | tr -d ' ')
    timestamp=$(date '+%H:%M:%S')
    log_file="$LOG_DIR/iter_$(printf '%02d' "$i").log"

    echo "[$timestamp] Iteration $i/$MAX_ITER — $FILE is $before lines"

    cat "$PROMPT_FILE" | claude --print --dangerously-skip-permissions \
        2>&1 | tee "$log_file"

    after=$(wc -l < "$FILE" | tr -d ' ')
    delta=$((before - after))
    timestamp_done=$(date '+%H:%M:%S')

    echo ""
    echo "[$timestamp_done] Iteration $i done: $before → $after lines (Δ$delta)"
    echo "---"
    echo "[$timestamp_done] Iter $i: $before → $after (Δ$delta)" >> "$LOG_DIR/summary.log"
done

final=$(wc -l < "$FILE" | tr -d ' ')
echo ""
echo "=== Complete ==="
echo "Result: $initial_lines → $final lines ($(($initial_lines - $final)) total removed)"
