#!/bin/bash
set -euo pipefail

# --- Configuration ---
MAX_ITER=${1:-10}
FILE="src/mechanism_enumeration.jl"
LOG_DIR=".ralph-logs"
PROMPT_FILE="PROMPT.md"

mkdir -p "$LOG_DIR"

# --- Live output formatting ---
# jq filter: shows Edit/Write/Bash/Task tool calls in color,
# assistant text as dimmed first-line summary
read -r -d '' JQ_FILTER << 'JQEOF' || true
if .type == "assistant" and .message.content != null then
  .message.content[] |
    if .type == "text" then
      (.text | length) as $len |
      (.text | split("\n") | map(select(length > 0)) | first // "") as $first |
      ($first | .[0:100]) as $headline |
      "\u001b[90m\u001b[3m" + $headline + (if $len > 120 then " (" + ($len | tostring) + " chars)" else "" end) + "\u001b[0m\n"
    elif .type == "tool_use" then
      if .name == "Edit" then
        "\u001b[33m[Edit]\u001b[0m " + (.input.file_path // "" | split("/") | last) + "\n"
      elif .name == "Write" then
        "\u001b[33m[Write]\u001b[0m " + (.input.file_path // "" | split("/") | last) + "\n"
      elif .name == "Bash" then
        "\u001b[36m[Bash]\u001b[0m " + (.input.command // "" | .[0:80]) + "\n"
      elif .name == "Task" then
        "\u001b[35m[Task]\u001b[0m " + (.input.description // "") + "\n"
      else empty end
    else empty end
else empty end
JQEOF

# --- Main loop ---
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

    # Retry loop: if claude hits a rate limit, sleep until reset and retry
    while true; do
        # Disable exit-on-error so we can inspect PIPESTATUS
        set +eo pipefail
        cat "$PROMPT_FILE" | claude -p --verbose --output-format stream-json --dangerously-skip-permissions \
            2>"$log_file.stderr" | tee "$log_file.jsonl" | \
            (grep --line-buffered '^{' || true) | \
            (jq --unbuffered -rj "$JQ_FILTER" || true)
        claude_exit=${PIPESTATUS[1]}  # exit status of `claude` (index 1 in the pipeline)
        set -eo pipefail

        [ "$claude_exit" -eq 0 ] && break

        # Claude failed — check if it was a rate limit by looking for the
        # last rate_limit_event in the stream-json output
        resets_at=$(grep '"rate_limit_event"' "$log_file.jsonl" 2>/dev/null \
            | tail -1 | jq -r '.rate_limit_info.resetsAt // empty' 2>/dev/null || true)

        # No rate limit info → unknown error, bail out
        if [ -z "$resets_at" ]; then
            echo ""
            echo "[$(date '+%H:%M:%S')] claude exited with status $claude_exit. Stopping."
            exit 1
        fi

        # Sleep until resetsAt + 60s buffer, minimum 60s
        now=$(date +%s)
        sleep_secs=$(( resets_at - now + 60 ))
        [ "$sleep_secs" -lt 60 ] && sleep_secs=60

        echo ""
        echo "[$(date '+%H:%M:%S')] Rate limited. Sleeping ${sleep_secs}s (~$((sleep_secs / 60))min) until $(date -r "$((resets_at + 60))" '+%H:%M:%S')..."
        sleep "$sleep_secs"
        echo "[$(date '+%H:%M:%S')] Resuming iteration $i"
    done

    # --- Iteration summary ---
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
