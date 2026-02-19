#!/bin/bash
# Ralph Loop: repeatedly runs `claude -p` on a prompt file to iteratively
# simplify a source file. Each iteration is a fresh claude session that reads
# the prompt, makes edits, and exits. The script tracks line-count progress
# across iterations.
#
# Error handling:
#   - Rate limit:    sleep until the API reset time + 1min buffer, then retry
#   - Transient 500: retry with backoff (1min, 5min, 30min), give up after 3
#   - Stale session: if no new log output for 15min, kill and retry as transient
#
# Logs: each attempt writes stream-json to .ralph-logs/iter_NNx.log.jsonl
# (where x = a, b, c... for retries within an iteration) and stderr to
# .ralph-logs/iter_NNx.log.stderr. A one-line-per-iteration summary is
# appended to .ralph-logs/summary.log.
#
#
# Alternatively, just run a command and see interactive claude sessions in real time without the error handling and logging:
# `while :; do cat PROMPT.md | claude --dangerously-skip-permissions; done`
# or
# `for i in $(seq 1 10); do cat PROMPT.md | claude --dangerously-skip-permissions; done`

set -euo pipefail

# --- Configuration ---
MAX_ITER=${1:-10}
FILE="src/mechanism_enumeration.jl"
LOG_DIR=".ralph-logs"
PROMPT_FILE="PROMPT.md"
STALE_TIMEOUT=1800  # seconds (30 min) — kill claude if no new turns completed

mkdir -p "$LOG_DIR"

# --- Live output formatting ---
# Parses claude's stream-json output and prints a compact live summary:
#   - Tool calls: colored [Edit]/[Write]/[Bash]/[Task] with filename or command
#   - Assistant text: dimmed italic first line (truncated to 100 chars)
#   - Everything else: suppressed
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

    # --- Retry loop ---
    # On failure, classifies the error and either sleeps+retries or gives up.
    # The pipeline runs in a background subshell so a foreground watchdog can
    # monitor log freshness and kill stale sessions.
    # Each attempt gets its own log file (iter_02a, iter_02b, ...) so failed
    # attempts are preserved for debugging. On success, the final attempt's
    # suffix letter is kept as-is (no renaming needed).
    attempt=0
    while true; do
        # Build log filename with attempt suffix: iter_02a.log.jsonl, iter_02b.log.jsonl, ...
        attempt_letter=$(printf "\\$(printf '%03o' $((97 + attempt)))")
        cur_log="${log_file}${attempt_letter}"
        rm -f "$cur_log.exit"

        # Run claude in a background subshell.
        # Pipeline: claude streams JSON → tee saves to log → grep filters to
        # valid JSON lines → jq formats live output.
        # `set +eo pipefail` inside the subshell so we can read PIPESTATUS
        # (otherwise the first non-zero exit aborts the subshell).
        # PIPESTATUS[1] is claude's exit code (index 0 is cat, 1 is claude).
        (
            set +eo pipefail
            cat "$PROMPT_FILE" | claude -p --verbose --output-format stream-json --dangerously-skip-permissions \
                2>"$cur_log.stderr" | tee "$cur_log.jsonl" | \
                (grep --line-buffered '^{' || true) | \
                (jq --unbuffered -rj "$JQ_FILTER" || true)
            echo "${PIPESTATUS[1]}" > "$cur_log.exit"
        ) &
        pipeline_pid=$!

        # Stale session watchdog: checks every 30s whether new turns have
        # completed. Counts "assistant" and "user" events in the log — these
        # represent actual API round-trips. Thinking deltas and system events
        # are ignored since they don't indicate real progress (claude can emit
        # thinking tokens for a long time while effectively stuck).
        last_turns=0
        last_turn_change=$(date +%s)
        while kill -0 "$pipeline_pid" 2>/dev/null; do
            sleep 30
            cur_turns=$(grep -c '"type":"assistant"\|"type":"user"' "$cur_log.jsonl" 2>/dev/null || echo 0)
            now=$(date +%s)
            if [ "$cur_turns" -ne "$last_turns" ]; then
                last_turns=$cur_turns
                last_turn_change=$now
            elif [ $((now - last_turn_change)) -gt "$STALE_TIMEOUT" ]; then
                echo ""
                echo "[$(date '+%H:%M:%S')] No new turns for $((STALE_TIMEOUT / 60))min ($last_turns turns completed). Killing stale session..."
                pkill -P "$pipeline_pid" 2>/dev/null || true
                kill "$pipeline_pid" 2>/dev/null || true
                break
            fi
        done

        wait "$pipeline_pid" 2>/dev/null || true
        # If the subshell completed normally, .exit has claude's exit code.
        # If it was killed by the watchdog, .exit doesn't exist → default to 1.
        claude_exit=$(cat "$cur_log.exit" 2>/dev/null || echo 1)

        # Success → move on to the next iteration
        if [ "$claude_exit" -eq 0 ]; then
            server_retries=0
            break
        fi

        # --- Failure classification ---

        # Check for rate limiting. Claude's stream-json emits rate_limit_event
        # objects throughout normal operation with status="allowed". An actual
        # rate limit produces status != "allowed" (e.g., "limited"). We only
        # look at non-"allowed" events to avoid false positives (a 500 error
        # would still have "allowed" events from before the crash).
        resets_at=$(grep '"rate_limit_event"' "$cur_log.jsonl" 2>/dev/null \
            | jq -r 'select(.rate_limit_info.status != "allowed") | .rate_limit_info.resetsAt // empty' \
            2>/dev/null | tail -1 || true)

        # Rate limit → sleep until the API quota resets, plus 1min buffer
        if [ -n "$resets_at" ]; then
            now=$(date +%s)
            sleep_secs=$(( resets_at - now + 60 ))
            [ "$sleep_secs" -lt 60 ] && sleep_secs=60

            echo ""
            echo "[$(date '+%H:%M:%S')] Rate limited. Sleeping ${sleep_secs}s (~$((sleep_secs / 60))min) until $(date -r "$((resets_at + 60))" '+%H:%M:%S')..."
            sleep "$sleep_secs"
            echo "[$(date '+%H:%M:%S')] Resuming iteration $i"
            attempt=$((attempt + 1))
            continue
        fi

        # Not a rate limit → transient error (API 500, stale timeout, etc.)
        # Retry with increasing backoff. Give up after 3 attempts.
        server_retries=$(( ${server_retries:-0} + 1 ))
        case "$server_retries" in
            1) wait_secs=60 ;;
            2) wait_secs=300 ;;
            3) wait_secs=1800 ;;
            *)
                echo ""
                echo "[$(date '+%H:%M:%S')] claude exited with status $claude_exit after 3 retries. Stopping."
                exit 1
                ;;
        esac
        echo ""
        echo "[$(date '+%H:%M:%S')] claude exited with status $claude_exit. Retrying in $((wait_secs / 60))min (attempt $server_retries/3)..."
        sleep "$wait_secs"
        attempt=$((attempt + 1))
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
