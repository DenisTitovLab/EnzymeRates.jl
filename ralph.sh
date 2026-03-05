#!/bin/bash
# Ralph Loop: repeatedly runs `claude -p` on a prompt file to sequentially
# implement a plan. Each iteration is a fresh claude session that reads
# the prompt, picks the next unfinished step, implements it, and exits.
#
# Error handling:
#   - Rate limit:    sleep until the API reset time + 1min buffer, then retry
#   - Transient 500: retry with backoff (1min, 5min, 30min), give up after 3
#                    Counter resets if claude made progress (>=5 completed turns)
#   - Stale session: if no new turns for 2hr, kill and retry as transient
#
# Logs: each attempt writes stream-json to .ralph-logs/iter_NNx.log.jsonl
# (where x = a, b, c... for retries within an iteration) and stderr to
# .ralph-logs/iter_NNx.log.stderr.
#
# Alternatively, just run a command and see interactive claude sessions in
# real time without the error handling and logging:
#   while :; do cat PROMPT.md | claude --dangerously-skip-permissions; done
# or:
#   for i in $(seq 1 10); do cat PROMPT.md | claude --dangerously-skip-permissions; done

# Exit on error, undefined vars, or pipe failures
set -euo pipefail

# --- Configuration ---
MAX_ITER=${1:-10}                # number of implementation iterations
LOG_DIR=".ralph-logs"
PROMPT_FILE="PLAN_IMPLEMENTATION_PROMPT.md"          # read fresh each attempt (edits take effect live)
FILES=(
    src/rate_eq_derivation.jl
    src/sym_poly_for_rate_eq_derivation.jl
    src/thermodynamic_constr_for_rate_eq_derivation.jl
)
STALE_TIMEOUT=7200               # seconds (2 hr) — kill claude if no new turns
PROGRESS_THRESHOLD=5             # turns needed to consider an attempt "productive"
LETTERS="abcdefghijklmnopqrstuvwxyz"  # attempt suffix lookup

mkdir -p "$LOG_DIR"
rm -f .ralph-done                    # clean sentinel from previous runs

# Count non-comment, non-docstring, non-blank lines in Julia files.
# Strips #-comments, tracks """...""" docstrings, skips blank lines.
count_code_lines() {
    awk '
    BEGIN { in_doc = 0 }
    {
        line = $0
        # Toggle docstring state on lines containing """
        while (match(line, /\x22\x22\x22/)) {
            in_doc = !in_doc
            line = substr(line, RSTART + 3)
        }
        if (in_doc) next
        # Strip inline comments (naive — ignores # inside strings)
        sub(/#.*/, "", $0)
        # Skip blank lines
        if ($0 ~ /^[[:space:]]*$/) next
        count++
    }
    END { print count+0 }
    ' "$@"
}

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
      ($first | .[0:300]) as $headline |
      "\u001b[90m\u001b[3m" + $headline + (if $len > 360 then " (" + ($len | tostring) + " chars)" else "" end) + "\u001b[0m\n"
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
initial_code=$(count_code_lines "${FILES[@]}")
echo "=== Ralph Loop ==="
echo "Files: ${FILES[*]}"
echo "Code lines (non-comment/doc/blank): $initial_code"
echo "Max iterations: $MAX_ITER"
echo "Logs: $LOG_DIR/"
echo ""

# Clean exit message on Ctrl+C
trap 'echo ""; echo "Interrupted at iteration $i. Check $LOG_DIR/ for logs."; exit 1' INT TERM

for i in $(seq 1 "$MAX_ITER"); do
    before=$(count_code_lines "${FILES[@]}")
    echo "[$(date '+%H:%M:%S')] Iteration $i/$MAX_ITER ($before code lines)"

    # Log file base: .ralph-logs/iter_02.log (suffixes added per attempt)
    log_file="$LOG_DIR/iter_$(printf '%02d' "$i").log"

    # --- Retry loop ---
    # On failure, classifies the error and either sleeps+retries or gives up.
    # The pipeline runs in a background subshell so a foreground watchdog can
    # monitor progress and kill stale sessions.
    # Each attempt gets its own log file (iter_02a, iter_02b, ...) so failed
    # attempts are preserved for debugging.
    attempt=0
    server_retries=0  # escalating backoff counter, resets on productive attempts
    while true; do
        # Log filename with attempt suffix: iter_02a.log.jsonl, iter_02b.log.jsonl
        cur_log="${log_file}${LETTERS:$attempt:1}"
        # Remove stale exit file — its absence after wait means the pipeline
        # was killed (by watchdog or signal) rather than exiting normally
        rm -f "$cur_log.exit"

        # Run claude in a background subshell so the watchdog can monitor it.
        # Pipeline stages:
        #   claude -p: runs in non-interactive piped mode, streams JSON to stdout
        #   tee: saves raw stream to log file AND passes it downstream
        #   grep '^{': filters to valid JSON lines (discards partial/non-JSON output)
        #   jq: formats selected events for live terminal display
        # `set +eo pipefail` inside subshell: lets pipeline finish so we can
        # read PIPESTATUS (otherwise first non-zero exit aborts the subshell).
        # `|| true` on grep/jq: prevents their exit codes from affecting PIPESTATUS.
        # PIPESTATUS[0] = claude's exit code (0=claude, 1=tee, 2=grep, 3=jq).
        (
            set +eo pipefail
            claude -p --verbose --output-format stream-json --dangerously-skip-permissions \
                < "$PROMPT_FILE" \
                2>"$cur_log.stderr" | tee "$cur_log.jsonl" | \
                (grep --line-buffered '^{' || true) | \
                (jq --unbuffered -rj "$JQ_FILTER" || true)
            echo "${PIPESTATUS[0]}" > "$cur_log.exit"
        ) &
        pipeline_pid=$!

        # --- Stale session watchdog ---
        # Checks every 30s whether new API round-trips have completed by counting
        # "assistant" and "user" events in the log. These represent actual turns.
        # Thinking deltas and system events are ignored — claude can emit thinking
        # tokens for a long time while effectively stuck in a loop.
        # `kill -0` checks if the subshell process is still alive (signal 0 = test only).
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
                # pkill -P: kill children of the subshell (claude, tee, grep, jq)
                # kill: kill the subshell itself (in case children already exited)
                pkill -P "$pipeline_pid" 2>/dev/null || true
                kill "$pipeline_pid" 2>/dev/null || true
                break
            fi
        done

        wait "$pipeline_pid" 2>/dev/null || true
        # If the subshell completed normally, .exit contains claude's exit code.
        # If it was killed by the watchdog, .exit doesn't exist → default to 1.
        claude_exit=$(cat "$cur_log.exit" 2>/dev/null || echo 1)

        # Success → move on to the next iteration
        if [ "$claude_exit" -eq 0 ]; then
            break
        fi

        # --- Failure classification ---

        # Check for rate limiting. Claude's stream-json emits rate_limit_event
        # objects throughout normal operation with status="allowed". An actual
        # rate limit produces status != "allowed" (e.g., "limited"). We filter
        # for non-"allowed" events to avoid false positives — a 500 error log
        # still contains "allowed" events from before the crash.
        resets_at=$(grep '"rate_limit_event"' "$cur_log.jsonl" 2>/dev/null \
            | jq -r 'select(.rate_limit_info.status != "allowed") | .rate_limit_info.resetsAt // empty' \
            2>/dev/null | tail -1 || true)

        # Rate limit → sleep until the API quota resets, plus 1min buffer
        if [ -n "$resets_at" ]; then
            now=$(date +%s)
            sleep_secs=$(( resets_at - now + 60 ))
            [ "$sleep_secs" -lt 60 ] && sleep_secs=60

            echo ""
            echo "[$(date '+%H:%M:%S')] Rate limited. Sleeping ${sleep_secs}s (~$((sleep_secs / 60))min) until $(date -d "@$((resets_at + 60))" '+%H:%M:%S')..."
            sleep "$sleep_secs"
            echo "[$(date '+%H:%M:%S')] Resuming iteration $i"
            attempt=$((attempt + 1))
            continue
        fi

        # Not a rate limit → transient error (API 500, stale timeout, etc.)
        # If claude completed enough turns before failing, it was making real
        # progress, so reset the backoff counter. This prevents a long productive
        # session that stalls near the end from burning through retry budget.
        # Only escalate backoff when claude fails repeatedly without useful work.
        completed_turns=$(grep -c '"type":"assistant"\|"type":"user"' "$cur_log.jsonl" 2>/dev/null || echo 0)
        if [ "$completed_turns" -ge "$PROGRESS_THRESHOLD" ]; then
            server_retries=0
        fi
        server_retries=$(( server_retries + 1 ))
        case "$server_retries" in
            1) wait_secs=60 ;;
            2) wait_secs=300 ;;
            3) wait_secs=1800 ;;
            *)
                echo ""
                echo "[$(date '+%H:%M:%S')] claude exited with status $claude_exit after 3 retries without progress. Stopping."
                exit 1
                ;;
        esac
        echo ""
        echo "[$(date '+%H:%M:%S')] claude exited with status $claude_exit ($completed_turns turns completed). Retrying in $((wait_secs / 60))min (attempt $server_retries/3)..."
        sleep "$wait_secs"
        attempt=$((attempt + 1))
    done

    # --- Iteration summary ---
    after=$(count_code_lines "${FILES[@]}")
    delta=$((before - after))
    echo ""
    echo "[$(date '+%H:%M:%S')] Iteration $i done: $before → $after code lines (Δ$delta)"
    echo "---"

    # --- Plan completion check ---
    # The agent creates .ralph-done when all plan steps are finished
    if [ -f ".ralph-done" ]; then
        echo ""
        echo "=== Plan fully implemented ==="
        cat .ralph-done
        break
    fi
done

final_code=$(count_code_lines "${FILES[@]}")
echo ""
echo "=== Complete ($i iterations): $initial_code → $final_code code lines (Δ$((initial_code - final_code))) ==="
