# Concrete-Types-Instead-of-Symbols Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace EnzymeRates' Symbol-string-juggling internals with one concrete struct family shared between enumeration and derivation, cutting src LOC ≥ 50% while preserving the `rate_equation` 0-alloc/<100ns invariant.

**Architecture:** Two struct flavors share field types. `Mechanism`/`AllostericMechanism` are non-parametric (one Julia type each) — what enumeration produces, cheap to instantiate millions of times. `EnzymeMechanism{Sig}`/`AllostericEnzymeMechanism{Sig}` are parametric — each unique `Sig` triggers ONE `@generated rate_equation` body-build. Conversion happens at the fitting boundary. Inside `@generated` bodies, a plain function `_mechanism_from_sig(Sig::Tuple) -> Mechanism` reconstructs concrete structs for body-build code to walk (no nested `@generated`).

**Tech Stack:** Julia 1.10+, Optimization.jl, Tables.jl, Test (stdlib).

**Spec:** `docs/superpowers/specs/2026-05-20-concrete-types-refactor-design.md`

---

## Operating Rules (apply to every task)

1. **TDD per CLAUDE.md.** Write failing test first; run to confirm failure; implement minimum; run to confirm pass; refactor if cleanup obvious.
2. **No test deletion.** See spec §2. Mechanical syntax adaptation only; assertion strength unchanged. No `@test_skip`, no `@test_broken`, no commenting out.
3. **`rate_equation` perf invariant.** `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl` must stay green at every commit. `allocs == 0`, `t < 100e-9` per call for every mechanism in `MECHANISM_TEST_SPECS`.
4. **Per-commit src LOC delta.** Every commit message ends with one line: `src delta: -X / +Y net Z, cumulative: -W`. Computed via `wc -l src/*.jl` before/after.
5. **No Symbol-string parsing.** After Stage 3, `Symbol("X$idx")`, `startswith(string(sym), …)`, `_form_name`, `_parse_bound`, `_rename_params_T` and friends should not exist in src outside the `name(p::Parameter, m)` chokepoint accessor.
6. **Mid-refactor checkpoint.** After Stage 4 (Task 4.6), STOP and review. Cumulative src delta must be ≤ −500; if not, revise spec before continuing.
7. **Commit cadence.** One git commit per task (not per step). Steps within a task are commits-in-progress.

---

# Stage 0 — Compile-budget + test-integrity gates (separate PR to main)

**Goal:** Land minimal CI gates first so the refactor branch can rebase against them. Two compile-budget gates (one trace-compile for `init_mechanisms`, one wall-clock for `rate_equation` body-build) + one test-integrity check script. Calibrate against current main with 2× headroom; recalibrate if Stage 3 trips the gate.

**Branch:** `add-compile-budget-tests` off `main` (new branch).

**Expected LOC delta:** +0 src, +~150 test (test + script infra).

### Task 0.1: Create the compile-budget test

**Files:**
- Create: `test/test_compile_budget.jl`

- [ ] **Step 1: Write the test file**

Create `test/test_compile_budget.jl`:

```julia
# ABOUTME: Compile-time regression gates for the EnzymeRates pipeline.
# ABOUTME: One trace-compile (init_mechanisms) + one wall-clock (rate_equation body-build).

using Test
using EnzymeRates

# Budgets calibrated against current main with 2× headroom (Task 0.3).
# Recalibrate if Stage 3 or 4 trips the gate.
const INIT_TRACE_BUDGET                  = 200   # calibrate in Task 0.3
const RATE_EQUATION_WALLCLOCK_BUDGET_S   = 5.0   # calibrate in Task 0.3

# Anchored to the EnzymeRates module prefix only. Counts every method
# specialization Julia compiles that touches our module — our functions,
# our types, Base methods specialized on our types (e.g.
# Base.hash(::EnzymeRates.Step, ...)), Core.kwcall plumbing, show/print.
# Filters out unrelated package precompiles (Optimization.jl, Tables.jl,
# etc.) which would create dependency-noise false positives.
#
# Intentionally NOT a suffix enumeration — that approach silently
# misses renamed or newly-introduced internal helpers during the
# refactor. The module-prefix anchor catches everything in our
# namespace automatically.
const RELEVANT_PRECOMPILE_PATTERN = r"EnzymeRates\."

function _count_relevant_precompiles(runner_script::String)
    trace_file = tempname()
    julia_exe = Base.julia_cmd().exec[1]
    cmd = Cmd([julia_exe, "--trace-compile=$(trace_file)",
               "--project=.", "-e", runner_script])
    try
        run(cmd; wait=true)
    catch e
        @warn "Subprocess failed: $e"
        return -1
    end
    isfile(trace_file) || (@info "trace-compile file missing"; return -1)
    n = try
        length(filter(line -> occursin(RELEVANT_PRECOMPILE_PATTERN, line) &&
                              !isempty(strip(line)),
                      collect(eachline(trace_file))))
    finally
        rm(trace_file, force=true)
    end
    n
end

@testset "compile-budget" begin
    # Trace-compile: init_mechanisms on a bi-bi reaction (non-trivial so
    # the gate is representative; uni-uni is too small to catch regressions).
    @testset "trace-compile: init_mechanisms (bi-bi)" begin
        script = """
            using EnzymeRates
            r = @enzyme_reaction begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
            end
            EnzymeRates.init_mechanisms(r)
            """
        n = _count_relevant_precompiles(script)
        @info "init_mechanisms trace-compile: $n (budget: $INIT_TRACE_BUDGET)"
        @test 0 <= n <= INIT_TRACE_BUDGET
    end

    # Wall-clock: rate_equation body-build (first call pays @generated cost).
    # The per-call runtime gate is separately enforced by
    # test_rate_eq_derivation.jl\'s test_rate_equation_performance
    # (0 allocs, <100ns per call) for every mechanism in MECHANISM_TEST_SPECS.
    @testset "wall-clock: rate_equation body-build (first call)" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        params = NamedTuple{Tuple(EnzymeRates.parameters(m))}(
            ntuple(_ -> 1.0, length(EnzymeRates.parameters(m))))
        concs = (S = 1.0, P = 0.5)
        t_first = @elapsed EnzymeRates.rate_equation(m, concs, params)
        @info "rate_equation first-call wall-clock: $(t_first)s " *
              "(budget: $RATE_EQUATION_WALLCLOCK_BUDGET_S s)"
        @test t_first < RATE_EQUATION_WALLCLOCK_BUDGET_S
    end

    # Warmup-reuse regression: on main today, calling init_mechanisms on
    # uni-uni FIRST makes the subsequent ter-ter init dramatically cheaper
    # — most enumeration specializations are shared across arity. If a
    # refactor commit accidentally introduces per-arity specialization
    # (e.g., a `where {N}` clause that triggers fresh @generated bodies
    # per substrate count), the ter-ter incremental cost approaches its
    # cold cost.
    #
    # Measurement: 3 subprocesses isolate the TER-TER INCREMENTAL cost
    # in both scenarios:
    #   n_uni_only:  using + uni-uni init                    (warmup-base)
    #   n_cold_ter:  using + ter-ter init                    (cold ter-ter)
    #   n_warm_ter:  using + uni-uni init + ter-ter init     (warm ter-ter)
    #
    # ter_cold_incremental = n_cold_ter - n_using_only      (~ cost of ter-ter alone)
    # ter_warm_incremental = n_warm_ter - n_uni_only        (what ter-ter ADDED after uni warmup)
    #
    # Gate: ter_warm_incremental << ter_cold_incremental (much less work
    # done because uni-uni already triggered the shared specializations).
    @testset "warmup-reuse: ter-ter incremental cost drops after uni-uni warmup" begin
        ter_block = """
            r_ter = @enzyme_reaction begin
                substrates: A[C], B[N], C[O]
                products:   P[C], Q[N], R[O]
            end
            EnzymeRates.init_mechanisms(r_ter)
            """
        uni_block = """
            r_uni = @enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
            end
            EnzymeRates.init_mechanisms(r_uni)
            """

        n_using_only  = _count_relevant_precompiles("using EnzymeRates")
        n_uni_only    = _count_relevant_precompiles("using EnzymeRates\n$uni_block")
        n_cold_ter    = _count_relevant_precompiles("using EnzymeRates\n$ter_block")
        n_warm_ter    = _count_relevant_precompiles(
            "using EnzymeRates\n$uni_block\n$ter_block")

        ter_cold_incr = n_cold_ter - n_using_only
        ter_warm_incr = n_warm_ter - n_uni_only

        @info "warmup-reuse: ter_cold_incremental=$ter_cold_incr, " *
              "ter_warm_incremental=$ter_warm_incr, " *
              "ratio=$(ter_cold_incr > 0 ? round(ter_warm_incr/ter_cold_incr; digits=3) : NaN)"

        # ter-ter alone should compile substantial new specializations.
        @test ter_cold_incr > 10

        # After uni-uni warmup, ter-ter should add MUCH less than its
        # cold cost (perfect sharing → near 0; tuned conservatively
        # to catch regressions where per-arity specialization snuck in).
        @test ter_warm_incr < ter_cold_incr ÷ 4
    end
end
```

- [ ] **Step 2: Verify it parses**

`julia --project -e \'include("test/test_compile_budget.jl")\''
Expected: tests run; the assertions may FAIL because budgets aren\'t calibrated yet — calibrate in Task 0.3.

### Task 0.2: Wire test_compile_budget.jl into runtests.jl

**Files:**
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add the include**

Add `include("test_compile_budget.jl")` in topological order with the other includes in `test/runtests.jl`.

- [ ] **Step 2: Run full suite**

`julia --project -e \'using Pkg; Pkg.test()\''
Expected: all existing tests pass; new compile-budget tests run (may FAIL until Task 0.3 calibrates).

### Task 0.3: Calibrate budgets against current main

**Files:**
- Modify: `test/test_compile_budget.jl`

- [ ] **Step 1: Measure on main**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates

# Trace-compile baseline
julia --project --trace-compile=/tmp/trace_init.log -e \'
    using EnzymeRates
    r = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end
    EnzymeRates.init_mechanisms(r)
\'
n=$(grep -E "EnzymeRates\." /tmp/trace_init.log | wc -l)
echo "init_mechanisms trace-compile baseline: $n"

# Wall-clock baseline
julia --project -e \'
    using EnzymeRates
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    params = NamedTuple{Tuple(EnzymeRates.parameters(m))}(
        ntuple(_ -> 1.0, length(EnzymeRates.parameters(m))))
    concs = (S = 1.0, P = 0.5)
    EnzymeRates.rate_equation(m, concs, params)   # warmup
    t = @elapsed EnzymeRates.rate_equation(m, concs, params)
    println("rate_equation wall-clock baseline: $t s")
\'
```

- [ ] **Step 2: Set budgets at ceil(2 × baseline)**

Update the constants in `test/test_compile_budget.jl`:

```julia
const INIT_TRACE_BUDGET                = ceil(Int, 2 * <baseline>)
const RATE_EQUATION_WALLCLOCK_BUDGET_S = 2 * <baseline_wc>
```

**NO RECALIBRATION. The budget is FIXED for the duration of the refactor.** If Stage 3's `@generated rate_equation` body-build rewrite (or any other stage) trips the gate, this is the gate doing its job — surfacing a compile-time regression for review. The response is to **REDESIGN the offending change to fit within budget**, NOT to raise the budget to accommodate it. Recalibration on trip would normalize regressions and defeat the gate.

Concrete example: if Stage 3's `_mechanism_from_sig` materialization pushes body-build trace-compile beyond 2×main, the executor must investigate WHY (recursive struct lifts? excessive `_to_sig`/`_from_sig` specializations? Mechanism construction overhead?) and FIX it — e.g., by walking `Sig` directly in the @generated body without materializing a concrete Mechanism, or by simplifying the Sig encoding to avoid per-Step Species sub-tuples. The gate's existence is what FORCES this redesign rather than letting the cost silently accumulate.

- [ ] **Step 3: Run with calibrated values**

`julia --project -e \'using Pkg; Pkg.test()\''
Expected: PASS.

### Task 0.4: Test-integrity check script

**Files:**
- Create: `scripts/check_test_integrity.sh`

**Rationale:** Per spec §2, no test may be deleted, commented out, or weakened. The script enforces the floor (catches deletion/skip/comment-out) AND emits WARN for any modified `@test`/`@test_throws` line for human review (per Reviewer C #4 — programmatic assertion weakening is hard to detect statically, so we surface it for eyeballs).

- [ ] **Step 1: Create the script**

```bash
mkdir -p scripts
cat > scripts/check_test_integrity.sh <<\'SCRIPT\'
#!/usr/bin/env bash
# ABOUTME: Test-integrity gate enforced at every stage closeout.
# ABOUTME: Hard checks: no file deletion, no @testset count drop,
# ABOUTME: no @test_skip/@test_broken, no commented-out @test lines.
# ABOUTME: Soft check: WARNs about any modified @test/@test_throws line
# ABOUTME: for human review (weakened assertions can\'t be caught statically).
# Per spec §2 (NON-NEGOTIABLE). Run from repo root.
set -uo pipefail

BASE_REF="${1:-main}"
fail=0

# Hard Check 1: no test file deleted (renames OK; deletions forbidden).
deleted=$(git diff --name-status "$BASE_REF"..HEAD -- test/ 2>/dev/null \
    | awk \'$1 == "D" { print $2 }\')
if [ -n "$deleted" ]; then
    echo "FAIL [Check 1]: test file(s) deleted vs $BASE_REF:"
    echo "$deleted" | sed \'s/^/    /\'
    fail=1
fi

# Hard Check 2: @testset count never decreases — EXCEPT for testsets
# documented in docs/superpowers/refactor-deleted-tests.md per spec §2.1
# narrow exception (tests of deleted helpers; conditions are strict).
n_base=0
for f in $(git ls-tree -r "$BASE_REF" --name-only test/ 2>/dev/null | grep "\.jl$"); do
    c=$(git show "$BASE_REF":"$f" 2>/dev/null | grep -c "@testset" || true)
    n_base=$((n_base + c))
done
n_head=$(grep -rh "@testset" test/ 2>/dev/null | grep -c "@testset" || true)
# Each "### test_..." heading in the deleted-tests log = one permitted deletion.
n_documented=$(grep -c "^### test_" docs/superpowers/refactor-deleted-tests.md 2>/dev/null || echo 0)
n_base_adj=$((n_base - n_documented))
echo "@testset count: $BASE_REF=$n_base, documented deletions=$n_documented, adjusted base=$n_base_adj, HEAD=$n_head"
if [ "$n_head" -lt "$n_base_adj" ]; then
    echo "FAIL [Check 2]: @testset count decreased beyond documented deletions"
    echo "    Add an entry to docs/superpowers/refactor-deleted-tests.md per spec §2.1,"
    echo "    OR restore the deleted testset(s) and re-apply as mechanical adaptation."
    fail=1
fi

# Hard Check 3: no @test_skip or @test_broken added.
forbidden=$(git diff "$BASE_REF"..HEAD -- test/ | grep "^+" | grep -v "^+++" \
    | grep -E "@test_skip|@test_broken" || true)
if [ -n "$forbidden" ]; then
    echo "FAIL [Check 3]: @test_skip or @test_broken added vs $BASE_REF:"
    echo "$forbidden" | sed \'s/^/    /\'
    fail=1
fi

# Hard Check 4: no @test lines commented out.
commented=$(git diff "$BASE_REF"..HEAD -- test/ | grep "^+" | grep -v "^+++" \
    | grep -E "^\+\s*#+\s*@test\b" || true)
if [ -n "$commented" ]; then
    echo "FAIL [Check 4]: @test line(s) commented out vs $BASE_REF:"
    echo "$commented" | sed \'s/^/    /\'
    fail=1
fi

# Soft Check 5: WARN on any modified @test/@test_throws line.
# Spec §2 forbids changing hardcoded values / weakening assertions.
# Static detection is unreliable — we emit a WARN list for human review.
modified_tests=$(git diff "$BASE_REF"..HEAD -U0 -- test/ \
    | grep -E "^[+-]\s*@test(_throws)?\b" \
    | grep -v "^---" | grep -v "^+++" \
    | grep -v "^+\s*@test\(_throws\)?\s*[^\s]" || true)
if [ -n "$modified_tests" ]; then
    echo ""
    echo "WARN [Check 5]: @test/@test_throws lines modified vs $BASE_REF."
    echo "  Human review REQUIRED — spec §2 forbids weakening assertions"
    echo "  (e.g., \'== 3\' → \'isa Int\', tolerance relaxation, swapped operators):"
    echo "$modified_tests" | sed \'s/^/    /\'
    echo "  If any change is a weakening, REVERT and re-apply as mechanical only."
fi

if [ "$fail" -eq 0 ]; then
    echo "Test-integrity hard checks PASSED vs $BASE_REF."
fi
exit "$fail"
SCRIPT
chmod +x scripts/check_test_integrity.sh
```

- [ ] **Step 2: Run it (no-op against fresh main)**

`bash scripts/check_test_integrity.sh main`
Expected: PASS.

- [ ] **Step 3: Create the test-timing-report script**

**Rationale:** Per-stage observability for test runtime. Not a CI-failing gate — purely informational. Surfaces any test file whose runtime grew ≥ 2× vs main baseline so the executor can investigate (most likely a real new mechanism added to MECHANISM_TEST_SPECS, but occasionally a refactor-introduced perf regression).

```bash
cat > scripts/test_timing_report.sh <<'SCRIPT'
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
    t=$(julia --project -e "
        using Pkg
        Pkg.activate(\"test\")
        try
            t = @elapsed include(\"$f\")
            print(round(t; digits=2))
        catch e
            print(\"FAIL\")
        end
    " 2>/dev/null)
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
SCRIPT
chmod +x scripts/test_timing_report.sh
```

- [ ] **Step 4: Record the main baseline**

Capture per-test-file timing on current `main`:

```bash
bash scripts/test_timing_report.sh /dev/null \
    docs/superpowers/refactor-test-timings-main-baseline.txt
```

The script's first arg (`/dev/null`) means "no baseline to compare against"; the second arg writes the output to the canonical baseline location. Resulting file looks like:

```
# Per-test-file runtime report (informational)
# Format: <test_file> <elapsed_seconds>
test_accessors.jl 1.23
test_aqua_jet.jl 12.45
...
```

The baseline is committed as part of Stage 0's PR so the refactor branch can compare against it for every stage closeout.

### Task 0.5: Commit Stage 0 and open PR

- [ ] **Step 1: Verify git status**

`git status`
Expected: `test/test_compile_budget.jl`, `test/runtests.jl`, `scripts/check_test_integrity.sh`, `scripts/test_timing_report.sh`, `docs/superpowers/refactor-test-timings-main-baseline.txt`.

- [ ] **Step 2: Commit**

```bash
git add test/test_compile_budget.jl test/runtests.jl \
        scripts/check_test_integrity.sh scripts/test_timing_report.sh \
        docs/superpowers/refactor-test-timings-main-baseline.txt
git commit -m "$(cat <<'EOF'
Add compile-time + test-integrity + test-timing gates for refactor

- test/test_compile_budget.jl: 3 gates (init_mechanisms trace-compile
  + rate_equation body-build wall-clock + warmup-reuse ter-ter
  incremental). Budget FIXED at 2x main; do NOT recalibrate on trip
  — redesign instead.
- scripts/check_test_integrity.sh: 4 hard checks + 1 WARN enforcing
  spec §2 (no test deletion / weakening / skip / broken). Honors
  spec §2.1 narrow exception via docs/superpowers/refactor-deleted-tests.md.
- scripts/test_timing_report.sh: per-test-file runtime report;
  flags >=2x growth vs main baseline for investigation (informational,
  not a failing gate).
- docs/superpowers/refactor-test-timings-main-baseline.txt: frozen
  per-test-file runtime on current main for stage closeouts to
  compare against.

Trace-compile regex anchored to EnzymeRates module prefix (no suffix
enumeration) so renamed/new internal helpers are caught automatically.

src delta: -0 / +0 net 0 (test + script + baseline infra only)
EOF
)"
```

- [ ] **Step 3: Push + open PR**

```bash
git push -u origin add-compile-budget-tests
gh pr create --title "Add compile-time + test-integrity gates" \
  --body "$(cat <<\'EOF\'
## Summary
- `test/test_compile_budget.jl`: 2 minimal gates (init_mechanisms trace-compile + rate_equation body-build wall-clock), calibrated against main with 2× headroom.
- `scripts/check_test_integrity.sh`: 4 hard checks + 1 WARN check enforcing spec §2.

## Motivation
Foundational infrastructure for the upcoming concrete-types refactor (spec: `docs/superpowers/specs/2026-05-20-concrete-types-refactor-design.md`). Lands first so the refactor branch can rebase and pick it up.

## Test plan
- [x] `julia --project -e \'using Pkg; Pkg.test()\'` passes locally
- [x] `bash scripts/check_test_integrity.sh main` PASSES (no-op against fresh main)
EOF
)"
```

- [ ] **Step 4: Wait for CI; merge to main**

After review and CI passes, merge. Then proceed to Stage 1 on `refactor-to-concrete-types-instead-of-symbols`, rebasing onto the new main.


# Stage 1 — New struct family + EnzymeMechanism{Sig} repack

**Goal:** Add the complete new concrete struct family with TDD-first tests, repack `EnzymeMechanism{M, R}` → `EnzymeMechanism{Sig}` (cosmetic; no behavior change), implement the `name(p, m)` chokepoint accessor, and add `Mechanism`/`AllostericMechanism` non-parametric forms with bidirectional converters.

**Expected LOC delta:** +700 src (foundation; recouped Stage 3+); negligible test loss. End of Stage 1: cumulative ~+700.

**Files in scope:** `src/types.jl` (heavy), `src/mechanism_enumeration.jl` (light — converter goes through `Mechanism`), `src/EnzymeRates.jl` (no change), every `src/*.jl` file that has `where {M, R}` or `where {S, P, R, N}` clauses (mechanical refactor in Task 1.13).

## Stage 1 commit consolidation (process simplification per Reviewer B #3)

The 17 Tasks below were drafted as 17 separate commits to keep TDD granular. **In practice, several can be combined without weakening the per-commit "tests green" invariant**, because Tasks 1.2–1.9 only ADD new types with no existing consumers — they cannot regress anything in the existing src or tests (which still consume the old parametric forms).

**Recommended commit grouping (5–6 commits for Stage 1):**

| Commit | Tasks | Why grouped |
|---|---|---|
| **1.A** | Task 1.1 (branch setup) | Standalone — branch + rebase |
| **1.B** | Tasks 1.2–1.9 (new type definitions: Metabolite hierarchy, Residual, Species, Step, RegulatorySite, Parameter family, ReactantAtoms, RegulatorMults) | All pure additions; no existing consumer. TDD progression within the commit is fine. Single `Pkg.test()` at the end. |
| **1.C** | Tasks 1.10–1.12 (concrete EnzymeReaction + Mechanism + AllostericMechanism) | The "transition" commit: rename old EnzymeReaction → EnzymeReactionLegacy, introduce new concrete EnzymeReaction, add Mechanism/AllostericMechanism. All in one atomic commit so the tests are coherent at every intermediate state. |
| **1.D** | Task 1.13 (Sig repack) | Standalone — affects every accessor in src; keep isolated for easy revert if anything breaks. |
| **1.E** | Tasks 1.14–1.15 (sig conversion functions + Mechanism↔EnzymeMechanism converters) | Both define the lift/sink layer; conceptually one piece of work. |
| **1.F** | Task 1.16 (name(p, m) chokepoint) + Task 1.17 (Stage 1 closeout/tag) | Final accessor + verification. |

**The TDD checkboxes inside each task remain — they're the per-step execution discipline (write test → run-fail → implement → run-pass).** The commit boundary is what changes: instead of committing after each task, the executor commits at the end of each grouping above. The per-commit "tests green" invariant is preserved because the grouped tasks don't have intra-group dependencies that would force tests to fail mid-commit.

The executor is free to use finer granularity (commit per task) if helpful for review or rollback, but the above is the recommended cadence.

### Task 1.1: Set up branch + rebase onto Stage 0

- [ ] **Step 1: Checkout the refactor branch**

```bash
git checkout refactor-to-concrete-types-instead-of-symbols
git fetch origin main
git rebase origin/main
```

- [ ] **Step 2: Verify Stage 0 tests are present**

Run: `ls test/test_compile_budget.jl && grep -q test_compile_budget test/runtests.jl && echo OK`
Expected: prints `OK`.

- [ ] **Step 3: Run full suite to confirm green baseline**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: all tests pass (including new compile-budget tests).

### Task 1.2: Add complete type hierarchy + chokepoint accessor (combined per Reviewer D #2)

**Files:**
- Modify: `src/types.jl` (add all structs from spec §5.1–5.7 + accessors from §6)
- Modify: `test/test_types.jl` (add full coverage)

**Rationale:** Tasks 1.2–1.9 in the prior plan iteration each added one struct per commit (Substrate, Product, Regulators, Residual, Species, Step, RegulatorySite, Parameter family, ReactantAtoms, RegulatorMults). None has any existing src consumer, so per-commit isolation buys nothing. This combined task ships them all in one atomic commit with one `Pkg.test()` run.

The new types are PURE ADDITIONS — they coexist with the existing parametric `EnzymeReaction{S,P,R,N}` and `EnzymeMechanism{M,R}` until later tasks rewire consumers. No regression risk.

**Critical implementation details (must be honored):**

- **`Step` does NOT carry `source_idx`.** Today's behavior (verified at `src/mechanism_enumeration.jl:2199-2208` + `rate_eq_derivation.jl:125`) is that `_canonicalize!` sorts steps + re-numbers groups, then `rep = first(steps_in_group)` uses POST-canonicalized step position as rep_idx. The new `name(p, m)` mirrors this: `_rep_idx_for_step` walks `steps(m)`, finds the group containing `p.step`, returns the step's position in the flattened canonicalized step list. No `source_idx` field needed; matches today exactly; one less field to propagate through every `_expand_*` move.

- **`ReactantAtoms`, `RegulatorMults`, and `EnzymeReaction` inner constructors canonicalize ordering** (Reviewer C #5). Two equivalent constructions (`Substrate(:ATP), atoms=[:C=>10, :H=>16]` vs `atoms=[:H=>16, :C=>10]`) must produce `==` results — otherwise dedup misses and the `EnzymeMechanism{Sig}` type universe explodes. The inner constructors sort: atoms by element name; regulators by name; reactants by metabolite name; allowed_multiplicities ascending; allo_states paired with their ligands.

- **`Species` inner constructor sorts `bound`** (already in spec §5.3). `Residual` inner constructor sorts `added` and `subtracted` (canonical comparison).

**Steps (TDD per CLAUDE.md):**

- [ ] **Step 1: Write the full test file**

Create `@testset` blocks in `test/test_types.jl` covering every struct from spec §5.1–5.7. Each struct gets:
- Constructor validation (rejection of invalid inputs)
- Accessor return-type and value tests
- `==` and `hash` consistency (especially: different field-order constructions are `==`)
- Canonicalization invariants (sorted output regardless of input order)
- `Step.source_idx` semantics: rep within a kinetic group = lowest source_idx (mirrors CLAUDE.md "Parameter naming convention")

Specific must-have tests:

@testset "Step has no source_idx field — rep_idx comes from position in mechanism" begin
    e   = Species([], :E)
    e_s = Species([Substrate(:S)], :E)
    e_p = Species([Product(:P)], :E)
    s1 = Step(e, e_s, Substrate(:S), true)
    s2 = Step(e_s, e_p, nothing, false)
    s3 = Step(e, e_p, Product(:P), true)
    @test fieldnames(Step) == (:from_species, :to_species,
                                :bound_metabolite, :is_equilibrium)
    @test !(:source_idx in fieldnames(Step))
    # rep_idx tests live alongside the Mechanism + name(p, m) tests:
    # the rep of a group is the position of its first step in the
    # flattened canonicalized step list (matches today's
    # `_canonicalize!` + `first(steps_in_group)` behavior).
end

@testset "ReactantAtoms canonicalizes atom ordering" begin
    ra1 = ReactantAtoms(Substrate(:ATP), [:C => 10, :H => 16, :N => 5])
    ra2 = ReactantAtoms(Substrate(:ATP), [:N => 5, :H => 16, :C => 10])
    @test ra1 == ra2                                # canonical sort
    @test hash(ra1) == hash(ra2)
end

@testset "EnzymeReaction canonicalizes reactant + regulator ordering" begin
    r1 = EnzymeReaction(
        [ReactantAtoms(Substrate(:B), [:C => 1]),
         ReactantAtoms(Substrate(:A), [:C => 1]),
         ReactantAtoms(Product(:P),   [:C => 2])],
        [RegulatorMults(AllostericRegulator(:Y), [2]),
         RegulatorMults(AllostericRegulator(:X), [2])],
        [3, 1, 2],
    )
    r2 = EnzymeReaction(
        [ReactantAtoms(Substrate(:A), [:C => 1]),
         ReactantAtoms(Substrate(:B), [:C => 1]),
         ReactantAtoms(Product(:P),   [:C => 2])],
        [RegulatorMults(AllostericRegulator(:X), [2]),
         RegulatorMults(AllostericRegulator(:Y), [2])],
        [1, 2, 3],
    )
    @test r1 == r2
    @test hash(r1) == hash(r2)
end
```

Plus the per-struct constructor/accessor/equality test patterns (one `@testset` per type — Substrate, Product, AllostericRegulator, CompetitiveInhibitor, Residual, Species, Step, RegulatorySite, each of the 10 Parameter subtypes, ReactantAtoms, RegulatorMults). See prior plan iteration (commit history) for specific examples; the test surface is mechanical.

- [ ] **Step 2: Run tests to verify all FAIL with UndefVarError**

`julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`

Expected: failures for each unknown identifier. Use the failure list as the implementation checklist.

- [ ] **Step 3: Implement the struct family in src/types.jl**

Add to `src/types.jl` in dependency order per spec §5:

1. Metabolite hierarchy (§5.1): `abstract type Metabolite end`, `Reactant <: Metabolite`, `Substrate`/`Product`, `Regulator <: Metabolite`, `AllostericRegulator`/`CompetitiveInhibitor`. `name(::Metabolite) -> Symbol` accessors.
2. Residual (§5.2): struct with `added::Vector{Substrate}`, `subtracted::Vector{Product}`, inner constructor sorts both; `Residual()` empty default; `==`/`hash`.
3. Species (§5.3): struct with `bound::Vector{Metabolite}`, `conformation::Symbol`, `residual::Residual`; inner constructor sorts `bound`; `name(::Species)` renderer (e.g., `:E_A_B`, `:Estar`); accessors.
4. RegulatorySite (§5.5): struct with `ligands::Vector{AllostericRegulator}`, `multiplicity::Int`, `allo_states::Vector{Symbol}`; inner constructor validates length match + multiplicity ≥ 1.
5. Step (§5.4): struct with `from_species`, `to_species`, `bound_metabolite`, `is_equilibrium`, **`source_idx::Int`**; inner constructor canonicalizes binding/iso direction (use metabolite-presence test per Reviewer C #4 — `met in bound(from)` vs `met in bound(to)`, NOT length comparison); accessors `from_species`/`to_species`/`bound_metabolite`/`is_equilibrium`/`is_binding`/`is_iso`/`direction`.
6. Parameter family (§5.6): `abstract type Parameter end`, then `Kd`/`Kiso`/`Kon`/`Koff`/`Kfor`/`Krev` (carry `step::Step`, `state::Symbol`), `Kreg` (carries `site::RegulatorySite`, `ligand::AllostericRegulator`, `state::Symbol`), `Keq`/`Etot`/`Lallo` (bare). Each step-bound subtype: `==`/`hash` on fields, `is_t_state(p) = p.state === :T`, `governing_step(p) = p.step`.
7. ReactantAtoms + RegulatorMults (§5.7): bundling structs with inner constructors that sort their respective collections.
8. EnzymeReaction (§5.7): concrete struct with `reactants::Vector{ReactantAtoms}`, `regulators::Vector{RegulatorMults}`, `allowed_catalytic_multiplicities::Vector{Int}`; inner constructor sorts all three.
9. Chokepoint `name(p::Parameter, m)` accessor (spec §6 + chokepoint scope note). Implementation uses `governing_step(p).source_idx` directly — NO linear scan into the mechanism (per Reviewer C #1). For `Kreg`: site index comes from `findfirst(==(p.site), regulatory_sites(m))` (small N — OK; or pre-compute via `_site_idx_table`).

- [ ] **Step 4: Run full test suite — all existing tests still pass + new tests pass**

`julia --project -e 'using Pkg; Pkg.test()'`

Expected: existing tests untouched (they consume the old `EnzymeReaction{S,P,R,N}` and `EnzymeMechanism{M,R}` still). New tests in `test_types.jl` all pass.

- [ ] **Step 5: Verify chokepoint discipline**

```bash
# No Symbol("K...") / Symbol("k...") string building outside name(p, m).
grep -nE 'Symbol\("[Kk][_0-9]' src/types.jl | grep -v "name(::Kd\|name(::Kiso\|name(::Kon\|name(::Koff\|name(::Kfor\|name(::Krev\|name(::Kreg"
```

Expected: empty (no matches). Per `[[feedback-chokepoint-accessors-for-future-migrations]]` — all parameter Symbol building routes through `name(p, m)`.

- [ ] **Step 6: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "$(cat <<'EOF'
Stage 1.2: add complete concrete type hierarchy + chokepoint accessor

Implements spec §5.1-5.7: Metabolite/Reactant/Regulator hierarchy,
Residual, Species, Step (with source_idx), RegulatorySite, Parameter
family (10 subtypes), ReactantAtoms/RegulatorMults bundling structs,
new concrete EnzymeReaction. All constructors canonicalize ordering.
Step.source_idx preserves today's positional parameter-naming
convention via the name(p::Parameter, m) chokepoint accessor.

No existing src code consumes these types yet (they coexist with the
existing parametric EnzymeReaction{S,P,R,N} / EnzymeMechanism{M,R}
until Stages 1.10-1.13 rewire consumers). Pure additions; zero
regression risk.

src delta: +0 / +~500 net +~500, cumulative: +~500
EOF
)"
```

### Task 1.10: Add new concrete EnzymeReaction (alongside old parametric)

**Files:**
- Modify: `src/types.jl` (rename old `EnzymeReaction` to `EnzymeReactionLegacy` to free the public name for the new struct)
- Modify: every src file referencing `EnzymeReaction` to use `EnzymeReactionLegacy` for the transitional period
- Modify: `test/test_types.jl`

**Rationale:** Tests, DSL macros, and downstream code all reference `EnzymeReaction` by name. We want the public name `EnzymeReaction` to refer to the new concrete struct from this task onward. The old parametric type gets renamed to `EnzymeReactionLegacy` in this commit and deleted in Stage 2 (after DSL rewrite).

- [ ] **Step 1: Write failing tests for the NEW EnzymeReaction**

Add to `test/test_types.jl`:

```julia
@testset "EnzymeReaction (new concrete)" begin
    r = EnzymeRates.EnzymeReaction(
        [
            EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:ATP), [:C => 10]),
            EnzymeRates.ReactantAtoms(EnzymeRates.Product(:ADP), [:C => 10]),
        ],
        [
            EnzymeRates.RegulatorMults(EnzymeRates.AllostericRegulator(:cAMP), [2]),
        ],
        [1, 2],
    )

    @test length(EnzymeRates.reactants(r)) == 2
    @test EnzymeRates.allowed_catalytic_multiplicities(r) == [1, 2]
    @test length(EnzymeRates.regulators(r)) == 1

    @test EnzymeRates.substrates(r) ==
          [EnzymeRates.Substrate(:ATP)]
    @test EnzymeRates.products(r) == [EnzymeRates.Product(:ADP)]
end
```

- [ ] **Step 2: Rename old EnzymeReaction → EnzymeReactionLegacy across src — EVERY reference, atomically**

In `src/types.jl`: rename the existing parametric `struct EnzymeReaction{S,P,R,N}` definition and the `EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=...)` outer constructor and the show method to use the name `EnzymeReactionLegacy`. Keep all internal logic identical for now. Also rename any accessor methods like `substrates(::EnzymeReaction{S,P,R,N})` → `substrates(::EnzymeReactionLegacy{...})`.

**CRITICAL — tests-green invariant.** This step must update EVERY caller of `EnzymeReaction(...)` in src AND test in the same commit, so the entire test suite still passes after the commit. The new EnzymeReaction concrete struct (Step 3) takes over the public name; existing callers explicitly call `EnzymeReactionLegacy(...)` until Stage 2 rewrites them. Below is the exhaustive call-site update list:

```bash
# Find ALL constructor / type-reference sites:
grep -rn "EnzymeReaction" src/ test/
```

Update every match using these rules:

| Site | Action |
|---|---|
| `src/dsl.jl` — `@enzyme_reaction` macro emits `:(EnzymeReaction($subs, $prods, ...))` | Change emitted form to `:(EnzymeReactionLegacy($subs, $prods, ...))`. Macro keeps existing grammar; will be rewritten in Stage 2. |
| `src/mechanism_enumeration.jl` — `EnzymeMechanism(spec::MechanismSpec)` reads `reaction::EnzymeReaction` field | Change field type to `EnzymeReactionLegacy`. |
| `src/identify_rate_equation.jl` — `IdentifyRateEquationProblem{R<:EnzymeReaction, D}` and field type | Change to `R<:EnzymeReactionLegacy`. |
| `src/fitting.jl` — no direct EnzymeReaction reference (consumes via accessors); typically unchanged. Verify with grep. |
| `src/EnzymeRates.jl` — `export EnzymeReaction` | LEAVE this export; the public name is taken over by the new concrete struct in Step 3. Users still get `EnzymeReaction` from the module. |
| `test/test_types.jl` — direct `EnzymeReaction(...)` construction tests | Change to `EnzymeReactionLegacy(...)` for the OLD-form tests. The NEW-form tests (in this task's Step 1) construct `EnzymeReaction(...)` of the new concrete struct. |
| `test/test_dsl.jl` — DSL output assertions | DSL still emits Legacy after this commit; assertions check Legacy field accessors. No change to assertion *shape*; only the type name swaps. |
| `test/test_accessors.jl` — accessor calls | Type-asserts on EnzymeReaction become EnzymeReactionLegacy. |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` — fixture constructors | Change `EnzymeReaction(...)` calls to `EnzymeReactionLegacy(...)`. |

After Step 2 the tests pass with the rename in place; after Step 3 the new EnzymeReaction is added and tested via direct construction; Stage 2.1 rewrites the DSL to emit the new EnzymeReaction; Stage 2.4 deletes EnzymeReactionLegacy.

- [ ] **Step 3: Implement the NEW EnzymeReaction concrete struct**

Add to `src/types.jl` (the new EnzymeReaction takes the public name):

```julia
struct EnzymeReaction
    reactants::Vector{ReactantAtoms}
    regulators::Vector{RegulatorMults}
    allowed_catalytic_multiplicities::Vector{Int}
end

reactants(r::EnzymeReaction)                       = r.reactants
regulators(r::EnzymeReaction)                      = r.regulators
allowed_catalytic_multiplicities(r::EnzymeReaction) =
    r.allowed_catalytic_multiplicities

# Derived filters.
function substrates(r::EnzymeReaction)
    Substrate[metabolite(ra) for ra in r.reactants
              if metabolite(ra) isa Substrate]
end
function products(r::EnzymeReaction)
    Product[metabolite(ra) for ra in r.reactants
            if metabolite(ra) isa Product]
end

function Base.:(==)(a::EnzymeReaction, b::EnzymeReaction)
    a.reactants == b.reactants && a.regulators == b.regulators &&
        a.allowed_catalytic_multiplicities == b.allowed_catalytic_multiplicities
end
Base.hash(r::EnzymeReaction, h::UInt) =
    hash(r.reactants, hash(r.regulators,
        hash(r.allowed_catalytic_multiplicities, h)))
```

- [ ] **Step 4: Run full test suite — verify no break + new test passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass, including the new `EnzymeReaction (new concrete)` testset AND existing tests that consume `EnzymeReactionLegacy`.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
Stage 1.10: add new concrete EnzymeReaction; rename old to EnzymeReactionLegacy

The new EnzymeReaction is a plain field struct holding ReactantAtoms +
RegulatorMults + allowed_catalytic_multiplicities. The old parametric
EnzymeReaction{S,P,R,N} is renamed to EnzymeReactionLegacy and remains
in use by DSL macros and downstream code until Stage 2 (DSL rewrite).

src delta: +50 / +70 net +20, cumulative: +343
EOF
)"
```

### Task 1.11: Add Mechanism non-parametric struct

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing tests**

Add to `test/test_types.jl`:

```julia
@testset "Mechanism (non-parametric)" begin
    r = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        EnzymeRates.RegulatorMults[],
        [1],
    )

    e   = EnzymeRates.Species([], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

    s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
    s_iso  = EnzymeRates.Step(e_s, e_p, nothing, false)
    s_rel  = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

    m = EnzymeRates.Mechanism(r, [[s_bind], [s_iso], [s_rel]])

    @test EnzymeRates.reaction(m) == r
    @test EnzymeRates.steps(m) == [[s_bind], [s_iso], [s_rel]]
    @test EnzymeRates.kinetic_groups(m) == 1:3
    @test EnzymeRates.n_steps(m) == 3
    @test EnzymeRates.rep_step(m, 2) == s_iso

    # equality + hash
    m2 = EnzymeRates.Mechanism(r, [[s_bind], [s_iso], [s_rel]])
    @test m == m2
    @test hash(m) == hash(m2)
end
```

- [ ] **Step 2: Run test to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: FAIL — `Mechanism not defined`.

- [ ] **Step 3: Implement**

Add to `src/types.jl`:

```julia
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}
end

reaction(m::Mechanism) = m.reaction
steps(m::Mechanism)    = m.steps
kinetic_groups(m::Mechanism) = 1:length(m.steps)
n_steps(m::Mechanism)        = sum(length, m.steps; init=0)
rep_step(m::Mechanism, g::Int) = first(m.steps[g])

function Base.:(==)(a::Mechanism, b::Mechanism)
    a.reaction == b.reaction && a.steps == b.steps
end
Base.hash(m::Mechanism, h::UInt) = hash(m.reaction, hash(m.steps, h))
```

- [ ] **Step 4: Run test to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Stage 1.11: add Mechanism non-parametric struct + accessors

src delta: +0 / +25 net +25, cumulative: +368"
```

### Task 1.12: Add AllostericMechanism non-parametric struct

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing tests**

Add to `test/test_types.jl`:

```julia
@testset "AllostericMechanism (non-parametric)" begin
    r = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        [EnzymeRates.RegulatorMults(EnzymeRates.AllostericRegulator(:A), [2])],
        [2],
    )

    e   = EnzymeRates.Species([], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

    cat_steps = [
        [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
        [EnzymeRates.Step(e_s, e_p, nothing, false)],
        [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
    ]
    cat_allo_states = [:EqualRT, :EqualRT, :EqualRT]

    site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A)], 2, [:OnlyR])

    am = EnzymeRates.AllostericMechanism(r, cat_steps, cat_allo_states, 2, [site])

    @test EnzymeRates.reaction(am) == r
    @test EnzymeRates.steps(am) == cat_steps
    @test EnzymeRates.cat_allo_state(am, 1) === :EqualRT
    @test EnzymeRates.cat_allo_state(am, 2) === :EqualRT
    @test EnzymeRates.cat_allo_state(am, 3) === :EqualRT
    @test EnzymeRates.catalytic_multiplicity(am) == 2
    @test EnzymeRates.regulatory_sites(am) == [site]
    @test EnzymeRates.allosteric_regulators(am) ==
          [EnzymeRates.AllostericRegulator(:A)]

    # :OnlyT catalytic groups must be rejected (R-state-active convention).
    @test_throws ErrorException EnzymeRates.AllostericMechanism(
        r, cat_steps, [:OnlyT, :EqualRT, :EqualRT], 2, [site])

    # Wrong-length cat_allo_states must be rejected.
    @test_throws ErrorException EnzymeRates.AllostericMechanism(
        r, cat_steps, [:EqualRT, :EqualRT], 2, [site])
end
```

- [ ] **Step 2: Run test to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: FAIL — `AllostericMechanism not defined`.

- [ ] **Step 3: Implement**

Add to `src/types.jl`:

```julia
const _VALID_CAT_ALLO_STATES = (:OnlyR, :EqualRT, :NonequalRT)

struct AllostericMechanism
    reaction::EnzymeReaction
    cat_steps::Vector{Vector{Step}}
    cat_allo_states::Vector{Symbol}    # parallel to cat_steps
    catalytic_multiplicity::Int
    regulatory_sites::Vector{RegulatorySite}

    function AllostericMechanism(reaction::EnzymeReaction,
                                  cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol},
                                  catalytic_multiplicity::Int,
                                  regulatory_sites::Vector{RegulatorySite})
        length(cat_allo_states) == length(cat_steps) ||
            error("AllostericMechanism: cat_allo_states length " *
                  "$(length(cat_allo_states)) must match cat_steps length " *
                  "$(length(cat_steps))")
        catalytic_multiplicity >= 1 ||
            error("AllostericMechanism: catalytic_multiplicity must be ≥ 1")
        for (g, tag) in enumerate(cat_allo_states)
            tag in _VALID_CAT_ALLO_STATES ||
                error("AllostericMechanism: catalytic group $g has " *
                      "invalid allo state $tag (must be one of " *
                      "$_VALID_CAT_ALLO_STATES); :OnlyT is rejected for " *
                      "catalytic groups (R-state-active convention)")
        end
        new(reaction, cat_steps, cat_allo_states,
            catalytic_multiplicity, regulatory_sites)
    end
end

reaction(m::AllostericMechanism)              = m.reaction
steps(m::AllostericMechanism)                 = m.cat_steps
cat_allo_state(m::AllostericMechanism, g::Int) = m.cat_allo_states[g]
catalytic_multiplicity(m::AllostericMechanism) = m.catalytic_multiplicity
regulatory_sites(m::AllostericMechanism)      = m.regulatory_sites
kinetic_groups(m::AllostericMechanism)        = 1:length(m.cat_steps)
n_steps(m::AllostericMechanism)               = sum(length, m.cat_steps; init=0)
rep_step(m::AllostericMechanism, g::Int)      = first(m.cat_steps[g])

function allosteric_regulators(m::AllostericMechanism)
    seen = AllostericRegulator[]
    for site in m.regulatory_sites
        for lig in site.ligands
            lig in seen || push!(seen, lig)
        end
    end
    seen
end

function competitive_inhibitors(m::AllostericMechanism)
    # CompetitiveInhibitors live on the reaction's regulators list,
    # filtered by subtype.
    [regulator(rm) for rm in m.reaction.regulators
     if regulator(rm) isa CompetitiveInhibitor]
end

function Base.:(==)(a::AllostericMechanism, b::AllostericMechanism)
    a.reaction == b.reaction && a.cat_steps == b.cat_steps &&
        a.cat_allo_states == b.cat_allo_states &&
        a.catalytic_multiplicity == b.catalytic_multiplicity &&
        a.regulatory_sites == b.regulatory_sites
end
Base.hash(m::AllostericMechanism, h::UInt) =
    hash(m.reaction, hash(m.cat_steps,
        hash(m.cat_allo_states,
            hash(m.catalytic_multiplicity, hash(m.regulatory_sites, h)))))
```

- [ ] **Step 4: Run test to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Stage 1.12: add AllostericMechanism non-parametric struct + accessors

src delta: +0 / +50 net +50, cumulative: +418"
```

### Task 1.13: Repack EnzymeMechanism{M,R} → EnzymeMechanism{Sig}

**Files:**
- Modify: `src/types.jl` (the EnzymeMechanism and AllostericEnzymeMechanism struct definitions; every accessor with `where {M, R}` etc.)
- Modify: all src files with `where {M, R}` or `where {CM, CS, RS}` clauses (`src/rate_eq_derivation.jl`, `src/mechanism_enumeration.jl`, `src/thermodynamic_constr_for_rate_eq_derivation.jl`)
- Modify: tests if any have type-asserts on `EnzymeMechanism{M, R}` (mechanical adaptation)

**Rationale:** Pure cosmetic refactor — collapse the existing two-param shape to single-Sig where `Sig === (Metabolites, Reactions)`. This gives later stages a single-Sig API to evolve. Zero behavior change.

- [ ] **Step 1: Add transition test asserting unchanged behavior**

Add to `test/test_types.jl`:

```julia
@testset "EnzymeMechanism Sig repack (Stage 1.13)" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end

    # m is still EnzymeMechanism; type parameter is now Sig (one param).
    @test m isa EnzymeMechanism
    Sig = typeof(m).parameters[1]
    @test Sig isa Tuple
    @test length(Sig) == 2   # (metabolites_sig, reactions_sig)

    # Accessors return the same values as before.
    @test substrates(m) == (:S,)
    @test products(m)   == (:P,)
    @test n_steps(m)    == 3
end
```

- [ ] **Step 2: Run test to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: FAIL — the inner type assertions fail because `typeof(m).parameters` still has 2 elements.

- [ ] **Step 3: Repack EnzymeMechanism in src/types.jl**

Locate the existing definition:

```julia
struct EnzymeMechanism{Metabolites, Reactions} <: AbstractEnzymeMechanism end
```

Replace with:

```julia
struct EnzymeMechanism{Sig} <: AbstractEnzymeMechanism end
```

Find every accessor with `where {M, R}` etc. — they're at lines ~526–650 in types.jl. Replace each `where {M, R}` with `where {Sig}` and destructure: e.g.,

```julia
# Before:
substrates(::EnzymeMechanism{M}) where {M} = M[1]
# After:
substrates(::EnzymeMechanism{Sig}) where {Sig} = Sig[1][1]
```

Apply the same pattern to `products`, `regulators`, `reactions`, `equilibrium_steps`, `n_steps`, `kinetic_group`, `kinetic_groups`, `steps_in_group`, `enzyme_forms`, `stoich_matrix`, `metabolites`.

Convert the constructor `EnzymeMechanism(mets, rxns)` at lines ~105–187 to return `EnzymeMechanism{(mets, rxns)}()` instead of `EnzymeMechanism{mets, rxns}()` after validation.

- [ ] **Step 4: Repack AllostericEnzymeMechanism similarly**

Locate:

```julia
struct AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites} <: AbstractEnzymeMechanism end
```

Replace with:

```julia
struct AllostericEnzymeMechanism{Sig} <: AbstractEnzymeMechanism end
```

Where `Sig === (CatalyticMech, CatSites, RegSites)`. Repack constructor + all accessors with `where {CM, CS, RS}` clauses.

- [ ] **Step 5: Update downstream src files**

Search for the old type parameter forms:

```bash
grep -rn "where {M, R}\|where {CM, CS, RS}\|EnzymeMechanism{M,\|AllostericEnzymeMechanism{CM," src/
```

For each match, change to `where {Sig}` form. Mechanical rewrite. Common cases:

- `where {M <: EnzymeMechanism}` (when M is the whole type): stays the same — Sig refactor doesn't change this pattern.
- `where {M, R}` in inner accessor: change to `where {Sig}` with destructuring inside.
- `@generated function f(::M, ...) where {M <: EnzymeMechanism}`: stays the same.

- [ ] **Step 6: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass. Numerical outputs unchanged; only type-parameter spelling differs.

- [ ] **Step 7: Run compile-budget tests specifically**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_compile_budget.jl"])'`
Expected: trace-compile counts ~same as main baseline (within calibrated 1.5× budget). If significantly higher, the Sig repack changed @generated body-build cost — investigate before committing.

- [ ] **Step 8: Commit**

```bash
git add src/ test/
git commit -m "$(cat <<'EOF'
Stage 1.13: repack EnzymeMechanism{M,R} → EnzymeMechanism{Sig}

Pure cosmetic refactor: collapse the existing two-param type
signature to single-Sig where Sig === (Metabolites, Reactions)
content-equivalent. Same for AllostericEnzymeMechanism{CM,CS,RS}
→ AllostericEnzymeMechanism{Sig} with Sig === (CM, CS, RS).

Zero behavior change. All accessors destructure Sig internally.
Sets up later stages to evolve a single Sig API.

src delta: -20 / +30 net +10, cumulative: +428
EOF
)"
```

### Task 1.14: Add _sig_of and _mechanism_from_sig conversion functions

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing tests**

Add to `test/test_types.jl`:

```julia
@testset "_sig_of / _mechanism_from_sig roundtrip" begin
    # Use REAL atom data (not just [:C => 1]) to catch type-parameter
    # validity issues. Per Julia spec, Pair{Symbol,Int} is NOT a valid
    # type-parameter value; encoding must use Tuple{Symbol,Int} leaves.
    r = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:ATP),
                                    [:C => 10, :H => 16, :N => 5,
                                     :O => 13, :P => 3]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:ADP),
                                    [:C => 10, :H => 15, :N => 5,
                                     :O => 10, :P => 2])],
        EnzymeRates.RegulatorMults[],
        [1],
    )
    e   = EnzymeRates.Species([], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:ATP)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:ADP)], :E)

    m = EnzymeRates.Mechanism(r, [
        [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:ATP), true)],
        [EnzymeRates.Step(e_s, e_p, nothing, false)],
        [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:ADP), true)],
    ])

    sig = EnzymeRates._sig_of(m)
    @test sig isa Tuple

    m_recon = EnzymeRates._mechanism_from_sig(sig)
    @test m_recon == m   # roundtrip

    # CRITICAL: sig MUST be usable as a type parameter (no Pairs, no
    # Vectors). This will throw TypeError if any leaf is invalid.
    # Catches Reviewer A's Finding #1 if it ever regresses.
    em_type = EnzymeRates.EnzymeMechanism{sig}
    @test em_type <: EnzymeRates.EnzymeMechanism
    em_inst = em_type()
    @test em_inst isa EnzymeRates.EnzymeMechanism

    # Roundtrip through the type-parameter form preserves everything.
    @test EnzymeRates.Mechanism(em_inst) == m
end
```

- [ ] **Step 2: Run test to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: FAIL — `_sig_of not defined`.

- [ ] **Step 3: Implement conversion functions via recursive dispatch**

Add to `src/types.jl`:

```julia
# ─── Mechanism ↔ Sig (parametric ↔ non-parametric) conversion ──
#
# IMPORTANT: every leaf in `sig` MUST be a valid Julia type-parameter
# value (isbits, Symbol, type, or Tuple of those). `Pair{Symbol,Int}` is
# NOT valid as a type parameter — encode pairs as `Tuple{Symbol,Int}`.
# Vectors are NEVER valid — always wrap in `Tuple(...)`.
#
# Strategy: one polymorphic `_to_sig` function with a method per source
# type; one matching `_from_sig` family taking a `Val{:Kind}` discriminator
# to dispatch reconstruction. Two functions per type (not ten paired
# helpers).

# --- _to_sig: type → Tuple ---
_to_sig(s::Substrate)            = (:Substrate, s.name)
_to_sig(p::Product)              = (:Product, p.name)
_to_sig(r::AllostericRegulator)  = (:AllostericRegulator, r.name)
_to_sig(c::CompetitiveInhibitor) = (:CompetitiveInhibitor, c.name)

_to_sig(r::Residual) = (
    Tuple(_to_sig(m) for m in r.added),
    Tuple(_to_sig(m) for m in r.subtracted),
)

_to_sig(s::Species) = (
    Tuple(_to_sig(m) for m in s.bound),
    s.conformation,
    _to_sig(s.residual),
)

_to_sig(s::Step) = (
    _to_sig(s.from_species),
    _to_sig(s.to_species),
    s.bound_metabolite === nothing ? nothing : _to_sig(s.bound_metabolite),
    s.is_equilibrium,
)

_to_sig(ra::ReactantAtoms) = (
    _to_sig(ra.metabolite),
    Tuple((p.first, p.second) for p in ra.atoms),    # (Symbol, Int) tuples — NOT Pairs
)

_to_sig(rm::RegulatorMults) = (
    _to_sig(rm.regulator),
    Tuple(rm.allowed_multiplicities),
)

_to_sig(r::EnzymeReaction) = (
    Tuple(_to_sig(ra) for ra in r.reactants),
    Tuple(_to_sig(rm) for rm in r.regulators),
    Tuple(r.allowed_catalytic_multiplicities),
)

# --- _from_sig: Tuple → type (dispatch on the leading kind Symbol) ---
function _metabolite_from_sig(sig::Tuple{Symbol, Symbol})
    kind, nm = sig
    kind === :Substrate             ? Substrate(nm)            :
    kind === :Product               ? Product(nm)              :
    kind === :AllostericRegulator   ? AllostericRegulator(nm)  :
    kind === :CompetitiveInhibitor  ? CompetitiveInhibitor(nm) :
    error("Unknown metabolite kind in sig: $kind")
end

function _residual_from_sig(sig::Tuple)
    added_sig, sub_sig = sig
    Residual(
        Substrate[_metabolite_from_sig(t) for t in added_sig if t[1] === :Substrate],
        Product[_metabolite_from_sig(t)   for t in sub_sig   if t[1] === :Product],
    )
end

function _species_from_sig(sig::Tuple)
    bound_sig, conformation, residual_sig = sig
    Species(
        Metabolite[_metabolite_from_sig(t) for t in bound_sig],
        conformation,
        _residual_from_sig(residual_sig),
    )
end

function _step_from_sig(sig::Tuple)
    from_sig, to_sig, met_sig, is_eq = sig
    met = met_sig === nothing ? nothing : _metabolite_from_sig(met_sig)
    Step(_species_from_sig(from_sig), _species_from_sig(to_sig), met, is_eq)
end

function _reactant_atoms_from_sig(sig::Tuple)
    met_sig, atoms_sig = sig
    ReactantAtoms(
        _metabolite_from_sig(met_sig)::Reactant,                 # type assert
        Pair{Symbol,Int}[s => c for (s, c) in atoms_sig],
    )
end

function _regulator_mults_from_sig(sig::Tuple)
    reg_sig, mults_sig = sig
    RegulatorMults(
        _metabolite_from_sig(reg_sig)::Regulator,                # type assert
        Int[m for m in mults_sig],
    )
end

function _reaction_from_sig(sig::Tuple)
    reactants_sig, regulators_sig, mults_sig = sig
    EnzymeReaction(
        ReactantAtoms[_reactant_atoms_from_sig(t) for t in reactants_sig],
        RegulatorMults[_regulator_mults_from_sig(t) for t in regulators_sig],
        Int[m for m in mults_sig],
    )
end

function _steps_from_sig(sig::Tuple)
    Vector{Step}[Step[_step_from_sig(s) for s in group] for group in sig]
end

# Public roundtrip functions.
_sig_of(m::Mechanism) = (_to_sig(m.reaction),
                          Tuple(Tuple(_to_sig(s) for s in g) for g in m.steps))

function _mechanism_from_sig(sig::Tuple)
    reaction_sig, steps_sig = sig
    Mechanism(_reaction_from_sig(reaction_sig), _steps_from_sig(steps_sig))
end
```

The recursive `_to_sig` collapses what was previously 6 named conversion
functions (`_reaction_to_sig`, `_species_to_sig`, `_step_to_sig`, etc.)
into one polymorphic dispatch family. Same for `_from_sig`. Code surface
~30% smaller than the original draft AND fixes the `Pair`-in-type-param
bug that would otherwise throw `TypeError` on the first `EnzymeMechanism`
construction with real atom data.

- [ ] **Step 4: Run test to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Stage 1.14: add _sig_of / _mechanism_from_sig roundtrip conversion

Plain Julia functions (NOT @generated) that bridge non-parametric
Mechanism ↔ Tuple sig used as EnzymeMechanism{Sig} type parameter.
Called from inside @generated rate_equation body at body-build time
in later stages.

src delta: +0 / +120 net +120, cumulative: +548"
```

### Task 1.15: Add EnzymeMechanism(::Mechanism) and Mechanism(::EnzymeMechanism) converters

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing tests**

Add to `test/test_types.jl`:

```julia
@testset "Mechanism ↔ EnzymeMechanism converters" begin
    r = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        EnzymeRates.RegulatorMults[],
        [1],
    )
    e   = EnzymeRates.Species([], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

    m = EnzymeRates.Mechanism(r, [
        [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
        [EnzymeRates.Step(e_s, e_p, nothing, false)],
        [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
    ])

    em = EnzymeMechanism(m)
    @test em isa EnzymeMechanism

    # Roundtrip via converter.
    m_back = EnzymeRates.Mechanism(em)
    @test m_back == m
end
```

- [ ] **Step 2: Run test to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: FAIL — `EnzymeMechanism(::Mechanism) not defined`.

- [ ] **Step 3: Implement converters**

Add to `src/types.jl`:

```julia
EnzymeMechanism(m::Mechanism) = EnzymeMechanism{_sig_of(m)}()
Mechanism(::EnzymeMechanism{Sig}) where {Sig} = _mechanism_from_sig(Sig)

# Allosteric counterparts use the same recursive _to_sig family.
_to_sig(s::RegulatorySite) = (
    Tuple(_to_sig(l) for l in s.ligands),
    s.multiplicity,
    Tuple(s.allo_states),
)

function _reg_site_from_sig(sig::Tuple)
    ligs_sig, mult, tags_sig = sig
    RegulatorySite(
        AllostericRegulator[_metabolite_from_sig(t)::AllostericRegulator
                            for t in ligs_sig],
        mult,
        Symbol[t for t in tags_sig],
    )
end

function AllostericEnzymeMechanism(am::AllostericMechanism)
    sig = (_to_sig(am.reaction),
           Tuple(Tuple(_to_sig(s) for s in g) for g in am.cat_steps),
           Tuple(am.cat_allo_states),
           am.catalytic_multiplicity,
           Tuple(_to_sig(rs) for rs in am.regulatory_sites))
    AllostericEnzymeMechanism{sig}()
end

function AllostericMechanism(::AllostericEnzymeMechanism{Sig}) where {Sig}
    reaction_sig, cat_steps_sig, cat_allo_states_sig,
        cat_mult, reg_sites_sig = Sig
    AllostericMechanism(
        _reaction_from_sig(reaction_sig),
        _steps_from_sig(cat_steps_sig),
        Symbol[s for s in cat_allo_states_sig],
        cat_mult,
        RegulatorySite[_reg_site_from_sig(t) for t in reg_sites_sig],
    )
end
```

- [ ] **Step 4: Run test to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Stage 1.15: add Mechanism ↔ EnzymeMechanism converters (both flavors)

src delta: +0 / +50 net +50, cumulative: +598"
```

### Task 1.16: Add name(p::Parameter, m) chokepoint implementations

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

**Rationale:** Single chokepoint for parameter Symbol production. Today's positional naming preserved. Future migration to path-based names is a single-function edit (see `feedback-chokepoint-accessors-for-future-migrations` memory note).

- [ ] **Step 1: Write failing tests**

Add to `test/test_types.jl`:

```julia
@testset "name(p::Parameter, m) chokepoint" begin
    r = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        EnzymeRates.RegulatorMults[],
        [1],
    )
    e   = EnzymeRates.Species([], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

    step1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
    step2 = EnzymeRates.Step(e_s, e_p, nothing, false)
    step3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

    m = EnzymeRates.Mechanism(r, [[step1], [step2], [step3]])

    # Positional naming: rep_idx for kinetic group g = position of first
    # step in steps(m)[g] within the flattened steps list.
    @test EnzymeRates.name(EnzymeRates.Kd(step1, :None), m) === :K1
    @test EnzymeRates.name(EnzymeRates.Kd(step1, :T),    m) === :K1_T
    @test EnzymeRates.name(EnzymeRates.Kon(step2, :None), m) === :k2f
    @test EnzymeRates.name(EnzymeRates.Koff(step2, :None), m) === :k2r
    @test EnzymeRates.name(EnzymeRates.Kfor(step2, :None), m) === :k2f
    @test EnzymeRates.name(EnzymeRates.Kd(step3, :None), m) === :K3

    @test EnzymeRates.name(EnzymeRates.Keq(),   m) === :Keq
    @test EnzymeRates.name(EnzymeRates.Etot(),  m) === :E_total
    @test EnzymeRates.name(EnzymeRates.Lallo(), m) === :L
end

@testset "name(p::Kreg, m) chokepoint" begin
    r = EnzymeRates.EnzymeReaction(
        [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
         EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
        [EnzymeRates.RegulatorMults(EnzymeRates.AllostericRegulator(:A), [2])],
        [2],
    )
    e   = EnzymeRates.Species([], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

    cat_steps = [
        [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
        [EnzymeRates.Step(e_s, e_p, nothing, false)],
        [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
    ]
    site_a = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:A)], 2, [:NonequalRT])

    am = EnzymeRates.AllostericMechanism(r, cat_steps, 2, [site_a])

    @test EnzymeRates.name(
            EnzymeRates.Kreg(site_a, EnzymeRates.AllostericRegulator(:A), :R),
            am) === :K_A_reg1
    @test EnzymeRates.name(
            EnzymeRates.Kreg(site_a, EnzymeRates.AllostericRegulator(:A), :T),
            am) === :K_A_T_reg1
end
```

- [ ] **Step 2: Run test to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: FAIL — `name(::Kd, ::Mechanism) method missing`.

- [ ] **Step 3: Implement**

Add to `src/types.jl`:

```julia
# rep_idx for a Parameter's Step = position of the group's first
# step in the flattened canonicalized step list. Matches today's
# `_canonicalize!` + `first(steps_in_group)` convention exactly
# (verified at src/mechanism_enumeration.jl:2199-2208 +
# rate_eq_derivation.jl:125).
#
# Linear scan: O(n_groups) to locate the group, then constant lookup.
# Called at @generated body-build only (not per call), and N is small
# (typically <30 steps). Memoize across all parameter names in a single
# body-build via a Dict{Step, Int} cache if the scan profile shows it
# matters — otherwise direct.
function _rep_idx_for_step(step::Step, m::Union{Mechanism, AllostericMechanism})
    groups = m isa Mechanism ? m.steps : m.cat_steps
    pos = 0
    for group in groups
        if step in group
            return pos + 1   # rep = first step in group; +1 since 1-indexed
        end
        pos += length(group)
    end
    error("Step not found in mechanism: $step")
end

function _rep_idx_for_step(step::Step, m::AbstractEnzymeMechanism)
    _rep_idx_for_step(step, _to_mechanism(m))
end

_to_mechanism(em::EnzymeMechanism)            = Mechanism(em)
_to_mechanism(aem::AllostericEnzymeMechanism) = AllostericMechanism(aem)

function _site_idx_of(site::RegulatorySite,
                      m::Union{AllostericMechanism, AllostericEnzymeMechanism})
    m_concrete = m isa AllostericMechanism ? m : AllostericMechanism(m)
    for (i, s) in enumerate(m_concrete.regulatory_sites)
        s == site && return i
    end
    error("RegulatorySite not found in mechanism")
end

# Chokepoint: name(p, m) renders today's positional Symbols using
# Step.source_idx-derived rep_idx. Future path-based migration:
# change these bodies only.
function name(p::Kd, m::Union{Mechanism, AbstractEnzymeMechanism})
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("K$(rep)_T") : Symbol("K$rep")
end
function name(p::Kiso, m::Union{Mechanism, AbstractEnzymeMechanism})
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("K$(rep)_T") : Symbol("K$rep")
end
function name(p::Kon, m::Union{Mechanism, AbstractEnzymeMechanism})
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("k$(rep)f_T") : Symbol("k$(rep)f")
end
function name(p::Koff, m::Union{Mechanism, AbstractEnzymeMechanism})
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("k$(rep)r_T") : Symbol("k$(rep)r")
end
function name(p::Kfor, m::Union{Mechanism, AbstractEnzymeMechanism})
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("k$(rep)f_T") : Symbol("k$(rep)f")
end
function name(p::Krev, m::Union{Mechanism, AbstractEnzymeMechanism})
    rep = _rep_idx_for_step(p.step, m)
    p.state === :T ? Symbol("k$(rep)r_T") : Symbol("k$(rep)r")
end
function name(p::Kreg, m::Union{AllostericMechanism, AllostericEnzymeMechanism})
    site_idx = _site_idx_of(p.site, m)
    lig_name = name(p.ligand)
    p.state === :T ? Symbol("K_$(lig_name)_T_reg$site_idx") :
                     Symbol("K_$(lig_name)_reg$site_idx")
end
name(::Keq,   _) = :Keq
name(::Etot,  _) = :E_total
name(::Lallo, _) = :L
```

- [ ] **Step 4: Run test to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_types.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "$(cat <<'EOF'
Stage 1.16: add name(p::Parameter, m) chokepoint accessor

Implements today's positional naming (:K1, :k1f, :K1_T, :K_<lig>_reg<i>)
through ONE function dispatching on Parameter subtype. Future
migration to path-based names is a single-function-body edit.

Per feedback-chokepoint-accessors-for-future-migrations memory note:
this is the central design enabler for the next refactor.

src delta: +0 / +90 net +90, cumulative: +688
EOF
)"
```

### Task 1.17: Stage 1 closeout — full test sweep + LOC delta check + test-integrity verification

**Files:** none (verification only)

- [ ] **Step 1: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass.

- [ ] **Step 2: Run compile-budget tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_compile_budget.jl"])'`
Expected: PASS within calibrated budgets for ALL five gates (init_mechanisms, expand_mechanisms, using EnzymeRates, rate_equation body-build, parameters body-build).

- [ ] **Step 3: Run rate_equation perf test specifically**

```bash
julia --project -e '
    include("test/test_rate_eq_derivation.jl")
    # The perf test runs over MECHANISM_TEST_SPECS; verify 0 allocs <100ns invariant.
'
```

Expected: PASS for every mechanism.

- [ ] **Step 4: Test-integrity verification (NON-NEGOTIABLE per spec §2)**

Spec §2: "Under NO circumstances may an existing test be deleted, commented out, or have its assertion logic changed. The ONLY change permitted to a test is mechanical syntax adaptation for the new struct surface." Enforced by:

```bash
bash scripts/check_test_integrity.sh main
```

The script (created in Stage 0 Task 0.4) runs four checks: no test file deleted, @testset count never decreases, no `@test_skip` / `@test_broken` introduced, no `@test` lines commented out.

If the script reports FAIL on any check: **STOP**. Revert the offending commit(s) and re-apply the change as mechanical syntax adaptation only (per spec §2). The stage is NOT complete until the script returns PASS. Do not proceed to the next stage.

- [ ] **Step N: Test-runtime report (informational — investigate if any file >=2x baseline)**

```bash
bash scripts/test_timing_report.sh
```

Per spec §3 implementor discipline: not a failing gate. The script reports per-test-file `@elapsed include()` time and flags any file whose runtime grew >=2x vs the main baseline (`docs/superpowers/refactor-test-timings-main-baseline.txt`, frozen at Stage 0). Investigate any 2x+ flag — most are legitimate (new mechanism added to `MECHANISM_TEST_SPECS`, new test added) but occasionally a refactor commit introduces an expensive helper. Document the cause in the stage commit message.


- [ ] **Step 5: Compute cumulative LOC delta**

```bash
wc -l src/*.jl
# Compare against main baseline: 7,136
```

Expected: cumulative around +688 (per task message accounting). If significantly different, investigate.

- [ ] **Step 6: Tag the stage**

```bash
git tag stage-1-complete
git log --oneline main..HEAD | head -20
```

Optional: push tag to remote for the mid-refactor checkpoint reference.

---

# Stage 2 — DSL rewrite

**Goal:** Replace the three DSL macros (`@enzyme_reaction`, `@enzyme_mechanism`, `@allosteric_mechanism`) with the prior branch's grammar, emitting constructors of the new concrete struct family directly (the new EnzymeReaction + Mechanism types from Stage 1). Delete the old DSL helpers. Delete `EnzymeReactionLegacy` (no longer used).

**Expected LOC delta:** −100 to −200 src. End of Stage 2: cumulative ~+500.

**Files in scope:** `src/dsl.jl` (full rewrite), `src/types.jl` (delete EnzymeReactionLegacy), `test/test_dsl.jl` (mechanical adaptation of assertions to new struct surface).

### Task 2.1: Rewrite @enzyme_reaction macro

**Files:**
- Modify: `src/dsl.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Read the prior branch's @enzyme_reaction implementation as reference**

```bash
git show refactor-to-use-structs-throughout:src/dsl.jl | sed -n '1,250p'
```

- [ ] **Step 2: Update test_dsl.jl assertions for new EnzymeReaction shape**

For each `@enzyme_reaction` testset in `test/test_dsl.jl`, find assertions that index into the old EnzymeReaction (e.g., `substrates(r)`, `regulators(r)`, etc.) and adapt them to expect the new return shape:

- Before: `substrates(r) == ((:S, ((:C, 6), (:H, 12), (:O, 6))),)` (tuple-of-tuples)
- After: `substrates(r) == [Substrate(:S)]` (Vector{Substrate})
- Before: `regulators(r) == ((:I, :dead_end),)` (tuple-of-pairs)
- After: `regulators(r) == [RegulatorMults(CompetitiveInhibitor(:I), [1])]` (Vector{RegulatorMults})

For each test, update the DSL input syntax to match the new prior-branch grammar:

- Old: `dead_end_inhibitors: I` → new: `competitive_inhibitors: I`
- Old: `oligomeric_state: 2` → new (still supported): `oligomeric_state: 2` (shorthand for `allowed_catalytic_multiplicities: (2,)`)
- New: per-regulator multiplicity required for allosteric regulators (`allosteric_regulators: A(1, 2)`)

**Do NOT delete or weaken any assertion**. Adapt only the input syntax and the expected output shape.

- [ ] **Step 3: Rewrite @enzyme_reaction macro implementation**

Replace the existing `@enzyme_reaction` macro and its helpers in `src/dsl.jl` with the implementation cribbed from the prior branch. Adapt to emit `EnzymeReaction(reactants::Vector{ReactantAtoms}, regulators::Vector{RegulatorMults}, mults::Vector{Int})` constructor calls of the NEW concrete struct (not `NonSingletonEnzymeReaction` from the prior branch).

Key code shape (from prior branch with adaptation):

```julia
macro enzyme_reaction(block)
    parsed = _parse_reaction_block(block)
    reactants_expr = _build_reactants_expr(parsed.subs, parsed.prods)
    regulators_expr = _build_regulators_expr(parsed.regs)
    mults_expr = _build_catalytic_mults_expr(parsed.mults)
    return esc(:(EnzymeReaction($reactants_expr, $regulators_expr, $mults_expr)))
end
```

`_parse_reaction_block`, `_parse_atom_bracket_entries`, `_parse_chemical_formula`, `_parse_regulator_entries`, `_parse_multiplicity_tuple`, `_parse_labeled_line` — port verbatim from prior branch's dsl.jl (lines ~54–200).

`_build_reactants_expr`, `_build_regulators_expr`, `_build_catalytic_mults_expr` — port and adapt to emit ReactantAtoms / RegulatorMults Expr nodes (lines ~210–260 of prior dsl.jl, replace the old `Substrate{Name}()` / `Product{Name}()` parametric calls with `Substrate(name)` / `Product(name)` concrete calls).

- [ ] **Step 4: Run test_dsl.jl**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_dsl.jl"])'`
Expected: PASS.

- [ ] **Step 5: Run full suite to catch downstream breakage**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: any test that uses `@enzyme_reaction` may need updating if its assertion looked at the old EnzymeReactionLegacy shape — adapt mechanically per spec §2.

- [ ] **Step 6: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Stage 2.1: rewrite @enzyme_reaction macro to prior-branch grammar

Emits EnzymeReaction(::Vector{ReactantAtoms}, ::Vector{RegulatorMults},
::Vector{Int}) constructor calls. Per-regulator multiplicities now
required for allosteric_regulators. dead_end_inhibitors: renamed to
competitive_inhibitors:.

src delta: -100 / +120 net +20, cumulative: +708"
```

### Task 2.2: Rewrite @enzyme_mechanism macro

**Files:**
- Modify: `src/dsl.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Read the prior branch's @enzyme_mechanism implementation**

```bash
git show refactor-to-use-structs-throughout:src/dsl.jl | sed -n '300,600p'
```

- [ ] **Step 2: Update test_dsl.jl assertions for new Mechanism shape**

For each `@enzyme_mechanism` testset, find assertions that index into the old EnzymeMechanism Sig (e.g., `reactions(m)`, `substrates(m)`) and adapt:

- Old: `substrates(m) == (:S,)` (Tuple{Symbol})
- After Stage 1: still `(:S,)` (accessor unchanged after Sig repack)
- For tests that should now look at the structural shape: `steps(Mechanism(m)) == [[Step(...)], [Step(...)]]`

Update the input syntax to the new grammar:

- Bare metabolite names in `substrates:` / `products:` / `regulators:`
- Function-call species notation `E(S)`, `Estar(B; residual = A - P)`
- Parenthesized step groups for shared kinetic parameters

- [ ] **Step 3: Rewrite @enzyme_mechanism macro implementation**

Replace the existing `@enzyme_mechanism` macro and its helpers in `src/dsl.jl`. Port from prior branch (lines ~310–600 of prior dsl.jl). Adapt to emit:

```julia
:($(Mechanism)($reaction_expr, $step_groups_expr))
```

where `step_groups_expr` is a `Vector{Vector{Step}}` constructed via the new `Step` outer constructor.

The macro needs to call `EnzymeReaction(...)` internally to build the reaction (using substrates/products/regulators declared at the top of the block, with default atom counts since `@enzyme_mechanism` doesn't take atom brackets). Use `[:C => 1]` as default atom signature for each metabolite in the mechanism-level macro (atoms only matter for `@enzyme_reaction`).

Then the macro emits `EnzymeMechanism(::Mechanism)` to deliver a parametric mechanism to the user — that's the user-facing entry point.

```julia
macro enzyme_mechanism(block)
    parsed = _parse_mechanism_block(block)
    reaction_expr = _build_default_reaction_expr(parsed)
    step_groups_expr = _build_step_groups_expr(parsed)
    return esc(:(EnzymeMechanism(Mechanism($reaction_expr, $step_groups_expr))))
end
```

`_parse_mechanism_block`, `_parse_steps_block`, `_parse_step_group`, `_build_step_expr`, `_split_side`, `_build_species_expr`, `_build_residual_expr`, `_walk_residual` — port from prior branch and adapt to emit `Step(...)` and `Species(...)` calls for the new concrete struct types.

- [ ] **Step 4: Run test_dsl.jl**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_dsl.jl"])'`
Expected: PASS.

- [ ] **Step 5: Run full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. Test fixtures that used the old `@enzyme_mechanism` grammar adapt mechanically.

- [ ] **Step 6: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Stage 2.2: rewrite @enzyme_mechanism to prior-branch grammar

Bare metabolite names + function-call species notation + parenthesized
step groups for shared kinetic parameters. Emits EnzymeMechanism(
Mechanism(reaction, step_groups)) constructor.

src delta: -200 / +180 net -20, cumulative: +688"
```

### Task 2.3: Rewrite @allosteric_mechanism macro

**Files:**
- Modify: `src/dsl.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Read prior branch's @allosteric_mechanism**

```bash
git show refactor-to-use-structs-throughout:src/dsl.jl | sed -n '600,900p'
```

- [ ] **Step 2: Update test_dsl.jl assertions**

For each `@allosteric_mechanism` testset, adapt input syntax (regulatory_site(multiplicity=N) blocks, ::AlloState annotations on cat steps) and expected output shape (AllostericMechanism + RegulatorySite vectors).

- [ ] **Step 3: Rewrite the macro**

Port from prior branch. Emit:

```julia
:($(AllostericEnzymeMechanism)($(AllostericMechanism)(
    $reaction_expr, $cat_steps_expr, $cat_mult, $reg_sites_expr)))
```

- [ ] **Step 4: Run test_dsl.jl + full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Stage 2.3: rewrite @allosteric_mechanism to prior-branch grammar

regulatory_site(multiplicity = N): begin ligands: ... end blocks
+ ::AlloState annotations on catalytic steps. Emits
AllostericEnzymeMechanism(AllostericMechanism(reaction, cat_steps,
cat_mult, reg_sites)).

src delta: -180 / +160 net -20, cumulative: +668"
```

### Task 2.4: Delete EnzymeReactionLegacy

**Files:**
- Modify: `src/types.jl` (delete the old parametric EnzymeReaction; rename or remove EnzymeReactionLegacy)
- Modify: any remaining src references

- [ ] **Step 1: Grep for remaining EnzymeReactionLegacy references**

```bash
grep -rn "EnzymeReactionLegacy" src/ test/
```

Expected: should only appear in `src/types.jl` (the struct itself) since DSL no longer constructs it.

- [ ] **Step 2: Delete EnzymeReactionLegacy struct + its accessors + its show method**

In `src/types.jl`, remove the `EnzymeReactionLegacy` struct definition, all its accessor method specializations, and the `Base.show(io::IO, ::EnzymeReactionLegacy)` method.

- [ ] **Step 3: Run full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/types.jl
git commit -m "Stage 2.4: delete EnzymeReactionLegacy (no longer used after DSL rewrite)

src delta: -180 / +0 net -180, cumulative: +488"
```

### Task 2.5: Stage 2 closeout

- [ ] **Step 1: Full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 2: Compile-budget check**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_compile_budget.jl"])'`
Expected: PASS within budgets for ALL 5 trace-compile gates (init_mechanisms, expand_mechanisms, using EnzymeRates, rate_equation body-build, parameters body-build) and ALL 4 wall-clock gates.

- [ ] **Step 3: Test-integrity verification (NON-NEGOTIABLE per spec §2)**

```bash
bash scripts/check_test_integrity.sh main
```

Spec §2: no test may be deleted, commented out, weakened, or `@test_skip`/`@test_broken`-tagged. Only mechanical syntax adaptation for new structs is permitted. If the script reports FAIL: **STOP**, revert the offending commit(s), re-apply as mechanical adaptation only. Stage is NOT complete until script returns PASS.

- [ ] **Step N: Test-runtime report (informational — investigate if any file >=2x baseline)**

```bash
bash scripts/test_timing_report.sh
```

Per spec §3 implementor discipline: not a failing gate. The script reports per-test-file `@elapsed include()` time and flags any file whose runtime grew >=2x vs the main baseline (`docs/superpowers/refactor-test-timings-main-baseline.txt`, frozen at Stage 0). Investigate any 2x+ flag — most are legitimate (new mechanism added to `MECHANISM_TEST_SPECS`, new test added) but occasionally a refactor commit introduces an expensive helper. Document the cause in the stage commit message.


- [ ] **Step 4: LOC delta**

```bash
wc -l src/*.jl
```

Expected: cumulative ~+488 (down from +688 after Stage 1). DSL grammar is more compact.

- [ ] **Step 5: Tag the stage**

```bash
git tag stage-2-complete
```

---

# Stage 3 — EnzymeMechanism derivation switch-over

**Goal:** Rewrite the `EnzymeMechanism` (non-allosteric) derivation pipeline to consume `Mechanism(em)` output and walk `Step`/`Species` structs directly. Delete Symbol classifier helpers (`is_k_parameter`, `_is_ss_rate_constant`, `_form_name`-driven code, `_split_reaction_side`'s Symbol version, etc.).

**Expected LOC delta:** −800 to −1200 src. End of Stage 3: cumulative target ≤ −300.

**Files in scope:** `src/rate_eq_derivation.jl` (heavy rewrite — non-allosteric sections), `src/thermodynamic_constr_for_rate_eq_derivation.jl` (rewrite), `src/sym_poly_for_rate_eq_derivation.jl` (replace `is_k_parameter`), `test/test_rate_eq_derivation.jl` (mechanical adaptation).

**Critical perf gate:** `test_rate_equation_performance` must stay green at every commit. The per-call body shape is unchanged (pure arithmetic Expr); only body-build code changes from "walk Symbol tuples" to "walk Step/Species structs".

### Task 3.1: Rewrite `parameters(::EnzymeMechanism, ::FullMode)` and `_raw_param_symbols` using Parameter struct walk

**Files:**
- Modify: `src/rate_eq_derivation.jl`

- [ ] **Step 1: Identify the code under replacement**

```bash
grep -n "parameters\|_raw_param_symbols\|_sorted_raw_param_symbols" src/rate_eq_derivation.jl
```

Read lines ~37–88 to understand the current implementation.

- [ ] **Step 2: Write the new helper using Parameter struct enumeration**

In `src/rate_eq_derivation.jl`, add (near the existing parameters definitions):

```julia
"""
Enumerate every Parameter for a mechanism, using positional naming.
Returns Vector{Parameter} in step+state order.
"""
function _enumerate_parameters_full(m::Mechanism)
    out = Parameter[]
    for (g_idx, group) in enumerate(m.steps)
        rep = first(group)
        if is_equilibrium(rep)
            if is_binding(rep)
                push!(out, Kd(rep, :None))
            else
                push!(out, Kiso(rep, :None))
            end
        else
            if is_binding(rep)
                push!(out, Kon(rep, :None))
                push!(out, Koff(rep, :None))
            else
                push!(out, Kfor(rep, :None))
                push!(out, Krev(rep, :None))
            end
        end
    end
    out
end
```

- [ ] **Step 3: Rewrite parameters(::EnzymeMechanism, ::FullMode)**

Replace the existing definition (around line 42):

```julia
@generated function parameters(::M, ::FullMode) where {M <: EnzymeMechanism}
    mech = _mechanism_from_sig(typeof(M()).parameters[1])
    params = _enumerate_parameters_full(mech)
    names = (name(p, mech) for p in params)
    Tuple((names..., :E_total))
end
```

- [ ] **Step 4: Run rate_eq_derivation tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_rate_eq_derivation.jl"])'`
Expected: PASS — same Symbols emitted in same order.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Stage 3.1: rewrite parameters(::EnzymeMechanism, Full) via Parameter struct walk

src delta: -30 / +25 net -5, cumulative: +483"
```

### Task 3.2: Rewrite `_dependent_param_exprs` to return Dict{Parameter, ...}

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl`

- [ ] **Step 1: Read current `_dependent_param_exprs` end-to-end**

```bash
cat src/thermodynamic_constr_for_rate_eq_derivation.jl
```

Identify how the function currently uses `Symbol("K$rep")` etc. to key the result Dict.

- [ ] **Step 2: Rewrite to use Parameter-keyed Dict + lifted Mechanism**

Replace the body of `_dependent_param_exprs(M::Type{<:EnzymeMechanism})`:

```julia
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    m = _mechanism_from_sig(typeof(M()).parameters[1])
    _dependent_param_exprs_for_mechanism(m)
end

function _dependent_param_exprs_for_mechanism(m::Mechanism)
    # ... rewrite using Step/Species walks and Parameter struct keys ...
    # Returns Dict{Parameter, Union{Symbol, Expr}}
end
```

(The body of `_dependent_param_exprs_for_mechanism` is a substantial rewrite of the Wegscheider/Haldane closure logic; preserve the algorithm and the cycle-detection structure but replace Symbol-keyed bookkeeping with Parameter struct keys. The Expr values remain Julia AST with leaves rendered via `name(p, m)`.)

- [ ] **Step 3: Update callers**

Find call sites of `_dependent_param_exprs`:

```bash
grep -rn "_dependent_param_exprs" src/
```

For each call site, adapt to the new Dict{Parameter, ...} return type. Callers that iterate over the result and previously did `Symbol → expr` lookup now do `Parameter → expr`, rendering the Symbol via `name(p, m)` when needed for splicing.

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_rate_eq_derivation.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/
git commit -m "Stage 3.2: rewrite _dependent_param_exprs to Dict{Parameter, ...}

Wegscheider/Haldane closure logic now keyed by Parameter structs;
Expr values unchanged in shape (Julia AST with rendered name leaves).

src delta: -200 / +150 net -50, cumulative: +433"
```

### Task 3.3: Rewrite `rate_equation(::EnzymeMechanism)` @generated body

**Files:**
- Modify: `src/rate_eq_derivation.jl` (rewrite `_build_rate_body`, `_raw_symbolic_rate_polys`, `_compute_alpha`, `_compute_re_groups`, `_split_reaction_side`)
- Modify: `test/test_rate_eq_derivation.jl` (mechanical adaptation only)

**Critical:** preserve the 0-alloc/<100ns invariant.

- [ ] **Step 1: Trace through existing rate_equation body-build**

Read `_build_rate_body` (~line 444 of rate_eq_derivation.jl), `_raw_rate_expr_and_symbols` (~line 491), `_raw_symbolic_rate_polys` (~line 294), `_compute_re_groups` (~line 183), `_compute_alpha` (~line 213), `_split_reaction_side` (~line 176), `_ss_contrib` (~line 279), `_compute_numerator` (look up).

Understand which use Symbol indexing.

- [ ] **Step 2: Rewrite `_split_reaction_side` to consume Step**

Replace:

```julia
function _split_reaction_side(side, enz_set)
    enzyme_sym = first(s for s in side if s in enz_set)
    met_syms = Symbol[s for s in side if s ∉ enz_set]
    return enzyme_sym, met_syms
end
```

with:

```julia
# Returns (enzyme_species, metabolites_in_side) for a step side.
# For our new Step struct, each side has exactly one Species and zero or
# one Metabolite — extract directly.
_lhs_species(s::Step)     = s.from_species
_lhs_metabolite(s::Step)  = s.bound_metabolite   # if iso, nothing
_rhs_species(s::Step)     = s.to_species
_rhs_metabolite(s::Step)  = nothing   # canonical form: met is on from side
```

- [ ] **Step 3: Rewrite `_compute_re_groups` to use Step**

```julia
"""Compute RE-connected groups via union-find over Species nodes."""
function _compute_re_groups(m::Mechanism)
    enz_species = _enumerate_species(m)
    N = length(enz_species)
    parent = collect(1:N)
    function find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
        x
    end
    for group in m.steps
        for s in group
            is_equilibrium(s) || continue
            i_from = findfirst(==(s.from_species), enz_species)
            i_to   = findfirst(==(s.to_species),   enz_species)
            ra, rb = find(i_from), find(i_to)
            ra != rb && (parent[ra] = rb)
        end
    end
    root_to_group = Dict{Int, Int}()
    groups = Vector{Vector{Int}}()
    form_to_group = zeros(Int, N)
    for i in 1:N
        r = find(i)
        g = get!(root_to_group, r) do; push!(groups, Int[]); length(groups) end
        push!(groups[g], i); form_to_group[i] = g
    end
    enz_species, groups, form_to_group
end

function _enumerate_species(m::Mechanism)
    seen = Species[]
    for group in m.steps
        for s in group
            s.from_species in seen || push!(seen, s.from_species)
            s.to_species   in seen || push!(seen, s.to_species)
        end
    end
    seen
end
```

- [ ] **Step 4: Rewrite `_compute_alpha` to use Step**

(Similar transformation — replace Symbol-form indexing with Species struct equality. Preserves the algorithm.)

- [ ] **Step 5: Rewrite `_raw_symbolic_rate_polys` and `_build_rate_body`**

Top-level rewrite to lift Sig → Mechanism once at body-build:

```julia
@generated function rate_equation(
    ::EnzymeMechanism{Sig}, concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {Sig}
    mech = _mechanism_from_sig(Sig)
    _build_rate_body_for_mechanism(mech, ReducedMode)
end

@generated function rate_equation(
    ::EnzymeMechanism{Sig}, concs::NamedTuple, params::NamedTuple, ::FullMode,
) where {Sig}
    mech = _mechanism_from_sig(Sig)
    _build_rate_body_for_mechanism(mech, FullMode)
end

function _build_rate_body_for_mechanism(m::Mechanism, mode)
    num, den = _raw_symbolic_rate_polys_for_mechanism(m)
    # ... build the arithmetic Expr identical in shape to today's body ...
    # leaves are param Symbols rendered via name(p, m); concs leaves use metabolite names.
end
```

Replace internal `Symbol("K$idx")` references with `name(Kd(step, :None), m)` etc. — the Symbol leaves of the emitted Expr are produced via the chokepoint accessor.

- [ ] **Step 6: Run rate_equation tests + perf gate**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_rate_eq_derivation.jl"])'`
Expected: PASS, including `test_rate_equation_performance` (0 allocs, <100ns per call).

If perf regresses, the issue is in the EMITTED body (per-call), not the body-BUILD. Inspect with `@code_llvm rate_equation(m, concs, params)` for a sample mechanism.

- [ ] **Step 7: Run compile-budget tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_compile_budget.jl"])'`
Expected: trace-compile counts within budget. Body-build cost may rise slightly (calling `_mechanism_from_sig` plus the struct walks) — that's fine as long as within budget.

- [ ] **Step 8: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Stage 3.3: rewrite rate_equation @generated body via Mechanism lift

@generated body now lifts Sig → Mechanism via plain function
_mechanism_from_sig at body-build time; walks Step/Species structs
directly. Per-call body shape unchanged (pure arithmetic Expr).

Deleted: _split_reaction_side (Symbol version), Symbol-keyed indexing
in _compute_re_groups, _compute_alpha, _compute_numerator.

Perf gate green: rate_equation 0 allocs, <100ns per call.

src delta: -400 / +200 net -200, cumulative: +233
EOF
)"
```

### Task 3.4: Rewrite `rate_equation_string`, `_kcat_forward`, kcat helpers

**Files:**
- Modify: `src/rate_eq_derivation.jl`

- [ ] **Step 1: Rewrite `rate_equation_string` to use Parameter walks**

Find `rate_equation_string(::M, ::FullMode)` and `rate_equation_string(::M, ::ReducedMode)` (~lines 549, 556). Adapt to lift Mechanism and use `name(p, m)` for header destructure lines instead of `_sorted_raw_param_symbols(M)`.

- [ ] **Step 2: Rewrite `_kcat_forward` and `_kcat_components`**

The `is_k_parameter`-based classifier in `_kcat_groups_from_polys` (~line 623) becomes a Parameter-type check. Adapt the monomial split:

```julia
function split_mono(mono::MONO, k_param_names::Set{Symbol})
    k_mono = MONO()
    met_mono = MONO()
    for (s, e) in mono
        if s in k_param_names || s == :Keq
            push!(k_mono, s => e)
        elseif s != :E_total
            push!(met_mono, s => e)
        end
    end
    sort!(k_mono; by=first), sort!(met_mono; by=first)
end
```

Pass `k_param_names = Set(name(p, m) for p in _enumerate_parameters_full(m))` from the caller.

- [ ] **Step 3: Delete `_is_ss_rate_constant`, `is_k_parameter` Symbol classifiers**

In `src/sym_poly_for_rate_eq_derivation.jl`, find `is_k_parameter` (~line 189) and delete it (move logic into the new Parameter-type-based classification).

In `src/rate_eq_derivation.jl`, find `_is_ss_rate_constant` (~line 613) and delete it.

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_rate_eq_derivation.jl"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/
git commit -m "Stage 3.4: rewrite rate_equation_string + _kcat_forward via Parameter walks

Deleted is_k_parameter (sym_poly_for_rate_eq_derivation.jl) and
_is_ss_rate_constant (rate_eq_derivation.jl) Symbol-string classifiers.
Classification now via Parameter subtype dispatch + name(p, m) chokepoint.

src delta: -250 / +100 net -150, cumulative: +83"
```

### Task 3.5: Rewrite `_build_kinetic_rename_map` or delete

**Files:**
- Modify: `src/rate_eq_derivation.jl`

**Rationale:** The kinetic-group rename map exists today because StepSpec uses `kinetic_group::Int` tags and the rep_idx might not equal the step index. With the new `Vector{Vector{Step}}` structural grouping, all steps in `steps(m)[g]` already share parameters — no rename map needed. The Pass-2 single-symbol Wegscheider absorption is also dead code per the memory note `project_dedup_pass2_dead_code`.

- [ ] **Step 1: Verify no production code calls `_build_kinetic_rename_map`**

```bash
grep -rn "_build_kinetic_rename_map" src/
```

If still used, refactor callers first.

- [ ] **Step 2: Delete `_build_kinetic_rename_map` and the Pass-1/Pass-2 rename code**

Remove the function definition (~line 119) and any internal callers that depended on it.

- [ ] **Step 3: Run full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Stage 3.5: delete _build_kinetic_rename_map (dead with structural groups)

src delta: -150 / +0 net -150, cumulative: -67"
```

### Task 3.6: Rewrite `src/thermodynamic_constr_for_rate_eq_derivation.jl`

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl`

- [ ] **Step 1: Identify the Symbol("K$idx") builders**

```bash
grep -n "Symbol(\"K\$\|Symbol(\"k\$" src/thermodynamic_constr_for_rate_eq_derivation.jl
```

Five locations: lines ~24, 27, 230, 239, 244, 269, 273.

- [ ] **Step 2: Rewrite each to use name(p, m) via Parameter struct**

For each Symbol-build site, replace with:

```julia
# Before:
Symbol("K$rep") → name(Kd(rep_step, :None), m)
Symbol("k$(rep)f") → name(Kon(rep_step, :None), m)   # or Kfor for iso steps
Symbol("k$(rep)r") → name(Koff(rep_step, :None), m)  # or Krev
```

Pass `m::Mechanism` through the call chain.

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "Stage 3.6: route Wegscheider Symbol building through name(p, m)

src delta: -100 / +60 net -40, cumulative: -107"
```

### Task 3.7: Stage 3 closeout

- [ ] **Step 1: Full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 2: Perf gate verification**

Specifically inspect `test_rate_equation_performance` results in the test output — `allocs == 0`, `t < 100e-9` for every mechanism in MECHANISM_TEST_SPECS.

- [ ] **Step 3: Compile-budget check (Stage 3 is the highest-risk stage for this)**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_compile_budget.jl"])'`

Expected: PASS within budget. This stage is the highest-risk for compile-time regression because it rewrites the @generated rate_equation/parameters bodies. If the gate trips:

**DO NOT RECALIBRATE THE BUDGET.** Per Task 0.3, the budget is fixed at 2×main for the entire refactor. A trip means this stage's implementation produces a real regression that must be **REDESIGNED**, not accepted. Likely causes and fixes:

- `_mechanism_from_sig` materializing the full Mechanism struct at every body-build adds N specializations per unique Sig. Fix: walk `Sig` directly in the @generated body without materializing concrete Mechanism — operate on the Tuple structure with explicit index access instead of struct field access.
- Per-step Species sub-tuples in Sig may be triggering more type specialization than necessary. Fix: flatten the Sig encoding so each step is one flat tuple rather than nested (species, species, metabolite, bool).
- The lift's `_from_sig` recursive dispatch may compile a fresh specialization for each unique tuple shape. Fix: use a single type-erased lift that returns boxed values.

If no design fix brings the gate within budget, ESCALATE — surface the issue at the mid-refactor checkpoint for Denis to decide whether the cost increase is genuinely worth the refactor benefit. Do NOT silently raise the budget.

- [ ] **Step 4: Test-integrity verification (NON-NEGOTIABLE per spec §2)**

```bash
bash scripts/check_test_integrity.sh main
```

Spec §2: no test may be deleted, commented out, weakened, or `@test_skip`/`@test_broken`-tagged. Only mechanical syntax adaptation for new structs is permitted. If the script reports FAIL: **STOP**, revert the offending commit(s), re-apply as mechanical adaptation only. Stage is NOT complete until script returns PASS.

- [ ] **Step N: Test-runtime report (informational — investigate if any file >=2x baseline)**

```bash
bash scripts/test_timing_report.sh
```

Per spec §3 implementor discipline: not a failing gate. The script reports per-test-file `@elapsed include()` time and flags any file whose runtime grew >=2x vs the main baseline (`docs/superpowers/refactor-test-timings-main-baseline.txt`, frozen at Stage 0). Investigate any 2x+ flag — most are legitimate (new mechanism added to `MECHANISM_TEST_SPECS`, new test added) but occasionally a refactor commit introduces an expensive helper. Document the cause in the stage commit message.


- [ ] **Step 5: LOC delta + tag**

```bash
wc -l src/*.jl
git tag stage-3-complete
```

Expected: cumulative around -107 (significant net deletion). Target was ≤ -300 by end of Stage 4 → on track.

---

# Stage 4 — AllostericEnzymeMechanism derivation switch-over (MID-REFACTOR CHECKPOINT)

**Goal:** Rewrite the AllostericEnzymeMechanism rate_equation, parameters, and T-state machinery to consume `AllostericMechanism(am)` output. Delete `_T_rename`, `_onlyR_syms`, `_all_t_state_names`, `_reg_param_name`, `_group_param_symbols`.

**Expected LOC delta:** −300 to −500 src. End of Stage 4: cumulative target ≤ −500.

**Files in scope:** `src/rate_eq_derivation.jl` (allosteric section, ~line 908 onward), `test/test_rate_eq_derivation.jl` (mechanical adaptation).

### Task 4.1: Rewrite `_T_rename` and friends via Parameter struct walks

**Files:**
- Modify: `src/rate_eq_derivation.jl`

- [ ] **Step 1: Identify allosteric Symbol helpers**

```bash
grep -n "_T_rename\|_onlyR_syms\|_all_t_state_names\|_reg_param_name\|_group_param_symbols\|_rename_params_T" src/rate_eq_derivation.jl
```

- [ ] **Step 2: Replace with Parameter struct enumerations**

For each helper, write a struct-based replacement. Examples:

```julia
"""Catalytic-cycle Parameters whose value is zeroed in the T-state branch
(any group with cat_allo_state == :OnlyR)."""
function _onlyR_parameters(am::AllostericMechanism)
    out = Parameter[]
    for (g, group) in enumerate(am.cat_steps)
        cat_allo_state(am, g) === :OnlyR || continue
        rep = first(group)
        if is_equilibrium(rep)
            push!(out, is_binding(rep) ? Kd(rep, :R) : Kiso(rep, :R))
        else
            if is_binding(rep)
                push!(out, Kon(rep, :R)); push!(out, Koff(rep, :R))
            else
                push!(out, Kfor(rep, :R)); push!(out, Krev(rep, :R))
            end
        end
    end
    out
end

"""R→T rename map for Parameters whose T-state name differs from R-state
(`:NonequalRT` groups)."""
function _T_rename_parameters(am::AllostericMechanism)
    rename = Dict{Parameter, Parameter}()
    for (g, group) in enumerate(am.cat_steps)
        cat_allo_state(am, g) === :NonequalRT || continue
        rep = first(group)
        if is_equilibrium(rep)
            T = is_binding(rep) ? Kd : Kiso
            rename[T(rep, :R)] = T(rep, :T)
        else
            if is_binding(rep)
                rename[Kon(rep, :R)]  = Kon(rep, :T)
                rename[Koff(rep, :R)] = Koff(rep, :T)
            else
                rename[Kfor(rep, :R)] = Kfor(rep, :T)
                rename[Krev(rep, :R)] = Krev(rep, :T)
            end
        end
    end
    rename
end

"""All T-state Parameters the rate-equation body emits as constraint LHSes."""
function _all_t_state_parameters(am::AllostericMechanism)
    out = Parameter[]
    for (g, group) in enumerate(am.cat_steps)
        tag = cat_allo_state(am, g)
        tag === :OnlyR && continue   # zeroed in T; not emitted
        rep = first(group)
        # ... per-step T-state Parameter emission, mirroring current
        # _all_t_state_names but returning struct instances ...
        if is_equilibrium(rep)
            T = is_binding(rep) ? Kd : Kiso
            push!(out, T(rep, :T))
        else
            if is_binding(rep)
                push!(out, Kon(rep, :T)); push!(out, Koff(rep, :T))
            else
                push!(out, Kfor(rep, :T)); push!(out, Krev(rep, :T))
            end
        end
    end
    # Regulator T-state Parameters
    for site in am.regulatory_sites
        for (lig, tag) in zip(site.ligands, site.allo_states)
            tag === :OnlyR && continue
            push!(out, Kreg(site, lig, :T))
        end
    end
    out
end
```

- [ ] **Step 3: Rewrite rate_equation(::AllostericEnzymeMechanism) body**

Mirror Stage 3.3 pattern: lift to AllostericMechanism, walk structs, emit body Expr.

- [ ] **Step 4: Run allosteric rate_equation tests + perf gate**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["test_rate_eq_derivation.jl"])'`
Expected: PASS, perf gate green.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Stage 4.1: rewrite allosteric rate_equation via AllostericMechanism lift

src delta: -350 / +200 net -150, cumulative: -257"
```

### Task 4.2: Rewrite `_dependent_param_exprs` for AllostericEnzymeMechanism

**Files:**
- Modify: `src/rate_eq_derivation.jl`

- [ ] **Step 1: Read current implementation**

Find `_dependent_param_exprs(::Type{AllostericEnzymeMechanism{...}})` around line 1090.

- [ ] **Step 2: Rewrite to use AllostericMechanism + Parameter structs**

Mirror Stage 3.2 pattern. T-state synthesized symbols become T-state Parameter struct instances.

- [ ] **Step 3: Test**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Stage 4.2: rewrite allosteric _dependent_param_exprs via Parameter structs

src delta: -150 / +80 net -70, cumulative: -327"
```

### Task 4.3: Rewrite `parameters(::AllostericEnzymeMechanism, ::FullMode)` and `_kcat_forward` for allosteric

**Files:**
- Modify: `src/rate_eq_derivation.jl`

- [ ] **Step 1: Rewrite `parameters(::AllostericEnzymeMechanism, ::FullMode)`**

Mirror Stage 3.1 pattern. Enumerate Parameters (Kd/Kiso/etc. + Kreg per ligand+site + T-state mirrors + Lallo + Etot).

- [ ] **Step 2: Rewrite `_kcat_forward(::AllostericEnzymeMechanism)`**

Find around line 800 area. Allosteric kcat has corner enumeration; preserve algorithm with struct-based name rendering.

- [ ] **Step 3: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/rate_eq_derivation.jl
git commit -m "Stage 4.3: rewrite allosteric parameters(Full) + _kcat_forward

src delta: -200 / +100 net -100, cumulative: -427"
```

### Task 4.4: Delete dead allosteric helpers

- [ ] **Step 1: Verify no remaining callers**

```bash
grep -rn "_rename_params_T\|_T_rename\|_onlyR_syms\|_all_t_state_names\|_reg_param_name\|_group_param_symbols" src/
```

- [ ] **Step 2: Delete each unused function**

- [ ] **Step 3: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/
git commit -m "Stage 4.4: delete dead allosteric Symbol helpers

Removed: _rename_params_T, _T_rename, _onlyR_syms, _all_t_state_names,
_reg_param_name, _group_param_symbols.

src delta: -180 / +0 net -180, cumulative: -607"
```

### Task 4.5: Rewrite `rate_equation_string` for allosteric

**Files:**
- Modify: `src/rate_eq_derivation.jl`

- [ ] **Step 1: Read current allosteric `rate_equation_string`**

Find around line 1379 area.

- [ ] **Step 2: Rewrite via Parameter walks**

Apply Stage 3.4 pattern to the allosteric variant.

- [ ] **Step 3: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/rate_eq_derivation.jl
git commit -m "Stage 4.5: rewrite allosteric rate_equation_string via Parameter walks

src delta: -80 / +40 net -40, cumulative: -647"
```

### Task 4.6: MID-REFACTOR CHECKPOINT

**STOP-AND-REVIEW gate per spec §11.** This is the hardest gate in the refactor; treat it accordingly.

- [ ] **Step 1: Cumulative LOC delta check**

```bash
wc -l src/*.jl
```

Expected: cumulative ≤ -500 (gate requirement). If not, STOP and revise the spec.

- [ ] **Step 2: Full test suite + perf gates**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: ALL tests pass; rate_equation perf gate green; compile-budget green for ALL 5 trace-compile + 4 wall-clock gates (init, expand, using, rate_equation body-build, parameters body-build).

- [ ] **Step 3: Test-integrity verification (NON-NEGOTIABLE per spec §2)**

```bash
bash scripts/check_test_integrity.sh main
```

Spec §2: no test may be deleted, commented out, weakened, or `@test_skip`/`@test_broken`-tagged. Only mechanical syntax adaptation for new structs is permitted. If the script reports FAIL: this is a serious violation given we're 4 stages in — STOP, revert the offending commit(s), re-apply as mechanical adaptation. Checkpoint is NOT cleared until script returns PASS.

- [ ] **Step N: Test-runtime report (informational — investigate if any file >=2x baseline)**

```bash
bash scripts/test_timing_report.sh
```

Per spec §3 implementor discipline: not a failing gate. The script reports per-test-file `@elapsed include()` time and flags any file whose runtime grew >=2x vs the main baseline (`docs/superpowers/refactor-test-timings-main-baseline.txt`, frozen at Stage 0). Investigate any 2x+ flag — most are legitimate (new mechanism added to `MECHANISM_TEST_SPECS`, new test added) but occasionally a refactor commit introduces an expensive helper. Document the cause in the stage commit message.


- [ ] **Step 4: Re-read spec end-to-end**

Open `docs/superpowers/specs/2026-05-20-concrete-types-refactor-design.md`. Verify Stages 5–7 plans still make sense given what we learned in Stages 1–4. If anything's wrong, AMEND the spec via a separate commit before proceeding.

- [ ] **Step 5: Decision point**

- If everything looks good: proceed to Stage 5.
- If LOC delta off-track: STOP, identify why, redesign Stages 5–7.
- If perf gate failing: STOP, fix root cause before proceeding.
- If test-integrity FAIL: STOP unconditionally, audit every stage commit for the violation, restore the deleted/weakened tests, do not proceed.

- [ ] **Step 6: Tag the checkpoint**

```bash
git tag stage-4-complete-checkpoint
```

---

# Stage 5 — Enumeration rewrite

**Goal:** Rewrite `init_mechanisms` / `expand_mechanisms` / `dedup!` to consume/produce `Vector{Mechanism}` directly using the new struct family. Delete `MechanismSpec`, `StepSpec`, `AllostericMechanismSpec`, and all the Symbol-form-name helpers (`_form_name`, `_parse_bound`, `_dead_end_form_name`, `_atoms_dict`, `_is_estar_form`).

**Expected LOC delta:** −1200 to −1500 src. End of Stage 5: cumulative target ≤ -1800.

**Files in scope:** `src/mechanism_enumeration.jl` (full rewrite, ~2,434 lines → target ~600–800 lines).

**Stage 5 stretch goal (per Reviewer D #5 — ~400 LOC additional savings):**

Beyond the deletion list above, the new struct surface enables three deeper simplifications that the per-task plan doesn't explicitly target. The executor should LOOK FOR these during Stage 5 implementation and pursue them when they're locally tractable; if all three land, Stage 5 ends ~400 LOC further negative than the +(-1200..-1500) baseline.

1. **`_parse_bound` / `_bound_metabolites_at_forms` chain collapses to a one-liner.** Lines ~892–953 of current `mechanism_enumeration.jl` greedily parse form names like `:E_A_B` back into `[:A, :B]`. With Species's `bound::Vector{Metabolite}` field, this is `bound.(unique_species_in_steps(m))` — one line. The plan's deletion list mentions `_parse_bound` but the saving is bigger than line count alone suggests because `_dead_end_form_name`, its callers, `_substrate_product_dead_end_opportunities`, and `_strip_reg_suffix` (lines 1287–1292) all exist as workarounds for Symbol-named-forms. Total: ~150 LOC additional.

2. **The 6 `_expand_*` moves share enormous structure.** Each looks up groups, copies, mutates, calls `_with_steps` (or a per-move variant). The per-move `_with_steps` divergence (`_dead_end_with_steps`, `_split_with_steps`) exists because each move has different tag-propagation rules to allosteric specs. With `Vector{Vector{Step}}` structural grouping and no `kinetic_group::Int` bookkeeping, most divergence collapses — there's no group renumbering. Look for unifying the 6 moves into one parametric `_expand(move::Symbol, spec)` or at least eliminating the per-move `_with_steps` variants. Estimated saving: ~200 LOC.

3. **`ANNOTATION_SUBSTITUTED` annotation / regex coupling between `rate_equation_string` (rate_eq_derivation.jl line 12) and the canonicalizer can be deleted in one pass** if Stage 5 and Stage 6 are sequenced carefully. Today the annotation exists to let the regex canonicalizer strip absorbed-tie display lines; the struct-based hash doesn't need that. Deleting the annotation removes the coupling between display + dedup that's currently spread across both files. Estimated saving: ~50 LOC.

These are STRETCH — pursue if convenient during the planned rewrites; don't add separate tasks for them unless they're nontrivial. Document any pursued in commit messages so cumulative LOC delta is clear at Stage 5 closeout.

### Task 5.1: Rewrite `_catalytic_topologies` to produce Vector{Vector{Step}}

**Files:**
- Modify: `src/mechanism_enumeration.jl`

- [ ] **Step 1: Read existing `_catalytic_topologies`**

```bash
grep -n "_catalytic_topologies\|backtrack" src/mechanism_enumeration.jl
```

Read the backtracking code (~line 150 onward).

- [ ] **Step 2: Adapt to produce Step instances directly**

Replace `_form_name(bound_subs, bound_prods, has_residual)` calls with direct `Species` construction:

```julia
# Before:
form = _form_name(bound_subs, bound_prods, has_residual)
# After:
species = Species(
    Metabolite[Substrate.(bound_subs)..., Product.(bound_prods)...],
    has_residual ? :Estar : :E,
)
```

Each catalytic topology generated becomes a `Vector{Step}` directly.

- [ ] **Step 3: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_mechanism_enumeration.jl"])'
git add src/mechanism_enumeration.jl
git commit -m "Stage 5.1: rewrite _catalytic_topologies to emit Step instances

src delta: -300 / +100 net -200, cumulative: -847"
```

### Task 5.2: Rewrite `init_mechanisms` to return Vector{Mechanism}

**Files:**
- Modify: `src/mechanism_enumeration.jl`

- [ ] **Step 1: Identify the wrapper**

Find `init_mechanisms(reaction)` and any helpers it calls (`_minimum_param_mechanisms`, etc.).

- [ ] **Step 2: Rewrite return type**

```julia
function init_mechanisms(reaction::EnzymeReaction)
    mechs = Mechanism[]
    for topology in _catalytic_topologies(reaction)
        # ... build Mechanism(reaction, [topology_grouped_steps]) ...
        push!(mechs, ...)
    end
    mechs
end
```

- [ ] **Step 3: Update return-type assertions in tests** (mechanical only)

- [ ] **Step 4: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_mechanism_enumeration.jl"])'
git add src/ test/test_mechanism_enumeration.jl
git commit -m "Stage 5.2: init_mechanisms now returns Vector{Mechanism}

src delta: -150 / +60 net -90, cumulative: -937"
```

### Task 5.3: Rewrite `expand_mechanisms` and each `_expand_*` move

**Files:**
- Modify: `src/mechanism_enumeration.jl`

The expand moves are:
- `_expand_re_to_ss(spec)`
- `_expand_split_kinetic_group(spec)` (formerly `_expand_remove_constraint`)
- `_expand_add_dead_end_regulator(spec)`
- `_expand_to_allosteric(spec)`
- `_expand_add_allosteric_regulator(spec)`
- `_expand_change_allo_state(spec)` (formerly `_expand_remove_tr_equiv`)

Each becomes a function `(::Mechanism) → Vector{Mechanism}` (or `(::AllostericMechanism) → Vector{AllostericMechanism}`).

- [ ] **Step 1: Rewrite `_expand_re_to_ss`**

(Specific code adaptation — read current implementation, replace StepSpec construction with Step construction.)

- [ ] **Step 2: Test the one move**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_mechanism_enumeration.jl"])'
```

- [ ] **Step 3: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Stage 5.3.1: rewrite _expand_re_to_ss for Mechanism

src delta: -80 / +30 net -50, cumulative: -987"
```

**Repeat Steps 1-3 for each of the remaining 5 expand moves**, committing per move. Expected per-move LOC delta: −50 to −100.

After all 6 moves rewritten, cumulative target: ~-1400.

### Task 5.4: Rewrite `dedup!` to use structural hash

**Files:**
- Modify: `src/mechanism_enumeration.jl`

- [ ] **Step 1: Identify current dedup**

`dedup!` at ~line 2260 uses `_canonicalize!` and `_dedup_key`. Both are Symbol-form-name based.

- [ ] **Step 2: Replace with struct-based hash**

```julia
function dedup!(cache::Dict{Int, Vector{Mechanism}})
    for (pc, specs) in cache
        for s in specs
            _canonicalize_mechanism!(s)
        end
        unique!(_dedup_key_mechanism, specs)
        isempty(specs) && delete!(cache, pc)
    end
    cache
end

_dedup_key_mechanism(m::Mechanism) = hash(_canonical_form(m))

function _canonicalize_mechanism!(m::Mechanism)
    # Sort steps within each group by struct-based key; sort groups by
    # canonical first-step key. Modifies in place.
    for group in m.steps
        sort!(group; by=_step_canonical_key)
    end
    sort!(m.steps; by=group -> _step_canonical_key(first(group)))
end

_step_canonical_key(s::Step) =
    (hash(s.from_species), hash(s.to_species),
     hash(s.bound_metabolite), s.is_equilibrium)
```

- [ ] **Step 3: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_mechanism_enumeration.jl"])'
git add src/mechanism_enumeration.jl
git commit -m "Stage 5.4: rewrite dedup! to struct-based hash

src delta: -150 / +40 net -110, cumulative: -1510"
```

### Task 5.5: Delete `MechanismSpec`, `StepSpec`, `AllostericMechanismSpec`

**Files:**
- Modify: `src/mechanism_enumeration.jl`

- [ ] **Step 1: Verify no remaining callers**

```bash
grep -rn "MechanismSpec\|StepSpec\|AllostericMechanismSpec" src/ test/
```

Expected: only the struct definitions themselves; tests use Mechanism directly.

- [ ] **Step 2: Delete the struct definitions + helpers**

In `src/mechanism_enumeration.jl`, delete the `StepSpec`, `MechanismSpec`, `AllostericMechanismSpec` struct definitions and their helpers (`all_form_names`, `step_metabolite`, `step_forms`, `_forms_with_binding_step`, etc.).

- [ ] **Step 3: Delete remaining Symbol-form-name helpers**

```bash
grep -n "_form_name\|_parse_bound\|_dead_end_form_name\|_atoms_dict\|_is_estar_form\|_can_pingpong\|_subtract_atoms" src/mechanism_enumeration.jl
```

Delete each function.

- [ ] **Step 4: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/mechanism_enumeration.jl
git commit -m "Stage 5.5: delete MechanismSpec/StepSpec/AllostericMechanismSpec + Symbol-form-name helpers

src delta: -500 / +0 net -500, cumulative: -2010"
```

### Task 5.6: Stage 5 closeout

- [ ] **Step 1: Full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 2: Compile-budget + perf gates**

Expected: PASS for all 5 trace-compile + 4 wall-clock gates AND `test_rate_equation_performance` 0-alloc/<100ns.

- [ ] **Step 3: Test-integrity verification (NON-NEGOTIABLE per spec §2)**

```bash
bash scripts/check_test_integrity.sh main
```

Spec §2: no test may be deleted, commented out, weakened, or `@test_skip`/`@test_broken`-tagged. Only mechanical syntax adaptation for new structs is permitted. If the script reports FAIL: **STOP**, revert the offending commit(s), re-apply as mechanical adaptation only. Stage is NOT complete until script returns PASS.

- [ ] **Step N: Test-runtime report (informational — investigate if any file >=2x baseline)**

```bash
bash scripts/test_timing_report.sh
```

Per spec §3 implementor discipline: not a failing gate. The script reports per-test-file `@elapsed include()` time and flags any file whose runtime grew >=2x vs the main baseline (`docs/superpowers/refactor-test-timings-main-baseline.txt`, frozen at Stage 0). Investigate any 2x+ flag — most are legitimate (new mechanism added to `MECHANISM_TEST_SPECS`, new test added) but occasionally a refactor commit introduces an expensive helper. Document the cause in the stage commit message.


- [ ] **Step 4: LOC delta + tag**

```bash
wc -l src/*.jl
git tag stage-5-complete
```

Expected: cumulative ≤ -1800.

---

# Stage 6 — identify_rate_equation canonicalizer

**Goal:** Replace the regex-over-string canonicalizer with a struct-based hash. Delete `_canonicalize_rate_eq_with_map`, `_sort_run_factors`, `_factor_sort_key`, regex pattern construction in `mechanism_enumeration.jl` and `identify_rate_equation.jl`.

**Expected LOC delta:** −300 to −500 src. End of Stage 6: cumulative target ≤ -2100.

**Files in scope:** `src/identify_rate_equation.jl`, `src/mechanism_enumeration.jl`, `test/test_identify_rate_equation.jl`.

### Task 6.0: Canonical-hash equivalence regression test (BLOCKING gate for Task 6.2)

**Files:**
- Create: `test/test_canonical_hash_equivalence.jl`
- Modify: `test/runtests.jl`

**Rationale (Reviewer A #6):** The existing regex canonicalizer strips single-symbol `(substituted into v)` ties so that mechanisms with different kinetic-group rename structures but identical `v` polynomials hash identically. If the new struct-based `_canonical_rate_eq_hash` (Task 6.1) doesn't produce the same equivalence classes, dedup will retain more candidates post-Stage 6 — blowing the runtime/compile budgets. This test runs both hashers in parallel and asserts equivalence-class agreement before the old one is deleted in Task 6.2.

- [ ] **Step 1: Write the equivalence test**

Create `test/test_canonical_hash_equivalence.jl`:

```julia
# ABOUTME: Verifies new struct-based _canonical_rate_eq_hash agrees with
# ABOUTME: the existing regex canonicalizer on equivalence-class assignment.

using Test
using EnzymeRates

@testset "canonical-hash equivalence (old regex vs new struct)" begin
    # Representative reactions covering the canonicalizer's edge cases:
    # - uni-uni: trivial structural equivalence
    # - bi-uni: kinetic-group merges of different orderings
    # - bi-bi: substituted-into-v ties
    # - ter-ter: largest enumeration; most opportunities for hash collision
    test_reactions = [
        @enzyme_reaction(begin
            substrates: S[C]
            products:   P[C]
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N], C[O]
            products:   P[C], Q[N], R[O]
        end),
    ]

    for reaction in test_reactions
        specs = EnzymeRates.init_mechanisms(reaction)
        cache = EnzymeRates.expand_mechanisms(specs, reaction)

        # Aggregate all candidates across param counts.
        all_specs = vcat(values(cache)...)

        # Bucket by OLD regex-based canonical hash.
        old_buckets = Dict{UInt64, Vector{Int}}()
        new_buckets = Dict{UInt64, Vector{Int}}()
        for (i, spec) in enumerate(all_specs)
            em = EnzymeRates.EnzymeMechanism(spec)  # or AllostericEnzymeMechanism
            old_h = EnzymeRates._canonical_rate_eq_hash_old(em)
            new_h = EnzymeRates._canonical_rate_eq_hash(em)
            push!(get!(old_buckets, old_h, Int[]), i)
            push!(get!(new_buckets, new_h, Int[]), i)
        end

        # Equivalence class agreement: same partition of specs into buckets
        # (modulo bucket key — what matters is which specs share a bucket).
        old_partition = Set(Set(b) for b in values(old_buckets))
        new_partition = Set(Set(b) for b in values(new_buckets))

        @test old_partition == new_partition
    end
end
```

- [ ] **Step 2: Extract the existing regex implementation into a private helper (per Reviewer C #6 — MANDATORY, not "either-or")**

The plan previously offered "duplicate body OR extract helper" — duplication risks drift if Task 6.1 modifies the body location. MANDATE the helper extraction:

```julia
# In src/mechanism_enumeration.jl (around line 2432):
# Step 2a: extract the existing implementation into a private helper.
function _canonical_rate_eq_hash_data_impl_regex(m::AbstractEnzymeMechanism)
    # ... existing body of _canonical_rate_eq_hash_data verbatim ...
end

# Step 2b: production name still points at the regex impl until Task 6.1.
function _canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)
    _canonical_rate_eq_hash_data_impl_regex(m)
end

function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)
    first(_canonical_rate_eq_hash_data(m))
end

# Step 2c: stable alias for the equivalence test to call.
function _canonical_rate_eq_hash_data_old(m::AbstractEnzymeMechanism)
    _canonical_rate_eq_hash_data_impl_regex(m)
end
function _canonical_rate_eq_hash_old(m::AbstractEnzymeMechanism)
    first(_canonical_rate_eq_hash_data_old(m))
end
```

After this refactor: `_canonical_rate_eq_hash_data_impl_regex` is the immutable witness of the old behavior. Task 6.1 introduces `_canonical_rate_eq_hash_data_impl_struct` and rewires `_canonical_rate_eq_hash_data` to call it — without touching the `_impl_regex` body. The equivalence test compares the two via the `_old` alias. Zero drift risk.

After Step 2 the equivalence test (Step 1) PASSES trivially (both names call the same impl). Task 6.1 then meaningfully diverges the two and the test does real work.

- [ ] **Step 3: Wire into runtests.jl**

```julia
# In test/runtests.jl, add:
include("test_canonical_hash_equivalence.jl")
```

- [ ] **Step 4: Run the test against current main (it should pass trivially since _canonical_rate_eq_hash and _canonical_rate_eq_hash_old reference the same body until Task 6.1)**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_canonical_hash_equivalence.jl"])'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_canonical_hash_equivalence.jl test/runtests.jl
git commit -m "$(cat <<'EOF'
Stage 6.0: add canonical-hash equivalence regression test

Preserves the old regex-based hash as _canonical_rate_eq_hash_old.
Adds a test asserting partition equivalence between old and new
canonical hashes across uni-uni / bi-bi / ter-ter representative
enumerations. Locked-in gate before Task 6.2 deletes the old hasher.

src delta: -0 / +30 net +30 (test infra; src unchanged), cumulative: <X>
EOF
)"
```

### Task 6.1: Write new struct-based `_canonical_rate_eq_hash`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (or move to a new home in src/identify_rate_equation.jl — choice during execution)
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write failing tests for the new hash + name_map**

Test that:
1. Two structurally-equivalent mechanisms (same `v` polynomial, different kinetic-group numbering) hash to the same UInt64 — same as the regex canonicalizer's behavior.
2. The hash function returns `(hash::UInt64, name_map::Dict{String,String})` matching the existing API shape (see `src/identify_rate_equation.jl:117,424-425,465-476,492-494`).
3. `_project_cached_params` continues to work: given cached params under representative mechanism `m_rep`, projecting to equivalent mechanism `m_eq` produces the right per-Parameter mapping.

```julia
@testset "_canonical_rate_eq_hash returns (hash, name_map) for param projection" begin
    # Two mechanisms with the same v polynomial but different
    # kinetic-group numbering should hash equal AND produce
    # name_maps that allow cached params to project across them.
    m_rep = ...   # canonical representative
    m_eq  = ...   # structurally equivalent variant

    h_rep, hex_rep, nm_rep = EnzymeRates._canonical_rate_eq_hash_data(m_rep)
    h_eq,  hex_eq,  nm_eq  = EnzymeRates._canonical_rate_eq_hash_data(m_eq)
    @test h_rep === h_eq
    @test hex_rep == hex_eq                                # hex display also stable
    @test length(hex_rep) == 16 && all(c -> c in "0123456789abcdef", hex_rep)
    @test nm_rep isa Dict{String, String}
    @test nm_eq  isa Dict{String, String}

    # Project cached params from m_rep onto m_eq via the maps.
    cached_params = (K1 = 1.0, K3 = 2.0, Keq = 0.5, E_total = 1.0)
    canon_to_rep = Dict(v => k for (k, v) in nm_rep)
    projected = EnzymeRates._project_cached_params(
        cached_params, canon_to_rep, nm_eq, EnzymeRates.fitted_params(m_eq))
    @test projected isa NamedTuple
    @test :Keq in keys(projected) && :E_total in keys(projected)
end
```

- [ ] **Step 2: Implement struct-based hash returning `(UInt64, name_map)`** (Reviewer C #2 — name_map is REQUIRED by identify_rate_equation.jl's param-projection contract)

The key insight: two mechanisms with the same algebraic rate equation but
different kinetic-group numbering (e.g., one has steps `[s1, s2, s3]` →
kinetic groups `[[s1, s2], [s3]]` and another has `[s3, s1, s2]` → `[[s3],
[s1, s2]]`) should hash identically. The canonical form must be invariant
under group renumbering.

The hash returns a tuple `(UInt64 hash, name_map::Dict{String,String})`
where `name_map` maps each parameter Symbol (as String, e.g., `"K1"`)
under this mechanism to its canonical token (e.g., `"p_1"`). Two
equivalent mechanisms produce the SAME canonical tokens for matching
roles — that's what enables `_project_cached_params` to relabel cached
fits across the equivalence class.

```julia
# Return contract is (UInt64, String, Dict{String,String}) — 3-tuple,
# matching the existing API at src/mechanism_enumeration.jl:2422
# (the h_short hex display is destructured into _CachedFitResult and
# surfaces as the user-facing eq_hash column in results).
function _canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)
    mech = _to_mechanism(m)
    canon, name_map = _canonicalize_for_hash(mech)
    h = hash(canon)
    (h, string(h, base=16, pad=16), name_map)
end

function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)
    first(_canonical_rate_eq_hash_data(m))
end

# Canonicalize: produce a canonical tuple representation + a
# String→String map from per-mechanism Symbol names to canonical
# tokens. The canonical tokens are stable across equivalent mechanisms.
function _canonicalize_for_hash(m::Mechanism)
    # Step 1: produce canonical Parameter ordering invariant to
    # position. Each Parameter gets a canonical key derived from its
    # role/Step structure (NOT its rendered Symbol).
    all_params = _enumerate_all_parameters_with_t_state(m)
    param_to_canon_key = Dict(
        p => _parameter_canonical_key(p) for p in all_params)

    # Step 2: assign canonical tokens by sorted canonical-key order.
    sorted_keys = sort(unique(values(param_to_canon_key));
                       by = repr)   # deterministic
    key_to_token = Dict(k => "p_$i" for (i, k) in enumerate(sorted_keys))

    # Step 3: name_map = per-mechanism Symbol name (as String) → token.
    name_map = Dict{String, String}()
    for p in all_params
        sym = name(p, m)
        token = key_to_token[param_to_canon_key[p]]
        name_map[String(sym)] = token
    end

    # Step 4: build the canonical form of the mechanism, using canonical
    # tokens for parameter leaves and structural hashes for steps/groups.
    groups_canon = Tuple(
        Tuple(sort([hash(s) for s in group]))
        for group in m.steps
    )
    groups_sorted = Tuple(sort(collect(groups_canon)))

    dep_exprs = _dependent_param_exprs_for_mechanism(m)
    dep_canon = sort([
        (key_to_token[_parameter_canonical_key(p)],
         _expr_canonical_via_name_map(expr, name_map))
        for (p, expr) in dep_exprs
    ])

    canon = (groups_sorted, Tuple(dep_canon))
    (canon, name_map)
end

# Per-Parameter canonical key independent of mechanism position.
function _parameter_canonical_key(p::Kd)   ; (:Kd,   hash(p.step), p.state) end
function _parameter_canonical_key(p::Kiso) ; (:Kiso, hash(p.step), p.state) end
function _parameter_canonical_key(p::Kon)  ; (:Kon,  hash(p.step), p.state) end
function _parameter_canonical_key(p::Koff) ; (:Koff, hash(p.step), p.state) end
function _parameter_canonical_key(p::Kfor) ; (:Kfor, hash(p.step), p.state) end
function _parameter_canonical_key(p::Krev) ; (:Krev, hash(p.step), p.state) end
function _parameter_canonical_key(p::Kreg)
    (:Kreg, hash(p.site), hash(p.ligand), p.state)
end
_parameter_canonical_key(::Keq)   = (:Keq,)
_parameter_canonical_key(::Etot)  = (:Etot,)
_parameter_canonical_key(::Lallo) = (:Lallo,)

# Walk Expr replacing Symbol leaves (parameter names) with their
# canonical tokens via name_map. Non-parameter Symbols (e.g., :Keq if
# unmapped, metabolite names) pass through unchanged.
function _expr_canonical_via_name_map(expr, name_map::Dict{String,String})
    if expr isa Symbol
        s = String(expr)
        return get(name_map, s, expr)
    end
    expr isa Expr || return expr
    Expr(expr.head,
         [_expr_canonical_via_name_map(a, name_map) for a in expr.args]...)
end
```

`_enumerate_all_parameters_with_t_state(m)` enumerates every Parameter the mechanism could emit names for — for non-allosteric: the catalytic step-bound family + Keq/Etot; for allosteric: catalytic + T-state mirrors + reg-site Kreg + Lallo. This is the same set covered by `parameters(m, Full)`.

**name_map equivalence proof obligation:** for two equivalent mechanisms `m_rep` and `m_eq`, the `name_map`s satisfy: there exists a bijection between their per-mechanism Symbols such that both maps assign the same canonical token to bijection-equivalent parameters. `_project_cached_params(cached, canon_to_rep, nm_eq, fitted_keys)` relies on this property to relabel cached fit values from `m_rep` onto `m_eq`. Task 6.0's equivalence test must verify this property on representative reaction fixtures (uni-uni, bi-bi, ter-ter) before Task 6.2's deletion.

- [ ] **Step 3: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_identify_rate_equation.jl"])'
git add src/
git commit -m "Stage 6.1: add struct-based _canonical_rate_eq_hash

src delta: -50 / +120 net +70, cumulative: -1940"
```

### Task 6.2: Delete regex canonicalizer (GATED on equivalence-test PASS)

- [ ] **Step 1: Confirm the equivalence test from Task 6.0 still passes**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["test_canonical_hash_equivalence.jl"])'
```

Expected: PASS. If FAIL, **DO NOT PROCEED** — the new hasher's partition diverges from old. Investigate which mechanisms hash differently and adjust Task 6.1's `_canonicalize_for_hash` before deleting the old code.

- [ ] **Step 2: Delete the regex canonicalizer**

In `src/mechanism_enumeration.jl`, delete `_canonicalize_rate_eq_with_map`, `_sort_run_factors`, `_factor_sort_key`, regex pattern construction, AND the `_canonical_rate_eq_hash_old` / `_canonical_rate_eq_hash_data_old` aliases preserved in Task 6.0.

- [ ] **Step 3: Repurpose the equivalence test as a permanent determinism + partition-stability regression test**

The equivalence test from Task 6.0 was a transient refactor gate (comparing old vs new hash). Once the old hash is gone, retarget the test to assert the new hash's properties that we care about long-term: (a) determinism (same mechanism → same hash across invocations) and (b) partition stability (the set of equivalence classes the new hash produces on the canonical reaction set doesn't change unexpectedly across future commits).

Rewrite `test/test_canonical_hash_equivalence.jl` (rename to `test/test_canonical_hash_partition.jl`) — same input data, but assert:

```julia
# Permanent regression test for _canonical_rate_eq_hash.
@testset "canonical-hash partition stability" begin
    test_reactions = [ ... same fixtures as Task 6.0 ... ]

    # Expected partition sizes per reaction. These are measured EMPIRICALLY
    # at Stage 6.1 (after the new hasher is wired in and the equivalence
    # test in Task 6.0 PASSED — meaning the new hash partition matches
    # the old hash partition exactly). The executor records the actual
    # numbers by running the test once and reading the printed values,
    # then commits the frozen literals here. If these change in a
    # future commit, the new hash has changed equivalence classes —
    # investigate before merging.
    #
    # MEASUREMENT PROCEDURE (Stage 6.1 closeout):
    #   1. Run the test with `expected_n_classes[label]` replaced by
    #      `@show length(new_buckets); 0` (forces test fail + prints count)
    #   2. Record the printed counts.
    #   3. Replace `0` with the actual count.
    #   4. Re-run; test must pass.
    expected_n_classes = Dict(
        "uni_uni"  => 0,  # ← measure at Stage 6.1, replace with literal
        "bi_bi"    => 0,  # ← measure at Stage 6.1, replace with literal
        "ter_ter"  => 0,  # ← measure at Stage 6.1, replace with literal
    )

    for (label, reaction) in test_reactions
        specs = EnzymeRates.init_mechanisms(reaction)
        cache = EnzymeRates.expand_mechanisms(specs, reaction)
        all_specs = vcat(values(cache)...)

        new_buckets = Dict{UInt64, Vector{Int}}()
        for (i, spec) in enumerate(all_specs)
            em = EnzymeRates.EnzymeMechanism(spec)
            h = EnzymeRates._canonical_rate_eq_hash(em)
            push!(get!(new_buckets, h, Int[]), i)
            # Determinism: same input, same hash.
            @test EnzymeRates._canonical_rate_eq_hash(em) === h
        end

        @test length(new_buckets) == expected_n_classes[label]
    end
end
```

The `expected_n_classes` numbers are FROZEN at Stage 6.1 (after the new hasher is wired in and the equivalence-test PASSED) — fill them in then. This test then catches any future regression of the partition (e.g., if someone changes `_canonicalize_for_hash` and accidentally narrows or broadens equivalence).

- [ ] **Step 4: Verify no remaining callers of the old hasher**

```bash
grep -rn "_canonicalize_rate_eq_with_map\|_sort_run_factors\|_factor_sort_key\|_canonical_rate_eq_hash_old\|_canonical_rate_eq_hash_data_old" src/ test/
```

Expected: no matches.

- [ ] **Step 5: Test + commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add src/ test/
git commit -m "Stage 6.2: delete regex canonicalizer; partition test now permanent regression

Old _canonical_rate_eq_hash_old + regex helpers deleted (gate cleared
in Task 6.0). test_canonical_hash_equivalence.jl renamed to
test_canonical_hash_partition.jl and repurposed as a permanent
determinism + partition-size regression test against frozen expected
class counts.

src delta: -400 / +0 net -400, cumulative: -2340"
```

### Task 6.3: Stage 6 closeout

- [ ] **Step 1: Full test suite + perf + compile-budget gates**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: PASS for all 5 trace-compile + 4 wall-clock gates AND `test_rate_equation_performance` 0-alloc/<100ns AND the new partition-stability regression test from Task 6.2.

- [ ] **Step 2: Test-integrity verification (NON-NEGOTIABLE per spec §2)**

```bash
bash scripts/check_test_integrity.sh main
```

Spec §2: no test may be deleted, commented out, weakened, or `@test_skip`/`@test_broken`-tagged. The repurposed `test_canonical_hash_partition.jl` (renamed in Task 6.2) is a permitted RENAME (not a deletion) — the script's Check 1 will correctly identify renames as not-deletions. Only mechanical syntax adaptation for new structs is permitted. If the script reports FAIL: **STOP**, revert the offending commit(s). Stage is NOT complete until script returns PASS.

- [ ] **Step N: Test-runtime report (informational — investigate if any file >=2x baseline)**

```bash
bash scripts/test_timing_report.sh
```

Per spec §3 implementor discipline: not a failing gate. The script reports per-test-file `@elapsed include()` time and flags any file whose runtime grew >=2x vs the main baseline (`docs/superpowers/refactor-test-timings-main-baseline.txt`, frozen at Stage 0). Investigate any 2x+ flag — most are legitimate (new mechanism added to `MECHANISM_TEST_SPECS`, new test added) but occasionally a refactor commit introduces an expensive helper. Document the cause in the stage commit message.


- [ ] **Step 3: LOC delta + tag**

```bash
wc -l src/*.jl
git tag stage-6-complete
```

---

# Stage 7 — Fitting touch-up + cleanup + docs

**Goal:** Final dead-code sweep, README + docstring updates, sanity passes. Cumulative LOC delta ≤ -3500.

**Expected LOC delta:** −300 to −500 src final cleanup.

### Task 7.1: Audit `src/fitting.jl` for any needed adaptations

- [ ] **Step 1: Re-read `src/fitting.jl`**

Most likely needs no changes — `loss!` uses `fitted_params`, `metabolites`, `rate_equation` which are public-API stable.

- [ ] **Step 2: If any adaptations needed, write test first + implement**

- [ ] **Step 3: Test + commit if changes**

### Task 7.2: Dead-code sweep across all src files

- [ ] **Step 1: Run unused-function detector**

```bash
for f in src/*.jl; do
    echo "=== $f ==="
    grep -oE '^function [a-zA-Z_][a-zA-Z_0-9!]*' "$f" | sed 's/function //' | \
    while read fn; do
        count=$(grep -r "\b$fn\b" src/ | wc -l)
        [ "$count" -le 1 ] && echo "POSSIBLY UNUSED: $fn"
    done
done
```

For each "POSSIBLY UNUSED" function, verify it's truly unused (check test files, external API) and delete if confirmed dead.

- [ ] **Step 2: Inline single-use private helpers**

For each private function (`_*`) called in exactly one place, consider inlining.

- [ ] **Step 3: Test after each round of deletions/inlinings**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 4: Commit each round**

```bash
git add src/
git commit -m "Stage 7.2: dead-code sweep — delete unused/inlined helpers

src delta: -X / +Y net Z, cumulative: -W"
```

### Task 7.3: Update README + docstrings

**Files:**
- Modify: `README.md`
- Modify: docstrings in `src/types.jl`, `src/dsl.jl`, `src/EnzymeRates.jl`

- [ ] **Step 1: Re-read current README**

```bash
cat README.md | head -100
```

- [ ] **Step 2: Update DSL examples to new grammar**

- [ ] **Step 3: Update architecture section to reflect new struct hierarchy**

- [ ] **Step 4: Commit**

```bash
git add README.md src/
git commit -m "Stage 7.3: update README + docstrings to new struct architecture

src delta: -0 / +0 net 0 (docs only)"
```

### Task 7.4: Update CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Update the "Source Layout" + "Key Architecture Decisions" sections**

- [ ] **Step 2: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Stage 7.4: update CLAUDE.md to reflect new struct architecture"
```

### Task 7.5: Final closeout — sanity sweep + PR preparation

- [ ] **Step 1: Full test suite — clean run**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: ALL tests pass; `test_rate_equation_performance` 0-alloc/<100ns gate green; ALL 5 trace-compile + 4 wall-clock compile-budget gates green.

- [ ] **Step 2: Final test-integrity verification (NON-NEGOTIABLE per spec §2)**

This is the FINAL check before opening the PR. The full no-deletion-since-main rule applies across all 7 stages:

```bash
bash scripts/check_test_integrity.sh main
```

Spec §2: across the entire refactor, no test may have been deleted, commented out, weakened, or `@test_skip`/`@test_broken`-tagged. Only mechanical syntax adaptation for new structs is permitted. If the script reports FAIL at this final stage: **STOP — DO NOT OPEN PR**. Audit every commit on the branch; restore the deleted/weakened tests; re-verify; only then proceed.

Also confirm the @testset count grew (you added new-type tests in Stage 1 plus the canonical-hash partition test in Stage 6):

```bash
n_main=$(git show main:test/runtests.jl 2>/dev/null | grep -c "@testset" || true)
for f in $(git ls-tree -r main --name-only test/ | grep "\.jl$"); do
    c=$(git show "main:$f" 2>/dev/null | grep -c "@testset" || true)
    n_main=$((n_main + c))
done
n_head=$(grep -rh "@testset" test/ | grep -c "@testset")
echo "@testset count: main=$n_main, HEAD=$n_head"
test "$n_head" -gt "$n_main" || { echo "FAIL: no new tests added across the entire refactor — suspicious"; exit 1; }
```

- [ ] **Step 3: Final LOC verification**

```bash
wc -l src/*.jl
```

Expected: total ≤ 3,600. If not, return to Task 7.2 for more aggressive cleanup or revise expected target.

- [ ] **Step 4: Generate PR description**

Based on `git log main..HEAD --oneline`, draft PR description covering:
- Motivation (spec link)
- Summary of changes (per-stage bullet points)
- Behavior changes (DSL grammar, structural kinetic groups)
- Migration path for users
- Perf gates passed
- LOC reduction achieved

- [ ] **Step 5: Push branch + open PR**

```bash
git push -u origin refactor-to-concrete-types-instead-of-symbols
gh pr create --title "Refactor to concrete types instead of Symbols" \
  --body "$(cat <<'EOF'
## Summary
- Replaced Symbol-string-juggling internals with one concrete struct family.
- Unified enumeration + derivation data structures (Step/Species/Parameter shared).
- src LOC reduced from 7,136 → <ACTUAL>, hitting the ≥50% reduction goal.

## Per-stage changes
[fill in from git log]

## Behavior changes
- `@enzyme_reaction`: per-regulator multiplicities; `competitive_inhibitors:` rename
- `@enzyme_mechanism`: bare metabolite names; parenthesized step groups; function-call species notation
- `@allosteric_mechanism`: `regulatory_site(multiplicity=N): begin ligands: ... end` blocks
- Parameter naming preserved (positional `:K1`, `:k1f`) via chokepoint `name(p, m)` accessor for future migration

## Perf gates passed
- `rate_equation` 0-alloc/<100ns per call: GREEN
- `loss!` runtime: no regression vs main baseline
- Compile-time budget (5 trace-compile + 4 wall-clock gates): GREEN

## Test integrity (spec §2 NON-NEGOTIABLE)
- [x] `bash scripts/check_test_integrity.sh main` PASSES
- [x] No test file deleted
- [x] @testset count grew vs main (new-type tests + canonical-hash partition test)
- [x] No `@test_skip` / `@test_broken` introduced
- [x] No `@test` lines commented out

## Test plan
- [x] Full test suite green
- [x] All perf gates passing

Spec: docs/superpowers/specs/2026-05-20-concrete-types-refactor-design.md
EOF
)"
```

- [ ] **Step 6: Tag final**

```bash
git tag refactor-complete
```

---

# Plan complete

Cumulative target met:
- src LOC: 7,136 → ≤3,600 (≥50% reduction)
- All tests preserved + adapted mechanically
- All perf gates green
- One large PR opened on `refactor-to-concrete-types-instead-of-symbols`
