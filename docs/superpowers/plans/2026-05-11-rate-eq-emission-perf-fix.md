# Rate equation emission perf fix: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore zero-allocation, sub-100-ns runtime for `rate_equation` by rebuilding the `+`/`*` Expr trees emitted by `_nest_binary` as balanced binary instead of flat n-ary. Lock the contract in with tightened perf-test thresholds, an Expr-shape regression test, a flat-printed-string regression test, and a CLAUDE.md non-negotiable.

**Architecture:** One TDD cycle fixes the bug (test threshold + `_nest_binary` body, committed together). Three regression-prevention commits follow — Expr shape, paren count, CLAUDE.md rule. Status doc closure note finishes. Five commits total.

**Tech Stack:** Julia 1.9+, EnzymeRates.jl package at `/home/denis.linux/.julia/dev/EnzymeRates`. Tests use Julia's `Test`.

**Reference spec:** `docs/superpowers/specs/2026-05-11-rate-eq-emission-perf-fix-design.md`.

---

## Pre-flight

### TDD development loop

Cold `Pkg.test()` runs are slow because of precompilation. Use a long-lived Julia REPL with Revise.jl for fast iteration:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project
```

In the REPL:

```julia
using Revise
using EnzymeRates
using Test
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
```

After source edits, re-`include` the affected test file. For the perf-test sweep specifically:

```julia
include("test/test_rate_eq_derivation.jl")
```

For final validation use `Pkg.test()`:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

### Witness mechanism

All regression-prevention tests use **Random-order Bi-Bi** as the single witness. It has 96 denominator monomials — the widest margin past the n-ary inlining cliff in the current test suite, so any regression of `_nest_binary` will show up here. Defined in `test/mechanism_definitions_for_test_enzyme_derivation.jl` with `name="Random-order Bi-Bi"`.

### Pre-fix baseline measurements

These are the numbers you should see at branch HEAD before Task 1. If they don't match, stop and investigate — something else has shifted.

| Mechanism             | den terms | allocs/call | ns/call (approx) |
|-----------------------|-----------|-------------|------------------|
| Segel Ordered Ter Bi  | 25        | 64 B        | 130 ns           |
| Segel Ordered Ter Ter | 36        | 960 B       | 575 ns           |
| RE Random Bi-Bi       | 42        | 1088 B      | 670 ns           |
| Random-order Bi-Bi    | 96        | 2432 B      | 1305 ns          |

Witness paren counts (current state, `_expr_to_string` already flattens):

- `rate_equation_string(m, Full)` → 5 open-parens
- `rate_equation_string(m, Reduced)` → 13 open-parens

The fix does not change these numbers — `_expr_to_string` already strips parens from balanced `+`/`*` nodes. The paren-count test guards against a future regression of that flattening behavior.

---

## Task 1: Fix `_nest_binary` + tighten perf-test thresholds

**One TDD cycle, two files committed together.** The threshold change makes the perf test fail; the `_nest_binary` change makes it pass.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl:542-547`
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:123-126`

### Step 1: Tighten the perf-test thresholds (failing test)

Replace lines 542-547 of `test/test_rate_eq_derivation.jl`. Strip the apologetic comment and snap thresholds back to the zero-allocation contract.

Before:

```julia
        allocs, t = test_rate_equation_performance(m, params, concs)
        # Rate equations are flat polynomial sums (up to ~100 monomials for
        # ter-ter mechanisms). Julia's optimizer sometimes spills these to
        # the heap; the bounds keep evaluations fast enough for fitting
        # loops while accommodating worst-case mechanisms.
        @test allocs < 4 * 1024
        @test t < 10e-6
```

After:

```julia
        allocs, t = test_rate_equation_performance(m, params, concs)
        @test allocs == 0
        @test t < 100e-9
```

- [ ] **Step 1a: Apply the threshold change.**

Use Edit on `test/test_rate_eq_derivation.jl:542-547` to replace the block above.

- [ ] **Step 1b: Run the perf sweep — confirm it fails.**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -60
```

Expected: `Pkg.test()` reports test failures in the `Performance` testset for every mid-and-large mechanism (Segel Ordered Ter Ter, RE Random Bi-Bi, Random-order Bi-Bi, …). Small mechanisms (≤25 den terms) still pass. This is the failing-test state — proceed to Step 2.

If the perf test passes anyway, stop — the bug already fixed itself somehow, and the rest of this plan is unnecessary. Investigate before continuing.

### Step 2: Apply the `_nest_binary` fix

Replace `_nest_binary` in `src/sym_poly_for_rate_eq_derivation.jl:123-126`.

Before:

```julia
"""Build flat n-ary expression: +(a, b, c, d). Single term returns unwrapped."""
function _nest_binary(op::Symbol, terms::Vector{Any})
    length(terms) == 1 ? terms[1] : Expr(:call, op, terms...)
end
```

After:

```julia
"""
Build a balanced binary `+`/`*` tree so every emitted call has exactly two
operands. Required for zero-allocation `rate_equation` runtime: Julia inlines
binary `+(::Float64, ::Float64)` into fused scalar arithmetic, but falls back
to a varargs path that boxes the operand tuple once the chain exceeds ~30
terms. See `test_rate_equation_performance` for the contract this enforces.
"""
function _nest_binary(op::Symbol, terms::Vector{Any})
    n = length(terms)
    n == 1 && return terms[1]
    n == 2 && return Expr(:call, op, terms[1], terms[2])
    mid = n >> 1
    Expr(:call, op,
        _nest_binary(op, terms[1:mid]),
        _nest_binary(op, terms[mid+1:end]))
end
```

- [ ] **Step 2a: Apply the source change.**

Use Edit on `src/sym_poly_for_rate_eq_derivation.jl:123-126` to replace the block.

- [ ] **Step 2b: Run the full test suite — confirm it now passes.**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

Expected: all testsets green, including `Performance` for every mechanism in `MECHANISM_TEST_SPECS`. Numerical tests (`Reference QSSA`, `Analytical Rate`, `Haldane`, `ODE Steady-State`, `Kcat …`) pass byte-identical to before — the change reorders Float64 sums by at most ~5000 ulps total, well below the tightest `rtol=1e-10` used anywhere in the suite.

If any non-perf test fails:
- A failing numerical test with `rtol=1e-10` and a delta below ~1e-12 means the reordering is responsible — relax the test's `rtol` only after confirming the underlying values differ by less than 1 ulp per summand. Document the rationale in the test comment.
- Any other failure is a real bug; stop and investigate.

If any `Performance` test fails:
- Find the failing mechanism. Run `_build_rate_body(typeof(m), Reduced)` in the REPL and walk it looking for `Expr(:call, op, args…)` with `op ∈ (:+, :*)` and `length(args) > 3`. Such a node bypasses `_nest_binary` and needs its own balancing pass. The most likely suspect is `_kcat_forward`'s `Expr(:call, :max, corner_exprs...)` — but `_kcat_forward` is not in `rate_equation`'s body, so this should not happen for the perf test. If it does, spec § "Out of scope" needs revision before fixing here.

### Step 3: Commit

- [ ] **Step 3a: Commit both files together.**

```bash
git add src/sym_poly_for_rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Fix rate_equation alloc regression: balanced binary _nest_binary

Emit balanced binary +/* trees so every call has exactly 2 operands.
Julia inlines binary Float64 +/* into fused scalar arithmetic; flat
n-ary above ~30 terms falls back to a varargs path that boxes the
operand tuple, turning 100ns/0B into 1µs/2KB per rate_equation call.

Same fix as commit 598a964 (Feb 2026), which was reverted in dec0707
(Mar) when factoring kept individual +/* chains short. Commit 28041fe
on this branch dropped factoring and re-exposed the codegen pathology.

Restore the zero-allocation, sub-100ns Performance thresholds.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Regression-prevention test — Expr shape

Adds a standalone test that walks the body of `_build_rate_body(typeof(m), Reduced)` for the witness mechanism and asserts every `+`/`*` call node has exactly two arguments. This is a guard test, not a feature test — it can't be written failing-first because Task 1 already made it pass. Per CLAUDE.md TDD rules this is permitted because the bugfix already had its own failing test (the perf assertions in Task 1).

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (append a new `@testset` after the main test loop, before line 754 `@testset "_is_ss_rate_constant"`)

### Step 1: Add the testset

- [ ] **Step 1a: Add the test.**

Scope the walk to the rate expression returned by `_raw_rate_expr_and_symbols` — that is exactly what `_nest_binary` controls. The full `_build_rate_body(ReducedMode)` block also contains dep-expr assignments emitted by `build_power_expr` (a separate codepath, `src/sym_poly_for_rate_eq_derivation.jl:230`) which uses small flat n-ary `*` (7-10 args) — those are harmless at runtime (the perf sweep proves it: the n-ary inlining cliff is at ~30, not 10) and out of scope for this test.

Insert this block immediately before line 754 (`@testset "_is_ss_rate_constant" begin`):

```julia
@testset "rate_equation polynomial body uses 2-arg +/* calls" begin
    # The fitter calls rate_equation millions of times per CV fold. The
    # polynomial body emitted by _poly_to_expr (via _nest_binary) MUST
    # have exactly 2 operands per +/* call so LLVM inlines the binary
    # Float64 path; n-ary varargs above ~30 terms boxes the argument
    # tuple and turns 100ns/0B into 1µs/2KB per call. See
    # docs/superpowers/specs/2026-05-11-rate-eq-emission-perf-fix-design.md.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    rate_expr, _, _ = EnzymeRates._raw_rate_expr_and_symbols(
        typeof(spec.mechanism))
    bad = Expr[]
    function walk!(e)
        if e isa Expr
            if e.head == :call && !isempty(e.args) &&
               e.args[1] isa Symbol && e.args[1] in (:+, :*) &&
               length(e.args) != 3
                push!(bad, e)
            end
            start = e.head == :call ? 2 : 1
            for i in start:length(e.args)
                walk!(e.args[i])
            end
        end
    end
    walk!(rate_expr)
    @test isempty(bad)
end
```

The `length(e.args) != 3` check: `e.args[1]` is the operator symbol, so a binary call has `length(e.args) == 3`.

### Step 2: Run and commit

- [ ] **Step 2a: Run the test, confirm pass.**

In the REPL:

```julia
include("test/test_rate_eq_derivation.jl")
```

Or from the shell, the targeted subset:

```bash
julia --project -e 'using Test; include("test/test_rate_eq_derivation.jl")' 2>&1 | grep -E "balanced binary|Test Summary" | head
```

Expected: `Test Summary: rate_equation Expr is balanced binary (+/* calls have 2 args) | Pass: 1 Total: 1`.

- [ ] **Step 2b: Commit.**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Regression test: rate_equation Expr body uses 2-arg +/* calls

Walks _build_rate_body output for the largest non-allosteric witness
(Random-order Bi-Bi, 96 den terms) and asserts every +/* call node
has exactly two operands. Guards against a future "simplification"
of _nest_binary back to flat n-ary slipping through unnoticed.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Regression-prevention test — flat printed string

Adds a paren-count assertion on `rate_equation_string(m, Full)` for the witness mechanism. Guards against a regression in `_expr_to_string`'s `+`/`*` flattening — without this guard, a printed body could explode into nested parens while still numerically correct.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (append a second `@testset` right after the Task 2 testset)

### Step 1: Add the testset

The current paren count for the witness in `Full` mode is exactly 5: `(; …) = params`, `(; …) = concs`, `(num)`, `(den)`, and one `(negatives)` wrap inside the numerator. Cap at 6 — strict enough that introducing one extra nested-`+` paren pair fails the test.

- [ ] **Step 1a: Add the test** immediately after the Task 2 testset, before line 754 `@testset "_is_ss_rate_constant"`:

```julia
@testset "rate_equation_string prints flat +/* sums (no nested parens)" begin
    # _expr_to_string is precedence-aware and flattens nested +/* nodes
    # transparently: balanced +(+(a,b), +(c,d)) prints as "a + b + c + d"
    # with no added parens. If that flattening regresses, the printed
    # output would gain nested parens but still evaluate correctly —
    # easy to miss without an explicit guard. Witness count for
    # Random-order Bi-Bi is 5 in Full mode (params destructure, concs
    # destructure, num wrap, den wrap, one negatives wrap inside num).
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    s = rate_equation_string(spec.mechanism, EnzymeRates.Full)
    @test count(==('('), s) <= 6
end
```

### Step 2: Run and commit

- [ ] **Step 2a: Run the test, confirm pass.**

```bash
julia --project -e 'using Test; include("test/test_rate_eq_derivation.jl")' 2>&1 | grep -E "flat.+sums|Test Summary" | head
```

Expected: `Test Summary: rate_equation_string prints flat +/* sums (no nested parens) | Pass: 1 Total: 1`.

- [ ] **Step 2b: Commit.**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Regression test: rate_equation_string keeps +/* chains flat-printed

Asserts that rate_equation_string(m, Full) on the largest witness
(Random-order Bi-Bi) has ≤6 open-parens — the structural ceiling
(2 destructurings + num wrap + den wrap + one negatives wrap).
Guards against _expr_to_string losing its +/* flattening behavior;
without this, a balanced-binary Expr body could explode into nested
parens in the printed form while still numerically correct.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add `rate_equation` perf rule to CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md` (append to the project-specific Testing section near the end of the file)

### Step 1: Find the insertion point

The project-specific CLAUDE.md ends with a `## Testing` section. Find the last line of that section — it should be the bullet starting `- kcat/rescaling tests …`. Insert the new rule as a new subsection immediately after the existing Testing bullets, before the end of the file.

- [ ] **Step 1a: Find the anchor line.**

```bash
grep -n "^## Testing\|^### \|kcat/rescaling tests" .claude/CLAUDE.md | tail -10
```

The bullet `- kcat/rescaling tests (scale invariance, rate proportionality, V≈1, custom target) run for ALL mechanism specs in the main \`run_all_tests\` loop — not in a separate file` is the last bullet of the section.

### Step 2: Insert the rule

- [ ] **Step 2a: Use Edit to add the new subsection** after that final bullet:

The new block to insert:

```markdown

### `rate_equation` runtime perf is non-negotiable

`rate_equation` MUST be allocation-free and sub-100-ns per call for every mechanism in `MECHANISM_TEST_SPECS`. Enforced by `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl` (`allocs == 0`, `t < 100e-9`) plus the Expr-shape and flat-string regression tests in the same file. The fitter evaluates `rate_equation` millions of times per cross-validation fold; any change that introduces allocations or microsecond-scale per-call time makes the package unusable in practice. If you find yourself considering a change that would force `rate_equation` to allocate or slow down — STOP and discuss with Denis first. This is one of the most important tests in the suite.
```

Use Edit with `old_string` being the existing final bullet (verbatim, including its leading `- ` and trailing newline if any) and `new_string` being the same bullet plus the block above.

### Step 3: Commit

- [ ] **Step 3a: Commit.**

```bash
git add .claude/CLAUDE.md
git commit -m "$(cat <<'EOF'
CLAUDE.md: rate_equation runtime perf is non-negotiable

Document the zero-allocation, sub-100ns contract as a project-level
rule. The fitter calls rate_equation millions of times per CV fold;
any regression that introduces allocations or microsecond-scale
per-call time makes the package unusable for parameter identification.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Closure note in `dedup_investigation/status_2026-05-11.md`

**Files:**
- Modify: `dedup_investigation/status_2026-05-11.md` (append to end)

### Step 1: Append the closure note

- [ ] **Step 1a: Read current end of file** to get the exact final lines for Edit:

```bash
tail -10 dedup_investigation/status_2026-05-11.md
```

The file currently ends with the `## Files of interest` section. Use Edit to append a new section after the last line (after the `src/identify_rate_equation.jl` bullet about the canonicalizer).

- [ ] **Step 1b: Append this block** to the end of the file:

```markdown

## Update — issue (1) closed

Root cause of the rate_equation alloc/slowdown regression: `_nest_binary` in `src/sym_poly_for_rate_eq_derivation.jl` was emitting flat n-ary `+(a, b, …, a_N)` Expr nodes. Above ~30 terms Julia's compiler stops inlining the varargs `+` path and falls back to a heap-boxed-argument-tuple codegen — turning the per-call cost from 100 ns / 0 B into ~1 µs / ~2 KB. Same regression as commit `598a964` (Feb), reverted in `dec0707` (Mar), re-exposed by commit `28041fe` on this branch when factoring was dropped. Fixed by restoring the balanced-binary emission. See `docs/superpowers/specs/2026-05-11-rate-eq-emission-perf-fix-design.md` and `docs/superpowers/plans/2026-05-11-rate-eq-emission-perf-fix.md`.

Issue (2) — residual same-hash-different-loss / same-loss-different-hash within `params_estimate_*.csv` — is deferred to a follow-up. The hypothesis in this status doc ("CMA-ES `maxtime=60` marginal under the slower per-call cost") is now testable: re-fit a same-hash mech pair against `LDH_data.csv` with the perf fix applied, at the default `maxtime`. If the within-hash spread collapses to fitter noise, the c8a3302 `fitted_params` filter fix is the only remaining substantive change and the dedup work is done. If anything residual remains, dump canonical strings for two non-collapsing mechanisms and diff to identify the next quirk.
```

### Step 2: Commit

- [ ] **Step 2a: Commit.**

```bash
git add dedup_investigation/status_2026-05-11.md
git commit -m "$(cat <<'EOF'
status doc: close issue (1), record issue (2) followup plan

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final validation

After all five tasks are committed:

- [ ] **Run the full test suite cold:**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: all green. The full sweep covers every spec in `MECHANISM_TEST_SPECS` with the tightened `allocs == 0` / `t < 100e-9` thresholds plus the two new regression-prevention testsets and the existing Aqua / JET quality gates.

- [ ] **Verify clean git state:**

```bash
git status && git log --oneline -7
```

Expected: clean working tree, 5 new commits at the tip of `improve-deduplication-of-ident-equations`.

- [ ] **Spot-check the printed body.**

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
spec = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi")
println(rate_equation_string(spec.mechanism, EnzymeRates.Full))
'
```

Eye-check: the `v = …` line should be a flat sum `a + b + c + …` with no nested parens around the `+` chain. The only parens should be: two `(; …)` destructurings, `(num)`/`(den)` outer wraps, and one `(neg_1 + … + neg_k)` wrap inside the numerator if there are negative monomials.
