# Rate equation emission: restore zero-allocation runtime

## Problem

On the `improve-deduplication-of-ident-equations` branch, `rate_equation`
calls allocate up to several KB per call and run in the microsecond range
for mid-sized mechanisms — versus the prior zero-allocation, sub-100-ns
target that makes `identify_rate_equation` practical for cross-validation
fitting. The `test_rate_equation_performance` thresholds in
`test/test_rate_eq_derivation.jl:546-547` were relaxed from
`allocs == 0` / `t < 100 ns` to `allocs < 4 KiB` / `t < 10 µs` in commit
`28041fe` to keep the test suite passing through the relaxation.

Empirical, current branch, real mechanisms:

| Mechanism             | den terms | allocs/call | ns/call |
|-----------------------|-----------|-------------|---------|
| Segel Ordered Ter Bi  | 25        | 64 B        | 130 ns  |
| Segel Ordered Ter Ter | 36        | 960 B       | 574 ns  |
| RE Random Bi-Bi       | 42        | 1088 B      | 670 ns  |
| Random-order Bi-Bi    | 96        | 2432 B      | 1305 ns |

The discontinuity between 25 and 36 terms is the LLVM-inlining boundary
for n-ary `+(::Float64, ::Float64, …)` codegen.

## Root cause

`_nest_binary` (`src/sym_poly_for_rate_eq_derivation.jl:123-126`) emits
the polynomial body as a single flat n-ary node:

```julia
function _nest_binary(op::Symbol, terms::Vector{Any})
    length(terms) == 1 ? terms[1] : Expr(:call, op, terms...)
end
```

When `terms` has more than ~30 entries, Julia's compiler stops inlining
the varargs `+`/`*` form and falls back to a path that allocates a
heap-boxed argument tuple. This had a known fix already shipped in
February (`598a964`, "Fix rate_equation StackOverflow"): emit a balanced
binary tree so every `+`/`*` call has exactly two operands. That fix
was reverted in commit `dec0707` (March, "Implement derived-eq
factorization") when the factoring pipeline kept individual `+` chains
short enough that flat n-ary was harmless. Commit `28041fe` on this
branch dropped factoring entirely and now emits the whole expanded
polynomial through one flat `_nest_binary(:+, …)` of up to several
thousand monomials — re-exposing the original codegen pathology.

Synthetic A/B over 36 identical monomials confirms the diagnosis:

```
flat  +(m1, …, m36):       880 B/call,  1264 ns
balanced +(+(…), +(…)):      0 B/call,    82 ns
```

## Approach

Restore the balanced-binary-tree implementation of `_nest_binary`. The
change is local — `_nest_binary` is the single helper through which
every `+`/`*` chain in every emission path passes:

- non-allosteric numerator / denominator (`_poly_to_expr`)
- allosteric R/T-state numerator / denominator
- regulatory-site saturation factors in `_kcat_forward`
- W-factor products in `_kcat_forward`
- power-product factors in `_build_allosteric_rate_body`

All callers benefit automatically. Print output is unchanged because
`_expr_to_string` already flattens nested `+`/`*` nodes transparently
(verified — balanced `+(+(a,b), +(c,d))` prints as `"a + b + c + d"`
with no added parentheses).

Approach was selected over (B) common-subexpression elimination in the
emitter and (C) reintroducing runtime-only factoring: (A) is a 3-line
behavior change, demonstrably hits the prior zero-allocation target,
and avoids re-introducing the factoring machinery that commit `28041fe`
deleted.

## Changes

### 1. `src/sym_poly_for_rate_eq_derivation.jl`

Replace `_nest_binary` body with the balanced-tree form from commit
`598a964`. Update docstring to record the codegen reason:

```julia
"""
Build a balanced binary `+`/`*` tree so every call has exactly two
operands. Required for zero-allocation runtime: Julia inlines binary
`+(::Float64, ::Float64)` into fused scalar arithmetic, but falls
back to a varargs path that boxes the operand list once the chain
exceeds ~30 terms. See `test_rate_equation_performance` for the
contract this enforces.
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

### 2. `test/test_rate_eq_derivation.jl`

Re-tighten the perf-test thresholds (lines 542-547). Remove the
apologetic comment; the budget is again real.

```julia
allocs, t = test_rate_equation_performance(m, params, concs)
@test allocs == 0
@test t < 100e-9
```

Add two regression-prevention tests in the same file, scoped to the
largest non-allosteric mechanism in `MECHANISM_TEST_SPECS` (currently
Random-order Bi-Bi, 96 den terms — picked because it crosses the
n-ary inlining cliff by the widest margin):

a. **Expr-shape assertion.** Walk the rate expression returned by
   `_raw_rate_expr_and_symbols(typeof(m))[1]` and assert every
   `Expr(:call, op, args…)` with `op ∈ (:+, :*)` has exactly two
   non-head args. This is exactly what `_nest_binary` controls,
   and what prevents a future "simplification" back to flat n-ary
   from slipping through unnoticed. Scope is narrowed to the rate
   expression — not the dep-expr assignments built by
   `build_power_expr` (a separate codepath using small flat n-ary
   `*` of 7-10 args, which the perf sweep proves does not allocate
   at runtime since the LLVM inlining cliff is at ~30 terms).

b. **Flat-string assertion.** After the fix, the only parens in
   `rate_equation_string(m)` come from a small, fixed set of sources:
   the destructuring lines `(; …) = params` and `(; …) = concs` (2),
   the `E_total * (num) / (den)` template (2), and at most one `(ne)`
   wrap per polynomial that has negative monomials (≤ 2 for the
   num+den pair). That's a hard ceiling of 6 paren pairs (12 chars)
   for any non-allosteric mechanism. Assert
   `count('(', rate_equation_string(m, Full)) ≤ 6` on the
   witness mechanism. The current count is 5 (2 destructurings + num
   wrap + den wrap + one negatives wrap inside the numerator); the
   cap is one above current, so any extra nested-`+` paren pair fails
   the test. `Full` mode is the right witness because `Reduced` adds
   dep-expr `1/(…)` divisor assignments (~13 parens) — structural
   noise that would mask a flattening regression.

### 3. `.claude/CLAUDE.md`

Add a paragraph under the project-specific Testing section
(end-of-file area, after the existing testing rules) calling out the
rate-equation runtime perf contract as a non-negotiable. Draft:

> ### `rate_equation` runtime perf is non-negotiable
>
> `rate_equation` MUST be allocation-free and sub-100-ns per call for
> every mechanism in `MECHANISM_TEST_SPECS`. Enforced by
> `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`.
> The fitter evaluates `rate_equation` millions of times per
> cross-validation fold; any change that introduces allocations or
> microsecond-scale per-call time makes the package unusable in
> practice. If you find yourself considering a change that would
> force `rate_equation` to allocate or slow down, STOP and discuss
> with Denis first — this is one of the most important tests in the
> suite.

### 4. `dedup_investigation/status_2026-05-11.md`

Append a short closure note at the end pointing at the fix commit:
issue (1) closed; issue (2) re-evaluation deferred to a follow-up
that re-fits LDH end-to-end against the perf-fixed branch. Leave the
existing diagnostic notes untouched as historical record.

## Validation

1. **Synthetic A/B** (already done): flat 36-term sum → 880 B / 1264 ns;
   balanced 36-term sum → 0 B / 82 ns. ~15× faster, zero allocations.
2. **Full `MECHANISM_TEST_SPECS` sweep** under tightened thresholds —
   every mechanism (including the 96-term Random-order Bi-Bi) must
   meet `allocs == 0` and `t < 100 ns`. If any mechanism fails,
   that's a real signal — most likely an `Expr(:call, op, …)` site
   outside `_nest_binary` (e.g., `Expr(:call, :max, corner_exprs...)`
   in `_kcat_forward:881`) needs the same treatment.
3. **Full `Pkg.test()`**: confirms no numerical / semantic
   regression. Reordering a sum of ≤5000 Float64 terms (the
   `MAX_RATE_EQUATION_TERMS` ceiling) accumulates at most ~5000 ulps
   ≈ 1.1e-12 relative — comfortably below the tightest `rtol=1e-10`
   used anywhere in the suite. Non-perf tests pass numerically
   indistinguishable from pre-fix; the perf assertions themselves
   move from passing the relaxed thresholds (`< 4 KiB` / `< 10 µs`)
   to passing the tightened thresholds (`== 0` / `< 100 ns`).
4. **Visual spot-check** on `rate_equation_string` for one mid-sized
   and one large mechanism: confirm the printed body is flat
   `a + b + c + …` with no nested parentheses around the `+` chain.

## Out of scope

- Re-running `identify_rate_equation` end-to-end on LDH to recompute
  `params_estimate_*.csv` and recount issue-(2) dedup residuals.
  Handled by a follow-up after this fix lands. The status note in
  `dedup_investigation/status_2026-05-11.md` already captures the
  hypothesis that CMA-ES `maxtime=60` was marginal because of the
  perf regression; the right next step is to remeasure with the
  regression gone.
- Touching `Expr(:call, :max, corner_exprs...)` in
  `_kcat_forward:881` or any other n-ary site not on the hot path.
  `_kcat_forward` runs once per fit, not per evaluation, so it is
  not perf-critical. Leave it alone unless the validation sweep
  surfaces a problem.

## Risk

Minimal. Balanced binary trees and flat n-ary nodes compute the same
floating-point value modulo summation order. Julia's compiler does
not assume associativity of FP `+` (we do not pass `--math-mode=fast`),
so the order changes deterministically per term position but not by
more than 1 ulp per addition. No downstream test compares
`rate_equation` output to a reference at tighter than ~1e-10 relative
tolerance, so the reordering is invisible.

The change is local to `_nest_binary`. No callers change. No types
change. No serialized state depends on `Expr` shape.
