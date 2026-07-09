# Equation-Complexity Filter — Design

**Goal:** Skip mechanisms whose rate equations are too complex to fit in practical time,
before any derivation runs — bounding fit cost and, we expect, preventing the codegen
segfault that killed the inhibitor run.

**Architecture:** One additive change to the beam search in `src/identify_rate_equation.jl`,
plus a complexity metric computed from the mechanism graph. It never touches the
`rate_equation` runtime path, so the 0-allocation / sub-120 ns contract is unaffected.

## Background

The LDH run with all four substrates/products also declared as competitive inhibitors
(`docs/ldh_hpc_results/2026_07_09_results`) crashed at iteration 6; the identical run
without inhibitors (`2026_07_08_results`) completed. A worker segfaulted compiling a giant
inhibitor `@generated rate_equation` body, and the uncaught `ProcessExitedException` killed
the whole distributed run.

The root cause is unbounded equation complexity: competitive inhibitors multiply enzyme
forms, producing rate equations with orders of magnitude more terms. `max_param_count` does
not bound this — parameter count and equation size are only loosely related. Such equations
are impractical to fit (µs-scale per call → hours per equation) and, at the extreme, blow up
derivation/codegen.

This spec adds a complexity filter that skips those mechanisms before deriving anything. We
expect it to remove the segfault at its source (the over-complex equations are never
compiled or fit).

## The metric — V×τ

The rate-equation denominator is the King–Altman enzyme-conservation sum. Its structure is
governed by the **segment graph**: rapid-equilibrium (RE) steps collapse the forms they
connect into a single node (a *segment*), and steady-state (SS) steps are the edges between
segments. Two quantities:

- **V** = number of segments (RE-connected components of the form graph).
- **τ** = number of spanning trees of the segment graph (Matrix–Tree theorem: the
  determinant of any cofactor of the graph Laplacian).

**V×τ is exactly the number of spanning-tree products in the denominator** — each of the V
segments contributes τ products, each a product of V−1 rate constants and the bound-metabolite
concentrations. This is the count of products the compiled equation evaluates on every call,
so it is the fit-cost driver. (The number of *distinct concentration monomials* — what one
would call the "terms" of the written equation — is smaller and does **not** predict cost:
ordered ter-ter has more monomials than random-order bi-bi yet is far cheaper.)

Verified on canonical all-SS mechanisms (exact denominator-product count matches V×τ):

| mechanism        | V | τ    | V×τ  | conc-monomials | per-call time |
|------------------|---|------|------|----------------|---------------|
| ordered bi-bi    | 5 | 5    | 25   | 11             | ~16 ns        |
| random-order bi-bi | 7 | 48 | 336  | 16             | ~125 ns       |
| ordered ter-ter  | 7 | 7    | 49   | 27             | ~25 ns        |
| random-order ter-ter | 15 | 393 216 | 5 898 240 | — | derivation blows up |

Connectivity, not size, dominates: ordered ter-ter (V×τ = 49) is cheaper than random-order
bi-bi (V×τ = 336) despite being a larger reaction. Allostericity is **not** a driver — the
`^N` multiplicity does not change the catalytic King–Altman complexity — so V×τ is computed
on the catalytic core. Each spanning-tree product multiplies exactly V−1 rate constants (a
spanning tree has V−1 edges) plus up to reaction-order concentrations; that depth is bounded
and folded into V, so V×τ alone is the guard metric.

## The threshold

Default `eq_complexity_filter = 336`, the V×τ of a random-order bi-bi. Rule of thumb: an
equation more complex than a random-order bi-bi is not worth the fitting effort. The next
tier up (random-order ter-ter, ~5.9 M) is a 17 000× gap, so the exact value is not delicate;
the cut lands in a wide empty band between "worth fitting" and "blows up."

## Computation

```julia
# V×τ of the catalytic segment graph. Pure struct walk + a small exact determinant;
# no compile, no fitted_params, no King–Altman derivation.
function _eq_complexity(m)
    cm = m isa AllostericMechanism ? catalytic_mechanism(m) : m
    # 1. union-find forms over RE edges (is_equilibrium == true) -> segments
    # 2. build the segment-graph Laplacian from SS edges (is_equilibrium == false)
    #    between distinct segments
    # 3. V = number of segments; τ = exact integer cofactor determinant (Matrix–Tree)
    return V * τ
end
```

- Compute on the **raw** mechanism (`Mechanism` / `AllostericMechanism`), before
  `compile_mechanism`/`fitted_params`, because those can themselves blow up on a monstrous
  mechanism.
- τ is an integer; compute the Matrix–Tree determinant with **exact** arithmetic
  (`Rational{BigInt}` or an integer Bareiss determinant) — τ can exceed `Int64` on dense
  graphs, and floating point would round.
- The filter is self-limiting: skipping an over-complex mechanism keeps it out of the beam
  frontier, so its children are never generated, so V stays small and the determinant stays
  cheap. Defensive guard: if V exceeds a generous bound (e.g. 40), skip immediately without
  the determinant.

## Placement and plumbing

- New keyword `eq_complexity_filter::Int = 336` on `identify_rate_equation`, threaded to
  `_beam_search` and `_process_batch` exactly as `max_param_count` is.
- Checked as the **first** line of the PASS-1 closure, before `compile_mechanism`. A cap
  skip identical in shape to the `max_param_count` cap, but earlier because it must precede
  the derivation-heavy calls.

## Skip tracking and reporting

The complexity skips are counted and printed in the progress log the same way
`max_param_count` skips already are. Today both a param-count skip and (with this change) a
complexity skip would return `nothing` from PASS-1, and `n_skipped` is derived as
`length(batch) − length(entries) − length(failures)` — so the reason is lost and everything
is mislabeled ">max_param_count params". The fix distinguishes them:

- PASS-1 returns a **distinct sentinel** for a complexity skip (e.g. `:complexity_skip`)
  versus `nothing` for a param-count skip. Downstream, both are treated as skips (dropped
  from the row loop), but they are counted separately.
- `_process_batch` returns the breakdown — `(entries, failures, n_param_skipped,
  n_complexity_skipped)` — instead of the callers re-deriving a single `n_skipped` from
  lengths.
- `_batch_summary` takes both counts and both thresholds, and prints both buckets, e.g.
  `… + N skipped (>20 params) + M skipped (>336 complexity) + …`. The buckets still sum to
  the child count.
- The "whole batch skipped" branch (`identify_rate_equation.jl:~735`) likewise reports the
  two skip reasons separately.

## Documentation

Add an `eq_complexity_filter` entry to the `identify_rate_equation` docstring, next to
`max_param_count`: what V×τ measures (spanning-tree product count of the catalytic segment
graph = number of products evaluated per call), the default (336 = random-order bi-bi), and
why it exists (skip equations too complex to fit in practical time; expected to prevent the
derivation-blow-up crash).

## Tests

- **Unit** — `_eq_complexity` returns the known V×τ for reference mechanisms: ordered bi-bi
  = 25, random-order bi-bi = 336, ordered ter-ter = 49 (all constructible via
  `@enzyme_mechanism`).
- **Filter** — `_process_batch` with a low `eq_complexity_filter` skips a mechanism above it
  (absent from `entries`, not fit) and keeps one below it.
- **Skip counts** — `_process_batch` reports the complexity skip in its own bucket
  (`n_complexity_skipped`), separate from `n_param_skipped`, and `_batch_summary` renders
  both; a mechanism dropped for complexity is not miscounted as a param-count skip.
- **Default** — `identify_rate_equation`'s default `eq_complexity_filter` is 336.

TDD per change (failing test → implement → green). Full `Pkg.test()` green gate before
finishing, including the `rate_equation` performance regression tests (unaffected — the
change is pre-fit control flow).

## Deferred: worker-death fault tolerance

A second fix — making the fitting `pmap` survive a worker crash — was scoped and then
deferred, because the complexity filter is expected to remove the crash at its source. Recorded
here so the investigation is not repeated if segfaults persist:

- `pmap`'s `on_error` does **not** catch a worker *process death*; the `ProcessExitedException`
  (connection EOF) propagates and crashes the run (verified). `retry_delays` + `retry_check`
  retries per-task but still throws and cascades on a deterministic death (verified — pool
  8→4 workers). This is a known, unresolved `pmap` limitation
  (JuliaLang/julia#36709, #44465).
- A dead worker is auto-removed from `workers()`; it is not restarted, and subsequent `pmap`s
  run on the survivors.
- Doing it right (preserve completed fits *and* skip the equation that killed the worker)
  requires a manual master–worker dispatch loop that tracks `worker → equation`, so a death
  marks that one equation errored (optionally after one retry) and continues on the survivors —
  roughly 50 lines replacing the two `pmap` calls in `_process_batch`.
- Revisit only if segfaults survive the complexity filter.
