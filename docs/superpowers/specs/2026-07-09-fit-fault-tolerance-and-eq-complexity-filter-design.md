# Fault-Tolerant Fitting and Equation-Complexity Filter — Design

**Goal:** Let a distributed `identify_rate_equation` run survive a worker crash, and
skip mechanisms whose rate equations are too complex to fit in practical time.

**Architecture:** Two independent changes to `src/identify_rate_equation.jl`, plus a
complexity metric computed from the mechanism graph in `src/rate_eq_derivation.jl` (or a
new small helper). Both are additive; neither touches the `rate_equation` runtime path,
so the 0-allocation / sub-120 ns contract is unaffected.

## Background

The LDH run with all four substrates/products also declared as competitive inhibitors
(`docs/ldh_hpc_results/2026_07_09_results`) crashed at iteration 6; the identical run
without inhibitors (`2026_07_08_results`) completed. Two independent root causes:

1. **Uncaught worker death.** A worker segfaulted compiling a giant inhibitor
   `@generated rate_equation` body. The master's fitting `pmap` had no fault tolerance,
   so the resulting `ProcessExitedException` propagated uncaught and killed the whole
   run. (slurm-35574432.out:563 segfault, :1022 uncaught `ProcessExitedException`.)

2. **No complexity bound.** Competitive inhibitors multiply enzyme forms, producing rate
   equations with orders of magnitude more terms. `max_param_count` does not bound this —
   parameter count and equation size are only loosely related. Such equations are
   impractical to fit (µs-scale per call → hours per equation) and, at the extreme, blow
   up derivation/codegen.

The two fixes are complementary: Fix 2 prevents most crashes by never fitting the
pathological mechanisms; Fix 1 survives any worker death that still slips through.

**Not in scope:** reducing the inhibitor-driven mechanism explosion at enumeration time
(a possible Fix 3). Deferred by request.

## Fix 1: Fault-tolerant fitting `pmap`

### Problem

`_process_batch` runs two `pmap` passes over the worker pool:
- PASS-1 (`src/identify_rate_equation.jl:512`): compile + cap-check + render each mechanism.
- PASS-2 (`:537`): fit one representative per distinct equation.

Each pass wraps its work in a `try/catch` **inside** the worker closure, which converts
Julia exceptions into a recorded failure. But a worker *process death* (segfault, OOM
kill) is not an exception inside the closure — `pmap` itself raises `ProcessExitedException`
on the master. That is uncaught, so `_process_batch`, `_beam_search`, and the whole run
crash.

### Fix

Add `on_error` to both `pmap` calls. On any task error (including `ProcessExitedException`),
`on_error(e)` returns a `WorkerCrashed(msg)` sentinel in that task's result slot. Because
`pmap` preserves input order, the crashed mechanism is recovered positionally
(`result[i]` ↔ `input[i]`) and folded into the existing failure paths:

- PASS-1: a `WorkerCrashed` at index `i` becomes `FitFailure(mechs[i], msg)` — recorded as
  a failure, like a caught exception.
- PASS-2: a `WorkerCrashed` at index `i` sets `fit_error[reps[i].eq_hash] = msg` — the same
  path a thrown representative fit already uses.

`pmap` drops the dead worker and continues the remaining tasks on survivors, so the batch
completes and the run proceeds.

```julia
struct WorkerCrashed
    msg::String
end

# PASS-1
raw = pmap(mechs; on_error = e -> WorkerCrashed(_exc_string(e))) do m
    ...  # unchanged body
end
compiled = Union{Nothing,FitFailure,NamedTuple}[
    c isa WorkerCrashed ? FitFailure(mechs[i], c.msg) : c
    for (i, c) in enumerate(raw)]

# PASS-2
rep_fits = pmap(reps; on_error = e -> WorkerCrashed(_exc_string(e))) do r
    ...  # unchanged body
end
for (i, r) in enumerate(rep_fits)
    if r isa WorkerCrashed
        fit_error[reps[i].eq_hash] = r.msg
    elseif r.error === nothing
        memo[r.eq_hash] = r.fit
    else
        fit_error[r.eq_hash] = r.error
    end
end
```

### Test

Distributed integration test (new, isolated): `addprocs(2)`, `@everywhere using EnzymeRates`,
run `_process_batch` on a small batch of valid mechanisms, kill a worker mid-run with
`rmprocs(w; waitfor = 0)`, and assert `_process_batch` returns `(entries, failures)`
(does not throw) with the affected mechanism recorded as errored and the others fit. This
reproduces the failure (before the fix it throws `ProcessExitedException`) and verifies the
fix (after, it returns).

## Fix 2: Equation-complexity filter

### The metric — V×τ

The rate-equation denominator is the King–Altman enzyme-conservation sum. Its structure is
governed by the **segment graph**: rapid-equilibrium (RE) steps collapse the forms they
connect into a single node (a *segment*), and steady-state (SS) steps are the edges
between segments. Two quantities:

- **V** = number of segments (RE-connected components of the form graph).
- **τ** = number of spanning trees of the segment graph (Matrix–Tree theorem: the
  determinant of any cofactor of the graph Laplacian).

**V×τ is exactly the number of spanning-tree products in the denominator** — each of the V
segments contributes τ products, each a product of V−1 rate constants and the bound-metabolite
concentrations. This is the count of terms the compiled equation evaluates on every call, so
it is the fit-cost driver. (The number of *distinct concentration monomials* — what one would
call "terms" of the written equation — is smaller and does **not** predict cost.)

Verified on canonical all-SS mechanisms (exact denominator-product count matches V×τ):

| mechanism        | V | τ    | V×τ  | per-call time |
|------------------|---|------|------|---------------|
| ordered bi-bi    | 5 | 5    | 25   | ~16 ns        |
| random-order bi-bi | 7 | 48 | 336  | ~125 ns       |
| ordered ter-ter  | 7 | 7    | 49   | ~25 ns        |
| random-order ter-ter | 15 | 393 216 | 5 898 240 | derivation blows up |

Connectivity, not size, dominates: ordered ter-ter (V×τ = 49) is cheaper than random-order
bi-bi (V×τ = 336) despite being a larger reaction. Allostericity is **not** a driver — the
`^N` multiplicity does not change the catalytic King–Altman complexity — so V×τ is computed
on the catalytic core.

### The threshold

Default `eq_complexity_filter = 336`, the V×τ of a random-order bi-bi. Rule of thumb: an
equation more complex than a random-order bi-bi is not worth the fitting effort. The next
tier up (random-order ter-ter, ~5.9 M) is a 17 000× gap, so the exact value is not delicate;
the cut lands in a wide empty band between "worth fitting" and "blows up."

### Computation

```julia
# V×τ of the catalytic segment graph. Pure struct walk + a small exact determinant;
# no compile, no fitted_params, no King–Altman derivation.
function _eq_complexity(m)
    cm = m isa AllostericMechanism ? catalytic_mechanism(m) : m
    # 1. union-find forms over RE edges -> segments
    # 2. build segment-graph Laplacian from SS edges (between distinct segments)
    # 3. V = number of segments; τ = exact integer cofactor determinant
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

### Placement and plumbing

- New keyword `eq_complexity_filter::Int = 336` on `identify_rate_equation`, threaded to
  `_beam_search` and `_process_batch` exactly as `max_param_count` is.
- Checked as the **first** line of the PASS-1 closure, before `compile_mechanism`:
  `_eq_complexity(m) > eq_complexity_filter && return nothing`. A cap skip (returns
  `nothing`), identical in shape to the `max_param_count` cap, but earlier because it must
  precede the derivation-heavy calls.
- Report the count in `_batch_summary` alongside the param-count skips (e.g.
  "N skipped (>eq_complexity_filter complexity)").

### Documentation

Add an `eq_complexity_filter` entry to the `identify_rate_equation` docstring, next to
`max_param_count`: what V×τ measures (spanning-tree product count of the catalytic segment
graph = number of terms evaluated per call), the default (336 = random-order bi-bi), and why
it exists (skip equations too complex to fit in practical time; also prevents the
derivation-blow-up crash).

### Tests

- **Unit** — `_eq_complexity` returns the known V×τ for reference mechanisms: ordered bi-bi
  = 25, random-order bi-bi = 336, ordered ter-ter = 49 (all constructible via
  `@enzyme_mechanism`).
- **Filter** — `_process_batch` with a low `eq_complexity_filter` skips a mechanism above it
  (absent from `entries`, not fit) and keeps one below it.
- **Default** — `identify_rate_equation`'s default `eq_complexity_filter` is 336.

## Testing strategy

TDD per fix (failing test → implement → green). Full `Pkg.test()` green gate before
finishing, including the `rate_equation` performance regression tests (unaffected — both
changes are pre-fit control flow).
