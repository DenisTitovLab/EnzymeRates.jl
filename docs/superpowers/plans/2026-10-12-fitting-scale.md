# Fitting Scale (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the analytic kcat scale of relative fits with the peak turnover over a
concentration range, apply measured kcats after the fit, select absolute mode with
`absolute = true`, and gate `rate_equation`'s speed above 400 denominator terms with a
bound proportional to the law's size, calibrated on CI.

**Architecture:** A deterministic multi-start coordinate search in log10 concentrations
(`_peak_turnover`, `src/fitting.jl`) finds a mechanism's largest net forward turnover
within `peak_range`. `fit_rate_equation` divides every steady-state rate constant of a
relative fit by it and reports whether the peak sits at the range's upper edge.
`rescale_parameter_values` applies the same scale, or a scale from measured turnovers,
after the fit. `_kcat_forward` and its helpers are deleted. The search loop takes the
turnover behind a `FunctionWrapper`, so it compiles once for every mechanism.

**Tech Stack:** Julia 1.10+ (CI runs the latest release), Optimization.jl,
FunctionWrappers.jl (already a dependency), Test, TestEnv for focused runs.

**Spec:** `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`, sections
11 (code generation and speed gates), 12 (fitting scale), 15 (phases) and 16 (docs and
ROADMAP). Read sections 12 and 18 before starting any task.

## Global Constraints

- Branch `fitting-scale` in the worktree
  `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/fitting_scale_wt`
  (below: `$WT`). It is stacked on PRs #85, #86, #87, #88 and #89, merged in locally;
  never rebase or rewrite those merge commits. Never touch the main checkout
  `/home/denis.linux/.julia/dev/EnzymeRates`.
- Run Julia only through the shared lock, one process at a time (7.7 GB, no swap), from
  inside `$WT`:
  `timeout -s KILL <seconds> flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --project=. --heap-size-hint=2500M <args>`
  (below: `$JL <args>`). Never park on a background monitor.
- Focused test run of one file (about 1–3 minutes):
  `$JL -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'`
- Full suite (about 6 minutes; precompile parallelism capped so precompilation fits in
  memory): `JULIA_NUM_PRECOMPILE_TASKS=2 $JL -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
  with `timeout -s KILL 3600`. Run it before every commit that changes `src/`.
- Code style: 92-character lines, 4-space indentation, match the surrounding code. Every
  new file starts with two `# ABOUTME: ` lines. Names say what code does; no "new",
  "old", "legacy", "improved" or history in names or comments. Comments are evergreen.
- Tests: every helper under `test/` carries the `_testhelper_` prefix (only a one-line
  closure local to a testset keeps a bare name). Test output must be pristine. Never
  delete a failing test; a test of removed behavior is rewritten to test the behavior
  that replaces it, keeping its fixture.
- `rate_equation` stays allocation-free and under 120 ns per call for every mechanism in
  `MECHANISM_TEST_SPECS` and the three largest enumerated bi-bi laws.
- Decided behavior (spec section 12, Denis's decisions E9): in relative mode
  (`absolute = false`, the default) `fit_rate_equation` and `identify_rate_equation`
  return constants scaled so the mechanism's peak forward turnover per enzyme, over every
  metabolite's concentration (products, inhibitors and regulators included) within
  `peak_range` (default `(1e-12, 1e3)`, in the data's concentration units), is 1. Every
  law is scaled; a law still rising at the upper edge by more than 10% over the range's
  last decade is flagged (`peak_at_range_edge`); the scale never turns a fit into a
  failure. Absolute mode never rescales. Measured kcats apply after the fit through
  `rescale_parameter_values(m, params, kcats)`, which returns `(params, spread)`.
  `scale_k_to_kcat` and `_kcat_forward` are removed; identify has no kcat keyword.
- Commits: never `git add -A` without `git status` first; never skip hooks; every
  commit message ends with exactly these two lines:
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`
  `Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE`
- PRs: never merge. A PR description never contains a claude.ai/code/session link; it
  ends with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`. When
  `gh pr edit` fails, use `gh api -X PATCH repos/DenisTitovLab/EnzymeRates.jl/pulls/<n>`.
- Write prose (docstrings, docs, comments, PR text) clearly and concisely: active voice,
  definite statements, no needless words.

## Review Focus

1. Data in µM or nM: the default range ends below the law's saturation, so most fits are
   flagged; they must still be scaled, flagged, and the docs must say how to widen the
   range (Task 2 tests a capped range; Task 7 documents it).
2. A law that needs a product or a modifier bound to turn over: the peak sits at an
   interior product concentration and must be found, not reported as zero (Task 1).
3. A law that is NaN or infinite somewhere in the range (extreme constants): the search
   must skip those points and never return NaN or Inf as the peak (Task 1 search-core
   test; Task 2 tests a fit whose turnover underflows to zero everywhere).
4. Anchors given with integer values, a reverse (negative) measured turnover, or a
   misspelled metabolite: integers and reverse anchors work, a misspelling is a clear
   error, never a silently ignored key (Task 4).
5. A fit with no positive turnover anywhere in the range: `fit_rate_equation` returns it
   unscaled without error, and an explicit `rescale_parameter_values(m, params)` call
   errors with a message naming the range (Tasks 2 and 4).

---

## File Structure

- `src/fitting.jl` (modify): `FittingProblem` gains `absolute` and `peak_range`;
  `fit_rate_equation` applies the peak scale; new: `DEFAULT_PEAK_RANGE`, `_peak_range`,
  `LogConcentrationTurnover`, `PEAK_SCAN_STEP`, `_maximize_on_box`,
  `_coordinate_ascent!`, `_line_maximum!`, `_finite_or_minus_inf`, `_peak_turnover`,
  `_scale_rate_constants`, `_check_has_keq`, and both `rescale_parameter_values` methods.
- `src/rate_eq_derivation.jl` (modify): delete `_kcat_forward` (both methods),
  `_kcat_groups_from_polys`, `_kcat_keys`, `_max_expr`, `_mwc_power_pair` and the keyword
  `rescale_parameter_values`; keep `_ss_rate_constant_names`.
- `src/identify_rate_equation.jl` (modify): `IdentifyRateEquationProblem` gains `absolute`
  and `peak_range`; the CSVs gain `peak_at_range_edge`; the `maxtime` docstring is fixed.
- `test/test_fitting_scale.jl` (create): the peak search and `rescale_parameter_values`.
- `test/runtests.jl`, `test/test_fitting.jl`, `test/test_identify_rate_equation.jl`,
  `test/test_compile_budget.jl`, `test/test_rate_eq_derivation.jl`,
  `test/allosteric_ground_truth.jl`, `test/test_mechanism_enumeration.jl`,
  `test/mechanism_definitions_for_test_enzyme_derivation.jl` (modify).
- `docs/src/fitting/normalized_vs_absolute.md`, `docs/src/fitting/tutorial.md`,
  `docs/src/fitting/loss_and_optimizers.md`, `docs/src/identify/model_selection.md`,
  `docs/src/developer.md`, `ROADMAP.md`, `Project.toml` (modify).

---

### Task 1: Find a mechanism's peak turnover over a concentration range

**Files:**
- Modify: `src/fitting.jl` (append after `fit_rate_equation`)
- Create: `test/test_fitting_scale.jl`
- Modify: `test/runtests.jl` (include the new file after `test_fitting.jl`)

**Interfaces:**
- Consumes: `metabolites(m)`, `rate_equation(m, concs, params)`,
  `FunctionWrapper` (already imported in `src/EnzymeRates.jl`).
- Produces (later tasks rely on these exact names):
  - `const DEFAULT_PEAK_RANGE = (1e-12, 1e3)`
  - `_peak_range(peak_range) -> Tuple{Float64,Float64}` (validates; error text below)
  - `_peak_turnover(m::AbstractEnzymeMechanism, params::NamedTuple,
    peak_range::Tuple{Float64,Float64}) -> (value::Float64, concs::NamedTuple,
    at_edge::Bool)`; `params` must hold `Keq`; its `E_total` is ignored; `value` is not
    positive (possibly `-Inf`) when the law has no forward turnover in the range.
  - `_maximize_on_box(f::FunctionWrapper{Float64,Tuple{Vector{Float64}}}, n::Int,
    a::Float64, b::Float64) -> (u::Vector{Float64}, value::Float64, at_edge::Bool)`

- [ ] **Step 1: Write the failing tests**

Create `test/test_fitting_scale.jl`:

```julia
# ABOUTME: Tests for the peak-turnover scale of relative fits: the search for a law's
# ABOUTME: peak turnover over concentrations, and rescale_parameter_values.
using Random

const _testhelper_TurnoverFunction = EnzymeRates.FunctionWrapper{Float64,
                                                                 Tuple{Vector{Float64}}}

# The largest turnover of `m` at `p` on the grid of whole decades from 10^lo_exp to
# 10^hi_exp for every metabolite; undefined points count as -Inf.
function _testhelper_grid_peak(m, p, lo_exp, hi_exp)
    mets = metabolites(m)
    turnover_params = merge(p, (E_total = 1.0,))
    maximum(Iterators.product(ntuple(_ -> lo_exp:hi_exp, length(mets))...)) do e
        v = rate_equation(m, NamedTuple{mets}(Tuple(10.0 .^ e)), turnover_params)
        isfinite(v) ? v : -Inf
    end
end

@testset "Peak turnover scale" begin

    # Uni-uni with a steady-state binding; at Keq = 2 these constants give kcat = 5 and
    # the Michaelis constant Km. Fitted: k_ES_to_E_P, k_ES_to_E_S, k_E_S_to_ES.
    mm = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    _testhelper_mm_params(Km) = (k_ES_to_E_P = 5.0, k_ES_to_E_S = 5.0,
                                 k_E_S_to_ES = 10.0 / Km, Keq = 2.0)

    # S also binds E(S) as a dead-end inhibitor. At zero product,
    # v = k·x/(1 + x + x²·K_S/K_I) with x = S/K_S, maximal at S = √(K_S·K_I), where
    # v = k/(1 + 2√(K_S/K_I)).
    si = @enzyme_mechanism begin
        substrates: S
        products:   P
        regulators: S
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            E(S) + S::Inh ⇌ E(S, S::Inh)
        end
    end
    si_params = (K_EP_to_E_P = 1.0, K_ESSinh_to_ES_Sinh = 4.0, K_ES_to_E_S = 1.0,
                 k_ES_to_EP = 3.0, Keq = 2.0)

    # Random-order steady-state bi-bi.
    rob = @enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(P) + Q
            E(P, Q) <--> E(Q) + P
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end

    # Allosteric uni-uni whose two regulators bind sites of multiplicity 2 and 1.
    two_site = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X::OnlyI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X
        end
        regulatory_site(multiplicity = 1): begin
            ligands: Y
        end
    end

    @testset "a saturating law peaks at the range's upper substrate corner" begin
        peak = EnzymeRates._peak_turnover(mm, _testhelper_mm_params(1.0), (1e-12, 1e3))
        # v = kcat·S/(Km + S) at S = 1e3, with the product at the range's lower bound.
        @test peak.value ≈ 5 * 1e3 / (1 + 1e3) rtol = 1e-12
        @test peak.concs.S ≈ 1e3 rtol = 1e-12
        @test peak.concs.P ≈ 1e-12 rtol = 1e-12
        @test peak.at_edge === false
        # The range bounds the search: capped at S = 10 the peak is the turnover there.
        capped = EnzymeRates._peak_turnover(mm, _testhelper_mm_params(1.0), (1e-12, 10.0))
        @test capped.value ≈ 5 * 10 / (1 + 10) rtol = 1e-12
        @test capped.concs.S ≈ 10.0 rtol = 1e-12
    end

    @testset "a law still rising at the range's upper edge is flagged" begin
        # Hyperbolic in S: v(1e3)/v(1e2) = 10·(Km + 1e2)/(Km + 1e3) exceeds 1.1 exactly
        # when Km > 1e3/89 ≈ 11.24.
        flagged, unflagged = _testhelper_mm_params(11.5), _testhelper_mm_params(10.5)
        @test EnzymeRates._peak_turnover(mm, flagged, (1e-12, 1e3)).at_edge
        @test !EnzymeRates._peak_turnover(mm, unflagged, (1e-12, 1e3)).at_edge
        # Capped at S = 10 with Km = 1: v(10)/v(1) = (10/11)/(1/2) ≈ 1.82.
        p = _testhelper_mm_params(1.0)
        @test EnzymeRates._peak_turnover(mm, p, (1e-12, 10.0)).at_edge
    end

    @testset "substrate inhibition peaks inside the range" begin
        peak = EnzymeRates._peak_turnover(si, si_params, (1e-12, 1e3))
        @test peak.value ≈ 3 / (1 + 2 * sqrt(1 / 4)) rtol = 1e-9
        @test peak.concs.S ≈ 2.0 rtol = 1e-6
        @test peak.at_edge === false
    end

    @testset "a product can raise the peak above the product-free turnover" begin
        # Every constant 1 except E(P) ⇌ E + P at 1e-3 both ways (Keq = 1): Q rebinding
        # E(P) reroutes the slow release of P through E(Q) + P, so the turnover peaks at
        # an intermediate Q, 125.62 times the turnover at saturating A and B without
        # product.
        pn = EnzymeRates.fitted_params(rob)
        p = merge(NamedTuple{pn}(Tuple(k === :k_E_P_to_EP ? 1e-3 : 1.0 for k in pn)),
                  (Keq = 1.0,))
        peak = EnzymeRates._peak_turnover(rob, p, (1e-12, 1e8))
        corner = rate_equation(rob, (A = 1e8, B = 1e8, P = 1e-12, Q = 1e-12),
                               merge(p, (E_total = 1.0,)))
        @test peak.value ≈ 0.249991171308 rtol = 1e-9
        @test peak.value / corner ≈ 125.62056353 rtol = 1e-8
        @test peak.concs.Q ≈ 14120.94 rtol = 1e-4
        @test peak.concs.A ≈ 1e8 rtol = 1e-12
        @test peak.concs.P ≈ 1e-12 rtol = 1e-12
    end

    @testset "an allosteric law can peak below saturating substrate" begin
        # The active conformation turns over fast but binds S weakly; the inactive one
        # binds S tightly and turns over slowly. At catalytic multiplicity 2, saturating
        # S drives the enzyme into the inactive conformation, so the turnover peaks near
        # S = 1 at about 5 and falls toward 0.2 at the range's upper edge.
        mwc = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        p = (K_A_ES_to_E_S = 100.0, K_EP_to_E_P = 1.0, K_I_ES_to_E_S = 0.01,
             k_A_ES_to_EP = 1000.0, k_I_ES_to_EP = 0.1, L = 1e-4, Keq = 2.0)
        peak = EnzymeRates._peak_turnover(mwc, p, (1e-12, 1e3))
        @test peak.value ≈ 5.00100010001 rtol = 1e-9
        @test 0.5 < peak.concs.S < 2
        @test peak.value > 20 * rate_equation(mwc, (S = 1e3, P = 1e-12),
                                              merge(p, (E_total = 1.0,)))
        @test peak.at_edge === false
    end

    @testset "a law that needs bound product peaks at a positive product level" begin
        # Catalysis runs only on the product-bound form E(P), so the forward rate needs P
        # and vanishes as P goes to zero.
        ep = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + P ⇌ E(P)
                E(P) + S ⇌ E(P, S)
                E(P, S) <--> E(P, P)
                E(P, P) ⇌ E(P) + P
            end
        end
        pn = EnzymeRates.fitted_params(ep)
        p = merge(NamedTuple{pn}(ntuple(_ -> 1.5, length(pn))), (Keq = 2.0,))
        peak = EnzymeRates._peak_turnover(ep, p, (1e-12, 1e3))
        @test peak.value ≈ 1.49227389554 rtol = 1e-9
        @test peak.concs.S ≈ 1e3 rtol = 1e-12
        @test 1e-3 < peak.concs.P < 1e3
    end

    @testset "the peak reaches every whole-decade grid point's turnover" begin
        for (m, n_draws) in ((rob, 5), (two_site, 3))
            pn = EnzymeRates.fitted_params(m)
            for seed in 1:n_draws
                rng = Random.MersenneTwister(seed)
                p = merge(NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn)),
                          (Keq = 2.0,))
                peak = EnzymeRates._peak_turnover(m, p, (1e-12, 1e3))
                @test peak.value >= _testhelper_grid_peak(m, p, -12, 3) * (1 - 1e-12)
                # The reported point reaches the reported value.
                @test rate_equation(m, peak.concs, merge(p, (E_total = 1.0,))) ==
                      peak.value
            end
        end
    end

    @testset "the peak is per enzyme: E_total in params is ignored" begin
        p = _testhelper_mm_params(1.0)
        @test EnzymeRates._peak_turnover(mm, merge(p, (E_total = 7.0,)),
                                         (1e-12, 1e3)).value ==
              EnzymeRates._peak_turnover(mm, p, (1e-12, 1e3)).value
    end

    @testset "the search is deterministic" begin
        rng = Random.MersenneTwister(7)
        pn = EnzymeRates.fitted_params(rob)
        p = merge(NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn)), (Keq = 2.0,))
        @test EnzymeRates._peak_turnover(rob, p, (1e-12, 1e3)) ===
              EnzymeRates._peak_turnover(rob, p, (1e-12, 1e3))
    end

    @testset "the search skips undefined points and flags a rising upper edge" begin
        # A concave bowl with its top at (-1, 2), undefined wherever u[1] > 0.
        bowl = _testhelper_TurnoverFunction(
            u -> u[1] > 0 ? NaN : -((u[1] + 1)^2 + (u[2] - 2)^2))
        u, v, at_edge = EnzymeRates._maximize_on_box(bowl, 2, -12.0, 3.0)
        @test v ≈ 0 atol = 1e-12
        @test u ≈ [-1.0, 2.0] atol = 1e-5
        @test at_edge === false
        # Rising up to the upper edge: f(3) = 3 exceeds 1.1·f(2).
        ramp = _testhelper_TurnoverFunction(u -> u[1])
        u, v, at_edge = EnzymeRates._maximize_on_box(ramp, 1, -12.0, 3.0)
        @test u == [3.0] && v == 3.0 && at_edge === true
        # A function undefined everywhere has no peak.
        nowhere = _testhelper_TurnoverFunction(u -> NaN)
        @test EnzymeRates._maximize_on_box(nowhere, 2, -12.0, 3.0)[2] == -Inf
    end

    @testset "peak_range must be a tuple (lo, hi) with 0 < lo < hi < Inf" begin
        @test EnzymeRates.DEFAULT_PEAK_RANGE === (1e-12, 1e3)
        @test EnzymeRates._peak_range((1, 1000)) === (1.0, 1000.0)
        @test EnzymeRates._peak_range((1e-12, 1e3)) === (1e-12, 1e3)
        for bad in ((0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (1e-12, Inf), (NaN, 1.0),
                    (-1.0, 1.0), 1.0, (1e-12, 1e3, 1e6), [1e-12, 1e3])
            @test_throws(ErrorException("peak_range must be a tuple (lo, hi) with " *
                                        "0 < lo < hi < Inf; got $(repr(bad))"),
                         EnzymeRates._peak_range(bad))
        end
    end

end
```

In `test/runtests.jl`, add `include("test_fitting_scale.jl")` on the line after
`include("test_fitting.jl")`.

- [ ] **Step 2: Run the tests to verify they fail**

Run the focused command for `test_fitting_scale`.
Expected: errors with `UndefVarError` for `_peak_turnover`, `_maximize_on_box`,
`DEFAULT_PEAK_RANGE` and `_peak_range` (the feature is missing).

- [ ] **Step 3: Implement the search**

Append to `src/fitting.jl`:

```julia
"""
The default `peak_range`: the concentrations, in the data's units, over which the peak
turnover that scales a relative fit is sought. It covers physiological concentrations in M
or in mM.
"""
const DEFAULT_PEAK_RANGE = (1e-12, 1e3)

"""
`peak_range` as a `Tuple{Float64,Float64}`, after checking that it is a tuple `(lo, hi)`
with `0 < lo < hi < Inf`.
"""
function _peak_range(peak_range)
    peak_range isa Tuple{Real,Real} && 0 < peak_range[1] < peak_range[2] < Inf ||
        error("peak_range must be a tuple (lo, hi) with 0 < lo < hi < Inf; got " *
              repr(peak_range))
    (Float64(peak_range[1]), Float64(peak_range[2]))
end

"""
The turnover `v/E_total` of `mechanism` at `params` (whose `E_total` is 1) as a function
of the log10 concentrations of `metabolites(mechanism)`, in that order; `N` counts them.
"""
struct LogConcentrationTurnover{M<:AbstractEnzymeMechanism, P<:NamedTuple, N}
    mechanism::M
    params::P
end

LogConcentrationTurnover(mechanism::AbstractEnzymeMechanism, params::NamedTuple) =
    LogConcentrationTurnover{typeof(mechanism), typeof(params),
                             length(metabolites(mechanism))}(mechanism, params)

(t::LogConcentrationTurnover{M, P, N})(u::Vector{Float64}) where {M, P, N} =
    rate_equation(t.mechanism,
        NamedTuple{metabolites(t.mechanism)}(ntuple(i -> exp10(u[i]), Val(N))), t.params)

"""Decades between the points of one coordinate's scan in the peak search."""
const PEAK_SCAN_STEP = 0.25

"""`v`, or `-Inf` when `v` is NaN or infinite, so the peak search never picks it."""
_finite_or_minus_inf(v) = isfinite(v) ? v : -Inf

"""
The largest value of `f` over the box `[a, b]^n` of log10 concentrations, the point that
reaches it, and whether `f` still rises at the box's upper edge. Points where `f` is NaN
or infinite never win.

The search is deterministic. It evaluates every corner of the box (the two extreme
corners only, above 12 coordinates) and its center, then climbs from the four best of
them one coordinate at a time (`_coordinate_ascent!`). `f` still rises at the upper edge
when some coordinate of the best point sits at `b` and lowering it by one decade (to `a`
at most) lowers `f` by more than 10%.
"""
function _maximize_on_box(f::FunctionWrapper{Float64, Tuple{Vector{Float64}}}, n::Int,
                          a::Float64, b::Float64)
    corners = n <= 12 ? (0:(2^n - 1)) : (0, 2^n - 1)
    starts = [Float64[(k >> (j - 1)) & 1 == 1 ? b : a for j in 1:n] for k in corners]
    push!(starts, fill((a + b) / 2, n))
    values = [_finite_or_minus_inf(f(u)) for u in starts]
    order = sort(eachindex(starts); by = i -> (-values[i], i))
    best_u, best_v = starts[order[1]], values[order[1]]
    for i in order[1:min(4, length(order))]
        u = copy(starts[i])
        v = _coordinate_ascent!(f, u, values[i], a, b)
        v > best_v && ((best_u, best_v) = (u, v))
    end
    at_edge = any(1:n) do j
        best_u[j] >= b - 1e-6 || return false
        lowered = copy(best_u)
        lowered[j] = max(a, b - 1)
        best_v > 1.1 * f(lowered)
    end
    best_u, best_v, at_edge
end

"""
Climb `f` from `u`, where its value is `v`, one coordinate at a time (`_line_maximum!`),
in place, until a full sweep gains less than one part in 10¹²; returns the final value.
"""
function _coordinate_ascent!(f::FunctionWrapper{Float64, Tuple{Vector{Float64}}},
                             u::Vector{Float64}, v::Float64, a::Float64, b::Float64)
    for _ in 1:100
        v_start = v
        for j in eachindex(u)
            v = _line_maximum!(f, u, j, v, a, b)
        end
        v > v_start && v - v_start > 1e-12 * abs(v) || break
    end
    v
end

"""
Move `u[j]` to the best point of `f` along coordinate `j` within `[a, b]`, in place, and
return the value there; `v` is the value at the current `u`. Scans the coordinate every
`PEAK_SCAN_STEP` decades, then refines the best scan point by golden-section search
between its neighbours. The current point stays unless a candidate beats it.
"""
function _line_maximum!(f::FunctionWrapper{Float64, Tuple{Vector{Float64}}},
                        u::Vector{Float64}, j::Int, v::Float64, a::Float64, b::Float64)
    n_steps = max(1, ceil(Int, (b - a) / PEAK_SCAN_STEP))
    h = (b - a) / n_steps
    best_t, best_v = u[j], v
    for k in 0:n_steps
        u[j] = k == n_steps ? b : a + k * h
        w = _finite_or_minus_inf(f(u))
        w > best_v && ((best_t, best_v) = (u[j], w))
    end
    lo, hi = max(a, best_t - h), min(b, best_t + h)
    φ = (sqrt(5) - 1) / 2
    x1, x2 = hi - φ * (hi - lo), lo + φ * (hi - lo)
    u[j] = x1; f1 = _finite_or_minus_inf(f(u))
    u[j] = x2; f2 = _finite_or_minus_inf(f(u))
    while hi - lo > 1e-9
        if f1 >= f2
            hi, x2, f2 = x2, x1, f1
            x1 = hi - φ * (hi - lo)
            u[j] = x1; f1 = _finite_or_minus_inf(f(u))
        else
            lo, x1, f1 = x1, x2, f2
            x2 = lo + φ * (hi - lo)
            u[j] = x2; f2 = _finite_or_minus_inf(f(u))
        end
    end
    f1 > best_v && ((best_t, best_v) = (x1, f1))
    f2 > best_v && ((best_t, best_v) = (x2, f2))
    u[j] = best_t
    best_v
end

"""
    _peak_turnover(m, params, peak_range) -> (; value, concs, at_edge)

The peak turnover of mechanism `m` at `params`, which must hold `Keq`: the largest net
forward rate per enzyme, `v/E_total`, over every metabolite's concentration within
`peak_range` — substrates, products, inhibitors and regulators alike
(`_maximize_on_box`). `concs` is the point that reaches it, and `at_edge` says whether
the turnover still rises at the range's upper edge. The `E_total` of `params` is
ignored. `value` is not positive when the mechanism has no forward turnover in the range.
The turnover sits behind a `FunctionWrapper`, so the search compiles once for every
mechanism.
"""
function _peak_turnover(@nospecialize(m::AbstractEnzymeMechanism),
                        @nospecialize(params::NamedTuple),
                        peak_range::Tuple{Float64, Float64})
    mets = metabolites(m)
    f = FunctionWrapper{Float64, Tuple{Vector{Float64}}}(
        LogConcentrationTurnover(m, merge(params, (E_total = 1.0,))))
    u, value, at_edge = _maximize_on_box(f, length(mets), log10(peak_range[1]),
                                         log10(peak_range[2]))
    (value = value, concs = NamedTuple{mets}(Tuple(exp10.(u))), at_edge = at_edge)
end
```

- [ ] **Step 4: Run the tests to verify they pass**

Run the focused command for `test_fitting_scale`.
Expected: every test passes. If the product-rescue, allosteric or E(P)-only values miss
their reference values, do not loosen the tolerance; find out why with
superpowers:systematic-debugging (the references come from a prototype of exactly this
algorithm).

- [ ] **Step 5: Run the full suite, then commit**

Run the full suite (Aqua and JET included). Expected: all pass. Then:

```bash
git status
git add src/fitting.jl test/test_fitting_scale.jl test/runtests.jl
git commit -m "Find a mechanism's peak turnover over a concentration range" \
  -m "A deterministic multi-start coordinate search in log10 concentrations finds the
largest net forward turnover per enzyme within a range, products and regulators
included, and reports whether it still rises at the range's upper edge." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 2: Select absolute mode with `absolute` and scale relative fits to peak turnover 1

**Files:**
- Modify: `src/fitting.jl` (`FittingProblem`, `_rate_table`, `loss!`,
  `fit_rate_equation`; add `_scale_rate_constants`)
- Modify: `src/identify_rate_equation.jl` (`IdentifyRateEquationProblem`, the
  `FittingProblem` call in `_process_batch`, `_fold_problem`, the `memo` comment in
  `_beam_search`)
- Modify: `test/test_fitting.jl`, `test/test_identify_rate_equation.jl`,
  `test/test_compile_budget.jl` (line 197–198)

**Interfaces:**
- Consumes: `_peak_turnover`, `_peak_range`, `DEFAULT_PEAK_RANGE` (Task 1);
  `_ss_rate_constant_names(m)` (`src/rate_eq_derivation.jl`).
- Produces:
  - `FittingProblem(m, table; Keq::Real, absolute::Bool = false,
    peak_range = DEFAULT_PEAK_RANGE)` with fields `absolute::Bool` and
    `peak_range::Tuple{Float64,Float64}` in place of `scale_k_to_kcat`.
  - `fit_rate_equation(...)` returns `(params, loss, retcode, peak_at_range_edge::Bool)`.
  - `IdentifyRateEquationProblem(reaction, table; Keq::Real, absolute::Bool = false,
    peak_range = DEFAULT_PEAK_RANGE)` with fields `absolute` and `peak_range`.
  - `_scale_rate_constants(m, params::NamedTuple, s::Real) -> NamedTuple` (same keys,
    every steady-state rate constant times `s`).
  - `_rate_table(table, mnames, Keq)` (three arguments).

- [ ] **Step 1: Write the failing tests**

In `test/test_fitting.jl`:

(a) Add at the top of the file, after `using Tables`:

```julia
using Optimization
using Optimization.SciMLBase: build_solution, ReturnCode, DefaultOptimizationCache

# An optimizer stub whose every solve returns `u` for each log-space coordinate, so a
# fit's raw constants are known exactly and the scale applied to them can be checked.
struct _testhelper_ConstantStubOpt
    u::Float64
end
Optimization.allowsbounds(::_testhelper_ConstantStubOpt) = true
function Optimization.SciMLBase.__solve(prob::Optimization.OptimizationProblem,
        opt::_testhelper_ConstantStubOpt; kwargs...)
    u = fill(opt.u, length(prob.u0))
    build_solution(DefaultOptimizationCache(prob.f, prob.p), opt, u, prob.f(u, prob.p);
                   retcode = ReturnCode.Success)
end
```

(b) In the testset "Absolute mode uncentered loss", replace both
`scale_k_to_kcat=nothing` with `absolute=true`, the trailing comment `# default 1.0`
with `# relative (default)`, and the section comment
`# ── Absolute mode: uncentered loss (scale_k_to_kcat=nothing) ──────────────` with
`# ── Absolute mode: uncentered loss (absolute=true) ─────────────────────────`. In the
testset "Zero allocations", replace `scale_k_to_kcat=nothing` with `absolute=true`.

(c) Replace the whole testset "scale_k_to_kcat validation" (and its section comment)
with:

```julia
    # ── absolute and peak_range keywords ──────────────────────────────────────
    @testset "absolute and peak_range keywords" begin
        ok_data = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1])
        fp = FittingProblem(uni_uni, ok_data; Keq=1.0)
        @test fp.absolute === false
        @test fp.peak_range === (1e-12, 1e3)
        @test FittingProblem(uni_uni, ok_data; Keq=1.0, absolute=true).absolute === true
        # Integer Keq and peak_range bounds convert to the Float64 fields.
        fp_int = FittingProblem(uni_uni, ok_data; Keq=2, peak_range=(1, 1000))
        @test fp_int.Keq === 2.0 && fp_int.peak_range === (1.0, 1000.0)
        # The concrete-mechanism constructor passes both keywords on.
        fp_conc = FittingProblem(EnzymeRates.Mechanism(uni_uni), ok_data; Keq=1.0,
                                 absolute=true, peak_range=(1e-6, 10.0))
        @test fp_conc.absolute === true && fp_conc.peak_range === (1e-6, 10.0)
        # There is no kcat keyword: measured kcats rescale a fit afterwards.
        @test_throws MethodError FittingProblem(uni_uni, ok_data; Keq=1.0,
                                                scale_k_to_kcat=1.0)
    end
```

(d) Replace the testset "scale_k_to_kcat normalization + retcode" (and its two comment
lines above it) with:

```julia
    # ── fit_rate_equation: relative fits come back at peak turnover 1 ─────────
    # The scale applies to every fitted point, so the assertions hold at any fit quality
    # and a half-second budget is enough.
    @testset "relative fits come back at peak turnover 1" begin
        using OptimizationBBO
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
            (S = 0.5, P = 0.5), (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        opt = BBO_adaptive_de_rand_1_bin_radiuslimited()

        # The default range (best of 2 restarts) and a custom one.
        for (fp_kwargs, n_restarts, range) in
                (((;), 2, (1e-12, 1e3)), ((; peak_range=(1e-6, 1e2)), 1, (1e-6, 1e2)))
            fp = FittingProblem(uni_uni, data; Keq=Keq_val, fp_kwargs...)
            res = _testhelper_fit_capturing_convergence() do
                fit_rate_equation(fp, opt; n_restarts, maxtime=0.5)
            end
            @test keys(res) == (:params, :loss, :retcode, :peak_at_range_edge)
            @test keys(res.params) == EnzymeRates.fitted_params(uni_uni)
            @test isfinite(res.loss)
            @test res.retcode isa Symbol
            @test res.peak_at_range_edge isa Bool
            peak = EnzymeRates._peak_turnover(uni_uni, merge(res.params, (Keq = Keq_val,)),
                                              range)
            @test peak.value ≈ 1.0 rtol=1e-9
        end
    end

    # ── fit_rate_equation: the scale on known raw constants ───────────────────
    @testset "peak scale, edge flag and absolute mode on known raw constants" begin
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        pn = EnzymeRates.fitted_params(uni_uni)
        # Every raw constant is exp(0) = 1. At Keq = 2 the peak lies at S = 1e3 and P at
        # the lower bound, v = 1e3/(2 + 1e3); every fitted constant here is a rate
        # constant, so each comes back as 1/v.
        rel = fit_rate_equation(FittingProblem(uni_uni, data; Keq=Keq_val),
                                _testhelper_ConstantStubOpt(0.0); n_restarts=1)
        @test all(v -> isapprox(v, (2 + 1e3) / 1e3; rtol=1e-12), values(rel.params))
        @test rel.peak_at_range_edge === false
        # Capped at S = 10 the law still rises at the upper edge,
        # v(10)/v(1) = (10/12)/(1/3) = 2.5, so the fit is flagged and scaled to the
        # turnover at the edge, v(10) = 10/12.
        edge = fit_rate_equation(
            FittingProblem(uni_uni, data; Keq=Keq_val, peak_range=(1e-12, 10.0)),
            _testhelper_ConstantStubOpt(0.0); n_restarts=1)
        @test edge.peak_at_range_edge === true
        @test all(v -> isapprox(v, 12 / 10; rtol=1e-12), values(edge.params))
        # Absolute mode returns the fitted constants as they are.
        absolute_fit = fit_rate_equation(
            FittingProblem(uni_uni, data; Keq=Keq_val, absolute=true),
            _testhelper_ConstantStubOpt(0.0); n_restarts=1)
        @test absolute_fit.params == NamedTuple{pn}((1.0, 1.0, 1.0))
        @test absolute_fit.peak_at_range_edge === false
        # At exp(-700) every product of two constants underflows to zero, so the law
        # has no positive turnover anywhere: the fit keeps its constants, unscaled.
        tiny = fit_rate_equation(FittingProblem(uni_uni, data; Keq=Keq_val),
                                 _testhelper_ConstantStubOpt(-700.0); n_restarts=1)
        @test all(==(exp(-700.0)), values(tiny.params))
        @test tiny.peak_at_range_edge === false
    end
```

(e) In the testset "Validation errors", replace the last block (comment
`# A Keq or a scale_k_to_kcat that is not positive and finite` and its loop) with:

```julia
        # A Keq that is not positive and finite, or a peak_range that is not a tuple
        # (lo, hi) with 0 < lo < hi < Inf
        data_ok = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1])
        for (bad, shown) in ((0, "0"), (-1, "-1"), (Inf, "Inf"), (NaN, "NaN"))
            @test_throws(ErrorException("Keq must be positive and finite; got $shown"),
                FittingProblem(uni_uni, data_ok; Keq=bad))
        end
        for bad in ((0.0, 1.0), (2.0, 1.0), (1e-12, Inf), 1.0)
            @test_throws(
                ErrorException("peak_range must be a tuple (lo, hi) with 0 < lo < hi < " *
                               "Inf; got $(repr(bad))"),
                FittingProblem(uni_uni, data_ok; Keq=1.0, peak_range=bad))
        end
```

In `test/test_identify_rate_equation.jl`:

(f) In the testset "construction", replace the block from
`# scale_k_to_kcat field: default 1.0, settable to nothing, validated.` through the
`@test_throws ErrorException IdentifyRateEquationProblem(... scale_k_to_kcat=0.0)` call
with:

```julia
        # absolute and peak_range: relative mode and the default range unless set.
        @test prob.absolute === false
        @test prob.peak_range === (1e-12, 1e3)
        prob_abs = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, absolute=true, peak_range=(1e-9, 1e2))
        @test prob_abs.absolute === true && prob_abs.peak_range === (1e-9, 1e2)
        @test_throws MethodError IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=1.0)
```

Replace the comment `# ≥ 0, or a Keq or scale_k_to_kcat that is not positive and
finite: each error` / `# names the column (with the row), Keq or scale_k_to_kcat.`
with `# ≥ 0, a Keq that is not positive and finite, or a peak_range out of order: each`
/ `# error names the column (with the row), Keq or peak_range.`, and replace the last
loop of that block (`for (bad, shown) in ((0, "0"), ...` with its two `@test_throws`)
with:

```julia
        for (bad, shown) in ((0, "0"), (-1, "-1"), (Inf, "Inf"), (NaN, "NaN"))
            @test_throws(ErrorException("Keq must be positive and finite; got $shown"),
                IdentifyRateEquationProblem(test_rxn, two_groups; Keq=bad))
        end
        @test_throws(
            ErrorException("peak_range must be a tuple (lo, hi) with 0 < lo < hi < Inf; " *
                           "got (1.0, 0.5)"),
            IdentifyRateEquationProblem(test_rxn, two_groups; Keq=1.0,
                                        peak_range=(1.0, 0.5)))
```

(g) In the testset "_fold_problem: the fold's rows in the full fit's problem type",
replace

```julia
        full = FittingProblem(em, prob.data; Keq=prob.Keq,
                              scale_k_to_kcat=prob.scale_k_to_kcat)
```

with

```julia
        full = FittingProblem(em, prob.data; Keq=prob.Keq, absolute=prob.absolute,
                              peak_range=prob.peak_range)
```

and the line `@test fold.Keq === 10.0 && fold.scale_k_to_kcat === 1.0` with:

```julia
        @test fold.Keq === 10.0 && fold.absolute === false &&
              fold.peak_range === (1e-12, 1e3)
        # A fold copies the problem's mode and range.
        prob_abs = IdentifyRateEquationProblem(prob.reaction, prob.data; Keq=prob.Keq,
                                               absolute=true, peak_range=(1e-9, 1e2))
        fold_abs = EnzymeRates._fold_problem(em, prob_abs,
                                             prob.data.group .== prob.data.group[3])
        @test fold_abs.absolute === true && fold_abs.peak_range === (1e-9, 1e2)
```

(h) In the testset "_cv_fold_loss: a fold compiles nothing for its mechanism", replace
`full = FittingProblem(em, prob.data; Keq=prob.Keq,` /
`scale_k_to_kcat=prob.scale_k_to_kcat)` with
`full = FittingProblem(em, prob.data; Keq=prob.Keq, absolute=prob.absolute,` /
`peak_range=prob.peak_range)`.

(i) In the `_process_batch` reuse test, replace the comment
`# scale_k_to_kcat=1.0 anchored kcat: the copied params are the rescaled fit,` /
`# not the raw 5.0 the stub optimizer returned.` with
`# The default peak scale: the copied params are the rescaled fit, not the raw` /
`# 5.0 the stub optimizer returned.`

In `test/test_compile_budget.jl` (the fit script, lines 197–198), replace
`Keq = PROB.Keq, scale_k_to_kcat = PROB.scale_k_to_kcat)` with
`Keq = PROB.Keq, absolute = PROB.absolute,` and a continuation line
`peak_range = PROB.peak_range)` aligned like its neighbours.

- [ ] **Step 2: Run the tests to verify they fail**

Run the focused commands for `test_fitting` and `test_identify_rate_equation`.
Expected: failures and errors such as `type FittingProblem has no field absolute`,
`MethodError: ... got unsupported keyword argument "absolute"` and
`keys(res) == (:params, :loss, :retcode, :peak_at_range_edge)` failing.

- [ ] **Step 3: Implement**

In `src/fitting.jl`:

1. Replace the `FittingProblem` struct docstring's `scale_k_to_kcat` field bullet with:

```
- `absolute`: `false` selects relative mode (per-group-centered loss, fit scaled to peak
  turnover 1); `true` selects absolute per-enzyme-turnover mode (uncentered loss, no
  scale)
- `peak_range`: the concentration range `(lo, hi)` of the peak-turnover scale
```

and the field `scale_k_to_kcat::Union{Float64,Nothing}` with the two fields
`absolute::Bool` and `peak_range::Tuple{Float64,Float64}` (in that order, after `Keq`).

2. Replace the outer constructor and its docstring with:

```julia
"""
    FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, absolute::Bool = false, peak_range = (1e-12, 1e3))

Construct a `FittingProblem` from an enzyme mechanism and tabular data.

The table must have columns: `group`, `Rate`, and one column per
metabolite matching `metabolites(mechanism)`, so no metabolite may be
named `group` or `Rate`. Uses `Tables.columntable` for conversion.

With `absolute = false` (the default) the data are relative: the loss centers each
group's log-ratios, and `fit_rate_equation` scales the fitted rate constants so the
mechanism's peak turnover — its largest net forward rate per enzyme over every
metabolite's concentration within `peak_range`, products and regulators included — is 1.
With `absolute = true` the data are absolute per-enzyme turnover: the loss is uncentered
and the fit is returned unscaled. `peak_range = (lo, hi)` is in the data's concentration
units; the default covers physiological concentrations in M or mM, and data in µM need a
higher upper bound.

Rate values must be finite and nonzero (the loss takes their log; a zero rate gives
`-Inf`), every concentration must be a finite number ≥ 0 (zero is valid), `Keq`
must be positive and finite, and `peak_range` must be a tuple `(lo, hi)` with
`0 < lo < hi < Inf`.
"""
function FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, absolute::Bool=false, peak_range=DEFAULT_PEAK_RANGE)
    data = _rate_table(table, metabolites(mechanism), Keq)
    group_map = Dict{eltype(data.group), Vector{Int}}()
    for (i, g) in enumerate(data.group)
        push!(get!(() -> Int[], group_map, g), i)
    end
    # The typed constructor converts Keq to the field type.
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, collect(values(group_map)), Keq, absolute,
        _peak_range(peak_range), log.(abs.(data.Rate)),
        Vector{Float64}(undef, length(data.Rate)))
end
```

3. Replace `_rate_table`'s docstring and first lines so it takes `(table, mnames, Keq)`:
delete the `scale_k_to_kcat` check (`scale_k_to_kcat === nothing || ... error(...)`),
and in the docstring delete "`scale_k_to_kcat` must be positive and finite or
`nothing`, " so it begins "Validate a rate table and return it as
`Tables.columntable(table)`: `Keq` must be positive and finite, the table needs ...".
Its signature line becomes `_rate_table(table, mnames, Keq) → NamedTuple`.

4. Replace the concrete-mechanism method with:

```julia
FittingProblem(mechanism::Union{Mechanism, AllostericMechanism}, table; kwargs...) =
    FittingProblem(compile_mechanism(mechanism), table; kwargs...)
```

(keep its docstring).

5. In `loss!`: in the docstring replace "When `fp.scale_k_to_kcat` is a `Real`, each
group's" with "In relative mode (`fp.absolute == false`) each group's", and "When it is
`nothing`, the log-ratios" with "In absolute mode (`fp.absolute`) the log-ratios". In
the body replace `if fp.scale_k_to_kcat === nothing` with `if fp.absolute`, and the
comment block above it with:

```julia
    # Pass 2: loss. Relative (fp.absolute false) → per-group mean-centered, removing
    # each group's arbitrary scale. Absolute → uncentered: the y-axis is absolute
    # per-enzyme turnover, so the absolute magnitude is meaningful.
```

6. In `fit_rate_equation`'s docstring replace the paragraph that begins "Rescaling is
driven by `fp.scale_k_to_kcat`" with:

```
In relative mode (`fp.absolute == false`) the returned rate constants are scaled so the
mechanism's peak turnover within `fp.peak_range` is 1 (see `FittingProblem`); a law with
no positive turnover in the range is returned unscaled. In absolute mode they are returned
as fitted.
```

and replace "Returns a NamedTuple `(params, loss, retcode)` where:" with "Returns a
NamedTuple `(params, loss, retcode, peak_at_range_edge)` where:", adding after the
`retcode` bullet:

```
- `peak_at_range_edge`: `true` when the scaled law still rises at the upper edge of
  `fp.peak_range` — by more than 10% over the range's last decade — so a wider range would
  change the scale; always `false` in absolute mode
```

7. Replace the end of `fit_rate_equation` (from `result_params = ...` to the `return`)
with:

```julia
    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(np)))
    peak_at_range_edge = false
    if !fp.absolute
        peak = _peak_turnover(fp.mechanism, merge(result_params, (Keq = fp.Keq,)),
                              fp.peak_range)
        peak_at_range_edge = peak.at_edge
        0 < peak.value < Inf && (result_params =
            _scale_rate_constants(fp.mechanism, result_params, 1 / peak.value))
    end
    return (params = result_params, loss = best_loss, retcode = best_retcode,
            peak_at_range_edge = peak_at_range_edge)
end

"""
`params` with every steady-state rate constant of `m` multiplied by `s`, every other
entry unchanged. A common factor on those constants multiplies every rate by the same
factor and leaves the binding constants, `Keq`, `E_total`, `L` and the regulatory
constants alone.
"""
function _scale_rate_constants(@nospecialize(m::AbstractEnzymeMechanism),
                               params::NamedTuple, s::Real)
    rate_constants = _ss_rate_constant_names(m)
    NamedTuple{keys(params)}(map((k, v) -> k in rate_constants ? v * s : v,
                                 keys(params), values(params)))
end
```

In `src/rate_eq_derivation.jl`, in the docstring of `_ss_rate_constant_names`, replace
"Used by `rescale_parameter_values` to scale only SS k's" with "Used by
`_scale_rate_constants` to scale only SS k's".

In `src/identify_rate_equation.jl`:

1. In the `IdentifyRateEquationProblem` docstring replace the `scale_k_to_kcat` bullet
with:

```
- `absolute`: `false` (the default) treats the data as relative and scales every fit to
  peak turnover 1; `true` treats it as absolute per-enzyme turnover (see `FittingProblem`)
- `peak_range`: the concentration range `(lo, hi)` of the peak-turnover scale (see
  `FittingProblem`)

Construct it with `IdentifyRateEquationProblem(reaction, table; Keq, absolute = false,
peak_range = (1e-12, 1e3))`.
```

2. Replace the field `scale_k_to_kcat::Union{Float64,Nothing}` with
`absolute::Bool` and `peak_range::Tuple{Float64,Float64}`; the constructor's keyword
`scale_k_to_kcat::Union{Real,Nothing}=1.0` with `absolute::Bool=false,
peak_range=DEFAULT_PEAK_RANGE`; `_rate_table(table, mnames, scale_k_to_kcat, Keq)` with
`_rate_table(table, mnames, Keq)`; and the final line with
`IdentifyRateEquationProblem{typeof(data)}(reaction, data, Keq, absolute,
_peak_range(peak_range))` (wrapped to 92 characters).

3. In `_process_batch` and in `_fold_problem`, replace
`Keq=prob.Keq, scale_k_to_kcat=prob.scale_k_to_kcat)` with
`Keq=prob.Keq, absolute=prob.absolute, peak_range=prob.peak_range)` (wrap to 92
characters).

4. In `_beam_search`, replace the comment line
`# Fits as `fit_rate_equation` returns them (rescaled when `scale_k_to_kcat` is` /
`# set), keyed by `eq_hash`, ...` with
`# Fits as `fit_rate_equation` returns them (scaled to peak turnover 1 in relative` /
`# mode), keyed by `eq_hash`, ...`, keeping the rest of the comment.

- [ ] **Step 4: Run the tests to verify they pass**

Run the focused commands for `test_fitting`, `test_fitting_scale` and
`test_identify_rate_equation`. Expected: all pass. Then grep:
`grep -rn "scale_k_to_kcat" src test` — expected: no match in `src/`; in `test/` only
the `@test_throws MethodError` lines from Steps 1(c) and 1(f) and the old-API keyword
uses inside `test/test_rate_eq_derivation.jl` and `test/allosteric_ground_truth.jl`
(`rescale_parameter_values(...; scale_k_to_kcat=...)`), which Task 4 rewrites.

- [ ] **Step 5: Run the full suite, then commit**

Run the full suite. Expected: all pass, including `test_compile_budget.jl`'s per-mechanism
fit budgets (`FIT_INSTANCE_BUDGET`, `FIT_NATIVE_BUDGET_KB`). If a budget fails, the peak
search compiles per mechanism somewhere; find where with
superpowers:systematic-debugging before changing anything.

```bash
git status
git add src/fitting.jl src/identify_rate_equation.jl src/rate_eq_derivation.jl \
  test/test_fitting.jl test/test_identify_rate_equation.jl test/test_compile_budget.jl
git commit -m "Scale relative fits to peak turnover 1; select absolute mode with absolute" \
  -m "FittingProblem and IdentifyRateEquationProblem take absolute and peak_range in
place of scale_k_to_kcat. In relative mode fit_rate_equation divides every steady-state
rate constant by the mechanism's peak turnover within peak_range and reports whether the
peak sits at the range's upper edge." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 3: Record the edge flag in identify's CSVs

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_process_batch` row, `_rows_to_dataframe`,
  `_failure_row`)
- Modify: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: `fit.peak_at_range_edge` (Task 2).
- Produces: a `peak_at_range_edge` column, after `fit_inherited`, in every identify CSV
  and in `cv_results`; `missing` for failed fits.

- [ ] **Step 1: Write the failing tests**

In the testset `_rows_to_dataframe`:
- add `peak_at_range_edge = true,` after `fit_inherited = false,` in the first `rows`
  entry, and after the existing `@test df.fit_inherited == [false]` add
  `@test df.peak_at_range_edge == [true]`;
- in `rows_with_failure`, add `peak_at_range_edge = false` to the fitted row and
  `peak_at_range_edge = missing` to the failure row (after their `fit_inherited`
  entries), and after `@test ismissing(df_with_failure.fit_inherited[2])` add:

```julia
        @test df_with_failure.peak_at_range_edge[1] == false
        @test ismissing(df_with_failure.peak_at_range_edge[2])
        # A row without the field records missing.
        @test ismissing(EnzymeRates._rows_to_dataframe([(n_params = 1, loss = 0.1,
            mechanism_type = "M", rate_equation = "v", retcode = "Success",
            error = missing, params = (;), eq_hash = "0", fit_inherited = false)
            ]).peak_at_range_edge[1])
```

In the testset "csv writer", add `peak_at_range_edge = true,` after
`fit_inherited = false,` in its row and, after `@test df.K_a == [2.0]`, add
`@test df.peak_at_range_edge == [true]`.

In the `_process_batch` reuse test, after
`@test [e.row.fit_inherited for e in entries] == [false, true]` add:

```julia
    # The edge flag belongs to the equation's fit, so the inherited row copies it.
    @test entries[1].row.peak_at_range_edge isa Bool
    @test entries[2].row.peak_at_range_edge === entries[1].row.peak_at_range_edge
```

In the testset "failure row preserves round-trippable mechanism", after
`@test row.error == "boom"` add `@test ismissing(row.peak_at_range_edge)`.

- [ ] **Step 2: Run the tests to verify they fail**

Run the focused command for `test_identify_rate_equation`.
Expected: failures on `peak_at_range_edge` (no such column / field).

- [ ] **Step 3: Implement**

In `_process_batch`, add `peak_at_range_edge = fit.peak_at_range_edge,` after
`fit_inherited = inherited,` in the row. In `_rows_to_dataframe`, after the
`fit_inherited = ...` column add:

```julia
        peak_at_range_edge =
            Union{Missing,Bool}[get(r, :peak_at_range_edge, missing) for r in rows],
```

In `_failure_row`, add `peak_at_range_edge = missing` after `fit_inherited = missing`
(move the closing parenthesis).

- [ ] **Step 4: Run the tests to verify they pass**

Run the focused command for `test_identify_rate_equation`. Expected: all pass.

- [ ] **Step 5: Run the full suite, then commit**

```bash
git status
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Record in identify's CSVs whether a fit's peak sits at the range edge" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 4: Rescale a fit to measured turnovers or to the peak scale after the fit

**Files:**
- Modify: `src/fitting.jl` (add both `rescale_parameter_values` methods and
  `_check_has_keq`)
- Modify: `src/rate_eq_derivation.jl` (delete the keyword method
  `rescale_parameter_values(m, params; scale_k_to_kcat)` and its heading comment
  `# ─── Public API: rescale_parameter_values ───...`)
- Modify: `test/test_fitting_scale.jl`, `test/test_rate_eq_derivation.jl`
  (`_testhelper_test_kcat_rescaling` and its call), `test/allosteric_ground_truth.jl`
  (lines 126–129)

**Interfaces:**
- Consumes: `_peak_turnover`, `_peak_range`, `DEFAULT_PEAK_RANGE` (Task 1);
  `_scale_rate_constants` (Task 2).
- Produces (public, exported):
  - `rescale_parameter_values(m::AbstractEnzymeMechanism, params::NamedTuple;
    peak_range = DEFAULT_PEAK_RANGE) -> (NamedTuple, Nothing)`
  - `rescale_parameter_values(m::AbstractEnzymeMechanism, params::NamedTuple,
    kcats::AbstractVector{<:Pair}) -> (NamedTuple, Float64)`

- [ ] **Step 1: Write the failing tests**

Append to `test/test_fitting_scale.jl`, inside the outer testset (before its final
`end`):

```julia
    @testset "rescale_parameter_values: peak scale" begin
        p = _testhelper_mm_params(1.0)
        q, spread = rescale_parameter_values(mm, p)
        @test spread === nothing
        @test keys(q) == keys(p)
        # Every fitted constant here is a rate constant; each is divided by the peak
        # turnover 5·1e3/(1 + 1e3).
        @test q.k_ES_to_E_P ≈ 1.001 rtol = 1e-12
        @test q.k_ES_to_E_S ≈ 1.001 rtol = 1e-12
        @test q.k_E_S_to_ES ≈ 2.002 rtol = 1e-12
        @test q.Keq == 2.0
        @test EnzymeRates._peak_turnover(mm, q, (1e-12, 1e3)).value ≈ 1 rtol = 1e-12
        # E_total passes through and does not enter the scale.
        qE, _ = rescale_parameter_values(mm, merge(p, (E_total = 3.0,)))
        @test qE.E_total == 3.0 && qE.k_ES_to_E_P == q.k_ES_to_E_P
        # peak_range sets the range: capped at S = 10 the peak turnover is 5·10/11.
        q10, _ = rescale_parameter_values(mm, p; peak_range = (1e-12, 10.0))
        @test q10.k_ES_to_E_P ≈ 1.1 rtol = 1e-12
        # Binding constants stay; only the rate constant of the substrate-inhibition law
        # scales, by 1/1.5.
        q_si, _ = rescale_parameter_values(si, si_params)
        @test q_si.k_ES_to_EP ≈ 2.0 rtol = 1e-9
        @test (q_si.K_EP_to_E_P, q_si.K_ESSinh_to_ES_Sinh, q_si.K_ES_to_E_S) ==
              (1.0, 4.0, 1.0)
    end

    @testset "rescale_parameter_values: measured turnovers" begin
        p = _testhelper_mm_params(1.0)
        _testhelper_mm_turnover(q, concs) = rate_equation(
            mm, merge((S = 0.0, P = 0.0), concs), merge(q, (E_total = 1.0,)))
        # One anchor: the rescaled law predicts the measured turnover there; P, unlisted,
        # is zero.
        q, spread = rescale_parameter_values(mm, p, [(S = 10.0,) => 42.0])
        @test _testhelper_mm_turnover(q, (S = 10.0,)) ≈ 42.0 rtol = 1e-12
        @test spread == 0.0
        @test q.Keq == 2.0
        # Two anchors combine by the geometric mean of measured over predicted turnover
        # (predicted 50/11 at S = 10 and 5/2 at S = 1); the spread is the standard
        # deviation of the two log ratios.
        r1, r2 = 42.0 / (50 / 11), 10.0 / (5 / 2)
        q2, spread2 = rescale_parameter_values(mm, p,
                                               [(S = 10.0,) => 42.0, (S = 1.0,) => 10.0])
        @test q2.k_ES_to_E_P ≈ 5.0 * sqrt(r1 * r2) rtol = 1e-12
        @test spread2 ≈ abs(log(r1) - log(r2)) / sqrt(2) rtol = 1e-12
        # Integer concentrations and turnovers work.
        @test rescale_parameter_values(mm, p, [(S = 10,) => 42])[1] == q
        # A reverse anchor: the product alone runs the reaction backwards.
        qr, _ = rescale_parameter_values(mm, p, [(P = 1.0,) => -2.0])
        @test _testhelper_mm_turnover(qr, (P = 1.0,)) ≈ -2.0 rtol = 1e-12
    end

    @testset "rescale_parameter_values: errors" begin
        p = _testhelper_mm_params(1.0)
        nokeq = (k_ES_to_E_P = 5.0, k_ES_to_E_S = 5.0, k_E_S_to_ES = 10.0)
        keq_message = "params must include Keq, which the dependent rate constants " *
                      "need; merge it into a fit's params"
        @test_throws ErrorException(keq_message) rescale_parameter_values(mm, nokeq)
        @test_throws ErrorException(keq_message) rescale_parameter_values(
            mm, nokeq, [(S = 1.0,) => 1.0])
        @test_throws(ErrorException("kcats is empty; pass at least one " *
                                    "(concentrations => turnover) anchor"),
                     rescale_parameter_values(mm, p, Pair{NamedTuple, Float64}[]))
        @test_throws(ErrorException("kcats: X is not a metabolite of the mechanism " *
                                    "(metabolites: S, P)"),
                     rescale_parameter_values(mm, p, [(X = 1.0,) => 1.0]))
        @test_throws(ErrorException("kcats: concentration of S must be a finite number " *
                                    "≥ 0; got -1.0"),
                     rescale_parameter_values(mm, p, [(S = -1.0,) => 1.0]))
        for (bad, shown) in ((0.0, "0.0"), (NaN, "NaN"), (Inf, "Inf"))
            @test_throws(ErrorException("kcats: turnover must be a finite, nonzero " *
                                        "number; got $shown"),
                         rescale_parameter_values(mm, p, [(S = 1.0,) => bad]))
        end
        # No substrate and no product: the predicted turnover is zero.
        @test_throws("an anchor needs a nonzero predicted turnover of the measured " *
                     "turnover's sign",
                     rescale_parameter_values(mm, p, [(S = 0.0,) => 1.0]))
        # The product alone runs the reaction backwards, against a forward measurement.
        @test_throws("an anchor needs a nonzero predicted turnover of the measured " *
                     "turnover's sign",
                     rescale_parameter_values(mm, p, [(P = 1.0,) => 1.0]))
        # Every product of two constants underflows to zero: no peak scale exists.
        tiny = (k_ES_to_E_P = 1e-200, k_ES_to_E_S = 1e-200, k_E_S_to_ES = 1e-200,
                Keq = 2.0)
        @test_throws(ErrorException("the mechanism has no positive forward turnover " *
                                    "within peak_range (1.0e-12, 1000.0); the peak scale " *
                                    "is undefined"),
                     rescale_parameter_values(mm, tiny))
    end
```

In `test/test_rate_eq_derivation.jl`, replace the whole function
`_testhelper_test_kcat_rescaling` with:

```julia
function _testhelper_test_peak_rescaling(spec::_testhelper_MechanismTestSpec; seed=100)
    m = spec.mechanism
    @testset "Peak rescaling" begin
        rng = Random.MersenneTwister(seed)
        params = _testhelper_random_reduced_params(m; rng)
        range = (1e-12, 1e3)
        peak = EnzymeRates._peak_turnover(m, params, range).value
        @test peak > 0

        # Scaling every steady-state rate constant by α scales the peak by α.
        ss_names = EnzymeRates._ss_rate_constant_names(m)
        α = 0.1 + 9.9 * rand(rng)
        scaled_params = NamedTuple{keys(params)}(Tuple(
            k in ss_names ? v * α : v
            for (k, v) in zip(keys(params), values(params))
        ))
        @test EnzymeRates._peak_turnover(m, scaled_params, range).value ≈
            α * peak rtol=1e-9

        # The peak scale makes the peak 1 and leaves every other constant alone.
        norm, spread = rescale_parameter_values(m, params)
        @test spread === nothing
        @test EnzymeRates._peak_turnover(m, norm, range).value ≈ 1.0 rtol=1e-9
        for k in keys(params)
            if !(k in ss_names)
                @test norm[k] == params[k]
            end
        end

        # Every rate changes by the same factor.
        met_names = metabolites(m)
        concs = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        v_orig = rate_equation(m, concs, params)
        @test rate_equation(m, concs, norm) / v_orig ≈ 1.0 / peak rtol=1e-8

        # The peak bounds the rescaled turnover at every corner of the range.
        norm_turnover = merge(norm, (E_total=1.0,))
        corners = Iterators.product(ntuple(_ -> range, length(met_names))...)
        @test all(c -> rate_equation(m, NamedTuple{met_names}(c), norm_turnover) <=
                       1 + 1e-9, corners)

        # An anchor at those concentrations sets the turnover there.
        target = sign(v_orig) * (0.1 + 9.9 * rand(rng))
        anchored, spread_one = rescale_parameter_values(m, params, [concs => target])
        @test spread_one == 0.0
        @test rate_equation(m, concs, merge(anchored, (E_total=1.0,))) ≈ target rtol=1e-9
    end
end
```

and in `_testhelper_run_all_tests` replace `_testhelper_test_kcat_rescaling(spec)` with
`_testhelper_test_peak_rescaling(spec)`.

In `test/allosteric_ground_truth.jl`, replace

```julia
    # kcat stays consistent with rate_equation under normalization.
    prm = NamedTuple{(fp..., :Keq, :E_total)}((1.3, 0.9, 2.1, 0.7, 3.0, 1.0))
    rescaled = ER.rescale_parameter_values(onlyA, prm; scale_k_to_kcat=5.0)
    @test isapprox(ER._kcat_forward(onlyA, rescaled), 5.0; rtol=1e-6)
```

with

```julia
    # The peak scale stays consistent with rate_equation.
    prm = NamedTuple{(fp..., :Keq, :E_total)}((1.3, 0.9, 2.1, 0.7, 3.0, 1.0))
    rescaled, _ = ER.rescale_parameter_values(onlyA, prm)
    @test isapprox(ER._peak_turnover(onlyA, rescaled, (1e-12, 1e3)).value, 1.0;
                   rtol=1e-9)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run the focused commands for `test_fitting_scale` and `test_rate_eq_derivation`.
Expected: `MethodError` for the three-argument `rescale_parameter_values`, and failures
because the two-argument call still returns a NamedTuple, not `(params, nothing)`.

- [ ] **Step 3: Implement**

Delete the keyword method `rescale_parameter_values(m, params; scale_k_to_kcat=1.0)`, its
docstring and its heading comment from `src/rate_eq_derivation.jl`. Append to
`src/fitting.jl`:

````julia
"""
    rescale_parameter_values(m, params; peak_range = (1e-12, 1e3)) -> (params, nothing)
    rescale_parameter_values(m, params, kcats) -> (params, spread)

Re-express the parameters `params` of mechanism `m` on another rate scale by multiplying
every steady-state rate constant (the `k_…` constants) by one factor. The binding
constants `K_…`, `Keq`, `E_total`, `L` and the regulatory constants stay as they are, and
every rate changes by that factor. `params` must include `Keq`, which the dependent rate
constants need; merge it into a fit's `params`.

With `kcats`, a vector of `concentrations => turnover` pairs, the factor matches measured
turnovers — a kcat, or a specific activity in turnover units — at stated concentrations,
given as a `NamedTuple` per anchor; metabolites an anchor leaves out are zero. Several
anchors combine by the geometric mean of measured over predicted turnover, the
least-squares scale in log space. `spread` is the standard deviation of the anchors'
`log(measured/predicted)`, 0 for one anchor; a large spread says the anchors disagree
with the model or with each other. An anchor where the model predicts zero turnover, or
turnover of the other sign, is an error.

Without `kcats`, the factor makes the mechanism's peak turnover — its largest net forward
rate per enzyme over every metabolite's concentration within `peak_range`, in the data's
concentration units — equal 1, the scale `fit_rate_equation` gives relative fits, and
`spread` is `nothing`. The default range covers physiological concentrations in M or mM;
data in µM need a higher upper bound.

```julia
fit = fit_rate_equation(FittingProblem(m, data; Keq = 2.0), optimizer)
kcats = [(A = 10.0, B = 5.0) => 42.0, (A = 20.0, B = 5.0) => 44.0]
params, spread = rescale_parameter_values(m, merge(fit.params, (; Keq = 2.0)), kcats)
```
"""
function rescale_parameter_values(@nospecialize(m::AbstractEnzymeMechanism),
                                  params::NamedTuple; peak_range=DEFAULT_PEAK_RANGE)
    _check_has_keq(params)
    range = _peak_range(peak_range)
    peak = _peak_turnover(m, params, range).value
    0 < peak < Inf || error("the mechanism has no positive forward turnover within " *
                            "peak_range $range; the peak scale is undefined")
    (_scale_rate_constants(m, params, 1 / peak), nothing)
end

function rescale_parameter_values(@nospecialize(m::AbstractEnzymeMechanism),
                                  params::NamedTuple, kcats::AbstractVector{<:Pair})
    _check_has_keq(params)
    isempty(kcats) &&
        error("kcats is empty; pass at least one (concentrations => turnover) anchor")
    mets = metabolites(m)
    turnover_params = merge(params, (E_total = 1.0,))
    log_ratios = map(kcats) do (concs, measured)
        concs isa NamedTuple ||
            error("kcats: concentrations must be a NamedTuple; got $(repr(concs))")
        for (k, c) in pairs(concs)
            k in mets || error("kcats: $k is not a metabolite of the mechanism " *
                               "(metabolites: $(join(mets, ", ")))")
            c isa Real && isfinite(c) && c >= 0 || error(
                "kcats: concentration of $k must be a finite number ≥ 0; got $(repr(c))")
        end
        measured isa Real && isfinite(measured) && measured != 0 || error(
            "kcats: turnover must be a finite, nonzero number; got $(repr(measured))")
        predicted = rate_equation(m,
            NamedTuple{mets}(Tuple(Float64(get(concs, k, 0.0)) for k in mets)),
            turnover_params)
        predicted != 0 && sign(predicted) == sign(measured) || error(
            "kcats: at $concs the mechanism predicts turnover $predicted; an anchor " *
            "needs a nonzero predicted turnover of the measured turnover's sign")
        log(measured / predicted)
    end
    center = sum(log_ratios) / length(log_ratios)
    spread = length(log_ratios) == 1 ? 0.0 :
        sqrt(sum(r -> (r - center)^2, log_ratios) / (length(log_ratios) - 1))
    (_scale_rate_constants(m, params, exp(center)), spread)
end

"""Errors unless `params` holds `Keq`, which the dependent rate constants need."""
_check_has_keq(params::NamedTuple) = haskey(params, :Keq) || error(
    "params must include Keq, which the dependent rate constants need; merge it into a " *
    "fit's params")
````

- [ ] **Step 4: Run the tests to verify they pass**

Run the focused commands for `test_fitting_scale`, `test_rate_eq_derivation` and
`allosteric_ground_truth`. Expected: all pass. `grep -rn "scale_k_to_kcat" src test`
now matches only the two `@test_throws MethodError` lines.

- [ ] **Step 5: Run the full suite, then commit**

```bash
git status
git add src/fitting.jl src/rate_eq_derivation.jl test/test_fitting_scale.jl \
  test/test_rate_eq_derivation.jl test/allosteric_ground_truth.jl
git commit -m "Rescale a fit to measured turnovers or to the peak scale after the fit" \
  -m "rescale_parameter_values(m, params; peak_range) applies the peak-turnover scale;
rescale_parameter_values(m, params, kcats) matches measured turnovers at stated
concentrations by their geometric mean and returns the spread of the anchors." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 5: Delete the analytic kcat

**Files:**
- Modify: `src/rate_eq_derivation.jl` (delete `_kcat_groups_from_polys`, both
  `_kcat_forward` methods, `_kcat_keys`, `_max_expr`, `_mwc_power_pair`, the heading
  `# ─── kcat Computation Helpers ───...` (rename it
  `# ─── Steady-state rate-constant names ───...` above `_ss_rate_constant_names`);
  update the ABOUTME line and two docstrings)
- Modify: `test/test_rate_eq_derivation.jl`, `test/allosteric_ground_truth.jl`,
  `test/test_mechanism_enumeration.jl`,
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `docs/src/developer.md`, `ROADMAP.md`

**Interfaces:**
- Consumes: `_peak_turnover` (Task 1).
- Produces: no `_kcat_forward` anywhere in `src/`, `test/` or `docs/src/`.

This task rewrites every test that calls `_kcat_forward` so it tests the peak turnover
on the same mechanism (CLAUDE.md: never delete a test; replace a removed behavior's test
with one of the behavior that replaces it).

- [ ] **Step 1: Rewrite the tests (they must still pass against the current code)**

In `test/test_rate_eq_derivation.jl`:

1. Line 2 of the file: replace `# ABOUTME: the parameters API, kcat/rescaling,
rate-equation strings, and the perf gate.` with `# ABOUTME: the parameters API, peak
rescaling, rate-equation strings, and the perf gate.`

2. Replace `_testhelper_test_analytical_kcat` with:

```julia
function _testhelper_test_analytical_kcat(spec::_testhelper_MechanismTestSpec; seed=42)
    spec.analytical_kcat_fn === nothing && return
    m = spec.mechanism
    @testset "Peak turnover approaches the analytical kcat" begin
        rng = Random.MersenneTwister(seed)
        params = _testhelper_random_reduced_params(m; rng)
        # Every constant lies in 0.1–10, so over a range reaching 1e12 the peak sits at
        # saturating substrates and zero products and equals kcat to within K/1e12.
        peak = EnzymeRates._peak_turnover(m, params, (1e-12, 1e12)).value
        # The oracle's positional formula may reference a forward rate that is
        # Haldane-DEPENDENT under the canonical step order (absent from the
        # reduced params), so bridge the FULL param set (dependent values
        # included) to positional names.
        p = merge(_testhelper_analytical_oracle_params(
                      m, _testhelper_compute_all_params(m, params);
                      source_steps=spec.source_steps,
                      source_reg_sites=spec.source_reg_sites),
                  (Et=params.E_total,))
        @test peak ≈ spec.analytical_kcat_fn(p) rtol=1e-9
    end
end
```

3. Line ~1448, in the testset "§5a I-state parameters are all defined (regression)",
replace the comment `# products = 0 (the kcat evaluation domain).` with
`# products = 0 (initial-rate conditions).`

4. Replace the heading `# ── Standalone kcat tests ──...` with
`# ── Standalone tests ──...` (keep the line's length).

5. Replace the testset "Fix B: _kcat_forward handles multiple saturating patterns" with
(same mechanism, now an oracle test of the peak search):

```julia
@testset "the peak of an allosteric random-order bi-bi reaches the grid maximum" begin
    m = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A <--> E(A)        :: NonequalAI
            E + B <--> E(B)        :: NonequalAI
            E(A) + B <--> E(A, B)  :: NonequalAI
            E(B) + A <--> E(A, B)  :: NonequalAI
            E(A, B) <--> E(P, Q)   :: NonequalAI
            E(P, Q) <--> E(P) + Q  :: NonequalAI
            E(P, Q) <--> E(Q) + P  :: NonequalAI
            E(P) <--> E + P        :: NonequalAI
            E(Q) <--> E + Q        :: NonequalAI
        end
    end
    rng = Random.MersenneTwister(1)
    pn = EnzymeRates.fitted_params(m)
    pv = NamedTuple{pn}(Tuple(0.2 + 2 * rand(rng) for _ in pn))
    peak = EnzymeRates._peak_turnover(m, merge(pv, (Keq = 1.0,)), (1e-12, 1e9))
    @test isfinite(peak.value)
    fp = merge(pv, (Keq = 1.0, E_total = 1.0))
    vgrid = maximum(rate_equation(m, (A = x, B = y, P = p, Q = q), fp)
                    for x in 10.0 .^ (-12:9), y in 10.0 .^ (-12:9),
                        p in 10.0 .^ (-12:3:9), q in 10.0 .^ (-12:3:9))
    @test peak.value >= vgrid * (1 - 1e-12)
end
```

6. In the testset "allosteric kcat: a one-sided regulator's corner is its
conformation's limit": rename it
`"the peak turnover covers every regulator corner at saturating substrate"`; in its
opening comment replace the last sentence ("kcat is the largest regulator corner at
saturating S.") with "The peak turnover is at least every regulator corner's turnover at
saturating S."; keep the four mechanisms (`only_a`, `only_i`, `only_i_dead`, `two_site`)
and their comments; and replace the loop from `BIG = 1e9` to the testset's last `end`
with:

```julia
    BIG = 1e9
    for m in (only_a, only_i, only_i_dead, two_site)
        pn = EnzymeRates.fitted_params(m)
        regs = Tuple(x for x in metabolites(m) if x ∉ (:S, :P))
        for seed in 1:5
            rng = Random.MersenneTwister(seed)
            p = merge(NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn)),
                      (Keq = 2.0, E_total = 1.0))
            # Corners of the searched range: each regulator at a bound, saturating S,
            # the product at the lower bound.
            vmax = maximum(0:(2^length(regs) - 1)) do mask
                r = Tuple((mask >> (i - 1)) & 1 == 1 ? BIG : 1e-12
                          for i in eachindex(regs))
                rate_equation(m, NamedTuple{(:S, :P, regs...)}((BIG, 1e-12, r...)), p)
            end
            @test EnzymeRates._peak_turnover(m, p, (1e-12, BIG)).value >=
                  vmax * (1 - 1e-12)
        end
    end
end
```

7. Replace the testset "Fix B: non-allosteric random-order bi-bi kcat = peak (contract
guard)" with:

```julia
@testset "the peak turnover of a random-order bi-bi reaches the grid maximum" begin
    # 40 parameter draws; the grid spans every metabolite, products included.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E(A) + B <--> E(A, B)
            E(B) + A <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(P) + Q
            E(P, Q) <--> E(Q) + P
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end
    pn = EnzymeRates.fitted_params(m)
    for seed in 1:40
        rng = Random.MersenneTwister(seed)
        pv = NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn))
        peak = EnzymeRates._peak_turnover(m, merge(pv, (Keq = 1.0,)), (1e-12, 1e3))
        fp = merge(pv, (Keq = 1.0, E_total = 1.0))
        vgrid = maximum(rate_equation(m, (A = x, B = y, P = p, Q = q), fp)
                        for x in 10.0 .^ (-12:3), y in 10.0 .^ (-12:3),
                            p in 10.0 .^ (-12:3:3), q in 10.0 .^ (-12:3:3))
        @test peak.value >= vgrid * (1 - 1e-12)
    end
end
```

8. Replace the testset "kcat with no product-free saturating pattern errors
descriptively" with (same mechanism):

```julia
@testset "the peak of a law that needs bound product lies at a positive product level" begin
    # Catalysis runs only on the product-bound form E(P), so every forward numerator
    # term carries P and the turnover vanishes at zero product.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
            E(P, S) <--> E(P, P)
            E(P, P) ⇌ E(P) + P
        end
    end
    pn = EnzymeRates.fitted_params(m)
    params = merge(NamedTuple{pn}(ntuple(_ -> 1.5, length(pn))), (Keq = 2.0,))
    peak = EnzymeRates._peak_turnover(m, params, (1e-12, 1e3))
    @test peak.value > 0
    @test 1e-3 < peak.concs.P < 1e3
end
```

In `test/allosteric_ground_truth.jl`, replace each of the five lines
`@test isfinite(ER._kcat_forward(X, prm))` (with `X` = `onlyA`, `multiA`, `metabD`
twice, `err1`) by `@test 0 < ER._peak_turnover(X, prm, (1e-12, 1e3)).value < Inf`; in the
testset "dead-inactive ping-pong :OnlyA derivation" replace

```julia
        kcat = ER._kcat_forward(dead, prm)
        @test isfinite(kcat) && kcat > 0
```

with

```julia
        @test 0 < ER._peak_turnover(dead, prm, (1e-12, 1e3)).value < Inf
```

and in the comment above that testset replace "and yields a finite positive kcat." with
"and has a finite positive peak turnover."

In `test/test_mechanism_enumeration.jl` (testset "change_allo_state drops mixed-iso
V-type relaxations (no :OnlyA binding)"), replace

```julia
        # end-to-end: every surviving child's saturating turnover is finite —
        # no relaxation reintroduces the covalent sink that crashes kcat.
        cm = EnzymeRates.compile_mechanism(k)
        fp = EnzymeRates.fitted_params(cm)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
        @test isfinite(EnzymeRates._kcat_forward(cm, prm))
```

with

```julia
        # end-to-end: every surviving child turns over without product, so no
        # relaxation reintroduces a covalent sink.
        cm = EnzymeRates.compile_mechanism(k)
        fp = EnzymeRates.fitted_params(cm)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
        @test EnzymeRates.rate_equation(cm, (ATP = 1.3, F6P = 1.3, ADP = 0.0,
                                             F16BP = 0.0), prm) > 0
```

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`:
- replace the field comment `# Optional analytical kcat function for kcat/rescaling
validation` with `# Optional analytical kcat: the turnover at saturating substrates
and zero products, which the peak turnover approaches on a wide range`;
- in the PFK-1 spec, replace the comment's last three lines (from
  `# for every regulator corner. Guards `_kcat_forward` for a saturating`) with
  `# for every regulator corner. Guards the peak turnover for a saturating` and
  `# pattern absent from the I-state.`; in the HK spec, make the same replacement (its
  comment says `kcat = k6f` on the line before).

- [ ] **Step 2: Run the rewritten tests against the current code**

Run the focused commands for `test_rate_eq_derivation`, `allosteric_ground_truth` and
`test_mechanism_enumeration` (the last takes about 3 minutes). Expected: all pass — they
no longer call `_kcat_forward`, which still exists.

- [ ] **Step 3: Delete the analytic kcat**

In `src/rate_eq_derivation.jl`:
- delete `_kcat_groups_from_polys`, both `_kcat_forward` methods with their docstrings,
  `_kcat_keys`, `_max_expr` and `_mwc_power_pair` with its docstring;
- rename the heading above `_ss_rate_constant_names` to
  `# ─── Steady-state rate-constant names ───` padded like its neighbours;
- ABOUTME line 2: `# ABOUTME: Parameters API and MWC allosteric rate assembly.`;
- `_mwc_combine`'s docstring: `"""MWC active + L·inactive state-combine `A + L * B` of
  the retained num/den sums (`_num_den_exprs`)."""`;
- `_mwc_state_polys`'s docstring: replace "that `rate_equation` (`_num_den_exprs`) and
  `_kcat_forward` share" with "that `rate_equation` (`_num_den_exprs`) uses", and delete
  its last two sentences ("The normalization is a common factor of kcat's
  saturating-limit ratio, ... consistent with `rate_equation`.").

Then run
`grep -rn "_kcat_forward\|_kcat_keys\|_kcat_groups\|_max_expr\|_mwc_power_pair" src test docs/src`.
Expected: no match (the test helper `_testhelper_test_analytical_kcat` and the spec field
`analytical_kcat_fn` stay).

In `docs/src/developer.md`, replace "Only the `@generated` methods (`rate_equation`,
`parameters`, `fitted_params`, `metabolites`, `_kcat_forward`) do work per singleton type:
`rate_equation` and `_kcat_forward` run the derivation," with "Only the `@generated`
methods (`rate_equation`, `parameters`, `fitted_params`, `metabolites`) do work per
singleton type: `rate_equation` runs the derivation," and rewrap the paragraph to 92
characters. At the end of the same paragraph (after "once per optimizer, not once per
mechanism."), add: "The peak search that scales a relative fit (`_peak_turnover`)
wraps the mechanism's turnover the same way, so its search loop compiles once."

In `ROADMAP.md`, delete the Known-issues entry "**`_kcat_forward` reports "no kcat
components" for mechanisms that never saturate.**" with its paragraph.

- [ ] **Step 4: Run the full suite**

Expected: all pass (Aqua and JET included).

- [ ] **Step 5: Commit**

```bash
git status
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl \
  test/allosteric_ground_truth.jl test/test_mechanism_enumeration.jl \
  test/mechanism_definitions_for_test_enzyme_derivation.jl docs/src/developer.md ROADMAP.md
git commit -m "Delete the analytic kcat; test the peak turnover in its place" \
  -m "Every test of _kcat_forward now tests the peak turnover on the same mechanism:
the analytic kcats against a wide range, the saturating-pattern cases against grid
maxima, and the guards as finite positive peaks." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 6: Gate laws above 400 denominator terms by a bound proportional to their size

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (after the testset "rate_equation stays
  allocation-free and under 120 ns on the largest laws")

**Interfaces:**
- Consumes: `_testhelper_test_rate_equation_performance`,
  `_testhelper_random_independent_params_concs` (same file),
  `EnzymeRates._raw_symbolic_rate_polys`.
- Produces: `_testhelper_large_ter_bi_mechanisms`,
  `_testhelper_seconds_per_denominator_term`, `_testhelper_denominator_terms`.

The spec (section 11, decision E1) keeps the 120-ns gate for the test specs and the three
largest bi-bi laws (396 denominator terms) and bounds larger laws by a constant times
their number of denominator terms, the constant calibrated on CI hardware in this phase.
Locally (aarch64) the three laws below take 134–188 ns, 0.19–0.24 ns per term.

- [ ] **Step 1: Write the gate with a provisional constant**

Add after the testset "rate_equation stays allocation-free and under 120 ns on the
largest laws":

```julia
# Three random-order seeds of A + B + C ⇌ P + Q, larger than the bi-bi laws above: 589,
# 799 and 858 denominator terms.
const _testhelper_large_ter_bi_mechanisms = [
    @enzyme_mechanism(begin substrates: A, B, C; products: P, Q; steps: begin
        (E + A <--> E(A), E(B) + A <--> E(A, B), E(B, C) + A <--> E(P, Q),
         E(C) + A <--> E(A, C), E(P) + A <--> E(A, P))
        (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(A, C) + B ⇌ E(P, Q), E(C) + B ⇌ E(B, C),
         E(C, Q) + B ⇌ E(B, C, Q), E(Q) + B ⇌ E(B, Q))
        (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(A, B) + C ⇌ E(P, Q), E(B) + C ⇌ E(B, C),
         E(B, Q) + C ⇌ E(B, C, Q), E(Q) + C ⇌ E(C, Q))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(B, C) + Q <--> E(B, C, Q),
         E(C) + Q <--> E(C, Q), E(P) + Q <--> E(P, Q)) end end),
    @enzyme_mechanism(begin substrates: A, B, C; products: P, Q; steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(B, C) + A ⇌ E(P, Q), E(C) + A ⇌ E(A, C),
         E(P) + A ⇌ E(A, P))
        (E + B <--> E(B), E(A) + B <--> E(A, B), E(A, C) + B <--> E(P, Q),
         E(C) + B <--> E(B, C), E(C, Q) + B <--> E(B, C, Q), E(Q) + B <--> E(B, Q))
        (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(A, B) + C ⇌ E(P, Q), E(B) + C ⇌ E(B, C),
         E(B, Q) + C ⇌ E(B, C, Q), E(Q) + C ⇌ E(C, Q))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(B, C) + Q <--> E(B, C, Q),
         E(C) + Q <--> E(C, Q), E(P) + Q <--> E(P, Q)) end end),
    @enzyme_mechanism(begin substrates: A, B, C; products: P, Q; steps: begin
        (E + A <--> E(A), E(C) + A <--> E(Q; residual = A + C - Q),
         E(C, P) + A <--> E(A, C, P), E(P) + A <--> E(A, P))
        (E + C <--> E(C), E(A) + C <--> E(Q; residual = A + C - Q),
         E(A, P) + C <--> E(A, C, P), E(P) + C <--> E(C, P))
        (E + P <--> E(P), E(A) + P <--> E(A, P), E(C) + P <--> E(C, P),
         E(Q; residual = A + C - Q) + P <--> E(A, C, P))
        (E(P) + Q ⇌ E(B, Q; residual = A + C - Q),
         E(; residual = A + C - Q) + Q ⇌ E(Q; residual = A + C - Q))
        (E(Q; residual = A + C - Q) + B ⇌ E(B, Q; residual = A + C - Q),
         E(; residual = A + C - Q) + B ⇌ E(P)) end end),
]

# Seconds per call allowed per denominator term above the 120-ns gate. Provisional:
# twice the worst local measurement; Step 3 calibrates it on the CI runners.
const _testhelper_seconds_per_denominator_term = 0.5e-9

"""The number of terms in the denominator of `m`'s derived rate law."""
_testhelper_denominator_terms(m) =
    length(EnzymeRates._raw_symbolic_rate_polys(typeof(m))[2])

@testset "rate_equation time per call grows at most in proportion to the law's size" begin
    @test _testhelper_denominator_terms.(_testhelper_large_ter_bi_mechanisms) ==
          [589, 799, 858]
    for m in _testhelper_large_ter_bi_mechanisms
        params, concs, _ = _testhelper_random_independent_params_concs(
            m, collect(metabolites(m)); rng = Random.MersenneTwister(42))
        allocs, t = _testhelper_test_rate_equation_performance(m, params, concs)
        n = _testhelper_denominator_terms(m)
        @info "rate_equation: $(round(1e9 * t; digits = 1)) ns per call, $n denominator " *
              "terms, $(round(1e9 * t / n; digits = 3)) ns per term"
        @test allocs == 0
        @test t < max(120e-9, _testhelper_seconds_per_denominator_term * n)
    end
end
```

- [ ] **Step 2: Run it locally**

Run the focused command for `test_rate_eq_derivation`. Expected: passes; the three
`@info` lines show about 0.2 ns per term. If the term counts differ from
`[589, 799, 858]`, the DSL text above does not match the seeds; stop and report.

- [ ] **Step 3: Calibrate the constant on the CI runners**

1. Commit (message "Gate rate_equation's time per call by the size of large laws", body
   "Laws above the three largest bi-bi laws get a bound proportional to their number of
   denominator terms." and the two trailer lines), then `git push -u origin fitting-scale`.
2. Trigger CI on the branch: `gh workflow run CI.yml --ref fitting-scale`; find the run
   with `gh run list --workflow CI.yml --branch fitting-scale --limit 1`; poll
   `gh run view <id> --json status,conclusion` in a bounded loop (sleep 60, at most 60
   minutes) until it completes.
3. Read the measurements: `gh run view <id> --log | grep "denominator terms"`. Expected:
   three lines per test job (ubuntu x64, windows x64, macOS aarch64). If a job failed on
   this gate, read its measured values; the calibration below fixes the constant.
4. Set `_testhelper_seconds_per_denominator_term` to twice the largest ns-per-term value
   across all jobs, rounded up to one significant figure (in seconds, e.g. `0.6e-9`), and
   replace its comment with: `# Seconds per call allowed per denominator term above the
   120-ns gate: twice the worst CI measurement, <x> ns per term (<runner>, <date>), so the
   bound stays above shared runners' noise.` Also record, in that comment, the worst
   value of each runner.
5. Run the focused command for `test_rate_eq_derivation` again (passes), commit
   ("Calibrate the size-scaled speed gate on the CI runners", trailers), push, and
   trigger CI again; confirm every test job passes.

---

### Task 7: Document the peak scale, absolute mode and rescaling; fix the `maxtime` docstring

**Files:**
- Modify: `docs/src/fitting/normalized_vs_absolute.md` (rewrite),
  `docs/src/fitting/tutorial.md`, `docs/src/fitting/loss_and_optimizers.md`,
  `docs/src/identify/model_selection.md`, `src/identify_rate_equation.jl` (docstring),
  `ROADMAP.md`

No test changes: documentation. Every `@example` block must run.

- [ ] **Step 1: Rewrite `docs/src/fitting/normalized_vs_absolute.md`**

Replace the whole file with:

````markdown
# Normalized vs absolute rate

Enzyme kinetic data is often reported on an arbitrary rate scale. The enzyme
concentration in the assay may be unknown, and figures from different
preparations may reflect different states of purity, so the absolute value of a
measured rate carries no reliable meaning. What survives is the *shape* of the
response: changing the amount of enzyme scales every rate by the same factor and
leaves that shape intact.

EnzymeRates.jl handles both cases through one `FittingProblem` keyword, `absolute`.
With the default `absolute = false` the data are relative: the fit ignores an
arbitrary factor on each group's rates, and the returned rate constants are scaled so
the mechanism's peak turnover is 1. With `absolute = true` the data are per-enzyme
turnover: the fit scores absolute magnitudes and returns the constants as fitted.
`IdentifyRateEquationProblem` takes the same keyword.

## The two modes

### Relative mode (the default)

- **Loss**: each group's log-ratios are **mean-centered before squaring**,
  which removes the arbitrary per-group `E_total` scale. The loss is invariant
  to rescaling all rates in a group by the same factor.
- **Output scale**: the loss cannot see a common factor on every steady-state rate
  constant, so the fitter fixes it after the fit: the returned constants make the
  mechanism's **peak turnover** equal 1. The peak turnover is the largest net forward
  rate per enzyme, `v/E_total`, over every metabolite's concentration — substrates,
  products, inhibitors and regulators — within `peak_range`. A deterministic search in
  log concentrations finds it, the range's corners included. Multiply the returned rate
  constants by a measured peak turnover to put them in absolute units, or rescale them
  to measured kcats (below).

Relative mode is the right choice when your data reports relative rates
(arbitrary units, percent activity, or values normalized per some `E_total`
estimate that you do not trust).

### Absolute mode

Pass `absolute = true`. In absolute mode:

- **Loss**: log-ratios are squared **without centering**. The absolute
  magnitude of the rate is scored, so the data must report per-enzyme turnover
  in consistent units.
- **Output**: the raw fitted parameters, never rescaled. The data fix the absolute
  scale.

Absolute mode is the right choice when your data is already on a per-enzyme
turnover basis.

!!! note "The fitter has no `E_total` input"
    Every prediction is evaluated at `E_total = 1`, so the loss compares your
    measured `Rate` against the per-enzyme turnover — in effect it fits
    `Rate / E_total`. The data table carries no `E_total` column, and the fitter
    ignores one if present. Relative mode makes this moot: per-group centering
    cancels any constant factor, `E_total` included. Absolute mode does not, so
    the measured `Rate` must already be per-enzyme — divide rates recorded at a
    known enzyme amount by that amount before fitting.

## The peak range

`peak_range = (lo, hi)` bounds the concentrations the peak search visits, in the units
of your data. The default `(1e-12, 1e3)` covers physiological concentrations whether
your data are in M, the unit [Fitting tutorial & data format](@ref) recommends, or in
mM. For data in µM, where `1e3` means 1 mM, raise the upper bound. `FittingProblem`,
`IdentifyRateEquationProblem` and `rescale_parameter_values` all take `peak_range`.

Most laws peak inside the range or saturate well below its upper edge. A law still
rising at the upper edge — by more than 10% over the range's last decade — is flagged:
the fit result's `peak_at_range_edge` is `true`, and identify's CSVs record the flag in
a column of the same name. The fit is still scaled, to the turnover at the edge; the
flag says that a wider range would change the scale.

## Rescaling after the fit

`rescale_parameter_values` re-expresses fitted parameters on another scale. It rescales
only the rate constants — the lowercase-`k` parameters every steady-state step
contributes, in both directions (`k_<from>_to_<to>`) — and leaves the binding constants
`K`, `Keq`, `E_total`, the allosteric `L`, and regulatory K values unchanged.

The rate constants are the only parameters that carry time in their units, and
the measured rate is the only thing that supplies a time scale, so it can pin
down only what lives in them. The binding constants, `Keq`, and `E_total` are
concentrations or concentration ratios, and `L` is a pure number — none carry a
time dimension, so the rate constrains them only as far as the response shape
already does. Multiplying every rate constant by a common factor multiplies every
rate by that factor, and the K's stay put: a rapid-equilibrium step stores its binding
constant directly, so the rescaling never touches it, while a steady-state binding
constant is the ratio of its two rate constants, in which the common factor cancels.

The parameters you pass must include `Keq`, which the dependent rate constants need; a
fit's `params` holds only the fitted constants, so merge it in.

- **Measured turnovers.** Pass `kcats`, a vector of `concentrations => turnover` pairs,
  each a measured turnover — a kcat or a specific activity in turnover units — at stated
  concentrations. Metabolites you leave out are zero. One anchor makes the predicted
  turnover equal the measured one there; several combine by the geometric mean of
  measured over predicted turnover. The second return value, the spread, is the standard
  deviation of the anchors' `log(measured/predicted)`: 0 for one anchor, and large when
  the anchors disagree with the model or with each other.
- **Peak scale.** Without `kcats`, the function applies the scale relative fits get:
  peak turnover 1 within `peak_range`. Its second return value is `nothing`.

```@example rescale
using EnzymeRates

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

params = (k_ES_to_E_P = 3.0, k_ES_to_E_S = 1.0, k_E_S_to_ES = 4.0, Keq = 2.0)

# A turnover of 42 measured at S = 1000 with no product
anchored, spread = rescale_parameter_values(uni_uni, params, [(S = 1e3,) => 42.0])
anchored
```

The rescaled law reproduces the anchor:

```@example rescale
rate_equation(uni_uni, (S = 1e3, P = 0.0), merge(anchored, (E_total = 1.0,)))
```

The peak scale, which relative fits return:

```@example rescale
peak_scaled, _ = rescale_parameter_values(uni_uni, params)
peak_scaled
```

Only the rate constants changed; `Keq` is the same in both.
````

- [ ] **Step 2: Update the other pages and the docstring**

1. `docs/src/fitting/tutorial.md`:
   - in the note "No `E_total` column", replace "For absolute mode (`scale_k_to_kcat =
     nothing`)," with "For absolute mode (`absolute = true`),";
   - in the paragraph listing construction errors, replace "a `Keq` that is not positive
     and finite, or a `scale_k_to_kcat` that is neither `nothing` nor positive and finite
     each raises" with "a `Keq` that is not positive and finite, or a `peak_range` that
     is not a tuple `(lo, hi)` with `0 < lo < hi < Inf` each raises" (rewrap);
   - replace the last paragraph (from "The returned `params` reflect kcat normalization"
     to the end) with:

     ```markdown
     The returned `params` are scaled so the mechanism's peak turnover — its largest net
     forward rate per enzyme over the concentrations in `peak_range` — is 1, and
     `result.peak_at_range_edge` says whether that peak sits at the range's upper edge.
     [Normalized vs absolute rate](@ref) explains the scale, absolute mode, and how to
     rescale a fit to a measured kcat; [Loss & optimizers](@ref) covers the loss
     function and how to supply your own optimizer.
     ```
2. `docs/src/fitting/loss_and_optimizers.md`, section "Centered vs uncentered": replace
   "The loss formula switches on `fp.scale_k_to_kcat`" with "The loss formula switches on
   `fp.absolute`", "**Relative mode** (`scale_k_to_kcat` is a `Real`)" with "**Relative
   mode** (`absolute = false`, the default)", and "**Absolute mode**
   (`scale_k_to_kcat === nothing`)" with "**Absolute mode** (`absolute = true`)".
3. `docs/src/identify/model_selection.md`, the `cv_results` columns table: after the
   `fit_inherited` row add
   `| `peak_at_range_edge` | `true` when the fitted law still rises at the upper edge of `peak_range` (by more than 10% over its last decade), so a wider range would change the scale of the reported rate constants; always `false` in absolute mode. |`
4. `src/identify_rate_equation.jl`, `identify_rate_equation`'s docstring: replace the
   `maxtime` bullet with:

   ```
   - `maxtime::Real = 600.0`: max time per optimizer run in seconds (common solver
     option, forwarded to `Optimization.solve`); each of a fit's `n_restarts` runs
     gets it, so a fit can take up to `n_restarts × maxtime`
   ```
5. `ROADMAP.md`: delete the Known-issues entry "**`identify_rate_equation`'s docstring
   calls `maxtime` the time per fit.**" with its paragraph.

- [ ] **Step 3: Run every `@example` block of the rewritten page**

Write the three `rescale` blocks into one script in the scratchpad and run it with
`$JL <script>`. Expected: it runs; `anchored` has `k_*` values scaled by the same factor
and `Keq = 2.0`; the second block prints `42.0` (to rounding); `peak_scaled` has `Keq =
2.0`. Then grep `grep -rn "scale_k_to_kcat\|_kcat_forward" docs/src src README.md`.
Expected: no match.

- [ ] **Step 4: Commit**

```bash
git status
git add docs/src/fitting/normalized_vs_absolute.md docs/src/fitting/tutorial.md \
  docs/src/fitting/loss_and_optimizers.md docs/src/identify/model_selection.md \
  src/identify_rate_equation.jl ROADMAP.md
git commit -m "Document the peak scale, absolute mode and rescaling after the fit" \
  -m "Also state that maxtime bounds each optimizer run, so a fit can take up to
n_restarts × maxtime." \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
```

---

### Task 8: Propose v0.10.0, verify, open the PR

**Files:**
- Modify: `Project.toml` (`version = "0.9.1"` → `version = "0.10.0"`)

- [ ] **Step 1: Bump the version**

Edit `Project.toml`, then run the full suite. Expected: all pass.

- [ ] **Step 2: Commit and push**

```bash
git status
git add Project.toml
git commit -m "Propose v0.10.0: the fitting scale changes the API" \
  -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE"
git push
```

- [ ] **Step 3: Open the PR (draft, stacked)**

Open a draft PR with base `main` and head `fitting-scale`, titled "Scale relative fits to
peak turnover and apply measured kcats after the fit (v0.10.0)". Body sections:
**Why** (the analytic kcat is 0, infinite or below the peak for valid mechanisms;
spec section 12), **What changes** (one bullet per Task 1–7), **Breaking changes**
(`scale_k_to_kcat` → `absolute`; `rescale_parameter_values` returns `(params, spread)`
and takes `kcats`; `fit_rate_equation` gains `peak_at_range_edge`; CSV column;
`_kcat_forward` removed), **Stacking** (based on the merge of #85, #86, #87, #88 and #89;
the commits of this PR are those after the merge commit
`Merge branch 'rewrite-spec-roadmap' into fitting-scale (stack base)`; it will be
rebased onto main once they merge), **Local HPC scripts** (the git-ignored
`docs/hpc_results/*/identify_*.jl` pass `scale_k_to_kcat=1.0`, the old default: delete
that keyword when running this version), **Version** (v0.10.0 proposed for Denis to
decide), **Tests** and **Verification** (full-suite count and time, CI results, the
calibrated speed-gate values). End with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`;
no session link.

- [ ] **Step 4: Wait for CI and report**

Poll `gh pr checks <n>` in a bounded loop (sleep 60, at most 60 minutes). Expected:
every job passes, the docs job included (it runs the `@example` blocks). If a job fails,
find the root cause with superpowers:systematic-debugging, fix it test-first, push, and
poll again.
