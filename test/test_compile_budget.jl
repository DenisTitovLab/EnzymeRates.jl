# ABOUTME: Compile-time regression gates for the EnzymeRates pipeline.
# ABOUTME: One trace-compile (init_mechanisms) + one wall-clock (rate_equation body-build).

using Test
using EnzymeRates

# Budgets calibrated against current main with 2× headroom (Step 4).
# Per spec §3: this budget is FIXED for the refactor's duration.
# Future stages REDESIGN on trip rather than recalibrate.
const INIT_TRACE_BUDGET                  = 58    # baseline 2026-05-20: 29; budget = 2× (rounded up)
const RATE_EQUATION_WALLCLOCK_BUDGET_S   = 2.1   # baseline 2026-05-20: 1.03s; budget = 2×

# Warmup-reuse: ratio + absolute ceiling on t_warm.
# Baseline 2026-05-20: t_cold=8.4s, t_warm=6.9s, ratio=0.82; budgets = 2×.
# baseline_ratio ≥ 0.5 → ratio gate effectively off (capped at 1.0);
# the absolute gate is what catches regressions in this case.
const WARMUP_TER_RATIO_MAX               = 1.0   # capped (2×0.82 > 1.0)
const WARMUP_T_WARM_TER_MAX_S            = 13.7  # 2× baseline t_warm

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

# Runs `script` in a fresh Julia subprocess; the script is expected to
# print a line `ELAPSED:<float>` (a single @elapsed measurement). Returns
# the parsed Float64, or NaN on any failure.
function _measure_elapsed_subprocess(script::String)
    julia_exe = Base.julia_cmd().exec[1]
    out_buf = IOBuffer()
    try
        run(pipeline(Cmd([julia_exe, "--project=.", "-e", script]);
                     stdout=out_buf, stderr=devnull); wait=true)
    catch e
        @warn "@elapsed subprocess failed: $e"
        return NaN
    end
    out = String(take!(out_buf))
    m = match(r"ELAPSED:([0-9.eE+-]+)", out)
    m === nothing ? NaN : parse(Float64, m.captures[1])
end

@testset "compile-budget" begin
    # Trace-compile: init_mechanisms on a bi-bi reaction (non-trivial so
    # the gate is representative; uni-uni is too small to catch regressions).
    @testset "trace-compile: init_mechanisms (bi-bi)" begin
        # Constructs EnzymeReactionLegacy directly so the init_mechanisms
        # trace-compile gate measures only the enumeration pipeline,
        # independent of DSL grammar changes.
        script = """
            using EnzymeRates
            r = EnzymeRates.EnzymeReactionLegacy(
                ((:A, ((:C, 1),)), (:B, ((:N, 1),))),
                ((:P, ((:C, 1),)), (:Q, ((:N, 1),))),
            )
            EnzymeRates.init_mechanisms(r)
            """
        n = _count_relevant_precompiles(script)
        @info "init_mechanisms trace-compile: $n (budget: $INIT_TRACE_BUDGET)"
        @test 0 <= n <= INIT_TRACE_BUDGET
    end

    # Wall-clock: rate_equation body-build (first call pays @generated cost).
    # MUST be measured in a fresh subprocess — runtests.jl runs test_dsl.jl
    # earlier, which JIT-builds the same EnzymeMechanism{...} body, so an
    # in-process measurement here would always be ~0.01s and the gate
    # would never trip on a regression. The per-call runtime gate is
    # separately enforced by test_rate_eq_derivation.jl's
    # test_rate_equation_performance (0 allocs, <100ns per call) for every
    # mechanism in MECHANISM_TEST_SPECS.
    @testset "wall-clock: rate_equation body-build (first call)" begin
        script = """
            using EnzymeRates
            m = @enzyme_mechanism begin
                substrates: S
                products:   P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                end
            end
            params = NamedTuple{Tuple(EnzymeRates.parameters(m))}(
                ntuple(_ -> 1.0, length(EnzymeRates.parameters(m))))
            concs = (S = 1.0, P = 0.5)
            t = @elapsed EnzymeRates.rate_equation(m, concs, params)
            println("ELAPSED:", t)
            """
        julia_exe = Base.julia_cmd().exec[1]
        out_buf = IOBuffer()
        try
            run(pipeline(Cmd([julia_exe, "--project=.", "-e", script]);
                         stdout=out_buf, stderr=devnull); wait=true)
        catch e
            @warn "Wall-clock subprocess failed: $e"
            @test false
            return
        end
        out = String(take!(out_buf))
        m_match = match(r"ELAPSED:([0-9.eE+-]+)", out)
        t_first = m_match === nothing ? NaN : parse(Float64, m_match.captures[1])
        @info "rate_equation first-call wall-clock: $(t_first)s " *
              "(budget: $RATE_EQUATION_WALLCLOCK_BUDGET_S s)"
        @test isfinite(t_first)
        @test t_first < RATE_EQUATION_WALLCLOCK_BUDGET_S
    end

    # Warmup-reuse regression: time init_mechanisms(r_ter) in two scenarios,
    # each in its own fresh Julia subprocess to avoid cross-test pollution:
    #
    #   t_cold: using EnzymeRates -> init_mechanisms(r_ter)
    #   t_warm: using EnzymeRates -> init_mechanisms(r_uni)  [warmup]
    #                              -> init_mechanisms(r_ter)  [@elapsed measured]
    #
    # If uni-uni warmup shares most specializations with ter-ter, t_warm
    # should be substantially less than t_cold. The parametric
    # EnzymeReactionLegacy{S,P,R,N} forces per-arity specialization;
    # the gate exists to catch a refactor commit that introduces NEW per-
    # arity specialization beyond today's baseline.
    #
    # Single-subprocess @elapsed per scenario avoids the subtraction-noise
    # problem the trace-compile-count approach had (small delta of two
    # larger noisy numbers).
    @testset "warmup-reuse: ter-ter post-warmup wall-clock bounded vs cold" begin
        # Constructs EnzymeReactionLegacy directly so the warmup-reuse
        # gate measures only the enumeration pipeline, independent of
        # DSL grammar changes.
        cold_script = """
            using EnzymeRates
            r_ter = EnzymeRates.EnzymeReactionLegacy(
                ((:A, ((:C, 1),)), (:B, ((:N, 1),)), (:C, ((:O, 1),))),
                ((:P, ((:C, 1),)), (:Q, ((:N, 1),)), (:R, ((:O, 1),))),
            )
            t = @elapsed EnzymeRates.init_mechanisms(r_ter)
            println("ELAPSED:", t)
            """
        warm_script = """
            using EnzymeRates
            r_uni = EnzymeRates.EnzymeReactionLegacy(
                ((:S, ((:C, 1),)),),
                ((:P, ((:C, 1),)),),
            )
            EnzymeRates.init_mechanisms(r_uni)   # warmup; not timed
            r_ter = EnzymeRates.EnzymeReactionLegacy(
                ((:A, ((:C, 1),)), (:B, ((:N, 1),)), (:C, ((:O, 1),))),
                ((:P, ((:C, 1),)), (:Q, ((:N, 1),)), (:R, ((:O, 1),))),
            )
            t = @elapsed EnzymeRates.init_mechanisms(r_ter)
            println("ELAPSED:", t)
            """

        t_cold = _measure_elapsed_subprocess(cold_script)
        t_warm = _measure_elapsed_subprocess(warm_script)

        @info "warmup-reuse: t_cold_ter=$(t_cold)s, t_warm_ter=$(t_warm)s, " *
              "ratio=$(round(t_warm/t_cold; digits=3))"

        # Ratio gate: t_warm should be a meaningful fraction less than t_cold
        # (uni-uni warmup paid for shared specializations).
        @test t_warm < t_cold * WARMUP_TER_RATIO_MAX

        # Absolute gate: post-warmup ter-ter wall-clock has a fixed ceiling.
        @test t_warm < WARMUP_T_WARM_TER_MAX_S
    end
end
