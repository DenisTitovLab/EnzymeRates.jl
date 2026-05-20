# ABOUTME: Compile-time regression gates for the EnzymeRates pipeline.
# ABOUTME: One trace-compile (init_mechanisms) + one wall-clock (rate_equation body-build).

using Test
using EnzymeRates

# Budgets calibrated against current main with 2× headroom (Step 4).
# Per spec §3: this budget is FIXED for the refactor's duration.
# Future stages REDESIGN on trip rather than recalibrate.
const INIT_TRACE_BUDGET                  = 58    # baseline 2026-05-20: 29; budget = 2× (rounded up)
const RATE_EQUATION_WALLCLOCK_BUDGET_S   = 2.1   # baseline 2026-05-20: 1.03s; budget = 2×

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
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
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

        # Two gates on warmup reuse: (a) the ratio of warm:cold must
        # stay below 0.7 (today's main is ~0.55; ~30% sharing demanded),
        # (b) an absolute ceiling catches a regression that explodes
        # arity-dependent specialization regardless of cold cost. Current
        # baseline: ter_warm_incr = 16; cap at 25 for 50% headroom.
        # See Stage 0 PR for empirical calibration rationale.
        @test ter_warm_incr <= 7 * ter_cold_incr ÷ 10
        @test ter_warm_incr <= 25
    end
end
