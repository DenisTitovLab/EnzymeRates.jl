# ABOUTME: Compile-time regression gates for the EnzymeRates pipeline:
# ABOUTME: init_mechanisms trace-compile, rate_equation body-build wall-clock, ter-ter→uni-uni compile reuse, dispatch identity.

using Test
using EnzymeRates

# Budgets calibrated against the current branch tip with 2× headroom.
# init_mechanisms trace-compile is dominated by Step / Species / Mechanism
# struct + @generated accessor specializations (EnzymeReaction is
# non-parametric, so there is no per-arity reaction-type specialization).
const INIT_TRACE_BUDGET                  = 100   # baseline 2026-05-27: 47; budget ≈ 2×
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

# Runs `script` in a fresh Julia subprocess; the script is expected to print
# `<label>:<float>` for each label in `labels`. Returns a Vector{Float64}
# parallel to `labels` (NaN for any label not found or on subprocess failure).
function _measure_labeled_subprocess(script::String, labels::Vector{String})
    julia_exe = Base.julia_cmd().exec[1]
    out_buf = IOBuffer()
    try
        run(pipeline(Cmd([julia_exe, "--project=.", "-e", script]);
                     stdout=out_buf, stderr=devnull); wait=true)
    catch e
        @warn "labeled subprocess failed: $e"
        return fill(NaN, length(labels))
    end
    out = String(take!(out_buf))
    map(labels) do label
        m = match(Regex("$(label):([0-9.eE+-]+)"), out)
        m === nothing ? NaN : parse(Float64, m.captures[1])
    end
end

@testset "compile-budget" begin
    # Trace-compile: init_mechanisms on a bi-bi reaction (non-trivial so
    # the gate is representative; uni-uni is too small to catch regressions).
    # Constructs EnzymeReaction via the direct constructor so the trace
    # measures only the enumeration pipeline, independent of the DSL
    # parser's macro-expansion cost.
    @testset "trace-compile: init_mechanisms (bi-bi)" begin
        script = """
            using EnzymeRates
            r = EnzymeRates.EnzymeReaction(
                [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:A), [:C => 1]),
                 EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:B), [:N => 1]),
                 EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P),   [:C => 1]),
                 EnzymeRates.ReactantAtoms(EnzymeRates.Product(:Q),   [:N => 1])],
                EnzymeRates.RegulatorMults[],
                Int[1],
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

    # Compile-reuse: ter-ter init_mechanisms compiles a superset of uni-uni's
    # machinery, so running uni-uni AFTER ter-ter in the same process is
    # essentially free. Measured in two fresh subprocesses (reactions built
    # via @enzyme_reaction, as the original main gate did):
    #   - cold:  uni-uni alone           → t_uni_cold ≈ 3 s
    #   - warm:  ter-ter, then uni-uni    → t_ter ≈ 15 s, t_uni_warm ≈ 0.1 ms
    # The warm/cold ratio (≈ 5e-5) is robust to machine speed (warm is ~0
    # regardless), unlike an absolute wall-clock ceiling on the cold time.
    @testset "compile reuse: ter-ter warms all of uni-uni" begin
        cold_script = """
            using EnzymeRates
            r_uni = @enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
            end
            t = @elapsed EnzymeRates.init_mechanisms(r_uni)
            println("UNI_COLD:", t)
            """
        warm_script = """
            using EnzymeRates
            r_ter = @enzyme_reaction begin
                substrates: A[C], B[N], C[O]
                products:   P[C], Q[N], R[O]
            end
            t_ter = @elapsed EnzymeRates.init_mechanisms(r_ter)
            r_uni = @enzyme_reaction begin
                substrates: S[C]
                products:   P[C]
            end
            t_uni = @elapsed EnzymeRates.init_mechanisms(r_uni)
            println("TER_COLD:", t_ter)
            println("UNI_WARM:", t_uni)
            """
        t_uni_cold = _measure_labeled_subprocess(cold_script, ["UNI_COLD"])[1]
        t_ter, t_uni_warm =
            _measure_labeled_subprocess(warm_script, ["TER_COLD", "UNI_WARM"])
        @info "compile reuse: ter_cold=$(round(t_ter; digits=2))s  " *
              "uni_cold=$(round(t_uni_cold; digits=2))s  " *
              "uni_warm=$(round(t_uni_warm * 1e6; digits=1))µs  " *
              "warm/cold=$(round(t_uni_warm / t_uni_cold; sigdigits=2))"
        # ter-ter cold-compile ceiling (~2× the 2026-05-27 baseline of ~15 s).
        @test isfinite(t_ter)
        @test t_ter < 30.0
        # Warm uni-uni must be near-instant relative to cold: ter-ter already
        # compiled the superset. Observed warm/cold ≈ 5e-5; the < 1e-3 gate
        # keeps a ~20× margin and is insensitive to machine speed.
        @test isfinite(t_uni_cold) && t_uni_cold > 0
        @test isfinite(t_uni_warm)
        @test t_uni_warm / t_uni_cold < 0.001
    end

    # Dispatch identity: EnzymeReaction is non-parametric, so uni-uni and
    # ter-ter are the same concrete type and `init_mechanisms(::EnzymeReaction)`
    # resolves to the same method instance — no per-arity specialization.
    @testset "dispatch identity: init_mechanisms uni-uni and ter-ter share method" begin
        r_uni = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P),   [:C => 1])],
            EnzymeRates.RegulatorMults[],
            Int[1],
        )
        r_ter = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:A), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:B), [:N => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:C), [:O => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P),   [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:Q),   [:N => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:R),   [:O => 1])],
            EnzymeRates.RegulatorMults[],
            Int[1],
        )
        @test typeof(r_uni) === typeof(r_ter)
        @test which(EnzymeRates.init_mechanisms, (typeof(r_uni),)) ===
              which(EnzymeRates.init_mechanisms, (typeof(r_ter),))
    end
end
