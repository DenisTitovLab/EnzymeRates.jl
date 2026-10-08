# ABOUTME: Compile-time regression gates: init_mechanisms trace-compile, rate_equation
# ABOUTME: body-build wall-clock, bi-bi→uni-uni compile reuse and dispatch identity.

using Test
using EnzymeRates

# Budgets calibrated against the current branch tip with 2× headroom.
# init_mechanisms trace-compile is dominated by Step / Species / Mechanism
# struct + @generated accessor specializations (EnzymeReaction is
# non-parametric, so there is no per-arity reaction-type specialization).
# The init_mechanisms baseline is 97 on Julia 1.12 and 101 on Julia 1.10.
const INIT_TRACE_BUDGET                  = 200   # baseline 2026-10-05: 97-101; budget ≈ 2×
const RATE_EQUATION_WALLCLOCK_BUDGET_S   = 6.0   # CI-runner baseline ~2.8s (local ~1.03s); budget = 2× CI

# Anchored to the EnzymeRates module prefix only. Counts every method
# specialization Julia compiles that touches our module — our functions,
# our types, Base methods specialized on our types (e.g.
# Base.hash(::EnzymeRates.Step, ...)), Core.kwcall plumbing, show/print.
# Filters out unrelated package precompiles (Optimization.jl, Tables.jl,
# etc.) which would create dependency-noise false positives.
#
# Intentionally NOT a suffix enumeration — that approach silently
# misses renamed or newly-introduced internal helpers. The module-prefix
# anchor catches everything in our namespace automatically.
const RELEVANT_PRECOMPILE_PATTERN = r"EnzymeRates\."

# Runs `runner_script` in a fresh Julia subprocess under --trace-compile. Returns
# (n, stdout): the number of EnzymeRates-prefixed compilations (-1 on failure) and
# whatever the script printed.
function _count_relevant_precompiles(runner_script::String)
    trace_file = tempname()
    julia_exe = Base.julia_cmd().exec[1]
    cmd = Cmd([julia_exe, "--trace-compile=$(trace_file)",
               "--project=.", "-e", runner_script])
    out_buf = IOBuffer()
    try
        run(pipeline(cmd; stdout=out_buf); wait=true)
    catch e
        @warn "Subprocess failed: $e"
        return -1, ""
    end
    isfile(trace_file) || (@info "trace-compile file missing"; return -1, "")
    n = try
        length(filter(line -> occursin(RELEVANT_PRECOMPILE_PATTERN, line) &&
                              !isempty(strip(line)),
                      collect(eachline(trace_file))))
    finally
        rm(trace_file, force=true)
    end
    n, String(take!(out_buf))
end

# Parses the `<label>:<float>` lines a subprocess printed. Returns a Vector{Float64}
# parallel to `labels` (NaN for any label not found).
function _testhelper_parse_labeled(out::String, labels::Vector{String})
    map(labels) do label
        m = match(Regex("$(label):([0-9.eE+-]+)"), out)
        m === nothing ? NaN : parse(Float64, m.captures[1])
    end
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
    _testhelper_parse_labeled(String(take!(out_buf)), labels)
end

@testset "compile-budget" begin
    # Two fresh subprocesses supply every measurement below. Each must be a fresh
    # process: the test process has already loaded the shared fixtures, and whenever
    # another test file runs before this one (a focused or reordered run)
    # init_mechanisms and the same EnzymeMechanism{...} body are already compiled, so
    # an in-process measurement could read near zero and no gate would trip on a
    # regression.
    #   - trace: init_mechanisms on a bi-bi reaction under --trace-compile, then
    #     uni-uni in the same process (the warm half of the compile-reuse gate).
    #     Reactions are built via the direct EnzymeReaction constructor so the trace
    #     measures only the enumeration pipeline, independent of the DSL parser's
    #     macro-expansion cost. bi-bi is non-trivial, so the gate is representative;
    #     uni-uni is too small to catch regressions.
    #   - cold: uni-uni init_mechanisms alone (the cold half of the compile-reuse
    #     gate), then the first rate_equation call. init_mechanisms compiles none of
    #     the @generated rate_equation body, so the first call costs what it costs in
    #     a process that never ran init_mechanisms.
    trace_script = """
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
        r_uni = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P),   [:C => 1])],
            EnzymeRates.RegulatorMults[],
            Int[1],
        )
        GC.gc()
        t_uni = @elapsed EnzymeRates.init_mechanisms(r_uni)
        println("UNI_WARM:", t_uni)
        """
    cold_script = """
        using EnzymeRates
        r_uni = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        t = @elapsed EnzymeRates.init_mechanisms(r_uni)
        println("UNI_COLD:", t)
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
    n, trace_out = _count_relevant_precompiles(trace_script)
    t_uni_warm = _testhelper_parse_labeled(trace_out, ["UNI_WARM"])[1]
    t_uni_cold, t_first = _measure_labeled_subprocess(cold_script, ["UNI_COLD", "ELAPSED"])

    # Trace-compile: the bi-bi init_mechanisms (the uni-uni that follows compiles
    # nothing new while reuse holds, so it adds no trace lines).
    @testset "trace-compile: init_mechanisms (bi-bi)" begin
        @info "init_mechanisms trace-compile: $n (budget: $INIT_TRACE_BUDGET)"
        @test 0 <= n <= INIT_TRACE_BUDGET
    end

    # Wall-clock: rate_equation body-build (first call pays @generated cost). The
    # per-call runtime gate is separately enforced by test_rate_eq_derivation.jl's
    # test_rate_equation_performance (0 allocs, <120ns per call) for every
    # mechanism in MECHANISM_TEST_SPECS.
    @testset "wall-clock: rate_equation body-build (first call)" begin
        @info "rate_equation first-call wall-clock: $(t_first)s " *
              "(budget: $RATE_EQUATION_WALLCLOCK_BUDGET_S s)"
        @test isfinite(t_first)
        @test t_first < RATE_EQUATION_WALLCLOCK_BUDGET_S
    end

    # Compile-reuse: bi-bi init_mechanisms compiles a superset of uni-uni's
    # machinery, so running uni-uni AFTER bi-bi in the same process is essentially
    # free:
    #   - cold:  uni-uni alone           → t_uni_cold ≈ 1-2 s
    #   - warm:  bi-bi, then uni-uni      → t_uni_warm ≈ 0.2-1.6 ms
    # The warm/cold ratio (≈ 1e-4 to 2e-3 on CI runners; macOS is the noisy high
    # end) is robust to machine speed, unlike an absolute wall-clock ceiling on
    # the cold time. The in-process ter-ter ceiling lives in
    # test_mechanism_enumeration.jl.
    @testset "compile reuse: bi-bi warms all of uni-uni" begin
        @info "compile reuse: " *
              "uni_cold=$(round(t_uni_cold; digits=2))s  " *
              "uni_warm=$(round(t_uni_warm * 1e6; digits=1))µs  " *
              "warm/cold=$(round(t_uni_warm / t_uni_cold; sigdigits=2))"
        # Warm uni-uni must be near-instant relative to cold: bi-bi already
        # compiled the superset. A lost reuse recompiles a sizeable fraction of
        # cold (warm/cold ≳ 0.1); the < 1e-2 gate sits ~6× above the noisiest
        # observed CI ratio and ~10× below a real failure.
        @test isfinite(t_uni_cold) && t_uni_cold > 0
        @test isfinite(t_uni_warm)
        @test t_uni_warm / t_uni_cold < 0.01
    end

    # Dispatch identity: EnzymeReaction is non-parametric, so uni-uni and
    # ter-ter are the same concrete type and `init_mechanisms(::EnzymeReaction)`
    # resolves to the same method instance — no per-arity specialization.
    @testset "dispatch identity: init_mechanisms uni-uni and ter-ter share method" begin
        r_uni = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        r_ter = @enzyme_reaction begin
            substrates: A[C], B[N], C[O]
            products:   P[C], Q[N], R[O]
        end
        @test typeof(r_uni) === typeof(r_ter)
        @test which(EnzymeRates.init_mechanisms, (typeof(r_uni),)) ===
              which(EnzymeRates.init_mechanisms, (typeof(r_ter),))
    end
end
