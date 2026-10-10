# ABOUTME: Compile-time regression gates: init_mechanisms trace-compile, rate_equation
# ABOUTME: body-build wall-clock, bi-bi→uni-uni reuse, dispatch identity, fit per mechanism.

using Test
using EnzymeRates

# Budgets calibrated against the current branch tip with 2× headroom.
# init_mechanisms trace-compile is dominated by Step / Species / Mechanism
# struct + @generated accessor specializations (EnzymeReaction is
# non-parametric, so there is no per-arity reaction-type specialization).
# The init_mechanisms baseline is 57 on Julia 1.13 (CI runners), 82 on Julia 1.12 and 87
# on Julia 1.10.
const INIT_TRACE_BUDGET                  = 200   # baseline 2026-10-09: 57-87; budget ≈ 2×
const RATE_EQUATION_WALLCLOCK_BUDGET_S   = 6.0   # CI 1.1-1.9 s, local 0.4 s
# Per fresh mechanism, over PASS1, the fit and one LOOCV fold, with baselines measured on
# Julia 1.12 (aarch64). A fit path that compiles the solver stack for each mechanism costs
# about 750 method instances and 650 KB. The budgets sit about 3× below that and leave
# room for the other Julia versions and architectures CI runs.
const FIT_INSTANCE_BUDGET                = 250   # baseline 2026-10-09: 64
const FIT_NATIVE_BUDGET_KB               = 300   # baseline 2026-10-09: 54
# Encoding the Sigs of the 239 bi-bi seeds, once per process: Base's collect widening
# compiles once per new mix of step shapes, so a search pays this early and once.
const SIG_INSTANCE_BUDGET                = 4000  # baseline 2026-10-09: 1888
const SIG_NATIVE_BUDGET_KB               = 1000  # baseline 2026-10-09: 508

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
function _testhelper_count_relevant_precompiles(runner_script::String)
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

# Runs `script` in a fresh Julia subprocess on the test environment, so the script can
# load test dependencies; the script is expected to print `<label>:<float>` for each
# label in `labels`. Returns a Vector{Float64} parallel to `labels` (NaN for any label
# not found or on subprocess failure).
function _testhelper_measure_labeled_subprocess(script::String, labels::Vector{String})
    julia_exe = Base.julia_cmd().exec[1]
    out_buf = IOBuffer()
    try
        run(pipeline(Cmd([julia_exe, "--project=$(Base.active_project())", "-e", script]);
                     stdout=out_buf, stderr=devnull); wait=true)
    catch e
        @warn "labeled subprocess failed: $e"
        return fill(NaN, length(labels))
    end
    _testhelper_parse_labeled(String(take!(out_buf)), labels)
end

@testset "compile-budget" begin
    # Three fresh subprocesses supply every measurement below. Each must be a fresh
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
    #   - fit: the identify path on bi-bi mechanisms, step for step as _process_batch
    #     and _cv_fold_loss take it: PASS1 (complexity, compile, fitted names, rate
    #     equation string, dedup key, mechanism_type string), PASS2 (FittingProblem and
    #     a CMA-ES fit) and one LOOCV fold. Encoding every seed's Sig and a full pass
    #     on ms[1] compile everything that does not depend on the mechanism, so each
    #     fresh mechanism after them compiles only what its own type needs. ms[1] and
    #     the fresh mechanisms differ in Sig shape and in fitted names. The census
    #     counts the method instances of every module, since a fit compiled per
    #     mechanism compiles Optimization, SciMLBase, CMAEvolutionStrategy and Base
    #     code; it also counts the show instances PASS1 adds, since Base's show prints
    #     a type by compiling show again for each new Sig shape.
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
    fit_script = """
        using EnzymeRates, OptimizationCMAEvolutionStrategy, Random
        const ER = EnzymeRates
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end
        const PROB = IdentifyRateEquationProblem(rxn,
            (group = repeat(["G1", "G2"], inner = 4),
             Rate = [0.3, 0.5, 0.6, 0.8, 0.2, 0.4, 0.7, 0.9],
             A = repeat([0.5, 1.0, 2.0, 4.0], 2), B = repeat([1.0, 2.0], 4),
             P = fill(0.1, 8), Q = fill(0.2, 8)); Keq = 2.0)
        const OPT = CMAEvolutionStrategyOpt()
        const FIT_KW = (; n_restarts = 1, maxtime = 10.0, maxiters = 20, abstol = nothing,
                        reltol = nothing, callback = nothing, solver_kwargs = (;))
        _testhelper_instances(meths) =
            sum(m -> count(_ -> true, Base.specializations(m)), meths)
        function _testhelper_all_instances()
            meths = Method[]
            Base.visit(m -> push!(meths, m), Core.methodtable)
            _testhelper_instances(meths)
        end
        # The method instances and native code bytes that PASS1, PASS2 and one LOOCV
        # fold of m add, their compile time, and the show instances PASS1 adds.
        function _testhelper_measure(m)
            n0, b0 = _testhelper_all_instances(), Base.jit_total_bytes()
            t0 = Base.cumulative_compile_time_ns()[1]
            s0 = _testhelper_instances(methods(show))
            ER._eq_complexity(m)
            em = ER.compile_mechanism(m)
            ER.fitted_params(em)
            eq_text = rate_equation_string(em)
            ER._rate_eq_dedup_key(eq_text)
            ER._mechanism_type_string(em)
            s1 = _testhelper_instances(methods(show))
            fp = FittingProblem(ER.compile_mechanism(m), PROB.data;
                                Keq = PROB.Keq, scale_k_to_kcat = PROB.scale_k_to_kcat)
            fit_rate_equation(fp, OPT; FIT_KW...)
            ER._cv_fold_loss(ER.compile_mechanism(m), PROB, "G1";
                             optimizer = OPT, FIT_KW...)
            (_testhelper_all_instances() - n0, Base.jit_total_bytes() - b0,
             (Base.cumulative_compile_time_ns()[1] - t0) / 1e9, s1 - s0)
        end
        ms = ER.init_mechanisms(rxn)
        # Encoding a Sig compiles Base's collect widening once per new mix of step
        # shapes, a cost a search pays early and once; encoding every seed first keeps
        # it out of the per-mechanism counts and measures its total for this process.
        n_sig, b_sig = _testhelper_all_instances(), Base.jit_total_bytes()
        foreach(ER._sig_of, ms)
        println("SIG_INSTANCES:", _testhelper_all_instances() - n_sig,
                " SIG_NATIVE:", Base.jit_total_bytes() - b_sig)
        fresh = [2, 6, 11, 18]
        Random.seed!(1)
        Base.cumulative_compile_timing(true)
        _testhelper_measure(ms[1])
        for (k, (n, b, t, s)) in enumerate(map(i -> _testhelper_measure(ms[i]), fresh))
            println("FIT_INSTANCES_", k, ":", n, " FIT_NATIVE_", k, ":", b,
                    " FIT_COMPILE_", k, ":", t, " PASS1_SHOW_", k, ":", s)
        end
        ems = ER.compile_mechanism.(ms[[1; fresh]])
        println("FIT_SHAPES:", length(unique(em -> typeof(typeof(em).parameters[1]), ems)))
        println("FIT_NAMES:", length(unique(ER.fitted_params, ems)))
        """
    n, trace_out = _testhelper_count_relevant_precompiles(trace_script)
    t_uni_warm = _testhelper_parse_labeled(trace_out, ["UNI_WARM"])[1]
    t_uni_cold, t_first = _testhelper_measure_labeled_subprocess(
        cold_script, ["UNI_COLD", "ELAPSED"])

    # Trace-compile: the bi-bi init_mechanisms (the uni-uni that follows compiles
    # nothing new while reuse holds, so it adds no trace lines).
    @testset "trace-compile: init_mechanisms (bi-bi)" begin
        @info "init_mechanisms trace-compile: $n (budget: $INIT_TRACE_BUDGET)"
        @test 0 <= n <= INIT_TRACE_BUDGET
    end

    # Wall-clock: rate_equation body-build (first call pays @generated cost). The
    # per-call runtime gate is separately enforced by test_rate_eq_derivation.jl's
    # _testhelper_test_rate_equation_performance (0 allocs, <120ns per call) for every
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
    #   - cold:  uni-uni alone           → t_uni_cold ≈ 1.6-2.3 s on CI runners
    #   - warm:  bi-bi, then uni-uni      → t_uni_warm ≈ 0.2-0.5 ms
    # The warm/cold ratio (≈ 1e-4 to 2.5e-4 on CI runners; macOS is the noisy high
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
        # cold (warm/cold ≳ 0.1); the < 1e-2 gate sits ~40× above the noisiest
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

    # The fit census walks Julia's global method table, which Julia has from 1.12 on.
    if isdefined(Core, :methodtable)
        fit_labels = [string(metric, "_", k)
                      for metric in ("FIT_INSTANCES", "FIT_NATIVE", "FIT_COMPILE",
                                     "PASS1_SHOW")
                      for k in 1:4]
        fit_vals = _testhelper_measure_labeled_subprocess(fit_script,
            [fit_labels; "FIT_SHAPES"; "FIT_NAMES"; "SIG_INSTANCES"; "SIG_NATIVE"])
        # One column per metric, one row per fresh mechanism.
        instances, native, compile_s, pass1_show = eachcol(reshape(fit_vals[1:16], 4, 4))
        n_shapes, n_name_tuples, sig_instances, sig_native = fit_vals[17:20]

        # Sig encoding: what encoding the 239 bi-bi seeds compiles, once per process. It
        # grows only with new mixes of step shapes, so a total far above the baseline
        # means encoding started compiling for each mechanism.
        @testset "Sig encoding: total compile for the bi-bi seeds" begin
            @info "Sig encoding of the bi-bi seeds: $sig_instances method instances " *
                  "(budget: $SIG_INSTANCE_BUDGET), $(round(sig_native / 1024; digits=1)) " *
                  "KB native code (budget: $SIG_NATIVE_BUDGET_KB KB)"
            @test sig_instances <= SIG_INSTANCE_BUDGET
            @test sig_native <= SIG_NATIVE_BUDGET_KB * 1024
        end

        # Fit compile: what PASS1, PASS2 and one LOOCV fold compile for each fresh
        # mechanism, counted in every module.
        @testset "fit compile: PASS1 + fit + LOOCV fold per fresh mechanism" begin
            @info "fit compile per fresh mechanism: method instances $instances " *
                  "(budget: $FIT_INSTANCE_BUDGET), native code " *
                  "$(round.(native ./ 1024; digits=1)) KB " *
                  "(budget: $FIT_NATIVE_BUDGET_KB KB), " *
                  "compile $(round.(compile_s; digits=3)) s"
            # ms[1] and the four fresh mechanisms have five Sig shapes and five
            # fitted-name tuples, so no fresh mechanism reuses another's per-type code.
            @test n_shapes == 5
            @test n_name_tuples == 5
            # Every fresh type compiles its own rate equation, so a census that counts
            # nothing is broken.
            @test minimum(instances) > 0 && minimum(native) > 0
            @test maximum(instances) <= FIT_INSTANCE_BUDGET
            @test maximum(native) <= FIT_NATIVE_BUDGET_KB * 1024
        end

        # mechanism_type: PASS1 prints the type of a mechanism with a fresh Sig shape
        # without compiling show for it.
        @testset "fit compile: PASS1 compiles no show for a fresh Sig shape" begin
            @info "PASS1 show instances per fresh mechanism: $pass1_show (budget: 0)"
            @test maximum(pass1_show) == 0
        end
    else
        @test_skip "fit compile: the census needs the global method table of Julia 1.12"
    end
end
