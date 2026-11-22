# ABOUTME: Tests of src/Fitting/fit_rate_equation.jl: FittingProblem, loss!, the scales
# ABOUTME: of a fit, its search around the data, and fit_rate_equation.

# ── FittingProblem, loss! and fit_rate_equation ───────────────────────────────────────

using Tables
using Optimization
using OptimizationCMAEvolutionStrategy: CMAEvolutionStrategyOpt
using Optimization.SciMLBase: build_solution, ReturnCode, DefaultOptimizationCache

# An optimizer stub whose every solve returns `u` for each log-space coordinate, so a
# fit's raw constants are known exactly and the scale applied to them can be checked.
# With `infinite = true` every solve reports an infinite objective, as a fit with no
# finite loss does.
struct _testhelper_ConstantStubOpt
    u::Float64
    infinite::Bool
end
_testhelper_ConstantStubOpt(u; infinite = false) = _testhelper_ConstantStubOpt(u, infinite)
Optimization.allowsbounds(::_testhelper_ConstantStubOpt) = true
function Optimization.SciMLBase.__solve(prob::Optimization.OptimizationProblem,
        opt::_testhelper_ConstantStubOpt; kwargs...)
    u = fill(opt.u, length(prob.u0))
    build_solution(DefaultOptimizationCache(prob.f, prob.p), opt, u,
                   opt.infinite ? Inf : prob.f(u, prob.p); retcode = ReturnCode.Success)
end

# An optimizer that hands the problem to `inner` and counts the objective's evaluations.
struct _testhelper_CountingOpt{O}
    inner::O
    evals::Base.RefValue{Int}
end
Optimization.allowsbounds(::_testhelper_CountingOpt) = true
function Optimization.SciMLBase.__solve(prob::Optimization.OptimizationProblem,
        opt::_testhelper_CountingOpt; kwargs...)
    objective = Optimization.OptimizationFunction((x, p) -> begin
        opt.evals[] += 1
        prob.f(x, p)
    end)
    counted = Optimization.SciMLBase.remake(prob; f = objective)
    Optimization.solve(counted, opt.inner; kwargs...)
end

@testset "Fitting" begin

    # ── Helper: build a Uni-Uni mechanism ─────────────────────────────────────
    # Fitted: k_E_S_to_ES, k_ES_to_E_S, k_ES_to_E_P. The true points below put
    # k_ES_to_E_P = 5 at Keq = 2, so the Haldane-derived product rebinding rate is
    # k_E_P_to_ES = k_E_S_to_ES·k_ES_to_E_P/(Keq·k_ES_to_E_S) = 1.
    uni_uni = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end

    # True point and 5-point concentration grid shared by the uni-uni testsets below.
    Keq_val = 2.0
    true_params = (k_E_S_to_ES = 10.0, k_ES_to_E_S = 25.0, k_ES_to_E_P = 5.0)
    concs5 = [
        (S = 1.0, P = 0.1),
        (S = 2.0, P = 0.1),
        (S = 5.0, P = 0.1),
        (S = 1.0, P = 0.5),
        (S = 2.0, P = 0.5),
    ]

    # ── Helper: an allosteric uni-uni with an inactive-state regulator R ──────
    allo = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    allo_fitted = parameters(allo)
    allo_params = NamedTuple{allo_fitted}(ntuple(i -> 1.0 + 0.1 * i, length(allo_fitted)))
    allo_concs = [(S = Float64(i), P = 0.1, R = 0.5) for i in 1:20]

    # ── Synthetic data generator ──────────────────────────────────────────────
    # The rates of `mechanism` at `true_params` and Keq = `Keq`, Keq_val by default.
    function _testhelper_make_synthetic_data(
            mechanism, true_params, concs_list;
            groups=fill("G1", length(concs_list)), Keq = Keq_val,
    )
        rates = Float64[]
        rate = compile_rate_equation(mechanism)
        for (i, concs) in enumerate(concs_list)
            r = rate(concs, true_params, Keq)
            push!(rates, r)
        end
        met_names = metabolites(mechanism)
        cols = Dict{Symbol, Vector}()
        cols[:group] = groups
        cols[:Rate] = rates
        for mn in met_names
            cols[mn] = [concs[mn] for concs in concs_list]
        end
        return (; (k => cols[k] for k in (:group, :Rate, met_names...))...)
    end

    # ── BlackBoxOptim fit: the only warning is its convergence ────────────────
    # BlackBoxOptim can converge on these fits before maxtime: an absolute loss has a
    # single optimum, and a centered population can collapse along its scale direction.
    # It then stops because the search has converged, a reason Optimization does not
    # recognize, so Optimization warns. Runs `fit`, checks that every warning it logs
    # is that one, and returns its result.
    function _testhelper_fit_capturing_convergence(fit)
        logs, result = Test.collect_test_logs(fit)
        @test all(l -> l.level < Base.CoreLogging.Warn ||
                       occursin("probably search has converged", string(l.message)), logs)
        result
    end

    # ── Mechanism-level accessors ─────────────────────────────────────────────
    @testset "Mechanism-level accessors" begin
        @test Set(parameters(uni_uni)) == Set(keys(true_params))
        @test metabolites(uni_uni) == [:S, :P]
    end

    # ── FittingProblem construction ───────────────────────────────────────────
    @testset "Construction" begin
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        fp = FittingProblem(uni_uni, data, Keq_val;
                            reference=_testhelper_row_reference(data))

        @test length(fp.log_abs_rates) == 5
        @test length(fp.group_point_indexes) == 1
        @test length(fp.group_point_indexes[1]) == 5
    end

    # ── Multi-group centering invariance ─────────────────────────────────────
    @testset "Multi-group centering invariance" begin
        # Two groups, each independently scaled
        data1 = _testhelper_make_synthetic_data(uni_uni, true_params, concs5;
            groups=["G1","G1","G1","G2","G2"])
        fp1 = FittingProblem(uni_uni, data1, Keq_val;
                             reference=_testhelper_row_reference(data1))

        # Scale group1 by 5x and group2 by 100x
        rates2 = copy(data1.Rate)
        rates2[1:3] .*= 5.0
        rates2[4:5] .*= 100.0
        data2 = merge(data1, (Rate = rates2,))
        fp2 = FittingProblem(uni_uni, data2, Keq_val;
                             reference=_testhelper_row_reference(data2))

        np = length(parameters(uni_uni))
        # A draw can fit these exact data closely: at a loss near 1e-8 the two losses differ
        # by round-off alone, 1e-11 relative but never above 1e-15 absolute over 5 million
        # draws, so a relative tolerance alone fails about one run in 900.
        @test all(1:10) do _
            x = randn(np) .* 2.0
            isapprox(EnzymeRates.Fitting.loss!(x, fp1), EnzymeRates.Fitting.loss!(x, fp2);
                     rtol=1e-12, atol=1e-14)
        end
    end

    # ── Absolute mode: uncentered loss (absolute=true) ─────────────────────────
    @testset "Absolute mode uncentered loss" begin
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        pn = parameters(uni_uni)
        x_true = [log(true_params[p]) for p in pn]

        fp_rel = FittingProblem(uni_uni, data, Keq_val;  # relative (default)
                                reference=_testhelper_row_reference(data))
        fp_abs = FittingProblem(uni_uni, data, Keq_val; absolute=true)

        # At true params both modes are ~0 (predictions match data exactly).
        @test EnzymeRates.Fitting.loss!(x_true, fp_rel) ≈ 0.0 atol=1e-20
        @test EnzymeRates.Fitting.loss!(x_true, fp_abs) ≈ 0.0 atol=1e-20

        # Scale every rate by 3 (a pure per-group offset). Relative loss is
        # invariant (centering removes it); absolute loss sees it: every
        # residual becomes log(3), so absolute loss = log(3)^2.
        data3 = merge(data, (Rate = data.Rate .* 3.0,))
        fp_rel3 = FittingProblem(uni_uni, data3, Keq_val;
                                 reference=_testhelper_row_reference(data3))
        fp_abs3 = FittingProblem(uni_uni, data3, Keq_val; absolute=true)
        @test EnzymeRates.Fitting.loss!(x_true, fp_rel3) ≈ 0.0 atol=1e-20
        @test EnzymeRates.Fitting.loss!(x_true, fp_abs3) ≈ log(3.0)^2 rtol=1e-8
    end

    # ── absolute and reference keywords ───────────────────────────────────────
    @testset "absolute and reference keywords" begin
        ok_data = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1])
        ref = _testhelper_row_reference(ok_data)
        fp = FittingProblem(uni_uni, ok_data, 1.0; reference=ref)
        @test fp.rate_data.absolute === false
        @test fp.rate_data.reference == [(S = 1.0, P = 0.1) => 1.0]
        fp_abs = FittingProblem(uni_uni, ok_data, 1.0; absolute=true)
        @test fp_abs.rate_data.absolute === true && isempty(fp_abs.rate_data.reference)
        # An integer Keq converts to the Float64 field.
        @test FittingProblem(uni_uni, ok_data, 2; reference=ref).rate_data.Keq === 2.0
        # The concrete-mechanism constructor passes the keywords on.
        fp_conc = FittingProblem(uni_uni, ok_data, 1.0;
                                 reference=ref)
        @test fp_conc.rate_data.absolute === false && fp_conc.rate_data.reference == ref
        @test FittingProblem(uni_uni, ok_data, 1.0;
                             absolute=true).rate_data.absolute === true
        # There is no kcat keyword: a `reference` sets the scale.
        @test_throws MethodError FittingProblem(uni_uni, ok_data, 1.0; reference=ref,
                                                scale_k_to_kcat=1.0)
    end

    # ── Sign-mismatch penalty ─────────────────────────────────────────────────
    # Regression test for all-mismatch groups: when every prediction in a
    # group is a sign mismatch, centering must not zero every deviation.
    # The loss must be nonzero to distinguish a bad mechanism from a perfect one.
    # Three sub-cases exercise each kind of mismatch that sets buf[i] = 10.0:
    #   (i)  pred == 0.0            (S=0, P=0)
    #   (ii) pred < 0, Rate > 0     (S=0, P>0: only the reverse term survives, so pred
    #        is always negative)
    #   (iii) pred > 0, Rate < 0    (S>0, P=0: only the forward term survives, so pred
    #        is always positive)
    # In all cases every point in the group is a mismatch so the expected
    # loss is: (0 from centering + 100.0 × n_mismatch) / n_data = 100.0
    @testset "All-mismatch group not cancelled by centering" begin
        pn = parameters(uni_uni)
        x = randn(length(pn))

        @testset "zero prediction (S=0, P=0)" begin
            data = (
                group = fill("G1", 5),
                Rate  = [1.0, 2.0, 3.0, 4.0, 5.0],
                S     = zeros(5),
                P     = zeros(5),
            )
            fp = FittingProblem(uni_uni, data, Keq_val;
                                reference=_testhelper_row_reference(data))
            l = EnzymeRates.Fitting.loss!(x, fp)
            @test l > 0.0
            @test l ≈ 100.0
        end

        @testset "sign mismatch: pred<0, Rate>0 (S=0, P>0)" begin
            # With S=0 the forward numerator term vanishes; only the reverse
            # term (∝ -P) remains, so pred < 0 for any positive parameters.
            data = (
                group = fill("G1", 5),
                Rate  = [1.0, 2.0, 3.0, 4.0, 5.0],
                S     = zeros(5),
                P     = [0.1, 0.2, 0.3, 0.4, 0.5],
            )
            fp = FittingProblem(uni_uni, data, Keq_val;
                                reference=_testhelper_row_reference(data))
            l = EnzymeRates.Fitting.loss!(x, fp)
            @test l > 0.0
            @test l ≈ 100.0
        end

        @testset "sign mismatch: pred>0, Rate<0 (S>0, P=0)" begin
            # With P=0 the reverse numerator term vanishes; only the forward
            # term (∝ S) remains, so pred > 0 for any positive parameters.
            data = (
                group = fill("G1", 5),
                Rate  = [-1.0, -2.0, -3.0, -4.0, -5.0],
                S     = [0.1, 0.2, 0.3, 0.4, 0.5],
                P     = zeros(5),
            )
            fp = FittingProblem(uni_uni, data, Keq_val;
                                reference=_testhelper_row_reference(data))
            l = EnzymeRates.Fitting.loss!(x, fp)
            @test l > 0.0
            @test l ≈ 100.0
        end
    end

    # ── Zero allocations ──────────────────────────────────────────────────────
    @testset "Zero allocations" begin
        concs_list = [(S = Float64(i), P = 0.1) for i in 1:20]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data, Keq_val;
                            reference=_testhelper_row_reference(data))

        x = randn(length(parameters(uni_uni)))
        _testhelper_loss_allocated(x, fp)  # warmup
        allocs = _testhelper_loss_allocated(x, fp)
        @test allocs == 0

        # Absolute mode (uncentered branch) is equally allocation-free.
        fp_abs = FittingProblem(uni_uni, data, Keq_val; absolute=true)
        _testhelper_loss_allocated(x, fp_abs)  # warmup
        allocs_abs = _testhelper_loss_allocated(x, fp_abs)
        @test allocs_abs == 0
    end

    # loss! builds each data point's concentration NamedTuple from metabolites(m), so an
    # allosteric mechanism's metabolite names must be a compile-time constant as well.
    @testset "Zero allocations: allosteric mechanism" begin
        data = _testhelper_make_synthetic_data(allo, allo_params, allo_concs)
        fp = FittingProblem(allo, data, Keq_val;
                            reference=_testhelper_row_reference(data))

        x = randn(length(allo_fitted))
        _testhelper_loss_allocated(x, fp)  # warmup
        allocs = _testhelper_loss_allocated(x, fp)
        @test allocs == 0
    end

    # ── loss! over table column types ─────────────────────────────────────────
    # The gate in test_rate_equation_derivation.jl measures loss! on every law;
    # these tests vary the table. A column of abstract element type makes every row
    # dispatch at run time while the rate law stays fast, and `@inferred loss!` passes
    # since the result is a Float64 either way, so they measure allocations.
    @testset "loss! allocates nothing whatever the table's column types" begin
        two_regulators = @enzyme_mechanism begin
            substrates: A, B
            products: P
            oligomeric_state: 4
            steps: begin
                E + A ⇌ EA :: EqualAI
                EA + B <--> EAB :: OnlyA
                EAB <--> EP :: OnlyA
                EP ⇌ E + P :: EqualAI
            end
            regulatory_site(4): R::OnlyI
            regulatory_site(4): Q::NonequalAI
        end
        no_regulators = @enzyme_mechanism begin
            substrates: A
            products: P
            oligomeric_state: 1
            steps: begin
                E + A ⇌ EA :: EqualAI
                EA <--> EP :: OnlyA
                EP ⇌ E + P :: EqualAI
            end
        end
        plain = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A ⇌ EA
                EA + B <--> EAB
                EAB <--> EP
                EP ⇌ E + P
            end
        end
        retype(table, T) = merge(table, (Rate = convert(Vector{T}, table.Rate),
                                         A = convert(Vector{T}, table.A)))
        mechanisms = (; two_regulators, no_regulators, plain)
        @testset "$(name)" for (name, m) in pairs(mechanisms)
            table = _testhelper_loss_table(metabolites(m))
            n = length(table.Rate)
            shapes = (
                union_missing = retype(table, Union{Missing, Float64}),
                any = retype(table, Any),
                float32 = retype(table, Float32),
                int_concentration = merge(table, (A = fill(2, n),)),
                views = map(c -> view(c, 1:n), table),
                extra_columns = merge(table, (Article = fill("a", n), Other = rand(n))),
                integer_group = merge(table, (group = [i % 4 for i in 1:n],)),
            )
            @testset "$(shape)" for (shape, data) in pairs(shapes)
                fp = FittingProblem(m, data, Keq_val)
                @test _testhelper_loss_allocations(fp) == 0
            end
            # A column's element type changes the cost of a row, never its value.
            x = 0.1 .* randn(length(parameters(m)))
            expected = EnzymeRates.Fitting.loss!(x, FittingProblem(m, table, Keq_val))
            for shape in (:union_missing, :any, :views)
                fp = FittingProblem(m, shapes[shape], Keq_val)
                @test EnzymeRates.Fitting.loss!(x, fp) == expected
            end
        end
    end

    @testset "loss! allocates nothing on a cross-validation fold" begin
        table = _testhelper_loss_table(metabolites(allo))
        fold = EnzymeRates._fold_problem(allo, table, table.group .!= "g0", Keq_val)
        @test _testhelper_loss_allocations(fold) == 0
    end

    # ── Optimizer boundary ────────────────────────────────────────────────────
    # fit_rate_equation hands Optimization.jl one problem type whatever the mechanism and
    # its data, so solve and the solver compile once per optimizer, not once per
    # mechanism. The objective must still return loss! and allocate nothing, for a Vector
    # and for the column view of its population that CMA-ES passes.
    function _testhelper_optimization_problem(fp)
        n = length(parameters(fp.mechanism))
        EnzymeRates.Fitting._optimization_problem(fp, zeros(n), fill(-15.0, n),
                                                  fill(15.0, n))
    end
    objective_allocs(prob, x) = @allocated prob.f(x, prob.p)

    @testset "Optimizer boundary" begin
        data_uni = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        fp_uni = FittingProblem(uni_uni, data_uni, Keq_val;
                                reference=_testhelper_row_reference(data_uni))
        data_allo = _testhelper_make_synthetic_data(allo, allo_params, allo_concs)
        fp_allo = FittingProblem(allo, data_allo, Keq_val;
                                 reference=_testhelper_row_reference(data_allo))

        @testset "one problem type for every mechanism" begin
            @test typeof(_testhelper_optimization_problem(fp_uni)) ===
                  typeof(_testhelper_optimization_problem(fp_allo))
        end

        @testset "objective returns loss! and allocates nothing" begin
            for fp in (fp_uni, fp_allo)
                prob = _testhelper_optimization_problem(fp)
                population = randn(length(prob.u0), 3)
                x = population[:, 2]
                for x_in in (x, view(population, :, 2))
                    @test prob.f(x_in, prob.p) == EnzymeRates.Fitting.loss!(x, fp)
                    objective_allocs(prob, x_in)  # warmup
                    @test objective_allocs(prob, x_in) == 0
                end
            end
        end
    end

    # ── Speed benchmark ───────────────────────────────────────────────────────
    @testset "Speed" begin
        # Build a larger mechanism: Ordered Bi-Bi
        # Decomposed ordered bi-bi: the central complex EAB↔EPQ becomes an
        # explicit iso step EAB <--> EPQ. 5 steps total (was 4 in
        # the lumped-central-complex form); fitter recovers k1-k5 params.
        ordered_bi_bi = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EAB
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end

        bb_Keq = 1.5
        bb_params_syms = parameters(ordered_bi_bi)
        # Generate random params
        param_vals = ntuple(i -> 1.0 + 9.0 * rand(), length(bb_params_syms))
        bb_true_params = NamedTuple{bb_params_syms}(param_vals)

        # 500 synthetic datapoints
        n_points = 500
        concs_list = [(A = 0.1 + 9.9*rand(), B = 0.1 + 9.9*rand(),
                       P = 0.1 + 9.9*rand(), Q = 0.1 + 9.9*rand()) for _ in 1:n_points]
        groups = [string("G", div(i-1, 50)+1) for i in 1:n_points]

        bb_rate = compile_rate_equation(ordered_bi_bi)
        rates = [bb_rate(c, bb_true_params, bb_Keq) for c in concs_list]
        met_names_bb = metabolites(ordered_bi_bi)
        data = (
            group = groups,
            Rate = rates,
            (mn => [c[mn] for c in concs_list] for mn in met_names_bb)...
        )
        fp = FittingProblem(ordered_bi_bi, data, bb_Keq;
                            reference=_testhelper_row_reference(data))

        x = randn(length(parameters(ordered_bi_bi)))
        EnzymeRates.Fitting.loss!(x, fp)  # warmup/compile

        # Minimum over several batches defeats the GC/scheduling inflation a
        # single mean suffers (matches _testhelper_test_rate_equation_performance); the
        # accumulator's finiteness check keeps the calls from being elided.
        best_us = Inf
        acc = 0.0
        for _ in 1:5
            acc = 0.0
            t = @elapsed for _ in 1:2000; acc += EnzymeRates.Fitting.loss!(x, fp); end
            best_us = min(best_us, t / 2000 * 1e6)
        end
        isfinite(acc) || error("loss! produced a non-finite result")
        # 500 datapoints; ~5 μs typical local, min-of-batches strips CI noise
        @test best_us < 70
    end

    # ── fit_rate_equation's time beyond loss! ─────────────────────────────────
    # A real fit is loss! times the solver's evaluations, plus the solver's own work and a
    # fixed cost per restart. Counting the evaluations makes the bound independent of how
    # long the solver searches, and timing both in one process cancels machine speed.
    # Locally a fit takes 1.1-1.25x its loss! calls on these tables, so a cost per evaluation
    # fails the bound only at about half of loss! itself. The allocation checks above
    # catch a smaller cost that allocates.
    @testset "a fit costs about its loss! calls" begin
        ordered_bi_bi = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EAB
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end
        allosteric_bi_bi = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            oligomeric_state: 2
            steps: begin
                E + A <--> EA            :: EqualAI
                EA + B <--> EAB      :: EqualAI
                EAB <--> EPQ       :: OnlyA
                EPQ <--> EQ + P      :: EqualAI
                EQ <--> E + Q            :: EqualAI
            end
            regulatory_site(2): R::OnlyI
        end
        n_restarts = 5
        for (name, m) in (("ordered bi-bi", ordered_bi_bi),
                          ("allosteric bi-bi", allosteric_bi_bi))
            rng = Random.MersenneTwister(5)
            law_params = NamedTuple{parameters(m)}(
                Tuple(1.0 .+ 9.0 .* rand(rng, length(parameters(m)))))
            mets = Tuple(metabolites(m))
            law_concs = [NamedTuple{mets}(Tuple(0.1 .+ 9.9 .* rand(rng, length(mets))))
                         for _ in 1:500]
            groups = [string("G", div(i - 1, 50) + 1) for i in 1:500]
            data = _testhelper_make_synthetic_data(m, law_params, law_concs; groups,
                                                   Keq = 1.5)
            fp = FittingProblem(m, data, 1.5;
                                reference = _testhelper_row_reference(data))
            counting = _testhelper_CountingOpt(CMAEvolutionStrategyOpt(), Ref(0))
            fit_rate_equation(fp, counting; n_restarts = 1)  # warmup/compile
            x = randn(rng, length(parameters(m)))
            EnzymeRates.Fitting.loss!(x, fp)  # warmup/compile
            t_loss = minimum(1:10) do _
                @elapsed(for _ in 1:200; EnzymeRates.Fitting.loss!(x, fp); end) / 200
            end
            attempts = map(1:5) do _
                counting.evals[] = 0
                Random.seed!(11)
                t = @elapsed fit_rate_equation(fp, counting; n_restarts)
                (; t, evals = counting.evals[])
            end
            (; t, evals) = argmin(a -> a.t / a.evals, attempts)
            # 3 ms is about 2-7x a restart's fixed cost on an aarch64 development machine.
            allowed = 1.5 * evals * t_loss + n_restarts * 3e-3
            @info "fit of the $name: $evals evaluations, " *
                  "$(round(1e6 * t_loss; digits = 2)) µs per loss!, " *
                  "$(round(t; digits = 3)) s, allowed $(round(allowed; digits = 3)) s"
            @test t <= allowed
        end
    end

    # ── fit_rate_equation: relative fits come back at the reference turnover ──
    # The scale applies to every fitted point, so the assertions hold at any fit quality
    # and a half-second budget is enough.
    @testset "relative fits come back at the reference turnover" begin
        using OptimizationBBO
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
            (S = 0.5, P = 0.5), (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        opt = BBO_adaptive_de_rand_1_bin_radiuslimited()

        # The first data row (best of 2 restarts) and the fifth.
        rate = compile_rate_equation(uni_uni)
        for (row, n_restarts) in ((1, 2), (5, 1))
            reference = _testhelper_row_reference(data; row)
            fp = FittingProblem(uni_uni, data, Keq_val; reference)
            res = _testhelper_fit_capturing_convergence() do
                fit_rate_equation(fp, opt; n_restarts, maxtime=0.5)
            end
            @test keys(res) == (:params, :loss, :retcode, :at_bound)
            @test keys(res.params) == parameters(uni_uni)
            @test isfinite(res.loss)
            @test res.retcode isa Symbol
            concs, turnover = only(reference)
            @test rate(concs, res.params, Keq_val) ≈ turnover rtol=1e-9
        end
    end

    @testset "a fit's params hold the fitted names" begin
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        fit = fit_rate_equation(FittingProblem(uni_uni, data, 2.0),
                                CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 50)
        @test keys(fit.params) == parameters(uni_uni)
        # Extra names in params, Keq among them, are ignored, and the Keq argument is
        # used.
        @test EnzymeRates.test_rate_equation(
                  uni_uni, merge(fit.params, (Keq = 99.0,)), data, 2.0) ==
              EnzymeRates.test_rate_equation(uni_uni, fit.params, data, 2.0)
    end

    # ── fit_rate_equation: the scale on known raw constants ───────────────────
    @testset "reference scale, failed reference and absolute mode on known constants" begin
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        pn = parameters(uni_uni)
        unscaled = NamedTuple{pn}((1.0, 1.0, 1.0))
        _testhelper_stub_fit(u; kwargs...) = fit_rate_equation(
            FittingProblem(uni_uni, data, Keq_val; kwargs...),
            _testhelper_ConstantStubOpt(u); n_restarts=1)
        # Every raw constant is exp(0) = 1, so at Keq = 2 the product rebinding rate is
        # 1/2 and v = (S − P/2)/(2 + S + P/2): 1/3 at S = 1 alone. A turnover of 2 there
        # scales every fitted constant, each a rate constant, by 6.
        rel = _testhelper_stub_fit(0.0; reference=[(S = 1.0,) => 2.0])
        @test keys(rel) == (:params, :loss, :retcode, :at_bound)
        @test keys(rel.params) == pn
        @test all(k -> isapprox(rel.params[k], 6.0; rtol=1e-12), pn)
        @test rel.retcode === :Success && isfinite(rel.loss)
        # The product alone runs the law backwards, against a forward reference turnover:
        # no factor matches it, so the fit fails without throwing and keeps its constants.
        wrong_sign = _testhelper_stub_fit(0.0; reference=[(P = 1.0,) => 1.0])
        @test wrong_sign.params == unscaled
        @test wrong_sign.loss == Inf && wrong_sign.retcode === :ReferenceFailed
        # At exp(-700) every product of two constants underflows to zero, so the law
        # predicts zero turnover at the reference: the fit fails the same way.
        tiny = _testhelper_stub_fit(-700.0; reference=[(S = 1.0,) => 2.0])
        @test all(k -> tiny.params[k] == exp(-700.0), pn)
        @test tiny.loss == Inf && tiny.retcode === :ReferenceFailed
        # With no finite restart the fit returns before the reference step, so its retcode
        # says so even where the reference would fail.
        no_finite = fit_rate_equation(
            FittingProblem(uni_uni, data, Keq_val; reference=[(P = 1.0,) => 1.0]),
            _testhelper_ConstantStubOpt(0.0; infinite=true); n_restarts=1)
        @test no_finite.params == unscaled
        @test no_finite.loss == Inf && no_finite.retcode === :NoFiniteLoss
        # Absolute mode returns the fitted constants as they are.
        absolute_fit = _testhelper_stub_fit(0.0; absolute=true)
        @test keys(absolute_fit) == (:params, :loss, :retcode, :at_bound)
        @test absolute_fit.params == unscaled
    end

    # ── solver-option forwarding (named commons + solver_kwargs) ──
    @testset "solver kwarg forwarding" begin
        using OptimizationCMAEvolutionStrategy
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
        ]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data, Keq_val;
                            reference=_testhelper_row_reference(data))

        # Default (empty) solver_kwargs runs on a solver that rejects unknown
        # options — no solver-specific option is force-injected.
        res = fit_rate_equation(fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=1.0)
        @test isfinite(res.loss)        # a finite loss means a fit actually ran

        # solver_kwargs is forwarded verbatim: an option no optimizer
        # recognizes surfaces as an error (proves the bag reaches `solve`);
        # a bogus name keeps this independent of any real option's support.
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=1.0,
            solver_kwargs=(; not_a_real_solver_option=1))

        # Merge semantics: the same key in both a named common option and
        # solver_kwargs does NOT raise a duplicate-keyword error (a naive
        # double-splat would). The override value taking effect is guaranteed
        # by `merge(common, solver_kwargs)` in fit_rate_equation; here we
        # assert the call succeeds and produces a finite-loss fit.
        res2 = fit_rate_equation(
            fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=60.0, solver_kwargs=(; maxtime=1.0))
        @test isfinite(res2.loss)

        # Clean break: popsize/verbose are no longer accepted named kwargs.
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt(); n_restarts=1, maxtime=1.0, popsize=200)
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt(); n_restarts=1, maxtime=1.0, verbose=-9)
    end

    # ── maxtime forwarded from fit_rate_equation to Optimization.solve (§6) ─
    @testset "maxtime forwarded to Optimization.solve" begin
        # A stub optimizer that records the `maxtime` kwarg it receives, so the
        # forwarding path (fit_rate_equation -> the `common = (; maxtime,
        # maxiters)` merge -> Optimization.solve) can be asserted end-to-end
        # without depending on a real solver's behavior.
        mutable struct _testhelper_MaxtimeStubOpt
            maxtime_seen::Union{Nothing, Real}
        end
        _testhelper_MaxtimeStubOpt() = _testhelper_MaxtimeStubOpt(nothing)
        Optimization.allowsbounds(::_testhelper_MaxtimeStubOpt) = true
        function Optimization.SciMLBase.__solve(
                prob::Optimization.OptimizationProblem,
                opt::_testhelper_MaxtimeStubOpt; kwargs...)
            opt.maxtime_seen = kwargs[:maxtime]
            u = zeros(length(prob.u0))
            cache = DefaultOptimizationCache(prob.f, prob.p)
            build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
        end

        concs_list = [(S = 1.0, P = 0.1), (S = 2.0, P = 0.1)]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data, Keq_val;
                            reference=_testhelper_row_reference(data))

        stub = _testhelper_MaxtimeStubOpt()
        fit_rate_equation(fp, stub; n_restarts=1, maxtime=1.23)
        @test stub.maxtime_seen == 1.23
    end

    # ── fit_rate_equation compiles no NamedTuple indexing per mechanism ───────
    # fit_rate_equation takes its FittingProblem @nospecialize and passes the fitted
    # parameters around as a NamedTuple keyed by their names. Indexing it by name would
    # dispatch at run time: each fit whose parameter names no earlier fit used compiles
    # a getindex method instance for them and infers Base's slow NamedTuple fallback.
    # Inferred only as a Tuple, the names also match the getindex method of Static.jl,
    # which the full suite loads before this file. Counts the getindex method instances,
    # inferred or compiled, on a NamedTuple keyed by any of `names`.
    function _testhelper_namedtuple_getindex_instances(names)
        keyed(T) = T isa DataType && T <: NamedTuple && !isdisjoint(fieldnames(T), names)
        count(mi -> mi.specTypes isa DataType && any(keyed, mi.specTypes.parameters),
              (mi for m in methods(getindex) for mi in Base.specializations(m)))
    end

    @testset "fit_rate_equation compiles no NamedTuple indexing per mechanism" begin
        using OptimizationCMAEvolutionStrategy
        optimizer = CMAEvolutionStrategyOpt()
        # A first fit compiles fit_rate_equation for this optimizer and these options.
        data_warm = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        fp_warm = FittingProblem(uni_uni, data_warm, Keq_val;
                                 reference=_testhelper_row_reference(data_warm))
        fit_rate_equation(fp_warm, optimizer; n_restarts=1, maxiters=10)
        # Metabolite names no other test uses, so no earlier test has fit this mechanism.
        unseen = @enzyme_mechanism begin
            substrates: Sidx
            products:   Pidx
            steps: begin
                E + Sidx <--> ESidx
                ESidx <--> E + Pidx
            end
        end
        data_unseen = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
                       Sidx = [1.0, 2.0, 3.0, 4.0], Pidx = [0.1, 0.2, 0.3, 0.4])
        fp_unseen = FittingProblem(unseen, data_unseen, Keq_val;
                                   reference=_testhelper_row_reference(data_unseen))
        unseen_names = parameters(unseen)
        n = _testhelper_namedtuple_getindex_instances(unseen_names)
        fit_rate_equation(fp_unseen, optimizer; n_restarts=1, maxiters=10)
        @test _testhelper_namedtuple_getindex_instances(unseen_names) == n
    end

    # ── Validation errors ──────────────────────────────────────────────────────
    # The table and reference checks are tested in test/Data/test_rate_data.jl.
    @testset "Validation errors reach the caller" begin
        data_zero = (group = ["G1"], Rate = [0.0], S = [1.0], P = [0.1])
        @test_throws(ErrorException("Zero rate at row 1: log(0) is undefined"),
                     FittingProblem(uni_uni, data_zero, 1.0;
                                    reference=[(S = 1.0, P = 0.1) => 1.0]))
        # A metabolite of the mechanism named after a required column
        rate_named = @enzyme_mechanism begin
            substrates: Rate
            products:   P
            steps: begin
                E + Rate <--> ERate
                ERate <--> E + P
            end
        end
        @test_throws(
            ErrorException("Metabolite Rate has the name of the required Rate column; " *
                           "rename the metabolite"),
            FittingProblem(rate_named, (group = ["G1"], Rate = [1.0], P = [0.1]), 1.0;
                           reference=[(P = 0.1,) => 1.0]))
    end

end

# ── The scales of a fit ───────────────────────────────────────────────────────────────

using EnzymeRates: test_rate_equation
using OptimizationCMAEvolutionStrategy

@testset "Reference turnover scale" begin

    # Uni-uni with a steady-state binding; at Keq = 2 these constants give kcat = 5 and
    # the Michaelis constant Km. Fitted: k_ES_to_E_P, k_ES_to_E_S, k_E_S_to_ES.
    mm = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    _testhelper_mm_params(Km) = (k_ES_to_E_P = 5.0, k_ES_to_E_S = 5.0,
                                 k_E_S_to_ES = 10.0 / Km)
    mm_concs = [(S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
                (S = 5.0, P = 0.1), (S = 10.0, P = 0.1), (S = 1.0, P = 0.5)]
    mm_rate = compile_rate_equation(mm)
    mm_data = (group = fill("G1", length(mm_concs)),
               Rate = [mm_rate(c, _testhelper_mm_params(1.0), 2.0) for c in mm_concs],
               S = [c.S for c in mm_concs], P = [c.P for c in mm_concs])

    @testset "relative mode accepts no reference" begin
        fp_none = FittingProblem(mm, mm_data, 2.0)
        @test !fp_none.rate_data.absolute && isempty(fp_none.rate_data.reference)
        # Absolute data set the scale themselves: no reference is needed, and one passed
        # is an error.
        fp_abs = FittingProblem(mm, mm_data, 2.0; absolute = true)
        @test fp_abs.rate_data.absolute && isempty(fp_abs.rate_data.reference)
        @test_throws(ErrorException("reference applies only to relative mode; absolute " *
                                    "turnover data set the rate scale themselves"),
                     FittingProblem(mm, mm_data, 2.0; absolute = true,
                                    reference = [(S = 10.0,) => 3.0]))
        # The problem stores the anchors with Float64 turnovers; integers convert.
        fp = FittingProblem(mm, mm_data, 2.0; reference = [(S = 10,) => 3])
        @test fp.rate_data.reference == [(S = 10,) => 3.0]
        @test fp.rate_data.reference isa Vector{Pair{NamedTuple, Float64}}
        # The concrete-mechanism constructor passes the reference on.
        fp_conc = FittingProblem(mm, mm_data, 2.0;
                                 reference = [(S = 10.0,) => 3.0])
        @test fp_conc.rate_data.reference == [(S = 10.0,) => 3.0]
    end

    @testset "reference anchors are validated at construction" begin
        _testhelper_problem(reference) =
            FittingProblem(mm, mm_data, 2.0; reference)
        @test_throws(ErrorException("reference must be a vector of (concentrations => " *
                                    "turnover) pairs; got (S = 10.0,) => 3.0"),
                     _testhelper_problem((S = 10.0,) => 3.0))
        @test_throws(ErrorException("reference is empty; pass at least one " *
                                    "(concentrations => turnover) anchor"),
                     _testhelper_problem(Pair{NamedTuple, Float64}[]))
        @test_throws(ErrorException("reference: concentrations must be a NamedTuple; " *
                                    "got (10.0,)"),
                     _testhelper_problem([(10.0,) => 3.0]))
        @test_throws(ErrorException("reference: X is not a metabolite of the mechanism " *
                                    "(metabolites: S, P)"),
                     _testhelper_problem([(X = 1.0,) => 3.0]))
        @test_throws(ErrorException("reference: concentration of S must be a finite " *
                                    "number ≥ 0; got -1.0"),
                     _testhelper_problem([(S = -1.0,) => 3.0]))
        for (bad, shown) in ((0.0, "0.0"), (NaN, "NaN"), (Inf, "Inf"))
            @test_throws(ErrorException("reference: turnover must be a finite, nonzero " *
                                        "number; got $shown"),
                         _testhelper_problem([(S = 10.0,) => bad]))
        end
    end

    @testset "relative fit matches the reference turnover" begin
        # The scale holds at any fit quality, so a short fit is enough.
        fp = FittingProblem(mm, mm_data, 2.0; reference = [(S = 10.0,) => 3.0])
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 50)
        @test keys(fit) == (:params, :loss, :retcode, :at_bound)
        @test isfinite(fit.loss)
        @test mm_rate((S = 10.0, P = 0.0), fit.params, 2.0) ≈ 3.0 rtol = 1e-8
        # Several anchors combine by the geometric mean of measured over predicted
        # turnover: after the scale, the log ratios average to zero.
        anchors = [(S = 10.0,) => 3.0, (S = 1.0,) => 2.0]
        fit2 = fit_rate_equation(
            FittingProblem(mm, mm_data, 2.0; reference = anchors),
            CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 50)
        log_ratios = [log(t / mm_rate(merge((S = 0.0, P = 0.0), c), fit2.params, 2.0))
                      for (c, t) in anchors]
        @test sum(log_ratios) / length(log_ratios) ≈ 0 atol = 1e-12
    end

    @testset "an unlisted metabolite sits at zero in an anchor" begin
        # Random-order bi-bi: with only A listed, B = P = Q = 0 and the turnover is 0,
        # which no scale reproduces, so the fit fails without throwing.
        bibi = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B <--> EAB
                EB + A <--> EAB
                EAB <--> EPQ
                EPQ <--> EP + Q
                EPQ <--> EQ + P
                EP <--> E + P
                EQ <--> E + Q
            end
        end
        data = (group = ["G1", "G1", "G1"], Rate = [0.5, 0.8, 1.0], A = [1.0, 2.0, 4.0],
                B = [1.0, 1.0, 1.0], P = [0.1, 0.1, 0.1], Q = [0.1, 0.1, 0.1])
        fp = FittingProblem(bibi, data, 2.0; reference = [(A = 1.0,) => 1.0])
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 20)
        @test fit.loss == Inf
        @test fit.retcode === :ReferenceFailed
        @test keys(fit.params) == parameters(bibi)
        @test all(v -> isfinite(v) && v > 0, values(fit.params))
    end

    @testset "reference scaling leaves L and binding constants" begin
        # MWC uni-uni, catalysis in the active conformation only; X binds the inactive
        # conformation, Y both with different constants.
        two_site = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::OnlyI
            regulatory_site(1): Y::NonequalAI
        end
        p, Keq = _testhelper_random_independent_params(two_site;
                                                       rng = Random.MersenneTwister(3))
        concs = (S = 2.0, P = 0.5, X = 0.3, Y = 1.5)
        predicted = compile_rate_equation(two_site)(concs, p, Keq)
        q = EnzymeRates.Fitting._scale_rate_constants(two_site, p, 7.0 / predicted)
        @test keys(q) == keys(p)
        rate_constants = [k for k in keys(p) if startswith(String(k), "k_")]
        @test !isempty(rate_constants)
        @test :L in keys(p)
        # Every k_… constant changes by one factor, measured over predicted turnover.
        for k in rate_constants
            @test q[k] / p[k] ≈ 7.0 / predicted rtol = 1e-12
        end
        # Every other entry — binding constants, L, regulator constants — stays exactly as
        # it was.
        for k in keys(p)
            k in rate_constants || @test q[k] === p[k]
        end
    end

    @testset "default scale sets the fastest data row to 1" begin
        fp = FittingProblem(mm, mm_data, 2.0)
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 50)
        @test isfinite(fit.loss) && fit.retcode !== :ReferenceFailed
        peak = maximum(abs(mm_rate(c, fit.params, 2.0)) for c in mm_concs)
        @test peak ≈ 1 rtol = 1e-10
    end

    @testset "default scale fails softly" begin
        zero_data = (group = ["G1", "G1", "G1"], Rate = [0.5, 0.8, 1.0],
                     S = [0.0, 0.0, 0.0], P = [0.0, 0.0, 0.0])
        fit = fit_rate_equation(FittingProblem(mm, zero_data, 2.0),
                                CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 20)
        @test fit.retcode === :ReferenceFailed
        @test fit.loss == Inf
        @test keys(fit.params) == parameters(mm)
        @test all(v -> isfinite(v) && v > 0, values(fit.params))
    end

    @testset "test_rate_equation is the fitting loss" begin
        fit = fit_rate_equation(FittingProblem(mm, mm_data, 2.0),
                                CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 50)
        @test keys(fit.params) == parameters(mm)
        @test test_rate_equation(mm, fit.params, mm_data, 2.0) ≈ fit.loss
        # Names beyond the mechanism's fitted parameters are ignored.
        @test test_rate_equation(mm, merge(fit.params, (Keq = 99.0, E_total = 1.0)),
                                 mm_data, 2.0) ≈ fit.loss
        # The concrete mechanism gives the same loss.
        @test test_rate_equation(mm, fit.params, mm_data, 2.0) ≈
              fit.loss
        # Absolute mode scores the magnitude: the true constants fit the data exactly,
        # constants ten times as fast do not.
        truth = _testhelper_mm_params(1.0)
        exact = test_rate_equation(mm, truth, mm_data, 2.0; absolute = true)
        @test exact ≈ 0 atol = 1e-20
        fast = NamedTuple{keys(truth)}(10 .* values(truth))
        @test test_rate_equation(mm, fast, mm_data, 2.0; absolute = true) > 1
        without = Base.structdiff(fit.params, (k_ES_to_E_S = 0.0,))
        @test_throws ErrorException("params is missing k_ES_to_E_S") test_rate_equation(
            mm, without, mm_data, 2.0)
    end

end

@testset "each constant starts at the typical value of its units" begin
    # Each spec's log centres at the typical rate 3 and concentrations 10^-j of
    # `_testhelper_typical_table`: the log of 3^rate × Π (10^-j)^power for the units the
    # rule reads from each constant's name (`_testhelper_units_from_name`).
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        fp = FittingProblem(m, _testhelper_typical_table(m), 2.0)
        log_s = [log(10.0^-j) for j in eachindex(metabolites(m))]
        expected = map(parameters(m)) do x
            rate, concentrations = _testhelper_units_from_name(x, metabolites(m))
            rate * log(3.0) + concentrations' * log_s
        end
        @test EnzymeRates.Fitting._log_centres(fp) ≈ collect(expected) atol = 1e-12
    end
end

@testset "a relative fit's scale multiplies exactly the rate constants" begin
    # The steady-state rate constants, named k_…, take the factor; every other fitted
    # constant (a K, a regulator's constant, L) keeps its value.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        ks = parameters(m)
        params = NamedTuple{ks}(Tuple(Float64.(eachindex(ks))))
        scaled = EnzymeRates.Fitting._scale_rate_constants(m, params, 10.0)
        @test keys(scaled) == ks
        @test [scaled[k] / params[k] for k in ks] ==
              [startswith(String(k), "k_") ? 10.0 : 1.0 for k in ks]
    end
end

# ── The fit's search around the data ──────────────────────────────────────────────────

# An optimizer stub that records each problem's start and box, and returns the start.
struct _testhelper_RecordingStubOpt
    seen::Vector{NamedTuple}
end
_testhelper_RecordingStubOpt() = _testhelper_RecordingStubOpt(NamedTuple[])
Optimization.allowsbounds(::_testhelper_RecordingStubOpt) = true
function Optimization.SciMLBase.__solve(prob::Optimization.OptimizationProblem,
        opt::_testhelper_RecordingStubOpt; kwargs...)
    push!(opt.seen, (u0 = copy(prob.u0), lb = copy(prob.lb), ub = copy(prob.ub)))
    build_solution(DefaultOptimizationCache(prob.f, prob.p), opt, prob.u0,
                   prob.f(prob.u0, prob.p); retcode = ReturnCode.Success)
end

# An optimizer stub that returns the fixed point `u` at a fixed objective, which may be Inf.
struct _testhelper_FixedPointStubOpt
    u::Vector{Float64}
    objective::Float64
end
Optimization.allowsbounds(::_testhelper_FixedPointStubOpt) = true
function Optimization.SciMLBase.__solve(prob::Optimization.OptimizationProblem,
        opt::_testhelper_FixedPointStubOpt; kwargs...)
    build_solution(DefaultOptimizationCache(prob.f, prob.p), opt, copy(opt.u),
                   opt.objective; retcode = ReturnCode.Success)
end

# A one-row table of the metabolites `mets` at the concentrations `s`, a Dict, and the rate
# `r`: its scales, the geometric means a fit's centres read, are `s` and `r`.
_testhelper_scales_table(mets, s, r) =
    (; group = ["g"], Rate = [r], (x => [s[x]] for x in mets)...)

# Tests that the centres of `mechanism` move with the law: scaling each metabolite by its
# own factor, `Keq` by those factors to their stoichiometry and each fitted parameter by the
# move of its centre multiplies the rate by exactly the factor on the rate reference.
function _testhelper_check_centres_match_law(mechanism, rng)
    rate = compile_rate_equation(mechanism)
    mets = Tuple(metabolites(mechanism))
    s = Dict(x => exp(randn(rng)) for x in mets)
    f = Dict(x => exp(2 * randn(rng)) for x in mets)
    sf = Dict(x => s[x] * f[x] for x in mets)
    r, g = exp(randn(rng)), exp(2 * randn(rng))
    keq_factor = prod(f[x] for x in _testhelper_products(mechanism)) /
                 prod(f[x] for x in _testhelper_substrates(mechanism))
    centres(scales, rate) = EnzymeRates.Fitting._log_centres(
        FittingProblem(mechanism, _testhelper_scales_table(mets, scales, rate), 1.0))
    move = centres(sf, g * r) .- centres(s, r)
    params, Keq = _testhelper_random_independent_params(mechanism; rng)
    fitted = parameters(mechanism)
    scaled = merge(params, NamedTuple{fitted}(Tuple(params[k] * exp(move[i])
                                                    for (i, k) in enumerate(fitted))))
    concs = NamedTuple{mets}(Tuple(exp(randn(rng)) for _ in mets))
    concs_scaled = NamedTuple{mets}(Tuple(concs[x] * f[x] for x in mets))
    @test rate(concs_scaled, scaled, Keq * keq_factor) ≈ g * rate(concs, params, Keq) rtol =
          1e-9
end

# A noiseless table of `mechanism` at the fitted constants `params` and `Keq`: one row, in
# one group, per point of the grid `axes`, a NamedTuple of each metabolite's concentrations.
function _testhelper_noiseless_table(mechanism, params, Keq, axes)
    mets = keys(axes)
    rows = vec(collect(Iterators.product(values(axes)...)))
    rate = compile_rate_equation(mechanism)
    (; group = fill("G1", length(rows)),
       Rate = [rate(NamedTuple{mets}(row), params, Keq) for row in rows],
       NamedTuple{mets}(Tuple(getindex.(rows, i) for i in eachindex(mets)))...)
end

@testset "Fit search around the data" begin
    log_centres = EnzymeRates.Fitting._log_centres
    concentration_scales = EnzymeRates.Fitting._concentration_scales
    rate_scale = EnzymeRates.Fitting._rate_scale
    optimizer = CMAEvolutionStrategyOpt()

    # Fitted: k_E_S_to_ES, k_ES_to_E_S, k_ES_to_E_P, all of which the data determine.
    steady_state = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    steady_state_params = (k_E_S_to_ES = 1e7, k_ES_to_E_S = 50.0, k_ES_to_E_P = 200.0)
    steady_state_axes = (S = [1e-6, 1e-5, 1e-4, 1e-3], P = [1e-5, 3e-4, 1e-3])

    # Fitted: K_ES_to_E_S, k_ES_to_E_P.
    rapid_binding = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end

    @testset "Facts match the law" begin
        rng = Random.MersenneTwister(97)
        @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
            _testhelper_check_centres_match_law(spec.mechanism, rng)
        end
        inline_laws = (
            two_in_steady_state = @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A + B <--> EAB
                    EAB <--> EPQ
                    EPQ <--> EQ + P
                    EQ <--> E + Q
                end
            end),
            two_in_rapid = @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A + B ⇌ EAB
                    EAB <--> EPQ
                    EPQ <--> EQ + P
                    EQ <--> E + Q
                end
            end),
            two_same_in_rapid = @enzyme_mechanism(begin
                substrates: A
                products: P
                steps: begin
                    E + A <--> F + P
                    F + A <--> E + P
                    F + A + A ⇌ FAA
                end
            end),
            two_out_steady_state = @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    EA + B <--> EAB
                    EAB <--> EPQ
                    EPQ <--> E + P + Q
                end
            end),
            modifier_copy = @enzyme_mechanism(begin
                substrates: S
                products: P
                steps: begin
                    E + S ⇌ ES
                    ES <--> EP
                    EP ⇌ E + P
                    EP + S ⇌ EPSmod
                end
            end),
            modifier_copy_steady_state = @enzyme_mechanism(begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> ES
                    ES <--> EP
                    EP <--> E + P
                    EP + S <--> EPSmod
                end
            end),
        )
        @testset "$label" for (label, law) in pairs(inline_laws)
            _testhelper_check_centres_match_law(law, rng)
        end
    end

    @testset "Each parameter's centre" begin
        # The centres at the table's scales `s` and `r`, by name.
        centres(m, s, r) = Dict(zip(parameters(m), log_centres(FittingProblem(
            m, _testhelper_scales_table(metabolites(m), s, r), 2.0))))
        s = Dict(:S => 10.0, :P => 100.0, :R => 1000.0)
        r = 4.0
        allosteric = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES     :: EqualAI
                ES <--> EP   :: OnlyA
                EP ⇌ E + P     :: EqualAI
            end
            regulatory_site(1): R::OnlyI
        end
        expected = Dict(:K_ES_to_E_S => log(10.0), :K_EP_to_E_P => log(100.0),
                        :k_A_ES_to_EP => log(4.0), :K_I_Rreg => log(1000.0), :L => 0.0)
        got = centres(allosteric, s, r)
        @test keys(got) == keys(expected)
        for (k, v) in expected
            @test got[k] ≈ v
        end

        # A Theorell–Chance step's rate constant reads what it takes up, and its
        # equilibrium constant, fitted where an isomerization closes the cycle, reads what
        # it gives off over what it takes up.
        s = Dict(:A => 2.0, :B => 10.0, :P => 100.0, :Q => 3.0)
        tc_steady = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B <--> EQ + P
                EQ ⇌ E + Q
            end
        end
        @test centres(tc_steady, s, r)[:k_EA_B_to_EQ_P] ≈ log(4.0) - log(10.0)
        tc_rapid = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> EA
                EA + B ⇌ EQ + P
                EQ <--> FQ
                FQ <--> F + Q
                F <--> E
            end
        end
        @test centres(tc_rapid, s, r)[:K_EA_B_to_EQ_P] ≈ log(100.0) - log(10.0)
    end

    @testset "Issue #97 regression" begin
        K_true = 0.2e-6
        table = _testhelper_noiseless_table(rapid_binding,
            (K_ES_to_E_S = K_true, k_ES_to_E_P = 5.0), 1000.0,
            (S = [2e-8, 1e-7, 2e-7, 5e-7, 2e-6, 1e-5], P = [1e-7, 1e-5]))
        for absolute in (false, true)
            Random.seed!(97)
            fp = FittingProblem(rapid_binding, table, 1000.0; absolute)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 5.0)
            @test fit.params.K_ES_to_E_S ≈ K_true rtol = 0.01
            @test abs(log(fit.params.K_ES_to_E_S) - (-15.0)) > 0.01
        end
    end

    @testset "Unit independence" begin
        # Each law, its true constants in M and the power of the concentration unit that
        # each fitted constant carries.
        laws = (
            (steady_state, steady_state_params,
             (k_E_S_to_ES = -1, k_ES_to_E_S = 0, k_ES_to_E_P = 0)),
            (rapid_binding, (K_ES_to_E_S = 2e-5, k_ES_to_E_P = 200.0),
             (K_ES_to_E_S = 1, k_ES_to_E_P = 0)),
        )
        for (law, params, powers) in laws, absolute in (true, false)
            table = _testhelper_noiseless_table(law, params, 2.0, steady_state_axes)
            fits = map((1.0, 1e3, 1e6)) do unit
                scaled = merge(table, (S = table.S .* unit, P = table.P .* unit))
                Random.seed!(11)
                fp = FittingProblem(law, scaled, 2.0; absolute)
                fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 10.0).params
            end
            for (fit, unit) in zip(fits, (1.0, 1e3, 1e6)), k in keys(powers)
                @test log(fit[k]) - powers[k] * log(unit) ≈ log(fits[1][k]) atol = 1e-3
            end
        end
    end

    @testset "Both modes" begin
        table = _testhelper_noiseless_table(steady_state, steady_state_params, 2.0,
                                            steady_state_axes)
        problems = (
            FittingProblem(steady_state, table, 2.0; absolute = true),
            FittingProblem(steady_state, table, 2.0;
                           reference = _testhelper_row_reference(table)),
        )
        for fp in problems
            Random.seed!(23)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 3, maxtime = 10.0)
            @test fit.loss < 1e-6
            for k in keys(steady_state_params)
                @test fit.params[k] ≈ steady_state_params[k] rtol = 0.01
            end
        end
    end

    @testset "Two-metabolite law" begin
        two_in = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A + B <--> EAB
                EAB <--> E + P + Q
            end
        end
        params = (k_E_A_B_to_EAB = 1e10, k_EAB_to_E_A_B = 100.0, k_EAB_to_E_P_Q = 300.0)
        @test issetequal(keys(params), parameters(two_in))
        table = _testhelper_noiseless_table(two_in, params, 20.0,
            (A = [1e-5, 3e-5, 1e-4, 3e-4], B = [1e-4, 1e-3], P = [1e-4, 1e-3],
             Q = [1e-4, 1e-3]))
        Random.seed!(29)
        fit = fit_rate_equation(FittingProblem(two_in, table, 20.0; absolute = true),
                                optimizer; n_restarts = 3, maxtime = 10.0)
        @test fit.loss < 1e-6
        for k in keys(params)
            @test fit.params[k] ≈ params[k] rtol = 0.01
        end
    end

    @testset "Always-zero metabolite" begin
        @testset "the fit runs through the positive concentrations' geometric mean" begin
            table = _testhelper_noiseless_table(rapid_binding,
                (K_ES_to_E_S = 1e-4, k_ES_to_E_P = 5.0), 2.0,
                (S = [1e-5, 1e-4, 1e-3], P = [0.0]))
            fp = FittingProblem(rapid_binding, table, 2.0)
            @test concentration_scales(fp)[:P] ≈ 1e-4 rtol = 1e-12
            Random.seed!(31)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 5.0)
            @test isfinite(fit.loss)
            @test fit.params.K_ES_to_E_S ≈ 1e-4 rtol = 0.01
        end

        @testset "the fallback pools the mechanism's metabolites' concentrations" begin
            two_in = @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A + B <--> EAB
                    EAB <--> E + P + Q
                end
            end
            # The other column, `Dose`, is no metabolite, so the fallback leaves it out.
            table = (group = fill("G1", 4), Rate = [1.0, 2.0, 3.0, 4.0],
                     A = [1.0, 0.0, 0.0, 0.0], B = [10.0, 10.0, 10.0, 10.0],
                     P = zeros(4), Q = zeros(4), Dose = fill(1e6, 4))
            s = concentration_scales(FittingProblem(two_in, table, 2.0))
            @test s[:A] == 1.0
            @test s[:B] ≈ 10.0 rtol = 1e-12
            @test s[:P] ≈ 10.0^0.8 rtol = 1e-12
            @test s[:Q] == s[:P]
        end

        @testset "a table with no positive concentration centres on 1" begin
            table = (group = fill("G1", 3), Rate = [0.5, 0.8, 1.0],
                     S = zeros(3), P = zeros(3))
            fp = FittingProblem(rapid_binding, table, 2.0)
            @test concentration_scales(fp) == Dict(:S => 1.0, :P => 1.0)
        end

        # A cross-validation fold builds its problem from its training rows, which can
        # hold a metabolite only at zero while the whole table does not.
        @testset "a training subset with a metabolite zero only in its rows" begin
            params = (K_ES_to_E_S = 1e-4, k_ES_to_E_P = 5.0)
            table = _testhelper_noiseless_table(rapid_binding, params, 2.0,
                (S = [1e-5, 1e-4, 1e-3], P = [0.0, 1e-2]))
            train = table.P .== 0.0
            @test any(>(0.0), table.P) && any(train)
            train_table = (group = table.group[train], Rate = table.Rate[train],
                           S = table.S[train], P = table.P[train])
            fp = FittingProblem(rapid_binding, train_table, 2.0)
            @test concentration_scales(fp)[:P] ≈ 1e-4 rtol = 1e-12
            whole = FittingProblem(rapid_binding, table, 2.0)
            @test concentration_scales(whole)[:P] ≈ 1e-2 rtol = 1e-12
            Random.seed!(37)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 5.0)
            @test isfinite(fit.loss)
            @test fit.params.K_ES_to_E_S ≈ 1e-4 rtol = 0.01
        end
    end

    @testset "Integer columns" begin
        integers = (group = fill("G1", 4), Rate = [3, 5, -7, 9],
                    S = [1, 2, 5, 10], P = [1, 1, 0, 2])
        floats = (group = integers.group, Rate = Float64.(integers.Rate),
                  S = Float64.(integers.S), P = Float64.(integers.P))
        fp_int = FittingProblem(rapid_binding, integers, 2.0)
        fp_float = FittingProblem(rapid_binding, floats, 2.0)
        @test log_centres(fp_int) == log_centres(fp_float)
        Random.seed!(41)
        fit = fit_rate_equation(fp_int, optimizer; n_restarts = 1, maxiters = 20)
        @test isfinite(fit.loss)
    end

    @testset "Search box" begin
        table = (group = ["G1", "G1", "G2", "G2"], Rate = [2.0, 8.0, -2.0, 8.0],
                 S = [1e-9, 1e-7, 1e-9, 1e-7], P = [1e-5, 1e-5, 1e-5, 1e-5])
        fp = FittingProblem(steady_state, table, 2.0)
        @test rate_scale(fp) ≈ 4.0 rtol = 1e-12
        s = concentration_scales(fp)
        @test s[:S] ≈ 1e-8 rtol = 1e-12
        @test s[:P] ≈ 1e-5 rtol = 1e-12
        names = parameters(steady_state)
        centre = log_centres(fp)
        expected = (k_E_S_to_ES = log(4.0) - log(1e-8), k_ES_to_E_S = log(4.0),
                    k_ES_to_E_P = log(4.0))
        @test centre ≈ [expected[k] for k in names]

        @testset "defaults sit 15 log units either side of the centre" begin
            stub = _testhelper_RecordingStubOpt()
            Random.seed!(43)
            fit_rate_equation(fp, stub; n_restarts = 300)
            @test length(stub.seen) == 300
            @test all(seen -> seen.lb == centre .- 15.0 && seen.ub == centre .+ 15.0,
                      stub.seen)
            @test all(seen -> all(seen.lb .<= seen.u0 .<= seen.ub), stub.seen)
            # Each restart draws anew, with a standard deviation of 2 around the centre.
            offsets = reduce(vcat, [seen.u0 .- centre for seen in stub.seen])
            @test allunique(offsets)
            @test abs(sum(offsets) / length(offsets)) < 0.3
            @test 1.8 < sqrt(sum(abs2, offsets) / length(offsets)) < 2.2
        end

        @testset "starts clamp into an explicit box, which holds absolute log values" begin
            stub = _testhelper_RecordingStubOpt()
            lb, ub = fill(-1.0, 3), fill(1.0, 3)
            Random.seed!(47)
            fit_rate_equation(fp, stub; n_restarts = 1, lb, ub)
            seen = only(stub.seen)
            @test seen.lb == lb && seen.ub == ub
            # The first constant's centre lies near 20, far above this box.
            @test seen.u0[findfirst(==(:k_E_S_to_ES), names)] == 1.0
            @test all(lb .<= seen.u0 .<= ub)
        end
    end

    @testset "At-bound flag" begin
        table = (group = ["G1", "G1", "G2", "G2"], Rate = [2.0, 8.0, -2.0, 8.0],
                 S = [1e-9, 1e-7, 1e-9, 1e-7], P = [1e-5, 1e-5, 1e-5, 1e-5])
        fp = FittingProblem(steady_state, table, 2.0)
        names = parameters(steady_state)
        centre = log_centres(fp)
        lb, ub = centre .- 15.0, centre .+ 15.0
        fit_at(u; objective = 1.0, problem = fp, kwargs...) = fit_rate_equation(
            problem, _testhelper_FixedPointStubOpt(u, objective); n_restarts = 1, kwargs...)

        @testset "within 0.01 log units of a default edge" begin
            for (edge, sign) in ((ub, -1), (lb, 1))
                u = copy(centre)
                u[2] = edge[2] + sign * 0.009
                @test fit_at(u).at_bound == (names[2],)
                u[2] = edge[2] + sign * 0.011
                @test fit_at(u).at_bound == ()
            end
            u = copy(centre)
            u[1], u[3] = ub[1], lb[3]
            @test fit_at(u).at_bound == (names[1], names[3])
        end

        @testset "an explicit box" begin
            lb_x, ub_x = fill(-1.0, 3), fill(1.0, 3)
            u = [0.0, 0.0, -1.0 + 0.009]
            @test fit_at(u; lb = lb_x, ub = ub_x).at_bound == (names[3],)
            u[3] = -1.0 + 0.011
            @test fit_at(u; lb = lb_x, ub = ub_x).at_bound == ()
        end

        @testset "a box that excludes the centre clamps the start to the edge" begin
            stub = _testhelper_RecordingStubOpt()
            lb_x, ub_x = fill(-1.0, 3), fill(1.0, 3)
            Random.seed!(47)
            fit_rate_equation(fp, stub; n_restarts = 1, lb = lb_x, ub = ub_x)
            @test all(lb_x .<= only(stub.seen).u0 .<= ub_x)
            start = only(stub.seen).u0
            fit = fit_at(start; lb = lb_x, ub = ub_x)
            # The start of k_E_S_to_ES, centred near 20, sits on `ub`.
            @test :k_E_S_to_ES in fit.at_bound
            @test fit.at_bound == Tuple(names[i] for i in 1:3 if start[i] == 1.0)
        end

        @testset "a point inside the box is not at a bound" begin
            @test fit_at(centre).at_bound == ()
        end

        @testset "relative mode reads the point before the scale" begin
            # The scale moves every fitted rate constant by a large factor, past the box.
            tight_lb, tight_ub = centre .- 1.0, centre .+ 1.0
            reference = [(S = 1e-6, P = 1e-9) => 1e6]
            fit = fit_at(centre; problem = FittingProblem(steady_state, table, 2.0;
                                                          reference),
                         lb = tight_lb, ub = tight_ub)
            @test fit.retcode === :Success
            @test any(k -> log(fit.params[k]) > tight_ub[findfirst(==(k), names)], names)
            @test fit.at_bound == ()
        end

        @testset "no finite restart gives ()" begin
            fit = fit_at(ub; objective = Inf)
            @test fit.retcode === :NoFiniteLoss
            @test fit.at_bound == ()
        end

        @testset "a failed reference scale still reports the best point" begin
            wrong_sign = FittingProblem(steady_state, table, 2.0;
                                        reference = [(P = 1e-5,) => 1.0])
            u = copy(centre)
            u[2] = ub[2]
            fit = fit_at(u; problem = wrong_sign)
            @test fit.retcode === :ReferenceFailed
            @test fit.at_bound == (names[2],)
        end

        @testset "destructures as params, loss, retcode" begin
            fitted, loss, retcode = fit_at(centre)
            @test fitted isa NamedTuple && loss == 1.0 && retcode === :Success
        end
    end
end

# ── The rate constants a fit scales ───────────────────────────────────────────────────

@testset "the fitted rate constants are what Fitting scales" begin
    # Every fitted steady-state rate constant of both conformations takes the factor;
    # scaling exactly these multiplies the rate by it.
    rng = Random.MersenneTwister(5)
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        names = _testhelper_steady_state_names(m)
        @test names isa Tuple{Vararg{Symbol}}
        @test allunique(names) && !isempty(names)
        @test issubset(names, parameters(m))
        @test isdisjoint(names, (:Keq, :L))
        params, Keq = _testhelper_random_independent_params(m; rng)
        met_names = metabolites(m)
        concs = NamedTuple{Tuple(met_names)}(Tuple(0.5 + rand(rng) for _ in met_names))
        rate = compile_rate_equation(m)
        v = rate(concs, params, Keq)
        scaled = EnzymeRates.Fitting._scale_rate_constants(m, params, 3.0)
        @test rate(concs, scaled, Keq) ≈ 3.0 * v rtol = 1e-9
    end
    # The exact fitted names: the set is the fitted steady-state constants of every
    # conformation's law (the inactive one's `I_` names included for MWC Dimer), and
    # rapid-equilibrium constants are absent (RE Uni-Uni).
    expected = Dict(
        "Uni-Uni" => (:k_ES_to_E_P, :k_ES_to_E_S, :k_E_S_to_ES),
        "RE Uni-Uni" => (:k_EA_to_E_P,),
        "MWC Dimer" => (:k_A_ES_to_EP, :k_I_ES_to_EP),
        "m_OnlyA_prod" => (:k_A_ES_to_EP,),
        "LDH i-state NonequalAI 5-group" =>
            (:k_A_E_NADH_to_ENADH, :k_A_ENADH_to_E_NADH,
             :k_A_ENADHPyruvate_to_ELactateNAD))
    for (name, names) in expected
        m = only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism
        @test Set(_testhelper_steady_state_names(m)) == Set(names)
    end
end

@testset "scaling the rate constants multiplies the steady-state ones only" begin
    for name in ("RE Uni-Uni", "MWC Dimer")
        m = only(s for s in MECHANISM_TEST_SPECS if s.name == name).mechanism
        ks = (parameters(m)...,)
        params = NamedTuple{ks}(Tuple(Float64.(eachindex(ks))))
        scaled = EnzymeRates.Fitting._scale_rate_constants(m, params, 10.0)
        ss = _testhelper_steady_state_names(m)
        @test keys(scaled) == ks
        @test !isempty(ss) && any(k -> !(k in ss), ks)
        @test all(scaled[k] == (k in ss ? 10 * params[k] : params[k]) for k in ks)
    end
end

# ── A 0/0 corner scored as a 0 ────────────────────────────────────────────────────────

@testset "a law the guard admits can be 0/0 at a mixed corner, scored as a 0 in a fit" begin
    # Census mechanism d7a09031d5dd2fc6 (ter-bi). At B = P = Q = 0 the rapid-equilibrium
    # segment {EBD, EBDP, EDP} has no form free of B and P, so the law is 0/0
    # there; the guard's faces never zero one substrate with every product (spec section
    # 8). At A = P = Q = 0, another such corner, the law is 0.
    am = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            (E + B ⇌ EB, EDP + B ⇌ EBDP, EP + B ⇌ EBP) :: OnlyA
            (E + P ⇌ EP, EB + P ⇌ EBP, EBD + P ⇌ EBDP,
             EQ + P ⇌ EPQ) :: EqualAI
            (E + Q <--> EQ, EP + Q <--> EPQ) :: EqualAI
            (EAB + D <--> EPQ, EB + D <--> EBD,
             EBP + D <--> EBDP, EP + D <--> EDP) :: OnlyA
            (EB + A <--> EAB, EBD + A <--> EPQ) :: OnlyA
            EQ + A ⇌ EAQ :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(am) ===
          EnzymeRates.Derivation._WELL_POSED
    params, Keq = _testhelper_random_independent_params(am; rng = Random.MersenneTwister(1))
    rate = compile_rate_equation(am)
    @test isnan(rate((A = 1.0, B = 0.0, D = 1.0, P = 0.0, Q = 0.0), params, Keq))
    @test rate((A = 0.0, B = 1.0, D = 1.0, P = 0.0, Q = 0.0), params, Keq) == 0.0
    # A fit scores the NaN as it scores the 0 a well-defined law predicts there: a zero rate
    # is refused, and a nonzero rate costs the same at the corner as at A = P = Q = 0, where
    # the law is 0.
    _testhelper_table(rate, A, B) =
        (group = [1, 1], Rate = [1.0, rate], A = [1.0, A], B = [1.0, B], D = [1.0, 1.0],
         P = [0.1, 0.0], Q = [0.1, 0.0])
    reference = _testhelper_row_reference(_testhelper_table(0.5, 1.0, 0.0))
    @test_throws(ErrorException("Zero rate at row 2: log(0) is undefined"),
                 FittingProblem(am, _testhelper_table(0.0, 1.0, 0.0), 1.0; reference))
    x = [log(params[p]) for p in parameters(am)]
    at_nan = FittingProblem(am, _testhelper_table(0.5, 1.0, 0.0), 1.0; reference)
    at_zero = FittingProblem(am, _testhelper_table(0.5, 0.0, 1.0), 1.0; reference)
    @test EnzymeRates.Fitting.loss!(x, at_nan) == EnzymeRates.Fitting.loss!(x, at_zero)
end
