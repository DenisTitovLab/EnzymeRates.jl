# ABOUTME: One row per checked call: those JET finds no run-time dispatch in, and those that
# ABOUTME: dispatch only for a stated reason. Covers every exported and public function.
using Test
using EnzymeRates
using OptimizationCMAEvolutionStrategy: CMAEvolutionStrategyOpt
include("type_stability_harness.jl")

@testset "type stability" begin
    Mech = EnzymeRates.Mechanisms
    Deriv = EnzymeRates.Derivation
    Enum = EnzymeRates.Enumeration
    Fit = EnzymeRates.Fitting

    # ── Why a row may dispatch ────────────────────────────────────────────────
    # A site is exempt when the function holding it is in a set its row lists. None of these
    # runs per data row: a site in `loss!` would run on every row of every evaluation.
    error_messages = (;
        reason = "builds an error message, so it runs only once something has failed",
        functions = [:_guard_message, :_face_message, :_atom_conserving_message,
                     :_exc_string])
    value_to_type = (;
        reason = "lifts a run-time value (a mechanism's steps, a table's or a row's " *
                 "names) to a type, which is what makes the code after it concrete",
        functions = [:rate_equation, :FittingProblem, :test_rate_equation, :RateData,
                     :_rate_table, :IdentifyRateEquationProblem, :_row_params,
                     :_fit_keywords])
    # `fit_rate_equation` hands the solver one problem type for every mechanism: a
    # FunctionWrapper around `loss!`, and a `@nospecialize` problem. `_fit_loss` pays one
    # dispatch per evaluation, not per row.
    solver_boundary = (;
        reason = "the solver's boundary, which hides the mechanism's type from the solver",
        functions = [:_fit_loss, :_optimization_problem, :fit_rate_equation,
                     :_concentration_scales, :_reference_scale, :FitObjective])
    pmap_and_tables = (;
        reason = "reads pmap's Vector{Any}, or table columns of unknown type",
        functions = [:_seed_node_children, :_pmap_to_vector, :_write_rows_csv!, :_write_csv,
                     :_rows_to_dataframe, :_single_row_dataframe, :_select_best_row])
    expr_walk = (;
        reason = "walks an Expr, whose args are a Vector{Any}",
        functions = [:parse_mechanism])
    prints_to_stdout = (;
        reason = "prints to stdout, whose type is the abstract IO",
        functions = [:_progress])
    checks_user_input = (;
        reason = "validates a user-supplied anchor of unknown type",
        functions = [:_check_anchors])
    exemption_sets = (error_messages, value_to_type, solver_boundary, pmap_and_tables,
                      expr_walk, prints_to_stdout, checks_user_input)

    # ── The values whose types the rows use; no row runs its call ────────────
    plain = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    allosteric = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): X::OnlyA
    end
    plain_concrete = EnzymeRates.Mechanism(plain)
    allosteric_concrete = EnzymeRates.AllostericMechanism(allosteric)
    Plain = typeof(plain_concrete)
    Allosteric = typeof(allosteric_concrete)
    Reaction = Mech.EnzymeReaction
    floats(names) = NamedTuple{Tuple(names)}(ntuple(_ -> 1.0, length(names)))
    table(m) = (; group = ["g"], Rate = [1.0], (x => [1.0] for x in metabolites(m))...)
    problem(m) = Fit.FittingProblem(m, table(m); Keq = 2.0)
    optimizer = CMAEvolutionStrategyOpt()
    plain_table = table(plain)

    concrete_forms = (; plain_concrete, allosteric_concrete)
    singleton_forms = (; plain, allosteric)
    every_form = (; concrete_forms..., singleton_forms...)

    @testset "Mechanisms" begin
        @_testhelper_no_dispatch Mech.check_conformation_label(:A, "label")
        @_testhelper_no_dispatch Mech.check_decodable_names(
            ::Vector{Symbol}, ::Vector{Symbol}, "label")
        @_testhelper_no_dispatch Mech.excludes(:A, :B)
        @_testhelper_no_dispatch Mech.form_text(plain_concrete, 1)
        @_testhelper_no_dispatch Mech.step_text(plain_concrete, 1)
        @_testhelper_no_dispatch Mech.has_turnover_cycle(plain_concrete, ::BitVector)
        @_testhelper_no_dispatch Mech.has_unbalanced_cycle(
            ::Vector{Int}, ::Vector{Int}, ::Vector{Int}, ::BitVector, ::Vector{Int},
            ::Vector{Int})
        @_testhelper_no_dispatch Mech.present_species(plain_concrete, ::BitVector)
        @_testhelper_no_dispatch Mech.rapid_equilibrium_segments(
            plain_concrete, ::BitVector)
        @_testhelper_no_dispatch Mech.rapid_equilibrium_segments(
            plain_concrete, ::BitVector, ::BitVector)
        @_testhelper_no_dispatch Mech.segment_walk(
            ::BitVector, ::Vector{Vector{Int}}, ::Vector{Int}, ::Matrix{Int})

        @_testhelper_no_dispatch EnzymeRates.Mechanism(
            ::Vector{Symbol}, ::Vector{Symbol}, ::Vector{Mech.EnzymeForm},
            ::Vector{Mech.Step})
        @_testhelper_no_dispatch EnzymeRates.Mechanism(plain)
        @_testhelper_no_dispatch EnzymeRates.AllostericMechanism(allosteric)
        @_testhelper_no_dispatch EnzymeRates.AllostericMechanism(
            ::Vector{Symbol}, ::Vector{Symbol}, ::Vector{Mech.EnzymeForm},
            ::Vector{Mech.Step}, ::Vector{Symbol}, 1, ::Vector{Mech.RegulatorySite})
        @_testhelper_no_dispatch EnzymeRates.AllostericMechanism(
            plain_concrete, ::Vector{Symbol}, 1, ::Vector{Mech.RegulatorySite})
        @_testhelper_no_dispatch EnzymeMechanism(plain_concrete)
        @_testhelper_no_dispatch AllostericEnzymeMechanism(allosteric_concrete)

        @testset "$form" for (form, m) in pairs(every_form)
            @_testhelper_no_dispatch metabolites(m)
        end
        @testset "$form" for (form, m) in pairs(concrete_forms)
            @_testhelper_no_dispatch Mech.compile_mechanism(m)
            @_testhelper_no_dispatch Mech.catalytic(m)
        end
    end

    @testset "Derivation" begin
        @testset "$form" for (form, m) in pairs(singleton_forms)
            # The hot path: each data row of each evaluation calls `rate_equation` here.
            @_testhelper_no_dispatch rate_equation(
                m, floats(metabolites(m)), floats(parameters(m)))
            @_testhelper_no_dispatch parameters(m)
            @_testhelper_no_dispatch parameters(m, FittedParameters())
            @_testhelper_no_dispatch parameters(m, IndependentParameters())
            @_testhelper_no_dispatch parameters(m, AllParameters())
        end
        @testset "$form" for (form, m) in pairs(concrete_forms)
            @_testhelper_dispatches_in value_to_type rate_equation(
                m, floats(metabolites(m)), floats(parameters(m)))
            @_testhelper_dispatches_in error_messages parameters(m)
            @_testhelper_dispatches_in error_messages parameters(m, FittedParameters())
            @_testhelper_dispatches_in error_messages parameters(m, IndependentParameters())
            @_testhelper_dispatches_in error_messages parameters(m, AllParameters())
            @_testhelper_dispatches_in error_messages Deriv.check_valid_rate_equation(m)
            @_testhelper_no_dispatch Deriv.check_valid_rate_equation(Bool, m)
            @_testhelper_no_dispatch Deriv.fitted_parameter_count(m)
            @_testhelper_no_dispatch Deriv.name(::Deriv.StepParameter, m)
            @_testhelper_dispatches_in error_messages Deriv.parameter_info(m)
            @_testhelper_no_dispatch Deriv.rate_equation_term_count(m)
            @_testhelper_dispatches_in error_messages Deriv.rate_equation_dedup_key(m)
        end
        @testset "$form" for (form, m) in pairs(every_form)
            @_testhelper_dispatches_in error_messages rate_equation_string(m)
            @_testhelper_dispatches_in error_messages rate_equation_string(
                m; coefficients = true)
            @_testhelper_dispatches_in error_messages rate_equation_string(
                m; coefficient_values = true)
            @_testhelper_dispatches_in error_messages rate_equation_string(
                m; coefficients = true, coefficient_values = true)
            @_testhelper_dispatches_in error_messages dependent_parameter_values(
                m, floats(parameters(m)))
        end
        @_testhelper_no_dispatch Deriv.name(::Deriv.RegulatorParameter, allosteric_concrete)
        @_testhelper_no_dispatch Deriv.name(
            ::Deriv.AllostericParameter, allosteric_concrete)
        @_testhelper_no_dispatch Deriv.rate_equation_term_count(plain_concrete, ::BitVector)
    end

    @testset "Enumeration" begin
        Names = Set{Symbol}
        @_testhelper_dispatches_in error_messages Enum.init_mechanisms(::Reaction)
        @_testhelper_dispatches_in error_messages pmap_and_tables Enum.seed_mechanisms(
            ::Reaction, ::Names, ::Names)
        @_testhelper_dispatches_in error_messages pmap_and_tables Enum.starting_mechanisms(
            ::Reaction, ::Names, ::Names)
        for Mechs in (Plain, Allosteric, Union{Plain, Allosteric})
            @_testhelper_dispatches_in error_messages Enum.expand_mechanisms(
                ::Vector{Mechs}, ::Reaction)
            @_testhelper_no_dispatch Enum._is_seed_node(
                ::Mechs, ::Reaction, ::Names, ::Names)
        end
    end

    @testset "Fitting" begin
        @testset "$form" for (form, m) in pairs(singleton_forms)
            # A site here would run on every data row of every evaluation of every fit.
            @_testhelper_no_dispatch Fit.loss!(::Vector{Float64}, problem(m))
            @_testhelper_dispatches_in error_messages solver_boundary fit_rate_equation(
                problem(m), optimizer)
            # `_scale_rate_constants` takes a `@nospecialize` mechanism, so it compiles for
            # `Any`, not for a concrete type.
            @_testhelper_dispatches_in error_messages Fit._scale_rate_constants(
                ::Any, floats(parameters(m)), 1.0)
        end
        # `fit_rate_equation` takes a `@nospecialize` problem, so this is what compiles.
        @_testhelper_dispatches_in error_messages solver_boundary fit_rate_equation(
            ::Fit.FittingProblem, optimizer)
        @testset "$form" for (form, m) in pairs(every_form)
            @_testhelper_dispatches_in value_to_type Fit.FittingProblem(
                m, table(m); Keq = 2.0)
            @_testhelper_dispatches_in value_to_type Fit.test_rate_equation(
                m, floats(parameters(m)), table(m); Keq = 2.0)
        end
    end

    @testset "Data and DSL" begin
        @_testhelper_dispatches_in value_to_type EnzymeRates.Data.RateData(
            plain_table, ::Vector{Symbol}; Keq = 2.0)
        @_testhelper_dispatches_in value_to_type EnzymeRates.Data.RateData(
            plain_table, ::Vector{Symbol}; Keq = 2.0, reference = [(S = 1.0,) => 1.0])
        @_testhelper_dispatches_in error_messages expr_walk EnzymeRates.DSL.parse_mechanism(
            "the mechanism's text")
    end

    @testset "identify_rate_equation and its workers" begin
        anchors = Pair{NamedTuple, Float64}[]
        # The keywords `_cv_model_selection`'s workers pass, with and without a reference.
        solver_options = (; optimizer, n_restarts = 1, maxtime = 1.0, maxiters = 1,
                          abstol = nothing, reltol = nothing, callback = nothing,
                          solver_kwargs = (;))
        without_reference = (; Keq = 2.0, absolute = false, solver_options...)
        with_reference = (; Keq = 2.0, absolute = false, reference = anchors,
                          solver_options...)

        @_testhelper_dispatches_in value_to_type IdentifyRateEquationProblem(
            ::Reaction, plain_table; Keq = 2.0)
        @_testhelper_no_dispatch EnzymeRates._fit_order(::EnzymeRates._Candidate)
        # The print, the tables and the CSV are what a search's master does between sweeps.
        @_testhelper_dispatches_in(
            error_messages, value_to_type, pmap_and_tables, prints_to_stdout,
            identify_rate_equation(
                ::IdentifyRateEquationProblem{typeof(plain_table)}; optimizer))
        @_testhelper_dispatches_in(
            error_messages, value_to_type, pmap_and_tables, prints_to_stdout,
            identify_rate_equation(
                ::IdentifyRateEquationProblem{typeof(plain_table)}; optimizer,
                save_dir = "results"))
        @testset "$form" for (form, m) in pairs(concrete_forms)
            @_testhelper_dispatches_in error_messages EnzymeRates._pass1_record(m, 1, 1)
        end
        @testset "$form" for (form, m) in pairs(singleton_forms)
            @_testhelper_dispatches_in(
                error_messages, value_to_type, solver_boundary,
                EnzymeRates._cross_validation_loss(
                    m, table(m), "g"; without_reference...))
            @_testhelper_dispatches_in(
                error_messages, value_to_type, solver_boundary, checks_user_input,
                EnzymeRates._cross_validation_loss(m, table(m), "g"; with_reference...))
            @_testhelper_dispatches_in value_to_type EnzymeRates._fit_keywords(
                m, ::IdentifyRateEquationProblem{typeof(table(m))})
            @_testhelper_dispatches_in(
                value_to_type, checks_user_input,
                EnzymeRates._fold_problem(m, table(m), ::BitVector; Keq = 2.0,
                                          absolute = false, reference = anchors))
        end
    end

    # ── What the rows must cover ──────────────────────────────────────────────
    # Only functions are required to have a row. A public type (`FittingProblem`,
    # `RateData`, `Mechanism`) is not a `Function`, so deleting its row would not fail.
    @testset "every public function has a row" begin
        for mod in (EnzymeRates, Mech, Deriv, Enum, Fit, EnzymeRates.DSL, EnzymeRates.Data),
                n in names(mod)
            startswith(string(n), "@") && continue  # macros run at parse time
            x = getfield(mod, n)
            x isa Function && x ∉ _TESTHELPER_CHECKED_FUNCTIONS &&
                @test "$mod.$n" == "has a row"
        end
    end

    @testset "every exempt function is reported by some row" begin
        exempt = reduce(union, (set.functions for set in exemption_sets))
        @test isempty(setdiff(exempt, _TESTHELPER_FUNCTIONS_WITH_SITES))
    end

    @testset "no exemption can cover loss!" begin
        @test all(set -> :loss! ∉ set.functions, exemption_sets)
    end
end
