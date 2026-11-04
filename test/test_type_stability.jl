# ABOUTME: Requires zero JET optimization reports (run-time dispatch, captured boxes) at
# ABOUTME: package code reached from every exported and public function, on Julia 1.13.
using Test
using EnzymeRates

if v"1.13.0-" <= VERSION < v"1.14.0-"
    using JET
    using OptimizationCMAEvolutionStrategy: CMAEvolutionStrategyOpt

    # A function name => why its site is allowed to dispatch at run time.
    const _TESTHELPER_DISPATCH_EXEMPTIONS = Pair{Symbol,String}[
        :_guard_message => "an error message",
        :_face_message => "an error message",
        :rate_equation => "a concrete mechanism is lifted to its type",
        :_assert_atom_conserving => "an error message",
        :_exc_string => "an error message: the text of a caught exception",
        :_seed_node_children => "a pmap boundary",
        :_fit_loss => "the solver barrier: a FunctionWrapper's call falls back on its " *
                      "callable of unknown type",
        :_optimization_problem => "the solver barrier: OptimizationFunction's constructor",
        :fit_rate_equation => "the solver barrier: solve, and the fit's values on the " *
                              "@nospecialize problem",
        :_reference_scale => "the solver barrier: rate_equation on the @nospecialize " *
                             "problem's mechanism",
        :FittingProblem => "compile_mechanism lifts a concrete mechanism and _rate_table " *
                           "the table's column names to the problem's type",
        :test_rate_equation => "compile_mechanism lifts a concrete mechanism and " *
                               "_rate_table the table's column names to the problem's type",
        :RateData => "the table's column names are lifted to its type",
        :_rate_table => "the table's column names are lifted to its type",
        :parse_mechanism => "Expr walking: it starts the walk on Meta.parse's result",
        :IdentifyRateEquationProblem => "the table's column names are lifted to its " *
                                        "RateData's type",
        :_row_params => "a row's parameter names are lifted to its type",
        :_pmap_to_vector => "a pmap boundary",
        :_write_rows_csv! => "a CSV boundary: the rows differ in type, so reading their " *
                             "fields dispatches",
        :_write_csv => "a CSV boundary",
        :_rows_to_dataframe => "a DataFrame boundary",
        :_single_row_dataframe => "a DataFrame boundary",
        :_select_best_row => "a DataFrame boundary: its columns' types are not known",
        :_progress => "the print to stdout, whose type is IO",
        :_default_save_dir => "Base's lock below Dates.format fails to optimize",
    ]

    # JET reports a function once per frame; the site is the innermost frame in package
    # code. A frame's linfo is a MethodInstance, or a CodeInstance holding one.
    _testhelper_frame_method(linfo::Core.MethodInstance) = linfo.def
    _testhelper_frame_method(linfo::Core.CodeInstance) =
        _testhelper_frame_method(linfo.def)
    _testhelper_frame_method(_) = nothing

    _testhelper_frame_compileable(linfo::Core.MethodInstance) =
        JET.is_compileable_mi(linfo)
    _testhelper_frame_compileable(linfo::Core.CodeInstance) =
        _testhelper_frame_compileable(linfo.def)
    _testhelper_frame_compileable(_) = true

    function _testhelper_in_package(frame)
        method = _testhelper_frame_method(frame.linfo)
        method isa Method || return false
        mod = method.module
        while true
            mod === EnzymeRates && return true
            parentmodule(mod) === mod && return false
            mod = parentmodule(mod)
        end
    end

    # `#_cv_model_selection#97` and `#_raw_step!##3` name `_cv_model_selection` and
    # `_raw_step!`.
    function _testhelper_undecorated_name(raw)
        raw = string(raw)
        stripped = replace(replace(raw, r"^#+" => ""), r"#+\d+$" => "")
        Symbol(isempty(stripped) ? raw : stripped)
    end

    function _testhelper_function_name(frame)
        method = _testhelper_frame_method(frame.linfo)
        _testhelper_undecorated_name(method isa Method ? method.name : frame.linfo)
    end

    # A row's `argtypes` for the call `f(argtypes...; kw...)` with keywords of type
    # `kwtype`. The row keeps `f` itself, so a check over the rows' functions finds it.
    _testhelper_keyword_call(kwtype, argtypes...) = (; kwtype, argtypes)

    _testhelper_report_opt(f, argtypes::Tuple) = JET.report_opt(f, argtypes)
    _testhelper_report_opt(f, call::NamedTuple) = JET.report_opt(
        Core.kwcall, (call.kwtype, f isa Type ? Type{f} : typeof(f), call.argtypes...))

    # `slowlock` is Base's path for a contended lock, not package code. Any package call that
    # takes a lock reaches it, such as `pkgversion`.
    function _testhelper_through_contended_lock(report)
        any(report.vst) do frame
            method = _testhelper_frame_method(frame.linfo)
            method isa Method && method.module === Base && method.name === :slowlock
        end
    end

    # The unique package sites at which `report_opt` finds a report on the row's call. A
    # report with no package frame counts its entry function, so none is dropped.
    function _testhelper_dispatch_sites(f, argtypes;
                                        exemptions = _TESTHELPER_DISPATCH_EXEMPTIONS)
        sites = NamedTuple[]
        for report in JET.get_reports(_testhelper_report_opt(f, argtypes))
            # A report below a call that cannot compile describes code that never runs at
            # those types: at run time the call dispatches on its arguments' concrete types.
            # JET skips such a frame but not a compileable callee of it. The entry always
            # counts, and the dispatch into such a call is reported at its caller.
            all(frame -> _testhelper_frame_compileable(frame.linfo), report.vst[2:end]) ||
                continue
            _testhelper_through_contended_lock(report) && continue
            i = findlast(_testhelper_in_package, report.vst)
            site = if i === nothing
                (; function_name = Symbol(nameof(f)), file = :none, line = 0)
            else
                frame = report.vst[i]
                (; function_name = _testhelper_function_name(frame),
                 file = frame.file, line = frame.line)
            end
            site in sites || push!(sites, site)
        end
        filter(s -> s.function_name ∉ first.(exemptions), sites)
    end

    const _testhelper_Mechs = EnzymeRates.Mechanisms
    const _testhelper_Deriv = EnzymeRates.Derivation
    const _testhelper_Enum = EnzymeRates.Enumeration
    const _testhelper_Mechanism = EnzymeRates.Mechanism
    const _testhelper_AllostericMechanism = EnzymeRates.AllostericMechanism

    const _testhelper_plain_singleton = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    const _testhelper_allosteric_singleton = @allosteric_mechanism begin
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
    const _testhelper_plain_concrete = _testhelper_Mechanism(_testhelper_plain_singleton)
    const _testhelper_allosteric_concrete =
        _testhelper_AllostericMechanism(_testhelper_allosteric_singleton)

    const _testhelper_PlainS = typeof(_testhelper_plain_singleton)
    const _testhelper_AlloS = typeof(_testhelper_allosteric_singleton)
    const _testhelper_PlainC = typeof(_testhelper_plain_concrete)
    const _testhelper_AlloC = typeof(_testhelper_allosteric_concrete)
    const _testhelper_Reaction = _testhelper_Mechs.EnzymeReaction

    # A NamedTuple type with the fields `names`, each a Float64.
    _testhelper_floats_type(names) =
        NamedTuple{Tuple(names), NTuple{length(names), Float64}}

    const _TESTHELPER_CONCRETE_FORMS = (("concrete", _testhelper_plain_concrete),
                                        ("concrete allosteric",
                                         _testhelper_allosteric_concrete))
    const _TESTHELPER_FORMS = (_TESTHELPER_CONCRETE_FORMS...,
                               ("singleton", _testhelper_plain_singleton),
                               ("singleton allosteric", _testhelper_allosteric_singleton))
    _testhelper_rate_equation_argtypes(m) = (typeof(m),
                                             _testhelper_floats_type(metabolites(m)),
                                             _testhelper_floats_type(parameters(m)))
    const _testhelper_Coefficients = NamedTuple{(:coefficients,), Tuple{Bool}}
    const _testhelper_CoefficientValues = NamedTuple{(:coefficient_values,), Tuple{Bool}}
    const _testhelper_BothCoefficientKeywords =
        NamedTuple{(:coefficients, :coefficient_values), Tuple{Bool, Bool}}

    const _testhelper_Fit = EnzymeRates.Fitting
    const _testhelper_KeqKeyword = NamedTuple{(:Keq,), Tuple{Float64}}
    _testhelper_rate_table(m) = (; group = ["g"], Rate = [1.0],
                                 (x => [1.0] for x in metabolites(m))...)
    _testhelper_problem_type(m) =
        typeof(_testhelper_Fit.FittingProblem(m, _testhelper_rate_table(m); Keq = 2.0))
    const _testhelper_PlainTable =
        typeof(_testhelper_rate_table(_testhelper_plain_singleton))
    const _testhelper_Optimizer = typeof(CMAEvolutionStrategyOpt())
    # A problem compiles a concrete mechanism, so it only ever holds a singleton.
    const _TESTHELPER_PROBLEM_FORMS = (("plain", _testhelper_plain_singleton),
                                       ("allosteric", _testhelper_allosteric_singleton))

    const _TESTHELPER_CHECKED_CALLS = Any[
        ("metabolites, concrete", metabolites, (_testhelper_PlainC,)),
        ("metabolites, concrete allosteric", metabolites, (_testhelper_AlloC,)),
        ("metabolites, singleton", metabolites, (_testhelper_PlainS,)),
        ("metabolites, singleton allosteric", metabolites, (_testhelper_AlloS,)),
        ("compile_mechanism", _testhelper_Mechs.compile_mechanism, (_testhelper_PlainC,)),
        ("compile_mechanism, allosteric", _testhelper_Mechs.compile_mechanism,
         (_testhelper_AlloC,)),
        ("catalytic", _testhelper_Mechs.catalytic, (_testhelper_PlainC,)),
        ("catalytic, allosteric", _testhelper_Mechs.catalytic, (_testhelper_AlloC,)),
        ("check_conformation_label", _testhelper_Mechs.check_conformation_label,
         (Symbol, String)),
        ("check_decodable_names", _testhelper_Mechs.check_decodable_names,
         (Vector{Symbol}, Vector{Symbol}, String)),
        ("excludes", _testhelper_Mechs.excludes, (Symbol, Symbol)),
        ("form_text", _testhelper_Mechs.form_text, (_testhelper_PlainC, Int)),
        ("step_text", _testhelper_Mechs.step_text, (_testhelper_PlainC, Int)),
        ("has_turnover_cycle", _testhelper_Mechs.has_turnover_cycle,
         (_testhelper_PlainC, BitVector)),
        ("has_unbalanced_cycle", _testhelper_Mechs.has_unbalanced_cycle,
         (Vector{Int}, Vector{Int}, Vector{Int}, BitVector, Vector{Int}, Vector{Int})),
        ("present_species", _testhelper_Mechs.present_species,
         (_testhelper_PlainC, BitVector)),
        ("rapid_equilibrium_segments", _testhelper_Mechs.rapid_equilibrium_segments,
         (_testhelper_PlainC, BitVector)),
        ("rapid_equilibrium_segments, with re",
         _testhelper_Mechs.rapid_equilibrium_segments,
         (_testhelper_PlainC, BitVector, BitVector)),
        ("segment_walk", _testhelper_Mechs.segment_walk,
         (BitVector, Vector{Vector{Int}}, Vector{Int}, Matrix{Int})),
        ("Mechanism constructor", _testhelper_Mechanism,
         (Vector{Symbol}, Vector{Symbol}, Vector{_testhelper_Mechs.EnzymeForm},
          Vector{_testhelper_Mechs.Step})),
        ("Mechanism from an EnzymeMechanism", _testhelper_Mechanism, (_testhelper_PlainS,)),
        ("AllostericMechanism from an AllostericEnzymeMechanism",
         _testhelper_AllostericMechanism, (_testhelper_AlloS,)),
        ("AllostericMechanism constructor from parts", _testhelper_AllostericMechanism,
         (Vector{Symbol}, Vector{Symbol}, Vector{_testhelper_Mechs.EnzymeForm},
          Vector{_testhelper_Mechs.Step}, Vector{Symbol}, Int,
          Vector{_testhelper_Mechs.RegulatorySite})),
        ("AllostericMechanism constructor", _testhelper_AllostericMechanism,
         (_testhelper_PlainC, Vector{Symbol}, Int,
          Vector{_testhelper_Mechs.RegulatorySite})),
        ("EnzymeMechanism from a Mechanism", EnzymeMechanism, (_testhelper_PlainC,)),
        ("AllostericEnzymeMechanism from an AllostericMechanism",
         AllostericEnzymeMechanism, (_testhelper_AlloC,)),
        (("rate_equation, $form", rate_equation, _testhelper_rate_equation_argtypes(m))
         for (form, m) in _TESTHELPER_FORMS)...,
        (("rate_equation_string, $form", rate_equation_string, (typeof(m),))
         for (form, m) in _TESTHELPER_FORMS)...,
        (("rate_equation_string with $kw, $form", rate_equation_string,
          _testhelper_keyword_call(Kw, typeof(m)))
         for (form, m) in _TESTHELPER_FORMS,
             (kw, Kw) in (("coefficients", _testhelper_Coefficients),
                          ("coefficient_values", _testhelper_CoefficientValues),
                          ("both", _testhelper_BothCoefficientKeywords)))...,
        (("parameters, $form", parameters, (typeof(m),))
         for (form, m) in _TESTHELPER_FORMS)...,
        (("parameters, $form, $(nameof(Kind))", parameters, (typeof(m), Kind))
         for (form, m) in _TESTHELPER_FORMS,
             Kind in (FittedParameters, IndependentParameters, AllParameters))...,
        (("dependent_parameter_values, $form", dependent_parameter_values,
          (typeof(m), _testhelper_floats_type(parameters(m))))
         for (form, m) in _TESTHELPER_FORMS)...,
        (("check_valid_rate_equation, $form", _testhelper_Deriv.check_valid_rate_equation,
          (typeof(m),)) for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        (("check_valid_rate_equation to Bool, $form",
          _testhelper_Deriv.check_valid_rate_equation, (Type{Bool}, typeof(m)))
         for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        (("fitted_parameter_count, $form", _testhelper_Deriv.fitted_parameter_count,
          (typeof(m),)) for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        (("name of a StepParameter, $form", _testhelper_Deriv.name,
          (_testhelper_Deriv.StepParameter, typeof(m)))
         for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        ("name of a RegulatorParameter", _testhelper_Deriv.name,
         (_testhelper_Deriv.RegulatorParameter, _testhelper_AlloC)),
        ("name of the AllostericParameter", _testhelper_Deriv.name,
         (_testhelper_Deriv.AllostericParameter, _testhelper_AlloC)),
        (("parameter_info, $form", _testhelper_Deriv.parameter_info, (typeof(m),))
         for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        (("rate_equation_term_count, $form", _testhelper_Deriv.rate_equation_term_count,
          (typeof(m),)) for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        ("rate_equation_term_count over a mask", _testhelper_Deriv.rate_equation_term_count,
         (_testhelper_PlainC, BitVector)),
        (("rate_equation_dedup_key, $form", _testhelper_Deriv.rate_equation_dedup_key,
          (typeof(m),)) for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
        ("init_mechanisms", _testhelper_Enum.init_mechanisms, (_testhelper_Reaction,)),
        ("seed_mechanisms", _testhelper_Enum.seed_mechanisms,
         (_testhelper_Reaction, Set{Symbol}, Set{Symbol})),
        ("starting_mechanisms", _testhelper_Enum.starting_mechanisms,
         (_testhelper_Reaction, Set{Symbol}, Set{Symbol})),
        (("expand_mechanisms on $(Mechs)", _testhelper_Enum.expand_mechanisms,
          (Mechs, _testhelper_Reaction))
         for Mechs in (Vector{_testhelper_PlainC}, Vector{_testhelper_AlloC},
                       Vector{Union{_testhelper_PlainC, _testhelper_AlloC}}))...,
        (("FittingProblem, $form", _testhelper_Fit.FittingProblem,
          _testhelper_keyword_call(_testhelper_KeqKeyword, typeof(m),
                                   typeof(_testhelper_rate_table(m))))
         for (form, m) in _TESTHELPER_FORMS)...,
        (("loss!, $form", _testhelper_Fit.loss!,
          (Vector{Float64}, _testhelper_problem_type(m)))
         for (form, m) in _TESTHELPER_PROBLEM_FORMS)...,
        (("fit_rate_equation, $form", fit_rate_equation,
          (_testhelper_problem_type(m), _testhelper_Optimizer))
         for (form, m) in _TESTHELPER_PROBLEM_FORMS)...,
        (("test_rate_equation, $form", _testhelper_Fit.test_rate_equation,
          _testhelper_keyword_call(_testhelper_KeqKeyword, typeof(m),
                                   _testhelper_floats_type(parameters(m)),
                                   typeof(_testhelper_rate_table(m))))
         for (form, m) in _TESTHELPER_FORMS)...,
        ("RateData", EnzymeRates.Data.RateData,
         _testhelper_keyword_call(_testhelper_KeqKeyword, _testhelper_PlainTable,
                                  Vector{Symbol})),
        ("RateData with a reference", EnzymeRates.Data.RateData,
         _testhelper_keyword_call(
             NamedTuple{(:Keq, :reference),
                        Tuple{Float64, Vector{Pair{@NamedTuple{S::Float64}, Float64}}}},
             _testhelper_PlainTable, Vector{Symbol})),
        ("parse_mechanism", EnzymeRates.DSL.parse_mechanism, (String,)),
        ("IdentifyRateEquationProblem", IdentifyRateEquationProblem,
         _testhelper_keyword_call(_testhelper_KeqKeyword, _testhelper_Reaction,
                                  _testhelper_PlainTable)),
        ("identify_rate_equation", identify_rate_equation,
         _testhelper_keyword_call(
             NamedTuple{(:optimizer,), Tuple{_testhelper_Optimizer}},
             IdentifyRateEquationProblem{_testhelper_PlainTable})),
        ("identify_rate_equation with save_dir", identify_rate_equation,
         _testhelper_keyword_call(
             NamedTuple{(:optimizer, :save_dir), Tuple{_testhelper_Optimizer, String}},
             IdentifyRateEquationProblem{_testhelper_PlainTable})),
        (("_pass1_record, $form", EnzymeRates._pass1_record, (typeof(m), Int, Int))
         for (form, m) in _TESTHELPER_CONCRETE_FORMS)...,
    ]

    @testset "type stability" begin
        for (name, f, argtypes) in _TESTHELPER_CHECKED_CALLS
            @testset "$name reports no run-time dispatch" begin
                sites = _testhelper_dispatch_sites(f, argtypes)
                isempty(sites) || foreach(s -> println(stderr, "  ", s), sites)
                @test isempty(sites)
            end
        end

        @testset "the harness finds a boxed variable and a call on a Ref{Any}" begin
            fixture = Module(:_testhelper_StabilityFixture)
            Core.eval(fixture, quote
                function ref_any_call()
                    r = Ref{Any}(1)
                    r[] + 1
                end
                function boxed_capture(n)
                    total = 0
                    add = k -> (total += k)
                    foreach(add, 1:n)
                    total = total + 0.5
                    total
                end
            end)
            for (f, argtypes) in ((fixture.ref_any_call, ()),
                                  (fixture.boxed_capture, (Int,)))
                @test !isempty(JET.get_reports(JET.report_opt(f, argtypes)))
                # None of their frames is in the package, so each counts its entry function.
                @test [s.function_name for s in _testhelper_dispatch_sites(f, argtypes)] ==
                      [nameof(f)]
            end
        end

        @testset "decorated names lose their decoration" begin
            for (raw, name) in (("#_cv_model_selection#97", :_cv_model_selection),
                                ("#_raw_step!##3", :_raw_step!),
                                ("#fit_rate_equation#8", :fit_rate_equation),
                                (:metabolites, :metabolites))
                @test _testhelper_undecorated_name(raw) == name
            end
        end

        @testset "a site in a package function is found, and exemptions remove it" begin
            # `_guard_message` builds an error message, so it dispatches for good.
            check = EnzymeRates.Derivation.check_valid_rate_equation
            none = Pair{Symbol,String}[]
            sites = _testhelper_dispatch_sites(check, (_testhelper_PlainC,);
                                               exemptions = none)
            guard = filter(s -> s.function_name == :_guard_message, sites)
            @test !isempty(guard)
            @test all(s -> endswith(string(s.file), "messages.jl"), guard)
            @test isempty(filter(s -> s.function_name == :_guard_message,
                                 _testhelper_dispatch_sites(check, (_testhelper_PlainC,);
                                     exemptions = [:_guard_message => "test"])))
        end

        @testset "a keyword row reaches the body of its function" begin
            # `rate_equation_string` runs the guard, whose message dispatches for good.
            sites = _testhelper_dispatch_sites(
                rate_equation_string,
                _testhelper_keyword_call(_testhelper_Coefficients, _testhelper_PlainC);
                exemptions = Pair{Symbol,String}[])
            @test any(s -> s.function_name == :_guard_message, sites)
        end

        # The exemption names `rate_equation`, whose singleton method is the hot path.
        @testset "rate_equation on a singleton needs no exemption" begin
            for m in (_testhelper_plain_singleton, _testhelper_allosteric_singleton)
                @test isempty(_testhelper_dispatch_sites(
                    rate_equation, _testhelper_rate_equation_argtypes(m);
                    exemptions = Pair{Symbol,String}[]))
            end
        end

        @testset "a report below a call that cannot compile counts for nothing" begin
            fixture = Module(:_testhelper_UncompileableFixture)
            Core.eval(fixture, quote
                struct _testhelper_Box{T}
                    x::T
                end
                _testhelper_dispatching_callee(r::Base.RefValue{Any}) = r[] + 1
                # Inlined, so the entry holds no dispatch into this frame of its own.
                @inline _testhelper_abstract_intermediate(b::_testhelper_Box) =
                    _testhelper_dispatching_callee(Ref{Any}(b.x))
                _testhelper_entry(boxes::Vector{_testhelper_Box}) =
                    _testhelper_abstract_intermediate(boxes[1])
            end)
            entry, argtypes = fixture._testhelper_entry, (Vector{fixture._testhelper_Box},)
            @test !isempty(JET.get_reports(JET.report_opt(entry, argtypes)))
            @test isempty(_testhelper_dispatch_sites(entry, argtypes))
            callee = fixture._testhelper_dispatching_callee
            @test !isempty(_testhelper_dispatch_sites(callee, (Base.RefValue{Any},)))
        end

        # A site here would run on every row of every evaluation of every fit.
        @testset "loss! needs no exemption" begin
            @test :loss! ∉ first.(_TESTHELPER_DISPATCH_EXEMPTIONS)
            for (_, m) in _TESTHELPER_PROBLEM_FORMS
                @test isempty(_testhelper_dispatch_sites(
                    _testhelper_Fit.loss!, (Vector{Float64}, _testhelper_problem_type(m));
                    exemptions = Pair{Symbol,String}[]))
            end
        end

        @testset "every public function has a row" begin
            checked = Set(f for (_, f, _) in _TESTHELPER_CHECKED_CALLS)
            for mod in (EnzymeRates, EnzymeRates.Mechanisms, EnzymeRates.Derivation,
                        EnzymeRates.Enumeration, EnzymeRates.Fitting, EnzymeRates.DSL,
                        EnzymeRates.Data), n in names(mod)
                startswith(string(n), "@") && continue  # macros run at parse time
                x = getfield(mod, n)
                x isa Function && x ∉ checked && @test "$mod.$n" == "has a row"
            end
        end

        @testset "every exemption is reported by some row" begin
            for (name, _) in _TESTHELPER_DISPATCH_EXEMPTIONS
                @test any(_TESTHELPER_CHECKED_CALLS) do (_, f, argtypes)
                    any(s -> s.function_name == name,
                        _testhelper_dispatch_sites(f, argtypes;
                                                   exemptions = Pair{Symbol,String}[]))
                end
            end
        end
    end
else
    @test_skip "zero JET reports is required on Julia 1.13 only"
end
