# ABOUTME: Runs JET's optimization analysis on one call and finds its dispatch sites.
# ABOUTME: Defines the row macros of test_type_stability.jl and tests this harness itself.
using Test
using EnzymeRates
using JET

# A frame's linfo is a MethodInstance, or a CodeInstance holding one.
_testhelper_frame_method(linfo::Core.MethodInstance) = linfo.def
_testhelper_frame_method(linfo::Core.CodeInstance) = _testhelper_frame_method(linfo.def)
_testhelper_frame_method(_) = nothing

_testhelper_frame_compileable(linfo::Core.MethodInstance) = JET.is_compileable_mi(linfo)
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

# `slowlock` is Base's path for a contended lock, not package code. Any package call that
# takes a lock reaches it, such as `pkgversion`.
function _testhelper_through_contended_lock(report)
    any(report.vst) do frame
        method = _testhelper_frame_method(frame.linfo)
        method isa Method && method.module === Base && method.name === :slowlock
    end
end

_testhelper_report_opt(f, argtypes, ::Nothing) = JET.report_opt(f, argtypes)
_testhelper_report_opt(f, argtypes, keywords::NamedTuple) =
    JET.report_opt(Core.kwcall, (typeof(keywords), Core.Typeof(f), argtypes...))

"""
The package sites at which JET reports on `f(argtypes...; keywords...)`, each once. A site's
function is the innermost frame in package code; JET gives a report in a shared callee once
per call, so an exempt function that reaches it first takes the site. A report with no
package frame counts its entry function, so none is dropped.
"""
function _testhelper_dispatch_sites(f, argtypes, keywords = nothing)
    sites = NamedTuple[]
    for report in JET.get_reports(_testhelper_report_opt(f, argtypes, keywords))
        # A report below a call that cannot compile mostly describes code that never runs at
        # those types: at run time the call dispatches on its arguments' concrete types. JET
        # skips such a frame but not a compileable callee of it. The exception is a frame
        # the compiler inlines: its compileable callee does run and dispatch, and its report
        # is dropped too. The entry always counts, and the dispatch into a call that is not
        # inlined is reported at its caller.
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
    sites
end

"""
What is wrong with a row's `sites` given the `exemption_sets` it lists, each a NamedTuple of
`reason` and `functions`: the sites no set allows, and the reasons of the sets no site uses.
"""
function _testhelper_row_problems(sites, exemption_sets)
    names = [s.function_name for s in sites]
    allowed = reduce(union, (set.functions for set in exemption_sets); init = Symbol[])
    (; unexpected_sites = filter(s -> s.function_name ∉ allowed, sites),
     unused_exemptions = [set.reason for set in exemption_sets
                          if isdisjoint(set.functions, names)])
end

# What the rows covered, for the checks that close test_type_stability.jl.
const _TESTHELPER_CHECKED_FUNCTIONS = Set{Any}()
const _TESTHELPER_FUNCTIONS_WITH_SITES = Set{Symbol}()

function _testhelper_check_row(label, f, argtypes, keywords, exemption_sets)
    push!(_TESTHELPER_CHECKED_FUNCTIONS, f)
    sites = _testhelper_dispatch_sites(f, argtypes, keywords)
    union!(_TESTHELPER_FUNCTIONS_WITH_SITES, (s.function_name for s in sites))
    (; unexpected_sites, unused_exemptions) =
        _testhelper_row_problems(sites, exemption_sets)
    @testset "$label" begin
        @test isempty(unexpected_sites)
        @test isempty(unused_exemptions)
    end
end

# An argument written `::T` stands for a value of type `T`, which can be abstract. Any other
# argument is evaluated for its type alone. The call itself never runs.
function _testhelper_row(call, exemption_sets)
    Meta.isexpr(call, :call) || error("a row is a call, got $call")
    f, args = call.args[1], call.args[2:end]
    any(a -> Meta.isexpr(a, :kw), args) && error("write the keywords after `;` in $call")
    keyword_blocks = filter(a -> Meta.isexpr(a, :parameters), args)
    positional = filter(a -> !Meta.isexpr(a, :parameters), args)
    argtype(a) = Meta.isexpr(a, :(::), 1) ? a.args[1] : :(Core.Typeof($a))
    argtypes = Expr(:tuple, map(argtype, positional)...)
    keywords = isempty(keyword_blocks) ? nothing : Expr(:tuple, only(keyword_blocks))
    esc(:(_testhelper_check_row($(string(call)), $f, $argtypes, $keywords,
                                $(Expr(:tuple, exemption_sets...)))))
end

"""
    @_testhelper_no_dispatch f(x, ::T; k = v)

Requires that JET reports nothing on the call in package code. The call does not run.
"""
macro _testhelper_no_dispatch(call)
    _testhelper_row(call, ())
end

"""
    @_testhelper_dispatches_in set1 set2 f(x, ::T; k = v)

Requires that JET reports on the call only inside functions of the listed exemption sets,
and at least once inside each set. The call does not run.
"""
macro _testhelper_dispatches_in(sets_and_call...)
    _testhelper_row(sets_and_call[end], sets_and_call[1:end-1])
end

@testset "the type-stability harness" begin
    plain_concrete = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    check = EnzymeRates.Derivation.check_valid_rate_equation

    @testset "it finds a boxed variable and a call on a Ref{Any}" begin
        fixture = Module(:_testhelper_StabilityFixture)
        Core.eval(fixture, quote
            function _testhelper_ref_any_call()
                r = Ref{Any}(1)
                r[] + 1
            end
            function _testhelper_boxed_capture(n)
                total = 0
                add = k -> (total += k)
                foreach(add, 1:n)
                total = total + 0.5
                total
            end
        end)
        for (f, argtypes) in ((fixture._testhelper_ref_any_call, ()),
                              (fixture._testhelper_boxed_capture, (Int,)))
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

    @testset "a site is found in the package function that dispatches" begin
        # `_guard_message` builds an error message, so it dispatches for good.
        guard = filter(s -> s.function_name == :_guard_message,
                       _testhelper_dispatch_sites(check, (typeof(plain_concrete),)))
        @test !isempty(guard)
        @test all(s -> endswith(string(s.file), "messages.jl"), guard)
    end

    @testset "a row's sets must allow every site, and each set must be used" begin
        sites = [(; function_name = :a, file = :none, line = 0)]
        allows_a = (; reason = "allows a", functions = [:a])
        allows_b = (; reason = "allows b", functions = [:b])
        @test _testhelper_row_problems(sites, (allows_a,)) ==
              (; unexpected_sites = NamedTuple[], unused_exemptions = String[])
        @test _testhelper_row_problems(sites, ()).unexpected_sites == sites
        @test _testhelper_row_problems(sites, (allows_a, allows_b)).unused_exemptions ==
              ["allows b"]
    end

    @testset "a keyword call reaches the body of its function" begin
        # `rate_equation_string` runs the guard, whose message dispatches for good.
        sites = _testhelper_dispatch_sites(rate_equation_string, (typeof(plain_concrete),),
                                           (; coefficients = true))
        @test any(s -> s.function_name == :_guard_message, sites)
    end

    @testset "a report below a call that cannot compile is dropped, even inlined" begin
        fixture = Module(:_testhelper_UncompileableFixture)
        Core.eval(fixture, quote
            struct _testhelper_Box{T}
                x::T
            end
            _testhelper_dispatching_callee(r::Base.RefValue{Any}) = r[] + 1
            # Inlined, so the entry holds no dispatch into this frame of its own and the
            # callee's dispatch does run; the filter drops its report all the same.
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
end
