# ABOUTME: Checks src/ by parsing it: public names, layers, underscore names, exports, file
# ABOUTME: names, one-way calls between files, ABOUTME lines, test file names and the docs.

import REPL  # defines Base.Docs.doc(::Binding), read by the docstring testset

include("module_map.jl")

"""The tables of a Markdown text, as `(; label, header, rows)`: `label` is the last
non-blank line above the table, `header` the cells of its first row and `rows` the cells
of each row below the separator."""
function _testhelper_doc_tables(text)
    cells(line) = String.(strip.(split(strip(line), '|')[2:(end - 1)]))
    tables, label, rows = [], "", Vector{String}[]
    for line in [split(text, '\n'); ""]
        if startswith(strip(line), "|")
            push!(rows, cells(line))
        else
            isempty(rows) ||
                push!(tables, (; label, header = rows[1], rows = rows[3:end]))
            rows = Vector{String}[]
            isempty(strip(line)) || (label = strip(line))
        end
    end
    tables
end

"""The names a table cell writes in backticks."""
_testhelper_cell_names(cell) = Set(Symbol(m[1]) for m in eachmatch(r"`([^`]+)`", cell))

"""The module a table cell or label names; the top level is `:EnzymeRates`."""
function _testhelper_doc_module(text)
    text = replace(text, r"[*`]" => "")
    text in ("top level", "The top level") ? :EnzymeRates : Symbol(text)
end

"""How the docs write a module."""
_testhelper_doc_label(m::Symbol) = m === :EnzymeRates ? "top level" : string(m)

"""One line per key where the documented names `doc` differ from the names `src` finds,
naming the key and the names the docs lack or add."""
function _testhelper_name_diffs(what, doc::Dict, src::Dict)
    show_names(s) = join(["`$n`" for n in sort!(_testhelper_doc_label.(collect(s)))], ", ")
    lines = String[]
    for k in sort!(collect(union(keys(doc), keys(src))))
        d, s = get(doc, k, Set{Symbol}()), get(src, k, Set{Symbol}())
        parts = String[]
        isempty(setdiff(s, d)) || push!(parts, "missing " * show_names(setdiff(s, d)))
        isempty(setdiff(d, s)) || push!(parts, "extra " * show_names(setdiff(d, s)))
        isempty(parts) ||
            push!(lines, "$what $(_testhelper_doc_label(k)): " * join(parts, "; "))
    end
    lines
end

"""Where the module map of the Markdown text `doc` disagrees with the tree at `root`: the
interface tables (names, "Used by" and ✓ marks) against each module's public names, the
cross-module uses and the exports."""
function _testhelper_module_map_problems(doc, root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    exports = Set(_testhelper_exports(joinpath(root, "EnzymeRates.jl")))
    owner(n) = something(findfirst(api -> n in api, apis), :EnzymeRates)
    apis[:EnzymeRates] = Set(n for n in exports if owner(n) === :EnzymeRates)
    names, ticks, used_by = (Dict{Any, Set{Symbol}}() for _ in 1:3)
    for t in _testhelper_doc_tables(doc)
        if t.header == ["Name", "What it does", "Used by", "Exported"]
            m, name = _testhelper_doc_module(t.label), nothing
            for r in t.rows
                # A row with an empty Name cell continues the row above.
                if !isempty(r[1])
                    name = only(_testhelper_cell_names(r[1]))
                    push!(get!(names, m, Set{Symbol}()), name)
                    r[4] == "✓" && push!(get!(ticks, m, Set{Symbol}()), name)
                end
                union!(get!(used_by, name, Set{Symbol}()),
                    (_testhelper_doc_module(strip(u)) for u in split(r[3], ',')
                        if !isempty(strip(u))))
            end
        end
    end
    users = Dict{Any, Set{Symbol}}()
    for u in _testhelper_uses(root)
        push!(get!(users, u.name, Set{Symbol}()), u.user)
    end
    exported_api = Dict(m => intersect(api, exports) for (m, api) in apis)
    [_testhelper_name_diffs("interface table", names, apis);
     _testhelper_name_diffs("✓ in interface table", ticks, exported_api);
     _testhelper_name_diffs("Used by, row", used_by, users)]
end

"""Runs the four rules over the tree at `root`."""
function _testhelper_check_boundaries(root)
    modules = sort!(collect(keys(_testhelper_module_sources(root))))
    @testset "every cross-module name is public in its module" begin
        @test isempty(_testhelper_api_violations(root))
    end
    @testset "modules respect the layers" begin
        @test isempty(_testhelper_layer_violations(root))
    end
    @testset "no cross-module underscore names" begin
        underscore = _testhelper_underscore_violations(root)
        for m in modules
            @test isempty(underscore[m])
        end
    end
    @testset "no module exports a name" begin
        @test isempty(_testhelper_export_violations(root))
    end
end

function _testhelper_with_tree(f, files::Pair...)
    mktempdir() do root
        for (path, text) in files
            mkpath(dirname(joinpath(root, path)))
            write(joinpath(root, path), text)
        end
        f(root)
    end
end

@testset "boundary scanner" begin
    @testset "cross_refs" begin
        @test _testhelper_cross_refs("Derivation.check_valid_rate_equation(m)") ==
              [(:Derivation, :check_valid_rate_equation, 1)]
        @test isempty(_testhelper_cross_refs("Base.show(io, x)"))
        @test _testhelper_cross_refs("\n\nf(x::Mechanisms.Mechanism) = 1") ==
              [(:Mechanisms, :Mechanism, 3)]
        @test _testhelper_cross_refs("function Mechanisms.show(io, m)\nend") ==
              [(:Mechanisms, :show, 1)]
        @test isempty(_testhelper_cross_refs("using ..Mechanisms: Mechanisms"))
        @test isempty(_testhelper_cross_refs("using LinearAlgebra: dot"))
        @test isempty(_testhelper_cross_refs("using ..Mechanisms"))
        @test _testhelper_cross_refs("\nimport ..Derivation: _x") ==
              [(:Derivation, :_x, 2)]
        @test _testhelper_cross_refs("using ..Derivation: foo, Derivation, bar as b") ==
              [(:Derivation, :foo, 1), (:Derivation, :bar, 1)]
        @test _testhelper_cross_refs("EnzymeRates.Mechanisms._hidden(m)") ==
              [(:Mechanisms, :_hidden, 1)]
        @test isempty(_testhelper_cross_refs("EnzymeRates.fit(m)"))
        @test isempty(_testhelper_cross_refs(":(EnzymeRates.Mechanisms.\$T(x))"))
        @test _testhelper_cross_refs("using .DSL: DSL, @enzyme_mechanism, parse") ==
              [(:DSL, Symbol("@enzyme_mechanism"), 1), (:DSL, :parse, 1)]
        @test _testhelper_cross_refs("DSL.@enzyme_mechanism begin end") ==
              [(:DSL, Symbol("@enzyme_mechanism"), 1)]
        @test Set(_testhelper_cross_refs("a = 1\nFitting.fit(\"\$(Derivation.f(x))\")")) ==
              Set([(:Derivation, :f, 2), (:Fitting, :fit, 2)])
    end

    cleanly_layered = (
        "Mechanisms/Mechanisms.jl" => "public Mechanism, compile_mechanism",
        "Mechanisms/types.jl" => "struct Mechanism\n    edges::Vector{Int}\nend",
        "Derivation/Derivation.jl" => "public check_valid_rate_equation",
        "Derivation/guard.jl" =>
            "f(m::Mechanisms.Mechanism) = Mechanisms.compile_mechanism(m)",
        "Fitting/Fitting.jl" => "public fit",
        "Fitting/fit.jl" => "g(m) = Derivation.check_valid_rate_equation(m)",
        "EnzymeRates.jl" => "h(m) = Fitting.fit(m)",
    )
    @testset "module_sources and api" begin
        _testhelper_with_tree(cleanly_layered...) do root
            @test _testhelper_module_sources(root) == Dict(
                :Mechanisms => [joinpath(root, "Mechanisms", "Mechanisms.jl"),
                    joinpath(root, "Mechanisms", "types.jl")],
                :Derivation => [joinpath(root, "Derivation", "Derivation.jl"),
                    joinpath(root, "Derivation", "guard.jl")],
                :Fitting => [joinpath(root, "Fitting", "Fitting.jl"),
                    joinpath(root, "Fitting", "fit.jl")])
            @test _testhelper_api(root, :Mechanisms) ==
                  Set([:Mechanism, :compile_mechanism])
            @test _testhelper_api(root, :Fitting) == Set([:fit])
        end
        _testhelper_with_tree(
            "DSL/DSL.jl" => "public parse, @enzyme_mechanism") do root
            @test _testhelper_api(root, :DSL) == Set([:parse, Symbol("@enzyme_mechanism")])
        end
        _testhelper_with_tree(
            "Fitting/Fitting.jl" => "public fit\npublic refit") do root
            @test _testhelper_api(root, :Fitting) == Set([:fit, :refit])
        end
        _testhelper_with_tree("Fitting/Fitting.jl" => "public fit",
            "Fitting/extra.jl" => "public refit") do root
            @test _testhelper_api(root, :Fitting) == Set([:fit, :refit])
        end
        _testhelper_with_tree("Fitting/Fitting.jl" => "fit() = 1") do root
            @test_throws "no public statement" _testhelper_api(root, :Fitting)
        end
    end
    @testset "a clean tree passes every rule" begin
        _testhelper_with_tree(cleanly_layered...) do root
            @test isempty(_testhelper_api_violations(root))
            @test isempty(_testhelper_layer_violations(root))
            @test all(isempty, values(_testhelper_underscore_violations(root)))
            @test isempty(_testhelper_export_violations(root))
        end
    end
    @testset "rule 1 finds a name that is not public, in a module and at the top" begin
        _testhelper_with_tree(cleanly_layered...,
            "Fitting/bad.jl" => "g(m) = Derivation.secret(m)",
            "identify.jl" => "h(m) = Mechanisms.hidden(m)\nMechanisms.Mechanism") do root
            found = _testhelper_api_violations(root)
            @test Set((u.user, u.file, u.mod, u.name, u.line) for u in found) == Set([
                (:Fitting, "Fitting/bad.jl", :Derivation, :secret, 1),
                (:EnzymeRates, "identify.jl", :Mechanisms, :hidden, 1)])
        end
    end
    @testset "rule 2 finds a sibling the layer forbids, not the top level" begin
        _testhelper_with_tree(cleanly_layered...,
            "Mechanisms/up.jl" => "f(m) = Derivation.check_valid_rate_equation(m)",
            "Derivation/across.jl" => "f(m) = Fitting.fit(m)",
            "Fitting/ok.jl" => "f(m) = Mechanisms.compile_mechanism(m)") do root
            found = _testhelper_layer_violations(root)
            @test Set((u.user, u.mod) for u in found) ==
                  Set([(:Mechanisms, :Derivation), (:Derivation, :Fitting)])
        end
    end
    @testset "rule 3 finds underscore names public or used across modules" begin
        _testhelper_with_tree(
            "Mechanisms/Mechanisms.jl" => "public Mechanism, _raw",
            "Derivation/Derivation.jl" => "public _guard, ok, _guard",
            "Fitting/f.jl" => "f() = Derivation._guard()",
            "Fitting/Fitting.jl" => "public fit") do root
            @test _testhelper_underscore_violations(root) == Dict(
                :Mechanisms => [:_raw], :Derivation => [:_guard], :Fitting => [])
        end
    end
    @testset "rule 4 finds an export in a module, not at the top level" begin
        _testhelper_with_tree(cleanly_layered...,
            "Derivation/api.jl" => "f() = 1\n\nexport f, g",
            "EnzymeRates.jl" => "export fit") do root
            @test _testhelper_export_violations(root) == [("Derivation/api.jl", 3)]
        end
    end
end

@testset "file scanner" begin
    problem_type = """
        # ABOUTME: Problem and _size, the data a
        # ABOUTME: fit reads.
        abstract type AbstractProblem end
        struct Problem <: AbstractProblem
            data::Vector{Float64}
        end
        @enum _Kind _PLAIN _ALLOSTERIC
        \"\"\"The size of `p`.\"\"\"
        _size(p::Problem) = length(p.data)
        """
    file_tree = (
        "EnzymeRates.jl" => """
            \"\"\"Does it all.\"\"\"
            module EnzymeRates
            export fit, identify
            include("Fitting/Fitting.jl")
            include("identify.jl")
            end
            """,
        "identify.jl" => """
            # ABOUTME: identify, the best
            # ABOUTME: law.
            identify(x) = Fitting.fit(x)
            """,
        "Fitting/Fitting.jl" => """
            \"\"\"Fits a law.\"\"\"
            module Fitting
            public fit, Problem
            include("problem_type.jl")
            include("loss.jl")
            include("fit.jl")
            include("report.jl")
            end
            """,
        "Fitting/problem_type.jl" => problem_type,
        "Fitting/loss.jl" => """
            # ABOUTME: The loss.
            # ABOUTME: Scaled.
            const _SCALE = 2.0
            macro _twice(x)
                :(2 * \$(esc(x)))
            end
            @inline _loss(p::Problem; kind::_Kind = _PLAIN) = _SCALE * _size(p) + p.fit
            """,
        "Fitting/fit.jl" => """
            # ABOUTME: fit, which minimizes
            # ABOUTME: the loss.
            \"\"\"Fits `p`, scaled by `_SCALE`.\"\"\"
            function fit(p::AbstractProblem)
                # _size(p) is not called here.
                @_twice _loss(p; _size = 1)
            end
            """,
        "Fitting/report.jl" => """
            # ABOUTME: The report of
            # ABOUTME: a fit.
            fit(p, n::AbstractProblem = p) = fit(p)
            report(p) = (_size = :_size, text = "_loss", code = :(f(\$_SCALE)),
                         value = "\$(_loss(p))")
            """,
    )
    @testset "file_scans: what each file defines and takes from its siblings" begin
        _testhelper_with_tree(file_tree...) do root
            scans = _testhelper_file_scans(root)
            # A keyword or named-tuple field, a field after a dot, a quoted symbol, a
            # string, a comment and a docstring refer to nothing; a quote's and a string's
            # interpolations do, and so does the type of an argument with a default. A
            # function both files define is taken by neither.
            @test [(f.mod, f.path, f.takes) for f in scans] == [
                (:EnzymeRates, "identify.jl", []),
                (:Fitting, "Fitting/problem_type.jl", []),
                (:Fitting, "Fitting/loss.jl",
                 ["Fitting/problem_type.jl" => [:Problem, :_Kind, :_PLAIN, :_size]]),
                (:Fitting, "Fitting/fit.jl",
                 ["Fitting/problem_type.jl" => [:AbstractProblem],
                  "Fitting/loss.jl" => [Symbol("@_twice"), :_loss]]),
                (:Fitting, "Fitting/report.jl",
                 ["Fitting/problem_type.jl" => [:AbstractProblem],
                  "Fitting/loss.jl" => [:_SCALE, :_loss]])]
            @test [f.defines for f in scans] == [Set([:identify]),
                Set([:AbstractProblem, :Problem, :_Kind, :_PLAIN, :_ALLOSTERIC, :_size]),
                Set([:_SCALE, Symbol("@_twice"), :_loss]), Set([:fit]),
                Set([:fit, :report])]
            @test !any(f -> f.includes, scans)
        end
    end

    named_tree = (
        "EnzymeRates.jl" => "module EnzymeRates\nexport fit, identify\nend",
        "identify.jl" => "identify(x) = 1",
        "Fitting/Fitting.jl" => "module Fitting\npublic fit, FitResult, @fitted\nend",
        "Fitting/fit_result_type.jl" => "struct FitResult{T}\n    x::T\nend",
        "Fitting/fitted_macro.jl" => "macro fitted(x)\n    x\nend",
        "Fitting/fit.jl" => "\"\"\"Fits.\"\"\"\nfit(x) = 1",
        "Fitting/loss.jl" => "_loss(x) = 1",
    )
    @testset "a clean tree passes the file rules" begin
        _testhelper_with_tree(named_tree...) do root
            @test _testhelper_file_name_and_include_problems(root) == String[]
        end
        _testhelper_with_tree(file_tree...) do root
            @test _testhelper_file_name_and_include_problems(root) == String[]
            @test _testhelper_one_way_problems(root, Set{Tuple{String, String}}()) ==
                  String[]
        end
    end
    @testset "the name and include rules find each fault and name the file" begin
        _testhelper_with_tree(named_tree...,
            "identify.jl" => "_identify(x) = 1",
            "Fitting/fit.jl" => "_fit(x) = 1",
            # The struct moved to another file; an outer constructor stayed.
            "Fitting/fit_result_type.jl" => "FitResult(x::Int) = FitResult{Int}(x)",
            "Fitting/fitted_macro.jl" => "macro fit_it(x)\n    x\nend",
            "Fitting/loss.jl" =>
                "include(\"more.jl\")\n_loss(x) = 1\nstruct FitResult{T}\n    x::T\nend",
            "Fitting/outcome_type.jl" => "struct FitOutcome end",
            "Fitting/utils.jl" => "_u(x) = 1") do root
            @test _testhelper_file_name_and_include_problems(root) == [
                "identify.jl is named for `identify`, a public name of EnzymeRates, but " *
                "does not define it",
                "Fitting/fit.jl is named for `fit`, a public name of Fitting, but does " *
                "not define it",
                "Fitting/fit_result_type.jl ends in `_type` but defines no type whose " *
                "name is `fit_result` in snake case",
                "Fitting/fitted_macro.jl ends in `_macro` but defines no macro `@fitted`",
                "Fitting/loss.jl calls `include`; only a module's header includes files",
                "Fitting/outcome_type.jl ends in `_type` but defines no type whose name " *
                "is `outcome` in snake case",
                "Fitting/utils.jl has a name that says nothing of what it holds"]
        end
    end
    @testset "the one-way rule finds a mutual pair and a stale exemption" begin
        pair = ("Fitting/loss.jl", "Fitting/problem_type.jl")
        _testhelper_with_tree(file_tree...,
            "Fitting/problem_type.jl" => problem_type * "_check(p) = _loss(p)\n") do root
            @test _testhelper_one_way_problems(root, Set{Tuple{String, String}}()) == [
                "Fitting/loss.jl and Fitting/problem_type.jl take names from each " *
                "other: loss.jl takes `Problem`, `_Kind`, `_PLAIN`, `_size`; " *
                "problem_type.jl takes `_loss`"]
            @test _testhelper_one_way_problems(root, Set([pair])) == String[]
        end
        _testhelper_with_tree(file_tree...) do root
            @test _testhelper_one_way_problems(root, Set([pair])) == [
                "Fitting/loss.jl and Fitting/problem_type.jl no longer take names from " *
                "each other; remove the pair from the exemption set in " *
                "test/test_boundaries.jl"]
        end
    end
    @testset "the ABOUTME rule finds a code file without its two lines" begin
        _testhelper_with_tree(
            "src/A.jl" => "# ABOUTME: One.\n# ABOUTME: Two.\nmodule A end",
            "src/A/b.jl" => "# ABOUTME: One.\nb() = 1",
            "ext/AExt.jl" => "# ABOUTME: One.\r\n# ABOUTME: Two.\r\n",
            "test/t.jl" => "t() = 1",
            "test/three.jl" => "# ABOUTME: One.\n# ABOUTME: Two.\n# ABOUTME: Three.\n",
            "docs/make.jl" => "# ABOUTME: One.\n# ABOUTME: Two.\n",
            "docs/src/index.md" => "text",
            # The built site and the local results are not the package's code.
            "docs/build/x.jl" => "x() = 1",
            "docs/hpc_results/r.jl" => "r() = 1",
            "scripts/s.jl" => "s() = 1") do root
            @test _testhelper_aboutme_problems(root) == [
                "src/A/b.jl does not start with exactly two `# ABOUTME: ` lines",
                "test/t.jl does not start with exactly two `# ABOUTME: ` lines",
                "test/three.jl does not start with exactly two `# ABOUTME: ` lines"]
            @test _testhelper_aboutme(joinpath(root, "ext", "AExt.jl")) == "One. Two."
            @test _testhelper_aboutme(joinpath(root, "src", "A", "b.jl")) === nothing
        end
    end
    @testset "file table" begin
        expected = join([
            "| Module | File | What it holds | Public names it defines | " *
            "Takes names from |",
            "|:--|:--|:--|:--|:--|",
            "| top level | `identify.jl` | identify, the best law. | ★ `identify` " *
            "|  |",
            "| Fitting | `problem_type.jl` | Problem and \\_size, the data a fit " *
            "reads. | `Problem` |  |",
            "| Fitting | `loss.jl` | The loss. Scaled. |  | `problem_type` (4) |",
            "| Fitting | `fit.jl` | fit, which minimizes the loss. | ★ `fit` | " *
            "`problem_type` (1), `loss` (2) |",
            "| Fitting | `report.jl` | The report of a fit. | ★ `fit` | " *
            "`problem_type` (1), `loss` (2) |", ""], "\n")
        _testhelper_with_tree(file_tree...) do root
            block = _testhelper_FILE_TABLE_BLOCK
            @test _testhelper_file_table(root) == expected
            doc = "text\n\n$(block.opening)$(expected)$(block.closing)\nmore\n"
            @test _testhelper_generated_block_problems(doc, block, root) == String[]
            stale = replace(doc, "(4)" => "(3)")
            @test _testhelper_generated_block_problems(stale, block, root) == [
                "the file table in docs/src/developer.md differs from the one generated " *
                "from src/; rewrite it with `julia test/update_module_map.jl`"]
            @test _testhelper_replace_generated_block(stale, block, expected) == doc
        end
    end
end

@testset "module boundaries" begin
    _testhelper_check_boundaries(_testhelper_SRC)
end

@testset "every code file starts with exactly two ABOUTME lines" begin
    @test _testhelper_aboutme_problems(joinpath(@__DIR__, "..")) == String[]
end

@testset "source files are named for what they give; only headers include files" begin
    @test _testhelper_file_name_and_include_problems(_testhelper_SRC) == String[]
end

# The pairs of files of a module that take names from each other today, each pair in
# alphabetical order. ROADMAP.md lists each as a known issue, with the names that cross
# and what would make the calls one-way; remove a pair here when they are.
const _testhelper_FILES_USING_EACH_OTHER = Set([
    ("Derivation/allosteric_rate_equation_correctness_guard.jl",
     "Derivation/rate_equation_correctness_guard.jl"),
    ("Derivation/rate_equation_correctness_guard.jl",
     "Derivation/thermodynamic_constraints.jl"),
])

@testset "the files of a module take names from each other one way" begin
    @test _testhelper_one_way_problems(_testhelper_SRC,
        _testhelper_FILES_USING_EACH_OTHER) == String[]
end

@testset "every public name is defined" begin
    for m in _testhelper_MODULES, n in _testhelper_api(_testhelper_SRC, m)
        @test isdefined(getfield(EnzymeRates, m), n)
    end
end

@testset "module map checker" begin
    tree = (
        "Mechanisms/Mechanisms.jl" => "public Mechanism, compile_mechanism",
        "Fitting/Fitting.jl" => "public fit",
        "Fitting/fit.jl" => "g(m) = Mechanisms.compile_mechanism(m)",
        "EnzymeRates.jl" => "export fit, identify\nh(m) = Fitting.fit(m)",
    )
    doc = """
        **`Mechanisms`**

        | Name | What it does | Used by | Exported |
        |:-----|:-------------|:--------|:---------|
        | `Mechanism` | A mechanism | | |
        | `compile_mechanism` | Lifts it | | |
        | | | Fitting | |

        **`Fitting`**

        | Name | What it does | Used by | Exported |
        |:-----|:-------------|:--------|:---------|
        | `fit` | Fits | top level | ✓ |

        **The top level**

        | Name | What it does | Used by | Exported |
        |:-----|:-------------|:--------|:---------|
        | `identify` | Identifies | | ✓ |
        """
    broken = replace(doc,
        "| `Mechanism` | A mechanism | | |\n" => "",
        "| `fit` | Fits | top level | ✓ |" =>
            "| `fit` | Fits | top level, Mechanisms | |\n| `fold` | Folds | | |",
        "| | | Fitting | |\n" => "")
    _testhelper_with_tree(tree...) do root
        @test _testhelper_module_map_problems(doc, root) == String[]
        @test _testhelper_module_map_problems(broken, root) == [
            "interface table Fitting: extra `fold`",
            "interface table Mechanisms: missing `Mechanism`",
            "✓ in interface table Fitting: missing `fit`",
            "Used by, row compile_mechanism: missing `Fitting`",
            "Used by, row fit: extra `Mechanisms`"]
    end
end

@testset "module map diagram and table" begin
    tree = (
        "Mechanisms/Mechanisms.jl" =>
            "\"\"\"Stores things. Hides more.\"\"\"\nmodule Mechanisms\n" *
            "public Mechanism, compile_mechanism, metabolites\nend",
        "Mechanisms/types.jl" => "struct Mechanism\n    x::Int\nend",
        "Fitting/Fitting.jl" =>
            "\"\"\"Fits a law\nto data. Hides more.\"\"\"\nmodule Fitting\n" *
            "public fit\nend",
        "Fitting/fit.jl" => "g(m) = Mechanisms.compile_mechanism(m)\n" *
            "Mechanisms.Mechanism\nMechanisms.metabolites(m)",
        "EnzymeRates.jl" => "\"\"\"Does it all.\"\"\"\nmodule EnzymeRates\n" *
            "export fit, identify, Result, metabolites\nh(m) = Fitting.fit(m)\n" *
            "struct Result end\nend",
    )
    expected = _testhelper_MAP_HEADER * join([
        "flowchart TB",
        "  top[\"<b><big>EnzymeRates</big></b><br/>Does it all.<br/><br/>" *
        "Defines:<br/><b>★#nbsp;identify()</b><br/><b>★#nbsp;Result</b>\"]",
        "  fit[\"<b><big>Fitting</big></b><br/>Fits a law<br/>to data.\"]",
        "  mech[\"<b><big>Mechanisms</big></b><br/>Stores things.<br/><br/>" *
        "Types every module but Data reads:<br/>Mechanism\"]",
        "  top -->|\"<b>★#nbsp;fit()</b>\"| fit",
        "  fit -->|\"<b>★#nbsp;metabolites()</b><br/>compile_mechanism()\"| mech",
        "  classDef mod fill:#eaf2ee,stroke:#666666,color:#222222;",
        "  class top,fit,mech mod;", ""], "\n")
    expected_table = join([
        "| Module | What it does | Public names (★ exported) | May call |",
        "|:--|:--|:--|:--|",
        "| EnzymeRates (top level) | Does it all. | ★ `identify`, ★ `Result` | all |",
        "| Fitting | Fits a law to data. | ★ `fit` | Mechanisms, Derivation, Data |",
        "| Mechanisms | Stores things. | ★ `metabolites`, `Mechanism`, " *
        "`compile_mechanism` | none |", ""], "\n")
    _testhelper_with_tree(tree...) do root
        map, table = _testhelper_MAP_BLOCK, _testhelper_TABLE_BLOCK
        @test _testhelper_module_map(root) == expected
        doc = "text\n\n```mermaid\n$expected```\n\nmore\n"
        @test _testhelper_generated_block_problems(doc, map, root) == String[]
        # A Windows checkout writes CRLF line endings.
        @test _testhelper_generated_block_problems(
            replace(doc, "\n" => "\r\n"), map, root) == String[]
        stale = replace(doc, "Stores things." => "Stores.")
        @test _testhelper_generated_block_problems(stale, map, root) == [
            "the Mermaid module map in docs/src/developer.md differs from the one " *
            "generated from src/; rewrite it with `julia test/update_module_map.jl`"]
        @test _testhelper_replace_generated_block(stale, map, expected) == doc
        @test _testhelper_generated_block_problems("no diagram", map, root) == [
            "the Mermaid module map in docs/src/developer.md differs from the one " *
            "generated from src/; rewrite it with `julia test/update_module_map.jl`"]

        @test _testhelper_module_table(root) == expected_table
        doc = "text\n\n$(table.opening)$(expected_table)$(table.closing)\nmore\n"
        @test _testhelper_generated_block_problems(doc, table, root) == String[]
        @test _testhelper_generated_block_problems(
            replace(doc, "\n" => "\r\n"), table, root) == String[]
        stale = replace(doc, "Stores things." => "Stores.")
        differs = ["the module table in .claude/CLAUDE.md differs from the one generated " *
                   "from src/; rewrite it with `julia test/update_module_map.jl`"]
        @test _testhelper_generated_block_problems(stale, table, root) == differs
        @test _testhelper_replace_generated_block(stale, table, expected_table) == doc
        @test _testhelper_replace_generated_block(
            replace(stale, "\n" => "\r\n"), table, expected_table) == doc
        @test _testhelper_generated_block_problems("no table", table, root) == differs
        # A character of several bytes just before the closing line.
        @test _testhelper_generated_block_problems(
            "$(table.opening)★$(table.closing)", table, root) == differs
        @test_throws "the module table in .claude/CLAUDE.md has no generated block" (
            _testhelper_replace_generated_block("no table", table, expected_table))
    end
end

@testset "the generated module map and tables match src/" begin
    root = joinpath(@__DIR__, "..")
    doc = read(joinpath(root, "docs", "src", "developer.md"), String)
    @test _testhelper_module_map_problems(doc, _testhelper_SRC) == String[]
    @test _testhelper_generated_block_problems(
        doc, _testhelper_MAP_BLOCK, _testhelper_SRC) == String[]
    @test _testhelper_generated_block_problems(
        doc, _testhelper_FILE_TABLE_BLOCK, _testhelper_SRC) == String[]
    @test _testhelper_generated_block_problems(
        read(joinpath(root, ".claude", "CLAUDE.md"), String), _testhelper_TABLE_BLOCK,
        _testhelper_SRC) == String[]
end

@testset "Mechanisms API is the data model" begin
    # The types of a mechanism and a reaction, and the functions that compute on them:
    # each conformation's steps, a step's text, and the rapid-equilibrium segments and
    # turnover cycles of a step mask.
    @test _testhelper_api(_testhelper_SRC, :Mechanisms) == Set([:Step, :Mechanism,
        :CatalyticScheme, :Allostery, :RegulatorySite, :Reactant, :Regulator,
        :EnzymeReaction, :metabolites, :step_text,
        :steps_running_in_allosteric_conformations, :rapid_equilibrium_segments,
        :has_turnover_cycle, :has_nonzero_weight_cycle])
end

@testset "documented internal paths resolve" begin
    @test isdefined(EnzymeRates, :test_rate_equation)
    @test isdefined(EnzymeRates.Enumeration, :candidate_mechanisms)
end

# REPL help (`?name`) finds each exported name's docstring through its top-level binding.
@testset "every exported name shows its docstring" begin
    for n in names(EnzymeRates)
        n === :EnzymeRates && continue
        text = string(Base.Docs.doc(Base.Docs.Binding(EnzymeRates, n)))
        @test !occursin("No documentation found", text)
    end
end

"""Where the files under `test/` of the tree at `root` break the test-file naming rule, one
line each:
- a folder of `test/` that mirrors no source folder and holds no data (`fixtures`,
  `reference`);
- in a module's test folder or `test/ext/`, a file not named `test_<file>.jl`, or a file
  in a folder below;
- a test file `test_<file>.jl` with no source file `<file>.jl` in the mirrored folder. A
  test file directly in `test/` may instead be named for a package-wide rule it checks,
  such as `test_boundaries.jl`;
- a source file, module headers aside, with no test file."""
function _testhelper_test_name_problems(root)
    package_wide_rules = ("boundaries", "type_stability", "compile_budget", "aqua_jet",
                          "oracles")
    data_folders = ("fixtures", "reference")
    entries(dir) = isdir(joinpath(root, dir)) ? readdir(joinpath(root, dir)) : String[]
    jl_files(dir) = filter(endswith(".jl"), entries(dir))
    folders(dir) = filter(name -> isdir(joinpath(root, dir, name)), entries(dir))
    modules = folders("src")
    problems = String[]
    for folder in folders("test")
        if !(folder in modules || folder == "ext" || folder in data_folders)
            push!(problems, "test/$folder mirrors no source folder")
        end
    end
    mirrored = [("src", "test"), ("ext", "test/ext")]
    append!(mirrored, [("src/$m", "test/$m") for m in modules])
    for (source_folder, test_folder) in mirrored
        header = basename(source_folder) * ".jl"
        sources = filter(f -> f != header && f != "EnzymeRates.jl", jl_files(source_folder))
        test_files = jl_files(test_folder)
        if test_folder != "test"
            for file in filter(!startswith("test_"), test_files)
                push!(problems, "$test_folder/$file is not named test_<file>.jl")
            end
            for below in folders(test_folder)
                for (dir, _, files) in walkdir(joinpath(root, test_folder, below))
                    for file in filter(endswith(".jl"), files)
                        path = relpath(joinpath(dir, file), root)
                        push!(problems, "$path is in a folder below $test_folder")
                    end
                end
            end
        end
        for file in filter(startswith("test_"), test_files)
            source = chopprefix(file, "test_")
            checks_a_rule = test_folder == "test" &&
                            chopsuffix(source, ".jl") in package_wide_rules
            if !(source in sources || checks_a_rule)
                push!(problems, "$test_folder/$file names no $source_folder/$source")
            end
        end
        for source in sources
            if !("test_$source" in test_files)
                push!(problems, "$source_folder/$source has no $test_folder/test_$source")
            end
        end
    end
    problems
end

# A test file is test_<x>.jl for the source file <x>.jl whose functions its testsets are
# about, in the mirrored folder; a test file directly in test/ may instead check a
# package-wide rule.
@testset "every test file is named for the source file it tests" begin
    @test _testhelper_test_name_problems(joinpath(@__DIR__, "..")) == String[]
end
