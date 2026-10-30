# ABOUTME: Checks the module boundaries of src/ by parsing it (API lists, layers, underscore
# ABOUTME: names, exports), the developer docs' module map, and the public names' docs.

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
interface tables (names, "Used by" and ✓ marks) against each `API` list, the cross-module
uses and the exports."""
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
    @testset "every cross-module name is in its API list" begin
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
        "Mechanisms/Mechanisms.jl" => "const API = (:Mechanism, :compile_mechanism)",
        "Mechanisms/types.jl" => "struct Mechanism\n    edges::Vector{Int}\nend",
        "Derivation/Derivation.jl" => "const API = (:check_valid_rate_equation,)",
        "Derivation/guard.jl" =>
            "f(m::Mechanisms.Mechanism) = Mechanisms.compile_mechanism(m)",
        "Fitting/Fitting.jl" => "const API = (:fit,)",
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
            "DSL/DSL.jl" => "const API = (:parse, Symbol(\"@enzyme_mechanism\"))") do root
            @test _testhelper_api(root, :DSL) == Set([:parse, Symbol("@enzyme_mechanism")])
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
    @testset "rule 1 finds a name outside the API list, in a module and at the top" begin
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
    @testset "rule 3 finds underscore names in an API list or used across modules" begin
        _testhelper_with_tree(
            "Mechanisms/Mechanisms.jl" => "const API = (:Mechanism, :_raw)",
            "Derivation/Derivation.jl" => "const API = (:_guard, :ok, :_guard)",
            "Fitting/f.jl" => "f() = Derivation._guard()",
            "Fitting/Fitting.jl" => "const API = (:fit,)") do root
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

@testset "module boundaries" begin
    _testhelper_check_boundaries(_testhelper_SRC)
end

@testset "module map checker" begin
    tree = (
        "Mechanisms/Mechanisms.jl" => "const API = (:Mechanism, :compile_mechanism)",
        "Fitting/Fitting.jl" => "const API = (:fit,)",
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

@testset "module map diagram" begin
    tree = (
        "Mechanisms/Mechanisms.jl" =>
            "\"\"\"Stores things. Hides more.\"\"\"\nmodule Mechanisms\n" *
            "const API = (:Mechanism, :compile_mechanism, :metabolites)\nend",
        "Mechanisms/types.jl" => "struct Mechanism\n    x::Int\nend",
        "Fitting/Fitting.jl" =>
            "\"\"\"Fits a law\nto data. Hides more.\"\"\"\nmodule Fitting\n" *
            "const API = (:fit,)\nend",
        "Fitting/fit.jl" => "g(m) = Mechanisms.compile_mechanism(m)\n" *
            "Mechanisms.Mechanism\nMechanisms.metabolites(m)",
        "EnzymeRates.jl" => "\"\"\"Does it all.\"\"\"\nmodule EnzymeRates\n" *
            "export fit, identify, Result, metabolites\nh(m) = Fitting.fit(m)\n" *
            "struct Result end\nend",
    )
    expected = _testhelper_MAP_HEADER * join([
        "flowchart TB",
        "  top[\"<b><big>EnzymeRates</big></b><br/>Does it all.<br/><br/>" *
        "Defines <b>★ identify()</b>, <b>★ Result</b>\"]",
        "  fit[\"<b><big>Fitting</big></b><br/>Fits a law<br/>to data.\"]",
        "  mech[\"<b><big>Mechanisms</big></b><br/>Stores things.<br/><br/>" *
        "Types every module but Data reads: Mechanism\"]",
        "  top -->|\"<b>★ fit()</b>\"| fit",
        "  fit -->|\"<b>★ metabolites()</b><br/>compile_mechanism()\"| mech",
        "  classDef mod fill:#eaf2ee,stroke:#666666,color:#222222;",
        "  class top,fit,mech mod;", ""], "\n")
    _testhelper_with_tree(tree...) do root
        @test _testhelper_module_map(root) == expected
        doc = "text\n\n```mermaid\n$expected```\n\nmore\n"
        @test _testhelper_module_map_block_problems(doc, root) == String[]
        stale = replace(doc, "Stores things." => "Stores.")
        @test _testhelper_module_map_block_problems(stale, root) == [
            "the Mermaid module map in docs/src/developer.md differs from the one " *
            "generated from src/; rewrite it with `julia test/update_module_map.jl`"]
        @test _testhelper_replace_module_map(stale, expected) == doc
        @test !isempty(_testhelper_module_map_block_problems("no diagram", root))
    end
    @testset "wrapping names into lines" begin
        @test _testhelper_wrap_names("A: ", ["aaa", "bbb", "ccc"], 10) ==
              "A: aaa,<br/>bbb, ccc"
    end
end

@testset "the developer docs' module map matches src/" begin
    doc = read(joinpath(@__DIR__, "..", "docs", "src", "developer.md"), String)
    @test _testhelper_module_map_problems(doc, _testhelper_SRC) == String[]
    @test _testhelper_module_map_block_problems(doc, _testhelper_SRC) == String[]
end

@testset "Mechanisms API is the data model" begin
    # The types of a mechanism and a reaction, the functions that compute on them, and the
    # two name checks the DSL runs while it reads a block, so that their errors name the
    # macro: a conformation label's shape and the decodability of the declared names.
    @test _testhelper_api(_testhelper_SRC, :Mechanisms) == Set([:Met, :Node, :Step,
        :Mechanism, :RegulatorySite, :AllostericMechanism, :Reactant, :Regulator,
        :EnzymeReaction, :EnzymeMechanism, :AllostericEnzymeMechanism,
        :AbstractEnzymeMechanism, :step_kind, :groups, :metabolites, :modifiers, :catalytic,
        :compile_mechanism, :form_text, :step_text, :step_sides, :runs_toward_products,
        :check_conformation_label, :check_decodable_names])
end

@testset "documented internal paths resolve" begin
    @test isdefined(EnzymeRates, :Mechanism)
    @test isdefined(EnzymeRates, :AllostericMechanism)
    @test isdefined(EnzymeRates, :compile_mechanism)
    @test isdefined(EnzymeRates.Enumeration, :init_mechanisms)
end

# REPL help (`?name`) finds each exported name's docstring through its top-level binding.
@testset "every exported name shows its docstring" begin
    for n in names(EnzymeRates)
        n === :EnzymeRates && continue
        text = string(Base.Docs.doc(Base.Docs.Binding(EnzymeRates, n)))
        @test !occursin("No documentation found", text)
    end
end
