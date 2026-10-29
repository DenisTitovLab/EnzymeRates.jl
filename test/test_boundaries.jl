# ABOUTME: Checks the module boundaries of src/ by parsing it (API lists, layers, underscore
# ABOUTME: names, exports), the developer docs' module map, and the public names' docs.

# The six submodules, each a folder src/<Module>/ with src/<Module>/<Module>.jl.
const _testhelper_MODULES = (:Data, :Mechanisms, :Derivation, :Enumeration, :Fitting,
    :DSL)

# The modules each module may name through `Module.name`. The top level (files outside
# src/<Module>/, owner :EnzymeRates) may name all six.
const _testhelper_LAYERS = Dict(
    :Data => Symbol[],
    :Mechanisms => Symbol[],
    :Derivation => [:Mechanisms],
    :Enumeration => [:Mechanisms, :Derivation],
    :Fitting => [:Mechanisms, :Derivation, :Data],
    :DSL => [:Mechanisms, :Derivation],
)

"""Calls `visit(ex, line)` on every `Expr` below `ex`, children before parents; `line` is
the source line of the statement holding `ex`."""
function _testhelper_walk(visit, ex, line = Ref(0))
    if ex isa LineNumberNode
        line[] = ex.line
    elseif ex isa Expr
        here = line[]
        foreach(a -> _testhelper_walk(visit, a, line), ex.args)
        visit(ex, here)
    end
end

"""The module `M` of a `M.name` or `EnzymeRates.M.name` node (one of the six modules), else
`nothing`. A node of one argument (`using X: name` lists) is no `M.name`, nor is a quoted
path that interpolates its name (`:(M.\$T)`)."""
function _testhelper_module_head(ex)
    (ex.head === :. && length(ex.args) == 2 && ex.args[2] isa QuoteNode &&
     ex.args[2].value isa Symbol) || return nothing
    head = ex.args[1]
    head isa Expr && head.head === :. && length(head.args) == 2 &&
        head.args[1] === :EnzymeRates && head.args[2] isa QuoteNode &&
        (head = head.args[2].value)
    head in _testhelper_MODULES ? head : nothing
end

"""`(module, name)` for each name that a `using Path.Mod: a, b` or `import Path.Mod: a`
statement takes from one of the six modules; `Mod` itself is no use."""
function _testhelper_imported_names(ex)
    (length(ex.args) == 1 && ex.args[1] isa Expr && ex.args[1].head === :(:)) || return []
    path, items... = ex.args[1].args
    mod = path.args[end]
    mod in _testhelper_MODULES || return []
    names = [(item.head === :as ? item.args[1] : item).args[end] for item in items]
    [(mod, name) for name in names if name !== mod]
end

"""`(module, name, line)` for each `Mod.name` in the source text, `Mod` one of the six
modules (also written `EnzymeRates.Mod.name`), and for each name imported from one. A type
annotation `x::Mod.T` and a method extension `function Mod.f(...)` are uses of `T` and
`f`."""
function _testhelper_cross_refs(src::AbstractString)
    refs = Tuple{Symbol, Symbol, Int}[]
    _testhelper_walk(Meta.parseall(src)) do ex, line
        mod = _testhelper_module_head(ex)
        if mod !== nothing
            push!(refs, (mod, ex.args[2].value, line))
        elseif ex.head in (:using, :import)
            for (mod, name) in _testhelper_imported_names(ex)
                push!(refs, (mod, name, line))
            end
        end
    end
    refs
end

_testhelper_source_files(dir) = sort!(String[joinpath(d, f)
    for (d, _, fs) in walkdir(dir) for f in fs if endswith(f, ".jl")])

"""Module name to the `.jl` files under `root/<Module>/`, for each module with a folder."""
_testhelper_module_sources(root) = Dict(
    m => _testhelper_source_files(joinpath(root, string(m)))
    for m in _testhelper_MODULES if isdir(joinpath(root, string(m))))

"""The name an `API` entry lists: `:name`, or `Symbol("@macro")` for a macro."""
_testhelper_api_name(q) = q isa QuoteNode ? q.value : Symbol(q.args[2])

"""The names in the `const API = (...)` tuple of `root/<mod>/<mod>.jl`."""
function _testhelper_api(root, mod::Symbol)
    path = joinpath(root, string(mod), string(mod) * ".jl")
    api = nothing
    _testhelper_walk(Meta.parseall(read(path, String))) do ex, _
        if ex.head === :const && ex.args[1] isa Expr && ex.args[1].head === :(=) &&
           ex.args[1].args[1] === :API
            api = Set{Symbol}(_testhelper_api_name(q) for q in ex.args[1].args[2].args)
        end
    end
    api === nothing && error("$path declares no `const API = (...)`")
    api
end

"""The path of `file` relative to `root`, with `/` as separator on every platform."""
_testhelper_relpath(root, file) = replace(relpath(file, root), '\\' => '/')

"""The module owning a file under `root`: the folder it is in, or `:EnzymeRates` for the
top level."""
function _testhelper_owner(root, file)
    first_part = Symbol(first(splitpath(_testhelper_relpath(root, file))))
    first_part in _testhelper_MODULES ? first_part : :EnzymeRates
end

"""Every cross-module use under `root`, as (user, file, module, name, line), leaving out
a module naming itself."""
function _testhelper_uses(root)
    uses = NamedTuple{(:user, :file, :mod, :name, :line),
        Tuple{Symbol, String, Symbol, Symbol, Int}}[]
    for file in _testhelper_source_files(root), (mod, name, line) in
        _testhelper_cross_refs(read(file, String))
        user = _testhelper_owner(root, file)
        user === mod || push!(uses,
            (; user, file = _testhelper_relpath(root, file), mod, name, line))
    end
    uses
end

"""Uses naming a module's name that its API list lacks (or naming a module with no
folder)."""
function _testhelper_api_violations(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    [u for u in _testhelper_uses(root) if !(u.name in get(apis, u.mod, Set{Symbol}()))]
end

"""Uses by a module of a sibling that its layer forbids; the top level is unrestricted."""
_testhelper_layer_violations(root) = [u for u in _testhelper_uses(root)
    if u.user !== :EnzymeRates && !(u.mod in _testhelper_LAYERS[u.user])]

"""For each module with a folder, the underscore names in its API list or used from
outside it."""
function _testhelper_underscore_violations(root)
    uses = _testhelper_uses(root)
    Dict(m => sort!(unique!(Symbol[n
        for n in Iterators.flatten((_testhelper_api(root, m),
            (u.name for u in uses if u.mod === m))) if startswith(string(n), "_")]))
        for m in keys(_testhelper_module_sources(root)))
end

"""`(file, line)` for each `export` in a module's folder: a module's names reach the others
only qualified, so an export would let `using ..Mod` bring bare names past rule 1."""
function _testhelper_export_violations(root)
    found = Tuple{String, Int}[]
    for files in values(_testhelper_module_sources(root)), file in files
        _testhelper_walk(Meta.parseall(read(file, String))) do ex, line
            ex.head === :export && push!(found, (_testhelper_relpath(root, file), line))
        end
    end
    found
end

const _testhelper_SRC = normpath(joinpath(@__DIR__, "..", "src"))

"""The names in the `export` lines of `file`."""
function _testhelper_exports(file)
    names = Set{Symbol}()
    _testhelper_walk(Meta.parseall(read(file, String))) do ex, _
        ex.head === :export && union!(names, ex.args)
    end
    names
end

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

"""How the docs write a module or a (user, used) pair of modules."""
_testhelper_doc_label(m::Symbol) = m === :EnzymeRates ? "top level" : string(m)
_testhelper_doc_label(pair::Tuple) = join(_testhelper_doc_label.(pair), " → ")

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
uses and the exports, the export table against the `export` lines, and the edge table
against the cross-module uses."""
function _testhelper_module_map_problems(doc, root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    exports = _testhelper_exports(joinpath(root, "EnzymeRates.jl"))
    owner(n) = something(findfirst(api -> n in api, apis), :EnzymeRates)
    apis[:EnzymeRates] = Set(n for n in exports if owner(n) === :EnzymeRates)
    names, ticks, used_by, exported, edges = (Dict{Any, Set{Symbol}}() for _ in 1:5)
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
        elseif t.header == ["Name", "Module", "What it does"]
            for r in t.rows
                union!(get!(exported, _testhelper_doc_module(r[2]), Set{Symbol}()),
                    _testhelper_cell_names(r[1]))
            end
        elseif t.header == ["From → To", "Names"]
            pair = nothing
            for r in t.rows
                isempty(r[1]) ||
                    (pair = Tuple(_testhelper_doc_module.(strip.(split(r[1], "→")))))
                union!(get!(edges, pair, Set{Symbol}()), _testhelper_cell_names(r[2]))
            end
        end
    end
    uses, users = Dict{Any, Set{Symbol}}(), Dict{Any, Set{Symbol}}()
    for u in _testhelper_uses(root)
        push!(get!(uses, (u.user, u.mod), Set{Symbol}()), u.name)
        push!(get!(users, u.name, Set{Symbol}()), u.user)
    end
    exported_api = Dict(m => intersect(api, exports) for (m, api) in apis)
    [_testhelper_name_diffs("interface table", names, apis);
     _testhelper_name_diffs("✓ in interface table", ticks, exported_api);
     _testhelper_name_diffs("Used by, row", used_by, users);
     _testhelper_name_diffs("export table, module", exported, exported_api);
     _testhelper_name_diffs("edge table, row", edges, uses)]
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
        ### What `using EnzymeRates` exports

        | Name | Module | What it does |
        |:-----|:-------|:-------------|
        | `fit` | Fitting | Fits |
        | `identify` | top level | Identifies |

        | From → To | Names |
        |:----------|:------|
        | top level → Fitting | `fit` |
        | Fitting → Mechanisms | |
        | | `compile_mechanism` |

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
        "| | | Fitting | |\n" => "",
        "| `identify` | top level | Identifies |" =>
            "| `identify` | Fitting | Identifies |",
        "| | `compile_mechanism` |\n" => "| | `Mechanism` |\n| DSL → Fitting | `fit` |\n")
    _testhelper_with_tree(tree...) do root
        @test _testhelper_module_map_problems(doc, root) == String[]
        @test _testhelper_module_map_problems(broken, root) == [
            "interface table Fitting: extra `fold`",
            "interface table Mechanisms: missing `Mechanism`",
            "✓ in interface table Fitting: missing `fit`",
            "Used by, row compile_mechanism: missing `Fitting`",
            "Used by, row fit: extra `Mechanisms`",
            "export table, module top level: missing `identify`",
            "export table, module Fitting: extra `identify`",
            "edge table, row DSL → Fitting: extra `fit`",
            "edge table, row Fitting → Mechanisms: missing `compile_mechanism`; " *
            "extra `Mechanism`"]
    end
end

@testset "the developer docs' module map matches src/" begin
    doc = read(joinpath(@__DIR__, "..", "docs", "src", "developer.md"), String)
    @test _testhelper_module_map_problems(doc, _testhelper_SRC) == String[]
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
