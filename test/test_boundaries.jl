# ABOUTME: Checks the module boundaries of src/ by parsing it (API lists, layering, field
# ABOUTME: access, underscore names), and that public names resolve and show their docs.

# The four submodules, each a folder src/<Module>/ with src/<Module>/<Module>.jl.
const _testhelper_MODULES = (:Mechanisms, :Derivation, :Enumeration, :Fitting)

# The modules each module may name through `Module.name`. The top level (files outside
# src/<Module>/, owner :EnzymeRates) may name all four.
const _testhelper_LAYERS = Dict(
    :Mechanisms => Symbol[],
    :Derivation => [:Mechanisms],
    :Enumeration => [:Mechanisms, :Derivation],
    :Fitting => [:Mechanisms, :Derivation],
)

# The modules whose API lists and cross-module uses may still name underscore names.
const _testhelper_UNDERSCORE_PENDING = Set{Symbol}([:Mechanisms])

# Per file (path relative to src/), the expressions allowed despite a banned name,
# such as "data.group". In Derivation they read its own values: the guard's networks
# (`_Network`'s `from`, `to` and `conformation`), its verdicts (`_Verdict.conformation`),
# the MWC law (`_AlloLaw.sites`) and the steps of `segment_graph`'s result, which
# Enumeration's filters read too.
const _testhelper_FIELD_ALLOWED = Dict{String, Set{String}}(
    "Derivation/admission.jl" => Set(["eps.from", "eps.to", "net.conformation", "net.from",
        "net.to", "v.conformation"]),
    "Derivation/codegen.jl" => Set(["A.sites"]),
    "Derivation/engine.jl" => Set(["a.sites", "b.sites", "s.steps"]),
    "Derivation/guard.jl" => Set(["net.conformation", "net.from", "net.to"]),
    "Derivation/key.jl" => Set(["A.sites"]),
    "Derivation/limit_test.jl" => Set(["net.conformation", "net.from", "net.to"]),
    "Derivation/messages.jl" => Set(["net.conformation", "net.from", "v.conformation"]),
    "Derivation/print.jl" => Set(["A.sites"]),
    "Enumeration/filters.jl" => Set(["s.steps"]),
    # The seed builder reads the written form (`RawStep`, `RawForm`, `RawMechanism`, `Met`)
    # that it assembles and transforms.
    "Enumeration/seeds.jl" => Set([
        "(raw.steps[j]).consumed", "(raw.steps[j]).released", "(raw.steps[j]).to",
        "(s.consumed[1]).name", "cur.bound", "entry.consumed", "entry.from", "f.added",
        "f.bound", "f.conformation", "f.subtracted", "f1.added", "f1.conformation",
        "f1.subtracted", "f2.added", "f2.conformation", "f2.subtracted", "joined.group",
        "joined.modifiers", "joined.products", "joined.steps", "joined.substrates",
        "raw.group", "raw.modifiers", "raw.products", "raw.steps", "raw.substrates",
        "release.consumed", "release.from", "s.consumed", "s.from", "s.rapid_equilibrium",
        "s.released", "s.to", "s0.from", "s0.to", "x.name", "y.modifier", "y.name"]),
)

const _testhelper_MECHANISM_TYPES = (:Mechanism, :AllostericMechanism, :Form, :Edge, :Met,
    :Metabolites, :RegulatorySite, :RawForm, :RawStep, :RawMechanism, :RawSite,
    :RawAllostericMechanism)

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

"""The module `M` of a `M.name` or `EnzymeRates.M.name` node (one of the four modules), else
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
statement takes from one of the four modules; `Mod` itself is no use."""
function _testhelper_imported_names(ex)
    (length(ex.args) == 1 && ex.args[1] isa Expr && ex.args[1].head === :(:)) || return []
    path, items... = ex.args[1].args
    mod = path.args[end]
    mod in _testhelper_MODULES || return []
    names = [(item.head === :as ? item.args[1] : item).args[end] for item in items]
    [(mod, name) for name in names if name !== mod]
end

"""`(module, name, line)` for each `Mod.name` in the source text, `Mod` one of the four
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

"""`(name, line, expression text)` for each `x.name` in the source text with
`name ∈ fields`, apart from `Mod.name` of a module. Destructuring (`(; edges) = m`) and
`getfield` are not seen."""
function _testhelper_field_accesses(src::AbstractString, fields)
    accesses = Tuple{Symbol, Int, String}[]
    _testhelper_walk(Meta.parseall(src)) do ex, line
        if ex.head === :. && length(ex.args) == 2 && ex.args[2] isa QuoteNode &&
           ex.args[2].value in fields && _testhelper_module_head(ex) === nothing
            push!(accesses, (ex.args[2].value, line, string(ex)))
        end
    end
    accesses
end

"""`(name, line)` for each `x.name` in the source text with `name ∈ fields`."""
_testhelper_field_refs(src::AbstractString, fields) =
    [(name, line) for (name, line, _) in _testhelper_field_accesses(src, fields)]

"""The field names of each `struct` in the source text, by struct name."""
function _testhelper_struct_fields(src::AbstractString)
    structs = Dict{Symbol, Vector{Symbol}}()
    _testhelper_walk(Meta.parseall(src)) do ex, _
        ex.head === :struct || return
        name = ex.args[2]
        while name isa Expr
            name = name.args[1]
        end
        fields = get!(structs, name, Symbol[])
        for item in ex.args[3].args
            item isa Expr && item.head === :const && (item = item.args[1])
            item isa Symbol && push!(fields, item)
            item isa Expr && item.head === :(::) && item.args[1] isa Symbol &&
                push!(fields, item.args[1])
        end
    end
    structs
end

_testhelper_source_files(dir) = sort!(String[joinpath(d, f)
    for (d, _, fs) in walkdir(dir) for f in fs if endswith(f, ".jl")])

"""Module name to the `.jl` files under `root/<Module>/`, for each module with a folder."""
_testhelper_module_sources(root) = Dict(
    m => _testhelper_source_files(joinpath(root, string(m)))
    for m in _testhelper_MODULES if isdir(joinpath(root, string(m))))

"""The names in the `const API = (...)` tuple of `root/<mod>/<mod>.jl`."""
function _testhelper_api(root, mod::Symbol)
    path = joinpath(root, string(mod), string(mod) * ".jl")
    api = nothing
    _testhelper_walk(Meta.parseall(read(path, String))) do ex, _
        if ex.head === :const && ex.args[1] isa Expr && ex.args[1].head === :(=) &&
           ex.args[1].args[1] === :API
            api = Set{Symbol}(q.value for q in ex.args[1].args[2].args)
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

"""Banned field accesses outside `Mechanisms` that `allowed` does not list."""
function _testhelper_field_violations(root, ban, allowed)
    found = NamedTuple{(:user, :file, :name, :line, :text),
        Tuple{Symbol, String, Symbol, Int, String}}[]
    for file in _testhelper_source_files(root)
        user = _testhelper_owner(root, file)
        banned = get(ban, user, Set{Symbol}())
        (user === :Mechanisms || isempty(banned)) && continue
        rel = _testhelper_relpath(root, file)
        for (name, line, text) in _testhelper_field_accesses(read(file, String), banned)
            text in get(allowed, rel, Set{String}()) ||
                push!(found, (; user, file = rel, name, line, text))
        end
    end
    found
end

"""The field names of the mechanism types, from their `struct` definitions under `root`."""
function _testhelper_mechanism_fields(root)
    structs = Dict{Symbol, Vector{Symbol}}()
    for file in _testhelper_source_files(root)
        merge!(structs, _testhelper_struct_fields(read(file, String)))
    end
    missing_types = [t for t in _testhelper_MECHANISM_TYPES if !haskey(structs, t)]
    isempty(missing_types) || error("no struct definition found for $missing_types")
    Set{Symbol}(f for t in _testhelper_MECHANISM_TYPES for f in structs[t])
end

const _testhelper_SRC = normpath(joinpath(@__DIR__, "..", "src"))
const _testhelper_MECHANISM_FIELDS = _testhelper_mechanism_fields(_testhelper_SRC)

# Per owner (a module or :EnzymeRates for the top level), the field names banned there.
const _testhelper_FIELD_BAN = Dict(:Derivation => _testhelper_MECHANISM_FIELDS,
    :Enumeration => _testhelper_MECHANISM_FIELDS)

"""Runs the four rules over the tree at `root`. A module in `pending` may still hold
underscore names: its rule 4 check is `@test_broken`."""
function _testhelper_check_boundaries(root, pending, ban, allowed)
    modules = sort!(collect(keys(_testhelper_module_sources(root))))
    @testset "every cross-module name is in its API list" begin
        @test isempty(_testhelper_api_violations(root))
    end
    @testset "modules respect the layers" begin
        @test isempty(_testhelper_layer_violations(root))
    end
    @testset "no mechanism field access outside Mechanisms" begin
        @test isempty(_testhelper_field_violations(root, ban, allowed))
    end
    @testset "no cross-module underscore names" begin
        @test issubset(pending, _testhelper_MODULES)
        underscore = _testhelper_underscore_violations(root)
        for m in modules
            if m in pending
                @test_broken isempty(underscore[m])
            else
                @test isempty(underscore[m])
            end
        end
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
        @test _testhelper_cross_refs("Derivation.has_valid_rate_equation(m)") ==
              [(:Derivation, :has_valid_rate_equation, 1)]
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
        @test Set(_testhelper_cross_refs("a = 1\nFitting.fit(\"\$(Derivation.f(x))\")")) ==
              Set([(:Derivation, :f, 2), (:Fitting, :fit, 2)])
    end
    @testset "field_refs" begin
        @test _testhelper_field_refs("x.edges[k].from", Set([:edges, :from])) ==
              [(:edges, 1), (:from, 1)]
        @test _testhelper_field_refs("a = 1\nx.edges", Set([:from])) == []
        @test _testhelper_field_refs("Mechanisms.edges", Set([:edges])) == []
        @test isempty(_testhelper_field_refs("using LinearAlgebra: dot", Set([:dot])))
        @test _testhelper_field_refs("EnzymeRates.Mechanisms.edges", Set([:edges])) == []
        @test _testhelper_field_accesses("data.group", Set([:group])) ==
              [(:group, 1, "data.group")]
    end
    @testset "struct_fields" begin
        src = """
        struct A
            x::Int
            const y::Vector{Int}
            z
            A(x) = new(x, [x], x)
        end
        struct B{T} <: Any
            w::T
        end
        """
        @test _testhelper_struct_fields(src) ==
              Dict(:A => [:x, :y, :z], :B => [:w])
    end
    @testset "mechanism fields come from the struct definitions" begin
        @test issubset([:substrates, :edges, :group, :bound, :from, :to, :consumed,
            :released, :ligands, :multiplicity, :catalytic, :name, :modifier,
            :steps, :sites, :regulatory_sites, :cat_allo_states, :counts],
            _testhelper_MECHANISM_FIELDS)
    end

    cleanly_layered = (
        "Mechanisms/Mechanisms.jl" => "const API = (:Mechanism, :compile_mechanism)",
        "Mechanisms/types.jl" => "struct Mechanism\n    edges::Vector{Int}\nend",
        "Derivation/Derivation.jl" => "const API = (:has_valid_rate_equation,)",
        "Derivation/guard.jl" =>
            "f(m::Mechanisms.Mechanism) = Mechanisms.compile_mechanism(m)",
        "Fitting/Fitting.jl" => "const API = (:fit,)",
        "Fitting/fit.jl" => "g(m) = Derivation.has_valid_rate_equation(m)",
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
    end
    @testset "a clean tree passes every rule" begin
        _testhelper_with_tree(cleanly_layered...) do root
            @test isempty(_testhelper_api_violations(root))
            @test isempty(_testhelper_layer_violations(root))
            @test isempty(_testhelper_field_violations(root,
                Dict(:Derivation => Set([:edges])), Dict{String, Set{String}}()))
            @test all(isempty, values(_testhelper_underscore_violations(root)))
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
            "Mechanisms/up.jl" => "f(m) = Derivation.has_valid_rate_equation(m)",
            "Derivation/across.jl" => "f(m) = Fitting.fit(m)",
            "Fitting/ok.jl" => "f(m) = Mechanisms.compile_mechanism(m)") do root
            found = _testhelper_layer_violations(root)
            @test Set((u.user, u.mod) for u in found) ==
                  Set([(:Mechanisms, :Derivation), (:Derivation, :Fitting)])
        end
    end
    @testset "rule 3 finds a banned field outside Mechanisms unless allowed" begin
        _testhelper_with_tree(cleanly_layered...,
            "Mechanisms/uses.jl" => "f(m) = m.edges",
            "Fitting/data.jl" => "f(m, data) = (m.edges, data.group)",
            "identify.jl" => "f(m) = m.edges") do root
            ban = Dict(:Fitting => Set([:edges, :group]), :EnzymeRates => Set([:edges]))
            allowed = Dict("Fitting/data.jl" => Set(["data.group"]))
            found = _testhelper_field_violations(root, ban, allowed)
            @test Set((u.file, u.text) for u in found) == Set([
                ("Fitting/data.jl", "m.edges"), ("identify.jl", "m.edges")])
        end
    end
    @testset "rule 4 finds underscore names in an API list or used across modules" begin
        _testhelper_with_tree(
            "Mechanisms/Mechanisms.jl" => "const API = (:Mechanism, :_raw)",
            "Derivation/Derivation.jl" => "const API = (:_guard, :ok, :_guard)",
            "Fitting/f.jl" => "f() = Derivation._guard()",
            "Fitting/Fitting.jl" => "const API = (:fit,)") do root
            @test _testhelper_underscore_violations(root) == Dict(
                :Mechanisms => [:_raw], :Derivation => [:_guard], :Fitting => [])
        end
    end
end

@testset "module boundaries" begin
    _testhelper_check_boundaries(_testhelper_SRC, _testhelper_UNDERSCORE_PENDING,
        _testhelper_FIELD_BAN, _testhelper_FIELD_ALLOWED)
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
