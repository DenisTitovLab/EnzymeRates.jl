# ABOUTME: Parses src/ for the boundary tests (cross-module uses, public names, exports,
# ABOUTME: the names each file defines and takes) and generates the map and the tables.

# The submodules, each a folder src/<Module>/ with src/<Module>/<Module>.jl.
const _testhelper_MODULES = (:Data, :Mechanisms, :Derivation, :Enumeration, :Fitting,
    :DSL, :Plotting)

# The modules each module may name through `Module.name`. The top level (files outside
# src/<Module>/, owner :EnzymeRates) may name all of them.
const _testhelper_LAYERS = Dict(
    :Data => Symbol[],
    :Mechanisms => Symbol[],
    :Derivation => [:Mechanisms],
    :Enumeration => [:Mechanisms, :Derivation],
    :Fitting => [:Mechanisms, :Derivation, :Data],
    :DSL => [:Mechanisms, :Derivation],
    :Plotting => [:Mechanisms, :Derivation, :Data],
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

"""The module `M` of a `M.name` or `EnzymeRates.M.name` node (one of the submodules), else
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
statement takes from one of the submodules; `Mod` itself is no use."""
function _testhelper_imported_names(ex)
    (length(ex.args) == 1 && ex.args[1] isa Expr && ex.args[1].head === :(:)) || return []
    path, items... = ex.args[1].args
    mod = path.args[end]
    mod in _testhelper_MODULES || return []
    names = [(item.head === :as ? item.args[1] : item).args[end] for item in items]
    [(mod, name) for name in names if name !== mod]
end

"""`(module, name, line)` for each `Mod.name` in the source text, `Mod` one of the
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

"""The names in the `public` statements of the files under `root/<mod>/`."""
function _testhelper_api(root, mod::Symbol)
    statements = Expr[]
    for file in _testhelper_source_files(joinpath(root, string(mod)))
        _testhelper_walk(Meta.parseall(read(file, String))) do ex, _
            ex.head === :public && push!(statements, ex)
        end
    end
    isempty(statements) && error("module $mod has no public statement")
    Set{Symbol}(n for ex in statements for n in ex.args)
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

"""Uses naming a module's name that is not public there (or naming a module with no
folder)."""
function _testhelper_api_violations(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    [u for u in _testhelper_uses(root) if !(u.name in get(apis, u.mod, Set{Symbol}()))]
end

"""Uses by a module of a sibling that its layer forbids; the top level is unrestricted."""
_testhelper_layer_violations(root) = [u for u in _testhelper_uses(root)
    if u.user !== :EnzymeRates && !(u.mod in _testhelper_LAYERS[u.user])]

"""For each module with a folder, the underscore names in its public statement or used from
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

"""The names in the `export` lines of `file`, in order."""
function _testhelper_exports(file)
    names = Symbol[]
    _testhelper_walk(Meta.parseall(read(file, String))) do ex, _
        ex.head === :export && append!(names, ex.args)
    end
    unique!(names)
end

# The diagram's styling header and the order of its boxes: the top level, then the
# submodules. An arrow runs from a box to a later one.
const _testhelper_MAP_HEADER = """
    %%{init: {"theme": "base", "flowchart": {"curve": "basis", "nodeSpacing": 40, "rankSpacing": 70, "htmlLabels": true},
      "themeVariables": {"fontFamily": "Helvetica, Arial, sans-serif", "fontSize": "13px",
      "primaryColor": "#ffffff", "primaryBorderColor": "#666666", "primaryTextColor": "#222222",
      "lineColor": "#777777", "edgeLabelBackground": "#ffffff", "clusterBkg": "#eaf2ee",
      "clusterBorder": "#666666", "titleColor": "#1f4f42"}}}%%
    """
const _testhelper_MAP_BOXES = (:EnzymeRates => "top", :DSL => "dsl", :Plotting => "plot",
    :Enumeration => "enum", :Fitting => "fit", :Derivation => "der", :Data => "data",
    :Mechanisms => "mech")

"""The names of the `struct` and `abstract type` definitions in `path`, a file or a folder
of files."""
function _testhelper_type_names(path)
    types = Set{Symbol}()
    typename(ex) = ex isa Symbol ? ex : typename(ex.args[1])
    for file in (isfile(path) ? [path] : _testhelper_source_files(path))
        _testhelper_walk(Meta.parseall(read(file, String))) do ex, _
            ex.head === :struct && push!(types, typename(ex.args[2]))
            ex.head === :abstract && push!(types, typename(ex.args[1]))
        end
    end
    types
end

"""The first sentence of the docstring of `module mod`, defined in `file`; a line break in
the docstring's sentence is a line break in the diagram's box."""
function _testhelper_box_sentence(file, mod::Symbol)
    docstring = nothing
    _testhelper_walk(Meta.parseall(read(file, String))) do ex, _
        ex.head === :macrocall || return
        target = ex.args[end]
        if target isa Expr && target.head === :module &&
           target.args[2] === mod && ex.args[end - 1] isa String
            docstring = ex.args[end - 1]
        end
    end
    docstring === nothing && error("$file has no docstring on `module $mod`")
    paragraph = first(split(strip(docstring), "\n\n"))
    sentence = match(r"^(.*?[.!?])(\s|$)"s, paragraph)[1]
    join(strip.(split(sentence, '\n')), "<br/>")
end

"""A name as the diagram writes it: a function as `name()`, a macro or a type bare, an
exported name in bold with a star. A non-breaking space (Mermaid's `#nbsp;`) joins the star
to the name, so a label that wraps never leaves the star on a line of its own."""
function _testhelper_map_name(name, types, exports)
    text = startswith(string(name), "@") || name in types ? string(name) : "$name()"
    name in exports ? "<b>★#nbsp;$text</b>" : text
end

"""The names in the order the diagram lists them: the exported ones first, as functions,
macros and types in the order of the `export` lines, then the others alphabetically."""
function _testhelper_map_order(names, types, exports)
    kind(n) = startswith(string(n), "@") ? 2 : n in types ? 3 : 1
    exported = sort!(collect(intersect(names, exports)); by = n -> (kind(n),
        findfirst(==(n), exports)))
    [exported; sort!(collect(setdiff(names, exports)); by = string)]
end

"""The Mermaid text of the module map of the tree at `root`: a box for the top level and
each module that has a folder, with the first sentence of its docstring, and an arrow for
each pair of modules where one takes names from the other, listing the names. The names
of `Mechanisms` types are listed in its box."""
function _testhelper_module_map(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    exports = _testhelper_exports(joinpath(root, "EnzymeRates.jl"))
    types = _testhelper_type_names(root)
    boxes = [(m, id) for (m, id) in _testhelper_MAP_BOXES
             if m === :EnzymeRates || haskey(apis, m)]
    ids = Dict(boxes)
    shown(names) = [_testhelper_map_name(n, types, exports)
        for n in _testhelper_map_order(names, types, exports)]
    lines = String["flowchart TB"]
    for (m, id) in boxes
        file = joinpath(root, m === :EnzymeRates ? "EnzymeRates.jl" : "$m/$m.jl")
        label = "<b><big>$m</big></b><br/>" * _testhelper_box_sentence(file, m)
        if m === :EnzymeRates
            defined = [n for n in exports if !any(api -> n in api, values(apis))]
            label *= "<br/><br/>Defines:<br/>" * join(shown(defined), "<br/>")
        elseif m === :Mechanisms
            label *= "<br/><br/>Types every module but Data reads:<br/>" *
                     join(shown([n for n in apis[m] if n in types]), "<br/>")
        end
        push!(lines, "  $id[\"$label\"]")
    end
    edges = Dict{Tuple{Symbol, Symbol}, Set{Symbol}}()
    for u in _testhelper_uses(root)
        (u.mod === :Mechanisms && u.name in types) ||
            push!(get!(edges, (u.user, u.mod), Set{Symbol}()), u.name)
    end
    order = first.(boxes)
    box_index((user, mod)) = (findfirst(==(user), order), findfirst(==(mod), order))
    for (user, mod) in sort!(collect(keys(edges)); by = box_index)
        label = join(shown(edges[(user, mod)]), "<br/>")
        push!(lines, "  $(ids[user]) -->|\"$label\"| $(ids[mod])")
    end
    push!(lines, "  classDef mod fill:#eaf2ee,stroke:#666666,color:#222222;",
        "  class $(join(last.(boxes), ',')) mod;")
    _testhelper_MAP_HEADER * join(lines, "\n") * "\n"
end

"""`names` as a table cell writes them: in the diagram's order, each in backticks, an
exported one after a ★."""
_testhelper_table_names(names, types, exports) = join([(n in exports ? "★ " : "") * "`$n`"
    for n in _testhelper_map_order(names, types, exports)], ", ")

"""The module table of the tree at `root`, as Markdown: a row for the top level and each
module that has a folder, with the first sentence of its docstring, its public names (the
top level's are the exports no module owns) and the modules its layer lets it call."""
function _testhelper_module_table(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    exports = _testhelper_exports(joinpath(root, "EnzymeRates.jl"))
    types = _testhelper_type_names(root)
    shown(names) = _testhelper_table_names(names, types, exports)
    callable(layer) = isempty(layer) ? "none" : join(layer, ", ")
    rows = ["| Module | What it does | Public names (★ exported) | May call |",
        "|:--|:--|:--|:--|"]
    for (m, _) in _testhelper_MAP_BOXES
        top = m === :EnzymeRates
        top || haskey(apis, m) || continue
        file = joinpath(root, top ? "EnzymeRates.jl" : "$m/$m.jl")
        sentence = replace(_testhelper_box_sentence(file, m), "<br/>" => " ")
        names = top ? [n for n in exports if !any(api -> n in api, values(apis))] : apis[m]
        push!(rows, "| $(top ? "$m (top level)" : m) | $sentence | $(shown(names)) | " *
                    "$(top ? "all" : callable(_testhelper_LAYERS[m])) |")
    end
    join(rows, "\n") * "\n"
end

"""The names that the top-level statement `ex` defines, pushed to `names`: a function's,
in long or short form; a macro's, as `Symbol("@name")`; a `const`'s; a type's; and an
`@enum`'s type and values; under a docstring or a macro such as `@inline` too. A method of
a qualified name (`Base.show`) or of a callable object defines no name. Not seen: a
definition inside a top-level `let`, `if`, `@static if` or loop (`for … @eval`), and the
names of `const A, B = …`."""
function _testhelper_defined_names!(names, ex)
    ex isa Expr || return names
    bare(x) = x isa Expr && x.head in (:where, :(::), :curly, :<:) ? bare(x.args[1]) : x
    add!(x) = x isa Symbol && push!(names, x)
    if ex.head === :macrocall && ex.args[1] === Symbol("@enum")
        for a in ex.args[3:end], x in (a isa Expr && a.head === :block ? a.args : [a])
            add!(bare(x))
        end
    elseif ex.head === :macrocall
        _testhelper_defined_names!(names, ex.args[end])
    elseif ex.head in (:function, :(=))
        sig = bare(ex.args[1])
        if sig isa Expr && sig.head === :call
            f = sig.args[1]
            add!(f isa Expr && f.head === :curly ? f.args[1] : f)
        elseif ex.head === :function
            add!(sig)  # `function f end`
        end
    elseif ex.head === :macro
        add!(Symbol("@", ex.args[1].args[1]))
    elseif ex.head === :const
        add!(bare(ex.args[1].args[1]))
    elseif ex.head === :struct
        add!(bare(ex.args[2]))
    elseif ex.head in (:abstract, :primitive)
        add!(bare(ex.args[1]))
    elseif ex.head in (:block, :toplevel)
        foreach(a -> _testhelper_defined_names!(names, a), ex.args)
    end
    names
end

"""The names the code of `ex` refers to, pushed to `found`. A quoted symbol refers to
nothing, though the code a quote interpolates does; nor do a field after a dot, a keyword
or named-tuple field name (`f(; a = 1)`, `(a = 1,)`) and the names of a `using`, `import`,
`export` or `public` statement. Strings, comments and docstrings hold no names, but a
string's interpolations do. A local variable or argument named like a sibling's definition
counts as a use of it, so rename the local. Code inside a quote, such as a macro's output
or a `@generated` body, is not seen."""
function _testhelper_referenced_names!(found, ex, quoted::Bool = false)
    if ex isa Symbol
        quoted || push!(found, ex)
    elseif ex isa Expr
        if ex.head === :$
            _testhelper_referenced_names!(found, ex.args[1])
        elseif quoted || ex.head === :quote
            foreach(a -> _testhelper_referenced_names!(found, a, true), ex.args)
        elseif ex.head === :kw
            # `k::T = d`, a keyword or an argument with a default, uses `T` and `d`.
            ex.args[1] isa Expr && ex.args[1].head === :(::) &&
                _testhelper_referenced_names!(found, ex.args[1].args[end])
            _testhelper_referenced_names!(found, ex.args[2])
        elseif ex.head === :tuple
            for a in ex.args
                field = a isa Expr && a.head === :(=) && a.args[1] isa Symbol
                _testhelper_referenced_names!(found, field ? a.args[2] : a)
            end
        elseif !(ex.head in (:using, :import, :export, :public))
            foreach(a -> _testhelper_referenced_names!(found, a), ex.args)
        end
    end
    found
end

"""Each source file of the tree at `root` but the module headers: the top level's, then
each module's, modules in the order of the module map and a module's files in the order
its header includes them (any it does not include last, by name). Each is
`(; mod, path, defines, takes, includes)`: its module (`:EnzymeRates` for the top level),
its path relative to `root`, the names it defines at top level, the names it takes from
each sibling file as `sibling path => names`, siblings in file order, and whether it calls
`include` (a call written `Base.include` is not seen). A file takes a name from a sibling
when the sibling defines it, the file does not, and the file's code refers to it
(`_testhelper_defined_names!`, `_testhelper_referenced_names!`)."""
function _testhelper_file_scans(root)
    scans = []
    for (mod, _) in _testhelper_MAP_BOXES
        top = mod === :EnzymeRates
        dir = top ? root : joinpath(root, string(mod))
        isdir(dir) || continue
        header = joinpath(dir, "$mod.jl")
        order = String[]
        _testhelper_walk(Meta.parseall(read(header, String))) do ex, _
            ex.head === :call && ex.args[1] === :include && ex.args[end] isa String &&
                push!(order, joinpath(dir, ex.args[end]))
        end
        files = top ? filter(endswith(".jl"), readdir(root; join = true)) :
                _testhelper_source_files(dir)
        sort!(filter!(!=(header), files);
            by = f -> (something(findfirst(==(f), order), length(order) + 1), f))
        parsed = [Meta.parseall(read(f, String)) for f in files]
        defines = [_testhelper_defined_names!(Set{Symbol}(), ex) for ex in parsed]
        for (i, file) in enumerate(files)
            uses = setdiff(_testhelper_referenced_names!(Set{Symbol}(), parsed[i]),
                defines[i])
            takes = [_testhelper_relpath(root, files[j]) =>
                         sort!(collect(intersect(uses, defines[j])))
                     for j in eachindex(files) if j != i]
            push!(scans, (; mod, path = _testhelper_relpath(root, file),
                defines = defines[i], takes = filter!(t -> !isempty(t.second), takes),
                includes = :include in uses))
        end
    end
    scans
end

# File names that say nothing of what a file gives the package.
const _testhelper_VAGUE_STEMS = ("utils", "util", "helpers", "helper", "common", "misc",
    "types", "core", "internal", "internals", "functions")

"""Where the source files of the tree at `root`, module headers aside, break the naming
and include rules, one line each: a file that calls `include`, which only a module's
header does; a file whose name says nothing of what it holds (`_testhelper_VAGUE_STEMS`);
a file named for a public name of its module (an export, at the top level) that it does
not define; a file `<x>_type.jl` that defines no type whose name is `<x>` in snake case;
and a file `<x>_macro.jl` that defines no macro `@<x>`. Any other name names the
calculation the file carries out and is not checked."""
function _testhelper_file_name_and_include_problems(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    apis[:EnzymeRates] = Set(_testhelper_exports(joinpath(root, "EnzymeRates.jl")))
    snake(name) = lowercase(replace(string(name), r"(?<=[a-z0-9])(?=[A-Z])" => "_"))
    problems = String[]
    for f in _testhelper_file_scans(root)
        stem = chopsuffix(basename(f.path), ".jl")
        type_name, macro_name = chopsuffix(stem, "_type"), chopsuffix(stem, "_macro")
        f.includes && push!(problems,
            "$(f.path) calls `include`; only a module's header includes files")
        stem in _testhelper_VAGUE_STEMS &&
            push!(problems, "$(f.path) has a name that says nothing of what it holds")
        Symbol(stem) in apis[f.mod] && !(Symbol(stem) in f.defines) &&
            push!(problems, "$(f.path) is named for `$stem`, a public name of $(f.mod), " *
                            "but does not define it")
        type_name == stem || any(n -> snake(n) == type_name,
                                 _testhelper_type_names(joinpath(root, f.path))) ||
            push!(problems, "$(f.path) ends in `_type` but defines no type whose name " *
                            "is `$type_name` in snake case")
        macro_name == stem || Symbol("@", macro_name) in f.defines ||
            push!(problems, "$(f.path) ends in `_macro` but defines no macro " *
                            "`@$macro_name`")
    end
    problems
end

"""The shortest route from file `from` to file `to` along `uses`, which maps each file to
the files it takes names from: the files in order, `from` first and `to` last, at least one
use long, so a loop when `from == to`; `nothing` when there is none."""
function _testhelper_shortest_route(uses, from, to)
    parent = Dict{String, String}()
    queue = [from]
    while !isempty(queue)
        file = popfirst!(queue)
        for next in uses[file]
            haskey(parent, next) && continue
            parent[next] = file
            if next == to
                route = [to]
                while length(route) == 1 || route[1] != from
                    pushfirst!(route, parent[route[1]])
                end
                return route
            end
            push!(queue, next)
        end
    end
    nothing
end

"""One line for each loop of files within a module of the tree at `root`, each file taking
names from the next and the last from the first, once the uses that `exempt` holds are set
aside: for each group of files that reach one another, the shortest loop through its first
file, naming the files in order and the names each takes from the next. And one line for
each use that `exempt` holds and that no longer closes a loop. A use is
`(file, sibling)`, two paths relative to `root`: `file` takes names from `sibling`."""
function _testhelper_file_loop_problems(root, exempt)
    scans = _testhelper_file_scans(root)
    takes = Dict((f.path, from) => names for f in scans for (from, names) in f.takes)
    uses = Dict(f.path => first.(f.takes) for f in scans)
    kept = Dict(file => filter(s -> !((file, s) in exempt), siblings)
        for (file, siblings) in uses)
    reaches(graph, a, b) = _testhelper_shortest_route(graph, a, b) !== nothing
    shown(names) = join(("`$n`" for n in names), ", ")
    problems, looped = String[], Set{String}()
    for file in (f.path for f in scans)
        loop = file in looped ? nothing : _testhelper_shortest_route(kept, file, file)
        loop === nothing && continue
        group = (g for g in keys(kept) if reaches(kept, file, g) && reaches(kept, g, file))
        union!(looped, group)
        push!(problems, "a loop of files, each taking names from the next: " *
            join(["$a takes $(shown(takes[(a, b)])) from $b"
                  for (a, b) in zip(loop, loop[2:end])], "; "))
    end
    for (file, sibling) in sort!(collect(exempt))
        sibling in get(uses, file, String[]) && reaches(uses, sibling, file) ||
            push!(problems, "$file taking names from $sibling no longer closes a loop; " *
                            "remove the use from the exemption set in " *
                            "test/test_boundaries.jl")
    end
    problems
end

"""The two `# ABOUTME: ` lines that `file` starts with, joined without the prefix, or
`nothing` unless it starts with exactly two such lines."""
function _testhelper_aboutme(file)
    prefix = "# ABOUTME: "
    lines = first(readlines(file), 3)
    count(startswith(prefix), lines) == 2 && all(startswith(prefix), lines[1:2]) ||
        return nothing
    join(chopprefix.(lines[1:2], prefix), " ")
end

"""One line for each `.jl` file under `src/`, `ext/`, `test/` and `docs/` of the
repository at `root` that does not start with exactly two `# ABOUTME: ` lines. The built
site `docs/build` and the local results `docs/hpc_results`, which git ignores, are
skipped."""
function _testhelper_aboutme_problems(root)
    skipped = (joinpath(root, "docs", "build"), joinpath(root, "docs", "hpc_results"))
    dirs = filter(isdir, [joinpath(root, d) for d in ("src", "ext", "test", "docs")])
    problems = String[]
    while !isempty(dirs)
        for path in readdir(pop!(dirs); join = true)
            if isdir(path)
                path in skipped || push!(dirs, path)
            elseif endswith(path, ".jl") && _testhelper_aboutme(path) === nothing
                push!(problems, "$(_testhelper_relpath(root, path)) does not start " *
                                "with exactly two `# ABOUTME: ` lines")
            end
        end
    end
    sort!(problems)
end

"""The table of the source files of the tree at `root`, as Markdown: a row for each file
in the order of `_testhelper_file_scans`, with its module, its name, its two ABOUTME lines
joined, the public names of its module that it defines (★ exported) and each sibling file
it takes names from, with the number of names."""
function _testhelper_file_table(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    exports = _testhelper_exports(joinpath(root, "EnzymeRates.jl"))
    apis[:EnzymeRates] = Set(exports)
    types = _testhelper_type_names(root)
    stem(path) = chopsuffix(basename(path), ".jl")
    rows = ["| Module | File | What it holds | Public names it defines | " *
            "Takes names from |", "|:--|:--|:--|:--|:--|"]
    for f in _testhelper_file_scans(root)
        # Documenter's Markdown reads an underscore inside a word as emphasis.
        about = replace(something(_testhelper_aboutme(joinpath(root, f.path)), ""),
            "_" => "\\_")
        takes = join(["`$(stem(from))` ($(length(names)))" for (from, names) in f.takes],
            ", ")
        cells = [f.mod === :EnzymeRates ? "top level" : string(f.mod),
            "`$(basename(f.path))`", about,
            _testhelper_table_names(intersect(apis[f.mod], f.defines), types, exports),
            takes]
        push!(rows, "| " * join(cells, " | ") * " |")
    end
    join(rows, "\n") * "\n"
end

# The three blocks generated from src/: the lines that open and close each, the function
# that generates its text from a source tree, and how a problem line describes it.
const _testhelper_MAP_BLOCK = (; opening = "```mermaid\n", closing = "```\n",
    generate = _testhelper_module_map,
    described = "the Mermaid module map in docs/src/developer.md")
const _testhelper_TABLE_BLOCK = (;
    opening = "<!-- module table, generated from src/ by test/update_module_map.jl -->\n",
    closing = "<!-- end of module table -->\n", generate = _testhelper_module_table,
    described = "the module table in .claude/CLAUDE.md")
# Documenter shows an HTML comment as text unless a raw block holds it.
const _testhelper_FILE_TABLE_BLOCK = (;
    opening = "```@raw html\n<!-- file table, generated from src/ by " *
              "test/update_module_map.jl -->\n```\n\n",
    closing = "\n```@raw html\n<!-- end of file table -->\n```\n",
    generate = _testhelper_file_table,
    described = "the file table in docs/src/developer.md")

"""The range of `doc` between `block`'s opening and closing lines, or `nothing` when `doc`
lacks the block."""
function _testhelper_generated_block_range(doc, block)
    opening = findfirst(block.opening, doc)
    opening === nothing && return nothing
    closing = findnext(block.closing, doc, last(opening) + 1)
    closing === nothing ? nothing : (last(opening) + 1):prevind(doc, first(closing))
end

"""One line when the generated `block` of the Markdown text `doc` is missing or differs
from the one generated from `root`, naming the command that rewrites it. CRLF line endings
(a Windows checkout) compare as LF."""
function _testhelper_generated_block_problems(doc, block, root)
    doc = replace(doc, "\r\n" => "\n")
    range = _testhelper_generated_block_range(doc, block)
    range !== nothing && doc[range] == block.generate(root) && return String[]
    ["$(block.described) differs from the one generated from src/; rewrite it with " *
     "`julia test/update_module_map.jl`"]
end

"""`doc`, with LF line endings, and the text of its generated `block` replaced by `text`.
Throws when `doc` lacks the block's opening or closing line."""
function _testhelper_replace_generated_block(doc, block, text)
    doc = replace(doc, "\r\n" => "\n")
    range = _testhelper_generated_block_range(doc, block)
    range === nothing && error("$(block.described) has no generated block: add the " *
                               "lines $(repr(block.opening)) and $(repr(block.closing))")
    doc[1:prevind(doc, first(range))] * text * doc[nextind(doc, last(range)):end]
end
