# ABOUTME: Parses src/ for the module boundary tests (cross-module uses, public names,
# ABOUTME: exports) and generates the module map and the module table from them.

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
const _testhelper_MAP_BOXES = (:EnzymeRates => "top", :DSL => "dsl",
    :Enumeration => "enum", :Fitting => "fit", :Derivation => "der", :Data => "data",
    :Mechanisms => "mech")

"""The names of the `struct` and `abstract type` definitions under `root`."""
function _testhelper_type_names(root)
    types = Set{Symbol}()
    typename(ex) = ex isa Symbol ? ex : typename(ex.args[1])
    for file in _testhelper_source_files(root)
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
        target = ex.args[end]
        if ex.head === :macrocall && target isa Expr && target.head === :module &&
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

"""The module table of the tree at `root`, as Markdown: a row for the top level and each
module that has a folder, with the first sentence of its docstring, its public names (the
top level's are the exports no module owns) and the modules its layer lets it call."""
function _testhelper_module_table(root)
    apis = Dict(m => _testhelper_api(root, m)
        for m in keys(_testhelper_module_sources(root)))
    exports = _testhelper_exports(joinpath(root, "EnzymeRates.jl"))
    types = _testhelper_type_names(root)
    shown(names) = join([(n in exports ? "★ " : "") * "`$n`"
        for n in _testhelper_map_order(names, types, exports)], ", ")
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

# The two blocks generated from src/: the lines that open and close each, the function that
# generates its text from a source tree, and how a problem line describes it.
const _testhelper_MAP_BLOCK = (; opening = "```mermaid\n", closing = "```\n",
    generate = _testhelper_module_map,
    described = "the Mermaid module map in docs/src/developer.md")
const _testhelper_TABLE_BLOCK = (;
    opening = "<!-- module table, generated from src/ by test/update_module_map.jl -->\n",
    closing = "<!-- end of module table -->\n", generate = _testhelper_module_table,
    described = "the module table in .claude/CLAUDE.md")

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
