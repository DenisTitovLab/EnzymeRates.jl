# ABOUTME: Slice-side helpers for the key/codegen study: census text to slice mechanisms, the key
# ABOUTME: with concentrations labelled by name and search counters, law sizes and timed codegen.
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
isdefined(Main, :Slice) || include(joinpath(SP, "slice", "src", "Slice.jl"))
const S = Slice
using Printf, Random, SHA, Serialization, Statistics

# ─── Census neutral text → slice mechanism (from research5/unit/ulib.jl:14-89, read
# and checked; the grammar is phase3/census/README.md "Neutral text") ───────────────────

function parse_label(lab::AbstractString, reactants)
    if endswith(lab, "mod") && Symbol(lab[1:(end - 3)]) in reactants
        return S.Met(Symbol(lab[1:(end - 3)]), true)
    end
    x = Symbol(lab)
    S.Met(x, !(x in reactants))
end

function parse_form(t::AbstractString, reactants)
    i = findfirst('(', t)
    conf, inner = Symbol(t[1:(i - 1)]), t[(i + 1):(end - 1)]
    res = ""
    occursin("; residual=", inner) && ((inner, res) = split(inner, "; residual="))
    bound = isempty(inner) ? S.Met[] : [parse_label(l, reactants) for l in split(inner, ",")]
    added, subtracted = Symbol[], Symbol[]
    for r in (isempty(res) ? SubString{String}[] : split(res, "_"))
        startswith(r, "plus") ? push!(added, Symbol(r[5:end])) :
            push!(subtracted, Symbol(r[6:end]))
    end
    S.RawForm(conf, bound, added, subtracted)
end

function parse_side(t::AbstractString, reactants)
    parts = split(t, " + ")
    parse_form(parts[1], reactants), S.Met[parse_label(l, reactants) for l in parts[2:end]]
end

"""The slice mechanism of a census neutral text over the given substrates and products."""
function parse_neutral(text::AbstractString, subs::Vector{Symbol}, prods::Vector{Symbol})
    reactants = Set([subs; prods])
    allo = startswith(text, "allosteric")
    mult, body, sitetext = 1, "", "none"
    if allo
        rest = text[length("allosteric n=") + 1:end]
        sp = findfirst(' ', rest)
        mult = parse(Int, rest[1:(sp - 1)])
        body, sitetext = split(rest[(sp + 1):end], " sites "; limit = 2)
    else
        body = text[length("plain ") + 1:end]
    end
    steps, group, tags = S.RawStep[], Int[], Symbol[]
    for (g, gt) in enumerate(split(body, " / "))
        tag = ""
        if !endswith(gt, "}")
            gt, tag = split(gt, "}:")
            gt *= "}"
        end
        allo && push!(tags, Symbol(tag))
        for st in split(gt[2:(end - 1)], " | ")
            flag, core = st[(end - 1):end], st[1:(end - 3)]
            a, b = split(core, " <-> ")
            (fa, ca), (fb, cb) = parse_side(a, reactants), parse_side(b, reactants)
            push!(steps, S.RawStep(fa, fb, ca, cb, flag == "RE"))
            push!(group, g)
        end
    end
    sites = S.Site[]
    if allo && sitetext != "none"
        for s in split(sitetext, " / ")
            mm = match(r"^site\((\d+): (.*)\)$", s)
            ligs = [split(l, "::") for l in split(mm[2], ", ")]
            push!(sites, S.Site(Symbol.(first.(ligs)), Symbol.(last.(ligs)),
                                parse(Int, mm[1])))
        end
    end
    allo ? S.Mechanism(subs, prods, steps, group; tags, multiplicity = mult, sites) :
           S.Mechanism(subs, prods, steps, group)
end

# ─── TSV rows (census_common.jl escapes) ────────────────────────────────────────────────

function unesc(s::AbstractString)
    io, i, cs = IOBuffer(), 1, collect(s)
    while i <= length(cs)
        c = cs[i]
        if c == '\\' && i < length(cs)
            n = cs[i + 1]
            print(io, n == 't' ? '\t' : n == 'n' ? '\n' : n == 'r' ? '\r' : n)
            i += 2
        else
            print(io, c)
            i += 1
        end
    end
    String(take!(io))
end
tsv(path) = [[unesc(f) for f in split(l, '\t')] for l in eachline(path) if !isempty(l)]
names_of(s) = isempty(s) ? Symbol[] : Symbol.(split(s, ","))

# ─── The key: slice/src/key.jl with concentrations labelled by name (spec line 537) ─────

"""
`canonical` of slice/src/key.jl:75-118 with two changes: a concentration is labelled
`c_<name>` and a term's fixed factors are listed by label, so the text depends on names,
not positions; and the search reports its leaves and nodes. `cap` is the leaf cap.
Returns `(canonical, leaves, nodes)`.
"""
function canonical_named(L; cap::Int = 100_000)
    comps, units, header = S.key_parts(L)
    vars = L.vars
    fixed(v) = vars.conc[v] || vars.names[v] in (:Keq, :L)
    terms = Tuple{Int, S.Mono, BigInt}[]
    for u in units
        s = S.lowest_shift(comps[u])
        g = foldl(gcd, (c for i in u for c in values(comps[i].terms)); init = big(0))
        append!(terms, ((i, S.mmul(mo, s), div(big(c), g)) for i in u
                        for (mo, c) in comps[i].terms))
    end
    params = [v for v in eachindex(vars.names) if !fixed(v)]
    held = Set(v for (_, mo, _) in terms for (v, _) in mo)
    present = filter(in(held), params)
    at = Dict(v => j for (j, v) in enumerate(present))
    name(v) = vars.conc[v] ? "c_" * string(vars.names[v]) : string(vars.names[v])
    labels = [string(i, "|", c, "|",
                     join(sort!(["$(name(v))^$e" for (v, e) in mo if fixed(v)]), "*"))
              for (i, mo, c) in terms]
    inc = [[(at[v], e) for (v, e) in mo if !fixed(v)] for (_, mo, _) in terms]
    occ = [Tuple{Int, S.Exp}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    best, order, leaves, nodes = nothing, Int[], 0, 0
    function search(pcol, tcol)
        nodes += 1
        pcol, tcol = S.refine(pcol, tcol, inc, occ)
        if allunique(pcol)
            (leaves += 1) > cap && error("canonical labelling: over $cap leaves")
            text = S.encode(pcol, inc, labels)
            (best === nothing || text < best) && ((best, order) = (text, sortperm(pcol)))
            return
        end
        c = minimum(x for x in pcol if count(==(x), pcol) > 1)
        for p in findall(==(c), pcol)
            q = 2 .* pcol
            q[p] -= 1
            search(q, tcol)
        end
    end
    search(ones(Int, length(present)), S.ranks(labels))
    absent = sort!(vars.names[setdiff(params, present)])
    text = string(header, "|P=", length(params), "|absent=", length(absent), "\n", best)
    S.Canonical([vars.names[present[order]]; absent], text, bytes2hex(sha256(text))), leaves,
        nodes
end

# ─── Law sizes ──────────────────────────────────────────────────────────────────────────

"""Terms of the polynomials the key reads (every unit of `key_parts`)."""
key_terms(L) = sum(S.nterms, S.key_parts(L)[1])

"""Denominator terms after substitution, both conformations for an MWC law (a flip-form-only
inactive conformation counts 1)."""
den_terms(L::S.Law) = S.nterms(L.D)
den_terms(A::S.AlloLaw) = S.nterms(A.active.D) + S.nterms(S.inactive_terms(A)[2])

"""V×τ of the active graph and of the inactive one (0 for a plain mechanism)."""
cx_pair(m) = (Int(S.complexity(m)),
              S.is_allosteric(m) ? (im = S.inactive_mask(m); any(im) ? Int(S.complexity(m, im)) : 0) : 0)

"""The calls a generated body makes (`+`, `-`, `*`, `/`, `inv`, `sqrt`, `^`), each binary."""
function body_ops(ex)
    n = 0
    walk(e) = e isa Expr && (e.head == :call && e.args[1] in (:+, :-, :*, :/, :inv, :sqrt, :^) &&
                             (n += max(1, length(e.args) - 2)); foreach(walk, e.args))
    walk(ex)
    n
end

# ─── Timing the generated rate_equation ────────────────────────────────────────────────

@noinline call_barrier(em, c, p) = S.rate_equation(em, c, p)

"""Random concentrations (log-uniform 0.1-10) and parameters (log-uniform 0.1-10, Keq 10,
E_total 1) of `m` as NamedTuples keyed by name."""
function random_inputs(m, rng)
    lu() = 10.0^(2 * rand(rng) - 1)
    concs = (; (x => lu() for x in S.metabolites(m))...)
    ps = S.parameters(m)
    params = (; (n => (n === :Keq ? 10.0 : n === :E_total ? 1.0 : lu()) for n in ps)...)
    concs, params
end

"""Seconds per call: the minimum over 5 batches of `n` calls through a `@noinline` barrier
(`barrier = true`) or inlined into the loop as the package's gate helper does; and the
bytes the second call allocates."""
function time_call(em, concs, params; n = 10_000, barrier = true)
    call_barrier(em, concs, params)
    allocs = @allocated call_barrier(em, concs, params)
    best, acc = barrier ? time_loop(call_barrier, em, concs, params, n) :
                          time_loop(S.rate_equation, em, concs, params, n)
    isfinite(acc) || error("non-finite rate")
    best, allocs
end

function time_loop(f::F, em, concs, params, n) where {F}
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:n
            acc += f(em, concs, params)
        end
        best = min(best, t / n)
    end
    best, acc
end
