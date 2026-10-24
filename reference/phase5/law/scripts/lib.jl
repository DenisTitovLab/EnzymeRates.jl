# ABOUTME: Shared helpers of the law-type study: census and hand mechanisms as slice-variant
# ABOUTME: mechanisms, laws as BigInt dictionaries, faces, class 5 and the RE-turnover check.

const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const R5 = joinpath(SCR, "phase45", "research5", "law")
using Printf, Random, Serialization, Statistics

module NetText
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research/oracle/net.jl")
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research/oracle/cases.jl")
end

"""Load law-type variant `v` (`"U"`, `"A"`, `"B"`, `"C"`) and return its module."""
function load_variant(v)
    Base.include(Main, joinpath(R5, "variants", v, "Slice$v.jl"))
    Base.invokelatest(getfield, Main, Symbol("Slice", v))
end

"""Roles `S=a,b;P=c;M=x` as `(substrates, products, modifiers)` name lists."""
function parse_roles(r::AbstractString)
    d = Dict{String, Vector{String}}()
    for kv in split(r, ';')
        k, v = split(kv, '=')
        d[String(k)] = String.(filter(!isempty, split(v, ',')))
    end
    d["S"], d["P"], d["M"]
end

"""The study's census sample: `(set, id, net)` rows of `data/sample.tsv`
(`scripts/extract_sample.py`)."""
function census_rows(path = joinpath(R5, "data", "sample.tsv"))
    out = NamedTuple[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        set, id, roles, regs, text = split(line, '\t')
        subs, prods, mods = parse_roles(roles)
        mods = [x for x in mods if !(x in subs || x in prods)]
        rg = String.(filter(!isempty, split(regs, ',')))
        net = NetText.census_net(text, (subs, prods, mods, rg); name = String(id))
        push!(out, (set = String(set), id = String(id), net = net))
    end
    out
end

"""Phase 4's hand cases (`oracle/cases.jl`) as `(set, id, net)` rows."""
hand_rows() = [(set = "hand", id = n.name, net = n)
               for n in [NetText.GUARD_CASES; NetText.ADMISSION_CASES; NetText.PINGPONG_CASES]]

"""The mechanism of `net` in slice-variant module `S`: `Xmod` is reactant X bound as a
modifier, a label that is no reactant is a modifier, a residual `plusA_minusP` adds A and
subtracts P."""
function slice_mech(S::Module, net)
    reactant(x) = x in net.subs || x in net.prods
    function met(l)
        endswith(l, "mod") && reactant(l[1:(end - 3)]) &&
            return S.Met(Symbol(l[1:(end - 3)]), true)
        S.Met(Symbol(l), !reactant(l))
    end
    function form(f)
        added, subtracted = Symbol[], Symbol[]
        for t in (isempty(f.res) ? String[] : split(f.res, '_'))
            startswith(t, "plus") ? push!(added, Symbol(t[5:end])) :
            startswith(t, "minus") ? push!(subtracted, Symbol(t[6:end])) :
            error("residual text $(f.res)")
        end
        S.RawForm(Symbol(f.conf), [met(x) for x in f.bound], added, subtracted)
    end
    steps = [S.RawStep(form(net.forms[s.a]), form(net.forms[s.b]), S.Met[met(x) for x in s.cons],
                       S.Met[met(x) for x in s.rel], s.re) for s in net.steps]
    group = [s.group for s in net.steps]
    isempty(net.tags) && return S.Mechanism(Symbol.(net.subs), Symbol.(net.prods), steps, group)
    sites = S.Site[S.Site([Symbol(l) for (l, _) in s.ligands], [t for (_, t) in s.ligands], s.mult)
             for s in net.sites]
    S.Mechanism(Symbol.(net.subs), Symbol.(net.prods), steps, group; tags = copy(net.tags),
                multiplicity = net.n, sites)
end

"""Census group of a sample set name: `spec`, or the reaction (`bibi`, `pingpong`, …)."""
group(set) = startswith(set, "spec:") ? "spec" : first(split(set, ':'))

# ─── Laws as plain data ─────────────────────────────────────────────────────────────────

const BMono = Vector{Pair{Int, Rational{Int}}}
"""A polynomial of any variant as `Dict{monomial, BigInt}`."""
bdict(p) = Dict{BMono, BigInt}(m => BigInt(c) for (m, c) in (p isa AbstractDict ? p : p.terms))

"""Every part of a law of any variant as comparable plain data."""
function lawdata(L)
    if hasproperty(L, :active)
        I = L.inactive === nothing ? nothing : lawdata(L.inactive)
        return (:allo, L.vars.names, bdict(L.H), lawdata(L.active), I, bdict(L.fA),
                bdict(L.fI), [(bdict(a), bdict(b), k) for (a, b, k) in L.sites], L.n)
    end
    (:plain, L.vars.names, L.vars.conc, bdict(L.H), bdict(L.Q), bdict(L.D),
     map(bdict, L.sigma), map(bdict, L.trees), L.shift)
end

"""Largest absolute coefficient of a polynomial of any variant."""
maxcoef(p) = maximum((abs(BigInt(c)) for c in values(p isa AbstractDict ? p : p.terms));
                     init = big(0))

# ─── Faces, class 5 and the RE-turnover check ───────────────────────────────────────────

"""Phase 4's faces (phase4_research.md section 1.2) of slice mechanism `m` as `(label,
zeroed names)`: the interior, every product zero, every substrate zero, each reactant alone
when its side has two or more; each again with every modifier-only name zero. A `::Mod`
copy is zeroed with its reactant (same concentration variable); site ligands never."""
function faces(m)
    subs, prods = m.substrates, m.products
    mods = sort!(unique!([x.name for x in m.mets if x.modifier && !(x.name in subs || x.name in prods)]))
    fs = Tuple{String, Vector{Symbol}}[("interior", Symbol[]), ("products=0", prods),
                                       ("substrates=0", subs)]
    for side in (subs, prods)
        length(side) > 1 && foreach(x -> push!(fs, ("$x=0", [x])), side)
    end
    isempty(mods) && return fs
    [fs; [(l == "interior" ? "modifiers=0" : l * ",modifiers=0", [z; mods]) for (l, z) in fs]]
end

"""Whether polynomial `p` holds a term free of every variable in `zs`."""
avoids(p, zs) = any(mo -> all(z -> all(pr -> pr.first != z || pr.second == 0, mo), zs),
                    keys(p isa AbstractDict ? p : p.terms))

"""
    class5(m, L) -> label of the first face whose denominator vanishes, or ""

Phase 4's support test (phase4_research.md section 1.6) on the new law: on every face some
term of the denominator survives. Plain: `L.D` in lowest terms. MWC: the denominator is
`(fI·DA)^n·WA + L·(fA·DI)^n·WI` with positive coefficients and sites with a constant term,
so it vanishes exactly when neither `fI·DA` nor `fA·DI` keeps a term, and a product keeps a
term exactly when each factor does.
"""
function class5(m, L)
    vars = L.vars
    for (label, names) in faces(m)
        zs = [vars.index[x] for x in names if haskey(vars.index, x)]
        ok = if hasproperty(L, :active)
            DI = L.inactive === nothing ? nothing : L.inactive.D
            (avoids(L.active.D, zs) && avoids(L.fI, zs)) ||
                ((DI === nothing || avoids(DI, zs)) && avoids(L.fA, zs))
        else
            avoids(L.D, zs)
        end
        ok || return label
    end
    ""
end

"""
    re_turnover_edge(S, m, mask) -> edge or 0

The RE-turnover check (today's `_compute_alpha` check, RD:169-180; phase 4's class 2) on
the slice's content-ordered forest: each rapid-equilibrium edge of `mask` must reproduce,
per metabolite name (a copy counts as its reactant), the difference its forms' potentials
have along the forest. The first edge that does not closes a cycle of rapid-equilibrium
steps that takes up some name net; 0 when every edge agrees.
"""
function re_turnover_edge(S::Module, m, mask)
    segment, parent, order = S.re_forest(m, mask)
    names = unique!([x.name for x in m.mets])
    row = [findfirst(==(x.name), names) for x in m.mets]
    ex = zeros(Int, length(names), length(m.forms))
    step(e, r) = sum((S.count_at(e.cons, j) - S.count_at(e.rel, j) for j in eachindex(m.mets)
                      if row[j] == r); init = 0)
    for i in order
        k = parent[i]
        k == 0 && continue
        e = m.edges[abs(k)]
        for r in eachindex(names)
            ex[r, i] = ex[r, k > 0 ? e.a : e.b] + sign(k) * step(e, r)
        end
    end
    for (k, e) in enumerate(m.edges)
        mask[k] && e.re || continue
        all(r -> ex[r, e.b] - ex[r, e.a] == step(e, r), eachindex(names)) || return k
    end
    0
end

"""`f()`'s best time over `reps` runs after one warm run, in µs."""
best(f; reps = 5) = (f(); minimum(@elapsed(f()) for _ in 1:reps) * 1e6)

"""Quantile `p` of `v` (nearest rank)."""
q(v, p) = isempty(v) ? NaN : sort(v)[clamp(round(Int, p * length(v)), 1, length(v))]

"""`v`'s median, p90 and max as text."""
stats(v) = @sprintf("n %4d median %8.1f p90 %8.1f max %9.1f", length(v), q(v, 0.5), q(v, 0.9),
                    isempty(v) ? NaN : maximum(v))
