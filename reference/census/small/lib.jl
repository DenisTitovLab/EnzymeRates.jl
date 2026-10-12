# ABOUTME: Census library: mirrors identify_rate_equation's candidate pipeline (seeds, base tier,
# ABOUTME: expansion, dedup, filters) and classifies candidates (guard, key, rank, saturation).
using EnzymeRates, Random, Serialization, Printf, LinearAlgebra, Statistics, SHA
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const CEN = joinpath(SCR, "census", "small")
include(joinpath(SCR, "round4", "guard", "lib.jl"))          # Pop, Slice as S, face_check, ...
include(joinpath(SCR, "proto", "common", "export", "export_lib.jl"))   # ER, _pspecies, _family

# ─── reactions ──────────────────────────────────────────────────────────────────────────

const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")
const RXN_NAMES = (:uni_uni_rxn, :uni_uni_with_reg, :uni_uni_allo, :uni_uni_allo_reg,
                   :uni_uni_allo_2reg, :uni_uni_reg_and_inhibitor, :uni_bi_rxn, :bi_bi_rxn,
                   :bi_bi_pp_rxn)
let ast = Meta.parseall(read(joinpath(TESTDIR, "test_mechanism_enumeration.jl"), String))
    for ex in ast.args
        ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
            ex.args[1].args[1] in RXN_NAMES && Core.eval(Main, ex)
    end
end

# ─── identify's pipeline, without the fits ──────────────────────────────────────────────

const CX_MAX = 337          # identify's eq_complexity_filter default
const NP_MAX = 20           # identify's max_param_count default

"""The seeds `_beam_search` starts from: `init_mechanisms`, or `seed_mechanisms` when the
reaction declares a regulator (every declared regulator required, the default)."""
function identify_seeds(rxn)
    req_allo, req_comp = ER._required_regulators(rxn, Symbol[], Symbol[])
    isempty(req_allo) && isempty(req_comp) ?
        (unique!(collect(ER.init_mechanisms(rxn))), "init_mechanisms") :
        (unique!(collect(ER.seed_mechanisms(rxn, req_allo, req_comp))),
         "seed_mechanisms(required allo $(sort!(collect(req_allo))), comp " *
         "$(sort!(collect(req_comp))))")
end

"""`expand_mechanisms([m], rxn)`, child by child with the move that made it, in the same
order; the same regulator-type filter and atom assertion."""
function expand_labeled(m, rxn)
    ER._assert_emission_rules(m)
    out = Tuple{Any, Symbol}[]
    add!(cs, mv) = foreach(c -> push!(out, (c, mv)), cs)
    add!(ER._expand_re_to_ss(m), :re_to_ss)
    add!(ER._expand_split_kinetic_group(m), :split)
    add!(ER._expand_add_dead_end_regulator(m, rxn), :dead_end)
    if m isa ER.Mechanism
        add!(ER._expand_to_allosteric(m, rxn), :to_allosteric)
    else
        add!(ER._expand_add_allosteric_regulator(m, rxn), :add_allo_reg)
        add!(ER._expand_change_allo_state(m), :change_allo_state)
        add!(ER._expand_merge_regulatory_sites(m), :merge_sites)
    end
    filter!(cm -> ER._respects_reg_type(cm[1], rxn), out)
    foreach(cm -> ER._assert_atom_conserving(cm[1]), out)
    out
end

"""The package mechanism as DSL-like text: kinetic groups in parentheses, tags for an
allosteric mechanism."""
function steps_text(m)
    gs = ER.steps(m)
    tags = m isa ER.AllostericMechanism ? collect(ER.cat_allo_states(m)) : nothing
    parts = String[]
    for (g, grp) in enumerate(gs)
        t = join((ER._step_text(s) for s in grp), ", ")
        length(grp) > 1 && (t = "(" * t * ")")
        tags === nothing || (t *= " ::" * string(tags[g]))
        push!(parts, t)
    end
    s = join(parts, "; ")
    if m isa ER.AllostericMechanism
        s *= " | n=$(ER.catalytic_multiplicity(m))"
        for site in ER.regulatory_sites(m)
            s *= " | site(" * join(("$(ER.name(l))::$t" for (l, t) in
                                   zip(ER.ligands(site), ER.allo_states(site))), ", ") *
                 "; n=$(ER.multiplicity(site)))"
        end
    end
    s
end

"""The structure of a package mechanism, plain or allosteric, as a `Pop.PMech` (no package
outputs), as `slice/integrate/export_terter_expand.jl` exports it."""
function structure(m, id, source)
    allo = m isa ER.AllostericMechanism
    groups = ER.steps(m)
    flat = [s for g in groups for s in g]
    species = ER._forms(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    rxn = ER.reaction(m)
    subs = Symbol[ER.name(x) for x in ER.substrates(rxn)]
    prods = Symbol[ER.name(x) for x in ER.products(rxn)]
    off = cumsum([0; length.(groups)])
    tags = allo ? collect(Symbol, ER.cat_allo_states(m)) : Symbol[]
    sites = allo ? ER.regulatory_sites(m) : ER.RegulatorySite[]
    Pop.PMech(; id, source, substrates = subs, products = prods,
              metabolites = [subs; prods], allosteric = allo, group_tags = tags,
              cat_multiplicity = allo ? ER.catalytic_multiplicity(m) : 1,
              sites = Pop.PSite[Pop.PSite(Symbol[ER.name(l) for l in ER.ligands(s)],
                                          ER.multiplicity(s), collect(Symbol, ER.allo_states(s)))
                                for s in sites],
              species = Pop.PSpecies[_pspecies(sp) for sp in species],
              steps = Pop.PStep[Pop.PStep(idx[ER.from_species(s)], idx[ER.to_species(s)],
                                          _pmets(ER.consumed(s)), _pmets(ER.released(s)),
                                          ER.is_equilibrium(s)) for s in flat],
              groups = [collect(off[g] + 1:off[g + 1]) for g in eachindex(groups)],
              family = _family(groups, species, flat, tags, sites, subs, prods))
end

"""
The census state of one run: every structure seen (identify's `seen`), every eq_hash
fitted (identify's `memo`), and the candidate records.
"""
mutable struct Run
    rxn::Any
    label::String
    seen::Set{UInt64}
    memo::Dict{UInt64, Int}            # eq_hash => gid of its fitted representative
    recs::Vector{Dict{Symbol, Any}}    # fresh candidates, in processing order
    mechs::Vector{Any}                 # package mechanism of each record
    funnel::Vector{Dict{Symbol, Any}}  # per round counts
end
Run(rxn, label) = Run(rxn, label, Set{UInt64}(), Dict{UInt64, Int}(), Dict{Symbol, Any}[],
                      Any[], Dict{Symbol, Any}[])

"""
identify's `_process_batch` without the fits: structural dedup against `seen`, then per
fresh mechanism the complexity cap, compile, the parameter cap and the rendered eq_hash;
one representative per eq_hash not already fitted. `items` are `(mechanism, parent gid,
move)`. Returns the funnel counts of the batch.
"""
function process_batch!(run::Run, items, round::Int)
    f = Dict{Symbol, Any}(:round => round, :generated => length(items))
    nseen = 0
    t0 = time()
    first_gid = length(run.recs) + 1
    for (m, par, mv) in items
        h = hash(m)
        if h in run.seen
            nseen += 1
            continue
        end
        push!(run.seen, h)
        rec = Dict{Symbol, Any}(:gid => length(run.recs) + 1, :round => round,
                                :parent => par, :move => mv,
                                :allo => m isa ER.AllostericMechanism)
        t = @elapsed try
            cx = ER._eq_complexity(m)
            rec[:cx] = Int(cx)
            if cx > CX_MAX
                rec[:status] = :cx_skip
            else
                em = ER.compile_mechanism(m)
                fk = ER.fitted_params(em)
                rec[:np] = length(fk)
                if length(fk) > NP_MAX
                    rec[:status] = :param_skip
                else
                    txt = ER.rate_equation_string(em)
                    rec[:eqh] = ER._rate_eq_dedup_key(txt)
                    rec[:status] = :ok
                end
            end
        catch e
            rec[:status] = :error
            rec[:err] = errmsg(e)
        end
        rec[:t_pkg] = t
        push!(run.recs, rec)
        push!(run.mechs, m)
    end
    batch = run.recs[first_gid:end]
    # one representative per eq_hash not already fitted, first in batch order
    nrep = ninh = 0
    for r in batch
        r[:status] === :ok || continue
        h = r[:eqh]
        if haskey(run.memo, h)
            r[:fitted] = false
            r[:rep] = run.memo[h]
            ninh += 1
        else
            run.memo[h] = r[:gid]
            r[:fitted] = true
            r[:rep] = r[:gid]
            nrep += 1
        end
    end
    f[:seen_skip] = nseen
    f[:fresh] = length(batch)
    f[:cx_skip] = count(r -> r[:status] === :cx_skip, batch)
    f[:param_skip] = count(r -> r[:status] === :param_skip, batch)
    f[:error] = count(r -> r[:status] === :error, batch)
    f[:ok] = count(r -> r[:status] === :ok, batch)
    f[:eqhash_dup] = ninh
    f[:fitted] = nrep
    f[:t_pipeline] = time() - t0
    f[:gids] = first_gid:length(run.recs)
    push!(run.funnel, f)
    f
end

"""Round 0: identify's seeds and base tier (`_base_tier`: each degenerate seed replaced by
its non-degenerate flip children), then the batch."""
function round0!(run::Run)
    rxn = run.rxn
    t = @elapsed seeds, source = identify_seeds(rxn)
    ndeg = count(ER._degenerate, seeds)
    base, fails = ER._base_tier(seeds, rxn)
    items = [(m, 0, m in seeds ? :seed : :seed_flip) for m in base]
    f = process_batch!(run, items, 0)
    f[:seed_source] = source
    f[:seeds] = length(seeds)
    f[:degenerate_seeds] = ndeg
    f[:flip_children_added] = count(m -> !(m in seeds), base)
    f[:base_failures] = length(fails)
    f[:t_seeds] = t
    f
end

"""One expansion round from `parents` (gids): `_expand_parents` (first parent wins on
structural dedup), then the batch."""
function expand_round!(run::Run, parents::Vector{Int}, round::Int)
    items = Tuple{Any, Int, Symbol}[]
    firstp = Dict{Any, Int}()
    nraw = 0
    expfail = 0
    t = @elapsed for g in parents
        kids = try
            expand_labeled(run.mechs[g], run.rxn)
        catch e
            expfail += 1
            run.recs[g][:expand_error] = errmsg(e)
            Tuple{Any, Symbol}[]
        end
        nraw += length(kids)
        for (c, mv) in kids
            haskey(firstp, c) && continue
            firstp[c] = g
            push!(items, (c, g, mv))
        end
    end
    f = process_batch!(run, items, round)
    f[:parents] = length(parents)
    f[:raw_children] = nraw
    f[:expand_failures] = expfail
    f[:t_expand] = t
    f
end

"""The frontier of a round: its candidates that compiled and passed the caps (identify
puts every such entry in the frontier, fitted or inherited, assuming a finite loss)."""
frontier(run::Run, round::Int) = [r[:gid] for r in run.recs
                                  if r[:round] == round && r[:status] === :ok]

# ─── the slice side: conversion, law, key ───────────────────────────────────────────────

"""The slice mechanism of a package mechanism (through its Pop export)."""
to_slice(m) = S.Converted(structure(m, "census", "census"), 1).mech

"""
`S.canonical` with concentrations labelled by name instead of by position (the skeptic's
note 1: by position, models that differ only in which modifier binds share a key).
"""
function canonical_named(L)
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
    nm(v) = vars.conc[v] ? "c_" * string(vars.names[v]) : string(vars.names[v])
    labels = [string(i, "|", c, "|", join(("$(nm(v))^$e" for (v, e) in mo if fixed(v)),
                                          "*")) for (i, mo, c) in terms]
    inc = [[(at[v], e) for (v, e) in mo if !fixed(v)] for (_, mo, _) in terms]
    occ = [Tuple{Int, S.Exp}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    best, order, leaves = nothing, Int[], 0
    function search(pcol, tcol)
        pcol, tcol = S.refine(pcol, tcol, inc, occ)
        if allunique(pcol)
            (leaves += 1) > 100_000 && error("canonical labelling: over 100000 leaves")
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
    (key = bytes2hex(sha256(text)), nabsent = length(absent), absent = absent,
     order = [vars.names[present[order]]; absent])
end

"""`N` and `D` of a slice law as sums of branches (products of `(polynomial, power)`), as
the guard probe's `reduced_branches` builds them."""
function law_branches(L)
    L isa S.Law && return (N = Branch[[(L.H, 1), (L.Q, 1)]], D = Branch[[(L.D, 1)]])
    A, T = L, eltype(values(L.H.terms))
    QI, DI = S.inactive_terms(A)
    Lp = S.term(T, S.mono(A.vars.index[:L] => 1))
    WA = [(s[1], s[3]) for s in A.sites]
    WI = [(s[2], s[3]) for s in A.sites]
    n = A.n
    N = Branch[[(A.H, 1); (A.fI, n); (A.active.Q, 1); (A.active.D, n - 1); WA],
               [(A.H, 1); (Lp, 1); (A.fA, n); (QI, 1); (DI, n - 1); WI]]
    D = Branch[[(A.fI, n); (A.active.D, n); WA], [(Lp, 1); (A.fA, n); (DI, n); WI]]
    (N = N, D = D)
end

"""The faces a census checks: every product at zero, every substrate at zero, and each
single reactant at zero when that differs from its side's face."""
function faces_of(sm)
    fs = Tuple{String, Vector{Symbol}}[("P=0", collect(sm.products)),
                                        ("S=0", collect(sm.substrates))]
    for side in (sm.substrates, sm.products)
        length(side) > 1 || continue
        foreach(x -> push!(fs, ("$x=0", [x])), side)
    end
    fs
end

"""The support test on every face: `N ≡ 0`, `D ≡ 0`, and per face whether `D` vanishes
there after the joint shift to lowest terms (0/0 or a pole)."""
function support_faces(N, D, vars, faces)
    concs = [i for i in eachindex(vars.names) if vars.conc[i]]
    Nz, Dz = all(zero_branch, N), all(zero_branch, D)
    live = [x for x in [N; D] if !zero_branch(x)]
    mu = Dict(v => isempty(live) ? 0 // 1 : minimum(bmin(x, v) for x in live) for v in concs)
    attains(f, Z) = (lo = Dict(z => minexp(f, z) for z in Z);
                     any(mo -> all(z -> S.mexp(mo, z) == lo[z], Z), keys(f.terms)))
    survives(Bs, Z) = any(Bs) do x
        zero_branch(x) && return false
        all(z -> bmin(x, z) == mu[z], Z) && all(((f, k),) -> k == 0 || attains(f, Z), x)
    end
    dvan = [!survives(D, [vars.index[x] for x in zs if haskey(vars.index, x)])
            for (_, zs) in faces]
    (Nz = Nz, Dz = Dz, dvan = dvan)
end

"""The guard's face check on one face, on the full graph and, for an allosteric mechanism,
on the inactive conformation's edges."""
function face_verdict(sm, zs)
    f = face_check(sm, zs)
    if f === nothing && S.is_allosteric(sm)
        im = S.inactive_mask(sm)
        any(im) && (f = face_check(sm, zs, im); f === nothing || (f = "inactive: " * f))
    end
    f
end
face_kind(f) = f === nothing ? :none : occursin("bottomless", f) ? :bottomless :
               occursin("split", f) ? :split : occursin("trap", f) ? :trap : :other

"""
Slice-side facts of an ok candidate: the solve, the law, its key, the parameter count and
law size, and the well-posedness guard (Denis's rule): RE turnover, v ≡ 0, 0/0 or a pole
on a face (bottomless segment or a vanishing D), face splits and traps, :OnlyA admission,
and the constructor-level boundary checks. Times each part.
"""
function slice_facts!(rec, m)
    t = time_ns()
    sm = to_slice(m)
    rec[:t_convert] = (time_ns() - t) / 1e9
    t = time_ns()
    sol = S.solve(sm)
    L = S.rate_law(sm, sol)
    rec[:t_law] = (time_ns() - t) / 1e9
    t = time_ns()
    c = canonical_named(L)
    rec[:t_key] = (time_ns() - t) / 1e9
    rec[:key] = c.key
    rec[:nabsent] = c.nabsent
    rec[:absent] = c.absent
    nreg = S.is_allosteric(sm) ? sum((length(ns) for (_, ns) in S.regulator_names(sm)); init = 0) + 1 : 0
    rec[:np_slice] = length(sol.indep) + nreg
    if L isa S.Law
        rec[:dterms] = S.nterms(L.D)
        rec[:qterms] = S.nterms(L.Q)
    else
        QI, DI = S.inactive_terms(L)
        rec[:dterms] = S.nterms(L.active.D) + S.nterms(DI)
        rec[:qterms] = S.nterms(L.active.Q) + S.nterms(QI)
        rec[:dterms_A], rec[:dterms_I] = S.nterms(L.active.D), S.nterms(DI)
        rec[:n] = L.n
    end
    # the guard
    t = time_ns()
    faces = faces_of(sm)
    rec[:faces] = first.(faces)
    br = law_branches(L)
    sup = support_faces(br.N, br.D, L.vars, faces)
    rec[:t_support] = (time_ns() - t) / 1e9
    rec[:v0] = sup.Nz
    rec[:D0] = sup.Dz
    rec[:dvan] = sup.dvan
    t = time_ns()
    rec[:re_turnover] = re_turnover(sm)
    fv = [face_verdict(sm, zs) for (_, zs) in faces]
    rec[:face] = fv
    rec[:face_kind] = face_kind.(fv)
    rec[:t_face] = (time_ns() - t) / 1e9
    t = time_ns()
    rec[:admit] = S.admit(sm)
    rec[:t_admit] = (time_ns() - t) / 1e9
    rec[:role] = (e = S.role_change(sm); e === nothing ? nothing : S.step_text(sm, e))
    rec[:names_err] = try_msg(() -> S.check_names(sm))
    rec[:unused] = unused_reactants(sm)
    sm, sol, L, br
end

# ─── identifiable rank of the slice law (exact log-gradients, BigFloat) ─────────────────

"""Value and log-gradient numerators `Σ_t e_tj·t` of polynomial `p` at `x` for the
variables `θpos` (variable => gradient slot)."""
function poly_grad(p, x::Vector{BigFloat}, θpos::Dict{Int, Int}, nθ)
    F = big(0.0)
    G = zeros(BigFloat, nθ)
    for (mo, c) in p.terms
        t = BigFloat(c)
        for (v, e) in mo
            t *= isinteger(e) ? x[v]^Int(e) : x[v]^(BigFloat(numerator(e)) / denominator(e))
        end
        F += t
        for (v, e) in mo
            j = get(θpos, v, 0)
            j == 0 || (G[j] += e * t)
        end
    end
    F, G
end

"""Value and log-gradient of a sum of branches at `x` (a zero branch contributes nothing)."""
function branches_grad(Bs, x, θpos, nθ, cache)
    V = big(0.0)
    G = zeros(BigFloat, nθ)
    for b in Bs
        zero_branch(b) && continue
        bv = big(1.0)
        bg = zeros(BigFloat, nθ)
        for (f, k) in b
            k == 0 && continue
            F, Gf = get!(() -> poly_grad(f, x, θpos, nθ), cache, objectid(f))
            bv *= F^k
            bg .+= k .* Gf ./ F
        end
        V += bv
        G .+= bv .* bg
    end
    V, G
end

"""Rank of `J` by Gaussian elimination with complete pivoting: pivots above `tol` times the
first; the smallest kept and the largest dropped pivot ratios."""
function rank_gap(J; tol = big(10)^-30)
    A = copy(J); m, n = size(A); piv = BigFloat[]
    for k in 1:min(m, n)
        _, idx = findmax(abs.(A[k:m, k:n]))
        i, j = idx[1] + k - 1, idx[2] + k - 1
        A[[k, i], :] = A[[i, k], :]; A[:, [k, j]] = A[:, [j, k]]
        push!(piv, abs(A[k, k]))
        iszero(A[k, k]) && break
        for r in (k + 1):m
            A[r, k:n] .-= (A[r, k] / A[k, k]) .* A[k, k:n]
        end
    end
    (isempty(piv) || iszero(piv[1]) || !isfinite(piv[1])) && return 0, 0.0, 0.0
    rk = count(q -> q > tol * piv[1], piv)
    rk, Float64(piv[rk] / piv[1]), rk < length(piv) ? Float64(piv[rk + 1] / piv[1]) : 0.0
end

"""
The identifiable rank of the slice law `L` (branches `br`): the rank of `∂ log v/∂ log θ`
over its fitted constants `θ` at `npts` random concentration points (θ, Keq random), exact
gradients in BigFloat. Also the null directions' involvement of each constant (Float64 SVD
of the same matrix) for the attribution.
"""
function slice_rank(L, br, rng; prec = 192)
    setprecision(BigFloat, prec) do
        vars = L.vars
        θv = [i for i in eachindex(vars.names) if !vars.conc[i] && vars.names[i] !== :Keq]
        nθ = length(θv)
        nθ == 0 && return (rank = 0, n = 0, kept = 1.0, dropped = 0.0, involved = Symbol[],
                           names = Symbol[], J = zeros(0, 0))
        θpos = Dict(v => j for (j, v) in enumerate(θv))
        npts = 2nθ + 6
        J = zeros(BigFloat, npts, nθ)
        x = zeros(BigFloat, length(vars.names))
        for i in eachindex(vars.names)
            vars.conc[i] && continue
            x[i] = vars.names[i] === :Keq ? logu(rng, -1, 1) : logu(rng, -1, 1)
        end
        for r in 1:npts
            for i in eachindex(vars.names)
                vars.conc[i] && (x[i] = logu(rng, -1, 1))
            end
            cache = Dict{UInt, Any}()
            Nv, Ng = branches_grad(br.N, x, θpos, nθ, cache)
            Dv_, Dg = branches_grad(br.D, x, θpos, nθ, cache)
            J[r, :] .= Ng ./ Nv .- Dg ./ Dv_
        end
        rk, kept, dropped = rank_gap(J)
        involved = Symbol[]
        if rk < nθ
            F = svd(Float64.(J))
            null = F.V[:, rk + 1:end]
            involved = [vars.names[θv[j]] for j in 1:nθ if norm(null[j, :]) > 1e-6]
        end
        (rank = rk, n = nθ, kept = kept, dropped = dropped, involved = involved,
         names = vars.names[θv], J = Float64.(J))
    end
end

# ─── the package law, evaluated from its own Exprs (no compile) ─────────────────────────

"""A value with its gradient over the log fitted constants (forward mode)."""
struct Dv{T}
    v::T
    g::Vector{T}
end
Base.:+(a::Dv, b::Dv) = Dv(a.v + b.v, a.g .+ b.g)
Base.:+(a::Dv, b::Real) = Dv(a.v + b, copy(a.g))
Base.:+(a::Real, b::Dv) = b + a
Base.:-(a::Dv, b::Dv) = Dv(a.v - b.v, a.g .- b.g)
Base.:-(a::Dv, b::Real) = Dv(a.v - b, copy(a.g))
Base.:-(a::Real, b::Dv) = Dv(a - b.v, .-b.g)
Base.:-(a::Dv) = Dv(-a.v, .-a.g)
Base.:*(a::Dv, b::Dv) = Dv(a.v * b.v, a.v .* b.g .+ b.v .* a.g)
Base.:*(a::Dv, b::Real) = Dv(a.v * b, a.g .* b)
Base.:*(a::Real, b::Dv) = b * a
Base.:/(a::Dv, b::Dv) = Dv(a.v / b.v, (a.g .* b.v .- a.v .* b.g) ./ b.v^2)
Base.:/(a::Dv, b::Real) = Dv(a.v / b, a.g ./ b)
Base.:/(a::Real, b::Dv) = Dv(a / b.v, (-a / b.v^2) .* b.g)
Base.:^(a::Dv, n::Integer) = Dv(a.v^n, (n * a.v^(n - 1)) .* a.g)
Base.:^(a::Dv, r::Real) = Dv(a.v^r, (r * a.v^(r - 1)) .* a.g)
Base.inv(a::Dv) = 1 / a
Base.sqrt(a::Dv) = a^0.5

"""Evaluate a rate-equation Expr (`+ - * / ^ inv sqrt`, numbers, symbols) in `env`."""
evx(x::Symbol, env) = env[x]
evx(x::Number, env) = x
function evx(ex::Expr, env)
    ex.head === :call || error("unexpected Expr head $(ex.head)")
    f = ex.args[1]
    if f === :^
        b = evx(ex.args[2], env)
        e = ex.args[3]
        e isa Number || (e = evx(e, env); e isa Dv && (e = e.v))
        return b^e
    end
    a = [evx(y, env) for y in @view ex.args[2:end]]
    f === :+ ? foldl(+, a) : f === :- ? (length(a) == 1 ? -a[1] : a[1] - a[2]) :
    f === :* ? foldl(*, a) : f === :/ ? a[1] / a[2] : f === :inv ? inv(a[1]) :
    f === :sqrt ? sqrt(a[1]) : error("unexpected operator $f")
end

"""The package's rate law of `em` as `rate_equation` generates it: dependents, fitted
constants, numerator and denominator Exprs, metabolite names."""
function pkg_law(em)
    M = typeof(em)
    dep, indep = ER._dependent_param_exprs(M)
    num, den = ER._num_den_exprs(M)
    (dep = sort!(collect(dep); by = first), indep = collect(indep), num = num, den = den,
     mets = collect(ER.metabolites(em)))
end

"""The package rate at fitted values `θ` (in `pl.indep` order), `keq` and concentrations
`conc` (a Dict by name, missing names at zero); values of type `T` (a `Dv` for gradients)."""
function pkg_rate(pl, θ, keq, conc, wrap)
    env = Dict{Symbol, Any}()
    for (j, s) in enumerate(pl.indep)
        env[s] = θ[j]
    end
    env[:Keq] = wrap(keq)
    env[:E_total] = wrap(one(keq))
    for x in pl.mets
        env[x] = wrap(get(conc, x, zero(keq)))
    end
    for (s, rhs) in pl.dep
        env[s] = evx(rhs, env)
    end
    n, d = evx(pl.num, env), evx(pl.den, env)
    n isa Dv || d isa Dv ? (n isa Dv ? n : wrap(n)) / d : n / d
end

"""
The identifiable rank of the package law (what identify fits): `∂v/∂ log θ` over its
fitted constants at `np + 6` random concentration points, exact forward-mode gradients in
Float64 and the SVD; an ambiguous gap is redone in BigFloat with complete pivoting.
"""
function pkg_rank(pl, rng)
    n = length(pl.indep)
    n == 0 && return (rank = 0, n = 0, kept = 1.0, dropped = 0.0, big = false)
    function jac(T)
        lu(lo, hi) = T(10)^(lo + (hi - lo) * T(rand(rng)))
        θv = [lu(-1, 1) for _ in 1:n]
        keq = lu(-1, 1)
        npts = n + 6
        J = zeros(T, npts, n)
        θ = [Dv(θv[j], [i == j ? θv[j] : zero(T) for i in 1:n]) for j in 1:n]
        wrap(x) = Dv(T(x), zeros(T, n))
        for r in 1:npts
            conc = Dict(x => lu(-1, 1) for x in pl.mets)
            v = pkg_rate(pl, θ, keq, conc, wrap)
            J[r, :] .= v.g
        end
        J
    end
    J = jac(Float64)
    all(isfinite, J) || return (rank = -1, n = n, kept = NaN, dropped = NaN, big = false)
    sv = svdvals(J)
    sv[1] == 0 && return (rank = 0, n = n, kept = 0.0, dropped = 0.0, big = false)
    rk = count(>(1e-9 * sv[1]), sv)
    kept, dropped = sv[rk] / sv[1], rk < length(sv) ? sv[rk + 1] / sv[1] : 0.0
    if kept < 1e-6 || dropped > 1e-12
        return setprecision(BigFloat, 192) do
            rkb, keptb, droppedb = rank_gap(jac(BigFloat))
            (rank = rkb, n = n, kept = keptb, dropped = droppedb, big = true)
        end
    end
    (rank = rk, n = n, kept = kept, dropped = dropped, big = false)
end

"""Whether the package law grows without bound as every substrate (`side = :fwd`, products
at zero) or every product (`:rev`, substrates at zero) scales from 1e6 to 1e9, every other
metabolite at zero, at two random parameter points (the test suite's degeneracy probe)."""
function pkg_unbounded(pl, subs, prods, rng)
    out = Dict(:fwd => false, :rev => false)
    for _ in 1:2
        θ = [10.0^(2rand(rng) - 1) for _ in pl.indep]
        keq = 10.0^(2rand(rng) - 1)
        for (side, names) in ((:fwd, subs), (:rev, prods))
            v(c) = pkg_rate(pl, θ, keq, Dict(x => c for x in names), identity)
            r = v(1e9) / v(1e6)
            isfinite(r) && abs(r) > 10 && (out[side] = true)
        end
    end
    out
end

"""The package's kcat precondition (`_kcat_keys`, which `rescale_parameter_values` reaches
after every fit with the default `scale_k_to_kcat`): its error message, or `nothing`."""
function kcat_error(m, em)
    try
        prods = Set(ER.name(x) for x in ER.products(ER.reaction(m)))
        if m isa ER.Mechanism
            num, den, _ = ER._raw_symbolic_rate_polys(typeof(em))
            ng, dg = ER._kcat_groups_from_polys(num, den, Set(ER.metabolites(em)))
            ER._kcat_keys(ng, dg, prods)
        else
            cat_mets = Set{Symbol}(ER.metabolites(ER.catalytic_mechanism(em)))
            nA, dA, _, _, _, dI = ER._mwc_state_polys(m, cat_mets)
            ng, dg = ER._kcat_groups_from_polys(ER.poly_mul(nA, dI), ER.poly_mul(dA, dI),
                                                cat_mets)
            ER._kcat_keys(ng, dg, prods)
        end
        nothing
    catch e
        errmsg(e)
    end
end

"""Package-side facts of an ok candidate: conditions V (no Vmax one way) and C, the kcat
precondition, the numeric unboundedness probe, and (when `rank`) the identifiable rank."""
function pkg_facts!(rec, m, rng; rank = false)
    rxn = ER.reaction(m)
    gs = ER.steps(m)
    t = time_ns()
    rec[:vmax_fwd] = ER._has_vmax(gs, rxn, ER.Substrate)
    rec[:vmax_rev] = ER._has_vmax(gs, rxn, ER.Product)
    rec[:condC] = ER._chemistry_equilibrates_both_sides(gs, rxn)
    rec[:t_vmax] = (time_ns() - t) / 1e9
    em = ER.compile_mechanism(m)
    rec[:kcat_err] = kcat_error(m, em)
    t = time_ns()
    pl = pkg_law(em)
    rec[:t_pkglaw] = (time_ns() - t) / 1e9
    subs = [ER.name(x) for x in ER.substrates(rxn)]
    prods = [ER.name(x) for x in ER.products(rxn)]
    t = time_ns()
    u = pkg_unbounded(pl, subs, prods, rng)
    rec[:unb_fwd], rec[:unb_rev] = u[:fwd], u[:rev]
    rec[:t_unbounded] = (time_ns() - t) / 1e9
    if rank
        t = time_ns()
        r = pkg_rank(pl, rng)
        rec[:rank_pkg] = r.rank
        rec[:rank_pkg_gap] = (r.kept, r.dropped, r.big)
        rec[:t_rank_pkg] = (time_ns() - t) / 1e9
    end
    pl
end

# ─── attribution of a slice rank deficit ────────────────────────────────────────────────

"""
The kind of each fitted constant `names` of the slice law: `:L_phantom_n1` (`L` at one
catalytic subunit), `:L_other`, `:regulator`, `:ss_noflux` / `:ss_flux` (a constant of a
steady-state edge without / with net flux in its conformation; an inactive constant counts
as flux-free when the inactive law turns nothing over), `:re_iso` (the `K` of a
rapid-equilibrium isomerization), `:re_binding`, `:re_other`.
"""
function kinds_of(sm, sol, L, names)
    flanks = Set(flank_groups(sm))
    cn = S.column_names(sm, sol)
    col = Dict(cn[i] => i for i in sol.indep)
    flux = S.edge_flux(sm)
    dead_I = L isa S.AlloLaw && iszero(first(S.inactive_terms(L)))
    regs = Set(n for (_, ns) in (S.is_allosteric(sm) ? S.regulator_names(sm) : []) for n in ns)
    map(names) do x
        x === :L && return sm.multiplicity == 1 ? :L_phantom_n1 : :L_other
        x in regs && return :regulator
        haskey(col, x) || return :unknown
        c = sol.cols[col[x]]
        e = sm.edges[c.edge]
        nonequal = S.is_allosteric(sm) && sm.tags[sm.group[c.edge]] === :NonequalAI
        if !e.re
            noflux = c.inactive ? (dead_I || !flux[c.edge]) : !flux[c.edge]
            noflux && return :ss_noflux
            sm.group[c.edge] in flanks && return :chain_flank
            return c.inactive && nonequal ? :nonequal_copy : :ss_flux
        end
        c.inactive && nonequal && return :nonequal_copy
        e.cons == 0 && e.rel == 0 ? :re_iso : S.is_binding(e) ? :re_binding : :re_other
    end
end

"""The causes a deficit is split by, in the order they are tried."""
const CAUSES = (:L_phantom_n1, :ss_noflux, :chain_flank, :nonequal_copy, :re_iso)

"""
The kinetic groups that flank a steady-state isomerization alone in its group, when each of
its forms has one other edge, a binding into that form alone in its group (the slice's
`chain_flanks`, here also for an allosteric mechanism): by the chain lemma, such a flank at
steady state adds a constant and no rate law in a plain mechanism.
"""
function flank_groups(m)
    alone = [count(==(g), m.group) == 1 for g in m.group]
    at = [Int[] for _ in m.forms]
    for (k, e) in enumerate(m.edges)
        push!(at[e.a], k)
        push!(at[e.b], k)
    end
    function flank(k0, x)
        rest = filter(!=(k0), at[x])
        length(rest) == 1 && S.is_binding(m.edges[rest[1]]) && m.edges[rest[1]].b == x &&
            alone[rest[1]] ? Int(m.group[rest[1]]) : 0
    end
    out = Int[]
    for (k, e) in enumerate(m.edges)
        alone[k] && e.cons == e.rel == 0 && !e.re || continue
        f1, f2 = flank(k, e.a), flank(k, e.b)
        f1 > 0 && f2 > 0 && push!(out, f1, f2)
    end
    out
end

"""Numeric rank of `J` (singular values above `tol` times the largest)."""
svdrank(J; tol = 1e-9) = (isempty(J) ? 0 : (sv = svdvals(J); sv[1] == 0 ? 0 :
                                              count(>(tol * sv[1]), sv)))

"""
The deficit `n - rk` split by cause: for each cause in turn (`:L_phantom_n1`, `:ss_noflux`,
`:re_iso`), drop each of its constants whose column the kept ones span (dropping it keeps
the rank); what is left is `:other`. Also whether the Float64 rank agrees with `rk`.
"""
function deficit_parts(J, kinds, rk)
    n = size(J, 2)
    parts = Dict{Symbol, Int}()
    agree = svdrank(J) == rk
    keep = trues(n)
    left = n - rk
    for cause in CAUSES, j in 1:n
        left == 0 && break
        kinds[j] === cause || continue
        keep[j] = false
        if svdrank(J[:, keep]) == rk
            left -= 1
            parts[cause] = get(parts, cause, 0) + 1
        else
            keep[j] = true
        end
    end
    left > 0 && (parts[:other] = left)
    parts, agree
end

"""
The kinds of the constants in the slice law's null directions: `:L` at one catalytic
subunit (the L phantom), a constant of a steady-state edge that carries no flux in its
conformation, the `K` of a rapid-equilibrium isomerization, a regulator constant, or
another constant (by edge kind). An inactive constant counts as flux-free when the
inactive law turns nothing over.
"""
function attribute(sm, sol, L, involved)
    kinds = Symbol[]
    isempty(involved) && return kinds
    names = S.column_names(sm, sol)
    col = Dict(names[i] => i for i in sol.indep)
    flux = S.edge_flux(sm)
    dead_I = L isa S.AlloLaw && iszero(first(S.inactive_terms(L)))
    regs = Set(n for (_, ns) in (S.is_allosteric(sm) ? S.regulator_names(sm) : []) for n in ns)
    for x in involved
        if x === :L
            push!(kinds, sm.multiplicity == 1 ? :L_phantom_n1 : :L_other)
        elseif x in regs
            push!(kinds, :regulator)
        elseif haskey(col, x)
            c = sol.cols[col[x]]
            e = sm.edges[c.edge]
            if !e.re
                noflux = c.inactive ? (dead_I || !flux[c.edge]) : !flux[c.edge]
                push!(kinds, noflux ? :ss_noflux : :ss_flux)
            elseif e.cons == 0 && e.rel == 0
                push!(kinds, :re_iso)
            else
                push!(kinds, S.is_binding(e) ? :re_binding : :re_other)
            end
        else
            push!(kinds, :unknown)
        end
    end
    kinds
end

"""The deficit's sub-class from the kinds of its null-direction constants."""
p_subclass(kinds) = :L_phantom_n1 in kinds ? :L_phantom_n1 :
                    :ss_noflux in kinds ? :ss_noflux :
                    :re_iso in kinds ? :re_iso_gauge : :other
