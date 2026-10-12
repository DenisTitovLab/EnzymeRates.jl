# ABOUTME: Census library: today's identify stages (seen, caps, eq_hash) on enumeration rounds,
# ABOUTME: and per-candidate slice classes (guard, rename-proof key, rank, saturation, size).

using EnzymeRates, Random, Serialization, Printf, Statistics, LinearAlgebra
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
include(joinpath(SCRATCH, "proto", "common", "export", "export_lib.jl"))   # defines const ER
include(joinpath(SCRATCH, "slice", "src", "Slice.jl"))
const S = Slice
const TESTDIR = joinpath(dirname(dirname(pathof(EnzymeRates))), "test")
const CX_FILTER = 337      # identify_rate_equation default eq_complexity_filter
const MAX_PARAMS = 20      # identify_rate_equation default max_param_count

errmsg(e) = replace(first(sprint(showerror, e), 300), '\n' => ' ')
try_msg(f) = try
    f()
    nothing
catch e
    errmsg(e)
end
secs(t0) = (time_ns() - t0) / 1e9

"""Top-level `const` statements of a test file whose names are in `names`."""
function test_consts(file, names)
    ast = Meta.parseall(read(joinpath(TESTDIR, file), String))
    [ex for ex in ast.args if ex isa Expr && ex.head === :const &&
     ex.args[1] isa Expr && ex.args[1].args[1] in names]
end
const RXNS = (:uni_uni_rxn, :uni_bi_rxn, :bi_bi_rxn, :bi_bi_pp_rxn, :ter_bi_rxn, :ter_ter_rxn,
              :pyruvate_carboxylase_rxn, :pyruvate_dehydrogenase_rxn)
foreach(ex -> Core.eval(Main, ex), test_consts("test_mechanism_enumeration.jl", RXNS))

# ─── package mechanism → slice mechanism (the verified Pop route) ────────────────────────

"""The structure of a package mechanism `m`, plain or allosteric, as a `Pop.PMech`
(as `slice/integrate/export_terter_expand.jl` exports it, without the family tags)."""
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
              groups = [collect(off[g] + 1:off[g + 1]) for g in eachindex(groups)])
end
to_slice(m) = S.Converted(structure(m, "c", "c"), 1).mech

"""The steps of a slice mechanism as DSL text, kinetic groups in parentheses, tags in
brackets."""
function mech_text(m::S.Mechanism)
    parts = String[]
    for g in 1:S.ngroups(m)
        t = [S.step_text(m, m.edges[k]) for k in S.members(m, g)]
        s = length(t) == 1 ? t[1] : "(" * join(t, ", ") * ")"
        S.is_allosteric(m) && (s *= " [" * string(m.tags[g]) * "]")
        push!(parts, s)
    end
    pre = S.is_allosteric(m) ? "allosteric n=$(m.multiplicity): " : ""
    pre * join(parts, "; ")
end

# ─── guard checks, copied from round4/guard/lib.jl (re_turnover, face_check) ─────────────

"""The first rapid-equilibrium edge (in `mask`) that closes a cycle of RE edges whose net
metabolite change is not zero (turnover inside one RE segment: no finite rate), as text."""
function re_turnover(m, mask = trues(length(m.edges)))
    segment, parent, order = S.re_forest(m, mask)
    tree = Set(abs(p) for p in parent if p != 0)
    _, names, ex = S.segment_offsets(m, mask)
    row = [searchsortedfirst(names, x.name) for x in m.mets]
    for (k, e) in enumerate(m.edges)
        mask[k] && e.re && !(k in tree) || continue
        d = ex[:, e.b] .- ex[:, e.a]
        for j in eachindex(m.mets)
            d[row[j]] -= S.count_at(e.cons, j) - S.count_at(e.rel, j)
        end
        any(!iszero, d) && return S.step_text(m, e)
    end
    nothing
end

"""
The rapid-equilibrium segments of `mask` at the face where the names `zs` are zero, as the
fast dynamics see them (round4/guard/lib.jl, verbatim): bottomless, a split of the
surviving forms, or a trap entered on the face; the first violation as text, or `nothing`.
"""
function face_check(m, zs, mask = trues(length(m.edges)))
    segment, names, ex = S.segment_offsets(m, mask)
    zr = [r for (r, x) in enumerate(names) if x in zs]
    nf = length(m.forms)
    survive = [all(r -> ex[r, i] == 0, zr) for i in 1:nf]
    takes(p) = any(j -> S.count_at(p, j) > 0 && m.mets[j].name in zs, eachindex(m.mets))
    fast = [Int[] for _ in 1:nf]
    alive = Tuple{Int, Int}[]
    touched = falses(nf)
    for (k, e) in enumerate(m.edges)
        mask[k] || continue
        touched[e.a] = touched[e.b] = true
        takes(e.cons) || (push!(alive, (e.a, e.b)); e.re && push!(fast[e.a], e.b))
        takes(e.rel) || (push!(alive, (e.b, e.a)); e.re && push!(fast[e.b], e.a))
    end
    reach(i) = (seen = falses(nf); seen[i] = true; todo = [i];
                while !isempty(todo); u = pop!(todo); for v in fast[u]
                    seen[v] || (seen[v] = true; push!(todo, v)); end; end; seen)
    for s in 1:maximum(segment; init = 0)
        fs = [i for i in findall(==(s), segment) if touched[i]]
        length(fs) < 2 && continue
        sv = filter(i -> survive[i], fs)
        isempty(sv) && return "bottomless {" * join((S.form_text(m, i) for i in fs), ", ") * "}"
        r = reach(sv[1])
        all(i -> r[i], sv) || return "surviving forms " * S.form_text(m, sv[1]) * " and " *
            S.form_text(m, sv[findfirst(i -> !r[i], sv)]) * " split on the face"
        for i in fs
            survive[i] && continue
            ri = reach(i)
            any(j -> ri[j] && survive[j], fs) && continue
            trap = findall(ri)                      # closed: no fast arc to a survivor
            entered = any(((a, b),) -> b in trap && !(a in trap), alive)
            entered && return "trap {" * join((S.form_text(m, j) for j in trap), ", ") *
                               "} entered on the face"
        end
    end
    nothing
end

"""The faces where data sit: every product zero, every substrate zero, each reactant zero."""
face_sets(m::S.Mechanism) =
    [("P=0", m.products); ("S=0", m.substrates);
     [(string(x, "=0"), [x]) for x in [m.substrates; m.products]]]

# ─── package stage: what `_process_batch` computes before fitting ────────────────────────

"""The package's pass-1 record of `m`: `:cx_skip` over the complexity filter, `:param_skip`
over the parameter cap, `:fail` when compiling or rendering throws, else `:ok` with the
fitted count and `eq_hash` (`_rate_eq_dedup_key` of `rate_equation_string`)."""
function pkg_stage(m)
    t0 = time_ns()
    cx = try
        ER._eq_complexity(m)
    catch e
        return (status = :fail, cx = -1, n_params = -1, eq_hash = UInt64(0), err = errmsg(e),
                t_cx = secs(t0), t_dep = 0.0, t_render = 0.0, eq_len = 0)
    end
    t_cx = secs(t0)
    cxi = cx > typemax(Int) ? typemax(Int) : Int(cx)
    cx > CX_FILTER && return (status = :cx_skip, cx = cxi, n_params = -1, eq_hash = UInt64(0),
                              err = "", t_cx, t_dep = 0.0, t_render = 0.0, eq_len = 0)
    t1 = time_ns()
    n = -1
    try
        em = ER.compile_mechanism(m)
        _, indep = ER._dependent_param_exprs(typeof(em))   # what `fitted_params(em)` returns
        n = length(indep)
        t_dep = secs(t1)
        n > MAX_PARAMS && return (status = :param_skip, cx = cxi, n_params = n,
                                  eq_hash = UInt64(0), err = "", t_cx, t_dep, t_render = 0.0,
                                  eq_len = 0)
        t2 = time_ns()
        txt = ER.rate_equation_string(em)
        h = ER._rate_eq_dedup_key(txt)
        return (status = :ok, cx = cxi, n_params = n, eq_hash = h, err = "", t_cx, t_dep,
                t_render = secs(t2), eq_len = length(txt))
    catch e
        return (status = :fail, cx = cxi, n_params = n, eq_hash = UInt64(0), err = errmsg(e),
                t_cx, t_dep = secs(t1), t_render = 0.0, eq_len = 0)
    end
end

# ─── slice stage: guard checks, law, key, law tests, rank ───────────────────────────────

const PolyT = Union{S.Poly{Int}, S.Poly{BigInt}}
big_poly(p) = convert(S.Poly{BigInt}, p)

"""The expanded numerator `N`, denominator `D` and `N/H` of law `L` (a `Law` or an `AlloLaw`),
the MWC combination expanded for an allosteric law; in `BigInt` after an overflow."""
function expand_law(L; big = false)
    cv(p) = big ? big_poly(p) : p
    if L isa S.Law
        H, Q, D = cv(L.H), cv(L.Q), cv(L.D)
        return L.vars, H * Q, D, Q
    end
    A = L
    T = big ? BigInt : eltype(values(A.H.terms))
    QI, DI = S.inactive_terms(A)
    QA, DA, fA, fI, H = cv(A.active.Q), cv(A.active.D), cv(A.fA), cv(A.fI), cv(A.H)
    QI, DI = cv(QI), cv(DI)
    Lp = S.term(T, S.mono(A.vars.index[:L] => 1))
    one_ = one(S.Poly{T})
    pw(p, k) = k == 0 ? one_ : foldl(*, (p for _ in 1:k))
    WA = foldl(*, (pw(cv(s[1]), s[3]) for s in A.sites); init = one_)
    WI = foldl(*, (pw(cv(s[2]), s[3]) for s in A.sites); init = one_)
    n = A.n
    Nb = pw(fI, n) * QA * pw(DA, n - 1) * WA + Lp * pw(fA, n) * QI * pw(DI, n - 1) * WI
    D = pw(fI, n) * pw(DA, n) * WA + Lp * pw(fA, n) * pw(DI, n) * WI
    A.vars, H * Nb, D, Nb
end
function expand_law_checked(L)
    try
        expand_law(L)
    catch e
        e isa OverflowError || rethrow()
        expand_law(L; big = true)
    end
end

"""A polynomial as Float64 terms for evaluation with gradients."""
struct CTerm
    coef::Float64
    vars::Vector{Int}
    exps::Vector{Float64}
end
cpoly(p) = [CTerm(Float64(c), [v for (v, _) in mo], [Float64(e) for (_, e) in mo])
            for (mo, c) in p.terms]

"""The value of `P` at `exp.(lx)` and, into `g`, its derivative in the log of each variable
`v` with `fpos[v] > 0` (column `fpos[v]`)."""
function evalgrad!(g, P::Vector{CTerm}, lx, fpos)
    fill!(g, 0.0)
    val = 0.0
    @inbounds for t in P
        s = 0.0
        for i in eachindex(t.vars)
            s += t.exps[i] * lx[t.vars[i]]
        end
        tv = t.coef * exp(s)
        val += tv
        for i in eachindex(t.vars)
            j = fpos[t.vars[i]]
            j > 0 && (g[j] += t.exps[i] * tv)
        end
    end
    val
end

"""Rank of `A` with singular values above `tol`."""
rank_tol(A, tol) = size(A, 2) == 0 ? 0 : count(>(tol), svdvals(A))

"""
The Jacobian `dv/dlog θ` (rows: random concentration points at one random parameter point,
each row scaled to unit max) of the law `N/D` over the fitted names `fnames`, its rank, and
the singular-value ratios kept/dropped around the rank.
"""
function law_jacobian(vars, N, D, fnames, rng; tol = 1e-9)
    nv = length(vars.names)
    fpos = zeros(Int, nv)
    for (j, nm) in enumerate(fnames)
        fpos[vars.index[nm]] = j
    end
    np = length(fnames)
    cN, cD = cpoly(N), cpoly(D)
    lx = zeros(nv)
    for v in 1:nv
        vars.conc[v] || (lx[v] = log(10.0) * (2 * rand(rng) - 1))
    end
    rows = 2np + 6
    J = zeros(rows, np)
    gN, gD = zeros(np), zeros(np)
    for r in 1:rows
        for v in 1:nv
            vars.conc[v] && (lx[v] = log(10.0) * (2 * rand(rng) - 1))
        end
        n = evalgrad!(gN, cN, lx, fpos)
        d = evalgrad!(gD, cD, lx, fpos)
        @. J[r, :] = (gN * d - n * gD) / d^2
        mx = maximum(abs, view(J, r, :); init = 0.0)
        mx > 0 && (J[r, :] ./= mx)
    end
    sv = np == 0 ? Float64[] : svdvals(J)
    s1 = isempty(sv) ? 1.0 : max(sv[1], eps())
    rk = count(>(tol * s1), sv)
    kept = rk == 0 ? 0.0 : sv[rk] / s1
    dropped = rk < length(sv) ? sv[rk + 1] / s1 : 0.0
    J, rk, kept, dropped, tol * s1
end

"""
Greedy attribution of the rank deficiency of `J` (absolute tolerance `atol`) over the
columns described by `info` (`kind`, `class`, `flux`): columns no term reaches (`absent`),
`L` when it adds no rank (`L_phantom`), deficiency inside one steady-state class's own
columns (`ss_local`; `ss_local_fluxfree` when that class carries no flux), a
rapid-equilibrium isomerization constant that adds no rank (`re_iso_gauge`), the rest
(`other`).
"""
function explain(J, atol, info)
    np = size(J, 2)
    r(cols) = rank_tol(J[:, cols], atol)
    rem = collect(1:np)
    out = Dict{Symbol, Int}()
    add!(k, n = 1) = (out[k] = get(out, k, 0) + n)
    for j in 1:np
        maximum(abs, view(J, :, j); init = 0.0) <= atol && (add!(:absent); filter!(!=(j), rem))
    end
    base = r(rem)
    Lc = findfirst(i -> info[i].kind === :L, rem)
    if Lc !== nothing
        c = rem[Lc]
        others = filter(!=(c), rem)
        if r(others) == base
            add!(:L_phantom)
            rem = others
        end
    end
    for cls in unique(info[i].class for i in rem if info[i].kind === :ss)
        cols = [i for i in rem if info[i].kind === :ss && info[i].class == cls]
        d = length(cols) - r(cols)
        d > 0 || continue
        add!(all(i -> !info[i].flux, cols) ? :ss_local_fluxfree : :ss_local, d)
        filter!(i -> !(i in cols[end - d + 1:end]), rem)
    end
    base = r(rem)
    for c in [i for i in rem if info[i].kind === :re_iso]
        others = filter(!=(c), rem)
        if r(others) == base
            add!(:re_iso_gauge)
            rem = others
        end
    end
    d = length(rem) - r(rem)
    d > 0 && add!(:other, d)
    out
end

"""The kind, kinetic class and flux of each fitted column of `m` under `sol`, then the
regulator constants and `L` of an allosteric mechanism."""
function column_info(m, sol)
    flux = S.edge_flux(m)
    info = NamedTuple[]
    for i in sol.indep
        c = sol.cols[i]
        e = m.edges[c.edge]
        kind = !e.re ? :ss : S.is_binding(e) ? :re_bind :
               (e.cons == 0 && e.rel == 0) ? :re_iso : :re_chem
        cls = S.class(m, c)
        f = any(k -> flux[k] && S.class(m, S.Col(k, c.slot, c.inactive)) == cls,
                eachindex(m.edges))
        push!(info, (kind = kind, class = (cls[1], cls[3], cls[4]), flux = f))
    end
    if S.is_allosteric(m)
        for (_, ns) in S.regulator_names(m), _ in ns
            push!(info, (kind = :reg, class = (0, false, false), flux = true))
        end
        push!(info, (kind = :L, class = (0, false, false), flux = true))
    end
    info
end

"""The fitted names of `m` under `sol` (`Slice.fitted_params` without a second solve)."""
function slice_fitted(m, sol)
    out = S.column_names(m, sol)[sol.indep]
    if S.is_allosteric(m)
        append!(out, (n for (_, ns) in S.regulator_names(m) for n in ns))
        push!(out, :L)
    end
    out
end

"""
Law tests on `N/D` after the joint shift to lowest terms in the concentrations: `v0` (N ≡ 0),
the faces where `D` vanishes (`dead`) or `N` vanishes (`nzero`), whether the rate grows
without bound as the substrates (`nosat_fwd`, products at zero) or the products
(`nosat_rev`, substrates at zero) grow together, and the reactants in which it grows without
bound alone (`nosat_single`, the other side at zero), and the concentration skeleton of
`N/H` and `D`.
"""
function law_tests(m, vars, N, D, Nb)
    sh = S.conc_shift((N, D), vars)
    N2, D2, Nb2 = N * sh, D * sh, Nb * sh
    idx(xs) = [vars.index[x] for x in xs]
    si, pi_ = idx(m.substrates), idx(m.products)
    faces = [(nm, idx(zs)) for (nm, zs) in face_sets(m)]
    survives(P, Z) = any(mo -> all(z -> S.mexp(mo, z) == 0, Z), keys(P.terms))
    dead = [nm for (nm, Z) in faces if !survives(D2, Z)]
    nzero = [nm for (nm, Z) in faces if !survives(N2, Z)]
    function topdeg(P, zero, grow)
        best = -1 // 1
        found = false
        for mo in keys(P.terms)
            all(z -> S.mexp(mo, z) == 0, zero) || continue
            d = sum((S.mexp(mo, x) for x in grow); init = 0 // 1)
            best = found ? max(best, d) : d
            found = true
        end
        found, best
    end
    function unbounded(grow, zero)
        fn, dn = topdeg(N2, zero, grow)
        fd, dd = topdeg(D2, zero, grow)
        fn && fd && dn > dd
    end
    single = [vars.names[x] for x in si if unbounded([x], pi_)]
    append!(single, [vars.names[x] for x in pi_ if unbounded([x], si)])
    concs(P) = sort!([string(k) for k in keys(S.split_conc(P, vars))])
    skel = hash((concs(Nb2), concs(D2)))
    (v0 = iszero(N2), dead = dead, nzero = nzero, nosat_fwd = unbounded(si, pi_),
     nosat_rev = unbounded(pi_, si), nosat_single = single, skeleton = skel,
     nN = S.nterms(N2), nD = S.nterms(D2))
end

"""Every slice-side verdict and measurement of package mechanism `mp`, with timings."""
function slice_stage(mp, rng)
    out = Dict{Symbol, Any}()
    t = time_ns()
    m = try
        to_slice(mp)
    catch e
        out[:conv_err] = errmsg(e)
        return out, nothing
    end
    out[:t_conv] = secs(t)
    out[:text] = mech_text(m)
    out[:allo] = S.is_allosteric(m)
    out[:mult] = m.multiplicity
    t = time_ns(); out[:boundary] = try_msg(() -> S.check_boundary(m)); out[:t_boundary] = secs(t)
    t = time_ns()
    out[:admit] = try
        S.admit(m)
    catch e
        "admit error: " * errmsg(e)
    end
    out[:t_admit] = secs(t)
    t = time_ns(); out[:re_turn] = re_turnover(m); out[:t_returnover] = secs(t)
    t = time_ns()
    fres = Pair{String, String}[]
    imask = S.is_allosteric(m) ? S.inactive_mask(m) : falses(length(m.edges))
    for (nm, zs) in face_sets(m)
        f = try
            face_check(m, zs)
        catch e
            "face_check error: " * errmsg(e)
        end
        if f === nothing && any(imask)
            f2 = face_check(m, zs, imask)
            f2 === nothing || (f = "inactive: " * f2)
        end
        f === nothing || push!(fres, nm => f)
    end
    out[:faces] = fres
    out[:t_faces] = secs(t)
    t = time_ns()
    sol, L = try
        sol = S.solve(m)
        sol, S.rate_law(m, sol)
    catch e
        out[:law_err] = errmsg(e)
        return out, m
    end
    out[:t_law] = secs(t)
    t = time_ns()
    out[:key] = try
        S.canonical(L).key
    catch e
        "key error: " * errmsg(e)
    end
    out[:t_key] = secs(t)
    t = time_ns()
    vars, N, D, Nb = expand_law_checked(L)
    out[:t_expand] = secs(t)
    t = time_ns()
    lt = law_tests(m, vars, N, D, Nb)
    foreach(k -> out[k] = getfield(lt, k), keys(lt))
    out[:t_tests] = secs(t)
    t = time_ns()
    fn = slice_fitted(m, sol)
    out[:n_slice] = length(fn)
    if lt.v0
        out[:rank] = -1
    else
        J, rk, kept, dropped, atol = law_jacobian(vars, N, D, fn, rng)
        out[:rank], out[:sv_kept], out[:sv_dropped] = rk, kept, dropped
        out[:t_rank] = secs(t)
        if rk < length(fn)
            t = time_ns()
            out[:explain] = explain(J, atol, column_info(m, sol))
            out[:t_explain] = secs(t)
        end
    end
    out, m
end

# ─── the package law's own rank: its Exprs on a tape, complex-step derivatives ──────────

"""A straight-line program of an Expr over `env` slots: op 1 constant, 2 slot, 3 +, 4 −,
5 ×, 6 ÷, 7 negation, 8 integer power, 9 real power."""
struct Tape
    op::Vector{Int8}
    a::Vector{Int32}
    b::Vector{Int32}
    c::Vector{Float64}
end
Tape() = Tape(Int8[], Int32[], Int32[], Float64[])
function emit!(t::Tape, op, a = 0, b = 0, c = 0.0)
    push!(t.op, op); push!(t.a, a); push!(t.b, b); push!(t.c, c)
    length(t.op)
end
function build!(t::Tape, ex, env::Dict{Symbol, Int})
    ex isa Symbol && return emit!(t, 2, env[ex])
    ex isa Number && return emit!(t, 1, 0, 0, Float64(ex))
    (ex isa Expr && ex.head === :call) || error("unexpected expression $ex")
    f, args = ex.args[1], ex.args[2:end]
    if f === :+ || f === :*
        i = build!(t, args[1], env)
        for a in args[2:end]
            i = emit!(t, f === :+ ? 3 : 5, i, build!(t, a, env))
        end
        return i
    elseif f === :-
        length(args) == 1 && return emit!(t, 7, build!(t, args[1], env))
        return emit!(t, 4, build!(t, args[1], env), build!(t, args[2], env))
    elseif f === :/
        return emit!(t, 6, build!(t, args[1], env), build!(t, args[2], env))
    elseif f === :^
        p = args[2]
        p isa Number || error("non-constant exponent $p")
        return emit!(t, isinteger(p) ? 8 : 9, build!(t, args[1], env), 0, Float64(p))
    end
    error("unexpected operator $f")
end
function run_tape(t::Tape, env::Vector{ComplexF64}, buf::Vector{ComplexF64})
    @inbounds for i in eachindex(t.op)
        op = t.op[i]
        buf[i] = op == 1 ? complex(t.c[i]) : op == 2 ? env[t.a[i]] :
                 op == 3 ? buf[t.a[i]] + buf[t.b[i]] : op == 4 ? buf[t.a[i]] - buf[t.b[i]] :
                 op == 5 ? buf[t.a[i]] * buf[t.b[i]] : op == 6 ? buf[t.a[i]] / buf[t.b[i]] :
                 op == 7 ? -buf[t.a[i]] : op == 8 ? buf[t.a[i]]^Int(t.c[i]) :
                 buf[t.a[i]]^t.c[i]
    end
    buf[length(t.op)]
end

"""The identifiable rank of the package's own reduced law of `m` (its dependent and
numerator/denominator Exprs, complex-step derivatives in the log of each fitted parameter)
and the fitted count."""
function pkg_rank(m, rng; tol = 1e-9, h = 1e-20)
    em = ER.compile_mechanism(m)
    M = typeof(em)
    dep, indep = ER._dependent_param_exprs(M)
    num, den = ER._num_den_exprs(M)
    mets = collect(Symbol, ER.metabolites(em))
    names = Symbol[indep..., :Keq, :E_total, mets..., keys(dep)...]
    env = Dict(s => i for (i, s) in enumerate(names))
    dtapes = [(env[s], (t = Tape(); build!(t, ex, env); t)) for (s, ex) in dep]
    tn, td = Tape(), Tape()
    build!(tn, num, env); build!(td, den, env)
    np = length(indep)
    x = zeros(ComplexF64, length(names))
    for i in 1:np
        x[i] = 10.0^(2 * rand(rng) - 1)
    end
    x[np + 1] = 10.0^(2 * rand(rng) - 1)
    x[np + 2] = 1.0
    buf = zeros(ComplexF64, max(length(tn.op), length(td.op),
                                maximum((length(t.op) for (_, t) in dtapes); init = 1)))
    function v(y)
        for (i, t) in dtapes
            y[i] = run_tape(t, y, buf)
        end
        run_tape(tn, y, buf) / run_tape(td, y, buf)
    end
    rows = 2np + 6
    J = zeros(rows, np)
    for r in 1:rows
        for (k, s) in enumerate(mets)
            x[np + 2 + k] = 10.0^(2 * rand(rng) - 1)
        end
        for j in 1:np
            y = copy(x)
            y[j] *= cis(h)
            J[r, j] = imag(v(y)) / h
        end
        mx = maximum(abs, view(J, r, :); init = 0.0)
        mx > 0 && (J[r, :] ./= mx)
    end
    sv = np == 0 ? Float64[] : svdvals(J)
    s1 = isempty(sv) ? 1.0 : max(sv[1], eps())
    rk = count(>(tol * s1), sv)
    rk, np, rk == 0 ? 0.0 : sv[rk] / s1, rk < length(sv) ? sv[rk + 1] / s1 : 0.0
end

# ─── the pipeline rounds ─────────────────────────────────────────────────────────────────

"""
One batch as `_process_batch` treats it: drop structures already in `seen` (counted), and
those in `full_seen` that the full run would have produced in an earlier batch (counted
apart), run the package stage on the rest, and mark one representative per `eq_hash` not
in `memo` as fitted (`memo` grows). Returns the records of the fresh mechanisms.
"""
function process_batch!(mechs, seen::Set{UInt64}, memo::Set{UInt64}, full_seen;
                        parent = nothing, rnd = 0)
    recs = Dict{Symbol, Any}[]
    n_seen = n_full = 0
    for (i, m) in enumerate(mechs)
        h = hash(m)
        if h in seen
            n_seen += 1
            continue
        end
        push!(seen, h)
        if full_seen !== nothing && h in full_seen
            n_full += 1
            continue
        end
        r = Dict{Symbol, Any}(:round => rnd, :mech => m, :hash => h,
                              :parent => parent === nothing ? 0 : parent[i])
        merge!(r, Dict(pairs(pkg_stage(m))))
        push!(recs, r)
    end
    reps = Dict{UInt64, Int}()
    n_memo = 0
    for (i, r) in enumerate(recs)
        r[:fitted] = false
        r[:status] === :ok || continue
        if r[:eq_hash] in memo
            n_memo += 1
            continue
        end
        get!(reps, r[:eq_hash], i)
    end
    for i in values(reps)
        recs[i][:fitted] = true
    end
    union!(memo, keys(reps))
    recs, (n_in = length(mechs), n_seen = n_seen, n_full = n_full, n_memo = n_memo)
end
