# ABOUTME: Solve-stage check (e) on allosteric mechanisms: phantom names, fold hazards, fitted
# ABOUTME: lists against the package's, package samples, two-state detailed balance, islands.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(PROTO, "p5a", "pkgvals.jl"))         # MWC oracle, package sample values

const POP1993 = ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
const POP848 = ("allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi")
records = [r for r in load_records() if S.is_allosteric(r.c.mech)]
const Λ = MWC.big1024(10)^40

"""Fundamental cycles of the subgraph of `keep` edges reachable from form `root`, as in
`Slice.cycles` (physical-forward signs, Keq multiples)."""
function subcycles(m, fwd, keep, root)
    E = length(m.edges)
    path = [zeros(Int, E) for _ in m.forms]
    seen, tree, grown = falses(length(m.forms)), falses(E), true
    seen[root] = true
    while grown
        grown = false
        for (k, e) in enumerate(m.edges)
            keep[k] && seen[e.a] != seen[e.b] || continue
            from, to, s = seen[e.a] ? (e.a, e.b, 1) : (e.b, e.a, -1)
            path[to] .= path[from]; path[to][k] += s
            seen[to] = tree[k] = grown = true
        end
    end
    names = unique!([x.name for x in m.mets])
    at = [searchsortedfirst(names, x.name) for x in m.mets]
    net = [x in m.substrates ? -1 : x in m.products ? 1 : 0 for x in names]
    i0 = findfirst(!=(0), net)
    zs, mults = Vector{Int}[], Int[]
    for (k, e) in enumerate(m.edges)
        keep[k] && !tree[k] && seen[e.a] || continue
        z = path[e.a] - path[e.b]; z[k] += 1
        tot = zeros(Int, length(names))
        for (j, f) in enumerate(m.edges), i in eachindex(m.mets)
            tot[at[i]] += z[j] * (S.count_at(f.rel, i) - S.count_at(f.cons, i))
        end
        push!(zs, [fwd[j] ? z[j] : -z[j] for j in 1:E])
        push!(mults, tot[i0] * net[i0])
    end
    zs, mults
end

"""Count of law independents when the inactive cycle rows come only from the law's own
inactive graph (`alt = :pruned`) or when `:OnlyA` members stay tied in the inactive state
(`alt = :tied`)."""
function alt_count(m, alt)
    cols, A, b = S.constraint_system(m)
    mask = S.inactive_mask(m)
    fwd = S.orient(m)
    index = Dict(c => i for (i, c) in enumerate(cols))
    ncyc = length(first(S.cycles(m, fwd)))
    rows = [A[i, :] for i in axes(A, 1)]
    rhs = copy(b)
    col(k, s) = index[S.Col(k, s, m.tags[m.group[k]] !== :EqualAI)]
    if alt === :pruned
        keep = ncyc + 1:2ncyc
        rows, rhs = rows[setdiff(eachindex(rows), keep)], rhs[setdiff(eachindex(rhs), keep)]
        zs, mults = subcycles(m, fwd, mask, findfirst(S.is_free, m.forms))
        for (z, c) in zip(zs, mults)
            row = zeros(Int, length(cols))
            for (k, zk) in enumerate(z)
                zk == 0 && continue
                o = S.oriented(m, k, fwd)
                row[col(k, 1)] += o.re && S.is_binding(o) ? -zk : zk
                o.re || (row[col(k, 2)] -= zk)
            end
            push!(rows, row); push!(rhs, c)
        end
    else
        lead = Dict{Tuple{Int, Int}, Int}()
        for (i, c) in enumerate(cols)
            c.inactive && m.tags[m.group[c.edge]] === :OnlyA || continue
            f = get!(lead, (Int(m.group[c.edge]), c.slot), i)
            f == i && continue
            row = zeros(Int, length(cols)); row[i], row[f] = 1, -1
            push!(rows, row); push!(rhs, 0)
        end
    end
    M = Rational{BigInt}.([permutedims(reduce(hcat, rows)) rhs])
    piv = S.rref!(M, length(cols))
    count(i -> !(i in piv) && (!cols[i].inactive || mask[cols[i].edge]), eachindex(cols))
end

stats = Dict{String, Dict{Symbol, Any}}()
st(set) = get!(() -> Dict{Symbol, Any}(:n => 0, :diff => Dict{Int, Int}(), :phantom => 0,
                                       :phantom_pkg => 0, :phantom_pkg_mechs => 0,
                                       :hazard => 0, :opp_pkg => 0, :onestate_pkg => 0,
                                       :opp_pkg_mechs => 0, :onestate_pkg_mechs => 0,
                                       :mask_ids => String[], :explain => Dict{Any, Int}(),
                                       :mask_diff => 0, :db => 0.0, :db_n => 0,
                                       :db_skip => 0, :same_eq => 0.0, :same_ne => Float64[],
                                       :pruned => 0, :tied => Dict{Int, Int}(),
                                       :errors => String[]), stats, set)
rng = Xoshiro(20261012)
for r in records
    m, p, c = r.c.mech, r.p, r.c
    t = st(r.set)
    t[:n] += 1
    sol = try
        S.solve(m)
    catch err
        push!(t[:errors], "$(p.id): $(sprint(showerror, err))")
        continue
    end
    names = S.column_names(m, sol)
    fitted = S.fitted_params(m)
    mask = S.inactive_mask(m)
    fwd = S.orient(m)
    d = length(fitted) - length(p.fitted_params)
    t[:diff][d] = get(t[:diff], d, 0) + 1
    # phantom names: fitted names that are not the full mechanism's declared constants
    act, ina = S.constant_names(m), S.constant_names(m; inactive = true)
    declared = Set{Symbol}([n for ns in act for n in ns])
    foreach(k -> m.tags[m.group[k]] === :EqualAI || union!(declared, ina[k]),
            eachindex(m.edges))
    union!(declared, (n for (_, ns) in S.regulator_names(m) for n in ns), (:L,))
    t[:phantom] += count(n -> !(n in declared), fitted) +
                   count(i -> sol.cols[i].inactive && !mask[sol.cols[i].edge], sol.indep)
    pdecl = Set{Symbol}([n for ns in [p.step_names_A; p.step_names_I] for n in ns])
    pph = count(n -> !(n in pdecl) && n !== :L && !endswith(String(n), "reg"),
                p.fitted_params)
    t[:phantom_pkg] += pph
    t[:phantom_pkg_mechs] += pph > 0
    # fold hazards: every name's columns are one independent with its folds, or dependents
    # with one right side
    law = [i for i in eachindex(sol.cols) if !sol.cols[i].inactive || mask[sol.cols[i].edge]]
    expr = Dict(i => ([i => 1 // 1], 0 // 1) for i in sol.indep)
    foreach(x -> expr[x.col] = (x.terms, x.keq), sol.deps)
    for n in unique(names[law])
        cs = [i for i in law if names[i] == n]
        length(unique(expr[i] for i in cs)) == 1 || (t[:hazard] += 1)
    end
    rA, rI = Dict(p.pass2_rename_A), Dict(p.pass2_rename_I)
    t[:opp_pkg] += count(((x, y),) -> get(rI, y, nothing) === x, collect(rA))
    syms(ps) = Set(s for q in ps for (mo, _) in q for (s, _) in mo)
    sA, sI = syms(p.state_polys_A), syms(p.state_polys_I)
    one = count(x -> x in sI, keys(rA)) + count(x -> x in sA, keys(rI))
    opp = any(((x, y),) -> get(rI, y, nothing) === x, collect(rA))
    t[:onestate_pkg] += one
    t[:onestate_pkg_mechs] += one > 0 && !opp
    t[:opp_pkg_mechs] += opp
    # what explains a count difference: package phantom names, package fold hazards
    key = (d, pph, opp || one > 0)
    t[:explain][key] = get(t[:explain], key, 0) + 1
    # the inactive law's edges against the package's pruning
    if any(j -> mask[c.edge_of[j]] != p.inactive_kept[j], eachindex(p.steps))
        t[:mask_diff] += 1
        push!(t[:mask_ids], "$(p.id) $(p.source) ok=$(p.ok)")
    end
    # the package's samples against my relations
    worst = big(0)
    for s in p.samples
        vals = pkg_values(p, s)
        kA, kI = pkg_step_constants(p, vals, :A), pkg_step_constants(p, vals, :I)
        kI = [kI[j] === nothing ? kA[j] : kI[j] for j in eachindex(kI)]
        xA, xI = S.edge_constants(c, kA), S.edge_constants(c, kI)
        x(i) = (q = sol.cols[i]; (q.inactive ? xI : xA)[q.edge][q.slot])
        for dp in sol.deps
            rhs = prod((x(i)^BigFloat(e) for (i, e) in dp.terms);
                       init = vals[:Keq]^BigFloat(dp.keq))
            worst = max(worst, abs(x(dp.col) / rhs - 1))
        end
        # an :EqualAI edge of the inactive law shares its active constants
        for k in eachindex(m.edges), q in eachindex(xA[k])
            mask[k] && m.tags[m.group[k]] === :EqualAI || continue
            worst = max(worst, abs(xI[k][q] / xA[k][q] - 1))
        end
    end
    d == 0 ? (t[:same_eq] = max(t[:same_eq], Float64(worst))) :
             push!(t[:same_ne], Float64(worst))
    # detailed balance of both states in the two-conformation oracle
    pars = random_params(m, rng)
    allp = S.all_parameters(m, pars)
    kA = S.pop_constants(c, edge_values(m, allp))
    kI = S.pop_constants(c, edge_values(m, allp; inactive = true, mask))
    act_I = [mask[c.edge_of[j]] for j in eachindex(p.steps)]
    keep = Set(f for (j, s) in enumerate(p.steps) if act_I[j] for f in (s.from, s.to))
    flipform = findfirst(S.is_free, m.forms)
    pf = MWC.free_form(p)
    j0 = findfirst(j -> c.edge_of[j] in findall(e -> flipform in (e.a, e.b), m.edges),
                   eachindex(p.steps))
    e0 = m.edges[c.edge_of[j0]]
    popflip = (e0.a == flipform) != c.flipped[j0] ? p.steps[j0].from : p.steps[j0].to
    if pf != popflip
        t[:db_skip] += 1
    else
        push!(keep, pf)
        regA, regI = Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
        rn = Dict(S.regulator_names(m))
        for (si, site) in enumerate(p.sites), (lig, state) in zip(site.ligands, site.states)
            ns = [getfield(allp, n) for n in rn[lig]]
            regA[(si, lig)] = state === :OnlyI ? nothing : ns[1]
            regI[(si, lig)] = state === :OnlyA ? nothing : ns[end]
        end
        reI = [s.re for s in p.steps]
        v(g) = MWC.v_factorized(p, kA, kI, reI, act_I, keep, allp.L, regA, regI,
                                Oracle.concs_at_ratio(p, rng, g);
                                n = p.cat_multiplicity, Λ = Λ)
        t[:db] = max(t[:db], Float64(abs(v(pars.Keq)) / abs(v(pars.Keq / 10))))
        t[:db_n] += 1
    end
    # island cycles and tied :OnlyA classes
    t[:pruned] += alt_count(m, :pruned) != length(sol.indep)
    dt = alt_count(m, :tied) - length(sol.indep)
    t[:tied][dt] = get(t[:tied], dt, 0) + 1
end

function report(sets, title)
    ts = [stats[s] for s in sets if haskey(stats, s)]
    isempty(ts) && return
    sm(k) = sum(t -> t[k], ts)
    diff = Dict{Int, Int}()
    foreach(t -> mergewith!(+, diff, t[:diff]), ts)
    tied = Dict{Int, Int}()
    foreach(t -> mergewith!(+, tied, t[:tied]), ts)
    ne = reduce(vcat, (t[:same_ne] for t in ts))
    @printf("== %s: %d mechanisms, %d solve errors\n", title, sm(:n),
            sum(t -> length(t[:errors]), ts))
    @printf("   fitted length minus the package's: %s\n", sort!(collect(diff)))
    @printf("   phantom names: mine %d; package %d names in %d mechanisms\n", sm(:phantom),
            sm(:phantom_pkg), sm(:phantom_pkg_mechs))
    @printf("   fold hazards: mine %d; package folds in opposite directions %d (%d \
            mechanisms), folded in one state but used by the other %d (%d further \
            mechanisms)\n",
            sm(:hazard), sm(:opp_pkg), sm(:opp_pkg_mechs), sm(:onestate_pkg),
            sm(:onestate_pkg_mechs))
    @printf("   inactive law edges differ from the package's pruning in %d\n",
            sm(:mask_diff))
    @printf("   package samples satisfy my relations: equal counts worst %.1e; unequal \
            counts %d, worst residual min %.1e max %.1e\n",
            maximum(t -> t[:same_eq], ts), length(ne), isempty(ne) ? NaN : minimum(ne),
            isempty(ne) ? NaN : maximum(ne))
    @printf("   two-state detailed balance: %d mechanisms, worst |v(Keq)|/|v(Keq/10)| \
            %.1e (%d skipped: flip form differs from the oracle's)\n", sm(:db_n),
            maximum(t -> t[:db], ts), sm(:db_skip))
    @printf("   island cycles change the count in %d; tying :OnlyA members changes it by \
            %s\n", sm(:pruned), sort!(collect(tied)))
    explain = Dict{Any, Int}()
    foreach(t -> mergewith!(+, explain, t[:explain]), ts)
    println("   (mine - package, package phantom names, package fold hazard) => mechanisms:")
    for (k, v) in sort!(collect(explain))
        k[1] == 0 && k[2] == 0 && !k[3] && continue
        println("      ", k, " => ", v)
    end
end
for s in sort!(collect(keys(stats)))
    report((s,), s)
    foreach(x -> println("      error ", x), first(stats[s][:errors], 5))
    foreach(x -> println("      inactive law differs from the package: ", x),
            first(stats[s][:mask_ids], 40))
end
report(POP1993, "the 1,993")
report(POP848, "the 848")
