# ABOUTME: Allosteric check of the admitted limit: a cone direction from component potentials,
# ABOUTME: certified on the joint system's rows, and the full oracle at small ε against deletion.
include(joinpath(@__DIR__, "lib.jl"))

"""Components of the forms joined by edges outside `:OnlyA` groups, and the arcs of the
`:OnlyA` plain bindings (unbound side → bound side), as `admit` builds them."""
function doors(m, fwd)
    kept = BitVector(m.tags[g] !== :OnlyA for g in m.group)
    comp = collect(eachindex(m.forms))
    top(i) = (while comp[i] != i; i = comp[i]; end; i)
    foreach(k -> kept[k] && (comp[top(m.edges[k].a)] = top(m.edges[k].b)), eachindex(kept))
    comp .= top.(eachindex(comp))
    arcs = [(comp[o.a], comp[o.b]) for o in (S.oriented(m, k, fwd) for k in findall(.!kept))
            if S.is_plain_binding(m, o)]
    comp, arcs
end

"""Potentials of the components: the longest arc path into each (sources at `src`), the flip
form's component at 0; `nothing` when the arcs close a cycle."""
function potentials(m, comp, arcs, src)
    cs = unique(comp)
    phi = Dict(c => (c == comp[1] ? 0 : src) for c in cs)
    for _ in cs, (u, v) in arcs
        phi[v] = max(phi[v], phi[u] + 1)
    end
    all(((u, v),) -> phi[v] > phi[u], arcs) && phi[comp[1]] == 0 ? phi : nothing
end

"""The ε-exponent of every column under the potentials `phi`, the speed exponent of each
:OnlyA RE chemistry edge, and whether the direction satisfies the joint system exactly."""
function direction(m, cols, A, fwd, comp, phi)
    a, speed = zeros(Rational{BigInt}, length(cols)), Dict{Int, Tuple{Int, Bool}}()
    for (i, c) in enumerate(cols)
        c.inactive && m.tags[m.group[c.edge]] === :OnlyA || continue
        o = S.oriented(m, c.edge, fwd)
        d = phi[comp[o.b]] - phi[comp[o.a]]
        bind = S.is_plain_binding(m, o)
        a[i] = o.re ? (S.is_binding(o) ? -d : d) : bind ? (c.slot == 1 ? d : 0) :
               1 + max(c.slot == 1 ? d : -d, 0)
        o.re && !bind && (speed[c.edge] = (1 + max(-d, 0), S.is_binding(o)))
    end
    a, speed, iszero(A * a)
end

records = [r for r in load_records() if S.is_allosteric(r.c.mech) &&
           any(==(:OnlyA), r.c.mech.tags)]
rng = Xoshiro(20261014)
rows = []
t0 = time()
for r in records
    m, p, c = r.c.mech, r.p, r.c
    msg = S.admit(m)
    msg !== nothing && !occursin("reach the free enzyme", msg) && continue   # flip or G1
    sol = S.solve(m)
    A = S.allo_law(m, sol)
    cols, Amat, b = S.constraint_system(m)
    n = length(cols)
    M = Rational{BigInt}.([Amat b])
    piv = S.rref!(M, n)
    logkeq = log(Oracle.logu(rng))
    ell = zeros(BigFloat, n)
    free = setdiff(1:n, piv)
    foreach(f -> ell[f] = log(Oracle.logu(rng)), free)
    for (row, q) in enumerate(piv)
        ell[q] = M[row, n + 1] * logkeq - sum((M[row, f] * ell[f] for f in free); init = big(0))
    end
    maximum(abs, Amat * ell - b * logkeq) < 1e-60 || error("$(p.id): base point")
    fwd = S.orient(m)
    comp, arcs = doors(m, fwd)
    names = S.column_names(m, sol)
    regs = [x for (_, ns) in S.regulator_names(m) for x in ns]
    pars = (; (names[i] => exp(ell[i]) for i in sol.indep)...,
            (x => Oracle.logu(rng) for x in regs)..., L = Oracle.logu(rng), Keq = exp(logkeq))
    conc = Oracle.concs_at_ratio(p, rng, pars.Keq * 30)
    vdel = allo_value(A, conc, pars)
    for src in (msg === nothing ? (0,) : (1, -1))
        phi = potentials(m, comp, arcs, src)
        phi === nothing && continue
        a, speed, certified = direction(m, cols, Amat, fwd, comp, phi)
        errs = map((big(10)^-15, big(10)^-30)) do eps
            val = [exp(ell[i] + a[i] * log(eps)) for i in 1:n]
            idx = Dict(q => i for (i, q) in enumerate(cols))
            col(k, s, ina) = val[idx[S.Col(k, s, ina && m.tags[m.group[k]] !== :EqualAI)]]
            edge(k, ina) = [col(k, s, ina) for s in 1:(m.edges[k].re ? 1 : 2)]
            kA = S.pop_constants(c, [edge(k, false) for k in eachindex(m.edges)])
            kI, reI = Vector{BigFloat}[], Bool[]
            for (j, s) in enumerate(p.steps)
                k = c.edge_of[j]
                x = !haskey(speed, k) ? edge(k, true) :
                    (q = col(k, 1, true); q = speed[k][2] ? 1 / q : q;
                     [q * eps^speed[k][1], eps^speed[k][1]])
                push!(kI, S.pop_flip(c, fwd, j, x))
                push!(reI, s.re && !haskey(speed, k))
            end
            rn = Dict(S.regulator_names(m))
            regA = Dict{Tuple{Int, Symbol}, Any}((si, l) => st === :OnlyI ? nothing :
                        getfield(pars, rn[l][1]) for (si, site) in enumerate(p.sites)
                        for (l, st) in zip(site.ligands, site.states))
            regI = Dict{Tuple{Int, Symbol}, Any}((si, l) => st === :OnlyA ? nothing :
                        getfield(pars, rn[l][end]) for (si, site) in enumerate(p.sites)
                        for (l, st) in zip(site.ligands, site.states))
            lo = minimum(min(v, 1 / v) for (k, e) in enumerate(m.edges) if e.re
                         for v in (edge(k, false)..., edge(k, true)...); init = big(1))
            lam = big(10)^40 * maximum(maximum, kA) * maximum(maximum, kI) / lo
            vl = v_flip(p, pop_forms(c, p)[1], kA, kI, reI, trues(length(p.steps)),
                        1:length(p.species), pars.L, regA, regI, conc;
                        n = p.cat_multiplicity, lam)
            relerr(vl, vdel)
        end
        push!(rows, (set = r.set, id = p.id, admitted = msg === nothing, src, certified,
                     e15 = errs[1], e30 = errs[2]))
    end
end
println("mechanisms with :OnlyA groups that pass the flip rule and G1: ",
        length(unique(x -> (x.set, x.id), rows)))
println("directions certified exactly on the joint system's rows: ",
        count(x -> x.certified, rows), " of ", length(rows))
foreach(x -> println("  not certified: ", x), first(filter(x -> !x.certified, rows), 10))
worst = sort(filter(x -> x.admitted, rows); by = x -> -x.e15)
foreach(x -> println("  slowest: ", x), first(worst, 5))
adm = filter(x -> x.admitted, rows)
@printf("admitted: %d; |v(ε) - v_deletion| / |v_deletion| worst at ε = 1e-15: %.1e, at ε = 1e-30: %.1e\n",
        length(adm), maximum(x -> x.e15, adm), maximum(x -> x.e30, adm))
rej = filter(x -> !x.admitted, rows)
for set in unique(x.set for x in rej)
    rs = filter(x -> x.set == set, rej)
    ids = unique(x.id for x in rs)
    differ = count(id -> any(x -> x.id == id && x.e30 > 1e-6, rs), ids)
    @printf("island-rejected %-16s %3d mechanisms: the limit differs from deletion (> 1e-6 at ε = 1e-30) along at least one of the two directions in %d\n",
            set, length(ids), differ)
    for id in ids
        any(x -> x.id == id && x.e30 > 1e-6, rs) && continue
        println("    no difference found: ", id, " ", [(x.src, x.e15, x.e30) for x in rs if x.id == id])
    end
end
for x in filter(x -> x.id in ("allo_terbi/131", "allo_terbi/130", "expand_bibi/308"), rows)
    println("  ", x)
end
@printf("time %.0f s\n", time() - t0)
