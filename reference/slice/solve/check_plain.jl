# ABOUTME: Solve-stage checks (a) and (b) on plain mechanisms: the relation space against the
# ABOUTME: package's, dependent choices, package-sample rebuilds and BigFloat detailed balance.
include(joinpath(@__DIR__, "lib.jl"))

records = [r for r in load_records() if !S.is_allosteric(r.c.mech)]
pkg = [r for r in records if r.p.ok && !isempty(r.p.fitted_params)]
const SHARED = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                "init_bibi_pp", "expand_bibi")
@printf("plain records: %d in all, %d package-built with package outputs\n",
        length(records), length(pkg))

# ─── (a) relation space against the package's ──────────────────────────────────────────
"""The package's relations as rows over the solve's columns plus a last `log Keq` column:
ties between steps that share a name, dependents, and Pass-2 folds."""
function pkg_rows(r, cols)
    c, p, m = r.c, r.p, r.c.mech
    fwd, n = S.orient(m), length(cols)
    idx = Dict(x => i for (i, x) in enumerate(cols))
    form, rows = Dict{Symbol, Vector{Rational{BigInt}}}(), Vector{Rational{BigInt}}[]
    for (j, ns) in enumerate(p.step_names), (q, nm) in enumerate(ns)
        slot, sg = pop_slot(c, fwd, j, q)
        v = zeros(Rational{BigInt}, n + 1)
        v[idx[S.Col(c.edge_of[j], slot, false)]] = sg
        haskey(form, nm) ? push!(rows, v - form[nm]) : (form[nm] = v)
    end
    for (d, mono) in p.dependent_monomials
        v = copy(form[d])
        for (s, e) in mono
            s === :Keq ? (v[n + 1] -= e) : (v .-= e .* form[s])
        end
        push!(rows, v)
    end
    foreach(((a, t),) -> push!(rows, form[a] - form[t]), p.pass2_rename)
    rows, form
end

"""The package's name of column `col` (through a Pop step on its edge)."""
function pkg_name(r, fwd, col)
    for (j, k) in enumerate(r.c.edge_of)
        k == col.edge || continue
        for q in eachindex(r.p.step_names[j])
            first(pop_slot(r.c, fwd, j, q)) == col.slot && return r.p.step_names[j][q]
        end
    end
    error("no Pop step on edge $(col.edge)")
end

nsame, ndiff, same_deps, diff_deps, diff_hal, diff_weg = 0, 0, 0, 0, 0, 0
tally = Dict{String, Vector{Int}}()     # set => [records, same space, same dependents]
count_mismatch, maxdepdiff = 0, 0
examples = String[]
for r in pkg
    m = r.c.mech
    cols, A, b = S.constraint_system(m)
    mine = Rational{BigInt}.([A -b])
    prow, _ = pkg_rows(r, cols)
    P = isempty(prow) ? zeros(Rational{BigInt}, 0, length(cols) + 1) :
        permutedims(reduce(hcat, prow))
    ra, rb, rab = rank_exact(mine), rank_exact(P), rank_exact([mine; P])
    t = get!(tally, r.set, [0, 0, 0])
    t[1] += 1
    if ra == rb == rab
        global nsame += 1
        t[2] += 1
    else
        global ndiff += 1
        push!(examples, "relation space differs: $(r.p.id) mine $ra pkg $rb joint $rab")
    end
    sol = S.solve(m)
    length(sol.indep) == length(r.p.fitted_params) || (global count_mismatch += 1)
    fwd = S.orient(m)
    mydeps = Set(pkg_name(r, fwd, sol.cols[d.col]) for d in sol.deps
                 if d.label !== :kinetic_group)
    pkdeps = Set([first.(r.p.dependent_monomials); first.(r.p.pass2_rename)])
    if mydeps == pkdeps
        global same_deps += 1
        t[3] += 1
    else
        global diff_deps += 1
        global maxdepdiff = max(maxdepdiff, length(setdiff(mydeps, pkdeps)))
        myhal = Set(pkg_name(r, fwd, sol.cols[d.col]) for d in sol.deps if d.label === :haldane)
        pkhal = Set(d for (d, mono) in r.p.dependent_monomials if any(x -> x[1] === :Keq, mono))
        myhal == pkhal || (global diff_hal += 1)
        setdiff(mydeps, myhal) == setdiff(pkdeps, pkhal) || (global diff_weg += 1)
    end
end
@printf("(a) relation space equal to the package's: %d of %d (differ %d)\n", nsame,
        length(pkg), ndiff)
@printf("    fitted count differs from length(fitted_params): %d\n", count_mismatch)
@printf("    dependent set (package names, folds included) same %d, differs %d \
        (Haldane choice differs %d, non-Haldane choice differs %d; at most %d names)\n",
        same_deps, diff_deps, diff_hal, diff_weg, maxdepdiff)
for set in sort!(collect(keys(tally)))
    t = tally[set]
    @printf("    %-14s %4d records: same relation space %4d, same dependent set %4d\n", set,
            t...)
end
sh = [tally[x] for x in SHARED if haskey(tally, x)]
@printf("    the 806 of the shared populations: %d records, same space %d, same dependents \
        %d, different dependents %d\n", sum(first, sh), sum(t -> t[2], sh),
        sum(t -> t[3], sh), sum(first, sh) - sum(t -> t[3], sh))
foreach(println, first(examples, 10))

# ─── (b+) package samples rebuilt through all_parameters ───────────────────────────────
"""BigFloat values of every package name at sample `s`: fitted, dependents from their
monomials, folds."""
function pkg_values(p, s)
    v = Dict{Symbol, BigFloat}(n => BigFloat(x) for (n, x) in s.params)
    for (d, mono) in p.dependent_monomials
        v[d] = prod((BigFloat(v[x])^BigFloat(e) for (x, e) in mono); init = big(1))
    end
    for (a, t) in p.pass2_rename
        v[a] = v[t]
    end
    v
end
worst_rebuild, nrebuild = big(0), 0
for r in pkg, s in r.p.samples
    m, c, p = r.c.mech, r.c, r.p
    vals = pkg_values(p, s)
    y = [BigFloat[vals[n] for n in ns] for ns in p.step_names]
    x = S.edge_constants(c, y)
    sol = S.solve(m)
    names = S.column_names(m, sol)
    pars = (; (names[i] => x[sol.cols[i].edge][sol.cols[i].slot] for i in sol.indep)...,
            Keq = BigFloat(vals[:Keq]), E_total = big(1))
    got = edge_values(m, S.all_parameters(m, pars))
    for k in eachindex(x), q in eachindex(x[k])
        global worst_rebuild = max(worst_rebuild, abs(got[k][q] - x[k][q]) / x[k][q])
    end
    global nrebuild += 1
end
@printf("(b+) package samples rebuilt from my independents through all_parameters: %d \
        samples, worst relative error of a per-edge constant %.1e\n", nrebuild,
        Float64(worst_rebuild))

# ─── (b) detailed balance through all_parameters ───────────────────────────────────────
rng = Xoshiro(20261010)
worst_rel, worst_abs, nbal, bad = big(0), big(0), 0, String[]
byset = Dict{String, Vector{Float64}}()
rejected = String[]
for r in records
    m = r.c.mech
    try
        S.solve(m)
    catch err
        push!(rejected, "$(r.set) $(r.p.id): $(sprint(showerror, err))")
        continue
    end
    for draw in 1:2
        pars = random_params(m, rng)
        x = edge_values(m, S.all_parameters(m, pars))
        k = S.pop_constants(r.c, x)
        v0 = Oracle.steady_state_v(r.p, k, Oracle.concs_at_ratio(r.p, rng, pars.Keq))
        v1 = Oracle.steady_state_v(r.p, k, Oracle.concs_at_ratio(r.p, rng, pars.Keq / 10))
        rel = abs(v0) / abs(v1)
        global worst_rel = max(worst_rel, rel)
        global worst_abs = max(worst_abs, abs(v0))
        global nbal += 1
        push!(get!(byset, r.set, Float64[]), Float64(rel))
        rel > 1e-30 && push!(bad, "$(r.set) $(r.p.id) rel $(Float64(rel))")
    end
end
@printf("(b) detailed balance: %d draws on %d solved plain records; worst |v(Keq)|/|v(Keq/10)| \
        %.1e, worst |v(Keq)| %.1e; above 1e-30: %d\n", nbal, length(records) - length(rejected),
        Float64(worst_rel), Float64(worst_abs), length(bad))
for set in sort!(collect(keys(byset)))
    @printf("    %-18s %4d draws, worst %.1e\n", set, length(byset[set]), maximum(byset[set]))
end
foreach(println, first(bad, 10))
println("solve rejected $(length(rejected)) plain records:")
foreach(x -> println("    ", x), rejected)
