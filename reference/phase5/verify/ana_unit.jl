# ABOUTME: Verifier unit-rule recheck on 200 laws the unit census did not use (new ter-ter init laws
# ABOUTME: and bi-bi round-2 allosteric laws): homogeneity, the spec rule, and N against Cleland's.
using Printf, Random
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
include(joinpath(SP, "slice/src/Slice.jl"))
const S = Slice
include(joinpath(VD, "sliceside.jl"))
include(joinpath(VD, "analib.jl"))
setprecision(BigFloat, 256)
T0 = time()
const U = Tuple{QQ, QQ}

"""Units (concentration, time exponents) of the law's variables, by name: each independent
from its step (slot 1 the physical-forward rate constant), concentrations, Keq."""
function var_units(m, sol)
    fwd, names = S.orient(m), S.column_names(m, sol)
    u = Dict{Symbol, U}(x => (QQ(1), QQ(0)) for x in S.metabolites(m))
    u[:Keq] = (QQ(length(m.products) - length(m.substrates)), QQ(0))
    for c in sol.indep
        col = sol.cols[c]
        o = S.oriented(m, col.edge, fwd)
        u[names[c]] = !o.re ? (QQ(-(col.slot == 1 ? S.nibsum(o.cons) : S.nibsum(o.rel))), QQ(-1)) :
                      S.is_binding(o) ? (QQ(1), QQ(0)) : (QQ(S.nibsum(o.rel) - S.nibsum(o.cons)), QQ(0))
    end
    u
end
munit(mo, vars, u) = foldl(((a, b), (v, e)) -> (a + e * u[vars.names[v]][1], b + e * u[vars.names[v]][2]),
                           mo; init = (QQ(0), QQ(0)))

"""Exact rank of the columns `cs` (vectors of rationals)."""
function qrank(cs)
    isempty(cs) && return 0
    M = reduce(hcat, cs)
    # rank via elimination: count pivots
    A = copy(M); r = 0
    for c in 1:size(A, 2)
        i = findfirst(k -> !iszero(A[k, c]), (r + 1):size(A, 1))
        i === nothing && continue
        i += r; r += 1
        A[[r, i], :] = A[[i, r], :]
        for k in (r + 1):size(A, 1)
            iszero(A[k, c]) || (A[k, :] .-= (A[k, c] // A[r, c]) .* A[r, :])
        end
    end
    r
end

"""The spec's unit rule (SPEC:490-509) on one conformation: `(:constant | :virtual_int |
:virtual_frac | :ratio, λ)` with λ as monomial => weight."""
function spec_rule(D, nc, subs)
    haskey(D, S.Mono()) && return :constant, [S.Mono() => QQ(1)]
    ex(mo) = [Int(S.mexp(mo, i)) for i in 1:nc]
    key(mo) = (sum(ex(mo)), count(i -> ex(mo)[i] != 0 && !(i in subs), 1:nc), -ex(mo))
    monos = sort!(collect(keys(D)); by = key)
    aff(mo) = QQ[1; ex(mo)]
    basis = S.Mono[]
    for mo in monos
        qrank(aff.([basis; [mo]])) > length(basis) && push!(basis, mo)
    end
    M = [reduce(hcat, aff.(basis)) QQ[1; zeros(Int, nc)]]
    λ = solve_exact(M)
    λ === nothing && return :ratio, [monos[1] => QQ(1)]
    (all(isinteger, λ) ? :virtual_int : :virtual_frac), [b => l for (b, l) in zip(basis, λ) if !iszero(l)]
end

function cleland(D, subs)
    whole = S.mono((i => 1 for i in subs)...)
    parts = [S.mmul(whole, S.mono(i => -1)) for i in subs]
    all(mo -> haskey(D, mo), [parts; [whole]]) || return nothing
    [[p => QQ(1) for p in parts]; whole => QQ(1 - length(subs))]
end

logN(λ, D, x) = sum(BigFloat(l) * log(S.evaluate(D[mo], x)) for (mo, l) in λ)

"""Check one conformation's (Q, D) over `vars`: homogeneity, the rule, Cleland."""
function check(m, Q, Dp, H, vars, u, rng)
    nc = length(S.metabolites(m))
    subs = [vars.index[x] for x in m.substrates]
    UD = unique([munit(mo, vars, u) for mo in keys(Dp.terms)])
    hom = length(UD) == 1
    UV = unique([munit(S.mmul(h, q), vars, u) .- UD[1] for h in keys(H.terms) for q in keys(Q.terms)])
    vok = isempty(Q.terms) || UV == [(QQ(0), QQ(-1))]
    D = S.split_conc(Dp, vars)
    rule, λ = spec_rule(D, nc, subs)
    # unit of N = Π c^λ (c the coefficient of mo: unit UD − a_mo), plus the ratio form's monomial
    uN = foldl((a, (mo, l)) -> a .+ l .* (UD[1] .- munit(mo, vars, u)), λ; init = (QQ(0), QQ(0)))
    rule === :ratio && (uN = uN .+ munit(λ[1].first, vars, u))
    dimless = (UD[1] .- uN) == (QQ(0), QQ(0))
    cle = rule === :constant ? nothing : cleland(D, subs)
    same = cle === nothing ? missing : begin
        all(1:2) do _
            x = draw(vars, rng)
            a, b = logN(λ, D, x), logN(cle, D, x)
            abs(a - b) <= big(10.0)^-60 * (1 + abs(a))
        end
    end
    (; hom, vok, rule, dimless, eligible = cle !== nothing, same, nterms = length(Dp.terms))
end

rows = Any[]
rng = Xoshiro(5)
for l in eachline(joinpath(VD, "data", get(ENV, "TT", "terter_new.tsv")))
    id, cx, raw = split(l, '\t'; limit = 3)
    m = mkmech(S, eval(Meta.parse(raw)))
    try
        sol = S.solve(m)
        L = S.law(m, sol)
        push!(rows, (:terter, id, check(m, L.Q, L.D, L.H, L.vars, var_units(m, sol), rng)))
    catch e
        push!(rows, (:terter, id, first(sprint(showerror, e), 120)))
    end
    length(rows) == parse(Int, get(ENV, "NTT", "150")) && break
end
@printf("[%4.0f s] ter-ter laws: %d\n", time() - T0, length(rows)); flush(stdout)
pl = [split(l, '\t') for l in eachline(joinpath(VD, "data", "phantom_pop.tsv"))]
pl = [f for f in pl if f[2] == "2" && f[5] == "1" && f[6] == "1"]
nallo = 0
for f in pl[randperm(Xoshiro(17), length(pl))]
    m = mkmech(S, eval(Meta.parse(f[9])))
    S.admit(m) === nothing || continue
    sol = S.solve(m)
    A = S.allo_law(m, sol)
    u = var_units(m, sol)
    push!(rows, (:alloA, f[1], check(m, A.active.Q, A.active.D, A.active.H, A.active.vars, u, rng)))
    A.inactive === nothing ||
        push!(rows, (:alloI, f[1], check(m, A.inactive.Q, A.inactive.D, A.inactive.H, A.inactive.vars, u, rng)))
    global nallo += 1
    nallo == 50 && break
end
for set in (:terter, :alloA, :alloI)
    rs = [r[3] for r in rows if r[1] === set]
    ok = [r for r in rs if !(r isa String)]
    @printf("%s: %d conformations, %d errors; homogeneous %d, v unit E/time %d, dimensionless %d; rules %s; Cleland-eligible %d, N = Cleland %d, N differs %d\n",
            set, length(rs), length(rs) - length(ok), count(r -> r.hom, ok), count(r -> r.vok, ok),
            count(r -> r.dimless, ok), string(sort(collect(Dict(k => count(r -> r.rule === k, ok) for k in unique(r.rule for r in ok))))),
            count(r -> r.eligible, ok), count(r -> r.same === true, ok), count(r -> r.same === false, ok))
    for r in rs
        r isa String && println("  error: ", r)
    end
end
@printf("done in %.0f s, maxrss %.0f MB\n", time() - T0, Sys.maxrss() / 2^20)
