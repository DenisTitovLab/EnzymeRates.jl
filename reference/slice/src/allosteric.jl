# ABOUTME: The allosteric rate law (MWC formulation 1): both conformations' laws from one
# ABOUTME: step list and one solve, cross-weighted at the flip form, regulator sites, kcat.

"""
An allosteric rate law over `vars`: `active.vars` (concentrations, independents of both
conformations, `:Keq`), then the regulator constants (`regulator_names` order) and `:L`.
With `QA, DA` from `active`, `QI, DI` from `inactive` (`0` and `1` when the inactive
conformation is the flip form alone), `WA = Π_s sites[s][1]^sites[s][3]` and `WI` alike:

    v = E_total·H·(fI^n·QA·DA^(n-1)·WA + L·fA^n·QI·DI^(n-1)·WI) /
                  (fI^n·DA^n·WA + L·fA^n·DI^n·WI)

where `fA`, `fI` are the flip form's weights in the two laws over their common factor.
"""
struct AlloLaw{T <: Integer}
    vars::Vars
    H::Poly{T}
    active::Law{T}
    inactive::Union{Nothing, Law{T}}
    fA::Poly{T}
    fI::Poly{T}
    sites::Vector{Tuple{Poly{T}, Poly{T}, Int}}
    n::Int
end

"""
    allo_law(m, sol = solve(m), mask = all edges) -> AlloLaw

The MWC law of the edges in `mask`: the active law of `mask`, the inactive law of the edges
of `inactive_mask(m)` in `mask` (deletion plus pruning; `admit` says when it is the limit),
the flip form's weight in each, each regulatory site's binding polynomial per conformation.
The flip form is `m.forms[1]`, its segment's root, so its weight is its node's tree sum.
"""
function allo_law(m::Mechanism, sol::Solution = solve(m), mask = trues(length(m.edges)))
    any(k -> mask[k] && m.edges[k].a == 1, eachindex(mask)) ||
        error("the conformation flips at $(form_text(m, 1)), which none of the steps holds")
    imask = inactive_mask(m) .& mask
    LA = law(m, sol, mask)
    LI = any(imask) ? law(m, sol, imask; inactive = true) : nothing
    if LI !== nothing && typeof(LA) != typeof(LI)
        LA, LI = law(BigInt, m, sol, mask, false), law(BigInt, m, sol, imask, true)
    end
    P, vars, regs = typeof(LA.H), deepcopy(LA.vars), regulator_names(m)
    foreach(((_, ns),) -> foreach(x -> var!(vars, x), ns), regs)
    var!(vars, :L)
    flip(L) = L.trees[1] * L.shift
    fA, fI = cross_weights(flip(LA), LI === nothing ? one(P) : flip(LI))
    K = Dict(regs)
    sites = map(m.sites) do s
        W = [one(P), one(P)]
        for (l, st) in zip(s.ligands, s.states), (i, absent) in ((1, :OnlyI), (2, :OnlyA))
            st === absent && continue
            k = vars.index[K[l][i == 1 ? 1 : end]]
            addto!(W[i], mono(vars.index[l] => 1, k => -1), 1)
        end
        (W[1], W[2], s.multiplicity)
    end
    AlloLaw(vars, LA.H, LA, LI, fA, fI, sites, m.multiplicity)
end

"""The flip form's weights `fA`, `fI` over their common factor: their integer and monomial
contents' common part, and their polynomial part when it is the same (then both are
monomials, and `1` when the weights are equal)."""
function cross_weights(fA::Poly{T}, fI::Poly{T}) where {T}
    (cA, uA, rA), (cI, uI, rI) = content(fA), content(fI)
    rA == rI && (rA = rI = nothing)
    g = gcd(cA, cI)
    u = mono((v => min(mexp(uA, v), mexp(uI, v)) for v in unique(first.([uA; uI])))...)
    w(c, mu, r) = (p = term(T, mmul(mu, mpow(u, -1)), div(c, g)); r === nothing ? p : p * r)
    w(cA, uA, rA), w(cI, uI, rI)
end

"""`(QI, DI)` of `A`: the inactive law's `Q` and `D`, or `0` and `1` when the inactive
conformation is the flip form alone."""
inactive_terms(A::AlloLaw{T}) where {T} =
    A.inactive === nothing ? (zero(Poly{T}), one(Poly{T})) : (A.inactive.Q, A.inactive.D)

"""The rate law of `m`: its `Law`, or its `AlloLaw` when it is allosteric."""
rate_law(m::Mechanism, sol::Solution = solve(m)) =
    is_allosteric(m) ? allo_law(m, sol) : law(m, sol)

"""The edges of an allosteric mechanism's kcat law: its catalytic core (`catalytic_core`)
and the rapid-equilibrium steps between free forms, so the flip form stays when it is a
free isomer on a dead-end branch (its weight vanishes at saturation all the same)."""
kcat_mask(m::Mechanism) = catalytic_core(m) .|
    [e.re && is_free(m.forms[e.a]) && is_free(m.forms[e.b]) for e in m.edges]

"""
    kcat_terms(m, A::AlloLaw) -> (vars, terms)

kcat's candidates `(numerator, denominator)` of an allosteric mechanism over `A.vars`, `A`
its law on `kcat_mask(m)`. For each corner, a set of the names that are no reactant
(inhibitors, regulators; at most one ligand per regulatory site, as a saturated site holds
one, and the rate between two such corners is monotone) at saturation and the others at
zero: the limit of `v/E_total` with the products at zero as the substrates grow, then the
corner's names. Each factor of the law keeps its leading face, a branch of numerator or
denominator multiplies its factors' faces, and the limit is the ratio of the leading
branches' faces. An error names the problem when the rate grows without bound, its limit
depends on how the concentrations grow, or the denominator vanishes.
"""
function kcat_terms(m::Mechanism, A::AlloLaw{T}) where {T}
    nc, names, n = length(metabolites(m)), A.vars.names, A.n
    subs = [i for i in 1:nc if names[i] in m.substrates]
    QI, DI = inactive_terms(A)
    L = term(T, mono(A.vars.index[:L] => 1))
    WA, WI = [(s[1], s[3]) for s in A.sites], [(s[2], s[3]) for s in A.sites]
    num = [[(A.H, 1), (A.fI, n), (A.active.Q, 1), (A.active.D, n - 1), WA...],
           [(A.H, 1), (L, 1), (A.fA, n), (QI, 1), (DI, n - 1), WI...]]
    den = [[(A.fI, n), (A.active.D, n), WA...], [(L, 1), (A.fA, n), (DI, n), WI...]]
    pure(sat) = all(s -> count(i -> names[i] in s.ligands, sat) <= 1, m.sites)
    corners = filter(pure, all_subsets([i for i in 1:nc if !is_reactant(m, names[i])]))
    terms = map(corners) do sat
        deg(mo) = (sum((e for (v, e) in mo if v in subs); init = zero(Exp)),
                   sum((e for (v, e) in mo if v in sat); init = zero(Exp)))
        function face(p)                       # (degree, top terms) with the rest at zero
            ts = filter(t -> all(((v, _),) -> v > nc || v in subs || v in sat, t.first),
                        p.terms)
            isempty(ts) && return nothing
            top = maximum(deg, keys(ts))
            top, Poly{T}(filter(t -> deg(t.first) == top, ts))
        end
        function lead(branches)                # the leading branches' degree and face
            fs = map(branches) do b
                parts = [(face(p), e) for (p, e) in b if e > 0]
                any(x -> x[1] === nothing, parts) && return nothing
                (foldl((x, y) -> x .+ y, (e .* f[1] for (f, e) in parts)),
                 prod(f[2] for (f, e) in parts for _ in 1:e))
            end
            fs = filter(!isnothing, fs)
            isempty(fs) && return nothing
            top = maximum(first, fs)
            top, sum(last(f) for f in fs if first(f) == top)
        end
        N, D = lead(num), lead(den)
        D === nothing && kcat_error("the rate is undefined (0/0) with the products at zero")
        N === nothing || N[1] <= D[1] ||
            kcat_error("the rate grows without bound with " *
                       join(names[N[1][1] > D[1][1] ? subs : sat], ", ") *
                       ", so the mechanism never saturates")
        (N === nothing || N[1] < D[1]) && return (zero(Poly{T}), one(Poly{T}))
        face_ratio(split_conc(N[2], A.vars), split_conc(D[2], A.vars), names[[subs; sat]])
    end
    A.vars, terms
end

"""kcat's candidates of `m` (`kcat_terms`): on its catalytic core, or for an allosteric
mechanism on the MWC law of `kcat_mask(m)`."""
kcat_candidates(m::Mechanism, sol::Solution = solve(m)) =
    is_allosteric(m) ? kcat_terms(m, allo_law(m, sol, kcat_mask(m))) : kcat_terms(m, sol)

"""kcat at the parameter values `params` (by name, as `parameters(m)`): the largest
candidate of `kcat_candidates`."""
function kcat(m::Mechanism, params::NamedTuple)
    vars, terms = kcat_candidates(m)
    T = typeof(float(params.Keq))
    x = T[vars.conc[i] ? one(T) : params[s] for (i, s) in enumerate(vars.names)]
    maximum(evaluate(n, x) / evaluate(d, x) for (n, d) in terms)
end
