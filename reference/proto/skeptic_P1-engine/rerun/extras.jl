# ABOUTME: Two design checks on the plain populations: (A) the package's reduced numerator divided exactly by
# ABOUTME: (1 − Γ/Keq) with a positive quotient; (B) speed-free SS steps as fitted phantoms at the group level.
const C = joinpath(@__DIR__, "..", "..", "common")
include(joinpath(C, "Oracle.jl")); include(joinpath(C, "Engine.jl"))
using Printf, Serialization
setprecision(BigFloat, Oracle.PREC)
const PR = Poly.P{Rational{Int}}
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
fromsym(pp, V) = PR(Dict{Poly.Mono, Rational{Int}}(
    Poly.mono((Poly.var!(V, s) => e for (s, e) in mono)...) => c for (mono, c) in pp))

"""Exact quotient of `p` by `1 − x` for a monomial `x` holding variable `kv` to the power −1
(chain division: terms equal up to powers of `x` form one Laurent polynomial in `x`), or
`nothing` when `1 − x` does not divide `p`."""
function divide_one_minus(p::PR, x::Poly.Mono, kv::Int)
    chains = Dict{Poly.Mono, Vector{Pair{Int, Rational{Int}}}}()
    for (t, c) in p.t
        e = Poly.mexp(t, kv)
        push!(get!(chains, Poly.mmul(t, Poly.mpow(x, e)), Pair{Int, Rational{Int}}[]), -e => c)
    end
    q = PR()
    for (key, terms) in chains
        sort!(terms; by = first)
        acc = 0 // 1
        for (k, (i, c)) in enumerate(terms)
            acc += c
            nxt = k < length(terms) ? terms[k + 1].first : i + 1
            for j in i:(nxt - 1)
                k == length(terms) && j == i && continue
                Poly.addto!(q.t, Poly.mmul(key, Poly.mpow(x, j)), acc)
            end
        end
        iszero(acc) || return nothing
    end
    q
end

rows = deserialize(joinpath(@__DIR__, "results.jls"))
byid = Dict(r[:id] => r for r in rows)
println("== (A) reduced numerator = (1 − Γ/Keq)·Q with Q > 0 (package pkg_num, dependents substituted)")
for pop in POPS
    ms = [m for m in Pop.load_population(pop) if !m.allosteric]
    ndiv = npos = 0; qterms = nterms = 0
    for m in ms
        V = Poly.Vars()
        num = fromsym(m.pkg_num, V)
        dep = Dict(Poly.var!(V, k) => Poly.mono((Poly.var!(V, s) => Int(e) for (s, e) in mono)...)
                   for (k, mono) in m.dependent_monomials)
        num = Poly.subst(num, v -> get(dep, v, nothing))
        kv = Poly.var!(V, :Keq)
        x = Poly.mono((Poly.var!(V, p) => 1 for p in m.products)...,
                      (Poly.var!(V, s) => -1 for s in m.substrates)..., kv => -1)
        q = divide_one_minus(num, x, kv)
        q === nothing && (println("   not divisible: ", m.id, "  ", m.source); continue)
        ndiv += 1
        npos += all(>(0), values(q.t))
        qterms += Poly.nterms(q); nterms += Poly.nterms(num)
        all(>(0), values(q.t)) || println("   quotient has a non-positive coefficient: ", m.id)
    end
    @printf("  %-13s n=%3d  divisible %3d  quotient all-positive %3d  terms: numerator %5d → quotient %5d\n",
            pop, length(ms), ndiv, npos, nterms, qterms)
end

println("\n== (B) speed-free SS steps (law unchanged when a step's kf and kr scale together) in the package's fitted set")
for pop in POPS
    ms = [m for m in Pop.load_population(pop) if !m.allosteric]
    nmech = nphantom = ngroups_checked = 0
    for m in ms
        r = byid[m.id]
        isempty(r[:speed_free_thermo]) && continue
        s = m.samples[1]
        base = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.params)
        concs = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in s.concs)
        function v_at(vals)
            for (k, mono) in m.dependent_monomials
                vals[k] = prod((vals[s]^e for (s, e) in mono); init = big(1))
            end
            all = merge(vals, concs)
            Oracle.eval_poly(m.pkg_num, all) / Oracle.eval_poly(m.pkg_den, all)
        end
        v0 = v_at(copy(base))
        found = 0
        for pair in unique(Tuple(m.step_names[j]) for j in r[:speed_free_thermo])
            ngroups_checked += 1
            vals = copy(base)
            fitted = [n for n in pair if n in m.fitted_params]
            length(fitted) == 2 || continue
            for n in fitted; vals[n] *= 3; end
            abs(v_at(vals) - v0) <= big(10)^-40 * abs(v0) && (found += 1)
        end
        found > 0 && (nmech += 1; nphantom += found)
    end
    @printf("  %-13s mechanisms with a fitted kf/kr pair the law cannot separate: %3d (%3d pairs; %3d speed-free groups checked)\n",
            pop, nmech, nphantom, ngroups_checked)
end
