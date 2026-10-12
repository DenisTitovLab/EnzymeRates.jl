# ABOUTME: P5b extra checks: the collapse-test assertions, determinism under step/group/form
# ABOUTME: permutations, and whether each free-enzyme weight stays positive at zero concentrations.
using Random, Printf
include(joinpath(@__DIR__, "Allo.jl"))

# ─── (e) the assertions of test/test_allosteric_collapse.jl on my outputs ──

hw = Dict(m.source => m for m in Pop.load_population("handwritten") if m.allosteric)
src(s) = "handwritten:test_allosteric_collapse.jl:" * s
"""Whether changing fitted `p` (×1.3 vs ×5) moves my rate at one random point."""
function moves(md, p)
    rng = Xoshiro(3)
    fit = Dict{Symbol, BigFloat}(n => big(0.4) + 2big(rand(rng)) for n in md.names.fitted)
    fit[:Keq] = big(3.0)
    c = Dict(x => big(0.4) + 2big(rand(rng)) for x in md.m.metabolites)
    v(x) = (f = copy(fit); f[p] = x; Allo.rate(md, Allo.all_values(md, f), c))
    !isapprox(v(big(1.3)), v(big(5.0)); rtol = 1e-12)
end
"""The right side of a name of `md` (`(Keq exponent, [name => exponent])`)."""
rhs(md, n) = md.names.value[n]
checks = Pair{String, Bool}[]
let md = Allo.derive(hw[src("uni([:NonequalAI,:EqualAI,:EqualAI])")]), t = Allo.law_text(md)
    push!(checks, "RE S:NonequalAI cat,P:EqualAI -> K_I collapses to K_A" =>
          !(:K_I_ES_to_E_S in md.names.fitted) &&
          get(md.names.folds, :K_I_ES_to_E_S, nothing) === :K_A_ES_to_E_S && !occursin("k_I_", t))
end
let md = Allo.derive(hw[src("uni([:NonequalAI,:EqualAI,:NonequalAI])")]), t = Allo.law_text(md)
    free = filter(p -> startswith(string(p), "K_I_"), md.names.fitted)
    push!(checks, "two NonequalAI bindings + EqualAI catalysis -> exactly one free K_I that moves v" =>
          length(free) == 1 && !occursin("k_I_", t) && moves(md, only(free)))
end
let md = Allo.derive(hw[src("uni([:NonequalAI,:NonequalAI,:EqualAI])")])
    push!(checks, "NonequalAI catalysis -> K_I and a k_I both free" =>
          :K_I_ES_to_E_S in md.names.fitted && any(p -> startswith(string(p), "k_I_"), md.names.fitted))
end
let md = Allo.derive(hw[src("ro_bibi(st[3]=:NonequalAI)")])
    push!(checks, "RO bi-bi one inner NonequalAI edge -> K_I not fitted and equal to K_A" =>
          !(:K_I_EAB_to_EB_A in md.names.fitted) && rhs(md, :K_I_EAB_to_EB_A) == rhs(md, :K_A_EAB_to_EB_A))
end
let md = Allo.derive(hw[src("uni([:NonequalAI,:EqualAI,:EqualAI]; re=(false,false,false))")]),
    t = replace(Allo.law_text(md), " " => "")
    push!(checks, "SS binding + EqualAI catalysis -> affinity collapses, speed free and moves v" =>
          !(:k_I_ES_to_E_S in md.names.fitted) && :k_I_E_S_to_ES in md.names.fitted &&
          occursin("k_I_ES_to_E_S=", t) && moves(md, :k_I_E_S_to_ES))
end
let md = Allo.derive(hw[src("uni([:NonequalAI,:NonequalAI,:EqualAI]; re=(false,false,false))")])
    push!(checks, "SS binding + NonequalAI catalysis -> both k_I free" =>
          :k_I_E_S_to_ES in md.names.fitted && :k_I_ES_to_E_S in md.names.fitted)
end
let md = Allo.derive(hw[src("uni([:NonequalAI,:OnlyA,:EqualAI])")]), t = replace(Allo.law_text(md), " " => "")
    push!(checks, "dead-I NonequalAI binding -> K_I free, not collapsed, moves v" =>
          :K_I_ES_to_E_S in md.names.fitted && !occursin("K_I_ES_to_E_S=K_A_ES_to_E_S", t) &&
          moves(md, :K_I_ES_to_E_S))
end
println("(e) collapse assertions:")
foreach(((k, v),) -> println("   ", v ? "PASS " : "FAIL ", k), checks)

# ─── determinism: permute forms, steps and groups ────────────────────────

"""`m` with its forms, flat steps and groups reordered at random (the content unchanged)."""
function permuted(m::Pop.PMech, rng)
    σ = randperm(rng, length(m.species)); inv = invperm(σ)        # new form k = old σ[k]
    gperm = randperm(rng, length(m.groups))
    order = reduce(vcat, [shuffle(rng, m.groups[g]) for g in gperm])  # new step i = old order[i]
    newidx = invperm(order)
    s(j) = (st = m.steps[j]; Pop.PStep(inv[st.from], inv[st.to], st.consumed, st.released, st.re))
    Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates, products = m.products,
              inhibitors = m.inhibitors, regulators = m.regulators, metabolites = m.metabolites,
              species = m.species[σ], steps = [s(j) for j in order],
              groups = [sort!([newidx[j] for j in m.groups[g]]) for g in gperm],
              allosteric = true, cat_multiplicity = m.cat_multiplicity,
              group_tags = m.group_tags[gperm], sites = m.sites,
              step_names_A = m.step_names_A[order], step_names_I = m.step_names_I[order])
end
lawkey(md) = [Allo.Poly.symbolic(p, md.T) for p in (md.A..., md.I...)]
"""All names reduced to independents: each law with every dependent substituted."""
function reduced(md)
    T2 = Allo.Vars()
    f(v) = (n = md.T.names[v]; haskey(md.names.value, n) ?
            (val = md.names.value[n]; Allo.mono([Allo.var!(T2, s) => Int(e) for (s, e) in
                                                [(:Keq, val[1]); val[2]] if e != 0]...)) :
            Allo.mono(Allo.var!(T2, n) => 1))
    [Allo.Poly.symbolic(Allo.subst(p, f), T2) for p in (md.A..., md.I...)]
end
rng = Xoshiro(11)
for name in (isempty(ARGS) ? ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi") : ARGS)
    pop = filter(m -> m.allosteric, Pop.load_population(name))
    nperm, fail_fit, fail_dep, fail_law, fail_red = 0, 0, 0, 0, 0
    examples = String[]
    for m in pop
        md = Allo.derive(m)
        for _ in 1:3
            mp = permuted(m, rng); mq = Allo.derive(mp); nperm += 1
            fail_fit += md.names.fitted != mq.names.fitted
            fail_dep += md.names.value != mq.names.value
            lk = lawkey(md) != lawkey(mq)
            fail_law += lk
            fail_red += reduced(md) != reduced(mq)
            lk && length(examples) < 3 && push!(examples, m.id)
        end
    end
    @printf("determinism %s: %d permutations; fitted differ %d, name values differ %d, named law differs %d, reduced law differs %d %s\n",
            name, nperm, fail_fit, fail_dep, fail_law, fail_red, isempty(examples) ? "" : "e.g. " * join(examples, ", "))
end

# ─── free-enzyme weights: positive at zero concentrations? ───────────────

for name in (isempty(ARGS) ? ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi") : ARGS)
    pop = filter(m -> m.allosteric, Pop.load_population(name))
    tab = Dict{Any, Int}()
    for m in pop
        md = Allo.derive(m)
        conc = Set(md.T.index[x] for x in m.metabolites if haskey(md.T.index, x))
        kind(d) = !any(v in conc for mo in keys(d.t) for (v, _) in mo) ?
                  (Allo.nterms(d) == 1 ? :mono : :poly) :
                  any(mo -> !any(v in conc for (v, _) in mo), keys(d.t)) ? :met_pos : :met_zero
        k = (Allo.norm_case(md), kind(md.A[3]), kind(md.I[3]))
        tab[k] = get(tab, k, 0) + 1
    end
    println("weights $name (case, d_A, d_I): ", join(("$k => $v" for (k, v) in sort!(collect(tab); by = string ∘ first)), "; "))
end
