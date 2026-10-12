# ABOUTME: Synthetic edge cases for Engine.bare_rate_law that the populations lack, checked against the oracle.
# ABOUTME: Disconnected graph, RE catalytic cycle, no free form, two reactive blocks, a cyclic dead-end branch.
const C = joinpath(@__DIR__, "..", "..", "common")
include(joinpath(C, "Oracle.jl")); include(joinpath(C, "Engine.jl"))
using Random, Printf
setprecision(BigFloat, Oracle.PREC)

S(name, mets...) = Pop.PSpecies([Pop.PMet(x, r) for (x, r) in mets], :E, Symbol[], Symbol[], name)
M(x, r = :S) = Pop.PMet(x, r)
step(f, t, c, r, re = false) = Pop.PStep(f, t, c, r, re)
function mech(id, subs, prods, species, steps; inh = Symbol[])
    Pop.PMech(; id, source = id, substrates = subs, products = prods, inhibitors = inh,
              metabolites = [subs; prods; inh], species, steps,
              groups = [[j] for j in eachindex(steps)])
end

cases = [
    # E ⇌ ES ⇌ EP ⇌ E all rapid: the RE steps close a catalytic cycle.
    mech("re_catalytic_cycle", [:S], [:P],
         [S(:E), S(:ES, (:S, :S)), S(:EP, (:P, :P))],
         [step(1, 2, [M(:S)], [], true), step(2, 3, [], [], true), step(1, 3, [M(:P, :P)], [], true)]),
    # Two components that share no form.
    mech("disconnected", [:S], [:P],
         [S(:E), S(:ES, (:S, :S)), S(:F), S(:FP, (:P, :P))],
         [step(1, 2, [M(:S)], []), step(3, 4, [M(:P, :P)], [])]),
    # Q stays bound throughout: no unbound residual-free form (dfree = 1).
    mech("no_free_form", [:A], [:P],
         [S(:EQ, (:Q, :Inh)), S(:EAQ, (:A, :S), (:Q, :Inh)), S(:EPQ, (:P, :P), (:Q, :Inh))],
         [step(1, 2, [M(:A)], []), step(2, 3, [], []), step(1, 3, [M(:P, :P)], [])]; inh = [:Q]),
    # Two catalytic routes through different complexes that share only E: blocks {E,ES}, {E,FS}.
    mech("two_reactive_blocks", [:S], [:P],
         [S(:E), S(:ES, (:S, :S)), S(:FS, (:S, :S))],
         [step(1, 2, [M(:S)], []), step(2, 1, [], [M(:P, :P)]),
          step(1, 3, [M(:S)], []), step(3, 1, [], [M(:P, :P)])]),
    # A square of SS inhibitor bindings hanging off E (a dead-end branch that holds a cycle).
    mech("cyclic_dead_end_branch", [:S], [:P],
         [S(:E), S(:ES, (:S, :S)), S(:EI, (:I, :Inh)), S(:EJ, (:J, :Inh)),
          S(:EIJ, (:I, :Inh), (:J, :Inh))],
         [step(1, 2, [M(:S)], []), step(2, 1, [], [M(:P, :P)]),
          step(1, 3, [M(:I, :Inh)], []), step(1, 4, [M(:J, :Inh)], []),
          step(3, 5, [M(:J, :Inh)], []), step(4, 5, [M(:I, :Inh)], [])]; inh = [:I, :J]),
    # One cycle that turns over twice: E → ES → F → FS → E takes up 2 S and gives off 2 P.
    mech("cycle_turns_twice", [:S], [:P],
         [S(:E), S(:ES, (:S, :S)), S(:F), S(:FS, (:S, :S))],
         [step(1, 2, [M(:S)], []), step(2, 3, [], [M(:P, :P)]),
          step(3, 4, [M(:S)], []), step(4, 1, [], [M(:P, :P)])]),
    # The same branch with rapid-equilibrium inhibitor bindings plus one SS dead-end leaf.
    mech("re_branch_plus_ss_leaf", [:S], [:P],
         [S(:E), S(:ES, (:S, :S)), S(:EI, (:I, :Inh)), S(:ESI, (:S, :S), (:I, :Inh))],
         [step(1, 2, [M(:S)], []), step(2, 1, [], [M(:P, :P)]),
          step(1, 3, [M(:I, :Inh)], [], true), step(2, 4, [M(:I, :Inh)], [])]; inh = [:I]),
]

vals_of(L, m, k, c) = (v = ones(BigFloat, length(L.vars.names));
    for (j, vs) in enumerate(L.stepvars), (i, x) in enumerate(vs); v[x] = k[j][i]; end;
    for x in m.metabolites; v[L.vars.index[x]] = c[x]; end; v)
lawv(L, vals) = Poly.evaluate(L.num, vals) / Poly.evaluate(L.den, vals)
function speed_free(L, m, k, c)
    v0 = lawv(L, vals_of(L, m, k, c))
    [j for (j, s) in enumerate(m.steps) if !s.re &&
     abs(lawv(L, vals_of(L, m, [i == j ? 3 .* q : q for (i, q) in enumerate(k)], c)) - v0) <=
     big(10)^-60 * abs(v0)]
end

rng = Random.Xoshiro(11)
for m in cases
    L = try
        Engine.bare_rate_law(m)
    catch e
        println(rpad(m.id, 24), " ERROR: ", sprint(showerror, e))
        continue
    end
    kc, _, _ = Oracle.consistent_constants(m, rng; Keq = 2.5)
    kb = [s.re ? kc[j] : [Oracle.logu(rng), Oracle.logu(rng)] for (j, s) in enumerate(m.steps)]
    c = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)
    err = Float64(abs(lawv(L, vals_of(L, m, kb, c)) / Oracle.steady_state_v(m, kb, c) - 1))
    s = L.stats
    @printf("%-24s oracle rel %.1e  blocks %s  cycles %d  max|n| %d  den %d  num %d  dfree %s\n", m.id,
            err, s[:block_sizes], s[:cycles], s[:max_n], s[:den_terms], s[:num_terms],
            Poly.nterms(L.dfree))
    println(" "^25, "dead-end prediction ", s[:dead_end], "  speed-free: thermo ",
            speed_free(L, m, kc, c), " bare ", speed_free(L, m, kb, c))
end
