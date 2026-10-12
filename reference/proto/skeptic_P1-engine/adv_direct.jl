# ABOUTME: Adversarial plain mechanisms built directly as Pop.PMech (no package): free form that is not its
# ABOUTME: segment root, speed-free predicate gaps, and the Hill engine's cycle-rank limit; checked against the oracle.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf
setprecision(BigFloat, Oracle.PREC)

"""A form `name` in conformation `conf` holding `(metabolite, role)` pairs."""
F(name, mets...; conf = :E) = Pop.PSpecies([Pop.PMet(x, r) for (x, r) in mets], conf, Symbol[],
                                           Symbol[], name)
M(x, r = :S) = Pop.PMet(x, r)
"""A mechanism from forms and steps `(from_name, to_name, consumed, released, re)`; one group per step."""
function mech(id, subs, prods, forms, steps; inh = Symbol[])
    idx = Dict(f.name => i for (i, f) in enumerate(forms))
    st = [Pop.PStep(idx[a], idx[b], c, r, re) for (a, b, c, r, re) in steps]
    Pop.PMech(; id, source = id, substrates = subs, products = prods, inhibitors = inh,
              metabolites = [subs; prods; inh], species = forms, steps = st,
              groups = [[j] for j in eachindex(st)],
              step_names = [s.re ? [Symbol("K", j)] : [Symbol("kf", j), Symbol("kr", j)]
                            for (j, s) in enumerate(st)])
end

cases = [
    # Free form E is first in species order; the RE segment {E, D} is rooted at D ("D" < "E").
    mech("D1_free_form_not_root", [:S], [:P],
         [F(:E), F(:ES, (:S, :S)), F(:D; conf = :D)],
         [(:E, :ES, [M(:S)], [], false), (:ES, :E, [], [M(:P, :P)], false),
          (:D, :E, [], [], true)]),
    # Two catalytic conformations joined only by an SS isomerization (a bridge).
    mech("D2_bridge_between_reactive_parts", [:S], [:P],
         [F(:E), F(:ES, (:S, :S)), F(:F; conf = :F), F(:FS, (:S, :S); conf = :F)],
         [(:E, :ES, [M(:S)], [], false), (:ES, :E, [], [M(:P, :P)], false),
          (:F, :FS, [M(:S)], [], false), (:FS, :F, [], [M(:P, :P)], false),
          (:E, :F, [], [], false)]),
    # All turnover inside one RE segment (SS chemistry self-loop) plus one SS dead-end leaf.
    mech("D3_single_segment_core_plus_ss_leaf", [:S], [:P],
         [F(:E), F(:ES, (:S, :S)), F(:EP, (:P, :P)), F(:EI, (:I, :Inh))],
         [(:E, :ES, [M(:S)], [], true), (:ES, :EP, [], [], false), (:EP, :E, [], [M(:P, :P)], true),
          (:E, :EI, [M(:I, :Inh)], [], false)]; inh = [:I]),
    # A non-reactive square of isomerizations between two catalytic conformations.
    mech("D4_nonreactive_square_between_reactive_parts", [:S], [:P],
         [F(:E), F(:ES, (:S, :S)), F(:F; conf = :F), F(:FS, (:S, :S); conf = :F),
          F(:G; conf = :G), F(:H; conf = :H)],
         [(:E, :ES, [M(:S)], [], false), (:ES, :E, [], [M(:P, :P)], false),
          (:F, :FS, [M(:S)], [], false), (:FS, :F, [], [M(:P, :P)], false),
          (:E, :G, [], [], false), (:G, :F, [], [], false), (:E, :H, [], [], false),
          (:H, :F, [], [], false)]),
    # 2 + 20 parallel SS edges between two RE segments: cycle rank 21.
    mech("D5_parallel_rank_21", [:S], [:P],
         [F(:E); F(:ES, (:S, :S)); [F(Symbol("EI", k), (Symbol("I", k), :Inh)) for k in 1:20];
          [F(Symbol("ESI", k), (Symbol("I", k), :Inh), (:S, :S)) for k in 1:20]],
         [(:E, :ES, [M(:S)], [], false); (:ES, :E, [], [M(:P, :P)], false);
          [(:E, Symbol("EI", k), [M(Symbol("I", k), :Inh)], [], true) for k in 1:20];
          [(:ES, Symbol("ESI", k), [M(Symbol("I", k), :Inh)], [], true) for k in 1:20];
          [(Symbol("EI", k), Symbol("ESI", k), [M(:S)], [], false) for k in 1:20]];
         inh = [Symbol("I", k) for k in 1:20]),
]

for (n, m) in enumerate(cases)
    println("\n=== ", m.id)
    rng = Random.Xoshiro(500 + n)
    t = @elapsed L = try
        Engine.bare_rate_law(m)
    catch e
        println("    Engine ERROR: ", first(sprint(showerror, e), 200)); nothing
    end
    tl = @elapsed lean = EngineLean.bare_rate_law(m)
    V, sv, ln, ld, lf = lean
    @printf("    Engine %.3f s, EngineLean %.3f s\n", t, tl)
    # Oracle on the lean law (always available).
    vals(k, c) = (v = ones(BigFloat, length(V.names));
                  for (j, vs) in enumerate(sv), (i, x) in enumerate(vs); v[x] = k[j][i]; end;
                  for x in m.metabolites; v[V.index[x]] = c[x]; end; v)
    kb = Sk.bare_constants(m, rng); c = Sk.random_concs(m, rng)
    vl = Poly.evaluate(ln, vals(kb, c)) / Poly.evaluate(ld, vals(kb, c))
    @printf("    lean vs oracle rel %.1e\n", Sk.relerr(vl, Oracle.steady_state_v(m, kb, c)))
    i, root = Sk.free_and_root(m)
    x = Oracle.steady_state_populations(m, kb, c)
    share = Poly.evaluate(lf, vals(kb, c)) / Poly.evaluate(ld, vals(kb, c))
    @printf("    free form %s (segment root %s): dfree/den = %.6g, oracle x[%s] = %.6g, x[%s] = %.6g\n",
            m.species[i].name, m.species[root].name, Float64(share), m.species[i].name,
            Float64(x[i]), m.species[root].name, Float64(x[root]))
    L === nothing && continue
    kc, _, _ = Oracle.consistent_constants(m, rng; Keq = 2.0)
    println("    probe predicate ", L.stats[:dead_end], " | exact thermo/bare ", Sk.speed_free_exact(m))
    println("    numeric thermo ", Sk.speed_free_numeric(L, m, kc, c), " oracle thermo ",
            Sk.speed_free_oracle(m, kc, c))
    println("    numeric bare   ", Sk.speed_free_numeric(L, m, kb, c), " oracle bare   ",
            Sk.speed_free_oracle(m, kb, c))
end
