# ABOUTME: Task 4: rate_equation and kcat values of deterministically enumerated allosteric
# ABOUTME: bi-bi mechanisms, written to a TSV for comparison between two package checkouts.
using EnzymeRates, Random, Printf
const ER = EnzymeRates

out = ARGS[1]
const N1, N2 = 200, 100            # picked from generation 1 and from generation 2

# The bi-bi reaction of test/test_mechanism_enumeration.jl (`bi_bi_rxn`).
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

key(m) = sprint(show, ER.compile_mechanism(m))
pick(v, k) = v[unique(round.(Int, range(1, length(v); length = min(k, length(v)))))]

t0 = time()
init = ER.init_mechanisms(rxn)
gen1 = ER.expand_mechanisms(init, rxn)
allo1 = sort!(unique(filter(m -> m isa ER.AllostericMechanism, gen1)); by = key)
parents = pick(allo1, 40)
gen2 = ER.expand_mechanisms(parents, rxn)
seen = Set(allo1)
allo2 = sort!(unique(filter(m -> m isa ER.AllostericMechanism && !(m in seen), gen2)); by = key)
println("init $(length(init)), generation 1: $(length(gen1)) children, $(length(allo1)) allosteric; ",
        "generation 2 from $(length(parents)) allosteric parents: $(length(gen2)) children, ",
        "$(length(allo2)) new allosteric; $(round(time() - t0; digits = 1)) s")
chosen = [pick(allo1, N1); pick(allo2, N2)]

nfree(am) = length(unique(ER.name(f) for g in ER.steps(am) for s in g
                          for f in (ER.from_species(s), ER.to_species(s))
                          if isempty(ER.bound(f)) && isempty(ER.residual(f))))
tags(am) = join(sort(unique(String.(ER.cat_allo_states(am)))), "+")

open(out, "w") do io
    println(io, join(["i", "key_hash", "n_free", "tags", "n_pruned_groups", "fitted_params",
                      "eqstring_hash", "v_generic", "v_products0", "v_low", "kcat"], '\t'))
    for (i, am) in enumerate(chosen)
        em = ER.compile_mechanism(am)
        fp = ER.fitted_params(em)
        rng = MersenneTwister(1000 + i)
        lu(lo, hi) = exp(log(lo) + (log(hi) - log(lo)) * rand(rng))
        vals = Tuple(p === :L ? lu(0.2, 5.0) : lu(0.1, 10.0) for p in fp)
        Keq = lu(0.5, 20.0)
        prm = NamedTuple{(fp..., :Keq, :E_total)}((vals..., Keq, 1.0))
        mets = Tuple(ER.metabolites(em))
        conc(f) = NamedTuple{mets}(Tuple(f(x) for x in mets))
        c1 = conc(_ -> lu(0.1, 10.0))
        c2 = conc(x -> x in (:P, :Q) ? 0.0 : lu(0.1, 10.0))
        c3 = conc(_ -> lu(0.01, 0.1))
        v = [ER.rate_equation(em, c, prm) for c in (c1, c2, c3)]
        kc = try ER._kcat_forward(em, prm) catch err; NaN end
        npruned = count(isempty, ER._inactive_groups(am)) -
                  count(==(:OnlyA), ER.cat_allo_states(am))
        println(io, join([string(i), string(hash(key(am))), string(nfree(am)), tags(am),
                          string(npruned), join(fp, ","),
                          string(hash(ER.rate_equation_string(em))),
                          repr.(v)..., repr(kc)], '\t'))
        i % 50 == 0 && (println("$i done, $(round(time() - t0; digits = 1)) s"); flush(stdout))
    end
end
println("wrote $(length(chosen)) mechanisms to $out in $(round(time() - t0; digits = 1)) s")
