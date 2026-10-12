# ABOUTME: Prints, as tagged step lists, a few enumerated allosteric ping-pong mechanisms that carry
# ABOUTME: residual forms and prune a non-:OnlyA group in the inactive state (census examples).
using EnzymeRates, Random, Printf
const ER = EnzymeRates

rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end
side(f, ms) = join([string(ER.name(f)); [string(ER.name(m)) for m in ms]], " + ")
steplines(am) = [side(ER.from_species(s), ER.consumed(s)) * (ER.is_equilibrium(s) ? " ⇌ " : " <--> ") *
                 side(ER.to_species(s), ER.released(s)) * " :: " * string(ER.cat_allo_states(am)[g])
                 for (g, grp) in enumerate(ER.steps(am)) for s in grp]
allo = unique(filter(m -> m isa ER.AllostericMechanism, ER.expand_mechanisms(ER.init_mechanisms(rxn), rxn)))
hasres(am) = any(!isempty(ER.residual(f)) for g in ER.steps(am) for s in g
                 for f in (ER.from_species(s), ER.to_species(s)))
npruned(am) = count(isempty, ER._inactive_groups(am)) - count(==(:OnlyA), ER.cat_allo_states(am))
pick = filter(am -> hasres(am) && npruned(am) > 0, allo)
println("allosteric ping-pong: $(length(allo)); with residual forms and a pruned non-:OnlyA group: $(length(pick))")
for am in pick[[1, end]]
    em = ER.compile_mechanism(am)
    fp = ER.fitted_params(em)
    prm = NamedTuple{(fp..., :Keq, :E_total)}((Tuple(p === :L ? 1.0 : 0.7 for p in fp)..., 5.0, 1.0))
    free = unique(ER.name(f) for g in ER.steps(am) for s in g for f in (ER.from_species(s), ER.to_species(s))
                  if isempty(ER.bound(f)) && isempty(ER.residual(f)))
    println("\nfree forms $(free); inactive forms kept: ",
            sort!(unique(string(ER.name(f)) for g in ER._inactive_groups(am) for s in g
                         for f in (ER.from_species(s), ER.to_species(s)))))
    foreach(l -> println("    ", l), steplines(am))
    @printf("  v(A=2, B=3, P=Q=0) = %.5g, kcat = %.5g\n",
            ER.rate_equation(em, (A = 2.0, B = 3.0, P = 0.0, Q = 0.0), prm), ER._kcat_forward(em, prm))
end
