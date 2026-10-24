# ABOUTME: Builds the bi-bi allosteric population (init + one expand round) and serializes it,
# ABOUTME: so a fresh process can read it back the way identify's workers do.
using EnzymeRates, Serialization
const rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
end
seeds = EnzymeRates.init_mechanisms(rxn)
kids = EnzymeRates.expand_mechanisms(seeds, rxn)
pop = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[seeds; kids]
serialize(ARGS[1], pop)
println("wrote ", length(pop), " mechanisms to ", ARGS[1])
