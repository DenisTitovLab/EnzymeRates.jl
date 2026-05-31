@testset "parameters API symmetry" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end
    m_mono = EnzymeRates.EnzymeMechanism(
        first(EnzymeRates.init_mechanisms(rxn)))
    full_mono = parameters(m_mono, Full)
    reduced_mono = parameters(m_mono, Reduced)
    @test :E_total in full_mono
    @test :E_total in reduced_mono
    @test :Keq in reduced_mono
    @test :Keq ∉ full_mono
    @test length(full_mono) >= length(reduced_mono)

    rxn_allo = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    base = first(EnzymeRates.init_mechanisms(rxn_allo))
    cat_allo_states = fill(:NonequalAI, length(EnzymeRates.kinetic_groups(base)))
    site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)),
        cat_allo_states, 2, [site])
    m_allo = EnzymeRates.AllostericEnzymeMechanism(am)
    full_allo = parameters(m_allo, Full)
    reduced_allo = parameters(m_allo, Reduced)
    @test :L in full_allo
    @test :E_total in full_allo
    @test :Keq ∉ full_allo
    @test any(occursin("_I_", string(p)) for p in full_allo)
    @test any(occursin("reg", string(p)) for p in full_allo)
    @test :L in reduced_allo
    @test :Keq in reduced_allo
    @test :E_total in reduced_allo
end

@testset "added field accessors: cat_allo_states" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    cas = fill(:NonequalAI, length(EnzymeRates.kinetic_groups(base)))
    site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)), cas, 2, [site])
    @test EnzymeRates.cat_allo_states(am) == cas
end
