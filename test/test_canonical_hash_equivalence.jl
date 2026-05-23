# ABOUTME: Verifies _canonical_rate_eq_hash agrees with the preserved
# ABOUTME: regex-based _canonical_rate_eq_hash_old on equivalence-class assignment.

using Test
using EnzymeRates

@testset "canonical-hash equivalence (old regex vs current)" begin
    # Representative reactions exercising the canonicalizer's edge cases:
    # - uni-uni: trivial structural equivalence
    # - bi-bi:   substituted-into-v ties across multiple kinetic groups
    # ter-ter intentionally omitted — `rate_equation` compilation is
    # extremely slow for mechanisms with >~30 enzyme forms (CLAUDE.md
    # "Known Issues"), and the canonical hasher invokes that path per
    # candidate. The bi-bi enumeration already covers every structural
    # symmetry the canonicalizer collapses.
    test_reactions = [
        @enzyme_reaction(begin
            substrates: S[C]
            products:   P[C]
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end),
    ]

    for reaction in test_reactions
        # init_mechanisms only — skip expand_mechanisms. The init level
        # already produces multiple structurally-equivalent variants
        # (mirror-step orderings, kinetic-group renumberings) that
        # exercise the canonicalizer's collapse rules. expand_mechanisms
        # adds variants at higher param counts whose canonicalizer
        # behavior is the same modulo size, at exponential compile cost.
        all_specs = EnzymeRates.init_mechanisms(reaction)

        old_buckets = Dict{UInt64, Vector{Int}}()
        new_buckets = Dict{UInt64, Vector{Int}}()
        for (i, spec) in enumerate(all_specs)
            em = EnzymeRates.compile_mechanism(spec)
            old_h = EnzymeRates._canonical_rate_eq_hash_old(em)
            new_h = EnzymeRates._canonical_rate_eq_hash(em)
            push!(get!(old_buckets, old_h, Int[]), i)
            push!(get!(new_buckets, new_h, Int[]), i)
        end

        # Equivalence-class agreement: same partition of specs into
        # buckets (modulo bucket key — what matters is which specs share
        # a bucket).
        old_partition = Set(Set(b) for b in values(old_buckets))
        new_partition = Set(Set(b) for b in values(new_buckets))

        @test old_partition == new_partition
    end
end
