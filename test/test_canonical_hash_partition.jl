# ABOUTME: Permanent regression test for _canonical_rate_eq_hash —
# ABOUTME: asserts determinism + partition stability vs frozen class counts.

using Test
using EnzymeRates

@testset "canonical-hash partition stability" begin
    # Representative reactions exercising the canonicalizer's edge cases:
    # - uni_uni: trivial structural equivalence
    # - bi_bi:   substituted-into-v ties across multiple kinetic groups
    # ter-ter intentionally omitted — `rate_equation` compilation is
    # extremely slow for mechanisms with >~30 enzyme forms (CLAUDE.md
    # "Known Issues"), and the canonical hasher invokes that path per
    # candidate. The bi-bi enumeration already covers every structural
    # symmetry the canonicalizer collapses.
    test_reactions = [
        ("uni_uni", @enzyme_reaction(begin
            substrates: S[C]
            products:   P[C]
        end)),
        ("bi_bi", @enzyme_reaction(begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end)),
    ]

    # Expected partition sizes per reaction. These literals were measured
    # empirically at Stage 6.2 (after the regex hasher was retired and the
    # struct-based hasher's partition matched the regex hasher's on these
    # fixtures per the Stage 6.0 equivalence test that gated the deletion).
    # If these counts ever change in a future commit, the canonical hasher's
    # equivalence classes have shifted — investigate before merging.
    expected_n_classes = Dict(
        "uni_uni" => 1,
        "bi_bi"   => 23,
    )

    for (label, reaction) in test_reactions
        # init_mechanisms only — skip expand_mechanisms. The init level
        # already produces multiple structurally-equivalent variants
        # (mirror-step orderings, kinetic-group renumberings) that
        # exercise the canonicalizer's collapse rules. expand_mechanisms
        # adds variants at higher param counts whose canonicalizer
        # behavior is the same modulo size, at exponential compile cost.
        all_mechs = EnzymeRates.init_mechanisms(reaction)

        new_buckets = Dict{UInt64, Vector{Int}}()
        for (i, m) in enumerate(all_mechs)
            em = EnzymeRates.compile_mechanism(m)
            h = EnzymeRates._canonical_rate_eq_hash(em)
            push!(get!(new_buckets, h, Int[]), i)
            # Determinism: same input, same hash across invocations.
            @test EnzymeRates._canonical_rate_eq_hash(em) === h
        end

        @test length(new_buckets) == expected_n_classes[label]
    end
end
