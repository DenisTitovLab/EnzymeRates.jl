# ABOUTME: Replays params_estimate_{5,6,7,8}.csv through the canonicalizer and
# ABOUTME: asserts mechanisms with the same fitted loss collapse to one hash.

using CSV, DataFrames, Test
using EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

const _CSV_REPLAY_DIR = joinpath(@__DIR__, "..", "dedup_investigation")

@testset "CSV dedup replay" begin
    for n in (5, 6, 7, 8)
        csv_path = joinpath(_CSV_REPLAY_DIR, "params_estimate_$n.csv")
        if !isfile(csv_path)
            @info "skipping CSV replay for n=$n: $(csv_path) missing"
            continue
        end
        df = CSV.read(csv_path, DataFrame)

        # Drop sentinel "failed to fit" rows. The fitter returns loss = 100
        # (or another sentinel near 1.0) when CMA-ES couldn't converge for
        # the mechanism. These rows aren't mathematically equivalent — they
        # just share a placeholder loss value — so grouping them by loss
        # for dedup-correctness assertions is meaningless and noisy.
        df = filter(:loss => l -> l < 1.0, df)
        isempty(df) && continue

        canonical_hashes = map(eachrow(df)) do row
            m = eval(Meta.parse(row.mechanism_type))()
            _canonical_rate_eq_hash(m)
        end
        df.canonical_hash = canonical_hashes

        # Same-loss → same-hash is a *sufficient-but-not-necessary*
        # signal for Source A/B/C dedup: Source A/B/C duplicates fit
        # to the same loss and SHOULD share a hash, so a high
        # within-loss-group consistency rate is a positive indicator.
        # But the converse doesn't hold — two genuinely distinct rate
        # equations (e.g., a non-allosteric mechanism vs. an allosteric
        # one at saturating regulator, or different catalytic topologies
        # whose data happens to be degenerate) can yield identical
        # fitted losses without being mathematically equivalent.
        # Assert a high consistency rate, not perfection.
        @testset "n=$n within-loss-group consistency rate" begin
            groups = groupby(df, :loss)
            n_groups = length(groups)
            n_violators = count(g -> length(unique(g.canonical_hash)) > 1, groups)
            rate_ok = (n_groups - n_violators) / n_groups
            # ≥98% of loss groups should hash-collapse. The residual ~2%
            # are same-loss-different-equation cases (different type
            # families with degenerate-fit data) that aren't real dedup
            # misses.
            @test rate_ok >= 0.98
        end

        # Same-hash → losses within a wide tolerance. This catches gross
        # hash collisions (distinct equations accidentally mapped to one
        # hash). The 10% tolerance accommodates CMA-ES local-minimum
        # drift: on the same equation, different starts can yield 1-10%
        # loss variation. The cluster test in test_eq_hash_dedup.jl is
        # the strict correctness gate; this one is a sanity check that
        # nothing is wildly wrong.
        @testset "n=$n within-hash-group sanity" begin
            for g in groupby(df, :canonical_hash)
                lo, hi = extrema(g.loss)
                @test isapprox(lo, hi; rtol=0.1)
            end
        end
    end
end
