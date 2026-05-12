# ABOUTME: Replays params_estimate_{5,6,7,8}.csv through the canonicalizer and
# ABOUTME: asserts mechanisms with the same fitted loss collapse to one hash.

using CSV, DataFrames, Test
using EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

"""
Classify a `fitted_params` symbol into its canonical kind for
shape-comparison within hash-equivalent mechanism buckets. The kind
is invariant under rep-step renaming (e.g. `:k9f`, `:k10f`, `:k11f`
all map to `:kf`); the underlying step-index is what varies across
topologically-distinct-but-equation-equivalent mechanisms.
"""
function _fp_kind(s::Symbol)
    str = string(s)
    is_T = endswith(str, "_T")
    base = is_T ? str[1:end-2] : str
    kind = if startswith(base, "K") && length(base) > 1 && isdigit(base[2])
        :K
    elseif startswith(base, "k") && length(base) > 1 && isdigit(base[2])
        endswith(base, "f") ? :kf : endswith(base, "r") ? :kr : :other
    elseif s == :L
        :L
    elseif startswith(str, "K_")
        :K_reg
    else
        :other
    end
    is_T ? Symbol(kind, :_T) : kind
end

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

        # Same-hash → same fitted_params SHAPE (count + kind multiset).
        # The strict assertion `length(unique(g.fp)) == 1` was originally
        # written to catch the c8a3302 Pass-2-absorbed-symbol leak, which
        # manifested as a fitted-params COUNT mismatch within a hash
        # bucket. After c8a3302, counts within each hash bucket are
        # consistent — but the NAMES vary because `parameters(m)` names
        # each kinetic group after its representative step (the first
        # source-order step in that group), and topologically-distinct-
        # but-equation-equivalent mechanisms have different rep-step
        # indices for the same group. Example: three mechs all hashing
        # to the same canonical equation can expose
        # `(:K1,:K2,:K3,:K4,:k9f)`, `(:K1,:K2,:K3,:K4,:k10f)`, and
        # `(:K1,:K2,:K3,:K4,:k11f)` — same shape, different rep-step
        # number. The fitter handles this via the param-projection logic
        # in `identify_rate_equation.jl:286` (cached params get mapped
        # onto each target spec's own `fitted_params` keys at fit time).
        # We assert the invariants the canonicalizer DOES guarantee:
        # same param count and same kind-multiset.
        @testset "n=$n hash-equivalent mechanisms share fitted_params shape" begin
            mechs_by_row = [eval(Meta.parse(row.mechanism_type))()
                            for row in eachrow(df)]
            df.fp = [EnzymeRates.fitted_params(m) for m in mechs_by_row]
            for g in groupby(df, :canonical_hash)
                # 1. Same param count — catches the c8a3302 leak
                #    regression (extra symbol slipping into indep).
                @test length(unique(length.(g.fp))) == 1
                # 2. Same multiset of param kinds — catches qualitative
                #    mismatches (e.g. one mech exposing kNr while
                #    another exposes only kNf for the same equation).
                kinds = [sort(_fp_kind.(fp)) for fp in g.fp]
                @test length(unique(kinds)) == 1
            end
        end
    end
end
