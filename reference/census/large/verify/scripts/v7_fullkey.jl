# ABOUTME: Ter-bi only: keys of the full base tier (every seed, not the 2,000 sample), and how many
# ABOUTME: fitted expansion children repeat a law of a base mechanism the sample left out.
include(joinpath(@__DIR__, "vlib.jl"))

rxn = ter_bi_rxn
seeds = unique!(collect(ER.init_mechanisms(rxn)))
base, _ = ER._base_tier(seeds, rxn)
t0 = time()
keys_full = Dict{String, Int}()
nskip = nerr = 0
for (i, m) in enumerate(base)
    ER._eq_complexity(m) > 337 && (global nskip += 1; continue)
    k = try
        S.canonical(S.rate_law(to_slice(m))).key
    catch
        global nerr += 1
        continue
    end
    get!(keys_full, k, i)
end
@printf("full base tier: %d mechanisms, %d over the complexity filter, %d errors, %d distinct keys (%.0f s)\n",
        length(base), nskip, nerr, length(keys_full), time() - t0)
data = load_run(:ter_bi_rxn)
sample_keys = Set(r[:key] for r in data.rounds[1] if r[:fitted])
println("distinct keys among the 2,023 fitted sampled base mechanisms: ", length(sample_keys),
        "; key duplicates inside the full base tier: ", length(base) - nskip - nerr - length(keys_full))
for (ri, recs) in enumerate(data.rounds)
    ri == 1 && continue
    fit = [r for r in recs if r[:fitted]]
    extra = [r for r in fit if !(:D in r[:labels]) && haskey(keys_full, r[:key]) &&
             !(r[:key] in sample_keys)]
    @printf("round %d: fitted %d; census D %d; not D but the key of an unsampled base mechanism %d (plain %d)\n",
            ri - 1, length(fit), count(r -> :D in r[:labels], fit), length(extra),
            count(r -> !get(r, :allo, false), extra))
    for r in extra[1:min(3, end)]
        println("   ", r[:id], ": ", first(r[:text], 300))
    end
end
