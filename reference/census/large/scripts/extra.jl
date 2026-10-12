# ABOUTME: Small follow-ups on saved runs: where round-0 rank-deficient candidates come from
# ABOUTME: (seed or flip child of a degenerate seed), and per-round W/D/P/Z overlaps by allostery.
include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))
for name in Symbol.(ARGS)
    data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
    info, rounds = data.info, data.rounds
    label!(rounds)
    r0 = info[:r0]
    nseed = r0.n_seeds - r0.n_degenerate          # base = these seeds, then the flip children
    fit0 = [r for r in rounds[1] if r[:fitted]]
    isflip(r) = r[:id][2] > nseed
    for l in (:P, :W_split, :W_trap)
        rs = [r for r in fit0 if l in r[:labels]]
        println("$name round 0 $l: $(length(rs)) fitted, $(count(isflip, rs)) of them flip children of degenerate seeds (fitted flip children: $(count(isflip, fit0)) of $(length(fit0)))")
    end
    for (ri, recs) in enumerate(rounds)
        fit = [r for r in recs if r[:fitted]]
        for allo in (false, true)
            fs = [r for r in fit if get(r, :allo, false) == allo]
            isempty(fs) && continue
            u = count(r -> any(l -> isW(l) || l in (:D, :P, :Ptoday, :Z), r[:labels]), fs)
            @printf("%s round %d %s: fitted %d, W %d, D %d, P %d, Ptoday %d, Z %d, any %d (%s)\n",
                    name, ri - 1, allo ? "allosteric" : "plain", length(fs),
                    count(r -> any(isW, r[:labels]), fs), count(r -> :D in r[:labels], fs),
                    count(r -> :P in r[:labels], fs), count(r -> :Ptoday in r[:labels], fs),
                    count(r -> :Z in r[:labels], fs), u, pct(u, length(fs)))
        end
    end
end
