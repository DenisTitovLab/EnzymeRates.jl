# ABOUTME: Compact summary table of the saved P8 runs: class counts per key, merges, key cost,
# ABOUTME: search and common-factor statistics, per population and mechanism kind.
using Serialization, Statistics, Printf
med(v) = isempty(v) ? NaN : median(v)
for tag in ARGS
    d = deserialize(joinpath(@__DIR__, "data", "$tag.jls"))
    for allo in (false, true)
        rs = [r for r in d.recs if r[:allosteric] == allo]
        isempty(rs) && continue
        @printf("%-26s %-5s n=%-5d eq_hash=%-5d\n", tag, allo ? "MWC" : "plain", length(rs),
                length(unique(r[:eq_hash] for r in rs)))
        for v in (:today, :redesign)
            ok = [r for r in rs if haskey(r, (v, :k1, :key))]
            cl(kk, f = :key) = length(unique(r[(v, kk, f)] for r in ok))
            kns = [r[(v, :k1c, :key_ns)] / 1e3 for r in ok]
            kns2 = [r[(v, :k2c, :key_ns)] / 1e3 for r in ok]
            canc = [r[(v, :k1c, :cancel_s)] * 1e6 for r in ok]
            fac = [r for r in ok if !isempty(r[(v, :k1c, :removed)])]
            t_before = sum(r[(v, :k1, :nterms)] for r in fac; init = 0)
            t_after = sum(r[(v, :k1c, :nterms)] for r in fac; init = 0)
            search = count(r -> any(r[(v, kk, :leaves)] > 1 for kk in (:k1, :k2, :k1c, :k2c)), ok)
            @printf("   %-8s exact=%-5d K1=%-5d K2=%-5d K1c=%-5d K2c=%-5d | key µs (K1c) med %.0f p99 %.0f max %.0f; (K2c) med %.0f max %.0f | cancel µs med %.0f p99 %.0f | factor laws %d terms %d→%d | search %d\n",
                    v, cl(:k1, :exact), cl(:k1), cl(:k2), cl(:k1c), cl(:k2c), med(kns),
                    quantile(kns, 0.99), maximum(kns), med(kns2), maximum(kns2), med(canc),
                    quantile(canc, 0.99), length(fac), t_before, t_after, search)
        end
    end
end
