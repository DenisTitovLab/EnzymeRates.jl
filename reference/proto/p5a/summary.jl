# ABOUTME: Final P5a numbers over every record set: path consistency, oracle agreement, cone
# ABOUTME: uniqueness, deletion agreement, sequestration, naive 0/0 and leading-face cancellations.
using Serialization, Printf
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
D = joinpath(@__DIR__, "data")
sets = Dict(n => deserialize(joinpath(D, "results_$n.jls")) for n in ("cases", "allo_uniuni", "allo_bibi_KV"))
ord(o) = o === nothing ? nothing : o[1]
for (name, recs) in sort!(collect(sets); by = first)
    println("== $name")
    for conv in (:i, :ii)
        rs = filter(r -> r[:conv] === conv && !haskey(r, :error), recs)
        adm = filter(r -> r[:admitted], rs)
        moved = filter(r -> haskey(r, :cone) ? r[:cone_dim] > 0 || length(r[:cone][1].a) > 0 : true, adm)
        @printf("  conv (%s): %d records (%d admitted); errors %d\n", conv, length(rs), length(adm),
                count(r -> r[:conv] === conv && haskey(r, :error), recs))
        @printf("    consistency along the RREF path: max |v(Γ=Keq)|/|v(Γ=Keq/10)| = %.1e\n",
                maximum(r[:consistency] for r in rs))
        @printf("    symbolic leading face vs oracle at ε=1e-30 (RREF path): max %.1e;  cone directions: max %.1e\n",
                maximum(r[:sym_vs_oracle] for r in rs),
                maximum((r[:cone_oracle_err] for r in rs if haskey(r, :cone_oracle_err)); init = 0.0))
        @printf("    oracle convergence at ε=1e-14 vs 1e-30: max %.1e\n", maximum(r[:conv_err][6] for r in rs))
        nz = count(r -> (o = ord(r[:ordI][3]); o !== nothing && o > 0), adm)
        nq = count(r -> (o = ord(r[:ordI][2]); o !== nothing && o > 0), adm)
        nseq = count(r -> (q = ord(r[:ordI][2]); d = ord(r[:ordI][3]); q !== nothing && d !== nothing && q < d), adm)
        @printf("    RREF path, admitted: D_I vanishes (naive substitution gives 0/0) %d; Q_I vanishes %d; Q_I of lower order than D_I (sequestration) %d; leading-face cancellation %d\n",
                nz, nq, nseq, count(r -> r[:cancel], adm))
        cz = count(r -> haskey(r, :cone) && any(c -> (o = ord(c.ordI[3]); o !== nothing && o > 0), r[:cone]), adm)
        @printf("    cone directions with D_I vanishing (naive 0/0): %d mechanisms\n", cz)
    end
end
for name in ("expand_bibi", "allo_bibi")
    f = joinpath(D, "pop_$name.jls")
    isfile(f) || continue
    recs = deserialize(f)
    println("== population $name (light study, no ε sweep)")
    for conv in (:i, :ii, :iu)
        rs = filter(r -> r[:conv] === conv, recs)
        ok = filter(r -> !haskey(r, :error), rs)
        feas = filter(r -> haskey(r, :cone_v), ok)
        uniq = filter(r -> r[:cone_spread] < 1e-20, feas)
        uq_del = filter(r -> r[:cone_vs_del] < 1e-20, uniq)
        nonu = filter(r -> r[:cone_spread] >= 1e-20, feas)
        @printf("  conv (%s): %d mechanisms, errors %d; cone empty %d; unique %d (= deletion %d, ≠ deletion %d); non-unique %d; sequestration possible %d; naive 0/0 %d; RREF order-dependent %d; RREF partner moves %d; oracle err max %.1e\n",
                conv, length(rs), length(rs) - length(ok), length(ok) - length(feas), length(uniq),
                length(uq_del), length(uniq) - length(uq_del), length(nonu),
                count(r -> r[:cone_seq], feas), count(r -> r[:zero_over_zero], feas),
                count(r -> r[:order_spread] > 1e-12, ok), count(r -> r[:npartner] > 0, ok),
                maximum((r[:oracle_err] for r in feas); init = 0.0))
        bad = [r for r in feas if !(r[:cone_spread] < 1e-20 && r[:cone_vs_del] < 1e-20)]
        isempty(bad) || println("     not (unique and = deletion): ", join([r[:id] for r in bad[1:min(end, 25)]], " "),
                                length(bad) > 25 ? " …" : "")
        empt = [r[:id] for r in ok if !haskey(r, :cone_v)]
        isempty(empt) || println("     cone empty: ", join(empt[1:min(end, 25)], " "))
        # step-kind profile of the mechanisms that are not unique-and-deletion
        prof = Dict{String, Int}()
        for r in bad
            t = r[:tags]
            k = join(sort(unique(string.(t))), "+")
            prof[k] = get(prof, k, 0) + 1
        end
        isempty(prof) || println("     tag profile of those: ", prof)
    end
end
