# ABOUTME: Recounts from the saved records: P sub-kinds per round, package-vs-slice count and rank
# ABOUTME: gaps, L-phantom shares, and the class union with and without the package-only P.
include(joinpath(@__DIR__, "vlib.jl"))

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
pct(a, b) = @sprintf("%.1f%%", 100a / b)

for name in RX
    data = load_run(name)
    println("="^100, "\n$name")
    for (ri, recs) in enumerate(data.rounds)
        fit = [r for r in recs if r[:fitted]]
        n = length(fit)
        P = [r for r in fit if :P in r[:labels]]
        kinds = Dict{String, Int}()
        for r in P, k in keys(get(r, :explain, Dict()))
            kinds[string(k)] = get(kinds, string(k), 0) + 1
        end
        Lph = count(r -> haskey(get(r, :explain, Dict()), :L_phantom), P)
        reonly = count(r -> haskey(get(r, :explain, Dict()), :re_iso_gauge) &&
                            !haskey(r[:explain], :L_phantom), P)
        Ponly = count(r -> :Ptoday in r[:labels] && !(:P in r[:labels]), fit)
        exc = n_params_gt = count(r -> r[:n_params] > r[:n_slice], fit)
        n_params_lt = count(r -> r[:n_params] < r[:n_slice], fit)
        X = [r for r in fit if :Xlaw in r[:labels]]
        gap = Dict{Int, Int}()
        foreach(r -> (d = r[:pkg_rank] - r[:rank]; gap[d] = get(gap, d, 0) + 1), X)
        isW(l) = startswith(string(l), "W")
        union_all = count(r -> any(l -> isW(l) || l in (:D, :P, :Ptoday, :N, :Z), r[:labels]), fit)
        union_slice = count(r -> any(l -> isW(l) || l in (:D, :P, :N, :Z), r[:labels]), fit)
        @printf("round %d: fitted %d; P %d (%s) by kind %s; L-phantom records %d (%s); RE-iso gauge without L phantom %d; package-only P %d\n",
                ri - 1, n, length(P), pct(length(P), n), kinds, Lph, pct(Lph, n), reonly, Ponly)
        @printf("   package count > slice count in %d fitted (of which P %d), < in %d; Xlaw %d, package rank - slice rank among them %s\n",
                n_params_gt, count(r -> r[:n_params] > r[:n_slice] && :P in r[:labels], fit),
                n_params_lt, length(X), sort!(collect(gap)))
        @printf("   W|D|P|Ptoday|N|Z %d (%s); without the package-only P (W|D|P|N|Z) %d (%s)\n",
                union_all, pct(union_all, n), union_slice, pct(union_slice, n))
    end
end
