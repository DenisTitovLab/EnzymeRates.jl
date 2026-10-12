# ABOUTME: Verifies class P with an independent rank (complex-step log-Jacobian of the compiled
# ABOUTME: rate_equation at new random points) and tests the census's same-law reading on package laws.
include(joinpath(@__DIR__, "vlib.jl"))
const SAMPLES = deserialize(joinpath(VER, "samples.jls"))

primary(r) = (dp = get(r, :deficit_parts, Dict{Symbol, Int}());
              isempty(dp) ? :package_only :
              first(k for k in (:L_phantom_n1, :ss_noflux, :chain_flank, :nonequal_copy, :re_iso, :other)
                    if haskey(dp, k)))

rng = Xoshiro(424242)
agree = Dict{String, Vector{Int}}()
bump!(k, ok) = (v = get!(agree, k, [0, 0]); v[1] += 1; v[2] += ok)
samelaw = Dict{Tuple{Symbol, String}, Int}()

for name in ["uni_uni_rxn", "uni_uni_with_reg", "uni_uni_allo", "uni_uni_allo_reg",
             "uni_uni_allo_2reg", "uni_uni_reg_and_inhibitor", "uni_bi_rxn", "bi_bi_rxn",
             "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = loadrun(name)
    s = SAMPLES[name]
    println("\n== $name")
    for (lab, gs) in (("P", s[:P]), ("unflagged", s[:P_unflagged]))
        for g in gs
            r = d.recs[g]
            pl = plaw(mechs[g])
            a = cs_rank(pl, rng)
            b = cs_rank(pl, rng; lo = -1.0, hi = 3.0)
            amb = a.kept < 1e-6 || a.dropped > 1e-11 || b.kept < 1e-6 || b.dropped > 1e-11
            ok = a.rank == b.rank == r[:rank_pkg] && length(pl.names) == r[:np]
            bump!("$lab rank agrees with census", ok)
            amb && bump!("$lab ambiguous gap", true)
            @printf("  [%s] #%d np %d/%d rank %d,%d census %d gaps kept %.1e,%.1e dropped %.1e,%.1e %s%s\n",
                    lab, g, length(pl.names), r[:np], a.rank, b.rank, r[:rank_pkg], a.kept, b.kept,
                    a.dropped, b.dropped, ok ? "" : "DISAGREE ", amb ? "AMBIGUOUS" : "")
            # same law as the parent (package laws), for the census's P reading
            if lab == "P" && r[:parent] != 0
                p = d.recs[r[:parent]]
                pp = plaw(mechs[r[:parent]])
                ab, ba = same_family(pl, pp, rng)
                verdict = ab < 1e-8 && ba < 1e-8 ? "same law as parent" :
                          ba < 1e-8 ? "contains parent's law" :
                          ab < 1e-8 ? "inside parent's law" : "different"
                k = (primary(r), verdict)
                samelaw[k] = get(samelaw, k, 0) + 1
                @printf("      %s vs parent #%d (np %d rank %d): child-in-parent %.1e, parent-in-child %.1e → %s\n",
                        primary(r), r[:parent], p[:np], p[:rank_pkg], ab, ba, verdict)
            end
        end
    end
    d = mechs = nothing
    GC.gc()
end
println("\nRank agreement (n, agree):")
for (k, v) in sort!(collect(agree))
    println("  $k: $(v[2]) of $(v[1])")
end
println("\nP candidate vs its parent, package laws (primary sub-class, verdict):")
for (k, v) in sort!(collect(samelaw))
    println("  $(k[1]) / $(k[2]): $v")
end
