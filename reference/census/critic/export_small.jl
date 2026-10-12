# ABOUTME: Exports the small census's per-candidate classes, sizes and keys to one TSV, using
# ABOUTME: the class rules of census/small/tables.jl, so the critic can weight them by fit cost.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small/lib.jl")
const OUT = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/critic/small.tsv"

function wlist(r)
    w = String[]
    get(r, :v0, false) && push!(w, "v=0")
    r[:re_turnover] === nothing || push!(w, "RE turnover")
    for (i, f) in enumerate(r[:faces])
        r[:face_kind][i] === :none || push!(w, string(r[:face_kind][i]))
        r[:dvan][i] && push!(w, "0/0")
    end
    r[:admit] === nothing || push!(w, "onlyA")
    unique!(w)
end
admit_kind(msg) = occursin("has none", msg) ? "nofree" :
                  occursin("flips at one free form", msg) ? "secondfree" :
                  occursin("still join the two sides", msg) ? "G1in" :
                  occursin("close a cycle", msg) ? "G1cycle" :
                  occursin("reach the free enzyme only", msg) ? "island" : "other"
function psub(r)
    get(r, :P, false) || return ""
    dp = get(r, :deficit_parts, nothing)
    dp === nothing && return "pkg_only"
    haskey(dp, :L_phantom_n1) && return "L_phantom"
    join(sort!(string.(collect(keys(dp)))), "+")
end
kh(r) = haskey(r, :key) ? string(hash(r[:key]); base = 16) : ""

const NAMES = ["uni_uni_rxn", "uni_uni_with_reg", "uni_uni_allo", "uni_uni_allo_reg",
               "uni_uni_allo_2reg", "uni_uni_reg_and_inhibitor", "uni_bi_rxn", "bi_bi_rxn",
               "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
open(OUT, "w") do io
    println(io, join(["run", "round", "gid", "parent", "move", "fitted", "allo", "allo_n",
                      "np", "np_slice", "dterms", "rank_pkg", "P", "Psub", "W", "Wkind",
                      "D", "Dof", "Z", "N", "key", "eqh"], '\t'))
    for name in NAMES
        d = deserialize(joinpath(CEN, "data", "$name.jls"))
        recs = d.recs
        byg = Dict(r[:gid] => r for r in recs)
        first_of = Dict{String, Int}()
        nfit = 0
        for r in recs
            r[:status] === :ok || continue
            fitted = get(r, :fitted, false)
            if fitted
                nfit += 1
                k = r[:key]
                haskey(first_of, k) ? (r[:D] = first_of[k]) : (first_of[k] = r[:gid])
                r[:W] = wlist(r)
                r[:N] = !r[:vmax_fwd] || !r[:vmax_rev]
                p = r[:parent]
                r[:Z] = p != 0 && get(byg[p], :key, nothing) == r[:key]
                r[:P] = r[:rank_pkg] < r[:np]
            end
            wk = !fitted ? "" : join([w == "onlyA" ? "onlyA:" * admit_kind(r[:admit]) : w
                                      for w in r[:W]], ",")
            println(io, join(Any[name, r[:round], r[:gid], r[:parent], r[:move], fitted,
                                 r[:allo], get(r, :allo_n, 1), r[:np],
                                 get(r, :np_slice, -1), get(r, :dterms, -1),
                                 get(r, :rank_pkg, -1), get(r, :P, false), psub(r),
                                 fitted && !isempty(r[:W]), wk, haskey(r, :D),
                                 get(r, :D, 0), get(r, :Z, false), get(r, :N, false), kh(r),
                                 string(get(r, :eqh, UInt64(0)); base = 16)], '\t'))
        end
        @printf("%s: %d records, %d fitted\n", name, length(recs), nfit)
        flush(stdout)
    end
end
println("maxrss MB ", round(Int, Sys.maxrss() / 2^20))
