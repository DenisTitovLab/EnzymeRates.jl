# ABOUTME: Summarizes the P5a case records: per uni-uni tag combination and binding kind, the
# ABOUTME: guard's verdict against cone feasibility, uniqueness and agreement with deletion.
using Serialization, Printf
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
const WHICH = isempty(ARGS) ? "cases" : ARGS[1]
recs = deserialize(joinpath(@__DIR__, "data", "results_$(WHICH).jls"))
byid = Dict{Tuple{String, Symbol}, Dict{Symbol, Any}}((r[:id], r[:conv]) => r for r in recs)
short(t) = t === :OnlyA ? "O" : t === :EqualAI ? "E" : "N"

"""One-word verdict of a cone record: empty, unique (=del or ≠del) or not unique."""
function verdict(r)
    haskey(r, :error) && return "ERR"
    r[:cone_dim] == 0 && !haskey(r, :cone) && return "empty"
    !haskey(r, :cone) && return "empty"
    u = r[:cone_spread] < 1e-20
    seq = r[:cone_seq] ? "+seq" : ""
    u ? (r[:cone_vs_del] < 1e-20 ? "=del" : @sprintf("≠del(%.0e)", r[:cone_vs_del])) * seq :
        @sprintf("multi(%.0e)", r[:cone_spread]) * seq
end

if WHICH == "cases"
    println("uni-uni grid: tags (S binding, chemistry, P binding); O = :OnlyA, E = :EqualAI, N = :NonequalAI")
    println("columns per binding kind S/P: guard | conv (i) cone | conv (ii) cone | RREF-path order spread (i)")
    combos = [(a, b, c) for a in (:OnlyA, :EqualAI, :NonequalAI) for b in (:OnlyA, :EqualAI, :NonequalAI)
              for c in (:OnlyA, :EqualAI, :NonequalAI)]
    for kind in ("RE/RE", "RE/SS", "SS/RE", "SS/SS")
        println("── $kind")
        for (tS, tC, tP) in combos
            id = "uni/$kind/$tS,$tC,$tP"
            ri, rii = byid[(id, :i)], byid[(id, :ii)]
            (tS !== :OnlyA && tC !== :OnlyA && tP !== :OnlyA) && continue
            @printf("   %s%s%s  %-4s | (i) %-22s | (ii) %-22s | %.1e%s\n", short(tS), short(tC),
                    short(tP), ri[:admitted] ? "adm" : "REJ", verdict(ri), verdict(rii),
                    ri[:order_spread], ri[:order_spread] > 1e-12 ? "  ← order-dependent" : "")
        end
    end
    # agreement of the guard with cone feasibility over the whole grid
    for conv in (:i, :ii)
        agree = 0; tot = 0; dis = String[]
        for r in recs
            startswith(r[:id], "uni/") && r[:conv] === conv || continue
            tot += 1
            feas = haskey(r, :cone)
            feas == r[:admitted] ? (agree += 1) : push!(dis, r[:id] * (feas ? " (cone feasible, guard rejects)" : " (guard admits, cone empty)"))
        end
        @printf("guard vs cone (%s) feasibility on the uni-uni grid: %d of %d agree\n", conv, agree, tot)
        foreach(d -> println("   ", d), dis)
    end
end

println("\nnon-grid cases:")
for r in recs
    startswith(r[:id], "uni/") && continue
    @printf("   %-20s %-3s %-4s n=%d cone(d=%d,n=%d): %-24s RREF in cone=%s order spread %.1e\n",
            r[:id], r[:conv], r[:admitted] ? "adm" : "REJ", r[:n], r[:cone_dim], r[:cone_n],
            verdict(r), r[:rref_in_cone], r[:order_spread])
end
# totals
for conv in (:i, :ii)
    rs = filter(r -> r[:conv] === conv && !haskey(r, :error) && r[:admitted], recs)
    withO = filter(r -> haskey(r, :cone) && r[:cone_dim] > 0, rs)
    nu = count(r -> r[:cone_spread] < 1e-20, withO)
    nd = count(r -> r[:cone_spread] < 1e-20 && r[:cone_vs_del] < 1e-20, withO)
    ns = count(r -> r[:cone_seq], withO)
    ne = count(r -> !haskey(r, :cone), rs)
    no = count(r -> r[:order_spread] > 1e-12, rs)
    nin = count(r -> r[:rref_in_cone], withO)
    @printf("conv (%s), admitted with a moved column: %d; cone empty %d; unique %d; unique and = deletion %d; sequestration possible %d; RREF path order-dependent %d of %d admitted; RREF(pkg) path inside the cone %d\n",
            conv, length(withO), ne, nu, nd, ns, no, length(rs), nin)
end
