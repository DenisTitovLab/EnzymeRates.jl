# ABOUTME: Runs the :OnlyA limit study (P5a) on the hand-written cases and selected population
# ABOUTME: members for conventions (i) and (ii); prints the case table and serializes every record.
include(joinpath(@__DIR__, "limits_lib.jl"))

const WHICH = isempty(ARGS) ? "cases" : ARGS[1]

function members()
    if WHICH == "cases"
        ms = MWC.Pop.PMech[retag(m) for m in deserialize(joinpath(@__DIR__, "data", "cases.jls"))]
        hw = Pop.load_population("handwritten")
        push!(ms, only(filter(m -> m.id == "handwritten/34", hw)))
        sp = filter(m -> m.allosteric && :OnlyA in m.group_tags, Pop.load_population("specs"))
        append!(ms, sp)
        return ms
    elseif WHICH == "allo_uniuni"
        return filter(m -> :OnlyA in m.group_tags, Pop.load_population("allo_uniuni"))
    elseif WHICH == "allo_bibi_KV"
        pop = Pop.load_population("allo_bibi")
        ch = Dict(m.id => MWC.is_chem(m) for m in pop)
        hasbind(m) = any(j -> MWC.step_tags(m)[j] === :OnlyA && !ch[m.id][j], eachindex(m.steps))
        V = filter(m -> !hasbind(m), pop)[1:20]
        Kpool = filter(hasbind, pop)
        K = Kpool[round.(Int, range(1, length(Kpool); length = 20))]
        return [K; V]
    end
    error("unknown set $WHICH")
end

ms = members()
recs = Dict{Symbol, Any}[]
fmt(x) = x isa AbstractFloat ? @sprintf("%.1e", x) : string(x)
for m in ms, conv in (:i, :ii)
    t = @elapsed r = try
        study(m, conv; seed = 7)
    catch e
        Dict{Symbol, Any}(:id => m.id, :conv => conv, :error => first(sprint(showerror, e), 200))
    end
    r[:time] = t; r[:source] = m.source; r[:tags] = m.group_tags
    push!(recs, r)
    if haskey(r, :error)
        @printf("%-36s %-3s ERROR %s\n", m.id, conv, r[:error]); continue
    end
    cone = haskey(r, :cone) ? @sprintf("cone d=%d n=%d spread=%.1e seq=%s Δdel=%.1e orc=%.0e v=%s",
            r[:cone_dim], r[:cone_n], r[:cone_spread], r[:cone_seq], r[:cone_vs_del],
            r[:cone_oracle_err], fmt(Float64(r[:cone_v]))) : @sprintf("cone d=%d EMPTY", r[:cone_dim])
    @printf("%-36s %-3s %s | rref: partner=%d wrong=%d inCone=%s cons=%.0e c14=%.0e sym=%.0e ordI(N,Q,D)=%s v=%s del=%s Δdel=%.1e orderSpread=%.1e | %s | %.1fs\n",
            m.id, conv, r[:admitted] ? "adm" : "REJ", length(r[:partner]), length(r[:wrong]),
            r[:rref_in_cone], r[:consistency], r[:conv_err][6], r[:sym_vs_oracle],
            join([o === nothing ? "-" : isinteger(o[1]) ? string(Int(o[1])) : string(o[1]) for o in r[:ordI]], ","),
            fmt(Float64(r[:v_lim_sym])), fmt(Float64(r[:v_del])), r[:lim_vs_del],
            r[:order_spread], cone, t)
    isempty(r[:partner]) || println("      partner moves: ", r[:partner])
    isempty(r[:wrong]) || println("      wrong-sign: ", r[:wrong])
end
serialize(joinpath(@__DIR__, "data", "results_$(WHICH).jls"), recs)
