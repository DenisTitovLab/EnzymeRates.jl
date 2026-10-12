# ABOUTME: Census of the non-monomial common factors the P8 driver flagged, per population and key.
# ABOUTME: Groups the flagged variable sets by shape (concentrations vs parameters) and lists examples.
using Serialization
for pop in ARGS
    d = deserialize(joinpath(@__DIR__, "data", "$pop.jls"))
    for variant in (:today, :redesign), kk in (:k1, :k2)
        shapes = Dict{Any, Vector{String}}()
        for r in d.recs
            cf = get(r, (variant, kk, :common_factor), nothing)
            (cf === nothing || isempty(cf)) && continue
            conc = sort([string(s) for s in cf if length(string(s)) == 1])
            par = sort([string(s) for s in cf if length(string(s)) > 1])
            push!(get!(shapes, (conc, length(par)), String[]), string(r[:id], " ", par))
        end
        tot = sum(length, values(shapes); init = 0)
        println(pop, " ", variant, " ", kk, ": ", tot, " laws")
        for (k, v) in sort!(collect(shapes); by = x -> -length(x[2]))
            println("   concs ", k[1], " + ", k[2], " params: ", length(v), "   e.g. ", first(v, 2))
        end
    end
end
