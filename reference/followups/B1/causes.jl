# ABOUTME: Follow-up B1: attributes each disagreement between the R_D kept-ratio class and P5a's (ii) verdict to
# ABOUTME: sampling (full (ii) agrees), the group tie (untied iiu agrees), the ratio semantics (Dka agrees) or other.
include(joinpath(@__DIR__, "b1_lib.jl"))
function load_full(name)
    d = Dict{Tuple{String, Symbol}, Dict{Symbol, Any}}()
    dir = joinpath(@__DIR__, "data")
    for f in readdir(dir)
        if (startswith(f, "full_$(name).jls") || startswith(f, "full_$(name)_")) && !endswith(f, "_ids.jls")
            for r in deserialize(joinpath(dir, f)); d[(r[:id], r[:conv])] = r; end
        elseif startswith(f, "dka_")
            for r in deserialize(joinpath(dir, f)); startswith(r[:id], name * "/") && (d[(r[:id], :Dka)] = r); end
        end
    end
    d
end
agrees(c, v) = c in ("no_island", "emptied") ? v == "U=D" : v == "U=D" ? false :
               c == "phantom" ? v == "U≠D" : c in ("undefined", "ambiguous") ? startswith(v, "NU") || v == "empty" :
               c == "sequestering" ? endswith(v, "+seq") : false
cause = Dict{String, Vector{String}}()
for name in POPS
    p5a = load_records(name); full = load_full(name)
    for m in load_allo(name)
        :OnlyA in m.group_tags || continue
        c = analyze(m, :D).cls
        v = verdict(get(p5a, (m.id, :ii), nothing))
        (v == "—" || agrees(c, v)) && continue
        fv(k) = verdict(get(full, (m.id, k), k === :Dka ? get(full, (m.id, :Dk), nothing) : nothing))
        k = agrees(c, fv(:ii)) ? "sampling (16 dirs)" : agrees(c, fv(:iiu)) ? "group tie" :
            agrees(c, fv(:Dka)) ? "ratio free vs kept" : "other"
        push!(get!(cause, "$k | $c", String[]), "$(m.id)[$v; ii=$(fv(:ii)) iiu=$(fv(:iiu)) Dka=$(fv(:Dka))]")
    end
end
for (k, ids) in sort!(collect(cause); by = first)
    println(rpad(k, 40), length(ids))
    println("    ", join(ids, " "))
end
