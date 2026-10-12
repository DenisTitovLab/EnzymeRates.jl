# ABOUTME: Follow-up B1: how often P5a's 16-direction verdict changes when every enumerated cone direction
# ABOUTME: is evaluated (conventions ii and iu), per population; lists the U=D → not-U=D flips.
include(joinpath(@__DIR__, "b1_lib.jl"))
function load_full(name)
    d = Dict{Tuple{String, Symbol}, Dict{Symbol, Any}}()
    dir = joinpath(@__DIR__, "data")
    for f in readdir(dir)
        (startswith(f, "full_$(name).jls") || startswith(f, "full_$(name)_")) && !endswith(f, "_ids.jls") || continue
        for r in deserialize(joinpath(dir, f)); d[(r[:id], r[:conv])] = r; end
    end
    d
end
for conv in (:ii, :iu)
    tot = 0; changed = 0; flips = String[]; adm = 0
    for name in POPS
        p5a = load_records(name); full = load_full(name)
        for ((id, c), r) in full
            c === conv || continue
            p = get(p5a, (id, conv), nothing); p === nothing && continue
            tot += 1
            vp, vf = verdict(p), verdict(r)
            vp == vf && continue
            changed += 1
            ok(v) = v in ("U=D", "U≠D")
            ok(vp) != ok(vf) && push!(flips, "$id:$vp→$vf")
        end
    end
    println("($conv): $tot mechanisms with both verdicts; verdict changed in $changed; admit/reject (unique, no seq) flipped in $(length(flips)):")
    println("   ", join(flips, " "))
end
