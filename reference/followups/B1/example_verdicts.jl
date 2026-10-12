# ABOUTME: Follow-up B1: verdicts of the report's example mechanisms under every source (P5a 16-direction,
# ABOUTME: full-direction ii/iu/iiu/Dk, Dka) plus a uni-uni no_island example. Usage: julia example_verdicts.jl id...
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
for id in ARGS
    name = String(split(id, "/")[1])
    p5a = load_records(name); full = load_full(name)
    m = only(filter(x -> x.id == id, load_allo(name)))
    a = analyze(m, :D)
    vs = ["P5a(ii)=" * verdict(get(p5a, (id, :ii), nothing)), "P5a(iu)=" * verdict(get(p5a, (id, :iu), nothing))]
    for c in (:ii, :iu, :iiu, :Dk, :Dka)
        r = get(full, (id, c), nothing)
        r === nothing && continue
        push!(vs, "$c=" * verdict(r) * (haskey(r, :cone_vs_del) ? @sprintf("(Δdel %.1e)", r[:cone_vs_del]) : ""))
    end
    println(rpad(id, 18), " R_D=", rpad(a.cls, 10), " R_iu=", rpad(analyze(m, :iu).cls, 10), " ", join(vs, " "))
end
um = filter(m -> analyze(m, :D).cls == "no_island" && :OnlyA in m.group_tags, load_allo("allo_uniuni"))
print_mech(stdout, um[1], :D)
