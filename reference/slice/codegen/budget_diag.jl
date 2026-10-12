# ABOUTME: Diagnostic for check (c): lists the method instances one fresh mechanism's
# ABOUTME: identify-like flow adds, grouped by method, to attribute them to codegen or harness.
include(joinpath(@__DIR__, "budget.jl"))
function mis()
    out = Set{Any}()
    Base.visit(Core.methodtable) do m
        for mi in Base.specializations(m)
            push!(out, mi)
        end
    end
    out
end
for (set, i) in (("init_bibi", 30), ("init_bibi", 31), ("expand_bibi", 600), ("allo_bibi", 700))
    m = load(set)[i]
    before = mis()
    measure(m, rng)
    new = setdiff(mis(), before)
    byname = Dict{String, Int}()
    for mi in new
        k = string(mi.def isa Method ? mi.def.module : "", ".", mi.def isa Method ? mi.def.name : mi.def)
        byname[k] = get(byname, k, 0) + 1
    end
    println("== $set/$i: $(length(new)) new instances")
    for (k, n) in sort!(collect(byname); by = x -> -x[2])
        println("   $n  $k")
    end
end
