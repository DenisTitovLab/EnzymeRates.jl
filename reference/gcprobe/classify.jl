# ABOUTME: Run at the end of a probe under the JuliaLICM detector: reads this process's own log,
# ABOUTME: resolves each hoisted gc_alloc_obj's type tag to a type, and says whether it holds pointers.
function classify_hoists(log)
    flush(stdout); flush(stderr)
    tags = Dict{UInt, Int}()
    for l in eachline(log)
        startswith(l, "+") && occursin("gc_alloc_obj", l) || continue
        m = match(r"gc_alloc_obj\(.*inttoptr \(i64 (\d+) to ptr\)", l)
        m === nothing && continue
        t = parse(UInt, m.captures[1]); tags[t] = get(tags, t, 0) + 1
    end
    println("HOISTED ALLOCATION TYPES (count, type, pointer-free):")
    for (t, n) in sort(collect(tags); by = x -> -x[2])
        T = unsafe_pointer_to_objref(Ptr{Cvoid}(t))
        pf = T isa DataType ? Base.datatype_pointerfree(T) : missing
        println("HOIST ", n, "  ", T, "  pointerfree=", pf)
    end
end
classify_hoists(ENV["CLASSIFY_LOG"])
