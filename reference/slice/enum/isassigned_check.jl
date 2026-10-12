# ABOUTME: Checks whether isassigned with an Int32 index takes Base's exception-based fallback,
# ABOUTME: the suspected cause of thermo.jl's slow cycles.
v = Vector{Vector{Int}}(undef, 8); v[1] = Int[]
f(v, i) = isassigned(v, i)
f(v, Int32(3)); f(v, 3)
t32 = minimum(@elapsed(for _ in 1:1000; f(v, Int32(3)); end) for _ in 1:5) / 1000
t64 = minimum(@elapsed(for _ in 1:1000; f(v, 3); end) for _ in 1:5) / 1000
println("isassigned unassigned slot: Int32 index ", round(1e9 * t32; digits = 1), " ns, Int index ",
        round(1e9 * t64; digits = 1), " ns")
println(which(isassigned, Tuple{Vector{Vector{Int}}, Int32}))
