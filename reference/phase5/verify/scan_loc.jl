# ABOUTME: Locates the verifier workload's hoisted allocations: loads scan_v.jl's functions for one
# ABOUTME: variant (ARGS[1]), then compiles and runs each loop after a marker line on stderr.
ENV["LOC"] = "1"
include(joinpath(@__DIR__, "scan_v.jl"))
ms = convert(Vector{S.Mechanism}, deserialize(IOBuffer(BYTES)))
ls = LT[S.rate_law(m) for m in pop]
for (name, f) in (("derive_untyped", () -> derive(BYTES)), ("derive_typed", () -> derive_typed(ms)),
                  ("merge_vec", () -> merge_vec(ls)), ("merge_nothing", () -> merge_nothing(ls)),
                  ("merge_same", () -> merge_same(ls)), ("pass1", () -> pass1(pop)))
    println(stderr, "=== MARK ", name); flush(stderr)
    f()
    println(stderr, "=== END ", name); flush(stderr)
end
