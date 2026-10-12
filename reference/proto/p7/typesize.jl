# ABOUTME: Memory a type-level law costs: creating RateLaw types for the measurement set under three
# ABOUTME: encodings (NTuple{N,Int64}, NTuple{N,Int16}, a Symbol) and the type's own summarysize.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Select.jl"))
using Printf, Base64
struct LawT{C} end
mem() = (GC.gc(true); GC.gc(true); (pool = Int(ccall(:jl_gc_pool_live_bytes, Int64, ())), live = Int(Base.gc_live_bytes())))
set = selected()
for variant in ("flat", "facth")
    codes = [code(variant == "flat" ? P7.flat_law(m) : P7.factored_law(m), (variant != "flat", variant == "facth", 0)) for (_, m) in set]
    n = sum(length, codes)
    @printf("%s: %d laws, %d Ints in total (median %d, max %d)\n", variant, length(codes), n,
            sort(length.(codes))[end ÷ 2], maximum(length, codes))
    for (enc, f) in (("NTuple{N,Int64}", c -> Tuple(c)), ("NTuple{N,Int16}", c -> Tuple(Int16.(c))),
                     ("Symbol", c -> Symbol(base64encode(reinterpret(UInt8, Int16.(c))))))
        m0 = mem()
        ts = [LawT{f(c)} for c in codes]
        m1 = mem()
        sz = sum(Base.summarysize(t) for t in ts)
        @printf("   %-16s pool +%8.1f KB  live +%8.1f KB  summarysize(types) %8.1f KB  per law %6.1f KB\n",
                enc, (m1.pool - m0.pool) / 1024, (m1.live - m0.live) / 1024, sz / 1024, (m1.live - m0.live) / 1024 / length(ts))
    end
end
