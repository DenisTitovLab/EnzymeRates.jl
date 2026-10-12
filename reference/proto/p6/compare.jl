# ABOUTME: Cross-checks P6's cached analyses against the package's (RE segments, exponent extras,
# ABOUTME: cycle space with Haldane multiples) and pairs the two timing runs mechanism by mechanism.
include(joinpath(@__DIR__, "harness.jl"))
using Printf, Serialization, Statistics

pkg = deserialize(joinpath(@__DIR__, "data", "pkg_measure.jls"))
gmt = deserialize(joinpath(@__DIR__, "data", "gm_timing.jls"))
ms = _testhelper_load_all()
@assert getfield.(pkg, :id) == getfield.(ms, :id) == getfield.(gmt, :id)

"""Rank of an integer matrix by exact elimination over the rationals."""
function _testhelper_rank(A::Matrix{Int})
    R = Rational{BigInt}.(A)
    rank, (nr, nc) = 0, size(R)
    for c in 1:nc
        p = findfirst(r -> R[r, c] != 0, rank+1:nr)
        p === nothing && continue
        p += rank
        R[[rank + 1, p], :] = R[[p, rank + 1], :]
        for r in 1:nr
            r != rank + 1 && R[r, c] != 0 && (R[r, :] .-= (R[r, c] / R[rank + 1, c]) .* R[rank + 1, :])
        end
        rank += 1
    end
    rank
end

bad = Dict(k => String[] for k in (:segments, :extras, :cycle_space, :haldane))
for (i, m) in enumerate(ms)
    g, edge_of, flipped = _testhelper_build(_testhelper_toraw(m))
    a = GM.analysis(g)
    names = [GM.form_name(g, f) for f in eachindex(g.forms)]
    mysegs = [sort([names[f] for f in eachindex(g.forms) if a.segment[f] == s]) for s in 1:maximum(a.segment)]
    Set(mysegs) == Set(pkg[i].segments) || push!(bad[:segments], m.id)
    myextras = Dict(names[f] => Dict(String(a.names[r]) => a.extras[r, f]
                                     for r in eachindex(a.names) if a.extras[r, f] > 0)
                    for f in eachindex(g.forms))
    myextras == pkg[i].extras || push!(bad[:extras], m.id)
    # My cycles over the package's flat steps in its stored direction.
    stepof = invperm(edge_of)
    M = zeros(Int, length(a.cycles), length(m.steps))
    for (r, c) in enumerate(a.cycles), sk in c
        j = stepof[abs(sk)]
        M[r, j] += sign(sk) * (flipped[j] ? -1 : 1)
    end
    C = pkg[i].cycles
    rC, rM = _testhelper_rank(C), _testhelper_rank(M)
    rC == rM == size(C, 1) == size(M, 1) == _testhelper_rank([C; M]) || push!(bad[:cycle_space], m.id)
    CH, MH = [C pkg[i].haldane], [M a.haldane]
    _testhelper_rank([CH; MH]) == rC || push!(bad[:haldane], m.id)
end
for (k, v) in sort!(collect(bad); by = first)
    @printf("%-12s mismatches %d of %d  %s\n", k, length(v), length(ms), join(first(v, 6), ", "))
end

q(v) = (v = filter(isfinite, Float64.(v)); isempty(v) ? (NaN, NaN, NaN) : (median(v), minimum(v), maximum(v)))
println("\nper-mechanism ratio package / GM — median, min, max")
for (label, sel) in (("plain", i -> !ms[i].allosteric), ("allosteric", i -> ms[i].allosteric))
    ix = filter(sel, eachindex(ms))
    println("  $label ($(length(ix)))")
    for (name, p, gfun) in (
            ("constructor (_with vs canonicalize)", r -> r.t_ctor, r -> r.t_ctor),
            ("canonical form only (_canonical_groups)", r -> r.t_canon, r -> r.t_ctor),
            ("allocated bytes", r -> r.alloc_bytes, r -> r.alloc_bytes),
            ("allocation count", r -> r.alloc_count, r -> r.alloc_count),
            ("hash", r -> r.t_hash, r -> r.t_hash),
            ("==", r -> r.t_eq, r -> r.t_eq),
            ("RE segments + cycle basis vs analyze", r -> r.t_segments + r.t_cycles, r -> r.t_analyze),
            ("bottomless check vs analyze + bottomless", r -> r.t_bottomless, r -> r.t_analyze + r.t_bottomless),
            ("summarysize", r -> r.size_fresh, r -> r.size_fresh),
            ("serialized bytes", r -> r.ser_fresh, r -> r.ser_fresh),
            ("inactive: _state_parts vs mask", r -> r.t_state_parts, r -> r.t_mask),
            ("inactive + its analyses", r -> r.t_state_parts + r.t_inactive_cycles, r -> r.t_mask + r.t_mask_analyze))
        @printf("    %-42s %8.1f %8.1f %8.1f\n", name, q([p(pkg[i]) / gfun(gmt[i]) for i in ix])...)
    end
end
