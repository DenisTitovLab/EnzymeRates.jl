# ABOUTME: Times the P6 graph type on every population (constructor under both orientation rules,
# ABOUTME: hash, ==, analyses, inactive mask) and records sizes; writes data/gm_timing.jls.
include(joinpath(@__DIR__, "harness.jl"))
using Printf, Serialization, Statistics

ms = _testhelper_load_all()
serbytes(x) = (io = IOBuffer(); serialize(io, x); position(io))
rows = NamedTuple[]
for (i, m) in enumerate(ms)
    raw = _testhelper_toraw(m)
    build() = GM.canonicalize(raw[1], raw[2], raw[3], raw[4]; raw[5]...)
    g = first(build())
    g2 = _testhelper_gm(_testhelper_toraw(m; rng = Xoshiro(i), flip = trues(length(m.steps))))
    @assert g == g2
    t_ctor = _testhelper_time_ns(build, 20)
    alloc_bytes = @allocated build()
    alloc_count = @allocations build()
    t_phys_index = _testhelper_time_ns(() -> GM.physical(first(build()); tier3 = :index), 20)
    t_phys_name = _testhelper_time_ns(() -> GM.physical(first(build()); tier3 = :name), 20)
    t_orient_index = _testhelper_time_ns(() -> GM.physical_forward(g; tier3 = :index), 20)
    t_hash = _testhelper_time_ns(() -> hash(g), 200)
    t_eq = _testhelper_time_ns(() -> g == g2, 200)
    t_analyze = _testhelper_time_ns(() -> GM.analyze(g), 20)
    size_fresh = Base.summarysize(g)
    ser_fresh = serbytes(g)
    a = GM.analysis(g)
    t_cached = _testhelper_time_ns(() -> GM.analysis(g), 1000)
    t_bottomless = _testhelper_time_ns(() -> GM.bottomless(g, a), 100)
    size_cached = Base.summarysize(g)
    ser_cached = serbytes(g)
    t_mask = NaN; t_mask_analyze = NaN
    if m.allosteric
        t_mask = _testhelper_time_ns(() -> GM.inactive_mask(g), 50)
        mask = GM.inactive_mask(g)
        t_mask_analyze = _testhelper_time_ns(() -> GM.analyze(g, mask), 20)
    end
    push!(rows, (; id = m.id, allosteric = m.allosteric, nforms = length(g.forms),
                 nsteps = length(g.edges), t_ctor, alloc_bytes, alloc_count, t_phys_index,
                 t_phys_name, t_orient_index, t_hash, t_eq, t_analyze, t_cached, t_bottomless,
                 t_mask, t_mask_analyze, size_fresh, size_cached, ser_fresh, ser_cached,
                 ser_analysis = serbytes(a), pkg_ctor_ns = m.ctor_ns))
end
serialize(joinpath(@__DIR__, "data", "gm_timing.jls"), rows)

q(v) = (v = filter(!isnan, Float64.(v)); isempty(v) ? (NaN, NaN, NaN) : (median(v), quantile(v, 0.9), maximum(v)))
for (label, sel) in (("plain", r -> !r.allosteric), ("allosteric", r -> r.allosteric), ("all", r -> true))
    rs = filter(sel, rows)
    println("== $label ($(length(rs)) mechanisms) — median / p90 / max")
    for f in (:t_ctor, :alloc_bytes, :alloc_count, :t_phys_index, :t_phys_name, :t_orient_index,
              :t_hash, :t_eq, :t_analyze, :t_cached, :t_bottomless, :t_mask, :t_mask_analyze,
              :size_fresh, :size_cached, :ser_fresh, :ser_cached, :ser_analysis, :pkg_ctor_ns)
        @printf("  %-15s %10.1f %10.1f %10.1f\n", f, q(getfield.(rs, f))...)
    end
end
