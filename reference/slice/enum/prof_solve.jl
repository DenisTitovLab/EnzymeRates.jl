# ABOUTME: Times the pieces of thermo.jl's solve on split candidates, to see what the split
# ABOUTME: move's gain test spends its time on.
include(joinpath(@__DIR__, "load.jl"))
d = deserialize(joinpath(S, "slice/enum/data/pkg_dead1.jls"))
ms = [Slice.Converted(r.parent, 0).mech for r in d.records[1:60]]
best(f, n = 20) = (f(); minimum(@elapsed(f()) for _ in 1:n))
for m in ms[1:6]
    fwd = Slice.orient(m)
    mask = Slice.inactive_mask(m)
    cols, A, b = Slice.constraint_system(m, mask)
    @printf("edges %2d forms %2d cols %3d rows %3d | solve %7.1f µs  cycles %7.1f  columns %7.1f  system %7.1f  rref %7.1f  orient %5.2f\n",
            length(m.edges), length(m.forms), length(cols), size(A, 1),
            1e6 * best(() -> Slice.solve(m)), 1e6 * best(() -> Slice.cycles(m, fwd)),
            1e6 * best(() -> Slice.columns(m, fwd, mask)),
            1e6 * best(() -> Slice.constraint_system(m, mask)),
            1e6 * best(() -> Slice.rref!(Slice.Exp.([A b]), length(cols))),
            1e6 * best(() -> Slice.orient(m)))
end
t = [best(() -> Slice.split_children(m), 3) for m in ms]
@printf("split median %.1f µs over %d parents\n", 1e6 * median(t), length(ms))
