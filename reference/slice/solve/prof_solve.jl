# ABOUTME: Profiles the pieces of the solve on a few converted mechanisms in a lean process
# ABOUTME: (only the converted mechanisms loaded), to find where its time goes.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "slice", "src", "Slice.jl"))
using Serialization, Printf, Statistics
const S = Slice
const DATA = joinpath(SCRATCH, "slice", "core", "data")
function cost(f)
    f()
    best = Inf
    for _ in 1:3
        n, t0 = 0, time_ns()
        while time_ns() - t0 < 2_000_000
            f(); n += 1
        end
        best = min(best, (time_ns() - t0) / n)
    end
    best
end
for set in ("init_uniuni", "specs", "init_bibi", "allo_bibi")
    conv = deserialize(joinpath(DATA, set * ".jls"))
    for c in conv[1:min(3, end)]
        m = c.mech
        fwd = S.orient(m)
        mask = S.inactive_mask(m)
        cols = S.columns(m, fwd, mask)
        al = @allocated S.solve(m)
        @printf("%-12s %-14s E=%2d cols=%3d  orient %7.2f  mask %7.2f  cycles %7.2f  columns %7.2f  system %7.2f  solve %7.2f µs  (solve allocates %d KB)\n",
                set, c.id, length(m.edges), length(cols),
                cost(() -> S.orient(m)) / 1e3, cost(() -> S.inactive_mask(m)) / 1e3,
                cost(() -> S.cycles(m, fwd)) / 1e3, cost(() -> S.columns(m, fwd, mask)) / 1e3,
                cost(() -> S.constraint_system(m)) / 1e3, cost(() -> S.solve(m)) / 1e3, al ÷ 1024)
    end
end
