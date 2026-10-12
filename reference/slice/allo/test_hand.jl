# ABOUTME: Hand-built admission cases, one per verdict: a directed cycle of :OnlyA bindings (G1),
# ABOUTME: a chord (G1), a kinetic island, a second free form, a plain mechanism and an admitted one.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
include(joinpath(SCRATCH, "slice", "src", "Slice.jl"))
using Test
const S = Slice
include(joinpath(@__DIR__, "g1_lp.jl"))

f(conf, bound...; res = Symbol[]) =
    S.RawForm(conf, [S.Met(x, false) for x in bound], [x for x in res if x in (:A, :B)],
              [x for x in res if x === :P])
step(a, b, cons, rel; re = true) =
    S.RawStep(a, b, [S.Met(x, false) for x in cons], [S.Met(x, false) for x in rel], re)
mech(steps, group, tags) = S.Mechanism([:A, :B], [:P], steps, group; tags,
                                        multiplicity = 2)

E, EA, EB, EAB, EP = f(:E), f(:E, :A), f(:E, :B), f(:E, :A, :B), f(:E, :P)
bi = [step(E, EA, [:A], []), step(EA, EAB, [:B], []), step(E, EB, [:B], []),
      step(EB, EAB, [:A], []), step(EAB, EP, [], []; re = false), step(EP, E, [], [:P])]
@testset "admission verdicts" begin
    # E + A → E(A) and E(A) + B → E(A, B) are :OnlyA; E(B) joins E and E(A, B) in the
    # inactive state, so the two arcs run E → E(A) → E: a directed cycle
    cyc = mech(bi, [1, 2, 3, 4, 5, 6], [:OnlyA, :OnlyA, :EqualAI, :EqualAI, :EqualAI,
                                       :EqualAI])
    msg = S.admit(cyc)
    @test occursin("close a cycle", msg) && !g1_lp(cyc)
    @test occursin("E + A ⇌ E(A)", msg) && occursin("E(A) + B ⇌ E(A, B)", msg)
    # the same with only E + A → E(A) :OnlyA: the kept path E – E(B) – E(A, B) – E(A) makes
    # it a chord
    chord = mech(bi, [1, 2, 3, 4, 5, 6], [:OnlyA, :EqualAI, :EqualAI, :EqualAI, :EqualAI,
                                         :EqualAI])
    @test occursin("still join the two sides of the :OnlyA binding E + A ⇌ E(A)",
                   S.admit(chord)) && !g1_lp(chord)
    # ping-pong with both half-reactions :OnlyA: the covalent forms are a kinetic island
    F = S.RawForm(:E, S.Met[], [:A], [:P])
    FB = S.RawForm(:E, [S.Met(:B, false)], [:A], [:P])
    pp = [step(E, EA, [:A], []), step(EA, F, [], [:P]; re = false), step(F, FB, [:B], []),
          step(FB, E, [], [:Q]; re = false)]
    island = S.Mechanism([:A, :B], [:P, :Q], pp, [1, 2, 3, 4];
                         tags = [:EqualAI, :OnlyA, :EqualAI, :OnlyA], multiplicity = 2)
    @test occursin("reach the free enzyme only through :OnlyA steps", S.admit(island))
    @test g1_lp(island)
    # bindings of A and B :OnlyA with chemistry and P release kept: E(A, B) stays joined to E
    # through E(P), so the bindings' arcs E → E(A) → E close a cycle (r_A + r_B is fixed)
    loop = mech(bi, [1, 2, 3, 4, 5, 6], [:OnlyA, :OnlyA, :OnlyA, :OnlyA, :EqualAI, :EqualAI])
    @test occursin("close a cycle", S.admit(loop)) && !g1_lp(loop)
    # chemistry :OnlyA as well: every island lies below E along bindings
    ok = mech(bi, [1, 2, 3, 4, 5, 6], [:OnlyA, :OnlyA, :OnlyA, :OnlyA, :OnlyA, :EqualAI])
    @test S.admit(ok) === nothing && g1_lp(ok)
    # a second free form F joined to E by a steady-state isomerization
    two = mech([bi; step(E, f(:F), [], []; re = false)], [1, 2, 3, 4, 5, 6, 7],
               fill(:EqualAI, 7))
    @test startswith(S.admit(two), "the conformation flips at one free form, E; the " *
                                   "free form F is not joined")
    @test S.admit(S.Mechanism([:A, :B], [:P], bi, [1, 2, 3, 4, 5, 6])) === nothing
end
