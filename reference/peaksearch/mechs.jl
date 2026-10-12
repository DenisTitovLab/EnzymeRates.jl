# ABOUTME: Mechanism variants for the peak search: textbook random-order bi-bi and structural
# ABOUTME: variants (lumped release, RE release, isomerization, parallel routes, abortive complexes).

const S = Symbol[]
ss(a, c, b, r) = (a, c, b, r, false)
re(a, c, b, r) = (a, c, b, r, true)

const BINDING = [ss(:E, [:A], :EA, S), ss(:E, [:B], :EB, S),
                 ss(:EA, [:B], :EAB, S), ss(:EB, [:A], :EAB, S)]
const RELEASE = [ss(:EPQ, S, :EQ, [:P]), ss(:EPQ, S, :EP, [:Q]),
                 ss(:EQ, S, :E, [:Q]), ss(:EP, S, :E, [:P])]

const MECHS = Dict{String, PeakSolver.Mech}(
    # E+A<->EA; E+B<->EB; EA+B<->EAB; EB+A<->EAB; EAB<->EPQ; random release, all SS
    "T" => mech("textbook", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ],
                vcat(BINDING, [ss(:EAB, S, :EPQ, S)], RELEASE), [5]),
    # (a) catalysis and release lumped: EAB <-> E + P + Q
    "a" => mech("lumped", [:E, :EA, :EB, :EAB],
                vcat(BINDING, [ss(:EAB, S, :E, [:P, :Q])]), [5]),
    # (b) rapid-equilibrium product release, steady-state binding and catalysis
    "b" => mech("REreleases", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ],
                vcat(BINDING, [ss(:EAB, S, :EPQ, S)],
                     [re(s[1], s[2], s[3], s[4]) for s in RELEASE]), [5]),
    # (c) ternary-complex isomerization before catalysis: EAB <-> EstarAB <-> EPQ
    "c" => mech("isomerization", [:E, :EA, :EB, :EAB, :EstarAB, :EPQ, :EP, :EQ],
                vcat(BINDING, [ss(:EAB, S, :EstarAB, S), ss(:EstarAB, S, :EPQ, S)], RELEASE),
                [6]),
    # (d1) parallel ternary complexes X1 = EAB (from EA+B), X2 = EstarAB (from EB+A),
    #      each catalyzing to a shared EPQ, random release
    "d1" => mech("parallel_sharedEPQ", [:E, :EA, :EB, :EAB, :EstarAB, :EPQ, :EP, :EQ],
                 vcat([ss(:E, [:A], :EA, S), ss(:E, [:B], :EB, S),
                       ss(:EA, [:B], :EAB, S), ss(:EB, [:A], :EstarAB, S),
                       ss(:EAB, S, :EPQ, S), ss(:EstarAB, S, :EPQ, S)], RELEASE), [5, 6]),
    # (d2) parallel ternary complexes, each catalyzing and releasing in one step
    "d2" => mech("parallel_lumped", [:E, :EA, :EB, :EAB, :EstarAB],
                 [ss(:E, [:A], :EA, S), ss(:E, [:B], :EB, S),
                  ss(:EA, [:B], :EAB, S), ss(:EB, [:A], :EstarAB, S),
                  ss(:EAB, S, :E, [:P, :Q]), ss(:EstarAB, S, :E, [:P, :Q])], [5, 6]),
    # (e) dead-end abortive complex EQ + A <-> EAQ
    "e" => mech("deadend_EAQ", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ, :EAQ],
                vcat(BINDING, [ss(:EAB, S, :EPQ, S)], RELEASE, [ss(:EQ, [:A], :EAQ, S)]), [5]),
    # (f) abortive complex that can release Q: EQ + A <-> EAQ <-> EA + Q
    "f" => mech("abortive_releasing", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ, :EAQ],
                vcat(BINDING, [ss(:EAB, S, :EPQ, S)], RELEASE,
                     [ss(:EQ, [:A], :EAQ, S), ss(:EAQ, S, :EA, [:Q])]), [5]),
    # (g) random binding, ordered release (the package's "Random-order Bi-Bi" fixture)
    "g" => mech("ordered_release", [:E, :EA, :EB, :EAB, :EPQ, :EQ],
                vcat(BINDING, [ss(:EAB, S, :EPQ, S), ss(:EPQ, S, :EQ, [:P]), ss(:EQ, S, :E, [:Q])]),
                [5]),
)
