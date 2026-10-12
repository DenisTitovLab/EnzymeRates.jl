# ABOUTME: An :OnlyA binding whose bound side lies in free E's component (F + P <--> E(P)): the
# ABOUTME: deletion law against three ε-readings of :OnlyA (away from E, kinetic block, K_I → ∞).
include(joinpath(@__DIR__, "vlib.jl"))

const DF = Dict("ES" => "E(S)", "EP" => "E(P)")
V(f, c, t, r, re, tag) = VStep(f, c, t, r, re, tag)
const EQ, NE, OA = :EqualAI, :NonequalAI, :OnlyA

for (label, re4) in (("F + P <--> E(P) :: OnlyA (SS)", false), ("F + P ⇌ E(P) :: OnlyA (RE)", true))
    steps = [V("E", ["S"], "ES", String[], true, EQ), V("ES", String[], "EP", String[], false, EQ),
             V("E", ["P"], "EP", String[], true, EQ), V("F", ["P"], "EP", String[], re4, OA)]
    println("\n## $label")
    m = build(dsltext(["S"], ["P"], steps, DF))
    m isa String && (println("NOT CONSTRUCTED: ", m); continue)
    println("fitted_params: ", ER.fitted_params(m))
    inact = sort!(unique([String(ER.name(f)) for grp in ER._inactive_groups(ER.AllostericMechanism(m))
                          for s in grp for f in (ER.from_species(s), ER.to_species(s))]))
    println("inactive forms kept by the package: ", inact)
    rng = MersenneTwister(5)
    for draw in 1:2
        fitted = Dict(String(p) => (p === :L ? lu(rng, 0.2, 5) : lu(rng, 0.1, 10))
                      for p in ER.fitted_params(m))
        Keq = lu(rng, 0.5, 20)
        vals = solve_values(steps, fitted, Keq, ["S"], ["P"])
        k4 = stepconsts(steps[4], :A, y -> haskey(vals, y))
        # (rate toward E: F + P → E(P), rate away from E: E(P) → F + P), per unit [P] and form.
        Λb = big(10)^30
        on, off = re4 ? (Λb * k4.ratio(y -> vals[y]), Λb) : k4.rates(y -> vals[y])
        readings = (("away", ε -> (on, ε * off)), ("block", ε -> (ε * on, ε * off)),
                    ("K_I→∞", ε -> (ε * on, off)))
        for P in (big(1), big(10)^-3, big(10)^-6, big(0))
            c = Dict("S" => big(2), "P" => P)
            vdel = pkg_rate(m, fitted, Keq, c)
            line = @sprintf("draw %d S=2 P=%-6.0e deletion law %+.5e |", draw, P, vdel)
            for (nm, f) in readings, ε in (big(10)^-6, big(10)^-12)
                r = try
                    oracle_rate(steps, vals, c, fitted["L"], "S"; blocked = Dict(4 => f(ε)))
                catch err
                    err isa LinearAlgebra.SingularException || rethrow()
                    nothing
                end
                line *= r === nothing ? " $nm ε=$(Float64(ε)): singular" :
                        @sprintf(" %s ε=%.0e: rel. to deletion %.1e (I-mass %.3f)", nm, ε,
                                 abs(r[1] - vdel) / abs(vdel), r[2])
            end
            println(line)
        end
    end
end
