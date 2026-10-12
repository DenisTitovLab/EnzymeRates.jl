# ABOUTME: Task 3: PR #86's rejected two-free-form fixtures under flip-at-free-E: what they
# ABOUTME: derive, and whether deleting the stranded inactive forms is the ε → 0 limit.
include(joinpath(@__DIR__, "lib.jl"))

const EPS = (big(10)^-4, big(10)^-6, big(10)^-8)
include(joinpath(@__DIR__, "fixtures_def.jl"))

fmt(x) = @sprintf("%.1e", x)
println("package: ", pathof(EnzymeRates))

"""Concentration points: generic, products at 0, and per fixture extras."""
function points(rng, subs, prods, inh, label)
    gen() = Dict{Symbol, BigFloat}(Symbol(x) => draw_conc(rng) for x in [subs; prods; inh])
    pts = Pair{String, Dict{Symbol, BigFloat}}[]
    push!(pts, "generic" => gen())
    c = gen(); for x in prods; c[Symbol(x)] = 0; end
    push!(pts, "products 0" => c)
    if startswith(label, "R2")
        for p in (big(10)^-3, big(10)^-6)
            c = gen(); c[:P] = p
            push!(pts, "P = $(fmt(p))" => c)
        end
    end
    if !isempty(inh) && !startswith(label, "S1")
        c = gen(); for x in inh; c[Symbol(x)] = 0; end
        push!(pts, "inhibitors 0" => c)
    end
    pts
end

for (label, subs, prods, inh, steps) in fixtures, n in unique((written_n(label), 1))
    println("\n=== $label  (catalytic_multiplicity $n)")
    m = build(dsl(subs, prods, steps; n, inhibitors = inh))
    if m isa String
        println("  NOT CONSTRUCTED: ", m)
        continue
    end
    am = ER.AllostericMechanism(m)
    fp = ER.fitted_params(m)
    allforms = sort!(unique([ER.name(f) for g in ER.steps(am) for s in g
                             for f in (ER.from_species(s), ER.to_species(s))]))
    keptI = inactive_forms(am)
    stranded = setdiff(allforms, keptI)
    ia = island_analysis(am)
    if n == written_n(label)
        println("  constructs; stored steps:")
        foreach(l -> println("    ", l), step_lines(am))
        println("  fitted_params: ", fp)
        println("  inactive forms kept (package): ", isempty(keptI) ? "(free E alone)" : join(keptI, ", "),
                "; stranded: ", isempty(stranded) ? "-" : join(stranded, ", "))
        println("  island rule: G1 ", ia.g1, "; admitted behind bindings: ",
                isempty(ia.admitted) ? "-" : join(join.(ia.admitted, ","), " | "),
                "; islands (rejected): ", isempty(ia.islands) ? "-" : join(join.(ia.islands, ","), " | "))
        s = ER.rate_equation_string(m)
        println("  rate_equation_string ($(length(s)) chars):")
        foreach(l -> println("    | ", l), split(s, '\n'))
    end
    x1 = Symbol(subs[1])
    rng = MersenneTwister(hash((label, n)) % UInt32)
    ndraw = n == 1 ? 4 : 2
    for draw in 1:ndraw
        fitted = draw_fitted(rng, fp)
        Keq = exp(BigFloat(log(0.5) + log(40) * rand(rng)))
        vals = all_values(am, fitted, Keq)
        sA, sI, miss = oracle_steps(am, vals)
        isempty(miss) || println("  draw $draw: missing inactive constants for $miss")
        for (pname, concs) in points(rng, subs, prods, inh, label)
            v = pkg_rate(m, vals, concs)
            vo, left, mI = oracle_v(sA, sI, concs, vals[:L], x1; n)
            line = "  draw $draw, $(rpad(pname, 13)): v = $(@sprintf("%+.6e", v)) | deletion oracle " *
                   "$(fmt(relerr(v, vo))) | I mass $(@sprintf("%.3f", mI))"
            if !isempty(stranded)
                for mode in (:literal, :phi)
                    errs = String[]
                    for ε in EPS
                        sAe, sIe, misse = oracle_steps(am, vals; mode, ε)
                        isempty(misse) || error("missing inactive constants $misse")
                        ve, lefte, mIe = oracle_v(sAe, sIe, concs, vals[:L], x1; n, strict = true)
                        push!(errs, fmt(relerr(ve, v)))
                    end
                    # The same limit with every ε-scaled step's ratio 7× its active value.
                    sA7, sI7, _ = oracle_steps(am, vals; mode, ε = EPS[end], base = big(7))
                    v7, _, _ = oracle_v(sA7, sI7, concs, vals[:L], x1; n, strict = true)
                    line *= " | $mode ε=1e-4,-6,-8: " * join(errs, " ") *
                            ", base×7: " * fmt(relerr(v7, v))
                end
            end
            println(line)
        end
        # Thermodynamic self-check of each ε network at equilibrium concentrations.
        if !isempty(stranded) && draw == 1
            ce = eq_concs(rng, am, Keq, Dict{Symbol, BigFloat}(Symbol(x) => draw_conc(rng) for x in inh))
            vref = abs(pkg_rate(m, vals, points(rng, subs, prods, inh, label)[1].second))
            for mode in (:literal, :phi)
                sAe, sIe, _ = oracle_steps(am, vals; mode, ε = EPS[1])
                ve, _, _ = oracle_v(sAe, sIe, ce, vals[:L], x1; n, strict = true)
                println("  $mode ε=1e-4 network at equilibrium: |v|/|v_typ| = ", fmt(abs(ve) / vref))
            end
        end
    end
    # How fast the stranded steps vanish matters at products 0: ε^a for the non-binding
    # :OnlyA steps, ε for the bindings.
    if !isempty(stranded) && n == 1
        rng2 = MersenneTwister(7)
        fitted = draw_fitted(rng2, fp)
        Keq = exp(BigFloat(log(0.5) + log(40) * rand(rng2)))
        vals = all_values(am, fitted, Keq)
        for (pname, concs) in points(rng2, subs, prods, inh, label)[1:2]
            v = pkg_rate(m, vals, concs)
            row = String[]
            for a in (1 // 2, 1, 2), ε in EPS
                sAe, sIe, _ = oracle_steps(am, vals; mode = :literal, ε, a)
                ve, _, _ = oracle_v(sAe, sIe, concs, vals[:L], x1; n, strict = true)
                push!(row, "a=$(a) ε=$(fmt(ε)): $(fmt(relerr(ve, v)))")
            end
            println("  vanishing exponent, $(pname): ", join(row, "; "))
        end
    end
end
println("\ndone")
