# ABOUTME: Arbitrates kcat with the independent oracle: the seven plain mechanisms where the slice
# ABOUTME: and the package disagree, the four "never saturates" ones, and ter-ter "depends" ones.
include(joinpath(@__DIR__, "spotlib.jl"))

function main()
    rng = Xoshiro(5)
    want = [("expand_bibi", "expand_bibi/$i") for i in (388, 442, 600, 609, 638, 663)]
    push!(want, ("p1_adv_pkg", "adv_pkg/12"))
    append!(want, [("handwritten", "handwritten/$i") for i in (3, 6, 10, 20)])
    cache = Dict{String, Any}()
    for (set, id) in want
        recs = get!(() -> skeptic_records(set), cache, set)
        r = only(x for x in recs if x.p.id == id)
        d = S.derive(r.c.mech)
        x = check(r, d, rng)
        @printf("%-14s %-18s kcat %-12s rel. diff %.1e (law vs oracle %.1e)\n", set, id, x.kc,
                x.kerr, x.errB)
    end
    rows = [x for k in 1:11
            for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls"))
            if x.ok && x.kcat isa String && occursin("depends on how", x.kcat)]
    println("ter-ter laws whose kcat 'depends on how the substrates grow': ", length(rows))
    for x in rows[randperm(rng, length(rows))[1:8]]
        ps = deserialize(joinpath(IDATA_IN, "pop_$(x.set).jls"))
        r = (set = x.set, p = ps[x.index], c = S.Converted(ps[x.index], x.index))
        d = S.derive(r.c.mech)
        y = check(r, d, rng)
        @printf("%-20s %-22s kcat %-14s law vs oracle %.1e\n", x.set, x.id, y.kc, y.errB)
    end
end

main()
