# ABOUTME: Validates P7's flat and factored laws on every population against the package's samples
# ABOUTME: (BigFloat) and, for plain mechanisms, by exact polynomial identity; prints size statistics.
include(joinpath(@__DIR__, "P7.jl"))
using Printf, Statistics
setprecision(BigFloat, 256)
const PI = P7.PI
pops = isempty(ARGS) ? Pop.populations() : ARGS
for pop in pops
    ms = Pop.load_population(pop)
    nbad = 0; maxerr_f = 0.0; maxerr_x = 0.0; nsplit = 0; nextract = 0; nexact = 0; nfail = 0
    tf = Int[]; tx = Int[]; failures = String[]
    for m in ms
        m.ok || continue
        local fl, fx
        try
            fl = P7.flat_law(m); fx = P7.factored_law(m; expand_den = get(ENV, "EXPAND_DEN", "0") == "1", blocks = get(ENV, "BLOCKS", "1") == "1")
        catch e
            nfail += 1; push!(failures, m.id * ": " * first(sprint(showerror, e), 200)); continue
        end
        for s in m.samples
            isnan(s.v) && continue
            v1 = P7.value(fl, P7.sample_vals(fl, s)); v2 = P7.value(fx, P7.sample_vals(fx, s))
            e1 = Float64(abs(v1 - s.v) / abs(s.v)); e2 = Float64(abs(v2 - s.v) / abs(s.v))
            maxerr_f = max(maxerr_f, e1); maxerr_x = max(maxerr_x, e2)
            (e1 > 1e-12 || e2 > 1e-12) && (nbad += 1)
        end
        push!(tf, P7.node_terms(fl)); push!(tx, P7.node_terms(fx))
        if !m.allosteric
            np, deps = fx.np, fx.deps
            sub(p) = P7.substdeps(p, np, deps)
            ok = P7.expand_node(fx, fx.den) == P7.topoly(m.pkg_den, fx.vars) &&
                 sub(P7.expand_node(fx, fx.num)) == sub(P7.topoly(m.pkg_num, fx.vars))
            ok ? (nexact += 1) : push!(failures, m.id * ": factored != package polys")
            sp = P7.haldane_split(sub(P7.topoly(m.pkg_num, fx.vars)), m, fx.vars)
            sp === nothing ? push!(failures, m.id * ": no Haldane split") : (nsplit += 1)
            hn = fx.nodes[fx.num][1]
            nextract += max(0, length(hn.facs) - 2)
        end
    end
    @printf("%-13s n=%4d fail=%d bad=%d maxerr flat %.1e fact %.1e exact %d Hsplit %d extracted %d terms flat %d fact %d\n",
            pop, length(ms), nfail, nbad, maxerr_f, maxerr_x, nexact, nsplit, nextract, sum(tf), sum(tx))
    foreach(f -> println("   ", f), first(failures, 8))
end
