# ABOUTME: Spot-checks random mechanisms of a population by the skeptic's independent route;
# ABOUTME: usage: spot.jl label count sets… (log on stdout, rows in skeptic/data/spot_<label>.jls)
include(joinpath(@__DIR__, "spotlib.jl"))

const LABEL, COUNT, SETS = ARGS[1], parse(Int, ARGS[2]), ARGS[3:end]
const OUT = joinpath(SCR, "slice", "skeptic", "data")

"""Checks up to `COUNT` derived laws of `SETS`, drawn at random; rows and rejections."""
function main()
    rng = Xoshiro(hash(LABEL) % UInt32)
    rows = NamedTuple[]
    t0 = time()
    recs = reduce(vcat, [skeptic_records(s) for s in SETS])
    order = randperm(rng, length(recs))
    derived, rejected = 0, String[]
    for i in order
        derived >= COUNT && break
        r = recs[i]
        m = r.c.mech
        d = try
            S.derive(m)
        catch err
            err isa ErrorException || rethrow()
            push!(rejected, "$(r.p.id): " * first(err.msg, 90))
            length(rejected) > 4 * COUNT && break
            continue
        end
        derived += 1
        x = try
            check(r, d, rng)
        catch err
            println("CHECK ERROR ", r.p.id, ": ", sprint(showerror, err)[1:min(end, 300)])
            continue
        end
        push!(rows, merge((id = r.p.id, set = r.set, allo = S.is_allosteric(m)), x))
        @printf("%-28s %s count %d/%d %s cond %.1e | B %.1e F %.1e P %.1e C %.1e dep %.1e | %s %s\n",
                r.p.id, S.is_allosteric(m) ? "A" : "P", x.mine, x.slice, x.count_ok ? "ok" : "DIFF",
                x.cond, x.errB, x.errF, x.errP, x.errC, x.worst_dep, x.kc,
                isnan(x.kerr) ? "" : @sprintf("%.1e", x.kerr))
        flush(stdout)
    end
    serialize(joinpath(OUT, "spot_$LABEL.jls"), (; rows, rejected))
    println("rejected while sampling: ", length(rejected))
    foreach(x -> println("  ", x), first(rejected, 10))
    fails = filter(x -> !x.count_ok || x.errB > 1e-25 || x.errP > 1e-25 || x.errC > 1e-25 ||
                        x.worst_dep > 1e-25 || x.errF > 1e-12 * max(1, x.cond) ||
                        occursin("REFUTED", x.kc) || x.kc == "kcat-DIFF", rows)
    @printf("%s: %d checked, %d failing; worst B %.1e P %.1e C %.1e dep %.1e F/cond %.1e; kcat %s; %.0f s\n",
            LABEL, length(rows), length(fails), maximum(x -> x.errB, rows; init = 0.0),
            maximum(x -> x.errP, rows; init = 0.0), maximum(x -> x.errC, rows; init = 0.0),
            maximum(x -> x.worst_dep, rows; init = 0.0),
            maximum(x -> x.errF / max(1, x.cond), rows; init = 0.0),
            join(("$k=$(count(x -> x.kc == k, rows))" for k in unique(x.kc for x in rows)), " "),
            time() - t0)
    foreach(x -> println("FAIL ", x), fails)

end

main()
