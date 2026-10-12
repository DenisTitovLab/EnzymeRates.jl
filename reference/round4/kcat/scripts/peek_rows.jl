# ABOUTME: Prints the stored per-point values (kcat, full-law limit, sup and its argmax) of a few
# ABOUTME: laws named on the command line as set/index.
using Serialization, Printf
const KDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/data"
rows = [x for l in ("core", "expand_terter", "terter_sample") for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok]
for a in ARGS
    set, idx = split(a, "/")
    x = only(filter(x -> x.set == set && x.index == parse(Int, idx), rows))
    println(x.id, " (", a, ") class ", x.cls, " full ", x.fcls, " ratio ", round(x.ratio; digits = 3))
    for p in x.pts
        @printf("    kcat %.4g full-law limit %.4g sup %.4g at %s\n", p.kc, p.kf, p.sup,
                sort!([k => round(v; sigdigits = 3) for (k, v) in p.at]))
    end
end
