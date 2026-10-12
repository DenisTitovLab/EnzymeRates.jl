# ABOUTME: Samples init_terter parts 5-11, which the probe's census did not reach: kcat classes,
# ABOUTME: exact positivity and per-substrate gaps of the kcat-zero and no-unique-limit laws.
include(joinpath(@__DIR__, "vlib.jl"))

rng = Xoshiro(511)
t0 = time()
taken = Set((s, i) for (s, i, _) in deserialize(joinpath(KDATA, "terter_sample.jls")))
out = Any[]
for k in 5:11
    set = "init_terter.part$k"
    rows = deserialize(joinpath(IDATA, "pipe_terter$k.jls"))
    okidx = [x.index for x in rows if x.ok && !((set, x.index) in taken)]
    pick = okidx[randperm(rng, length(okidx))[1:min(400, end)]]
    for r in pop_records(set, pick)
        ok, why, sol, L = derive_law(r.m)
        ok || (push!(out, (part = k, index = r.index, ok = false)); continue)
        cls, _, _ = kcat_class(r.m, sol)
        zf = cls in (:zero, :nounique) ? zero_facts(r.m, L) : nothing
        zc = cls in (:zero, :nounique) ? zero_facts(r.m, S.law(r.m, sol, S.catalytic_core(r.m))) : nothing
        push!(out, (part = k, index = r.index, ok = true, cls, zf, zc, nok = length(okidx)))
    end
    @printf("part %d: %d derivable rows, %d sampled, %.0f s\n", k, length(okidx), length(pick), time() - t0)
    flush(stdout)
end
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (y = f(x); d[y] = get(d, y, 0) + 1), xs);
              sort!(collect(d); by = string))
ok = filter(x -> x.ok, out)
println("sampled laws ", length(out), ", derivable ", length(ok))
println("class per part: ", cnt(x -> (x.part, x.cls), ok))
z = filter(x -> x.cls === :zero, ok)
println("zero: ", length(z), "; exact v>0 at P=0: ", count(x -> x.zf.pos && x.zf.num0 && x.zf.den0, z),
        "; any gap > 0: ", count(x -> any(g -> g !== nothing && g > 0, x.zf.gap), z),
        "; core-law culprits: ", cnt(x -> count(g -> g !== nothing && g < 0, x.zc.gap), z),
        "; full-law culprits: ", cnt(x -> count(g -> g !== nothing && g < 0, x.zf.gap), z))
nu = filter(x -> x.cls === :nounique, ok)
println("nounique: ", length(nu), "; core-law gap patterns: ", cnt(x -> sort([g === nothing ? 99 : Int(g) for g in x.zc.gap]), nu))
# weighted estimate of the zero share over parts 5-11 (weights: derivable rows per part)
w = Dict(x.part => x.nok for x in ok)
est = sum(w[k] * count(x -> x.part == k && x.cls === :zero, ok) / count(x -> x.part == k, ok) for k in keys(w)) / sum(values(w))
@printf("estimated kcat-zero share over parts 5-11 (by part weights): %.3f; done in %.0f s\n", est, time() - t0)
