# ABOUTME: Verifier phantom recheck on a fresh population: pre-filter, BigFloat ranks of the variant
# ABOUTME: and its L→0 law, the certificate, keys, fit sources and fits removed. Args: <lo> <hi> <out>.
using Printf, Random, Serialization
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
include(joinpath(SP, "slice/src/Slice.jl"))
const S = Slice
include(joinpath(VD, "sliceside.jl"))
include(joinpath(VD, "analib.jl"))
setprecision(BigFloat, 256)
const LO, HI, OUT = parse(Int, ARGS[1]), parse(Int, ARGS[2]), ARGS[3]
T0 = time()

lines = readlines(joinpath(VD, "data", "phantom_pop.tsv"))
res = Dict{Int, Any}()
nrow = 0
for (li, l) in enumerate(lines)
    LO <= li <= HI || continue
    f = split(l, '\t')
    gid, round, parent, move, allo, n = parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]),
                                         Symbol(f[4]), f[5] == "1", parse(Int, f[6])
    r = eval(Meta.parse(f[9]))
    out = Dict{Symbol, Any}(:gid => gid, :round => round, :parent => parent, :move => move,
                            :allo => allo, :n => n)
    try
        sm = mkmech(S, r)
        if !allo
            out[:key] = S.canonical(S.law(sm)).key
        else
            sol = S.solve(sm)
            A = S.allo_law(sm, sol)
            out[:key] = S.canonical(A).key
            if n == 1 && round == 2
                rng = Xoshiro(1000 + gid)
                pm = mkmech(S, r; plain = true)
                out[:admit] = S.admit(sm) === nothing
                im = S.inactive_mask(sm)
                out[:pf_empty] = !any(im)
                out[:pf_bind] = !any(k -> im[k] && (sm.edges[k].cons != 0 || sm.edges[k].rel != 0),
                                     eachindex(im))
                out[:pf_sites] = all(s -> all(==(:OnlyA), s.states), sm.sites)
                fv = collect(S.fitted_params(sm))
                LP = S.law(pm)
                fp = collect(S.fitted_params(pm))
                out[:np], out[:np0] = length(fv), length(fp)
                rv, okv, kv, dv = law_rank(A.vars, fv, allo_parts(S, A), rng)
                r0, ok0, k0, d0 = law_rank(LP.vars, fp, plain_parts(LP), rng)
                out[:r_var], out[:r_0] = rv, r0
                out[:draws_agree] = okv && ok0
                out[:kept] = Float64(min(kv, k0))
                out[:dropped] = Float64(max(dv, d0))
                cA, c0 = S.canonical(A.active), S.canonical(LP)
                out[:key_act], out[:key0] = cA.key, c0.key
                c = certificate(S, A)
                out[:cert] = c isa Dict
                out[:cert_why] = c isa Dict ? :ok : c
                if c isa Dict && cA.key == c0.key
                    # the map v(θ, L) = v0(θ ⊙ (1+L)^w) at two points, names matched through
                    # the two laws' canonical orders
                    nmap = Dict(zip(c0.order, cA.order))
                    worst = 0.0
                    for _ in 1:2
                        x = draw(A.vars, rng)
                        Lv = x[A.vars.index[:L]]
                        num, den = allo_parts(S, A)(x, log.(x), Dict{Int, Int}(), 0)
                        y = draw(LP.vars, rng)
                        for (i, nm) in enumerate(LP.vars.names)
                            LP.vars.conc[i] || nm === :Keq || haskey(nmap, nm) || continue
                            vn = haskey(nmap, nm) ? nmap[nm] : nm
                            j = A.vars.index[vn]
                            y[i] = x[j] * (haskey(c, vn) ? (1 + Lv)^BigFloat(c[vn]) : big(1.0))
                        end
                        q0, d0p = plain_parts(LP)(y, log.(y), Dict{Int, Int}(), 0)
                        worst = max(worst, Float64(abs((num.v / den.v) / (q0.v / d0p.v) - 1)))
                    end
                    out[:cert_err] = worst
                end
            end
        end
        out[:status] = :ok
    catch e
        out[:status] = :error
        out[:err] = first(sprint(showerror, e), 200)
    end
    res[gid] = out
    global nrow += 1
    nrow % 500 == 0 && (@printf("[%4.0f s, %4.0f MB] %d rows (line %d)\n", time() - T0,
                                Sys.maxrss() / 2^20, nrow, li); flush(stdout))
    time() - T0 > 540 && (println("stopping for time at line $li"); break)
end
serialize(OUT, res)
open(replace(OUT, ".jls" => ".tsv"), "w") do io
    for g in sort!(collect(keys(res)))
        println(io, join((string(k, "=", replace(string(v), r"[\t\n]" => " ")) for (k, v) in sort!(collect(res[g]); by = first)), "\t"))
    end
end
@printf("done %d rows in %.0f s, maxrss %.0f MB\n", nrow, time() - T0, Sys.maxrss() / 2^20)
