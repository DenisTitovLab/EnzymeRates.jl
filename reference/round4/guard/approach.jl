# ABOUTME: For records whose law is 0/0 on a face (or disagrees with the oracle there): the law
# ABOUTME: along random directions into the face, and the oracle just off the face.
include(joinpath(@__DIR__, "lib.jl"))

"""The slice mechanism of a plain-data spec."""
function from_spec(sp)
    met(x) = S.Met(x[1], x[2])
    form(f) = S.RawForm(f[1], S.Met[met(x) for x in f[2]], collect(Symbol, f[3]), collect(Symbol, f[4]))
    steps = [S.RawStep(form(s[1]), form(s[2]), S.Met[met(x) for x in s[3]],
                       S.Met[met(x) for x in s[4]], s[5]) for s in sp.steps]
    S.Mechanism(sp.substrates, sp.products, steps, sp.group; tags = sp.tags,
                multiplicity = sp.multiplicity,
                sites = [S.Site(s[1], s[2], s[3]) for s in sp.sites])
end

sl = deserialize(joinpath(R, "data", "slice_full.jls"))
sp = Dict(r.key => r.spec for r in deserialize(joinpath(R, "data", "specs_full.jls")))
lawv(s) = haskey(s, :law) && haskey(s.law, :reduced) ? s.law[:reduced] : nothing
picks = String[]
for s in sl
    t = lawv(s)
    t === nothing && continue
    (t.dfwd || t.drev) && push!(picks, s.key)
end
append!(picks, ["p2_fuzz_dump/5", "p2_fuzz_dump_flip/9"])
println("records: ", length(picks))
rng = Xoshiro(5)
setprecision(BigFloat, 256)
summary = Dict{String, Int}()
for key in picks
    m = from_spec(sp[key])
    sol = S.solve(m)
    rb = reduced_branches(m, sol)
    p, vals = random_params(m, rng)
    merge!(vals, Dict(pairs(p)))
    base = Dict{Symbol, BigFloat}(x => logu(rng, -1, 1) for x in S.metabolites(m))
    foreach(l -> base[l] = big(0.0), S.regulators(m))
    t = lawv(only(filter(s -> s.key == key, sl)))
    faces = [(f, zs) for (f, zs, bad) in ((:fwd, m.products, t.dfwd), (:rev, m.substrates, t.drev))
             if bad || startswith(key, "p2_fuzz")]
    for (f, zs) in faces
        N1, D1 = n1_branches(m, rb)
        vs = BigFloat[]
        out = String[]
        for d in 1:3
            u = Dict(z => logu(rng, -1, 1) for z in zs)
            row = String[]
            for ε in (big(1e-4), big(1e-8), big(1e-12))
                c = copy(base)
                foreach(z -> c[z] = ε * u[z], zs)
                x = varvals(rb.vars, c, vals)
                v = evalbranches(N1, x) / evalbranches(D1, x)
                push!(row, @sprintf("%.6g", Float64(v)))
                ε == big(1e-12) && push!(vs, v)
            end
            push!(out, join(row, " → "))
        end
        c0 = copy(base); foreach(z -> c0[z] = big(0.0), zs)
        o0 = try oracle_v(m, vals, c0, big(10)^20, ones(BigFloat, length(m.edges) + 1)) catch e; big(NaN) end
        c1 = copy(base); foreach(z -> c1[z] = big(1e-30), zs)
        o1 = try oracle_v(m, vals, c1, big(10)^60, ones(BigFloat, length(m.edges) + 1)) catch e; big(NaN) end
        spread = maximum(vs) - minimum(vs)
        kind = spread > 1e-6 * maximum(abs.(vs)) + 1e-12 ? "direction-dependent" : "removable (one limit)"
        summary["$kind"] = get(summary, "$kind", 0) + 1
        println("$key $f: $kind; law along 3 directions (ε=1e-4,1e-8,1e-12): ", join(out, " | "),
                @sprintf("; oracle at face %.6g, at 1e-30 with Λ=1e60 %.6g", Float64(o0), Float64(o1)))
    end
end
println("summary: ", summary)
