# ABOUTME: Codegen check (d) on the plain records of the sets named in ARGS: generated code
# ABOUTME: against the BigFloat law and the oracle, BigFloat inputs, zeros, half exponents.
include(joinpath(@__DIR__, "lib.jl"))

const NPTS = 4
rows = NamedTuple[]
notes = String[]
note(s) = (push!(notes, s); length(notes) <= 60 && println(s))
for set in ARGS
    t0 = time()
    for r in load_records((set,))
        m, p = r.c.mech, r.p
        S.is_allosteric(m) && continue
        rng = Xoshiro(hash(p.id, hash(set)))
        L = try
            S.law(m)
        catch err
            note("$set $(p.id): no law: $(sprint(showerror, err))")
            continue
        end
        em = S.compile(m)
        # Float64 inputs against the BigFloat law, and BigFloat inputs through the same code
        e64, ebig, econd = 0.0, big(0), 0.0
        for _ in 1:NPTS
            pr, c = rand_params(m, rng), rand_concs(m, rng)
            v, vb = S.rate_equation(em, c, pr), law_big(L, c, pr)
            vB = S.rate_equation(em, map(BigFloat, c), map(BigFloat, pr))
            sp = prod(c[x] for x in m.substrates)
            pp = prod(c[x] for x in m.products) / pr.Keq
            e = Float64(abs(v - vb) / abs(vb))
            e64 = max(e64, e)
            econd = max(econd, e / ((sp + pp) / abs(sp - pp)))
            ebig = max(ebig, abs(vB - vb) / abs(vb))
        end
        # each metabolite at zero, then every product at zero
        pr, c = rand_params(m, rng), rand_concs(m, rng)
        nz, nfin, ezero, zbad = 0, 0, 0.0, String[]
        for z in [[[x] for x in S.metabolites(m)]; [m.products]]
            cz = merge(c, (; (x => 0.0 for x in z)...))
            v, vb = S.rate_equation(em, cz, pr), law_big(L, cz, pr)
            nz += 1
            nfin += isfinite(v)
            ok = isfinite(vb) ? (iszero(vb) ? iszero(v) || abs(v) < 1e-300 :
                                 abs(v - vb) / abs(vb) < 1e-12) : !isfinite(v)
            ok || push!(zbad, "$(join(z, ",")): $v vs $(Float64(vb))")
            isfinite(vb) && !iszero(vb) && (ezero = max(ezero, Float64(abs(v - vb) / abs(vb))))
        end
        isempty(zbad) || note("$set $(p.id): zeros: " * join(zbad, "; "))
        # the oracle at consistent points (BigFloat constants, Float64 code)
        eor = 0.0
        for _ in 1:2
            pars = merge(random_params(m, rng), (E_total = big(1),))
            k = S.pop_constants(r.c, edge_values(m, S.all_parameters(m, pars)))
            u = (0.3 + 1.7 * rand(rng)) * (rand(rng, Bool) ? 1 : -1)
            conc = Oracle.concs_at_ratio(p, rng, pars.Keq * big(10)^u)
            vo = Oracle.steady_state_v(p, k, conc)
            v = S.rate_equation(em, (; (x => Float64(conc[x]) for x in S.metabolites(m))...),
                                map(Float64, pars))
            eor = max(eor, Float64(abs(v - vo) / abs(vo)))
        end
        b = S.body(S.program(m))
        half = any(!isinteger(e) for P in (L.H, L.Q, L.D, L.sigma..., L.trees...)
                   for mo in keys(P.terms) for (_, e) in mo)
        sq = occursin("sqrt(", string(b))
        half == sq || note("$set $(p.id): half exponents $half but sqrt in body $sq")
        max(e64, eor) > 1e-12 && note("$set $(p.id): Float64 rel $e64, oracle rel $eor")
        ebig > 1e-60 && note("$set $(p.id): BigFloat inputs rel $(Float64(ebig))")
        push!(rows, (set, id = p.id, e64, econd, ebig = Float64(ebig), nz, nfin, ezero,
                     nzbad = length(zbad), eor, half, sq, ops = countops(b)))
    end
    @printf("%s: %d rows so far, %.0f s\n", set, length(rows), time() - t0)
    flush(stdout)
end
serialize(joinpath(@__DIR__, "data", "values_$(join(ARGS, "+")).jls"), rows)
@printf("%d plain laws: Float64 vs BigFloat law worst %.1e (above 1e-12: %d; worst / \
        condition %.1e); BigFloat inputs worst %.1e; oracle worst %.1e (above 1e-12: %d)\n",
        length(rows), maximum(r -> r.e64, rows), count(r -> r.e64 > 1e-12, rows),
        maximum(r -> r.econd, rows), maximum(r -> r.ebig, rows), maximum(r -> r.eor, rows),
        count(r -> r.eor > 1e-12, rows))
@printf("zeros: %d evaluations, finite %d, matching the law %d, worst rel %.1e\n",
        sum(r -> r.nz, rows), sum(r -> r.nfin, rows), sum(r -> r.nz - r.nzbad, rows),
        maximum(r -> r.ezero, rows))
hv = filter(r -> r.half, rows)
println("half exponents: ", length(hv), " laws, sqrt in body: ", count(r -> r.sq, hv), ": ",
        join(("$(r.set) $(r.id) (Float64 $(r.e64), oracle $(r.eor))" for r in hv), ", "))
println("notes: ", length(notes))
println("maxrss ", round(Sys.maxrss() / 2^20), " MB")
