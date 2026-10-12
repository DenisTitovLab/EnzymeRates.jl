# ABOUTME: Codegen check (d) on the allosteric records of the sets in ARGS: generated MWC code
# ABOUTME: against the BigFloat AlloLaw and the two-conformation oracle, BigFloat inputs, zeros.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(SCRATCH, "slice", "allo", "lib.jl"))      # allo_value, oracle_v (MWC.jl)

const NPTS = 4
rows = NamedTuple[]
notes = String[]
note(s) = (push!(notes, s); length(notes) <= 60 && println(s))
for set in ARGS
    t0 = time()
    for r in load_records((set,))
        m, p = r.c.mech, r.p
        S.is_allosteric(m) || continue
        rng = Xoshiro(hash(p.id, hash(set)))
        A = try
            S.allo_law(m)
        catch err
            note("$set $(p.id): no law: $(sprint(showerror, err))")
            continue
        end
        em = S.compile(m)
        e64, econd, ebig = 0.0, 0.0, big(0)
        for _ in 1:NPTS
            pr, c = rand_params(m, rng), rand_concs(m, rng)
            v, vb = S.rate_equation(em, c, pr), allo_value(A, Dict(pairs(c)), pr)
            vB = S.rate_equation(em, map(BigFloat, c), map(BigFloat, pr))
            sp = prod(c[x] for x in m.substrates)
            pp = prod(c[x] for x in m.products) / pr.Keq
            e = Float64(abs(v - vb) / abs(vb))
            e64 = max(e64, e)
            econd = max(econd, e / ((sp + pp) / abs(sp - pp)))
            ebig = max(ebig, abs(vB - vb) / abs(vb))
        end
        pr, c = rand_params(m, rng), rand_concs(m, rng)
        nz, nfin, zbad = 0, 0, String[]
        for z in [[[x] for x in S.metabolites(m)]; [m.products]]
            cz = merge(c, (; (x => 0.0 for x in z)...))
            v, vb = S.rate_equation(em, cz, pr), allo_value(A, Dict(pairs(cz)), pr)
            nz += 1
            nfin += isfinite(v)
            ok = isfinite(vb) ? (iszero(vb) ? abs(v) < 1e-300 : abs(v - vb) / abs(vb) < 1e-12) :
                 !isfinite(v)
            ok || push!(zbad, "$(join(z, ",")): $v vs $(Float64(vb))")
        end
        isempty(zbad) || note("$set $(p.id): zeros: " * join(zbad, "; "))
        eor = 0.0
        for _ in 1:2
            pars = random_params(m, rng)
            allp = S.all_parameters(m, pars)
            u = (0.3 + 1.7 * rand(rng)) * (rand(rng, Bool) ? 1 : -1)
            conc = Oracle.concs_at_ratio(p, rng, pars.Keq * big(10)^u)
            vo = oracle_v(r, allp, conc)
            v = S.rate_equation(em, (; (x => Float64(conc[x]) for x in S.metabolites(m))...),
                                merge(map(Float64, pars), (E_total = 1.0,)))
            eor = max(eor, Float64(abs(v - vo) / abs(vo)))
        end
        b = S.body(S.program(m))
        half = any(!isinteger(e) for P in (A.H, A.active.Q, A.active.D, A.fA, A.fI)
                   for mo in keys(P.terms) for (_, e) in mo) ||
               (A.inactive !== nothing && any(!isinteger(e) for P in (A.inactive.Q, A.inactive.D)
                                              for mo in keys(P.terms) for (_, e) in mo))
        sq = occursin("sqrt(", string(b))
        half == sq || note("$set $(p.id): half exponents $half but sqrt in body $sq")
        max(e64, eor) > 1e-12 && note("$set $(p.id): Float64 rel $e64, oracle rel $eor")
        ebig > 1e-60 && note("$set $(p.id): BigFloat inputs rel $(Float64(ebig))")
        push!(rows, (set, id = p.id, e64, econd, ebig = Float64(ebig), nz, nfin,
                     nzbad = length(zbad), eor, half, sq, admitted = S.admit(m) === nothing,
                     n = m.multiplicity,
                     ops = countops(b), inactive = A.inactive !== nothing,
                     cross = !(S.nterms(A.fA) == 1 && S.nterms(A.fI) == 1)))
    end
    @printf("%s: %d rows so far, %.0f s\n", set, length(rows), time() - t0)
    flush(stdout)
end
serialize(joinpath(@__DIR__, "data", "allo_values_$(join(ARGS, "+")).jls"), rows)
@printf("%d MWC laws (%d admitted): Float64 vs BigFloat law worst %.1e (above 1e-12: %d; \
        worst / condition %.1e); BigFloat inputs worst %.1e; oracle worst %.1e (above 1e-12: \
        %d)\n", length(rows), count(r -> r.admitted, rows), maximum(r -> r.e64, rows),
        count(r -> r.e64 > 1e-12, rows), maximum(r -> r.econd, rows), maximum(r -> r.ebig, rows),
        maximum(r -> r.eor, rows), count(r -> r.eor > 1e-12, rows))
@printf("zeros: %d evaluations, finite %d, matching the law %d\n", sum(r -> r.nz, rows),
        sum(r -> r.nfin, rows), sum(r -> r.nz - r.nzbad, rows))
@printf("multiplicities %s; empty inactive state %d; polynomial cross weights %d; half \
        exponents %d (sqrt in body %d)\n", sort!(unique(r.n for r in rows)),
        count(r -> !r.inactive, rows), count(r -> r.cross, rows), count(r -> r.half, rows),
        count(r -> r.sq, rows))
println("notes: ", length(notes))
println("maxrss ", round(Sys.maxrss() / 2^20), " MB")
