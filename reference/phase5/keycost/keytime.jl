# ABOUTME: Times the slice law and the by-name key (and, with "code", codegen) per exported row;
# ABOUTME: resumable chunks. Args: <in.tsv> <first row> <budget s> <out.tsv> <key|code>.
include(joinpath(@__DIR__, "klib.jl"))

const IN, FIRST, BUDGET, OUT, MODE = ARGS[1], parse(Int, ARGS[2]), parse(Float64, ARGS[3]),
                                     ARGS[4], ARGS[5]

"""Largest exponent denominator among the solve's dependents (1 when all are integers)."""
maxden(sol) = maximum((denominator(e) for d in sol.deps for e in [d.keq; last.(d.terms)]);
                      init = 1)

"""One row's measurements as strings; a function barrier, so no loop reads a field of the
slice's `Union{Law{Int}, Law{BigInt}}` (RULINGS.md R18)."""
@noinline function measure(r)
    set, label, id, subs, prods, pa, pi_ = r[1], r[2], r[3], names_of(r[4]), names_of(r[5]),
                                          r[6], r[7]
    m = parse_neutral(r[end], subs, prods)
    a, i = cx_pair(m)
    sol = S.solve(m)
    tl = @elapsed L = try
        S.rate_law(m, sol)
    catch e
        return [set, label, id, pa, pi_, string(a), string(i), "error: " *
                first(replace(sprint(showerror, e), '\n' => ' '), 120)]
    end
    tk = @elapsed (c, leaves, nodes) = canonical_named(L)
    if tk > 0.005                                     # a second timing for the larger keys
        tk = min(tk, @elapsed (c2, _, _) = canonical_named(L))
        c2.key == c.key || error("key differs between two calls: $id")
    end
    np = count(n -> n !== :Keq && n !== :E_total, S.parameters(m))
    out = [set, label, id, pa, pi_, string(a), string(i), "ok", string(den_terms(L)),
           string(key_terms(L)), string(np), @sprintf("%.6f", tl), @sprintf("%.6f", tk),
           string(leaves), string(nodes), c.key[1:16], string(maxden(sol)),
           string(S.is_allosteric(m) ? m.multiplicity : 0)]
    if MODE == "code"
        tp = @elapsed ex = S.body(S.program(m, sol, L))
        em = S.compile(m)
        concs, params = random_inputs(m, Xoshiro(7))
        tc = @elapsed S.rate_equation(em, concs, params)
        tb, al = time_call(em, concs, params; barrier = true)
        ti, _ = time_call(em, concs, params; barrier = false)
        append!(out, [@sprintf("%.6f", tp), @sprintf("%.4f", tc), @sprintf("%.2f", 1e9 * tb),
                      @sprintf("%.2f", 1e9 * ti), string(al), string(body_ops(ex)),
                      bytes2hex(sha256(string(ex)))[1:16]])
    end
    out
end

rows = tsv(IN)
t0 = time()
measure(rows[FIRST])                                   # compile the measuring code once
i = FIRST
open(OUT, "a") do io
    global i
    while i <= length(rows) && time() - t0 < BUDGET
        println(io, join(measure(rows[i]), '\t'))
        flush(io)
        i += 1
    end
end
@printf("rows %d-%d of %d in %.0f s; maxrss %.0f MB\n", FIRST, i - 1, length(rows),
        time() - t0, Sys.maxrss() / 2^20)
println(i <= length(rows) ? "NEXT $i" : "DONE")
