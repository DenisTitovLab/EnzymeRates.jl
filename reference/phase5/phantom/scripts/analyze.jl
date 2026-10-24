# ABOUTME: Per-candidate L-phantom facts on the slice types: plain keys, and for each one-subunit
# ABOUTME: allosteric variant the pre-filter, modular and BigFloat ranks, certificate, slot map, L→0 key.
# usage: analyze.jl <data dir> <reaction> <first row> <last row> <out tsv> [bigfloat: 1|0]
include(joinpath(@__DIR__, "plib.jl"))
include(joinpath(SP, "phase3", "census", "census_common.jl"))

const DATA, RXN = ARGS[1], ARGS[2]
const LO, HI = parse(Int, ARGS[3]), parse(Int, ARGS[4])
const OUTF = ARGS[5]
const BIG = length(ARGS) >= 6 ? ARGS[6] == "1" : true
const SUBS, PRODS = startswith(RXN, "uniuni") ? ([:S], [:P]) : ([:A, :B], [:P, :Q])
const T0 = time()

rows = read_rows(joinpath(DATA, "recs_$RXN.tsv"))
const HEAD = rows[1]
col(r, f) = r[findfirst(==(f), HEAD)]
done = Set{Int}()
isfile(OUTF) && foreach(l -> push!(done, parse(Int, split(l, '\t')[1])), eachline(OUTF))

errtext(e) = replace(first(sprint(showerror, e), 160), r"[\t\n]" => " ")

"""Whether the L → 0 law `LA` (variant names) equals the plain law `LP` under the edge-slot
map `pmap` (plain name => variant name), tested at three random points of each field."""
function ident_slot(LP, LA, pmap, rng)
    BP, BA = blaw(LP), blaw(LA)
    Q = lcm(denlcm(BP), denlcm(BA))
    for F in (F1, F2), _ in 1:3
        tP = Dict{Symbol, F}()
        for x in BP.vars.names
            tP[x] = F(rand(rng, 1:typemax(Int)))
        end
        tA = Dict{Symbol, F}()
        for x in BA.vars.names
            BA.vars.conc[BA.vars.index[x]] || x === :Keq || continue
            tA[x] = get(tP, x, F(rand(rng, 1:typemax(Int))))   # a site ligand: not in the plain law
        end
        for (p, a) in pmap
            tA[a] = tP[p]
        end
        for x in BA.vars.names
            haskey(tA, x) || return false        # a variant constant with no plain partner
        end
        value_fp(BP, tP, Q) == value_fp(BA, tA, Q) || return false
    end
    true
end

"""The largest relative error of `v(θ, L) = v₀(θ ⊙ (1 + L)^w)` over five random points in
256-bit BigFloat."""
function cert_error(A, LP, w, pmapinv, rng)
    BV, B0 = blaw(A), blaw(A.active)
    setprecision(BigFloat, 256) do
        worst = big(0.0)
        for _ in 1:5
            vals = Dict{Symbol, BigFloat}()
            for x in BV.vars.names
                vals[x] = big(10)^(2 * rand(rng) - 1)
            end
            u = 1 + vals[:L]
            v0vals = copy(vals)
            for (x, wx) in w
                v0vals[x] = vals[x] * u^(BigFloat(numerator(wx)) / denominator(wx))
            end
            a, b = value(BV, vals), value(B0, v0vals)
            worst = max(worst, abs(a / b - 1))
        end
        Float64(worst)
    end
end

function analyze(r)
    gid, round, move, parent = parse(Int, col(r, "gid")), col(r, "round"), col(r, "move"),
                               col(r, "parent")
    allo, n, text = col(r, "allo") == "1", parse(Int, col(r, "n")), col(r, "text")
    base = [string(gid), "", round, move, parent, col(r, "np")]
    t = time()
    if !allo
        base[2] = "plain"
        try
            pm = parse_neutral(text, SUBS, PRODS)
            L = S.law(pm)
            return [base; "ok"; string(length(S.fitted_params(pm))); canonical_named(L).key;
                    @sprintf("%.1f", 1e3 * (time() - t))]
        catch e
            return [base; "error: " * errtext(e); ""; ""; ""]
        end
    end
    n == 1 || return [base[1]; "allo_n$n"; base[3:end]]
    base[2] = "allo1"
    out = Dict{String, String}()
    rng = Xoshiro(gid)
    try
        sm = parse_neutral(text, SUBS, PRODS)
        pm = parse_neutral(text, SUBS, PRODS; plain = true)
        adm = S.admit(sm)
        out["admit"] = adm === nothing ? "ok" : first(adm, 60)
        im = S.inactive_mask(sm)
        out["pf_empty"] = string(Int(!any(im)))
        out["pf_bind"] = string(Int(!any(k -> im[k] && (sm.edges[k].cons != 0 ||
                                                         sm.edges[k].rel != 0), eachindex(im))))
        out["pf_sites"] = string(Int(all(s -> all(==(:OnlyA), s.states), sm.sites)))
        out["nsites"] = string(length(sm.sites))
        sol = S.solve(sm)
        A = S.allo_law(sm, sol)
        fv = collect(S.fitted_params(sm))
        out["np_slice"] = string(length(fv))
        BV = blaw(A)
        cv = colsof(BV, fv)
        out["r_var"] = string(rank_mod(BV, cv, rng))
        out["r_noL"] = string(rank_mod(BV, colsof(BV, filter(!=(:L), fv)), rng))
        out["r_var_rel"] = string(rank_mod(BV, cv, rng; ones = true))
        LP = S.law(pm)
        fp = collect(S.fitted_params(pm))
        out["np0"] = string(length(fp))
        BP = blaw(LP)
        cp = colsof(BP, fp)
        out["r_0"] = string(rank_mod(BP, cp, rng))
        out["r_0_rel"] = string(rank_mod(BP, cp, rng; ones = true))
        out["t_mod"] = @sprintf("%.1f", 1e3 * (time() - t))
        if BIG
            tb = time()
            rb, kept, dropped = rank_big(BV, cv, rng)
            out["rb_var"], out["rb_kept"], out["rb_drop"] = string(rb), @sprintf("%.1e", kept),
                                                           @sprintf("%.1e", dropped)
            rb0, kept0, dropped0 = rank_big(BP, cp, rng)
            out["rb_0"], out["rb0_kept"], out["rb0_drop"] = string(rb0), @sprintf("%.1e", kept0),
                                                           @sprintf("%.1e", dropped0)
            out["t_big"] = @sprintf("%.1f", 1e3 * (time() - tb))
        end
        # the edge-slot map between the plain and the allosteric solves (key.md Q3 option b)
        solP = S.solve(pm)
        out["edges_same"] = string(Int(pm.forms == sm.forms && pm.edges == sm.edges))
        slotP = Dict((solP.cols[i].edge, solP.cols[i].slot) => i for i in solP.indep)
        slotV = Dict((sol.cols[i].edge, sol.cols[i].slot) => i for i in sol.indep
                     if !sol.cols[i].inactive)
        out["n_inact_indep"] = string(count(i -> sol.cols[i].inactive, sol.indep))
        out["slot_same"] = string(Int(Set(keys(slotP)) == Set(keys(slotV))))
        nP, nV = S.column_names(pm, solP), S.column_names(sm, sol)
        pmap = Dict(nP[i] => nV[slotV[k]] for (k, i) in slotP if haskey(slotV, k))
        out["ident_slot"] = string(Int(ident_slot(LP, A.active, pmap, rng)))
        out["key0"] = canonical_named(LP).key
        out["key_act"] = canonical_named(A.active).key
        w = certificate(A)
        out["cert"] = w === nothing ? "0" : "1"
        if w !== nothing
            out["cert_err"] = @sprintf("%.1e", cert_error(A, LP, w, nothing, rng))
            out["cert_w"] = join(("$(k)=$(v)" for (k, v) in sort!(collect(w); by = first)), ",")
        end
        out["status"] = "ok"
    catch e
        out["status"] = "error: " * errtext(e)
    end
    out["t_ms"] = @sprintf("%.1f", 1e3 * (time() - t))
    [base; [string(k, "=", v) for (k, v) in sort!(collect(out); by = first)]]
end

open(OUTF, "a") do io
    nrow = 0
    for i in max(LO, 2):min(HI, length(rows))
        r = rows[i]
        gid = parse(Int, col(r, "gid"))
        (gid in done || col(r, "status") != "ok") && continue
        println(io, join(analyze(r), '\t'))
        flush(io)
        nrow += 1
        if nrow % 200 == 0
            @printf("%s rows to %d: %d done (%.0f s, maxrss %.0f MB)\n", RXN, i, nrow,
                    time() - T0, Sys.maxrss() / 2^20)
            flush(stdout)
        end
        time() - T0 > 540 && (println("stopping for time at row $i"); break)
    end
    @printf("%s finished %d rows in %.0f s, maxrss %.0f MB\n", RXN, nrow, time() - T0,
            Sys.maxrss() / 2^20)
end
