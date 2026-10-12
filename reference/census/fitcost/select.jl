# ABOUTME: Fit-cost calibration, step 1: identify-like candidate pools for bi_bi_rxn and ter_bi_rxn
# ABOUTME: (seeds, two expansion rounds), their classes (rank, saturation, guard) and a stratified sample.
const FCD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small/lib.jl")
include(joinpath(FCD, "fitlib.jl"))

let ast = Meta.parseall(read(joinpath(TESTDIR, "test_mechanism_enumeration.jl"), String))
    for ex in ast.args
        ex isa Expr && ex.head === :const && ex.args[1] isa Expr &&
            ex.args[1].args[1] === :ter_bi_rxn && Core.eval(Main, ex)
    end
end

# The census's well-posedness verdict (census/small/report.jl `wclasses`): the sub-classes a
# record fails, empty when it passes.
function wclasses(r)
    w = String[]
    haskey(r, :slice_err) && return ["no slice law: " * first(r[:slice_err], 60)]
    get(r, :v0, false) && push!(w, "v=0")
    r[:re_turnover] === nothing || push!(w, "RE turnover")
    fk = r[:face_kind]
    for (i, f) in enumerate(r[:faces])
        fk[i] === :bottomless && push!(w, "0/0 bottomless @" * f)
        r[:dvan][i] && push!(w, "0/0 D vanishes @" * f)
        fk[i] === :split && push!(w, "split @" * f)
        fk[i] === :trap && push!(w, "trap @" * f)
        fk[i] === :other && push!(w, "face other @" * f)
    end
    r[:admit] === nothing || push!(w, "no :OnlyA limit")
    r[:role] === nothing || push!(w, "role change")
    r[:names_err] === nothing || push!(w, "names")
    r[:unused] === nothing || push!(w, "unused reactant")
    unique!(w)
end
nclass(r) = haskey(r, :vmax_fwd) && (!r[:vmax_fwd] || !r[:vmax_rev])
nnum(r) = get(r, :unb_fwd, false) || get(r, :unb_rev, false)

"""Slice and package facts (census lib) of the records `gids`, plus the package law's
arithmetic size (`pkg_ops`: call nodes of numerator, denominator and dependents)."""
function facts_for!(run, gids, rng)
    for g in gids
        rec, m = run.recs[g], run.mechs[g]
        haskey(rec, :facts_done) && continue
        rec[:facts_done] = true
        rec[:text] = steps_text(m)
        rec[:status] === :ok || continue
        try
            slice_facts!(rec, m)
        catch e
            rec[:slice_err] = errmsg(e)
        end
        try
            pl = pkg_facts!(rec, m, rng; rank = true)
            rec[:pkg_ops] = ncalls(pl.num) + ncalls(pl.den) +
                            sum((ncalls(r) for (_, r) in pl.dep); init = 0)
        catch e
            rec[:pkg_err] = errmsg(e)
        end
    end
end

"""Census classes of the fitted records among `gids`: D (rename-proof key of an earlier
record), W, N, Z (key of its parent), P (package-law rank below the fitted count)."""
function classify!(run, gids)
    first_with = Dict{String, Int}()
    for g in sort(gids)
        r = run.recs[g]
        (r[:status] === :ok && get(r, :fitted, false) && haskey(r, :key)) || continue
        haskey(first_with, r[:key]) ? (r[:D] = first_with[r[:key]]) :
                                       (first_with[r[:key]] = g)
    end
    for g in gids
        r = run.recs[g]
        (r[:status] === :ok && get(r, :facts_done, false)) || continue
        r[:W] = wclasses(r)
        r[:N] = nclass(r) || nnum(r)
        p = r[:parent]
        r[:Z] = p != 0 && haskey(r, :key) && haskey(run.recs[p], :key) &&
                run.recs[p][:key] == r[:key]
        r[:P] = haskey(r, :rank_pkg) && r[:rank_pkg] >= 0 && r[:rank_pkg] < r[:np]
    end
end

eligible(r) = r[:status] === :ok && get(r, :fitted, false) && haskey(r, :W) &&
              haskey(r, :pkg_ops) && haskey(r, :rank_pkg) && r[:rank_pkg] >= 0 &&
              isempty(r[:W]) && !haskey(r, :D) && !r[:Z] && r[:kcat_err] === nothing

"""`v` scaled to [0, 1] (all zeros when constant)."""
function unit_scale(v)
    lo, hi = extrema(v)
    hi > lo ? (v .- lo) ./ (hi - lo) : zero(v)
end

"""Farthest-point sample of `k` records in (parameter count, log law size), both scaled to
[0, 1], starting from the largest law; then random records up to `k`."""
function pick(recs, k, rng)
    k = min(k, length(recs))
    k == 0 && return eltype(recs)[]
    X = hcat(unit_scale([Float64(r[:np]) for r in recs]),
             unit_scale([log(r[:pkg_ops]) for r in recs]))
    chosen = [argmax(X[:, 2])]
    d = [norm(X[i, :] .- X[chosen[1], :]) for i in axes(X, 1)]
    while length(chosen) < cld(k, 2)
        i = argmax(d)
        d[i] <= 0 && break
        push!(chosen, i)
        d .= min.(d, [norm(X[j, :] .- X[i, :]) for j in axes(X, 1)])
    end
    rest = setdiff(eachindex(recs), chosen)
    append!(chosen, rest[randperm(rng, length(rest))[1:min(k - length(chosen), length(rest))]])
    recs[chosen]
end

quant(v) = isempty(v) ? "-" : join((@sprintf("%g", quantile(v, q)) for q in (0, 0.25, 0.5, 0.75, 0.9, 1)), "/")

function build(name; npar1, npar2, nfacts0, nP = 12, ntotal = 40, seed = 20261011,
               pick_seed = nothing)
    rxn = getfield(Main, name)
    run = Run(rxn, string(name))
    rng = Xoshiro(seed)
    println("== $name"); flush(stdout)
    t0 = time()
    f0 = round0!(run)
    fr0 = frontier(run, 0)
    par1 = length(fr0) <= npar1 ? fr0 : fr0[sort!(randperm(rng, length(fr0))[1:npar1])]
    f1 = expand_round!(run, par1, 1)
    fr1 = frontier(run, 1)
    par2 = length(fr1) <= npar2 ? fr1 : fr1[sort!(randperm(rng, length(fr1))[1:npar2])]
    f2 = expand_round!(run, par2, 2)
    t = time() - t0
    for f in (f0, f1, f2)
        println("  round $(f[:round]): parents=$(get(f, :parents, "-")) fresh=$(f[:fresh]) ",
                "cx_skip=$(f[:cx_skip]) param_skip=$(f[:param_skip]) error=$(f[:error]) ",
                "ok=$(f[:ok]) eqhash_dup=$(f[:eqhash_dup]) fitted=$(f[:fitted])")
    end
    @printf("  pools built in %.1f s\n", t); flush(stdout)
    fitted(round) = [r[:gid] for r in run.recs if r[:round] == round && r[:status] === :ok &&
                     get(r, :fitted, false)]
    g0, g1, g2 = fitted(0), fitted(1), fitted(2)
    s0 = length(g0) <= nfacts0 ? g0 :
         sort!(union(g0[randperm(rng, length(g0))[1:nfacts0]], intersect(g0, par1)))
    fg = sort!(union(s0, g1, g2, par2))
    t0 = time()
    facts_for!(run, fg, rng)
    t = time() - t0
    classify!(run, fg)
    @printf("  facts of %d records (%d of %d round-0 fitted) in %.1f s\n", length(fg),
            length(s0), length(g0), t); flush(stdout)
    for (round, gs) in ((0, s0), (1, g1), (2, g2))
        rs = [run.recs[g] for g in gs if haskey(run.recs[g], :W)]
        n = length(rs)
        c(f) = count(f, rs)
        println("  round $round with facts: $n; W $(c(r -> !isempty(r[:W]))), D $(c(r -> haskey(r, :D))), ",
                "Z $(c(r -> r[:Z])), P $(c(r -> r[:P])), N $(c(r -> r[:N])), ",
                "eligible $(c(eligible)) (P $(c(r -> eligible(r) && r[:P])))")
        println("    np q0/25/50/75/90/100: ", quant([r[:np] for r in rs]),
                "; pkg_ops: ", quant([r[:pkg_ops] for r in rs if haskey(r, :pkg_ops)]),
                "; slice D terms: ", quant([r[:dterms] for r in rs if haskey(r, :dterms)]))
    end
    E = [run.recs[g] for g in fg if eligible(run.recs[g])]
    # the reference: a plain, full-rank round-0 law one parameter above the round-0 median,
    # median size among those
    r0 = [r for r in E if r[:round] == 0 && !r[:allo] && !r[:P] && !r[:N]]
    nps = sort([r[:np] for r in r0])
    target = nps[cld(length(nps), 2)] + 1
    cands = [r for r in r0 if r[:np] == target]
    isempty(cands) && (cands = r0)
    sort!(cands; by = r -> (r[:pkg_ops], r[:gid]))
    ref = cands[cld(length(cands), 2)]
    println("  reference #$(ref[:gid]) np=$(ref[:np]) pkg_ops=$(ref[:pkg_ops]) dterms=$(ref[:dterms]): $(ref[:text])")
    # the sample (its own random stream when `pick_seed` is given; the pools and facts
    # above do not depend on it)
    pick_seed === nothing || (rng = Xoshiro(pick_seed))
    Pset = [r for r in E if r[:P]]
    Nset = [r for r in E if r[:N] && !r[:P]]
    Oset = [r for r in E if !r[:P] && !r[:N] && r[:gid] != ref[:gid]]
    Ppl, Pal = [r for r in Pset if !r[:allo]], [r for r in Pset if r[:allo]]
    kpl = min(length(Ppl), max(nP ÷ 2, nP - length(Pal)))
    selP = vcat(pick(Ppl, kpl, rng), pick(Pal, nP - kpl, rng))
    selN = pick(Nset, 5, rng)
    selO = pick(Oset, ntotal - 1 - length(selP) - length(selN), rng)
    println("  sample: OK $(length(selO)) + reference, P $(length(selP)) ($(length(Ppl)) plain, ",
            "$(length(Pal)) allosteric in the pool), N $(length(selN)) (pool $(length(Nset)))")
    row(r, cls) = (rxn = name, gid = r[:gid], cls = cls, round = r[:round], move = r[:move],
                   parent = r[:parent], allo = r[:allo],
                   allo_n = r[:allo] ? ER.catalytic_multiplicity(run.mechs[r[:gid]]) : 1,
                   np = r[:np], np_slice = get(r, :np_slice, -1), cx = r[:cx],
                   dterms = get(r, :dterms, -1), qterms = get(r, :qterms, -1),
                   pkg_ops = r[:pkg_ops], rank = r[:rank_pkg],
                   gap_kept = r[:rank_pkg_gap][1], gap_dropped = r[:rank_pkg_gap][2],
                   text = r[:text])
    sel = vcat([(row(ref, "REF"), run.mechs[ref[:gid]])],
               [(row(r, "OK"), run.mechs[r[:gid]]) for r in selO],
               [(row(r, "P"), run.mechs[r[:gid]]) for r in selP],
               [(row(r, "N"), run.mechs[r[:gid]]) for r in selN])
    # a plain record outside the sample, fitted first to compile the fitter's shared code
    taken = Set(r.gid for (r, _) in sel)
    warm = first(r for r in Oset if !r[:allo] && !(r[:gid] in taken))
    # the population the sample stands for: every eligible record's np and size
    popn = [(round = r[:round], allo = r[:allo], P = r[:P], np = r[:np], pkg_ops = r[:pkg_ops],
             dterms = get(r, :dterms, -1), cx = r[:cx]) for r in E]
    (name = name, ref = (row(ref, "REF"), run.mechs[ref[:gid]]), sel = sel,
     warm = (row(warm, "WARM"), run.mechs[warm[:gid]]), pop = popn)
end

function main_select(; file = "selection", ntotal = 40, nP = 12, pick_seed = nothing)
    out = [build(:bi_bi_rxn; npar1 = 10^6, npar2 = 30, nfacts0 = 10^6, ntotal, nP, pick_seed),
           build(:ter_bi_rxn; npar1 = 25, npar2 = 15, nfacts0 = 400, ntotal, nP, pick_seed)]
    serialize(joinpath(FCD, "data", "$file.jls"), out)
    open(joinpath(FCD, "data", "$file.tsv"), "w") do io
        cols = (:rxn, :gid, :cls, :round, :move, :parent, :allo, :allo_n, :np, :np_slice, :cx,
                :dterms, :qterms, :pkg_ops, :rank, :gap_kept, :gap_dropped, :text)
        println(io, join(cols, '\t'))
        for s in out, (r, _) in s.sel
            println(io, join((getfield(r, c) for c in cols), '\t'))
        end
    end
    open(joinpath(FCD, "data", "population_$file.tsv"), "w") do io
        println(io, "rxn\tround\tallo\tP\tnp\tpkg_ops\tdterms\tcx")
        for s in out, p in s.pop
            println(io, join((s.name, p.round, p.allo, p.P, p.np, p.pkg_ops, p.dterms, p.cx), '\t'))
        end
    end
    println("selection saved: ", sum(length(s.sel) for s in out), " candidates")
    flush(stdout)
end

# usage: select.jl <file stem> <candidates per reaction> <rank-deficient per reaction> <pick seed>
abspath(PROGRAM_FILE) == @__FILE__() &&
    main_select(; file = ARGS[1], ntotal = parse(Int, ARGS[2]), nP = parse(Int, ARGS[3]),
                pick_seed = parse(Int, ARGS[4]))
