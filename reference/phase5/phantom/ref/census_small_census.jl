# ABOUTME: Runs the census for one reaction: round 0 (identify's base tier), round 1 (every
# ABOUTME: round-0 frontier entry expanded), round 2 (≤60 sampled round-1 parents); saves records.
# usage: census.jl <reaction const names, comma-separated> [rank sample cap per round] [round-2 parents]
include(joinpath(@__DIR__, "lib.jl"))

const NAMES = Symbol.(split(ARGS[1], ","))
const RANK_CAP = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : typemax(Int)
const NPAR2 = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 60
const T0 = time()

"""Slice and package facts of every ok candidate of `round`; ranks for the fitted ones (a
seeded sample of at most `RANK_CAP` of them, the rest marked unranked)."""
function facts!(run::Run, round::Int, rng)
    gids = [r[:gid] for r in run.recs if r[:round] == round]
    fitted = [g for g in gids if get(run.recs[g], :fitted, false)]
    ranked = Set(length(fitted) <= RANK_CAP ? fitted :
                 fitted[sort!(randperm(Xoshiro(round + 7), length(fitted))[1:RANK_CAP])])
    tl = time()
    for (k, g) in enumerate(gids)
        rec, m = run.recs[g], run.mechs[g]
        rec[:text] = steps_text(m)
        rec[:status] === :ok || continue
        sm = sol = L = br = nothing
        try
            sm, sol, L, br = slice_facts!(rec, m)
        catch e
            rec[:slice_err] = errmsg(e)
        end
        try
            pkg_facts!(rec, m, rng; rank = g in ranked)
        catch e
            rec[:pkg_err] = errmsg(e)
        end
        if g in ranked && L !== nothing
            try
                t = time_ns()
                r = slice_rank(L, br, rng)
                rec[:t_rank_slice] = (time_ns() - t) / 1e9
                rec[:rank_slice] = r.rank
                rec[:rank_slice_n] = r.n
                rec[:rank_slice_gap] = (r.kept, r.dropped)
                rec[:involved] = r.involved
                if r.rank < r.n
                    kinds = kinds_of(sm, sol, L, r.names)
                    rec[:kinds] = Dict(zip(r.names, kinds))
                    parts, agree = deficit_parts(r.J, kinds, r.rank)
                    rec[:deficit_parts] = parts
                    rec[:f64_rank_agrees] = agree
                end
            catch e
                rec[:rank_err] = errmsg(e)
            end
        end
        rec[:ranked] = g in ranked
        if k % 500 == 0
            @printf("    facts %d/%d of round %d (%.0f s)\n", k, length(gids), round,
                    time() - tl)
            flush(stdout)
        end
    end
end

function show_funnel(f)
    println("  round $(f[:round]): ", join(("$k=$(f[k])" for k in
        (:seeds, :degenerate_seeds, :flip_children_added, :parents, :raw_children,
         :generated, :seen_skip, :fresh, :cx_skip, :param_skip, :error, :ok, :eqhash_dup,
         :fitted, :expand_failures) if haskey(f, k)), ", "),
        @sprintf(" (pipeline %.1f s)", f[:t_pipeline]))
    flush(stdout)
end

function census(name)
    rxn = getfield(Main, name)
    run = Run(rxn, string(name))
    rng = Xoshiro(20261011)
    println("reaction $name"); flush(stdout)
    f0 = round0!(run)
    show_funnel(f0)
    t = @elapsed facts!(run, 0, rng)
    @printf("  facts round 0: %.1f s\n", t)
    f1 = expand_round!(run, frontier(run, 0), 1)
    show_funnel(f1)
    t = @elapsed facts!(run, 1, rng)
    @printf("  facts round 1: %.1f s\n", t)
    par1 = frontier(run, 1)
    pick = length(par1) <= NPAR2 ? par1 :
           par1[sort!(randperm(Xoshiro(60), length(par1))[1:NPAR2])]
    f2 = expand_round!(run, pick, 2)
    f2[:parent_pool] = length(par1)
    show_funnel(f2)
    t = @elapsed facts!(run, 2, rng)
    @printf("  facts round 2: %.1f s\n", t)
    for r in run.recs
        r[:allo_n] = r[:allo] ? ER.catalytic_multiplicity(run.mechs[r[:gid]]) : 1
    end
    out = NPAR2 == 60 ? string(name) : string(name, ".r2_", NPAR2)
    serialize(joinpath(CEN, "data", "$out.jls"),
              (label = run.label, funnel = run.funnel, recs = run.recs, rank_cap = RANK_CAP,
               npar2 = NPAR2))
    serialize(joinpath(CEN, "data", "$out.mechs.jls"), run.mechs)
    @printf("done %s in %.0f s; maxrss %.0f MB\n", name, time() - T0, Sys.maxrss() / 2^20)
end

foreach(census, NAMES)
