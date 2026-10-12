# ABOUTME: Fit-cost calibration, step 2: fits each sampled candidate as identify does (20 CMA-ES
# ABOUTME: restarts, maxtime 600 s each, relative loss) on synthetic data, logging per-fit and per-restart cost.
using EnzymeRates, Serialization, Random, Printf, Statistics, LinearAlgebra
isdefined(Main, :ER) || Core.eval(Main, :(const ER = EnzymeRates))
isdefined(Main, :FCD) || Core.eval(Main, :(const FCD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost"))
isdefined(Main, :fit_instrumented) || include(joinpath(FCD, "fitlib.jl"))

const KEQ = 10.0
fc_err(e) = replace(first(sprint(showerror, e), 200), '\n' => ' ', '\t' => ' ')

const FIT_COLS = (:rxn, :gid, :cls, :allo, :allo_n, :round, :move, :np, :np_slice, :rank,
    :cx, :dterms, :qterms, :pkg_ops, :t_pass1, :t_problem, :t_first_loss, :t_loss_med, :t_loss_min,
    :t_fit, :fit_bytes, :fit_gctime, :fit_compile, :n_restarts, :truncated, :gens_total,
    :gens_med, :gens_min, :gens_max, :fevals_total, :lambda, :t_restart1, :t_restart_med,
    :t_restart_max, :n_success, :n_stalled, :n_maxtime, :n_maxiters, :n_other, :best_loss,
    :best_ret, :n_near_best, :t_rescale, :rescale_err, :t_fit_pkg, :pkg_compile, :pkg_gctime,
    :pkg_bytes, :loss_pkg, :ret_pkg, :error)

"""Fit one candidate the way identify does and return its row (a Dict over FIT_COLS):
identify's pass-1 calls, the `FittingProblem`, the first `loss!` call, `fit_rate_equation`
itself (cold), then the same restarts timed one by one; append the restarts to `io_rs`."""
function measure(name, row, m, data, opt, io_rs; maxtime, n_restarts, fit_cap)
    out = Dict{Symbol, Any}(c => "" for c in FIT_COLS)
    for c in (:rxn, :gid, :cls, :allo, :allo_n, :round, :move, :np, :np_slice, :rank, :cx,
              :dterms, :qterms, :pkg_ops)
        out[c] = getfield(row, c)
    end
    seed = row.gid + (name === :bi_bi_rxn ? 0 : 1_000_000)
    try
        Random.seed!(seed)
        # identify's PASS 1 for one fresh mechanism (`_process_batch`)
        t = time_ns()
        ER._eq_complexity(m)
        em = ER.compile_mechanism(m)
        pn = ER.fitted_params(em)
        eq_text = ER.rate_equation_string(em)
        ER._mechanism_type_string(em)
        ER._rate_eq_dedup_key(eq_text)
        out[:t_pass1] = (time_ns() - t) / 1e9
        np = length(pn)
        np == row.np || error("fitted count $np differs from the selection's $(row.np)")
        t = time_ns()
        fp = ER.FittingProblem(ER.compile_mechanism(m), data; Keq = KEQ, scale_k_to_kcat = 1.0)
        out[:t_problem] = (time_ns() - t) / 1e9
        x = clamp.(randn(np) .* 2.0, -15.0, 15.0)
        t = time_ns()
        ER.loss!(x, fp)
        out[:t_first_loss] = (time_ns() - t) / 1e9
        out[:t_loss_med], out[:t_loss_min], _ = time_loss(fp, np, Xoshiro(seed))
        # identify's PASS 2 fit of this candidate, cold (first fit of its problem type)
        Random.seed!(seed)
        sp = @timed ER.fit_rate_equation(fp, opt; n_restarts, maxtime)
        out[:t_fit_pkg] = sp.time
        out[:pkg_compile] = hasproperty(sp, :compile_time) ? sp.compile_time : -1.0
        out[:pkg_gctime] = sp.gctime
        out[:pkg_bytes] = sp.bytes
        out[:loss_pkg] = sp.value.loss
        out[:ret_pkg] = sp.value.retcode
        # the same restarts again (same seed), timed one by one
        Random.seed!(seed)
        st = @timed fit_instrumented(fp, opt; n_restarts, maxtime, fit_cap)
        f = st.value
        rs = f.restarts
        out[:t_fit] = st.time
        out[:fit_bytes] = st.bytes
        out[:fit_gctime] = st.gctime
        out[:fit_compile] = hasproperty(st, :compile_time) ? st.compile_time : -1.0
        out[:n_restarts] = length(rs)
        out[:truncated] = f.truncated
        gens = [r.gens for r in rs]
        ts = [r.t for r in rs]
        out[:gens_total] = sum(gens)
        out[:gens_med] = median(gens)
        out[:gens_min], out[:gens_max] = extrema(gens)
        out[:fevals_total] = sum(r.fevals for r in rs)
        out[:lambda] = rs[1].lambda
        out[:t_restart1] = ts[1]
        out[:t_restart_med] = median(ts)
        out[:t_restart_max] = maximum(ts)
        out[:n_success] = count(r -> r.retcode === :Success, rs)
        out[:n_stalled] = count(r -> r.retcode === :Stalled, rs)
        out[:n_maxtime] = count(r -> r.retcode === :MaxTime, rs)
        out[:n_maxiters] = count(r -> r.retcode === :MaxIters, rs)
        out[:n_other] = length(rs) - out[:n_success] - out[:n_stalled] - out[:n_maxtime] -
                        out[:n_maxiters]
        out[:best_loss] = f.loss
        out[:best_ret] = f.retcode
        out[:n_near_best] = count(r -> r.loss <= f.loss * 1.01 + 1e-12, rs)
        for r in rs
            println(io_rs, join((name, row.gid, row.cls, r.restart, r.t, r.gens, r.lambda,
                                 r.fevals, r.retcode, r.reason, r.loss), '\t'))
        end
        flush(io_rs)
        rp = NamedTuple{pn}(Tuple(exp.(f.x)))
        t = time_ns()
        try
            ER.rescale_parameter_values(fp.mechanism, merge(rp, (Keq = fp.Keq, E_total = 1.0));
                                        scale_k_to_kcat = fp.scale_k_to_kcat)
        catch e
            out[:rescale_err] = fc_err(e)
        end
        out[:t_rescale] = (time_ns() - t) / 1e9
    catch e
        out[:error] = fc_err(e)
    end
    out
end

function main_fit(; budget = 4200.0, maxtime = 600.0, n_restarts = 20, fit_cap = 1200.0,
                  files = ["selection.jls"], tag = "",
                  ref_gid = Dict{Symbol, Int}(), data_seed = Dict(:bi_bi_rxn => 101, :ter_bi_rxn => 202))
    T0 = time()
    sels = [deserialize(joinpath(FCD, "data", f)) for f in files]
    sel = sels[1]            # its references make the data; every file adds candidates
    opt = CMAEvolutionStrategyOpt()
    @printf("fit phase: budget %.0f s, maxtime %.0f s per restart, %d restarts, cap %.0f s per fit\n",
            budget, maxtime, n_restarts, fit_cap)
    out_name(base) = joinpath(FCD, "data", isempty(tag) ? "$base.tsv" : "$(base)_$tag.tsv")
    datas = Dict{Symbol, Any}()
    # a candidate of any file can stand in for a reaction's reference (`ref_gid`)
    cand = Dict((s.name, r.gid) => (r, m) for sk in sels for s in sk for (r, m) in s.sel)
    for s in sel
        (rrow, rm) = haskey(ref_gid, s.name) ? cand[(s.name, ref_gid[s.name])] : s.ref
        em = ER.compile_mechanism(rm)
        pn = ER.fitted_params(em)
        rng = Xoshiro(data_seed[s.name])
        θ = NamedTuple{pn}(Tuple(10.0 .^ (rand(rng, length(pn)) .- 0.5)))
        data, V, scale = make_data(s.name, em, θ, KEQ, rng)
        datas[s.name] = data
        open(joinpath(FCD, "data", isempty(tag) ? "data_$(s.name).tsv" : "data_$(s.name)_$tag.tsv"), "w") do io
            cols = collect(keys(data))
            println(io, join([cols; :v_true], '\t'))
            for i in eachindex(data.Rate)
                println(io, join([[data[c][i] for c in cols]; V[i]], '\t'))
            end
        end
        println("data $(s.name): $(length(data.Rate)) points in $(length(unique(data.group))) groups; ",
                "reference #$(rrow.gid) θ = $θ; group scales $scale")
        println("  rate range ", extrema(data.Rate))
    end
    flush(stdout)
    # compile the fitter's shared code on a mechanism outside the sample
    for s in sel
        (wrow, wm) = s.warm
        t = time_ns()
        fp = ER.FittingProblem(ER.compile_mechanism(wm), datas[s.name]; Keq = KEQ,
                               scale_k_to_kcat = 1.0)
        Random.seed!(1)
        w = ER.fit_rate_equation(fp, opt; n_restarts = 2, maxtime)
        fit_instrumented(fp, opt; n_restarts = 1, maxtime)
        @printf("warm-up %s #%d (np %d): %.1f s, loss %.4g\n", s.name, wrow.gid, wrow.np,
                (time_ns() - t) / 1e9, w.loss)
    end
    t = time_ns()
    ov = cma_overhead(opt)
    @printf("CMA-ES per-generation overhead measured in %.1f s\n", (time_ns() - t) / 1e9)
    open(out_name("cma_overhead"), "w") do io
        println(io, "n\trep\tt\tgens\tlambda\treason")
        for o in ov
            println(io, join((o.n, o.rep, o.t, o.gens, o.lambda, o.reason), '\t'))
        end
    end
    flush(stdout)
    # per file, interleave the reactions; each starts with its reference, then a seeded
    # shuffle; a candidate already listed by an earlier file is skipped
    order = Any[]
    taken = Set{Tuple{Symbol, Int}}()
    for (k, sk) in enumerate(sels)
        lists = [vcat([(s.name, s.sel[1]...)],
                      [(s.name, x...) for x in s.sel[2:end][randperm(Xoshiro(4 + k), length(s.sel) - 1)]])
                 for s in sk]
        for i in 1:maximum(length, lists), l in lists
            i <= length(l) || continue
            key = (l[i][1], l[i][2].gid)
            key in taken && continue
            push!(taken, key)
            push!(order, l[i])
        end
    end
    println("candidates to fit: $(length(order)) from $(length(files)) selection file(s)")
    io_fit = open(out_name("fits"), "w")
    io_rs = open(out_name("restarts"), "w")
    println(io_fit, join(FIT_COLS, '\t'))
    println(io_rs, "rxn\tgid\tcls\trestart\tt\tgens\tlambda\tfevals\tretcode\treason\tloss")
    tf = time()
    done = 0
    for (name, row, m) in order
        if time() - tf > budget
            println("budget reached after $done fits; $(length(order) - done) not fitted")
            break
        end
        r = measure(name, row, m, datas[name], opt, io_rs; maxtime, n_restarts, fit_cap)
        println(io_fit, join((r[c] for c in FIT_COLS), '\t'))
        flush(io_fit)
        done += 1
        @printf("%3d %-10s #%-5d %-3s np=%-2d rank=%-2d ops=%-5d fit %7.2f s gens %7d (med %6.0f) ret S%d/St%d/T%d best %.4g %s%s\n",
                done, name, row.gid, row.cls, row.np, row.rank, row.pkg_ops,
                r[:t_fit] == "" ? NaN : r[:t_fit], r[:gens_total] == "" ? 0 : r[:gens_total],
                r[:gens_med] == "" ? NaN : r[:gens_med],
                r[:n_success] == "" ? 0 : r[:n_success], r[:n_stalled] == "" ? 0 : r[:n_stalled],
                r[:n_maxtime] == "" ? 0 : r[:n_maxtime],
                r[:best_loss] == "" ? NaN : r[:best_loss],
                r[:t_fit_pkg] == "" ? "" : @sprintf("[pkg fit %.2f s loss %.4g]", r[:t_fit_pkg], r[:loss_pkg]),
                r[:error] == "" ? "" : " ERROR " * r[:error])
        flush(stdout)
    end
    close(io_fit)
    close(io_rs)
    @printf("fit phase done: %d fits in %.0f s (phase total %.0f s); maxrss %.0f MB\n", done,
            time() - tf, time() - T0, Sys.maxrss() / 2^20)
end

# usage: fit.jl <budget s> <maxtime s> <fit cap s> <tag> <selection files, comma-separated>
#               [reference gids, e.g. bi_bi_rxn=49,ter_bi_rxn=296] [data seeds, same form]
pairs_arg(a) = Dict(Symbol(first(split(x, "="))) => parse(Int, last(split(x, "="))) for x in split(a, ","))
if abspath(PROGRAM_FILE) == @__FILE__()
    main_fit(; budget = parse(Float64, ARGS[1]), maxtime = parse(Float64, ARGS[2]),
             fit_cap = parse(Float64, ARGS[3]), tag = ARGS[4],
             files = String.(split(ARGS[5], ",")),
             ref_gid = length(ARGS) >= 6 ? pairs_arg(ARGS[6]) : Dict{Symbol, Int}(),
             data_seed = length(ARGS) >= 7 ? pairs_arg(ARGS[7]) :
                         Dict(:bi_bi_rxn => 101, :ter_bi_rxn => 202))
end
