# ABOUTME: Converts EnzymeRates mechanisms into Pop.PMech records holding the package's outputs.
# ABOUTME: Included by export_pop.jl after `using EnzymeRates` and Population.jl (module Pop).

const ER = EnzymeRates

_role(::ER.Substrate) = :S
_role(::ER.Product) = :P
_role(::ER.CompetitiveInhibitor) = :Inh
_role(::ER.AllostericRegulator) = :Reg
_pmet(x) = Pop.PMet(ER.name(x), _role(x))
_pmets(xs) = Pop.PMet[_pmet(x) for x in xs]

_pspecies(sp::ER.Species) = Pop.PSpecies(
    _pmets(ER.bound(sp)), ER.conformation(sp),
    Symbol[ER.name(a) for a in ER.added(ER.residual(sp))],
    Symbol[ER.name(p) for p in ER.subtracted(ER.residual(sp))], ER.name(sp))

_ppoly(p) = sort!(Pop.PPoly([Pop.PMono(collect(mono)) => c for (mono, c) in p]);
                  by = x -> string(x.first))

_errtext(e) = first(sprint(showerror, e), 600)

"""Log-uniform draw in [1e-2, 1e2]."""
_logu(rng) = 10.0^(4 * rand(rng) - 2)

"""Evaluate a dependent-parameter Expr (`*`, `/`, `^`, `+`, `-`, numbers, symbols)."""
_evalx(x::Number, vals) = Float64(x)
_evalx(x::Symbol, vals) = vals[x]
function _evalx(ex::Expr, vals)
    ex.head === :call || error("unexpected Expr head $(ex.head)")
    op = ex.args[1]
    a = [_evalx(y, vals) for y in ex.args[2:end]]
    op === :* ? prod(a) : op === :+ ? sum(a) : op === :/ ? a[1] / a[2] :
    op === :^ ? a[1]^a[2] : op === :- ? (length(a) == 1 ? -a[1] : a[1] - a[2]) :
    error("unexpected operator $op")
end

"""A dependent's power-product Expr (`build_power_expr`: `Keq`, symbols, `*`, `1 / x`,
`x ^ n`, a Float64 exponent for a non-integer power, or the literal `1.0`) as
`symbol => exponent` pairs sorted by symbol."""
function _expr_mono(ex)
    acc = Dict{Symbol, Rational{Int}}()
    function walk(x, f)
        if x isa Symbol
            acc[x] = get(acc, x, 0) + f
        elseif x isa Number
            isone(x) || error("constant factor $x in a dependent expression")
        elseif x isa Expr && x.head === :call && x.args[1] === :*
            foreach(a -> walk(a, f), x.args[2:end])
        elseif x isa Expr && x.head === :call && x.args[1] === :/
            walk(x.args[2], f); walk(x.args[3], -f)
        elseif x isa Expr && x.head === :call && x.args[1] === :^
            e = x.args[3]
            walk(x.args[2], f * (e isa Integer ? Rational{Int}(e) : rationalize(Int, e)))
        else
            error("unexpected dependent expression $x")
        end
    end
    walk(ex, 1 // 1)
    sort!([k => v for (k, v) in acc if v != 0]; by = first)
end

"""Per-call time of `rate_equation(em, c, p)` as the package's performance test measures
it (minimum over 5 batches of 10_000 calls), with the allocations of one call."""
function _bench(em, c, p)
    ER.rate_equation(em, c, p)
    allocs = @allocated ER.rate_equation(em, c, p)
    best = Inf
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += ER.rate_equation(em, c, p)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("rate_equation produced a non-finite result")
    best * 1e9, allocs
end

"""Feature tags of a mechanism for population summaries."""
function _family(groups, species, flat, allo_tags, sites, rxn_subs, rxn_prods)
    words = ("uni", "bi", "ter", "quad", "quint", "hexa")
    tags = Symbol[Symbol(words[length(rxn_subs)], "_", words[length(rxn_prods)])]
    any(sp -> !isempty(ER.residual(sp)), species) && push!(tags, :pingpong)
    any(s -> !isempty(ER.consumed(s)) && !isempty(ER.released(s)), flat) &&
        push!(tags, :theorell_chance)
    any(s -> ER.is_binding(s) && ER._is_chemistry(s), flat) && push!(tags, :fused_binding)
    any(s -> length(ER.consumed(s)) > 1 || length(ER.released(s)) > 1, flat) &&
        push!(tags, :multi_metabolite_side)
    any(s -> any(x -> x isa ER.CompetitiveInhibitor,
                 [ER.consumed(s); ER.released(s)]), flat) && push!(tags, :inhibitor)
    deg = Dict{ER.Species, Int}()
    for s in flat, sp in (ER.from_species(s), ER.to_species(s))
        deg[sp] = get(deg, sp, 0) + 1
    end
    any(==(1), values(deg)) && push!(tags, :dead_end_leaf)
    nre = count(ER.is_equilibrium, flat)
    push!(tags, nre == 0 ? :all_ss : nre == length(flat) ? :all_re : :mixed_re_ss)
    any(g -> length(g) > 1, groups) && push!(tags, :multi_step_group)
    any(sp -> ER.conformation(sp) !== :E, species) && push!(tags, :conformations)
    for t in unique(allo_tags); push!(tags, Symbol("tag_", t)); end
    isempty(sites) || push!(tags, :reg_sites)
    tags
end

"""Median wall time in ns of rebuilding `m` through its constructor (`_with`)."""
function _ctor_ns(m; reps = 20)
    ER._with(m)
    ts = Float64[]
    for _ in 1:reps
        t0 = time_ns(); ER._with(m); push!(ts, Float64(time_ns() - t0))
    end
    sort!(ts)[(reps + 1) ÷ 2]
end

"""
`m` (a `Mechanism` or `AllostericMechanism`) as a `Pop.PMech`. Structure comes from the
mechanism the derivation sees (the lift of `compile_mechanism(m)`); every package output is
computed stage by stage, and a stage that throws is recorded in `stage_errors` (the first
fatal one also in `error`, with `ok = false`). `seed` fixes the numeric samples.
"""
function to_pmech(m, id::String, source::String; seed::Integer, nsamples = 4,
                  do_kcat = true, do_bench = true)
    d = Dict{Symbol, Any}(:id => id, :source => source)
    errs = Pair{Symbol, String}[]
    fail!(stage, e) = (push!(errs, stage => _errtext(e)); nothing)
    em = ER.compile_mechanism(m)
    lm = ER._concrete(em)
    allo = lm isa ER.AllostericMechanism
    groups = ER.steps(lm)
    flat = [s for g in groups for s in g]
    gof = [g for (g, grp) in enumerate(groups) for _ in grp]
    species = ER._forms(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    rxn = ER.reaction(lm)
    subs = Symbol[ER.name(x) for x in ER.substrates(rxn)]
    prods = Symbol[ER.name(x) for x in ER.products(rxn)]
    d[:substrates], d[:products] = subs, prods
    d[:inhibitors] = Symbol[ER.name(ER.regulator(r)) for r in ER.regulators(rxn)
                            if ER.regulator(r) isa ER.CompetitiveInhibitor]
    d[:regulators] = allo ? Symbol[ER.name(l) for l in ER.allosteric_regulators(lm)] :
                     Symbol[]
    d[:metabolites] = collect(Symbol, ER.metabolites(em))
    d[:species] = Pop.PSpecies[_pspecies(sp) for sp in species]
    d[:steps] = Pop.PStep[Pop.PStep(idx[ER.from_species(s)], idx[ER.to_species(s)],
                                    _pmets(ER.consumed(s)), _pmets(ER.released(s)),
                                    ER.is_equilibrium(s)) for s in flat]
    gl = Vector{Int}[]
    for g in eachindex(groups); push!(gl, findall(==(g), gof)); end
    d[:groups] = gl
    d[:allosteric] = allo
    tags = allo ? collect(Symbol, ER.cat_allo_states(lm)) : Symbol[]
    sites = allo ? ER.regulatory_sites(lm) : ER.RegulatorySite[]
    d[:group_tags] = tags
    d[:cat_multiplicity] = allo ? ER.catalytic_multiplicity(lm) : 1
    d[:sites] = Pop.PSite[Pop.PSite(Symbol[ER.name(l) for l in ER.ligands(s)],
                                    ER.multiplicity(s), collect(Symbol, ER.allo_states(s)))
                          for s in sites]
    d[:family] = _family(groups, species, flat, tags, sites, subs, prods)

    # Parameter names per flat step.
    if allo
        tA(g) = tags[g] === :EqualAI ? :EqualAI : :A
        tI(g) = tags[g] === :EqualAI ? :EqualAI : :I
        d[:step_names_A] = [Symbol[ER.name(p, lm) for p in ER._step_constants(s, tA(g))]
                            for (s, g) in zip(flat, gof)]
        d[:step_names_I] = [Symbol[ER.name(p, lm) for p in ER._step_constants(s, tI(g))]
                            for (s, g) in zip(flat, gof)]
        try
            ig = ER._inactive_groups(lm)
            d[:inactive_kept] = Bool[s in ig[g] for (s, g) in zip(flat, gof)]
        catch e
            fail!(:inactive_groups, e)
        end
        try
            cm, sp = ER._state_parts(lm, :I)
            names_pkg = [Symbol[] for _ in flat]
            rev = falses(length(flat))
            for (j, s) in enumerate(s for g in ER.steps(cm) for s in g)
                k = findfirst(t -> t == s || t == ER._reversed(s), flat)
                k === nothing && error("inactive step $(ER._step_text(s)) not in the mechanism")
                names_pkg[k] = Symbol[ER.name(p, cm) for p in sp[j]]
                rev[k] = flat[k] != s
            end
            d[:step_names_I_pkg] = names_pkg
            d[:inactive_reversed] = collect(rev)
        catch e
            fail!(:state_parts_I, e)
        end
    else
        d[:step_names] = [Symbol[ER.name(p, lm) for p in ER._step_constants(s, :None)]
                          for s in flat]
    end

    # The derivation: dependents, fitted parameters, the equation string and its key.
    dep = nothing
    try
        dep, indep = ER._dependent_param_exprs(lm)
        d[:dependents] = Pair{Symbol, String}[k => ER._expr_to_string(v)
                                              for (k, v) in sort!(collect(dep); by = first)]
        d[:dependent_monomials] = [k => _expr_mono(v)
                                   for (k, v) in sort!(collect(dep); by = first)]
        d[:fitted_params] = collect(Symbol, ER.fitted_params(em))
        collect(indep) == d[:fitted_params] ||
            push!(errs, :fitted_mismatch => "fitted_params != _dependent_param_exprs indep")
        d[:independent_param_count] = ER._independent_param_count(lm)
        s = ER.rate_equation_string(em)
        d[:rate_equation_string] = s
        d[:eq_hash] = ER._rate_eq_dedup_key(s)
    catch e
        fail!(:derivation, e)
        d[:ok] = false
        d[:error] = _errtext(e)
    end

    # Polynomials: plain mechanisms with no Wegscheider rename and with the package's own
    # default (Pass-2) rename; allosteric per-state polynomials as the package derives them.
    if !allo && dep !== nothing
        try
            num, den, dfree = ER._raw_symbolic_rate_polys(lm, ER._step_parameters(lm),
                                                          Dict{Symbol, Symbol}())
            pnum, pden, pdfree = ER._raw_symbolic_rate_polys(lm)
            d[:raw_num], d[:raw_den], d[:raw_dfree] = _ppoly(num), _ppoly(den), _ppoly(dfree)
            d[:pkg_num], d[:pkg_den], d[:pkg_dfree] = _ppoly(pnum), _ppoly(pden), _ppoly(pdfree)
            d[:pass2_changes_polys] = (num, den, dfree) != (pnum, pden, pdfree)
            d[:pass2_rename] = sort!(collect(ER._build_wegscheider_rename_map(lm)); by = first)
            d[:full_params] = collect(Symbol, ER.parameters(em, ER.Full))
        catch e
            fail!(:polys, e)
        end
    elseif allo && dep !== nothing
        for (st, key) in ((:A, :state_polys_A), (:I, :state_polys_I))
            try
                d[key] = Pop.PPoly[_ppoly(p) for p in ER._state_rate_polys(lm, st)]
                cm, sp = ER._state_parts(lm, st)
                d[st === :A ? :pass2_rename_A : :pass2_rename_I] =
                    sort!(collect(ER._build_wegscheider_rename_map(cm; step_params = sp));
                          by = first)
            catch e
                fail!(Symbol(:state_polys_, st), e)
            end
        end
    end
    try
        d[:eq_complexity] = Int(ER._eq_complexity(lm))
    catch e
        fail!(:eq_complexity, e)
    end

    # Numeric samples.
    if dep !== nothing
        rng = Random.Xoshiro(seed)
        mets = Tuple(d[:metabolites])
        fitted = Tuple(d[:fitted_params])
        fullp = Tuple(get(d, :full_params, Symbol[]))
        samples = Pop.PSample[]
        for _ in 1:nsamples
            c = NamedTuple{mets}(Tuple(_logu(rng) for _ in mets))
            p = NamedTuple{(fitted..., :Keq, :E_total)}(
                (Tuple(_logu(rng) for _ in fitted)..., _logu(rng), 1.0))
            fp = NamedTuple{fullp}(Tuple(k === :E_total ? 1.0 : _logu(rng) for k in fullp))
            v = try
                t0 = time()
                r = Float64(ER.rate_equation(em, c, p))
                isempty(samples) && (d[:first_call_s] = time() - t0)
                r
            catch e
                fail!(:rate_equation, e); NaN
            end
            vals = Dict{Symbol, Float64}(pairs(p))
            depv = Pair{Symbol, Float64}[]
            try
                for (k, ex) in sort!(collect(dep); by = first)
                    push!(depv, k => _evalx(ex, vals))
                end
            catch e
                fail!(:dep_values, e)
            end
            vf = NaN
            if !allo && !isempty(fullp)
                vf = try
                    Float64(ER.rate_equation(em, c, fp, ER.Full))
                catch e
                    fail!(:rate_equation_full, e); NaN
                end
            end
            kc = NaN
            if do_kcat
                kc = try
                    Float64(ER._kcat_forward(em, p))
                catch e
                    fail!(:kcat, e); NaN
                end
            end
            push!(samples, Pop.PSample(collect(pairs(c)), collect(pairs(p)), depv, v,
                                       collect(pairs(fp)), vf, kc))
        end
        d[:samples] = samples
        if do_bench && !isempty(samples) && !isnan(samples[1].v)
            try
                c1 = NamedTuple{mets}(Tuple(last.(samples[1].concs)))
                p1 = NamedTuple{(fitted..., :Keq, :E_total)}(Tuple(last.(samples[1].params)))
                d[:eval_ns], d[:eval_allocs] = _bench(em, c1, p1)
            catch e
                fail!(:bench, e)
            end
        end
    end
    bad = findfirst(x -> first(x) in (:rate_equation, :dep_values), errs)
    if bad !== nothing && get(d, :ok, true)
        d[:ok] = false
        d[:error] = last(errs[bad])
    end
    try
        d[:ctor_ns] = _ctor_ns(m)
    catch e
        fail!(:constructor, e)
    end
    d[:stage_errors] = unique(errs)
    Pop.PMech(; d...)
end
