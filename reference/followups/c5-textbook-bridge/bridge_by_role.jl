# ABOUTME: C5 follow-up: today's positional Segel bridge against a bridge by role (each k_i named by the
# ABOUTME: content and written direction of step i), on every spec with an analytical oracle. Read-only.
using EnzymeRates, LinearAlgebra, Random
const ER = EnzymeRates
const ORACLE = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/er_oracle"
const FIX = joinpath(ORACLE, "test", "mechanism_definitions_for_test_enzyme_derivation.jl")
const TRD = joinpath(ORACLE, "test", "test_rate_eq_derivation.jl")
include(FIX)

# ── Today's bridge, evaluated from the test file's AST without running its testsets ──
const WANT = Set([:_testhelper_positional_flat_idx, :_testhelper_positional_site_idx,
    :_testhelper_positional_params, :_testhelper_analytical_oracle_params,
    :_testhelper_compute_all_params, :_testhelper_eval_dep_expr,
    :_testhelper_random_independent_params_concs])
function defname(ex)
    ex isa Expr || return nothing
    ex.head === :macrocall && ex.args[1] isa GlobalRef && ex.args[1].name === Symbol("@doc") &&
        return defname(ex.args[end])
    ex.head in (:function, :(=)) || return nothing
    sig = ex.args[1]
    sig isa Expr && sig.head === :where && (sig = sig.args[1])
    sig isa Expr && sig.head === :call ? sig.args[1] : nothing
end
for ex in Meta.parseall(read(TRD, String)).args
    defname(ex) in WANT && Core.eval(Main, ex)
end

# ── The written blocks of the fixtures file, by spec name ──
walk(f, ex) = (f(ex); ex isa Expr && foreach(a -> walk(f, a), ex.args))
const SRC_MACROS = (Symbol("@_testhelper_enzyme_mechanism_src"),
                    Symbol("@_testhelper_allosteric_mechanism_src"))
function spec_blocks()
    out = Dict{String, Expr}()
    for blk in Meta.parseall(read(FIX, String)).args
        blk isa Expr && blk.head === :let || continue
        mac, nm = nothing, nothing
        walk(blk) do ex
            ex isa Expr || return
            ex.head === :macrocall && ex.args[1] in SRC_MACROS && mac === nothing &&
                (mac = ex.args[end])
            if ex.head === :call && ex.args[1] === :_testhelper_MechanismTestSpec
                for a in ex.args[2:end]
                    a isa Expr && a.head === :kw && a.args[1] === :name && (nm = a.args[2])
                end
            end
        end
        mac === nothing || nm === nothing || (out[nm] = mac)
    end
    out
end

# ── The bridge by role (this is the sketch; everything above is scaffolding) ──
# Labeled lines of a mechanism block: label => values (`substrates: A, B` => [:A, :B]).
function _testhelper_labeled(block)
    out = Pair{Any, Vector{Any}}[]
    for a in block.args
        a isa Expr || continue
        head, rest = a.head === :tuple ? (a.args[1], a.args[2:end]) : (a, Any[])
        push!(out, head.args[2] => Any[head.args[3]; rest])
    end
    out
end
# The steps of a block as written: per step `(lhs, rhs)`, each side `(conformation, bound, free)`,
# a metabolite `(name, inhibitor?)`; a `::Tag` other than `::Inh` is dropped.
function _testhelper_written_steps(block)
    lab = _testhelper_labeled(block)
    vals(l) = reduce(vcat, (v for (k, v) in lab if k in l); init = Any[])
    inh = setdiff(vals((:regulators, :catalytic_inhibitors)), vals((:substrates, :products)))
    untag(x) = x isa Expr && x.head === :(::) && x.args[2] !== :Inh ? x.args[1] : x
    met(x) = (x = untag(x); x isa Symbol ? (x, x in inh) : (x.args[1], true))
    form(f) = f isa Symbol ? (f, Tuple{Symbol, Bool}[]) : (f.args[1], sort(met.(f.args[2:end])))
    side(x) = (x = untag(x); x isa Expr && x.args[1] === :+ ?
               (form(x.args[2])..., sort(met.(x.args[3:end]))) : (form(x)..., Tuple{Symbol, Bool}[]))
    lines = [untag(l) for l in only(vals((:steps, :catalytic_steps))).args if l isa Expr]
    [(side(s.args[2]), side(s.args[3])) for l in lines for s in (l.head === :tuple ? l.args : [l])]
end
_testhelper_side(sp, mets) = (ER.conformation(sp),
    sort([(ER.name(b), b isa ER.CompetitiveInhibitor) for b in ER.bound(sp)]),
    sort([(ER.name(x), x isa ER.CompetitiveInhibitor) for x in mets]))
# k{i}f/k{i}r (SS) or K{i} (RE) of written step i, `_T` for the inactive state: the package
# constant of the stored step with the same content, read in the written direction.
function _testhelper_written_constants(m, written, nt)
    mech = _testhelper_lift(m)
    groups, fes = ER.steps(mech), ER._free_enz_set(mech)
    out = Pair{Symbol, Any}[]
    for (i, (lhs, rhs)) in enumerate(written)
        g, s = only((g, s) for (g, grp) in enumerate(groups) for s in grp if
            Set([_testhelper_side(ER.from_species(s), ER.consumed(s)),
                 _testhelper_side(ER.to_species(s), ER.released(s))]) == Set([lhs, rhs]))
        rep = ER._group_rep(groups[g], fes)
        same = _testhelper_side(ER.from_species(s), ER.consumed(s)) == lhs
        tag = mech isa ER.AllostericMechanism ? ER.cat_allo_state(mech, g) : :None
        states = tag === :NonequalAI ? ((:A, ""), (:I, "_T")) :
                 ((tag in (:None, :EqualAI) ? tag : :A, ""),)
        for (st, sfx) in states
            if ER.is_equilibrium(rep)
                push!(out, Symbol("K", i, sfx) => ER.name(ER.Krapid(rep, st), mech))
            else
                f, r = ER.name(ER.Kfor(rep, st), mech), ER.name(ER.Krev(rep, st), mech)
                same || ((f, r) = (r, f))
                push!(out, Symbol("k", i, "f", sfx) => f, Symbol("k", i, "r", sfx) => r)
            end
        end
    end
    (; (k => nt[v] for (k, v) in out if haskey(nt, v))...)
end
# Ligand sets of the regulatory sites in written order: explicit `regulatory_site` blocks,
# then one site per remaining allosteric regulator (the DSL's own rule, dsl.jl:644).
function _testhelper_written_sites(block)
    lab = _testhelper_labeled(block)
    explicit = [Set(Symbol[v for (k, vs) in _testhelper_labeled(only(v)) if k === :ligands for v in vs])
                for (k, v) in lab if k isa Expr && k.args[1] === :regulatory_site]
    regs = [x.args[1] for (k, vs) in lab if k === :allosteric_regulators for x in vs]
    [explicit; [Set([r]) for r in regs if !any(s -> r in s, explicit)]]
end

# ── Walk-through: Segel Ordered Bi Bi ──
blocks = spec_blocks()
spec = only(filter(s -> s.name == "Segel Ordered Bi Bi", MECHANISM_TEST_SPECS))
m = spec.mechanism
mech = ER.Mechanism(m)
fmt(sp, mets) = string(ER.name(sp), isempty(mets) ? "" : " + " * join(ER.name.(mets), " + "))
println("== Segel Ordered Bi Bi ==")
println("stored steps (group order, stored direction) and their package constants:")
for (g, grp) in enumerate(ER.steps(mech)), s in grp
    println("  group $g: ", fmt(ER.from_species(s), ER.consumed(s)), " -> ",
            fmt(ER.to_species(s), ER.released(s)), "   Kfor = ", ER.name(ER.Kfor(s, :None), mech),
            ", Krev = ", ER.name(ER.Krev(s, :None), mech))
end
println("source_steps (as the DSL returns them, after the Step constructor):")
for (i, s) in enumerate(Iterators.flatten(spec.source_steps))
    println("  $i: ", fmt(ER.from_species(s), ER.consumed(s)), " -> ", fmt(ER.to_species(s), ER.released(s)))
end
println("flat_idx (stored group -> written index): ", _testhelper_positional_flat_idx(mech, spec.source_steps))
dep, _ = ER._dependent_param_exprs(typeof(m))
println("fitted: ", ER.fitted_params(m))
println("dependents: ", [k => e for (k, e) in dep])
allnames = (ER.fitted_params(m)..., first.(collect(dep))..., :Keq, :E_total)
ntsym = NamedTuple{allnames}(allnames)                    # each value is its own name
pos = _testhelper_positional_params(m, ntsym; source_steps = spec.source_steps)
today = _testhelper_analytical_oracle_params(m, ntsym; source_steps = spec.source_steps)
role = _testhelper_written_constants(m, _testhelper_written_steps(blocks[spec.name]), ntsym)
for k in sort(collect(keys(role)); by = string)
    println("  ", rpad(k, 5), " before swap: ", rpad(pos[k], 16), " after swap: ", rpad(today[k], 16),
            " by role: ", role[k])
end

# ── Every spec with an analytical oracle: same values, and the oracle still matches ──
println("\n== all specs with an analytical oracle ==")
poskey(k) = occursin(r"^(k\d+[fr]|K\d+)(_T)?$", String(k))
nspec = 0; nsame = 0; npass = 0; nsites = 0; nallo = 0
for spec in MECHANISM_TEST_SPECS
    spec.analytical_rate_fn === nothing && continue
    global nspec += 1
    m = spec.mechanism
    blk = blocks[spec.name]
    written = _testhelper_written_steps(blk)
    nsrc = sum(length, spec.source_steps)
    rng = Random.MersenneTwister(1001)
    same_all = true; pass_all = true
    for _ in 1:5
        new_params, concs, all_params =
            _testhelper_random_independent_params_concs(m, spec.metabolite_names; rng = rng)
        Et = 0.1 + 9.9 * rand(rng)
        today = _testhelper_analytical_oracle_params(m, all_params; source_steps = spec.source_steps,
                                                     source_reg_sites = spec.source_reg_sites)
        role = _testhelper_written_constants(m, written, all_params)
        tpos = Dict(k => v for (k, v) in pairs(today) if poskey(k))
        same_all &= tpos == Dict(pairs(role))
        rest = (; (k => v for (k, v) in pairs(today) if !poskey(k))...)
        p = merge(rest, role, (Et = Et,))
        v_pkg = rate_equation(m, concs, merge(new_params, (E_total = Et,)))
        pass_all &= isapprox(v_pkg, spec.analytical_rate_fn(p, concs); rtol = 1e-10)
    end
    sites_ok = spec.source_reg_sites === nothing ? missing :
        _testhelper_written_sites(blk) ==
        [Set(ER.name.(ER.ligands(s))) for s in spec.source_reg_sites]
    global nsame += same_all; global npass += pass_all
    sites_ok === missing || (global nallo += 1; global nsites += sites_ok)
    println(rpad(spec.name, 66), " written steps ", lpad(length(written), 2), " (source ", lpad(nsrc, 2),
            ")  same k-values: ", same_all, "  oracle matches: ", pass_all,
            sites_ok === missing ? "" : "  site order: $(sites_ok)")
end
println("specs: $nspec; identical positional values: $nsame; oracle matches by role: $npass; ",
        "allosteric site order read from the block: $nsites of $nallo")
