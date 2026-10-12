# ABOUTME: Allosteric printing check: both printed forms of every admitted MWC law evaluated against
# ABOUTME: the law, per-conformation units, artifacts, and examples written to logs/show_allo.log.
include(joinpath(@__DIR__, "..", "engine", "lib.jl"))
include(joinpath(@__DIR__, "value.jl"))

"""`v` of a printed text: its lines evaluated, each once every name it reads is defined."""
function eval_lines(text, pars, conc)
    mod = Module()
    Core.eval(mod, :(params = $pars))
    Core.eval(mod, :(concs = $conc))
    foreach(x -> Core.eval(mod, :($x = $(getfield(pars, x)))), keys(pars))
    foreach(x -> Core.eval(mod, :($x = $(getfield(conc, x)))), keys(conc))
    todo = [Meta.parse(l) for l in split(text, "\n") if !startswith(l, "#")]
    while !isempty(todo)
        left = filter(todo) do ex
            try
                Core.eval(mod, ex); false
            catch err
                err isa UndefVarError || rethrow(); true
            end
        end
        length(left) == length(todo) && error("printed law: undefined names")
        todo = left
    end
    Core.eval(mod, :v)
end

"""Whether every term of conformation law `L`'s printed denominator, and of its flip form's
weight over the same normalizer, is dimensionless, and every numerator term times ΠS a rate."""
function units_ok(m, L, sol)
    f = L.trees[1] * L.shift
    u = var_units(m, L, sol)
    rule, lambda = S.normalizer(m, L)
    rule === :coefficient && return true
    Dc = S.split_conc(L.D, L.vars)
    punits(p) = p === nothing ? (0 // 1, 0 // 1) : mono_units(first(keys(p.terms)), u)
    apart = reduce((a, b) -> a .+ b, (e .* punits(S.content(Dc[mo])[3]) for (mo, e) in lambda);
                   init = (0 // 1, 0 // 1))
    tu(cm, x) = mono_units(S.mmul(cm, x[2]), u) .+ punits(x[3]) .- (x[4] ? apart : (0, 0))
    P = S.printed(m, L)
    F = S.printed(m, S.Law(L.vars, L.H, f, L.D, L.sigma, L.trees, L.shift)).q
    su = mono_units(S.conc_mono(L.vars, m.substrates), u)
    all(tu(g...) == (0, 0) for g in [P.d; F]) && all(tu(g...) .+ su == (0, -1) for g in P.q)
end

const SHOW = ("specs/27", "specs/35", "specs/40", "handwritten/30", "allo_uniuni/12",
              "allo_bibi/1460", "AD3d_copy_square_SS_split")
records = [r for r in load_records() if S.is_allosteric(r.c.mech) && S.admit(r.c.mech) === nothing]
rng = Xoshiro(20261015)
worst, units_bad, artifacts, wrong, n = big(0), String[], String[], String[], 0
show = open(joinpath(@__DIR__, "logs", "show_allo.log"), "w")
for r in records
    m = r.c.mech
    sol = S.solve(m)
    A = S.allo_law(m, sol)
    full, coef = S.rate_equation_string(m, A, sol), S.coefficient_form(m, A)
    pars = random_params(m, rng)
    cd = Oracle.concs_at_ratio(r.p, rng, pars.Keq * 7)
    conc = (; (x => cd[x] for x in S.metabolites(m))...)
    v = pars.E_total * allo_value(A, cd, pars)
    for t in (full, coef)
        e = abs(eval_lines(t, pars, conc) - v) / abs(v)
        e > 1e-100 && push!(wrong, r.p.id)
        global worst = max(worst, e)
        occursin(r"\^1\b|\b1 \^|\* 1\b|\b1 \*|L \* 1\b|\(1\)", t) && push!(artifacts, r.p.id)
    end
    for L in (A.active, A.inactive)
        L === nothing || units_ok(m, L, sol) || push!(units_bad, r.p.id)
    end
    global n += 1
    if r.p.id in SHOW
        println(show, "== ", r.set, " ", r.p.id, " ", r.p.source, " (n = ", m.multiplicity,
                ", tags ", m.tags, ")\n", full, "\n-- coefficient form:\n", coef, "\n")
    end
end
close(show)
@printf("admitted MWC laws printed: %d; both forms evaluate to the law, worst relative %.1e\n",
        n, worst)
println("printed value off by more than 1e-100: ", length(unique(wrong)), " ", first(unique(wrong), 10))
println("units: every printed denominator term and flip weight dimensionless, numerator ",
        "terms rates, except in ", length(units_bad), " ", first(unique(units_bad), 10))
println("artifacts (^1, 1 ^, * 1, 1 *, L * 1, (1)): ", length(unique(artifacts)), " ",
        first(unique(artifacts), 10))
