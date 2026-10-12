# ABOUTME: Runs the pipeline (derive) on every record of the given sets: verdicts, stage times,
# ABOUTME: sizes, keys, the law against its oracle; usage: [RANGE=a:b] run_pipeline.jl label sets…
include(joinpath(@__DIR__, "lib.jl"))

const LABEL, SETS = ARGS[1], ARGS[2:end]
const NPOINTS = parse(Int, get(ENV, "NPOINTS", "2"))
const RANGE = eval(Meta.parse(get(ENV, "RANGE", "1:typemax(Int)")))

"""The pipeline's artifacts and measurements for record `r` (a NamedTuple row)."""
function measure(r, rng)
    m = r.c.mech
    base = (set = r.set, index = r.c.index, id = r.p.id, allo = S.is_allosteric(m),
            nforms = length(m.forms), nedges = length(m.edges), ngroups = S.ngroups(m))
    t_derive = @elapsed d = try
        S.derive(m)
    catch err
        err isa ErrorException || rethrow()
        return merge(base, (ok = false, msg = err.msg))
    end
    sol, L = d.solution, d.law
    # each stage timed again on its own; for init_terter on every tenth record only
    timed = !startswith(r.set, "init_terter") || r.c.index % 10 == 0
    t_solve, t_law, t_key, t_print, t_code, t_kcat = fill(NaN, 6)
    if timed
        t_solve = @elapsed S.solve(m)
        t_law = @elapsed S.rate_law(m, sol)
        t_key = @elapsed S.canonical(L)
        t_print = @elapsed (S.rate_equation_string(m, L, sol); S.coefficient_form(m, L))
        t_code = @elapsed S.body(S.program(m, sol, L))
        t_kcat = @elapsed try S.kcat_candidates(m, sol) catch end
    end
    states = L isa S.Law ? [L] : filter(!isnothing, [L.active, L.inactive])
    cols, A, _ = S.constraint_system(m)
    masks = L isa S.Law ? [trues(length(m.edges))] :
            [trues(length(m.edges)), S.inactive_mask(m)]
    worst = 0.0
    t_oracle = @elapsed for _ in 1:NPOINTS
        x, cv, conc = random_point(r, sol, rng)
        worst = max(worst, relerr(law_value(L, conc, x), oracle_value(r, sol, x, cv, conc)))
    end
    merge(base, (ok = true, msg = "", key = d.key.key, order = d.key.order,
                 nfit = length(S.fitted_params(m)), rref = size(A), ncols = length(cols),
                 trees = [any(k) ? tree_count(m, k) : 0 for k in masks],
                 complexity = sum(S.complexity(m, k) for k in masks if any(k)),
                 nQ = sum(S.nterms(s.Q) for s in states), nD = sum(S.nterms(s.D) for s in states),
                 rules = [first(n) for n in d.normalizers], ops = countops(d.code),
                 printed = length(d.printed), coefficients = length(d.coefficients),
                 kcat = d.kcat isa String ? d.kcat : length(d.kcat[2]),
                 halves = any(!isinteger(e) for s in states for P in (s.Q, s.D)
                              for mo in keys(P.terms) for (_, e) in mo),
                 big = !(states[1] isa S.Law{Int}), oracle = worst, t_derive, t_solve, t_law,
                 t_key, t_print, t_code, t_kcat, t_oracle))
end

rng = Xoshiro(hash(LABEL) % UInt32)
rows = NamedTuple[]
t0 = time()
first_set = records_of(SETS[1])
measure(first_set[1], rng)                       # compiles the pipeline before timing
for set in SETS
    recs = set == SETS[1] ? first_set : records_of(set)
    for r in recs[intersect(RANGE, eachindex(recs))]
        push!(rows, measure(r, rng))
    end
    @printf("%-20s done at %.0f s\n", set, time() - t0)
    flush(stdout)
end
serialize(joinpath(IDATA, "pipe_$LABEL.jls"), rows)
ok = filter(x -> x.ok, rows)
@printf("%s: %d records, %d derived, %d rejected; oracle worst %.1e (above 1e-25: %d); \
        derive %s ms; %.0f s, maxrss %.0f MB\n", LABEL, length(rows), length(ok),
        length(rows) - length(ok), maximum(x -> x.oracle, ok; init = 0.0),
        count(x -> x.oracle > 1e-25, ok), q3([1e3 * x.t_derive for x in ok]), time() - t0,
        Sys.maxrss() / 2^20)
