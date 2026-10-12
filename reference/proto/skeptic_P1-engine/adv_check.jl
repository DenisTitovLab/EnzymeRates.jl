# ABOUTME: Runs the probe's Engine and EngineLean on the adversarial package-built mechanisms and checks every
# ABOUTME: claim: exact package agreement, oracle, detailed balance, speed-free predicate, dfree, Haldane division, rank.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf, Serialization
setprecision(BigFloat, Oracle.PREC)

pop = deserialize(joinpath(@__DIR__, "data", length(ARGS) >= 1 ? ARGS[1] : "adv_pkg.jls"))
rows = Dict{Symbol, Any}[]
for (n, m) in enumerate(pop)
    r = Dict{Symbol, Any}(:id => m.id, :src => m.source, :ok => m.ok, :err => m.stage_errors)
    push!(rows, r)
    println("\n=== ", m.id, "  ", m.source, "  pkg ok=", m.ok,
            isempty(m.stage_errors) ? "" : "  pkg stage errors: " * join(("$(a): $(first(b, 160))" for (a, b) in m.stage_errors), " | "))
    for (j, s) in enumerate(m.steps)
        println("    step $j: ", Pop.step_text(m, s), "   ", m.step_names[j])
    end
    println("    species: ", [s.name for s in m.species], "  fitted: ", m.fitted_params)
    L = try
        Engine.bare_rate_law(m)
    catch e
        r[:engine_err] = sprint(showerror, e); println("    Engine ERROR: ", first(r[:engine_err], 300)); nothing
    end
    lean = try
        EngineLean.bare_rate_law(m)
    catch e
        r[:lean_err] = sprint(showerror, e); println("    EngineLean ERROR: ", first(r[:lean_err], 300)); nothing
    end
    L === nothing && lean === nothing && continue
    if L !== nothing && lean !== nothing
        V, sv, ln, ld, lf = lean
        r[:lean_same] = (Sk.sym(ln, V), Sk.sym(ld, V), Sk.sym(lf, V)) ==
                        (Sk.sym(L.num, L.vars), Sk.sym(L.den, L.vars), Sk.sym(L.dfree, L.vars))
        println("    lean == full engine: ", r[:lean_same])
    end
    L === nothing && continue
    st = L.stats
    println("    segments ", st[:segments], " blocks ", st[:block_sizes], " cycles ", st[:cycles],
            " max|n| ", st[:max_n], " loops ", st[:loops], " re_closing ", st[:re_closing],
            " den ", st[:den_terms], " num ", st[:num_terms])
    # (1) exact agreement with the package's raw polynomials
    if !isempty(m.raw_den)
        num, den, dfree = Engine.grouped(L, m)
        r[:num_eq] = Sk.sym(num, L.vars) == Sk.topoly(m.raw_num)
        r[:den_eq] = Sk.sym(den, L.vars) == Sk.topoly(m.raw_den)
        r[:dfree_eq] = Sk.sym(dfree, L.vars) == Sk.topoly(m.raw_dfree)
        println("    package raw polys identical: num ", r[:num_eq], " den ", r[:den_eq], " dfree ", r[:dfree_eq])
    end
    # (2) oracle at bare constants; detailed balance at consistent constants
    rng = Random.Xoshiro(4242 + n)
    errs = Float64[]
    for _ in 1:3
        kb = Sk.bare_constants(m, rng); c = Sk.random_concs(m, rng)
        push!(errs, Sk.relerr(Sk.lawv(L, Sk.vals_of(L, m, kb, c)), Oracle.steady_state_v(m, kb, c)))
    end
    r[:oracle] = maximum(errs)
    Keq = Oracle.logu(rng)
    kc, _, _ = Oracle.consistent_constants(m, rng; Keq)
    c = Oracle.concs_at_ratio(m, rng, Keq)
    neq = Poly.evaluate(L.num, Sk.vals_of(L, m, kc, c)); c[m.products[1]] /= 10
    r[:db] = Float64(abs(neq) / abs(Poly.evaluate(L.num, Sk.vals_of(L, m, kc, c))))
    @printf("    oracle max rel %.1e   detailed balance %.1e\n", r[:oracle], r[:db])
    # (3) speed-free: engine numeric, oracle numeric, probe predicate, exact predicates
    kb = Sk.bare_constants(m, rng); c2 = Sk.random_concs(m, rng)
    kc2, _, _ = Oracle.consistent_constants(m, rng; Keq = Oracle.logu(rng))
    r[:sf_thermo] = Sk.speed_free_numeric(L, m, kc2, c2)
    r[:sf_bare] = Sk.speed_free_numeric(L, m, kb, c2)
    r[:sf_thermo_oracle] = Sk.speed_free_oracle(m, kc2, c2)
    r[:sf_bare_oracle] = Sk.speed_free_oracle(m, kb, c2)
    r[:pred_probe] = st[:dead_end]
    r[:pred_thermo], r[:pred_bare] = Sk.speed_free_exact(m)
    println("    speed-free thermo: engine ", r[:sf_thermo], " oracle ", r[:sf_thermo_oracle],
            " | probe predicate ", r[:pred_probe], " exact predicate ", r[:pred_thermo])
    println("    speed-free bare:   engine ", r[:sf_bare], " oracle ", r[:sf_bare_oracle],
            " | exact bare predicate ", r[:pred_bare])
    # (4) dfree vs the oracle's free-form population
    i, root = Sk.free_and_root(m)
    r[:dfree_pop] = Sk.dfree_vs_population(L, m, kb, c2)
    r[:free_is_root] = i == root
    println("    free form ", i === nothing ? "-" : m.species[i].name, " segment root ",
            root === nothing ? "-" : m.species[root].name, "; free forms ",
            [s.name for s in m.species if Pop.is_free(s)],
            @sprintf("; dfree/den vs oracle x_free rel %.1e", r[:dfree_pop]))
    # (5) Reduced: package dependents, non-integer exponents, (1 − Γ/Keq) division, group phantoms
    m.ok && !isempty(m.pkg_den) || continue
    nonint = [k for (k, mono) in m.dependent_monomials if any(!isinteger(e) for (_, e) in mono)]
    r[:nonint] = nonint
    println("    dependents: ", m.dependents)
    isempty(nonint) || println("    NON-INTEGER dependent exponents: ", nonint)
    if isempty(nonint)
        V = Poly.Vars()
        fromsym(pp) = Sk.PR(Dict{Poly.Mono, Rational{Int}}(
            Poly.mono((Poly.var!(V, s) => e for (s, e) in mono)...) => c for (mono, c) in pp))
        num = fromsym(m.pkg_num)
        dep = Dict(Poly.var!(V, k) => Poly.mono((Poly.var!(V, s) => Int(e) for (s, e) in mono)...)
                   for (k, mono) in m.dependent_monomials)
        num = Poly.subst(num, v -> get(dep, v, nothing))
        kv = Poly.var!(V, :Keq)
        x = Poly.mono((Poly.var!(V, p) => 1 for p in m.products)...,
                      (Poly.var!(V, s) => -1 for s in m.substrates)..., kv => -1)
        q = EngineLean.divide_one_minus(num, x, kv)
        r[:divisible] = q !== nothing
        r[:q_positive] = q !== nothing && all(>(0), values(q.t))
        println("    (1 − Γ/Keq) divides the reduced numerator: ", r[:divisible],
                q === nothing ? "" : "; quotient all positive: $(r[:q_positive]); terms $(Poly.nterms(num)) → $(Poly.nterms(q))")
    end
    # identifiable rank of the package's Reduced law in its fitted parameters
    try
        nθ, rk = Sk.identifiable_rank(m, Random.Xoshiro(99 + n))
        r[:rank] = (nθ, rk)
        println("    fitted params ", nθ, ", identifiable rank ", rk)
    catch e
        println("    rank ERROR: ", first(sprint(showerror, e), 200))
    end
end
serialize(joinpath(@__DIR__, "data", "adv_check_rows.jls"), rows)
