# ABOUTME: Q11 on the package (graph_mechanism_wt): a structured joint-solve output (Keq exponent and
# ABOUTME: Rational{Int} factor exponents per dependent) rebuilds today's Exprs; plus a Q10 Bareiss stress test.
using EnzymeRates, Random, Printf, LinearAlgebra, Test
const ER = EnzymeRates
const WT = dirname(dirname(pathof(EnzymeRates)))
println("julia ", VERSION, "; EnzymeRates from ", WT)

"""The 42 `MECHANISM_TEST_SPECS`, evaluated as census/corpus.jl's `load_specs` does."""
function load_specs()
    path = joinpath(WT, "test", "mechanism_definitions_for_test_enzyme_derivation.jl")
    mod = Module(:Q11Specs)
    Core.eval(mod, :(using EnzymeRates, LinearAlgebra, Random, Test))
    Core.eval(mod, quote
        macro _testhelper_enzyme_mechanism_src(block)
            esc(Expr(:tuple, Expr(:macrocall, Symbol("@enzyme_mechanism"), __source__, block), nothing))
        end
        macro _testhelper_allosteric_mechanism_src(block)
            esc(Expr(:tuple, Expr(:macrocall, Symbol("@allosteric_mechanism"), __source__, block),
                     nothing, nothing))
        end
        macro _testhelper_with_raw(call)
            esc(Expr(:tuple, call, nothing))
        end
    end)
    srcmacro(ex) = ex isa Expr && (ex.head === :macro && ex.args[1].args[1] in
        (:_testhelper_enzyme_mechanism_src, :_testhelper_allosteric_mechanism_src,
         :_testhelper_with_raw) ||
        ex.head === :macrocall && ex.args[1] isa GlobalRef &&
        ex.args[1].name === Symbol("@doc") && srcmacro(ex.args[end]))
    for ex in Meta.parseall(read(path, String); filename = path).args
        ex isa LineNumberNode && continue
        srcmacro(ex) && continue
        try Core.eval(mod, ex) catch end
    end
    [(s.name, s.mechanism) for s in Core.eval(mod, :MECHANISM_TEST_SPECS)]
end

"""Structured solve output: per dependent name, `(keq, [independent => exponent])` with
`Rational{Int}` exponents, from the same joint solve as `_dependent_param_exprs`; also
whether the row reduction ran in `Rational{BigInt}`."""
function structured(m)
    sys = ER._joint_system(m)
    M, pivots = ER._reduce_joint_system(sys.A, sys.rhs)
    n = length(sys.columns)
    names = [ER.name(ER._column_parameter(m, sys.flat, c), m) for c in sys.columns]
    class(c) = ER._class(c, sys.flat[c.step][2], sys.tags[sys.flat[c.step][2]])
    last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
    kept_in_law(i) = (c = sys.columns[i];
                      (!c.inactive || sys.law[c.step]) && last_of[class(c)] == i)
    free = setdiff(1:n, pivots)
    deps = Dict{Symbol, Tuple{Rational{Int}, Vector{Pair{Symbol, Rational{Int}}}}}()
    for (r, p) in enumerate(pivots)
        kept_in_law(p) || continue
        deps[names[p]] = (Rational{Int}(M[r, n + 1]),
                          [names[f] => Rational{Int}(-M[r, f]) for f in free if !iszero(M[r, f])])
    end
    deps, eltype(M) != Rational{Int}
end

value(dep, params) = prod((params[x]^Float64(e) for (x, e) in dep[2]); init = params.Keq^Float64(dep[1]))

# specs come compiled; the solve reads the concrete mechanism, as the Type method does
mechs = Tuple{String, Any}[("spec:" * n, m isa ER.AbstractEnzymeMechanism ? ER._concrete(m) : m)
                           for (n, m) in load_specs()]
for (label, text) in (
    ("bibi", "@enzyme_reaction begin\n substrates: A[C], B[N]\n products: P[C], Q[N]\n competitive_inhibitors: A, I\nend"),
    ("pingpong", "@enzyme_reaction begin\n substrates: A[CX], B[N]\n products: P[C], Q[NX]\n allosteric_regulators: R\n oligomeric_state: 2\n shared_catalytic_site: (A, P)\nend"),
    ("terbi", "@enzyme_reaction begin\n substrates: A[C], B[N], D[X]\n products: P[CN], Q[X]\nend"))
    rxn = Core.eval(Main, Meta.parse(text))
    init = ER.init_mechanisms(rxn)
    label == "terbi" && (init = init[1:20:end])
    append!(mechs, [(label * ":init", m) for m in init])
    label == "terbi" && continue
    for p in init[1:min(end, 25)]
        append!(mechs, [(label * ":child", c) for c in ER.expand_mechanisms([p], rxn)])
    end
end
println("mechanisms: ", length(mechs), " (", count(x -> x[2] isa ER.AllostericMechanism, mechs), " allosteric)")
same, nbig, ndep, worst, errs, nexp, nhalf = 0, 0, 0, 0.0, 0, 0, 0
rng = Xoshiro(11)
for (label, m) in mechs
    try
        dep, indep = ER._dependent_param_exprs(m)
        st, big = structured(m)
        global nbig += big
        rebuilt = Dict(k => ER.build_power_expr(v[1], [(x, e) for (x, e) in v[2]]) for (k, v) in st)
        global same += rebuilt == dep
        global ndep += length(st)
        global nhalf += count(v -> any(e -> denominator(e) != 1, last.(v[2])) || denominator(v[1]) != 1, values(st))
        names = unique([collect(indep); [x for v in values(st) for (x, _) in v[2]]])
        params = (; (x => exp(randn(rng)) for x in names if x !== :L)..., Keq = exp(randn(rng)))
        for (k, v) in st
            a, b = value(v, params), ER._evaluate_dependent(dep[k], params)
            global worst = max(worst, abs(a - b) / abs(b))
        end
    catch e
        global errs += 1
        errs <= 3 && println("  error on ", label, ": ", first(sprint(showerror, e), 160))
    end
end
@printf("Q11: structured output rebuilds _dependent_param_exprs byte-identically (==) in %d of %d mechanisms; %d dependents (%d with a non-integer exponent); row reduction in Rational{BigInt} in %d; worst relative value difference %.1e; errors %d\n",
        same, length(mechs), ndep, nhalf, nbig, worst, errs)

# ── Q10 stress: checked-Int Bareiss on random multigraph Laplacians ─────────────────────
"""Determinant of `A` by fraction-free (Bareiss) elimination in checked `Int`
(checks.jl's `det_bareiss`); an `OverflowError` when a product leaves `Int`'s range."""
function det_bareiss(A::Matrix{Int})
    M, n = copy(A), size(A, 1)
    n == 0 && return 1
    sgn, prev = 1, 1
    for k in 1:(n - 1)
        if M[k, k] == 0
            p = findfirst(i -> M[i, k] != 0, (k + 1):n)
            p === nothing && return 0
            p += k
            for j in 1:n
                M[k, j], M[p, j] = M[p, j], M[k, j]
            end
            sgn = -sgn
        end
        for i in (k + 1):n, j in (k + 1):n
            M[i, j] = div(Base.checked_sub(Base.checked_mul(M[i, j], M[k, k]),
                                           Base.checked_mul(M[i, k], M[k, j])), prev)
        end
        prev = M[k, k]
    end
    sgn * M[n, n]
end
rng = Xoshiro(7)
exact, wrong, over, small_over = 0, 0, 0, typemax(Int128)
for n in 2:14, trial in 1:30
    μ = rand(rng, 1:3)
    A = zeros(Int, n, n)
    for a in 1:n, b in (a + 1):n
        rand(rng) < 0.5 || continue
        k = rand(rng, 1:μ)
        A[a, a] += k; A[b, b] += k; A[a, b] -= k; A[b, a] -= k
    end
    Lp = A[2:end, 2:end]
    d = det(BigInt.(Lp))
    try
        det_bareiss(Lp) == d ? (global exact += 1) : (global wrong += 1)
    catch err
        err isa OverflowError || rethrow()
        global over += 1
        global small_over = min(small_over, Int128(min(d, big(typemax(Int128)))))
    end
end
println("Q10 stress (checks.jl's 390 random multigraph Laplacians, same seed): exact $exact, wrong $wrong, overflow $over (smallest tree count that overflowed: $small_over)")
