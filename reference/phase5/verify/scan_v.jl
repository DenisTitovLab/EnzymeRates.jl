# ABOUTME: Verifier GC workload for one law-type variant (ARGS[1] = U or C): caller loops over a
# ABOUTME: population reading law fields, printing, keys, codegen, plus adversarial Poly-merging loops.
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
const VAR = ARGS[1]
include(joinpath(SP, "phase45/research5/law/variants", VAR, "Slice$VAR.jl"))
const S = getfield(Main, Symbol("Slice", VAR))
include(joinpath(VD, "sliceside.jl"))
using Serialization
t0 = time()
pop = S.Mechanism[]
for l in eachline(joinpath(VD, "data", "scan_pop.tsv"))
    try
        push!(pop, mkmech(S, eval(Meta.parse(l))))
    catch
    end
end
println("julia ", VERSION, " variant ", VAR, ": ", length(pop), " mechanisms, ",
        count(S.is_allosteric, pop), " allosteric; rate_law returns ",
        Base.return_types(S.rate_law, (S.Mechanism, S.Solution))); flush(stdout)
const BYTES = let io = IOBuffer(); serialize(io, pop); take!(io) end

# the crash probe's caller pattern: a loop reading a field of what rate_law returns
function derive(bytes)
    tot = zeros(Int, 6)
    for m in deserialize(IOBuffer(bytes))
        sol = S.solve(m)
        L = S.rate_law(m, sol)
        nD = L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D)
        ns = 0
        if L isa S.AlloLaw
            for (WA, WI, k) in L.sites
                ns += S.nterms(WA) + S.nterms(WI) + k
            end
            L.inactive === nothing || (ns += S.nterms(L.inactive.D))
        end
        s = S.rate_equation_string(m, L, sol)
        c = S.coefficient_form(m, L)
        k = S.canonical(L)
        b = S.body(S.program(m, sol, L))
        tot .+= (1, nD, ns, length(s) + length(c), length(k.order), length(string(b)))
    end
    tot
end

# the same loop over a typed population (m::Mechanism, so rate_law's return is a small Union)
function derive_typed(ms::Vector{S.Mechanism})
    tot = zeros(Int, 3)
    for m in ms
        sol = S.solve(m)
        L = S.rate_law(m, sol)
        nD = L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D)
        ns = 0
        if L isa S.AlloLaw
            for (WA, WI, k) in L.sites
                ns += S.nterms(WA) + S.nterms(WI) + k
            end
            L.inactive === nothing || (ns += S.nterms(L.inactive.D))
        end
        tot .+= (1, nD, ns)
    end
    tot
end
const LT = Union{Base.return_types(S.rate_law, (S.Mechanism, S.Solution))...}

# adversarial: merge a Poly with another type inside a loop
function merge_vec(ls)
    t = 0
    for L in ls
        y = L isa S.Law ? L.H : L.sites          # Union{Poly, Vector}
        t += y isa S.Poly ? S.nterms(y) : length(y)
    end
    t
end
function merge_nothing(ls)
    t = 0
    for (i, L) in enumerate(ls)
        z = isodd(i) ? L.H : nothing             # Union{Nothing, Poly}
        t += z === nothing ? 0 : S.nterms(z)
    end
    t
end
function merge_same(ls)
    t = 0
    for L in ls
        x = L isa S.Law ? L.D : L.fA             # Poly in both branches
        t += S.nterms(x)
    end
    t
end

# identify-like: filter, key, names; then compiled rate_equation in a loss loop
function pass1(ms)
    out = Any[]
    for m in ms
        Int(S.complexity(m)) > 337 && continue
        sol = S.solve(m)
        L = S.rate_law(m, sol)
        push!(out, (S.canonical(L).key, length(S.fitted_params(m))))
    end
    out
end
function loss(em, rows, p)
    s = 0.0
    for c in rows
        s += abs(S.rate_equation(em, c, p) - 1.0)
    end
    s
end

ref = nothing
for rd in 1:(get(ENV, "LOC", "0") == "1" ? 0 : 3)
    tot = derive(BYTES)
    tot2 = derive_typed(convert(Vector{S.Mechanism}, deserialize(IOBuffer(BYTES))))
    ls = LT[S.rate_law(m) for m in pop]
    a, b, c = merge_vec(ls), merge_nothing(ls), merge_same(ls)
    keys1 = pass1(pop)
    lsum = 0.0
    for m in pop[1:10]
        em = S.compile(m)
        names = S.parameters(em)
        p = NamedTuple{names}(Tuple(1.0 + 0.1 * i for i in 1:length(names)))
        mets = S.metabolites(em)
        rows = [NamedTuple{mets}(Tuple(0.5 + 0.1 * j + 0.01 * i for j in 1:length(mets))) for i in 1:50]
        lsum += loss(em, rows, p)
    end
    sig = (tot, a, b, c, length(keys1), hash(keys1), round(lsum; digits = 6), tot2)
    global ref = ref === nothing ? sig : ref
    GC.gc(true)
    println("round $rd: ", sig == ref ? "same" : "DIFFERENT", " ", sig[1], " ", sig[2:4], " keys ", sig[5]); flush(stdout)
end
println("done in ", round(time() - t0; digits = 1), " s")
