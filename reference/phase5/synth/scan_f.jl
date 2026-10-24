# ABOUTME: Typed GC workload for one law-type variant (ARGS[1] = C, E or F): the verifier's loops (scan_v.jl)
# ABOUTME: plus loops that merge a freshly derived law with another type, each also run with full GCs inside.
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
const VAR = ARGS[1]
const VPATH = Dict("C" => joinpath(SP, "phase45/research5/law/variants/C/SliceC.jl"),
                   "U" => joinpath(SP, "phase45/research5/law/variants/U/SliceU.jl"),
                   "E" => joinpath(VD, "variantE/SliceE.jl"),
                   "F" => joinpath(SP, "phase45/research5/synth/variantF/SliceF.jl"))
include(VPATH[VAR])
const S = getfield(Main, Symbol("Slice", VAR))
include(joinpath(VD, "sliceside.jl"))
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
        Base.return_types(S.rate_law, (S.Mechanism, S.Solution)), "; ismutabletype Poly/Law/AlloLaw: ",
        ismutabletype(S.Poly{Int}), "/", ismutabletype(S.Law{Int}), "/", ismutabletype(S.AlloLaw{Int}))
flush(stdout)

nI(L) = hasproperty(L, :inactive) && L.inactive !== nothing ? S.nterms(L.inactive.D) : 0

# the verifier's typed derivation loop (scan_v.jl derive_typed), inactive read through nI
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
            ns += nI(L)
        end
        tot .+= (1, nD, ns)
    end
    tot
end
const LT = Union{Base.return_types(S.rate_law, (S.Mechanism, S.Solution))...}

# the verifier's adversarial loops over stored laws (scan_v.jl)
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

# loops that merge a law derived in the loop with another type; `full` adds a full GC per step
function fresh_or_nothing(ms::Vector{S.Mechanism}, full::Bool)
    t = 0
    for (i, m) in enumerate(ms)
        z = isodd(i) ? S.rate_law(m) : nothing  # Union{Nothing, Law, AlloLaw}
        t += z === nothing ? 0 : z isa S.Law ? S.nterms(z.D) : S.nterms(z.active.D)
        full && GC.gc(true)
    end
    t
end
function fresh_or_vec(ms::Vector{S.Mechanism}, full::Bool)
    t = 0
    for m in ms
        L = S.rate_law(m)
        y = L isa S.Law ? L : L.sites            # Union{Law, Vector}
        t += y isa S.Law ? S.nterms(y.Q) : length(y)
        full && GC.gc(true)
    end
    t
end
function fresh_active_or_poly(ms::Vector{S.Mechanism}, full::Bool)
    t = 0
    for m in ms
        L = S.rate_law(m)
        y = L isa S.AlloLaw ? L.active : L.H     # Union{Law, Poly}
        t += y isa S.Law ? S.nterms(y.D) : S.nterms(y)
        full && GC.gc(true)
    end
    t
end
function fresh_poly_or_vec(ms::Vector{S.Mechanism}, full::Bool)
    t = 0
    for m in ms
        L = S.rate_law(m)
        y = L isa S.Law ? L.H : L.sites          # crash_c.jl's loop
        t += y isa S.Poly ? S.nterms(y) : length(y)
        full && GC.gc(true)
    end
    t
end

ref = nothing
for rd in 1:3
    tot2 = derive_typed(pop)
    ls = LT[S.rate_law(m) for m in pop]
    a, b = merge_vec(ls), merge_nothing(ls)
    f = (fresh_or_nothing(pop, false), fresh_or_vec(pop, false),
         fresh_active_or_poly(pop, false), fresh_poly_or_vec(pop, false))
    g = (fresh_or_nothing(pop[1:60], true), fresh_or_vec(pop[1:60], true),
         fresh_active_or_poly(pop[1:60], true), fresh_poly_or_vec(pop[1:60], true))
    sig = (tot2, a, b, f, g)
    global ref = ref === nothing ? sig : ref
    GC.gc(true)
    println("round $rd: ", sig == ref ? "same" : "DIFFERENT", " ", sig, ", full gc ok"); flush(stdout)
end
println("done in ", round(time() - t0; digits = 1), " s")
