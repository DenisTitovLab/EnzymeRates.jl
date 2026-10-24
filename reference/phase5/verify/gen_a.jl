# ABOUTME: Verifier run A (package + stock slice): R5.5 round-2 counts rebuilt independently, the
# ABOUTME: phantom population (bi-bi, round 2 from parents outside the study's), and a scan population.
using EnzymeRates, Random, Printf
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
include(joinpath(VD, "pkgside.jl"))
include(joinpath(SP, "slice/src/Slice.jl"))
include(joinpath(VD, "sliceside.jl"))
const S = Slice
T0 = time()
stamp(msg) = (@printf("[%5.0f s, %4.0f MB] %s\n", time() - T0, Sys.maxrss() / 2^20, msg); flush(stdout))

# ── Part 1: R5.5, census bi-bi (competitive A and I), rounds 0-2, uniform ──────────────────
rxnI = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    competitive_inhibitors: A, I
end
function r55()
    seeds = unique!(collect(ER.init_mechanisms(rxnI)))
    base, _ = ER._base_tier(seeds, rxnI)
    seen = Set{UInt64}(hash.(base))
    pool0 = [m for m in base if cxA(m) <= 337]
    next1 = Any[]
    n1 = 0
    for p in pool0, c in ER.expand_mechanisms([p], rxnI)
        h = hash(c)
        h in seen && continue
        push!(seen, h)
        n1 += 1
        cxA(c) <= 337 && push!(next1, c)
    end
    stamp("r55 round 0: $(length(seeds)) seeds, $(length(base)) base; round 1: $n1 new, $(length(next1)) within 337")
    for (label, rng) in (("study seed Xoshiro(20261014)", Xoshiro(20261014)),
                         ("fresh seed Xoshiro(7)", Xoshiro(7)))
        parents = next1[sort(randperm(rng, length(next1))[1:1500])]
        sn = copy(seen)
        ncand = nadm = nallo = nover = nmax = 0
        cmpn = cmpbad = slicefail = 0
        overs = Int[]
        for p in parents, c in ER.expand_mechanisms([p], rxnI)
            h = hash(c)
            h in sn && continue
            push!(sn, h)
            ncand += 1
            a = cxA(c)
            a <= 337 || continue
            nadm += 1
            c isa ER.AllostericMechanism || continue
            nallo += 1
            i = cxI(c)
            a + i > 337 && (nover += 1; push!(overs, a + i))
            max(a, i) > 337 && (nmax += 1)
            # independent recount on the slice's graph: V×τ of all edges and of the inactive mask
            try
                sm = mkmech(S, rawrec(c))
                sa, si = Int(S.complexity(sm)), Int(S.complexity(sm, S.inactive_mask(sm)))
                cmpn += 1
                (sa, si) == (a, i) || (cmpbad += 1; cmpbad <= 5 &&
                    println("  slice/package V×τ mismatch: package ($a, $i) slice ($sa, $si)"))
            catch e
                slicefail += 1
                slicefail <= 3 && println("  slice failed: ", first(sprint(showerror, e), 150))
            end
        end
        stamp("r55 round 2 [$label]: 1500 parents of $(length(next1)), $ncand new, $nadm within 337 (active), " *
              "$nallo allosteric within; newly over (sum) $nover, newly over (max) $nmax; sums $(sort(overs)); " *
              "slice recount $cmpn compared, $cmpbad mismatches, $slicefail slice failures")
    end
end
r55()
GC.gc()

# ── Part 2: phantom population: bi-bi without inhibitors, identify's base tier ─────────────
rxn2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
mutable struct Pop
    mechs::Vector{Any}
    rows::Vector{Any}
    gid::Dict{Any, Int}
end
function register!(pop, m, round, par, mv)
    g = get(pop.gid, m, 0)
    g == 0 || return g
    g = length(pop.mechs) + 1
    pop.gid[m] = g
    push!(pop.mechs, m)
    allo = m isa ER.AllostericMechanism
    n = allo ? ER.catalytic_multiplicity(m) : 0
    cx, np, status = -1, -1, "ok"
    try
        cx = cxA(m)
        np = ER._independent_param_count(m)
        status = cx > 337 ? "cx_skip" : np > 20 ? "np_skip" : "ok"
    catch e
        status = "error"
    end
    push!(pop.rows, (; gid = g, round, parent = par, move = mv, allo, n, cx, np, status))
    g
end
okg(pop, r) = [x.gid for x in pop.rows if x.round == r && x.status == "ok"]
function expand!(pop, parents, round)
    for p in parents, (c, mv) in children_labeled(pop.mechs[p], rxn2)
        register!(pop, c, round, p, mv)
    end
end
pop = Pop(Any[], Any[], Dict{Any, Int}())
seeds2 = unique!(collect(ER.init_mechanisms(rxn2)))
base2, _ = ER._base_tier(seeds2, rxn2)
foreach(m -> register!(pop, m, 0, 0, m in seeds2 ? :seed : :seed_flip), base2)
expand!(pop, okg(pop, 0), 1)
pool1 = okg(pop, 1)
theirs = pool1[sort!(randperm(Xoshiro(60), length(pool1))[1:1000])]
rest = setdiff(pool1, theirs)
mine = rest[sort!(randperm(Xoshiro(424242), length(rest))[1:300])]
n_before = length(pop.rows)
expand!(pop, mine, 2)
c1 = count(x -> x.round == 1, pop.rows)
v1 = count(x -> x.round == 1 && x.status == "ok" && x.allo && x.n == 1, pop.rows)
c2 = count(x -> x.round == 2, pop.rows)
v2 = count(x -> x.round == 2 && x.status == "ok" && x.allo && x.n == 1, pop.rows)
stamp("phantom pop: round 0 $(length(base2)) base; round 1 $c1 candidates, $v1 one-subunit ok " *
      "(study: 2,752 and 1,833); round 2 from $(length(mine)) parents outside the study's 1000: " *
      "$c2 candidates, $v2 one-subunit ok; overlap with study parents $(length(intersect(mine, theirs)))")
open(joinpath(VD, "data", "phantom_pop.tsv"), "w") do io
    for x in pop.rows
        x.status == "ok" || continue
        println(io, join((x.gid, x.round, x.parent, x.move, Int(x.allo), x.n, x.cx, x.np,
                          repr(rawrec(pop.mechs[x.gid]))), '\t'))
    end
end
stamp("wrote phantom_pop.tsv")

# ── Part 3: scan population (raw data only): bi-bi init and round 1, regulated bi-bi seeds,
# ping-pong init ─────────────────────────────────────────────────────────────────────────────
rxnR = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
end
rxnPP = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end
scan = Any[]
append!(scan, (pop.mechs[g] for g in okg(pop, 0)))
r1 = okg(pop, 1)
append!(scan, (pop.mechs[g] for g in r1[sort!(randperm(Xoshiro(11), length(r1))[1:400])]))
segR = unique!(collect(ER.seed_mechanisms(rxnR, Set([:R]), Set{Symbol}())))
append!(scan, segR[sort!(randperm(Xoshiro(12), length(segR))[1:200])])
ppi = [m for m in unique!(collect(ER.init_mechanisms(rxnPP))) if cxA(m) <= 337]
append!(scan, ppi[sort!(randperm(Xoshiro(13), length(ppi))[1:min(150, length(ppi))])])
open(joinpath(VD, "data", "scan_pop.tsv"), "w") do io
    foreach(m -> println(io, repr(rawrec(m))), scan)
end
stamp("wrote scan_pop.tsv: $(length(scan)) mechanisms, $(count(m -> m isa ER.AllostericMechanism, scan)) allosteric")
