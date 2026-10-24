# ABOUTME: R5.5 populations: identify's seed tier and expansion rounds of a reaction, each candidate's
# ABOUTME: V×τ per conformation. Args: <rxn> <mode uniform|largest> <parent quota per round>...
using EnzymeRates, Random, Printf, SHA
const ER = EnzymeRates
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/census"
include(joinpath(CEN, "census_common.jl"))
include(joinpath(CEN, "neutral_steps.jl"))
const WT = dirname(dirname(pathof(EnzymeRates)))
println("EnzymeRates ", pkgversion(EnzymeRates), " at ",
        readchomp(`git -C $WT rev-parse --short HEAD`), "; Julia ", VERSION)

# bibiR: the census bi-bi without inhibitors plus a regulator (two subunits allowed);
# pingpong: the census ping-pong reaction (phase3/census/corpus.jl:40-48); bibiR1 and
# pingpong1: the same with one subunit; bibi: the census bi-bi (no regulator).
const RXN = Dict(
    "bibiR" => """
@enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
    oligomeric_state: 2
end""",
    "pingpong" => """
@enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
    shared_catalytic_site: (A, P)
end""",
    "bibiR1" => """
@enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
end""",
    "pingpong1" => """
@enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    shared_catalytic_site: (A, P)
end""",
    "bibi" => """
@enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    competitive_inhibitors: A, I
end""")

const FILTER = 337
cxA(m) = Int(ER._eq_complexity(m))
cxI(m) = m isa ER.AllostericMechanism ? Int(ER._eq_complexity(ER._inactive_groups(m))) : 0
"""Whether the inactive conformation's law is the active one's: every step kept and every
group `:EqualAI`, so both conformations read the same constants on the same graph."""
same_law(m) = m isa ER.AllostericMechanism &&
    all(==(:EqualAI), ER.cat_allo_states(m)) &&
    sum(length, ER._inactive_groups(m)) == sum(length, ER.steps(m))

function record(io, round, m)
    t = census_text(m)
    a, i = cxA(m), cxI(m)
    s, p, _ = census_roles(m)
    allo = m isa ER.AllostericMechanism
    write_row(io, string(round), mech_id(t), allo ? "1" : "0",
              string(allo ? ER.catalytic_multiplicity(m) : 0), string(a), string(i),
              same_law(m) ? "1" : "0", join(s, ","), join(p, ","), t)
    a
end

rxname, mode, quotas = ARGS[1], ARGS[2], parse.(Int, ARGS[3:end])
rxn = Core.eval(Main, Meta.parse(RXN[rxname]))
out = joinpath(@__DIR__, "data", "r55_$(rxname)_$mode.tsv")
rng = Xoshiro(20261014)
t0 = time()
open(out, "w") do io
    # identify's seeds: fully regulated when the reaction declares a regulator; the bi-bi
    # seeds as if both competitive inhibitors were optional (init_mechanisms)
    seeds = rxname == "bibi" ? unique!(collect(ER.init_mechanisms(rxn))) :
            unique!(collect(ER.seed_mechanisms(rxn, Set([:R]), Set{Symbol}())))
    base, fails = ER._base_tier(seeds, rxn)
    seen = Set{UInt64}(hash.(base))
    pool = [m for m in base if record(io, 0, m) <= FILTER]
    @printf("%s round 0: %d seeds, %d base, %d base-tier failures, %d within %d; %.0f s, maxrss %.0f MB\n",
            rxname, length(seeds), length(base), length(fails), length(pool), FILTER,
            time() - t0, Sys.maxrss() / 2^20)
    for (round, quota) in enumerate(quotas)
        # uniform: a random sample; largest: the quota with the largest V×τ (ties by order)
        parents = length(pool) <= quota ? pool :
                  mode == "uniform" ? pool[sort(randperm(rng, length(pool))[1:quota])] :
                  pool[sort(sortperm(cxA.(pool); rev = true, alg = MergeSort)[1:quota])]
        next, nkids, errs = Any[], 0, 0
        for p in parents
            kids = try
                ER.expand_mechanisms([p], rxn)
            catch e
                errs += 1
                continue
            end
            for c in kids
                h = hash(c)
                h in seen && continue
                push!(seen, h)
                nkids += 1
                record(io, round, c) <= FILTER && round < length(quotas) && push!(next, c)
            end
        end
        @printf("%s round %d: %d parents (of %d), %d new children, %d parent errors; %.0f s, maxrss %.0f MB\n",
                rxname, round, length(parents), length(pool), nkids, errs, time() - t0,
                Sys.maxrss() / 2^20)
        flush(stdout)
        pool = next
    end
end
println("wrote ", out)
