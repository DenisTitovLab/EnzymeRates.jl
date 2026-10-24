# ABOUTME: Regenerates the L-phantom study's populations with the phase-2 package: identify's base
# ABOUTME: tier (round 0), every entry expanded (round 1), a sampled round 2; neutral text, moves, parents.
# usage: gen.jl <out dir> <reaction> <last round> <round-2 parents> [round-1 parents]
using EnzymeRates, Random, Printf
const ER = EnzymeRates
const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/census"
include(joinpath(CEN, "census_common.jl"))
include(joinpath(CEN, "neutral_steps.jl"))

const OUT, RXN, LAST, NPAR2 = ARGS[1], ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4])
const NPAR1 = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : typemax(Int)
const CX_MAX = 337          # identify's eq_complexity_filter default
const NP_MAX = 20           # identify's max_param_count default

# the census's bi_bi_rxn and bi_bi_pp_rxn (WT test/test_mechanism_enumeration.jl:205-213)
const REACTIONS = Dict(
    "bibi" => @enzyme_reaction(begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end),
    "pingpong" => @enzyme_reaction(begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end),
    # regulator populations for the pre-filter's site scope (key.md item 9)
    "uniuni_R" => @enzyme_reaction(begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
    end),
    "bibi_R" => @enzyme_reaction(begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        allosteric_regulators: R
    end))

"""identify's seeds: `init_mechanisms`, or `seed_mechanisms` with every declared regulator
required (WT src/identify_rate_equation.jl `_required_regulators`)."""
function identify_seeds(rxn)
    req_allo, req_comp = ER._required_regulators(rxn, Symbol[], Symbol[])
    isempty(req_allo) && isempty(req_comp) ? unique!(collect(ER.init_mechanisms(rxn))) :
        unique!(collect(ER.seed_mechanisms(rxn, req_allo, req_comp)))
end

"""`expand_mechanisms([m], rxn)` child by child with its move, in the same order, with the
same regulator-type filter and atom check (WT src/mechanism_enumeration.jl:1920-1940)."""
function expand_labeled(m, rxn)
    ER._assert_emission_rules(m)
    out = Tuple{Any, Symbol}[]
    add!(cs, mv) = foreach(c -> push!(out, (c, mv)), cs)
    add!(ER._expand_re_to_ss(m), :re_to_ss)
    add!(ER._expand_split_kinetic_group(m), :split)
    add!(ER._expand_add_dead_end_regulator(m, rxn), :dead_end)
    if m isa ER.Mechanism
        add!(ER._expand_to_allosteric(m, rxn), :to_allosteric)
    else
        add!(ER._expand_add_allosteric_regulator(m, rxn), :add_allo_reg)
        add!(ER._expand_change_allo_state(m), :change_allo_state)
        add!(ER._expand_merge_regulatory_sites(m), :merge_sites)
    end
    filter!(cm -> ER._respects_reg_type(cm[1], rxn), out)
    foreach(cm -> ER._assert_atom_conserving(cm[1]), out)
    out
end

mutable struct Pop
    mechs::Vector{Any}
    recs::Vector{NamedTuple}
    gid::Dict{Any, Int}
    edges::Vector{Tuple{Int, Int, Symbol}}
end

"""Register `m` (first parent `par`, move `mv`) unless seen; returns its gid."""
function register!(pop::Pop, m, round, par, mv)
    g = get(pop.gid, m, 0)
    g == 0 || return g
    g = length(pop.mechs) + 1
    pop.gid[m] = g
    push!(pop.mechs, m)
    allo = m isa ER.AllostericMechanism
    n = allo ? ER.catalytic_multiplicity(m) : 0
    cx = np = -1
    status = "ok"
    try
        cx = Int(ER._eq_complexity(m))
        np = ER._independent_param_count(m)
        status = cx > CX_MAX ? "cx_skip" : np > NP_MAX ? "np_skip" : "ok"
    catch e
        status = "error: " * first(sprint(showerror, e), 120)
    end
    push!(pop.recs, (; gid = g, round, parent = par, move = mv, allo, n, cx, np, status,
                     text = census_text(m)))
    g
end

ok(pop, round) = [r.gid for r in pop.recs if r.round == round && r.status == "ok"]

function expand!(pop::Pop, rxn, parents, round)
    nfail = 0
    for p in parents
        kids = try
            expand_labeled(pop.mechs[p], rxn)
        catch e
            nfail += 1
            println("  expand error at gid $p: ", first(sprint(showerror, e), 200))
            continue
        end
        for (c, mv) in kids
            g = register!(pop, c, round, p, mv)
            push!(pop.edges, (p, g, mv))
        end
    end
    nfail
end

function main()
    rxn = REACTIONS[RXN]
    pop = Pop(Any[], NamedTuple[], Dict{Any, Int}(), Tuple{Int, Int, Symbol}[])
    t = @elapsed begin
        seeds = identify_seeds(rxn)
        base, fails = ER._base_tier(seeds, rxn)
        for m in base
            register!(pop, m, 0, 0, m in seeds ? :seed : :seed_flip)
        end
    end
    @printf("%s round 0: %d seeds, %d base entries, %d ok, %d base failures (%.1f s)\n", RXN,
            length(seeds), length(base), length(ok(pop, 0)), length(fails), t)
    for round in 1:LAST
        pool = ok(pop, round - 1)
        cap = round == 1 ? NPAR1 : NPAR2
        parents = length(pool) <= cap ? pool :
                  pool[sort!(randperm(Xoshiro(60), length(pool))[1:cap])]
        n0 = length(pop.recs)
        t = @elapsed nfail = expand!(pop, rxn, parents, round)
        new = pop.recs[(n0 + 1):end]
        @printf("%s round %d: %d parents of %d, %d fresh, %d ok, %d allosteric n=1 ok, %d expand errors (%.1f s, maxrss %.0f MB)\n",
                RXN, round, length(parents), length(pool), length(new),
                count(r -> r.status == "ok", new),
                count(r -> r.status == "ok" && r.allo && r.n == 1, new), nfail, t,
                Sys.maxrss() / 2^20)
        flush(stdout)
    end
    open(joinpath(OUT, "recs_$RXN.tsv"), "w") do io
        write_row(io, "gid", "round", "parent", "move", "allo", "n", "cx", "np", "status", "text")
        for r in pop.recs
            write_row(io, string(r.gid), string(r.round), string(r.parent), string(r.move),
                      string(Int(r.allo)), string(r.n), string(r.cx), string(r.np), r.status,
                      r.text)
        end
    end
    open(joinpath(OUT, "edges_$RXN.tsv"), "w") do io
        for (p, c, mv) in pop.edges
            write_row(io, string(p), string(c), string(mv))
        end
    end
    @printf("done %s: %d records, %d edges, maxrss %.0f MB\n", RXN, length(pop.recs),
            length(pop.edges), Sys.maxrss() / 2^20)
end

main()
