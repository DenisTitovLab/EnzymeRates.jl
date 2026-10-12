# ABOUTME: Skeptic of probe P2: random plain mechanisms (two conformations, an inhibitor, fused,
# ABOUTME: Theorell-Chance, exchange and iso steps, random groups) checked against Thermo.solve.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf
setprecision(BigFloat, 512)

const S_ = Pop.PMet(:S, :S); const P_ = Pop.PMet(:P, :P); const I_ = Pop.PMet(:I, :Inh)
const NONE = Pop.PMet[]
const FLIP = get(ARGS, 2, "") == "flip"     # store each step in a random direction
const DUMP = get(ARGS, 3, "fuzz_dump")
key(c, b) = (c, sort(b; by = x -> x.name))

"""Candidate steps from form (c, b): (to form, consumed, released, kind, chemistry?)."""
function moves(c, b)
    out = Tuple{Any, Vector{Pop.PMet}, Vector{Pop.PMet}, Symbol, Bool}[]
    if length(b) < 2
        for M in (S_, P_, I_)
            push!(out, (key(c, [b; M]), [M], NONE, Symbol(:bind, M.name), false))
            push!(out, (key(c == :E ? :F : :E, [b; M]), [M], NONE, Symbol(:bindconf, M.name),
                        false))
        end
    end
    i = findfirst(==(S_), b)
    if i !== nothing
        b2 = copy(b); b2[i] = P_
        push!(out, (key(c, b2), NONE, NONE, :chem, true))
        push!(out, (key(c == :E ? :F : :E, b2), NONE, NONE, :chemconf, true))
        push!(out, (key(c, deleteat!(copy(b), i)), NONE, [P_], :fusedrelease, true))
    end
    push!(out, (key(c == :E ? :F : :E, b), NONE, NONE, :conf, false))
    length(b) < 2 && push!(out, (key(c, [b; P_]), [S_], NONE, :fusedbind, true))
    push!(out, (key(c == :E ? :F : :E, b), [S_], [P_], :tc, true))
    j = findfirst(==(I_), b)
    if j !== nothing
        b3 = copy(b); b3[j] = S_
        push!(out, (key(c, b3), [S_], [I_], :exchange, false))
    end
    out
end

function random_mech(rng, id)
    forms = Any[key(:E, Pop.PMet[])]
    fidx = Dict(forms[1] => 1)
    steps = Tuple{Int, Int, Vector{Pop.PMet}, Vector{Pop.PMet}, Symbol, Bool}[]
    seen = Set{Tuple}()
    function add!(f, t, cons, rel, kind, chem)
        t in keys(fidx) || (push!(forms, t); fidx[t] = length(forms))
        a, b = fidx[f], fidx[t]
        a == b && return false
        sides = minmax((a, sort(string.(getfield.(cons, :name)))), (b, sort(string.(getfield.(rel, :name)))))
        sides in seen && return false
        push!(seen, sides); push!(steps, (a, b, cons, rel, kind, chem)); true
    end
    # a catalytic core
    core = rand(rng, 1:4)
    E0 = key(:E, Pop.PMet[])
    if core == 1                      # ordered uni-uni
        add!(E0, key(:E, [S_]), [S_], NONE, :bindS, false)
        add!(key(:E, [S_]), key(:E, [P_]), NONE, NONE, :chem, true)
        add!(key(:E, [P_]), E0, NONE, [P_], :bindP, false)
    elseif core == 2                  # Theorell-Chance-like through F
        add!(E0, key(:F, Pop.PMet[]), [S_], [P_], :tc, true)
        add!(key(:F, Pop.PMet[]), E0, NONE, NONE, :conf, false)
    elseif core == 3                  # iso uni-uni
        add!(E0, key(:E, [S_]), [S_], NONE, :bindS, false)
        add!(key(:E, [S_]), key(:F, [P_]), NONE, NONE, :chemconf, true)
        add!(key(:F, [P_]), key(:F, Pop.PMet[]), NONE, [P_], :bindP, false)
        add!(key(:F, Pop.PMet[]), E0, NONE, NONE, :conf, false)
    else                              # fused binding then release
        add!(E0, key(:E, [P_]), [S_], NONE, :fusedbind, true)
        add!(key(:E, [P_]), E0, NONE, [P_], :bindP, false)
    end
    for _ in 1:rand(rng, 1:6)
        f = rand(rng, forms)
        mv = moves(f...)
        isempty(mv) && continue
        t, cons, rel, kind, chem = rand(rng, mv)
        add!(f, t, cons, rel, kind, chem)
    end
    species = [Pop.PSpecies(b, c, Symbol[], Symbol[], Symbol(c, join(string.(getfield.(b, :name)))))
               for (c, b) in forms]
    re = [!chem && rand(rng) < 0.5 for (_, _, _, _, _, chem) in steps]
    # stored in the generated forward direction, as the package's canonical form stores
    # these kinds (the scrambles below test every other direction)
    psteps = [!FLIP || rand(rng, Bool) ? Pop.PStep(a, b, cons, rel, r) : Pop.PStep(b, a, rel, cons, r)
              for ((a, b, cons, rel, _, _), r) in zip(steps, re)]
    # random groups of same-kind non-iso steps with the same flag
    groups = [[j] for j in eachindex(steps)]
    for _ in 1:3
        i, j = rand(rng, eachindex(steps)), rand(rng, eachindex(steps))
        ki, kj = steps[i][5], steps[j][5]
        (i != j && ki == kj && re[i] == re[j] && !(ki in (:chem, :chemconf, :conf))) || continue
        gi = findfirst(g -> i in g, groups); gj = findfirst(g -> j in g, groups)
        gi == gj && continue
        append!(groups[gi], groups[gj]); deleteat!(groups, gj)
    end
    mets = any(s -> any(x -> x.name === :I, [s[3]; s[4]; forms[s[1]][2]; forms[s[2]][2]]), steps) ?
           [:S, :P, :I] : [:S, :P]
    Pop.PMech(; id, source = id, substrates = [:S], products = [:P],
              inhibitors = :I in mets ? [:I] : Symbol[], metabolites = mets,
              species, steps = psteps, groups)
end

rng = Xoshiro(31337)
N = parse(Int, get(ARGS, 1, "400"))
stats = Dict{Symbol, Int}()
bump(k) = (stats[k] = get(stats, k, 0) + 1)
fails = String[]
gaps = Dict{Int, Int}()
dump = Pop.PMech[]
t0 = time()
for n in 1:N
    m = random_mech(rng, "fuzz/$n")
    sol = try
        Thermo.solve(m)
    catch e
        bump(:solve_error); push!(fails, "SOLVE-ERR $(m.id): " * first(sprint(showerror, e), 200)); continue
    end
    bump(:solved)
    tc = Adv.true_count(m)
    tc.contradiction && bump(:contradiction)
    tc.count == length(sol.indep) || (bump(:count_mismatch); push!(dump, m);
        push!(fails, "COUNT $(m.id) potentials $(tc.count) probe $(length(sol.indep))"))
    k, keq, _ = Adv.random_constants(sol, rng)
    Adv.ties_hold(m, k) || (bump(:tie_fail); push!(fails, "TIES $(m.id)"); push!(dump, m))
    v0 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq))
    v1 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq / 10))
    abs(v0) < big(10)^-30 * abs(v1) || (bump(:db_fail); push!(fails, @sprintf("DB %s %.1e", m.id, Float64(abs(v0) / abs(v1)))))
    sig = Adv.signature(sol)
    for _ in 1:5
        Adv.signature(Thermo.solve(Adv.scramble(m, rng))) == sig || (bump(:scramble_fail); push!(fails, "SCRAMBLE $(m.id)"); break)
    end
    es = Adv.exponent_stats(sol)
    es.nonint > 0 && bump(:nonint_exponents)
    es.maxint > 2 && bump(:rref_entry_gt2)
    z = Thermo.zero_flux_ss(m)
    J, v = Adv.step_fluxes(m, k, Oracle.concs_at_ratio(m, rng, 10keq))
    sc = maximum(abs(x) for x in J if x !== nothing; init = abs(v))
    for (j, x) in enumerate(J)
        x === nothing && continue
        if z[j] && abs(x) > big(10)^-30 * sc
            bump(:flagged_but_flux); push!(fails, "FALSE-FLUX-FREE $(m.id) step $j")
        elseif !z[j] && abs(x) < big(10)^-30 * sc
            bump(:zero_flux_unflagged); push!(fails, "UNFLAGGED-ZERO $(m.id) step $j: " * Pop.step_text(m, m.steps[j]))
            push!(dump, m)
        end
    end
    any(z) && bump(:has_flux_free)
    if n % 4 == 0
        m2, nconv = Thermo.with_flux_free_as_re(m)
        s2 = nconv > 0 ? Thermo.solve(m2) : sol
        r = Adv.fd_rank(m2, s2, rng)[1]
        g = length(s2.indep) - r
        gaps[g] = get(gaps, g, 0) + 1
    end
end
@printf("fuzz: %d mechanisms in %.0f s\n", N, time() - t0)
println("stats: ", sort!(collect(stats); by = first))
println("count after flux-free conversion minus FD rank (every 4th): ", sort!(collect(gaps)))
foreach(println, fails[1:min(end, 60)])
Pop.save_population(DUMP, unique(dump); path = joinpath(PROTO, "skeptic_P2-solve", "data", DUMP * ".jls"))
