# ABOUTME: P2 extra checks: arithmetic size, pinned constants, Keq-bearing dependents, the one-Haldane-row
# ABOUTME: claim, AllParameters with RE cycles closed, and the solve stressed on allosteric catalytic graphs.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
using Random
const T = Thermo
include(joinpath(@__DIR__, "determinism_helpers.jl"))

function one_haldane_row(m, sol)
    # RREF of [rhs | cycle rows] with rhs first: exactly one row keeps a nonzero rhs
    cs = sol.cs
    colidx = Dict((c.j, c.r) => i for (i, c) in enumerate(sol.cols))
    cyc = T.fundamental_cycles(length(m.species), cs)
    deltas = T.step_delta(cs)
    M = zeros(T.Q, length(cyc), length(sol.cols) + 1)
    for (i, z) in enumerate(cyc)
        M[i, 1] = T.cycle_multiple(m, deltas, z); M[i, 2:end] .= T.cycle_row(cs, colidx, z)
    end
    piv, _ = T.rref!(M)
    count(r -> M[r, 1] != 0, axes(M, 1))
end

function stats(label, ms; perms = 0, rng = Xoshiro(1))
    maxint = 0; maxexp = 0 // 1; nonint = 0; pin1 = 0; pinK = 0; multiK = 0; nm = 0
    onerow = Dict{Int, Int}(); permfail = 0; errs = 0; allre_diff = 0; allre_mechs = 0
    for m in ms
        nm += 1
        sol = try
            T.solve(m)
        catch e
            errs += 1; println("ERR ", m.id, ": ", sprint(showerror, e)[1:min(end, 160)]); continue
        end
        maxint = max(maxint, sol.maxint)
        for d in sol.deps
            for (_, e) in d.terms
                abs(e) > maxexp && (maxexp = abs(e)); isinteger(e) || (nonint += 1)
            end
            isempty(d.terms) && (d.keq == 0 ? (pin1 += 1) : (pinK += 1))
        end
        count(d -> d.keq != 0, sol.deps) > 1 && (multiK += 1)
        k = one_haldane_row(m, sol); onerow[k] = get(onerow, k, 0) + 1
        for _ in 1:perms
            signature(T.solve(permuted(m, rng))) == signature(sol) || (permfail += 1)
        end
    end
    println(label, ": mechanisms $nm, errors $errs, max |num|,|den| in RREF $maxint, max |exponent| $maxexp, ",
            "non-integer exponents $nonint, pinned to 1: $pin1, pinned to Keq^c: $pinK, ",
            "mechanisms with >1 Keq-bearing dependent: $multiK, Haldane rows after RREF with rhs first: ",
            sort(collect(onerow)), perms > 0 ? ", permutation failures $permfail of $(perms * nm)" : "")
end

plain = [m for p in Pop.populations() for m in Pop.load_population(p) if !m.allosteric]
allo = [m for p in Pop.populations() for m in Pop.load_population(p) if m.allosteric]
stats("plain", plain)
stats("allosteric catalytic graphs as plain", allo; perms = 3)
# AllParameters with the RE-only cycle rows added (closes open RE cycles)
nd = 0; extra_deps = Dict{Int, Int}()
for m in plain
    a = T.solve(m; mode = :all)
    cs = a.cs
    colidx = Dict((c.j, c.r) => i for (i, c) in enumerate(a.cols))
    reidx = [j for (j, c) in enumerate(cs) if c.re]
    cyc = T.fundamental_cycles(length(m.species), cs[reidx])
    al = T.group_alignment(m, cs, m.groups)
    rows = T.group_rows(cs, colidx, al)
    for z0 in cyc
        z = zeros(Int, length(cs)); z[reidx] .= z0
        push!(rows, T.cycle_row(cs, colidx, z))
    end
    M = zeros(T.Q, length(rows), length(a.cols)); for (i, r) in enumerate(rows); M[i, :] .= r; end
    r1 = length(T.rref!(M)[1])
    d = length(a.cols) - r1
    if d != length(a.indep)
        global nd += 1; k = length(a.indep) - d; extra_deps[k] = get(extra_deps, k, 0) + 1
    end
end
println("AllParameters + RE-cycle rows: mechanisms whose count drops $nd, drop sizes ", sort(collect(extra_deps)))
