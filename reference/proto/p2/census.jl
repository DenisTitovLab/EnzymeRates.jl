# ABOUTME: P2 census over every plain mechanism: cycle classes, group orientation vs the package,
# ABOUTME: flux-free SS steps, open RE cycles, AllParameters counts, dependent choices, sizes, time.
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
using Printf, Statistics
const T = Thermo

"""(group index, stored slot) of a package parameter name: the first step that renders it."""
function pkg_slot(m, name)
    for (gi, g) in enumerate(m.groups), j in g, (r, n) in enumerate(m.step_names[j])
        n === name && return (gi, r)
    end
    nothing
end

"""(group, stored slot) of one of my columns: slot in the stored orientation."""
function my_slot(sol, i)
    c = sol.cols[i]; s = sol.cs[c.j]
    g = findfirst(gr -> c.j in gr, sol.m.groups)
    (g, s.re ? 1 : (s.flipped ? 3 - c.r : c.r))
end

"""A dependent's chemistry description for the textbook check."""
function describe(sol, i)
    c = sol.cols[i]; s = sol.cs[c.j]
    kind = ("bind", "chem", "iso")[s.kind]
    s.re && return "RE $(kind) K"
    against = c.prio[4] == 1
    "SS $(kind) $(against ? "reverse" : "forward")"
end

"""Whether class (group g, stored slot r) is a chemistry member's constant against its physical
forward direction, or an RE chemistry constant."""
function textbook_class(sol, cl)
    g, r = cl
    any(sol.m.groups[g]) do j
        s = sol.cs[j]
        s.kind >= 2 || return false
        s.re && return true
        i = findfirst(c -> c.j == j && c.r == (s.flipped ? 3 - r : r), sol.cols)
        sol.cols[i].prio[4] == 1
    end
end

"""Short description of a class: the most chemical member's kind, RE/SS, against/forward."""
function class_desc(sol, cl)
    g, r = cl
    j = sol.m.groups[g][argmax([sol.cs[j].kind for j in sol.m.groups[g]])]
    s = sol.cs[j]
    kind = ("bind", "chem", "iso")[s.kind]
    s.re && return "RE $(kind) K"
    i = findfirst(c -> c.j == j && c.r == (s.flipped ? 3 - r : r), sol.cols)
    "SS $(kind) $(sol.cols[i].prio[4] == 1 ? "reverse" : "forward")"
end

hal_diff = 0; weg_diff = 0
rows = String[]
nm = 0; mu_hist = Dict{Int, Int}(); hfund_hist = Dict{Int, Int}(); g_hist = Dict{Int, Int}()
orient_bad = 0; orient_groups = 0
zf_mechs = String[]; open_re = String[]
allp_diff = 0; same_depset = 0; diff_depset = 0; examples = String[]
my_textbook = 0; pkg_textbook = 0; nhal_mechs = 0
label_counts = Dict{Symbol, Int}(); fold_counts = Dict{Symbol, Int}()
maxint = 0; times = Float64[]
mine_desc = Dict{String, Int}(); pkg_desc = Dict{String, Int}()
for pop in Pop.populations()
    for m in Pop.load_population(pop)
        m.allosteric && continue
        global nm += 1
        sol = T.solve(m)
        t = minimum(@elapsed(T.solve(m)) for _ in 1:5); push!(times, t)
        global maxint = max(maxint, sol.maxint)
        mu_hist[sol.ncycles] = get(mu_hist, sol.ncycles, 0) + 1
        hfund_hist[sol.nhaldane_fund] = get(hfund_hist, sol.nhaldane_fund, 0) + 1
        g_hist[sol.haldane_gcd] = get(g_hist, sol.haldane_gcd, 0) + 1
        # group orientation: aligned-vs-stored flip must be constant inside a group
        al = T.group_alignment(m, sol.cs, m.groups)
        for g in m.groups
            length(g) < 2 && continue
            # only orientation-sensitive members matter: SS steps and non-binding RE steps
            sens = [j for j in g if !(sol.cs[j].re && sol.cs[j].bindK)]
            isempty(sens) && continue
            global orient_groups += 1
            fl = [sol.cs[j].flipped ⊻ al[j][2] for j in sens]
            allequal(fl) || (global orient_bad += 1; println("ORIENT ", m.id, " group ", g))
        end
        # flux-free SS steps
        z = T.zero_flux_ss(m)
        any(z) && push!(zf_mechs, string(m.id, " [", m.source, "] steps ", findall(z)))
        # open RE cycles: RE-only fundamental cycles not implied by the group rows
        cs = sol.cs
        colidx = Dict((c.j, c.r) => i for (i, c) in enumerate(sol.cols))
        G = T.group_rows(cs, colidx, al)
        re_cs = [c.re ? c : nothing for c in cs]
        reidx = [j for (j, c) in enumerate(cs) if c.re]
        if !isempty(reidx)
            sub = cs[reidx]
            cyc = T.fundamental_cycles(length(m.species), sub)
            nopen = 0
            for z0 in cyc
                zfull = zeros(Int, length(cs)); zfull[reidx] .= z0
                row = T.cycle_row(cs, colidx, zfull)
                A = zeros(T.Q, length(G) + 1, length(colidx))
                for (k, r) in enumerate(G); A[k, :] .= r; end
                rG = length(T.rref!(copy(A[1:length(G), :]))[1])
                A[end, :] .= row
                rA = length(T.rref!(A)[1])
                rA > rG && (nopen += 1)
            end
            nopen > 0 && push!(open_re, "$(m.id) open=$(nopen) of $(length(cyc))")
        end
        # AllParameters count vs the package's Full list
        alls = T.solve(m; mode = :all)
        nall = length(alls.indep)
        nall == length(m.full_params) - 1 || (global allp_diff += 1;
            println("ALLP ", m.id, " mine=", nall, " pkg=", length(m.full_params) - 1))
        # labels and folds
        for d in sol.deps
            label_counts[d.label] = get(label_counts, d.label, 0) + 1
            d.fold && (fold_counts[d.label] = get(fold_counts, d.label, 0) + 1)
        end
        # dependent choice vs the package, as (group, stored slot) of non-kinetic-group dependents
        mine = Set(my_slot(sol, d.col) for d in sol.deps if d.label !== :kinetic_group)
        pk = Set{Tuple{Int, Int}}()
        for (n, _) in m.dependents; push!(pk, pkg_slot(m, n)); end
        for (a, _) in m.pass2_rename; push!(pk, pkg_slot(m, a)); end
        if mine == pk
            global same_depset += 1
        else
            global diff_depset += 1
            length(examples) < 400 && push!(examples, string(m.id, " [", m.source, "]\n      mine: ",
                join(sort([string(sol.cols[d.col].name, " (", d.label, ", ", describe(sol, d.col), ")")
                           for d in sol.deps if d.label !== :kinetic_group]), "; "),
                "\n      pkg:  ", join(sort(string.(first.(m.dependents))), "; "),
                isempty(m.pass2_rename) ? "" : "  pass2: " * join(string.(first.(m.pass2_rename)), ", ")))
        end
        # textbook, by class (kinetic group + stored slot): a Haldane class is textbook when a
        # chemistry member (fused, Theorell-Chance, isomerization) of the group runs this slot
        # against its physical forward direction (or the class is an RE chemistry K)
        hal = [my_slot(sol, d.col) for d in sol.deps if d.label === :haldane]
        pkh = [pkg_slot(m, n) for (n, mono) in m.dependent_monomials if any(p -> p.first === :Keq, mono)]
        if !isempty(hal)
            global nhal_mechs += 1
            all(cl -> textbook_class(sol, cl), hal) && (global my_textbook += 1)
            all(cl -> textbook_class(sol, cl), pkh) && (global pkg_textbook += 1)
            for cl in hal; k = class_desc(sol, cl); mine_desc[k] = get(mine_desc, k, 0) + 1; end
            for cl in pkh; k = class_desc(sol, cl); pkg_desc[k] = get(pkg_desc, k, 0) + 1; end
            Set(hal) == Set(pkh) || (global hal_diff += 1)
        end
        wmine = Set(my_slot(sol, d.col) for d in sol.deps if d.label === :wegscheider)
        wpkg = setdiff(pk, Set(pkh))
        wmine == wpkg || (global weg_diff += 1)
    end
end
println("plain mechanisms: ", nm)
println("cycle-space size mu: ", sort(collect(mu_hist)))
println("fundamental Haldane cycles per mechanism: ", sort(collect(hfund_hist)))
println("gcd of Haldane multiples: ", sort(collect(g_hist)))
println("groups with orientation-sensitive members: ", orient_groups, "; aligned-vs-stored inconsistent: ", orient_bad)
println("AllParameters count != package Full count: ", allp_diff)
println("mechanisms with flux-free SS steps: ", length(zf_mechs)); foreach(x -> println("   ", x), zf_mechs)
println("mechanisms with open RE cycles: ", length(open_re)); foreach(x -> println("   ", x), open_re)
println("dependent labels: ", label_counts, "  folds by label: ", fold_counts)
println("dependent set (group, stored slot) same as package: ", same_depset, ", different: ", diff_depset)
println("Haldane class differs from the package: ", hal_diff, "; Wegscheider classes differ: ", weg_diff)
println("Haldane dependents textbook (reverse of a chemistry step): mine ", my_textbook, " / ", nhal_mechs,
        ", package ", pkg_textbook, " / ", nhal_mechs)
println("my Haldane dependent kinds: ", sort(collect(mine_desc)))
println("package Haldane dependent kinds: ", sort(collect(pkg_desc)))
println("max |num|,|den| in RREF: ", maxint)
@printf("solve time: median %.1f us, max %.1f us, total %.3f s\n", 1e6 * median(times), 1e6 * maximum(times), sum(times))
open(joinpath(@__DIR__, "depset_examples.txt"), "w") do io
    foreach(x -> println(io, x), examples)
end
