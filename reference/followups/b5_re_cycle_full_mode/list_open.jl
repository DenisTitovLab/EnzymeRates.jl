# ABOUTME: Lists the plain mechanisms whose RE cycles the kinetic grouping leaves open,
# ABOUTME: with their size and the Full-mode constants absent from the package's Full polynomials.
const C = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common"
include(joinpath(C, "Population.jl"))
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "handwritten"]
syms(p) = Set(s for (mono, _) in p for (s, _) in mono)
rows = []
for pop in POPS, m in Pop.load_population(pop)
    m.allosteric && continue
    cyc, open = Pop.re_cycle_closure(m)
    open == 0 && continue
    used = union(syms(m.pkg_num), syms(m.pkg_den))
    missing_full = [p for p in m.full_params if p != :E_total && p ∉ used]
    push!(rows, (m.id, m.source, length(m.steps), length(m.species), length(m.groups), cyc, open,
                 missing_full, m.pass2_rename))
end
sort!(rows; by = r -> (r[3], r[4], r[1]))
println("open-RE-cycle plain mechanisms: ", length(rows))
for r in rows
    println(rpad(r[1], 18), " steps=", r[3], " forms=", r[4], " groups=", r[5], " cycles=", r[6],
            " open=", r[7], " missing_in_Full_v=", r[8], " pass2=", r[9], "  [", r[2], "]")
end
