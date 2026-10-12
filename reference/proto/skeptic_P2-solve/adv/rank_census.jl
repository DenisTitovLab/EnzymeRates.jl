# ABOUTME: Skeptic of probe P2: count vs FD identifiable rank on EVERY plain mechanism (the probe
# ABOUTME: sampled 40 enumerated ones), before and after the flux-free conversion; stdlib only.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/AdvLib.jl")
using Random, Printf

rng = Xoshiro(4242)
hist = Dict{String, Dict{Int, Int}}()
examples = Dict{String, Vector{String}}()
zfbridge = String[]
t0 = time()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    sol = Thermo.solve(m)
    m2, nconv = Thermo.with_flux_free_as_re(m)
    s2 = nconv > 0 ? Thermo.solve(m2) : sol
    r, kept, dropped = Adv.fd_rank(m2, s2, rng)
    gap = length(s2.indep) - r
    h = get!(hist, pop, Dict{Int, Int}()); h[gap] = get(h, gap, 0) + 1
    if gap > 0
        ex = get!(examples, pop, String[])
        length(ex) < 12 && push!(ex, @sprintf("%s count %d (after flux-free %d) rank %d kept %.0e dropped %.0e  %s",
                                              m.id, length(sol.indep), length(s2.indep), r, kept, dropped,
                                              join(Pop.step_text.(Ref(m), m.steps), " | ")))
    end
    z = Thermo.zero_flux_ss(m)
    if any(z)
        b = Adv.plain_bridges(m)
        whole = [g for g in m.groups if all(j -> !m.steps[j].re && z[j], g)]
        isempty(whole) || push!(zfbridge, string(m.id, " whole flux-free groups ", whole,
                                                 " plain bridges among them ", [j for g in whole for j in g if b[j]]))
    end
end
@printf("plain mechanisms, count after the flux-free conversion minus FD rank (%.0f s):\n", time() - t0)
for (pop, h) in sort!(collect(hist); by = first)
    println("  ", rpad(pop, 14), sort!(collect(h)))
end
for (pop, ex) in sort!(collect(examples); by = first)
    println("== ", pop); foreach(x -> println("   ", x), ex)
end
println("whole flux-free SS groups and whether they are plain bridges:"); foreach(println, zfbridge)
