# ABOUTME: Census of mixed kinetic groups (a plain binding and a chemistry step in one group)
# ABOUTME: and of rapid-equilibrium chemistry steps over every exported population; stdlib only.
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
using Printf

metkey(x::Pop.PMet) = (x.name, x.role)

"""Step class: :plain (one metabolite consumed, none released, `to` = `from` plus it with
the same residual), :fused (a binding that is not plain), :iso, :tc (one taken up, one given
off), :other."""
function step_class(m::Pop.PMech, s::Pop.PStep)
    nc, nr = length(s.consumed), length(s.released)
    nc == 0 && nr == 0 && return :iso
    if nc == 1 && nr == 0
        a, b = m.species[s.from], m.species[s.to]
        plain = sort(metkey.([a.bound; s.consumed])) == sort(metkey.(b.bound)) &&
                sort(a.res_added) == sort(b.res_added) &&
                sort(a.res_subtracted) == sort(b.res_subtracted)
        return plain ? :plain : :fused
    end
    nc == 1 && nr == 1 && return :tc
    :other
end

"""Whether a form holds both a substrate and a product (a dead-end complex)."""
deadend(sp::Pop.PSpecies) = any(x -> x.role === :S, sp.bound) && any(x -> x.role === :P, sp.bound)

group_text(m, g) = join([Pop.step_text(m, m.steps[j]) for j in m.groups[g]], ", ")

const POPS = ["init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "allo_uniuni", "allo_bibi", "specs", "handwritten"]

shapes = Dict{String, Dict{String, Int}}()      # population => signature => count
examples = Dict{String, Vector{String}}()
for name in POPS
    mechs = Pop.load_population(name)
    nmech = length(mechs)
    nonuniform = 0
    n_mixed_mech = 0; n_mech_re = 0; n_mech_ss = 0
    n_mixed_groups = 0; n_groups_re = 0; n_groups_ss = 0
    tagcount = Dict{Symbol, Int}()
    mech_onlya_mixed = 0
    chem_member_classes = Dict{Symbol, Int}()
    plain_member_deadend = 0; plain_member_catalytic = 0
    # RE chemistry census (every step, mixed or not)
    re_chem_mech = Dict{Symbol, Int}()
    re_fused_in_mixed_mech = 0
    sh = get!(shapes, name, Dict{String, Int}())
    for m in mechs
        cls = [step_class(m, s) for s in m.steps]
        for k in unique(cls[j] for j in eachindex(cls) if m.steps[j].re && cls[j] !== :plain)
            re_chem_mech[k] = get(re_chem_mech, k, 0) + 1
        end
        has_mixed = false; has_re = false; has_ss = false; has_onlya = false; re_fused_mixed = false
        for (g, grp) in enumerate(m.groups)
            flags = unique(m.steps[j].re for j in grp)
            length(flags) == 1 || (nonuniform += 1)
            any(j -> cls[j] === :plain, grp) && any(j -> cls[j] !== :plain, grp) || continue
            has_mixed = true; n_mixed_groups += 1
            re = only(flags)
            re ? (n_groups_re += 1; has_re = true) : (n_groups_ss += 1; has_ss = true)
            for j in grp
                if cls[j] === :plain
                    deadend(m.species[m.steps[j].to]) ? (plain_member_deadend += 1) :
                                                        (plain_member_catalytic += 1)
                else
                    chem_member_classes[cls[j]] = get(chem_member_classes, cls[j], 0) + 1
                    re && cls[j] === :fused && (re_fused_mixed = true)
                end
            end
            tag = m.allosteric ? m.group_tags[g] : :none
            tagcount[tag] = get(tagcount, tag, 0) + 1
            tag === :OnlyA && (has_onlya = true)
            sig = (re ? "RE  " : "SS  ") * (m.allosteric ? "[$(tag)] " : "") * group_text(m, g)
            sh[sig] = get(sh, sig, 0) + 1
            ex = get!(examples, name * (re ? "/RE" : "/SS"), String[])
            length(ex) < 3 && push!(ex, m.id * "  " * sig)
        end
        has_mixed && (n_mixed_mech += 1)
        has_re && (n_mech_re += 1)
        has_ss && (n_mech_ss += 1)
        has_onlya && (mech_onlya_mixed += 1)
        re_fused_mixed && (re_fused_in_mixed_mech += 1)
    end
    @printf("== %-13s %5d mechanisms | groups with two flags: %d\n", name, nmech, nonuniform)
    @printf("   mechanisms with a mixed group: %d (with an RE mixed group %d, with an SS mixed group %d, both %d)\n",
            n_mixed_mech, n_mech_re, n_mech_ss, n_mech_re + n_mech_ss - n_mixed_mech)
    @printf("   mixed groups: %d (RE %d, SS %d); chemistry members by class %s; plain members dead-end %d, catalytic %d\n",
            n_mixed_groups, n_groups_re, n_groups_ss, sort!(collect(chem_member_classes)),
            plain_member_deadend, plain_member_catalytic)
    any(m -> m.allosteric, mechs) &&
        @printf("   mixed groups by tag %s; mechanisms with an :OnlyA mixed group %d\n",
                sort!(collect(tagcount)), mech_onlya_mixed)
    @printf("   mechanisms with an RE chemistry step, by class: %s; with an RE fused binding inside a mixed group: %d\n",
            sort!(collect(re_chem_mech)), re_fused_in_mixed_mech)
end

println("\n== distinct mixed-group shapes (flag, tag, steps) per population, most frequent first")
for name in POPS
    sh = shapes[name]
    isempty(sh) && continue
    println("-- $name: $(length(sh)) shapes")
    for (sig, c) in first(sort!(collect(sh); by = x -> (-x[2], x[1])), 12)
        @printf("   %4d  %s\n", c, sig)
    end
end
println("\n== first examples")
for k in sort!(collect(keys(examples))), e in examples[k]
    println("   ", k, ": ", e)
end
