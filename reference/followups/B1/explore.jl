# ABOUTME: First pass of follow-up B1: class counts and configuration flags per population and reading,
# ABOUTME: plus a check of the plain-binding test against Thermo's kind; stdlib only.
include(joinpath(@__DIR__, "b1_lib.jl"))
include(joinpath(S, "proto/common/Thermo.jl"))
for name in POPS
    pop = load_allo(name)
    println("== $name: ", length(pop), " allosteric, ", count(m -> :OnlyA in m.group_tags, pop), " with :OnlyA")
    # plain test vs Thermo kind 1 (no conformation check)
    nd = 0
    for m in pop
        cs = Thermo.canonical_steps(m)
        for (j, s) in enumerate(m.steps)
            ((step_class(m, s) === :plain) != (cs[j].kind == 1)) && (nd += 1)
        end
    end
    cens = Dict{Symbol, Int}()
    for m in pop, (j, s) in enumerate(m.steps)
        blocked_flags(m)[j] || continue
        k = Symbol(step_class(m, s), s.re ? "_RE" : "_SS")
        cens[k] = get(cens, k, 0) + 1
    end
    println("   plain-test mismatches vs Thermo kind: $nd; blocked step classes: ", sort!(collect(cens); by = first))
    for rd in (:D, :ii, :iu, :i)
        as = [analyze(m, rd) for m in pop]
        cnt = Dict{String, Int}()
        for a in as; cnt[a.cls] = get(cnt, a.cls, 0) + 1; end
        @printf("   %-3s %s\n", rd, join(("$k=$v" for (k, v) in sort!(collect(cnt); by = first)), " "))
        @printf("       nfree≠1 %d | unbalanced %d | directed cycle %d | aff chords %d | kin chords %d | forced-ratio self-loops %d | cyc through main %d | hanging cyc %d | ambiguous signs %d | lit_seq %d | lit_emptied&!no_island %d\n",
                count(a -> a.nfree != 1, as), count(a -> a.unbalanced, as), count(a -> a.directed_cycle, as),
                count(a -> !isempty(a.aff_chords), as), count(a -> !isempty(a.kchords), as),
                count(a -> !isempty(a.forced), as), count(a -> a.cyc_main, as), count(a -> a.cyc_hang, as),
                count(a -> 2 in a.isl_sign, as), count(a -> a.lit_seq, as),
                count(a -> a.lit_emptied && a.nisl > 0, as))
        for c in ("phantom", "undefined", "sequestering", "emptied", "no_island", "ambiguous", "infeasible", "free_forms")
            sel = filter(a -> a.cls == c, as)
            isempty(sel) && continue
            @printf("       %-12s kinetic doors 0:%d 1:%d ≥2:%d | islands 1:%d ≥2:%d | affinity doors 0:%d ≥1:%d\n", c,
                    count(a -> a.nkdoor == 0, sel), count(a -> a.nkdoor == 1, sel), count(a -> a.nkdoor >= 2, sel),
                    count(a -> a.nisl == 1, sel), count(a -> a.nisl >= 2, sel),
                    count(a -> a.nadoor == 0, sel), count(a -> a.nadoor >= 1, sel))
        end
    end
end
