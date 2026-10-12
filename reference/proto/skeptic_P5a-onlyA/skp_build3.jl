# ABOUTME: Skeptic P5a: builds the multi-member :OnlyA chemistry-group cases (third batch) with EnzymeRates
# ABOUTME: and exports Pop.PMech records with the guard verdict, as skp_build.jl does.
using EnzymeRates, Random, Serialization, Printf, Statistics
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "common", "export", "export_lib.jl"))

const SB = Module(:SkpCases)
Core.eval(SB, :(using EnzymeRates))

head(subs, prods; inh = String[], n = 2) = """
    substrates: $(join(subs, ", "))
    products: $(join(prods, ", "))
""" * (isempty(inh) ? "" : "    catalytic_inhibitors: $(join(inh, ", "))\n") *
"    catalytic_multiplicity: $n\n"
dsl(subs, prods, lines; inh = String[], n = 2) = "@allosteric_mechanism begin\n" *
    head(subs, prods; inh, n) * "    catalytic_steps: begin\n        " *
    join(lines, "\n        ") * "\n    end\nend"

C = Tuple{String, String, String}[]
add!(id, what, txt) = push!(C, (id, what, txt))
# AD14: a two-member :OnlyA chemistry group whose members need different ratio exponents.
add!("AD14_chem_group_split", "OnlyA chemistry in two contexts (free and I-bound), S binds E(I) OnlyA",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "(E(S) <--> E(P), E(S, I::Inh) <--> E(P, I::Inh)) :: OnlyA",
        "E(I::Inh) + S ⇌ E(S, I::Inh) :: OnlyA",
        "E(I::Inh) + P ⇌ E(P, I::Inh) :: EqualAI",
        "E + I::Inh ⇌ E(I::Inh) :: EqualAI"]; inh = ["I"]))
add!("AD14b_chem_group_split_SS", "AD14 with the I-context S binding SS",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "(E(S) <--> E(P), E(S, I::Inh) <--> E(P, I::Inh)) :: OnlyA",
        "E(I::Inh) + S <--> E(S, I::Inh) :: OnlyA",
        "E(I::Inh) + P ⇌ E(P, I::Inh) :: EqualAI",
        "E + I::Inh ⇌ E(I::Inh) :: EqualAI"]; inh = ["I"]))
add!("AD14c_chem_group_ok", "control: same groups, the I-context S binding EqualAI",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "(E(S) <--> E(P), E(S, I::Inh) <--> E(P, I::Inh)) :: OnlyA",
        "E(I::Inh) + S ⇌ E(S, I::Inh) :: EqualAI",
        "E(I::Inh) + P ⇌ E(P, I::Inh) :: EqualAI",
        "E + I::Inh ⇌ E(I::Inh) :: EqualAI"]; inh = ["I"]))

out = Pop.PMech[]
for (i, (id, what, txt)) in enumerate(C)
    try
        ex = Meta.parse(txt)
        block = ex.args[end]
        mech_expr, _, _ = ER._parse_mechanism_body(block, true)
        inner = mech_expr.args[2]
        rxn = Core.eval(SB, inner.args[2]); groups = Core.eval(SB, inner.args[3])
        tags = Core.eval(SB, inner.args[4]); n = inner.args[5]; sites = Core.eval(SB, inner.args[6])
        gs, perm = ER._canonical_groups(rxn, groups)
        ctags = tags[perm]
        viol = ER._onlya_haldane_violation(rxn, gs, ctags)
        if viol === nothing
            am = ER.AllostericMechanism(rxn, groups, tags, n, sites)
            ER.cat_allo_states(am) == ctags || error("tag permutation mismatch")
            p = to_pmech(am, id, "skp:" * what; seed = 7000 + i, do_bench = false)
            push!(out, p)
            @printf("%-30s ADMITTED ok=%s fitted=%d stage_errors=%s\n", id, p.ok,
                    length(p.fitted_params), isempty(p.stage_errors) ? "-" : string(p.stage_errors))
        else
            tw = ER.AllostericMechanism(rxn, groups, fill(:NonequalAI, length(groups)), n, sites)
            p = to_pmech(tw, id, "skp:" * what; seed = 7000 + i, do_bench = false, do_kcat = false)
            d = Dict(f => getfield(p, f) for f in fieldnames(Pop.PMech))
            d[:ok] = false; d[:error] = "rejected: " * viol; d[:group_tags] = ctags
            push!(out, Pop.PMech(; d...))
            @printf("%-30s REJECTED %s\n", id, first(viol, 150))
        end
    catch e
        @printf("%-30s BUILD ERROR %s\n", id, first(sprint(showerror, e isa LoadError ? e.error : e), 300))
    end
end
serialize(joinpath(@__DIR__, "data", "adv_cases3.jls"), out)
println("wrote $(length(out)) cases")
for m in out
    println("── ", m.id, "  tags=", m.group_tags, " n=", m.cat_multiplicity, " ok=", m.ok)
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("     g", g[j], " [", m.group_tags[g[j]], "] ", Pop.step_text(m, s),
                isempty(m.inactive_kept) ? "" : (m.inactive_kept[j] ? "" : "   (pruned in I)"))
    end
    println("   species: ", [sp.name for sp in m.species], "  free: ",
            [sp.name for sp in m.species if Pop.is_free(sp)])
    m.ok && println("   v samples: ", [s.v for s in m.samples])
end
