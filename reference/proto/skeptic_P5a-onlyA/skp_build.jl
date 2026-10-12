# ABOUTME: Skeptic P5a: builds adversarial allosteric :OnlyA mechanisms with EnzymeRates, records the
# ABOUTME: guard verdict, and exports Pop.PMech records (twins with restored tags for rejected ones).
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

# AD1: a mixed :OnlyA group (fused binding + plain binding of B) whose plain member sits
# on a binding square of :EqualAI steps.
add!("AD1_mixed_square_SS", "mixed OnlyA group, plain member on an EqualAI square (SS)",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A ⇌ E(A) :: EqualAI",
        "(E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q)) :: OnlyA",
        "E(Q) + P ⇌ E(P, Q) :: EqualAI",
        "E + Q ⇌ E(Q) :: EqualAI",
        "E + B ⇌ E(B) :: EqualAI",
        "E(B) + Q ⇌ E(B, Q) :: EqualAI"]))
add!("AD1_mixed_square_RE", "mixed OnlyA group, plain member on an EqualAI square (RE)",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A ⇌ E(A) :: EqualAI",
        "(E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q)) :: OnlyA",
        "E(Q) + P <--> E(P, Q) :: EqualAI",
        "E + Q ⇌ E(Q) :: EqualAI",
        "E + B ⇌ E(B) :: EqualAI",
        "E(B) + Q ⇌ E(B, Q) :: EqualAI"]))
add!("AD1_mixed_square_SS_n1", "AD1 SS at multiplicity 1",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A ⇌ E(A) :: EqualAI",
        "(E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q)) :: OnlyA",
        "E(Q) + P ⇌ E(P, Q) :: EqualAI",
        "E + Q ⇌ E(Q) :: EqualAI",
        "E + B ⇌ E(B) :: EqualAI",
        "E(B) + Q ⇌ E(B, Q) :: EqualAI"]; n = 1))
# AD2: two free forms (conformations E and F)
add!("AD2a_iso_chemOnlyA", "iso uni-uni, two free forms, chemistry OnlyA",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: OnlyA",
        "F + P ⇌ E(P) :: EqualAI",
        "F <--> E :: EqualAI"]))
add!("AD2a_iso_chemNonequal", "iso uni-uni, two free forms, no OnlyA (V1-type check)",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: NonequalAI",
        "F + P ⇌ E(P) :: EqualAI",
        "F <--> E :: NonequalAI"]))
add!("AD2b_iso_bindOnlyA", "iso uni-uni, two free forms, P binding to F and iso OnlyA",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: EqualAI",
        "F + P ⇌ E(P) :: OnlyA",
        "F <--> E :: OnlyA"]))
add!("AD2c_two_free_disconnected", "second free form F joined only by an OnlyA iso, F binds an inhibitor",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "E <--> F :: OnlyA",
        "F + I::Inh ⇌ F(I::Inh) :: EqualAI"]; inh = ["I"]))
add!("AD2d_iso_isoOnlyA", "iso uni-uni, two free forms, only the free-form iso OnlyA",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: EqualAI",
        "F + P ⇌ E(P) :: EqualAI",
        "F <--> E :: OnlyA"]))
# AD3: inhibitors and ::Inh copies
add!("AD3a_inhibitor_OnlyA", "competitive inhibitor on free E tagged OnlyA",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: NonequalAI",
        "E + P ⇌ E(P) :: EqualAI",
        "E + I::Inh ⇌ E(I::Inh) :: OnlyA"]; inh = ["I"]))
add!("AD3b_copy_square_tied", "substrate-inhibition copy square, both copy bindings one OnlyA group",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "(E + S::Inh ⇌ E(S::Inh), E(S) + S::Inh ⇌ E(S, S::Inh)) :: OnlyA",
        "E(S::Inh) + S ⇌ E(S, S::Inh) :: EqualAI"]; inh = ["S"]))
add!("AD3c_copy_square_onesided", "copy square with one copy binding OnlyA (Wegscheider-only violation)",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: OnlyA",
        "E + P ⇌ E(P) :: EqualAI",
        "E + S::Inh ⇌ E(S::Inh) :: OnlyA",
        "E(S) + S::Inh ⇌ E(S, S::Inh) :: EqualAI",
        "E(S::Inh) + S ⇌ E(S, S::Inh) :: EqualAI"]; inh = ["S"]))
add!("AD3d_copy_square_SS_split", "copy square, the two copy bindings separate OnlyA SS groups",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "E + S::Inh <--> E(S::Inh) :: OnlyA",
        "E(S) + S::Inh <--> E(S, S::Inh) :: OnlyA",
        "E(S::Inh) + S ⇌ E(S, S::Inh) :: EqualAI"]; inh = ["S"]))
# AD4: several metabolites on one side
add!("AD4_multi_release", "uni-bi with a two-product release step tagged OnlyA",
    dsl(["S"], ["P", "Q"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P, Q) :: EqualAI",
        "E(P, Q) <--> E + P + Q :: OnlyA"]))
add!("AD4b_multi_release_chem", "uni-bi, chemistry OnlyA, two-product release RE",
    dsl(["S"], ["P", "Q"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) <--> E(P, Q) :: OnlyA",
        "E(P, Q) ⇌ E + P + Q :: EqualAI"]))
# AD5: a group whose members bind different conformations
add!("AD5_conf_selection", "S binds E and F in one OnlyA group; conformational selection",
    dsl(["S"], ["P"], [
        "E <--> F :: EqualAI",
        "(E + S ⇌ E(S), F + S ⇌ F(S)) :: OnlyA",
        "E(S) <--> F(S) :: EqualAI",
        "F(S) <--> F(P) :: OnlyA",
        "F + P ⇌ F(P) :: EqualAI"]))
# AD6: parallel steps between the same two forms
add!("AD6_parallel", "parallel RE and SS isomerization between E(S) and E(P)",
    dsl(["S"], ["P"], [
        "E + S ⇌ E(S) :: EqualAI",
        "E(S) ⇌ E(P) :: EqualAI",
        "E(S) <--> E(P) :: OnlyA",
        "E + P <--> E(P) :: EqualAI"]))
# AD7: thermodynamics pins a constant to 1 (handwritten/24 made allosteric)
add!("AD7a_pin_bindOnlyA", "pins-to-1 iso: inhibitor group OnlyA",
    dsl(["A"], ["P"], [
        "E + A <--> E(A) :: EqualAI",
        "(E + I::Inh ⇌ E(I::Inh), F + I::Inh ⇌ E(I::Inh)) :: OnlyA",
        "E + P <--> E(P) :: EqualAI",
        "E ⇌ F :: EqualAI",
        "E(A) <--> E(P) :: EqualAI"]; inh = ["I"]))
add!("AD7b_pin_isoOnlyA", "pins-to-1 iso: E ⇌ F OnlyA",
    dsl(["A"], ["P"], [
        "E + A <--> E(A) :: EqualAI",
        "(E + I::Inh ⇌ E(I::Inh), F + I::Inh ⇌ E(I::Inh)) :: EqualAI",
        "E + P <--> E(P) :: EqualAI",
        "E ⇌ F :: OnlyA",
        "E(A) <--> E(P) :: EqualAI"]; inh = ["I"]))
add!("AD7c_pin_chemOnlyA", "pins-to-1 iso: chemistry OnlyA, F a second free form",
    dsl(["A"], ["P"], [
        "E + A <--> E(A) :: EqualAI",
        "(E + I::Inh ⇌ E(I::Inh), F + I::Inh ⇌ E(I::Inh)) :: EqualAI",
        "E + P <--> E(P) :: EqualAI",
        "E ⇌ F :: EqualAI",
        "E(A) <--> E(P) :: OnlyA"]; inh = ["I"]))
# AD8: Theorell-Chance OnlyA with a dead end
add!("AD8_tc_deadend", "Theorell-Chance step OnlyA, dead-end B on E(Q)",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A ⇌ E(A) :: EqualAI",
        "E(A) + B <--> E(Q) + P :: OnlyA",
        "E + Q ⇌ E(Q) :: EqualAI",
        "E(Q) + B ⇌ E(B, Q) :: EqualAI"]))
# AD9: RE chemistry OnlyA inside an RE segment with SS releases, plus an SS dead end
add!("AD9_re_chem_ss_deadend", "RE chemistry OnlyA, SS releases, SS dead-end OnlyA",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A ⇌ E(A) :: EqualAI",
        "E(A) + B ⇌ E(A, B) :: EqualAI",
        "E(A, B) ⇌ E(P, Q) :: OnlyA",
        "E(Q) + P <--> E(P, Q) :: EqualAI",
        "E + Q <--> E(Q) :: EqualAI",
        "E(A) + Q <--> E(A, Q) :: OnlyA"]))
# AD10: a binding square with two separate OnlyA SS bindings and OnlyA chemistry
add!("AD10_square_two_ss", "random-order A: two OnlyA SS A-bindings, chemistry OnlyA",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A <--> E(A) :: OnlyA",
        "E + B ⇌ E(B) :: EqualAI",
        "E(A) + B ⇌ E(A, B) :: EqualAI",
        "E(B) + A <--> E(A, B) :: OnlyA",
        "E(A, B) <--> E(P, Q) :: OnlyA",
        "E(Q) + P ⇌ E(P, Q) :: EqualAI",
        "E + Q ⇌ E(Q) :: EqualAI"]))
# AD11: mixed group whose plain member is a dead end AND a second plain member on a square
add!("AD11_mixed_two_plain", "mixed OnlyA group: fused B binding + dead-end plain B + plain B on a square",
    dsl(["A", "B"], ["P", "Q"], [
        "E + A ⇌ E(A) :: EqualAI",
        "(E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q), E(P) + B <--> E(B, P)) :: OnlyA",
        "E(Q) + P ⇌ E(P, Q) :: EqualAI",
        "E + Q ⇌ E(Q) :: EqualAI",
        "E + P ⇌ E(P) :: EqualAI",
        "E + B ⇌ E(B) :: EqualAI",
        "E(B) + Q ⇌ E(B, Q) :: EqualAI"]))

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
            p = to_pmech(am, id, "skp:" * what; seed = 5000 + i, do_bench = false)
            push!(out, p)
            @printf("%-30s ADMITTED ok=%s fitted=%d stage_errors=%s\n", id, p.ok,
                    length(p.fitted_params), isempty(p.stage_errors) ? "-" : string(p.stage_errors))
        else
            tw = ER.AllostericMechanism(rxn, groups, fill(:NonequalAI, length(groups)), n, sites)
            p = to_pmech(tw, id, "skp:" * what; seed = 5000 + i, do_bench = false, do_kcat = false)
            d = Dict(f => getfield(p, f) for f in fieldnames(Pop.PMech))
            d[:ok] = false; d[:error] = "rejected: " * viol; d[:group_tags] = ctags
            push!(out, Pop.PMech(; d...))
            @printf("%-30s REJECTED %s\n", id, first(viol, 150))
        end
    catch e
        @printf("%-30s BUILD ERROR %s\n", id, first(sprint(showerror, e isa LoadError ? e.error : e), 300))
    end
end
serialize(joinpath(@__DIR__, "data", "adv_cases.jls"), out)
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
