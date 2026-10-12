# ABOUTME: Builds adversarial plain mechanisms with the package DSL and exports them as Pop.PMech records.
# ABOUTME: Needs EnzymeRates (oracle worktree); writes skeptic_P1-engine/data/adv_pkg.jls and a log of package errors.
using EnzymeRates, Random, Serialization, Printf, Statistics
const C = joinpath(@__DIR__, "..", "common")
include(joinpath(C, "Population.jl"))
include(joinpath(C, "export", "export_lib.jl"))

# Each case: (id, DSL source). Parsed and evaluated one by one so a DSL or constructor
# error is recorded instead of aborting the run.
cases = Pair{String, String}[]
mech(body; subs = "S", prods = "P", regs = "") =
    "@enzyme_mechanism begin\n substrates: $subs\n products: $prods\n" *
    (isempty(regs) ? "" : " regulators: $regs\n") * " steps: begin\n$body\n end\nend"

push!(cases, "A1_iso_parallel_ss" => mech("""
    E + S <--> F + P
    F <--> E"""))
push!(cases, "A2_iso_parallel_re_iso" => mech("""
    E + S <--> F + P
    F ⇌ E"""))
push!(cases, "A3_two_reactive_conformations_ss_bridge" => mech("""
    E + S <--> E(S)
    E(S) <--> E + P
    F + S <--> F(S)
    F(S) <--> F + P
    E <--> F"""))
push!(cases, "A3re_two_reactive_conformations_re_bindings" => mech("""
    E + S ⇌ E(S)
    E(S) <--> E + P
    F + S ⇌ F(S)
    F(S) <--> F + P
    E <--> F"""))
push!(cases, "A3g_bridge_grouped_bindings" => mech("""
    (E + S <--> E(S), F + S <--> F(S))
    E(S) <--> E + P
    F(S) <--> F + P
    E <--> F"""))
push!(cases, "A4_single_segment_core_ss_leaf" => mech("""
    E + S ⇌ E(S)
    E(S) <--> E(P)
    E(P) ⇌ E + P
    E + I <--> E(I)"""; regs = "I"))
push!(cases, "A4b_single_segment_core_ss_chain" => mech("""
    E + S ⇌ E(S)
    E(S) <--> E(P)
    E(P) ⇌ E + P
    E + I <--> E(I)
    E(I) <--> F(I)"""; regs = "I"))
push!(cases, "A5_cyclic_dead_end_ss_square" => mech("""
    E + S <--> E(S)
    E(S) <--> E + P
    E + I <--> E(I)
    E + J <--> E(J)
    E(I) + J <--> E(I, J)
    E(J) + I <--> E(I, J)"""; regs = "I, J"))
push!(cases, "A6_free_forms_one_re_segment_D_E" => mech("""
    E + S <--> E(S)
    E(S) <--> E + P
    D ⇌ E"""))
push!(cases, "A6b_free_forms_ss_iso_D_E" => mech("""
    E + S ⇌ E(S)
    E(S) <--> E + P
    D <--> E"""))
push!(cases, "A6c_free_forms_one_re_segment_written_D_first" => mech("""
    D ⇌ E
    D + S <--> D(S)
    D(S) <--> E + P"""))
push!(cases, "A7_substrate_copy_in_reactive_cycles" => mech("""
    E + S ⇌ E(S)
    E(S) <--> E + P
    E + S::Inh ⇌ E(S::Inh)
    E(S::Inh) + S ⇌ E(S, S::Inh)
    E(S) + S::Inh ⇌ E(S, S::Inh)
    E(S, S::Inh) <--> E(S::Inh) + P"""; regs = "S"))
push!(cases, "A8_fused_two_metabolites_re" => mech("""
    E + A + B ⇌ E(A, B)
    E(A, B) <--> E + P + Q"""; subs = "A, B", prods = "P, Q"))
push!(cases, "A8ss_fused_two_metabolites_ss" => mech("""
    E + A + B <--> E(A, B)
    E(A, B) <--> E(Q) + P
    E(Q) <--> E + Q"""; subs = "A, B", prods = "P, Q"))
push!(cases, "A9_slip_cycle" => mech("""
    E + A <--> E(A)
    E(A) + B <--> E(A, B)
    E(A, B) <--> E(P, Q)
    E(P, Q) <--> E(Q) + P
    E(Q) <--> E + Q
    E(A) <--> E(P)
    E(P) <--> E + P"""; subs = "A, B", prods = "P, Q"))
push!(cases, "A10_cycle_turns_twice" => mech("""
    E + S <--> E(S)
    E(S) <--> F + P
    F + S <--> F(S)
    F(S) <--> E + P"""))
push!(cases, "A10g_cycle_turns_twice_grouped_bindings" => mech("""
    (E + S <--> E(S), F + S <--> F(S))
    E(S) <--> F + P
    F(S) <--> E + P"""))
push!(cases, "A10gg_cycle_turns_twice_grouped_all" => mech("""
    (E + S <--> E(S), F + S <--> F(S))
    (E(S) <--> F + P, F(S) <--> E + P)"""))
push!(cases, "A11_parallel_rank_12" => mech(join([
    "    E + S <--> E(S)", "    E(S) <--> E + P",
    ["    E + I$k ⇌ E(I$k)\n    E(S) + I$k ⇌ E(S, I$k)\n    E(I$k) + S <--> E(S, I$k)"
     for k in 1:11]...], "\n"); regs = join(["I$k" for k in 1:11], ", ")))
push!(cases, "A12_parallel_rank_21" => mech(join([
    "    E + S <--> E(S)", "    E(S) <--> E + P",
    ["    E + I$k ⇌ E(I$k)\n    E(S) + I$k ⇌ E(S, I$k)\n    E(I$k) + S <--> E(S, I$k)"
     for k in 1:20]...], "\n"); regs = join(["I$k" for k in 1:20], ", ")))
push!(cases, "A13_group_across_conformations_tc" => mech("""
    E + S <--> E(S)
    (E(S) <--> E(P), F(S) <--> F(P))
    E(P) <--> E + P
    F + S <--> F(S)
    F(P) <--> F + P
    E <--> F"""))
push!(cases, "A14_uniuni_two_intermediates_all_ss" => mech("""
    E + S <--> E(S)
    E(S) <--> E(P)
    E(P) <--> E + P"""))
push!(cases, "A16_nonreactive_square_between_reactive_conformations" => mech("""
    E + S <--> E(S)
    E(S) <--> E + P
    F + S <--> F(S)
    F(S) <--> F + P
    E <--> G
    G <--> F
    E <--> H
    H <--> F"""))
push!(cases, "A17_dead_end_on_ping_pong_residual_form" => mech("""
    E + A <--> E(A)
    E(A) <--> E(; residual = A - P) + P
    E(; residual = A - P) + B <--> E(B; residual = A - P)
    E(B; residual = A - P) <--> E + Q
    E(; residual = A - P) + I <--> E(I; residual = A - P)"""; subs = "A, B", prods = "P, Q", regs = "I"))
push!(cases, "A15_dead_end_branch_with_re_inside" => mech("""
    E + S <--> E(S)
    E(S) <--> E + P
    E + I <--> E(I)
    E(I) + J ⇌ E(I, J)"""; regs = "I, J"))

# Second set (ARGS[1] == "extra"): the first free form in species order is not its segment's root.
const EXTRA = Pair{String, String}[
    "A18_free_form_E_seen_before_root_D" => mech("""
    E + S <--> C(S)
    C(S) + B <--> E + P
    D ⇌ E"""; subs = "B, S", prods = "P"),
    "A18b_free_form_E_seen_before_root_D_ss_iso" => mech("""
    E + S ⇌ C(S)
    C(S) + B <--> E + P
    D ⇌ E"""; subs = "B, S", prods = "P"),
]
const SET = isempty(ARGS) ? "main" : ARGS[1]
SET == "extra" && (cases = EXTRA)

out = Pop.PMech[]
errs = String[]
for (k, (id, src)) in enumerate(cases)
    t0 = time()
    m = try
        x = Core.eval(Main, Meta.parse(src))
        x isa EnzymeRates.AbstractEnzymeMechanism ? ER._concrete(x) : x
    catch e
        push!(errs, "$id: construction: " * first(sprint(showerror, e), 400))
        println(errs[end]); continue
    end
    pm = try
        to_pmech(m, "adv_pkg/$k", id; seed = 777_000 + k)
    catch e
        push!(errs, "$id: to_pmech: " * first(sprint(showerror, e), 400))
        println(errs[end]); continue
    end
    push!(out, pm)
    @printf("%-48s ok=%d %.2fs  fitted=%d  errors: %s\n", id, pm.ok, time() - t0,
            length(pm.fitted_params),
            join(("$(a): $(first(b, 200))" for (a, b) in pm.stage_errors), " | "))
end
mkpath(joinpath(@__DIR__, "data"))
serialize(joinpath(@__DIR__, "data", SET == "extra" ? "adv_pkg2.jls" : "adv_pkg.jls"), out)
println("wrote $(length(out)) mechanisms; construction/export errors: $(length(errs))")
