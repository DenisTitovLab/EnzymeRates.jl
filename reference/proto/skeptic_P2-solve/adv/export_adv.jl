# ABOUTME: Skeptic of probe P2: builds adversarial plain mechanisms with EnzymeRates and exports
# ABOUTME: them as Pop.PMech (package outputs included); also re-checks one same-model split.
using EnzymeRates, Random, Serialization, Printf, Statistics
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const OUT = joinpath(PROTO, "skeptic_P2-solve", "data", "adv_pkg.jls")
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(PROTO, "common", "export", "export_lib.jl"))

mech(subs, prods, regs, body) = """
@enzyme_mechanism begin
    substrates: $subs
    products: $prods
    $(isempty(regs) ? "" : "regulators: " * regs)
    steps: begin
$body
    end
end"""

const DEFS = [
    # ROADMAP known issue 3 exactly, and the same step at RE
    "exchange_bridge" => mech("S", "P", "I, J", """
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I ⇌ E(I)
        E(I) + J <--> E(J) + I"""),
    "exchange_bridge_RE" => mech("S", "P", "I, J", """
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I ⇌ E(I)
        E(I) + J ⇌ E(J) + I"""),
    "exchange_in_wegscheider_cycle" => mech("S", "P", "I, J", """
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I ⇌ E(I)
        E(I) + J <--> E(J) + I
        E + J ⇌ E(J)"""),
    "dead_end_ss_binding" => mech("S", "P", "I", """
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
        E + I <--> E(I)"""),
    # dead-end SS chain and SS triangle (a non-reactive block)
    "dead_end_ss_chain" => mech("S", "P", "I", """
        E + S <--> E(S)
        E(S) <--> E(P)
        E(P) <--> E + P
        E + I <--> E(I)
        E(I) <--> F(I)"""),
    "dead_end_ss_triangle" => mech("S", "P", "I", """
        E + S <--> E(S)
        E(S) <--> E(P)
        E(P) <--> E + P
        E + I <--> E(I)
        E(I) <--> F(I)
        E + I <--> F(I)"""),
    # an SS step inside an RE segment, parallel to the RE path (zero-change cycle)
    "ss_in_re_segment" => mech("S", "P", "", """
        E + S ⇌ E(S)
        E(S) ⇌ F(S)
        E + S <--> F(S)
        F(S) <--> F(P)
        F(P) ⇌ F + P
        F <--> E"""),
    # a flux-free SS binding grouped with the productive SS binding
    "mixed_flux_free_group" => mech("S", "P", "J", """
        (E + S <--> E(S), E(J) + S <--> E(J, S))
        E(S) <--> E(P)
        E(P) <--> E + P
        E + J ⇌ E(J)"""),
    # chemistry isomerizations that change conformation in opposite directions, grouped
    "iso_group_opposite_chem" => mech("S", "P", "", """
        E + S <--> E(S)
        F + S <--> F(S)
        (E(S) <--> F(P), F(S) <--> E(P))
        F(P) <--> F + P
        E(P) <--> E + P"""),
    # a chemistry isomerization (Tier 1) grouped with a conformational one (Tier 3)
    "iso_group_mixed_tier" => mech("S", "P", "I", """
        F + S <--> F(S)
        (F(S) <--> E(P), E(I) <--> F(I))
        E(P) <--> E + P
        E <--> F
        E + I ⇌ E(I)
        F + I ⇌ F(I)"""),
    # bindings in two conformations grouped; conformational isos grouped
    "conf_group_bindings" => mech("S", "P", "", """
        (E + S <--> E(S), F + S <--> F(S))
        E <--> F
        E(S) <--> F(S)
        E(S) <--> E(P)
        F(S) <--> F(P)
        E(P) <--> E + P
        F(P) <--> F + P"""),
    "conf_group_isos" => mech("S", "P", "", """
        E + S <--> E(S)
        F + S <--> F(S)
        (E <--> F, E(S) <--> F(S))
        E(S) <--> E(P)
        E(P) <--> E + P"""),
    # a group of two Theorell-Chance-like steps between two free forms: k_r = k_f / Keq
    "tc_pair_group" => mech("S", "P", "I", """
        (E + S <--> F + P, F + S <--> E + P)
        E + I <--> E(I)"""),
    # thermodynamics pins K_TC to Keq through a grouped dead-end binding
    "keq_pin" => mech("S", "P", "I", """
        E + S ⇌ F + P
        E(S) <--> E(P)
        E + S ⇌ E(S)
        E(P) ⇌ E + P
        (E + I <--> E(I), F + I <--> E(I))"""),
    # a group whose class enters a Wegscheider cycle twice: half-integer exponents
    "half_exponent" => mech("A", "Z", "", """
        E + A ⇌ E(A)
        E(A) + A ⇌ E(A, A)
        (E(A, A) <--> E(A, Z), E(A, Z) <--> D(Z, Z))
        E(A, A) <--> D(Z, Z)
        D(Z, Z) ⇌ D(Z) + Z
        D(Z) ⇌ D + Z
        D <--> E"""),
    # an inhibitor copy of the substrate that isomerizes into the productive complex
    "inh_copy_iso_cycle" => mech("A", "P", "A", """
        E + A ⇌ E(A)
        E + A::Inh ⇌ E(A::Inh)
        E(A::Inh) <--> E(A)
        E(A) <--> E(P)
        E(P) ⇌ E + P"""),
    # an SS exchange of the substrate for its inhibitor copy (zero net change by name)
    "inh_copy_exchange" => mech("A", "P", "A", """
        E + A ⇌ E(A)
        E + A::Inh ⇌ E(A::Inh)
        E(A) + A::Inh <--> E(A::Inh) + A
        E(A) <--> E(P)
        E(P) ⇌ E + P"""),
    # Theorell-Chance step in a Wegscheider cycle with the ternary route
    "tc_ternary_wegscheider" => mech("A, B", "P, Q", "", """
        E + A ⇌ E(A)
        E(A) + B <--> E(Q) + P
        E(A) + B ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) ⇌ E(Q) + P
        E(Q) ⇌ E + Q"""),
    # ping-pong with SS dead ends at the residual form and at E
    "pingpong_dead_ends" => mech("A, B", "P, Q", "I, B", """
        E + A <--> E(A)
        E(A) <--> E(; residual = A - P) + P
        E(; residual = A - P) + B <--> E(B; residual = A - P)
        E(B; residual = A - P) <--> E + Q
        E(; residual = A - P) + I <--> E(I; residual = A - P)
        E + B::Inh <--> E(B::Inh)"""),
    # iso ordered bi-bi: free forms E and F in different SS segments, grouped dead ends
    "iso_bibi_two_free" => mech("A, B", "P, Q", "I", """
        E + A <--> E(A)
        E(A) + B <--> E(A, B)
        E(A, B) <--> F(P, Q)
        F(P, Q) <--> F(Q) + P
        F(Q) <--> F + Q
        F <--> E
        (E + I ⇌ E(I), F + I ⇌ F(I))"""),
    # two metabolites on one side, in a Wegscheider cycle with the sequential route
    "multi_metabolite_side" => mech("A, B", "P, Q", "", """
        E + A + B ⇌ E(A, B)
        E + A ⇌ E(A)
        E(A) + B ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) ⇌ E + P + Q"""),
    # mirror-symmetric routes in E and F, K(E -> F) pinned to 1: the bridge E(S) <--> F(S)
    # and E <--> F carry no flux by symmetry though neither is structurally flux-free
    "symmetric_routes_pinned" => mech("S", "P", "I", """
        (E + S <--> E(S), F + S <--> F(S))
        (E(S) <--> E(P), F(S) <--> F(P))
        (E(P) <--> E + P, F(P) <--> F + P)
        E(S) <--> F(S)
        E <--> F
        (E + I ⇌ E(I), F + I ⇌ E(I))"""),
    # package-rejected inputs (recorded as rejections)
    "parallel_identical_steps" => mech("S", "P", "", """
        E + S ⇌ E(S)
        E(S) <--> E(P)
        (E(S) <--> E(P))
        E(P) ⇌ E + P"""),
    "group_different_conformational_changes" => mech("S", "P", "", """
        E + S <--> E(S)
        E(S) <--> E(P)
        E(P) <--> E + P
        (E <--> F, F <--> G)
        G <--> E"""),
]

out = Pop.PMech[]
for (i, (id, src)) in enumerate(DEFS)
    t = time()
    em = try
        Core.eval(Main, Meta.parse(src))
    catch e
        println(rpad(id, 40), "REJECTED by the package: ", first(sprint(showerror, e), 300))
        continue
    end
    pm = try
        to_pmech(ER._concrete(em), "adv/$i", id; seed = 900_000 + i, do_kcat = false,
                 do_bench = false)
    catch e
        println(rpad(id, 40), "EXPORT FAILED: ", first(sprint(showerror, e), 300))
        continue
    end
    push!(out, pm)
    @printf("%-40s ok=%d fitted=%d ipc=%d %.1fs %s\n", id, pm.ok, length(pm.fitted_params),
            pm.independent_param_count, time() - t,
            isempty(pm.stage_errors) ? "" : string(pm.stage_errors))
    println("    fitted: ", pm.fitted_params)
    for d in pm.dependents; println("    dep ", d); end
    isempty(pm.pass2_rename) || println("    pass2: ", pm.pass2_rename)
    flush(stdout)
end
Pop.save_population("adv", out; path = OUT)
println("wrote ", OUT, " (", length(out), " mechanisms)")
