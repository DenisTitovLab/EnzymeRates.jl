# ABOUTME: Self-contained reproducers for the two derivation/enumeration bugs behind the LDH futile
# ABOUTME: cycle: A (non-idempotent _canonical_mechanism) and B2 (Wegscheider tie enforced in child, not
# ABOUTME: parent → a split that yields the identical rate function). Run: julia --project=. <this file>
using EnzymeRates, Random
recon(s) = EnzymeRates.AllostericMechanism(Core.eval(EnzymeRates, Meta.parse(s))())
C1 = EnzymeRates._canonical_mechanism
ehash(em) = string(EnzymeRates._rate_eq_dedup_key(EnzymeRates.rate_equation_string(em)), base=16, pad=16)

# Parent mechanisms from docs/ldh_hpc_results/2026_07_09_results_2 (pre-#64 run), as compiled-type Sig
# strings (the `mechanism_type` CSV column). eq_hash of A parent = 5b270db6065fc543, B2 parent = 7f09e180b56580aa.
const A_SIG  = "AllostericEnzymeMechanism{EnzymeMechanism{(((((:Product, :Lactate), ((:C, 3), (:H, 6), (:O, 3))), ((:Product, :NAD), ((:C, 21), (:H, 27), (:N, 7), (:O, 14), (:P, 2))), ((:Substrate, :NADH), ((:C, 21), (:H, 29), (:N, 7), (:O, 14), (:P, 2))), ((:Substrate, :Pyruvate), ((:C, 3), (:H, 4), (:O, 3)))), (((:CompetitiveInhibitor, :Lactate), (4,)), ((:CompetitiveInhibitor, :NAD), (4,)), ((:CompetitiveInhibitor, :NADH), (4,)), ((:CompetitiveInhibitor, :Pyruvate), (4,))), (4,)), (((((), :E, ((), ())), (((:Product, :Lactate),), :E, ((), ())), (:Product, :Lactate), true), ((((:CompetitiveInhibitor, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:CompetitiveInhibitor, :Lactate)), :E, ((), ())), (:Product, :Lactate), true), ((((:Product, :NAD),), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), (:Product, :Lactate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:Product, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Product, :Lactate), true), ((((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (((:Product, :Lactate), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:Product, :Lactate), true)), ((((), :E, ((), ())), (((:CompetitiveInhibitor, :Lactate),), :E, ((), ())), (:CompetitiveInhibitor, :Lactate), true), ((((:CompetitiveInhibitor, :NADH),), :E, ((), ())), (((:CompetitiveInhibitor, :Lactate), (:CompetitiveInhibitor, :NADH)), :E, ((), ())), (:CompetitiveInhibitor, :Lactate), true)), ((((), :E, ((), ())), (((:Product, :NAD),), :E, ((), ())), (:Product, :NAD), true), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), (:Product, :NAD), true), ((((:CompetitiveInhibitor, :NADH),), :E, ((), ())), (((:Product, :NAD), (:CompetitiveInhibitor, :NADH)), :E, ((), ())), (:Product, :NAD), true), ((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E, ((), ())), (:Product, :NAD), true)), ((((), :E, ((), ())), (((:Substrate, :NADH),), :E, ((), ())), (:Substrate, :NADH), true), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Substrate, :NADH), true), ((((:CompetitiveInhibitor, :Lactate),), :E, ((), ())), (((:CompetitiveInhibitor, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Substrate, :NADH), true), ((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :NADH), true), ((((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (((:Substrate, :NADH), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:Substrate, :NADH), true)), ((((), :E, ((), ())), (((:CompetitiveInhibitor, :NADH),), :E, ((), ())), (:CompetitiveInhibitor, :NADH), true), ((((:CompetitiveInhibitor, :Lactate),), :E, ((), ())), (((:CompetitiveInhibitor, :Lactate), (:CompetitiveInhibitor, :NADH)), :E, ((), ())), (:CompetitiveInhibitor, :NADH), true), ((((:Product, :NAD),), :E, ((), ())), (((:Product, :NAD), (:CompetitiveInhibitor, :NADH)), :E, ((), ())), (:CompetitiveInhibitor, :NADH), true), ((((:Substrate, :Pyruvate),), :E, ((), ())), (((:CompetitiveInhibitor, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :NADH), true), ((((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (((:CompetitiveInhibitor, :NADH), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :NADH), true)), ((((), :E, ((), ())), (((:Substrate, :Pyruvate),), :E, ((), ())), (:Substrate, :Pyruvate), false), ((((:Product, :NAD),), :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), false), ((((:Substrate, :NADH),), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), false), ((((:CompetitiveInhibitor, :NADH),), :E, ((), ())), (((:CompetitiveInhibitor, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), false)), ((((), :E, ((), ())), (((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), true), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:Substrate, :NADH), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), true), ((((:CompetitiveInhibitor, :NADH),), :E, ((), ())), (((:CompetitiveInhibitor, :NADH), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), true)), (((((:Product, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:CompetitiveInhibitor, :Lactate)), :E, ((), ())), (:CompetitiveInhibitor, :Lactate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:CompetitiveInhibitor, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:CompetitiveInhibitor, :Lactate), true)), (((((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), nothing, false),)))}, (4, (:EqualAI, :NonequalAI, :EqualAI, :OnlyA, :EqualAI, :EqualAI, :EqualAI, :NonequalAI, :EqualAI)), ()}"
const B2_SIG = "AllostericEnzymeMechanism{EnzymeMechanism{(((((:Product, :Lactate), ((:C, 3), (:H, 6), (:O, 3))), ((:Product, :NAD), ((:C, 21), (:H, 27), (:N, 7), (:O, 14), (:P, 2))), ((:Substrate, :NADH), ((:C, 21), (:H, 29), (:N, 7), (:O, 14), (:P, 2))), ((:Substrate, :Pyruvate), ((:C, 3), (:H, 4), (:O, 3)))), (((:CompetitiveInhibitor, :Lactate), (4,)), ((:CompetitiveInhibitor, :NAD), (4,)), ((:CompetitiveInhibitor, :NADH), (4,)), ((:CompetitiveInhibitor, :Pyruvate), (4,))), (4,)), (((((), :E, ((), ())), (((:Product, :NAD),), :E, ((), ())), (:Product, :NAD), true), ((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E, ((), ())), (:Product, :NAD), true), ((((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (((:Product, :NAD), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:Product, :NAD), true)), ((((), :E, ((), ())), (((:Substrate, :NADH),), :E, ((), ())), (:Substrate, :NADH), true), ((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :NADH), true)), ((((), :E, ((), ())), (((:CompetitiveInhibitor, :NADH),), :E, ((), ())), (:CompetitiveInhibitor, :NADH), true),), ((((), :E, ((), ())), (((:Substrate, :Pyruvate),), :E, ((), ())), (:Substrate, :Pyruvate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), true)), ((((), :E, ((), ())), (((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), false), ((((:Product, :NAD),), :E, ((), ())), (((:Product, :NAD), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), false), ((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Substrate, :Pyruvate), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:CompetitiveInhibitor, :Pyruvate), false)), (((((:Product, :NAD),), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), (:Product, :Lactate), true), ((((:Substrate, :NADH),), :E, ((), ())), (((:Product, :Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Product, :Lactate), true)), (((((:Product, :NAD),), :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), true),), (((((:Substrate, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, ((), ())), nothing, false),), (((((:CompetitiveInhibitor, :Pyruvate),), :E, ((), ())), (((:Substrate, :Pyruvate), (:CompetitiveInhibitor, :Pyruvate)), :E, ((), ())), (:Substrate, :Pyruvate), false),)))}, (4, (:EqualAI, :EqualAI, :EqualAI, :NonequalAI, :OnlyA, :EqualAI, :NonequalAI, :EqualAI, :NonequalAI)), ()}"

println("========== A: _canonical_mechanism is NOT idempotent ==========")
amA = recon(A_SIG); mc = C1(amA)
foundA = false
for g in EnzymeRates.kinetic_groups(amA), idx in 1:length(EnzymeRates.steps(amA)[g])
    length(EnzymeRates.steps(amA)[g]) >= 2 || continue
    ng = EnzymeRates._split_one_step(EnzymeRates.steps(amA), g, idx)
    ns = vcat(EnzymeRates.cat_allo_states(amA), [EnzymeRates.cat_allo_states(amA)[g]])
    raw = EnzymeRates._with_steps_and_cat_states(amA, ng, ns)
    c1 = C1(raw); c2 = C1(c1)
    if c1 != mc && c2 == mc
        global foundA = true
        println("  split(group=$g, member=$idx): C(raw)=$(length(EnzymeRates.steps(c1))) groups != parent's $(length(EnzymeRates.steps(mc))),")
        println("  so `_expand_split_kinetic_group`'s guard `child == C(parent)` KEEPS it — but C(C(raw))=$(length(EnzymeRates.steps(c2))) groups == parent.")
        println("  => a second canonicalization pass would drop this no-op split. FIX 1 = iterate the merge to a fixed point.")
        break
    end
end
foundA || println("  (not reproduced under this build)")

println("\n========== B2: split fuses two Wegscheider-tied Pyruvate groups the parent left separate ==========")
amB = recon(B2_SIG)
emp = EnzymeRates.compile_mechanism(amB); pn = EnzymeRates.fitted_params(emp); mets = EnzymeRates.metabolites(emp); ph = ehash(emp)
println("  parent: np=$(length(pn)) eqhash=$ph  (identifiable rank 7 of 10 — over-parameterized)")
for ch in EnzymeRates._expand_split_kinetic_group(amB)
    emc = EnzymeRates.compile_mechanism(ch)
    EnzymeRates.fitted_params(emc) == pn || continue
    hc = ehash(emc); hc == ph && continue
    Random.seed!(1); m = 0.0
    for _ in 1:3000
        theta = [exp(2.5*randn()) for _ in pn]; c = [exp(3.0*randn()) for _ in mets]
        p = merge(NamedTuple{pn}(Tuple(theta)), (Keq=20000.0, E_total=1.0)); cc = NamedTuple{mets}(Tuple(c))
        vp = EnzymeRates.rate_equation(emp, cc, p); vc = EnzymeRates.rate_equation(emc, cc, p)
        m = max(m, abs(vp-vc)/max(abs(vp), abs(vc), 1e-300))
    end
    if m < 1e-9
        println("  split child: eqhash=$hc (!= parent) yet max rel |v_parent - v_child| = $m over 3000 points")
        println("  => IDENTICAL rate function, different rendered text. The parent keeps two Pyruvate-RE groups")
        println("     (K_A_Pyruvate_E, K_A_Pyruvate_ENAD) SEPARATE even though a Wegscheider box ties them equal;")
        println("     the tie is resolved as two multi-symbol deps (both = the dead-end-SS ratio koff/kon), which the")
        println("     single-symbol rename map misses. The split perturbs the graph so canonicalization DOES fuse them.")
        println("     FIX 2 (dead-end SS binding -> RE) removes the ratio so the tie is single-symbol and the parent fuses too.")
        break
    end
end
println("DONE reproducers")
