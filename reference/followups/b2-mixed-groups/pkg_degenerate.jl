# ABOUTME: Counts init mechanisms with rapid-equilibrium chemistry (fused binding, isomerization,
# ABOUTME: mixed group) and how many of them the beam fits directly (not `_degenerate`).
using EnzymeRates
const ER = EnzymeRates

const bi_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
const bi_bi_pp_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

flat(m) = [s for g in ER.steps(m) for s in g]
re_fused(m) = any(s -> ER.is_equilibrium(s) && ER.is_binding(s) && ER._is_chemistry(s), flat(m))
re_iso(m) = any(s -> ER.is_equilibrium(s) && ER.is_iso(s), flat(m))
re_tc(m) = any(s -> ER.is_equilibrium(s) && !ER.is_iso(s) && !ER.is_binding(s), flat(m))
re_mixed(m) = any(g -> ER.is_equilibrium(first(g)) && any(ER._is_chemistry, g) &&
                       any(!ER._is_chemistry, g), ER.steps(m))
ss_mixed(m) = any(g -> !ER.is_equilibrium(first(g)) && any(ER._is_chemistry, g) &&
                       any(!ER._is_chemistry, g), ER.steps(m))

for (label, rxn) in (("bi_bi", bi_bi_rxn), ("bi_bi_pp", bi_bi_pp_rxn))
    ms = unique(ER.init_mechanisms(rxn))
    println("== $label: $(length(ms)) init mechanisms, $(count(ER._degenerate, ms)) degenerate")
    for (what, f) in (("RE fused binding", re_fused), ("RE isomerization", re_iso),
                      ("RE Theorell-Chance or multi-metabolite step", re_tc),
                      ("RE mixed group", re_mixed), ("SS mixed group", ss_mixed))
        sel = filter(f, ms)
        println("   ", rpad(what, 45), lpad(length(sel), 4), " mechanisms, ",
                lpad(count(!ER._degenerate, sel), 4), " not degenerate (fitted as they are)")
    end
end
