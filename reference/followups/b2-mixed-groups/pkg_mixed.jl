# ABOUTME: Package-side checks for mixed kinetic groups: which group/flag combinations the
# ABOUTME: constructor accepts, and why the n = 2 allosteric bi-bi children keep so few of them.
using EnzymeRates
const ER = EnzymeRates

const bi_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
const bi_bi_allo_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    oligomeric_state: 2
end

println("== constructor verdicts")
for (label, body) in (
        ("mixed group, both RE",
         "E + A <--> E(A); (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q)); E(Q) + P <--> E(P, Q); E + Q ⇌ E(Q)"),
        ("mixed group, both SS",
         "E + A ⇌ E(A); (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q)); E(Q) + P ⇌ E(P, Q); E + Q <--> E(Q)"),
        ("mixed group, fused RE + plain SS",
         "E + A <--> E(A); (E(A) + B ⇌ E(P, Q), E(Q) + B <--> E(B, Q)); E(Q) + P <--> E(P, Q); E + Q ⇌ E(Q)"),
        ("Theorell-Chance step grouped with a plain binding",
         "E + A ⇌ E(A); (E(A) + B <--> E(Q) + P, E(Q) + B <--> E(B, Q)); E + Q ⇌ E(Q)"))
    src = "@enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin " *
          replace(body, "; " => "\n") * " end end)"
    try
        m = Core.eval(Main, Meta.parse(src))
        println("  accepted: ", label)
    catch e
        msg = sprint(showerror, e isa LoadError ? e.error : e)
        println("  rejected: ", label, "\n      ", first(split(msg, '\n')))
    end
end

is_plain(s) = !ER._is_chemistry(s)
mixed(g) = any(is_plain, g) && any(ER._is_chemistry, g)
dead_end(sp) = any(x -> x isa ER.Substrate, ER.bound(sp)) && any(x -> x isa ER.Product, ER.bound(sp))

for (label, rxn) in (("bi_bi (multiplicity 1)", bi_bi_rxn), ("bi_bi oligomeric_state 2", bi_bi_allo_rxn))
    ms = unique(ER.init_mechanisms(rxn))
    rows = Dict{Tuple{String, String, Bool}, Vector{Int}}()   # key => [mechanisms, children]
    for m in ms
        mg = [g for g in ER.steps(m) if mixed(g)]
        isempty(mg) && continue
        flag = any(g -> ER.is_equilibrium(first(g)), mg) ? "has an RE mixed group" : "SS mixed only"
        member = any(g -> any(s -> is_plain(s) && !dead_end(ER.to_species(s)), g), mg) ?
                 "a catalytic plain member" : "dead-end plain members only"
        key = (flag, member, ER._hyperbolic_catalysis(m))
        r = get!(rows, key, [0, 0])
        r[1] += 1; r[2] += length(ER._expand_to_allosteric(m, rxn))
    end
    println("== $label: $(length(ms)) init mechanisms; those with a mixed group:")
    for (k, r) in sort!(collect(rows))
        println("   ", rpad(k[1], 22), rpad(k[2], 28), "hyperbolic=", rpad(k[3], 6),
                " mechanisms ", rpad(r[1], 4), " _expand_to_allosteric children ", r[2])
    end
end
