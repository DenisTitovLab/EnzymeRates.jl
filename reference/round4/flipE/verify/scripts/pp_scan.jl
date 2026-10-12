# ABOUTME: Scans every allosteric ping-pong bi-bi mechanism of one expansion for laws that vanish
# ABOUTME: at zero products and for _kcat_forward errors (context for "let the loss decide").
using EnzymeRates, Random, Printf
const ER = EnzymeRates

rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end
t0 = time()
init = ER.init_mechanisms(rxn)
gen1 = ER.expand_mechanisms(init, rxn)
allo = unique(filter(m -> m isa ER.AllostericMechanism, gen1))
println("bi_bi_pp: init $(length(init)), gen1 $(length(gen1)), allosteric $(length(allo)); $(round(time() - t0)) s")
nzero, nkerr, nerr, ex = 0, 0, 0, String[]
tagline(am) = join(String.(ER.cat_allo_states(am)), ",")
for (i, am) in enumerate(allo)
    try
        em = ER.compile_mechanism(am)
        fp = ER.fitted_params(em)
        rng = MersenneTwister(i)
        prm = NamedTuple{(fp..., :Keq, :E_total)}(
            (Tuple(p === :L ? 1.0 : exp(log(0.1) + log(100) * rand(rng)) for p in fp)..., 5.0, 1.0))
        v0 = ER.rate_equation(em, (A = 2.0, B = 3.0, P = 0.0, Q = 0.0), prm)
        vP = ER.rate_equation(em, (A = 2.0, B = 3.0, P = 1.0, Q = 1.0), prm)
        kerr = try ER._kcat_forward(em, prm); false catch; true end
        global nzero += v0 == 0
        global nkerr += kerr
        if (v0 == 0 || kerr) && length(ex) < 3
            push!(ex, "v(products 0) = $v0, v(products 1) = $vP, kcat error = $kerr; tags " *
                      tagline(am) * "\n" * sprint(show, MIME("text/plain"), em))
        end
    catch err
        global nerr += 1
    end
    i % 500 == 0 && (println("$i done, $(round(time() - t0)) s"); flush(stdout))
end
println("allosteric ping-pong laws: $(length(allo)); v = 0 at products 0: $nzero; kcat errors: $nkerr; derivation errors: $nerr")
foreach(e -> println("\n", e), ex)
