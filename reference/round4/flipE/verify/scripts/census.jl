# ABOUTME: Free-form census of enumerated mechanisms (ping-pong, regulator, ter-bi reactions) and
# ABOUTME: a TSV of laws and kcat for a deterministic allosteric sample, to diff two checkouts.
using EnzymeRates, Random, Printf
const ER = EnzymeRates

out = ARGS[1]
rxns = [
    ("bi_bi_pp", @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end),
    ("uni_uni_allo_reg", @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end),
    ("uni_uni_reg_and_inhibitor", @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R2
        competitive_inhibitors: R1
        oligomeric_state: 2
    end),
    ("uni_bi", @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
    end),
    ("ter_bi", @enzyme_reaction begin
        substrates: A[C], B[N], D[X]
        products: P[CN], Q[X]
    end),
]

key(m) = sprint(show, ER.compile_mechanism(m))
pick(v, k) = v[unique(round.(Int, range(1, length(v); length = min(k, length(v)))))]
forms(m) = unique(f for g in ER.steps(m) for s in g for f in (ER.from_species(s), ER.to_species(s)))
isfree(f) = isempty(ER.bound(f)) && isempty(ER.residual(f))

t0 = time()
io = open(out, length(ARGS) > 1 ? "a" : "w")
length(ARGS) > 1 || println(io, join(["rxn", "i", "key_hash", "fitted_params", "eqstring_hash", "v_generic",
                  "v_products0", "v_low", "kcat"], '\t'))
for (rname, rxn) in rxns
    length(ARGS) > 1 && !(rname in ARGS[2:end]) && continue
    init = ER.init_mechanisms(rxn)
    seeds = rname == "ter_bi" ? pick(sort(init; by = key), 30) : init
    gen1 = ER.expand_mechanisms(seeds, rxn)
    allo1 = sort!(unique(filter(m -> m isa ER.AllostericMechanism, gen1)); by = key)
    parents = pick(allo1, 25)
    # Expand each parent alone: expand_mechanisms can throw on a parent whose law is too
    # large to derive; such parents are counted and skipped.
    gen2, nthrow = similar(gen1, 0), 0
    for p in parents
        try
            append!(gen2, ER.expand_mechanisms([p], rxn))
        catch err
            nthrow += 1
            nthrow == 1 && println("  expand threw on a parent: ", first(sprint(showerror, err), 150))
        end
    end
    nthrow > 0 && println("  $rname: expand threw on $nthrow of $(length(parents)) parents")
    seen = Set(allo1)
    allo2 = sort!(unique(filter(m -> m isa ER.AllostericMechanism && !(m in seen), gen2)); by = key)
    everything = unique([init; gen1; gen2])
    nfree = [count(isfree, forms(m)) for m in everything]
    conf = unique(ER.conformation(f) for m in everything for f in forms(m))
    nres = count(m -> any(f -> !isempty(ER.residual(f)), forms(m)), everything)
    @printf("%-26s init %d, gen1 %d (%d allosteric), gen2 %d (%d new allosteric); census over %d distinct: free forms per mechanism %s; conformations %s; with residual forms %d; %.0f s\n",
            rname, length(init), length(gen1), length(allo1), length(gen2), length(allo2),
            length(everything), string(sort(collect(Set(nfree)))), string(conf), nres, time() - t0)
    flush(stdout)
    chosen = [pick(allo1, 40); pick(allo2, 20)]
    for (i, am) in enumerate(chosen)
        em = ER.compile_mechanism(am)
        fp = ER.fitted_params(em)
        rng = MersenneTwister(2000 + i)
        lu(lo, hi) = exp(log(lo) + (log(hi) - log(lo)) * rand(rng))
        vals = Tuple(p === :L ? lu(0.2, 5.0) : lu(0.1, 10.0) for p in fp)
        prm = NamedTuple{(fp..., :Keq, :E_total)}((vals..., lu(0.5, 20.0), 1.0))
        mets = Tuple(ER.metabolites(em))
        prods = Set(ER.name.(ER.products(rxn)))
        conc(f) = NamedTuple{mets}(Tuple(f(x) for x in mets))
        cs = (conc(_ -> lu(0.1, 10.0)), conc(x -> x in prods ? 0.0 : lu(0.1, 10.0)),
              conc(_ -> lu(0.01, 0.1)))
        row = try
            v = [ER.rate_equation(em, c, prm) for c in cs]
            kc = try repr(ER._kcat_forward(em, prm)) catch err "ERR:" * first(sprint(showerror, err), 60) end
            [string(hash(ER.rate_equation_string(em))); repr.(v); kc]
        catch err
            fill("ERR:" * first(sprint(showerror, err), 60), 5)
        end
        println(io, join([rname, string(i), string(hash(key(am))), join(fp, ","), row...], '\t'))
    end
    flush(io)
end
close(io)
println("done in $(round(time() - t0; digits = 1)) s")
