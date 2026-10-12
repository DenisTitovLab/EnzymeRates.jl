# ABOUTME: Checks the package's existing degeneracy screens (no Vmax, chemistry at equilibrium
# ABOUTME: with both sides) that the probe's guard proposal leaves out: counts and one example law.
using EnzymeRates, Printf
const ER = EnzymeRates

const NoVmax = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)
        E(S) ⇌ E(P)
        E + P ⇌ E(P)
    end
end

function main()
    m = ER.Mechanism(NoVmax)
    println("## a no-Vmax uni-uni (E + S <--> E(S) SS; E(S) ⇌ E(P), E + P ⇌ E(P) RE)")
    println("  _degenerate: ", ER._degenerate(m))
    s = ER.rate_equation_string(NoVmax)
    println("  ", strip(last(split(s, '\n'; keepempty = false))))
    names = ER.fitted_params(NoVmax)
    params = (; (x => 1.0 for x in names)..., Keq = 10.0, E_total = 1.0)
    for S in (1.0, 1e2, 1e4, 1e6)
        @printf("  v(S = %.0e, P = 0) = %.6g\n", S, ER.rate_equation(NoVmax, (S = S, P = 0.0), params))
    end
    for (label, rxn) in (("bi-bi", @enzyme_reaction(begin
                                substrates: A[C], B[N]
                                products: P[C], Q[N]
                            end)),
                         ("bi-bi ping-pong", @enzyme_reaction(begin
                                substrates: A[CX], B[N]
                                products: P[C], Q[NX]
                            end)),
                         ("uni-bi", @enzyme_reaction(begin
                                substrates: S[AB]
                                products: P[A], Q[B]
                            end)))
        ms = ER.init_mechanisms(rxn)
        nd = count(ER._degenerate, ms)
        nv = count(m -> !(ER._has_vmax(ER.steps(m), rxn, ER.Substrate) &&
                          ER._has_vmax(ER.steps(m), rxn, ER.Product)), ms)
        @printf("  init_mechanisms(%s): %d mechanisms, %d degenerate (%d without a Vmax in some direction)\n",
                label, length(ms), nd, nv)
    end
end

main()
