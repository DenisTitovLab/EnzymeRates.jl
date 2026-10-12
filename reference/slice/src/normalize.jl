# ABOUTME: The normalizer of a printed law: its constant term's coefficient, else a virtual
# ABOUTME: constant term (Cleland's, then an affine basis), else one term's coefficient.

"""
    normalizer(m, L) -> (rule, lambda)

The virtual constant term `N = Π c(mono)^e` over `lambda = [mono => e]`, `c(mono)` being the
coefficient of the concentration monomial `mono` in `L.D`, by which printed laws divide
their numerator and denominator. `rule` names the construction, the first that applies:
`:constant`, the constant term's coefficient; `:cleland`, `Π_X c(ΠS/X) / c(ΠS)^(n-1)` when
every factor is present; `:affine`, the coefficients of the first affinely independent
monomials (fewest metabolites, then fewest non-substrates, then substrates first) raised to
the integer exponents, summing to 1, that combine their exponent vectors to zero;
`:coefficient`, the first such monomial's coefficient. Each denominator term divided by `N`
is dimensionless except under `:coefficient`.
"""
function normalizer(m::Mechanism, L::Law)
    D, nc = split_conc(L.D, L.vars), length(metabolites(m))
    haskey(D, Mono()) && return :constant, [Mono() => 1]
    S = conc_mono(L.vars, m.substrates)
    cle = [[mmul(S, mono(L.vars.index[x] => -1)) => 1 for x in m.substrates];
           S => 1 - length(m.substrates)]
    all(p -> haskey(D, p.first), cle) && return :cleland, cle
    ex(mo) = [Int(mexp(mo, i)) for i in 1:nc]
    nonsub(mo) = count(i -> mexp(mo, i) != 0 && !(L.vars.names[i] in m.substrates), 1:nc)
    monos = sort!(collect(keys(D)); by = mo -> (sum(ex(mo)), nonsub(mo), -ex(mo)))
    aff(mo) = Exp[1; ex(mo)]
    basis = Mono[]
    for mo in monos
        length(rref!(hcat(aff.([basis; [mo]])...), length(basis) + 1)) > length(basis) &&
            push!(basis, mo)
    end
    nb = length(basis)
    M = hcat(aff.(basis)..., Exp[1; zeros(Int, nc)])
    rref!(M, nb)
    lam = M[1:nb, end]
    all(iszero, M[(nb + 1):end, end]) && all(isinteger, lam) &&
        return :affine, [b => Int(l) for (b, l) in zip(basis, lam) if !iszero(l)]
    :coefficient, [monos[1] => 1]
end
