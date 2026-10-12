# ABOUTME: Evaluates an AlloLaw (MWC formula over its polynomials) at named values in any
# ABOUTME: number type; shared by the allosteric-stage check scripts.

"""`v/E_total` of the AlloLaw `A` at concentrations `conc` and parameters `pars` (by name), in
the number type `T`."""
function allo_value(A, conc, pars, ::Type{T} = BigFloat) where {T}
    x = [A.vars.conc[i] ? T(conc[n]) : T(getfield(pars, n)) for (i, n) in enumerate(A.vars.names)]
    ev(p) = S.evaluate(p, x)
    QA, DA = ev(A.active.Q), ev(A.active.D)
    QI, DI = A.inactive === nothing ? (zero(T), one(T)) : (ev(A.inactive.Q), ev(A.inactive.D))
    fA, fI, n, L = ev(A.fA), ev(A.fI), A.n, x[A.vars.index[:L]]
    WA = prod((ev(s[1])^s[3] for s in A.sites); init = one(T))
    WI = prod((ev(s[2])^s[3] for s in A.sites); init = one(T))
    ev(A.H) * (fI^n * QA * DA^(n - 1) * WA + L * fA^n * QI * DI^(n - 1) * WI) /
        (fI^n * DA^n * WA + L * fA^n * DI^n * WI)
end
