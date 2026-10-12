# ABOUTME: Printed rate laws: the full law (parameters, dependent relations, v divided by
# ABOUTME: the normalizer) and the coefficient form, with one legend line per coefficient.

"""`x^e` as text: `x`, `x^2`, `sqrt(x)`, `x^(3/2)`."""
power_text(x, e) = e == 1 ? string(x) : e == 1 // 2 ? "sqrt($x)" :
    isinteger(e) ? "$x^$(Int(e))" : "$x^($(numerator(e))/$(denominator(e)))"

"""Factors as a ratio: `a * b / (c * d)`, `1` when none is above the bar."""
ratio_text(up, dn) = string(isempty(up) ? "1" : join(up, " * "), isempty(dn) ? "" :
    " / " * (length(dn) == 1 ? only(dn) : "(" * join(dn, " * ") * ")"))

"""`p`'s terms as `(c, mo)`, in content order."""
terms_of(p::Poly) = sort!([(c, mo) for (mo, c) in p.terms]; by = last)

"""
    term_text(c, mo, names; rest = nothing, over = false) -> String

The term `c·mo` (`c` positive) as text, times the polynomial `rest` and over `N` when
`over`: `2 * A * k1 * (k2 + k3) / (k4 * N)`; `1` for an empty term.
"""
function term_text(c, mo::Mono, names; rest = nothing, over = false)
    up = [power_text(names[v], e) for (v, e) in mo if e > 0]
    dn = [power_text(names[v], -e) for (v, e) in mo if e < 0]
    numerator(c) == 1 || pushfirst!(up, string(numerator(c)))
    denominator(c) == 1 || pushfirst!(dn, string(denominator(c)))
    rest === nothing || push!(up, "(" * sum_text(rest, names) * ")")
    over && push!(dn, "N")
    ratio_text(up, dn)
end

"""The polynomial `p` as a sum, `a + 2 * b - c`."""
function sum_text(p::Poly, names)
    s = join(string(c < 0 ? " - " : " + ", term_text(abs(c), mo, names))
             for (c, mo) in terms_of(p))
    (startswith(s, " + ") ? "" : "-") * chop(s; head = 3, tail = 0)
end

"""`p` as `c·mu·rest`: its integer content `c`, its monomial content `mu` (each variable's
lowest exponent over the terms) and the rest, `nothing` when that is 1."""
function content(p::Poly{T}) where {T}
    c = foldl(gcd, values(p.terms))
    mu = mono((v => minimum(mo -> mexp(mo, v), keys(p.terms))
               for v in unique(first(t) for mo in keys(p.terms) for t in mo))...)
    rest = Poly{T}(Dict(mmul(mo, mpow(mu, -1)) => div(x, c) for (mo, x) in p.terms))
    c, mu, nterms(rest) == 1 ? nothing : rest
end

"""Print order of concentration monomials: fewest metabolites, then `metabolites` order."""
conc_key(mo, nc) = (sum(last, mo; init = zero(Exp)), [-mexp(mo, i) for i in 1:nc])

"""
    printed(m, L) -> (; rule, N, q, d)

`L` as printed: numerator and denominator divided by the normalizer (`normalizer`, which
names the `rule`), whose integer and monomial contents divide every coefficient while its
other factors, if any, print apart as `N` (their text, else `nothing`). `q` and `d` list
each concentration monomial of `Q` and `D`, in print order, with its coefficient
`(c, mo, rest, over)`: the term `c·mo` times the polynomial `rest` (`nothing` for 1), over
`N` when `over`, which is false when `rest` is `N`'s one factor and cancels.
"""
function printed(m::Mechanism, L::Law)
    rule, lambda = normalizer(m, L)
    D, names = split_conc(L.D, L.vars), L.vars.names
    k, mu, apart = 1 // 1, Mono(), Tuple{Poly, Int}[]
    for (mo, e) in lambda
        c, u, rest = content(D[mo])
        k, mu = k * (c // 1)^e, mmul(mu, mpow(u, e))
        rest === nothing || push!(apart, (rest, e))
    end
    paren(p) = "(" * sum_text(p, names) * ")"
    N = isempty(apart) ? nothing :
        length(apart) == 1 && apart[1][2] == 1 ? sum_text(apart[1][1], names) :
        ratio_text([power_text(paren(p), e) for (p, e) in apart if e > 0],
                   [power_text(paren(p), -e) for (p, e) in apart if e < 0])
    function coef(q)
        c, u, rest = content(q)
        whole = apart == [(rest, 1)]
        (c // k, mmul(u, mpow(mu, -1)), whole ? nothing : rest, !isempty(apart) && !whole)
    end
    groups(p) = sort!([(mo, coef(q)) for (mo, q) in split_conc(p, L.vars)];
                      by = g -> conc_key(g[1], length(metabolites(m))))
    (; rule, N, q = groups(L.Q), d = groups(L.D))
end

"""The printed groups `gs` of `printed` as a sum, each concentration monomial times its
coefficient."""
group_sum(gs, names) = join((term_text(c, mmul(cm, mo), names; rest, over)
                             for (cm, (c, mo, rest, over)) in gs), " + ")

"""The `v` line of a printed law: `ΠS − ΠP/Keq` times the numerator sum `q` of `nq` groups
over the denominator sum `d`."""
v_text(L::Law, q, nq, d) =
    "v = E_total * ($(sum_text(L.H, L.vars.names))) * $(nq == 1 ? q : "($q)") / ($d)"

"""
    rate_equation_string(m, L::Law, sol = solve(m)) -> String

The printed law: the parameters and concentrations it reads, its Wegscheider and Haldane
relations (each dependent constant in the independents), the normalizer's factor `N` printed
apart, if any, and `v = E_total·(ΠS − ΠP/Keq)·Q/D` in the independents with numerator and
denominator divided by the normalizer (`printed`), so every denominator term is
dimensionless.
"""
function rate_equation_string(m::Mechanism, L::Law, sol::Solution = solve(m))
    P, names, cn = printed(m, L), L.vars.names, column_names(m, sol)
    lines = ["(; $(join(parameters(m), ", "))) = params",
             "(; $(join(metabolites(m), ", "))) = concs"]
    for (label, title) in ((:wegscheider, "Wegscheider"), (:haldane, "Haldane"))
        ds = sort!(filter(d -> d.label === label, sol.deps); by = d -> cn[d.col])
        isempty(ds) || push!(lines, "# $title constraints:")
        for d in ds
            mo = mono((L.vars.index[cn[f]] => e for (f, e) in d.terms)...,
                      L.vars.index[:Keq] => d.keq)
            push!(lines, "$(cn[d.col]) = $(term_text(1, mo, names))")
        end
    end
    P.N === nothing || push!(lines, "N = $(P.N)")
    push!(lines, v_text(L, group_sum(P.q, names), length(P.q), group_sum(P.d, names)))
    join(lines, "\n")
end

const SUPERSCRIPTS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")

"""A coefficient's name for the concentration monomial `mo`: `A_B`, `A²_B`, `0` for none."""
conc_label(mo, names) = isempty(mo) ? "0" :
    join((string(names[v], e == 1 ? "" : String([SUPERSCRIPTS[d - '0' + 1]
                                                  for d in string(Int(e))]))
          for (v, e) in mo), "_")

"""
    coefficient_sum(prefix, gs, names) -> (sum, legend)

The printed groups `gs` (`printed`) as a sum of named coefficients, `prefix_A_B * A * B`
(`prefix_0` alone for the constant term, `1` when that is 1), and one legend line per named
coefficient giving its value in the fitted constants.
"""
function coefficient_sum(prefix, gs, names)
    unit(cm, x) = isempty(cm) && x == (1 // 1, Mono(), nothing, false)
    label(cm) = "$(prefix)_$(conc_label(cm, names))"
    join((unit(cm, x) ? "1" : isempty(cm) ? label(cm) :
          "$(label(cm)) * $(term_text(1, cm, names))" for (cm, x) in gs), " + "),
    ["$(label(cm)) = $(term_text(c, mo, names; rest, over))"
     for (cm, (c, mo, rest, over)) in gs if !unit(cm, (c, mo, rest, over))]
end

"""
    coefficient_form(m, L::Law) -> String

The law as `v = E_total·(ΠS − ΠP/Keq)·(n_0 + n_A·A + …) / (1 + d_A·A + d_A_B·A·B + …)`,
each coefficient named by the concentrations it multiplies (`d_A²_B` for `A^2·B`) and
divided by the normalizer (`printed`), then one legend line per coefficient in the fitted
constants, after the normalizer's factor `N` when it prints apart.
"""
function coefficient_form(m::Mechanism, L::Law)
    P, names = printed(m, L), L.vars.names
    (q, lq), (d, ld) = coefficient_sum("n", P.q, names), coefficient_sum("d", P.d, names)
    lines = [v_text(L, q, length(P.q), d),
             "# The coefficients are functions of the fitted constants, not parameters:"]
    P.N === nothing || push!(lines, "N = $(P.N)")
    join([lines; lq; ld], "\n")
end
