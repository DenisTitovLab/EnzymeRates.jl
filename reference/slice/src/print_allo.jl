# ABOUTME: Printed allosteric (MWC) laws per conformation, in full and in coefficients: each
# ABOUTME: conformation over its own normalizer, with regulatory sites and cross weights.

"""
Conformation `tag`'s law `L` as printed: its normalizer's text `N` (`nothing` for a
monomial) and the groups (`printed`) of `L.Q`, `L.D` and the flip form's weight, each over
the normalizer, so every term of `D` and of the weight is dimensionless. The normalizer's
factor printed apart is the extra last variable of `names`, `N_tag`.
"""
function conformation(m::Mechanism, L::Law, tag)
    names = [L.vars.names; Symbol("N_", tag)]
    over(gs) = [(cm, (c, o ? mmul(mo, mono(length(names) => -1)) : mo, r, false))
                for (cm, (c, mo, r, o)) in gs]
    P = printed(m, L)
    flip = Law(L.vars, L.H, L.trees[1] * L.shift, L.D, L.sigma, L.trees, L.shift)
    (; tag, names, N = P.N, q = over(P.q), d = over(P.d), f = over(printed(m, flip).q))
end

"""
    mwc_lines(m, A, coefficients) -> (definitions, legend, v line)

The MWC law `A` of `m` as printed after its header. Each conformation `s` (`A`, then `I`
unless it is the flip form alone) defines `Q_s` (unless it does not turn over) and `D_s`
over its normalizer `N_s`; the flip weights over the normalizers, `F_A` and `F_I`, only
when they differ and are not 1; `W_s`, the product of its regulatory-site factors. With
`coefficients`, the sums read in coefficients named by conformation and concentrations
(`dA_S`, `nI_0`) that the legend defines in the fitted constants.
"""
function mwc_lines(m::Mechanism, A::AlloLaw, coefficients::Bool)
    cs = Any[conformation(m, A.active, "A")]
    A.inactive === nothing || push!(cs, conformation(m, A.inactive, "I"))
    defs, legend = String[], String[]
    function define!(name, gs, c)
        coefficients || return push!(defs, "$name = $(group_sum(gs, c.names))")
        prefix = string(name[1] == 'Q' ? 'n' : lowercase(name[1]), c.tag)
        text, lines = coefficient_sum(prefix, gs, c.names)
        push!(defs, "$name = $text")
        append!(legend, lines)
    end
    for c in cs
        c.N === nothing || push!(coefficients ? legend : defs, "N_$(c.tag) = $(c.N)")
        isempty(c.q) || define!("Q_$(c.tag)", c.q, c)
        define!("D_$(c.tag)", c.d, c)
    end
    flips = [group_sum(c.f, c.names) for c in cs]
    flips[1] == get(flips, 2, "1") ||
        foreach(((c, t),) -> t == "1" || define!("F_$(c.tag)", c.f, c), zip(cs, flips))
    for (i, tag) in ((1, "A"), (2, "I"))
        ws = [string("(", sum_text(s[i], A.vars.names), ")", s[3] == 1 ? "" : "^$(s[3])")
              for s in A.sites if nterms(s[i]) > 1]
        isempty(ws) || push!(defs, "W_$tag = " * join(ws, " * "))
    end
    has(x) = any(l -> startswith(l, x * " = "), defs)
    pw(x, e) = !has(x) || e == 0 ? nothing : e == 1 ? x : "$x^$e"
    product(fs...) = (f = filter(!isnothing, [fs...]); isempty(f) ? "1" : join(f, " * "))
    branch(s, o, e) = product(pw("F_$o", A.n), pw("Q_$s", e), pw("D_$s", A.n - e),
                              pw("W_$s", 1))
    num = filter(!isnothing, [has("Q_A") ? branch("A", "I", 1) : nothing,
                              has("Q_I") ? "L * " * branch("I", "A", 1) : nothing])
    rest = branch("I", "A", 0)
    den = branch("A", "I", 0) * " + L" * (rest == "1" ? "" : " * " * rest)
    defs, legend, "v = E_total * ($(sum_text(A.H, A.vars.names))) * " *
                  "($(isempty(num) ? "0" : join(num, " + "))) / ($den)"
end

"""
    rate_equation_string(m, A::AlloLaw, sol = solve(m)) -> String

The printed MWC law, valid Julia read in order: the active law's header (parameters,
concentrations, Wegscheider and Haldane relations), each conformation's definitions
(`mwc_lines`), and `v = E_total·(ΠS − ΠP/Keq)·(F_I^n·Q_A·D_A^(n−1)·W_A +
L·F_A^n·Q_I·D_I^(n−1)·W_I) / (F_I^n·D_A^n·W_A + L·F_A^n·D_I^n·W_I)` without factors of 1.
"""
function rate_equation_string(m::Mechanism, A::AlloLaw, sol::Solution = solve(m))
    plain = split(rate_equation_string(m, A.active, sol), "\n")
    stop = findfirst(l -> startswith(l, "N = ") || startswith(l, "v = "), plain)
    defs, _, v = mwc_lines(m, A, false)
    join([plain[1:stop - 1]; defs; v], "\n")
end

"""
    coefficient_form(m, A::AlloLaw) -> String

The MWC law's `v` line, each conformation's sums in coefficients (`D_A = 1 + dA_S * S`),
then the normalizers printed apart and one legend line per coefficient.
"""
function coefficient_form(m::Mechanism, A::AlloLaw)
    defs, legend, v = mwc_lines(m, A, true)
    join([v; defs; "# The coefficients are functions of the fitted constants, not " *
                   "parameters:"; legend], "\n")
end

"""The printed law of `m`, plain or MWC (`rate_law`), and its coefficient form."""
rate_equation_string(m::Mechanism) =
    (sol = solve(m); rate_equation_string(m, rate_law(m, sol), sol))
coefficient_form(m::Mechanism) = coefficient_form(m, rate_law(m))
