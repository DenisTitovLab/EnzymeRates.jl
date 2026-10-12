# ABOUTME: The derivation pipeline: one entry point that takes a mechanism through every
# ABOUTME: stage (checks, solve, law, normalization, printing, key, code generation, kcat).

"""
    derive(m) -> NamedTuple

Every artifact of `m`'s rate law, in pipeline order: the joint solve (`solution`); the law
(`law`: a `Law`, or an `AlloLaw` when `m` is allosteric); the normalizer of each printed
conformation (`normalizers`, `(rule, lambda)` as `normalizer` gives them); the printed law
(`printed`) and its coefficient form (`coefficients`); the dedup key (`key`, a `Canonical`
of the law, before any code); the body that `rate_equation` generates (`code`) and the
compiled singleton it belongs to (`compiled`); and kcat's candidates (`kcat`, as
`kcat_candidates` gives them, or the message of the error saying why there is no kcat).
Errors when `m` fails a boundary check (`check_boundary`) or admission (`admit`).
"""
function derive(m::Mechanism)
    check_boundary(m)
    rejection = admit(m)
    rejection === nothing || error(rejection)
    sol = solve(m)
    L = rate_law(m, sol)
    states = L isa Law ? [L] : filter(!isnothing, [L.active, L.inactive])
    printed, coefficients = rate_equation_string(m, L, sol), coefficient_form(m, L)
    key, code = canonical(L), body(program(m, sol, L))
    kc = try
        kcat_candidates(m, sol)
    catch err
        err isa ErrorException && startswith(err.msg, "kcat: ") || rethrow()
        err.msg
    end
    (; mechanism = m, solution = sol, law = L,
     normalizers = [normalizer(m, s) for s in states], printed, coefficients, key, code,
     compiled = compile(m), kcat = kc)
end
