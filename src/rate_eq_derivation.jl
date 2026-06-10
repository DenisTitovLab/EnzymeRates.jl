# ABOUTME: King-Altman/Cha rate equation derivation via @generated functions.
# ABOUTME: Parameters API, kcat computation, and MWC allosteric rate assembly.

# ─── Parameters API ─────────────────────────────────────────

const _AnyMechanism = AbstractEnzymeMechanism

"""
Suffix appended to single-symbol equality lines whose LHS got folded
into the kinetic-group rename map. Both display sites (User defined
kinetic-group merges and absorbed single-symbol Wegscheider ties) emit
this exact string so the rate-equation dedup key in
`identify_rate_equation.jl` can strip these provenance lines.
"""
const ANNOTATION_SUBSTITUTED = "  (substituted into v)"

"""
    parameters(m::EnzymeMechanism, [mode])
    parameters(m::AllostericEnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `Reduced` (default): independent k's + Keq + E_total. The set of
  symbols the user supplies to evaluate the Haldane-reduced rate
  equation. Returned for both `EnzymeMechanism` and
  `AllostericEnzymeMechanism`.
- `Full`: all raw rate-constant symbols + E_total. For
  `EnzymeMechanism` this is "all 2N k's + E_total." For
  `AllostericEnzymeMechanism` it composes the catalytic raw A-state
  symbols + every I-state mirror (catalytic + regulatory + synthesized
  dep) + reg-site A-state K's (skipping `:OnlyI` ligands) + `:L` +
  `:E_total`. The allosteric Full mode enumerates the complete raw-symbol
  set; no `rate_equation` method is defined for
  `(::AllostericEnzymeMechanism, ::FullMode)`, so this mode is for
  symbol enumeration, not runtime evaluation.
"""
function parameters end

parameters(m::_AnyMechanism) = parameters(m, Reduced)

# ── EnzymeMechanism ───────────────────────────────────────────
@generated function parameters(
    ::EnzymeMechanism{Sig}, ::FullMode,
) where {Sig}
    mech = Mechanism(EnzymeMechanism{Sig}())
    params = _enumerate_parameters_full(mech)
    names = Tuple(name(p, mech) for p in params)
    Tuple((names..., :E_total))
end

@generated function parameters(::M, ::ReducedMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

# ── AllostericEnzymeMechanism ────────────────────────────────
@generated function parameters(
    ::M, ::FullMode,
) where {M <: AllostericEnzymeMechanism}
    aem = M()
    am = AllostericMechanism(aem)
    params = _enumerate_parameters_full_allosteric(am)
    names = Symbol[name(p, am) for p in params]

    # Synthesized-dep I-symbols are constraint LHSes the rate-equation
    # body emits when a non-`:NonequalAI`-tagged dep's RHS references a
    # `:NonequalAI` symbol. They have no Parameter-struct rep (the LHS
    # is a derived dep name, not a base k/K/Kreg), so they're produced
    # Symbol-level. Spliced just before the first Kreg / Lallo so they
    # land between catalytic-I mirrors and reg-I mirrors in the output
    # tuple.
    synth_names = _synthesized_dep_i_names(typeof(catalytic_mechanism(aem)),
                                            am)
    if !isempty(synth_names)
        insert_pos = findfirst(p -> p isa Union{Kreg, Lallo}, params)
        idx = insert_pos === nothing ? length(names) + 1 : insert_pos
        splice!(names, idx:(idx - 1), synth_names)
    end

    Tuple((names..., :E_total))
end

@generated function parameters(
    ::M, ::ReducedMode,
) where {M <: AllostericEnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

"""Independent rate constant names for fitting (excludes Keq, E_total)."""
@generated function fitted_params(::M) where {M <: _AnyMechanism}
    _, indep = _dependent_param_exprs(M)
    indep
end

# Convenience methods for the concrete Mechanism / AllostericMechanism
# working representation: lift to the compiled singleton via
# `compile_mechanism`, then dispatch to the singleton method. Callers that
# hold a `Mechanism` (e.g. the `identify_rate_equation` pipeline) use these
# without an explicit lift. `compile_mechanism` does Sig conversion (it
# allocates) — it is NOT the rate_equation hot path. Pass an
# `EnzymeMechanism` / `AllostericEnzymeMechanism` for the 0-alloc/<100ns
# `rate_equation` guarantee.
parameters(m::Union{Mechanism, AllostericMechanism},
           mode::AbstractRateEquationMode) =
    parameters(compile_mechanism(m), mode)
parameters(m::Union{Mechanism, AllostericMechanism}) = parameters(m, Reduced)
fitted_params(m::Union{Mechanism, AllostericMechanism}) =
    fitted_params(compile_mechanism(m))

"""
Build a renaming map for single-symbol Wegscheider RE ties between two
binding K's. Calls `_dependent_param_exprs_kernel` to discover
binding-K-to-binding-K Wegscheider closures of the form `K_a = K_b`
(RHS is a bare Symbol). Both sides must be binding K's (RE step with
metabolite on LHS) — absorbing a binding-K-to-iso-K tie would produce
inconsistent sign-flips when the kernel runs with the full rename,
since the binding-K column is sign-flipped (Kd convention) but the
iso-K column is not.

The rename means the polynomial in `v` uses the representative symbol
directly, so Source-C duplicates (split kinetic groups that
Wegscheider ties back together) collapse at hash time.
"""
function _build_wegscheider_rename_map(M::Type{<:EnzymeMechanism})
    mech = Mechanism(M())
    rename = Dict{Symbol, Symbol}()
    step_params = _step_parameters(mech)
    # binding-K set: value-context rep name of each RE binding step. Walk
    # Mechanism.steps directly — an RE step carrying a bound metabolite is
    # a binding step; step_params is indexed in the same flat order.
    binding_set = Set{Symbol}()
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_set, name(step_params[idx][1], mech))
    end
    # Pass 2: single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _dependent_param_exprs_kernel(M, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

_build_wegscheider_rename_map(m::EnzymeMechanism) = _build_wegscheider_rename_map(typeof(m))

# ─── RE Group Helpers ───────────────────────────────────────

"""
Per-step side breakdown for the rate equation derivation. Returns
`(from_species_sym, to_species_sym, m_lhs_syms, m_rhs_syms)` for a
single `Step`. This is the canonical metabolite-on-which-side
projection: it reads Step fields directly, placing the bound
metabolite on the binding (m_lhs) or release (m_rhs) side from the
canonical metabolite-on-`to_species` placement plus the SS-dissociation
rule. The five-branch logic is load-bearing: SS catalytic-release steps
where the bound metabolite is a Product that doesn't appear in either
bound list put the metabolite on m_rhs, not m_lhs.
"""
function _step_sides(s::Step)
    e_lhs = name(from_species(s))
    e_rhs = name(to_species(s))
    is_iso(s) && return (e_lhs, e_rhs, Symbol[], Symbol[])
    bm = bound_metabolite(s)
    bm_name = name(bm)
    from_bound_names = Set{Symbol}(name(m) for m in bound(from_species(s)))
    to_bound_names = Set{Symbol}(name(m) for m in bound(to_species(s)))
    if bm_name in to_bound_names
        # Canonical binding: bound met on to_species → emit on m_lhs
        return (e_lhs, e_rhs, Symbol[bm_name], Symbol[])
    elseif bm_name in from_bound_names
        # Reverse-canonical (defensive: ctor swap should have prevented this)
        return (e_lhs, e_rhs, Symbol[], Symbol[bm_name])
    elseif !is_equilibrium(s) && bm isa Product &&
           !(isempty(from_bound_names) && isempty(to_bound_names))
        # SS dissociation rule: bound_metabolite is a Product released
        # in an SS catalytic step; not in either bound list (because the
        # Species canonicalization moved it). Emit on m_rhs.
        return (e_lhs, e_rhs, Symbol[], Symbol[bm_name])
    elseif length(from_bound_names) > length(to_bound_names)
        # Bound-list-size fallback: the side with fewer bound metabolites
        # is the release side, so the metabolite goes on m_rhs.
        return (e_lhs, e_rhs, Symbol[], Symbol[bm_name])
    else
        return (e_lhs, e_rhs, Symbol[bm_name], Symbol[])
    end
end

"""
Compute RE-connected groups via union-find over enzyme Species. Walks
`m.steps` directly; nodes are Species values (compared by `==`). Returns
`(enz_species, groups, form_to_group)` with `groups[g] :: Vector{Int}`
indexing into `enz_species`.
"""
function _compute_re_groups(mech::Mechanism)
    enz_species = _enumerate_species(mech)
    N = length(enz_species)
    parent = collect(1:N)
    function find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
        x
    end
    for group in steps(mech)
        for s in group
            is_equilibrium(s) || continue
            i_from = findfirst(==(from_species(s)), enz_species)
            i_to   = findfirst(==(to_species(s)),   enz_species)
            ra, rb = find(i_from), find(i_to)
            ra != rb && (parent[ra] = rb)
        end
    end
    root_to_group = Dict{Int, Int}()
    groups = Vector{Vector{Int}}()
    form_to_group = zeros(Int, N)
    for i in 1:N
        r = find(i)
        g = get!(root_to_group, r) do; push!(groups, Int[]); length(groups) end
        push!(groups[g], i); form_to_group[i] = g
    end
    enz_species, groups, form_to_group
end

"""Distinct enzyme Species in `m.steps`, in step-walk order."""
function _enumerate_species(m::Mechanism)
    seen = Species[]
    for group in steps(m)
        for s in group
            from_species(s) in seen || push!(seen, from_species(s))
            to_species(s)   in seen || push!(seen, to_species(s))
        end
    end
    seen
end

"""
Concentration symbols (substrates ∪ products ∪ regulators) that may appear in
the rate-equation polynomials — the same set `_poly_to_expr` treats as
concentrations and the only symbols `_reduce_conc_lowest_terms` is allowed to
shift. Parameter symbols are never in this set, so the concentration-GCD can
never drop a fitted parameter.
"""
function _concentration_symbols(mech::Mechanism)
    rxn = reaction(mech)
    cs = Set{Symbol}(name(s) for s in substrates(rxn))
    for p in products(rxn);    push!(cs, name(p)); end
    for rm in regulators(rxn); push!(cs, name(regulator(rm))); end
    cs
end

"""
Free enzyme of a rapid-equilibrium segment: the form with the fewest bound
metabolites, tie-broken toward no covalent residual, then a deterministic
name. Referencing each segment's alphas to this form makes the derivation
independent of the step order and yields the readable `1 + [S]/K + …` form.
"""
function _segment_root(group, enz_species)
    argmin(i -> (length(bound(enz_species[i])),
                 has_residual(enz_species[i]) ? 1 : 0,
                 string(name(enz_species[i]))), group)
end

"""
Compute alpha factors (relative concentrations within RE groups) as POLY
values. Iterates `mech.steps` directly. Binding steps' direction comes
from the canonical Step form (metabolite-on-`to_species`); iso steps
are physical-forward (canonicalized in the Mechanism constructor).
`step_to_K[idx]` is the parameter Symbol for the RE step at flat
position `idx` (rep-renamed via the `name(p, m)` chokepoint).
"""
function _compute_alpha(mech::Mechanism, enz_species,
                        enz_name_to_form, groups, step_to_K)
    N = length(enz_species)
    alpha = Vector{POLY}(fill(poly_one(), N))
    flat = _flat_steps(mech)

    for group in groups
        length(group) == 1 && continue
        root = _segment_root(group, enz_species)
        visited = Set{Int}([root])
        queue = [root]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (idx, (s, _)) in enumerate(flat)
                is_equilibrium(s) || continue
                e_l, e_r, m_l, m_r = _step_sides(s)
                i_f = enz_name_to_form[e_l]
                j_f = enz_name_to_form[e_r]
                Ksym = step_to_K[idx]
                Kp = poly_sym(Ksym)
                Kinv = POLY(_mono(Ksym => -1) => 1)
                is_iso = isempty(m_l) && isempty(m_r)
                if i_f == cur && j_f ∉ visited
                    if is_iso
                        alpha[j_f] = poly_mul(alpha[cur], Kp)
                    else
                        f = poly_mul(poly_sym(m_l[1]), Kinv)
                        isempty(m_r) ||
                            (f = poly_mul(f, POLY(_mono(m_r[1] => -1) => 1)))
                        alpha[j_f] = poly_mul(alpha[cur], f)
                    end
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    if is_iso
                        alpha[i_f] = poly_mul(alpha[cur], Kinv)
                    else
                        f = poly_mul(Kp, POLY(_mono(m_l[1] => -1) => 1))
                        alpha[i_f] = poly_mul(alpha[cur], f)
                    end
                    push!(visited, i_f); push!(queue, i_f)
                end
            end
        end
    end
    alpha
end

"""Build rate poly for one SS step direction in Laurent (fractional) form."""
function _ss_contrib(k_poly, mets, i_form, alpha)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(mets)))
    poly_mul(r, alpha[i_form])
end

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───

"""
Build raw numerator and denominator POLYs for the rate equation by
walking the lifted `Mechanism`. Parameter Symbols on the leaves of
`num`/`den` are produced via the `name(p, mech)` chokepoint (which
collapses kinetic-group members to their rep's name). `rename_map` then
applies any single-symbol Wegscheider ties as a post-pass.
"""
function _raw_symbolic_rate_polys(mech::Mechanism, step_params, rename_map,
                                  subs_species, prods_species)
    enz_species, groups, form_to_group = _compute_re_groups(mech)
    enz_name_to_form = Dict{Symbol, Int}(
        name(es) => i for (i, es) in enumerate(enz_species))
    flat = _flat_steps(mech)
    step_to_K = Dict{Int, Symbol}(
        i => name(step_params[i][1], mech)
        for i in eachindex(flat) if is_equilibrium(flat[i][1]))
    alpha = _compute_alpha(mech, enz_species,
                           enz_name_to_form, groups, step_to_K)
    G = length(groups)

    R = [poly_zero() for _ in 1:G, _ in 1:G]
    for (idx, (s, _)) in enumerate(flat)
        is_equilibrium(s) && continue
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(s)
        i_form = enz_name_to_form[e_lhs]
        j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        kf_poly = poly_sym(name(step_params[idx][1], mech))
        kr_poly = poly_sym(name(step_params[idx][2], mech))
        R[g1, g2] = poly_add(R[g1, g2],
            _ss_contrib(kf_poly, m_lhs, i_form, alpha))
        R[g2, g1] = poly_add(R[g2, g1],
            _ss_contrib(kr_poly, m_rhs, j_form, alpha))
    end

    L = [i == j ? poly_zero() : poly_neg(R[i,j])
         for i in 1:G, j in 1:G]
    for i in 1:G
        L[i, i] = reduce(poly_add, R[i, j] for j in 1:G if j != i; init=poly_zero())
    end
    D = [begin
        idx = [r for r in 1:G if r != root]
        isempty(idx) ? poly_one() : sym_det(L[idx, idx], G - 1)
    end for root in 1:G]

    den = poly_zero()
    for g in 1:G
        sigma = reduce(poly_add, (alpha[i] for i in groups[g]); init=poly_zero())
        csigma = _rename_symbols(sigma, rename_map)
        den = poly_add(den, poly_mul(csigma, D[g]))
    end

    num = _compute_numerator(
        mech, enz_name_to_form, step_params,
        alpha, form_to_group,
        D, subs_species, prods_species)

    num = _rename_symbols(num, rename_map)
    den = _rename_symbols(den, rename_map)
    conc_set = _concentration_symbols(mech)
    num, den = _reduce_conc_lowest_terms(num, den, conc_set)
    num, den
end

function _raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})
    mech = Mechanism(M())
    step_params = _step_parameters(mech)
    rename_map = _build_wegscheider_rename_map(M)
    # substrates(::EnzymeReaction) returns Vector{Substrate} (concrete metabolite
    # structs); _compute_numerator compares them as Symbols (`name(...) == S`,
    # `S in subs_in(...)`). Wrap explicitly.
    subs_syms = Symbol[name(s) for s in substrates(mech.reaction)]
    prods_syms = Symbol[name(p) for p in products(mech.reaction)]
    _raw_symbolic_rate_polys(mech, step_params, rename_map,
                              subs_syms, prods_syms)
end

"""
Numerator = net flux across one complete steady-state reaction-cut. Each
per-turnover-conserved "event" is a candidate cut whose SS-step fluxes sum to v:
metabolite cuts (bind a substrate / release a product / iso-convert a substrate /
iso-produce a product) and central-species cuts (produce / consume an iso-step
endpoint form). Dead-end binding/release steps touching a substrate-product mixed
complex are excluded (a chemistry step producing such a complex — a ping-pong
covalent intermediate — is kept). A candidate is usable iff all its steps are SS;
NUM = oriented-flux sum over the chosen cut (prefer a metabolite cut — always
complete — over a central-species cut, then fewest steps, then a chemistry cut,
then sorted indices). No usable candidate ⇒ a complete all-RE catalytic cycle ⇒
no finite rate ⇒ raise. A central-species cut is trusted only when its form is
unique up to bound regulators; a regulator-variant sibling is a parallel route the
single-form cut would undercount, so raise rather than return a wrong rate.
"""
function _compute_numerator(
    mech::Mechanism, enz_name_to_form, step_params,
    alpha, form_to_group, D, subs_species, prods_species,
)
    is_mixed_substrate_product_complex(f) =
        any(m -> m isa Substrate, bound(f)) && any(m -> m isa Product, bound(f))
    subs_in(f)  = Set(name(m) for m in bound(f) if m isa Substrate)
    prods_in(f) = Set(name(m) for m in bound(f) if m isa Product)

    # Forward-oriented reaction steps (skip inhibitor/regulator binding, pure
    # conformational isos, and product-rebinding dead-ends). For each: type ∈
    # {:bind,:chem,:release}; (ff,ft) = forward (from,to). A product release stored
    # as product-binding (`E+P→EP`, product on m_lhs) is the reverse reaction, so
    # its forward direction swaps the endpoints; an SS-dissociation release
    # (`EA→E+P`, product on m_rhs) is already forward.
    rsteps = NamedTuple[]
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        bm = bound_metabolite(s)
        local ff, ft, typ
        if bm === nothing                                   # iso step
            f, t = from_species(s), to_species(s)
            (Set(name(m) for m in bound(f)) == Set(name(m) for m in bound(t)) &&
             residual(f) == residual(t)) && continue        # pure conformational iso
            ff, ft, typ = f, t, :chem
        elseif bm isa Substrate
            ff, ft, typ = from_species(s), to_species(s), :bind
        elseif bm isa Product
            _, _, m_lhs, _ = _step_sides(s)
            rev = !isempty(m_lhs)                           # product-binding storage
            ff = rev ? to_species(s) : from_species(s)
            ft = rev ? from_species(s) : to_species(s)
            typ = :release
        else
            continue                                        # inhibitor / regulator
        end
        typ !== :chem &&
            (is_mixed_substrate_product_complex(ff) ||
             is_mixed_substrate_product_complex(ft)) && continue
        push!(rsteps, (idx = idx, ff = ff, ft = ft, typ = typ, s = s))
    end

    # Candidate cuts: (steps into rsteps, central form or `nothing` for metabolite).
    cands = Tuple{Vector{Int}, Union{Species, Nothing}}[]
    add_cand!(central, pred) = (g = [k for k in eachindex(rsteps) if pred(rsteps[k])];
                                isempty(g) || push!(cands, (g, central)))
    for S in subs_species
        add_cand!(nothing, r -> r.typ === :bind && name(bound_metabolite(r.s)) == S)
        add_cand!(nothing, r -> r.typ === :chem && S in subs_in(r.ff) && !(S in subs_in(r.ft)))
    end
    for P in prods_species
        add_cand!(nothing, r -> r.typ === :release && name(bound_metabolite(r.s)) == P)
        add_cand!(nothing, r -> r.typ === :chem && P in prods_in(r.ft) && !(P in prods_in(r.ff)))
    end
    central_forms = Set{Species}()      # iso-step endpoints
    for r in rsteps
        r.typ === :chem && (push!(central_forms, r.ff); push!(central_forms, r.ft))
    end
    for X in sort(collect(central_forms); by = x -> string(name(x)))  # sorted: precompile-stable
        add_cand!(X, r -> r.ft == X)    # produce X
        add_cand!(X, r -> r.ff == X)    # consume X
    end

    usable = [c for c in cands if all(!is_equilibrium(rsteps[k].s) for k in c[1])]
    isempty(usable) && error(
        "rate_equation: no rapid-equilibrium-consistent reaction cut — a complete " *
        "all-RE catalytic cycle exists, so the mechanism has no finite rate.")

    # Prefer a metabolite cut (central === nothing) over a central-species cut, then
    # fewest steps, then a chemistry cut, then sorted indices (precompile-stable).
    has_chem(c) = any(rsteps[k].typ === :chem for k in c[1])
    steps, central = usable[argmin(
        i -> (usable[i][2] === nothing ? 0 : 1, length(usable[i][1]),
              has_chem(usable[i]) ? 0 : 1, sort([rsteps[k].idx for k in usable[i][1]])),
        eachindex(usable))]

    # A central-species cut is complete only when its form is unique up to bound
    # regulators; a regulator-variant sibling is a parallel route this cut would
    # undercount, so raise rather than return a wrong rate.
    if central !== nothing
        xs = sort!([name(m) for m in bound(central) if m isa Substrate])
        xp = sort!([name(m) for m in bound(central) if m isa Product])
        for Y in _enumerate_species(mech)
            Y == central && continue
            sort!([name(m) for m in bound(Y) if m isa Substrate]) == xs &&
                sort!([name(m) for m in bound(Y) if m isa Product]) == xp &&
                residual(Y) == residual(central) && error(
                "rate_equation: ambiguous central-complex cut on $(name(central)) — " *
                "the regulator-variant form $(name(Y)) is a parallel route this cut " *
                "would undercount; cannot derive a unique rate.")
        end
    end

    num = poly_zero()
    for k in steps
        r = rsteps[k]; idx = r.idx
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(r.s)
        i_form = enz_name_to_form[e_lhs]; j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        rf = _ss_contrib(poly_sym(name(step_params[idx][1], mech)), m_lhs, i_form, alpha)
        rr = _ss_contrib(poly_sym(name(step_params[idx][2], mech)), m_rhs, j_form, alpha)
        canon = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        num = poly_add(num, (r.typ === :release && !isempty(m_lhs)) ? poly_neg(canon) : canon)
    end
    num
end

# ─── Expr generation from POLY ──────────────────────────────

"""
Compute the raw rate expression (bare symbols) and sorted parameter/concentration symbols.
Returns `(expr, all_params, sorted_concs)`.
"""
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    num_expr = _poly_to_expr(num, param_syms, conc_syms)
    den_expr = _poly_to_expr(den, param_syms, conc_syms)
    expr = :(E_total * ($num_expr) / ($den_expr))
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end

# ─── Mode-dispatched rate_equation ────────────────────────────

"""
    rate_equation(m::EnzymeMechanism, concs, params, [mode])

Compute the QSSA steady-state rate. The body is generated at compile time
as a single arithmetic expression with no allocations, loops, or matrix ops.
"""
function rate_equation end

rate_equation(m::_AnyMechanism, concs, params) = rate_equation(m, concs, params, Reduced)

# Concrete-mechanism convenience: lift to the singleton (see the note on
# the parameters/fitted_params methods above — allocates, not the hot path).
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params,
              mode::AbstractRateEquationMode) =
    rate_equation(compile_mechanism(m), concs, params, mode)
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params) =
    rate_equation(m, concs, params, Reduced)

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple, ::FullMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, FullMode)
end

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple,
    ::ReducedMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, ReducedMode)
end

# ─── String Representation ────────────────────────────────────

"""
    rate_equation_string(m, [mode]) -> String

Return the symbolic rate equation for mechanism `m` as a multi-line
`String` (it returns the text — it does not print). `mode` is `Reduced`
(default) or `Full`; pass a concrete `Mechanism` / `AllostericMechanism`
or its compiled [`EnzymeMechanism`](@ref) singleton.

The string is a runnable transcript of how [`rate_equation`](@ref)
evaluates: a `(; …) = params` destructure line, a `(; …) = concs`
destructure line, then the `v = E_total * (num) / (den)` line. In
`Reduced` mode, dependent rate constants are listed first under
`# Wegscheider constraints:` and `# Haldane constraints:` headers — the
thermodynamic identities that eliminate parameters — and only the
independent set appears in the `params` destructure. In `Full` mode every
rate constant is independent, so there is no constraint section. `Full`
mode is defined for `EnzymeMechanism` only; an `AllostericEnzymeMechanism`
supports `Reduced` mode only.

Use `print` on the result to see the multi-line layout without escaped
newlines.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> print(rate_equation_string(m))
(; K_P_E, K_S_E, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E) / (1 + P / K_P_E + S / K_S_E)
```
"""
function rate_equation_string end

rate_equation_string(m::_AnyMechanism) = rate_equation_string(m, Reduced)

# Concrete-mechanism convenience: lift to the singleton.
rate_equation_string(m::Union{Mechanism, AllostericMechanism},
                     mode::AbstractRateEquationMode) =
    rate_equation_string(compile_mechanism(m), mode)
rate_equation_string(m::Union{Mechanism, AllostericMechanism}) =
    rate_equation_string(m, Reduced)

"""Build the `v = E_total * (num) / (den)` line from the raw symbolic rate polys."""
function _rate_v_line(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(m))
    cs = Set{Symbol}(metabolites(m))
    "v = E_total * ($(_expr_to_string(_poly_to_expr(num, ps, cs)))) / " *
        "($(_expr_to_string(_poly_to_expr(den, ps, cs))))"
end

function rate_equation_string(::M, ::FullMode) where {M<:EnzymeMechanism}
    mech = Mechanism(M())
    param_names = Symbol[name(p, mech) for p in _enumerate_parameters_full(mech)]
    lines = ["(; $(join((param_names..., :E_total), ", "))) = params",
             "(; $(join(metabolites(M()), ", "))) = concs"]
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end

function rate_equation_string(::M, ::ReducedMode) where {M<:EnzymeMechanism}
    m = M()
    _, indep = _dependent_param_exprs(M)

    # Wegscheider/Haldane: single-symbol entries get the substituted-into-v
    # annotation; multi-symbol RHSes get runtime assignment in
    # `_build_rate_body` (no annotation).
    dep_raw, _ = _dependent_param_exprs_kernel(M, Dict{Symbol, Symbol}())
    keq_set = Set([:Keq])
    weg_lines, hal_lines = String[], String[]
    for (sym, expr) in sort(collect(dep_raw); by=p -> string(p[1]))
        is_haldane = _expr_references_any(expr, keq_set)
        suffix = expr isa Symbol ? ANNOTATION_SUBSTITUTED : ""
        line = "$sym = $(string(expr))$suffix"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    lines = ["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
             "(; $(join(metabolites(m), ", "))) = concs"]
    isempty(weg_lines)  ||
        (push!(lines, "# Wegscheider constraints:");  append!(lines, weg_lines))
    isempty(hal_lines)  ||
        (push!(lines, "# Haldane constraints:");      append!(lines, hal_lines))
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end

# ─── kcat Computation Helpers ──────────────────────────────────

"""
Set of Symbol names for SS rate-constant parameters (Kon, Koff, Kfor,
Krev) of `em`. For `AllostericEnzymeMechanism`, also includes the
`_T`-suffixed names of every SS rate constant that lives in the
inactive state polynomial. Routes Symbol production through the
`name(p, m)` chokepoint via Parameter-subtype dispatch. Used by
`rescale_parameter_values` to scale only SS k's without touching RE
Kd's, Keq, E_total, L, or regulatory K's.
"""
function _ss_rate_constant_names(em::EnzymeMechanism)
    mech = Mechanism(em)
    Set{Symbol}(name(p, mech) for p in _enumerate_parameters_full(mech)
                if p isa Union{Kon, Koff, Kfor, Krev})
end

function _ss_rate_constant_names(em::AllostericEnzymeMechanism)
    am = AllostericMechanism(em)
    a_names = Set{Symbol}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        rep = _group_rep(group, fes)
        is_equilibrium(rep) && continue
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        for p in _emit_cat_params_for_rep(rep, st)
            push!(a_names, name(p, am))
        end
    end
    i_names = Set{Symbol}(name(p, am) for p in _all_i_state_parameters(am)
                          if p isa Union{Kon, Koff, Kfor, Krev})
    union(a_names, i_names)
end

"""Group `num` and `den` POLYs by metabolite monomial pattern, using
`k_param_names` to classify each symbol. Returns `(num_groups,
den_groups)` where each value is a POLY of k-monomials sharing the
same met-monomial. Reverse (negative-coefficient) terms are dropped
from the numerator. Used by `_kcat_forward` to compare saturating
metabolite patterns across A/I states."""
function _kcat_groups_from_polys(num::POLY, den::POLY,
                                  k_param_names::Set{Symbol})
    function split_mono(mono::MONO)
        k_mono = MONO()
        met_mono = MONO()
        for (s, e) in mono
            if s in k_param_names || s == :Keq
                push!(k_mono, s => e)
            elseif s != :E_total
                push!(met_mono, s => e)
            end
        end
        sort!(k_mono; by=first), sort!(met_mono; by=first)
    end

    num_groups = Dict{MONO, POLY}()
    for (mono, coeff) in num
        coeff > 0 || continue
        k_part, met_part = split_mono(mono)
        p = get!(num_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    den_groups = Dict{MONO, POLY}()
    for (mono, coeff) in den
        k_part, met_part = split_mono(mono)
        p = get!(den_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    num_groups, den_groups
end

"""
    _kcat_forward(m::EnzymeMechanism, params) → Float64

Compute kcat (forward) analytically from the polynomial structure.
kcat is the maximum rate at saturating substrates, zero products,
and E_total=1. For mechanisms with multiple catalytic paths
(e.g., non-essential activator), returns the max over all paths.

Groups the rate-equation numerator (forward terms only, positive coefficients)
and denominator by metabolite pattern. For each matching pair, builds k-only
(num_k_expr / den_k_expr) candidates and emits `max(nk/dk, ...)`. Equilibrium
constants cancel between num and den at matching metabolite levels.

Multiple candidates arise for mechanisms with alternative catalytic pathways
(e.g., non-essential activator with/without activator bound).
"""
@generated function _kcat_forward(
    ::M, params::NamedTuple,
) where {M <: EnzymeMechanism}
    num, den = _raw_symbolic_rate_polys(M)
    k_param_names = Set{Symbol}(_raw_param_symbols(M()))
    num_groups, den_groups = _kcat_groups_from_polys(num, den, k_param_names)

    # Build kcat candidates: for each forward numerator metabolite group
    # with a matching denominator group, create (num_k_expr, den_k_expr)
    empty_set = Set{Symbol}()
    components = Tuple{Any, Any}[]
    for (met_key, num_k) in sort!(collect(num_groups); by=first)
        den_k = get(den_groups, met_key, nothing)
        den_k === nothing && continue
        num_expr = _poly_to_expr(num_k, empty_set, empty_set)
        den_expr = _poly_to_expr(den_k, empty_set, empty_set)
        push!(components, (num_expr, den_expr))
    end

    dep_exprs, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    candidates = [:($nk / $dk) for (nk, dk) in components]
    result = length(candidates) == 1 ?
        candidates[1] : Expr(:call, :max, candidates...)
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        assignments...,
        :(return $result))
end

"""
    _kcat_forward(m::AllostericEnzymeMechanism, params) → Float64

Compute kcat (forward) for an MWC allosteric enzyme.

kcat is the maximum rate at saturating substrates, zero products,
E_total=1, over all regulator concentration corners (each regulator
either 0 or saturating).

With regulatory sites: regulators shift A/I balance,
so kcat depends on the regulator corner. We enumerate all 2^n_lig
corners and return the max.
"""
@generated function _kcat_forward(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    params::NamedTuple,
) where {CM,CS,RS}
    M_type = AllostericEnzymeMechanism{CM,CS,RS}
    aem = M_type()
    am  = AllostericMechanism(aem)
    CatN = catalytic_multiplicity(am)

    # Build A-state and I-state polynomials separately so the saturating
    # metabolite pattern can be matched across conformations. I-state zeros
    # `:OnlyA` k symbols at the polynomial level — patterns that only exist
    # via an `:OnlyA` step (e.g. an `:OnlyA` substrate-binding K) drop out
    # of I-state, so their `(num_k_I, den_k_I)` are 0 and the I-state
    # contribution at saturation vanishes.
    num_A_poly, den_A_poly = _raw_symbolic_rate_polys_allosteric(am)
    a_only_syms = _a_only_syms(am)
    rename_I = _a_to_i_rename(am)
    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol pick up an
    # I-state name. Mirrors `_dependent_param_exprs` Pass 2 below.
    dep_A_all, _ = _dependent_param_exprs_allosteric(am)
    _add_case_b_renames!(rename_I, dep_A_all, am)
    num_I_poly = _rename_symbols(
        _zero_symbols_in_poly(num_A_poly, a_only_syms),
        rename_I)
    den_I_poly = _rename_symbols(
        _zero_symbols_in_poly(den_A_poly, a_only_syms),
        rename_I)
    rename_A = _A_rename_parameters(am)
    a_param_names = Set{Symbol}(get(rename_A, s, s) for s in _raw_param_symbols(CM()))
    # Inactive-state param-name set: every active-state name plus its I mirror,
    # sourced via Parameter-struct dispatch through `name(p, am)`. Covers
    # catalytic rate-constant Parameters; synth-dep I-mirrors aren't
    # classifier-relevant since polynomial monomials only contain base k/K
    # symbols.
    i_param_names = union(
        a_param_names,
        Set{Symbol}(name(p, am) for p in _all_i_state_parameters(am)
                    if !(p isa Kreg)))
    num_A_groups, den_A_groups =
        _kcat_groups_from_polys(num_A_poly, den_A_poly, a_param_names)
    num_I_groups, den_I_groups =
        _kcat_groups_from_polys(num_I_poly, den_I_poly, i_param_names)

    # Choose the saturating active-state met pattern (single component for
    # mechanisms exercised here; assert keeps that constraint visible).
    a_keys = sort!([k for k in keys(num_A_groups) if haskey(den_A_groups, k)])
    isempty(a_keys) &&
        error("_kcat_forward: AllostericEnzymeMechanism produced no kcat " *
              "components — saturating-substrate pattern not found in numerator")
    length(a_keys) == 1 ||
        error("_kcat_forward: AllostericEnzymeMechanism with multiple " *
              "saturating-substrate kcat components ($(length(a_keys)) found) " *
              "is unsupported")
    met_key = a_keys[1]
    empty_set = Set{Symbol}()
    num_k_A_expr = _poly_to_expr(num_A_groups[met_key], empty_set, empty_set)
    den_k_A_expr = _poly_to_expr(den_A_groups[met_key], empty_set, empty_set)

    # Inactive state at the same metabolite pattern. Missing → set 0 (the
    # saturating pattern is unreachable in the inactive state, so I contributes
    # neither flux nor enzyme mass at saturation along that pattern).
    num_I_p = get(num_I_groups, met_key, nothing)
    den_I_p = get(den_I_groups, met_key, nothing)
    num_k_I_expr = num_I_p === nothing ? 0 :
        _poly_to_expr(num_I_p, empty_set, empty_set)
    den_k_I_expr = den_I_p === nothing ? 0 :
        _poly_to_expr(den_I_p, empty_set, empty_set)
    i_pattern_dead = den_I_p === nothing

    # Build dependent param assignments for active state
    a_assignments, i_assignments_ =
        _build_dep_assignments(M_type)

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq)

    i_state_dead = _i_state_dead(aem)

    # A_c = num_k_c * den_k_c^(CatN-1), B_c = den_k_c^CatN
    A_A = CatN == 1 ? num_k_A_expr :
        :($(num_k_A_expr) * $(den_k_A_expr)^$(CatN - 1))
    B_A = :($(den_k_A_expr)^$(CatN))

    # Inactive-state catalytic A/B. `A_I` is zero whenever the catalytic cycle
    # is broken; `B_I` additionally drops to 0 if the saturating pattern is
    # absent in the inactive state (its mass is lower-order in subs and vanishes
    # relative to the active state at saturation).
    if i_pattern_dead
        A_I = 0
        B_I = 0
    else
        A_I = i_state_dead ? 0 :
              CatN == 1 ? num_k_I_expr :
                          :($(num_k_I_expr) * $(den_k_I_expr)^$(CatN - 1))
        B_I = :($(den_k_I_expr)^$(CatN))
    end

    # Skip inactive-state dep-assignments when the inactive state is dead (they
    # may reference zero-valued params and cause Inf via 1/K substitutions).
    i_assignments = (i_state_dead || i_pattern_dead) ? Expr[] : i_assignments_

    if isempty(RS)
        # No regulatory sites: single kcat value
        result = :($(CatN) * ($(A_A) + L * $(A_I)) /
                   ($(B_A) + L * $(B_I)))
        return Expr(:block,
            _destructuring_expr(hw_params, :params),
            a_assignments...,
            i_assignments...,
            :(return $result))
    end

    # Enumerate all unique ligands across reg sites
    all_ligs = AllostericRegulator[]
    for site in am.regulatory_sites
        for lig in site.ligands
            lig in all_ligs || push!(all_ligs, lig)
        end
    end
    n_ligs = length(all_ligs)
    lig_idx = Dict(lig => i - 1 for (i, lig) in enumerate(all_ligs))

    # For each corner (2^n_ligs), compute kcat expression
    # At corner: each lig is either 0 or ∞.
    # At ∞: Q_reg_c_i → sum(inv(K_j_c_i) for saturating j)
    # At 0: Q_reg_c_i = 1
    corner_exprs = Any[]
    for mask in 0:(2^n_ligs - 1)
        W_A_factors = Any[]
        W_I_factors = Any[]
        for (site_idx, site) in enumerate(am.regulatory_sites)
            n_reg = site.multiplicity
            sat_terms_A = Any[]
            sat_terms_I = Any[]
            for (lig, tag) in zip(site.ligands, site.allo_states)
                if (mask >> lig_idx[lig]) & 1 == 1
                    if tag !== :OnlyI
                        K_A_sym = name(Kreg(site, lig, :A), am)
                        push!(sat_terms_A, :(inv($K_A_sym)))
                    end
                    if tag !== :OnlyA
                        # `:EqualAI` ligands share the A-state symbol; the
                        # body emits an `:EqualAI` ligand's I-state slot
                        # via the A-state name (no `_T` rename).
                        K_I_state = tag === :EqualAI ? :A : :I
                        K_I_sym = name(Kreg(site, lig, K_I_state), am)
                        push!(sat_terms_I, :(inv($K_I_sym)))
                    end
                end
            end
            if !isempty(sat_terms_A)
                q_A = length(sat_terms_A) == 1 ?
                    sat_terms_A[1] :
                    _nest_binary(:+, sat_terms_A)
                push!(W_A_factors, _power_expr(q_A, n_reg))
            end
            if !isempty(sat_terms_I)
                q_I = length(sat_terms_I) == 1 ?
                    sat_terms_I[1] :
                    _nest_binary(:+, sat_terms_I)
                push!(W_I_factors, _power_expr(q_I, n_reg))
            end
        end
        # Build kcat at this corner. Empty W_A_factors / W_I_factors mean
        # no saturating regulator at that conformation — they default to 1.
        if isempty(W_A_factors) && isempty(W_I_factors)
            # Zero regulators corner
            kcat_expr = :($(CatN) * ($(A_A) + L * $(A_I)) /
                          ($(B_A) + L * $(B_I)))
        else
            W_A = isempty(W_A_factors) ? 1 :
                length(W_A_factors) == 1 ?
                    W_A_factors[1] :
                    _nest_binary(:*, W_A_factors)
            W_I = isempty(W_I_factors) ? 1 :
                length(W_I_factors) == 1 ?
                    W_I_factors[1] :
                    _nest_binary(:*, W_I_factors)
            kcat_expr = :($(CatN) *
                ($(A_A) * $(W_A) + L * $(A_I) * $(W_I)) /
                ($(B_A) * $(W_A) + L * $(B_I) * $(W_I)))
        end
        push!(corner_exprs, kcat_expr)
    end

    result = length(corner_exprs) == 1 ?
        corner_exprs[1] : Expr(:call, :max, corner_exprs...)
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        a_assignments...,
        i_assignments...,
        :(return $result))
end

# ─── Public API: rescale_parameter_values ──────────────────────────

"""
    rescale_parameter_values(m, params::NamedTuple; scale_k_to_kcat=1.0)

Rescale SS rate constants so that `_kcat_forward(m, result) ≈ scale_k_to_kcat`.
Non-SS parameters (K's, Keq, E_total, L, regulatory K's) are unchanged.
"""
function rescale_parameter_values(
    m::_AnyMechanism, params::NamedTuple; scale_k_to_kcat=1.0,
)
    kcat_current = _kcat_forward(m, params)
    scale = scale_k_to_kcat / kcat_current
    ss_names = _ss_rate_constant_names(m)
    NamedTuple{keys(params)}(Tuple(
        k in ss_names ? v * scale : v
        for (k, v) in zip(keys(params), values(params))
    ))
end

# ═══════════════════════════════════════════════════════════════════
# AllostericEnzymeMechanism rate equations (MWC)
#
# MWC rate formula (per conformation c, summed over conformations).
# Let cat_n = catalytic_multiplicity(m):
#   num = cat_n * sum_c( L_c * N_cat_c * Q_cat_c^(cat_n - 1)
#             * prod(Q_reg_i_c^n_reg_i for all regulatory sites i) )
#   den = sum_c( L_c * Q_cat_c^cat_n * prod(Q_reg_i_c^n_reg_i for all regulatory sites i) )
#   v = E_total * num / den
#
# Regulatory sites contribute to BOTH numerator and denominator at their
# multiplicity, regardless of whether n_reg_i matches cat_n.
# ═══════════════════════════════════════════════════════════════════

"""
The I-state catalytic cycle cannot close — and therefore both forward
and reverse net flux vanish — when any `:OnlyA` kinetic group is
present. The Cha polynomial-zeroing approach kills only one half of
the catalytic flux (forward for substrate-OnlyA, reverse for
product-OnlyA), leaving the other half non-zero at chemical
equilibrium. Forcing `N_I = 0` ensures Haldane consistency. Used by
both `rate_equation` (via `_allosteric_num_den_exprs`) and
`_kcat_forward`.
"""
function _i_state_dead(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    any(cat_allo_state(m, g) == :OnlyA for g in kinetic_groups(cm))
end

"""
Set of A-state catalytic parameter Symbol names for an
`AllostericMechanism`. Cached helper for the four call sites in
`_kcat_forward`, `_dependent_param_exprs`, `_build_dep_assignments`,
and `_allosteric_num_den_exprs` so the set is assembled in one place.
"""
_a_only_syms(am::AllostericMechanism) =
    Set{Symbol}(name(p, am) for p in _onlyA_parameters(am))

"""
A → I rename map (Symbol → Symbol) for `:NonequalAI` catalytic-group
parameters. Routes through `name(p, am)`; both keys and values are the
rendered Symbol names of `Kd/Kiso/Kon/Koff/Kfor/Krev` parameters.
"""
_a_to_i_rename(am::AllostericMechanism) =
    Dict{Symbol, Symbol}(
        name(p_A, am) => name(p_I, am)
        for (p_A, p_I) in _I_rename_parameters(am))

"""
Catalytic-cycle `Parameter`s zeroed in the I-state branch (one entry per
`:OnlyA` kinetic group). Per kinetic group, the representative step's
binding/iso × equilibrium/steady-state pair determines the emitted
parameter type(s).
"""
function _onlyA_parameters(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :OnlyA || continue
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :A))
    end
    out
end

"""
A→I rename map for `:NonequalAI` kinetic groups — `Parameter` form. Maps
each A-state catalytic Parameter to its I-state counterpart. Synthesized
dep I-symbols (whose RHS references a renamed `:NonequalAI` symbol) are
emitted Symbol-level by the dep-assignment machinery; they have no
Parameter representation.
"""
function _I_rename_parameters(am::AllostericMechanism)
    rename = Dict{Parameter, Parameter}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :NonequalAI || continue
        rep = _group_rep(group, fes)
        if is_equilibrium(rep)
            if is_binding(rep)
                rename[Kd(rep, :A)] = Kd(rep, :I)
            else
                rename[Kiso(rep, :A)] = Kiso(rep, :I)
            end
        else
            if is_binding(rep)
                rename[Kon(rep, :A)]  = Kon(rep, :I)
                rename[Koff(rep, :A)] = Koff(rep, :I)
            else
                rename[Kfor(rep, :A)] = Kfor(rep, :I)
                rename[Krev(rep, :A)] = Krev(rep, :I)
            end
        end
    end
    rename
end

"""
Rename map converting the catalytic polynomial's `:None`-rendered Symbols to
their `:A`-rendered counterparts for non-`:EqualAI` catalytic groups.

The catalytic polynomial is derived from a non-allosteric `EnzymeMechanism`
whose parameters carry `state=:None`, producing Symbols like `:K_ATP_E`.
Allosteric parameter enumeration uses `state=:A`, producing `:K_A_ATP_E`.
This map bridges the gap so polynomial Symbols align with the allosteric
machinery's expectations.

`:EqualAI` groups are skipped: both `:None` and `:EqualAI` render with an
empty state token, so `name(p_None, am) == name(p_EqualAI, am)` — the rename
is an identity and can be omitted.
"""
function _A_rename_parameters(am::AllostericMechanism)
    rename = Dict{Symbol, Symbol}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :EqualAI && continue
        rep = _group_rep(group, fes)
        if is_equilibrium(rep)
            P = is_binding(rep) ? Kd : Kiso
            rename[name(P(rep, :None), am)] = name(P(rep, :A), am)
        else
            if is_binding(rep)
                rename[name(Kon(rep, :None), am)]  = name(Kon(rep, :A), am)
                rename[name(Koff(rep, :None), am)] = name(Koff(rep, :A), am)
            else
                rename[name(Kfor(rep, :None), am)] = name(Kfor(rep, :A), am)
                rename[name(Krev(rep, :None), am)] = name(Krev(rep, :A), am)
            end
        end
    end
    rename
end

"""
Apply the `:None`→`:A` rename to a raw catalytic polynomial for use in
allosteric context. `:EqualAI` groups need no rename (both states share
an empty state token and thus identical symbol rendering).
"""
function _raw_symbolic_rate_polys_allosteric(am::AllostericMechanism)
    CM = typeof(compile_mechanism(Mechanism(reaction(am), steps(am))))
    num_poly, den_poly = _raw_symbolic_rate_polys(CM)
    rename_A = _A_rename_parameters(am)
    isempty(rename_A) && return num_poly, den_poly
    _rename_symbols(num_poly, rename_A), _rename_symbols(den_poly, rename_A)
end

"""
Return `(dep_exprs, indep_params)` for the catalytic sub-mechanism in
allosteric context, with dep keys and RHS Symbols renamed from `:None`-state
to `:A`-state. `:EqualAI` groups are skipped (identity rename).
"""
function _dependent_param_exprs_allosteric(am::AllostericMechanism)
    CM = typeof(compile_mechanism(Mechanism(reaction(am), steps(am))))
    dep, indep = _dependent_param_exprs(CM)
    rename_A = _A_rename_parameters(am)
    isempty(rename_A) && return dep, indep
    renamed_dep = Dict{Symbol, Union{Symbol, Expr}}()
    for (k, v) in dep
        new_k = get(rename_A, k, k)
        new_v = substitute_params_expr(v, rename_A)
        renamed_dep[new_k] = new_v
    end
    renamed_indep = Symbol[get(rename_A, p, p) for p in indep]
    renamed_dep, renamed_indep
end

"""
All I-state `Parameter`s the rate-equation body emits as constraint LHSes
— `Parameter` form. Catalytic groups: every non-`:OnlyA` group
contributes I-state Parameter(s) for its rep step (`Kd`/`Kiso`/`Kon`+
`Koff`/`Kfor`+`Krev`). Regulator sites: every non-`:OnlyA` ligand
contributes an I-state `Kreg`. Synthesized-dep I-mirrors (deps whose RHS
references a `:NonequalAI` symbol) belong to dep-parameter machinery and
are emitted Symbol-level by the dep-assignment builder.
"""
function _all_i_state_parameters(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :OnlyA && continue
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :I))
    end
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyA && continue
            push!(out, Kreg(site, lig, :I))
        end
    end
    out
end

"""
Enumerate every raw rate-constant `Parameter` for an `AllostericMechanism`
(the complete Full-mode symbol set). Order is:

1. Catalytic A-state Parameter per kinetic group (every group). `:OnlyA`
   and `:NonequalAI` use `state = :A`; `:EqualAI` uses `state = :EqualAI`
   because the symbol is shared with the I-state branch (the chokepoint
   `name(p, m)` renders both to the same `Symbol`).
2. Catalytic I-state mirrors via `_all_i_state_parameters` (skips
   `:OnlyA` groups).
3. Reg-site `Kreg(site, lig, :I)` and `Kreg(site, lig, :A)` per ligand
   (the I-state set skips `:OnlyA` ligands; the A-state set skips
   `:OnlyI` ligands).
4. `Lallo()` for the MWC coupling `L`.

Synthesized-dep I-symbols (derived deps whose RHS references a
`:NonequalAI` symbol) are NOT included — they have no Parameter
representation and are emitted Symbol-level by the caller.

Non-appearing names (e.g., a catalytic I-state mirror that's elided in
a `t_state_dead` mechanism) are harmless to include: unused names are
inert — nothing consumes a name that does not appear in the
rate-equation Exprs. This enumeration intentionally over-emits.
"""
function _enumerate_parameters_full_allosteric(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        rep = _group_rep(group, fes)
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        append!(out, _emit_cat_params_for_rep(rep, st))
    end
    append!(out, _all_i_state_parameters(am))
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyI && continue
            push!(out, Kreg(site, lig, :A))
        end
    end
    push!(out, Lallo())
    out
end

"""
Synthesized I-state dep-symbol names emitted as constraint LHSes by the
rate-equation body whenever a non-`:NonequalAI`-tagged dep symbol's RHS
references a `:NonequalAI` catalytic symbol. These symbols have no
Parameter struct representation (they're derived deps from
Wegscheider/Haldane elimination, not base k/K/Kreg), so they are
produced Symbol-level. Used by `parameters(Full)` to ensure the
complete symbol set covers them.

When the I-state cycle is dead (any `:OnlyA` catalytic group present)
the rate-equation body elides every I-state constraint assignment, so
synthesized dep I-names never appear in the body and need not be
emitted here.
"""
function _synthesized_dep_i_names(::Type{CM}, am::AllostericMechanism,
                                  ) where {CM <: EnzymeMechanism}
    any(cat_allo_state(am, g) === :OnlyA for g in 1:length(steps(am))) &&
        return Symbol[]
    dep_A_all, _ = _dependent_param_exprs_allosteric(am)
    nonequalai_set = Set{Symbol}(
        name(p, am) for (p, _) in _I_rename_parameters(am))
    isempty(nonequalai_set) && return Symbol[]
    names = Symbol[]
    for (k, v) in sort(collect(dep_A_all); by=first)
        k in nonequalai_set && continue
        _expr_references_any(v, nonequalai_set) || continue
        push!(names, _dep_inactive_name(am, k))
    end
    names
end

# ─── Dependent parameter expressions ─────────────────────────────

# Distinct inactive-state name for a *dependent* parameter being promoted to
# per-state (Case B: an `:EqualAI` dep whose Haldane/Wegscheider RHS references
# a `:NonequalAI` symbol). For a `:NonequalAI`/`:A` dep, `_flip_to_inactive`
# already yields a distinct `:I` name; for an `:EqualAI` dep it is a no-op, so
# fall back to the forced `:I` variant to avoid a self-map.
function _dep_inactive_name(am, k::Symbol)
    p = _param_for_symbol(am, k)
    nm = name(_flip_to_inactive(p), am)
    nm == k ? name(_force_inactive(p), am) : nm
end

# Pass 2 of the I-rename construction, shared across the four synth-dep sites so
# the synthesized inactive names are identical everywhere. Returns the Pass-1
# key set (callers that need to distinguish synthesized entries use it).
function _add_case_b_renames!(rename_I::Dict{Symbol, Symbol}, deps, am)
    renamed_set = Set{Symbol}(keys(rename_I))
    for (k, v) in deps
        haskey(rename_I, k) && continue
        _expr_references_any(v, renamed_set) || continue
        rename_I[k] = _dep_inactive_name(am, k)
    end
    return renamed_set
end

"""
    _dependent_param_exprs(M::Type{<:AllostericEnzymeMechanism})

Return `(dep_exprs, indep_params)` for an AllostericEnzymeMechanism.
A-state expressions come from `_dependent_param_exprs(CM)`; I-state entries
copy A with `_I` rename, with per-group filtering (`:OnlyA` zeroed in I,
`:EqualAI` shares A symbol).
Adds reg site params and L to indep.

Parameter Symbol production routes through `name(p, am)` (chokepoint),
walking `Kd/Kiso/Kon/Koff/Kfor/Krev`/`Kreg` structs via the
`_onlyA_parameters` / `_I_rename_parameters` helpers. RHS expression
substitution stays Symbol-keyed because dep RHSes are `Expr` trees the
substitution machinery walks at Symbol level.
"""
function _dependent_param_exprs(
    ::Type{AllostericEnzymeMechanism{CM,CS,RS}},
) where {CM,CS,RS}
    aem = AllostericEnzymeMechanism{CM,CS,RS}()
    am  = AllostericMechanism(aem)
    dep_A_all, indep_A_all = _dependent_param_exprs_allosteric(am)

    a_only_syms = _a_only_syms(am)

    dep_A = Dict{Symbol, Union{Symbol, Expr}}(dep_A_all)
    indep_A = collect(indep_A_all)

    # Pass 1: base I-rename for `:NonequalAI` catalytic-group Parameters.
    rename_I = _a_to_i_rename(am)
    # Pass 2: dep RHSes whose expression references a `:NonequalAI`
    # symbol need their own I-state name. Mirrors the second pass in
    # `_build_dep_assignments`: after Gaussian elimination, dep RHSes
    # reference only independent params, so a single non-iterating pass
    # suffices. Synthesized dep I-names are produced through the chokepoint
    # via _flip_to_inactive + _param_for_symbol so no string surgery is needed.
    _add_case_b_renames!(rename_I, dep_A_all, am)

    t_state_dead_flag = _i_state_dead(aem)

    dep_I = Dict{Symbol, Union{Symbol, Expr}}()
    indep_I_list = Symbol[]

    # Generate I-state dep entries for every A-state dep that has an
    # I-state version per `rename_I`. Covers both Case A (dep symbol's
    # catalytic group is `:NonequalAI` — the symbol itself is in
    # rename_I) and Case B (dep symbol is `:EqualAI`-tagged but its RHS
    # references a `:NonequalAI` symbol — Pass 2 above added the
    # synthesized mapping).
    for (k, v) in dep_A_all
        _expr_references_any(v, a_only_syms) && continue
        i_k = get(rename_I, k, nothing)
        i_k === nothing && continue
        dep_I[i_k] = substitute_params_expr(v, rename_I)
    end

    # When the I-state cycle is dead (any `:OnlyA` group), skip
    # generating `:EqualAI` mirror entries (K1_T = K1, k5f_T = k5f,
    # etc.). They're already elided from the rate equation body in
    # `_build_allosteric_rate_body`, so producing them here only
    # inflates `length(dep_exprs)`.
    #
    # Additionally: even when `t_state_dead_flag` is false, a
    # `:NonequalAI` symbol `p` can be "phantom" — declared as
    # `:NonequalAI` but whose underlying `p` only appears in A-state
    # monomials that ALSO contain a `:OnlyA` symbol. Those monomials
    # get zeroed in the I-state polynomial (via
    # `_zero_symbols_in_poly`), so `p` never survives I-state masking
    # and its `_T`-suffixed name doesn't appear anywhere in the rate
    # equation body. Adding it to indep would expose a fittable
    # parameter the optimizer searches over but that has no effect on
    # the loss — pure dimension bloat. Filter against the set of
    # symbols that actually survive into the I-state polynomial.
    num_A, den_A = _raw_symbolic_rate_polys_allosteric(am)
    i_state_survivors = _i_state_surviving_syms(num_A, den_A, a_only_syms)
    for p in indep_A_all
        p ∈ a_only_syms && continue
        if haskey(rename_I, p)
            p ∈ i_state_survivors || continue
            push!(indep_I_list, rename_I[p])
        elseif !t_state_dead_flag
            p ∈ i_state_survivors || continue
            p_inactive = _flip_to_inactive(_param_for_symbol(am, p))
            i_name = name(p_inactive, am)
            i_name == p && continue  # EqualAI: A and I share the same symbol
            dep_I[i_name] = p
        end
    end

    # Reg-site Parameters via `Kreg` structs + the `name(::Kreg, am)`
    # chokepoint.
    reg_params_a = Symbol[]
    reg_params_i_indep = Symbol[]
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            K_A = name(Kreg(site, lig, :A), am)
            K_I = name(Kreg(site, lig, :I), am)
            tag === :OnlyI || push!(reg_params_a, K_A)
            if tag === :EqualAI
                dep_I[K_I] = K_A
            elseif tag === :NonequalAI || tag === :OnlyI
                push!(reg_params_i_indep, K_I)
            end
        end
    end

    merged_dep = merge(dep_A, dep_I)
    # No final name-sort here (unlike the non-allosteric path): each component
    # is already content-canonical — `indep_A` inherits the sorted catalytic
    # `_dependent_param_exprs`, and the reg params follow the canonical
    # `regulatory_sites` order — so rate-equivalent allosteric mechanisms still
    # render identical strings.
    merged_indep = (indep_A..., indep_I_list...,
                    reg_params_a..., reg_params_i_indep..., :L)
    return merged_dep, merged_indep
end

# `parameters` and `fitted_params` for `AllostericEnzymeMechanism`
# dispatch on explicit per-type methods at the top of this file.

# ─── Rate body building helpers ───────────────────────────────────

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ...
Skips ligands absent from the given conformation (`:OnlyA` in I-state, `:OnlyI`
in A-state). Uses the A-state K symbol when the ligand tag is `:EqualAI`.
Renders K-names via the `name(::Kreg, am)` chokepoint."""
function _reg_site_expr(am::AllostericMechanism, site_idx::Int, inactive::Bool)
    site = regulatory_sites(am)[site_idx]
    terms = Any[1]
    for (lig, tag) in zip(ligands(site), allo_states(site))
        if inactive
            tag === :OnlyA && continue
        else
            tag === :OnlyI && continue
        end
        # `:EqualAI` ligands share the A-state symbol in both conformations;
        # `:NonequalAI` / `:OnlyI` ligands carry a distinct I-state K name.
        state = (inactive && tag in (:NonequalAI, :OnlyI)) ? :I : :A
        K_sym = name(Kreg(site, lig, state), am)
        push!(terms, :($(name(lig)) / $K_sym))
    end
    _nest_binary(:+, terms)
end

"""Raise an expression to an integer power (returns 1 for n=0, expr for n=1)."""
function _power_expr(expr, n::Int)
    n == 0 && return 1
    n == 1 && return expr
    :(($expr)^$n)
end

"""
Build active-state and inactive-state dep-param assignment Exprs.
Returns `(a_assignments::Vector{Expr}, i_assignments::Vector{Expr})`.
Shared by `_build_allosteric_rate_body` and `rate_equation_string`.

Routes all parameter Symbol production through the `name(p, am)`
chokepoint via `_onlyA_parameters`, `_I_rename_parameters`, and
`_all_i_state_parameters`. The set of I-state assignment LHSes that get
emitted for catalytic dep symbols mirrors the catalytic side of
`_all_i_state_parameters` plus synthesized dep I-names (derived deps
referencing a `:NonequalAI` symbol). When the I-state cycle is dead,
synthesized dep I-names are elided here for the same reason the caller
elides `i_assignments` entirely in that case.
"""
function _build_dep_assignments(
    M_type::Type{<:AllostericEnzymeMechanism},
)
    m = M_type()
    am = AllostericMechanism(m)

    dep_A, indep_A = _dependent_param_exprs_allosteric(am)
    sorted_deps = sort(collect(dep_A); by=first)

    a_only_syms = _a_only_syms(am)
    # Pass 1: base I-rename for `:NonequalAI` catalytic-group Parameters.
    rename_I = _a_to_i_rename(am)
    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol need their own
    # I-state name. Mirrors `_dependent_param_exprs` Pass 2; the
    # synthesized entries here are the same ones the rate-equation body
    # consumes via `i_names_set`.
    renamed_set = _add_case_b_renames!(rename_I, dep_A, am)
    I_subs = rename_I

    # Catalytic I-state dep-symbol filter: a catalytic dep gets an I-state
    # assignment if its I-mirror name appears among the catalytic I-state
    # Parameters (every non-`:OnlyA` group when alive; only `:NonequalAI`
    # groups when i_dead — `:EqualAI` groups share the A-state symbol so
    # emit no constraint mirror) OR was synthesized in Pass 2 above
    # (covers `:EqualAI` deps whose RHS references a `:NonequalAI`
    # symbol). When i_dead, synthesized entries are elided since the
    # caller elides `i_assignments` entirely in that case.
    i_dead = _i_state_dead(m)
    i_names_set = Set{Symbol}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        tag = cat_allo_state(am, g)
        tag === :OnlyA && continue
        i_dead && tag !== :NonequalAI && continue
        rep = _group_rep(group, fes)
        if is_equilibrium(rep)
            push!(i_names_set,
                  name(is_binding(rep) ? Kd(rep, :I) : Kiso(rep, :I), am))
        else
            if is_binding(rep)
                push!(i_names_set, name(Kon(rep, :I), am))
                push!(i_names_set, name(Koff(rep, :I), am))
            else
                push!(i_names_set, name(Kfor(rep, :I), am))
                push!(i_names_set, name(Krev(rep, :I), am))
            end
        end
    end
    if !i_dead
        # Synthesized dep I-names: catalytic deps whose RHS references a
        # `:NonequalAI` symbol but whose own catalytic group is NOT
        # `:NonequalAI`. They are entries `rename_I[k]` where `k` was
        # added in the Pass-2 loop above (i.e., not in the original
        # `renamed_set`).
        for (k, v) in rename_I
            k in renamed_set && continue
            push!(i_names_set, v)
        end
    end

    a_assignments = Expr[Expr(:(=), sym, expr_kd) for (sym, expr_kd) in sorted_deps]

    i_assignments = Expr[]

    # `:EqualAI` independent catalytic params: p_I = p_A
    # (must come before dep assignments that reference them)
    for p in indep_A
        p ∈ a_only_syms && continue
        haskey(rename_I, p) && continue  # NonequalAI — handled below
        p_inactive = _flip_to_inactive(_param_for_symbol(am, p))
        i_name = name(p_inactive, am)
        i_name == p && continue  # EqualAI: A and I share the same symbol
        push!(i_assignments, Expr(:(=), i_name, p))
    end

    # `:EqualAI` reg params: K_I_reg = K_A_reg. Routes through `Kreg`
    # chokepoint so the symbol scheme stays in one place.
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :EqualAI || continue
            push!(i_assignments,
                  Expr(:(=), name(Kreg(site, lig, :I), am),
                             name(Kreg(site, lig, :A), am)))
        end
    end

    # Emit an I-state assignment for every dep with an I-state name in
    # `i_names_set`. Covers both `:NonequalAI`-tagged dep symbols (Case A)
    # and synthesized I-names for `:EqualAI`-tagged derived deps whose
    # RHS references a `:NonequalAI` symbol (Case B).
    for (sym, expr_kd) in sorted_deps
        i_sym = get(rename_I, sym, nothing)
        i_sym === nothing && continue
        i_sym in i_names_set || continue
        if _expr_references_any(expr_kd, a_only_syms)
            push!(i_assignments, Expr(:(=), i_sym, 0))
        else
            push!(i_assignments, Expr(:(=), i_sym,
                substitute_params_expr(expr_kd, I_subs)))
        end
    end

    return a_assignments, i_assignments
end

"""
Assemble the MWC numerator and denominator Exprs.
Returns `(full_num, full_den)` where the numerator already includes the
`catalytic_multiplicity` factor.
"""
function _allosteric_num_den_exprs(M_type::Type{<:AllostericEnzymeMechanism})
    m = M_type()
    am = AllostericMechanism(m)
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(am)

    num_A_poly, den_A_poly = _raw_symbolic_rate_polys_allosteric(am)
    rename_A = _A_rename_parameters(am)
    cat_params = Set{Symbol}(get(rename_A, s, s) for s in _raw_param_symbols(CM()))
    cat_mets = Set{Symbol}(metabolites(CM()))

    a_only_syms = _a_only_syms(am)
    # Pass 1: base I-rename for `:NonequalAI` catalytic-group Parameters.
    rename_I = _a_to_i_rename(am)
    # Pass 2: dep RHSes referencing a `:NonequalAI` symbol need their own
    # I-state name so the polynomial rename covers synthesized deps.
    # Mirrors the second pass in `_dependent_param_exprs`.
    dep_A_all, _ = _dependent_param_exprs_allosteric(am)
    _add_case_b_renames!(rename_I, dep_A_all, am)

    N_A = _poly_to_expr(num_A_poly, cat_params, cat_mets)
    Q_A = _poly_to_expr(den_A_poly, cat_params, cat_mets)

    # I-state catalytic Exprs.
    # Zero `:OnlyA` symbols at POLY level, then rename `:NonequalAI`
    # symbols to I-suffixed counterparts. `:EqualAI` symbols pass through
    # unchanged (A-state binding) and resolve through dep-param assignments.
    # When the I-state cycle is broken, force N_I = 0: the Cha framework
    # otherwise produces a non-physical reverse flux from products that
    # have nowhere to go.
    num_i_poly = _rename_symbols(
        _zero_symbols_in_poly(num_A_poly, a_only_syms),
        rename_I)
    den_i_poly = _rename_symbols(
        _zero_symbols_in_poly(den_A_poly, a_only_syms),
        rename_I)
    N_I = _i_state_dead(m) ? 0 :
          _poly_to_expr(num_i_poly, cat_params, cat_mets)
    Q_I = _poly_to_expr(den_i_poly, cat_params, cat_mets)

    reg_Q_A = Any[_reg_site_expr(am, i, false) for i in eachindex(RS)]
    reg_Q_I = Any[_reg_site_expr(am, i, true) for i in eachindex(RS)]

    # Numerator: N × Q_cat^(CatN-1) × all reg-site factors at multiplicity.
    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i], multiplicity(RS[i])))
        end
        _nest_binary(:*, factors)
    end

    # Denominator: Q_cat^CatN × all reg-site factors at multiplicity.
    function make_den_term(Q, reg_Qs)
        factors = Any[_power_expr(Q, CatN)]
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i], multiplicity(RS[i])))
        end
        _nest_binary(:*, factors)
    end

    num_A = make_num_term(N_A, Q_A, reg_Q_A)
    den_A = make_den_term(Q_A, reg_Q_A)
    den_I = make_den_term(Q_I, reg_Q_I)

    if _i_state_dead(m)
        # I-state cycle broken: N_I = 0, so drop the L*num_I term
        # entirely (skip dead numerator branch). Q_I still contributes
        # to denominator as enzyme mass.
        :($(CatN) * $(num_A)), :($(den_A) + L * $(den_I))
    else
        num_I = make_num_term(N_I, Q_I, reg_Q_I)
        :($(CatN) * ($(num_A) + L * $(num_I))), :($(den_A) + L * $(den_I))
    end
end

"""Build the MWC rate equation body as an Expr block."""
function _build_allosteric_rate_body(M_type::Type{<:AllostericEnzymeMechanism})
    full_num, full_den = _allosteric_num_den_exprs(M_type)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # When the I-state cycle is broken, i_assignments (I-state Haldanes
    # and :EqualAI catalytic mirrors K_I = K) become dead code — they're
    # only referenced from the L*num_I branch, which is now elided.
    i_assignments = _i_state_dead(M_type()) ? Expr[] : i_assignments_

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M_type())

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(mets, :concs),
        a_assignments...,
        i_assignments...,
        :(return $rate_expr))
end

# ─── Rate equation dispatch ───────────────────────────────────────

@generated function rate_equation(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {CM,CS,RS}
    _build_allosteric_rate_body(AllostericEnzymeMechanism{CM,CS,RS})
end

# ─── String representation ────────────────────────────────────────

function rate_equation_string(
    ::AllostericEnzymeMechanism{CM,CS,RS}, ::ReducedMode,
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    m = M()
    cm = catalytic_mechanism(m)
    CMT = typeof(cm)
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(m)

    # Active-state Wegscheider/Haldane: single-symbol entries get the
    # substituted-into-v annotation; multi-symbol RHSes get runtime
    # assignment in `_build_rate_body` (no annotation).
    dep_A_raw, _ = _dependent_param_exprs_kernel(CMT, Dict{Symbol, Symbol}())
    keq_set = Set([:Keq])
    weg_lines, hal_lines = String[], String[]
    for (sym, expr) in sort(collect(dep_A_raw); by=p -> string(p[1]))
        is_haldane = _expr_references_any(expr, keq_set)
        suffix = expr isa Symbol ? ANNOTATION_SUBSTITUTED : ""
        line = "$sym = $(string(expr))$suffix"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    # Inactive-state assignments — partitioned by Keq-reference predicate so
    # they fold into the same two sections as active-state lines. Allosteric
    # mechanisms thus use the same three-section structure as
    # non-allosteric.
    _, i_assignments_ = _build_dep_assignments(M)
    i_assignments = _i_state_dead(m) ? Expr[] : i_assignments_
    for a in i_assignments
        sym = a.args[1]
        expr = a.args[2]
        is_haldane = _expr_references_any(expr, keq_set)
        line = "$sym = $(_expr_to_string(expr))"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    # Sort each section lexicographically — load-bearing for eq_hash
    # dedup of allosteric Source-C clusters, since inactive-state lines are
    # appended in iteration order rather than the lexicographic active-state
    # order.
    sort!(weg_lines)
    sort!(hal_lines)

    full_num, full_den = _allosteric_num_den_exprs(M)
    v_line = "v = E_total * ($(_expr_to_string(full_num))) / ($(_expr_to_string(full_den)))"

    lines = ["(; $(join(hw_params, ", "))) = params",
             "(; $(join(mets, ", "))) = concs"]
    isempty(weg_lines)  ||
        (push!(lines, "# Wegscheider constraints:");  append!(lines, weg_lines))
    isempty(hal_lines)  ||
        (push!(lines, "# Haldane constraints:");      append!(lines, hal_lines))
    push!(lines, v_line)
    join(lines, "\n")
end
