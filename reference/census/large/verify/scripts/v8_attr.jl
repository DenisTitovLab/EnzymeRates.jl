# ABOUTME: The P sub-kinds checked by column deletion at a new point: whether dropping L (L phantom)
# ABOUTME: or a rapid-equilibrium isomerization constant (gauge) keeps the Jacobian's rank.
include(joinpath(@__DIR__, "vlib.jl"))

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)

"""The BigFloat Jacobian `dv/dlog θ` of the slice law of `m` at a random point (rows scaled
to unit max), with the fitted names."""
function law_jac(m, rng; h = big(10.0)^-40)
    L = S.rate_law(m)
    names = collect(S.fitted_params(m))
    keq = lu10(rng, -1, 1)
    θ0 = Dict{Symbol, BigFloat}(n => lu10(rng, -1.5, 1.5) for n in names)
    np = length(names)
    J = zeros(BigFloat, 2np + 6, np)
    for r in axes(J, 1)
        conc = random_conc(m, rng; lo = -1.5, hi = 1.5)
        f(θ) = law_value(L, lawx(L, merge(θ, Dict(:Keq => keq)), conc))
        for (j, n) in enumerate(names)
            θp, θm = copy(θ0), copy(θ0)
            θp[n] *= exp(h); θm[n] *= exp(-h)
            J[r, j] = (f(θp) - f(θm)) / (2h)
        end
        J[r, :] ./= maximum(abs, view(J, r, :))
    end
    J, names
end
rk(J) = count(>(1e-30), pivoted_norms(J))

rng = Xoshiro(20261011 + 8)
tally = Dict{String, Int}()
for name in RX
    data = load_run(name)
    fit = [r for recs in data.rounds for r in recs if r[:fitted] && :P in r[:labels]]
    for kind in (:L_phantom, :re_iso_gauge, :ss_local, :other)
        pool = [r for r in fit if haskey(get(r, :explain, Dict()), kind) && !(:W_admit in r[:labels])]
        isempty(pool) && continue
        for r in pool[sort!(randperm(rng, length(pool))[1:min(5, length(pool))])]
            m = r[:smech]
            J, names = law_jac(m, rng)
            full = rk(J)
            fwd = S.orient(m)
            cn = S.constant_names(m)
            iso = Set(cn[k][1] for (k, e) in enumerate(m.edges) if e.re && e.cons == 0 && e.rel == 0)
            dropL = :L in names ? rk(J[:, findall(!=(:L), names)]) : -1
            isoc = [j for (j, n) in enumerate(names) if n in iso]
            dropiso = [rk(J[:, setdiff(1:length(names), [j])]) for j in isoc]
            ok = kind === :L_phantom ? dropL == full :
                 kind === :re_iso_gauge ? any(==(full), dropiso) : true
            k = string(name, " ", kind, ok ? " confirmed" : " NOT confirmed")
            tally[k] = get(tally, k, 0) + 1
            @printf("   %s %s %s: rank %d of %d; without L %d; without each RE-isomerization constant %s\n",
                    name, r[:id], kind, full, length(names), dropL, dropiso)
        end
    end
end
println("tally: ", sort!(collect(tally)))
