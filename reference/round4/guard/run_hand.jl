# ABOUTME: Hand-built guard cases (entered bottomless segment, an LDH-like face trap, RE turnover
# ABOUTME: beside an SS route, no turnover, SS dead end at a face, odd-but-fine mechanisms).
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "advbuild.jl"))
const OUT = joinpath(R, "data")

H(id, subs, prods, steps; kw...) = (id, build(id, subs, prods, steps; kw...))
const HAND = [
    H("bottomless_entered", [:A, :B], [:P, :Q], [
        st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(:A, :B); c = [:B], re = true, g = 2),
        st(E(:A, :B), E(:P, :Q); g = 3),
        st(E(:P), E(:P, :Q); c = [:Q], re = true, g = 4), st(E(:Q), E(:P, :Q); c = [:P], re = true, g = 5),
        st(E(), E(:P); c = [:P], g = 6), st(E(), E(:Q); c = [:Q], g = 7)]),
    H("face_trap_ldh_like", [:A, :B], [:P, :Q], [
        st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(:A, :B); c = [:B], g = 2),
        st(E(:A, :B), E(:P, :Q); g = 3),
        st(E(:P), E(:P, :Q); c = [:Q], re = true, g = 4), st(E(:Q), E(:P, :Q); c = [:P], re = true, g = 5),
        st(E(), E(:P); c = [:P], g = 6), st(E(), E(:Q); c = [:Q], re = true, g = 7)]),
    H("re_turnover_beside_ss_route", [:S], [:P], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); re = true, g = 2),
        st(E(), E(:P); c = [:P], re = true, g = 3),
        st(E(:S), Fc(:P); g = 4), st(Fc(:P), E(); r = [:P], g = 5)]),
    H("no_turnover_path", [:S], [:P], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); g = 2),
        st(E(:P), Fc(); r = [:P], g = 3)]),
    H("ss_dead_end_at_face", [:S], [:P, :Q], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P, :Q); g = 2),
        st(E(:Q), E(:P, :Q); c = [:P], re = true, g = 3), st(E(), E(:Q); c = [:Q], g = 4),
        st(E(:P), E(:P, :Q); c = [:Q], g = 5)]),
    H("abortive_complex", [:A, :B], [:P, :Q], [
        st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(:A, :B); c = [:B], g = 2),
        st(E(:A, :B), E(:P, :Q); g = 3), st(E(:Q), E(:P, :Q); c = [:P], g = 4),
        st(E(), E(:Q); c = [:Q], g = 5), st(E(:Q), E(:A, :Q); c = [:A], re = true, g = 6)]),
    H("product_activated", [:S], [:P], [
        st(E(), E(:P); c = [:P], re = true, g = 1), st(E(:P), E(:P, :S); c = [:S], re = true, g = 2),
        st(E(:P, :S), E(:P, :P); g = 3), st(E(:P), E(:P, :P); c = [:P], re = true, g = 4)]),
    H("ss_flux_free_dead_end", [:S], [:P], [
        st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
        st(E(), E(:I); c = [:I], g = 3), st(E(:S), Fc(:S); g = 4)]),
]

include(joinpath(@__DIR__, "approach_lib.jl"))

function main()
    recs = NamedTuple[]
    for (id, p) in HAND
        mech, err = try
            S.Converted(p, 0).mech, nothing
        catch e
            nothing, errmsg(e)
        end
        push!(recs, (key = "hand/$id", set = "hand", prov = "hand:$id", mech = mech,
                     spec = spec(S.pop_raw(p)), ctor_err = err))
    end
    serialize(joinpath(OUT, "specs_hand.jls"), [(key = r.key, set = r.set, spec = r.spec) for r in recs])
    rng = Xoshiro(3)
    out = Any[]
    for r in recs
        m = r.mech
        if m === nothing
            push!(out, (key = r.key, set = r.set, prov = r.prov, ctor_err = r.ctor_err))
            println(r.key, ": slice ctor error ", r.ctor_err)
            continue
        end
        st_ = slice_struct(m)
        lt = law_test(m)
        num = haskey(lt, :law) ? numeric_checks(m, lt[:law], rng) : nothing
        println("\n", r.key, "\n  struct: ", st_, "\n  law reduced: ", get(lt, :reduced, get(lt, :reduced_err, "-")),
                "\n  law bare: ", get(lt, :bare, get(lt, :bare_err, "-")))
        if num !== nothing
            for f in (:pos, :fwd, :rev)
                n = num[f]
                @printf("  oracle %s: %s law=%.6g o=%s\n", f, n.verdict, n.vlaw, string(round.(n.o; sigdigits = 7)))
            end
            t = lt[:reduced]
            for (f, zs, bad) in ((:fwd, m.products, t.dfwd), (:rev, m.substrates, t.drev))
                bad && println("  approach $f: ", approach(m, lt[:law], zs, rng))
            end
        end
        delete!(lt, :law)
        push!(out, (key = r.key, set = r.set, prov = r.prov, ctor_err = r.ctor_err,
                    allo = S.is_allosteric(m), struct_ = st_, law = lt, num = num))
    end
    serialize(joinpath(OUT, "slice_hand.jls"), out)
end

main()
