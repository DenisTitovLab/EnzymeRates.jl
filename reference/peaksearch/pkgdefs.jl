# ABOUTME: Package (EnzymeRates) versions of every peak-search mechanism and the map from this
# ABOUTME: solver's rate constants to the package's parameter names, shared by the cross-checks.

const PKG = Dict{String, Any}()
PKG["T"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(P, Q) <--> E(P) + Q
        E(Q) <--> E + Q
        E(P) <--> E + P
    end
end
PKG["a"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E + P + Q
    end
end
PKG["b"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) ⇌ E(Q) + P
        E(P, Q) ⇌ E(P) + Q
        E(Q) ⇌ E + Q
        E(P) ⇌ E + P
    end
end
PKG["c"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> Estar(A, B)
        Estar(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(P, Q) <--> E(P) + Q
        E(Q) <--> E + Q
        E(P) <--> E + P
    end
end
PKG["d1"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> Estar(A, B)
        E(A, B) <--> E(P, Q)
        Estar(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(P, Q) <--> E(P) + Q
        E(Q) <--> E + Q
        E(P) <--> E + P
    end
end
PKG["d2"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> Estar(A, B)
        E(A, B) <--> E + P + Q
        Estar(A, B) <--> E + P + Q
    end
end
PKG["e"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(P, Q) <--> E(P) + Q
        E(Q) <--> E + Q
        E(P) <--> E + P
        E(Q) + A <--> E(A, Q)
    end
end
PKG["f"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(P, Q) <--> E(P) + Q
        E(Q) <--> E + Q
        E(P) <--> E + P
        E(Q) + A <--> E(A, Q)
        E(A, Q) <--> E(A) + Q
    end
end
PKG["g"] = () -> @enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)
        E + B <--> E(B)
        E(A) + B <--> E(A, B)
        E(B) + A <--> E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(Q) <--> E + Q
    end
end

"""
The package's parameter values for my rate constants: each package step is matched to my
step with the same forms and ligands, in either direction.
"""
function pkg_values(em, m, kf, kr)
    mech = ER.Mechanism(em)
    sp = ER._step_parameters(mech)
    vals = Dict{Symbol, BigFloat}()
    for (idx, (s, _)) in enumerate(ER._flat_steps(mech))
        key = (ER.name(ER.from_species(s)), ER.name(ER.to_species(s)),
               sort([ER.name(x) for x in ER.consumed(s)]),
               sort([ER.name(x) for x in ER.released(s)]))
        hits = Tuple{BigFloat, BigFloat}[]
        for (i, st) in enumerate(m.steps)
            mine = (m.forms[st.from], m.forms[st.to],
                    sort([LIGANDS[l] for l in st.cons]), sort([LIGANDS[l] for l in st.rel]))
            flipped = (mine[2], mine[1], mine[4], mine[3])
            key == mine && push!(hits, (kf[i], kr[i]))
            key == flipped && push!(hits, (kr[i], kf[i]))
        end
        length(hits) == 1 || error("step $key matched $(length(hits)) of my steps")
        f, r = only(hits)
        if ER.is_equilibrium(s)
            Ka = f / r
            vals[ER.name(sp[idx][1], mech)] = ER.is_binding(s) ? 1 / Ka : Ka
        else
            vals[ER.name(sp[idx][1], mech)] = f
            vals[ER.name(sp[idx][2], mech)] = r
        end
    end
    vals
end

