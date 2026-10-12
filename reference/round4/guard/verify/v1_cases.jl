# ABOUTME: The hand-written mechanisms of the guard verification as plain step lists,
# ABOUTME: shared by the oracle scripts (each case lists its faces and reference substrate).
S(f, t, c, r, re) = (f, t, Symbol[c...], Symbol[r...], re)
RE, SS = true, false

const CASES = [
    (name = "A re_turnover_beside_ss_route [S -> P]", subs = [:S], prods = [:P], ref = :S,
     faces = [], steps = [
        S("E", "E(S)", [:S], [], RE), S("E(S)", "E(P)", [], [], RE), S("E", "E(P)", [:P], [], RE),
        S("E(S)", "F(P)", [], [], SS), S("F(P)", "E", [], [:P], SS)]),
    (name = "B mut/expand_bibi/453/flip2 [A + B -> P + Q]", subs = [:A, :B], prods = [:P, :Q],
     ref = :A, faces = [], steps = [
        S("E", "E(A)", [:A], [], RE), S("E", "E(P)", [:P], [], RE), S("E", "E(Q)", [:Q], [], RE),
        S("E(A)", "E(P,Q)", [:B], [], RE), S("E(P)", "E(P,Q)", [:Q], [], RE),
        S("E(Q)", "E(P,Q)", [:P], [], SS)]),
    (name = "C hand LDH-like trap [A + B -> P + Q]", subs = [:A, :B], prods = [:P, :Q], ref = :A,
     faces = [[:P, :Q], [:Q], [:P]], steps = [
        S("E", "E(A)", [:A], [], RE), S("E(A)", "E(A,B)", [:B], [], SS),
        S("E(A,B)", "E(P,Q)", [], [], SS), S("E(P)", "E(P,Q)", [:Q], [], RE),
        S("E(Q)", "E(P,Q)", [:P], [], RE), S("E", "E(P)", [:P], [], SS), S("E", "E(Q)", [:Q], [], RE)]),
    (name = "D enumerator flip child 1 (split E/E(B)) [A + B -> P + Q]", subs = [:A, :B],
     prods = [:P, :Q], ref = :A, faces = [[:P, :Q], [:Q], [:P]], steps = [
        S("E", "E(A)", [:A], [], SS), S("E", "E(B)", [:B], [], SS), S("E", "E(Q)", [:Q], [], RE),
        S("E(A)", "E(P,Q)", [:B], [], SS), S("E(B)", "E(B,Q)", [:Q], [], RE),
        S("E(B)", "E(P,Q)", [:A], [], SS), S("E(Q)", "E(B,Q)", [:B], [], RE),
        S("E(Q)", "E(P,Q)", [:P], [], SS)]),
    (name = "E spec 40 LDH, active conformation as a plain network", subs = [:NADH, :Pyr],
     prods = [:Lac, :NAD], ref = :NADH, faces = [[:Lac, :NAD], [:NAD], [:Lac]], steps = [
        S("E", "E(Lac)", [:Lac], [], SS), S("E", "E(NAD)", [:NAD], [], RE),
        S("E", "E(NADH)", [:NADH], [], SS), S("E", "E(Pyr)", [:Pyr], [], RE),
        S("E(Lac)", "E(Lac,NAD)", [:NAD], [], RE), S("E(Lac)", "E(Lac,NADH)", [:NADH], [], SS),
        S("E(NAD)", "E(Lac,NAD)", [:Lac], [], RE), S("E(NAD)", "E(NAD,Pyr)", [:Pyr], [], RE),
        S("E(NADH)", "E(Lac,NADH)", [:Lac], [], SS), S("E(NADH)", "E(NADH,Pyr)", [:Pyr], [], RE),
        S("E(Pyr)", "E(NAD,Pyr)", [:NAD], [], RE), S("E(Pyr)", "E(NADH,Pyr)", [:NADH], [], SS),
        S("E(Lac,NAD)", "E(NADH,Pyr)", [], [], SS)]),
    (name = "F spec 14 [S1 + S2 -> P] (reverse face)", subs = [:S1, :S2], prods = [:P], ref = :S1,
     faces = [[:S1, :S2]], steps = [
        S("E", "E(P)", [:P], [], RE), S("E", "E(S1)", [:S1], [], SS), S("E", "E(S2)", [:S2], [], RE),
        S("E(P)", "E(S1,S2)", [], [], SS), S("E(S1)", "E(S1,S2)", [:S2], [], RE),
        S("E(S2)", "E(S1,S2)", [:S1], [], RE)]),
    (name = "G hand ss_dead_end_at_face [S -> P + Q]", subs = [:S], prods = [:P, :Q], ref = :S,
     faces = [[:P, :Q]], steps = [
        S("E", "E(S)", [:S], [], RE), S("E(S)", "E(P,Q)", [], [], SS),
        S("E(Q)", "E(P,Q)", [:P], [], RE), S("E", "E(Q)", [:Q], [], SS),
        S("E(P)", "E(P,Q)", [:Q], [], SS)]),
    (name = "H hand no_turnover_path [S -> P]", subs = [:S], prods = [:P], ref = :S,
     faces = [], steps = [
        S("E", "E(S)", [:S], [], RE), S("E(S)", "E(P)", [], [], SS), S("E(P)", "F", [], [:P], SS)]),
    (name = "I enumerator flip child 9 (trap E(Q)) [A + B -> P + Q]", subs = [:A, :B],
     prods = [:P, :Q], ref = :A, faces = [[:P, :Q]], steps = [
        S("E", "E(A)", [:A], [], SS), S("E", "E(P)", [:P], [], RE), S("E", "E(Q)", [:Q], [], SS),
        S("E(A)", "E(A,B)", [:B], [], RE), S("E(A)", "E(A,Q)", [:Q], [], SS),
        S("E(P)", "E(B,P)", [:B], [], RE), S("E(P)", "E(P,Q)", [:Q], [], RE),
        S("E(Q)", "E(A,Q)", [:A], [], SS), S("E(Q)", "E(P,Q)", [:P], [], RE),
        S("E(A,B)", "E(P,Q)", [], [], SS)]),
]
