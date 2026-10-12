# ABOUTME: The skeptic's adversarial mechanisms as specs (40 cases), each with the outcome it
# ABOUTME: targets: a derived law, a rejection by derive, or a constructor error.

E(b...; kw...) = F(:E, b...; kw...)
Fc(b...; kw...) = F(:F, b...; kw...)
site(ls, sts, n) = (ligands = ls, states = sts, multiplicity = n)

const CASES = Any[]
case(name, expect, thunk) = push!(CASES, (name = name, expect = expect, thunk = thunk))

# ── rational exponents ─────────────────────────────────────────────────────────────────
case("half_keq_conf", :derive, () -> build("half_keq_conf", [:S], [:P], [
    st(E(), Fc(); g = 1), st(E(:S), Fc(:S); g = 1),
    st(Fc(), E(:S); c = [:S], g = 2), st(Fc(:S), E(); r = [:P], g = 3)]))
# three and four traversals of one conformational class per turnover, every other step a
# binding or fused chemistry (kind 1 or 2), so the class is eliminated first: Keq^(1/3), ^(1/4)
cube(re) = [st(E(), Fc(); g = 1), st(E(:A), Fc(:A); g = 1), st(E(:A, :B), Fc(:A, :B); g = 1),
            st(Fc(), E(:A); c = [:A], re = re, g = 2), st(Fc(:A), E(:A, :B); c = [:B], re = re, g = 3),
            st(Fc(:A, :B), E(); r = [:P], g = 4)]
case("cube_root", :derive, () -> build("cube_root", [:A, :B], [:P], cube(false)))
case("cube_root_rapid", :derive, () -> build("cube_root_rapid", [:A, :B], [:P], cube(true)))
case("quarter_root", :derive, () -> build("quarter_root", [:A, :B], [:P, :Q], [
    st(E(), Fc(); g = 1), st(E(:A), Fc(:A); g = 1), st(E(:A, :B), Fc(:A, :B); g = 1),
    st(E(:Q), Fc(:Q); g = 1),
    st(Fc(), E(:A); c = [:A], g = 2), st(Fc(:A), E(:A, :B); c = [:B], g = 3),
    st(Fc(:A, :B), E(:Q); r = [:P], g = 4), st(Fc(:Q), E(); r = [:Q], g = 5)]))
# ── several Haldane dependents, Keq in Q ───────────────────────────────────────────────
two_routes(tag...) = [st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); g = 2),
                      st(E(:S), Fc(:P); g = 3), st(E(:P), E(); r = [:P], re = true, g = 4),
                      st(Fc(:P), E(); r = [:P], g = 5)]
case("two_haldane", :derive, () -> build("two_haldane", [:S], [:P], two_routes()))
case("three_haldane_ss", :derive, () -> build("three_haldane_ss", [:A], [:P, :Q], [
    st(E(), E(:A); c = [:A], g = 1), st(E(:A), E(:P, :Q); g = 2),
    st(E(:A), Fc(:P, :Q); g = 3), st(E(:A), E(:Q); r = [:P], g = 4),
    st(E(:P, :Q), E(:Q); r = [:P], g = 5), st(Fc(:P, :Q), E(:P); r = [:Q], g = 6),
    st(E(:Q), E(); r = [:Q], g = 7), st(E(:P), E(); r = [:P], g = 8)]))
# ── modifier copies, inhibitors, copies of one metabolite ──────────────────────────────
case("substrate_inhibition_mod", :derive, () -> build("substrate_inhibition_mod", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); g = 2),
    st(E(:P), E(); r = [:P], re = true, g = 3),
    st(E(:S), E(:S, Mod(:S)); c = [Mod(:S)], re = true, g = 4),
    st(E(), E(Mod(:P)); c = [Mod(:P)], re = true, g = 5)]))
case("inhibitor_two_sites", :derive, () -> build("inhibitor_two_sites", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(), E(:I); c = [:I], re = true, g = 3), st(E(:S), E(:S, :I); c = [:I], re = true, g = 4)]))
case("two_copies_inhibitor", :derive, () -> build("two_copies_inhibitor", [:S], [:P], [
    st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(), E(:I); c = [:I], re = true, g = 3), st(E(:I), E(:I, :I); c = [:I], re = true, g = 4)]))
case("copies_cap_16", :error => "15 copies", () -> build("copies_cap_16", [:S], [:P], [
    st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(), E(fill(:I, 16)...); c = fill(:I, 16), re = true, g = 3)]))
case("entries_cap_33", :error => "32 metabolite entries", () -> build("entries_cap_33",
    [:S], [:P], [[st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(); r = [:P], g = 2)];
                 [st(E(), E(Symbol("I", i)); c = [Symbol("I", i)], re = true, g = 2 + i)
                  for i in 1:31]]))
case("role_change", :reject => "role", () -> build("role_change", [:S], [:P], [
    st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(:P); g = 2), st(E(:P), E(); r = [:P], g = 3),
    st(E(), E(:P); c = [Mod(:P)], re = true, g = 4)]))
# ── rapid-equilibrium segments sharing forms, closed cycles, all RE ────────────────────
case("re_square_closed", :derive, () -> build("re_square_closed", [:A, :B], [:P], [
    st(E(), E(:A); c = [:A], re = true, g = 1), st(E(), E(:B); c = [:B], re = true, g = 2),
    st(E(:A), E(:A, :B); c = [:B], re = true, g = 3),
    st(E(:B), E(:A, :B); c = [:A], re = true, g = 4), st(E(:A, :B), E(); r = [:P], g = 5)]))
case("re_square_tied", :derive, () -> build("re_square_tied", [:A, :B], [:P], [
    st(E(), E(:A); c = [:A], re = true, g = 1), st(E(), E(:B); c = [:B], re = true, g = 2),
    st(E(:A), E(:A, :B); c = [:B], re = true, g = 2),
    st(E(:B), E(:A, :B); c = [:A], re = true, g = 1), st(E(:A, :B), E(); r = [:P], g = 3)]))
case("re_star_shared", :derive, () -> build("re_star_shared", [:A, :B], [:P], [
    st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(:A, :B); c = [:B], re = true, g = 2),
    st(E(:A), E(:A, :I); c = [:I], re = true, g = 3), st(E(:A), Fc(:A); re = true, g = 4),
    st(Fc(:A), Fc(:A, :B); c = [:B], re = true, g = 5), st(E(:A, :B), E(); r = [:P], g = 6),
    st(Fc(:A, :B), E(); r = [:P], g = 7)]))
case("all_re", :reject => "", () -> build("all_re", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); re = true, g = 2),
    st(E(:P), E(); r = [:P], re = true, g = 3)]))
case("missing_substrate", :reject => "", () -> build("missing_substrate", [:A, :B], [:P], [
    st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(); r = [:P], g = 2)]))
case("theorell_chance", :derive, () -> build("theorell_chance", [:A, :B], [:P, :Q], [
    st(E(), E(:A); c = [:A], re = true, g = 1), st(E(:A), E(:Q); c = [:B], r = [:P], g = 2),
    st(E(:Q), E(); r = [:Q], re = true, g = 3)]))
case("deadend_isomer_saturated", :derive, () -> build("deadend_isomer_saturated", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(:S), Fc(:S); re = true, g = 3)]))
case("ss_noflux_deadend", :derive, () -> build("ss_noflux_deadend", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(), E(:I); c = [:I], g = 3)]))
# ── printed names that collide ─────────────────────────────────────────────────────────
case("coef_name_collision", :derive, () -> build("coef_name_collision", [:A, :B], [:P], [
    st(E(), E(:A); c = [:A], g = 1), st(E(:A), E(:A, :B); c = [:B], g = 2),
    st(E(:A, :B), E(); r = [:P], g = 3), st(E(), E(:A_B); c = [:A_B], re = true, g = 4)]))
case("reserved_N", :derive, () -> build("reserved_N", [:S], [:P], [
    st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(:P); g = 2), st(E(:P), E(); r = [:P], g = 3),
    st(E(), E(:N); c = [:N], re = true, g = 4)]))
case("reserved_Keq", :reject => "", () -> build("reserved_Keq", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(), E(:Keq); c = [:Keq], re = true, g = 3)]))
case("reserved_E_total", :reject => "", () -> build("reserved_E_total", [:S], [:P], [
    st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(); r = [:P], g = 2),
    st(E(), E(:E_total); c = [:E_total], re = true, g = 3)]))
# ── allosteric: sites at n = 3 and 4, flip forms, islands, mixed :OnlyA groups ────────
uni3(t1, t2, t3) = [st(E(), E(:S); c = [:S], re = true, g = 1), st(E(:S), E(:P); g = 2),
                    st(E(:P), E(); r = [:P], re = true, g = 3)]
case("allo_n3_sites", :derive, () -> build("allo_n3_sites", [:S], [:P], uni3(0, 0, 0);
    tags = [:NonequalAI, :OnlyA, :EqualAI], mult = 3,
    sites = [site([:R1, :R2], [:OnlyA, :OnlyI], 4), site([:R3], [:NonequalAI], 3),
             site([:R4], [:EqualAI], 2)]))
case("allo_n4_island_binding", :derive, () -> build("allo_n4_island_binding", [:S], [:P],
    uni3(0, 0, 0); tags = [:OnlyA, :OnlyA, :EqualAI], mult = 4,
    sites = [site([:R], [:NonequalAI], 4)]))
case("allo_n1_island_binding", :derive, () -> build("allo_n1_island_binding", [:S], [:P],
    uni3(0, 0, 0); tags = [:OnlyA, :OnlyA, :EqualAI], mult = 1))
case("allo_flip_re_other_free", :derive, () -> build("allo_flip_re_other_free", [:S], [:P], [
    st(E(), Fc(); re = true, g = 1), st(Fc(), Fc(:S); c = [:S], re = true, g = 2),
    st(Fc(:S), Fc(); r = [:P], g = 3)]; tags = [:EqualAI, :NonequalAI, :NonequalAI], mult = 1))
case("allo_flip_ss_other_free", :reject => "free form", () -> build("allo_flip_ss_other_free",
    [:S], [:P], [st(E(), Fc(); g = 1), st(Fc(), Fc(:S); c = [:S], re = true, g = 2),
                 st(Fc(:S), Fc(); r = [:P], g = 3)]; tags = [:EqualAI, :NonequalAI, :NonequalAI]))
case("allo_flip_re_onlyA", :reject => "free form", () -> build("allo_flip_re_onlyA",
    [:S], [:P], [st(E(), Fc(); re = true, g = 1), st(Fc(), Fc(:S); c = [:S], re = true, g = 2),
                 st(Fc(:S), Fc(); r = [:P], g = 3)]; tags = [:OnlyA, :NonequalAI, :NonequalAI]))
pp(t) = [st(E(), E(:A); c = [:A], re = true, g = 1),
         st(E(:A), E(; add = [:A], sub = [:P]); r = [:P], g = 2),
         st(E(; add = [:A], sub = [:P]), E(:B; add = [:A], sub = [:P]); c = [:B], re = true, g = 3),
         st(E(:B; add = [:A], sub = [:P]), E(); r = [:Q], g = 4)]
case("allo_pingpong_island", :reject => "reach the free enzyme", () -> build(
    "allo_pingpong_island", [:A, :B], [:P, :Q], pp(0);
    tags = [:OnlyA, :OnlyA, :EqualAI, :OnlyA], mult = 2))
case("allo_pingpong_behind_binding", :derive, () -> build("allo_pingpong_behind_binding",
    [:A, :B], [:P, :Q], pp(0); tags = [:NonequalAI, :NonequalAI, :OnlyA, :OnlyA], mult = 1))
rand_bibi(tc) = [st(E(), E(:A); c = [:A], re = true, g = 1),
                 st(E(:B), E(:A, :B); c = [:A], re = true, g = 1),
                 st(E(), E(:B); c = [:B], re = true, g = 2),
                 st(E(:A), E(:A, :B); c = [:B], re = true, g = 2),
                 st(E(:A, :B), E(); r = [:P], g = 3)]
case("allo_mixed_onlyA_group", :derive, () -> build("allo_mixed_onlyA_group", [:A, :B], [:P],
    rand_bibi(0); tags = [:NonequalAI, :OnlyA, :OnlyA], mult = 1))
case("allo_mixed_chord", :reject => "still join", () -> build("allo_mixed_chord", [:A, :B], [:P],
    rand_bibi(0); tags = [:NonequalAI, :OnlyA, :NonequalAI], mult = 1))
case("allo_flip_alone", :derive, () -> build("allo_flip_alone", [:S], [:P], uni3(0, 0, 0);
    tags = [:OnlyA, :EqualAI, :OnlyA], mult = 2))
case("allo_equal_all_regulator", :derive, () -> build("allo_equal_all_regulator", [:S], [:P],
    uni3(0, 0, 0); tags = [:EqualAI, :EqualAI, :EqualAI], mult = 2,
    sites = [site([:R], [:OnlyA], 2)]))
case("allo_nonequal_two_haldane", :derive, () -> build("allo_nonequal_two_haldane", [:S], [:P],
    two_routes(); tags = fill(:NonequalAI, 5), mult = 1))
case("allo_onlyA_ss_binding", :derive, () -> build("allo_onlyA_ss_binding", [:S], [:P], [
    st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(:P); g = 2), st(E(:P), E(); r = [:P], g = 3)];
    tags = [:OnlyA, :OnlyA, :NonequalAI], mult = 1))
case("allo_cube_root_nonequal", :derive, () -> build("allo_cube_root_nonequal", [:A, :B],
    [:P], cube(false); tags = fill(:NonequalAI, 4), mult = 2,
    sites = [site([:R], [:NonequalAI], 3)]))
case("allo_reserved_L", :reject => "", () -> build("allo_reserved_L", [:S], [:P], uni3(0, 0, 0);
    tags = [:NonequalAI, :NonequalAI, :NonequalAI], sites = [site([:L], [:OnlyA], 1)]))
