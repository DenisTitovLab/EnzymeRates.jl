# ABOUTME: Cleland's textbook rate equations and Haldane relations as Julia functions of named constants.
# ABOUTME: Shared by the P4 bridge and census scripts.
# ─── Textbook equations (Cleland 1963; Segel ch. IX) ───────────────────
uniuni(k, x) = (k.V1 * x.A / k.Ka - k.V2 * x.P / k.Kp) / (1 + x.A / k.Ka + x.P / k.Kp)
uniuni_v1v2(k, x) = k.V1 * k.V2 * (x.A - x.P / k.Keq) / (k.V2 * k.Ka + k.V2 * x.A + k.V1 * x.P / k.Keq)
ordered_core(k, x) = k.V2 * k.Kia * k.Kb + k.V2 * k.Kb * x.A + k.V2 * k.Ka * x.B +
    k.V1 * k.Kq * x.P / k.Keq + k.V1 * k.Kp * x.Q / k.Keq + k.V2 * x.A * x.B +
    k.V1 * k.Kq * x.A * x.P / (k.Kia * k.Keq) + k.V2 * k.Ka * x.B * x.Q / k.Kiq + k.V1 * x.P * x.Q / k.Keq
ordered(k, x) = k.V1 * k.V2 * (x.A * x.B - x.P * x.Q / k.Keq) /
    (ordered_core(k, x) + k.V2 * x.A * x.B * x.P / k.Kip + k.V1 * x.B * x.P * x.Q / (k.Kib * k.Keq))
theorell_chance(k, x) = k.V1 * k.V2 * (x.A * x.B - x.P * x.Q / k.Keq) / ordered_core(k, x)
pingpong(k, x) = k.V1 * k.V2 * (x.A * x.B - x.P * x.Q / k.Keq) /
    (k.V2 * k.Kb * x.A + k.V2 * k.Ka * x.B + k.V1 * k.Kq * x.P / k.Keq + k.V1 * k.Kp * x.Q / k.Keq +
     k.V2 * x.A * x.B + k.V1 * k.Kq * x.A * x.P / (k.Kia * k.Keq) + k.V2 * k.Ka * x.B * x.Q / k.Kiq +
     k.V1 * x.P * x.Q / k.Keq)
re_random(k, x) = (k.V1 * x.A * x.B / (k.Kia * k.Kb) - k.V2 * x.P * x.Q / (k.Kip * k.Kq)) /
    (1 + x.A / k.Kia + x.B / k.Kib + x.A * x.B / (k.Kia * k.Kb) + x.P / k.Kip + x.Q / k.Kiq +
     x.P * x.Q / (k.Kip * k.Kq))
haldanes = [
    "V1*Kp*Kiq/(V2*Kia*Kb)" => k -> k.V1 * k.Kp * k.Kiq / (k.V2 * k.Kia * k.Kb),
    "(V1/V2)^2*Kip*Kq/(Ka*Kib)" => k -> (k.V1 / k.V2)^2 * k.Kip * k.Kq / (k.Ka * k.Kib),
    "V1*Kip*Kq/(V2*Kia*Kb)" => k -> k.V1 * k.Kip * k.Kq / (k.V2 * k.Kia * k.Kb),
    "(V1/V2)^2*Kip*Kiq/(Kia*Kib)" => k -> (k.V1 / k.V2)^2 * k.Kip * k.Kiq / (k.Kia * k.Kib),
    "V1*Kp*Kiq/(V2*Ka*Kib)" => k -> k.V1 * k.Kp * k.Kiq / (k.V2 * k.Ka * k.Kib),
    "(V1/V2)^2*Kp*Kq/(Ka*Kb)" => k -> (k.V1 / k.V2)^2 * k.Kp * k.Kq / (k.Ka * k.Kb),
    "V1*Kp/(V2*Ka) (uni-uni)" => k -> k.V1 * k.Kp / (k.V2 * k.Ka)]

