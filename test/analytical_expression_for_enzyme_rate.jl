# ============================================================================
# Steady-State Rate Equations from King-Altman Method
# Reference: Segel, Enzyme Kinetics (Chapter IX)
# ============================================================================

# --------------------------------------------------------------------------
# 1. Uni Uni Mechanism (Segel, Eq. IX-8)
#
#        k1f       k2f       k3f
# E + A ⇌   EA  ⇌   EP  ⇌   E + P
#        k1r       k2r       k3r
# --------------------------------------------------------------------------
"""
    rate_uni_uni(params, concs)

Steady-state rate for a Uni Uni mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, Etotal)`
- `concs`:  NamedTuple `(A, P)`
"""
function rate_uni_uni(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
    (; A, P) = concs

    num = k1f * k2f * k3f * A - k1r * k2r * k3r * P

    denom = (k1r * k2r + k1r * k3f + k2f * k3f) +
            k1f * (k2r + k2f + k3f) * A +
            k3r * (k1r + k2r + k2f) * P

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 2. Iso Uni Uni Mechanism (Segel, Eq. IX-45)
#
#        k1f       k2f        k3f
# E + A ⇌   EA  ⇌   E'P  ⇌   E' + P
#        k1r       k2r        k3r
#
#               k4f
#         E' ⇌ E
#               k4r
#
# Four enzyme species: E, EA, E'P, E'
# King-Altman figure is a square: E → EA → E'P → E' → E
# --------------------------------------------------------------------------
"""
    rate_iso_uni_uni(params, concs)

Steady-state rate for an Iso Uni Uni mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal)`
- `concs`:  NamedTuple `(A, P)`
"""
function rate_iso_uni_uni(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
    (; A, P) = concs

    num = k1f * k2f * k3f * k4f * A - k1r * k2r * k3r * k4r * P

    denom = (k4f + k4r) * (k1r * k3f + k1r * k2r + k2f * k3f) +
            k1f * (k2f * k3f + k2f * k4f + k2r * k4f + k3f * k4f) * A +
            k3r * (k1r * k2r + k1r * k4r + k2f * k4r + k2r * k4r) * P +
            k1f * k3r * (k2f + k2r) * A * P

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 3. Ordered Uni Bi Mechanism (Segel, Eq. IX-60)
#
#        k1f                k2f            k3f
# E + A ⇌   (EA ⇌ EPQ)  ⇌   P + EQ  ⇌   E + Q
#        k1r                k2r            k3r
#
# Three enzyme species: E, (EA≡EPQ), EQ
# King-Altman figure is a triangle
# --------------------------------------------------------------------------
"""
    rate_ordered_uni_bi(params, concs)

Steady-state rate for an Ordered Uni Bi mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, Etotal)`
- `concs`:  NamedTuple `(A, P, Q)`
"""
function rate_ordered_uni_bi(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
    (; A, P, Q) = concs

    num = k1f * k2f * k3f * A - k1r * k2r * k3r * P * Q

    denom = k3f * (k2f + k1r) +
            k1f * (k2f + k3f) * A +
            k1r * k2r * P +
            k3r * (k2f + k1r) * Q +
            k1f * k2r * A * P +
            k2r * k3r * P * Q

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 4. Ordered Bi Bi Mechanism (Segel, Eq. IX-87)
#
#        k1f         k2f               k3f            k4f
# E + A ⇌   EA + B ⇌   (EAB ⇌ EPQ) ⇌   EQ + P  ⇌   E + Q
#        k1r         k2r               k3r            k4r
#
# Four enzyme species: E, EA, (EAB≡EPQ), EQ
# King-Altman figure is a square: E → EA → (EAB≡EPQ) → EQ → E
# --------------------------------------------------------------------------
"""
    rate_ordered_bi_bi(params, concs)

Steady-state rate for an Ordered Bi Bi mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal)`
- `concs`:  NamedTuple `(A, B, P, Q)`
"""
function rate_ordered_bi_bi(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
    (; A, B, P, Q) = concs

    num = k1f * k2f * k3f * k4f * A * B - k1r * k2r * k3r * k4r * P * Q

    denom = k1r * k4f * (k2r + k3f) +
            k1f * k4f * (k2r + k3f) * A +
            k2f * k3f * k4f * B +
            k1r * k2r * k3r * P +
            k1r * k4r * (k2r + k3f) * Q +
            k1f * k2f * (k3f + k4f) * A * B +
            k1f * k2r * k3r * A * P +
            k2f * k3f * k4r * B * Q +
            k3r * k4r * (k1r + k2r) * P * Q +
            k1f * k2f * k3r * A * B * P +
            k2f * k3r * k4r * B * P * Q

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 5. Theorell-Chance Bi Bi Mechanism (Segel, Eq. IX-122)
#
#        k1f          k2f           k3f
# E + A ⇌   EA + B ⇌   EQ + P  ⇌   E + Q
#        k1r          k2r           k3r
#
# No ternary central complex. Three enzyme species: E, EA, EQ.
# King-Altman figure is a triangle.
# --------------------------------------------------------------------------
"""
    rate_theorell_chance_bi_bi(params, concs)

Steady-state rate for a Theorell-Chance Bi Bi mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, Etotal)`
- `concs`:  NamedTuple `(A, B, P, Q)`
"""
function rate_theorell_chance_bi_bi(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
    (; A, B, P, Q) = concs

    num = k1f * k2f * k3f * A * B - k1r * k2r * k3r * P * Q

    denom = k1r * k3f +
            k1f * k3f * A +
            k2f * k3f * B +
            k1r * k2r * P +
            k1r * k3r * Q +
            k1f * k2f * A * B +
            k1f * k2r * A * P +
            k2f * k3r * B * Q +
            k2r * k3r * P * Q

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 6. Ping Pong Bi Bi Mechanism (Segel, Eq. IX-140)
#
#        k1f             k2f          k3f             k4f
# E + A ⇌   (EA ⇌ FP) ⇌   F + P, F + B ⇌   (FB ⇌ EQ) ⇌   E + Q
#        k1r             k2r          k3r             k4r
#
# E = free enzyme, F = modified (amino) enzyme form
# Four enzyme species: E, (EA≡FP), F, (FB≡EQ)
# King-Altman figure is a square: E → (EA≡FP) → F → (FB≡EQ) → E
# --------------------------------------------------------------------------
"""
    rate_ping_pong_bi_bi(params, concs)

Steady-state rate for a Ping Pong Bi Bi mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal)`
- `concs`:  NamedTuple `(A, B, P, Q)`
"""
function rate_ping_pong_bi_bi(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
    (; A, B, P, Q) = concs

    num = k1f * k2f * k3f * k4f * A * B - k1r * k2r * k3r * k4r * P * Q

    denom = k1f * k2f * (k3r + k4f) * A +
            k3f * k4f * (k1r + k2f) * B +
            k1r * k2r * (k3r + k4f) * P +
            k3r * k4r * (k1r + k2f) * Q +
            k1f * k3f * (k2f + k4f) * A * B +
            k1f * k2r * (k3r + k4f) * A * P +
            k3f * k4r * (k1r + k2f) * B * Q +
            k2r * k4r * (k1r + k3r) * P * Q

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 7. Ordered Ter Bi Mechanism (Segel, Eq. IX-195)
#
#        k1f          k2f          k3f                k4f           k5f
# E + A ⇌   EA + B ⇌   EAB + C ⇌   (EABC ⇌ EPQ) ⇌   EQ + P  ⇌   E + Q
#        k1r          k2r          k3r                k4r           k5r
#
# Five enzyme species: E, EA, EAB, (EABC≡EPQ), EQ
# King-Altman figure is a pentagon.
# --------------------------------------------------------------------------
"""
    rate_ordered_ter_bi(params, concs)

Steady-state rate for an Ordered Ter Bi mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal)`
- `concs`:  NamedTuple `(A, B, C, P, Q)`
"""
function rate_ordered_ter_bi(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal) = params
    (; A, B, C, P, Q) = concs

    num = k1f * k2f * k3f * k4f * k5f * A * B * C -
          k1r * k2r * k3r * k4r * k5r * P * Q

    denom = # constant
            k1r * k2r * k5f * (k3r + k4f) +
            # single substrate/product
            k1f * k2r * k5f * (k3r + k4f) * A +
            k1r * k3f * k4f * k5f * C +
            k1r * k2r * k3r * k4r * P +
            k1r * k2r * k5r * (k3r + k4f) * Q +
            # two substrates/products
            k1f * k2f * k5f * (k3r + k4f) * A * B +
            k1f * k3f * k4f * k5f * A * C +
            k2f * k3f * k4f * k5f * B * C +
            k1f * k2r * k3r * k4r * A * P +
            k1r * k3f * k4f * k5r * C * Q +
            k4r * k5r * (k1r * k2r + k1r * k3r + k2r * k3r) * P * Q +
            # three substrates/products
            k1f * k2f * k3f * (k4f + k5f) * A * B * C +
            k1f * k2f * k3r * k4r * A * B * P +
            k2f * k3f * k4f * k5r * B * C * Q +
            k1r * k3f * k4r * k5r * C * P * Q +
            k2f * k3r * k4r * k5r * B * P * Q +
            # four substrates/products
            k1f * k2f * k3f * k4r * A * B * C * P +
            k2f * k3f * k4r * k5r * B * C * P * Q

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 8. Bi Uni Uni Uni Ping Pong Ter Bi Mechanism (Segel, Eq. IX-228)
#
#        k1f          k2f              k3f          k4f              k5f
# E + A ⇌   EA + B ⇌   (EAB ⇌ FP) ⇌   F + P, F + C ⇌   (FC ⇌ EQ) ⇌   E + Q
#        k1r          k2r              k3r          k4r              k5r
#
# E = free enzyme, F = modified enzyme form after first half-reaction
# Five enzyme species: E, EA, (EAB≡FP), F, (FC≡EQ)
# King-Altman figure is a pentagon: E → EA → (EAB≡FP) → F → (FC≡EQ) → E
# --------------------------------------------------------------------------
"""
    rate_bi_uni_uni_uni_ping_pong_ter_bi(params, concs)

Steady-state rate for a Bi Uni Uni Uni Ping Pong Ter Bi mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal)`
- `concs`:  NamedTuple `(A, B, C, P, Q)`
"""
function rate_bi_uni_uni_uni_ping_pong_ter_bi(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal) = params
    (; A, B, C, P, Q) = concs

    num = k1f * k2f * k3f * k4f * k5f * A * B * C -
          k1r * k2r * k3r * k4r * k5r * P * Q

    denom = # single substrate/product
            k1r * k4f * k5f * (k2r + k3f) * C +
            k1r * k2r * k3r * (k4r + k5f) * P +
            k1r * k4r * k5r * (k2r + k3f) * Q +
            # two substrates/products
            k1f * k2f * k3f * (k4r + k5f) * A * B +
            k1f * k4f * k5f * (k2r + k3f) * A * C +
            k1f * k2r * k3r * (k4r + k5f) * A * P +
            k2f * k3f * k4f * k5f * B * C +
            k2f * k3f * k4r * k5r * B * Q +
            k1r * k4f * k5r * (k2r + k3f) * C * Q +
            k3r * k5r * (k1r * k2r + k1r * k4r + k2r * k4r) * P * Q +
            # three substrates/products
            k1f * k2f * k4f * (k3f + k5f) * A * B * C +
            k1f * k2f * k3r * (k4r + k5f) * A * B * P +
            k2f * k3f * k4f * k5r * B * C * Q +
            k2f * k3r * k4r * k5r * B * P * Q

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 9. Ordered Ter Ter Mechanism (Segel, Eq. IX-261)
#
#        k1f          k2f          k3f                  k4f           k5f           k6f
# E + A ⇌   EA + B ⇌   EAB + C ⇌   (EABC ⇌ EPQR) ⇌   EQR + P  ⇌   ER + Q  ⇌   E + R
#        k1r          k2r          k3r                  k4r           k5r           k6r
#
# Six enzyme species: E, EA, EAB, (EABC≡EPQR), EQR, ER
# King-Altman figure is a hexagon.
# Denominator contains 27 grouped terms.
# --------------------------------------------------------------------------
"""
    rate_ordered_ter_ter(params, concs)

Steady-state rate for an Ordered Ter Ter mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal)`
- `concs`:  NamedTuple `(A, B, C, P, Q, R)`
"""
function rate_ordered_ter_ter(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
    (; A, B, C, P, Q, R) = concs

    num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
          k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R

    denom = # constant
            k1r * k2r * k5f * k6f * (k3r + k4f) +
            # single substrate/product
            k1f * k2r * k5f * k6f * (k3r + k4f) * A +
            k1r * k3f * k4f * k5f * k6f * C +
            k1r * k2r * k3r * k4r * k6f * P +
            k1r * k2r * k5f * k6r * (k3r + k4f) * R +
            # two substrates/products
            k1f * k2f * k5f * k6f * (k3r + k4f) * A * B +
            k1f * k3f * k4f * k5f * k6f * A * C +
            k1f * k2r * k3r * k4r * k6f * A * P +
            k2f * k3f * k4f * k5f * k6f * B * C +
            k1r * k3f * k4f * k5f * k6r * C * R +
            k1r * k2r * k3r * k4r * k5r * P * Q +
            k1r * k2r * k3r * k4r * k6r * P * R +
            k1r * k2r * k5r * k6r * (k3r + k4f) * Q * R +
            # three substrates/products
            k1f * k2f * k3f * (k4f * k5f + k4f * k6f + k5f * k6f) * A * B * C +
            k1f * k2f * k3r * k4r * k6f * A * B * P +
            k1f * k2r * k3r * k4r * k5r * A * P * Q +
            k2f * k3f * k4f * k5f * k6r * B * C * R +
            k1r * k3f * k4f * k5r * k6r * C * Q * R +
            k4r * k5r * k6r * (k1r * k2r + k1r * k3r + k2r * k3r) * P * Q * R +
            # four substrates/products
            k1f * k2f * k3f * k4r * k6f * A * B * C * P +
            k1f * k2f * k3f * k4f * k5r * A * B * C * Q +
            k1f * k2f * k3r * k4r * k5r * A * B * P * Q +
            k2f * k3f * k4f * k5r * k6r * B * C * Q * R +
            k2f * k3r * k4r * k5r * k6r * B * P * Q * R +
            k1r * k3f * k4r * k5r * k6r * C * P * Q * R +
            # five substrates/products
            k1f * k2f * k3f * k4r * k5r * A * B * C * P * Q +
            k2f * k3f * k4r * k5r * k6r * B * C * P * Q * R

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 10. Bi Uni Uni Bi Ping Pong Ter Ter Mechanism (Segel, Eq. IX-278)
#
#        k1f          k2f              k3f           k4f              k5f           k6f
# E + A ⇌   EA + B ⇌   (EAB ⇌ FP) ⇌   F + P, F + C ⇌   (FC ⇌ EQR) ⇌   ER + Q  ⇌   E + R
#        k1r          k2r              k3r           k4r              k5r           k6r
#
# Six enzyme species: E, EA, (EAB≡FP), F, (FC≡EQR), ER
# King-Altman figure is a hexagon.
# --------------------------------------------------------------------------
"""
    rate_bi_uni_uni_bi_ping_pong_ter_ter(params, concs)

Steady-state rate for a Bi Uni Uni Bi Ping Pong Ter Ter mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal)`
- `concs`:  NamedTuple `(A, B, C, P, Q, R)`
"""
function rate_bi_uni_uni_bi_ping_pong_ter_ter(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
    (; A, B, C, P, Q, R) = concs

    num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
          k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R

    denom = # single substrate/product
            k1r * k4f * k5f * k6f * (k2r + k3f) * C +
            k1r * k2r * k3r * k6f * (k4r + k5f) * P +
            # two substrates/products
            k1f * k2f * k3f * k6f * (k4r + k5f) * A * B +
            k1f * k4f * k5f * k6f * (k2r + k3f) * A * C +
            k1f * k2r * k3r * k6f * (k4r + k5f) * A * P +
            k2f * k3f * k4f * k5f * k6f * B * C +
            k1r * k4f * k5f * k6r * (k2r + k3f) * C * R +
            k1r * k2r * k3r * k4r * k5r * P * Q +
            k1r * k2r * k3r * k6r * (k4r + k5f) * P * R +
            k1r * k4r * k5r * k6r * (k2r + k3f) * Q * R +
            # three substrates/products
            k1f * k2f * k4f * (k3f * k5f + k3f * k6f + k5f * k6f) * A * B * C +
            k1f * k2f * k3r * k6f * (k4r + k5f) * A * B * P +
            k1f * k2f * k3f * k4r * k5r * A * B * Q +
            k1f * k2r * k3r * k4r * k5r * A * P * Q +
            k2f * k3f * k4f * k5f * k6r * B * C * R +
            k2f * k3f * k4r * k5r * k6r * B * Q * R +
            k1r * k4f * k5r * k6r * (k2r + k3f) * C * Q * R +
            k3r * k5r * k6r * (k1r * k2r + k1r * k4r + k2r * k4r) * P * Q * R +
            # four substrates/products
            k1f * k2f * k3f * k4f * k5r * A * B * C * Q +
            k1f * k2f * k3r * k4r * k5r * A * B * P * Q +
            k2f * k3f * k4f * k5r * k6r * B * C * Q * R +
            k2f * k3r * k4r * k5r * k6r * B * P * Q * R

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 11. Bi Bi Uni Uni Ping Pong Ter Ter Mechanism (Segel, Eq. IX-288)
#
#        k1f          k2f               k3f           k4f          k5f             k6f
# E + A ⇌   EA + B ⇌   (EAB ⇌ FPQ) ⇌   FQ + P, FQ ⇌   F + Q, F + C ⇌   (FC ⇌ ER) ⇌   E + R
#        k1r          k2r               k3r           k4r          k5r             k6r
#
# Six enzyme species: E, EA, (EAB≡FPQ), FQ, F, (FC≡ER)
# King-Altman figure is a hexagon.
# --------------------------------------------------------------------------
"""
    rate_bi_bi_uni_uni_ping_pong_ter_ter(params, concs)

Steady-state rate for a Bi Bi Uni Uni Ping Pong Ter Ter mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal)`
- `concs`:  NamedTuple `(A, B, C, P, Q, R)`
"""
function rate_bi_bi_uni_uni_ping_pong_ter_ter(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
    (; A, B, C, P, Q, R) = concs

    num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
          k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R

    denom = # single substrate/product
            k1r * k4f * k5f * k6f * (k2r + k3f) * C +
            k1r * k4f * k5r * k6r * (k2r + k3f) * R +
            # two substrates/products
            k1f * k2f * k3f * k4f * (k5r + k6f) * A * B +
            k1f * k4f * k5f * k6f * (k2r + k3f) * A * C +
            k2f * k3f * k4f * k5f * k6f * B * C +
            k2f * k3f * k4f * k5r * k6r * B * R +
            k1r * k2r * k3r * k5f * k6f * C * P +
            k1r * k4f * k5f * k6r * (k2r + k3f) * C * R +
            k1r * k2r * k3r * k4r * (k5r + k6f) * P * Q +
            k1r * k2r * k3r * k5r * k6r * P * R +
            k1r * k4r * k5r * k6r * (k2r + k3f) * Q * R +
            # three substrates/products
            k1f * k2f * k5f * (k3f * k4f + k3f * k6f + k4f * k6f) * A * B * C +
            k1f * k2f * k3f * k4r * (k5r + k6f) * A * B * Q +
            k1f * k2r * k3r * k5f * k6f * A * C * P +
            k1f * k2r * k3r * k4r * (k5r + k6f) * A * P * Q +
            k2f * k3f * k4f * k5f * k6r * B * C * R +
            k2f * k3f * k4r * k5r * k6r * B * Q * R +
            k1r * k2r * k3r * k5f * k6r * C * P * R +
            k3r * k4r * k6r * (k1r * k2r + k1r * k5r + k2r * k5r) * P * Q * R +
            # four substrates/products
            k1f * k2f * k3r * k5f * k6f * A * B * C * P +
            k1f * k2f * k3r * k4r * (k5r + k6f) * A * B * P * Q +
            k2f * k3r * k4r * k5r * k6r * B * P * Q * R

    return Etotal * num / denom
end

# --------------------------------------------------------------------------
# 12. Hexa Uni Ping Pong Mechanism (Segel, Eq. IX-308)
#
#        k1f              k2f          k3f              k4f          k5f              k6f
# E + A ⇌   (EA ⇌ FP) ⇌   F + P, F + B ⇌   (FB ⇌ GQ) ⇌   G + Q, G + C ⇌   (GC ⇌ ER) ⇌   E + R
#        k1r              k2r          k3r              k4r          k5r              k6r
#
# E = free enzyme, F = first modified form, G = second modified form
# Six enzyme species: E, (EA≡FP), F, (FB≡GQ), G, (GC≡ER)
# King-Altman figure is a hexagon. Denominator has 17 grouped terms.
# --------------------------------------------------------------------------
"""
    rate_hexa_uni_ping_pong(params, concs)

Steady-state rate for a Hexa Uni Ping Pong mechanism.

# Arguments
- `params`: NamedTuple `(k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal)`
- `concs`:  NamedTuple `(A, B, C, P, Q, R)`
"""
function rate_hexa_uni_ping_pong(params, concs)
    (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
    (; A, B, C, P, Q, R) = concs

    num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
          k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R

    denom = # two substrates/products (9 classes, each with 2 terms → factored)
            k1f * k2f * k3f * k4f * (k5r + k6f) * A * B +
            k1f * k2f * k5f * k6f * (k3r + k4f) * A * C +
            k1f * k2f * k3r * k4r * (k5r + k6f) * A * Q +
            k3f * k4f * k5f * k6f * (k1r + k2f) * B * C +
            k3f * k4f * k5r * k6r * (k1r + k2f) * B * R +
            k1r * k2r * k5f * k6f * (k3r + k4f) * C * P +
            k1r * k2r * k3r * k4r * (k5r + k6f) * P * Q +
            k1r * k2r * k5r * k6r * (k3r + k4f) * P * R +
            k3r * k4r * k5r * k6r * (k1r + k2f) * Q * R +
            # three substrates/products (8 classes)
            k1f * k3f * k5f * (k2f * k4f + k2f * k6f + k4f * k6f) * A * B * C +
            k1f * k2f * k3f * k4r * (k5r + k6f) * A * B * Q +
            k1f * k2r * k5f * k6f * (k3r + k4f) * A * C * P +
            k1f * k2r * k3r * k4r * (k5r + k6f) * A * P * Q +
            k3f * k4f * k5f * k6r * (k1r + k2f) * B * C * R +
            k3f * k4r * k5r * k6r * (k1r + k2f) * B * Q * R +
            k1r * k2r * k5f * k6r * (k3r + k4f) * C * P * R +
            k2r * k4r * k6r * (k1r * k3r + k1r * k5r + k3r * k5r) * P * Q * R

    return Etotal * num / denom
end
