# Unit tests for symbolic polynomial types (FactoredPoly, FactoredSigma, DenomTerm)

@testset "Factored denominator types" begin
    # Aliases for internal constructors/functions
    FP = EnzymeRates.FactoredPoly
    FS = EnzymeRates.FactoredSigma
    DT = EnzymeRates.DenomTerm

    @testset "unfactored_denom_term wrapping" begin
        # A simple sigma = 1 + S/K1 (2 terms), cofactor = poly_one()
        sigma = EnzymeRates.poly_add(
            EnzymeRates.poly_one(),
            EnzymeRates.poly_sym(:S),
        )
        dt = EnzymeRates.unfactored_denom_term(sigma, EnzymeRates.poly_one())
        @test dt isa DT
        @test length(dt.sigma.coefficients) == 1
        @test dt.sigma.coefficients[1] == EnzymeRates.poly_one()
        @test length(dt.sigma.products) == 1
        @test length(dt.sigma.products[1].factors) == 1
        @test dt.sigma.products[1].factors[1] == sigma
        @test dt.sigma.products[1].exponents == [1]
        @test dt.cofactor == EnzymeRates.poly_one()
    end

    @testset "_expand_to_poly matches manual expansion" begin
        # Test: (1 + S)^2 * (1 + R) = 1 + 2S + S^2 + R + 2SR + S^2*R
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)
        f1 = EnzymeRates.poly_add(one, pS)  # 1 + S
        f2 = EnzymeRates.poly_add(one, pR)  # 1 + R

        fp = FP([f1, f2], [2, 1])
        fs = FS([one], [fp])
        dt = DT(fs, one)

        expanded = EnzymeRates._expand_to_poly([dt])

        # Manual expansion: (1+S)^2 * (1+R)
        manual = EnzymeRates.poly_mul(
            EnzymeRates.poly_mul(f1, f1),
            f2,
        )
        @test expanded == manual
        @test length(expanded) == 6  # 1, S, S^2, R, SR, S^2*R
    end

    @testset "_expand_to_poly with multiple DenomTerms" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)
        pK = EnzymeRates.poly_sym(:K)

        # Term 1: 1 * (1+S) with cofactor (1+R)
        fp1 = FP([EnzymeRates.poly_add(one, pS)], [1])
        fs1 = FS([one], [fp1])
        dt1 = DT(fs1, EnzymeRates.poly_add(one, pR))

        # Term 2: K * (1+S) with cofactor poly_one()
        fp2 = FP([EnzymeRates.poly_add(one, pS)], [1])
        fs2 = FS([pK], [fp2])
        dt2 = DT(fs2, one)

        expanded = EnzymeRates._expand_to_poly([dt1, dt2])

        # Manual: (1+S)*(1+R) + K*(1+S)
        manual = EnzymeRates.poly_add(
            EnzymeRates.poly_mul(
                EnzymeRates.poly_add(one, pS),
                EnzymeRates.poly_add(one, pR),
            ),
            EnzymeRates.poly_mul(pK, EnzymeRates.poly_add(one, pS)),
        )
        @test expanded == manual
    end

    @testset "_estimate_expanded_term_count" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)

        # (1+S)^2 * (1+R)^1 with cofactor 1
        # estimate = 1 * (2^2 * 2^1) * 1 = 8
        f1 = EnzymeRates.poly_add(one, pS)  # 2 terms
        f2 = EnzymeRates.poly_add(one, pR)  # 2 terms
        fp = FP([f1, f2], [2, 1])
        fs = FS([one], [fp])
        dt = DT(fs, one)
        @test EnzymeRates._estimate_expanded_term_count([dt]) == 8

        # Actual expanded = 6 (due to combining S*S = S^2)
        actual = EnzymeRates._expand_to_poly([dt])
        @test length(actual) == 6
        # Estimate >= actual (upper bound)
        @test EnzymeRates._estimate_expanded_term_count([dt]) >= length(actual)
    end

    @testset "_estimate_expanded_term_count with sum-of-products" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)
        pK = EnzymeRates.poly_sym(:K)

        # FactoredSigma: 1*(1+S)^2 + K*(1+R)^2, cofactor = 1
        f1 = EnzymeRates.poly_add(one, pS)
        f2 = EnzymeRates.poly_add(one, pR)
        fp1 = FP([f1], [2])
        fp2 = FP([f2], [2])
        fs = FS([one, pK], [fp1, fp2])
        dt = DT(fs, one)

        # Estimate: 1*2^2 + 1*2^2 = 4 + 4 = 8
        @test EnzymeRates._estimate_expanded_term_count([dt]) == 8
    end

    @testset "Expr generation for factored types" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)
        pK1 = EnzymeRates.poly_sym(:K1)
        ps = Set([:K1])
        cs = Set([:S, :R])
        inv = Set{Symbol}()

        # Simple factored: (1 + K1*S) with exp 1
        f = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK1, pS))
        fp = FP([f], [1])

        expr = EnzymeRates._factored_poly_to_expr(fp, ps, cs, inv)
        @test expr isa Union{Symbol, Expr, Int}

        # With exponent 2: ((1 + K1*S))^2
        fp2 = FP([f], [2])
        expr2 = EnzymeRates._factored_poly_to_expr(fp2, ps, cs, inv)
        @test expr2 isa Expr

        # FactoredSigma with 1 sub-group
        fs = FS([one], [fp])
        sigma_expr = EnzymeRates._factored_sigma_to_expr(fs, ps, cs, inv)
        @test sigma_expr isa Union{Symbol, Expr, Int}

        # DenomTerm with cofactor = poly_one()
        dt = DT(fs, one)
        den_expr = EnzymeRates._denom_terms_to_expr([dt], ps, cs, inv)
        @test den_expr isa Union{Symbol, Expr, Int}
    end

    @testset "Expr generation: numerical equivalence" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)
        pK1 = EnzymeRates.poly_sym(:K1)
        pK2 = EnzymeRates.poly_sym(:K2)
        ps = Set([:K1, :K2])
        cs = Set([:S, :R])
        inv = Set{Symbol}()

        # Build: (1 + K1*S)^2 * (1 + K2*R)
        f1 = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK1, pS))
        f2 = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK2, pR))
        fp = FP([f1, f2], [2, 1])
        fs = FS([one], [fp])
        dt = DT(fs, one)

        # Factored Expr
        den_expr = EnzymeRates._denom_terms_to_expr([dt], ps, cs, inv)

        # Expanded Expr
        expanded = EnzymeRates._expand_to_poly([dt])
        flat_expr = EnzymeRates._poly_to_expr(expanded, ps, cs, inv)

        # Evaluate both at specific values
        vals = Dict(:K1 => 2.0, :K2 => 3.0, :S => 0.5, :R => 1.5)
        eval_expr = (e) -> begin
            bindings = ["$k = $(vals[k])" for k in keys(vals)]
            code = "let $(join(bindings, ", "))\n  $e\nend"
            eval(Meta.parse(code))
        end

        v_factored = eval_expr(den_expr)
        v_flat = eval_expr(flat_expr)
        @test isapprox(v_factored, v_flat; rtol=1e-12)

        # Also verify against manual calculation:
        # (1 + 2*0.5)^2 * (1 + 3*1.5) = (2)^2 * (5.5) = 22.0
        @test isapprox(v_factored, 22.0; rtol=1e-12)
    end

    @testset "Expr generation: sum-of-products numerical equivalence" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pR = EnzymeRates.poly_sym(:R)
        pK1 = EnzymeRates.poly_sym(:K1)
        pK2 = EnzymeRates.poly_sym(:K2)
        pKc = EnzymeRates.poly_sym(:Kc)
        ps = Set([:K1, :K2, :Kc])
        cs = Set([:S, :R])
        inv = Set{Symbol}()

        # Build: 1*(1+K1*S)*(1+K2*R) + Kc*(1+K1*S)
        f1 = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK1, pS))
        f2 = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK2, pR))

        fp1 = FP([f1, f2], [1, 1])
        fp2 = FP([f1], [1])
        fs = FS([one, pKc], [fp1, fp2])
        dt = DT(fs, one)

        den_expr = EnzymeRates._denom_terms_to_expr([dt], ps, cs, inv)
        expanded = EnzymeRates._expand_to_poly([dt])
        flat_expr = EnzymeRates._poly_to_expr(expanded, ps, cs, inv)

        vals = Dict(:K1 => 2.0, :K2 => 3.0, :Kc => 0.7, :S => 0.5, :R => 1.5)
        eval_expr = (e) -> begin
            bindings = ["$k = $(vals[k])" for k in keys(vals)]
            code = "let $(join(bindings, ", "))\n  $e\nend"
            eval(Meta.parse(code))
        end

        v_factored = eval_expr(den_expr)
        v_flat = eval_expr(flat_expr)
        @test isapprox(v_factored, v_flat; rtol=1e-12)

        # Manual: 1*(1+2*0.5)*(1+3*1.5) + 0.7*(1+2*0.5) = 2*5.5 + 0.7*2 = 12.4
        @test isapprox(v_factored, 12.4; rtol=1e-12)
    end

    @testset "to_rate_expr overload with DenomTerms" begin
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pK1 = EnzymeRates.poly_sym(:K1)
        ps = Set([:K1])
        cs = Set([:S])
        inv = Set{Symbol}()

        num = EnzymeRates.poly_mul(pK1, pS)  # K1*S
        f = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK1, pS))
        fp = FP([f], [1])
        fs = FS([one], [fp])
        dt = DT(fs, one)

        expr = EnzymeRates.to_rate_expr(num, [dt], ps, cs, inv)
        @test expr isa Expr
        @test string(expr) |> s -> occursin("E_total", s)
    end

    @testset "Symbol renaming on factored types" begin
        one = EnzymeRates.poly_one()
        pK1 = EnzymeRates.poly_sym(:K1)
        pK2 = EnzymeRates.poly_sym(:K2)
        pS = EnzymeRates.poly_sym(:S)

        # Rename map: K2 → K1 (kinetic-group alias)
        rename = Dict(:K2 => :K1)

        # FactoredPoly with K2*S in a factor
        f = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK2, pS))
        fp = FP([f], [2])
        fp_c = EnzymeRates._rename_symbols(fp, rename)
        # After rename: factor should contain K1*S instead of K2*S
        expected_f = EnzymeRates.poly_add(one, EnzymeRates.poly_mul(pK1, pS))
        @test fp_c.factors[1] == expected_f
        @test fp_c.exponents == [2]

        # FactoredSigma with K2 as coefficient
        fs = FS([pK2], [fp])
        fs_c = EnzymeRates._rename_symbols(fs, rename)
        @test fs_c.coefficients[1] == pK1

        # DenomTerm: check cofactor also gets renamed
        cofactor = EnzymeRates.poly_mul(pK2, pS)
        dt = DT(fs, cofactor)
        dt_c = EnzymeRates._rename_symbols(dt, rename)
        expected_cof = EnzymeRates.poly_mul(pK1, pS)
        @test dt_c.cofactor == expected_cof
    end

    @testset "unfactored_denom_term roundtrip: expand matches original" begin
        # Create a sigma poly, wrap it, expand it back, verify identity
        one = EnzymeRates.poly_one()
        pS = EnzymeRates.poly_sym(:S)
        pP = EnzymeRates.poly_sym(:P)
        pK1 = EnzymeRates.poly_sym(:K1)
        sigma = EnzymeRates.poly_add(
            one,
            EnzymeRates.poly_add(
                EnzymeRates.poly_mul(pK1, pS),
                EnzymeRates.poly_mul(pK1, pP),
            ),
        )
        cofactor = EnzymeRates.poly_add(one, pS)

        dt = EnzymeRates.unfactored_denom_term(sigma, cofactor)
        expanded = EnzymeRates._expand_to_poly([dt])
        manual = EnzymeRates.poly_mul(sigma, cofactor)
        @test expanded == manual
    end
end
