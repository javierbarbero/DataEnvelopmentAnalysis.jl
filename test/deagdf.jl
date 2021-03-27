# Tests for Generalized DF DEA Models
@testset "GeneralizedDFDEAModel" begin

    ## Test Generalized DF DEA Model with Zofío and Prieto (2006) data
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]

    # alpha = 0.5 CRS equals Input Oriented CRS
    deaio = dea(X, Y, orient = :Input, rts = :CRS)
    deagdf05crs = deagdf(X, Y, alpha = 0.5, rts = :CRS)

    @test typeof(deagdf05crs) == GeneralizedDFDEAModel

    @test efficiency(deaio) ≈ efficiency(deagdf05crs) atol = 1e-7

    @test nobs(deagdf05crs) == 5
    @test ninputs(deagdf05crs) == 2
    @test noutputs(deagdf05crs) == 2

    @test efficiency(deagdf(targets(deagdf05crs, :X), targets(deagdf05crs, :Y), alpha = 0.5, rts = :CRS, slack = false)) ≈ ones(5)
    @test efficiency(deaadd(targets(deagdf05crs, :X), targets(deagdf05crs, :Y))) ≈ zeros(5) atol=1e-8

    # alphpa = 0.5 VRS
    deagdf05vrs = deagdf(X, Y, alpha = 0.5, rts = :VRS)

    @test nobs(deagdf05vrs) == 5
    @test ninputs(deagdf05vrs) == 2
    @test noutputs(deagdf05vrs) == 2
    @test efficiency(deagdf05vrs) ≈ [0.682;
                               1;
                               1;
                               0.250;
                               0.360] atol = 1e-3
    @test slacks(deagdf05vrs, :X) ≈ [0.605935 0;
                                     0 0;
                                     0 0;
                                     0 0;
                                     0.2 3.4] atol = 1e-5

    @test slacks(deagdf05vrs, :Y) ≈ [0 4.67865;
                                     0 0;
                                     0 0;
                                     0 0;
                                     3.0 0] atol = 1e-5

    @test efficiency(deagdf(targets(deagdf05vrs, :X), targets(deagdf05vrs, :Y), alpha = 0.5, rts = :CRS, slack = false)) ≈ ones(5)
    @test efficiency(deaadd(targets(deagdf05vrs, :X), targets(deagdf05vrs, :Y))) ≈ zeros(5) atol=1e-6

    # Test no slacks
    deagdfnoslack = deagdf(X, Y, alpha = 0.5, rts = :VRS, slack = false)
    @test efficiency(deagdfnoslack) == efficiency(deagdf05vrs)
    @test isempty(slacks(deagdfnoslack, :X)) == 1
    @test isempty(slacks(deagdfnoslack, :Y)) == 1

    @test efficiency(deagdf(targets(deagdfnoslack, :X), targets(deagdfnoslack, :Y), alpha = 0.5, rts = :VRS, slack = false)) ≈ ones(5)
    @test efficiency(deaadd(targets(deagdfnoslack, :X), targets(deagdfnoslack, :Y))) != zeros(5) # Different as there is no slacks in first model

    # Test default alpha is 0.5
    @test efficiency(deagdf(X, Y, rts = :CRS)) == efficiency(deagdf(X, Y, alpha = 0.5, rts = :CRS))

    ## Test if one-by-one DEA using evaluation and reference sets match initial results
    deagdf05crs_ref_eff = zeros(size(X, 1))

    deagdf05vrs_ref_eff = zeros(size(X, 1))

    Xref = X[:,:]
    Yref = Y[:,:]

    for i = 1:size(X, 1)
        Xeval = X[i:i,:]
        Xeval = Xeval[:,:]
        Yeval = Y[i:i,:]
        Yeval = Yeval[:,:]

        deagdf05crs_ref_eff[i] = efficiency(deagdf(Xeval, Yeval, alpha = 0.5, rts = :CRS, Xref = Xref, Yref = Yref, slack = false))[1]

        deagdf05vrs_ref_eff[i] = efficiency(deagdf(Xeval, Yeval, alpha = 0.5, rts = :VRS, Xref = Xref, Yref = Yref, slack = false))[1]
    end

    @test deagdf05crs_ref_eff ≈ efficiency(deagdf05crs)
    @test deagdf05vrs_ref_eff ≈ efficiency(deagdf05vrs)

    # Test with another data
    X = [2; 4; 8; 12; 6; 14; 14; 9.412];
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

    gdf2 = deagdf(X, Y, alpha = 0.7, rts = :VRS)

    @test efficiency(gdf2) ≈ [1; 1; 1; 1; 0.423216; 0.69836; 1; 0.2315] atol = 1e-4
    @test slacks(gdf2, :X) ≈ [0; 0; 0; 0; 0; 0.570479; 2.0; 0]  atol = 1e-5
    @test slacks(gdf2, :Y) ≈ [0; 0; 0; 0; 0; 0; 0; 0] atol = 1e-6

    # Print
    show(IOBuffer(), deagdf05crs)
    show(IOBuffer(), deagdfnoslack)

    # Test errors
    @test_throws DimensionMismatch deagdf([1; 2 ; 3], [4 ; 5], alpha = 0.5) #  Different number of observations
    @test_throws DimensionMismatch deagdf([1; 2], [4 ; 5], alpha = 0.5, Xref = [1; 2; 3; 4]) # Different number of observations in reference sets
    @test_throws DimensionMismatch deagdf([1 1; 2 2], [4 4; 5 5], alpha = 0.5, Xref = [1 1 1; 2 2 2]) # Different number of inputs
    @test_throws DimensionMismatch deagdf([1 1; 2 2], [4 4; 5 5], alpha = 0.5, Yref = [4 4 4; 5 5 5]) # Different number of inputs
    @test_throws ArgumentError deagdf([1; 2; 3], [4; 5; 6], alpha = 0.5, rts = :Error) # Invalid returns to scale

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deagdf(X, Y, rts = :VRS, slack = false)) ≈ [1; 1; 1; 1; 1; 1; 1; 1] atol = 1e-5

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]

    @test efficiency(deagdf(X, Y, rts = :VRS, slack = false)) ≈ [1; 1; 1; 1; 1; 1; 1; 1] atol = 1e-5

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]

    @test efficiency(deagdf(X, Y, rts = :VRS)) ≈ [1; 1; 1; 1; 0.410097; 0.634489; 1; 0.20504] atol = 1e-5

end
