# Tests for Revenue DEA Models
@testset "RevenueDEAModel" begin

    ## Test Revenue DEA Model with Cooper et al. (2007)
    # Test agains book results
    X = [3 2; 1 3; 4 6]
    Y = [3; 5; 6]
    P = [6; 6; 6]

    dearevenuecooper = dearevenue(X, Y, P, rts = :CRS)

    @test typeof(dearevenuecooper) == RevenueDEAModel
    
    @test nobs(dearevenuecooper) == 3
    @test ninputs(dearevenuecooper) == 2
    @test noutputs(dearevenuecooper) == 1
    @test ismonetary(dearevenuecooper) == false

    @test efficiency(dearevenuecooper, :Economic)   ≈ [0.9; 1; 0.6] atol = 1e-3
    @test efficiency(dearevenuecooper, :Technical)  ≈ [0.9; 1; 0.6  ] atol = 1e-3
    @test efficiency(dearevenuecooper, :Allocative) ≈ [1; 1; 1] atol = 1e-3

    ## Test Revenue DEA Model with Zofío and Prieto (2006) data.
    # Test agains results in R
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]
    P = [3 2; 3 2; 3 2; 3 2; 3 2.0]

    # Revnue CRS
    dearevenuecrs = dearevenue(X, Y, P, rts = :CRS)
    @test efficiency(dearevenuecrs) ≈ [0.4915254237;
                                    1.000;
                                    1.000;
                                    0.250;
                                    0.1735537190]
    @test efficiency(dearevenuecrs, :Technical) ≈ [0.6364;
                                1.000;
                                1.000;
                                0.250;
                                0.2609]  atol = 1e-3
    @test efficiency(dearevenuecrs, :Allocative) ≈ [0.7723970944;
                                1.000;
                                1.000;
                                1.000;
                                0.6652892562]

    @test efficiency(dearevenue(targets(dearevenuecrs, :X), targets(dearevenuecrs, :Y), P, rts = :CRS)) ≈ ones(5)

    # Revenue VRS
    dearevenuevrs = dearevenue(X, Y, P, rts = :VRS)
    @test efficiency(dearevenuevrs) ≈ [0.6444444444 ;
                                    1.000;
                                    1.000;
                                    0.500;
                                    0.4565217391]
    @test efficiency(dearevenuevrs, :Technical) ≈ [0.7777777778 ;
                                1.000;
                                1.000;
                                0.500 ;
                                0.6000000000]
    @test efficiency(dearevenuevrs, :Allocative) ≈ [0.8285714286;
                                1.000;
                                1.000;
                                1.000;
                                0.7608695652]

    @test efficiency(dearevenue(targets(dearevenuevrs, :X), targets(dearevenuevrs, :Y), P, rts = :VRS)) ≈ ones(5)

    # Revnue FDH
    dearevenuefdh = dearevenue(X, Y, P, rts = :FDH)
    @test efficiency(dearevenuefdh) ≈ [0.6590909091;
                                1.0000000000;
                                1.0000000000;
                                0.5000000000;
                                0.4565217391]
    @test efficiency(dearevenuefdh, :Technical) ≈ [0.875;
                                1.000;
                                1.000;
                                0.500;
                                0.600]  atol = 1e-3
    @test efficiency(dearevenuefdh, :Allocative) ≈ [0.7532467532;
                                1.0000000000;
                                1.0000000000;
                                1.0000000000;
                                0.7608695652]

    @test efficiency(dearevenue(targets(dearevenuefdh, :X), targets(dearevenuefdh, :Y), P, rts = :FDH)) ≈ ones(5)

    # Check defaults
    @test efficiency(dearevenue(X, Y, P)) == efficiency(dearevenuevrs)
    @test efficiency(dearevenuevrs, :Economic) == efficiency(dearevenuevrs)

    # Print
    show(IOBuffer(), dearevenuecooper)

    # Test errors
    @test_throws DimensionMismatch dearevenue([1; 2 ; 3], [4 ; 5], [1; 1; 1]) #  Different number of observations
    @test_throws DimensionMismatch dearevenue([1; 2; 3], [4; 5; 6], [1; 2; 3; 4]) # Different number of observation in prices
    @test_throws DimensionMismatch dearevenue([1; 2; 3], [4 4; 5 5; 6 6], [4 4 4; 5 5 5; 6 6 6]) # Different number of output prices and outputs
    @test_throws ArgumentError dearevenue([1; 2; 3], [4; 5; 6], [1; 2; 3], rts = :Error) # Invalid returns to scale
    @test_throws ArgumentError dearevenue([1; 2; 3], [4; 5; 6], [1; 2; 3], dispos = :Error) # Invalid disposability
    @test_throws ArgumentError normfactor(dearevenue(X, Y, P)) # ERROR: RevenueDEAModel has no normalization factor

    # ------------------
    # Weak Disposability Tests
    # ------------------

    X = [1; 2; 3; 2; 4]
    Y = [2; 3; 4; 1; 3]
    P = [1; 1; 1; 1; 1]

    dearevenueStrong = dearevenue(X, Y, P, dispos = :Strong)
    @test efficiency(dearevenueStrong, :Economic) ≈ [1.0; 1.0; 1.0; 0.3333333333333333; 0.75]
    @test efficiency(dearevenueStrong, :Technical) ≈ [1.0; 1.0; 1.0; 0.3333333333333333; 0.75]
    @test efficiency(dearevenueStrong, :Allocative) ≈ [1.0; 1.0; 1.0; 1.0; 1.0]

    dearevenueWeak = dearevenue(X, Y, P, dispos = :Weak)
    @test efficiency(dearevenueWeak, :Economic) ≈ [1.0; 1.0; 1.0; 0.3333333333333333; 1.0]
    @test efficiency(dearevenueWeak, :Technical) ≈ [1.0; 1.0; 1.0; 0.3333333333333333; 1.0]
    @test efficiency(dearevenueWeak, :Allocative) ≈ [1.0; 1.0; 1.0; 1.0; 1.0]

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]
    P = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(dearevenue(X, Y, P)) ≈ [1; 1; 1; 1; 1; 1; 1; 1]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]
    P = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]

    @test efficiency(dearevenue(X, Y, P)) ≈ [1; 0.8571428571; 0.8571428571; 0.5714285714; 0.4285714286; 0.7142857143; 0.7142857143; 0.4642857143]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]
    P = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(dearevenue(X, Y, P)) ≈ [1; 1; 1; 1; 0.4615384615; 0.7777777778; 1; 0.2816951993]

end
