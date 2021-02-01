# Tests for Cost DEA Models
@testset "CostDEAModel" begin

    ## Test Cost DEA Model with Cooper et al. (2007)
    # Test agains book results
    X = [3 2; 1 3; 4 6]
    Y = [3; 5; 6]
    W = [4 2; 4 2; 4 2]

    deacostcooper = deacost(X, Y, W, rts = :CRS)

    @test typeof(deacostcooper) == CostDEAModel
    @test ismonetary(deacostcooper) == false

    @test efficiency(deacostcooper, :Economic)   ≈ [0.375; 1; 0.429] atol = 1e-3
    @test efficiency(deacostcooper, :Technical)  ≈ [0.9  ; 1; 0.6  ] atol = 1e-3
    @test efficiency(deacostcooper, :Allocative) ≈ [0.417; 1; 0.714] atol = 1e-3

    ## Test Cost DEA Model with Zofío and Prieto (2006) data.
    # Test agains results in R
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]
    W = [2 1; 2 1; 2 1; 2 1; 2 1.0]

    # Cost CRS
    deacostcrs = deacost(X, Y, W, rts = :CRS)
    @test efficiency(deacostcrs) ≈ [0.4307692308;
                                    1.000;
                                    1.000;
                                    0.250;
                                    0.2608695652]
    @test efficiency(deacostcrs, :Technical) ≈ [0.6364;
                                1.000;
                                1.000;
                                0.250;
                                0.2609]  atol = 1e-3
    @test efficiency(deacostcrs, :Allocative) ≈ [0.6769230769;
                                1.000;
                                1.000;
                                1.000;
                                1.000]

    @test efficiency(deacost(targets(deacostcrs, :X), targets(deacostcrs, :Y), W, rts = :CRS)) ≈ ones(5)

    # Cost VRS
    deacostvrs = deacost(X, Y, W, rts = :VRS)
    @test efficiency(deacostvrs) ≈ [0.6153846154;
                                    1.000;
                                    1.000;
                                    0.500;
                                    0.3478260870]
    @test efficiency(deacostvrs, :Technical) ≈ [0.750;
                                1.000;
                                1.000;
                                0.500 ;
                                0.375]  atol = 1e-3
    @test efficiency(deacostvrs, :Allocative) ≈ [0.8205128205;
                                1.000;
                                1.000;
                                1.000;
                                0.9275362319]

    @test efficiency(deacost(targets(deacostvrs, :X), targets(deacostvrs, :Y), W, rts = :VRS)) ≈ ones(5)

    # Check defaults
    @test efficiency(deacost(X, Y, W)) == efficiency(deacostvrs)
    @test efficiency(deacostvrs, :Economic) == efficiency(deacostvrs)

    # Print
    show(IOBuffer(), deacostcooper)

    # Test errors
    @test_throws ErrorException deacost([1; 2 ; 3], [4 ; 5], [1; 1; 1]) #  Different number of observations
    @test_throws ErrorException deacost([1; 2; 3], [4; 5; 6], [1; 2; 3], rts = :Error) # Invalid returns to scale
    @test_throws ErrorException deacost([1; 2; 3], [4; 5; 6], [1; 2; 3; 4]) # Different number of observation in prices
    @test_throws ErrorException deacost([1 1; 2 2; 3 3 ], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3]) # Different number of input prices and inputs
    @test_throws ErrorException deacost([1; 2; 3], [4; 5; 6], [1; 2; 3], dispos = :Error) # Invalid disposability
    @test_throws ErrorException normfactor(deacost(X, Y, W)) # CostDEAModel has no normalization factor

    # ------------------
    # Weak Disposability Tests
    # ------------------

    X = [1; 2; 3; 2; 4]
    Y = [2; 3; 4; 1; 3]
    W = [1; 1; 1; 1; 1]

    deacostStrong = deacost(X, Y, W, dispos = :Strong)
    @test efficiency(deacostStrong, :Economic) ≈ [1.0; 1.0; 1.0; 0.5; 0.5]
    @test efficiency(deacostStrong, :Technical) ≈ [1.0; 1.0; 1.0; 0.5; 0.5]
    @test efficiency(deacostStrong, :Allocative) ≈ [1.0; 1.0; 1.0; 1.0; 1.0]

    deacostWeak = deacost(X, Y, W, dispos = :Weak)
    @test efficiency(deacostWeak, :Economic) ≈ [1.0; 1.0; 1.0; 1.0; 0.5]
    @test efficiency(deacostWeak, :Technical) ≈ [1.0; 1.0; 1.0; 1.0; 0.5]
    @test efficiency(deacostWeak, :Allocative) ≈ [1.0; 1.0; 1.0; 1.0; 1.0]

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]
    W = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]

    @test efficiency(deacost(X, Y, W)) ≈ [1; 0.8; 0.8; 0.5714285714; 0.4; 0.5714285714; 0.5714285714; 0.4166666667]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]
    W = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deacost(X, Y, W)) ≈ [1; 1; 1; 1; 1; 1; 1; 1]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]
    W = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deacost(X, Y, W)) ≈ [1; 1; 1; 1; 0.5; 0.4761904762; 0.8571428571; 0.2843710157]

end
