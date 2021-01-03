# Tests for Profitability DEA Models
@testset "ProfitabilityDEAModel" begin

    ## Test Profitability DEA Model with Zofío and Prieto (2006) data
    X = [5 3; 2 4; 4 2; 4 8; 7 9]
    Y = [7 4; 10 8; 8 10; 5 4; 3 6]
    W = [2 1; 2 1; 2 1; 2 1; 2 1.0]
    P = [3 2; 3 2; 3 2; 3 2; 3 2.0]

    # alpha = 0.5 CRS equals Input Oriented CRS
    deaprofbl = deaprofitability(X, Y, W, P)

    @test typeof(deaprofbl) == ProfitabilityDEAModel

    @test efficiency(deaprofbl) ≈ [0.388;
                                   1.000;
                                   0.765;
                                   0.250;
                                   0.159] atol = 1e-3
    @test efficiency(deaprofbl, :CRS) ≈ [0.636;
                                1.000;
                                1.000;
                                0.250;
                                0.261] atol = 1e-3
    @test efficiency(deaprofbl, :VRS) ≈ [0.682;
                                1.000;
                                1.000;
                                0.250;
                                0.360] atol = 1e-3
    @test efficiency(deaprofbl, :Scale) ≈ [0.933;
                                1.000;
                                1.000;
                                1.000;
                                0.725] atol = 1e-3
    @test efficiency(deaprofbl, :Allocative) ≈ [0.610;
                                1.000;
                                0.765;
                                1.000;
                                0.609] atol = 1e-3

    @test efficiency(deaprofitability(targets(deaprofbl, :X), targets(deaprofbl, :Y), W, P)) ≈ ones(5) atol=1e-7

    # Check defaults
    @test efficiency(deaprofitability(X, Y, W, P, alpha = 0.5)) == efficiency(deaprofbl)
    @test efficiency(deaprofbl, :Economic) == efficiency(deaprofbl)


    # Print
    show(IOBuffer(), deaprofbl)

    # Test errors
    @test_throws ErrorException deaprofitability([1; 2 ; 3], [4 ; 5], [1; 1; 1], [4; 5]) #  Different number of observations
    @test_throws ErrorException deaprofitability([1; 2; 3], [4; 5; 6], [1; 2; 3; 4], [4; 5; 6]) # Different number of observation in input prices
    @test_throws ErrorException deaprofitability([1; 2; 3], [4; 5; 6], [1; 2; 3], [4; 5; 6; 7]) # Different number of observation in output prices
    @test_throws ErrorException deaprofitability([1 1; 2 2; 3 3], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3], [4; 5; 6]) # Different number of input prices and inputs
    @test_throws ErrorException deaprofitability([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6]) # Different number of oputput prices and outputs
    @test_throws ErrorException efficiency(deaprofbl, :Error)
    @test_throws ErrorException normfactor(deaprofitability(X, Y, W, P)) # ERROR: ProfitabilityDEAModel has no normalization factor

    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]
    W = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]
    P = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deaprofitability(X, Y, W, P)) ≈ ( sum(Y .* P, dims = 2) ./ sum(X .* W, dims = 2) )  / 0.25 atol = 1e-5

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]
    W = [1; 1; 1; 1; 1; 1; 1; 1]
    P = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]

    @test efficiency(deaprofitability(X, Y, W, P)) ≈ ( sum(Y .* P, dims = 2) ./ sum(X .* W, dims = 2) )  / 14 atol = 1e-5

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]
    W = [1; 1; 1; 1; 1; 1; 1; 1]
    P = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deaprofitability(X, Y, W, P)) ≈ ( sum(Y .* P, dims = 2) ./ sum(X .* W, dims = 2) )  / 1.25 atol = 1e-5

end
