# Tests for Profit DEA Models

# Struct for testing DEA model without monetary option
struct wrongProfitDEAmodel <: AbstractProfitDEAModel
    n::Int64
    eff::Vector
end

@testset "ProfitDEAModel" begin

    ## Test Profit DEA Model
    # Test with Zofio, Pastor and Aparicio (2013) data
    X = [1 1; 1 1; 0.75 1.5; 0.5 2; 0.5 2; 2 2; 2.75 3.5; 1.375 1.75]
    Y = [1 11; 5 3; 5 5; 2 9; 4 5; 4 2; 3 3; 4.5 3.5]
    P = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1]
    W = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1]

    GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2))
    GxGydollar = repeat(GxGydollar, 1, 2)

    deaprofitdollar = deaprofit(X, Y, W, P, Gx = GxGydollar, Gy = GxGydollar)

    @test typeof(deaprofitdollar) == ProfitDEAModel

    @test efficiency(deaprofitdollar, :Economic)   ≈ [2; 2; 0; 2; 2; 8; 12; 4] atol = 1e-3
    @test efficiency(deaprofitdollar, :Technical)  ≈ [0; 0; 0; 0; 0; 6; 12; 3] atol = 1e-3
    @test efficiency(deaprofitdollar, :Allocative) ≈ [2; 2; 0; 2; 2; 2; 0; 1] atol = 1e-3

    @test efficiency(deaprofit(X, Y, W, P, Gx = :Monetary, Gy = :Monetary)) == efficiency(deaprofitdollar)

    @test efficiency(deaprofit(targets(deaprofitdollar, :X), targets(deaprofitdollar, :Y), W, P, Gx = :Monetary, Gy = :Monetary)) ≈ zeros(8)

    @test normfactor(deaprofitdollar) == ones(8);

    @test peersmatrix(deaprofitdollar) == deaprofitdollar.lambda

    # Check directions checking technical efficiency
    @test efficiency(deaprofit(X, Y, W, P, Gx = :Zeros, Gy = :Ones), :Technical) == efficiency(deaddf(X, Y, Gx = :Zeros, Gy = :Ones, rts = :VRS))
    @test efficiency(deaprofit(X, Y, W, P, Gx = :Ones, Gy = :Zeros), :Technical) == efficiency(deaddf(X, Y, Gx = :Ones, Gy = :Zeros, rts = :VRS))
    @test efficiency(deaprofit(X, Y, W, P, Gx = :Observed, Gy = :Observed), :Technical) == efficiency(deaddf(X, Y, Gx = :Observed, Gy = :Observed, rts = :VRS))
    @test efficiency(deaprofit(X, Y, W, P, Gx = :Mean, Gy = :Mean), :Technical) == efficiency(deaddf(X, Y, Gx = :Mean, Gy = :Mean, rts = :VRS))

    # Print
    show(IOBuffer(), deaprofitdollar)

    # Check normalization factor with different direction
    @test normfactor(deaprofit(X, Y, W, P, Gx = :Ones, Gy = :Ones)) == vec(sum(P .* ones(size(Y)), dims = 2) .+ sum(W .* ones(size(X)), dims = 2))

    # Check monetary option
    @test ismonetary(deaprofitdollar) == false

    deaprofitmonetary = deaprofit(X, Y, W, P, Gx = :Ones, Gy = :Ones, monetary = true)
    @test efficiency(deaprofitmonetary, :Economic)   ≈ [2; 2; 0; 2; 2; 8; 12; 4] atol = 1e-3
    @test efficiency(deaprofitmonetary, :Technical)  ≈ [0; 0; 0; 0; 0; 6; 12; 3] atol = 1e-3
    @test efficiency(deaprofitmonetary, :Allocative) ≈ [2; 2; 0; 2; 2; 2; 0; 1] atol = 1e-3
    @test ismonetary(deaprofitmonetary) == true

    # Test errors
    @test_throws ErrorException deaprofit([1; 2 ; 3], [4 ; 5], [1; 1; 1], [4; 5], Gx = [1; 2 ; 3], Gy = [4 ; 5]) #  Different number of observations
    @test_throws ErrorException deaprofit([1; 2; 3], [4; 5; 6], [1; 2; 3; 4], [4; 5; 6], Gx = [1; 2; 3], Gy = [4; 5; 6]) # Different number of observation in input prices
    @test_throws ErrorException deaprofit([1; 2; 3], [4; 5; 6], [1; 2; 3], [4; 5; 6; 7], Gx = [1; 2; 3], Gy = [4; 5; 6]) # Different number of observation in output prices
    @test_throws ErrorException deaprofit([1 1; 2 2; 3 3], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3], [4; 5; 6], Gx = [1 1; 2 2; 3 3], Gy = [4; 5; 6]) # Different number of input prices and inputs
    @test_throws ErrorException deaprofit([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6], Gx = [1; 2; 3], Gy = [4 4; 5 5; 6 6]) # Different number of oputput prices and outputs
    @test_throws ErrorException deaprofit([1 1; 2 2; 3 3], [4; 5; 6], [1 1; 2 2; 3 3], [4; 5; 6], Gx = [1 1 1; 2 2 2; 3 3 3], Gy = [4; 5; 6]) # Different size of inputs direction
    @test_throws ErrorException deaprofit([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4; 5 5; 6 6], Gx = [1; 2; 3], Gy = [4 4 4; 5 5 5; 6 6 6]) # Different size of outputs direction
    @test_throws ErrorException deaprofit([1; 2; 3], [1; 2; 3], [1; 1; 1], [1; 1; 1], Gx = :Error, Gy = :Ones) # Invalid inputs direction
    @test_throws ErrorException deaprofit([1; 2; 3], [1; 2; 3], [1; 1; 1], [1; 1; 1], Gx = :Ones, Gy = :Error) # Invalid outputs direction

    # Test struct defaults and errors
    @test_throws ErrorException ismonetary(wrongProfitDEAmodel(3, [1; 2; 3])) # Model does not have info on peers
    
    # ------------------
    # Test Vector and Matrix inputs and outputs
    # ------------------
    # Tests against results in R

    # Inputs is Matrix, Outputs is Vector
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]
    W = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]
    P = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deaprofit(X, Y, W, P, Gx = X, Gy = Y)) ≈ [0; 0.1666666667; 0.1666666667; 0.375; 0.5454545455; 0.375; 0.375; 0.5283018868]

    # Inputs is Vector, Output is Matrix
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]
    W = [1; 1; 1; 1; 1; 1; 1; 1]
    P = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]

    @test efficiency(deaprofit(X, Y, W, P, Gx = X, Gy = Y)) ≈ [0; 0.1538461538; 0.1538461538; 0.6666666667; 1.142857143; 0.3636363636; 0.3636363636; 1]

    # Inputs is Vector, Output is Vector
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]
    W = [1; 1; 1; 1; 1; 1; 1; 1]
    P = [1; 1; 1; 1; 1; 1; 1; 1]

    @test efficiency(deaprofit(X, Y, W, P, Gx = X, Gy = Y)) ≈ [0.6666666667; 0; 0.0625; 0.1904761905; 0.4444444444; 0.3809523810; 0.2608695652; 0.6849978751]

end
