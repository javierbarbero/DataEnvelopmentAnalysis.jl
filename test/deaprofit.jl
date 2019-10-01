# Tests for Profit DEA Models
@testset "ProfitDEAModel" begin

    ## Test Profit DEA Model
    # Test with Zofio, Pastor and Aparicio (2013) data
    X = [1 1; 1 1; 0.75 1.5; 0.5 2; 0.5 2; 2 2; 2.75 3.5; 1.375 1.75]
    Y = [1 11; 5 3; 5 5; 2 9; 4 5; 4 2; 3 3; 4.5 3.5]
    P = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1]
    W = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1]

    GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2))
    GxGydollar = repeat(GxGydollar, 1, 2)

    deaprofitdollar = deaprofit(X, Y, W, P, GxGydollar, GxGydollar)
    @test efficiency(deaprofitdollar, :Economic)   ≈ [2; 2; 0; 2; 2; 8; 12; 4] atol = 1e-3
    @test efficiency(deaprofitdollar, :Technical)  ≈ [0; 0; 0; 0; 0; 6; 12; 3] atol = 1e-3
    @test efficiency(deaprofitdollar, :Allocative) ≈ [2; 2; 0; 2; 2; 2; 0; 1] atol = 1e-3

    # Print
    show(IOBuffer(), deaprofitdollar)

    # Test errors
    @test_throws ErrorException deaprofit([1; 2 ; 3], [4 ; 5], [1; 1; 1], [4; 5], [1; 2 ; 3], [4 ; 5]) #  Different number of observations
    @test_throws ErrorException deaprofit([1; 2; 3], [4; 5; 6], [1; 2; 3; 4], [4; 5; 6], [1; 2; 3], [4; 5; 6]) # Different number of observation in input prices
    @test_throws ErrorException deaprofit([1; 2; 3], [4; 5; 6], [1; 2; 3], [4; 5; 6; 7], [1; 2; 3], [4; 5; 6]) # Different number of observation in output prices
    @test_throws ErrorException deaprofit([1 1; 2 2; 3 3], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3], [4; 5; 6], [1 1; 2 2; 3 3], [4; 5; 6]) # Different number of input prices and inputs
    @test_throws ErrorException deaprofit([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6], [1; 2; 3], [4 4; 5 5; 6 6]) # Different number of oputput prices and outputs
    @test_throws ErrorException deaprofit([1 1; 2 2; 3 3], [4; 5; 6], [1 1; 2 2; 3 3], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3], [4; 5; 6]) # Different size of inputs direction
    @test_throws ErrorException deaprofit([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6]) # Different size of inputs direction

end
