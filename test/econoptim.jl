# Tests for economic optimization problems. Tet only basic functionality and errors as other test are in specific DEA models.
@testset "EconOptimProblems" begin

    # ------------------
    # Cost Minimization
    # ------------------
    X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6	8]
    Y = [1; 1; 1; 1; 1; 1; 1; 1]
    W = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]

    Xtarget, clambda = deamincost(X, Y, W)

    @test Xtarget == 2 * ones(8,2)
    @test all(clambda[:,1] .== 1)

    @test_throws DimensionMismatch deamincost([1; 2 ; 3], [4 ; 5], [1; 1; 1]) #  Different number of observations
    @test_throws DimensionMismatch deamincost([1; 2; 3], [4; 5; 6], [1; 2; 3; 4]) # Different number of observation in prices
    @test_throws DimensionMismatch deamincost([1 1; 2 2; 3 3 ], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3]) # Different number of input prices and inputs
    @test_throws ArgumentError deamincost([1; 2; 3], [4; 5; 6], [1; 2; 3], rts = :Error) # Invalid returns to scale
    @test_throws ArgumentError deamincost([1; 2; 3], [4; 5; 6], [1; 2; 3], dispos = :Error) # Invalid disposability

    # ------------------
    # Revenue Maximization
    # ------------------
    X = [1; 1; 1; 1; 1; 1; 1; 1]
    Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5]
    P = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]

    Ytarget, rlambda = deamaxrevenue(X, Y, P)

    @test Ytarget == 7 * ones(8,2)
    @test all(clambda[:,1] .== 1)

    @test_throws DimensionMismatch  deamaxrevenue([1; 2 ; 3], [4 ; 5], [1; 1; 1]) #  Different number of observations
    @test_throws DimensionMismatch  deamaxrevenue([1; 2; 3], [4; 5; 6], [1; 2; 3; 4]) # Different number of observation in prices
    @test_throws DimensionMismatch  deamaxrevenue([1; 2; 3], [4 4; 5 5; 6 6], [4 4 4; 5 5 5; 6 6 6]) # Different number of output prices and outputs
    @test_throws ArgumentError deamaxrevenue([1; 2; 3], [4; 5; 6], [1; 2; 3], rts = :Error) # Invalid returns to scale
    @test_throws ArgumentError deamaxrevenue([1; 2; 3], [4; 5; 6], [1; 2; 3], dispos = :Error) # Invalid disposability

    # ------------------
    # Profit Maximization
    # ------------------
    X = [2; 4; 8; 12; 6; 14; 14; 9.412]
    Y = [1; 5; 8; 9; 3; 7; 9; 2.353]
    W = [1; 1; 1; 1; 1; 1; 1; 1]
    P = [2; 2; 2; 2; 2; 2; 2; 2]

    Xtarget, Ytarget, plambda = deamaxprofit(X, Y, W, P)

    @test Xtarget == 8 * ones(8,1)
    @test Ytarget == 8 * ones(8,1)
    @test all(plambda[:,3] .== 1)

    @test_throws DimensionMismatch deamaxprofit([1; 2 ; 3], [4 ; 5], [1; 1; 1], [4; 5]) #  Different number of observations
    @test_throws DimensionMismatch deamaxprofit([1; 2; 3], [4; 5; 6], [1; 2; 3; 4], [4; 5; 6]) # Different number of observation in input prices
    @test_throws DimensionMismatch deamaxprofit([1; 2; 3], [4; 5; 6], [1; 2; 3], [4; 5; 6; 7]) # Different number of observation in output prices
    @test_throws DimensionMismatch deamaxprofit([1 1; 2 2; 3 3], [4; 5; 6], [1 1 1; 2 2 2; 3 3 3], [4; 5; 6]) # Different number of input prices and inputs
    @test_throws DimensionMismatch deamaxprofit([1; 2; 3], [4 4; 5 5; 6 6], [1; 2; 3], [4 4 4; 5 5 5; 6 6 6]) # Different number of oputput prices and outputs

end
