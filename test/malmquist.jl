# Tests for Mamlmquist DEA Model
@testset "MalmquistDEAModel" begin

    ## Test Mamlmquist DEA Model with 1 input and 1 output
    X = Array{Float64,3}(undef, 5, 1, 2)
    X[:, :, 1] = [2; 3; 5; 4; 4];
    X[:, :, 2] = [1; 2; 4; 3; 4];

    Y = Array{Float64,3}(undef, 5, 1, 2)
    Y[:, :, 1] = [1; 4; 6; 3; 5];
    Y[:, :, 2] = [1; 4; 6; 3; 3];

    # Default Malmquist Productivity Index
    mprod = malmquist(X, Y)

    @test nobs(mprod) == 5
    @test ninputs(mprod) == 1
    @test noutputs(mprod) == 1
    @test nperiods(mprod) == 2
    @test prodchange(mprod) ≈ [2.0000000000;
                               1.5000000000;
                               1.2500000000;
                               1.3333333333;
                               0.6000000000]
    @test prodchange(mprod, :Prod) == prodchange(mprod)
    @test prodchange(mprod, :EC) ≈ [1.3333333333;
                                    1.0000000000;
                                    0.8333333333;
                                    0.8888888889;
                                    0.4000000000];
    @test prodchange(mprod, :TC) ≈ [1.5; 1.5; 1.5; 1.5; 1.5];

    # Default output oriented
    mprodoo = malmquist(X, Y, orient = :Output)

    @test prodchange(mprodoo) == prodchange(mprod)
    @test prodchange(mprodoo, :Prod) == prodchange(mprod, :Prod)
    @test prodchange(mprodoo, :EC) == prodchange(mprod, :EC)
    @test prodchange(mprodoo, :TC) == prodchange(mprod, :TC)

    # Test geomean is the geometric mean of TC
    mprodbase = malmquist(X, Y, refperiod = :Base)
    mprodcomparison = malmquist(X, Y, refperiod = :Comparison)

    @test prodchange(mprod, :TC) == sqrt.( prodchange(mprodbase, :TC) .* prodchange(mprodcomparison, :TC) )

    ## Test Mamlmquist DEA Model with 1 input and 1 output; and 3 years
    X = Array{Float64,3}(undef, 5, 1, 3)
    X[:, :, 1] = [2; 3; 5; 4; 4];
    X[:, :, 2] = [1; 2; 4; 3; 4];
    X[:, :, 3] = [0.5; 1.5; 3; 2; 4]

    Y = Array{Float64,3}(undef, 5, 1, 3)
    Y[:, :, 1] = [1; 4; 6; 3; 5];
    Y[:, :, 2] = [1; 4; 6; 3; 3];
    Y[:, :, 3] = [2; 4; 6; 3; 1];

    # Default Malmquist Productivity Index
    mprod3 = malmquist(X, Y)

    @test nobs(mprod3) == 5
    @test ninputs(mprod3) == 1
    @test noutputs(mprod3) == 1
    @test nperiods(mprod3) == 3
    @test prodchange(mprod3) ≈ [2.0000000000 4.0000000;
                               1.5000000000 1.3333333;
                               1.2500000000 1.3333333;
                               1.3333333333 1.5000000;
                               0.6000000000 0.3333333]
    @test prodchange(mprod3, :Prod) == prodchange(mprod3)
    @test prodchange(mprod3, :EC) ≈ [1.3333333333 2.0000000;
                                    1.0000000000 0.6666667;
                                    0.8333333333 0.6666667;
                                    0.8888888889 0.7500000;
                                    0.4000000000 0.1666667] atol = 1e-7;
    @test prodchange(mprod3, :TC) ≈ [1.5 2.0; 1.5 2.0; 1.5 2.0; 1.5 2.0; 1.5 2.0];

    # Print
    show(IOBuffer(), mprod)
    show(IOBuffer(), mprod3)

    # Test errors
    @test_throws ErrorException malmquist(X[1:4,:,:], X[1:5,:,:]) # Different number of observations in inputs and outputs
    @test_throws ErrorException malmquist(X[:,:,1:2], X[:,:,1:3]) # Different number of time periods in inputs and outputs
    @test_throws ErrorException prodchange(mprod, :Error)

end
