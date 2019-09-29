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

    # Print
    show(IOBuffer(), mprod)

end
