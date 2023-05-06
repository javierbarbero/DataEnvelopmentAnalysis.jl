# Tests for BootstrapRadial DEA Returns To Scale Test
@testset "ReturnsToScaleTest" begin

    X = [2, 4, 3, 5, 6]
    Y = [1, 2, 3, 4, 5]

    # Input oriented
    rtstest = deartstest(X, Y, orient = :Input, rng = StableRNG(1234567))

    @test rtstest.S ≈ 0.802947 atol = 1e-5
    @test rtstest.p ≈ 0.335 atol = 1e-2
    @test rtstest.h ≈ 0.309757 atol = 1e-5
    @test criticalvalue(rtstest, 0.05) ≈ 0.674615 atol = 1e-5

    # Output oriented
    rtstest = deartstest(X, Y, orient = :Output, rng = StableRNG(1234567))

    @test rtstest.S ≈ 0.813093 atol = 1e-5
    @test rtstest.p ≈ 0.235 atol = 1e-2
    @test rtstest.h ≈ 0.259230 atol = 1e-5
    @test criticalvalue(rtstest, 0.05) ≈ 0.724159 atol = 1e-5

    # Print
    show(IOBuffer(), rtstest)
    
end

